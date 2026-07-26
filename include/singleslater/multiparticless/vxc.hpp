/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *
 */
#pragma once

#include <singleslater/multiparticless.hpp>
#include <singleslater/kohnsham.hpp>

#include <dft.hpp>
#include <dft/util.hpp>
#include <dft/epc.hpp>
#include <grid/integrator.hpp>
#include <util/threads.hpp>
#include <util/mpi.hpp>
#include <util/timer.hpp>
#include <cqlinalg.hpp>

#include <chrono>
#include <type_traits>

namespace ChronusQ {

  /**
   *  \brief Unified Kohn-Sham DFT driver for MultiParticleSS.
   *
   *  Forms the XC contributions for target subsystems (default = all subsystems)
   *  in a SINGLE traversal of the shared molecular grid. 
   *
   *  Energy accounting: each subsystem's intra-XC energy is written back to its
   *  KohnSham::XCEnergy (picked up by KohnSham::computeEnergy); the total inter-
   *  particle correlation energies are cached per functional pair and added
   *  once by MultiParticleSS::computeEnergy.
   *
  */
  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::formXC(EMPerturbation& pert,
      const XCTerms& xcTerms) {

    if( xcTerms.subsystems.empty() ) return;

    if( intParam.useGauXC )
      formXCGauXC(pert, xcTerms);
    else
      formXCInHouse(pert, xcTerms);

    needsFullXCEnergyUpdate = false;

  }

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::formXCInHouse(EMPerturbation& pert,
      const XCTerms& xcTerms) {

    const auto& parts = xcTerms.subsystems;
    const auto& pairs = xcTerms.interTerms;
    using EPCPair = typename XCTerms::InterTerm;

    ProgramTimer::tick("Form VXC");

    const size_t nP = parts.size();

    // ------------------------------------------------------------------------
    //  Grid / parallelism setup
    // ------------------------------------------------------------------------
    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(this->comm);
    size_t mpiSize   = MPISize(this->comm);
    size_t nAtoms    = this->molecule().nAtoms;

    int color = ((mpiSize < nAtoms) or (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL )
#endif
    {

      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
      size_t NTNPPB = nthreads * NPtsMaxPerBatch;

      SetLAThreads(1);

      // ----------------------------------------------------------------------
      //  Per-subsystem static info (basis, intra-XC functionals, grad needs)
      // ----------------------------------------------------------------------
      std::vector<SubSSPtr>   ss(nP);
      std::vector<BasisSet*>  bases(nP);
      std::vector<size_t>     NB(nP);
      std::vector<bool>       hasZ(nP);
      std::vector<std::shared_ptr<KohnSham<MatsT,IntsT>>> ks(nP);  // intra-XC owner (or null)
      std::vector<bool>       intraGGA(nP, false);
      std::vector<bool>       needGrad(nP, false);
      size_t NBmax = 0;

      for(size_t p = 0; p < nP; p++) {
        ss[p]    = subsystems.at(parts[p].label);
        bases[p] = &ss[p]->basisSet();
        NB[p]    = bases[p]->nBasis;
        hasZ[p]  = ss[p]->onePDM->hasZ();
        NBmax    = std::max(NBmax, NB[p]);
        if( ss[p]->onePDM->hasXY() )
          CErr("Relativistic (2C/4C) inter-particle EPC NYI for MultiParticleSS!");

        auto k = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss[p]);
        if( parts[p].formIntraXC ) {
          ks[p] = k;
          intraGGA[p] = std::any_of(k->functionals.begin(),k->functionals.end(),
                          [](const std::shared_ptr<DFTFunctional>& x){ return x->isGGA(); });
          if( intraGGA[p] ) needGrad[p] = true;
        }
      }

      // Inter-particle functional properties are derived from the registered
      //   model rather than duplicated in the active XC terms.
      std::vector<bool> interGGA(interFunctionals.size(), false);
      bool multiFunc = false;
      for(const auto& pair : pairs) {
        const auto& functionals = interFunctionals[pair.functionalIndex].functionals;
        interGGA[pair.functionalIndex] = std::any_of(functionals.begin(), functionals.end(),
          [](const std::shared_ptr<DFTFunctional>& functional) {
            return functional->isGGA();
          });
        if( interGGA[pair.functionalIndex] ) {
          needGrad[pair.firstIndex] = true;
          needGrad[pair.secondIndex] = true;
        }
        if( functionals.size() > 1 ) multiFunc = true;
      }

      bool anyGGA = false;
      for(size_t p = 0; p < nP; p++) anyGGA = anyGGA || needGrad[p];
      bool anyIntraGGA = false;
      for(size_t p = 0; p < nP; p++) anyIntraGGA = anyIntraGGA || intraGGA[p];

      for(size_t p = 0; p < nP; p++) if( ks[p] && ks[p]->functionals.size() > 1 ) multiFunc = true;

      // ----------------------------------------------------------------------
      //  Per-subsystem VXC accumulators (output + per-thread scratch)
      // ----------------------------------------------------------------------
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<IntsT>>> VXC(nP);
      std::vector<std::vector<IntsT*>> VXC_SZYX(nP);
      std::vector<std::vector<std::vector<IntsT*>>> integrateVXC(nP);
      std::vector<IntsT*> intVXC_RAW(nP, nullptr);

      for(size_t p = 0; p < nP; p++) {
        if( !parts[p].formVXC ) continue;
        VXC[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(
          NB[p], ss[p]->onePDM->hasXY(), ss[p]->onePDM->hasZ());
        VXC[p]->clear();
        VXC_SZYX[p] = VXC[p]->SZYXPointers();

        size_t NB2 = NB[p]*NB[p];
        intVXC_RAW[p] = (nthreads == 1) ? nullptr :
          CQMemManager::get().malloc<IntsT>(VXC_SZYX[p].size() * nthreads * NB2);
        for(size_t k = 0; k < VXC_SZYX[p].size(); k++) {
          integrateVXC[p].emplace_back();
          if( nthreads != 1 )
            for(size_t ith = 0; ith < nthreads; ith++)
              integrateVXC[p].back().emplace_back(intVXC_RAW[p] + (ith + k*nthreads) * NB2);
          else
            integrateVXC[p].back().emplace_back(VXC_SZYX[p][k]);
        }
        for(auto &X : integrateVXC[p]) for(auto &Y : X) std::fill_n(Y, NB2, IntsT(0.));
      }

      // Per-thread energy accumulators
      std::vector<std::vector<double>> integrateInterEnergy(interFunctionals.size(),
        std::vector<double>(nthreads, 0.));                                    // inter, per pair
      std::vector<std::vector<double>> integrateXCEnergy(nP,
        std::vector<double>(nthreads, 0.));                                    // intra, per subsystem

      // Per-thread wall-time accumulators (core-seconds) for a coarse
      //   density-eval / intra-XC / inter-XC breakdown of the single grid pass.
      std::vector<double> tDen(nthreads, 0.), tIntra(nthreads, 0.), tInter(nthreads, 0.),
                          tAsm(nthreads, 0.);
      using dclock = std::chrono::high_resolution_clock;

      // ----------------------------------------------------------------------
      //  Per-subsystem density / U-variable scratch (evaluated once per batch)
      // ----------------------------------------------------------------------
      std::vector<double*> DenS(nP, nullptr), DenZ(nP, nullptr);
      std::vector<double*> GDenS(nP, nullptr), GDenZ(nP, nullptr);
      std::vector<double*> U_n(nP, nullptr), U_gamma(nP, nullptr);
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>> Re(nP);

      for(size_t p = 0; p < nP; p++) {
        DenS[p] = CQMemManager::get().malloc<double>(NTNPPB);
        if( hasZ[p] ) DenZ[p] = CQMemManager::get().malloc<double>(NTNPPB);
        if( needGrad[p] ) {
          GDenS[p] = CQMemManager::get().malloc<double>(3*NTNPPB);
          if( hasZ[p] ) GDenZ[p] = CQMemManager::get().malloc<double>(3*NTNPPB);
        }
        U_n[p] = CQMemManager::get().malloc<double>(2*NTNPPB);
        if( needGrad[p] ) U_gamma[p] = CQMemManager::get().malloc<double>(3*NTNPPB);
        Re[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(ss[p]->onePDM->real_part());
      }

      // ----------------------------------------------------------------------
      //  Shared kernel / Z-vector scratch (per-thread, reused across
      //  subsystems and pairs)
      // ----------------------------------------------------------------------
      size_t NTNBmax2 = nthreads * NBmax * NBmax;
      double *SCRATCHNBNB = CQMemManager::get().malloc<double>(NTNBmax2);
      double *SCRATCHNBNP = CQMemManager::get().malloc<double>(NTNPPB*NBmax);

      double *epsEval   = CQMemManager::get().malloc<double>(NTNPPB);
      double *epcEval   = CQMemManager::get().malloc<double>(NTNPPB);
      double *dVU_n     = CQMemManager::get().malloc<double>(2*NTNPPB);
      double *dVU_gamma = anyGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      double *ZrhoVar1  = CQMemManager::get().malloc<double>(NTNPPB);
      double *ZgammaVar1= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *ZgammaVar2= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *ZgammaVar3= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *cross_U_gamma   = anyGGA ? CQMemManager::get().malloc<double>(4*NTNPPB) : nullptr;
      double *cross_dVU_gamma = anyGGA ? CQMemManager::get().malloc<double>(4*NTNPPB) : nullptr;
      // Intra GGA scratch used by formZ_vxc (Kx/Ky/Kz, Hx/Hy/Hz)
      double *KScratch = anyIntraGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      double *HScratch = anyIntraGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      IntsT  *ZMAT = CQMemManager::get().malloc<IntsT>(NTNPPB*NBmax);

      // Per-subsystem accumulated Z (one buffer per spin component, per thread).
      //   Every intra-XC and EPC contribution for a subsystem adds its Z here
      //   Then afterwards we perform a single syr2k operation per subsystem
      std::vector<std::vector<IntsT*>> ZACC(nP);
      std::vector<IntsT*> ZACC_RAW(nP, nullptr);
      for(size_t p = 0; p < nP; p++) {
        if( !parts[p].formVXC ) continue;
        size_t nComp = VXC_SZYX[p].size();
        size_t stride = nthreads * NB[p] * NPtsMaxPerBatch;
        ZACC_RAW[p] = CQMemManager::get().malloc<IntsT>(nComp * stride);
        for(size_t k = 0; k < nComp; k++)
          ZACC[p].emplace_back(ZACC_RAW[p] + k * stride);
      }

      double *epsSCR = nullptr, *dVU_n_SCR = nullptr, *dVU_gamma_SCR = nullptr;
      if( multiFunc ) {
        epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        dVU_n_SCR = CQMemManager::get().malloc<double>(2*NTNPPB);
        if( anyGGA ) dVU_gamma_SCR = CQMemManager::get().malloc<double>(3*NTNPPB);
      }

      auto epsScreenOf = [&]() {
        double e = intParam.epsilon / nAtoms / intParam.nAng / intParam.nRad;
        return std::max(e, std::numeric_limits<double>::epsilon());
      };

      // ----------------------------------------------------------------------
      //  Intra-particle XC assembly for subsystem p (currently only applies to electronic kohn-sham)
      // ----------------------------------------------------------------------
      auto intraVXC = [&](size_t p, size_t NPts, size_t thread_id,
                          const std::vector<size_t>& NBE_vec,
                          const std::vector<IntsT*>& BasisEval_vec,
                          std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec,
                          std::vector<double>& weights) {

        if( !ks[p] ) return;

        bool gga = intraGGA[p];
        const auto& func = ks[p]->functionals;

        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;
        auto   onePDM  = ss[p]->onePDM;

        size_t NBE_main  = NBE_vec[p];
        IntsT* BasisEval = BasisEval_vec[p];
        size_t IOff = NBE_main*NPts;

        double* DenS_loc  = DenS[p] + TIDNPPB;
        double* DenZ_loc  = DenZ[p]  ? DenZ[p]  + TIDNPPB : nullptr;
        double* GDenS_loc = GDenS[p] ? GDenS[p] + 3*TIDNPPB : nullptr;
        double* GDenZ_loc = GDenZ[p] ? GDenZ[p] + 3*TIDNPPB : nullptr;
        double* U_n_loc   = U_n[p] + 2*TIDNPPB;
        double* U_gamma_loc = U_gamma[p] ? U_gamma[p] + 3*TIDNPPB : nullptr;

        double* epsEval_loc   = epsEval + TIDNPPB;
        double* dVU_n_loc     = dVU_n + 2*TIDNPPB;
        double* dVU_gamma_loc = dVU_gamma ? dVU_gamma + 3*TIDNPPB : nullptr;
        double* ZrhoVar1_loc  = ZrhoVar1 + TIDNPPB;
        double* ZgammaVar1_loc= ZgammaVar1 ? ZgammaVar1 + TIDNPPB : nullptr;
        double* ZgammaVar2_loc= ZgammaVar2 ? ZgammaVar2 + TIDNPPB : nullptr;
        double* epsSCR_loc        = epsSCR        ? epsSCR + TIDNPPB : nullptr;
        double* dVU_n_SCR_loc     = dVU_n_SCR     ? dVU_n_SCR + 2*TIDNPPB : nullptr;
        double* dVU_gamma_SCR_loc = dVU_gamma_SCR ? dVU_gamma_SCR + 3*TIDNPPB : nullptr;
        double* KScratch_loc = KScratch ? KScratch + 3*TIDNPPB : nullptr;
        double* HScratch_loc = HScratch ? HScratch + 3*TIDNPPB : nullptr;
        IntsT*  ZMAT_loc = ZMAT + NBmax * TIDNPPB;

        double epsScreen = epsScreenOf();

        // Kernel derivatives wrt U variables (overwrites eps / dVU)
        loadVXCder(func, NPts, U_n_loc, U_gamma_loc, epsEval_loc, dVU_n_loc,
          dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc, dVU_gamma_SCR_loc);

        integrateXCEnergy[p][thread_id] += energy_vxc(NPts, weights, epsEval_loc, DenS_loc);

        if( !parts[p].formVXC ) return;

        IntsT* ZACC_S = ZACC[p][SCALAR] + thread_id * NB[p] * NPtsMaxPerBatch;

        // ---- SCALAR component ----------------------------------------------
        constructZVars(onePDM, SCALAR, gga, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);
        formZ_vxc(onePDM, SCALAR, gga, NPts, NBE_main, IOff, epsScreen, weights,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, DenZ_loc, nullptr, nullptr,
          GDenS_loc, GDenZ_loc, nullptr, nullptr,
          KScratch_loc, KScratch_loc ? KScratch_loc + NPts : nullptr,
          KScratch_loc ? KScratch_loc + 2*NPts : nullptr,
          HScratch_loc, HScratch_loc ? HScratch_loc + NPts : nullptr,
          HScratch_loc ? HScratch_loc + 2*NPts : nullptr, BasisEval, ZMAT_loc);

        // Accumulate into the shared per-subsystem Z (single syr2k done later)
        blas::axpy(NBE_main*NPts, IntsT(1.), ZMAT_loc, 1, ZACC_S, 1);

        if( not onePDM->hasZ() ) return;

        // ---- MZ component (UKS) --------------------------------------------
        IntsT* ZACC_Z = ZACC[p][MZ] + thread_id * NB[p] * NPtsMaxPerBatch;

        constructZVars(onePDM, MZ, gga, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);
        formZ_vxc(onePDM, MZ, gga, NPts, NBE_main, IOff, epsScreen, weights,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, DenZ_loc, nullptr, nullptr,
          GDenS_loc, GDenZ_loc, nullptr, nullptr,
          KScratch_loc, KScratch_loc ? KScratch_loc + NPts : nullptr,
          KScratch_loc ? KScratch_loc + 2*NPts : nullptr,
          HScratch_loc, HScratch_loc ? HScratch_loc + NPts : nullptr,
          HScratch_loc ? HScratch_loc + 2*NPts : nullptr, BasisEval, ZMAT_loc);

        blas::axpy(NBE_main*NPts, IntsT(1.), ZMAT_loc, 1, ZACC_Z, 1);
      }; // intraVXC

      // ----------------------------------------------------------------------
      //  Per-side EPC potential assembly for a given pair.
      //    mainIsFirst == true  => main = pr.firstIndex, aux = pr.secondIndex
      //    mainIsFirst == false => main = pr.secondIndex, aux = pr.firstIndex
      // ----------------------------------------------------------------------
      auto sideVXC = [&](const EPCPair& pr, bool mainIsFirst, size_t NPts, size_t thread_id,
                         const std::vector<size_t>& NBE_vec,
                         const std::vector<IntsT*>& BasisEval_vec,
                         std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec,
                         std::vector<double>& weights) {

        size_t m = mainIsFirst ? pr.firstIndex : pr.secondIndex;
        size_t x = mainIsFirst ? pr.secondIndex : pr.firstIndex;
        bool pairGGA = interGGA[pr.functionalIndex];
        const auto& epc_functionals = interFunctionals[pr.functionalIndex].functionals;

        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;

        auto  onePDM_main = ss[m]->onePDM;
        auto  onePDM_aux  = ss[x]->onePDM;
        bool  mainIsElectron = ss[m]->particle.charge < 0;

        size_t NBE_main  = NBE_vec[m];
        IntsT* BasisEval_main = BasisEval_vec[m];
        size_t IOff = NBE_main*NPts;

        double* DenS_main_loc  = DenS[m] + TIDNPPB;
        double* DenZ_main_loc  = DenZ[m] ? DenZ[m] + TIDNPPB : nullptr;
        double* DenS_aux_loc   = DenS[x] + TIDNPPB;
        double* DenZ_aux_loc   = DenZ[x] ? DenZ[x] + TIDNPPB : nullptr;
        double* GDenS_main_loc = GDenS[m] ? GDenS[m] + 3*TIDNPPB : nullptr;
        double* GDenZ_main_loc = GDenZ[m] ? GDenZ[m] + 3*TIDNPPB : nullptr;
        double* GDenS_aux_loc  = GDenS[x] ? GDenS[x] + 3*TIDNPPB : nullptr;
        double* GDenZ_aux_loc  = GDenZ[x] ? GDenZ[x] + 3*TIDNPPB : nullptr;
        double* U_n_main_loc   = U_n[m] + 2*TIDNPPB;
        double* U_n_aux_loc    = U_n[x] + 2*TIDNPPB;
        double* U_gamma_main_loc = U_gamma[m] ? U_gamma[m] + 3*TIDNPPB : nullptr;
        double* U_gamma_aux_loc  = U_gamma[x] ? U_gamma[x] + 3*TIDNPPB : nullptr;

        double* epsEval_loc   = epsEval + TIDNPPB;
        double* epcEval_loc   = epcEval + TIDNPPB;
        double* dVU_n_loc     = dVU_n + 2*TIDNPPB;
        double* dVU_gamma_loc = dVU_gamma ? dVU_gamma + 3*TIDNPPB : nullptr;
        double* ZrhoVar1_loc  = ZrhoVar1 + TIDNPPB;
        double* ZgammaVar1_loc= ZgammaVar1 ? ZgammaVar1 + TIDNPPB : nullptr;
        double* ZgammaVar2_loc= ZgammaVar2 ? ZgammaVar2 + TIDNPPB : nullptr;
        double* ZgammaVar3_loc= ZgammaVar3 ? ZgammaVar3 + TIDNPPB : nullptr;
        double* cross_U_gamma_loc   = cross_U_gamma   ? cross_U_gamma   + 4*TIDNPPB : nullptr;
        double* cross_dVU_gamma_loc = cross_dVU_gamma ? cross_dVU_gamma + 4*TIDNPPB : nullptr;
        double* epsSCR_loc        = epsSCR        ? epsSCR + TIDNPPB : nullptr;
        double* dVU_n_SCR_loc     = dVU_n_SCR     ? dVU_n_SCR + 2*TIDNPPB : nullptr;
        double* dVU_gamma_SCR_loc = dVU_gamma_SCR ? dVU_gamma_SCR + 3*TIDNPPB : nullptr;
        IntsT*  ZMAT_loc = ZMAT + NBmax * TIDNPPB;

        // Zero kernel-derivative accumulators (EPC kernel accumulates into them)
        std::fill_n(epsEval_loc, NPts, 0.);
        std::fill_n(epcEval_loc, NPts, 0.);
        std::fill_n(dVU_n_loc, 2*NPts, 0.);
        if( pairGGA ) std::fill_n(dVU_gamma_loc, 3*NPts, 0.);

        double epsScreen = epsScreenOf();

        // Cross V -> U variables for the GGA (EPC-19) kernel
        if( pairGGA )
          mkCrossAuxVar(false, mainIsElectron,
            onePDM_main, onePDM_aux, epsScreen, NPts,
            GDenS_main_loc, GDenS_main_loc + NPts, GDenS_main_loc + 2*NPts,
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            GDenZ_main_loc, GDenZ_main_loc ? GDenZ_main_loc + NPts : nullptr,
            GDenZ_main_loc ? GDenZ_main_loc + 2*NPts : nullptr,
            GDenS_aux_loc, GDenS_aux_loc + NPts, GDenS_aux_loc + 2*NPts,
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            GDenZ_aux_loc, GDenZ_aux_loc ? GDenZ_aux_loc + NPts : nullptr,
            GDenZ_aux_loc ? GDenZ_aux_loc + 2*NPts : nullptr,
            cross_U_gamma_loc);

        // EPC kernel derivatives wrt U variables (adds into epsEval / dVU)
        loadEPCVXCder(mainIsElectron,
          epc_functionals,
          NPts, U_n_main_loc, U_gamma_main_loc, U_n_aux_loc,
          U_gamma_aux_loc, cross_U_gamma_loc, epsEval_loc, dVU_n_loc,
          dVU_gamma_loc, cross_dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc,
          dVU_gamma_SCR_loc, cross_dVU_gamma_loc, epcEval_loc);

        // EPC energy: counted once, on the electron side. epsEval holds the
        //   pure EPC energy-per-particle (no intra-XC here), so
        //   integral(eps_epc * rho_e) is exactly this pair's EPC energy.
        if( mainIsElectron )
          integrateInterEnergy[pr.functionalIndex][thread_id] +=
            energy_vxc(NPts, weights, epsEval_loc, DenS_main_loc);

        if( !parts[m].formVXC ) return;

        // ---- SCALAR component ------------------------------------------------
        constructZVars(onePDM_main, SCALAR, pairGGA, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

        if( pairGGA ) {
          constructEPCZVars(mainIsElectron, SCALAR, NPts, cross_dVU_gamma_loc, ZgammaVar3_loc);
          formZ_vxc_epc(onePDM_main, onePDM_aux, SCALAR, pairGGA, NPts, NBE_main, IOff,
            epsScreen, weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_main_loc, DenZ_main_loc, nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc,
            nullptr, nullptr, DenS_aux_loc, DenZ_aux_loc, nullptr, nullptr,
            GDenS_aux_loc, GDenZ_aux_loc, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        } else {
          formZ_vxc(onePDM_main, SCALAR, pairGGA, NPts, NBE_main, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_main_loc, DenZ_main_loc,
            nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        }

        // Accumulate into the shared per-subsystem Z (single syr2k done later)
        blas::axpy(NBE_main*NPts, IntsT(1.), ZMAT_loc, 1,
          ZACC[m][SCALAR] + thread_id * NB[m] * NPtsMaxPerBatch, 1);

        if( not onePDM_main->hasZ() ) return;

        // ---- MZ component (UKS) ---------------------------------------------
        constructZVars(onePDM_main, MZ, pairGGA, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

        if( pairGGA ) {
          constructEPCZVars(mainIsElectron, MZ, NPts, cross_dVU_gamma_loc, ZgammaVar3_loc);
          formZ_vxc_epc(onePDM_main, onePDM_aux, MZ, pairGGA, NPts, NBE_main, IOff,
            epsScreen, weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_main_loc, DenZ_main_loc, nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc,
            nullptr, nullptr, DenS_aux_loc, DenZ_aux_loc, nullptr, nullptr,
            GDenS_aux_loc, GDenZ_aux_loc, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        } else {
          formZ_vxc(onePDM_main, MZ, pairGGA, NPts, NBE_main, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_main_loc, DenZ_main_loc,
            nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        }

        blas::axpy(NBE_main*NPts, IntsT(1.), ZMAT_loc, 1,
          ZACC[m][MZ] + thread_id * NB[m] * NPtsMaxPerBatch, 1);
      }; // sideVXC

      // Evaluate a pair side when its VXC is requested, or when it is the
      //   electronic side used to integrate the pair energy.
      auto needsPairSide = [&](size_t p) {
        return parts[p].formVXC or ss[p]->particle.charge < 0;
      };

      // ----------------------------------------------------------------------
      //  Batch callback: evaluate every active subsystem density once, then
      //  form the requested intra- and inter-particle DFT contributions.
      // ----------------------------------------------------------------------
      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch,
        std::vector<double> &weights, std::vector<size_t> NBE_vec,
        std::vector<IntsT*> BasisEval_vec,
        std::vector<std::vector<size_t>> & batchEvalShells_vec,
        std::vector<std::vector<std::pair<size_t,size_t>>> & subMatCut_vec) {

        size_t NPts = batch.size();
        size_t thread_id = GetThreadID();
        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;

        double epsScreen = epsScreenOf();

        double* SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NBmax*NBmax;
        double* SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NBmax*NPtsMaxPerBatch;

        // ---- Evaluate densities + U variables for each subsystem ----------
        auto _t0 = dclock::now();
        for(size_t p = 0; p < nP; p++) {
          SHELL_EVAL_TYPE denTyp = needGrad[p] ? GRADIENT : NOGRAD;
          size_t NBEp = NBE_vec[p];
          IntsT* BEp  = BasisEval_vec[p];
          auto&  subMatp = subMatCut_vec[p];
          double* gS = GDenS[p] ? GDenS[p] + 3*TIDNPPB : nullptr;

          evalDen(denTyp, NPts, NBEp, NB[p], subMatp, SCRATCHNBNB_loc, SCRATCHNBNP_loc,
            Re[p]->S().pointer(), DenS[p] + TIDNPPB,
            gS, gS ? gS + NPts : nullptr, gS ? gS + 2*NPts : nullptr, BEp);

          if( hasZ[p] ) {
            double* gZ = GDenZ[p] ? GDenZ[p] + 3*TIDNPPB : nullptr;
            evalDen(denTyp, NPts, NBEp, NB[p], subMatp, SCRATCHNBNB_loc, SCRATCHNBNP_loc,
              Re[p]->Z().pointer(), DenZ[p] + TIDNPPB,
              gZ, gZ ? gZ + NPts : nullptr, gZ ? gZ + 2*NPts : nullptr, BEp);
          }

          double* gSu = GDenS[p] ? GDenS[p] + 3*TIDNPPB : nullptr;
          double* gZu = GDenZ[p] ? GDenZ[p] + 3*TIDNPPB : nullptr;
          mkAuxVar(ss[p]->onePDM, needGrad[p], epsScreen, NPts,
            DenS[p] + TIDNPPB, DenZ[p] ? DenZ[p] + TIDNPPB : nullptr, nullptr, nullptr,
            gSu, gSu ? gSu + NPts : nullptr, gSu ? gSu + 2*NPts : nullptr,
            gZu, gZu ? gZu + NPts : nullptr, gZu ? gZu + 2*NPts : nullptr,
            nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, nullptr,
            nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, nullptr,nullptr,nullptr,
            U_n[p] + 2*TIDNPPB, U_gamma[p] ? U_gamma[p] + 3*TIDNPPB : nullptr);
        }

        // Zero the per-subsystem Z accumulators for this batch
        for(size_t p = 0; p < nP; p++)
          for(size_t k = 0; k < ZACC[p].size(); k++)
            std::fill_n(ZACC[p][k] + thread_id * NB[p] * NPtsMaxPerBatch,
                        NBE_vec[p]*NPts, IntsT(0.));

        auto _t1 = dclock::now();
        tDen[thread_id] += std::chrono::duration<double>(_t1 - _t0).count();

        // ---- Intra-particle XC (electron B3LYP, etc.) --------------------
        for(size_t p = 0; p < nP; p++)
          intraVXC(p, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);

        auto _t2 = dclock::now();
        tIntra[thread_id] += std::chrono::duration<double>(_t2 - _t1).count();

        // ---- Inter-particle correlation (EPC) for each required side ------
        for(auto& pr : pairs) {
          if( needsPairSide(pr.firstIndex) )
            sideVXC(pr, true, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);
          if( needsPairSide(pr.secondIndex) )
            sideVXC(pr, false, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);
        }

        auto _t3 = dclock::now();
        tInter[thread_id] += std::chrono::duration<double>(_t3 - _t2).count();

        // ---- Assemble VXC: ONE basis syr2k per subsystem/component -------
        for(size_t p = 0; p < nP; p++) {
          size_t NBEp = NBE_vec[p];
          IntsT* BEp  = BasisEval_vec[p];
          auto&  subMatp = subMatCut_vec[p];
          for(size_t k = 0; k < ZACC[p].size(); k++) {
            IntsT* ZACC_loc = ZACC[p][k] + thread_id * NB[p] * NPtsMaxPerBatch;
            blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,
              NBEp,NPts,IntsT(1.),BEp,NBEp,ZACC_loc,NBEp,IntsT(0.),
              SCRATCHNBNB_loc,NBEp);
            IncBySubMat(NB[p],NB[p],NBEp,NBEp,integrateVXC[p][k][thread_id],NB[p],
              SCRATCHNBNB_loc,NBEp,subMatp);
          }
        }

        auto _t4 = dclock::now();
        tAsm[thread_id] += std::chrono::duration<double>(_t4 - _t3).count();

      }; // vxcbuild

      // ----------------------------------------------------------------------
      //  Single N-basis integration over the shared molecular grid
      // ----------------------------------------------------------------------
      std::vector<SHELL_EVAL_TYPE> typsN(nP);
      for(size_t p = 0; p < nP; p++) typsN[p] = needGrad[p] ? GRADIENT : NOGRAD;

      BeckeIntegrator<EulerMac> integrator(intComm, this->molecule(), bases, typsN,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch, intParam.epsilon);

      size_t nGridPts = 0;
      integrator.integrateN<size_t>(nGridPts, vxcbuild);

      // ----------------------------------------------------------------------
      //  Finish VXC per subsystem: 4 pi factor, thread reduction, hermitize
      // ----------------------------------------------------------------------
      for(size_t p = 0; p < nP; p++) {
        size_t NB2 = NB[p]*NB[p];
        for(size_t k = 0; k < VXC_SZYX[p].size(); k++) {
          if( nthreads == 1 )
            blas::scal(NB2, IntsT(4*M_PI), VXC_SZYX[p][k], 1);
          else
            for(size_t ithread = 0; ithread < nthreads; ithread++)
              MatAdd('N','N',NB[p],NB[p],IntsT((ithread == 0) ? 0. : 1.),VXC_SZYX[p][k],NB[p],
                IntsT(4*M_PI),integrateVXC[p][k][ithread],NB[p], VXC_SZYX[p][k],NB[p]);
          HerMat('L',NB[p],VXC_SZYX[p][k],NB[p]);
        }
      }

      std::vector<double> interEnergy(interFunctionals.size(), 0.);
      for(const auto& pair : pairs)
        for(auto& e : integrateInterEnergy[pair.functionalIndex])
          interEnergy[pair.functionalIndex] += 4*M_PI*e;

      std::vector<double> intraEnergy(nP, 0.);
      for(size_t p = 0; p < nP; p++)
        for(auto& e : integrateXCEnergy[p]) intraEnergy[p] += 4*M_PI*e;

#ifdef CQ_ENABLE_MPI
      // Reduce VXC + energies across the integration communicator
      {
        IntsT* mpiScr = (mpiRank == 0) ? CQMemManager::get().malloc<IntsT>(NBmax*NBmax) : nullptr;
        for(size_t p = 0; p < nP; p++)
          for(auto &V : VXC_SZYX[p]) {
            MPIReduce(V, NB[p]*NB[p], mpiScr, 0, intComm);
            if( mpiRank == 0 ) std::copy_n(mpiScr, NB[p]*NB[p], V);
          }
        if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);
        for(const auto& pair : pairs)
          interEnergy[pair.functionalIndex] =
            MPIReduce(interEnergy[pair.functionalIndex], 0, intComm);
        for(size_t p = 0; p < nP; p++) intraEnergy[p] = MPIReduce(intraEnergy[p], 0, intComm);
      }
#endif

      // Add target VXC matrices and update the active energy caches.
      if( mpiRank == 0 ) {
        for(size_t p = 0; p < nP; p++) {
          if( parts[p].formVXC ) *ss[p]->fockMatrix += *VXC[p];
          if( parts[p].formIntraXC and ks[p] )
            ks[p]->XCEnergy = intraEnergy[p];
        }
        for(const auto& pair : pairs)
          interFunctionals[pair.functionalIndex].energy = interEnergy[pair.functionalIndex];
      }

      // ----------------------------------------------------------------------
      //  Free scratch
      // ----------------------------------------------------------------------
      CQMemManager::get().free(SCRATCHNBNB, SCRATCHNBNP, epsEval, epcEval, dVU_n,
        ZrhoVar1, ZMAT);
      for(size_t p = 0; p < nP; p++)
        if( ZACC_RAW[p] ) CQMemManager::get().free(ZACC_RAW[p]);
      if( anyGGA )
        CQMemManager::get().free(dVU_gamma, ZgammaVar1, ZgammaVar2, ZgammaVar3,
          cross_U_gamma, cross_dVU_gamma);
      if( anyIntraGGA )
        CQMemManager::get().free(KScratch, HScratch);
      if( multiFunc ) {
        CQMemManager::get().free(epsSCR, dVU_n_SCR);
        if( anyGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }
      for(size_t p = 0; p < nP; p++) {
        CQMemManager::get().free(DenS[p], U_n[p]);
        if( DenZ[p] )    CQMemManager::get().free(DenZ[p]);
        if( GDenS[p] )   CQMemManager::get().free(GDenS[p]);
        if( GDenZ[p] )   CQMemManager::get().free(GDenZ[p]);
        if( U_gamma[p] ) CQMemManager::get().free(U_gamma[p]);
        if( intVXC_RAW[p] ) CQMemManager::get().free(intVXC_RAW[p]);
      }

      SetLAThreads(LAThreads);
      MPICommFree(intComm);

    } // valid intComm

    MPI_Barrier(this->comm);

    ProgramTimer::tock("Form VXC");

  }; // MultiParticleSS<MatsT,IntsT>::formXCInHouse


  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::formXCGauXC(EMPerturbation&,
      const XCTerms& terms) {

    if constexpr (std::is_same_v<IntsT,dcomplex>) {
      CErr("Complex (dcomplex,dcomplex) GauXC MultiParticleSS XC NYI");
    } else {

      if(!this->gauxcUtils)
        CErr("GauXC MultiParticleSS requested without GauXCUtils");
      if(!this->gauxcUtils->integrator_pointer)
        CErr("GauXC MultiParticleSS requested without GauXC integrator");

      ProgramTimer::tick("Form VXC");

      const size_t nP = order_.size();
      auto globalIndexOf = [&](const std::string& label) {
        auto iter = std::find(order_.begin(), order_.end(), label);
        if( iter == order_.end() ) CErr("Unknown MultiParticleSS XC label " + label);
        return static_cast<size_t>(std::distance(order_.begin(), iter));
      };

      const auto& functionalSpec = this->gauxcUtils->multiparticle_functional_spec;
      if(functionalSpec.inter_functionals.size() != interFunctionals.size())
        CErr("CQ/GauXC inter-particle functional topology mismatch");
      for(size_t i = 0; i < interFunctionals.size(); ++i) {
        const auto& cqPair = interFunctionals[i];
        const auto& gxPair = functionalSpec.inter_functionals[i];
        if(gxPair.electron != globalIndexOf(cqPair.first) or
           gxPair.particle != globalIndexOf(cqPair.second))
          CErr("CQ/GauXC inter-particle functional ordering mismatch");
      }

      GauXC::MultiParticleXCTerms xcTerms;
      std::vector<bool> formVXC(nP, false);
      for(const auto& subsystem : terms.subsystems) {
        size_t p = globalIndexOf(subsystem.label);
        if(subsystem.formIntraXC) xcTerms.active_intra.push_back(p);
        if(subsystem.formVXC) {
          xcTerms.vxc_targets.push_back(p);
          formVXC[p] = true;
        }
      }
      for(const auto& pair : terms.interTerms) {
        xcTerms.active_inter.push_back(pair.functionalIndex);
      }

      std::vector<SubSSPtr> ss(nP);
      std::vector<Eigen::MatrixXd> Ps(nP);
      std::vector<Eigen::MatrixXd> Pz(nP);
      std::vector<GauXC::XCIntegrator<Eigen::MatrixXd>::multiparticle_density> densities(nP);

      // CQ stores a total scalar density for RKS. GauXC's RKS path expects the
      // spin density, so pass half the scalar density and double the VXC below.
      for(size_t p = 0; p < nP; ++p) {
        ss[p] = subsystems.at(order_[p]);
        if(ss[p]->onePDM->hasXY())
          CErr("Relativistic (2C/4C) GauXC MultiParticleSS XC NYI");

        const size_t NB = ss[p]->basisSet().nBasis;
        Ps[p] = Eigen::Map<Eigen::Matrix<double,-1,-1>>(
          ss[p]->onePDM->real_part().S().pointer(), NB, NB);

        if(ss[p]->onePDM->hasZ()) {
          Pz[p] = Eigen::Map<Eigen::Matrix<double,-1,-1>>(
            ss[p]->onePDM->real_part().Z().pointer(), NB, NB);
          densities[p] = {&Ps[p], &Pz[p]};
        } else {
          Ps[p] /= 2.0;
          densities[p] = {&Ps[p], nullptr};
        }
      }

      auto result = this->gauxcUtils->integrator_pointer->eval_exc_vxc(
        densities, functionalSpec, xcTerms);

      if(result.VXCs.size() != nP or result.VXCz.size() != nP)
        CErr("GauXC MultiParticleSS returned inconsistent VXC dimensions");

      // Only requested blocks are allocated and returned by GauXC.
      for(size_t p = 0; p < nP; ++p) {
        if( not formVXC[p] ) continue;
        Eigen::MatrixXd VXCs = result.VXCs[p];
        VXCs *= 2.0;
        ss[p]->fockMatrix->S() += VXCs;

        if(ss[p]->onePDM->hasZ()) {
          Eigen::MatrixXd VXCz = result.VXCz[p];
          VXCz *= 2.0;
          ss[p]->fockMatrix->Z() += VXCz;
        }
      }

      for(const auto& subsystem : terms.subsystems) {
        if( not subsystem.formIntraXC ) continue;
        size_t p = globalIndexOf(subsystem.label);
        auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss[p]);
        if( not ks or p >= result.intra_exc.size() )
          CErr("GauXC returned inconsistent intra-XC energy data");
        ks->XCEnergy = result.intra_exc[p];
      }
      if(result.inter_pair_exc.size() != interFunctionals.size())
        CErr("GauXC returned inconsistent inter-XC energy data");
      for(const auto& pair : terms.interTerms) {
        interFunctionals[pair.functionalIndex].energy =
          result.inter_pair_exc[pair.functionalIndex];
      }
      ProgramTimer::tock("Form VXC");

    }
  }


  // Complex GIAO in-house multiparticle XC driver.
  template <>
  inline void MultiParticleSS<dcomplex,dcomplex>::formXCInHouse(EMPerturbation& pert,
      const XCTerms& xcTerms) {

    const auto& parts = xcTerms.subsystems;
    const auto& pairs = xcTerms.interTerms;
    using EPCPair = XCTerms::InterTerm;

    ProgramTimer::tick("Form VXC");

    const size_t nP = parts.size();

    // ------------------------------------------------------------------------
    //  Grid / parallelism setup
    // ------------------------------------------------------------------------
    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(this->comm);
    size_t mpiSize   = MPISize(this->comm);
    size_t nAtoms    = this->molecule().nAtoms;

    int color = ((mpiSize < nAtoms) or (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL )
#endif
    {

      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
      size_t NTNPPB = nthreads * NPtsMaxPerBatch;
      size_t nEval = 2*nP;

      SetLAThreads(1);

      // ----------------------------------------------------------------------
      //  Per-subsystem static info (basis, intra-XC functionals, grad needs)
      // ----------------------------------------------------------------------
      std::vector<SubSSPtr>   ss(nP);
      std::vector<BasisSet*>  bases(nP);
      std::vector<size_t>     NB(nP);
      std::vector<bool>       hasZ(nP);
      std::vector<std::shared_ptr<KohnSham<dcomplex,dcomplex>>> ks(nP);
      std::vector<bool>       intraGGA(nP, false);
      std::vector<bool>       needGrad(nP, false);
      size_t NBmax = 0;

      for(size_t p = 0; p < nP; p++) {
        ss[p]    = subsystems.at(parts[p].label);
        bases[p] = &ss[p]->basisSet();
        NB[p]    = bases[p]->nBasis;
        hasZ[p]  = ss[p]->onePDM->hasZ();
        NBmax    = std::max(NBmax, NB[p]);
        if( ss[p]->onePDM->hasXY() )
          CErr("Relativistic (2C/4C) inter-particle EPC NYI for MultiParticleSS!");

        auto k = std::dynamic_pointer_cast<KohnSham<dcomplex,dcomplex>>(ss[p]);
        if( parts[p].formIntraXC ) {
          ks[p] = k;
          intraGGA[p] = std::any_of(k->functionals.begin(),k->functionals.end(),
                          [](const std::shared_ptr<DFTFunctional>& x){ return x->isGGA(); });
          if( intraGGA[p] ) needGrad[p] = true;
        }
      }

      // Inter-particle functional properties are derived from the registered
      //   model rather than duplicated in the active XC terms.
      std::vector<bool> interGGA(interFunctionals.size(), false);
      bool multiFunc = false;
      for(const auto& pair : pairs) {
        const auto& functionals = interFunctionals[pair.functionalIndex].functionals;
        interGGA[pair.functionalIndex] = std::any_of(functionals.begin(), functionals.end(),
          [](const std::shared_ptr<DFTFunctional>& functional) {
            return functional->isGGA();
          });
        if( interGGA[pair.functionalIndex] ) {
          needGrad[pair.firstIndex] = true;
          needGrad[pair.secondIndex] = true;
        }
        if( functionals.size() > 1 ) multiFunc = true;
      }

      bool anyGGA = false;
      for(size_t p = 0; p < nP; p++) anyGGA = anyGGA || needGrad[p];
      bool anyIntraGGA = false;
      for(size_t p = 0; p < nP; p++) anyIntraGGA = anyIntraGGA || intraGGA[p];

      for(size_t p = 0; p < nP; p++) if( ks[p] && ks[p]->functionals.size() > 1 ) multiFunc = true;


      // ----------------------------------------------------------------------
      //  Per-subsystem VXC accumulators (output + per-thread scratch)
      // ----------------------------------------------------------------------
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> VXC(nP);
      std::vector<std::vector<dcomplex*>> VXC_SZYX(nP);
      std::vector<std::vector<std::vector<dcomplex*>>> integrateVXC(nP);
      std::vector<dcomplex*> intVXC_RAW(nP, nullptr);

      for(size_t p = 0; p < nP; p++) {
        if( !parts[p].formVXC ) continue;
        VXC[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(
          NB[p], ss[p]->onePDM->hasXY(), ss[p]->onePDM->hasZ());
        VXC[p]->clear();
        VXC_SZYX[p] = VXC[p]->SZYXPointers();

        size_t NB2 = NB[p]*NB[p];
        intVXC_RAW[p] = (nthreads == 1) ? nullptr :
          CQMemManager::get().malloc<dcomplex>(VXC_SZYX[p].size() * nthreads * NB2);
        for(size_t k = 0; k < VXC_SZYX[p].size(); k++) {
          integrateVXC[p].emplace_back();
          if( nthreads != 1 )
            for(size_t ith = 0; ith < nthreads; ith++)
              integrateVXC[p].back().emplace_back(intVXC_RAW[p] + (ith + k*nthreads) * NB2);
          else
            integrateVXC[p].back().emplace_back(VXC_SZYX[p][k]);
        }
        for(auto &X : integrateVXC[p]) for(auto &Y : X) std::fill_n(Y, NB2, dcomplex(0.));
      }

      // Per-thread energy accumulators
      std::vector<std::vector<double>> integrateInterEnergy(interFunctionals.size(),
        std::vector<double>(nthreads, 0.));                                    // inter, per pair
      std::vector<std::vector<double>> integrateXCEnergy(nP,
        std::vector<double>(nthreads, 0.));                                    // intra, per subsystem

      // ----------------------------------------------------------------------
      //  Per-subsystem density / U-variable scratch (evaluated once per batch).
      //  GIAO keeps two evaluated copies of each subsystem basis:
      //    slot p     = main basis, old two-basis GIAO phase convention (+)
      //    slot p+nP  = auxiliary basis, old two-basis GIAO convention (-)
      // ----------------------------------------------------------------------
      std::vector<double*> DenS(nEval, nullptr), DenZ(nEval, nullptr);
      std::vector<double*> GDenS(nEval, nullptr), GDenZ(nEval, nullptr);
      std::vector<double*> U_n(nEval, nullptr), U_gamma(nEval, nullptr);
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> Re(nP);

      for(size_t p = 0; p < nP; p++) {
        Re[p] = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<dcomplex>>(ss[p]->onePDM);
        if( !Re[p] ) CErr("GIAO MultiParticleSS XC requested with non-complex density matrix");
      }
      for(size_t e = 0; e < nEval; e++) {
        size_t p = (e < nP) ? e : e - nP;
        DenS[e] = CQMemManager::get().malloc<double>(NTNPPB);
        if( hasZ[p] ) DenZ[e] = CQMemManager::get().malloc<double>(NTNPPB);
        if( needGrad[p] ) {
          GDenS[e] = CQMemManager::get().malloc<double>(3*NTNPPB);
          if( hasZ[p] ) GDenZ[e] = CQMemManager::get().malloc<double>(3*NTNPPB);
        }
        U_n[e] = CQMemManager::get().malloc<double>(2*NTNPPB);
        if( needGrad[p] ) U_gamma[e] = CQMemManager::get().malloc<double>(3*NTNPPB);
      }

      // ----------------------------------------------------------------------
      //  Shared kernel / Z-vector scratch (per-thread, reused across
      //  subsystems and pairs)
      // ----------------------------------------------------------------------
      size_t NTNBmax2 = nthreads * NBmax * NBmax;
      dcomplex *SCRATCHNBNB = CQMemManager::get().malloc<dcomplex>(NTNBmax2);
      dcomplex *SCRATCHNBNP = CQMemManager::get().malloc<dcomplex>(NTNPPB*NBmax);

      double *epsEval   = CQMemManager::get().malloc<double>(NTNPPB);
      double *epcEval   = CQMemManager::get().malloc<double>(NTNPPB);
      double *dVU_n     = CQMemManager::get().malloc<double>(2*NTNPPB);
      double *dVU_gamma = anyGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      double *ZrhoVar1  = CQMemManager::get().malloc<double>(NTNPPB);
      double *ZgammaVar1= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *ZgammaVar2= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *ZgammaVar3= anyGGA ? CQMemManager::get().malloc<double>(NTNPPB) : nullptr;
      double *cross_U_gamma   = anyGGA ? CQMemManager::get().malloc<double>(4*NTNPPB) : nullptr;
      double *cross_dVU_gamma = anyGGA ? CQMemManager::get().malloc<double>(4*NTNPPB) : nullptr;
      double *KScratch = anyIntraGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      double *HScratch = anyIntraGGA ? CQMemManager::get().malloc<double>(3*NTNPPB) : nullptr;
      dcomplex *ZMAT = CQMemManager::get().malloc<dcomplex>(NTNPPB*NBmax);

      // Per-subsystem accumulated Z (one buffer per spin component, per thread).
      //   Every intra-XC and EPC contribution for a subsystem adds its Z here.
      //   Then afterwards we perform a single GIAO Z -> VXC transform.
      std::vector<std::vector<dcomplex*>> ZACC(nP);
      std::vector<dcomplex*> ZACC_RAW(nP, nullptr);
      for(size_t p = 0; p < nP; p++) {
        if( !parts[p].formVXC ) continue;
        size_t nComp = VXC_SZYX[p].size();
        size_t stride = nthreads * NB[p] * NPtsMaxPerBatch;
        ZACC_RAW[p] = CQMemManager::get().malloc<dcomplex>(nComp * stride);
        for(size_t k = 0; k < nComp; k++)
          ZACC[p].emplace_back(ZACC_RAW[p] + k * stride);
      }

      double *epsSCR = nullptr, *dVU_n_SCR = nullptr, *dVU_gamma_SCR = nullptr;
      if( multiFunc ) {
        epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        dVU_n_SCR = CQMemManager::get().malloc<double>(2*NTNPPB);
        if( anyGGA ) dVU_gamma_SCR = CQMemManager::get().malloc<double>(3*NTNPPB);
      }

      auto epsScreenOf = [&]() {
        double e = intParam.epsilon / nAtoms / intParam.nAng / intParam.nRad;
        return std::max(e, std::numeric_limits<double>::epsilon());
      };

      auto assembleGIAOVXC = [&](size_t NBfull, size_t NBE, size_t NPts,
        dcomplex* BasisEval, dcomplex* ZBUF, dcomplex* SCRATCH,
        dcomplex* target, std::vector<std::pair<size_t,size_t>>& subMatCut) {

        // Creating according to J. Chem. Theory Comput. 2011, 7, 3097-3104
        //   Eq. 14, with the explicit conjugations needed for GIAO basis values.
        size_t NBNP = NBE * NPts;
        for(size_t i = 0; i < NBNP; i++) ZBUF[i] = std::conj(ZBUF[i]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,
          NBE,NBE,NPts,dcomplex(1.),ZBUF,NBE,BasisEval,NBE,
          dcomplex(0.),SCRATCH,NBE);

        for(size_t i = 0; i < NBNP; i++) ZBUF[i] = std::conj(ZBUF[i]);
        for(size_t i = 0; i < NBNP; i++) BasisEval[i] = std::conj(BasisEval[i]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,
          NBE,NBE,NPts,dcomplex(1.),BasisEval,NBE,ZBUF,NBE,
          dcomplex(1.),SCRATCH,NBE);

        for(size_t i = 0; i < NBNP; i++) BasisEval[i] = std::conj(BasisEval[i]);

        IncBySubMat(NBfull,NBfull,NBE,NBE,target,NBfull,SCRATCH,NBE,subMatCut);
      };

      // ----------------------------------------------------------------------
      //  Intra-particle XC assembly for subsystem p
      // ----------------------------------------------------------------------
      auto intraVXC = [&](size_t p, size_t NPts, size_t thread_id,
                          const std::vector<size_t>& NBE_vec,
                          const std::vector<dcomplex*>& BasisEval_vec,
                          std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec,
                          std::vector<double>& weights) {

        if( !ks[p] ) return;

        bool gga = intraGGA[p];
        const auto& func = ks[p]->functionals;

        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;
        auto   onePDM  = ss[p]->onePDM;

        size_t NBE_main  = NBE_vec[p];
        dcomplex* BasisEval = BasisEval_vec[p];
        size_t IOff = NBE_main*NPts;
        size_t mainSlot = p;

        double* DenS_loc  = DenS[mainSlot] + TIDNPPB;
        double* DenZ_loc  = DenZ[mainSlot]  ? DenZ[mainSlot]  + TIDNPPB : nullptr;
        double* GDenS_loc = GDenS[mainSlot] ? GDenS[mainSlot] + 3*TIDNPPB : nullptr;
        double* GDenZ_loc = GDenZ[mainSlot] ? GDenZ[mainSlot] + 3*TIDNPPB : nullptr;
        double* U_n_loc   = U_n[mainSlot] + 2*TIDNPPB;
        double* U_gamma_loc = U_gamma[mainSlot] ? U_gamma[mainSlot] + 3*TIDNPPB : nullptr;

        double* epsEval_loc   = epsEval + TIDNPPB;
        double* dVU_n_loc     = dVU_n + 2*TIDNPPB;
        double* dVU_gamma_loc = dVU_gamma ? dVU_gamma + 3*TIDNPPB : nullptr;
        double* ZrhoVar1_loc  = ZrhoVar1 + TIDNPPB;
        double* ZgammaVar1_loc= ZgammaVar1 ? ZgammaVar1 + TIDNPPB : nullptr;
        double* ZgammaVar2_loc= ZgammaVar2 ? ZgammaVar2 + TIDNPPB : nullptr;
        double* epsSCR_loc        = epsSCR        ? epsSCR + TIDNPPB : nullptr;
        double* dVU_n_SCR_loc     = dVU_n_SCR     ? dVU_n_SCR + 2*TIDNPPB : nullptr;
        double* dVU_gamma_SCR_loc = dVU_gamma_SCR ? dVU_gamma_SCR + 3*TIDNPPB : nullptr;
        double* KScratch_loc = KScratch ? KScratch + 3*TIDNPPB : nullptr;
        double* HScratch_loc = HScratch ? HScratch + 3*TIDNPPB : nullptr;
        dcomplex* ZMAT_loc = ZMAT + NBmax * TIDNPPB;

        double epsScreen = epsScreenOf();

        // Kernel derivatives wrt U variables (overwrites eps / dVU)
        loadVXCder(func, NPts, U_n_loc, U_gamma_loc, epsEval_loc, dVU_n_loc,
          dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc, dVU_gamma_SCR_loc);

        integrateXCEnergy[p][thread_id] += energy_vxc(NPts, weights, epsEval_loc, DenS_loc);

        if( !parts[p].formVXC ) return;

        dcomplex* ZACC_S = ZACC[p][SCALAR] + thread_id * NB[p] * NPtsMaxPerBatch;

        // ---- SCALAR component ----------------------------------------------
        constructZVars(onePDM, SCALAR, gga, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);
        formZ_vxc(onePDM, SCALAR, gga, NPts, NBE_main, IOff, epsScreen, weights,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, DenZ_loc, nullptr, nullptr,
          GDenS_loc, GDenZ_loc, nullptr, nullptr,
          KScratch_loc, KScratch_loc ? KScratch_loc + NPts : nullptr,
          KScratch_loc ? KScratch_loc + 2*NPts : nullptr,
          HScratch_loc, HScratch_loc ? HScratch_loc + NPts : nullptr,
          HScratch_loc ? HScratch_loc + 2*NPts : nullptr, BasisEval, ZMAT_loc);

        // Accumulate into the shared per-subsystem Z (single transform done later)
        blas::axpy(NBE_main*NPts, dcomplex(1.), ZMAT_loc, 1, ZACC_S, 1);

        if( not onePDM->hasZ() ) return;

        // ---- MZ component (UKS) --------------------------------------------
        dcomplex* ZACC_Z = ZACC[p][MZ] + thread_id * NB[p] * NPtsMaxPerBatch;

        constructZVars(onePDM, MZ, gga, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);
        formZ_vxc(onePDM, MZ, gga, NPts, NBE_main, IOff, epsScreen, weights,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, DenZ_loc, nullptr, nullptr,
          GDenS_loc, GDenZ_loc, nullptr, nullptr,
          KScratch_loc, KScratch_loc ? KScratch_loc + NPts : nullptr,
          KScratch_loc ? KScratch_loc + 2*NPts : nullptr,
          HScratch_loc, HScratch_loc ? HScratch_loc + NPts : nullptr,
          HScratch_loc ? HScratch_loc + 2*NPts : nullptr, BasisEval, ZMAT_loc);

        blas::axpy(NBE_main*NPts, dcomplex(1.), ZMAT_loc, 1, ZACC_Z, 1);
      }; // intraVXC

      // ----------------------------------------------------------------------
      //  Per-side EPC potential assembly for a given pair.
      //    mainIsFirst == true  => main = pr.firstIndex, aux = pr.secondIndex
      //    mainIsFirst == false => main = pr.secondIndex, aux = pr.firstIndex
      // ----------------------------------------------------------------------
      auto sideVXC = [&](const EPCPair& pr, bool mainIsFirst, size_t NPts, size_t thread_id,
                         const std::vector<size_t>& NBE_vec,
                         const std::vector<dcomplex*>& BasisEval_vec,
                         std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec,
                         std::vector<double>& weights) {

        size_t m = mainIsFirst ? pr.firstIndex : pr.secondIndex;
        size_t x = mainIsFirst ? pr.secondIndex : pr.firstIndex;
        bool pairGGA = interGGA[pr.functionalIndex];
        const auto& epc_functionals = interFunctionals[pr.functionalIndex].functionals;

        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;

        auto  onePDM_main = ss[m]->onePDM;
        auto  onePDM_aux  = ss[x]->onePDM;
        bool  mainIsElectron = ss[m]->particle.charge < 0;

        size_t NBE_main  = NBE_vec[m];
        dcomplex* BasisEval_main = BasisEval_vec[m];
        size_t IOff = NBE_main*NPts;
        size_t mainSlot = m;       // + phase copy
        size_t auxSlot  = x + nP;  // - phase copy, matching old aux basis path

        double* DenS_main_loc  = DenS[mainSlot] + TIDNPPB;
        double* DenZ_main_loc  = DenZ[mainSlot] ? DenZ[mainSlot] + TIDNPPB : nullptr;
        double* DenS_aux_loc   = DenS[auxSlot] + TIDNPPB;
        double* DenZ_aux_loc   = DenZ[auxSlot] ? DenZ[auxSlot] + TIDNPPB : nullptr;
        double* GDenS_main_loc = GDenS[mainSlot] ? GDenS[mainSlot] + 3*TIDNPPB : nullptr;
        double* GDenZ_main_loc = GDenZ[mainSlot] ? GDenZ[mainSlot] + 3*TIDNPPB : nullptr;
        double* GDenS_aux_loc  = GDenS[auxSlot] ? GDenS[auxSlot] + 3*TIDNPPB : nullptr;
        double* GDenZ_aux_loc  = GDenZ[auxSlot] ? GDenZ[auxSlot] + 3*TIDNPPB : nullptr;
        double* U_n_main_loc   = U_n[mainSlot] + 2*TIDNPPB;
        double* U_n_aux_loc    = U_n[auxSlot] + 2*TIDNPPB;
        double* U_gamma_main_loc = U_gamma[mainSlot] ? U_gamma[mainSlot] + 3*TIDNPPB : nullptr;
        double* U_gamma_aux_loc  = U_gamma[auxSlot] ? U_gamma[auxSlot] + 3*TIDNPPB : nullptr;

        double* epsEval_loc   = epsEval + TIDNPPB;
        double* epcEval_loc   = epcEval + TIDNPPB;
        double* dVU_n_loc     = dVU_n + 2*TIDNPPB;
        double* dVU_gamma_loc = dVU_gamma ? dVU_gamma + 3*TIDNPPB : nullptr;
        double* ZrhoVar1_loc  = ZrhoVar1 + TIDNPPB;
        double* ZgammaVar1_loc= ZgammaVar1 ? ZgammaVar1 + TIDNPPB : nullptr;
        double* ZgammaVar2_loc= ZgammaVar2 ? ZgammaVar2 + TIDNPPB : nullptr;
        double* ZgammaVar3_loc= ZgammaVar3 ? ZgammaVar3 + TIDNPPB : nullptr;
        double* cross_U_gamma_loc   = cross_U_gamma   ? cross_U_gamma   + 4*TIDNPPB : nullptr;
        double* cross_dVU_gamma_loc = cross_dVU_gamma ? cross_dVU_gamma + 4*TIDNPPB : nullptr;
        double* epsSCR_loc        = epsSCR        ? epsSCR + TIDNPPB : nullptr;
        double* dVU_n_SCR_loc     = dVU_n_SCR     ? dVU_n_SCR + 2*TIDNPPB : nullptr;
        double* dVU_gamma_SCR_loc = dVU_gamma_SCR ? dVU_gamma_SCR + 3*TIDNPPB : nullptr;
        dcomplex* ZMAT_loc = ZMAT + NBmax * TIDNPPB;

        // Zero kernel-derivative accumulators (EPC kernel accumulates into them)
        std::fill_n(epsEval_loc, NPts, 0.);
        std::fill_n(epcEval_loc, NPts, 0.);
        std::fill_n(dVU_n_loc, 2*NPts, 0.);
        if( pairGGA ) std::fill_n(dVU_gamma_loc, 3*NPts, 0.);

        double epsScreen = epsScreenOf();

        // Cross V -> U variables for the GGA (EPC-19) kernel
        if( pairGGA )
          mkCrossAuxVar(false, mainIsElectron,
            onePDM_main, onePDM_aux, epsScreen, NPts,
            GDenS_main_loc, GDenS_main_loc + NPts, GDenS_main_loc + 2*NPts,
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            GDenZ_main_loc, GDenZ_main_loc ? GDenZ_main_loc + NPts : nullptr,
            GDenZ_main_loc ? GDenZ_main_loc + 2*NPts : nullptr,
            GDenS_aux_loc, GDenS_aux_loc + NPts, GDenS_aux_loc + 2*NPts,
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            GDenZ_aux_loc, GDenZ_aux_loc ? GDenZ_aux_loc + NPts : nullptr,
            GDenZ_aux_loc ? GDenZ_aux_loc + 2*NPts : nullptr,
            cross_U_gamma_loc);

        // EPC kernel derivatives wrt U variables (adds into epsEval / dVU)
        loadEPCVXCder(mainIsElectron,
          epc_functionals,
          NPts, U_n_main_loc, U_gamma_main_loc, U_n_aux_loc,
          U_gamma_aux_loc, cross_U_gamma_loc, epsEval_loc, dVU_n_loc,
          dVU_gamma_loc, cross_dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc,
          dVU_gamma_SCR_loc, cross_dVU_gamma_loc, epcEval_loc);

        // EPC energy: counted once, on the electron side. epsEval holds the
        //   pure EPC energy-per-particle (no intra-XC here), so
        //   integral(eps_epc * rho_e) is exactly this pair's EPC energy.
        if( mainIsElectron )
          integrateInterEnergy[pr.functionalIndex][thread_id] +=
            energy_vxc(NPts, weights, epsEval_loc, DenS_main_loc);

        if( !parts[m].formVXC ) return;

        // ---- SCALAR component ----------------------------------------------
        constructZVars(onePDM_main, SCALAR, pairGGA, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

        if( pairGGA ) {
          constructEPCZVars(mainIsElectron, SCALAR, NPts, cross_dVU_gamma_loc, ZgammaVar3_loc);
          formZ_vxc_epc(onePDM_main, onePDM_aux, SCALAR, pairGGA, NPts, NBE_main, IOff,
            epsScreen, weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_main_loc, DenZ_main_loc, nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc,
            nullptr, nullptr, DenS_aux_loc, DenZ_aux_loc, nullptr, nullptr,
            GDenS_aux_loc, GDenZ_aux_loc, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        } else {
          formZ_vxc(onePDM_main, SCALAR, pairGGA, NPts, NBE_main, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_main_loc, DenZ_main_loc,
            nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        }

        // Accumulate into the shared per-subsystem Z (single transform done later)
        blas::axpy(NBE_main*NPts, dcomplex(1.), ZMAT_loc, 1,
          ZACC[m][SCALAR] + thread_id * NB[m] * NPtsMaxPerBatch, 1);

        if( not onePDM_main->hasZ() ) return;

        // ---- MZ component (UKS) ---------------------------------------------
        constructZVars(onePDM_main, MZ, pairGGA, NPts, dVU_n_loc, dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

        if( pairGGA ) {
          constructEPCZVars(mainIsElectron, MZ, NPts, cross_dVU_gamma_loc, ZgammaVar3_loc);
          formZ_vxc_epc(onePDM_main, onePDM_aux, MZ, pairGGA, NPts, NBE_main, IOff,
            epsScreen, weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_main_loc, DenZ_main_loc, nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc,
            nullptr, nullptr, DenS_aux_loc, DenZ_aux_loc, nullptr, nullptr,
            GDenS_aux_loc, GDenZ_aux_loc, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        } else {
          formZ_vxc(onePDM_main, MZ, pairGGA, NPts, NBE_main, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_main_loc, DenZ_main_loc,
            nullptr, nullptr, GDenS_main_loc, GDenZ_main_loc, nullptr, nullptr,
            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, BasisEval_main, ZMAT_loc);
        }

        blas::axpy(NBE_main*NPts, dcomplex(1.), ZMAT_loc, 1,
          ZACC[m][MZ] + thread_id * NB[m] * NPtsMaxPerBatch, 1);
      }; // sideVXC

      // Evaluate a pair side when its VXC is requested, or when it is the
      //   electronic side used to integrate the pair energy.
      auto needsPairSide = [&](size_t p) {
        return parts[p].formVXC or ss[p]->particle.charge < 0;
      };

      // ----------------------------------------------------------------------
      //  Batch callback: evaluate every active subsystem density once in each
      //  GIAO phase convention, then form the requested DFT terms.
      // ----------------------------------------------------------------------
      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch,
        std::vector<double> &weights, std::vector<size_t> NBE_vec,
        std::vector<dcomplex*> BasisEval_vec,
        std::vector<std::vector<size_t>> & batchEvalShells_vec,
        std::vector<std::vector<std::pair<size_t,size_t>>> & subMatCut_vec) {

        size_t NPts = batch.size();
        size_t thread_id = GetThreadID();
        size_t TIDNPPB = thread_id * NPtsMaxPerBatch;

        double epsScreen = epsScreenOf();

        dcomplex* SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NBmax*NBmax;
        dcomplex* SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NBmax*NPtsMaxPerBatch;

        // ---- Evaluate densities + U variables for each subsystem ----------
        for(size_t e = 0; e < nEval; e++) {
          size_t p = (e < nP) ? e : e - nP;
          SHELL_EVAL_TYPE denTyp = needGrad[p] ? GRADIENT : NOGRAD;
          size_t NBEp = NBE_vec[e];
          dcomplex* BEp  = BasisEval_vec[e];
          auto&  subMatp = subMatCut_vec[e];
          double* gS = GDenS[e] ? GDenS[e] + 3*TIDNPPB : nullptr;

          evalDen(denTyp, NPts, NBEp, NB[p], subMatp, SCRATCHNBNB_loc, SCRATCHNBNP_loc,
            Re[p]->S().pointer(), DenS[e] + TIDNPPB,
            gS, gS ? gS + NPts : nullptr, gS ? gS + 2*NPts : nullptr, BEp);

          if( hasZ[p] ) {
            double* gZ = GDenZ[e] ? GDenZ[e] + 3*TIDNPPB : nullptr;
            evalDen(denTyp, NPts, NBEp, NB[p], subMatp, SCRATCHNBNB_loc, SCRATCHNBNP_loc,
              Re[p]->Z().pointer(), DenZ[e] + TIDNPPB,
              gZ, gZ ? gZ + NPts : nullptr, gZ ? gZ + 2*NPts : nullptr, BEp);
          }

          double* gSu = GDenS[e] ? GDenS[e] + 3*TIDNPPB : nullptr;
          double* gZu = GDenZ[e] ? GDenZ[e] + 3*TIDNPPB : nullptr;
          mkAuxVar(ss[p]->onePDM, needGrad[p], epsScreen, NPts,
            DenS[e] + TIDNPPB, DenZ[e] ? DenZ[e] + TIDNPPB : nullptr, nullptr, nullptr,
            gSu, gSu ? gSu + NPts : nullptr, gSu ? gSu + 2*NPts : nullptr,
            gZu, gZu ? gZu + NPts : nullptr, gZu ? gZu + 2*NPts : nullptr,
            nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, nullptr,
            nullptr,nullptr,nullptr, nullptr,nullptr,nullptr, nullptr,nullptr,nullptr,
            U_n[e] + 2*TIDNPPB, U_gamma[e] ? U_gamma[e] + 3*TIDNPPB : nullptr);
        }

        // Zero the per-subsystem Z accumulators for this batch
        for(size_t p = 0; p < nP; p++)
          for(size_t k = 0; k < ZACC[p].size(); k++)
            std::fill_n(ZACC[p][k] + thread_id * NB[p] * NPtsMaxPerBatch,
                        NBE_vec[p]*NPts, dcomplex(0.));

        // ---- Intra-particle XC (electron B3LYP, etc.) --------------------
        for(size_t p = 0; p < nP; p++)
          intraVXC(p, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);

        // ---- Inter-particle correlation (EPC) for each required side ------
        for(auto& pr : pairs) {
          if( needsPairSide(pr.firstIndex) )
            sideVXC(pr, true, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);
          if( needsPairSide(pr.secondIndex) )
            sideVXC(pr, false, NPts, thread_id, NBE_vec, BasisEval_vec, subMatCut_vec, weights);
        }

        // ---- Assemble VXC: ONE GIAO Z -> VXC transform per subsystem/component
        for(size_t p = 0; p < nP; p++) {
          size_t NBEp = NBE_vec[p];
          dcomplex* BEp  = BasisEval_vec[p];
          auto&  subMatp = subMatCut_vec[p];
          for(size_t k = 0; k < ZACC[p].size(); k++) {
            dcomplex* ZACC_loc = ZACC[p][k] + thread_id * NB[p] * NPtsMaxPerBatch;
            assembleGIAOVXC(NB[p], NBEp, NPts, BEp, ZACC_loc, SCRATCHNBNB_loc,
              integrateVXC[p][k][thread_id], subMatp);
          }
        }

      }; // vxcbuild

      // ----------------------------------------------------------------------
      //  Single N-basis GIAO integration over the shared molecular grid.
      // ----------------------------------------------------------------------
      std::vector<BasisSet*> evalBases(nEval);
      std::vector<SHELL_EVAL_TYPE> typsN(nEval);
      std::vector<double> phaseScales(nEval);
      for(size_t p = 0; p < nP; p++) {
        typsN[p] = needGrad[p] ? GRADIENT : NOGRAD;
        typsN[p+nP] = typsN[p];
        evalBases[p] = bases[p];
        evalBases[p+nP] = bases[p];
        phaseScales[p] = 1.0;     // first basis
        phaseScales[p+nP] = -1.0; // auxiliary basis
      }

      BeckeIntegrator<EulerMac> integrator(intComm, this->molecule(), evalBases, typsN,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch, intParam.epsilon);

      size_t nGridPts = 0;
      integrator.integrateN<size_t>(nGridPts, vxcbuild, pert, phaseScales);

      // ----------------------------------------------------------------------
      //  Finish VXC per subsystem: 4 pi factor and thread reduction.
      //  HerMat is only for GTO; GIAO forms the full complex matrix explicitly
      //  in the assembly above.
      // ----------------------------------------------------------------------
      for(size_t p = 0; p < nP; p++) {
        size_t NB2 = NB[p]*NB[p];
        for(size_t k = 0; k < VXC_SZYX[p].size(); k++) {
          if( nthreads == 1 )
            blas::scal(NB2, dcomplex(4*M_PI), VXC_SZYX[p][k], 1);
          else
            for(size_t ithread = 0; ithread < nthreads; ithread++)
              MatAdd('N','N',NB[p],NB[p],dcomplex((ithread == 0) ? 0. : 1.),VXC_SZYX[p][k],NB[p],
                dcomplex(4*M_PI),integrateVXC[p][k][ithread],NB[p], VXC_SZYX[p][k],NB[p]);
        }
      }

      std::vector<double> interEnergy(interFunctionals.size(), 0.);
      for(const auto& pair : pairs)
        for(auto& e : integrateInterEnergy[pair.functionalIndex])
          interEnergy[pair.functionalIndex] += 4*M_PI*e;

      std::vector<double> intraEnergy(nP, 0.);
      for(size_t p = 0; p < nP; p++)
        for(auto& e : integrateXCEnergy[p]) intraEnergy[p] += 4*M_PI*e;

#ifdef CQ_ENABLE_MPI
      {
        dcomplex* mpiScr = (mpiRank == 0) ? CQMemManager::get().malloc<dcomplex>(NBmax*NBmax) : nullptr;
        for(size_t p = 0; p < nP; p++)
          for(auto &V : VXC_SZYX[p]) {
            MPIReduce(V, NB[p]*NB[p], mpiScr, 0, intComm);
            if( mpiRank == 0 ) std::copy_n(mpiScr, NB[p]*NB[p], V);
          }
        if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);
        for(const auto& pair : pairs)
          interEnergy[pair.functionalIndex] =
            MPIReduce(interEnergy[pair.functionalIndex], 0, intComm);
        for(size_t p = 0; p < nP; p++) intraEnergy[p] = MPIReduce(intraEnergy[p], 0, intComm);
      }
#endif

      // Add target VXC matrices and update the active energy caches.
      if( mpiRank == 0 ) {
        for(size_t p = 0; p < nP; p++) {
          if( parts[p].formVXC ) *ss[p]->fockMatrix += *VXC[p];
          if( parts[p].formIntraXC and ks[p] )
            ks[p]->XCEnergy = intraEnergy[p];
        }
        for(const auto& pair : pairs)
          interFunctionals[pair.functionalIndex].energy = interEnergy[pair.functionalIndex];
      }

      // ----------------------------------------------------------------------
      //  Free scratch
      // ----------------------------------------------------------------------
      CQMemManager::get().free(SCRATCHNBNB, SCRATCHNBNP, epsEval, epcEval, dVU_n,
        ZrhoVar1, ZMAT);
      for(size_t p = 0; p < nP; p++)
        if( ZACC_RAW[p] ) CQMemManager::get().free(ZACC_RAW[p]);
      if( anyGGA )
        CQMemManager::get().free(dVU_gamma, ZgammaVar1, ZgammaVar2, ZgammaVar3,
          cross_U_gamma, cross_dVU_gamma);
      if( anyIntraGGA )
        CQMemManager::get().free(KScratch, HScratch);
      if( multiFunc ) {
        CQMemManager::get().free(epsSCR, dVU_n_SCR);
        if( anyGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }
      for(size_t e = 0; e < nEval; e++) {
        CQMemManager::get().free(DenS[e], U_n[e]);
        if( DenZ[e] )    CQMemManager::get().free(DenZ[e]);
        if( GDenS[e] )   CQMemManager::get().free(GDenS[e]);
        if( GDenZ[e] )   CQMemManager::get().free(GDenZ[e]);
        if( U_gamma[e] ) CQMemManager::get().free(U_gamma[e]);
      }
      for(size_t p = 0; p < nP; p++)
        if( intVXC_RAW[p] ) CQMemManager::get().free(intVXC_RAW[p]);

      SetLAThreads(LAThreads);
      MPICommFree(intComm);

    } // valid intComm

    MPI_Barrier(this->comm);

    ProgramTimer::tock("Form VXC");

  }; // MultiParticleSS<dcomplex,dcomplex>::formXCInHouse

} // namespace ChronusQ
