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

#include <singleslater/kohnsham.hpp>

#include <dft/epc.hpp>
#include <dft/util.hpp>
#include <grid/integrator.hpp>
#include <util/mpi.hpp>
#include <util/threads.hpp>
#include <util/timer.hpp>
#include <array>
#include <limits>
#include <type_traits>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  std::vector<double> MultiParticleSS<MatsT,IntsT>::formXCGradient(
      EMPerturbation& pert, const XCTerms& xcTerms) {

    if( intParam.useGauXC )
      return formXCGradientGauXC(pert, xcTerms);
    else
      return formXCGradientInHouse(pert, xcTerms);
  }


  /**
   *  \brief EXC nuclear gradient via Inhouse engine.
   */
  template <typename MatsT, typename IntsT>
  std::vector<double> MultiParticleSS<MatsT,IntsT>::formXCGradientInHouse(
      EMPerturbation& pert, const XCTerms& xcTerms) {

    const auto& parts = xcTerms.subsystems;
    const auto& pairs = xcTerms.interTerms;
    const size_t nAtoms = this->molecule().nAtoms;
    const size_t nGradient = 3*nAtoms;
    std::vector<double> gradient(nGradient,0.);

    if( parts.empty() ) return gradient;

    ProgramTimer::tick("Form MultiParticle XC Gradient");

    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    const size_t nthreads = GetNumThreads();
    const size_t LAThreads = GetLAThreads();
    const size_t mpiRank = MPIRank(this->comm);
    const size_t mpiSize = MPISize(this->comm);
    const size_t nP = parts.size();
    const size_t nEval = nP;

    int color = ((mpiSize < nAtoms) or (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL )
#endif
    {
      SetLAThreads(1);

      const size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
      const size_t NTNPPB = nthreads * NPtsMaxPerBatch;

      // Gather the subsystem data needed by the active functionals.
      std::vector<SubSSPtr> ss(nP);
      std::vector<BasisSet*> bases(nP);
      std::vector<size_t> NB(nP);
      std::vector<bool> hasZ(nP);
      std::vector<bool> intraGGA(nP,false);
      std::vector<std::shared_ptr<KohnSham<MatsT,IntsT>>> ks(nP);
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<IntsT>>> gridPDM(nP);
      size_t NBmax = 0;
      bool multipleIntraFunctionals = false;

      for(size_t p = 0; p < nP; ++p) {
        ss[p] = subsystems.at(parts[p].label);
        bases[p] = &ss[p]->basisSet();
        NB[p] = bases[p]->nBasis;
        NBmax = std::max(NBmax,NB[p]);
        hasZ[p] = ss[p]->onePDM->hasZ();

        if( ss[p]->onePDM->hasXY() )
          CErr("Nuclear gradient is not implemented for 2-component systems (GHF or X2C)");

        if( parts[p].formIntraXC ) {
          ks[p] = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss[p]);
          if( not ks[p] )
            CErr("Intra-particle XC gradient requested for a non-Kohn-Sham subsystem");
          intraGGA[p] = std::any_of(ks[p]->functionals.begin(),ks[p]->functionals.end(),
            [](const std::shared_ptr<DFTFunctional>& functional) {
              return functional->isGGA();
            });
          multipleIntraFunctionals = multipleIntraFunctionals or
            ks[p]->functionals.size() > 1;
        }

        if constexpr (std::is_same_v<MatsT,IntsT>)
          gridPDM[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(*ss[p]->onePDM);
        else
          gridPDM[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(
            ss[p]->onePDM->real_part());
      }

      // EPC gradients currently support one LDA EPC17 functional per pair.
      for(const auto& pair : pairs) {
        const auto& functionals = interFunctionals[pair.functionalIndex].functionals;
        if( functionals.size() > 1 )
          CErr("Multiple EPC functionals in an XC gradient are NYI");
        if( std::any_of(functionals.begin(),functionals.end(),
              [](const std::shared_ptr<DFTFunctional>& functional) {
                return functional->isGGA();
              }) )
          CErr("GGA EPC gradient NYI");

        const bool firstElectronic = ss[pair.firstIndex]->particle.charge < 0.;
        const bool secondElectronic = ss[pair.secondIndex]->particle.charge < 0.;
        if( firstElectronic == secondElectronic )
          CErr("EPC gradient requires exactly one electronic subsystem");
      }

      struct DensityScratch {
        std::vector<double> denS, denZ;
        std::vector<double> gDenS, gDenZ;
        std::array<std::vector<double>,3> nuclearDenS, nuclearDenZ;
        std::array<std::array<std::vector<double>,3>,3> nuclearGDenS, nuclearGDenZ;
        std::vector<double> uN, uGamma;
        std::array<std::vector<double>,3> gradUN, gradUGamma;

        DensityScratch(size_t nAtoms, size_t NTNPPB, bool hasZ, bool isGGA) :
          denS(NTNPPB), uN(2*NTNPPB) {

          if( hasZ ) denZ.resize(NTNPPB);
          if( isGGA ) {
            gDenS.resize(3*NTNPPB);
            uGamma.resize(3*NTNPPB);
            if( hasZ ) gDenZ.resize(3*NTNPPB);
          }

          for(size_t xyz = 0; xyz < 3; ++xyz) {
            nuclearDenS[xyz].resize(nAtoms*NTNPPB);
            gradUN[xyz].resize(nAtoms*2*NTNPPB);
            if( hasZ ) nuclearDenZ[xyz].resize(nAtoms*NTNPPB);
            if( isGGA ) {
              gradUGamma[xyz].resize(nAtoms*3*NTNPPB);
              for(size_t uvw = 0; uvw < 3; ++uvw) {
                nuclearGDenS[xyz][uvw].resize(nAtoms*NTNPPB);
                if( hasZ ) nuclearGDenZ[xyz][uvw].resize(nAtoms*NTNPPB);
              }
            }
          }
        }
      };

      // GIAO uses one + phase copy for intra/electronic densities and one -
      //   phase copy for the non-electronic side of each inter-particle term.
      std::vector<DensityScratch> density;
      density.reserve(nEval);
      for(size_t e = 0; e < nEval; ++e) {
        const size_t p = e % nP;
        density.emplace_back(nAtoms,NTNPPB,hasZ[p],intraGGA[p]);
      }

      // Basis-contraction scratch is shared by subsystems because their
      //   densities are evaluated sequentially within each grid batch.
      std::vector<IntsT> scratchNBNB(nthreads*NBmax*NBmax);
      std::vector<IntsT> scratchNBNP(nthreads*NBmax*NPtsMaxPerBatch);
      std::array<std::vector<IntsT>,3> scratchNBNP2;
      for(auto& scratch : scratchNBNP2)
        scratch.resize(nthreads*NBmax*NPtsMaxPerBatch);

      std::vector<double> epsEval(nthreads*NPtsMaxPerBatch);
      std::vector<double> epsEvalAux(nthreads*NPtsMaxPerBatch);
      std::vector<double> dVUn(nthreads*2*NPtsMaxPerBatch);
      std::vector<double> dVUnAux(nthreads*2*NPtsMaxPerBatch);
      std::vector<double> dVUgamma(nthreads*3*NPtsMaxPerBatch);
      std::vector<double> epsScratch;
      std::vector<double> dVUnScratch;
      std::vector<double> dVUgammaScratch;
      if( multipleIntraFunctionals ) {
        epsScratch.resize(nthreads*NPtsMaxPerBatch);
        dVUnScratch.resize(nthreads*2*NPtsMaxPerBatch);
        dVUgammaScratch.resize(nthreads*3*NPtsMaxPerBatch);
      }

      std::vector<std::vector<double>> threadGradient(
        nthreads,std::vector<double>(nGradient,0.));

      auto epsScreenOf = [&]() {
        double epsilon = intParam.epsilon / nAtoms / intParam.nAng / intParam.nRad;
        return std::max(epsilon,std::numeric_limits<double>::epsilon());
      };

      auto gradientBuild = [&](size_t& res, std::vector<cart_t>& batch,
        std::vector<double>& weights, std::vector<size_t> NBE,
        std::vector<IntsT*> BasisEval,
        std::vector<std::vector<size_t>>& batchEvalShells,
        std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut) {

        const size_t NPts = batch.size();
        const size_t thread = GetThreadID();
        const size_t TIDNPPB = thread*NPtsMaxPerBatch;
        const double epsScreen = epsScreenOf();

        IntsT* scratchNBNBLoc = scratchNBNB.data() + thread*NBmax*NBmax;
        IntsT* scratchNBNPLoc = scratchNBNP.data() + thread*NBmax*NPtsMaxPerBatch;
        std::vector<IntsT*> scratchNBNP2Loc(3);
        for(size_t xyz = 0; xyz < 3; ++xyz)
          scratchNBNP2Loc[xyz] = scratchNBNP2[xyz].data() +
            thread*NBmax*NPtsMaxPerBatch;

        // Evaluate each active density and its nuclear derivative once.
        for(size_t e = 0; e < nEval; ++e) {
          const size_t p = e % nP;
          auto& den = density[e];
          const bool isGGA = intraGGA[p];
          const SHELL_EVAL_TYPE typ = isGGA ? GRADIENT : NOGRAD;
          const size_t NDer = isGGA ? 4 : 1;
          IntsT* basisGradEval = BasisEval[e] + NDer*NPts*NBE[e];

          std::array<std::vector<double*>,3> nuclearDenSPtr, nuclearDenZPtr;
          std::array<std::array<std::vector<double*>,3>,3>
            nuclearGDenSPtr, nuclearGDenZPtr;
          std::array<std::vector<double*>,3> gradUNPtr, gradUGammaPtr;

          for(size_t xyz = 0; xyz < 3; ++xyz) {
            nuclearDenSPtr[xyz].resize(nAtoms);
            if( hasZ[p] ) nuclearDenZPtr[xyz].resize(nAtoms);
            gradUNPtr[xyz].resize(nAtoms);
            if( isGGA ) gradUGammaPtr[xyz].resize(nAtoms);

            for(size_t atom = 0; atom < nAtoms; ++atom) {
              nuclearDenSPtr[xyz][atom] = den.nuclearDenS[xyz].data() +
                atom*NTNPPB + TIDNPPB;
              gradUNPtr[xyz][atom] = den.gradUN[xyz].data() +
                atom*2*NTNPPB + 2*TIDNPPB;
              if( hasZ[p] )
                nuclearDenZPtr[xyz][atom] = den.nuclearDenZ[xyz].data() +
                  atom*NTNPPB + TIDNPPB;
              if( isGGA ) {
                gradUGammaPtr[xyz][atom] = den.gradUGamma[xyz].data() +
                  atom*3*NTNPPB + 3*TIDNPPB;
                for(size_t uvw = 0; uvw < 3; ++uvw) {
                  nuclearGDenSPtr[xyz][uvw].resize(nAtoms);
                  nuclearGDenSPtr[xyz][uvw][atom] =
                    den.nuclearGDenS[xyz][uvw].data() + atom*NTNPPB + TIDNPPB;
                  if( hasZ[p] ) {
                    nuclearGDenZPtr[xyz][uvw].resize(nAtoms);
                    nuclearGDenZPtr[xyz][uvw][atom] =
                      den.nuclearGDenZ[xyz][uvw].data() + atom*NTNPPB + TIDNPPB;
                  }
                }
              }
            }
          }

          std::vector<std::vector<IntsT*>> unusedMatrixScratch(nAtoms);
          std::vector<std::vector<IntsT*>> unusedPDMGradient(nAtoms);

          double* gDenS = isGGA ? den.gDenS.data() + 3*TIDNPPB : nullptr;
          evalDen(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
            scratchNBNPLoc,gridPDM[p]->S().pointer(),den.denS.data() + TIDNPPB,
            gDenS,gDenS ? gDenS + NPts : nullptr,gDenS ? gDenS + 2*NPts : nullptr,
            BasisEval[e]);

          evalDenGrad(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
            scratchNBNPLoc,unusedMatrixScratch,unusedMatrixScratch,scratchNBNP2Loc,
            unusedPDMGradient,gridPDM[p]->S().pointer(),nuclearDenSPtr[0],
            nuclearDenSPtr[1],nuclearDenSPtr[2],nuclearGDenSPtr[0][0],
            nuclearGDenSPtr[0][1],nuclearGDenSPtr[0][2],nuclearGDenSPtr[1][0],
            nuclearGDenSPtr[1][1],nuclearGDenSPtr[1][2],nuclearGDenSPtr[2][0],
            nuclearGDenSPtr[2][1],nuclearGDenSPtr[2][2],BasisEval[e],
            basisGradEval,this->molecule(),*bases[p]);

          double* gDenZ = nullptr;
          if( hasZ[p] ) {
            gDenZ = isGGA ? den.gDenZ.data() + 3*TIDNPPB : nullptr;
            evalDen(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
              scratchNBNPLoc,gridPDM[p]->Z().pointer(),den.denZ.data() + TIDNPPB,
              gDenZ,gDenZ ? gDenZ + NPts : nullptr,gDenZ ? gDenZ + 2*NPts : nullptr,
              BasisEval[e]);

            evalDenGrad(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
              scratchNBNPLoc,unusedMatrixScratch,unusedMatrixScratch,scratchNBNP2Loc,
              unusedPDMGradient,gridPDM[p]->Z().pointer(),nuclearDenZPtr[0],
              nuclearDenZPtr[1],nuclearDenZPtr[2],nuclearGDenZPtr[0][0],
              nuclearGDenZPtr[0][1],nuclearGDenZPtr[0][2],nuclearGDenZPtr[1][0],
              nuclearGDenZPtr[1][1],nuclearGDenZPtr[1][2],nuclearGDenZPtr[2][0],
              nuclearGDenZPtr[2][1],nuclearGDenZPtr[2][2],BasisEval[e],
              basisGradEval,this->molecule(),*bases[p]);
          }

          mkAuxVar(ss[p]->onePDM,isGGA,epsScreen,NPts,
            den.denS.data() + TIDNPPB,hasZ[p] ? den.denZ.data() + TIDNPPB : nullptr,
            nullptr,nullptr,gDenS,gDenS ? gDenS + NPts : nullptr,
            gDenS ? gDenS + 2*NPts : nullptr,gDenZ,gDenZ ? gDenZ + NPts : nullptr,
            gDenZ ? gDenZ + 2*NPts : nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,
            nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,
            nullptr,nullptr,den.uN.data() + 2*TIDNPPB,
            isGGA ? den.uGamma.data() + 3*TIDNPPB : nullptr);

          mkAuxVarGrad(ss[p]->onePDM,isGGA,epsScreen,NPts,nuclearDenSPtr[0],
            nuclearDenSPtr[1],nuclearDenSPtr[2],nuclearDenZPtr[0],
            nuclearDenZPtr[1],nuclearDenZPtr[2],gDenS,
            gDenS ? gDenS + NPts : nullptr,gDenS ? gDenS + 2*NPts : nullptr,
            gDenZ,gDenZ ? gDenZ + NPts : nullptr,gDenZ ? gDenZ + 2*NPts : nullptr,
            nuclearGDenSPtr[0][0],nuclearGDenSPtr[0][1],nuclearGDenSPtr[0][2],
            nuclearGDenSPtr[1][0],nuclearGDenSPtr[1][1],nuclearGDenSPtr[1][2],
            nuclearGDenSPtr[2][0],nuclearGDenSPtr[2][1],nuclearGDenSPtr[2][2],
            nuclearGDenZPtr[0][0],nuclearGDenZPtr[0][1],nuclearGDenZPtr[0][2],
            nuclearGDenZPtr[1][0],nuclearGDenZPtr[1][1],nuclearGDenZPtr[1][2],
            nuclearGDenZPtr[2][0],nuclearGDenZPtr[2][1],nuclearGDenZPtr[2][2],
            gradUNPtr[0],gradUNPtr[1],gradUNPtr[2],gradUGammaPtr[0],
            gradUGammaPtr[1],gradUGammaPtr[2],nAtoms);
        }

        double* epsEvalLoc = epsEval.data() + TIDNPPB;
        double* epsEvalAuxLoc = epsEvalAux.data() + TIDNPPB;
        double* dVUnLoc = dVUn.data() + 2*TIDNPPB;
        double* dVUnAuxLoc = dVUnAux.data() + 2*TIDNPPB;
        double* dVUgammaLoc = dVUgamma.data() + 3*TIDNPPB;
        double* epsScratchLoc = multipleIntraFunctionals ?
          epsScratch.data() + TIDNPPB : nullptr;
        double* dVUnScratchLoc = multipleIntraFunctionals ?
          dVUnScratch.data() + 2*TIDNPPB : nullptr;
        double* dVUgammaScratchLoc = multipleIntraFunctionals ?
          dVUgammaScratch.data() + 3*TIDNPPB : nullptr;

        // Add each subsystem's intra-particle XC derivative.
        for(size_t p = 0; p < nP; ++p) {
          if( not ks[p] ) continue;

          auto& den = density[p];
          loadVXCder(ks[p]->functionals,NPts,den.uN.data() + 2*TIDNPPB,
            intraGGA[p] ? den.uGamma.data() + 3*TIDNPPB : nullptr,epsEvalLoc,
            dVUnLoc,intraGGA[p] ? dVUgammaLoc : nullptr,epsScratchLoc,
            dVUnScratchLoc,intraGGA[p] ? dVUgammaScratchLoc : nullptr);

          for(size_t atom = 0; atom < nAtoms; ++atom)
            for(size_t xyz = 0; xyz < 3; ++xyz)
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(
                intraGGA[p],NPts,weights,dVUnLoc,
                intraGGA[p] ? dVUgammaLoc : nullptr,
                den.gradUN[xyz].data() + atom*2*NTNPPB + 2*TIDNPPB,
                intraGGA[p] ? den.gradUGamma[xyz].data() +
                  atom*3*NTNPPB + 3*TIDNPPB : nullptr);
        }

        // Add each inter-particle derivative once, with the electronic and
        //   non-electronic density derivatives contributing to the same force.
        for(const auto& pair : pairs) {
          const bool firstElectronic = ss[pair.firstIndex]->particle.charge < 0.;
          const size_t electronic = firstElectronic ? pair.firstIndex : pair.secondIndex;
          const size_t other = firstElectronic ? pair.secondIndex : pair.firstIndex;
          const size_t electronicSlot = electronic;
          const size_t otherSlot = other;
          auto& electronicDen = density[electronicSlot];
          auto& otherDen = density[otherSlot];

          std::fill_n(epsEvalLoc,NPts,0.);
          std::fill_n(epsEvalAuxLoc,NPts,0.);
          loadEPCGradder(interFunctionals[pair.functionalIndex].functionals,NPts,
            electronicDen.uN.data() + 2*TIDNPPB,
            otherDen.uN.data() + 2*TIDNPPB,epsEvalLoc,epsEvalAuxLoc,
            dVUnLoc,dVUnAuxLoc);

          for(size_t atom = 0; atom < nAtoms; ++atom)
            for(size_t xyz = 0; xyz < 3; ++xyz) {
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(false,
                NPts,weights,dVUnLoc,nullptr,electronicDen.gradUN[xyz].data() +
                  atom*2*NTNPPB + 2*TIDNPPB,nullptr);
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(false,
                NPts,weights,dVUnAuxLoc,nullptr,otherDen.gradUN[xyz].data() +
                  atom*2*NTNPPB + 2*TIDNPPB,nullptr);
            }
        }
      };

      // Evaluate every active basis on the same molecular grid.
      std::vector<BasisSet*> evalBases(nEval);
      std::vector<SHELL_EVAL_TYPE> typs(nEval);
      for(size_t e = 0; e < nEval; ++e) {
        const size_t p = e % nP;
        evalBases[e] = bases[p];
        typs[e] = intraGGA[p] ? GRADIENT : NOGRAD;
      }

      BeckeIntegrator<EulerMac> integrator(intComm,this->molecule(),evalBases,typs,
        EulerMac(intParam.nRad),intParam.nAng,intParam.nRadPerBatch,intParam.epsilon);
      integrator.turn_on_grad();

      size_t nGridPts = 0;
      integrator.integrateN<size_t>(nGridPts,gradientBuild);

      for(const auto& localGradient : threadGradient)
        for(size_t i = 0; i < nGradient; ++i)
          gradient[i] += 4*M_PI*localGradient[i];

      SetLAThreads(LAThreads);
      MPICommFree(intComm);
    }

#ifdef CQ_ENABLE_MPI
    std::vector<double> reducedGradient(nGradient);
    MPIAllReduce(gradient.data(),nGradient,reducedGradient.data(),this->comm);
    gradient.swap(reducedGradient);
#endif

    MPI_Barrier(this->comm);
    ProgramTimer::tock("Form MultiParticle XC Gradient");
    return gradient;
  }

  /**
   *  \brief Complex GIAO EXC nuclear gradient via the in-house engine.
   */
  template <>
  inline std::vector<double> MultiParticleSS<dcomplex,dcomplex>::formXCGradientInHouse(
      EMPerturbation& pert, const XCTerms& xcTerms) {

    const auto& parts = xcTerms.subsystems;
    const auto& pairs = xcTerms.interTerms;
    const size_t nAtoms = this->molecule().nAtoms;
    const size_t nGradient = 3*nAtoms;
    std::vector<double> gradient(nGradient,0.);

    if( parts.empty() ) return gradient;

    ProgramTimer::tick("Form MultiParticle XC Gradient");

    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    const size_t nthreads = GetNumThreads();
    const size_t LAThreads = GetLAThreads();
    const size_t mpiRank = MPIRank(this->comm);
    const size_t mpiSize = MPISize(this->comm);
    const size_t nP = parts.size();
    const size_t nEval = 2*nP;

    int color = ((mpiSize < nAtoms) or (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL )
#endif
    {
      SetLAThreads(1);

      const size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
      const size_t NTNPPB = nthreads * NPtsMaxPerBatch;

      // Gather the subsystem data needed by the active functionals.
      std::vector<SubSSPtr> ss(nP);
      std::vector<BasisSet*> bases(nP);
      std::vector<size_t> NB(nP);
      std::vector<bool> hasZ(nP);
      std::vector<bool> intraGGA(nP,false);
      std::vector<std::shared_ptr<KohnSham<dcomplex,dcomplex>>> ks(nP);
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> gridPDM(nP);
      size_t NBmax = 0;
      bool multipleIntraFunctionals = false;

      for(size_t p = 0; p < nP; ++p) {
        ss[p] = subsystems.at(parts[p].label);
        bases[p] = &ss[p]->basisSet();
        NB[p] = bases[p]->nBasis;
        NBmax = std::max(NBmax,NB[p]);
        hasZ[p] = ss[p]->onePDM->hasZ();

        if( ss[p]->onePDM->hasXY() )
          CErr("Nuclear gradient is not implemented for 2-component systems (GHF or X2C)");

        if( parts[p].formIntraXC ) {
          ks[p] = std::dynamic_pointer_cast<KohnSham<dcomplex,dcomplex>>(ss[p]);
          if( not ks[p] )
            CErr("Intra-particle XC gradient requested for a non-Kohn-Sham subsystem");
          intraGGA[p] = std::any_of(ks[p]->functionals.begin(),ks[p]->functionals.end(),
            [](const std::shared_ptr<DFTFunctional>& functional) {
              return functional->isGGA();
            });
          multipleIntraFunctionals = multipleIntraFunctionals or
            ks[p]->functionals.size() > 1;
        }

        gridPDM[p] = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(*ss[p]->onePDM);
      }

      // EPC gradients currently support one LDA EPC17 functional per pair.
      for(const auto& pair : pairs) {
        const auto& functionals = interFunctionals[pair.functionalIndex].functionals;
        if( functionals.size() > 1 )
          CErr("Multiple EPC functionals in an XC gradient are NYI");
        if( std::any_of(functionals.begin(),functionals.end(),
              [](const std::shared_ptr<DFTFunctional>& functional) {
                return functional->isGGA();
              }) )
          CErr("GGA EPC gradient NYI");

        const bool firstElectronic = ss[pair.firstIndex]->particle.charge < 0.;
        const bool secondElectronic = ss[pair.secondIndex]->particle.charge < 0.;
        if( firstElectronic == secondElectronic )
          CErr("EPC gradient requires exactly one electronic subsystem");
      }

      struct DensityScratch {
        std::vector<double> denS, denZ;
        std::vector<double> gDenS, gDenZ;
        std::array<std::vector<double>,3> nuclearDenS, nuclearDenZ;
        std::array<std::array<std::vector<double>,3>,3> nuclearGDenS, nuclearGDenZ;
        std::vector<double> uN, uGamma;
        std::array<std::vector<double>,3> gradUN, gradUGamma;

        DensityScratch(size_t nAtoms, size_t NTNPPB, bool hasZ, bool isGGA) :
          denS(NTNPPB), uN(2*NTNPPB) {

          if( hasZ ) denZ.resize(NTNPPB);
          if( isGGA ) {
            gDenS.resize(3*NTNPPB);
            uGamma.resize(3*NTNPPB);
            if( hasZ ) gDenZ.resize(3*NTNPPB);
          }

          for(size_t xyz = 0; xyz < 3; ++xyz) {
            nuclearDenS[xyz].resize(nAtoms*NTNPPB);
            gradUN[xyz].resize(nAtoms*2*NTNPPB);
            if( hasZ ) nuclearDenZ[xyz].resize(nAtoms*NTNPPB);
            if( isGGA ) {
              gradUGamma[xyz].resize(nAtoms*3*NTNPPB);
              for(size_t uvw = 0; uvw < 3; ++uvw) {
                nuclearGDenS[xyz][uvw].resize(nAtoms*NTNPPB);
                if( hasZ ) nuclearGDenZ[xyz][uvw].resize(nAtoms*NTNPPB);
              }
            }
          }
        }
      };

      // GIAO uses one + phase copy for intra/electronic densities and one -
      //   phase copy for the non-electronic side of each inter-particle term.
      std::vector<DensityScratch> density;
      density.reserve(nEval);
      for(size_t e = 0; e < nEval; ++e) {
        const size_t p = e % nP;
        density.emplace_back(nAtoms,NTNPPB,hasZ[p],intraGGA[p]);
      }

      // Basis-contraction scratch is shared by subsystems because their
      //   densities are evaluated sequentially within each grid batch.
      std::vector<dcomplex> scratchNBNB(nthreads*NBmax*NBmax);
      std::vector<dcomplex> scratchNBNP(nthreads*NBmax*NPtsMaxPerBatch);
      std::array<std::vector<dcomplex>,3> scratchNBNP2;
      for(auto& scratch : scratchNBNP2)
        scratch.resize(nthreads*NBmax*NPtsMaxPerBatch);

      std::vector<double> epsEval(nthreads*NPtsMaxPerBatch);
      std::vector<double> epsEvalAux(nthreads*NPtsMaxPerBatch);
      std::vector<double> dVUn(nthreads*2*NPtsMaxPerBatch);
      std::vector<double> dVUnAux(nthreads*2*NPtsMaxPerBatch);
      std::vector<double> dVUgamma(nthreads*3*NPtsMaxPerBatch);
      std::vector<double> epsScratch;
      std::vector<double> dVUnScratch;
      std::vector<double> dVUgammaScratch;
      if( multipleIntraFunctionals ) {
        epsScratch.resize(nthreads*NPtsMaxPerBatch);
        dVUnScratch.resize(nthreads*2*NPtsMaxPerBatch);
        dVUgammaScratch.resize(nthreads*3*NPtsMaxPerBatch);
      }

      std::vector<std::vector<double>> threadGradient(
        nthreads,std::vector<double>(nGradient,0.));

      auto epsScreenOf = [&]() {
        double epsilon = intParam.epsilon / nAtoms / intParam.nAng / intParam.nRad;
        return std::max(epsilon,std::numeric_limits<double>::epsilon());
      };

      auto gradientBuild = [&](size_t& res, std::vector<cart_t>& batch,
        std::vector<double>& weights, std::vector<size_t> NBE,
        std::vector<dcomplex*> BasisEval,
        std::vector<std::vector<size_t>>& batchEvalShells,
        std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut) {

        const size_t NPts = batch.size();
        const size_t thread = GetThreadID();
        const size_t TIDNPPB = thread*NPtsMaxPerBatch;
        const double epsScreen = epsScreenOf();

        dcomplex* scratchNBNBLoc = scratchNBNB.data() + thread*NBmax*NBmax;
        dcomplex* scratchNBNPLoc = scratchNBNP.data() + thread*NBmax*NPtsMaxPerBatch;
        std::vector<dcomplex*> scratchNBNP2Loc(3);
        for(size_t xyz = 0; xyz < 3; ++xyz)
          scratchNBNP2Loc[xyz] = scratchNBNP2[xyz].data() +
            thread*NBmax*NPtsMaxPerBatch;

        // Evaluate each active density and its nuclear derivative once.
        for(size_t e = 0; e < nEval; ++e) {
          const size_t p = e % nP;
          auto& den = density[e];
          const bool isGGA = intraGGA[p];
          const SHELL_EVAL_TYPE typ = isGGA ? GRADIENT : NOGRAD;
          const size_t NDer = isGGA ? 4 : 1;
          dcomplex* basisGradEval = BasisEval[e] + NDer*NPts*NBE[e];

          std::array<std::vector<double*>,3> nuclearDenSPtr, nuclearDenZPtr;
          std::array<std::array<std::vector<double*>,3>,3>
            nuclearGDenSPtr, nuclearGDenZPtr;
          std::array<std::vector<double*>,3> gradUNPtr, gradUGammaPtr;

          for(size_t xyz = 0; xyz < 3; ++xyz) {
            nuclearDenSPtr[xyz].resize(nAtoms);
            if( hasZ[p] ) nuclearDenZPtr[xyz].resize(nAtoms);
            gradUNPtr[xyz].resize(nAtoms);
            if( isGGA ) gradUGammaPtr[xyz].resize(nAtoms);

            for(size_t atom = 0; atom < nAtoms; ++atom) {
              nuclearDenSPtr[xyz][atom] = den.nuclearDenS[xyz].data() +
                atom*NTNPPB + TIDNPPB;
              gradUNPtr[xyz][atom] = den.gradUN[xyz].data() +
                atom*2*NTNPPB + 2*TIDNPPB;
              if( hasZ[p] )
                nuclearDenZPtr[xyz][atom] = den.nuclearDenZ[xyz].data() +
                  atom*NTNPPB + TIDNPPB;
              if( isGGA ) {
                gradUGammaPtr[xyz][atom] = den.gradUGamma[xyz].data() +
                  atom*3*NTNPPB + 3*TIDNPPB;
                for(size_t uvw = 0; uvw < 3; ++uvw) {
                  nuclearGDenSPtr[xyz][uvw].resize(nAtoms);
                  nuclearGDenSPtr[xyz][uvw][atom] =
                    den.nuclearGDenS[xyz][uvw].data() + atom*NTNPPB + TIDNPPB;
                  if( hasZ[p] ) {
                    nuclearGDenZPtr[xyz][uvw].resize(nAtoms);
                    nuclearGDenZPtr[xyz][uvw][atom] =
                      den.nuclearGDenZ[xyz][uvw].data() + atom*NTNPPB + TIDNPPB;
                  }
                }
              }
            }
          }

          std::vector<std::vector<dcomplex*>> unusedMatrixScratch(nAtoms);
          std::vector<std::vector<dcomplex*>> unusedPDMGradient(nAtoms);

          double* gDenS = isGGA ? den.gDenS.data() + 3*TIDNPPB : nullptr;
          evalDen(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
            scratchNBNPLoc,gridPDM[p]->S().pointer(),den.denS.data() + TIDNPPB,
            gDenS,gDenS ? gDenS + NPts : nullptr,gDenS ? gDenS + 2*NPts : nullptr,
            BasisEval[e]);

          evalDenGrad(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
            scratchNBNPLoc,unusedMatrixScratch,unusedMatrixScratch,scratchNBNP2Loc,
            unusedPDMGradient,gridPDM[p]->S().pointer(),nuclearDenSPtr[0],
            nuclearDenSPtr[1],nuclearDenSPtr[2],nuclearGDenSPtr[0][0],
            nuclearGDenSPtr[0][1],nuclearGDenSPtr[0][2],nuclearGDenSPtr[1][0],
            nuclearGDenSPtr[1][1],nuclearGDenSPtr[1][2],nuclearGDenSPtr[2][0],
            nuclearGDenSPtr[2][1],nuclearGDenSPtr[2][2],BasisEval[e],
            basisGradEval,this->molecule(),*bases[p]);

          double* gDenZ = nullptr;
          if( hasZ[p] ) {
            gDenZ = isGGA ? den.gDenZ.data() + 3*TIDNPPB : nullptr;
            evalDen(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
              scratchNBNPLoc,gridPDM[p]->Z().pointer(),den.denZ.data() + TIDNPPB,
              gDenZ,gDenZ ? gDenZ + NPts : nullptr,gDenZ ? gDenZ + 2*NPts : nullptr,
              BasisEval[e]);

            evalDenGrad(typ,NPts,NBE[e],NB[p],subMatCut[e],scratchNBNBLoc,
              scratchNBNPLoc,unusedMatrixScratch,unusedMatrixScratch,scratchNBNP2Loc,
              unusedPDMGradient,gridPDM[p]->Z().pointer(),nuclearDenZPtr[0],
              nuclearDenZPtr[1],nuclearDenZPtr[2],nuclearGDenZPtr[0][0],
              nuclearGDenZPtr[0][1],nuclearGDenZPtr[0][2],nuclearGDenZPtr[1][0],
              nuclearGDenZPtr[1][1],nuclearGDenZPtr[1][2],nuclearGDenZPtr[2][0],
              nuclearGDenZPtr[2][1],nuclearGDenZPtr[2][2],BasisEval[e],
              basisGradEval,this->molecule(),*bases[p]);
          }

          mkAuxVar(ss[p]->onePDM,isGGA,epsScreen,NPts,
            den.denS.data() + TIDNPPB,hasZ[p] ? den.denZ.data() + TIDNPPB : nullptr,
            nullptr,nullptr,gDenS,gDenS ? gDenS + NPts : nullptr,
            gDenS ? gDenS + 2*NPts : nullptr,gDenZ,gDenZ ? gDenZ + NPts : nullptr,
            gDenZ ? gDenZ + 2*NPts : nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,
            nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,
            nullptr,nullptr,den.uN.data() + 2*TIDNPPB,
            isGGA ? den.uGamma.data() + 3*TIDNPPB : nullptr);

          mkAuxVarGrad(ss[p]->onePDM,isGGA,epsScreen,NPts,nuclearDenSPtr[0],
            nuclearDenSPtr[1],nuclearDenSPtr[2],nuclearDenZPtr[0],
            nuclearDenZPtr[1],nuclearDenZPtr[2],gDenS,
            gDenS ? gDenS + NPts : nullptr,gDenS ? gDenS + 2*NPts : nullptr,
            gDenZ,gDenZ ? gDenZ + NPts : nullptr,gDenZ ? gDenZ + 2*NPts : nullptr,
            nuclearGDenSPtr[0][0],nuclearGDenSPtr[0][1],nuclearGDenSPtr[0][2],
            nuclearGDenSPtr[1][0],nuclearGDenSPtr[1][1],nuclearGDenSPtr[1][2],
            nuclearGDenSPtr[2][0],nuclearGDenSPtr[2][1],nuclearGDenSPtr[2][2],
            nuclearGDenZPtr[0][0],nuclearGDenZPtr[0][1],nuclearGDenZPtr[0][2],
            nuclearGDenZPtr[1][0],nuclearGDenZPtr[1][1],nuclearGDenZPtr[1][2],
            nuclearGDenZPtr[2][0],nuclearGDenZPtr[2][1],nuclearGDenZPtr[2][2],
            gradUNPtr[0],gradUNPtr[1],gradUNPtr[2],gradUGammaPtr[0],
            gradUGammaPtr[1],gradUGammaPtr[2],nAtoms);
        }

        double* epsEvalLoc = epsEval.data() + TIDNPPB;
        double* epsEvalAuxLoc = epsEvalAux.data() + TIDNPPB;
        double* dVUnLoc = dVUn.data() + 2*TIDNPPB;
        double* dVUnAuxLoc = dVUnAux.data() + 2*TIDNPPB;
        double* dVUgammaLoc = dVUgamma.data() + 3*TIDNPPB;
        double* epsScratchLoc = multipleIntraFunctionals ?
          epsScratch.data() + TIDNPPB : nullptr;
        double* dVUnScratchLoc = multipleIntraFunctionals ?
          dVUnScratch.data() + 2*TIDNPPB : nullptr;
        double* dVUgammaScratchLoc = multipleIntraFunctionals ?
          dVUgammaScratch.data() + 3*TIDNPPB : nullptr;

        // Add each subsystem's intra-particle XC derivative.
        for(size_t p = 0; p < nP; ++p) {
          if( not ks[p] ) continue;

          auto& den = density[p];
          loadVXCder(ks[p]->functionals,NPts,den.uN.data() + 2*TIDNPPB,
            intraGGA[p] ? den.uGamma.data() + 3*TIDNPPB : nullptr,epsEvalLoc,
            dVUnLoc,intraGGA[p] ? dVUgammaLoc : nullptr,epsScratchLoc,
            dVUnScratchLoc,intraGGA[p] ? dVUgammaScratchLoc : nullptr);

          for(size_t atom = 0; atom < nAtoms; ++atom)
            for(size_t xyz = 0; xyz < 3; ++xyz)
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(
                intraGGA[p],NPts,weights,dVUnLoc,
                intraGGA[p] ? dVUgammaLoc : nullptr,
                den.gradUN[xyz].data() + atom*2*NTNPPB + 2*TIDNPPB,
                intraGGA[p] ? den.gradUGamma[xyz].data() +
                  atom*3*NTNPPB + 3*TIDNPPB : nullptr);
        }

        // Add each inter-particle derivative once, with the electronic and
        //   non-electronic density derivatives contributing to the same force.
        for(const auto& pair : pairs) {
          const bool firstElectronic = ss[pair.firstIndex]->particle.charge < 0.;
          const size_t electronic = firstElectronic ? pair.firstIndex : pair.secondIndex;
          const size_t other = firstElectronic ? pair.secondIndex : pair.firstIndex;
          const size_t electronicSlot = electronic;
          const size_t otherSlot = other + nP;
          auto& electronicDen = density[electronicSlot];
          auto& otherDen = density[otherSlot];

          std::fill_n(epsEvalLoc,NPts,0.);
          std::fill_n(epsEvalAuxLoc,NPts,0.);
          loadEPCGradder(interFunctionals[pair.functionalIndex].functionals,NPts,
            electronicDen.uN.data() + 2*TIDNPPB,
            otherDen.uN.data() + 2*TIDNPPB,epsEvalLoc,epsEvalAuxLoc,
            dVUnLoc,dVUnAuxLoc);

          for(size_t atom = 0; atom < nAtoms; ++atom)
            for(size_t xyz = 0; xyz < 3; ++xyz) {
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(false,
                NPts,weights,dVUnLoc,nullptr,electronicDen.gradUN[xyz].data() +
                  atom*2*NTNPPB + 2*TIDNPPB,nullptr);
              threadGradient[thread][3*atom + xyz] += energy_vxc_grad(false,
                NPts,weights,dVUnAuxLoc,nullptr,otherDen.gradUN[xyz].data() +
                  atom*2*NTNPPB + 2*TIDNPPB,nullptr);
            }
        }
      };

      // Evaluate every active basis on the same molecular grid.
      std::vector<BasisSet*> evalBases(nEval);
      std::vector<SHELL_EVAL_TYPE> typs(nEval);
      for(size_t e = 0; e < nEval; ++e) {
        const size_t p = e % nP;
        evalBases[e] = bases[p];
        typs[e] = intraGGA[p] ? GRADIENT : NOGRAD;
      }

      BeckeIntegrator<EulerMac> integrator(intComm,this->molecule(),evalBases,typs,
        EulerMac(intParam.nRad),intParam.nAng,intParam.nRadPerBatch,intParam.epsilon);
      integrator.turn_on_grad();

      size_t nGridPts = 0;
      std::vector<double> phaseScales(nEval,1.);
      for(size_t p = 0; p < nP; ++p) phaseScales[p+nP] = -1.;
      integrator.integrateN<size_t>(nGridPts,gradientBuild,pert,phaseScales);

      for(const auto& localGradient : threadGradient)
        for(size_t i = 0; i < nGradient; ++i)
          gradient[i] += 4*M_PI*localGradient[i];

      SetLAThreads(LAThreads);
      MPICommFree(intComm);
    }

#ifdef CQ_ENABLE_MPI
    std::vector<double> reducedGradient(nGradient);
    MPIAllReduce(gradient.data(),nGradient,reducedGradient.data(),this->comm);
    gradient.swap(reducedGradient);
#endif

    MPI_Barrier(this->comm);
    ProgramTimer::tock("Form MultiParticle XC Gradient");
    return gradient;
  }

  /**
   *  \brief EXC nuclear gradient via GauXC.
  */
  template <typename MatsT, typename IntsT>
  std::vector<double> MultiParticleSS<MatsT,IntsT>::formXCGradientGauXC(
      EMPerturbation&, const XCTerms& terms) {

    const size_t nAtoms = this->molecule().nAtoms;
    std::vector<double> gradient(3*nAtoms, 0.);

    if( terms.subsystems.empty() ) return gradient;

    if constexpr (std::is_same_v<IntsT,dcomplex>) {
      CErr("Complex (dcomplex,dcomplex) GauXC MultiParticleSS XC gradient NYI");
    } else {

      if(!this->gauxcUtils)
        CErr("GauXC MultiParticleSS gradient requested without GauXCUtils");
      if(!this->gauxcUtils->integrator_pointer)
        CErr("GauXC MultiParticleSS gradient requested without GauXC integrator");

      ProgramTimer::tick("Form MultiParticle XC Gradient");

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

      // The gradient plan activates the same functionals as the energy pass.
      // vxc_targets is unused by the gradient driver.
      GauXC::MultiParticleXCTerms xcTerms;
      for(const auto& subsystem : terms.subsystems)
        if(subsystem.formIntraXC)
          xcTerms.active_intra.push_back(globalIndexOf(subsystem.label));
      for(const auto& pair : terms.interTerms)
        xcTerms.active_inter.push_back(pair.functionalIndex);

      std::vector<SubSSPtr> ss(nP);
      std::vector<Eigen::MatrixXd> Ps(nP);
      std::vector<Eigen::MatrixXd> Pz(nP);
      std::vector<GauXC::XCIntegrator<Eigen::MatrixXd>::multiparticle_density> densities(nP);

      // Same density convention as formXCGauXC: RKS passes half the scalar
      //   density, UKS passes scalar + magnetization. GauXC's internal
      //   xmat_fac reconstructs the physical density, so the returned gradient
      //   needs no rescaling.
      for(size_t p = 0; p < nP; ++p) {
        ss[p] = subsystems.at(order_[p]);
        if(ss[p]->onePDM->hasXY())
          CErr("Relativistic (2C/4C) GauXC MultiParticleSS XC gradient NYI");

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

      auto gxGradient = this->gauxcUtils->integrator_pointer->eval_exc_grad(
        densities, functionalSpec, xcTerms);

      if(gxGradient.size() != 3*nAtoms)
        CErr("GauXC MultiParticleSS returned inconsistent gradient dimensions");
      for(size_t i = 0; i < 3*nAtoms; ++i) gradient[i] = gxGradient[i];

      ProgramTimer::tock("Form MultiParticle XC Gradient");
    }

    return gradient;
  }

} // namespace ChronusQ
