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

#include <fockbuilder/neofock.hpp>

#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/gradints/direct.hpp>
#include <particleintegrals/gradints/incore.hpp>

#include <dft.hpp>
#include <dft/util.hpp>
#include <grid/integrator.hpp>
#include <util/threads.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void NEOFockBuilder<MatsT,IntsT>::formepJ(SingleSlater<MatsT,IntsT>& ss,
    bool increment)
  {
    // Validate internal pointers
    if( this->aux_ss == nullptr )
      CErr("aux_ss uninitialized in formepJ!");
    if( this->outMat == nullptr )
      CErr("outMat uninitialized in formepJ!");
    if( contraction == nullptr )
      CErr("contraction uninitialized in formepJ!");

    // Decide onePDM to use
    cqmatrix::PauliSpinorMatrices<MatsT>& contract1PDM
        = increment ? *this->aux_ss->deltaOnePDM : *this->aux_ss->onePDM;

    size_t NB = ss.basisSet().nBasis;

    // Zero out J and K[i]
    if(not increment)
      this->outMat->clear();
  
    std::vector<TwoBodyContraction<MatsT>> contract =
      { {contract1PDM.S().pointer(), this->outMat->pointer(), true, COULOMB} };

    EMPerturbation pert;

    if(contraction->printContractionTiming)
        std::cout << "      " << std::string(ss.particle.charge>0? "Protonic" : "Electronic") << " Subsystem Asymm-Contraction Timing: " << std::endl;
    auto beginEPJContract = tick();
    contraction->twoBodyContract(ss.comm, contract, pert);

///***
    // Print asymm CD contraction total timing (note this is with other density)
    if(contraction->printContractionTiming)
        std::cout << "        " << std::left << std::setw(38) << "Cholesky-Asymm-Contraction duration = " << tock(beginEPJContract) << " s \n" << std::endl;
//***/

    // Get the EPJ Energy (by tracing with own density)
    //double epjEnergy = ss.template computeOBProperty<SCALAR>( this->outMat->pointer());
    //std::cout << "  " << std::string(ss.particle.charge>0? "ProtDensity" : "ElecDensity")
    //    << " Contracted EPJ Energy  = " << std::setprecision(12) <<  epjEnergy  << std::endl;
  }

  template <typename MatsT, typename IntsT>
  std::vector<double> NEOFockBuilder<MatsT,IntsT>::formepJGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& pert, double xHFX) {

    if( this->aux_ss == nullptr )
      CErr("aux_ss uninitialized in formepJGrad!");
    if( gradTPI == nullptr )
      CErr("gradTPI uninitialized in formepJGrad!");

    size_t NB = ss.basisSet().nBasis;
    size_t nGrad = 3*ss.molecule().nAtoms;

    // Form contraction
    std::unique_ptr<GradContractions<MatsT,IntsT>> contract = nullptr;
    if ( std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>((*gradTPI)[0]) ) {
      contract = std::make_unique<InCore4indexGradContraction<MatsT,IntsT>>(*gradTPI);
    }
    else if ( std::dynamic_pointer_cast<DirectTPI<IntsT>>((*gradTPI)[0]) ) {
      contract = std::make_unique<DirectGradContraction<MatsT,IntsT>>(*gradTPI);
    }
    else
      CErr("Gradients of RI NYI!");
    // Assume that the order of the TPI is the same between the regular and
    //   gradient integrals
    contract->contractSecond = contraction->contractSecond;

    // Create contraction list
    std::vector<std::vector<TwoBodyContraction<MatsT>>> cList;

    std::vector<cqmatrix::Matrix<MatsT>> JList;
    JList.reserve(nGrad);

    for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {
      std::vector<TwoBodyContraction<MatsT>> tempCont;

      // Coulomb
      JList.emplace_back(NB);
      JList.back().clear();
      tempCont.push_back(
         {this->aux_ss->onePDM->S().pointer(), JList.back().pointer(), true, COULOMB}
      );

      cList.push_back(tempCont);
    }

    // Contract to J/K
    contract->gradTwoBodyContract(MPI_COMM_WORLD, true, cList, pert);

    // Contract to gradient
    std::vector<double> gradient;
    cqmatrix::PauliSpinorMatrices<MatsT> twoEGrad(NB, false, false);

    for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {

      // G[S] = -2 * J
      // TODO: The negative here is implicitly taking care of the
      //   electron/proton charges.
      twoEGrad.S() = -2. * JList[iGrad];

      double gradVal = ss.template computeOBProperty<SCALAR>(
        twoEGrad.S().pointer()
      );
      gradient.push_back(0.25*gradVal);

    }

    return gradient; 

  }

  template <typename MatsT, typename IntsT>
  void NEOFockBuilder<MatsT,IntsT>::formFock(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& empert, bool increment,
    double xHFX)
  {
    if( this->upstream == nullptr )
      CErr("Upstream FockBuilder uninitialized in formFock!");

    // Call all upstream FockBuilders
    this->upstream->formFock(ss, empert, increment, xHFX);

    formepJ(ss, increment);

    *ss.twoeH -= 2. * *this->outMat;
    *ss.fockMatrix -= 2. * *this->outMat;
  }

  template <typename MatsT, typename IntsT>
  std::vector<double> NEOFockBuilder<MatsT,IntsT>::getGDGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& pert, double xHFX) {
    if( this->upstream == nullptr )
      CErr("Upstream FockBuilder uninitialized in getGDGrad!");

    size_t nGrad = 3*ss.molecule().nAtoms;

    std::vector<double> gradient = this->upstream->getGDGrad(ss, pert, xHFX);
    std::vector<double> epjGrad = formepJGrad(ss, pert, xHFX);

    /***
    // Add additional printout for debugging NEO-Ehrenfest
    auto printGrad = [&](std::string name, std::vector<double>& vecgrad) {
      std::cout << name << std::endl;
      std::cout << std::setprecision(12);
      for( auto iAt = 0; iAt < ss.molecule().nAtoms; iAt++ ) {
        std::cout << " Gradient@I = " << iAt << ":";
        for( auto iCart = 0; iCart < 3; iCart++ ) {
          std::cout << "  " << vecgrad[iAt*3 + iCart];
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    };
    printGrad("EPJ Gradient:", epjGrad);
    ***/

    std::transform(gradient.begin(), gradient.end(), epjGrad.begin(),
                   gradient.begin(), std::plus<double>());

    return gradient;

  };


  /*
   * NEO KS METHODS
   */


  template <typename MatsT, typename IntsT>
  void NEOKohnShamBuilder<MatsT,IntsT>::formVXC(SingleSlater<MatsT,IntsT>& ss) {

    KohnSham<MatsT,IntsT>* ks = dynamic_cast<KohnSham<MatsT,IntsT>*>(&ss);
    bool isKS(ks);

    // TODO: initialize integration parameter storage
    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    // Parallelism

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(ss.comm);
    size_t mpiSize   = MPISize(ss.comm);

    // MPI Communicator for numerical integration
    // *** Assumes that MPI will be done only across atoms in integration

    size_t nAtoms = ss.molecule().nAtoms;
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(ss.comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) 
#endif
    {
  
      // Define several useful quantities for later on
      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
  
      bool isGGA = false;
      if( isKS )
        isGGA = std::any_of(ks->functionals.begin(),ks->functionals.end(),
                     [](std::shared_ptr<DFTFunctional> &x) {
                       return x->isGGA(); 
                     }); 

      bool epcisGGA = std::any_of(epc_functionals.begin(),epc_functionals.end(),
                       [](std::shared_ptr<DFTFunctional> &x) {
                       return x->isGGA(); 
                     }); 
  
      // Turn off LA threads
      SetLAThreads(1);
  
      BasisSet &basis = ss.basisSet();
      size_t NB     = basis.nBasis;
      size_t NB2    = NB*NB;
      size_t NTNB2  = nthreads * NB2;
      size_t NTNPPB = nthreads*NPtsMaxPerBatch;

      // Auxiliary basis set information for NEO-DFT
      BasisSet &aux_basis = this->aux_ss->basisSet();
      size_t aux_NB    = aux_basis.nBasis;
      size_t aux_NB2   = aux_NB*aux_NB;
      size_t aux_NTNB2 = nthreads * aux_NB2;

      // Make sure we've allocated VXC if we need to
      if( !VXC ) {
        if( ks ) {
          VXC = ks->VXC;
        } else {
          VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
            NB, ss.nC > 1, !ss.iCS
          );
        }
      }

      // Clean up all VXC components for a the evaluation for a new 
      // batch of points
  
      std::vector<std::vector<double*>> integrateVXC;
      double* intVXC_RAW = (nthreads == 1) ? nullptr :
        CQMemManager::get().malloc<double>(VXC->nComponent() * NTNB2);
  
      std::vector<double*> VXC_SZYX = VXC->SZYXPointers(); 
      for(auto k = 0; k < VXC_SZYX.size(); k++) {
        integrateVXC.emplace_back();

        if( nthreads != 1 ) 
          for(auto ith = 0; ith < nthreads; ith++)
            integrateVXC.back().emplace_back(
              intVXC_RAW + (ith + k*nthreads) * NB2
            );

        else 
          integrateVXC.back().emplace_back(VXC_SZYX[k]);

      }
      
      for(auto &X : integrateVXC) for(auto &Y : X) std::fill_n(Y,NB2,0.);

      std::vector<double> integrateXCEnergy(nthreads,0.);
      std::vector<double> integrateEPCEnergy(nthreads,0.);
  
      // Allocating Memory
      // ----------------------------------------------------------//
      if( ks ) ks->XCEnergy = 0.;
      XCEnergy = 0.;
      double *SCRATCHNBNB = 
        CQMemManager::get().malloc<double>(NTNB2); 
      double *SCRATCHNBNP = 
        CQMemManager::get().malloc<double>(NTNPPB*NB); 

      double *DenS, *DenZ, *DenX ;
      DenS = CQMemManager::get().malloc<double>(NTNPPB);

      if( ss.onePDM->hasZ() )
        DenZ = CQMemManager::get().malloc<double>(NTNPPB);

      if( ss.onePDM->hasXY() )
        CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

      double *epsEval = CQMemManager::get().malloc<double>(NTNPPB);
      double *epcEval = CQMemManager::get().malloc<double>(NTNPPB);
      std::fill_n(epcEval, NTNPPB, 0.);

      // Density U-Variables
      double *U_n   = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 
      double *dVU_n = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 

      double *ZrhoVar1, *ZgammaVar1, *ZgammaVar2;

      ZrhoVar1 = CQMemManager::get().malloc<double>(NTNPPB);
      
      // These quantities are only used for GGA functionals
      double *GDenS, *GDenZ, *U_gamma, *dVU_gamma;
      if( isGGA || epcisGGA ) {
        GDenS = CQMemManager::get().malloc<double>(3*NTNPPB);
        ZgammaVar1 = CQMemManager::get().malloc<double>(NTNPPB);
        ZgammaVar2 = CQMemManager::get().malloc<double>(NTNPPB);

        if( ss.onePDM->hasZ() )
          GDenZ = CQMemManager::get().malloc<double>(3*NTNPPB);

        if( ss.onePDM->hasXY() )
          CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

        // Gamma U-Variables
        U_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
        dVU_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
      }

      // These quantities are used when there are multiple xc functionals
      double *epsSCR, *dVU_n_SCR, *dVU_gamma_SCR;
      if( ks && ks->functionals.size() > 1 ) {
        epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        dVU_n_SCR    = CQMemManager::get().malloc<double>(2*NTNPPB);
        if(isGGA || epcisGGA) 
          dVU_gamma_SCR = 
            CQMemManager::get().malloc<double>(3*NTNPPB);
      }

      // ZMatrix
      double *ZMAT = CQMemManager::get().malloc<double>(NTNPPB*NB);
 
      // Decide if we need to allocate space for real part of the 
      // densities and copy over the real parts
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> Re1PDM
          = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
              ss.onePDM->real_part());


      // Scratch pointers for auxiliary systems 
      // ---------------NEO------------------------------------------------
      double *AUX_SCRATCHNBNB = 
        CQMemManager::get().malloc<double>(aux_NTNB2);
      double *AUX_SCRATCHNBNP = 
        CQMemManager::get().malloc<double>(NTNPPB*aux_NB);        

      // Density pointers for auxiliary system
      double *aux_DenS, *aux_DenZ;
      aux_DenS = CQMemManager::get().malloc<double>(NTNPPB);

        if ( this->aux_ss->onePDM->hasZ() )
          aux_DenZ = CQMemManager::get().malloc<double>(NTNPPB);

        if ( this->aux_ss->onePDM->hasXY() )
          CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

      // NEO auxiliary Density U-Variables
      double *aux_U_n   = 
        CQMemManager::get().malloc<double>(2*NTNPPB);
      double *aux_dVU_n  = 
        CQMemManager::get().malloc<double>(2*NTNPPB);

      double *aux_ZrhoVar1, *aux_ZgammaVar1, *aux_ZgammaVar2;
      // These quantities are only used for GGA functionals
      double *aux_GDenS, *aux_GDenZ, *aux_U_gamma, *aux_dVU_gamma;

      aux_ZrhoVar1 = CQMemManager::get().malloc<double>(NTNPPB);

      if( epcisGGA ) {
        aux_GDenS = CQMemManager::get().malloc<double>(3*NTNPPB);
        aux_ZgammaVar1 = CQMemManager::get().malloc<double>(NTNPPB);
        aux_ZgammaVar2 = CQMemManager::get().malloc<double>(NTNPPB);

        if( this->aux_ss->onePDM->hasZ() )
          aux_GDenZ = CQMemManager::get().malloc<double>(3*NTNPPB);

        if( this->aux_ss->onePDM->hasXY() )
          CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

        // Gamma U-Variables
        aux_U_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
        aux_dVU_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
      }

      double *aux_epsSCR, *aux_dVU_n_SCR, *aux_dVU_gamma_SCR;
      if( epc_functionals.size() > 1 ) {
        aux_epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        aux_dVU_n_SCR    = CQMemManager::get().malloc<double>(2*NTNPPB);
        if(epcisGGA) 
          aux_dVU_gamma_SCR = 
            CQMemManager::get().malloc<double>(3*NTNPPB);
      }

      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> aux_Re1PDM
          = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
              this->aux_ss->onePDM->real_part());


      // ---------------end NEO------------------------------------------------

      // --------------Cross Memory--------------------------------------------
      double *cross_U_gamma, *cross_dVU_gamma, *ZgammaVar3;
      if (epcisGGA) {
        cross_U_gamma = CQMemManager::get().malloc<double>(4*NTNPPB); 
        cross_dVU_gamma = CQMemManager::get().malloc<double>(4*NTNPPB);
        ZgammaVar3 = CQMemManager::get().malloc<double>(NTNPPB);
      }
      // --------------end cross-----------------------------------------------
 

      // -------------------------------------------------------------//
      // End allocating Memory

      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
        std::vector<double> &weights, std::vector<size_t> NBE_vec, 
        std::vector<double*> BasisEval_vec, 
        std::vector<std::vector<size_t>> & batchEvalShells_vec, 
        std::vector<std::vector<std::pair<size_t,size_t>>> & subMatCut_vec) {

        // main system
        size_t NBE = NBE_vec[0];
        double * BasisEval = BasisEval_vec[0];
        std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
        std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

        // auxiliary system
        size_t aux_NBE = NBE_vec[1];
        double * aux_BasisEval = BasisEval_vec[1];
        std::vector<size_t> & aux_batchEvalShells = batchEvalShells_vec[1];
        std::vector<std::pair<size_t,size_t>> & aux_subMatCut = subMatCut_vec[1];


        // intParam.epsilon / ntotalpts (NANG * NRAD * NATOMS)
        double epsScreen = intParam.epsilon / nAtoms /
          intParam.nAng / intParam.nRad;

        epsScreen = std::max(epsScreen,
                             std::numeric_limits<double>::epsilon());

        size_t NPts = batch.size();
        size_t IOff = NBE*NPts;

        size_t thread_id = GetThreadID();
        size_t TIDNPPB   = thread_id * NPtsMaxPerBatch;

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto topevalDen = std::chrono::high_resolution_clock::now();
#endif

        // --------------Main System------------------------------------------
        // Setup local pointers
        double * SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NB2;
        double * SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NB*NPtsMaxPerBatch;

        // local density pointers 
        double * DenS_loc = DenS + TIDNPPB;
        double * DenZ_loc = DenZ + TIDNPPB;

        // local density gradient pointers
        double * GDenS_loc = GDenS + 3*TIDNPPB;
        double * GDenZ_loc = GDenZ + 3*TIDNPPB;

        // local vxc energy and U, V pointer
        double * epsEval_loc   = epsEval   +  TIDNPPB;
        double * epcEval_loc   = epcEval   +  TIDNPPB;
        double * U_n_loc       = U_n       + 2*TIDNPPB;
        double * dVU_n_loc     = dVU_n     + 2*TIDNPPB;
        double * U_gamma_loc   = U_gamma   + 3*TIDNPPB;
        double * dVU_gamma_loc = dVU_gamma + 3*TIDNPPB;

        // local gradient pointer
        double * ZrhoVar1_loc   = ZrhoVar1   + TIDNPPB;
        double * ZgammaVar1_loc = ZgammaVar1 + TIDNPPB;
        double * ZgammaVar2_loc = ZgammaVar2 + TIDNPPB;

        // local multiple functional pointer
        double * epsSCR_loc        = epsSCR        +   TIDNPPB;
        double * dVU_n_SCR_loc     = dVU_n_SCR     + 2*TIDNPPB;
        double * dVU_gamma_SCR_loc = dVU_gamma_SCR + 3*TIDNPPB;

        // local Zmat
        double *ZMAT_loc = ZMAT + NB * TIDNPPB;

        // ---------------End Main System------------------------------------


        // ---------------Auxiliary System-----------------------------------
        // aux local scratch pointers 
        double * AUX_SCRATCHNBNB_loc, * AUX_SCRATCHNBNP_loc;
        AUX_SCRATCHNBNB_loc = AUX_SCRATCHNBNB + thread_id * aux_NB2;
        AUX_SCRATCHNBNP_loc = AUX_SCRATCHNBNP + thread_id * aux_NB*NPtsMaxPerBatch;         

        // aux local density pointers
        double * aux_DenS_loc, * aux_DenZ_loc;
        aux_DenS_loc = aux_DenS + TIDNPPB;
        aux_DenZ_loc = aux_DenZ + TIDNPPB;

        // aux local density gradient pointers
        double *aux_GDenS_loc, *aux_GDenZ_loc;
        aux_GDenS_loc = aux_GDenS + 3*TIDNPPB;
        aux_GDenZ_loc = aux_GDenZ + 3*TIDNPPB;

        // aux local U and V pointers
        double * aux_U_n_loc, * aux_dVU_n_loc, * aux_U_gamma_loc, * aux_dVU_gamma_loc;
        aux_U_n_loc       = aux_U_n       + 2*TIDNPPB;
        aux_dVU_n_loc     = aux_dVU_n     + 2*TIDNPPB;
        aux_U_gamma_loc   = aux_U_gamma   + 3*TIDNPPB;
        aux_dVU_gamma_loc = aux_dVU_gamma + 3*TIDNPPB;

        // aux local Z pointers
        double * aux_ZrhoVar1_loc, * aux_ZgammaVar1_loc, * aux_ZgammaVar2_loc;        
        aux_ZrhoVar1_loc   = aux_ZrhoVar1   + TIDNPPB;
        aux_ZgammaVar1_loc = aux_ZgammaVar1 + TIDNPPB;
        aux_ZgammaVar2_loc = aux_ZgammaVar2 + TIDNPPB;

        // aux local multiple functional pointers
        double * aux_epsSCR_loc, * aux_dVU_n_SCR_loc, * aux_dVU_gamma_SCR_loc;
        aux_epsSCR_loc        = aux_epsSCR        +   TIDNPPB;
        aux_dVU_n_SCR_loc     = aux_dVU_n_SCR     + 2*TIDNPPB;
        aux_dVU_gamma_SCR_loc = aux_dVU_gamma_SCR + 3*TIDNPPB;

        // ---------------End Auxiliary System--------------------------------

        // ---------------Cross Between Main and Auxiliary system-------------
        double * cross_U_gamma_loc, *cross_dVU_gamma_loc, *ZgammaVar3_loc;
        if (epcisGGA) {
          cross_U_gamma_loc = cross_U_gamma + 4*TIDNPPB;
          cross_dVU_gamma_loc = cross_dVU_gamma + 4*TIDNPPB;
          ZgammaVar3_loc = ZgammaVar3 + TIDNPPB;
        }
        // ---------------End Cross-------------------------------------------
        
        // This evaluates the V variables for all components of the main system
        // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
        evalDen((isGGA || epcisGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc, SCRATCHNBNP_loc, Re1PDM->S().pointer(), DenS_loc, 
            GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

        // evaluate density for auxiliary components
        evalDen((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, 
          aux_NB, aux_subMatCut, 
          AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, aux_Re1PDM->S().pointer(), aux_DenS_loc,
          aux_GDenS_loc, aux_GDenS_loc + NPts, aux_GDenS_loc + 2*NPts, aux_BasisEval);


        if ( ss.onePDM->hasZ() )
          evalDen((isGGA || epcisGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
                SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Z().pointer(), DenZ_loc, 
                GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

        if ( this->aux_ss->onePDM->hasZ() )
          evalDen((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, 
            aux_NB, aux_subMatCut, 
            AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, aux_Re1PDM->Z().pointer(), 
            aux_DenZ_loc,
            aux_GDenZ_loc, aux_GDenZ_loc + NPts, aux_GDenZ_loc + 2*NPts, aux_BasisEval);

        if ( ss.onePDM->hasXY() )
          CErr("Relativistic NEO-Kohn-Sham NYI!");

        if ( this->aux_ss->onePDM->hasXY() )
          CErr("Relativistic NEO-Kohn-Sham NYI!");

        
        // V -> U variables for NEO-DFT kernal derivatives 
        mkAuxVar(this->aux_ss->onePDM,
          epcisGGA,epsScreen,NPts,
          aux_DenS_loc,aux_DenZ_loc,nullptr,nullptr,
          aux_GDenS_loc,aux_GDenS_loc + NPts,aux_GDenS_loc + 2*NPts,
          aux_GDenZ_loc,aux_GDenZ_loc + NPts,aux_GDenZ_loc + 2*NPts,
          nullptr,nullptr,nullptr,
          nullptr,nullptr,nullptr,
          nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr, 
          aux_U_n_loc, aux_U_gamma_loc
        );


        // V -> U variables for evaluating the kernel derivatives.
        mkAuxVar(ss.onePDM,
          isGGA || epcisGGA,epsScreen,NPts,
          DenS_loc,DenZ_loc,nullptr,nullptr,
          GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
          GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
          nullptr,nullptr,nullptr,
          nullptr,nullptr,nullptr,
          nullptr, 
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr, 
          U_n_loc,U_gamma_loc
        );


        // Cross V -> U variables
        if (epcisGGA)
          mkCrossAuxVar(false,ss.particle.charge < 0,
            ss.onePDM, this->aux_ss->onePDM, epsScreen,NPts,
            GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
            nullptr,nullptr,nullptr,
            nullptr,nullptr,nullptr,
            GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
            aux_GDenS_loc,aux_GDenS_loc + NPts,aux_GDenS_loc + 2*NPts,
            nullptr,nullptr,nullptr,
            nullptr,nullptr,nullptr,
            aux_GDenZ_loc,aux_GDenZ_loc + NPts,aux_GDenZ_loc + 2*NPts,
            cross_U_gamma_loc);


        // Get NEO-DFT Energy derivatives wrt U variables
        if( ks )
          loadVXCder(
            ks->functionals,
            NPts, U_n_loc, U_gamma_loc,
            epsEval_loc, dVU_n_loc, dVU_gamma_loc,
            epsSCR_loc, dVU_n_SCR_loc, dVU_gamma_SCR_loc);
        else {
          std::fill_n(epsEval_loc, NPts, 0.);
          std::fill_n(dVU_n_loc, NPts, 0.);
          if( isGGA || epcisGGA ) std::fill_n(dVU_gamma_loc, NPts, 0.);
          if( ks && ks->functionals.size() > 1 ) {
            std::fill_n(epsSCR_loc, NPts, 0.);
            std::fill_n(dVU_n_SCR_loc, NPts, 0.);
            if( isGGA || epcisGGA ) std::fill_n(dVU_gamma_SCR_loc, NPts, 0.);
          }
        }


        loadEPCVXCder(ss.particle.charge < 0, 
          epc_functionals,
          NPts, U_n_loc, U_gamma_loc, aux_U_n_loc,
          aux_U_gamma_loc, cross_U_gamma_loc, epsEval_loc, dVU_n_loc,
          dVU_gamma_loc, cross_dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc, 
          dVU_gamma_SCR_loc, cross_dVU_gamma_loc, epcEval_loc);


        // Compute for the current batch the XC energy and increment the 
        // total XC energy.
        integrateXCEnergy[thread_id] += 
          energy_vxc(NPts, weights, epsEval_loc, DenS_loc);
        //integrateEPCEnergy[thread_id] += 
        //  energy_vxc(NPts, weights, epcEval_loc, DenS_loc);


        // Construct the required quantities for the formation of the Z 
        // vector (SCALAR) given the kernel derivatives wrt U variables. 
        constructZVars(ss.onePDM, SCALAR, isGGA || epcisGGA, 
          NPts,dVU_n_loc,dVU_gamma_loc, ZrhoVar1_loc, ZgammaVar1_loc,
          ZgammaVar2_loc);

        // Construct the required quantities for the formation of the Z
        // vector for EPC-19 functional
        if (epcisGGA)
          constructEPCZVars(ss.particle.charge < 0,
            SCALAR,NPts,cross_dVU_gamma_loc, ZgammaVar3_loc);

        // Creating ZMAT (SCALAR) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        if ( epcisGGA )
          formZ_vxc_epc(ss.onePDM, this->aux_ss->onePDM,
            SCALAR, epcisGGA, NPts, NBE, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_loc, DenZ_loc, nullptr, nullptr, GDenS_loc, GDenZ_loc, nullptr,
            nullptr, aux_DenS_loc, aux_DenZ_loc, nullptr, nullptr,
            aux_GDenS_loc, aux_GDenZ_loc, nullptr, nullptr,
            BasisEval, ZMAT_loc);
        else
          formZ_vxc(ss.onePDM, SCALAR,isGGA, NPts, NBE, IOff, epsScreen, weights, 
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
            DenZ_loc, nullptr, nullptr, GDenS_loc, GDenZ_loc, nullptr, 
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            nullptr, BasisEval, ZMAT_loc);

        bool evalZ = true;

        if (evalZ) {
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat - SCALAR)
          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,1.,BasisEval,NBE,ZMAT_loc,NBE,0.,
            SCRATCHNBNB_loc,NBE);

          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[SCALAR][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);
        }


        if( not ss.onePDM->hasZ() ) return;

        //
        // ---------------   UKS or 2C ------------- Mz ---------------
        //    See J. Chem. Theory Comput. 2017, 13, 2591-2603  
        //

        // Construct the required quantities for the formation of 
        // the Z vector (Mz) given the kernel derivatives wrt U 
        // variables.
        constructZVars(ss.onePDM, MZ,isGGA || epcisGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Construct the required quantities for the formation of the Z
        // vector for EPC-19 functional
        if (epcisGGA)
          constructEPCZVars(ss.particle.charge < 0, 
            MZ,NPts,cross_dVU_gamma_loc, ZgammaVar3_loc);

        // Creating ZMAT (Mz) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        if (epcisGGA)
          formZ_vxc_epc(ss.onePDM, this->aux_ss->onePDM,
            MZ, epcisGGA, NPts, NBE, IOff, epsScreen, weights,
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc,
            DenS_loc, DenZ_loc, nullptr, nullptr, GDenS_loc, GDenZ_loc, nullptr,
            nullptr, aux_DenS_loc, aux_DenZ_loc, nullptr, nullptr,
            aux_GDenS_loc, aux_GDenZ_loc, nullptr, nullptr,
            BasisEval, ZMAT_loc);
        else
          formZ_vxc(ss.onePDM,MZ,isGGA, NPts, NBE, IOff, epsScreen, weights, 
            ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
            DenZ_loc, nullptr, nullptr, GDenS_loc, GDenZ_loc, nullptr, 
            nullptr, nullptr, nullptr, 
            nullptr, nullptr, nullptr, 
            nullptr, BasisEval, ZMAT_loc);


        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,1.,BasisEval,NBE,ZMAT_loc,NBE,0.,
            SCRATCHNBNB_loc,NBE);
    
    
          // Locating the submatrix in the right position given 
          // the subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MZ][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           

        }
 


        if( not ss.onePDM->hasXY() ) return;
        
        CErr("Relativistic NEO-Kohn-Sham NYI!",std::cout);

      }; // VXC integrate



      // Create the BeckeIntegrator object
      if( ks ) intParam = ks->intParam;
      else if( auto aux_ks = dynamic_cast<KohnSham<MatsT,IntsT>*>( this->aux_ss ) )
        intParam = aux_ks->intParam;
      BeckeIntegrator<EulerMac> 
        integrator(intComm,ss.molecule(),basis,aux_basis,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch,
          (isGGA ? GRADIENT : NOGRAD), (epcisGGA ? GRADIENT : NOGRAD), 
          intParam.epsilon);

      // Integrate the VXC
      integrator.integrate<size_t>(vxcbuild);



      // Finishing up the VXC
      // factor in the 4 pi (Lebedev) and built the upper triagolar part
      // since we create only the lower triangular. For all components
      for(auto k = 0; k < VXC_SZYX.size(); k++) {
        if( nthreads == 1 )
          blas::scal(NB2,4*M_PI,VXC_SZYX[k],1);
        else
          for(auto ithread = 0; ithread < nthreads; ithread++)
            MatAdd('N','N',NB,NB,((ithread == 0) ? 0. : 1.),VXC_SZYX[k],NB,
              4*M_PI,integrateVXC[k][ithread],NB, VXC_SZYX[k],NB);
        
        HerMat('L',NB,VXC_SZYX[k],NB);
      }

      for(auto &X : integrateXCEnergy) {
        XCEnergy += 4*M_PI*X;
      }

      if( ks )
        ks->XCEnergy = XCEnergy;
      else
        ss.extraEnergy = XCEnergy;


// Combine MPI results
#ifdef CQ_ENABLE_MPI

      double* mpiScr;
      if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<double>(NB*NB);

      for(auto &V : VXC_SZYX) {

        MPIReduce(V, NB*NB, mpiScr, 0, intComm);

        if( mpiRank == 0 ) std::copy_n(mpiScr,NB*NB,V);

      }

      if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

      if(ks) ks->XCEnergy = MPIReduce(ks->XCEnergy,0,intComm);

#endif


      // Freeing the memory
      // ----------------Main System-------------------------------------------- //
      CQMemManager::get().free(SCRATCHNBNB,SCRATCHNBNP,DenS,epsEval,epcEval,U_n,
        dVU_n, ZrhoVar1,ZMAT);
      if( isGGA || epcisGGA )  
        CQMemManager::get().free(ZgammaVar1,ZgammaVar2,GDenS,U_gamma,
          dVU_gamma);

      if( ss.onePDM->hasZ() ) {
        CQMemManager::get().free(DenZ);
        if( isGGA || epcisGGA )  CQMemManager::get().free(GDenZ);
      }

      if( ss.onePDM->hasXY() )
        CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

      if( ks && ks->functionals.size() > 1 ) {
        CQMemManager::get().free(epsSCR,dVU_n_SCR);
        if( isGGA || epcisGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }

      if( nthreads != 1 ) CQMemManager::get().free(intVXC_RAW);

      Re1PDM = nullptr;

      // -----------------End Main System-------------------------------------------- //

      // ----------------Auxiliary System-------------------------------------------- //
      CQMemManager::get().free(AUX_SCRATCHNBNB,AUX_SCRATCHNBNP,aux_DenS,aux_U_n,aux_dVU_n,aux_ZrhoVar1);
      if(epcisGGA )  
        CQMemManager::get().free(aux_ZgammaVar1,aux_ZgammaVar2,aux_GDenS,aux_U_gamma,
          aux_dVU_gamma);

      if( this->aux_ss->onePDM->hasZ() ) {
        CQMemManager::get().free(aux_DenZ);
        if ( epcisGGA ) CQMemManager::get().free(aux_GDenZ);
      }

      if( this->aux_ss->onePDM->hasXY() )
        CErr("Relativistic NEO-Kohn-Sham NYI!", std::cout);

      if( epc_functionals.size() > 1 ) {
        CQMemManager::get().free(aux_epsSCR,aux_dVU_n_SCR);
        if( epcisGGA ) CQMemManager::get().free(aux_dVU_gamma_SCR);
      }

      aux_Re1PDM = nullptr;

      // -----------------End Aux-------------------------------------------- //
      // -----------------Cross---------------------------------------------- //
      if(epcisGGA)
        CQMemManager::get().free(cross_U_gamma, cross_dVU_gamma, ZgammaVar3);
      // -----------------End Cross------------------------------------------ //
      // End freeing the memory

      // Turn back on LA threads
      SetLAThreads(LAThreads);

      MPICommFree(intComm); // Free communicator

    } // Valid intComm

    MPI_Barrier(ss.comm); // Syncronize the MPI processes

  }

  template <>
  inline void NEOKohnShamBuilder<dcomplex,dcomplex>::formVXC(SingleSlater<dcomplex,dcomplex>&) {
    CErr("GIAO for NEO NYI!");
  };


  template <typename MatsT, typename IntsT>
  void NEOKohnShamBuilder<MatsT,IntsT>::formFock(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& empert, bool increment,
    double xHFX)
  {
    if( this->upstream == nullptr )
      CErr("Upstream FockBuilder uninitialized in formFock!");

    // Call all upstream FockBuilders
    this->upstream->formFock(ss, empert, increment, xHFX);

    std::string sys = ss.particle.charge > 0 ? "Protonic" : "Electronic";
    auto vxcBegin = tick();
    if(not this->intParam.useGauXC){
      // InHouse NEO-DFT:
      formVXC(ss);
      *ss.fockMatrix += *VXC;
      //KohnSham<MatsT,IntsT>* ks = dynamic_cast<KohnSham<MatsT,IntsT>*>(&ss);
      //if(ss.particle.charge < 0) std::cout << "TOTAL Electronic EXC: " << ks->XCEnergy << std::endl;
      //if(ss.particle.charge > 0) std::cout << "TOTAL Protonic   EXC: " << ks->XCEnergy << std::endl;
    }else{
      // GauXC NEO-DFT:
      // EPC will be done only in Electronic formVXC() call to avoid evaluating rho_e and rho_p twice
      if(ss.particle.charge > 0){
        // When ss is protonic, do nothing
        return;
      } else{
        // When ss is electronic, evaluate EPC for both systems 
        bool is_uks = ss.onePDM->hasZ();
        bool is_rks = not is_uks;
        size_t elec_NB = ss.basisSet().nBasis; 
        size_t prot_NB = this->aux_ss->basisSet().nBasis;
        
        // Convert CQ matrices to be Eigen matrices to feed into GauXC
        Eigen::Matrix<double, -1, -1> elec_Ps, elec_Pz, prot_Ps, prot_Pz;
        elec_Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(ss.onePDM->real_part().S().pointer(), elec_NB, elec_NB);
        prot_Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->aux_ss->onePDM->real_part().S().pointer(), prot_NB, prot_NB);
        prot_Pz = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->aux_ss->onePDM->real_part().Z().pointer(), prot_NB, prot_NB);

        // Initialize return values
        double  elec_EXC = 0.0, prot_EXC = 0.0;
        Eigen::MatrixXd elec_VXCs, elec_VXCz, prot_VXCs, prot_VXCz;

        // Call corresonding epc evaluation functions 
        if (is_rks){ 
          // RKS Electronic + UKS Protonic
          elec_Ps /= 2.0; // Need to scale by 0.5 due to GauXC's RKS logic
          std::tie(elec_EXC, prot_EXC, elec_VXCs, prot_VXCs, prot_VXCz)
              = ss.gauxcUtils->integrator_pointer->neo_eval_exc_vxc( elec_Ps, prot_Ps, prot_Pz );
        } else{      
          // UKS Electronic + UKS Protonic
          elec_Pz  = Eigen::Map<Eigen::Matrix<double, -1, -1>>(ss.onePDM->real_part().Z().pointer(), elec_NB, elec_NB);
          std::tie(elec_EXC, prot_EXC, elec_VXCs, elec_VXCz, prot_VXCs, prot_VXCz)
              = ss.gauxcUtils->integrator_pointer->neo_eval_exc_vxc( elec_Ps, elec_Pz, prot_Ps, prot_Pz );
        }

        // Update Energy
        KohnSham<MatsT,IntsT>* elec_ks = dynamic_cast<KohnSham<MatsT,IntsT>*>( &ss );
        KohnSham<MatsT,IntsT>* prot_ks = dynamic_cast<KohnSham<MatsT,IntsT>*>( &(*this->aux_ss) );
        elec_ks->XCEnergy = elec_EXC;
        prot_ks->XCEnergy = prot_EXC;
        
        // Update Electronic VXC (with a scaling factor of 2)
        elec_VXCs *= 2.0;
        ss.fockMatrix->S() +=  elec_VXCs;
        if (is_uks)  {
          elec_VXCz *= 2.0;
          ss.fockMatrix->Z() +=  elec_VXCz;
        }

        // Update Electronic VXC (with a scaling factor of 2)
        prot_VXCs *= 2.0; prot_VXCz *= 2.0;
        this->aux_ss->fockMatrix->S() +=  prot_VXCs;         
        this->aux_ss->fockMatrix->Z() +=  prot_VXCz;         

        //std::cout << "TOTAL Protonic   EXC: " << prot_EXC << std::endl;
        //std::cout << "TOTAL Electronic EXC: " << elec_EXC << std::endl;

      } // End Electronic GauXC
    } // End GauXC DFT


    if(ss.TPI->printContractionTiming)
      if(this->intParam.useGauXC)
        std::cout << "      Total NEO VXC duration using GauXC DFT Engine: " << tock(vxcBegin) << " s \n" << std::endl;
      else
        std::cout << "      " << (ss.particle.charge > 0 ? "Protonic" : "Electronic") 
            << " VXC duration using In-house DFT Engine: " << tock(vxcBegin) << " s \n" << std::endl;
            
} //NEOKohnShamBuilder<MatsT,IntsT>::formFock


  template <typename MatsT, typename IntsT>
  std::vector<double> NEOKohnShamBuilder<MatsT,IntsT>::getGDGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& empert, double xHFX)
  {

    // Call all upstream FockBuilders
    return this->upstream->getGDGrad(ss, empert, xHFX);

  }
}
