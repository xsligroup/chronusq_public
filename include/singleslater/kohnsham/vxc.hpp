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

#include <grid/integrator.hpp>
#include <basisset/basisset_util.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blasext.hpp>

#include <dft/util.hpp>

#include <util/threads.hpp>

// VXC_DEBUG_LEVEL == 1 - Timing
// VXC_DEBUG_LEVEL == 2 - VXC/rho/gamma + Timing
// VXC_DEBUG_LEVEL == 3 - Debug 2 + Overlap + no screening
// VXC_DEBUG_LEVEL  > 3 - Debug 3 + print everthing
#ifndef VXC_DEBUG_LEVEL
#  define VXC_DEBUG_LEVEL 0
#endif

namespace ChronusQ {
  
  /**
   *  \brief assemble the VXC for all density componet over batch
   *  of points. 
   *
   *  It handles submatrix of the VXC (for a given 
   *  subset of shell to be evaluated for the provided batch) 
   *  and assemble them. 
   *
   *  This function is integrated by the BeckeIntegrator
   *  object.
   */  
  template <typename MatsT, typename IntsT>
  void KohnSham<MatsT,IntsT>::formVXC(EMPerturbation &pert) {

#if VXC_DEBUG_LEVEL >= 1
    // TIMING 
    auto topMem = std::chrono::high_resolution_clock::now();
#endif
    ProgramTimer::tick("Form VXC");

    assert( intParam.nRad % intParam.nRadPerBatch == 0 );




    // Parallelism

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(this->comm);
    size_t mpiSize   = MPISize(this->comm);

    // MPI Communicator for numerical integration
    // *** Assumes that MPI will be done only across atoms in integration

    size_t nAtoms = this->molecule().nAtoms;
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);







#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) 
#endif
    {
  
  
  
  
  
  
      // Define several useful quantities for later on
      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
  
      bool isGGA = std::any_of(functionals.begin(),functionals.end(),
                     [](std::shared_ptr<DFTFunctional> &x) {
                       return x->isGGA(); 
                     }); 
  
      // Turn off LA threads
      SetLAThreads(1);
  
      BasisSet &basis = this->basisSet();
      size_t NB     = basis.nBasis;
      size_t NB2    = NB*NB;
      size_t NTNB2  = nthreads * NB2;
      size_t NTNPPB = nthreads*NPtsMaxPerBatch;

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
  
      // Start Debug quantities
#if VXC_DEBUG_LEVEL >= 3
      // tmp VXC submat
      double *tmpS      = CQMemManager::get().malloc<double>(NB2); 
      std::fill_n(tmpS,basis.nBasis*basis.nBasis,0.);
#endif

#if VXC_DEBUG_LEVEL >= 2
      double sumrho   = 0.;
      double sumspin  = 0.;
      double sumgamma = 0.;
#endif
      // END Debug quantities
 
      // Allocating Memory
      // ----------------------------------------------------------//
    
      XCEnergy = 0.;
      double *SCRATCHNBNB = 
        CQMemManager::get().malloc<double>(NTNB2); 
      double *SCRATCHNBNP = 
        CQMemManager::get().malloc<double>(NTNPPB*NB); 

      double *DenS, *DenZ, *DenX, *DenY, *Mnorm ;
      double *KScratch;
      double *HScratch;
      bool   *Msmall;
      DenS = CQMemManager::get().malloc<double>(NTNPPB);

      if( this->onePDM->hasZ() )
        DenZ = CQMemManager::get().malloc<double>(NTNPPB);

      if( this->onePDM->hasXY() ) {
        DenY = CQMemManager::get().malloc<double>(NTNPPB);
        DenX = CQMemManager::get().malloc<double>(NTNPPB);

        Mnorm    = CQMemManager::get().malloc<double>(NTNPPB);
        KScratch = CQMemManager::get().malloc<double>(3*NTNPPB);
        Msmall   = CQMemManager::get().malloc<bool>(NTNPPB);
      }

      double *epsEval = CQMemManager::get().malloc<double>(NTNPPB);

      // Density U-Variables
      double *U_n   = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 
      double *dVU_n = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 

      double *ZrhoVar1, *ZgammaVar1, *ZgammaVar2;

      ZrhoVar1 = CQMemManager::get().malloc<double>(NTNPPB);
      
      double *GDenS, *GDenZ, *GDenX, *GDenY, *U_gamma, *dVU_gamma;
      if( isGGA ) {
        GDenS = CQMemManager::get().malloc<double>(3*NTNPPB);
        ZgammaVar1 = CQMemManager::get().malloc<double>(NTNPPB);
        ZgammaVar2 = CQMemManager::get().malloc<double>(NTNPPB);

        if( this->onePDM->hasZ() )
          GDenZ = CQMemManager::get().malloc<double>(3*NTNPPB);

        if( this->onePDM->hasXY() ) {
          GDenY = CQMemManager::get().malloc<double>(3*NTNPPB);
          GDenX = CQMemManager::get().malloc<double>(3*NTNPPB);
          HScratch = CQMemManager::get().malloc<double>(3*NTNPPB);
        }

        // Gamma U-Variables
        U_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
        dVU_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
      }

      double *epsSCR, *dVU_n_SCR, *dVU_gamma_SCR;
      if( functionals.size() > 1 ) {
        epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        dVU_n_SCR    = CQMemManager::get().malloc<double>(2*NTNPPB);
        if(isGGA) 
          dVU_gamma_SCR = 
            CQMemManager::get().malloc<double>(3*NTNPPB);
      }
 
      // ZMatrix
      double *ZMAT = CQMemManager::get().malloc<double>(NTNPPB*NB);

      // Decide if we need to allocate space for real part of the
      // densities and copy over the real parts
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> Re1PDM
          = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
              this->onePDM->real_part());


      // -------------------------------------------------------------//
      // End allocating Memory

#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto botMem = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> durevalDen(0.)  ;
      std::chrono::duration<double> durmkAuxVar(0.) ;
      std::chrono::duration<double> durloadVXCder(0.) ;
      std::chrono::duration<double> durenergy_vxc(0.) ;
      std::chrono::duration<double> durconstructZVars(0.) ;
      std::chrono::duration<double> durformZ_vxc(0.) ;
      std::chrono::duration<double> durDSYR2K(0.) ;
      std::chrono::duration<double> durIncBySubMat(0.) ;
#endif

      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
        std::vector<double> &weights, std::vector<size_t> NBE_vec, 
        std::vector<double*> BasisEval_vec, 
        std::vector<std::vector<size_t>>& batchEvalShells_vec, 
        std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec) {

#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cerr,"BASIS SCR",BasisEval,NBE,
          4*batch.size(),NBE);
#endif

        size_t NBE = NBE_vec[0];
        double* BasisEval = BasisEval_vec[0];
        std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
        std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

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

        // Setup local pointers
        double * SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NB2;
        double * SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NB*NPtsMaxPerBatch;

        double * DenS_loc = DenS + TIDNPPB;
        double * DenZ_loc = DenZ + TIDNPPB;
        double * DenY_loc = DenY + TIDNPPB;
        double * DenX_loc = DenX + TIDNPPB;

        double * GDenS_loc = GDenS + 3*TIDNPPB;
        double * GDenZ_loc = GDenZ + 3*TIDNPPB;
        double * GDenY_loc = GDenY + 3*TIDNPPB;
        double * GDenX_loc = GDenX + 3*TIDNPPB;
        

        double * epsEval_loc   = epsEval   +  TIDNPPB;
        double * U_n_loc       = U_n       + 2*TIDNPPB;
        double * dVU_n_loc     = dVU_n     + 2*TIDNPPB;
        double * U_gamma_loc   = U_gamma   + 3*TIDNPPB;
        double * dVU_gamma_loc = dVU_gamma + 3*TIDNPPB;


        double * ZrhoVar1_loc   = ZrhoVar1   + TIDNPPB;
        double * ZgammaVar1_loc = ZgammaVar1 + TIDNPPB;
        double * ZgammaVar2_loc = ZgammaVar2 + TIDNPPB;

        double * epsSCR_loc        = epsSCR        +   TIDNPPB;
        double * dVU_n_SCR_loc     = dVU_n_SCR     + 2*TIDNPPB;
        double * dVU_gamma_SCR_loc = dVU_gamma_SCR + 3*TIDNPPB;

        double *ZMAT_loc = ZMAT + NB * TIDNPPB;

        //2C
        double * Mnorm_loc    = Mnorm        +   TIDNPPB;
        double * KScratch_loc = KScratch     + 3*TIDNPPB;
        bool   * Msmall_loc   = Msmall       +   TIDNPPB;
        double * HScratch_loc = HScratch     + 3*TIDNPPB;

        // This evaluates the V variables for all components
        // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, Re1PDM->S().pointer(), DenS_loc,
          GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

#if VXC_DEBUG_LEVEL < 3
        // Coarse screen on Density
        double MaxDenS_loc = *std::max_element(DenS_loc,DenS_loc+NPts);
        if (MaxDenS_loc < epsScreen) { return; }
#endif

        if( this->onePDM->hasZ() )
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Z().pointer(), DenZ_loc,
            GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

        if( this->onePDM->hasXY() ) {
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Y().pointer(), DenY_loc,
            GDenY_loc, GDenY_loc + NPts, GDenY_loc + 2*NPts, BasisEval);
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->X().pointer(), DenX_loc,
            GDenX_loc, GDenX_loc + NPts, GDenX_loc + 2*NPts, BasisEval);
        }

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botevalDen  = std::chrono::high_resolution_clock::now();
        auto topmkAuxVar = std::chrono::high_resolution_clock::now();
#endif

        // V -> U variables for evaluating the kernel derivatives.
        mkAuxVar(this->onePDM,isGGA,epsScreen,NPts,
          DenS_loc,DenZ_loc,DenY_loc,DenX_loc,
          GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
          GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
          GDenY_loc,GDenY_loc + NPts,GDenY_loc + 2*NPts,
          GDenX_loc,GDenX_loc + NPts,GDenX_loc + 2*NPts,
          Mnorm_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          nullptr, nullptr, 
          Msmall_loc,U_n_loc,U_gamma_loc
        );

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botmkAuxVar = std::chrono::high_resolution_clock::now();
#endif


#if VXC_DEBUG_LEVEL >= 2
        assert(nthreads == 1);
        // Debug int
        for(auto iPt = 0; iPt < NPts; iPt++) { 
          sumrho  += weights[iPt] * (U_n[2*iPt] + U_n[2*iPt + 1]);
          sumspin += weights[iPt] * (U_n[2*iPt] - U_n[2*iPt + 1]);
          if(isGGA) 
            sumgamma += weights[iPt] * 
              ( U_gamma[3*iPt] + U_gamma[3*iPt+1] + U_gamma[3*iPt+2]);
        };
        // end debug
#endif
      
#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto toploadVXCder = std::chrono::high_resolution_clock::now();
#endif

        // Get DFT Energy derivatives wrt U variables
        loadVXCder(functionals, NPts, U_n_loc, U_gamma_loc, epsEval_loc,
          dVU_n_loc, dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc,
          dVU_gamma_SCR_loc); 
  
#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botloadVXCder = std::chrono::high_resolution_clock::now();
        auto topenergy_vxc = std::chrono::high_resolution_clock::now();
#endif

        // Compute for the current batch the XC energy and increment the 
        // total XC energy.
        integrateXCEnergy[thread_id] += 
          energy_vxc(NPts, weights, epsEval_loc, DenS_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botenergy_vxc     = 
          std::chrono::high_resolution_clock::now();
        auto topconstructZVars = 
          std::chrono::high_resolution_clock::now();
#endif
   
        // Construct the required quantities for the formation of the Z 
        // vector (SCALAR) given the kernel derivatives wrt U variables. 

        constructZVars(this->onePDM,SCALAR,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botconstructZVars = 
          std::chrono::high_resolution_clock::now();
        auto topformZ_vxc      = 
          std::chrono::high_resolution_clock::now();
#endif

        // Creating ZMAT (SCALAR) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, SCALAR, isGGA, NPts, NBE, IOff, epsScreen, 
          weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botformZ_vxc = std::chrono::high_resolution_clock::now();
#endif

        bool evalZ = true;

#if VXC_DEBUG_LEVEL < 3
        // Coarse screen on ZMat
        double MaxBasis = *std::max_element(BasisEval,BasisEval+IOff);
        double MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif

        if (evalZ) {

 #if VXC_DEBUG_LEVEL >= 1
          auto topDSYR2K    = std::chrono::high_resolution_clock::now();
 #endif
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat - SCALAR)
          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,1.,BasisEval,NBE,ZMAT_loc,NBE,0.,
            SCRATCHNBNB_loc,NBE);

 #if VXC_DEBUG_LEVEL >= 1
          // TIMING
          auto botDSYR2K      = std::chrono::high_resolution_clock::now();
          durDSYR2K += botDSYR2K - topDSYR2K;
          auto topIncBySubMat = std::chrono::high_resolution_clock::now();
 #endif

          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[SCALAR][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);

 #if VXC_DEBUG_LEVEL >= 1
          // TIMING
          auto botIncBySubMat = std::chrono::high_resolution_clock::now();
          durIncBySubMat += botIncBySubMat - topIncBySubMat;
 #endif
        }



#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cerr,"Basis   ",BasisEval,NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisX  ",BasisEval+NBE*NPts,NBE,
          NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisY  ",BasisEval+2*NBE*NPts,
          NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisZ  ",BasisEval+3*NBE*NPts,
          NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"ZMAT  ",ZMAT_loc,NBE,NPts,NBE);
#endif

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        durevalDen += botevalDen - topevalDen;
        durmkAuxVar += botmkAuxVar - topmkAuxVar;
        durloadVXCder += botloadVXCder - toploadVXCder;
        durenergy_vxc += botenergy_vxc - topenergy_vxc;
        durconstructZVars += botconstructZVars - topconstructZVars;
        durformZ_vxc += botformZ_vxc - topformZ_vxc;
#endif

#if VXC_DEBUG_LEVEL >= 3
        // Create Numerical Overlap
        for(auto iPt = 0; iPt < NPts; iPt++)
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,1,weights[iPt],BasisEval + iPt*NB,NB, 
            BasisEval + iPt*NB,NB, 1.,tmpS,NB);
#endif

        if( not this->onePDM->hasZ() ) return;

        //
        // ---------------   UKS or 2C ------------- Mz ---------------
        //    See J. Chem. Theory Comput. 2017, 13, 2591-2603  
        //

        // Construct the required quantities for the formation of 
        // the Z vector (Mz) given the kernel derivatives wrt U 
        // variables.
        constructZVars(this->onePDM,MZ,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (Mz) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MZ, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
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
 


        if( not this->onePDM->hasXY() ) return;

        //
        // ---------------  2C ------------- My ----------------------
        //

        // Construct the required quantities for the formation of the 
        // Z vector (Mz) given the kernel derivatives wrt U variables. 
        constructZVars(this->onePDM,MY,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (My) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MY, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,1.,BasisEval,NBE,ZMAT_loc,NBE,0.,
            SCRATCHNBNB_loc,NBE);
    
    
          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MY][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           

        }

        //
        // ---------------  2C ------------- Mx ----------------------
        //

        // Construct the required quantities for the formation of the 
        // Z vector (Mx) given the kernel derivatives wrt U variables. 
        constructZVars(this->onePDM,MX,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (Mz) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MX, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,1.,BasisEval,NBE,ZMAT_loc,NBE,0.,
            SCRATCHNBNB_loc,NBE);
    
    
          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MX][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           
        }

      }; // VXC integrate


      // Create the BeckeIntegrator object
      BeckeIntegrator<EulerMac> 
        integrator(intComm,this->molecule(),basis,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch,
          (isGGA ? GRADIENT : NOGRAD), intParam.epsilon);

      // Integrate the VXC
      integrator.integrate<size_t>(vxcbuild);


#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto toptransform    = std::chrono::high_resolution_clock::now();
#endif

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

      for(auto &X : integrateXCEnergy)
        XCEnergy += 4*M_PI*X;



// Combine MPI results
#ifdef CQ_ENABLE_MPI

      double* mpiScr;
      if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<double>(NB*NB);

      for(auto &V : VXC_SZYX) {

        MPIReduce(V,NB*NB,mpiScr,0,intComm);

        if( mpiRank == 0 ) std::copy_n(mpiScr,NB*NB,V);

      }

      if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

      XCEnergy = MPIReduce(XCEnergy,0,intComm);

#endif

#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto bottransform = std::chrono::high_resolution_clock::now();
#endif

#if VXC_DEBUG_LEVEL >= 3
      // DebugPrint
      std::cerr << std::endl;
      blas::scal(NB2,4*M_PI,tmpS,1);
      prettyPrintSmart(std::cerr,"Analytic  Overlap",
        this->aoints.overlap,NB,NB,NB);
      prettyPrintSmart(std::cerr,"Numeric  Overlap",tmpS,NB,NB,NB);
      std::cerr << std::endl;
      std::cerr << std::endl;
      for(auto i = 0; i < NB2; i++)
        tmpS[i] = std::abs(tmpS[i] - this->aoints.overlap[i]);
      std::cerr << "MAX DIFF OVERLAP = " << 
        *std::max_element(tmpS,tmpS+NB2) << std::endl;
#endif



#if VXC_DEBUG_LEVEL >= 2
      // DEBUG
      std::cerr << std::scientific << "\n";
      std::cerr << "N     electrons      = " << 4*M_PI*sumrho << "\n";
      std::cerr << "N unp electrons      = " << 4*M_PI*sumspin << "\n";
      std::cerr << "sum gamma        = " << 4*M_PI*sumgamma << "\n";
      std::cerr << "EXC              = " << XCEnergy << "\n";

      prettyPrintSmart(std::cerr,"onePDM Scalar",this->onePDM[SCALAR],
        NB,NB,NB);
      prettyPrintSmart(std::cerr,"Numerical Scalar VXC ",
        integrateVXC[SCALAR][0],NB,NB,NB);

      if( not this->iCS or this->nC > 1 ) { 
        prettyPrintSmart(std::cerr,"onePDM Mz",this->onePDM[MZ],
          NB,NB,NB);
        prettyPrintSmart(std::cerr,"Numerical Mz VXC",
          integrateVXC[MZ][0],NB,NB,NB);

        if( this->onePDM->hasXY() ) {
          prettyPrintSmart(std::cerr,"onePDM My",this->onePDM[MY],
            NB,NB,NB);
          prettyPrintSmart(std::cerr,"Numerical My VXC",
            integrateVXC[MY][0],NB,NB,NB);

          prettyPrintSmart(std::cerr,"onePDM Mx",this->onePDM[MX],
            NB,NB,NB);
          prettyPrintSmart(std::cerr,"Numerical Mx VXC",
            integrateVXC[MX][0],NB,NB,NB);
        }
      }
#endif

      // Freeing the memory
      // ------------------------------------------------------------ //
      CQMemManager::get().free(SCRATCHNBNB,SCRATCHNBNP,DenS,epsEval,U_n,
        dVU_n, ZrhoVar1,ZMAT);
      if( isGGA )  
        CQMemManager::get().free(ZgammaVar1,ZgammaVar2,GDenS,U_gamma,
          dVU_gamma);

      if( this->onePDM->hasZ() ) {
        CQMemManager::get().free(DenZ);
        if( isGGA )  CQMemManager::get().free(GDenZ);
      }

      if( this->onePDM->hasXY() ) {
        CQMemManager::get().free(DenX,DenY,Mnorm,KScratch,Msmall);
        if( isGGA )  CQMemManager::get().free(GDenX,GDenY,HScratch);
      }

      if( functionals.size() > 1 ) {
        CQMemManager::get().free(epsSCR,dVU_n_SCR);
        if( isGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }


      if( nthreads != 1 ) CQMemManager::get().free(intVXC_RAW);

      Re1PDM = nullptr;
      // ------------------------------------------------------------- //
      // End freeing the memory


#if VXC_DEBUG_LEVEL >= 1
     // TIMING
     double d_batch = this->molecule().nAtoms * 
                        intParam.nRad / intParam.nRadPerBatch;

     std::chrono::duration<double> durMem = botMem - topMem;
     std::chrono::duration<double> durtransform = 
       bottransform - toptransform;

     std::cerr << std::scientific << "\n";
     std::cerr << "Mem " << durMem.count()/d_batch << "\n";
     std::cerr << "transform " << durtransform.count()/d_batch << "\n";
     std::cerr << "evalDen " << durevalDen.count()/d_batch << "\n";
     std::cerr << "mkAuxVar " << durmkAuxVar.count()/d_batch << "\n";
     std::cerr << "loadVXCder " << durloadVXCder.count()/d_batch << "\n";
     std::cerr << "energy_vxc " << durenergy_vxc.count()/d_batch << "\n";
     std::cerr << "constructZVars " << durconstructZVars.count()/d_batch 
               << "\n";
     std::cerr << "formZ_vxc " << durformZ_vxc.count()/d_batch << "\n";
     std::cerr << "DSYR2K " << durDSYR2K.count()/d_batch << "\n";
     std::cerr << "IncBySubMat " << durIncBySubMat.count()/d_batch;
     std::cerr <<  "\n\n\n";
#endif


  
      // Turn back on LA threads
      SetLAThreads(LAThreads);

      MPICommFree(intComm); // Free communicator

    } // Valid intComm

    MPI_Barrier(this->comm); // Syncronize the MPI processes

    ProgramTimer::tock("Form VXC");

  }; // KohnSham::formVXC

// TangDD: FormVXC for GIAO
  template <>
  void KohnSham<dcomplex,dcomplex>::formVXC(EMPerturbation &pert) {

#if VXC_DEBUG_LEVEL >= 1
    // TIMING 
    auto topMem = std::chrono::high_resolution_clock::now();
#endif
    ProgramTimer::tick("Form VXC");

    assert( intParam.nRad % intParam.nRadPerBatch == 0 );




    // Parallelism

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(this->comm);
    size_t mpiSize   = MPISize(this->comm);

    // MPI Communicator for numerical integration
    // *** Assumes that MPI will be done only across atoms in integration

    size_t nAtoms = this->molecule().nAtoms;
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(this->comm,color,mpiRank);







#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) 
#endif
    {
  
  
  
  
  
  
      // Define several useful quantities for later on
      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
  
      bool isGGA = std::any_of(functionals.begin(),functionals.end(),
                     [](std::shared_ptr<DFTFunctional> &x) {
                       return x->isGGA(); 
                     }); 
  
      // Turn off LA threads
      SetLAThreads(1);
  
      BasisSet &basis = this->basisSet();
      size_t NB     = basis.nBasis;
      size_t NB2    = NB*NB;
      size_t NTNB2  = nthreads * NB2;
      size_t NTNPPB = nthreads*NPtsMaxPerBatch;

      // Clean up all VXC components for a the evaluation for a new 
      // batch of points
  
      std::vector<std::vector<dcomplex*>> integrateVXC;
      dcomplex* intVXC_RAW = (nthreads == 1) ? nullptr :
        CQMemManager::get().malloc<dcomplex>(VXC->nComponent() * NTNB2);
  
      std::vector<dcomplex*> VXC_SZYX = VXC->SZYXPointers();
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
  
      // Start Debug quantities
#if VXC_DEBUG_LEVEL >= 3
      // tmp VXC submat
      double *tmpS      = CQMemManager::get().malloc<double>(NB2); 
      std::fill_n(tmpS,basis.nBasis*basis.nBasis,0.);
#endif

#if VXC_DEBUG_LEVEL >= 2
      double sumrho   = 0.;
      double sumspin  = 0.;
      double sumgamma = 0.;
#endif
      // END Debug quantities
 
      // Allocating Memory
      // ----------------------------------------------------------//
    
      XCEnergy = 0.;
      dcomplex *SCRATCHNBNB = 
        CQMemManager::get().malloc<dcomplex>(NTNB2); 
      dcomplex *SCRATCHNBNP = 
        CQMemManager::get().malloc<dcomplex>(NTNPPB*NB); 

      double *DenS, *DenZ, *DenX, *DenY, *Mnorm ;
      double *KScratch;
      double *HScratch;
      bool   *Msmall;
      DenS = CQMemManager::get().malloc<double>(NTNPPB);

      if( this->onePDM->hasZ() )
        DenZ = CQMemManager::get().malloc<double>(NTNPPB);

      if( this->onePDM->hasXY() ) {
        DenY = CQMemManager::get().malloc<double>(NTNPPB);
        DenX = CQMemManager::get().malloc<double>(NTNPPB);

        Mnorm    = CQMemManager::get().malloc<double>(NTNPPB);
        KScratch = CQMemManager::get().malloc<double>(3*NTNPPB);
        Msmall   = CQMemManager::get().malloc<bool>(NTNPPB);
      }

      double *epsEval = CQMemManager::get().malloc<double>(NTNPPB);

      // Density U-Variables
      double *U_n   = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 
      double *dVU_n = 
        CQMemManager::get().malloc<double>(2*NTNPPB); 

      double *ZrhoVar1, *ZgammaVar1, *ZgammaVar2;

      ZrhoVar1 = CQMemManager::get().malloc<double>(NTNPPB);
      
      double *GDenS, *GDenZ, *GDenX, *GDenY, *U_gamma, *dVU_gamma;
      if( isGGA ) {
        GDenS = CQMemManager::get().malloc<double>(3*NTNPPB);
        ZgammaVar1 = CQMemManager::get().malloc<double>(NTNPPB);
        ZgammaVar2 = CQMemManager::get().malloc<double>(NTNPPB);

        if( this->onePDM->hasZ() )
          GDenZ = CQMemManager::get().malloc<double>(3*NTNPPB);

        if( this->onePDM->hasXY() ) {
          GDenY = CQMemManager::get().malloc<double>(3*NTNPPB);
          GDenX = CQMemManager::get().malloc<double>(3*NTNPPB);
          HScratch = CQMemManager::get().malloc<double>(3*NTNPPB);
        }

        // Gamma U-Variables
        U_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
        dVU_gamma = CQMemManager::get().malloc<double>(3*NTNPPB); 
      }

      double *epsSCR, *dVU_n_SCR, *dVU_gamma_SCR;
      if( functionals.size() > 1 ) {
        epsSCR    = CQMemManager::get().malloc<double>(NTNPPB);
        dVU_n_SCR    = CQMemManager::get().malloc<double>(2*NTNPPB);
        if(isGGA) 
          dVU_gamma_SCR = 
            CQMemManager::get().malloc<double>(3*NTNPPB);
      }
 
      // ZMatrix
      dcomplex *ZMAT = CQMemManager::get().malloc<dcomplex>(NTNPPB*NB);

      // Decide if we need to allocate space for real part of the
      // densities and copy over the real parts
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> Re1PDM;
      Re1PDM = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<dcomplex>>(
        this->onePDM);

/*
#ifdef DebugPrint
      prettyPrintSmart(std::cout,"Re1PDM",Re1PDM->S().pointer(),NB,NB,NB);
      if( this->onePDM->hasZ() )
        prettyPrintSmart(std::cout,"Re1PDM",Re1PDM->Z().pointer(),NB,NB,NB);
      if( this->onePDM->hasXY() ) {
        prettyPrintSmart(std::cout,"Re1PDM",Re1PDM->X().pointer(),NB,NB,NB);
        prettyPrintSmart(std::cout,"Re1PDM",Re1PDM->Y().pointer(),NB,NB,NB);
      }
#endif 
*/
      // -------------------------------------------------------------//
      // End allocating Memory

#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto botMem = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> durevalDen(0.)  ;
      std::chrono::duration<double> durmkAuxVar(0.) ;
      std::chrono::duration<double> durloadVXCder(0.) ;
      std::chrono::duration<double> durenergy_vxc(0.) ;
      std::chrono::duration<double> durconstructZVars(0.) ;
      std::chrono::duration<double> durformZ_vxc(0.) ;
      std::chrono::duration<double> durDSYR2K(0.) ;
      std::chrono::duration<double> durIncBySubMat(0.) ;
#endif

      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
        std::vector<double> &weights, std::vector<size_t> NBE_vec, 
        std::vector<dcomplex*> BasisEval_vec, 
        std::vector<std::vector<size_t>>& batchEvalShells_vec, 
        std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec) {

#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cerr,"BASIS SCR",BasisEval,NBE,
          4*batch.size(),NBE);
#endif

        size_t NBE = NBE_vec[0];
        dcomplex* BasisEval = BasisEval_vec[0];
        std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
        std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

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

        // Setup local pointers
        dcomplex * SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NB2;
        dcomplex * SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NB*NPtsMaxPerBatch;

        double * DenS_loc = DenS + TIDNPPB;
        double * DenZ_loc = DenZ + TIDNPPB;
        double * DenY_loc = DenY + TIDNPPB;
        double * DenX_loc = DenX + TIDNPPB;

        double * GDenS_loc = GDenS + 3*TIDNPPB;
        double * GDenZ_loc = GDenZ + 3*TIDNPPB;
        double * GDenY_loc = GDenY + 3*TIDNPPB;
        double * GDenX_loc = GDenX + 3*TIDNPPB;
        

        double * epsEval_loc   = epsEval   +  TIDNPPB;
        double * U_n_loc       = U_n       + 2*TIDNPPB;
        double * dVU_n_loc     = dVU_n     + 2*TIDNPPB;
        double * U_gamma_loc   = U_gamma   + 3*TIDNPPB;
        double * dVU_gamma_loc = dVU_gamma + 3*TIDNPPB;


        double * ZrhoVar1_loc   = ZrhoVar1   + TIDNPPB;
        double * ZgammaVar1_loc = ZgammaVar1 + TIDNPPB;
        double * ZgammaVar2_loc = ZgammaVar2 + TIDNPPB;

        double * epsSCR_loc        = epsSCR        +   TIDNPPB;
        double * dVU_n_SCR_loc     = dVU_n_SCR     + 2*TIDNPPB;
        double * dVU_gamma_SCR_loc = dVU_gamma_SCR + 3*TIDNPPB;

        dcomplex *ZMAT_loc = ZMAT + NB * TIDNPPB;

        //2C
        double * Mnorm_loc    = Mnorm        +   TIDNPPB;
        double * KScratch_loc = KScratch     + 3*TIDNPPB;
        bool   * Msmall_loc   = Msmall       +   TIDNPPB;
        double * HScratch_loc = HScratch     + 3*TIDNPPB;

        // This evaluates the V variables for all components
        // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, Re1PDM->S().pointer(), DenS_loc,
          GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

#if VXC_DEBUG_LEVEL < 3
        // Coarse screen on Density
        double MaxDenS_loc = *std::max_element(DenS_loc,DenS_loc+NPts);
        if (MaxDenS_loc < epsScreen) { return; }
#endif

        if( this->onePDM->hasZ() )
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Z().pointer(), DenZ_loc,
            GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

        if( this->onePDM->hasXY() ) {
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Y().pointer(), DenY_loc,
            GDenY_loc, GDenY_loc + NPts, GDenY_loc + 2*NPts, BasisEval);
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->X().pointer(), DenX_loc,
            GDenX_loc, GDenX_loc + NPts, GDenX_loc + 2*NPts, BasisEval);
        }

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botevalDen  = std::chrono::high_resolution_clock::now();
        auto topmkAuxVar = std::chrono::high_resolution_clock::now();
#endif

        // V -> U variables for evaluating the kernel derivatives.
        mkAuxVar(this->onePDM,isGGA,epsScreen,NPts,
          DenS_loc,DenZ_loc,DenY_loc,DenX_loc,
          GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
          GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
          GDenY_loc,GDenY_loc + NPts,GDenY_loc + 2*NPts,
          GDenX_loc,GDenX_loc + NPts,GDenX_loc + 2*NPts,
          Mnorm_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          nullptr, nullptr, 
          Msmall_loc,U_n_loc,U_gamma_loc
        );

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botmkAuxVar = std::chrono::high_resolution_clock::now();
#endif


#if VXC_DEBUG_LEVEL >= 2
        assert(nthreads == 1);
        // Debug int
        for(auto iPt = 0; iPt < NPts; iPt++) { 
          sumrho  += weights[iPt] * (U_n[2*iPt] + U_n[2*iPt + 1]);
          sumspin += weights[iPt] * (U_n[2*iPt] - U_n[2*iPt + 1]);
          if(isGGA) 
            sumgamma += weights[iPt] * 
              ( U_gamma[3*iPt] + U_gamma[3*iPt+1] + U_gamma[3*iPt+2]);
        };
        // end debug
#endif
      
#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto toploadVXCder = std::chrono::high_resolution_clock::now();
#endif

        // Get DFT Energy derivatives wrt U variables
        loadVXCder(functionals, NPts, U_n_loc, U_gamma_loc, epsEval_loc,
          dVU_n_loc, dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc,
          dVU_gamma_SCR_loc); 
  
#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botloadVXCder = std::chrono::high_resolution_clock::now();
        auto topenergy_vxc = std::chrono::high_resolution_clock::now();
#endif

        // Compute for the current batch the XC energy and increment the 
        // total XC energy.
        integrateXCEnergy[thread_id] += 
          energy_vxc(NPts, weights, epsEval_loc, DenS_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botenergy_vxc     = 
          std::chrono::high_resolution_clock::now();
        auto topconstructZVars = 
          std::chrono::high_resolution_clock::now();
#endif
   
        // Construct the required quantities for the formation of the Z 
        // vector (SCALAR) given the kernel derivatives wrt U variables. 

        constructZVars(this->onePDM,SCALAR,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botconstructZVars = 
          std::chrono::high_resolution_clock::now();
        auto topformZ_vxc      = 
          std::chrono::high_resolution_clock::now();
#endif

        // Creating ZMAT (SCALAR) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, SCALAR, isGGA, NPts, NBE, IOff, epsScreen, 
          weights, ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botformZ_vxc = std::chrono::high_resolution_clock::now();
#endif

        bool evalZ = true;

#if VXC_DEBUG_LEVEL < 3
        // Coarse screen on ZMat
        // dcomplex MaxBasis = *std::max_element(BasisEval,BasisEval+IOff);
        // dcomplex MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        // evalZ = ( std::abs(dcomplex(2.0) * MaxBasis * MaxZ) > epsScreen); 
#endif

        if (evalZ) {

 #if VXC_DEBUG_LEVEL >= 1
          auto topDSYR2K    = std::chrono::high_resolution_clock::now();
 #endif
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat - SCALAR)
// SS+TangDD: 
// here all the DSYR2K need to be replaced
// BasisEval ^* * ZMAT^T +ZMAT^* * BasisEval^T 
// FIXME
          //blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,dcomplex(1.),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(0.),
          //  SCRATCHNBNB_loc,NBE);
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),ZMAT_loc,NBE,BasisEval,NBE,dcomplex(0.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(1.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
            BasisEval[icount] = std::conj(BasisEval[icount]);          

 #if VXC_DEBUG_LEVEL >= 1
          // TIMING
          auto botDSYR2K      = std::chrono::high_resolution_clock::now();
          durDSYR2K += botDSYR2K - topDSYR2K;
          auto topIncBySubMat = std::chrono::high_resolution_clock::now();
 #endif

          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[SCALAR][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);

 #if VXC_DEBUG_LEVEL >= 1
          // TIMING
          auto botIncBySubMat = std::chrono::high_resolution_clock::now();
          durIncBySubMat += botIncBySubMat - topIncBySubMat;
 #endif
        }



#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cerr,"Basis   ",BasisEval,NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisX  ",BasisEval+NBE*NPts,NBE,
          NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisY  ",BasisEval+2*NBE*NPts,
          NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"BasisZ  ",BasisEval+3*NBE*NPts,
          NBE,NPts,NBE);
        prettyPrintSmart(std::cerr,"ZMAT  ",ZMAT_loc,NBE,NPts,NBE);
#endif

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        durevalDen += botevalDen - topevalDen;
        durmkAuxVar += botmkAuxVar - topmkAuxVar;
        durloadVXCder += botloadVXCder - toploadVXCder;
        durenergy_vxc += botenergy_vxc - topenergy_vxc;
        durconstructZVars += botconstructZVars - topconstructZVars;
        durformZ_vxc += botformZ_vxc - topformZ_vxc;
#endif

#if VXC_DEBUG_LEVEL >= 3
        // Create Numerical Overlap
        for(auto iPt = 0; iPt < NPts; iPt++)
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,1,weights[iPt],BasisEval + iPt*NB,NB, 
            BasisEval + iPt*NB,NB, 1.,tmpS,NB);
#endif

        if( not this->onePDM->hasZ() ) return;

        //
        // ---------------   UKS or 2C ------------- Mz ---------------
        //    See J. Chem. Theory Comput. 2017, 13, 2591-2603  
        //

        // Construct the required quantities for the formation of 
        // the Z vector (Mz) given the kernel derivatives wrt U 
        // variables.
        constructZVars(this->onePDM,MZ,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (Mz) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MZ, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        //MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        //evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
// SS+TangDD: 
// here all the DSYR2K need to be replaced
// BasisEval ^* * ZMAT^T +ZMAT^* * BasisEval^T 
// FIXME
          //blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,dcomplex(1.),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(0.),
          //  SCRATCHNBNB_loc,NBE);
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),ZMAT_loc,NBE,BasisEval,NBE,dcomplex(0.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(1.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
            BasisEval[icount] = std::conj(BasisEval[icount]);      
    
          // Locating the submatrix in the right position given 
          // the subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MZ][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           

        }
 


        if( not this->onePDM->hasXY() ) return;

        //
        // ---------------  2C ------------- My ----------------------
        //

        // Construct the required quantities for the formation of the 
        // Z vector (Mz) given the kernel derivatives wrt U variables. 
        constructZVars(this->onePDM,MY,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (My) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MY, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        // MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        // evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
          //blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,dcomplex(1.),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(0.),
          //  SCRATCHNBNB_loc,NBE);
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),ZMAT_loc,NBE,BasisEval,NBE,dcomplex(0.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(1.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
            BasisEval[icount] = std::conj(BasisEval[icount]);  

    
          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MY][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           

        }

        //
        // ---------------  2C ------------- Mx ----------------------
        //

        // Construct the required quantities for the formation of the 
        // Z vector (Mx) given the kernel derivatives wrt U variables. 
        constructZVars(this->onePDM,MX,isGGA,NPts,dVU_n_loc,dVU_gamma_loc,ZrhoVar1_loc,
          ZgammaVar1_loc, ZgammaVar2_loc);

        // Creating ZMAT (Mz) according to 
        //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 15 
        formZ_vxc(this->onePDM, MX, isGGA, NPts, NBE, IOff, epsScreen, weights, 
          ZrhoVar1_loc, ZgammaVar1_loc, ZgammaVar2_loc, DenS_loc, 
          DenZ_loc, DenY_loc, DenX_loc, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, KScratch_loc, KScratch_loc + NPts, 
          KScratch_loc + 2* NPts, HScratch_loc, HScratch_loc + NPts, 
          HScratch_loc + 2* NPts, BasisEval, ZMAT_loc);


#if VXC_DEBUG_LEVEL < 3
        // MaxZ     = *std::max_element(ZMAT_loc,ZMAT_loc+IOff);
        // evalZ = ( std::abs(2 * MaxBasis * MaxZ) > epsScreen); 
#endif
        // Coarse screen on ZMat
        if(evalZ) {
  
          // Creating according to 
          //   J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 14 
          //
          // Z -> VXC (submat)
          //blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,dcomplex(1.),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(0.),
          //  SCRATCHNBNB_loc,NBE);
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),ZMAT_loc,NBE,BasisEval,NBE,dcomplex(0.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            ZMAT_loc[icount] = std::conj(ZMAT_loc[icount]);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,
            dcomplex(1.0),BasisEval,NBE,ZMAT_loc,NBE,dcomplex(1.0),SCRATCHNBNB_loc,NBE);

          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
            BasisEval[icount] = std::conj(BasisEval[icount]);  

    
          // Locating the submatrix in the right position given the 
          // subset of shells for the given batch.
          IncBySubMat(NB,NB,NBE,NBE,integrateVXC[MX][thread_id],NB,
            SCRATCHNBNB_loc,NBE,subMatCut);           
        }

      }; // VXC integrate


      // Create the BeckeIntegrator object
      BeckeIntegrator<EulerMac> 
        integrator(intComm,this->molecule(),basis,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch,
          (isGGA ? GRADIENT : NOGRAD), intParam.epsilon);

      // Integrate the VXC
      integrator.integrate<size_t>(vxcbuild,pert);


#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto toptransform    = std::chrono::high_resolution_clock::now();
#endif

      // Finishing up the VXC
      // factor in the 4 pi (Lebedev) and built the upper triagolar part
      // since we create only the lower triangular. For all components
      for(auto k = 0; k < VXC_SZYX.size(); k++) {
        if( nthreads == 1 )
          blas::scal(NB2,4*M_PI,VXC_SZYX[k],1);
        else
          for(auto ithread = 0; ithread < nthreads; ithread++)
            MatAdd('N','N',NB,NB,dcomplex((ithread == 0) ? 0. : 1.),VXC_SZYX[k],NB,
              dcomplex(4*M_PI),integrateVXC[k][ithread],NB, VXC_SZYX[k],NB);

        // TangDD: HerMat only for GTO        
        //  HerMat('L',NB,VXC_SZYX[k],NB);
      }

      for(auto &X : integrateXCEnergy)
        XCEnergy += 4*M_PI*X;



// Combine MPI results
#ifdef CQ_ENABLE_MPI

      dcomplex* mpiScr;
      if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<dcomplex>(NB*NB);

      for(auto &V : VXC_SZYX) {

        MPIReduce(V,NB*NB,mpiScr,0,intComm);

        if( mpiRank == 0 ) std::copy_n(mpiScr,NB*NB,V);

      }

      if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

      XCEnergy = MPIReduce(XCEnergy,0,intComm);

#endif

#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto bottransform = std::chrono::high_resolution_clock::now();
#endif

#if VXC_DEBUG_LEVEL >= 3
      // DebugPrint
      std::cerr << std::endl;
      blas::scal(NB2,4*M_PI,tmpS,1);
      prettyPrintSmart(std::cerr,"Analytic  Overlap",
        this->aoints.overlap,NB,NB,NB);
      prettyPrintSmart(std::cerr,"Numeric  Overlap",tmpS,NB,NB,NB);
      std::cerr << std::endl;
      std::cerr << std::endl;
      for(auto i = 0; i < NB2; i++)
        tmpS[i] = std::abs(tmpS[i] - this->aoints.overlap[i]);
      std::cerr << "MAX DIFF OVERLAP = " << 
        *std::max_element(tmpS,tmpS+NB2) << std::endl;
#endif



#if VXC_DEBUG_LEVEL >= 2
      // DEBUG
      std::cerr << std::scientific << "\n";
      std::cerr << "N     electrons      = " << 4*M_PI*sumrho << "\n";
      std::cerr << "N unp electrons      = " << 4*M_PI*sumspin << "\n";
      std::cerr << "sum gamma        = " << 4*M_PI*sumgamma << "\n";
      std::cerr << "EXC              = " << XCEnergy << "\n";

      prettyPrintSmart(std::cerr,"onePDM Scalar",this->onePDM[SCALAR],
        NB,NB,NB);
      prettyPrintSmart(std::cerr,"Numerical Scalar VXC ",
        integrateVXC[SCALAR][0],NB,NB,NB);

      if( not this->iCS or this->nC > 1 ) { 
        prettyPrintSmart(std::cerr,"onePDM Mz",this->onePDM[MZ],
          NB,NB,NB);
        prettyPrintSmart(std::cerr,"Numerical Mz VXC",
          integrateVXC[MZ][0],NB,NB,NB);

        if( this->onePDM->hasXY() ) {
          prettyPrintSmart(std::cerr,"onePDM My",this->onePDM[MY],
            NB,NB,NB);
          prettyPrintSmart(std::cerr,"Numerical My VXC",
            integrateVXC[MY][0],NB,NB,NB);

          prettyPrintSmart(std::cerr,"onePDM Mx",this->onePDM[MX],
            NB,NB,NB);
          prettyPrintSmart(std::cerr,"Numerical Mx VXC",
            integrateVXC[MX][0],NB,NB,NB);
        }
      }
#endif

      // Freeing the memory
      // ------------------------------------------------------------ //
      CQMemManager::get().free(SCRATCHNBNB,SCRATCHNBNP,DenS,epsEval,U_n,
        dVU_n, ZrhoVar1,ZMAT);
      if( isGGA )  
        CQMemManager::get().free(ZgammaVar1,ZgammaVar2,GDenS,U_gamma,
          dVU_gamma);

      if( this->onePDM->hasZ() ) {
        CQMemManager::get().free(DenZ);
        if( isGGA )  CQMemManager::get().free(GDenZ);
      }

      if( this->onePDM->hasXY() ) {
        CQMemManager::get().free(DenX,DenY,Mnorm,KScratch,Msmall);
        if( isGGA )  CQMemManager::get().free(GDenX,GDenY,HScratch);
      }

      if( functionals.size() > 1 ) {
        CQMemManager::get().free(epsSCR,dVU_n_SCR);
        if( isGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }


      if( nthreads != 1 ) CQMemManager::get().free(intVXC_RAW);

      Re1PDM = nullptr;
      // ------------------------------------------------------------- //
      // End freeing the memory


#if VXC_DEBUG_LEVEL >= 1
     // TIMING
     double d_batch = this->molecule().nAtoms * 
                        intParam.nRad / intParam.nRadPerBatch;

     std::chrono::duration<double> durMem = botMem - topMem;
     std::chrono::duration<double> durtransform = 
       bottransform - toptransform;

     std::cerr << std::scientific << "\n";
     std::cerr << "Mem " << durMem.count()/d_batch << "\n";
     std::cerr << "transform " << durtransform.count()/d_batch << "\n";
     std::cerr << "evalDen " << durevalDen.count()/d_batch << "\n";
     std::cerr << "mkAuxVar " << durmkAuxVar.count()/d_batch << "\n";
     std::cerr << "loadVXCder " << durloadVXCder.count()/d_batch << "\n";
     std::cerr << "energy_vxc " << durenergy_vxc.count()/d_batch << "\n";
     std::cerr << "constructZVars " << durconstructZVars.count()/d_batch 
               << "\n";
     std::cerr << "formZ_vxc " << durformZ_vxc.count()/d_batch << "\n";
     std::cerr << "DSYR2K " << durDSYR2K.count()/d_batch << "\n";
     std::cerr << "IncBySubMat " << durIncBySubMat.count()/d_batch;
     std::cerr <<  "\n\n\n";
#endif


  
      // Turn back on LA threads
      SetLAThreads(LAThreads);

      MPICommFree(intComm); // Free communicator

    } // Valid intComm

    MPI_Barrier(this->comm); // Syncronize the MPI processes

    ProgramTimer::tock("Form VXC");

  }; // KohnSham::formVXC GIAO

}; // namespace ChronusQ

