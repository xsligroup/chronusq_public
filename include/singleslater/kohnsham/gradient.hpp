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
   *  \brief assemble the nuclear gradient of EXC for all density componet over batch
   *  of points. 
   *
   *
   *  This function is integrated by the BeckeIntegrator
   *  object.
   */  
  template <typename MatsT, typename IntsT>
  void KohnSham<MatsT,IntsT>::formEXCGradient() {
#if VXC_DEBUG_LEVEL >= 1
    // TIMING 
    auto topMem = std::chrono::high_resolution_clock::now();
#endif
    ProgramTimer::tick("Form EXC Gradient");

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
      size_t nAtoms = this->molecule_.atoms.size();

      // integrated XC Energy per thread
      std::vector<double> integrateXCEnergy(nthreads,0.);

      // integrated XC Energy Gradient per thread
      std::vector<std::vector<double>> integrateXCEnergyGradX(nthreads);
      std::vector<std::vector<double>> integrateXCEnergyGradY(nthreads);
      std::vector<std::vector<double>> integrateXCEnergyGradZ(nthreads);
      for(size_t it = 0; it < nthreads; it++) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          integrateXCEnergyGradX[it].push_back(0.);
          integrateXCEnergyGradY[it].push_back(0.);
          integrateXCEnergyGradZ[it].push_back(0.);
        }
      }
 
      // Allocating Memory
      // ----------------------------------------------------------//
    
      XCEnergy = 0.;
      double *SCRATCHNBNB = 
        CQMemManager::get().template malloc<double>(NTNB2); 
      double *SCRATCHNBNP = 
        CQMemManager::get().template malloc<double>(NTNPPB*NB); 

      std::vector<std::vector<double*>> VEC_SCRATCHNBNB;
      for (size_t ic = 0; ic < nAtoms; ic++) {
        std::vector<double*> tmp = 
        {CQMemManager::get().template malloc<double>(NTNB2), CQMemManager::get().template malloc<double>(NTNB2), CQMemManager::get().template malloc<double>(NTNB2)};
        VEC_SCRATCHNBNB.emplace_back(tmp);
      }

      std::vector<std::vector<double*>> VEC_SCRATCHNBNP;
      for (size_t ic = 0; ic < nAtoms; ic++) {
        std::vector<double*> tmp = 
        {CQMemManager::get().template malloc<double>(NTNPPB*NB), CQMemManager::get().template malloc<double>(NTNPPB*NB), CQMemManager::get().template malloc<double>(NTNPPB*NB)};
        VEC_SCRATCHNBNP.emplace_back(tmp);
      }

      std::vector<double*> VEC_SCRATCHNBNP2 = 
        {CQMemManager::get().template malloc<double>(NTNPPB*NB), CQMemManager::get().template malloc<double>(NTNPPB*NB), CQMemManager::get().template malloc<double>(NTNPPB*NB)};


      double *DenS, *DenZ;
      std::vector<double*> nGDenS_dX, nGDenS_dY, nGDenS_dZ;
      std::vector<double*> nGDenZ_dX, nGDenZ_dY, nGDenZ_dZ;
      DenS = CQMemManager::get().template malloc<double>(NTNPPB);

      for (size_t ic = 0; ic < nAtoms; ic++) {
        nGDenS_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        nGDenS_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        nGDenS_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
      }


      if( this->onePDM->hasZ() ) {
        DenZ = CQMemManager::get().template malloc<double>(NTNPPB);

        for (size_t ic = 0; ic < nAtoms; ic++) {
          nGDenZ_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          nGDenZ_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          nGDenZ_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        }
      }

      if( this->onePDM->hasXY() )
        CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

      double *epsEval = CQMemManager::get().template malloc<double>(NTNPPB);

      // Density U-Variables
      double *U_n   = 
        CQMemManager::get().template malloc<double>(2*NTNPPB); 
      double *dVU_n = 
        CQMemManager::get().template malloc<double>(2*NTNPPB);

      // Density gradient U-Variables
      std::vector<double*> GradU_nX, GradU_nY, GradU_nZ;
      for(size_t ic = 0; ic < nAtoms; ic++) {
        GradU_nX.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
        GradU_nY.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
        GradU_nZ.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
      }

      // These quantities are only used for GGA functionals
      double *GDenS, *GDenZ, *U_gamma, *dVU_gamma;
      std::vector<double*> GGDenS_dxX, GGDenS_dxY, GGDenS_dxZ;
      std::vector<double*> GGDenS_dyX, GGDenS_dyY, GGDenS_dyZ;
      std::vector<double*> GGDenS_dzX, GGDenS_dzY, GGDenS_dzZ;
      std::vector<double*> GGDenZ_dxX, GGDenZ_dxY, GGDenZ_dxZ;
      std::vector<double*> GGDenZ_dyX, GGDenZ_dyY, GGDenZ_dyZ;
      std::vector<double*> GGDenZ_dzX, GGDenZ_dzY, GGDenZ_dzZ;
      std::vector<double*> GradU_gamma_X, GradU_gamma_Y, GradU_gamma_Z;
      if( isGGA ) {
        GDenS = CQMemManager::get().template malloc<double>(3*NTNPPB);
        for (size_t ic = 0; ic < nAtoms; ic++) {
          GGDenS_dxX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dxY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dxZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dyX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dyY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dyZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dzX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dzY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          GGDenS_dzZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        }
        if( this->onePDM->hasZ() ) {
          GDenZ = CQMemManager::get().template malloc<double>(3*NTNPPB);
          for (size_t ic = 0; ic < nAtoms; ic++) {
            GGDenZ_dxX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dxY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dxZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dyX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dyY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dyZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dzX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dzY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            GGDenZ_dzZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          }
        }
        if( this->onePDM->hasXY() )
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

        // Gamma U-Variables
        U_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 
        dVU_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 

        // Gradients of Gamma U-Variables
        for(size_t ic = 0; ic < nAtoms; ic++) {
          GradU_gamma_X.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
          GradU_gamma_Y.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
          GradU_gamma_Z.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
        }
      }

      double *epsSCR, *dVU_n_SCR, *dVU_gamma_SCR;
      if( functionals.size() > 1 ) {
        epsSCR    = CQMemManager::get().template malloc<double>(NTNPPB);
        dVU_n_SCR = CQMemManager::get().template malloc<double>(2*NTNPPB);
        if (isGGA)
          dVU_gamma_SCR = 
            CQMemManager::get().template malloc<double>(3*NTNPPB);
      }

      // Decide if we need to allocate space for real part of the 
      // densities and copy over the real parts
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> Re1PDM;
      if (std::is_same<MatsT,double>::value)
        Re1PDM = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<double>>(
            this->onePDM);
      else
        Re1PDM = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
            this->onePDM->real_part());


      // Obtain the 1PDM Gradient (stored in SingleSlater), which is calculated during Pulay Gradient
      std::vector<std::vector<double*>> Re1PDMGrad_scalar(nAtoms); 
      std::vector<std::vector<double*>> Re1PDMGrad_mz(nAtoms); 
      for(auto ic = 0; ic < nAtoms; ic++) {
        for (auto xyz = 0; xyz < 3; xyz++ )//{
          if( std::is_same<MatsT,double>::value )
            Re1PDMGrad_scalar[ic].push_back(reinterpret_cast<double*>(this->onePDMGrad[3*ic+xyz]->S().pointer()));
          else{
            Re1PDMGrad_scalar[ic].push_back(CQMemManager::get().template malloc<double>(NB2));
            GetMatRE('N',NB,NB,1.,this->onePDMGrad[3*ic+xyz]->S().pointer(),NB,Re1PDMGrad_scalar[ic].back(),NB);
          }
          //prettyPrintSmart(std::cout,"1PDMGRAD DFT "+std::to_string(3*ic+xyz),Re1PDMGrad_scalar[ic][xyz],NB,NB,NB);}
      }
      
      if( this->onePDM->hasZ() ) {
        for(auto ic = 0; ic < nAtoms; ic++) {
          for(auto xyz = 0; xyz < 3; xyz++)
            if( std::is_same<MatsT,double>::value )
              Re1PDMGrad_mz[ic].push_back(reinterpret_cast<double*>(this->onePDMGrad[3*ic+xyz]->Z().pointer()));
            else{
              Re1PDMGrad_mz[ic].push_back(CQMemManager::get().template malloc<double>(NB2));
              GetMatRE('N',NB,NB,1.,this->onePDMGrad[3*ic+xyz]->Z().pointer(),NB,Re1PDMGrad_mz[ic].back(),NB);
            }
        }
      }

      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
        std::vector<double> &weights, std::vector<size_t> NBE_vec, 
        std::vector<double*> BasisEval_vec, 
        std::vector<std::vector<size_t>> & batchEvalShells_vec, 
        std::vector<std::vector<std::pair<size_t,size_t>>> & subMatCut_vec) {

#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cout,"BASIS SCR",BasisEval,NBE,
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

        size_t NDer = isGGA ? 4 : 1;
        double *BasisGradEval = BasisEval + NDer * NPts * NBE;

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto topevalDen = std::chrono::high_resolution_clock::now();
#endif

        // --------------Main System------------------------------------------
        // Setup local pointers
        double * SCRATCHNBNB_loc = SCRATCHNBNB + thread_id * NB2;
        double * SCRATCHNBNP_loc = SCRATCHNBNP + thread_id * NB*NPtsMaxPerBatch;

        // local vector scratch pointers
        std::vector<std::vector<double*>> VEC_SCRATCHNBNB_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          std::vector<double*> tmp = 
            {VEC_SCRATCHNBNB[ic][0] + thread_id * NB2, VEC_SCRATCHNBNB[ic][1] + thread_id * NB2, VEC_SCRATCHNBNB[ic][2] + thread_id * NB2};
          VEC_SCRATCHNBNB_loc.emplace_back(tmp);
        }

        std::vector<std::vector<double*>> VEC_SCRATCHNBNP_loc;
        for (size_t ic = 0; ic < nAtoms; ic ++) {
          std::vector<double*> tmp = 
            {VEC_SCRATCHNBNP[ic][0] + thread_id * NB*NPtsMaxPerBatch, VEC_SCRATCHNBNP[ic][1] + thread_id * NB*NPtsMaxPerBatch, VEC_SCRATCHNBNP[ic][2] + thread_id * NB*NPtsMaxPerBatch};
          VEC_SCRATCHNBNP_loc.emplace_back(tmp);
        }

        std::vector<double*> VEC_SCRATCHNBNP2_loc = 
          {VEC_SCRATCHNBNP2[0] + thread_id * NB*NPtsMaxPerBatch, VEC_SCRATCHNBNP2[1] + thread_id * NB*NPtsMaxPerBatch, VEC_SCRATCHNBNP2[2] + thread_id * NB*NPtsMaxPerBatch};

        // local density pointers 
        double * DenS_loc = DenS + TIDNPPB;
        double * DenZ_loc = DenZ + TIDNPPB;

        // local density nuclear gradient pointers
        std::vector<double*> nGDenS_dX_loc, nGDenS_dY_loc, nGDenS_dZ_loc;
        std::vector<double*> nGDenZ_dX_loc, nGDenZ_dY_loc, nGDenZ_dZ_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          nGDenS_dX_loc.emplace_back(nGDenS_dX[ic] + TIDNPPB);
          nGDenS_dY_loc.emplace_back(nGDenS_dY[ic] + TIDNPPB);
          nGDenS_dZ_loc.emplace_back(nGDenS_dZ[ic] + TIDNPPB);

          if (this->onePDM->hasZ() ) {
            nGDenZ_dX_loc.emplace_back(nGDenZ_dX[ic] + TIDNPPB);
            nGDenZ_dY_loc.emplace_back(nGDenZ_dY[ic] + TIDNPPB);
            nGDenZ_dZ_loc.emplace_back(nGDenZ_dZ[ic] + TIDNPPB);
          }
        }

        // local density gradient pointers
        double * GDenS_loc = GDenS + 3*TIDNPPB;
        double * GDenZ_loc = GDenZ + 3*TIDNPPB;

        // local density gradient nuclear gradient pointers
        std::vector<double*> GGDenS_dxX_loc, GGDenS_dxY_loc, GGDenS_dxZ_loc;
        std::vector<double*> GGDenS_dyX_loc, GGDenS_dyY_loc, GGDenS_dyZ_loc;
        std::vector<double*> GGDenS_dzX_loc, GGDenS_dzY_loc, GGDenS_dzZ_loc;
        std::vector<double*> GGDenZ_dxX_loc, GGDenZ_dxY_loc, GGDenZ_dxZ_loc;
        std::vector<double*> GGDenZ_dyX_loc, GGDenZ_dyY_loc, GGDenZ_dyZ_loc;
        std::vector<double*> GGDenZ_dzX_loc, GGDenZ_dzY_loc, GGDenZ_dzZ_loc;
        if ( isGGA ) {
          for (size_t ic = 0; ic < nAtoms; ic++) {
            GGDenS_dxX_loc.emplace_back(GGDenS_dxX[ic] + TIDNPPB);
            GGDenS_dxY_loc.emplace_back(GGDenS_dxY[ic] + TIDNPPB);
            GGDenS_dxZ_loc.emplace_back(GGDenS_dxZ[ic] + TIDNPPB);

            GGDenS_dyX_loc.emplace_back(GGDenS_dyX[ic] + TIDNPPB);
            GGDenS_dyY_loc.emplace_back(GGDenS_dyY[ic] + TIDNPPB);
            GGDenS_dyZ_loc.emplace_back(GGDenS_dyZ[ic] + TIDNPPB);

            GGDenS_dzX_loc.emplace_back(GGDenS_dzX[ic] + TIDNPPB);
            GGDenS_dzY_loc.emplace_back(GGDenS_dzY[ic] + TIDNPPB);
            GGDenS_dzZ_loc.emplace_back(GGDenS_dzZ[ic] + TIDNPPB);

            if (this->onePDM->hasZ() ) {

              GGDenZ_dxX_loc.emplace_back(GGDenZ_dxX[ic] + TIDNPPB);
              GGDenZ_dxY_loc.emplace_back(GGDenZ_dxY[ic] + TIDNPPB);
              GGDenZ_dxZ_loc.emplace_back(GGDenZ_dxZ[ic] + TIDNPPB);

              GGDenZ_dyX_loc.emplace_back(GGDenZ_dyX[ic] + TIDNPPB);
              GGDenZ_dyY_loc.emplace_back(GGDenZ_dyY[ic] + TIDNPPB);
              GGDenZ_dyZ_loc.emplace_back(GGDenZ_dyZ[ic] + TIDNPPB);

              GGDenZ_dzX_loc.emplace_back(GGDenZ_dzX[ic] + TIDNPPB);
              GGDenZ_dzY_loc.emplace_back(GGDenZ_dzY[ic] + TIDNPPB);
              GGDenZ_dzZ_loc.emplace_back(GGDenZ_dzZ[ic] + TIDNPPB);
            }
          }
        }

        // local vxc energy and U, V pointer
        double * epsEval_loc   = epsEval   +  TIDNPPB;
        double * U_n_loc       = U_n       + 2*TIDNPPB;
        double * dVU_n_loc     = dVU_n     + 2*TIDNPPB;
        double * U_gamma_loc   = U_gamma   + 3*TIDNPPB;
        double * dVU_gamma_loc = dVU_gamma + 3*TIDNPPB;

        // local U gradient ponter
        std::vector<double*> GradU_nX_loc, GradU_nY_loc, GradU_nZ_loc;
        std::vector<double*> GradU_gamma_X_loc, GradU_gamma_Y_loc, GradU_gamma_Z_loc;
        for(size_t ic = 0; ic < nAtoms; ic++) {
          GradU_nX_loc.emplace_back(GradU_nX[ic] + 2*TIDNPPB);
          GradU_nY_loc.emplace_back(GradU_nY[ic] + 2*TIDNPPB);
          GradU_nZ_loc.emplace_back(GradU_nZ[ic] + 2*TIDNPPB);
          if (isGGA) {
            GradU_gamma_X_loc.emplace_back(GradU_gamma_X[ic] + 3*TIDNPPB);
            GradU_gamma_Y_loc.emplace_back(GradU_gamma_Y[ic] + 3*TIDNPPB);
            GradU_gamma_Z_loc.emplace_back(GradU_gamma_Z[ic] + 3*TIDNPPB);
          }
        }
  
        // local multiple functional pointer
        double * epsSCR_loc        = epsSCR        +   TIDNPPB;
        double * dVU_n_SCR_loc     = dVU_n_SCR     + 2*TIDNPPB;
        double * dVU_gamma_SCR_loc = dVU_gamma_SCR + 3*TIDNPPB;

        // ---------------End Main System------------------------------------
        

        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, Re1PDM->S().pointer(), DenS_loc, 
          GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

        // This evaluates the Gradient of V Variables
        evalDenGrad((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut,
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, VEC_SCRATCHNBNB_loc, VEC_SCRATCHNBNP_loc,
          VEC_SCRATCHNBNP2_loc, Re1PDMGrad_scalar, Re1PDM->S().pointer(), nGDenS_dX_loc,
          nGDenS_dY_loc, nGDenS_dZ_loc, GGDenS_dxX_loc, GGDenS_dxY_loc, 
          GGDenS_dxZ_loc, GGDenS_dyX_loc, GGDenS_dyY_loc,
          GGDenS_dyZ_loc, GGDenS_dzX_loc, GGDenS_dzY_loc, 
          GGDenS_dzZ_loc, BasisEval, BasisGradEval, nAtoms, basis);

          //std::cout << "after evalDen " << std::endl;
          //for(auto iPt = 0; iPt < NPts; iPt++) {
          //    std::cout << std::fixed << std::setprecision(16) << "DenS: " << DenS_loc[iPt] << std::endl;
          //}

          //std::cout << "after evalDenGrad " << std::endl;
          //for(auto iPt = 0; iPt < NPts; iPt++) {
          //  for(auto ic = 0; ic < nAtoms; ic++) {
          //    std::cout << std::fixed << std::setprecision(16) << "GDenX: " << nGDenS_dX_loc[ic][iPt] << std::endl;
          //    std::cout << std::fixed << std::setprecision(16) << "GDenY: " << nGDenS_dY_loc[ic][iPt] << std::endl;
          //    std::cout << std::fixed << std::setprecision(16) << "GDenZ: " << nGDenS_dZ_loc[ic][iPt] << std::endl;
          //    }
          //}

#if VXC_DEBUG_LEVEL < 3
        // Coarse screen on Density
        double MaxDenS_loc = *std::max_element(DenS_loc,DenS_loc+NPts);
        //if (MaxDenS_loc < epsScreen) { return; }
#endif

        if( this->onePDM->hasZ() ) {
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Z().pointer(), DenZ_loc, 
            GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

          evalDenGrad((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut,
            SCRATCHNBNB_loc, SCRATCHNBNP_loc, VEC_SCRATCHNBNB_loc, VEC_SCRATCHNBNP_loc,
            VEC_SCRATCHNBNP2_loc, Re1PDMGrad_mz, Re1PDM->Z().pointer(), nGDenZ_dX_loc,
            nGDenZ_dY_loc, nGDenZ_dZ_loc, GGDenZ_dxX_loc, GGDenZ_dxY_loc, 
            GGDenZ_dxZ_loc, GGDenZ_dyX_loc, GGDenZ_dyY_loc,
            GGDenZ_dyZ_loc, GGDenZ_dzX_loc, GGDenZ_dzY_loc, 
            GGDenZ_dzZ_loc, BasisEval, BasisGradEval, nAtoms, basis);
        }

        if( this->onePDM->hasXY() )
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botevalDen  = std::chrono::high_resolution_clock::now();
        auto topmkAuxVar = std::chrono::high_resolution_clock::now();
#endif

        // V -> U variables for evaluating the kernel derivatives, not needed for energy only evaluation
        mkAuxVar(this->onePDM,isGGA,epsScreen,NPts,
          DenS_loc,DenZ_loc,nullptr,nullptr,
          GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
          GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
          nullptr,nullptr,nullptr,
          nullptr,nullptr,nullptr,
          nullptr, 
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr,
          nullptr, U_n_loc,U_gamma_loc
        );

        //grad V -> grad U variables for evaluateing the kernel gradients
        mkAuxVarGrad(this->onePDM,isGGA,epsScreen,NPts,
          nGDenS_dX_loc, nGDenS_dY_loc, nGDenS_dZ_loc,
          nGDenZ_dX_loc, nGDenZ_dY_loc, nGDenZ_dZ_loc,
          GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
          GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
          GGDenS_dxX_loc, GGDenS_dxY_loc, GGDenS_dxZ_loc,
          GGDenS_dyX_loc, GGDenS_dyY_loc, GGDenS_dyZ_loc,
          GGDenS_dzX_loc, GGDenS_dzY_loc, GGDenS_dzZ_loc,
          GGDenZ_dxX_loc, GGDenZ_dxY_loc, GGDenZ_dxZ_loc,
          GGDenZ_dyX_loc, GGDenZ_dyY_loc, GGDenZ_dyZ_loc,
          GGDenZ_dzX_loc, GGDenZ_dzY_loc, GGDenZ_dzZ_loc, 
          GradU_nX_loc, GradU_nY_loc, GradU_nZ_loc,
          GradU_gamma_X_loc, GradU_gamma_Y_loc, GradU_gamma_Z_loc, nAtoms);

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

        bool isGGA = std::any_of(functionals.begin(),functionals.end(),
                     [](std::shared_ptr<DFTFunctional> &x) {
                       return x->isGGA(); 
                     }); 

        // XC energy gradient
        for(size_t ic = 0; ic < nAtoms; ic++) {
          if (isGGA) {
            integrateXCEnergyGradX[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nX_loc[ic], GradU_gamma_X_loc[ic]);
            integrateXCEnergyGradY[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nY_loc[ic], GradU_gamma_Y_loc[ic]);
            integrateXCEnergyGradZ[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nZ_loc[ic], GradU_gamma_Z_loc[ic]);
          } else {
            integrateXCEnergyGradX[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nX_loc[ic], nullptr);
            integrateXCEnergyGradY[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nY_loc[ic], nullptr);
            integrateXCEnergyGradZ[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nZ_loc[ic], nullptr);
          }
        }

#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto botenergy_vxc     = 
          std::chrono::high_resolution_clock::now();
        auto topconstructZVars = 
          std::chrono::high_resolution_clock::now();
#endif
   
      }; // EXC gradient integrate

      // Create the BeckeIntegrator object
      BeckeIntegrator<EulerMac> 
        integrator(intComm,this->molecule(),basis,
        EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch,
          (isGGA ? GRADIENT : NOGRAD), intParam.epsilon);

      integrator.turn_on_grad();

      // Integrate the VXC
      integrator.integrate<size_t>(vxcbuild);


#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto toptransform    = std::chrono::high_resolution_clock::now();
#endif

      // Finishing up the EXC
      // factor in the 4 pi (Lebedev) and built the upper triagolar part
      // since we create only the lower triangular. For all components
      for(auto &X : integrateXCEnergy)
        XCEnergy += 4*M_PI*X;

      //std::cout << "XCEnergy in fromEXCGrad " << XCEnergy << std::endl;
      //std::cout << "XCGradient size " << XCGradient.size() << std::endl;
      // Finishing up the EXC Gradient
      for(size_t it = 0; it < nthreads; it++) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          this->XCGradient[ic][0] += 4*M_PI*integrateXCEnergyGradX[it][ic];
          this->XCGradient[ic][1] += 4*M_PI*integrateXCEnergyGradY[it][ic];
          this->XCGradient[ic][2] += 4*M_PI*integrateXCEnergyGradZ[it][ic];
        }
      }

      std::cout << "XCEnergy gradient" << std::endl;
      for(size_t ic = 0; ic < nAtoms; ic++) {
        std::cout << this->XCGradient[ic][0] << "  " << this->XCGradient[ic][1] << "  " << this->XCGradient[ic][2] << std::endl;
      }



// Combine MPI results
#ifdef CQ_ENABLE_MPI

      XCEnergy = MPIReduce(XCEnergy,0,intComm);

      for(size_t it = 0; it < nthreads; it++) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          XCGradient[ic][0] = MPIReduce(XCGradient[ic][0],0,intComm);
          XCGradient[ic][1] = MPIReduce(XCGradient[ic][1],0,intComm);
          XCGradient[ic][2] = MPIReduce(XCGradient[ic][2],0,intComm);
        }
      }

#endif

#if VXC_DEBUG_LEVEL >= 1
      // TIMING
      auto bottransform = std::chrono::high_resolution_clock::now();
#endif

#if VXC_DEBUG_LEVEL >= 3
      // DebugPrint
      std::cerr << std::endl;
      Scale(NB2,4*M_PI,tmpS,1);
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

#endif

      // Freeing the memory
      // ----------------Main System-------------------------------------------- //
      CQMemManager::get().free(SCRATCHNBNB,SCRATCHNBNP,DenS,epsEval,U_n);
      for(size_t ic = 0; ic < nAtoms; ic++) {
        CQMemManager::get().free(VEC_SCRATCHNBNB[ic][0], VEC_SCRATCHNBNB[ic][1], VEC_SCRATCHNBNB[ic][2]);
        CQMemManager::get().free(VEC_SCRATCHNBNP[ic][0], VEC_SCRATCHNBNP[ic][1], VEC_SCRATCHNBNP[ic][2]);
      }
      CQMemManager::get().free(VEC_SCRATCHNBNP2[0], VEC_SCRATCHNBNP2[1], VEC_SCRATCHNBNP2[2]);
      for(size_t ic = 0; ic < nAtoms; ic++) {
        CQMemManager::get().free(nGDenS_dX[ic],nGDenS_dY[ic],nGDenS_dZ[ic]);
        CQMemManager::get().free(GradU_nX[ic],GradU_nY[ic],GradU_nZ[ic]);
      }
      if( isGGA ) {
        CQMemManager::get().free(GDenS,U_gamma);
        for(size_t ic = 0; ic < nAtoms; ic++) {
          CQMemManager::get().free(GGDenS_dxX[ic],GGDenS_dxY[ic],GGDenS_dxZ[ic]);
          CQMemManager::get().free(GGDenS_dyX[ic],GGDenS_dyY[ic],GGDenS_dyZ[ic]);
          CQMemManager::get().free(GGDenS_dzX[ic],GGDenS_dzY[ic],GGDenS_dzZ[ic]);
          CQMemManager::get().free(GradU_gamma_X[ic],GradU_gamma_Y[ic],GradU_gamma_Z[ic]);
        }
      }

      if( this->onePDM->hasZ() ) {
        CQMemManager::get().free(DenZ);
        for(size_t ic = 0; ic < nAtoms; ic++) {
          CQMemManager::get().free(nGDenZ_dX[ic],nGDenZ_dY[ic],nGDenZ_dZ[ic]);
        }
        if( isGGA ) { 
          CQMemManager::get().free(GDenZ);
          for(size_t ic = 0; ic < nAtoms; ic++) {
            CQMemManager::get().free(GGDenZ_dxX[ic],GGDenZ_dxY[ic],GGDenZ_dxZ[ic]);
            CQMemManager::get().free(GGDenZ_dyX[ic],GGDenZ_dyY[ic],GGDenZ_dyZ[ic]);
            CQMemManager::get().free(GGDenZ_dzX[ic],GGDenZ_dzY[ic],GGDenZ_dzZ[ic]);
          }
        }
      }


      if( this->onePDM->hasXY() ) {
        CErr("Nuclear gradient is not implemented for 2-component systems (GHF or X2C)", std::cout);
      }

      // These quantities are used when there are multiple xc functionals
      if( functionals.size() > 1 ) {
        CQMemManager::get().free(epsSCR, dVU_n_SCR);
        if( isGGA ) CQMemManager::get().free(dVU_gamma_SCR);
      }

      if( not std::is_same<MatsT,double>::value ) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          for(auto &X : Re1PDMGrad_scalar[ic]) CQMemManager::get().free(X);
          for(auto &X : Re1PDMGrad_mz[ic]) CQMemManager::get().free(X);
        }
      }

      // -----------------End Main System-------------------------------------------- //
      // End freeing the memory


#if VXC_DEBUG_LEVEL >= 1
     // TIMING
     double d_batch = this->aoints.molecule().nAtoms * 
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

  }; // KohnSham::formEXCGradient



}; // namespace ChronusQ

