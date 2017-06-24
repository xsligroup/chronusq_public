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
#include <singleslater/neoss.hpp>

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
   *  \brief assemble NEO-DFT Gradient (exc gradient + epc gradient)
   *
   *
   *  This function is integrated by the BeckeIntegrator
   *  object.
   */  
  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::formEXCGradient() {

    // Obtaining subsystem as KS objects
    auto ess = getSubsystem<KohnSham>("Electronic");
    auto pss = getSubsystem<KohnSham>("Protonic");
    auto prot_ks = std::dynamic_pointer_cast<NEOKohnShamBuilder<MatsT,IntsT>>(this->fockBuilders["Protonic"].back());
    if(!ess or !pss or !prot_ks) CErr("Unsuccessful cast in NEOSS<MatsT,IntsT>::formEXCGradient()");
    
    auto epc_functionals = prot_ks->getFunctionals();
    auto intParam = ess->intParam;

#if VXC_DEBUG_LEVEL >= 1
    // TIMING 
    auto topMem = std::chrono::high_resolution_clock::now();
#endif
    ProgramTimer::tick("Form NEO-DFT Gradient");

    assert( intParam.nRad % intParam.nRadPerBatch == 0 );

    // Parallelism

    size_t nthreads = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(this->comm);
    size_t mpiSize   = MPISize(this->comm);

    // MPI Communicator for numerical integration
    // *** Assumes that MPI will be done only across atoms in integration

    size_t nAtoms = ess->molecule().nAtoms;
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(ess->comm,color,mpiRank);

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) 
#endif
    {  
  
      // Define several useful quantities for later on
      size_t NPtsMaxPerBatch = intParam.nRadPerBatch * intParam.nAng;
  
      bool isGGA = std::any_of(ess->functionals.begin(),ess->functionals.end(),
                      [](std::shared_ptr<DFTFunctional> &x) {
                        return x->isGGA(); 
                      });

      // NOTE: GGA epc gradient NYI!
      bool epcisGGA = std::any_of(epc_functionals.begin(),epc_functionals.end(),
                          [](std::shared_ptr<DFTFunctional> &x) {
                              return x->isGGA(); 
                          }); 

      if(epcisGGA) CErr(" GGA epc gradient NYI!");
  
      // Turn off LA threads
      SetLAThreads(1);

      // Electronic basis information
      BasisSet &basis = ess->basisSet();
      size_t NB       = basis.nBasis;
      size_t NB2      = NB*NB;
      size_t NTNB2    = nthreads * NB2;
      size_t NTNPPB   = nthreads*NPtsMaxPerBatch;
      size_t nAtoms   = this->molecule_.atoms.size();

      // Protonic basis information
      BasisSet &aux_basis = pss->basisSet();
      size_t aux_NB       = aux_basis.nBasis;
      size_t aux_NB2      = aux_NB*aux_NB;
      size_t aux_NTNB2    = nthreads * aux_NB2;

    
      // integrated EXC Energy per thread
      std::vector<double> integrateEXCEnergy(nthreads,0.);
      // integrated EPC Energy per thread (for benchmarking only)
      std::vector<double> integrateEPCEnergy(nthreads,0.);
      
      // integrated EXC Energy Gradient per thread
      std::vector<std::vector<double>> integrateEXCEnergyGradX(nthreads);
      std::vector<std::vector<double>> integrateEXCEnergyGradY(nthreads);
      std::vector<std::vector<double>> integrateEXCEnergyGradZ(nthreads);
      // integrated EPC Energy Gradient per thread
      std::vector<std::vector<double>> integrateEPCEnergyGradX(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradY(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradZ(nthreads);

      // integrated EPC Energy Gradient w.r.t electronic, per thread
      std::vector<std::vector<double>> integrateEPCEnergyGradXE(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradYE(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradZE(nthreads);

      // integrated EPC Energy Gradient w.r.t protonic, per thread
      std::vector<std::vector<double>> integrateEPCEnergyGradXP(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradYP(nthreads);
      std::vector<std::vector<double>> integrateEPCEnergyGradZP(nthreads);

      for(size_t it = 0; it < nthreads; it++) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          integrateEXCEnergyGradX[it].push_back(0.);
          integrateEXCEnergyGradY[it].push_back(0.);
          integrateEXCEnergyGradZ[it].push_back(0.);
          integrateEPCEnergyGradX[it].push_back(0.);
          integrateEPCEnergyGradY[it].push_back(0.);
          integrateEPCEnergyGradZ[it].push_back(0.);
          integrateEPCEnergyGradXE[it].push_back(0.);
          integrateEPCEnergyGradYE[it].push_back(0.);
          integrateEPCEnergyGradZE[it].push_back(0.);
          integrateEPCEnergyGradXP[it].push_back(0.);
          integrateEPCEnergyGradYP[it].push_back(0.);
          integrateEPCEnergyGradZP[it].push_back(0.);
        }
      }

 
      // Allocating Memory
      // --------------------------------------------------------------------------------//

      ess->XCEnergy = 0.;
      pss->XCEnergy = 0.;
      
      // --------------------------------------------------------------------------------//
      // --------------------Start Allocating Memory for Electronic System---------------//
      // --------------------------------------------------------------------------------//
      
      // Allocate scratch matrices
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

      // Allocate density and density gradient matrices
      double *DenS, *DenZ;
      std::vector<double*> nGDenS_dX, nGDenS_dY, nGDenS_dZ;
      std::vector<double*> nGDenZ_dX, nGDenZ_dY, nGDenZ_dZ;
      DenS = CQMemManager::get().template malloc<double>(NTNPPB);
      for (size_t ic = 0; ic < nAtoms; ic++) {
        nGDenS_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        nGDenS_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        nGDenS_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
      }
      if( ess->onePDM->hasZ() ) {
        DenZ = CQMemManager::get().template malloc<double>(NTNPPB);
        for (size_t ic = 0; ic < nAtoms; ic++) {
          nGDenZ_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          nGDenZ_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          nGDenZ_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        }
      }
      if( ess->onePDM->hasXY() )
        CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

      double *epsEval_elec = CQMemManager::get().template malloc<double>(NTNPPB);
      double *epsEval_prot = CQMemManager::get().template malloc<double>(NTNPPB);

      double *epcEval_elec = CQMemManager::get().template malloc<double>(NTNPPB);
      double *epcEval_prot = CQMemManager::get().template malloc<double>(NTNPPB);
      std::fill_n(epcEval_elec, NTNPPB, 0.);
      std::fill_n(epcEval_prot, NTNPPB, 0.);

      // Allocate Density U-Variables
      double *U_n   = 
        CQMemManager::get().template malloc<double>(2*NTNPPB); 
      double *dVU_n = 
        CQMemManager::get().template malloc<double>(2*NTNPPB);

      // EPC derivative d \epsilon / d rho^e
      double *dVU_n_elec = 
        CQMemManager::get().template malloc<double>(2*NTNPPB);
      // EPC derivative d \epsilon / d rho^p
      double *dVU_n_prot = 
        CQMemManager::get().template malloc<double>(2*NTNPPB);

      // Allocate Density gradient U-Variables
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
        if( ess->onePDM->hasZ() ) {
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
        if( ess->onePDM->hasXY() )
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

        // Allocate Gamma U-Variables
        U_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 
        dVU_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 

        // Allocate Gradients of Gamma U-Variables
        for(size_t ic = 0; ic < nAtoms; ic++) {
          GradU_gamma_X.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
          GradU_gamma_Y.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
          GradU_gamma_Z.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB));
        }
      }

      // These quantities are used when there are multiple xc functionals
      double *epsSCR, *dVU_n_SCR, *dVU_gamma_SCR;
      if( ess->functionals.size() > 1 ) {
        epsSCR    = CQMemManager::get().template malloc<double>(NTNPPB);
        dVU_n_SCR = CQMemManager::get().template malloc<double>(2*NTNPPB);
        if (isGGA)
          dVU_gamma_SCR = CQMemManager::get().template malloc<double>(3*NTNPPB);
      }

      // Decide if we need to allocate space for real part of the 
      // densities and copy over the real parts
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> Re1PDM;
      if (std::is_same<MatsT,double>::value)
        Re1PDM = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<double>>(ess->onePDM);
      else
        Re1PDM = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(ess->onePDM->real_part());

      // Obtain the 1PDM Gradient (stored in SingleSlater), which is calculated during Pulay Gradient
      // Note: Currently the 1PDM gradient is NOT used in this function. 
      //       We include this part of contribution in total pulay force
      std::vector<std::vector<double*>> Re1PDMGrad_scalar(nAtoms); 
      std::vector<std::vector<double*>> Re1PDMGrad_mz(nAtoms); 
      for(auto ic = 0; ic < nAtoms; ic++) {
        for (auto xyz = 0; xyz < 3; xyz++ )
          if( std::is_same<MatsT,double>::value )
            Re1PDMGrad_scalar[ic].push_back(reinterpret_cast<double*>(ess->onePDMGrad[3*ic+xyz]->S().pointer()));
          else{
            Re1PDMGrad_scalar[ic].push_back(CQMemManager::get().template malloc<double>(NB2));
            GetMatRE('N',NB,NB,1.,ess->onePDMGrad[3*ic+xyz]->S().pointer(),NB,Re1PDMGrad_scalar[ic].back(),NB);
          }
      }
      if( ess->onePDM->hasZ() ) {
        for(auto ic = 0; ic < nAtoms; ic++) {
          for(auto xyz = 0; xyz < 3; xyz++)
            if( std::is_same<MatsT,double>::value )
              Re1PDMGrad_mz[ic].push_back(reinterpret_cast<double*>(ess->onePDMGrad[3*ic+xyz]->Z().pointer()));
            else{
              Re1PDMGrad_mz[ic].push_back(CQMemManager::get().template malloc<double>(NB2));
              GetMatRE('N',NB,NB,1.,ess->onePDMGrad[3*ic+xyz]->Z().pointer(),NB,Re1PDMGrad_mz[ic].back(),NB);
            }
        }
      }
      // --------------------------------------------------------------------------------//
      // --------------------Finish Allocating Memory for Electronic System---------------//
      // --------------------------------------------------------------------------------//
    


      // --------------------------------------------------------------------------------//
      // --------------------Start Allocating Memory for Protonic System---------------//
      // --------------------------------------------------------------------------------//

      // Allocate scratch matrices
      double *AUX_SCRATCHNBNB =
        CQMemManager::get().template malloc<double>(aux_NTNB2);
      double *AUX_SCRATCHNBNP =
        CQMemManager::get().template malloc<double>(NTNPPB*aux_NB); 
      std::vector<std::vector<double*>> AUX_VEC_SCRATCHNBNB;
      for (size_t ic = 0; ic < nAtoms; ic++) {
        std::vector<double*> tmp =
          {CQMemManager::get().template malloc<double>(aux_NTNB2), CQMemManager::get().template malloc<double>(aux_NTNB2), CQMemManager::get().template malloc<double>(aux_NTNB2)};
        AUX_VEC_SCRATCHNBNB.emplace_back(tmp);
      }
      std::vector<std::vector<double*>> AUX_VEC_SCRATCHNBNP;
      for (size_t ic = 0; ic < nAtoms; ic++) {
        std::vector<double*> tmp = 
          {CQMemManager::get().template malloc<double>(NTNPPB*aux_NB), CQMemManager::get().template malloc<double>(NTNPPB*aux_NB), CQMemManager::get().template malloc<double>(NTNPPB*aux_NB)};
        AUX_VEC_SCRATCHNBNP.emplace_back(tmp);
      }
      std::vector<double*> AUX_VEC_SCRATCHNBNP2 = 
          {CQMemManager::get().template malloc<double>(NTNPPB*aux_NB), CQMemManager::get().template malloc<double>(NTNPPB*aux_NB), CQMemManager::get().template malloc<double>(NTNPPB*aux_NB)};

      // Density pointers for auxiliary system
      double *aux_DenS, *aux_DenZ;
      std::vector<double*> aux_nGDenS_dX, aux_nGDenS_dY, aux_nGDenS_dZ; // pointers for nuclear gradients
      std::vector<double*> aux_nGDenZ_dX, aux_nGDenZ_dY, aux_nGDenZ_dZ; // pointers for nuclear gradients
      aux_DenS = CQMemManager::get().template malloc<double>(NTNPPB);
      for(size_t ic = 0; ic < nAtoms; ic++) {
        aux_nGDenS_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        aux_nGDenS_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        aux_nGDenS_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
      }
      if( pss->onePDM->hasZ() ) {
        aux_DenZ = CQMemManager::get().template malloc<double>(NTNPPB);
        for(size_t ic = 0; ic < nAtoms; ic++) {
          aux_nGDenZ_dX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_nGDenZ_dY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_nGDenZ_dZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        }
      }
      if( pss->onePDM->hasXY() )
        CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

      // NEO auxiliary Density U-Variables
      double *aux_U_n =
        CQMemManager::get().template malloc<double>(2*NTNPPB);
      double *aux_dVU_n = 
        CQMemManager::get().template malloc<double>(2*NTNPPB);

      // NEO Density gradient U-Variables
      std::vector<double*> aux_GradU_nX, aux_GradU_nY, aux_GradU_nZ;
      for(size_t ic = 0; ic < nAtoms; ic++) {
        aux_GradU_nX.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
        aux_GradU_nY.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
        aux_GradU_nZ.emplace_back(CQMemManager::get().template malloc<double>(2*NTNPPB));
      }

      // These quantities are only used for GGA functionals
      double *aux_GDenS, *aux_GDenZ, *aux_U_gamma, *aux_dVU_gamma;
      std::vector<double*> aux_GGDenS_dxX, aux_GGDenS_dxY, aux_GGDenS_dxZ; 
      std::vector<double*> aux_GGDenS_dyX, aux_GGDenS_dyY, aux_GGDenS_dyZ; 
      std::vector<double*> aux_GGDenS_dzX, aux_GGDenS_dzY, aux_GGDenS_dzZ; 
      std::vector<double*> aux_GGDenZ_dxX, aux_GGDenZ_dxY, aux_GGDenZ_dxZ; 
      std::vector<double*> aux_GGDenZ_dyX, aux_GGDenZ_dyY, aux_GGDenZ_dyZ; 
      std::vector<double*> aux_GGDenZ_dzX, aux_GGDenZ_dzY, aux_GGDenZ_dzZ; 
      std::vector<double*> aux_GradU_gamma_X, aux_GradU_gamma_Y, aux_GradU_gamma_Z;
      if( epcisGGA ) {
        aux_GDenS = CQMemManager::get().template malloc<double>(3*NTNPPB);
        for (size_t ic = 0; ic < nAtoms; ic++) {
          aux_GGDenS_dxX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dxY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dxZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dyX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dyY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dyZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dzX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dzY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          aux_GGDenS_dzZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
        }
        if( pss->onePDM->hasZ() ) {
          aux_GDenZ = CQMemManager::get().template malloc<double>(3*NTNPPB);
          for (size_t ic = 0; ic < nAtoms; ic++) {
            aux_GGDenZ_dxX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dxY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dxZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dyX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dyY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dyZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dzX.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dzY.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
            aux_GGDenZ_dzZ.emplace_back(CQMemManager::get().template malloc<double>(NTNPPB));
          }
        }
        if( pss->onePDM->hasXY() ) 
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

        // Gamma U-Variables
        aux_U_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 
        aux_dVU_gamma = CQMemManager::get().template malloc<double>(3*NTNPPB); 

        // gradients of Gamma U-Variables
        for(size_t ic = 0; ic < nAtoms; ic++) {
          aux_GradU_gamma_X.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB)); 
          aux_GradU_gamma_Y.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB)); 
          aux_GradU_gamma_Z.emplace_back(CQMemManager::get().template malloc<double>(3*NTNPPB)); 
        }
      }

      double *aux_epsSCR, *aux_dVU_n_SCR, *aux_dVU_gamma_SCR;
      if( epc_functionals.size() > 1 ) {
        aux_epsSCR    = CQMemManager::get().template malloc<double>(NTNPPB);
        aux_dVU_n_SCR    = CQMemManager::get().template malloc<double>(2*NTNPPB);
        if(epcisGGA) 
          aux_dVU_gamma_SCR = CQMemManager::get().template malloc<double>(3*NTNPPB);
      }

      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> aux_Re1PDM;
      if (std::is_same<MatsT,double>::value)
        aux_Re1PDM = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<double>>(pss->onePDM);
      else
        aux_Re1PDM = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(pss->onePDM->real_part());

      // Note: Currently the 1PDM gradient is NOT used in this function. 
      //       We include this part of contribution in total pulay force
      std::vector<std::vector<double*>> aux_Re1PDMGrad_scalar(nAtoms); 
      std::vector<std::vector<double*>> aux_Re1PDMGrad_mz(nAtoms); 
      for(auto ic = 0; ic < nAtoms; ic++) {
        for (auto xyz = 0; xyz < 3; xyz++ )
          if( std::is_same<MatsT,double>::value )
            aux_Re1PDMGrad_scalar[ic].push_back(reinterpret_cast<double*>(pss->onePDMGrad[3*ic+xyz]->S().pointer()));
          else{
            aux_Re1PDMGrad_scalar[ic].push_back(CQMemManager::get().template malloc<double>(aux_NB2));
            GetMatRE('N',aux_NB,aux_NB,1.,pss->onePDMGrad[3*ic+xyz]->S().pointer(),aux_NB,aux_Re1PDMGrad_scalar[ic].back(),aux_NB);
          }
      }
      if( pss->onePDM->hasZ() ) {
        for(auto ic = 0; ic < nAtoms; ic++) {
          for(auto xyz = 0; xyz < 3; xyz++)
            if( std::is_same<MatsT,double>::value )
              aux_Re1PDMGrad_mz[ic].push_back(reinterpret_cast<double*>(pss->onePDMGrad[3*ic+xyz]->Z().pointer()));
            else{
              aux_Re1PDMGrad_mz[ic].push_back(CQMemManager::get().template malloc<double>(aux_NB2));
              GetMatRE('N',aux_NB,aux_NB,1.,pss->onePDMGrad[3*ic+xyz]->Z().pointer(),aux_NB,aux_Re1PDMGrad_mz[ic].back(),aux_NB);
            }
        }
      }
      // --------------------------------------------------------------------------------//
      // --------------------Finish Allocating Memory for Protonic System---------------//
      // --------------------------------------------------------------------------------//

      
      auto vxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
        std::vector<double> &weights, std::vector<size_t> NBE_vec, 
        std::vector<double*> BasisEval_vec, 
        std::vector<std::vector<size_t>> & batchEvalShells_vec, 
        std::vector<std::vector<std::pair<size_t,size_t>>> & subMatCut_vec) {

        // Electronic system
        size_t NBE = NBE_vec[0];
        double* BasisEval = BasisEval_vec[0];
        std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
        std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

        // Protonic system
        size_t aux_NBE = NBE_vec[1];
        double * aux_BasisEval = BasisEval_vec[1];
        std::vector<size_t> & aux_batchEvalShells = batchEvalShells_vec[1];
        std::vector<std::pair<size_t,size_t>> & aux_subMatCut = subMatCut_vec[1];

#if VXC_DEBUG_LEVEL > 3
        prettyPrintSmart(std::cout,"BASIS SCR",BasisEval,NBE,
          4*batch.size(),NBE);
        prettyPrintSmart(std::cout,"PROTONIC BASIS SCR",aux_BasisEval,NBE,
          4*batch.size(),NBE);
#endif 

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

        size_t aux_NDer = epcisGGA ? 4 : 1;
        double *aux_BasisGradEval = aux_BasisEval + aux_NDer * NPts * aux_NBE;
#if VXC_DEBUG_LEVEL >= 1
        // TIMING
        auto topevalDen = std::chrono::high_resolution_clock::now();
#endif

        // --------------------------------------------------------------------------------//
        // --------------------Start Set Up Local Memory for Electronic System------------//
        // --------------------------------------------------------------------------------//
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
          if (ess->onePDM->hasZ()) {
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
        if (isGGA) {
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

            if (ess->onePDM->hasZ()) {

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
        double * epsEval_elec_loc   = epsEval_elec   +  TIDNPPB;
        double * epsEval_prot_loc   = epsEval_prot   +  TIDNPPB;
        double * epcEval_elec_loc   = epcEval_elec   +  TIDNPPB;
        double * epcEval_prot_loc   = epcEval_prot   +  TIDNPPB;
        double * U_n_loc       = U_n       + 2*TIDNPPB;
        double * dVU_n_loc     = dVU_n     + 2*TIDNPPB;
        double * U_gamma_loc   = U_gamma   + 3*TIDNPPB;
        double * dVU_gamma_loc = dVU_gamma + 3*TIDNPPB;
        double * dVU_n_elec_loc     = dVU_n_elec     + 2*TIDNPPB;
        double * dVU_n_prot_loc     = dVU_n_prot     + 2*TIDNPPB;

        // local U gradient ponter
        std::vector<double*> GradU_nX_loc, GradU_nY_loc, GradU_nZ_loc;
        std::vector<double*> GradU_gamma_X_loc, GradU_gamma_Y_loc, GradU_gamma_Z_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
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
        // --------------------------------------------------------------------------------//
        // --------------------Finish Set Up Local Memory for Electronic System------------//
        // --------------------------------------------------------------------------------//


        // --------------------------------------------------------------------------------//
        // --------------------Start Set Up Local Memory for Protonic System------------//
        // --------------------------------------------------------------------------------//
        // aux local scratch pointers 
        double *  AUX_SCRATCHNBNB_loc = AUX_SCRATCHNBNB + thread_id * aux_NB2;
        double *  AUX_SCRATCHNBNP_loc = AUX_SCRATCHNBNP + thread_id * aux_NB*NPtsMaxPerBatch;         

        // aux local vector scratch pointers
        std::vector<std::vector<double*>> AUX_VEC_SCRATCHNBNB_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          std::vector<double*> tmp = 
            {AUX_VEC_SCRATCHNBNB[ic][0] + thread_id * aux_NB2, AUX_VEC_SCRATCHNBNB[ic][1] + thread_id * aux_NB2, AUX_VEC_SCRATCHNBNB[ic][2] + thread_id * aux_NB2};
          AUX_VEC_SCRATCHNBNB_loc.emplace_back(tmp);
        }

        std::vector<std::vector<double*>> AUX_VEC_SCRATCHNBNP_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          std::vector<double*> tmp = 
            {AUX_VEC_SCRATCHNBNP[ic][0] + thread_id * aux_NB*NPtsMaxPerBatch, AUX_VEC_SCRATCHNBNP[ic][1] + thread_id * aux_NB*NPtsMaxPerBatch, AUX_VEC_SCRATCHNBNP[ic][2] + thread_id * aux_NB*NPtsMaxPerBatch};
          AUX_VEC_SCRATCHNBNP_loc.emplace_back(tmp);
        }
        
        std::vector<double*> AUX_VEC_SCRATCHNBNP2_loc =
          {AUX_VEC_SCRATCHNBNP2[0] + thread_id * aux_NB*NPtsMaxPerBatch, AUX_VEC_SCRATCHNBNP2[1] + thread_id * aux_NB*NPtsMaxPerBatch, AUX_VEC_SCRATCHNBNP2[2] + thread_id * aux_NB*NPtsMaxPerBatch};

        // aux local density pointers
        double * aux_DenS_loc = aux_DenS + TIDNPPB;
        double * aux_DenZ_loc = aux_DenZ + TIDNPPB;

        std::vector<double*> aux_nGDenS_dX_loc, aux_nGDenS_dY_loc, aux_nGDenS_dZ_loc;
        std::vector<double*> aux_nGDenZ_dX_loc, aux_nGDenZ_dY_loc, aux_nGDenZ_dZ_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          aux_nGDenS_dX_loc.emplace_back(aux_nGDenS_dX[ic] + TIDNPPB);
          aux_nGDenS_dY_loc.emplace_back(aux_nGDenS_dY[ic] + TIDNPPB);
          aux_nGDenS_dZ_loc.emplace_back(aux_nGDenS_dZ[ic] + TIDNPPB);
          if (pss->onePDM->hasZ()) {
            aux_nGDenZ_dX_loc.emplace_back(aux_nGDenZ_dX[ic] + TIDNPPB);
            aux_nGDenZ_dY_loc.emplace_back(aux_nGDenZ_dY[ic] + TIDNPPB);
            aux_nGDenZ_dZ_loc.emplace_back(aux_nGDenZ_dZ[ic] + TIDNPPB);
          }
        }

        // aux local density gradient pointers
        double * aux_GDenS_loc = aux_GDenS + 3*TIDNPPB;
        double * aux_GDenZ_loc = aux_GDenZ + 3*TIDNPPB;

        // aux local density gradient nuclear gradient pointers
        std::vector<double*> aux_GGDenS_dxX_loc, aux_GGDenS_dxY_loc, aux_GGDenS_dxZ_loc;
        std::vector<double*> aux_GGDenS_dyX_loc, aux_GGDenS_dyY_loc, aux_GGDenS_dyZ_loc;
        std::vector<double*> aux_GGDenS_dzX_loc, aux_GGDenS_dzY_loc, aux_GGDenS_dzZ_loc;
        std::vector<double*> aux_GGDenZ_dxX_loc, aux_GGDenZ_dxY_loc, aux_GGDenZ_dxZ_loc;
        std::vector<double*> aux_GGDenZ_dyX_loc, aux_GGDenZ_dyY_loc, aux_GGDenZ_dyZ_loc;
        std::vector<double*> aux_GGDenZ_dzX_loc, aux_GGDenZ_dzY_loc, aux_GGDenZ_dzZ_loc;
        if (epcisGGA) {
          for (size_t ic = 0; ic < nAtoms; ic++) {
            aux_GGDenS_dxX_loc.emplace_back(aux_GGDenS_dxX[ic] + TIDNPPB);
            aux_GGDenS_dxY_loc.emplace_back(aux_GGDenS_dxY[ic] + TIDNPPB);
            aux_GGDenS_dxZ_loc.emplace_back(aux_GGDenS_dxZ[ic] + TIDNPPB);

            aux_GGDenS_dyX_loc.emplace_back(aux_GGDenS_dyX[ic] + TIDNPPB);
            aux_GGDenS_dyY_loc.emplace_back(aux_GGDenS_dyY[ic] + TIDNPPB);
            aux_GGDenS_dyZ_loc.emplace_back(aux_GGDenS_dyZ[ic] + TIDNPPB);

            aux_GGDenS_dzX_loc.emplace_back(aux_GGDenS_dzX[ic] + TIDNPPB);
            aux_GGDenS_dzY_loc.emplace_back(aux_GGDenS_dzY[ic] + TIDNPPB);
            aux_GGDenS_dzZ_loc.emplace_back(aux_GGDenS_dzZ[ic] + TIDNPPB);

            if (pss->onePDM->hasZ()) {
              aux_GGDenZ_dxX_loc.emplace_back(aux_GGDenZ_dxX[ic] + TIDNPPB);
              aux_GGDenZ_dxY_loc.emplace_back(aux_GGDenZ_dxY[ic] + TIDNPPB);
              aux_GGDenZ_dxZ_loc.emplace_back(aux_GGDenZ_dxZ[ic] + TIDNPPB);

              aux_GGDenZ_dyX_loc.emplace_back(aux_GGDenZ_dyX[ic] + TIDNPPB);
              aux_GGDenZ_dyY_loc.emplace_back(aux_GGDenZ_dyY[ic] + TIDNPPB);
              aux_GGDenZ_dyZ_loc.emplace_back(aux_GGDenZ_dyZ[ic] + TIDNPPB);

              aux_GGDenZ_dzX_loc.emplace_back(aux_GGDenZ_dzX[ic] + TIDNPPB);
              aux_GGDenZ_dzY_loc.emplace_back(aux_GGDenZ_dzY[ic] + TIDNPPB);
              aux_GGDenZ_dzZ_loc.emplace_back(aux_GGDenZ_dzZ[ic] + TIDNPPB);
            }
          }
        }

        // aux local U and V pointers
        double * aux_U_n_loc       = aux_U_n       + 2*TIDNPPB;
        double * aux_dVU_n_loc     = aux_dVU_n     + 2*TIDNPPB;
        double * aux_U_gamma_loc   = aux_U_gamma   + 3*TIDNPPB;
        double * aux_dVU_gamma_loc = aux_dVU_gamma + 3*TIDNPPB;

        // aux local U gradient ponter
        std::vector<double*> aux_GradU_nX_loc, aux_GradU_nY_loc, aux_GradU_nZ_loc;
        std::vector<double*> aux_GradU_gamma_X_loc, aux_GradU_gamma_Y_loc, aux_GradU_gamma_Z_loc;
        for (size_t ic = 0; ic < nAtoms; ic++) {
          aux_GradU_nX_loc.emplace_back(aux_GradU_nX[ic] + 2*TIDNPPB);
          aux_GradU_nY_loc.emplace_back(aux_GradU_nY[ic] + 2*TIDNPPB);
          aux_GradU_nZ_loc.emplace_back(aux_GradU_nZ[ic] + 2*TIDNPPB);
          if (epcisGGA) {
            aux_GradU_gamma_X_loc.emplace_back(aux_GradU_gamma_X[ic] + 3*TIDNPPB);
            aux_GradU_gamma_Y_loc.emplace_back(aux_GradU_gamma_Y[ic] + 3*TIDNPPB);
            aux_GradU_gamma_Z_loc.emplace_back(aux_GradU_gamma_Z[ic] + 3*TIDNPPB);
          }
        }

        // aux local multiple functional pointers
        double * aux_epsSCR_loc        = aux_epsSCR        +   TIDNPPB;
        double * aux_dVU_n_SCR_loc     = aux_dVU_n_SCR     + 2*TIDNPPB;
        double * aux_dVU_gamma_SCR_loc = aux_dVU_gamma_SCR + 3*TIDNPPB;
        // --------------------------------------------------------------------------------//
        // --------------------Finish Set Up Local Memory for Protonic System------------//
        // --------------------------------------------------------------------------------//

        

        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, Re1PDM->S().pointer(), DenS_loc, 
          GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

        // This evaluates the Gradient of V Variables
        evalDenGrad((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut,
          SCRATCHNBNB_loc, SCRATCHNBNP_loc, VEC_SCRATCHNBNB_loc, VEC_SCRATCHNBNP_loc,
          VEC_SCRATCHNBNP2_loc, Re1PDMGrad_scalar, Re1PDM->S().pointer(), 
          nGDenS_dX_loc,  nGDenS_dY_loc,  nGDenS_dZ_loc, 
          GGDenS_dxX_loc, GGDenS_dxY_loc, GGDenS_dxZ_loc, 
          GGDenS_dyX_loc, GGDenS_dyY_loc, GGDenS_dyZ_loc, 
          GGDenS_dzX_loc, GGDenS_dzY_loc, GGDenS_dzZ_loc, 
          BasisEval, BasisGradEval, nAtoms, basis);
        
        if( ess->onePDM->hasZ() ) {
          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            SCRATCHNBNB_loc ,SCRATCHNBNP_loc, Re1PDM->Z().pointer(), DenZ_loc, 
            GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

          evalDenGrad((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut,
            SCRATCHNBNB_loc, SCRATCHNBNP_loc, VEC_SCRATCHNBNB_loc, VEC_SCRATCHNBNP_loc,
            VEC_SCRATCHNBNP2_loc, Re1PDMGrad_mz, Re1PDM->Z().pointer(), 
            nGDenZ_dX_loc,  nGDenZ_dY_loc,  nGDenZ_dZ_loc, 
            GGDenZ_dxX_loc, GGDenZ_dxY_loc, GGDenZ_dxZ_loc, 
            GGDenZ_dyX_loc, GGDenZ_dyY_loc, GGDenZ_dyZ_loc, 
            GGDenZ_dzX_loc, GGDenZ_dzY_loc, GGDenZ_dzZ_loc, 
            BasisEval, BasisGradEval, nAtoms, basis);
        }

        if( ess->onePDM->hasXY() )
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);

        evalDen((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, aux_NB, aux_subMatCut, 
          AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, aux_Re1PDM->S().pointer(), aux_DenS_loc,
          aux_GDenS_loc, aux_GDenS_loc + NPts, aux_GDenS_loc + 2*NPts, aux_BasisEval);

        evalDenGrad((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, aux_NB, aux_subMatCut,
          AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, AUX_VEC_SCRATCHNBNB_loc, AUX_VEC_SCRATCHNBNP_loc,
          AUX_VEC_SCRATCHNBNP2_loc, aux_Re1PDMGrad_scalar, aux_Re1PDM->S().pointer(), 
          aux_nGDenS_dX_loc,  aux_nGDenS_dY_loc,  aux_nGDenS_dZ_loc, 
          aux_GGDenS_dxX_loc, aux_GGDenS_dxY_loc, aux_GGDenS_dxZ_loc, 
          aux_GGDenS_dyX_loc, aux_GGDenS_dyY_loc, aux_GGDenS_dyZ_loc, 
          aux_GGDenS_dzX_loc, aux_GGDenS_dzY_loc, aux_GGDenS_dzZ_loc, 
          aux_BasisEval, aux_BasisGradEval, nAtoms, aux_basis);

        if( pss->onePDM->hasZ() ) {
          evalDen((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, aux_NB, aux_subMatCut, 
            AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, aux_Re1PDM->Z().pointer(), aux_DenZ_loc,
            aux_GDenZ_loc, aux_GDenZ_loc + NPts, aux_GDenZ_loc + 2*NPts, aux_BasisEval);

          evalDenGrad((epcisGGA ? GRADIENT : NOGRAD), NPts, aux_NBE, aux_NB, aux_subMatCut,
            AUX_SCRATCHNBNB_loc, AUX_SCRATCHNBNP_loc, 
            AUX_VEC_SCRATCHNBNB_loc, AUX_VEC_SCRATCHNBNP_loc,
            AUX_VEC_SCRATCHNBNP2_loc, aux_Re1PDMGrad_mz, aux_Re1PDM->Z().pointer(),
            aux_nGDenZ_dX_loc,  aux_nGDenZ_dY_loc,  aux_nGDenZ_dZ_loc,  
            aux_GGDenZ_dxX_loc, aux_GGDenZ_dxY_loc, aux_GGDenZ_dxZ_loc, 
            aux_GGDenZ_dyX_loc, aux_GGDenZ_dyY_loc, aux_GGDenZ_dyZ_loc, 
            aux_GGDenZ_dzX_loc, aux_GGDenZ_dzY_loc, aux_GGDenZ_dzZ_loc, 
            aux_BasisEval, aux_BasisGradEval, nAtoms, aux_basis);
        }

        if( pss->onePDM->hasXY() )
          CErr("Nuclear Gradient is not implemented for 2-component systems (GHF or X2C)!", std::cout);        


        // V -> U variables for evaluating the kernel derivatives, not needed for energy only evaluation
        mkAuxVar(ess->onePDM,isGGA,epsScreen,NPts,
          DenS_loc,DenZ_loc,nullptr,nullptr,
          GDenS_loc, GDenS_loc+NPts, GDenS_loc+2*NPts,
          GDenZ_loc, GDenZ_loc+NPts, GDenZ_loc+2*NPts,
          nullptr,nullptr,nullptr,
          nullptr,nullptr,nullptr,
          nullptr, 
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr,
          nullptr, U_n_loc,U_gamma_loc
        );

        //grad V -> grad U variables for evaluateing the kernel gradients
        mkAuxVarGrad(ess->onePDM,isGGA,epsScreen,NPts,
          nGDenS_dX_loc, nGDenS_dY_loc, nGDenS_dZ_loc,
          nGDenZ_dX_loc, nGDenZ_dY_loc, nGDenZ_dZ_loc,
          GDenS_loc, GDenS_loc+NPts, GDenS_loc+2*NPts,
          GDenZ_loc, GDenZ_loc+NPts, GDenZ_loc+2*NPts,
          GGDenS_dxX_loc, GGDenS_dxY_loc, GGDenS_dxZ_loc,
          GGDenS_dyX_loc, GGDenS_dyY_loc, GGDenS_dyZ_loc,
          GGDenS_dzX_loc, GGDenS_dzY_loc, GGDenS_dzZ_loc,
          GGDenZ_dxX_loc, GGDenZ_dxY_loc, GGDenZ_dxZ_loc,
          GGDenZ_dyX_loc, GGDenZ_dyY_loc, GGDenZ_dyZ_loc,
          GGDenZ_dzX_loc, GGDenZ_dzY_loc, GGDenZ_dzZ_loc, 
          GradU_nX_loc, GradU_nY_loc, GradU_nZ_loc,
          GradU_gamma_X_loc, GradU_gamma_Y_loc, GradU_gamma_Z_loc, nAtoms
        );
        

        mkAuxVar(pss->onePDM,epcisGGA,epsScreen,NPts,
          aux_DenS_loc,aux_DenZ_loc,nullptr,nullptr,
          aux_GDenS_loc, aux_GDenS_loc+NPts, aux_GDenS_loc+2*NPts,
          aux_GDenZ_loc, aux_GDenZ_loc+NPts, aux_GDenZ_loc+2*NPts,
          nullptr,nullptr,nullptr,
          nullptr,nullptr,nullptr,
          nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr,
          nullptr, nullptr, nullptr, 
          aux_U_n_loc, aux_U_gamma_loc
        );

        //grad V -> grad U variables for evaluating the kernel gradients
        mkAuxVarGrad(pss->onePDM,epcisGGA,epsScreen,NPts,
          aux_nGDenS_dX_loc, aux_nGDenS_dY_loc, aux_nGDenS_dZ_loc,
          aux_nGDenZ_dX_loc, aux_nGDenZ_dY_loc, aux_nGDenZ_dZ_loc,
          aux_GDenS_loc, aux_GDenS_loc+NPts, aux_GDenS_loc+2*NPts,
          aux_GDenZ_loc, aux_GDenZ_loc+NPts, aux_GDenZ_loc+2*NPts,
          aux_GGDenS_dxX_loc, aux_GGDenS_dxY_loc, aux_GGDenS_dxZ_loc,
          aux_GGDenS_dyX_loc, aux_GGDenS_dyY_loc, aux_GGDenS_dyZ_loc,
          aux_GGDenS_dzX_loc, aux_GGDenS_dzY_loc, aux_GGDenS_dzZ_loc,
          aux_GGDenZ_dxX_loc, aux_GGDenZ_dxY_loc, aux_GGDenZ_dxZ_loc,
          aux_GGDenZ_dyX_loc, aux_GGDenZ_dyY_loc, aux_GGDenZ_dyZ_loc,
          aux_GGDenZ_dzX_loc, aux_GGDenZ_dzY_loc, aux_GGDenZ_dzZ_loc, 
          aux_GradU_nX_loc, aux_GradU_nY_loc, aux_GradU_nZ_loc,
          aux_GradU_gamma_X_loc, aux_GradU_gamma_Y_loc, aux_GradU_gamma_Z_loc, nAtoms);


        // Get EXC Energy derivatives wrt electronic U variables
        loadVXCder(ess->functionals, NPts, U_n_loc, U_gamma_loc, epsEval_elec_loc, 
          dVU_n_loc, dVU_gamma_loc, epsSCR_loc, dVU_n_SCR_loc, 
          dVU_gamma_SCR_loc); 

        // Get EPC Energy derivatives wrt electronic and protonic U variables
        loadEPCGradder(epc_functionals, NPts, 
          U_n_loc, aux_U_n_loc, 
          epsEval_elec_loc, epsEval_prot_loc, 
          dVU_n_elec_loc,dVU_n_prot_loc);

        // Compute for the current batch the XC energy and increment the 
        // total XC energy.
        integrateEXCEnergy[thread_id] += 
          energy_vxc(NPts, weights, epsEval_elec_loc, DenS_loc);

        integrateEPCEnergy[thread_id] += 
          energy_vxc(NPts, weights, epsEval_prot_loc, aux_DenS_loc);
        
        // EXC energy gradient
        for(size_t ic = 0; ic < nAtoms; ic++) {
          if (isGGA) {
            integrateEXCEnergyGradX[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nX_loc[ic], GradU_gamma_X_loc[ic]);
            integrateEXCEnergyGradY[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nY_loc[ic], GradU_gamma_Y_loc[ic]);
            integrateEXCEnergyGradZ[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nZ_loc[ic], GradU_gamma_Z_loc[ic]);
          } else {
            integrateEXCEnergyGradX[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nX_loc[ic], nullptr);
            integrateEXCEnergyGradY[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nY_loc[ic], nullptr);
            integrateEXCEnergyGradZ[thread_id][ic] += 
              energy_vxc_grad(isGGA, NPts, weights, dVU_n_loc, dVU_gamma_loc, GradU_nZ_loc[ic], nullptr);
          }
        }


        for(size_t ic = 0; ic < nAtoms; ic++) {
          integrateEPCEnergyGradXE[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_elec_loc, nullptr, GradU_nX_loc[ic], nullptr);
          integrateEPCEnergyGradYE[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_elec_loc, nullptr, GradU_nY_loc[ic], nullptr);
          integrateEPCEnergyGradZE[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_elec_loc, nullptr, GradU_nZ_loc[ic], nullptr);
          integrateEPCEnergyGradXP[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_prot_loc, nullptr, aux_GradU_nX_loc[ic], nullptr);
          integrateEPCEnergyGradYP[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_prot_loc, nullptr, aux_GradU_nY_loc[ic], nullptr);
          integrateEPCEnergyGradZP[thread_id][ic] += 
            energy_vxc_grad(false, NPts, weights, dVU_n_prot_loc, nullptr, aux_GradU_nZ_loc[ic], nullptr);
        }
        
        }; // VXC integrate


        
        // Create the BeckeIntegrator object
        BeckeIntegrator<EulerMac> 
          integrator(intComm,ess->molecule(),basis,aux_basis,
          EulerMac(intParam.nRad), intParam.nAng, intParam.nRadPerBatch,
            (isGGA ? GRADIENT : NOGRAD), (epcisGGA ? GRADIENT : NOGRAD), 
            intParam.epsilon);

        integrator.turn_on_grad();

        // Integrate the VXC
        integrator.integrate<size_t>(vxcbuild);

        double EXCEnergy = 0.0;
        double EPCEnergy = 0.0;
        for(auto &X : integrateEXCEnergy)  EXCEnergy += 4*M_PI*X;
        for(auto &X : integrateEPCEnergy)  EPCEnergy += 4*M_PI*X;

        //std::cout << std::setprecision(16) << "EXC Energy in Gradient: " << EXCEnergy << std::endl;
        //std::cout << std::setprecision(16) << "EPC Energy in Gradient: " << EPCEnergy << std::endl;

        // Finishing up the EXC Gradient
        for(size_t it = 0; it < nthreads; it++) {
          for(size_t ic = 0; ic < nAtoms; ic++) {
            this->EXCGradient[ic][0]  += 4*M_PI*integrateEXCEnergyGradX[it][ic];
            this->EXCGradient[ic][1]  += 4*M_PI*integrateEXCEnergyGradY[it][ic];
            this->EXCGradient[ic][2]  += 4*M_PI*integrateEXCEnergyGradZ[it][ic];
            this->EPCGradientE[ic][0] += 4*M_PI*integrateEPCEnergyGradXE[it][ic];
            this->EPCGradientE[ic][1] += 4*M_PI*integrateEPCEnergyGradYE[it][ic];
            this->EPCGradientE[ic][2] += 4*M_PI*integrateEPCEnergyGradZE[it][ic];
            this->EPCGradientP[ic][0] += 4*M_PI*integrateEPCEnergyGradXP[it][ic];
            this->EPCGradientP[ic][1] += 4*M_PI*integrateEPCEnergyGradYP[it][ic];
            this->EPCGradientP[ic][2] += 4*M_PI*integrateEPCEnergyGradZP[it][ic];
            this->EPCGradient[ic][0]  += this->EPCGradientE[ic][0] + this->EPCGradientP[ic][0];
            this->EPCGradient[ic][1]  += this->EPCGradientE[ic][1] + this->EPCGradientP[ic][1];
            this->EPCGradient[ic][2]  += this->EPCGradientE[ic][2] + this->EPCGradientP[ic][2];
          }
        }
        
        /***
        std::cout << "EXCEnergy gradient" << std::endl;

        std::cout << "EXCEnergy gradient" << std::endl;
        for(size_t ic = 0; ic < nAtoms; ic++) {
          std::cout << this->EXCGradient[ic][0] << "  " << this->EXCGradient[ic][1] << "  " << this->EXCGradient[ic][2] << std::endl;
        }

        std::cout << "EPCEnergy gradient" << std::endl;
        for(size_t ic = 0; ic < nAtoms; ic++) {
          std::cout << this->EPCGradient[ic][0] << "  " << this->EPCGradient[ic][1] << "  " << this->EPCGradient[ic][2] << std::endl;
        }

        std::cout << "EPCEnergy gradient E" << std::endl;
        for(size_t ic = 0; ic < nAtoms; ic++) {
          std::cout << this->EPCGradientE[ic][0] << "  " << this->EPCGradientE[ic][1] << "  " << this->EPCGradientE[ic][2] << std::endl;
        }

        std::cout << "EPCEnergy gradient P" << std::endl;
        for(size_t ic = 0; ic < nAtoms; ic++) {
          std::cout << this->EPCGradientP[ic][0] << "  " << this->EPCGradientP[ic][1] << "  " << this->EPCGradientP[ic][2] << std::endl;
        }
        ***/



// Combine MPI results
#ifdef CQ_ENABLE_MPI

      EXCEnergy = MPIReduce(EXCEnergy,0,intComm);
      EPCEnergy = MPIReduce(EPCEnergy,0,intComm);

      for(size_t it = 0; it < nthreads; it++) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          EXCGradient[ic][0] = MPIReduce(EXCGradient[ic][0],0,intComm);
          EXCGradient[ic][1] = MPIReduce(EXCGradient[ic][1],0,intComm);
          EXCGradient[ic][2] = MPIReduce(EXCGradient[ic][2],0,intComm);
          EPCGradient[ic][0] = MPIReduce(EPCGradient[ic][0],0,intComm);
          EPCGradient[ic][1] = MPIReduce(EPCGradient[ic][1],0,intComm);
          EPCGradient[ic][2] = MPIReduce(EPCGradient[ic][2],0,intComm);
          EPCGradientE[ic][0] = MPIReduce(EPCGradientE[ic][0],0,intComm);
          EPCGradientE[ic][1] = MPIReduce(EPCGradientE[ic][1],0,intComm);
          EPCGradientE[ic][2] = MPIReduce(EPCGradientE[ic][2],0,intComm);
          EPCGradientP[ic][0] = MPIReduce(EPCGradientP[ic][0],0,intComm);
          EPCGradientP[ic][1] = MPIReduce(EPCGradientP[ic][1],0,intComm);
          EPCGradientP[ic][2] = MPIReduce(EPCGradientP[ic][2],0,intComm);
        }
      }

#endif


      // Freeing the memory
      // ----------------Main System-------------------------------------------- //
      CQMemManager::get().free(SCRATCHNBNB,SCRATCHNBNP,DenS,epsEval_elec,epsEval_prot,U_n);
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

      if( ess->onePDM->hasZ() ) {
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


      if( ess->onePDM->hasXY() ) {
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
      
      // ----------------Auxiliary System-------------------------------------------- //
      CQMemManager::get().free(AUX_SCRATCHNBNB,AUX_SCRATCHNBNP,aux_DenS,aux_U_n);
      for(size_t ic = 0; ic < nAtoms; ic++) {
        CQMemManager::get().free(AUX_VEC_SCRATCHNBNB[ic][0], AUX_VEC_SCRATCHNBNB[ic][1], AUX_VEC_SCRATCHNBNB[ic][2]);
        CQMemManager::get().free(AUX_VEC_SCRATCHNBNP[ic][0], AUX_VEC_SCRATCHNBNP[ic][1], AUX_VEC_SCRATCHNBNP[ic][2]);
      }
      CQMemManager::get().free(AUX_VEC_SCRATCHNBNP2[0], AUX_VEC_SCRATCHNBNP2[1], AUX_VEC_SCRATCHNBNP2[2]);
      for(size_t ic = 0; ic < nAtoms; ic++) {
        CQMemManager::get().free(aux_nGDenS_dX[ic],aux_nGDenS_dY[ic],aux_nGDenS_dZ[ic]);
        CQMemManager::get().free(aux_GradU_nX[ic], aux_GradU_nY[ic], aux_GradU_nZ[ic]);
      }
      if( epcisGGA ) {
        CQMemManager::get().free(aux_GDenS,aux_U_gamma);
        for(size_t ic = 0; ic < nAtoms; ic++) {
          CQMemManager::get().free(aux_GGDenS_dxX[ic],aux_GGDenS_dxY[ic],aux_GGDenS_dxZ[ic]);
          CQMemManager::get().free(aux_GGDenS_dyX[ic],aux_GGDenS_dyY[ic],aux_GGDenS_dyZ[ic]);
          CQMemManager::get().free(aux_GGDenS_dzX[ic],aux_GGDenS_dzY[ic],aux_GGDenS_dzZ[ic]);
          CQMemManager::get().free(aux_GradU_gamma_X[ic],aux_GradU_gamma_Y[ic],aux_GradU_gamma_Z[ic]);
        }
      }

      if( pss->onePDM->hasZ() ) {
        CQMemManager::get().free(aux_DenZ);
        for(size_t ic = 0; ic < nAtoms; ic++) {
          CQMemManager::get().free(aux_nGDenZ_dX[ic],aux_nGDenZ_dY[ic],aux_nGDenZ_dZ[ic]);
        }
        if ( epcisGGA ) {
          CQMemManager::get().free(aux_GDenZ);
          for(size_t ic = 0; ic < nAtoms; ic++) {
            CQMemManager::get().free(aux_GGDenZ_dxX[ic],aux_GGDenZ_dxY[ic],aux_GGDenZ_dxZ[ic]);
            CQMemManager::get().free(aux_GGDenZ_dyX[ic],aux_GGDenZ_dyY[ic],aux_GGDenZ_dyZ[ic]);
            CQMemManager::get().free(aux_GGDenZ_dzX[ic],aux_GGDenZ_dzY[ic],aux_GGDenZ_dzZ[ic]);
          }
        }
      }

      if( pss->onePDM->hasXY() ) {
        CErr("Nuclear gradient is not implemented for 2-component systems (GHF or X2C)", std::cout);
      }

      if( epc_functionals.size() > 1 ) {
        CQMemManager::get().free(aux_epsSCR, aux_dVU_n_SCR);
        if( isGGA ) CQMemManager::get().free(aux_dVU_gamma_SCR);
      }

      if( not std::is_same<MatsT,double>::value ) {
        for(size_t ic = 0; ic < nAtoms; ic++) {
          for(auto &X : aux_Re1PDMGrad_scalar[ic]) CQMemManager::get().free(X);
          for(auto &X : aux_Re1PDMGrad_mz[ic]) CQMemManager::get().free(X);
        }
      }
      // -----------------End Aux-------------------------------------------- //
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

