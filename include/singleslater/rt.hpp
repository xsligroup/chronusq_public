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

#include <singleslater.hpp>
#include <corehbuilder.hpp>
#include <fockbuilder.hpp>
#include <physcon.hpp>

#include <util/timer.hpp>
#include <cqlinalg/blasext.hpp>

#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>


//#define _DEBUGORTHO
//#define __DEBUGTPB__

namespace ChronusQ {

/**
 *  \brief Compute tau matrix for traveling proton basis
 * 
 *  τ_{QP} = \sum_{a = x,y,z} S_{QP}^{a} * \dot{R_P}
 * 
 */
template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeTau() {

    // Sanity checks
    if(this->particle.charge < 0)     CErr("No tau term for electronic subsystem");
    
    // Initialize / clear tau matrix
    size_t NB = basisSet().nBasis;
    if(!tau)
      tau = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
    std::fill_n(tau->pointer(), NB*NB, MatsT(0.));

    size_t numQProt  = this->molecule_.atomsQ.size();
    size_t NBPerProt = NB/numQProt;
    for (size_t iProt = 0; iProt < numQProt; iProt++) {
      for (size_t XYZ = 0; XYZ < 3; XYZ++) {
        size_t indQ = this->molecule_.atomsQ[iProt];

        // Multi-proton tau methods:
        // Option 1: Use velocity of the proton that the basis function center on 
        MatAdd('N','N', NB, NBPerProt, 
            MatsT(1.), tau->pointer()+iProt*NBPerProt*NB, NB, 
            MatsT(this->molecule_.atoms[indQ].velocity[XYZ]), (*this->aoints_->S0a)[XYZ]->matrix().pointer()+iProt*NBPerProt*NB, NB,
            tau->pointer()+iProt*NBPerProt*NB, NB);
        
        // Option 2: Diagonal approximation. Only generate 𝜏 within each proton
        //MatAdd('N','N', NBPerProt, NBPerProt, 
        //    MatsT(1.), tau->pointer()+iProt*NBPerProt*NB+iProt*NBPerProt, NB, 
        //    MatsT(this->molecule_.atoms[indQ].velocity[XYZ]), (*this->aoints_->S0a)[XYZ]->matrix().pointer()+iProt*NBPerProt*NB+iProt*NBPerProt, NB,
        //    tau->pointer()+iProt*NBPerProt*NB+iProt*NBPerProt, NB);

      //std::cout << "Proton velocity "+std::to_string(XYZ) << MatsT(this->molecule_.atoms[indQ].velocity[XYZ]) << std::endl;
      }
    }
    
#ifdef __DEBUGTPB__
    prettyPrintSmart(std::cout,"Tau Matrix", tau->pointer(),NB,NB,NB);
#endif
  }



  template <typename MatsT, typename IntsT> 
  void SingleSlater<MatsT,IntsT>::checkIdempotency(std::string system) {

    size_t NB  = basisSet().nBasis;
    size_t nSQ = NB*NB;

    MatsT* SCR = CQMemManager::get().malloc<MatsT>(nSQ);

    // compute P * P
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),
        onePDMOrtho->S().pointer(),NB,onePDMOrtho->S().pointer(),NB,MatsT(0.),SCR,NB);

    // compute P * P - P
    MatAdd('N','N',NB,NB,static_cast<MatsT>(1.),SCR,NB,static_cast<MatsT>(-1.),onePDMOrtho->S().pointer(), NB, SCR, NB);

    MatsT idem = 0.0;
    for(size_t p = 0; p < NB; p++)
      idem += SCR[p*NB+p];

    std::cout << std::scientific <<  std::setprecision(16) << system << " Density Idempotency: " << idem << std::endl;

    // free up scratch space 
    CQMemManager::get().free(SCR);

  }; // checkIdempotency



  template <typename MatsT, typename IntsT> 
  void SingleSlater<MatsT,IntsT>::mcWeenyPurification() {

    size_t NB  = basisSet().nBasis;
    size_t nSQ = NB*NB;

    MatsT* P2 = CQMemManager::get().malloc<MatsT>(nSQ);
    MatsT* P3 = CQMemManager::get().malloc<MatsT>(nSQ);

    std::vector<cqmatrix::Matrix<MatsT>> onePDMAB = this->onePDM->template spinGatherToBlocks<MatsT>(false);
    std::vector<cqmatrix::Matrix<MatsT>> purifiedOnePDMAB = {};

    for( size_t i = 0; i < onePDMAB.size(); i++ ) {
      purifiedOnePDMAB.emplace_back(NB);
      // compute P * P
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),
          onePDMAB[i].pointer(),NB,onePDMAB[i].pointer(),NB,MatsT(0.),P2,NB);
      // compute P * P * P
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),
          P2,NB,onePDMAB[i].pointer(),NB,MatsT(0.),P3,NB);
      // compute 3 P * P - 2 P * P * P
      MatAdd('N','N',NB,NB,MatsT(3.),P2,NB,MatsT(-2.), P3, NB, purifiedOnePDMAB[i].pointer(), NB);
    }

    setOnePDMAO(purifiedOnePDMAB.data());

    CQMemManager::get().free(P2, P3);

  }; // checkIdempotency



  template <typename MatsT, typename IntsT>
  template <typename M, enable_if_dcomplex<M>>
  void SingleSlater<MatsT, IntsT>::addTauToFock() {

    // Sanity checks
    if(this->particle.charge < 0)     CErr("No tau term for electronic subsystem");
    if(!tau)                          CErr("Tau should be initialized and computed");
      
    size_t NB = basisSet().nBasis;
#ifdef __DEBUGTPB__
    fockMatrix->output(std::cout, "Fock before adding Tau", true);
#endif
    MatAdd('N','N',NB, NB, dcomplex(1.), fockMatrix->S().pointer(), NB, dcomplex(0,-1.0),tau->pointer(), NB, fockMatrix->S().pointer(), NB);
    MatAdd('N','N',NB, NB, dcomplex(1.), fockMatrix->Z().pointer(), NB, dcomplex(0,-1.0),tau->pointer(), NB, fockMatrix->Z().pointer(), NB);
#ifdef __DEBUGTPB__
    fockMatrix->output(std::cout, "Fock after adding Tau", true);
#endif
  }




/**
 *  \brief Compute and returns time derivative of onePDM
 * 
 *  dP/dt = -i * ( i dP/dt )
 *        = -i * ( [F,P] - i( τ P + P τ^*) )
 *        = -i [F,P] - ( τ P + P τ^*)
 * 
 */
  template <typename MatsT, typename IntsT>
  template <typename M, enable_if_dcomplex<M>>
  cqmatrix::PauliSpinorMatrices<MatsT> SingleSlater<MatsT, IntsT>::getTimeDerDen(bool includeTau) {

    if(includeTau and !tau)           CErr("Tau should be initialized and computed");
      
    size_t NB = basisSet().nBasis;
      
    // Transfroma from spinor form (S/Z) to spin-block form (alpha/beta)
    std::vector<cqmatrix::Matrix<dcomplex>> onePDMOrthoAB     = onePDMOrtho->template spinGatherToBlocks<dcomplex>(false);
    std::vector<cqmatrix::Matrix<dcomplex>> fockMatrixOrthoAB = fockMatrixOrtho->template spinGatherToBlocks<dcomplex>(false);
    std::vector<cqmatrix::Matrix<dcomplex>> timeDerDenOrthoAB = {};

    for( size_t i = 0; i < onePDMOrthoAB.size(); i++ ) {
      timeDerDenOrthoAB.emplace_back(NB);
      // F P 
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,dcomplex(1.),
          fockMatrixOrthoAB[i].pointer(),NB,onePDMOrthoAB[i].pointer(),NB,dcomplex(0.),timeDerDenOrthoAB[i].pointer(),NB);
      // F P - P F
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,dcomplex(-1.),
          onePDMOrthoAB[i].pointer(),NB,fockMatrixOrthoAB[i].pointer(),NB,dcomplex(1.),timeDerDenOrthoAB[i].pointer(),NB);      
        
#ifdef __DEBUGTPB__
      timeDerDenOrthoAB[i].output(std::cout,"F P - P F for " + std::to_string(i), true);
#endif

      if(includeTau){
        cqmatrix::Matrix<dcomplex> SCR(NB);
        // τ * P 
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,dcomplex(1.),
            tauOrtho->pointer(),NB,onePDMOrthoAB[i].pointer(),NB,dcomplex(0.),SCR.pointer(),NB);
        // τ P + P τ^*
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,dcomplex(1.),
            onePDMOrthoAB[i].pointer(),NB,tauOrtho->pointer(),NB,dcomplex(1.),SCR.pointer(),NB); 
        // -i * [F,P] - ( τ P + P τ^*) 
        MatAdd('N','N',NB, NB, dcomplex(0.,-1.), timeDerDenOrthoAB[i].pointer(), NB, dcomplex(-1.), SCR.pointer(), NB, timeDerDenOrthoAB[i].pointer(), NB);
#ifdef __DEBUGTPB__
        SCR.output(std::cout,"τ P + P τ^* for " + std::to_string(i), true);
        timeDerDenOrthoAB[i].output(std::cout,"dP/dt for " + std::to_string(i), true);
#endif
      } else {
        // -i * ( [F,P]
        timeDerDenOrthoAB[i] *= dcomplex(0.,-1.);
      }
    }

    // Transform from spin-block form (alpha/beta) to spinor form (S/Z) 
    cqmatrix::PauliSpinorMatrices<MatsT> timeDerDenOrtho(NB, onePDMOrtho->hasXY(), onePDMOrtho->hasZ());
    if(nC == 1) {
      if(iCS) {
        timeDerDenOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(timeDerDenOrthoAB[0]);
      } else {
        timeDerDenOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(timeDerDenOrthoAB[0],timeDerDenOrthoAB[1]);
      }
    } else {
      timeDerDenOrtho = timeDerDenOrthoAB[0].template spinScatter<MatsT>();
    }

    return timeDerDenOrtho;

  }



/**
 *  \brief Perform Runge-kutta 4th order method to propagate the density to next timestep
 * 
 *  P(t+1) = P(t) + Δt/6 (k1 + 2*k1 + 2*k3 + k4)
 *
 *  Upon entry: onePDMOrtho stores the density matrix that need to be updated
 *  Upon exit:  onePDMOrtho stores the density matrix that have been propagated with RK4
 * 
 */
  template <typename MatsT, typename IntsT>
  template <typename M, enable_if_dcomplex<M>>
  void SingleSlater<MatsT,IntsT>::RK4Propagation(bool includeTau, double dt, bool increment, EMPerturbation& pert_tp5, EMPerturbation& pert_t1){

    size_t NB  = basisSet().nBasis;
    cqmatrix::PauliSpinorMatrices<MatsT> onePDMOrthoSave = *onePDMOrtho;

    if(includeTau){
      if(!tau) CErr("Tau should be initialized and computed");
      // Convert tau matrix from ao basis to orthonormal basis 
      if(!tauOrtho) tauOrtho = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
      *(tauOrtho) = orthoSpinor->nonortho2ortho(*tau);
    }
#ifdef __DEBUGTPB__
    tauOrtho->output(std::cout,"Ortho Tau Matrix",true);
    onePDMOrtho->output(std::cout, "initial P(t)", true);
#endif

    // =========================================================================================
    // Compute K1:
    // =========================================================================================
  
    // Compute k1 = dP(t)/dt = -i[F,P] - ( τ P + P τ^*) in spinor form (S/Z)
    cqmatrix::PauliSpinorMatrices<MatsT> k1 = getTimeDerDen(includeTau);
#ifdef __DEBUGTPB__
    k1.output(std::cout, "k1", true);
#endif

  

    // =========================================================================================
    // Compute K2:
    // =========================================================================================

    // Compute P^(k1) = P(t) + 0.5 * Δt * k1
    *onePDMOrtho += k1 * MatsT(0.5*dt);
    // Obtain new fock matrix at P^(k1)
    ortho2aoDen();
    formFock(pert_tp5, increment);
    ao2orthoFock();
    // Compute k2 = dP^(k1)/dt = -i[F,P] - ( τ P + P τ^*) in spinor form (S/Z)
    cqmatrix::PauliSpinorMatrices<MatsT> k2 = getTimeDerDen(includeTau);
#ifdef __DEBUGTPB__
    onePDMOrtho->output(std::cout, "P^(k1)", true);
    k2.output(std::cout, "k2", true);
#endif
  


    // =========================================================================================
    // Compute K3:
    // =========================================================================================
  
    // Restore P(t)
    *onePDMOrtho = onePDMOrthoSave;
    // Compute P^(k2) = P(t) + 0.5 * Δt * k2
    *onePDMOrtho += k2 * MatsT(0.5*dt);
    // Obtain new fock matrix at P^(k2)
    ortho2aoDen();
    formFock(pert_tp5, increment);
    ao2orthoFock();
    // Compute k3 = dP^(k2)/dt = -i[F,P] - ( τ P + P τ^*) in spinor form (S/Z)
    cqmatrix::PauliSpinorMatrices<MatsT> k3 = getTimeDerDen(includeTau);
#ifdef __DEBUGTPB__ 
    onePDMOrtho->output(std::cout, "P^(k2)", true);
    k3.output(std::cout, "k3", true);
#endif



    // =========================================================================================
    // Compute K4:
    // =========================================================================================

    // Restore P(t)
    
    // Compute P^(k3) = P(t) + Δt * k3
    *onePDMOrtho += k3 * MatsT(dt);
    // Obtain new fock matrix at P^(k3)
    ortho2aoDen();
    formFock(pert_t1, increment);
    ao2orthoFock();
    // Compute k4 = dP^(k3)/dt = -i[F,P] - ( τ P + P τ^*) in spinor form (S/Z)
    cqmatrix::PauliSpinorMatrices<MatsT> k4 = getTimeDerDen(includeTau);
#ifdef __DEBUGTPB__ 
    onePDMOrtho->output(std::cout, "P^(k3)", true);
    k4.output(std::cout, "k4", true);
#endif


  
    // =========================================================================================
    // Compute final P(t+1):
    // =========================================================================================

    // Restore P(t)
    *onePDMOrtho = onePDMOrthoSave;
    // P(t+1) = P(t) + Δt/6 (k1 + 2*k1 + 2*k3 + k4)
    *onePDMOrtho += dt / MatsT(6.0) * (k1 + MatsT(2.0)*k2 + MatsT(2.0)*k3 + k4);
#ifdef __DEBUGTPB__ 
    onePDMOrtho->output(std::cout, "P(t+1)", true);
#endif

    // Normalize progagated density
    computeNaturalOrbitals();
    formDensity();
  }



/**
 *  \brief Perform unitary propagation to propagate the density to next timestep,
 *         possibly perform magnus2 propagation.
 *         Propagation is done in alpha/beta basis
 * 
 *  P(t+1) = U P(t) U*
 *       U = e^(-i Δt F)
 * 
 *  Upon entry: onePDMOrtho stores the density matrix that need to be updated
 *  Upon exit:  onePDMOrtho stores the density matrix that have been propagated with a unitary propagator 
 * 
 */
template <typename MatsT, typename IntsT>
template <typename M, enable_if_dcomplex<M>>
void SingleSlater<MatsT, IntsT>::unitaryPropagation(bool includeTau, double dt, bool doMagnus2, EMPerturbation& pert_t1) {
  
  size_t NB  = basisSet().nBasis;

  // Gather orthonormal density from S/Z to A/B blocks
  std::vector<cqmatrix::Matrix<dcomplex>> onePDMOrthoAB  = getOnePDMOrtho();
  // Save a copy of orthonormal density before propagation if doing magnus2
  std::vector<cqmatrix::Matrix<dcomplex>> onePDMOrthoABSave  = {};
  if (doMagnus2) onePDMOrthoABSave = getOnePDMOrtho();

  // For traveling basis, add tau term to protonic fock matrix
  if(includeTau) addTauToFock();

  // Gather AO fock matrices from S/Z to A/B blocks
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fock_k = getFock();
  // Gather tranformation matrices
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho = getOrtho();
  // Convert AO fock matrices to orthonormal basis (in A/B blocks)
  std::vector<cqmatrix::Matrix<dcomplex>> fockMatrixOrthoAB = {};
  for( size_t i = 0; i < fock_k.size(); i++ ) {
    fockMatrixOrthoAB.emplace_back(fock_k[i]->dimension());
    if( ortho[i]->hasOverlap() )
      fockMatrixOrthoAB[i] = ortho[i]->nonortho2ortho(*fock_k[i]);
    else 
      fockMatrixOrthoAB[i] = *(fock_k[i]);
  }

  // Do propagation
  for (size_t i = 0; i < fockMatrixOrthoAB.size(); i++) {
    cqmatrix::Matrix<dcomplex> U(NB);
    cqmatrix::Matrix<dcomplex> SCR(NB);

    // U = e^(-i Δt F)
    MatExp('D', NB, dcomplex(0., -dt), fockMatrixOrthoAB[i].pointer(), NB, U.pointer(), NB);

    // P(t+1) = U P(t) U*
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, dcomplex(1.),
               U.pointer(), NB, onePDMOrthoAB[i].pointer(), NB, dcomplex(0.), SCR.pointer(), NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, dcomplex(1.),
               SCR.pointer(), NB, U.pointer(), NB, dcomplex(0.), onePDMOrthoAB[i].pointer(), NB);
  }

  if (doMagnus2) {
    // Convert orthonormal P(t+1) in A/B blocks to S/Z in ss.onePDMOrtho
    setOnePDMOrtho(onePDMOrthoAB.data());
    // Form a new fock matrix using AO P(t+1)
    ortho2aoDen();
    formFock(pert_t1, false);
    
    // For traveling basis, add tau term to protonic fock matrix
    if(includeTau) addTauToFock();
    
    // Gather AO fock matrices from S/Z to A/B blocks
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fock_k1 = getFock();
    // Compute 0.5 * (F(k) + F(k+1))
    for( size_t i = 0; i < fock_k.size(); i++ )
      *fock_k[i] = 0.5 * (*fock_k[i] + *fock_k1[i]); 
    // Convert AO fock matrices to orthonormal basis (in A/B blocks)
    for( size_t i = 0; i < fock_k.size(); i++ ) {
      if( ortho[i]->hasOverlap() )
        fockMatrixOrthoAB[i] = ortho[i]->nonortho2ortho(*fock_k[i]);
      else 
        fockMatrixOrthoAB[i] = *(fock_k[i]);
    }

    // Restore density before propagation
    onePDMOrthoAB = onePDMOrthoABSave;

    // Do progatation using new fock matrix
    for (size_t i = 0; i < fockMatrixOrthoAB.size(); i++) {
      cqmatrix::Matrix<dcomplex> U(NB);
      cqmatrix::Matrix<dcomplex> SCR(NB);

      // U = e^(-i Δt F)
      MatExp('D', NB, dcomplex(0., -dt), fockMatrixOrthoAB[i].pointer(), NB, U.pointer(), NB);

      // P(t+1) = U P(t) U*
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, dcomplex(1.),
                U.pointer(), NB, onePDMOrthoAB[i].pointer(), NB, dcomplex(0.), SCR.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, dcomplex(1.),
                SCR.pointer(), NB, U.pointer(), NB, dcomplex(0.), onePDMOrthoAB[i].pointer(), NB);
    }
  }

  // Convert orthonormal P(t+1) in A/B blocks to S/Z in ss.onePDMOrtho
  setOnePDMOrtho(onePDMOrthoAB.data());

}
  

}; // namespace ChronusQ

