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

namespace ChronusQ{
/*
 * Brief: Saves the initial set of orbitals to be used as the
 *        reference for the unitary transformation
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::saveRefMOs(vecMORef<MatsT>& mo) {

  // Store Reference MO's from which to compute rotations
  size_t nMO = mo.size();

  // Allocate Matrices
  if( not refMOsAllocated ) {
    refMO.reserve(nMO);
    for( size_t i = 0; i < nMO; i++ ) {
      size_t NBC = mo[i].get().dimension();
      refMO.emplace_back(NBC);
    }
  }

  // Copy MO's
  for( size_t i = 0; i < nMO; i++ ) {
    size_t NBC = mo[i].get().dimension();
    refMO[i] = mo[i].get();
  }

  storedRef = true;
}

/*
 *     Brief: Rotate the vector of MO's using the computed OrbRot parameters
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::rotateMOs(vecMORef<MatsT>& mo) {

  std::vector<cqmatrix::Matrix<MatsT>> U = computeUnitary();

  size_t disp = 0;
  for( size_t i = 0; i < mo.size(); i++ ) {
    size_t NBC   = mo[i].get().dimension();

    // Rotate MO's
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        NBC, NBC, NBC, 
        MatsT(1.), refMO[i].pointer(), NBC, 
        U[i].pointer(), NBC, 
        MatsT(0.), mo[i].get().pointer(), NBC);

#ifdef _NRSCF_PRINT_MOS
    prettyPrintSmart(std::cout, "MO " + std::to_string(i), mo[i].get().pointer(), NBC, NBC, NBC);
#endif
  }
}

template<typename MatsT>
std::vector<cqmatrix::Matrix<MatsT>> NewtonRaphsonSCF<MatsT>::computeUnitary(){

  vecShrdPtrMat<MatsT> den = this->orbitalModifierDrivers.getOnePDM();
  size_t nMat = den.size();

  // Initialize Anti-hermitian matrix
  std::vector<cqmatrix::Matrix<MatsT>> A;
  A.reserve(nMat);
  for( size_t i=0; i<nMat; ++i )
    A.emplace_back(den[i]->dimension());
  for( auto& a : A ) a.clear();

  // Parse parameters to make antisymmetric matrices
  size_t disp = 0;
  for( auto& rO : rotOpt ){
    for( auto& iO : rO.rotIndices.first ){
      for( auto& iV : rO.rotIndices.second ){
        A[rO.spaceIndex](iO, iV) += -orbRot[disp];
        A[rO.spaceIndex](iV, iO) += SmartConj(orbRot[disp]);
        ++disp;
      }
    }
  }

  // Compute the Unitary matrices
  std::vector<cqmatrix::Matrix<MatsT>> U;
  U.reserve(nMat);
  for( size_t i=0; i<nMat; ++i )
    U.emplace_back(den[i]->dimension());
  for( auto& u : U ) u.clear();

  // Compute either exp(A) or Cayley transform
  for( size_t i=0; i<nMat; ++i ){
    size_t N = U[i].dimension();
    // Try to compute exp(A) or use Cayley transform if failed
    try{
      MatExp(N,A[i].pointer(),N,U[i].pointer(),N);
    }
    catch(...){
        CErr("Matrix Exponential failed to converge");
    }
#ifdef _NRSCF_DEBUG_UNITARY
    prettyPrintSmart(std::cout, "Unitary Matrix", U[i].pointer(),N,N,N);
    cqmatrix::Matrix<MatsT> SCR(N);
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
        N, N, N, 
        MatsT(1.), U[i].pointer(), N, 
        U[i].pointer(), N, 
        MatsT(0.), SCR.pointer(), N);
    for( size_t p = 0; p < N; p++ )
      SCR(p, p) -= MatsT(1.);
    double sum = 0.;
    for( size_t p = 0; p < N * N; p++ )
      sum += std::abs(SCR.pointer()[p]);
    std::cerr << "Unitary Test:  U^H * U = 1 -> Error =  " << sum << std::endl;
#endif
  }
  return U;
}


/*
 *  Brief: Computes the orbital gradients and the diagonal Hessian and stores
 *         them in their vectors. The total gradient is the concatonation of
 *         orbital gradients for each Fock matrix. e.g. UHF->[Grad_alpha,Grad_beta]
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::computeGradient(vecMORef<MatsT>& mo) {

  // Compute the Gradient of orbital rotations
  if( this->orbitalModifierDrivers.computeNROrbGrad ) {
    // Compute user-defined gradient
    this->orbitalModifierDrivers.computeNROrbGrad(orbGrad);
  } else {
    // Compute default gradient from fock and mo's
      
      
    // Compute MO Fock
    size_t nMat = mo.size();
    vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
    std::vector<cqmatrix::Matrix<MatsT>> moFock;
    moFock.reserve(nMat);
    for( size_t i=0; i<nMat; ++i ) {
      size_t NB = mo[i].get().dimension();
      moFock.emplace_back( fock[i]->transform('N',mo[i].get().pointer(),NB,NB) );
    }

    // Compute Gradient for parameters
    size_t disp = 0;
    for( auto& rO : rotOpt ){
      for( auto& iO : rO.rotIndices.first ){
        for( auto& iV : rO.rotIndices.second ){
           orbGrad[disp] = moFock[rO.spaceIndex](iO,iV);
          ++disp;
        }
      }
    }

    // Conjugate the gradient for complex wave functions.
    // this is necessary to make quasi-newton and gradient descent
    // consistent with the full NR implementation
    /*
    if( std::is_same<MatsT, dcomplex>::value ){
      for(size_t iP=0; iP<nParam; ++iP)
          orbGrad[iP] = SmartConj(orbGrad[iP]);
    }
    */
  } // if( computeNROrbGrad ) else
#ifdef _NRSCF_DEBUG_GRAD
  printParamVec("Gradient",orbGrad);
#endif
}; // NewtonRaphsonSCF<MatsT> :: computeGradient

template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::computeDiagHess(vecEPtr& eps){

  size_t disp = 0;
  for( auto& rO : rotOpt ){
    for( auto& iO : rO.rotIndices.first ){
      double eO = eps[rO.spaceIndex][iO] + this->scfControls.nrLevelShift;
      for( auto& iV : rO.rotIndices.second ){
        orbDiagHess[disp] = eps[rO.spaceIndex][iV] - eO;
        ++disp;
      }
    }
  }
#ifdef _NRSCF_DEBUG_GRAD
  printParamVec("Diagonal Orbital Hessian",orbDiagHess);
#endif
}

/*
 * Brief: This computes the new set of parameters from 
 *        a given search direction (dx)
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::takeStep(MatsT* dx) {

  // Compute norm of step
  double norm = blas::nrm2(nParam, dx, 1);
#ifdef _NRSCF_DEBUG
  std::cout << "Search Direction Norm = " << norm << std::endl;
  printParamVec("Search Direction", dx);
#endif

  double scale = norm > this->scfControls.nrTrust ? this->scfControls.nrTrust / norm : 1.;

  // Perform Damping
  if( this->doingDamp )
    scale *= double(1. - this->scfControls.dampParam);

  // take Step
  blas::axpy(nParam, MatsT(-scale), dx, 1, orbRot, 1);
};

}// namespace :: ChronusQ
