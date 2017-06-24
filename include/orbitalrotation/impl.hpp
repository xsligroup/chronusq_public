/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <orbitalrotation.hpp>
#include <util/print.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/ortho.hpp>
#include <cqlinalg/matfunc.hpp>

//#define DEBUG_ORBITALROTATION_IMPL

namespace ChronusQ {
  
  /*
   * \brief rotate MO for one step size using Newton-Raphson:
   * 
   * steps:
   *   1. compute Gradient if not computed
   *   2. compute (approximated) Hessian and inverse it
   *   3. compute rotation paramemter X = - g / H
   *   4. compute rotation matrix by matrix exponential U = exp(X)
   *   5. rotate MO for one step size as C' = C U
   */ 

  template <typename MatsT, typename IntsT>  
  void OrbitalRotation<MatsT, IntsT>::rotateMO(EMPerturbation & pert, 
    cqmatrix::Matrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM) {
    
    auto & mopart = mcwfn_.MOPartition;
    auto & mo     = mcwfn_.reference().mo[0];
    size_t nTOrb  = mopart.nMO;
    size_t nTOrb2 = nTOrb * nTOrb;
    size_t nAO    = mo.dimension(); 
    
    // offset it by inactive core index
    auto mo_pointer = mo.pointer(); 
    
    // get gradient
    if (not orbitalGradient_) computeOrbGradient(pert, oneRDM, twoRDM); 
    MatsT * G = orbitalGradient_->pointer();

    // allocate memory    
    MatsT * X = CQMemManager::get().malloc<MatsT>(nTOrb2);
    MatsT * U = CQMemManager::get().malloc<MatsT>(nTOrb2);

    MatsT * H = nullptr; 
    if (settings.alg == ORB_ROT_2ND_ORDER) {
      H = CQMemManager::get().malloc<MatsT>(nTOrb2*nTOrb2); 
    } else {
      H = CQMemManager::get().malloc<MatsT>(nTOrb2); 
    }
    
    ProgramTimer::tick("Form Hessian");
    // X = - H ^{-1} * G
    if (settings.alg == ORB_ROT_2ND_ORDER) {
      // this->computeOrbOrbHessian(H);
      CErr("2nd Order orbital rotation not implemented");
    } else {
      this->computeOrbOrbHessianDiag(pert, oneRDM, twoRDM, H);
      for (auto i = 0ul; i < nTOrb2; i++) X[i] = - G[i] / H[i];
    }
    ProgramTimer::tock("Form Hessian");
    
    ProgramTimer::tick("Rotate MO");
    
    // Damp X to avoid Large Rotation
    auto X_max = std::max_element(X, X+nTOrb2, 
      [&] (MatsT a, MatsT b) { return std::norm(a) < std::norm(b);} );
    
    if(std::abs(*X_max) > settings.XDampTol) {
      double XDampFc = settings.XDampTol / std::abs(*X_max);
      auto pos = std::distance(X, X_max);
      auto j   = int(pos / nTOrb);
      auto i   = pos - j*nTOrb;
      std::cout << "    WARNING!: Large Rotation! XMax at X[" 
                << i << ", " << j << "] = " << *X_max
                << ", step scaled by " << XDampFc
                << "\n" << std::endl;
      
      blas::scal(nTOrb2, MatsT(XDampFc), X, 1);
    }
    
    // U = exp(X)
    MatExp(nTOrb, X, nTOrb, U, nTOrb);
    
#ifdef DEBUG_ORBITALROTATION_IMPL
    double HNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, H, nTOrb); 
    prettyPrintSmart(std::cout, " OR H ", H, nTOrb, nTOrb, nTOrb);
    double XNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, X, nTOrb); 
    prettyPrintSmart(std::cout, " OR X ", X, nTOrb, nTOrb, nTOrb);
    double UNorm = lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, U, nTOrb); 
    prettyPrintSmart(std::cout, " OR U ", U, nTOrb, nTOrb, nTOrb);
    std::cout << "OR H Norm = " << HNorm << std::endl;
    std::cout << "OR X Norm = " << XNorm << std::endl;
    std::cout << "OR U Norm = " << UNorm << std::endl;
#endif

    // Orthonormalized U and disable GramSchmidt printining
    std::cout.setstate(std::ios_base::failbit);
    size_t NUOrtho = GramSchmidt(nTOrb, 0, nTOrb, U, nTOrb);   
    std::cout.clear();
    
    if(NUOrtho != nTOrb) CErr("Failed at Orthonormalizing Rotation U Matix.");
    
    // Rotate MO Coefficient C' = C * U
    // use X as SCR
    
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
      nAO, nTOrb, nTOrb, MatsT(1.), mo_pointer, 
      nAO, U, nTOrb, MatsT(0.), X, nAO);
    SetMat('N', nAO, nTOrb, MatsT(1.), X, nAO, mo_pointer, nAO);  
    
    ProgramTimer::tock("Rotate MO");
    
    // free memory
    if(orbitalGradient_) orbitalGradient_ = nullptr;
    if(H) CQMemManager::get().free(H);
    if(X) CQMemManager::get().free(X);
    if(U) CQMemManager::get().free(U);

  }; // OrbitalRotation::rotateOrbitals()

}; // namespace ChronusQ

// Other headers
#include <orbitalrotation/fock.hpp>      // fock matrix implementation 
#include <orbitalrotation/gradient.hpp>  // gradient implementation
#include <orbitalrotation/hessian.hpp>   // hessian implementation
#include <orbitalrotation/ivo.hpp>   // hessian implementation
