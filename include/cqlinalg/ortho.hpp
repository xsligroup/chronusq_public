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

#include <cqlinalg/cqlinalg_config.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

  template <typename F, typename _VecNorm, typename _MatInner>
  size_t GramSchmidt(size_t N, size_t Mold, size_t Mnew, F *V, size_t LDV, 
    _VecNorm vecNorm, _MatInner matInner, size_t NRe = 0,
    double eps = 1e-12) {

    F * SCR = CQMemManager::get().malloc<F>(Mold + Mnew);

    if( Mold == 0 ) {
      // Normalize the first vector
      F inner = vecNorm(V);
      if(std::abs(inner) < std::sqrt(N)*eps) CErr("Zero inner product incurred!");
      blas::scal(N, 1./inner, V, 1);
    }

    // Orthonormalize the rest of the matrix using GS
    size_t iOrtho = (Mold == 0) ? Mold + 1: Mold;
    for(auto k = iOrtho; k < (Mold + Mnew); k++) {

      F* V_c = V + k      * LDV;
      F* V_p = V + iOrtho * LDV;

      if( k != iOrtho ) std::copy_n(V_c,N,V_p);

      // Project out the inner products
      for(auto iRe = 0; iRe < (NRe+1); iRe++) {
        matInner(iOrtho,1,V,LDV,V_p,LDV,SCR);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,1,iOrtho,F(-1.),V,LDV,SCR,iOrtho,F(1.),V_p,LDV);
      }

      // Normalize the new vector
      F inner = vecNorm(V_p);
      std::cout << k << " " << inner << std::endl;
      if(std::abs(inner) < std::sqrt(N)*eps) {
        std::cout << "Zero inner product incurred! " << k << "\n";
        blas::scal(N,F(0.),V_p,1);
      } else {
        blas::scal(N, 1./inner, V_p, 1);
        iOrtho++;
      }
    }

    CQMemManager::get().free(SCR);


#if 0

    size_t M = iOrtho;
    F* tInner = CQMemManager::get().malloc<F>(M*M);
    matInner(M,M,V,LDV,V,LDV,tInner);
    for(auto k = 0ul; k < M; k++) tInner[k*(M+1)] -= 1.;
    std::cerr << "Error after " << blas::nrm2(M*M,tInner,1) 
              << std::endl;

    CQMemManager::get().free(tInner);

#endif

    return iOrtho;


  };


  template <typename F>
  size_t GramSchmidt(size_t N, size_t Mold, size_t Mnew, F *V, size_t LDV, 
    size_t NRe = 0, double eps = 1e-12) {

    return
    GramSchmidt(N,Mold,Mnew,V,LDV,
      [&](F* Vc){ return std::sqrt(std::abs(blas::dot(N,Vc,1,Vc,1))); },
      [&](size_t i, size_t j, F* Vi, size_t LDVi, F* Vj, size_t LDVj, 
        F* inner){
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N,F(1.),Vi,LDVi,Vj,LDVj,F(0.),inner,i);
      },
      NRe,eps);

  };







  template <typename F, typename _MatInner>
  void GramSchmidt(size_t N, size_t Mold, size_t Mnew, F *VL, size_t LDVL, 
    F *VR, size_t LDVR, _MatInner matInner, size_t NRe = 0) {

    F * SCR = CQMemManager::get().malloc<F>(Mold + Mnew);

    if( Mold == 0 ) {
      // Normalize the first vector of each set
      F inner; matInner(1,1,VL,LDVL,VR,LDVR,&inner);
      if( std::abs(inner) < 1e-12 ) CErr("Zero Inner product incurred");

      blas::scal(N, F(1.) / std::sqrt(std::abs(inner)), VL, 1);
      blas::scal(N, F(1.) / std::sqrt(std::abs(inner)), VR, 1);
    }

    // Orthonormalize the rest of the matrix using GS
    for(auto k = Mold + 1; k < (Mold + Mnew); k++) {

      F* VR_c = VR + k*LDVR;
      F* VL_c = VL + k*LDVL;

      // Project out the inner products
      for(auto iRe = 0; iRe < (NRe+1); iRe++) {
        matInner(k,1,VR,LDVR,VL_c,LDVL,SCR);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,1,k,F(-1.),VL,LDVL,SCR,k,F(1.),VL_c,LDVL);
      }

      for(auto iRe = 0; iRe < (NRe+1); iRe++) {
        matInner(k,1,VL,LDVL,VR_c,LDVR,SCR);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,1,k,F(-1.),VR,LDVR,SCR,k,F(1.),VR_c,LDVR);
      }

      // Normalize the new vector
      F inner; matInner(1,1,VL_c,LDVL,VR_c,LDVR,&inner);
      if( std::abs(inner) < 1e-12 ) CErr("Zero Inner product incurred");

      blas::scal(N, F(1.) / std::sqrt(std::abs(inner)), VR_c, 1);
      blas::scal(N, F(1.) / std::sqrt(std::abs(inner)), VL_c, 1);

    }

    CQMemManager::get().free(SCR);


#if 0

    size_t M = Mold + Mnew;
    F* tInner = CQMemManager::get().malloc<F>(M*M);
    matInner(M,M,VL,LDVL,VR,LDVR,tInner);
    for(auto k = 0ul; k < M; k++) tInner[k*(M+1)] -= 1.;
    std::cerr << "Error after " << blas::nrm2(M*M,tInner,1) 
              << std::endl;

    CQMemManager::get().free(tInner);

#endif


  };


} // namespace ChronusQ

