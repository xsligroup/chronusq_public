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

#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <util/matout.hpp>
#include <util/timer.hpp>
#include <util/threads.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>

// Use stupid but bullet proof incore contraction for debug
//#define _BULLET_PROOF_INCORE_REL_RI

namespace ChronusQ {

  /**
   *  \brief Perform various tensor contractions of the full ERI
   *  tensor in core. Wraps other helper functions and provides
   *  loop structure
   *
   *  Currently supports
   *    - Coulomb-type (34,12) contractions
   *    - Exchange-type (23,12) contractions
   *
   *  Works with both real and complex matricies
   *
   *  \param [in/out] list Contains information pertinent to the 
   *    matricies to be contracted with. See TwoBodyContraction
   *    for details
   */
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::twoBodyContract3Index(
    MPI_Comm comm, std::vector<TwoBodyContraction<MatsT>> &list) const {

    ROOT_ONLY(comm);

    auto topIncore = tick();

    // Loop over matricies to contract with
    for(auto &C : list) {

      // Coulomb-type (34,12) ERI contraction
      // AX(mn) = (mn | kl) X(kl)
      if( C.contType == COULOMB ) {
        JContract3Index(comm,C);
      // Exchange-type (23,12) ERI contraction
      // AX(mn) = (mk |ln) X(kl)
      } else if( C.contType == EXCHANGE ) {
        KContract3Index(comm,C);
      }

    } // loop over matricies

    auto durIncore = tock(topIncore);


  }; // GTODirectRelERIContraction::twoBodyContractIncore


  /**
   *  \brief Perform a Coulomb-type (34,12) ERI contraction with
   *  a one-body operator.
   */
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::JContract3Index(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C) const {

    size_t NB  = this->ints()->nBasis();
    size_t NB2 = NB*NB;
    size_t NB3 = NB*NB2;

    // XSLI: This line breaks GIAO
//    if( C.ERI4 == nullptr) C.ERI4 = reinterpret_cast<double*>(ERI);


#ifdef _BULLET_PROOF_INCORE_REL_RI

    memset(C.AX,0.,NB*sizeof(MatsT));

    if( C.intTrans == TRANS_KL ) {

      // D(μν) = D(λκ)(μν|[κλ]^T) = D(λκ)(μν|λκ)
      #pragma omp parallel for
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[n] += C.ERI4[   n + l*NB + k*NB2] * C.X[l + k*NB];

      }
    } else if( C.intTrans == TRANS_MNKL ) {

      // D(μν) = D(λκ)(μν|κλ)^T = D(λκ)(κλ|μν)
      CErr("Invalid C.intTrans in 3-Index AO Direct",std::cout);
      
    } else if( C.intTrans == TRANS_NONE ) {

      // D(μν) = D(λκ)(μν|κλ)
      #pragma omp parallel for
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[n] += C.ERI4[    n + k*NB + l*NB2] * C.X[l + k*NB];

      }

    }
 

#else //_BULLET_PROOF_INCORE_REL_RI


    if( std::is_same<IntsT,dcomplex>::value )
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,NB2,MatsT(1.),C.ERI4,NB2,C.X,NB2,MatsT(0.),C.AX,NB2);
    else
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB2,1,NB2,MatsT(1.),C.ERI4,NB2,C.X,NB2,MatsT(0.),C.AX,NB2);

    // if Complex ints + Hermitian, conjugate
    if( std::is_same<IntsT,dcomplex>::value and C.HER )
      IMatCopy('R',NB,NB,MatsT(1.),C.AX,NB,NB);

    // If non-hermetian, transpose
    if( not C.HER )  {

      IMatCopy('T',NB,NB,MatsT(1.),C.AX,NB,NB);

    }

#endif //_BULLET_PROOF_INCORE_REL_RI

  }; // GTODirectRelERIContraction::JContract3Index


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::KContract3Index(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C) const {

    size_t NB  = this->ints()->nBasis();
    size_t NB2 = NB*NB;
    size_t NB3 = NB*NB2;

    // XSLI: This line breaks GIAO
//    if( C.ERI4 == nullptr) C.ERI4 = reinterpret_cast<double*>(ERI);

#ifdef _BULLET_PROOF_INCORE_REL_RI

    memset(C.AX,0.,NB*sizeof(MatsT));

    if( C.intTrans == TRANS_KL ) {
      
      // D(μν) = D(λκ)(μλ|[κν]^T) = D(λκ)(μλ|νκ)
      #pragma omp parallel for
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[n] += C.ERI4[   l + n*NB + k*NB2] * C.X[l + k*NB];

      }

    } else if( C.intTrans == TRANS_MNKL ) {
      
      // D(μν) = D(λκ)(μλ|κν)^T = D(λκ)(κν|μλ)
      CErr("Invalid C.intTrans in 3-Index AO Direct",std::cout);

    } else if( C.intTrans == TRANS_NONE ){

      // D(μν) = D(λκ)(μλ|κν)
      #pragma omp parallel for
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[n] += C.ERI4[  l + k*NB + n*NB2] * C.X[l + k*NB];

      }
    }

#else //_BULLET_PROOF_INCORE_REL_RI

    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);

    #pragma omp parallel for
    for(auto nu = 0; nu < NB; nu++)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,1,NB2,MatsT(1.),C.ERI4+nu*NB3,NB,C.X,NB2,MatsT(0.),C.AX+nu*NB,NB);

    SetLAThreads(LAThreads);

#endif //_BULLET_PROOF_INCORE_REL_RI

  }; // GTODirectRelERIContraction::KContract3Index

}; // namespace ChronusQ


