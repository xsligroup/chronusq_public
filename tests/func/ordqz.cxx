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

#include <func.hpp>
#include <memmanager.hpp>

#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>


using namespace ChronusQ;

#ifndef _CQ_GENERATE_TESTS

/**
 *  \brief Generate a randon number of different
 *  types using the STL random number generator
 */
template <typename T>
inline T RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis);

template <>
inline double RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dis(e);

}; // RAND_NUMBER<double>

template <>
inline dcomplex RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dcomplex(dis(e),dis(e));

}; // RAND_NUMBER<dcomplex>


// OrdQZ test suite


  template <typename T>
  void ordqz_test(size_t N, double SIGMA) {

    std::random_device r; 
    std::default_random_engine e(r());
    std::uniform_real_distribution<> dis(-50,68);


    CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,256e6,256);

    T* A   = CQMemManager::get().template malloc<T>(N*N);
    T* B   = CQMemManager::get().template malloc<T>(N*N);
    T* VSL = CQMemManager::get().template malloc<T>(N*N);
    T* VSR = CQMemManager::get().template malloc<T>(N*N);



    T* AC = CQMemManager::get().template malloc<T>(N*N);
    T* BC = CQMemManager::get().template malloc<T>(N*N);

    for(auto j = 0ul; j < N*N; j++) {

      A[j] = RAND_NUMBER<T>(e,dis);
      B[j] = RAND_NUMBER<T>(e,dis);

    }

    std::copy_n(A,N*N,AC);
    std::copy_n(B,N*N,BC);


    dcomplex *ALPHA = CQMemManager::get().template malloc<dcomplex>(N);
    T        *BETA  = CQMemManager::get().template malloc<T>(N);

    // Compute QZ factorization
    OrdQZ('V','V',N,A,N,B,N,ALPHA,BETA,SIGMA,VSL,N,VSR,N);


    // Overwrite ALPHA with W - SIGMA
    for(auto j = 0; j < N; j++) 
      if( std::abs(BETA[j]) < 1e-13 )
        ALPHA[j] = std::numeric_limits<double>::infinity();
      else
        ALPHA[j] = std::abs(ALPHA[j] / BETA[j] - SIGMA);



    std::vector<size_t> indx(N);
    std::iota(indx.begin(),indx.end(),0);
    std::vector<size_t> indx_c(indx);

    // Sort based on W - SIGMA
    std::stable_sort(indx.begin(),indx.end(),
      [&](size_t i, size_t j) { 
        return ALPHA[i].real() < ALPHA[j].real(); }
    );

    // Check that the sort didn't do anything
    for(auto j = 0; j < N; j++) 
      EXPECT_EQ(indx[j],indx_c[j]);


    T* SCR = CQMemManager::get().template malloc<T>(N*N);


    // Compute AC = AC - VSL * S * VSR**H
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,N,N,N,T(1.) ,A  ,N,VSR,N,T(0.),SCR,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,T(-1.),VSL,N,SCR,N,T(1.),AC ,N);

    // Compute BC = BC - VSL * T * VSR**H
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,N,N,N,T(1.) ,B  ,N,VSR,N,T(0.),SCR,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,T(-1.),VSL,N,SCR,N,T(1.),BC ,N);

    auto abs_comp = []( T a, T b ) {

      return (std::abs(a) < std::abs(b));

    };

    // Check that the returned moities are indeed QZ
    double maxA = std::abs(*std::max_element(AC,AC+N*N,abs_comp));
    double maxB = std::abs(*std::max_element(BC,BC+N*N,abs_comp));


    CQMemManager::get().free(A, B, VSL, VSR, AC, BC, ALPHA, BETA, SCR);


    EXPECT_TRUE(maxA < 1e-10) <<  maxA;
    EXPECT_TRUE(maxB < 1e-10) <<  maxB;

  };


  TEST( ORDQZ, REAL_ORDQZ ) {
  
    ordqz_test<double>(50,-2);
  
  }
  
  TEST( ORDQZ, COMPLEX_ORDQZ ) {
  
    ordqz_test<dcomplex>(50,-2);
  
  }



#endif
