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

#include <chronusq_sys.hpp>

// Includes typedefs 
#ifndef LAPACK_COMPLEX_CPP
  #define LAPACK_COMPLEX_CPP
#endif
//#include <lapack/config.h>

#ifdef CQ_HAS_TA
  #include <tiledarray.h>
#endif

#include <lapack.hh>
#include <blas.hh>


// Choose linear algebra headers
#ifdef _CQ_MKL
#if 0
  #define MKL_Complex16 dcomplex // Redefine MKL complex type
  #define MKL_Complex8  std::complex<float> // Redefine MKL complex type 

  #include <mkl.h> // MKL
  //#include <blas/fortran.h>

  #ifdef CQ_ENABLE_MPI
    #include <mkl_scalapack.h>  
    #include <mkl_pblas.h>

    //#define CXXBLACS_BLACS_Complex16 double
    //#define CXXBLACS_BLACS_Complex8  float
    //
    //#define CXXBLACS_HAS_BLAS
    //#define CXXBLACS_HAS_LAPACK
    //#define CXXBLACS_HAS_PBLAS
    //#define CXXBLACS_HAS_SCALAPACK
  #endif
#endif
  extern "C" {
    int  MKL_Get_Max_Threads();
    void MKL_Set_Num_Threads(int);
    #define mkl_set_num_threads         MKL_Set_Num_Threads
    #define mkl_get_max_threads         MKL_Get_Max_Threads
  }
#else

  #ifdef CQ_ENABLE_MPI
//  #error CXXBLAS + nonMKL Not Tested!
  #endif

  //#define CXXBLACS_HAS_BLAS
  //#define CXXBLACS_HAS_LAPACK

  //// Why?
  //#ifndef CQ_HAS_TA
  //  #define CXXBLACS_BLAS_Complex16 std::complex<double>
  //  #define CXXBLACS_BLAS_Complex8  std::complex<float>
  //  //#define CXXBLACS_LAPACK_Complex16 double
  //  //#define CXXBLACS_LAPACK_Complex8  float
  //#endif

  //#include <blas/fortran.h>
  //#include <lapack/fortran.h>

  extern "C" {
    int openblas_get_num_threads();
    void openblas_set_num_threads(int*);
  }


#ifdef CQ_HAS_TA
  #define TA_Complex16 dcomplex // Redefine MKL complex type
  #define TA_Complex8  std::complex<float> // Redefine MKL complex type
  #define TA_INT int 
#endif


#endif



#include <blas.hh>
#include <lapack.hh>

//#ifdef CQ_ENABLE_MPI
//  #include <cxxblacs.hpp>
//#else
//  #define CB_INT int32_t
//#endif


#include <memmanager.hpp>
#include <Eigen/Core>

extern "C" void openblas_get_num_threads_(int*);
extern "C" void openblas_set_num_threads_(int*);

