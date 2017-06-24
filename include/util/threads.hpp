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
#include <cqlinalg/cqlinalg_config.hpp>
#ifdef __has_include
#  if __has_include(<openblas_config.h>)
#    define _CQ_OPENBLAS 1
#  endif  // __has_include(<openblas_config.h>)
#endif  // defined(__has_include)
#if !defined(_CQ_MKL) && !defined(_CQ_OPENBLAS)
#  include <thread>
#endif

namespace ChronusQ {

  inline void SetLAThreads(size_t n) {
#ifdef _CQ_MKL
    mkl_set_num_threads(n);
#else
#ifdef _CQ_OPENBLAS
    int N = n;
    openblas_set_num_threads_(&N);
#endif  // defined(_CQ_OPENBLAS)
#endif
    Eigen::setNbThreads(n);
  };


  inline size_t GetLAThreads() {
#ifdef _CQ_MKL
    return mkl_get_max_threads();
#else
#ifdef _CQ_OPENBLAS
    return openblas_get_num_threads();
#else
    return std::thread::hardware_concurrency();
#endif  // defined(_CQ_OPENBLAS)
#endif
  };


  inline void SetNumThreads(size_t n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
    SetLAThreads(n);
  };

  inline size_t GetNumThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  };

  inline size_t GetThreadID() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
  };


}; // namespace ChronusQ

