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

#ifdef CQ_ENABLE_MPI
#include <scalapackpp/linear_systems/gesv.hpp>
#endif

namespace ChronusQ {

#ifdef CQ_ENABLE_MPI


  template <typename _F>
  inline int64_t LinSolve(const int64_t N, const int64_t NRHS, _F *A, 
    const int64_t IA, const int64_t JA, const scalapackpp::scalapack_desc DESCA, 
    _F *B, const int64_t IB, const int64_t JB,
    const scalapackpp::scalapack_desc DESCB, int64_t *IPIV) {

    return scalapackpp::pgesv(N,NRHS,A,IA,JA,DESCA,IPIV,B,IB,JB,DESCB);

  }

  template <typename _F>
  inline int64_t LinSolve(const int64_t N, const int64_t NRHS, _F *A, 
    const int64_t IA, const int64_t JA, const scalapackpp::scalapack_desc DESCA, 
    _F *B, const int64_t IB, const int64_t JB,
    const scalapackpp::scalapack_desc DESCB) {

    int64_t* iPIV = CQMemManager::get().calloc<int64_t>(DESCA[8] + DESCA[4]); // LLD + MB

    int64_t INFO = LinSolve(N,NRHS,A,IA,JA,DESCA,B,IB,JB,DESCB,iPIV);

    CQMemManager::get().free(iPIV);

    return INFO;
  }


#endif

}; // namespace ChronusQ

