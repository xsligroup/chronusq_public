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
#include <scalapackpp/pblas/gemm.hpp>
#include <scalapackpp/util/type_conversions.hpp>
#endif

namespace ChronusQ {

#ifdef CQ_ENABLE_MPI

  template <typename _F>
  void Gemm_MPI(char TRANSA, char TRANSB, int64_t M, int64_t N, int64_t K, _F ALPHA,
    _F *A, int64_t IA, int64_t JA, scalapackpp::scalapack_desc DESCA, 
    _F *B, int64_t IB, int64_t JB, scalapackpp::scalapack_desc DESCB, 
    _F BETA, _F *C, int64_t IC, int64_t JC, scalapackpp::scalapack_desc DESCC) {

    auto opa = scalapackpp::detail::to_op(TRANSA);
    auto opb = scalapackpp::detail::to_op(TRANSB);
    scalapackpp::pgemm(opa,opb,M,N,K,ALPHA,A,IA,JA,DESCA,B,IB,JB,DESCB,
      BETA,C,IC,JC,DESCC);

  }

#endif

}; // namespace ChronusQ


