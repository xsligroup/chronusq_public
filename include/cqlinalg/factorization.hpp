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

namespace ChronusQ {

  /**
   *  \brief Computes the inverse of a non-singular matrix A. Initially
   *  computes the LU factorization then wraps DGETRI / ZGETRI for the
   *  matrix inversion depending on context.
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dgetri.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zgetri.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int LUInv(int N, _F *A, int LDA);

  /**
   *  \brief Computes the QR factorization of a general matrix A:
   *  returns Q in place of A and returns an upper triangular R.
   *  Smart wrapper around DGEQRF / DORGQR or ZGEQRF / ZUNGQR
   *  depending on the context
   */
  template <typename _F>
  int QR(int M, int N, _F *A, int LDA, _F *R, int LDR);


  /**
   *  \brief Computes the QR factorization of a general matrix A:
   *  returns Q in place of A and discards R
   *  Smart wrapper around DGEQRF / DORGQR or ZGEQRF / ZUNGQR
   *  depending on the context
   */
  template <typename _F>
  inline int QR(int M, int N, _F *A, int LDA) {

    int LDR = std::min(M,N);
    _F *R = CQMemManager::get().malloc<_F>(LDR*LDR);

    int INFO = QR(M,N,A,LDA,R,LDR);

    CQMemManager::get().free(R);

    return INFO;

  }


  template <typename _F>
  int OrdQZ(char JOBVSL, char JOBVSR, int N, _F *A, int LDA, _F *B, 
    int LDB, dcomplex *ALPHA, _F *BETA, double SIMGA, _F *VSL, 
    int LDVSL, _F *VSR, int LDVSR);

  template <typename _F>
  int OrdQZ2(char JOBVSL, char JOBVSR, int N, _F *A, int LDA, _F *B, 
    int LDB, dcomplex *ALPHA, _F *BETA, double hLim, double SIMGA, _F *VSL, 
    int LDVSL, _F *VSR, int LDVSR);


  /*
   *   Brief: This function computes an approximate matrix inverse. If the matrix
   *          is singular then it approximates the inverse by removing
   *          the singular columns from the matrix. Result is stored
   *          in place (A). num is the threshold to consider a singular
   *          value to be zero.
   */
  template<typename MatsT>
  void SVDInverse(const size_t N, MatsT* A, const size_t LDA, const double num);
  
}; // namespace ChronusQ

