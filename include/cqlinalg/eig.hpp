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

  // General (Non-Hermetian) Eigenproblem
    
  /**
   *  \brief Smart wrapper around DGEEV and ZGEEV depending on context.
   *
   *  Templated function which diagonalizes a general (non-hermetian)
   *  matrix which is general to real and complex. Handles memory managment
   *  internally through a CQMemManager object. Always returns complex
   *  eigenvalues for generality.
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dgeev.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zgeev.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int GeneralEigen(char JOBVL, char JOBVR, int N, _F *A, int LDA, dcomplex *W,
                   _F *VL, int LDVL, _F *VR, int LDVR);


  // Hermetian Eigenproblem

  /**
   *  \brief Smart wrapper around DSYEV and ZHEEV depending on context.
   *
   *  Templated function which diagonalizes a hermetian matrix which 
   *  is general to real and complex. Handles memory managment
   *  internally through a CQMemManager object. Returns real
   *  eigenvalues.
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dsyev.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zheev.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int HermetianEigen(char JOBZ, char UPLO, int N, _F *A, int LDA, double *W);

  /**
   *  \brief Smart wrapper around DSYEV and ZHEEV depending on context.
   *
   *  Templated function which diagonalizes a hermetian matrix which 
   *  is general to real and complex. Handles memory managment
   *  internally through a CQMemManager object. Returns complex
   *  eigenvalues.
   *
   *  See http://www.netlib.org/lapack/lapack-3.1.1/html/dsyev.f.html or
   *      http://www.netlib.org/lapack/lapack-3.1.1/html/zheev.f.html for
   *  parameter documentation.
   */ 
  template <typename _F>
  int HermetianEigen(char JOBZ, char UPLO, int N, _F *A, int LDA, dcomplex *W);

  template <typename _F>
  int HermetianEigen(char JOBZ, char UPLO, int N, double *A, int LDA, dcomplex *W);


}; // namespace ChronusQ

