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
   *  \brief Helper function to obtain optimal LAPACK workspace
   *  dimension.
   *
   *  \param [in] func Function which returns the optimal workspace
   *  \returns         Optimal workspace dimension
   *
   *  z.B.
   *  
   *  For DSYEV, the optimal workspace may be queried by passing LWORK = -1
   *  (like all LAPACK functions). The parameter list for DSYEV is given as
   *
   *  void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, 
   *    double *W, double *WORK, int *LWORK, int *INFO);
   *
   *  To use this function to query the workspace, we may bind a set of 
   *  arguments to the funcion (i.e. the parameters from a DSYEV envocation)
   *
   *  auto test = std::bind(dsyev,&JOBZ,&UPLO,&N,A,&LDA,W,
   *    std::placeholders::_1, std::placeholders::_2,&INFO);
   *
   *  And then get the worksapce 
   *
   *  int LWORK = getLWork<double>(test);
   */ 
  template <typename T>
  int getLWork(std::function<void(T*,int*)> func) {
    int LWORK = -1;
    T LWorkQuery;
    func(&LWorkQuery,&LWORK);
  
    return (int)std::real(LWorkQuery);
  }; // getLWork

}; // namespace ChronusQ

