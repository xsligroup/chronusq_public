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
#include <cqlinalg/solve.hpp>
#include <matrix.hpp>

namespace ChronusQ {

  /**
   *   \brief The DIIS class. A class to perform a DIIS extrapolation 
   *    based on a series of error metrics stored in core. 
   *
   */

  template <typename T>
  class DIIS {

  protected:

    // Useful typedefs
    typedef T*                        oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

  public:

    size_t         nExtrap;     ///< Size of extrapolation space
    std::vector<T> coeffs;      ///< Vector of extrapolation coeficients
    const std::vector<std::vector<cqmatrix::Matrix<T>>> &errorMetric; ///< Vector of vectors containing error metrics

    // Constructor
      
    /**
     *  DIIS Constructor. Constructs a DIIS object
     *
     *  \param [in]  nExtrap     Size of extrapolation space
     *  \param [in]  nMat        Number of matrices to trace for each element of B
     *  \param [in]  OSize       Size of the error metrics used to construct B
     *  \param [in]  errorMetric Vector of vectors containing error metrics
     *  \param [out] InvFail     Boolean of whether matrix inversion failed
     */ 
    DIIS(size_t nExtrap, const std::vector<std::vector<cqmatrix::Matrix<T>>> &errorMetric) :
      nExtrap(nExtrap), errorMetric(errorMetric) {

      coeffs.resize(nExtrap+1);

    };


    // Constructors for default, copy, and move
    DIIS() = delete; 
    DIIS(const DIIS &) = delete; 
    DIIS(DIIS &&) = delete;

    // Public Member functions
    bool extrapolate();

  }; // class DIIS



  /**
   *  \brief Performs a DIIS extrapolation using the vectors stored 
   *  in errorMetric
   *
   */ 
  template<typename T>
  bool DIIS<T>::extrapolate(){

    int N          = nExtrap + 1;
    int NRHS       = 1;
    bool InvFail   = false;
    cqmatrix::Matrix<T> B(N);
    B.clear();

    // Build the B matrix
    for( auto l=0ul; l<errorMetric[0].size(); l++){
      size_t NB = errorMetric[0][l].dimension();
      size_t OSize = NB*NB;
      for(auto j = 0ul; j < nExtrap; j++){
        for(auto k = 0ul; k <= j; k++){
          B(k,j) += blas::dot(OSize,errorMetric[k][l].pointer(),1,
                                 errorMetric[j][l].pointer(),1);
          /*
          if (errorMetric[k].hasZ())
            B(k,j) += blas::dot(OSize,errorMetric[k][l].Z().pointer(),1,
                                   errorMetric[j][l].Z().pointer(),1);
          if (errorMetric[k].hasXY()) {
            B(k,j) += blas::dot(OSize,errorMetric[k][l].Y().pointer(),1,
                                   errorMetric[j][l].Y().pointer(),1);
            B(k,j) += blas::dot(OSize,errorMetric[k][l].X().pointer(),1,
                                 errorMetric[j][l].X().pointer(),1);
           */
        }
      }
    }
    for(auto j = 0ul; j < nExtrap; j++){
      for(auto k = 0ul; k < j; k++){
         B(j,k) = B(k,j);
      }
    }
    for(auto l = 0ul; l < nExtrap; l++){
      B(nExtrap,l) = -1.0;
      B(l,nExtrap) = -1.0;
    }
    B(nExtrap,nExtrap) = 0.0;
  //prettyPrintSmart(std::cout,"B Matrix",&B[0],N,N,N);

    // Initialize LHS of the linear problem
    std::fill_n(&coeffs[0],N,0.);
    coeffs[nExtrap] = -1.0;
 
    int64_t* IPIV = CQMemManager::get().malloc<int64_t>(N);
    int INFO = lapack::gesv(N, 1, B.pointer(), N, IPIV, &coeffs[0], N);
    CQMemManager::get().free(IPIV);

//  for(auto i = 0ul; i < N; i++)
//    std::cout << "coeff = " << coeffs[i] << std::endl;

  //prettyPrintSmart(std::cout,"B Matrix Factored",&B[0],N,N,N);
    InvFail = (INFO != 0);
    return !InvFail;
  }; // function DIIS


}; // namespace ChronusQ

