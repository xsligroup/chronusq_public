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

#include <cqlinalg/svd.hpp>
#include <cqlinalg/util.hpp>
#include <cerr.hpp>
#include <lapack.hh>
#include <limits>

namespace ChronusQ {


  template <typename T>
  struct real_type {
    using type = T;
  };

  template <typename T>
  struct real_type<std::complex<T>> {
    using type = T;
  };

  template <typename _F>
  size_t ORTH(int M, int N, _F *A, int LDA, double *S,
    _F *U, int LDU) {

    if (not U) U = A;

    int INFO = lapack::gesvd(A == U ? lapack::Job::OverwriteVec : lapack::Job::SomeVec, 
                 lapack::Job::NoVec, M, N, A, LDA, S, U, LDU, reinterpret_cast<_F*>(NULL),LDA);

    if (INFO < 0)
      CErr(std::to_string(-INFO) + "-th argument had an illegal value.");
    else if (INFO > 0)
      CErr("SVD did not converge.");

    double tol = std::max(M,N) * S[0]
      * std::numeric_limits<typename real_type<_F>::type>::epsilon();

    return std::distance(S, std::find_if(S,S+std::min(M,N),
      [tol](double x){ return x < tol; }));

  }; // ORTH

  template size_t ORTH(int, int, double*, int, double*, double*, int);
  template size_t ORTH(int, int, dcomplex*, int, double*, dcomplex*, int);

}; // namespace ChronusQ

