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

#include <matrix.hpp>

namespace ChronusQ {
namespace cqmatrix {

/**
 *  \brief < p |O| q> = T(mu, p)^H @ < mu |O| nu > @ T(nu, q)
 *
 *  \param [in]  TRANS     Whether transpose/adjoint T
 *  \param [in]  T         Transformation matrix
 *  \param [in]  LDT       Leading dimension of T
 *  \param [in]  off_sizes Vector of 2 pairs,
 *                         a pair of offset and size for each index.
 *  \param [out] out       Return the contraction result.
 *  \param [in]  increment Perform out += result if true
 */
template <typename MatsT>
template <typename TransT, typename OutT>
void Matrix<MatsT>::subsetTransform(
    char TRANS, const TransT* T, int LDT,
    const std::vector<std::pair<size_t,size_t>> &off_sizes,
    OutT* out, bool increment) const {
  typedef typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type ResultsT;
  
  if (not this->isSquareMatrix()) CErr("transform only supported for square matrix");
  
  size_t N_ = nRow_;
  ResultsT* SCR = CQMemManager::get().malloc<ResultsT>(N_ * off_sizes[0].second);
  std::fill_n(SCR, N_ * off_sizes[0].second, ResultsT(0.0));
  MatsT * dummy = nullptr;

  // SCR(nu, p) = < mu |O| nu >^H @ T(mu, p)
  // < p |O| q> = SCR(nu, p)^H @ T(nu, q)
  PairTransformation(TRANS, T, LDT, off_sizes[0].first, off_sizes[1].first,
    'N', pointer(), N_, N_, 1, 'T', out, off_sizes[0].second, off_sizes[1].second,
    dummy, SCR, increment); 
  
  CQMemManager::get().free(SCR);
}
template void Matrix<double>::subsetTransform(
    char TRANS, const double* T, int LDT,
    const std::vector<std::pair<size_t,size_t>> &off_sizes,
    double* out, bool increment) const;
template void Matrix<double>::subsetTransform(
    char TRANS, const dcomplex* T, int LDT,
    const std::vector<std::pair<size_t,size_t>> &off_sizes,
    dcomplex* out, bool increment) const;
template void Matrix<dcomplex>::subsetTransform(
    char TRANS, const dcomplex* T, int LDT,
    const std::vector<std::pair<size_t,size_t>> &off_sizes,
    dcomplex* out, bool increment) const;
template void Matrix<dcomplex>::subsetTransform(
    char TRANS, const double* T, int LDT,
    const std::vector<std::pair<size_t,size_t>> &off_sizes,
    dcomplex* out, bool increment) const;


//
//  template <>
//  template <>
//  void Matrix<dcomplex>::subsetTransform(
//      char TRANS, const double* T, int LDT,
//      const std::vector<std::pair<size_t,size_t>> &off_sizes,
//      dcomplex* out, bool increment) const {
//    std::vector<size_t> offs;
//    for (const auto &off_size : off_sizes) {
//      if (TRANS == 'T' or TRANS == 'C')
//        offs.push_back(off_size.first);
//      else if (TRANS == 'N')
//        offs.push_back(off_size.first * LDT);
//    }
//    dcomplex* SCR = CQMemManager::get().malloc<dcomplex>(off_sizes[1].second * N_);
//    if (TRANS == 'T' or TRANS == 'C')
//      TRANS = 'N';
//    else if (TRANS == 'N')
//      TRANS = 'C';
//
//
//    blas::Op OP_TRANS;
//    if (TRANS == 'T') {
//      OP_TRANS = blas::Op::Trans;
//    } else if (TRANS == 'C') {
//      OP_TRANS = blas::Op::ConjTrans;
//    } else if (TRANS == 'N') {
//      OP_TRANS = blas::Op::NoTrans;
//    }
//
//
//    // SCR(q, mu) = T(nu, q)^H @ < mu |O| nu >^H
//    blas::gemm(blas::Layout::ColMajor,OP_TRANS, blas::Op::ConjTrans, off_sizes[1].second, N_, N_,
//        dcomplex(1.), T+offs[1], LDT, pointer(), N_,
//        dcomplex(0.), SCR, off_sizes[1].second);
//    // < p |O| q> = T(mu, p)^H @ SCR(q, mu)^H
//    //            = T(mu, p)^H @ < mu |O| nu > @ T(nu, q)
//    dcomplex outFactor = increment ? 1.0 : 0.0;
//    blas::gemm(blas::Layout::ColMajor,OP_TRANS, blas::Op::ConjTrans, off_sizes[0].second, off_sizes[1].second, N_,
//        dcomplex(1.), T+offs[0], LDT, SCR, off_sizes[1].second,
//        outFactor, out, off_sizes[0].second);
//    CQMemManager::get().free(SCR);
//  }

/**
 *  \brief < p |O| q> = T(mu, p)^H @ < mu |O| nu > @ T(nu, q)
 *
 *  \param [in] TRANS Whether transpose/adjoint T
 *  \param [in] T     Transformation matrix
 *  \param [in] NT    Number of columns for T
 *  \param [in] LDT   Leading dimension of T
 *
 *  \return Matrix object with element type derived from
 *          MatsT and TransT.
 */
template <typename MatsT>
template <typename TransT>
Matrix<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<TransT, dcomplex>::value),
dcomplex, double>::type> Matrix<MatsT>::transform(
    char TRANS, const TransT* T, int NT, int LDT) const {
  Matrix<typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(NT);
  
  if (not this->isSquareMatrix()) CErr("transform only supported for square matrix");
  transInts.clear();
  subsetTransform(TRANS,T,LDT,{{0,NT},{0,NT}},transInts.pointer(),false);
  return transInts;
}

template Matrix<double> Matrix<double>::transform(
    char TRANS, const double* T, int NT, int LDT) const;
template Matrix<dcomplex> Matrix<double>::transform(
    char TRANS, const dcomplex* T, int NT, int LDT) const;
template Matrix<dcomplex> Matrix<dcomplex>::transform(
    char TRANS, const double* T, int NT, int LDT) const;
template Matrix<dcomplex> Matrix<dcomplex>::transform(
    char TRANS, const dcomplex* T, int NT, int LDT) const;

} // namespace cqmatrix
} // namespace ChronusQ
