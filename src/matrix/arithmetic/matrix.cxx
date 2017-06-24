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
#include <cqlinalg.hpp>

namespace ChronusQ {
namespace cqmatrix {

template <typename MatsT>
template <typename ScalarT, typename MatsU>
Matrix<MatsT>::Matrix( const ScaledMatrix<ScalarT, MatsU> &scaled ):
    Matrix(scaled.getScalarMatrix().nRows(), scaled.getScalarMatrix().nColumns()) {
  if (scaled.isPauli() and scaled.getPauliSpinorMatrices().hasZ())
    CErr("Cannot create a Matrix from a PauliSpinorMatrices with XYZ components.");
  SetMat('N',nRows(),nColumns(),scaled.scalar(),scaled.getScalarMatrix().pointer(),nRows(),pointer(),nRows());
}
template Matrix<double>::Matrix( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>::Matrix( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>::Matrix( const ScaledMatrix<dcomplex, double>& );
template Matrix<dcomplex>::Matrix( const ScaledMatrix<double, dcomplex>& );
template Matrix<dcomplex>::Matrix( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
Matrix<MatsT>&
Matrix<MatsT>::operator=( const ScaledMatrix<ScalarT, MatsU> &scaled ) {
  if (not isSameDimension(scaled.getScalarMatrix()))
    CErr("Cannot assign Matrix of different size.");
  if (std::is_same<MatsT, double>::value and
      std::is_same<MatsU, dcomplex>::value)
    CErr("Cannot assign a complex Matrix object to a real one.");
  if (scaled.isPauli() and scaled.hasZ())
    CErr("Cannot assign a Matrix from a PauliSpinorMatrices with XYZ components.");
  SetMat('N',nRows(),nColumns(),scaled.scalar(),scaled.getScalarMatrix().pointer(),nRows(),pointer(),nRows());
  return *this;
}
template Matrix<double>&
Matrix<double>::operator=( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>&
Matrix<dcomplex>::operator=( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>&
Matrix<dcomplex>::operator=( const ScaledMatrix<dcomplex, double>& );
template Matrix<dcomplex>&
Matrix<dcomplex>::operator=( const ScaledMatrix<double, dcomplex>& );
template Matrix<dcomplex>&
Matrix<dcomplex>::operator=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
Matrix<MatsT>& Matrix<MatsT>::operator=( const Matrix<MatsT> &other ) {
  if (this != &other) { // self-assignment check expected
    if (not isSameDimension(other))
      CErr("Cannot assign Matrix of different size.");
    std::copy_n(other.ptr_, nRow_ * nCol_, ptr_);
  }
  return *this;
}
template Matrix<double>& Matrix<double>::operator=( const Matrix<double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator=( const Matrix<dcomplex>& );

template <typename MatsT>
Matrix<MatsT>& Matrix<MatsT>::operator=( Matrix<MatsT> &&other ) {
  if (this != &other) { // self-assignment check expected
    if (not isSameDimension(other))
      CErr("Cannot assign Matrix of different size.");
    CQMemManager::get().free(ptr_);
    ptr_ = other.ptr_;
    other.ptr_ = nullptr;
  }
  return *this;
}
template Matrix<double>& Matrix<double>::operator=( Matrix<double>&& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator=( Matrix<dcomplex>&& );

template <typename MatsT>
Matrix<MatsT>& Matrix<MatsT>::operator*=( MatsT scalar ) {
  blas::scal(nRows() * nColumns(), scalar, pointer(), 1);
  return *this;
}
template Matrix<double>& Matrix<double>::operator*=( double scalar );
template Matrix<dcomplex>& Matrix<dcomplex>::operator*=( dcomplex scalar );

template <typename MatsT>
template <typename MatsU>
Matrix<MatsT>& Matrix<MatsT>::operator+=( const Matrix<MatsU> &other ) {
  if (not isSameDimension(other))
    CErr("Cannot add two Matrix of different size.");
  if (std::is_same<MatsT, double>::value and
      std::is_same<MatsU, dcomplex>::value)
    CErr("Cannot assign a complex Matrix object to a real one.");
  blas::axpy(nRows() * nColumns(), 1.0, other.pointer(), 1, pointer(), 1);
  return *this;
}
template Matrix<double>& Matrix<double>::operator+=( const Matrix<double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const Matrix<double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const Matrix<dcomplex>& );

// operator+= for Eigen::Matrix
template <typename MatsT>
template <typename MatsU>
Matrix<MatsT>& Matrix<MatsT>::operator+=(const Eigen::Matrix<MatsU, Eigen::Dynamic, Eigen::Dynamic>& eigen_mat) {


  (*this) += Matrix(eigen_mat);
  return *this;
}
template Matrix<double>&   Matrix<double>::operator+=(   const Eigen::Matrix<double,   Eigen::Dynamic, Eigen::Dynamic>&);
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const Eigen::Matrix<double,   Eigen::Dynamic, Eigen::Dynamic>&);
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const Eigen::Matrix<dcomplex, Eigen::Dynamic, Eigen::Dynamic>&);

template <typename MatsT>
template <typename MatsU>
Matrix<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type> Matrix<MatsT>::operator+( const Matrix<MatsU> &other ) const {
  Matrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> sum(*this);
  sum += other;
  return sum;
}
template Matrix<double> Matrix<double>::operator+( const Matrix<double>& ) const;
template Matrix<dcomplex> Matrix<dcomplex>::operator+( const Matrix<double>& ) const;
template Matrix<dcomplex> Matrix<double>::operator+( const Matrix<dcomplex>& ) const;
template Matrix<dcomplex> Matrix<dcomplex>::operator+( const Matrix<dcomplex>& ) const;

template <typename MatsT>
template <typename MatsU>
Matrix<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type> Matrix<MatsT>::operator-( const Matrix<MatsU> &other ) const {
  return *this + (-other);
}
template Matrix<double> Matrix<double>::operator-( const Matrix<double>& ) const;
template Matrix<dcomplex> Matrix<dcomplex>::operator-( const Matrix<double>& ) const;
template Matrix<dcomplex> Matrix<double>::operator-( const Matrix<dcomplex>& ) const;
template Matrix<dcomplex> Matrix<dcomplex>::operator-( const Matrix<dcomplex>& ) const;

template <typename MatsT>
template <typename ScalarT, typename MatsU>
Matrix<MatsT>& Matrix<MatsT>::operator+=( const ScaledMatrix<ScalarT, MatsU> &scaled ) {
  if (not isSameDimension(scaled.getScalarMatrix()))
    CErr("Cannot add two Matrix of different size.");
  if (scaled.isPauli() and scaled.hasZ())
    CErr("Cannot assign a Matrix from a PauliSpinorMatrices with XYZ components.");
  blas::axpy(nRows() * nColumns(), scaled.scalar(), scaled.getScalarMatrix().pointer(), 1, pointer(), 1);
  return *this;
}
template Matrix<double>& Matrix<double>::operator+=( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const ScaledMatrix<dcomplex, double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const ScaledMatrix<double, dcomplex>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator+=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
Matrix<MatsT>& Matrix<MatsT>::operator-=( const ScaledMatrix<ScalarT, MatsU> &scaled ) {
  return operator+=(-scaled);
}
template Matrix<double>& Matrix<double>::operator-=( const ScaledMatrix<double, double>& );
template <> template <> Matrix<double>& Matrix<double>::operator-=( const ScaledMatrix<double, dcomplex>& ) {
  CErr("Cannot subtract a complex ScaledMatrix from a real Matrix.");
  return *this;
}
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const ScaledMatrix<double, double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const ScaledMatrix<dcomplex, double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const ScaledMatrix<double, dcomplex>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type>
Matrix<MatsT>::operator+( const ScaledMatrix<ScalarT, MatsU> &scaled ) const {
  PauliSpinorMatrices<typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<MatsU, dcomplex>::value or
       std::is_same<ScalarT, dcomplex>::value),
      dcomplex, double>::type> sum(*this);
  sum += scaled;
  return sum;
}
template PauliSpinorMatrices<double> Matrix<double>::operator+( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator+( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator+( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator+( const ScaledMatrix<dcomplex, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator+( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator+( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator+( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator+( const ScaledMatrix<dcomplex, dcomplex>& ) const;

template PauliSpinorMatrices<double>
operator+( const ScaledMatrix<double, double>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, double>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, dcomplex>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, dcomplex>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, double>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, double>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, dcomplex>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, dcomplex>&, const Matrix<dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type>
Matrix<MatsT>::operator-( const ScaledMatrix<ScalarT, MatsU> &scaled ) const {
  return *this + (-scaled);
}
template PauliSpinorMatrices<double> Matrix<double>::operator-( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator-( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator-( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::operator-( const ScaledMatrix<dcomplex, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator-( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator-( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator-( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::operator-( const ScaledMatrix<dcomplex, dcomplex>& ) const;

template PauliSpinorMatrices<double>
operator-( const ScaledMatrix<double, double>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, double>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, dcomplex>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, dcomplex>&, const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, double>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, double>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, dcomplex>&, const Matrix<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, dcomplex>&, const Matrix<dcomplex>& );

template <typename MatsT>
template <typename MatsU>
Matrix<MatsT>& Matrix<MatsT>::operator-=( const Matrix<MatsU> &other ) {
  return operator+=(-other);
}
template Matrix<double>& Matrix<double>::operator-=( const Matrix<double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const Matrix<double>& );
template Matrix<dcomplex>& Matrix<dcomplex>::operator-=( const Matrix<dcomplex>& );

template <typename MatsT>
template <typename _FScale>
Matrix<MatsT> Matrix<MatsT>::scaleT(_FScale scale, char TRANS) const {
  size_t outNRow = nRows();
  size_t outNCol = nColumns();
  if (TRANS == 'T' or TRANS == 'C') std::swap(outNRow, outNCol);
  Matrix<MatsT> out(outNRow, outNCol);
  SetMat(TRANS, nRows(), nColumns(), scale, pointer(), nRows(), out.pointer(), outNRow);
  return out;
}
template Matrix<double>   Matrix<double>::scaleT(double scale, char TRANS) const;
template Matrix<dcomplex> Matrix<dcomplex>::scaleT(double scale, char TRANS) const;
template Matrix<dcomplex> Matrix<dcomplex>::scaleT(dcomplex scale, char TRANS) const;

template <typename MatsT>
template <typename _FScale>
void Matrix<MatsT>::inplace_scaleT(_FScale scale, char TRANS) {
  size_t LDA_trans = (TRANS == 'T' or TRANS == 'C') ? nColumns(): nRows();
  IMatCopy(TRANS, nRows(), nColumns(), scale, pointer(), nRows(), LDA_trans);
  if (TRANS == 'T' or TRANS == 'C') resize(nColumns(), nRows());
  return;  
}
template void Matrix<double>::inplace_scaleT(double scale, char TRANS);
template void Matrix<dcomplex>::inplace_scaleT(double scale, char TRANS);
template void Matrix<dcomplex>::inplace_scaleT(dcomplex scale, char TRANS);

template <typename MatsT>
void Matrix<MatsT>::setTriangle(blas::Uplo upLo, MatsT value, bool setDiag, MatsT diagValue) {
  
  if (not this->isSquareMatrix()) CErr("setTriangle only supported for square matrix");
  size_t N_ = nRow_;
  switch (upLo) {
    case blas::Uplo::Upper:
    for (size_t j = 1; j < N_; j++)
      for (size_t i = 0; i < j; i++)
        operator()(i,j) = value;
    break;
  case blas::Uplo::Lower:
    for (size_t j = 0; j < N_; j++)
      for (size_t i = j + 1; i < N_; i++)
        operator()(i,j) = value;
    break;
    default:
      break;
  }

  if (setDiag)
    for (size_t i = 0; i < N_; i++)
      operator()(i,i) = diagValue;
}
template void Matrix<double>::setTriangle(blas::Uplo upLo, double value, bool setDiag, double diagValue);
template void Matrix<dcomplex>::setTriangle(blas::Uplo upLo, dcomplex value, bool setDiag, dcomplex diagValue);






} // namespace cqmatrix
} // namespace ChronusQ
