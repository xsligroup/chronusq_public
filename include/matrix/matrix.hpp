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
#include <matrix/ndarray.hpp>
#include <util/matout.hpp>

namespace ChronusQ {
namespace cqmatrix {

template <typename ScalarT, typename MatsT>
class ScaledMatrix;

template <typename MatsT>
class PauliSpinorMatrices;

template <typename MatsT>
class Matrix {

  template <typename MatsU>
  friend class Matrix;

protected:
  std::shared_ptr<NDArray<MatsT>> array_ = nullptr;     ///< Raw NDArray storage (2 index)

public:

  // Constructor
  Matrix() = delete;
  /**
   * @brief Construct a new Matrix object
   *
   * @param nRow Number of rows
   * @param nCol Number of columns
   */
  Matrix(size_t nRow, size_t nCol):
      array_(std::make_shared<NDArray<MatsT>>(std::vector<size_t>{nRow, nCol})) {}
  Matrix(const std::vector<size_t> &dims):
      array_(std::make_shared<NDArray<MatsT>>(dims)) {
    if (dims.size() != 2)
      CErr("Cannot create a Matrix with more or less than 2 dimensions.");
  }
  Matrix(std::shared_ptr<NDArray<MatsT>> array):
      array_(array) {}
  Matrix(const NDArray<MatsT> &ndarray):
      array_(std::make_shared<NDArray<MatsT>>(ndarray)) {
    if (not array_->isMatrix())
      CErr("Cannot create a Matrix with more or less than 2 dimensions.");
  }
  Matrix(NDArray<MatsT> &&ndarray):
      array_(std::make_shared<NDArray<MatsT>>(std::move(ndarray))) {
    if (not array_->isMatrix())
      CErr("Cannot create a Matrix with more or less than 2 dimensions.");
  }
  Matrix(size_t n):
      Matrix(n, n) {}

  Matrix( const Matrix &other ):
      Matrix(*other.array_){}
  Matrix( Matrix &&other ) = default;
  template <typename MatsU>
  Matrix( const Matrix<MatsU> &other, int = 0 ):
      array_(std::make_shared<NDArray<MatsT>>(*other.array_, 0)) {}
  template <typename MatsU>
  Matrix( const PauliSpinorMatrices<MatsU> &other ):
      Matrix(other.S(), 0) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real Matrix from a Complex one.");
    if (other.hasZ())
      CErr("Cannot create a Matrix from a PauliSpinorMatrices"
           " with XYZ components.");
  }
  Matrix( PauliSpinorMatrices<MatsT> &&other ):
      Matrix(std::move(other.S())) {
    if (other.hasZ())
      CErr("Cannot create a Matrix from a PauliSpinorMatrices"
           " with XYZ components.");
  }
  template <typename ScalarT, typename MatsU>
  Matrix( const ScaledMatrix<ScalarT, MatsU>& );

  // constructors that take an Eigen::Matrix
  Matrix( const Eigen::Matrix<MatsT, Eigen::Dynamic, Eigen::Dynamic>& eigen_mat ):
      Matrix(eigen_mat.rows(), eigen_mat.cols()) {
    std::copy_n(eigen_mat.data(), eigen_mat.rows()*eigen_mat.cols(), pointer());
  }
  template <typename MatsU>
  Matrix( const Eigen::Matrix<MatsU, Eigen::Dynamic, Eigen::Dynamic>& eigen_mat ):
      Matrix(eigen_mat.rows(), eigen_mat.cols()) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real Matrix from a Complex Eigen matrix.");
    std::copy_n(eigen_mat.data(), eigen_mat.rows()*eigen_mat.cols(), pointer());
  }

  Matrix& operator=( const Matrix &other );
  Matrix& operator=( Matrix &&other );

  template <typename ScalarT, typename MatsU>
  Matrix& operator=( const ScaledMatrix<ScalarT, MatsU>& );

  const std::vector<size_t>& dimensions() const{ return array_->dimensions(); }
  size_t nRows() const { return dimensions()[0]; }
  size_t nColumns() const { return dimensions()[1]; }

  std::shared_ptr<NDArray<MatsT>> getNDArray() const { return array_; }

  bool isSquareMatrix() const { return nRows() == nColumns(); }

  template <typename MatU>
  bool isSameDimension(const Matrix<MatU>& other) const {
    return nRows() == other.nRows() and nColumns() == other.nColumns();
  }
  
  void resize(size_t nRow, size_t nCol) {
    array_->resize({nRow, nCol});
  }

  void swap(Matrix& other) noexcept {
    array_.swap(other.array_);
  }

  Matrix& operator*=( MatsT );
  ScaledMatrix<double, MatsT> operator-() const {
    return ScaledMatrix<double, MatsT>(-1.0, *this);
  }

  template <typename MatsU>
  Matrix& operator+=( const Matrix<MatsU>& );
  template <typename MatsU>
  Matrix& operator+=(const Eigen::Matrix<MatsU, Eigen::Dynamic, Eigen::Dynamic>&);
  template <typename MatsU>
  Matrix& operator-=( const Matrix<MatsU>& );
  template <typename MatsU>
  Matrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator+( const Matrix<MatsU>& ) const;
  template <typename MatsU>
  Matrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator-( const Matrix<MatsU>& ) const;

  template <typename ScalarT, typename MatsU>
  Matrix& operator+=( const ScaledMatrix<ScalarT, MatsU>& );
  template <typename ScalarT, typename MatsU>
  Matrix& operator-=( const ScaledMatrix<ScalarT, MatsU>& );
  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  operator+( const ScaledMatrix<ScalarT, MatsU>& ) const;
  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  operator-( const ScaledMatrix<ScalarT, MatsU>& ) const;

  MatsT& operator()(size_t p, size_t q) {
    return (*array_)(p,q);
  }
  const MatsT& operator()(size_t p, size_t q) const {
    return (*array_)(p,q);
  }

  // Matrix direct access
  MatsT* pointer() { return array_->pointer(); }
  const MatsT* pointer() const { return array_->pointer(); }

  Matrix<double> real_part() {
    return array_->real_part();
  }

  Matrix<double> imag_part() {
    return array_->imag_part();
  }
  
  // transform and return the transformed matrix
  Matrix<MatsT> T(char TRANS = 'T') const { return scaleT(1.0, TRANS); }
  
  template <typename _FScale>
  Matrix<MatsT> scaleT(_FScale scale, char TRANS) const;
  
  void inplace_T(char TRANS = 'T') { inplace_scaleT(1.0, TRANS); }

  template <typename _FScale>
  void inplace_scaleT(_FScale scale, char TRANS);
  
  void setTriangle(blas::Uplo upLo, MatsT value, bool setDiag, MatsT diagValue = 1.0);
  
  void clear() {
    array_->clear();
  }

  void output(std::ostream &out, const std::string &s = "",
                      bool printFull = false) const {
    std::string matStr;
    if (s == "")
      matStr = "Matrix";
    else
      matStr = "Matrix[" + s + "]";
    if (printFull)
      prettyPrintSmart(out, matStr, pointer(), nRows(), nColumns(), nRows());
    else {
      out << matStr << std::endl;
    }
  }

  void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

#ifdef CQ_ENABLE_MPI
    // BCast matrix to all MPI processes
    array_->broadcast(comm, root);
#endif

  }

  template <typename MatsU>
  void spinScatter(PauliSpinorMatrices<MatsU>& pauli,
      bool hasXY = true, bool hasZ = true) const;
  
  template <typename MatsU>
  PauliSpinorMatrices<MatsU> spinScatter(
      bool hasXY = true, bool hasZ = true) const;

  template <typename MatsU>
  Matrix<MatsU> spatialToSpinBlock() const;
  
  template <typename MatsU>
  void componentScatter(Matrix<MatsU> & LL,
                        Matrix<MatsU> & LS,
                        Matrix<MatsU> & SL,
                        Matrix<MatsU> & SS,
                        bool increment = false) const;
   
  template <typename MatsU>
  void componentGather(const Matrix<MatsU> & LL,
                       const Matrix<MatsU> & LS,
                       const Matrix<MatsU> & SL,
                       const Matrix<MatsU> & SS,
                       bool increment = false);
   
  template <typename MatsU>
  static Matrix<MatsT>
  componentGatherBuild(const Matrix<MatsU> & LL,
                       const Matrix<MatsU> & LS,
                       const Matrix<MatsU> & SL,
                       const Matrix<MatsU> & SS);
  
  
  
  template <typename TransT>
  Matrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<TransT, dcomplex>::value),
  dcomplex, double>::type> transform(
      char TRANS, const TransT* T, int NT, int LDT) const;

  template <typename TransT, typename OutT>
  void subsetTransform(
      char TRANS, const TransT* T, int LDT,
      const std::vector<std::pair<size_t,size_t>> &off_size,
      OutT* out, bool increment = false) const;

  
  double norm(lapack::Norm norm) const {
    return lapack::lange(norm, nRows(), nColumns(), pointer(), nRows());
  }

  virtual bool hasNaN() const {
    return array_->hasNaN();
  }

  ~Matrix() {}

}; // class Matrix

template <typename _F1, typename _F2, typename _FScale1>
void MatrixAXPY(char TRANS, _FScale1 ALPHA, const Matrix<_F1>& X, Matrix<_F2>& Y) {
  size_t nOpXRow = X.nRows();
  size_t nYRow = Y.nRows();
  size_t nOpXCol = X.nColumns();
  size_t nYCol = Y.nColumns();
  if (TRANS == 'T' or TRANS == 'C') std::swap(nOpXRow, nOpXCol);
  if (nOpXRow != nYRow or nOpXCol != nYCol)
    CErr("X and Y must have matching dimensions in MatrixAXPY");
   
  MatAdd('N', TRANS, nYRow, nYCol, _F2(1.), Y.pointer(), Y.nRows(), ALPHA, X.pointer(), X.nRows(),
      Y.pointer(), Y.nRows());
}


template <typename ScalarT, typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type> operator+(
    const ScaledMatrix<ScalarT, MatsT> &lhs, const Matrix<MatsU> &rhs ) {
  return rhs + lhs;
}

template <typename ScalarT, typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type> operator-(
    const ScaledMatrix<ScalarT, MatsT> &lhs, const Matrix<MatsU> &rhs ) {
  return lhs + (-rhs);
}

template <typename MatsT>
std::ostream& operator<<(std::ostream&, const Matrix<MatsT>&);

} // namespace cqmatrix
} // namespace ChronusQ
