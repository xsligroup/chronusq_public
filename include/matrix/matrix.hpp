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
#include <memmanager.hpp>
#include <cerr.hpp>
#include <util/matout.hpp>
#include <cqlinalg.hpp>

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
  size_t nRow_;
  size_t nCol_;
  MatsT *ptr_ = nullptr;     ///< Raw matrix storage (2 index)

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
      nRow_(nRow), nCol_(nCol) {
    malloc();
  }
  Matrix(size_t n):
      Matrix(n, n) { }

  Matrix( const Matrix &other ):
      Matrix(other.nRow_, other.nCol_) {
    std::copy_n(other.ptr_, nRow_ * nCol_, ptr_);
  }
  template <typename MatsU>
  Matrix( const Matrix<MatsU> &other, int = 0 ):
      Matrix(other.nRow_, other.nCol_) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real Matrix from a Complex one.");
    std::copy_n(other.ptr_, nRow_ * nCol_, ptr_);
  }
  Matrix( Matrix &&other ):
      nRow_(other.nRow_), nCol_(other.nCol_),
      ptr_(other.ptr_) { other.ptr_ = nullptr; }
  template <typename MatsU>
  Matrix( const PauliSpinorMatrices<MatsU> &other ):
      Matrix(other.nRows(), other.nColumns()) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real Matrix from a Complex one.");
    if (other.hasZ())
      CErr("Cannot create a Matrix from a PauliSpinorMatrices"
           " with XYZ components.");
    std::copy_n(other.S().pointer(), nRow_ * nCol_, ptr_);
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
    std::copy_n(eigen_mat.data(), eigen_mat.rows()*eigen_mat.cols(), ptr_);
  }
  template <typename MatsU>
  Matrix( const Eigen::Matrix<MatsU, Eigen::Dynamic, Eigen::Dynamic>& eigen_mat ): 
      Matrix(eigen_mat.rows(), eigen_mat.cols()) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real Matrix from a Complex Eigen matrix.");
    std::copy_n(eigen_mat.data(), eigen_mat.rows()*eigen_mat.cols(), ptr_);
  }
  
  Matrix& operator=( const Matrix &other );
  Matrix& operator=( Matrix &&other );

  template <typename ScalarT, typename MatsU>
  Matrix& operator=( const ScaledMatrix<ScalarT, MatsU>& );

  size_t dimension() const{ return nRow_; }
  size_t nColumns() const { return nCol_; }
  size_t nRows() const { return nRow_; }
   
  bool isSquareMatrix() const { return nRow_ == nCol_; }

  template <typename MatU>
  bool isSameDimension(const Matrix<MatU>& other) const {
    return nRow_ == other.nRows() and nCol_ == other.nColumns();
  }
  
  void resize(size_t nRow, size_t nCol) {
    if (nRow * nCol != nRow_ * nCol_) {
      nRow_ = nRow;
      nCol_ = nCol;
      #pragma omp critical
      {
        malloc();
      }
    } else {
      nRow_ = nRow;
      nCol_ = nCol;
    }
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
    return ptr_[p + q * nRow_];
  }
  MatsT operator()(size_t p, size_t q) const {
    return ptr_[p + q * nRow_];
  }

  // Matrix direct access
  MatsT* pointer() { return ptr_; }
  const MatsT* pointer() const { return ptr_; }

  Matrix<double> real_part() {
    Matrix<double> realMat(nRow_, nCol_);
    GetMatRE('N', nRow_, nCol_, 1., pointer(), nRow_, realMat.pointer(), nRow_);
    return realMat;
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
    std::fill_n(ptr_, nRow_ * nCol_, MatsT(0.));
  }

  void output(std::ostream &out, const std::string &s = "",
                      bool printFull = false) const {
    std::string matStr;
    if (s == "")
      matStr = "Square Matrix";
    else
      matStr = "Square Matrix[" + s + "]";
    if (printFull)
      prettyPrintSmart(out, matStr, pointer(), nRow_, nCol_, nRow_);
    else {
      out << matStr << std::endl;
    }
  }

  void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

#ifdef CQ_ENABLE_MPI
    // BCast matrix to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering a matrix ***\n";
      size_t nRow_bcast = nRow_;
      size_t nCol_bcast = nCol_;
      MPIBCast(nRow_bcast,root,comm);
      MPIBCast(nCol_bcast,root,comm);

      if (nRow_bcast != nRow_ or nCol_bcast != nCol_) {
        nRow_ = nRow_bcast;
        nCol_ = nCol_bcast;
        malloc();
      }

      MPIBCast(ptr_, nRow_ * nCol_, root, comm);
    }
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
    return lapack::lange(norm, nRow_, nCol_, ptr_, nRow_);
  }

  virtual bool hasNaN() const {
    for (size_t i = 0; i < nRow_ * nCol_; i++) {
      if (std::isnan(std::real(ptr_[i]))) return true;
      if (std::isnan(std::imag(ptr_[i]))) return true;
    }
    return false;
  }

  void malloc() {
    if (ptr_) CQMemManager::get().free(ptr_);
    size_t N = nRow_ * nCol_;
    if (N != 0) {
      try { ptr_ = CQMemManager::get().malloc<MatsT>(N); }
      catch(...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for the full INTS matrix ("
                  << (N /1e9) * sizeof(double) << " GB)" << std::endl;
        std::cout << std::endl << CQMemManager::get() << std::endl;
        throw std::bad_alloc();
      }
    }
  }

  ~Matrix() {
    if(ptr_) CQMemManager::get().free(ptr_);
  }

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
template <typename MatsT>
bool hasNaN(MatsT * ptr, size_t N) {
  for (size_t i = 0; i < N; i++) {
    if (std::isnan(std::real(ptr[i]))) return true;
    if (std::isnan(std::imag(ptr[i]))) return true;
  }
  return false;
}
} // namespace ChronusQ
