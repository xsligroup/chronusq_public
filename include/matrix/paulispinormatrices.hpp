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
#include <matrix.hpp>

namespace ChronusQ {
namespace cqmatrix {

enum class PAULI_SPINOR_COMPS : size_t {
  S = 0, Z = 1, Y = 2, X = 3
};

/**
 *  \brief Templated class to handle the storage of
 *  one-electron two-spinor integral matrices
 *  Scalar, X, Y, and Z in Pauli representation.
 *
 */
template <typename MatsT>
class PauliSpinorMatrices {

  template <typename MatsU>
  friend class PauliSpinorMatrices;
  template <typename MatsU>
  friend class Matrix;

protected:
  std::vector<Matrix<MatsT>> components_;

public:

  // Constructor
  PauliSpinorMatrices() = delete;
  PauliSpinorMatrices( const PauliSpinorMatrices & ) = default;
  PauliSpinorMatrices( PauliSpinorMatrices && ) = default;
  /**
   *  \brief Constructor for PauliSpinorMatrices
   *         If memory provided, it will be split evenly into components
   *
   *  \param nRow Number of rows
   *  \param nCol Number of columns
   *  \param hasXY Flag for XY components
   *  \param hasZ Flag for Z component
   */
  PauliSpinorMatrices(size_t nRow, size_t nCol,
      bool hasXY = true, bool hasZ = true) {
    size_t nComp = hasXY ? 4 : hasZ ? 2 : 1;
    components_.reserve(nComp);
    for (size_t i = 0; i < nComp; i++)
      components_.emplace_back(nRow, nCol);
  }
  PauliSpinorMatrices(size_t n, bool hasXY = true, bool hasZ = true):
      PauliSpinorMatrices(n, n, hasXY, hasZ) { }
  template <typename MatsU>
  PauliSpinorMatrices(const Matrix<MatsU> &other,
                      bool addXY = false, bool addZ = false ) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real PauliSpinorMatrices from a Complex one.");
    size_t nComp = addXY ? 4 : addZ ? 2 : 1;
    components_.reserve(nComp--);
    components_.emplace_back(other, 0);
    while (nComp > 0) {
      components_.emplace_back(this->nRows(), this->nColumns());
      components_.back().clear();
      nComp--;
    }
  }
  template <typename MatsU>
  PauliSpinorMatrices(const PauliSpinorMatrices<MatsU> &other,
                      bool addXY = false, bool addZ = false ) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real PauliSpinorMatrices from a Complex one.");
    size_t nComp = std::max(addXY ? 4ul : addZ ? 2ul : 1ul, other.components_.size());
    components_.reserve(nComp);
    for (auto &p : other.components_)
      components_.emplace_back(p);
    while (nComp > other.components_.size()) {
      components_.emplace_back(this->nRows(), this->nColumns());
      components_.back().clear();
      nComp--;
    }
  }
  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices( const ScaledMatrix<ScalarT, MatsU>&,
                             bool addXY = false, bool addZ = false );

  PauliSpinorMatrices& operator=( const PauliSpinorMatrices &other );
  PauliSpinorMatrices& operator=( PauliSpinorMatrices &&other );
  template <typename MatsU>
  PauliSpinorMatrices& operator=( const Matrix<MatsU> &other ) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real PauliSpinorMatrices from a Complex Matrix.");
    if (not S().isSameDimension(other))
      CErr("Cannot assign a Matrix of different size to PauliSpinorMatrices.");
    S() = other;
    for (size_t i = 1; i < components_.size(); i++) components_[i].clear();
    return *this;
  }
  template <typename MatsU>
  PauliSpinorMatrices& operator=( Matrix<MatsU> &&other ) {
    if (std::is_same<MatsU, dcomplex>::value
        and std::is_same<MatsT, double>::value)
      CErr("Cannot create a Real PauliSpinorMatrices from a Complex Matrix.");
    if (not S().isSameDimension(other))
      CErr("Cannot assign a Matrix of different size to PauliSpinorMatrices.");
    S() = std::move(other);
    for (size_t i = 1; i < components_.size(); i++) components_[i].clear();
    return *this;
  }
  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices& operator=( const ScaledMatrix<ScalarT, MatsU>& );

  PauliSpinorMatrices& operator*=( MatsT );

  ScaledMatrix<double, MatsT> operator-() const {
    return ScaledMatrix<double, MatsT>(-1.0, *this);
  }

  template <typename MatsU>
  PauliSpinorMatrices& operator+=( const PauliSpinorMatrices<MatsU>& );
  template <typename MatsU>
  PauliSpinorMatrices& operator-=( const PauliSpinorMatrices<MatsU>& );
  template <typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator+( const PauliSpinorMatrices<MatsU>& ) const;
  template <typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator-( const PauliSpinorMatrices<MatsU>& ) const;

  template <typename MatsU>
  PauliSpinorMatrices& operator+=( const Matrix<MatsU>& );
  template <typename MatsU>
  PauliSpinorMatrices& operator-=( const Matrix<MatsU>& );
  template <typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator+( const Matrix<MatsU>& ) const;
  template <typename MatsU>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator-( const Matrix<MatsU>& ) const;

  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices& operator+=( const ScaledMatrix<ScalarT, MatsU>& );
  template <typename ScalarT, typename MatsU>
  PauliSpinorMatrices& operator-=( const ScaledMatrix<ScalarT, MatsU>& );
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

  bool hasXY() const { return components_.size() == 4; }
  bool hasZ() const { return components_.size() >= 2; }
  size_t nComponent() const { return components_.size(); }

  size_t dimension() const{ return S().dimension(); }
  size_t nColumns() const { return S().nColumns(); }
  size_t nRows() const { return S().nRows(); }

  bool isSquareMatrix() const { return S().isSquareMatrix(); }
  
  void resize(size_t nRow, size_t nCol) {
    S().resize(nRow, nCol);
    if (hasZ()) {
      Z().resize(nRow, nCol);
    }
    if (hasXY()) {
      X().resize(nRow, nCol);
      Y().resize(nRow, nCol);
    }
  }

  Matrix<MatsT>& operator[](PAULI_SPINOR_COMPS comp) {
    size_t i = static_cast<size_t>(comp);
    if (i >= components_.size())
      CErr("Requested component is NOT in this PauliSpinorMatrices object.");
    return components_[i];
  }

  const Matrix<MatsT>& operator[](PAULI_SPINOR_COMPS comp) const {
    size_t i = static_cast<size_t>(comp);
    if (i >= components_.size())
      CErr("Requested component is NOT in this PauliSpinorMatrices object.");
    return components_[i];
  }

  Matrix<MatsT>& S() { return operator[](PAULI_SPINOR_COMPS::S); }
  const Matrix<MatsT>& S() const { return operator[](PAULI_SPINOR_COMPS::S); }
  Matrix<MatsT>& X() { return operator[](PAULI_SPINOR_COMPS::X); }
  const Matrix<MatsT>& X() const { return operator[](PAULI_SPINOR_COMPS::X); }
  Matrix<MatsT>& Y() { return operator[](PAULI_SPINOR_COMPS::Y); }
  const Matrix<MatsT>& Y() const { return operator[](PAULI_SPINOR_COMPS::Y); }
  Matrix<MatsT>& Z() { return operator[](PAULI_SPINOR_COMPS::Z); }
  const Matrix<MatsT>& Z() const { return operator[](PAULI_SPINOR_COMPS::Z); }

  std::vector<MatsT*> SZYXPointers() {
    if (hasXY())
      return { S().pointer(), Z().pointer(), Y().pointer(), X().pointer() };
    if (hasZ())
      return { S().pointer(), Z().pointer() };
    return { S().pointer() };
  }
  
  PauliSpinorMatrices<double> real_part() {
    PauliSpinorMatrices<double> realMats(S().real_part(), false, false);
    for (size_t i = 1; i < components_.size(); i++)
      realMats.components_.emplace_back(components_[i].real_part());
    return realMats;
  }

  void clear() {
    for (Matrix<MatsT>& c : components_)
      c.clear();
  }

  void output(std::ostream &out, const std::string &s = "",
                      bool printFull = false) const {
    if (printFull) {
      std::string oeiStr;
      if (s == "")
        oeiStr = "PauliSpinorOEI";
      else
        oeiStr = "PauliSpinorOEI[" + s + "]";
      prettyPrintSmart(out, oeiStr+".S", S().pointer(),
          this->nRows(), this->nColumns(), this->nRows());
      if(hasZ())
        prettyPrintSmart(out, oeiStr+".Z", Z().pointer(),
            this->nRows(), this->nColumns(), this->nRows());
      if(hasXY()) {
        prettyPrintSmart(out, oeiStr+".Y", Y().pointer(),
            this->nRows(), this->nColumns(), this->nRows());
        prettyPrintSmart(out, oeiStr+".X", X().pointer(),
            this->nRows(), this->nColumns(), this->nRows());
      }
    } else {
      std::string oeiStr;
      if (s == "")
        oeiStr = "Pauli spinor one-electron integrals";
      else
        oeiStr = "PauliSpinorOEI[" + s + "]";
      out << oeiStr;
      if(not hasZ())
        out << " without XYZ components";
      else if(not hasXY())
        out << " without XY components";
      out << std::endl;
    }
  }

  void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
    for (Matrix<MatsT>& comp : components_)
      comp.broadcast(comm, root);
  }
  
  template <typename MatsU>
  void spinGather(Matrix<MatsU>& mat) const;

  template <typename MatsU>
  Matrix<MatsU> spinGather() const;

  template <typename MatsU>
  std::vector<Matrix<MatsU>> spinGatherToBlocks(
      bool genABBA = true, bool genBB = true) const;

  template <typename MatsU>
  static PauliSpinorMatrices<MatsT>
  spinBlockScatterBuild(const Matrix<MatsU> &AA, bool hasXY = false, bool hasZ = false);

  template <typename MatsU>
  static PauliSpinorMatrices<MatsT>
  spinBlockScatterBuild(const Matrix<MatsU> &AA, const Matrix<MatsU> &BB,
                        bool hasXY = false, bool hasZ = true);

  template <typename MatsU>
  static PauliSpinorMatrices<MatsT>
  spinBlockScatterBuild(const Matrix<MatsU> &AA, const Matrix<MatsU> &AB,
                        const Matrix<MatsU> &BA, const Matrix<MatsU> &BB,
                        bool hasXY = true, bool hasZ = true);
  
  template <typename MatsU>
  void componentScatter(PauliSpinorMatrices<MatsU> & LL, 
                        PauliSpinorMatrices<MatsU> & LS,
                        PauliSpinorMatrices<MatsU> & SL,
                        PauliSpinorMatrices<MatsU> & SS,
                        bool increment = false) const; 
  
  template <typename MatsU>
  void componentGather(const PauliSpinorMatrices<MatsU> & LL, 
                       const PauliSpinorMatrices<MatsU> & LS,
                       const PauliSpinorMatrices<MatsU> & SL,
                       const PauliSpinorMatrices<MatsU> & SS,
                       bool increment = false); 

  template <typename MatsU>
  static PauliSpinorMatrices<MatsT> 
  componentGatherBuild(const PauliSpinorMatrices<MatsU> & LL, 
                       const PauliSpinorMatrices<MatsU> & LS,
                       const PauliSpinorMatrices<MatsU> & SL,
                       const PauliSpinorMatrices<MatsU> & SS); 
  
  template <typename MatsU>
  void componentAdd(const char TRANS, MatsU scale, const std::string & comp,
    const PauliSpinorMatrices<MatsU> & pauli);

  void symmetrizeLSSL(char TRANS, bool get_SL_from_LS = true);
  
  template <typename TransT>
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<TransT, dcomplex>::value),
  dcomplex, double>::type> transform(
      char TRANS, const TransT* T, int NT, int LDT) const {
    
    if (not this->isSquareMatrix()) CErr("transform only supported for square matrix");
          
    PauliSpinorMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transMats(NT, hasXY(), hasZ());
    transMats.S() = S().transform(TRANS, T, NT, LDT);
    if (hasZ())
      transMats.Z() = Z().transform(TRANS, T, NT, LDT);
    if (hasXY()) {
      transMats.Y() = Y().transform(TRANS, T, NT, LDT);
      transMats.X() = X().transform(TRANS, T, NT, LDT);
    }
    return transMats;
  }

  double norm(lapack::Norm type) const {
    if(type == lapack::Norm::Fro)
      CErr("There is a possible error in PauliSpinorMatrix::Norm for Frobenius norms.");
    double result = S().norm(type);
    if (hasZ()) result = std::max(result, Z().norm(type));
    if (hasXY()) {
       result = std::max(result, Y().norm(type));
       result = std::max(result, X().norm(type));
    }
    return result; 
  }

  // Check if the matrix has NaN
  virtual bool hasNaN() const {
    if (S().hasNaN()) return true;
    if (hasZ() and Z().hasNaN()) return true;
    if (hasXY() and (Y().hasNaN() or X().hasNaN())) return true;
    return false;
  }
  
  ~PauliSpinorMatrices() {}

}; // class PauliSpinorMatrices

template <typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type> operator+(
    const Matrix<MatsT> &lhs, const PauliSpinorMatrices<MatsU> &rhs ) {
  return rhs + lhs;
}

template <typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type> operator-(
    const Matrix<MatsT> &lhs, const PauliSpinorMatrices<MatsU> &rhs ) {
  return lhs + (-rhs);
}

template <typename ScalarT, typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type> operator+(
    const ScaledMatrix<ScalarT, MatsT> &lhs, const PauliSpinorMatrices<MatsU> &rhs ) {
  return rhs + lhs;
}

template <typename ScalarT, typename MatsT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type> operator-(
    const ScaledMatrix<ScalarT, MatsT> &lhs, const PauliSpinorMatrices<MatsU> &rhs ) {
  return lhs + (-rhs);
}

} // namespace cqmatrix
} // namespace ChronusQ
