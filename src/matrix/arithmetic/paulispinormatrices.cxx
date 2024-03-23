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

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<MatsT>::PauliSpinorMatrices(
    const ScaledMatrix<ScalarT, MatsU> &scaled, bool addXY, bool addZ ):
    Matrix<MatsT>(scaled.matrix().memManager(), scaled.matrix().nRows(), 
        scaled.matrix().nColumns()) {
  size_t nZYX = addXY ? 3 : addZ ? 1 : 0;
  if (scaled.isPauli()) {
    auto scaledMat = dynamic_cast<const PauliSpinorMatrices<MatsU>&>(scaled.matrix());
    Matrix<MatsT>::operator=(scaled.scalar() * scaledMat.S());
    nZYX = std::max(nZYX, scaledMat.components_.size());
    components_.reserve(nZYX);
    for (auto &mat : scaledMat.components_)
      components_.emplace_back(scaled.scalar() * mat);
    while (nZYX > scaledMat.components_.size()) {
      components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
      components_.back().clear();
      nZYX--;
    }
  } else {
    Matrix<MatsT>::operator=(scaled);
    components_.reserve(nZYX);
    while (nZYX > 0) {
      components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
      components_.back().clear();
      nZYX--;
    }
  }
}
template PauliSpinorMatrices<double>::PauliSpinorMatrices(
    const ScaledMatrix<double, double>&, bool, bool );
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(
    const ScaledMatrix<double, double>&, bool, bool );
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(
    const ScaledMatrix<dcomplex, double>&, bool, bool );
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(
    const ScaledMatrix<double, dcomplex>&, bool, bool );
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(
    const ScaledMatrix<dcomplex, dcomplex>&, bool, bool );

template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator=(const Matrix<double>&);
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=(const Matrix<double>&);
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=(const Matrix<dcomplex>&);
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator=(Matrix<double>&&);
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=(Matrix<double>&&);
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=(Matrix<dcomplex>&&);

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator=( const ScaledMatrix<ScalarT, MatsU> &rhs ) {
  if (rhs.isPauli()) {
    const PauliSpinorMatrices<MatsU>& r =
        dynamic_cast<const PauliSpinorMatrices<MatsU>&>(rhs.matrix());
    S() = rhs.scalar() * r.S();
    if (r.hasZ()) {
      if (not hasZ())
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
      Z() = rhs.scalar() * r.Z();
    } else if (hasZ()) {
      Z().clear();
    }
    if (r.hasXY()) {
      if (not hasXY()) {
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
      }
      Y() = rhs.scalar() * r.Y();
      X() = rhs.scalar() * r.X();
    } else if (hasXY()) {
      Y().clear();
      X().clear();
    }
  } else {
    S() = rhs;
    for (auto &mat : components_) mat.clear();
  }
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( const ScaledMatrix<dcomplex, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( const ScaledMatrix<double, dcomplex>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator=( const PauliSpinorMatrices<MatsT> &other ) {
  if (this != &other) { // self-assignment check expected
    return operator=(1.0 * other);
  }
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator=( const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( const PauliSpinorMatrices<dcomplex>& );

template <typename MatsT>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator=( PauliSpinorMatrices<MatsT> &&other ) {
  if (this != &other) { // self-assignment check expected
    S() = std::move(other.S());
    if (other.hasZ()) {
      if (hasZ())
        Z() = std::move(other.Z());
      else
        components_.emplace_back(std::move(other.Z()));
    } else if (hasZ()) {
      Z().clear();
    }
    if (other.hasXY()) {
      if (hasXY()) {
        Y() = std::move(other.Y());
        X() = std::move(other.X());
      } else {
        components_.emplace_back(std::move(other.Y()));
        components_.emplace_back(std::move(other.X()));
      }
    } else if (hasXY()) {
      Y().clear();
      X().clear();
    }
  }
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator=( PauliSpinorMatrices<double>&& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator=( PauliSpinorMatrices<dcomplex>&& );

template <typename MatsT>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator*=( MatsT scalar ) {
  S() *= scalar;
  for (Matrix<MatsT> &mat : components_)
    mat *= scalar;
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator*=( double scalar );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator*=( dcomplex scalar );

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator+=( const PauliSpinorMatrices<MatsU> &other ) {
  S() += other.S();
  if (hasZ() and other.hasZ())
    Z() += other.Z();
  else if (other.hasZ())
    components_.emplace_back(other.Z());
  if (hasXY() and other.hasXY()) {
    Y() += other.Y();
    X() += other.X();
  } else if (other.hasXY()) {
    components_.emplace_back(other.Y());
    components_.emplace_back(other.X());
  }
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator+=( const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const PauliSpinorMatrices<dcomplex>& );

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator+=( const Matrix<MatsU> &other ) {
  S() += other;
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator+=( const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const Matrix<dcomplex>& );

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator-=( const Matrix<MatsU> &other ) {
  S() = other;
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator-=( const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const Matrix<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const Matrix<dcomplex>& );

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator+( const Matrix<MatsU> &other ) const {
  PauliSpinorMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> sum(*this);
  sum += other;
  return sum;
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator+( const Matrix<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator+( const Matrix<dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const Matrix<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const Matrix<dcomplex>& ) const;

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value),
dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator-( const Matrix<MatsU> &other ) const {
  return *this + (-other);
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator-( const Matrix<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator-( const Matrix<dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const Matrix<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const Matrix<dcomplex>& ) const;

template PauliSpinorMatrices<double>
operator+( const Matrix<double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const Matrix<dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const Matrix<double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const Matrix<dcomplex>&, const PauliSpinorMatrices<dcomplex>& );

template PauliSpinorMatrices<double>
operator-( const Matrix<double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const Matrix<dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const Matrix<double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const Matrix<dcomplex>&, const PauliSpinorMatrices<dcomplex>& );

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value), dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator+( const PauliSpinorMatrices<MatsU> &other ) const {
  typedef typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<MatsU, dcomplex>::value),
      dcomplex, double>::type ResultsT;
  if (other.components_.size() > components_.size()) {
    PauliSpinorMatrices<ResultsT> sum(other);
    sum += *this;
    return sum;
  }
  PauliSpinorMatrices<ResultsT> sum(*this);
  sum += other;
  return sum;
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator+( const PauliSpinorMatrices<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const PauliSpinorMatrices<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator+( const PauliSpinorMatrices<dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const PauliSpinorMatrices<dcomplex>& ) const;

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value), dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator-( const PauliSpinorMatrices<MatsU> &other ) const {
  return *this + (-other);
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator-( const PauliSpinorMatrices<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const PauliSpinorMatrices<double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator-( const PauliSpinorMatrices<dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const PauliSpinorMatrices<dcomplex>& ) const;

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator-=( const PauliSpinorMatrices<MatsU> &other ) {
  return operator+=(-other);
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator-=( const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const PauliSpinorMatrices<dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator+=( const ScaledMatrix<ScalarT, MatsU> &scaled ) {
  if (scaled.isPauli()) {
    const PauliSpinorMatrices<MatsU>& scaledMat =
        dynamic_cast<const PauliSpinorMatrices<MatsU>&>(scaled.matrix());
    S() += scaled.scalar() * scaledMat.S();
    if (scaledMat.hasZ()) {
      if (hasZ())
        Z() += scaled.scalar() * scaledMat.Z();
      else {
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
        Z() = scaled.scalar() * scaledMat.Z();
      }
    }
    if (scaledMat.hasXY()) {
      if (hasXY()) {
        Y() += scaled.scalar() * scaledMat.Y();
        X() += scaled.scalar() * scaledMat.X();
      } else {
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
        Y() = scaled.scalar() * scaledMat.Y();
        components_.emplace_back(this->memManager(), this->nRows(), this->nColumns());
        X() = scaled.scalar() * scaledMat.X();
      }
    }
  } else {
    S() += scaled;
  }
  return *this;
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator+=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const ScaledMatrix<dcomplex, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const ScaledMatrix<double, dcomplex>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator+=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<MatsT>&
PauliSpinorMatrices<MatsT>::operator-=( const ScaledMatrix<ScalarT, MatsU> &scaled ) {
  return operator+=(-scaled);
}
template PauliSpinorMatrices<double>&
PauliSpinorMatrices<double>::operator-=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const ScaledMatrix<double, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const ScaledMatrix<dcomplex, double>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const ScaledMatrix<double, dcomplex>& );
template PauliSpinorMatrices<dcomplex>&
PauliSpinorMatrices<dcomplex>::operator-=( const ScaledMatrix<dcomplex, dcomplex>& );

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator+( const ScaledMatrix<ScalarT, MatsU> &scaled ) const {
  PauliSpinorMatrices<typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<MatsU, dcomplex>::value or
       std::is_same<ScalarT, dcomplex>::value),
      dcomplex, double>::type> sum(*this);
  sum += scaled;
  return sum;
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator+( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator+( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator+( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator+( const ScaledMatrix<dcomplex, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator+( const ScaledMatrix<dcomplex, dcomplex>& ) const;

template <typename MatsT>
template <typename ScalarT, typename MatsU>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<MatsT, dcomplex>::value or
 std::is_same<MatsU, dcomplex>::value or
 std::is_same<ScalarT, dcomplex>::value),
dcomplex, double>::type>
PauliSpinorMatrices<MatsT>::operator-( const ScaledMatrix<ScalarT, MatsU> &scaled ) const {
  return *this + (-scaled);
}
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::operator-( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator-( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator-( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<double>::operator-( const ScaledMatrix<dcomplex, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const ScaledMatrix<double, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const ScaledMatrix<dcomplex, double>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const ScaledMatrix<double, dcomplex>& ) const;
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::operator-( const ScaledMatrix<dcomplex, dcomplex>& ) const;

template PauliSpinorMatrices<double>
operator+( const ScaledMatrix<double, double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<double, dcomplex>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator+( const ScaledMatrix<dcomplex, dcomplex>&, const PauliSpinorMatrices<dcomplex>& );

template PauliSpinorMatrices<double>
operator-( const ScaledMatrix<double, double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<double, dcomplex>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, double>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, dcomplex>&, const PauliSpinorMatrices<double>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, double>&, const PauliSpinorMatrices<dcomplex>& );
template PauliSpinorMatrices<dcomplex>
operator-( const ScaledMatrix<dcomplex, dcomplex>&, const PauliSpinorMatrices<dcomplex>& );

} // namespace cqmatrix
} // namespace ChronusQ
