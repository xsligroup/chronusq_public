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

template <typename ScalarT, typename MatsT>
class ScaledMatrix {

  template <typename ScalarT2, typename MatsT2>
  friend class ScaledMatrix;

protected:
  bool isPauli_; // Keep track of declaration type for static binding effect
  ScalarT scalar_;
  const Matrix<MatsT> &mat_;

public:

  ScaledMatrix() = delete;
  ScaledMatrix( const ScaledMatrix& ) = default;
  ScaledMatrix( ScaledMatrix&& ) = default;
  ScaledMatrix(ScalarT scalar, const Matrix<MatsT> &ints, bool isPauli = false):
      isPauli_(isPauli), scalar_(scalar), mat_(ints) {}
  ScaledMatrix(ScalarT scalar, const PauliSpinorMatrices<MatsT> &ints):
      isPauli_(true), scalar_(scalar), mat_(ints) {}
  template <typename ScalarT1, typename ScalarT2>
  ScaledMatrix(ScalarT1 scalar,
      const ScaledMatrix<ScalarT2, MatsT> &scaled):
      isPauli_(scaled.isPauli_), scalar_(scalar * scaled.scalar_), mat_(scaled.mat_) {}

  bool isPauli() const { return isPauli_; }
  ScalarT scalar() const { return scalar_; }
  const Matrix<MatsT>& matrix() const { return mat_; }

  ScaledMatrix operator-() const {
    return ScaledMatrix(-1.0, *this);
  }

}; // class ScaledMatrix

template <typename ScalarT, typename MatsT>
ScaledMatrix<ScalarT, MatsT>
operator*(ScalarT a, const Matrix<MatsT> &ints) {
  return ScaledMatrix<ScalarT, MatsT>(a, ints);
}
template <typename ScalarT, typename MatsT>
ScaledMatrix<ScalarT, MatsT>
operator*(const Matrix<MatsT> &ints, ScalarT a) {
  return ScaledMatrix<ScalarT, MatsT>(a, ints);
}
template <typename ScalarT, typename MatsT>
ScaledMatrix<ScalarT, MatsT>
operator*(ScalarT a, const PauliSpinorMatrices<MatsT> &ints) {
  return ScaledMatrix<ScalarT, MatsT>(a, ints);
}
template <typename ScalarT, typename MatsT>
ScaledMatrix<ScalarT, MatsT>
operator*(const PauliSpinorMatrices<MatsT> &ints, ScalarT a) {
  return ScaledMatrix<ScalarT, MatsT>(a, ints);
}
template <typename ScalarT1, typename ScalarT2, typename MatsT>
ScaledMatrix<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value),
dcomplex, double>::type, MatsT>
operator*(ScalarT1 a, const ScaledMatrix<ScalarT2, MatsT> &ints) {
  return ScaledMatrix<typename std::conditional<
      (std::is_same<ScalarT1, dcomplex>::value or
       std::is_same<ScalarT2, dcomplex>::value),
      dcomplex, double>::type, MatsT>(a, ints);
}
template <typename ScalarT1, typename ScalarT2, typename MatsT>
ScaledMatrix<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value),
dcomplex, double>::type, MatsT>
operator*(const ScaledMatrix<ScalarT1, MatsT> &ints, ScalarT2 a) {
  return ScaledMatrix<typename std::conditional<
      (std::is_same<ScalarT1, dcomplex>::value or
       std::is_same<ScalarT2, dcomplex>::value),
      dcomplex, double>::type, MatsT>(a, ints);
}

template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value or
 std::is_same<MatsT1, dcomplex>::value or
 std::is_same<MatsT2, dcomplex>::value),
dcomplex, double>::type>
operator+(const ScaledMatrix<ScalarT1, MatsT1>&, const ScaledMatrix<ScalarT2, MatsT2>&);

template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value or
 std::is_same<MatsT1, dcomplex>::value or
 std::is_same<MatsT2, dcomplex>::value),
dcomplex, double>::type>
operator-(const ScaledMatrix<ScalarT1, MatsT1>&, const ScaledMatrix<ScalarT2, MatsT2>&);

} // namespace cqmatrix
} // namespace ChronusQ
