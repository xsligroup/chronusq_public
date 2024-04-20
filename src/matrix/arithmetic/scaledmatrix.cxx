/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute itnd/or modify
 *  it under the terms of the GNU General Public Licenses published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option)ny later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUTNY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received copy of the GNU General Public Licenselong
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

template ScaledMatrix<double, double>::
    ScaledMatrix(double, const ScaledMatrix<double, double>&);
template ScaledMatrix<dcomplex, double>::
    ScaledMatrix(dcomplex, const ScaledMatrix<dcomplex, double>&);
template ScaledMatrix<dcomplex, double>::
    ScaledMatrix(dcomplex, const ScaledMatrix<double, double>&);
template ScaledMatrix<dcomplex, double>::
    ScaledMatrix(double, const ScaledMatrix<dcomplex, double>&);
template ScaledMatrix<double, dcomplex>::
    ScaledMatrix(double, const ScaledMatrix<double, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>::
    ScaledMatrix(dcomplex, const ScaledMatrix<dcomplex, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>::
    ScaledMatrix(dcomplex, const ScaledMatrix<double, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>::
    ScaledMatrix(double, const ScaledMatrix<dcomplex, dcomplex>&);

template ScaledMatrix<double, double>
operator*(double, const Matrix<double>&);
template ScaledMatrix<dcomplex, double>
operator*(dcomplex, const Matrix<double>&);
template ScaledMatrix<double, dcomplex>
operator*(double, const Matrix<dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>
operator*(dcomplex, const Matrix<dcomplex>&);
template ScaledMatrix<double, double>
operator*(const Matrix<double>&, double);
template ScaledMatrix<dcomplex, double>
operator*(const Matrix<double>&, dcomplex);
template ScaledMatrix<double, dcomplex>
operator*(const Matrix<dcomplex>&, double);
template ScaledMatrix<dcomplex, dcomplex>
operator*(const Matrix<dcomplex>&, dcomplex);

template ScaledMatrix<double, double>
operator*(double, const PauliSpinorMatrices<double>&);
template ScaledMatrix<dcomplex, double>
operator*(dcomplex, const PauliSpinorMatrices<double>&);
template ScaledMatrix<double, dcomplex>
operator*(double, const PauliSpinorMatrices<dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>
operator*(dcomplex, const PauliSpinorMatrices<dcomplex>&);
template ScaledMatrix<double, double>
operator*(const PauliSpinorMatrices<double>&, double);
template ScaledMatrix<dcomplex, double>
operator*(const PauliSpinorMatrices<double>&, dcomplex);
template ScaledMatrix<double, dcomplex>
operator*(const PauliSpinorMatrices<dcomplex>&, double);
template ScaledMatrix<dcomplex, dcomplex>
operator*(const PauliSpinorMatrices<dcomplex>&, dcomplex);

template ScaledMatrix<double, double>
operator*(double, const ScaledMatrix<double, double>&);
template ScaledMatrix<dcomplex, double>
operator*(dcomplex, const ScaledMatrix<double, double>&);
template ScaledMatrix<dcomplex, double>
operator*(double, const ScaledMatrix<dcomplex, double>&);
template ScaledMatrix<dcomplex, double>
operator*(dcomplex, const ScaledMatrix<dcomplex, double>&);
template ScaledMatrix<double, dcomplex>
operator*(double, const ScaledMatrix<double, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>
operator*(dcomplex, const ScaledMatrix<double, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>
operator*(double, const ScaledMatrix<dcomplex, dcomplex>&);
template ScaledMatrix<dcomplex, dcomplex>
operator*(dcomplex, const ScaledMatrix<dcomplex, dcomplex>&);
template ScaledMatrix<double, double>
operator*(const ScaledMatrix<double, double>&, double);
template ScaledMatrix<dcomplex, double>
operator*(const ScaledMatrix<double, double>&, dcomplex);
template ScaledMatrix<dcomplex, double>
operator*(const ScaledMatrix<dcomplex, double>&, double);
template ScaledMatrix<dcomplex, double>
operator*(const ScaledMatrix<dcomplex, double>&, dcomplex);
template ScaledMatrix<double, dcomplex>
operator*(const ScaledMatrix<double, dcomplex>&, double);
template ScaledMatrix<dcomplex, dcomplex>
operator*(const ScaledMatrix<double, dcomplex>&, dcomplex);
template ScaledMatrix<dcomplex, dcomplex>
operator*(const ScaledMatrix<dcomplex, dcomplex>&, double);
template ScaledMatrix<dcomplex, dcomplex>
operator*(const ScaledMatrix<dcomplex, dcomplex>&, dcomplex);

template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value or
 std::is_same<MatsT1, dcomplex>::value or
 std::is_same<MatsT2, dcomplex>::value),
dcomplex, double>::type>
operator+(const ScaledMatrix<ScalarT1, MatsT1> &lhs, const ScaledMatrix<ScalarT2, MatsT2> &rhs) {
  if (not lhs.getScalarMatrix().isSameDimension(rhs.getScalarMatrix()))
    CErr("Cannot add two Matrix of different size.");
  typedef typename std::conditional<
      (std::is_same<ScalarT1, dcomplex>::value or
       std::is_same<ScalarT2, dcomplex>::value or
       std::is_same<MatsT1, dcomplex>::value or
       std::is_same<MatsT2, dcomplex>::value),
      dcomplex, double>::type RetT;
  if (rhs.isPauli()) {
    PauliSpinorMatrices<RetT> result(rhs);
    result += lhs;
    return result;
  }
  PauliSpinorMatrices<RetT> result(lhs);
  result += rhs;
  return result;
}
template PauliSpinorMatrices<double>
operator+(const ScaledMatrix<double, double>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, double>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, double>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, double>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator+(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<dcomplex, dcomplex>&);

template <typename ScalarT1, typename ScalarT2, typename MatsT1, typename MatsT2>
PauliSpinorMatrices<typename std::conditional<
(std::is_same<ScalarT1, dcomplex>::value or
 std::is_same<ScalarT2, dcomplex>::value or
 std::is_same<MatsT1, dcomplex>::value or
 std::is_same<MatsT2, dcomplex>::value),
dcomplex, double>::type>
operator-(const ScaledMatrix<ScalarT1, MatsT1>& lhs, const ScaledMatrix<ScalarT2, MatsT2>& rhs) {
  return lhs + (-rhs);
}
template PauliSpinorMatrices<double>
operator-(const ScaledMatrix<double, double>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, double>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, double>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, double>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<double, dcomplex>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<double, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<double, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<dcomplex, double>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, double>&, const ScaledMatrix<dcomplex, dcomplex>&);
template PauliSpinorMatrices<dcomplex>
operator-(const ScaledMatrix<dcomplex, dcomplex>&, const ScaledMatrix<dcomplex, dcomplex>&);

} // namespace cqmatrix
} // namespace ChronusQ
