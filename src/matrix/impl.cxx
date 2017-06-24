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

template class Matrix<double>;
template class Matrix<dcomplex>;

template class PauliSpinorMatrices<double>;
template class PauliSpinorMatrices<dcomplex>;

template <typename MatsT>
std::ostream& operator<<(std::ostream &out, const Matrix<MatsT> &mat) {
  mat.output(out);
  return out;
}

template std::ostream& operator<<(std::ostream&, const Matrix<double>&);
template std::ostream& operator<<(std::ostream&, const Matrix<dcomplex>&);

template class ScaledMatrix<double, double>;
template class ScaledMatrix<dcomplex, double>;
template class ScaledMatrix<double, dcomplex>;
template class ScaledMatrix<dcomplex, dcomplex>;

template Matrix<double>::Matrix(const PauliSpinorMatrices<double>&);
template Matrix<dcomplex>::Matrix(const PauliSpinorMatrices<double>&);
template Matrix<dcomplex>::Matrix(const PauliSpinorMatrices<dcomplex>&);

template PauliSpinorMatrices<double>::PauliSpinorMatrices(const Matrix<double>&, bool, bool);
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(const Matrix<double>&, bool, bool);
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(const Matrix<dcomplex>&, bool, bool);

template PauliSpinorMatrices<double>::PauliSpinorMatrices(const PauliSpinorMatrices<double>&, bool, bool);
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(const PauliSpinorMatrices<double>&, bool, bool);
template PauliSpinorMatrices<dcomplex>::PauliSpinorMatrices(const PauliSpinorMatrices<dcomplex>&, bool, bool);

} // namespace cqmatrix 
} // namespace ChronusQ
