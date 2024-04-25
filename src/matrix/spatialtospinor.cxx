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
template <typename MatsU>
void Matrix<MatsT>::spinScatter(PauliSpinorMatrices<MatsU>& pauli, bool hasXY, bool hasZ) const {
  size_t nRow = nRows() / 2;
  size_t nCol = nColumns() / 2;
  assert(pauli.nRows() == nRow and pauli.nColumns() == nCol);

  MatsU *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
  if (hasZ) Z = pauli.Z().pointer();
  if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }

  SpinScatter(nRow, nCol, pointer(), nRows(), S, nRow, Z, nRow, Y, nRow, X, nRow);
}

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsU>
Matrix<MatsT>::spinScatter(bool hasXY, bool hasZ) const {
  size_t nRow = nRows() / 2;
  size_t nCol = nColumns() / 2;
  PauliSpinorMatrices<MatsU> pauli(nRow, nCol, hasXY, hasZ);
  spinScatter(pauli, hasXY, hasZ);
  return pauli;
}

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>
PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(
    const Matrix<MatsU> &AA, bool hasXY, bool hasZ) {
  size_t nRow = AA.nRows();
  size_t nCol = AA.nColumns();
  PauliSpinorMatrices<MatsT> pauli(nRow, nCol, hasXY, hasZ);
  MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
  if (hasZ) Z = pauli.Z().pointer();
  if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
  SpinScatter(nRow, nCol, AA.pointer(), nRow, reinterpret_cast<MatsU*>(NULL), nRow,
      reinterpret_cast<MatsU*>(NULL), nRow, reinterpret_cast<MatsU*>(NULL), nRow,
      S, nRow, Z, nRow, Y, nRow, X, nRow, true, true);
  return pauli;
}

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>
PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(
    const Matrix<MatsU> &AA, const Matrix<MatsU> &BB,
    bool hasXY, bool hasZ) {
  size_t nRow = AA.nRows();
  size_t nCol = AA.nColumns();
  PauliSpinorMatrices<MatsT> pauli(nRow, nCol, hasXY, hasZ);
  MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
  if (hasZ) Z = pauli.Z().pointer();
  if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
  SpinScatter(nRow, nCol, AA.pointer(), nRow, reinterpret_cast<MatsU*>(NULL), nRow,
      reinterpret_cast<MatsU*>(NULL), nRow, BB.pointer(), nRow,
      S, nRow, Z, nRow, Y, nRow, X, nRow, true, false);
  return pauli;
}

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT>
PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(
    const Matrix<MatsU> &AA, const Matrix<MatsU> &AB,
    const Matrix<MatsU> &BA, const Matrix<MatsU> &BB,
    bool hasXY, bool hasZ) {
  size_t nRow = AA.nRows();
  size_t nCol = AA.nColumns();
  PauliSpinorMatrices<MatsT> pauli(nRow, nCol, hasXY, hasZ);
  MatsT *S = pauli.S().pointer(), *Z = nullptr, *Y = nullptr, *X = nullptr;
  if (hasZ) Z = pauli.Z().pointer();
  if (hasXY) { Y = pauli.Y().pointer(); X = pauli.X().pointer(); }
  SpinScatter(nRow, nCol, AA.pointer(), nRow, AB.pointer(), nRow,
                    BA.pointer(), nRow, BB.pointer(), nRow,
                    S, nRow, Z, nRow, Y, nRow, X, nRow, false, false);
  return pauli;
}


template <typename MatsT>
template <typename MatsU>
void PauliSpinorMatrices<MatsT>::spinGather(Matrix<MatsU>& mat) const {
  size_t nRow = this->nRows();
  size_t nCol = this->nColumns();
  assert(nRow * 2 == mat.nRows() and nCol * 2 == mat.nColumns());

  const MatsT *AS = S().pointer(), *AZ = nullptr, *AY = nullptr, *AX = nullptr;
  if (hasZ()) AZ = Z().pointer();
  if (hasXY()) { AY = Y().pointer(); AX = X().pointer(); }
  
  SpinGather(nRow, nCol, mat.pointer(), 2 * nRow, AS, nRow, AZ, nRow, AY, nRow, AX, nRow, not hasXY(), not hasZ());
}

template <typename MatsT>
template <typename MatsU>
Matrix<MatsU> PauliSpinorMatrices<MatsT>::spinGather() const {
  size_t nRow = this->nRows();
  size_t nCol = this->nColumns();
  Matrix<MatsU> mat(2 * nRow, 2 * nCol);
  spinGather(mat);
  return mat;
}

template <typename MatsT>
template <typename MatsU>
std::vector<Matrix<MatsU>>
PauliSpinorMatrices<MatsT>::spinGatherToBlocks(
    bool genABBA, bool genBB) const {
  size_t nRow = this->nRows();
  size_t nCol = this->nColumns();
  std::vector<Matrix<MatsU>> blocks;
  blocks.reserve(1 + (genABBA ? 2 : 0) + (genBB ? 1 : 0));
  blocks.emplace_back(nRow, nCol);
  MatsU *AA = nullptr, *AB = nullptr, *BA = nullptr, *BB = nullptr;
  if (genABBA) {
    blocks.emplace_back(nRow, nCol);
    blocks.emplace_back(nRow, nCol);
  }
  if (genBB) { blocks.emplace_back(nRow, nCol); BB = blocks.back().pointer(); }
  AA = blocks[0].pointer();
  if (genABBA) { AB = blocks[1].pointer(); BA = blocks[2].pointer(); }

  const MatsT *AS = S().pointer(), *AZ = nullptr, *AY = nullptr, *AX = nullptr;
  if (hasZ()) AZ = Z().pointer();
  if (hasXY()) { AY = Y().pointer(); AX = X().pointer(); }

  SpinGather(nRow, nCol, AA, nRow, AB, nRow, BA, nRow, BB, nRow,
      AS, nRow, AZ, nRow, AY, nRow, AX, nRow, not hasXY(), not hasZ());
  return blocks;
}

template <typename MatsT>
template <typename MatsU>
Matrix<MatsU> Matrix<MatsT>::spatialToSpinBlock() const {
  Matrix<MatsU> spinor(2 * nRow_, 2 * nCol_);
/*
    for ( auto sp = 0ul; sp < 2; sp++)
    for ( auto nu = 0ul; nu < N_; nu++)
    for ( auto mu = 0ul; mu < N_; mu++) {
      spinor(sp*N_ + mu, sp*N_ + nu) = (*this)(mu, nu);
    }
*/
  SetMatDiag(nRow_, nCol_, pointer(), nRow_, spinor.pointer(), 2 * nRow_);
  return spinor;
}

template void Matrix<double>::spinScatter(PauliSpinorMatrices<double>&,bool,bool) const;
template void Matrix<double>::spinScatter(PauliSpinorMatrices<dcomplex>&,bool,bool) const;
template void Matrix<dcomplex>::spinScatter(PauliSpinorMatrices<dcomplex>&,bool,bool) const;
template PauliSpinorMatrices<double> Matrix<double>::spinScatter(bool,bool) const;
template PauliSpinorMatrices<dcomplex> Matrix<double>::spinScatter(bool,bool) const;
template PauliSpinorMatrices<dcomplex> Matrix<dcomplex>::spinScatter(bool,bool) const;

template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::spinBlockScatterBuild(const Matrix<double> &AA, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(const Matrix<double> &AA, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(const Matrix<dcomplex> &AA, bool, bool);
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::spinBlockScatterBuild(
    const Matrix<double> &AA, const Matrix<double> &BB, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(
    const Matrix<double> &AA, const Matrix<double> &BB, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(
    const Matrix<dcomplex> &AA, const Matrix<dcomplex> &BB, bool, bool);
template PauliSpinorMatrices<double>
PauliSpinorMatrices<double>::spinBlockScatterBuild(
    const Matrix<double> &AA, const Matrix<double> &AB,
    const Matrix<double> &BA, const Matrix<double> &BB, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(
    const Matrix<double> &AA, const Matrix<double> &AB,
    const Matrix<double> &BA, const Matrix<double> &BB, bool, bool);
template PauliSpinorMatrices<dcomplex>
PauliSpinorMatrices<dcomplex>::spinBlockScatterBuild(
    const Matrix<dcomplex> &AA, const Matrix<dcomplex> &AB,
    const Matrix<dcomplex> &BA, const Matrix<dcomplex> &BB, bool, bool);

template void PauliSpinorMatrices<double>::spinGather(Matrix<double>&) const;
template void PauliSpinorMatrices<double>::spinGather(Matrix<dcomplex>&) const;
template void PauliSpinorMatrices<dcomplex>::spinGather(Matrix<dcomplex>&) const;
template Matrix<double> PauliSpinorMatrices<double>::spinGather() const;
template Matrix<dcomplex> PauliSpinorMatrices<double>::spinGather() const;
template Matrix<dcomplex> PauliSpinorMatrices<dcomplex>::spinGather() const;

template std::vector<Matrix<double>>
PauliSpinorMatrices<double>::spinGatherToBlocks(bool genABBA, bool genBB) const;
template std::vector<Matrix<dcomplex>>
PauliSpinorMatrices<double>::spinGatherToBlocks(bool genABBA, bool genBB) const;
template std::vector<Matrix<dcomplex>>
PauliSpinorMatrices<dcomplex>::spinGatherToBlocks(bool genABBA, bool genBB) const;

template Matrix<double> Matrix<double>::spatialToSpinBlock() const;
template Matrix<dcomplex> Matrix<double>::spatialToSpinBlock() const;
template Matrix<dcomplex> Matrix<dcomplex>::spatialToSpinBlock() const;

} // namespace cqmatrix
} // namespace ChronusQ
