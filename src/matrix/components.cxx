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
void Matrix<MatsT>::componentScatter(Matrix<MatsU> & LL,
    Matrix<MatsU> & LS, Matrix<MatsU> & SL, Matrix<MatsU> & SS, 
    bool increment) const {
    
  // dimension check
  size_t nRow = this->nRows();
  size_t nCol = this->nColumns();
  if (nRow % 2 == 1 or nCol % 2 == 1) {
    CErr("componentScatter is only supported for SqaureMatrix with even size"); 
  }
  nRow /= 2; 
  nCol /= 2;

  if( (LL.pointer() and (nRow != LL.nRows() or nCol != LL.nColumns())) or
      (LS.pointer() and (nRow != LS.nRows() or nCol != LS.nColumns())) or 
      (SL.pointer() and (nRow != SL.nRows() or nCol != SL.nColumns())) or
      (SS.pointer() and (nRow != SS.nRows() or nCol != SS.nColumns())) )
    CErr("Mismatched dimensioins in componentScatter!");
  
  ComponentScatter(nRow, nCol, this->pointer(), this->nRows(),
    LL.pointer(), nRow, LS.pointer(), nRow, SL.pointer(), nRow, SS.pointer(), nRow, increment); 
}

template <typename MatsT>
template <typename MatsU>
void Matrix<MatsT>::componentGather(const Matrix<MatsU> & LL,
    const Matrix<MatsU> & LS, const Matrix<MatsU> & SL, const Matrix<MatsU> & SS, 
    bool increment) {
    
  // dimension check
  size_t nRow = this->nRows();
  size_t nCol = this->nColumns();
  if (nRow % 2 == 1 or nCol % 2 == 1) {
    CErr("componentGather is only supported for SqaureMatrix with even size"); 
  }
  nRow /= 2; 
  nCol /= 2;

  if( (LL.pointer() and (nRow != LL.nRows() or nCol != LL.nColumns())) or
      (LS.pointer() and (nRow != LS.nRows() or nCol != LS.nColumns())) or 
      (SL.pointer() and (nRow != SL.nRows() or nCol != SL.nColumns())) or
      (SS.pointer() and (nRow != SS.nRows() or nCol != SS.nColumns())))
    CErr("Mismatched dimensioins in componentGather!");

  ComponentGather(nRow, nCol, this->pointer(), this->nRows(),
      LL.pointer(), nRow, LS.pointer(), nRow, SL.pointer(), nRow, SS.pointer(), nRow, increment); 
}

template <typename MatsT>
template <typename MatsU>
Matrix<MatsT> Matrix<MatsT>::componentGatherBuild(const Matrix<MatsU> & LL,
    const Matrix<MatsU> & LS, const Matrix<MatsU> & SL, const Matrix<MatsU> & SS) { 
  
  // get dimension
  size_t nRow, nCol;
  if (LL.pointer()) {
    nRow = LL.nRows();
    nCol = LL.nColumns();
  } else if (LS.pointer()) {
    nRow = LS.nRows();
    nCol = LS.nColumns();
  } else if (SL.pointer()) {
    nRow = SL.nRows();
    nCol = SL.nColumns();
  } else if (SS.pointer()) {
    nRow = SS.nRows();
    nCol = SS.nColumns();
  } else {
    CErr("Nothing to build in componentGatherBuild");
  }
  Matrix<MatsT> mat(nRow * 2, nCol * 2);
  mat.componentGather(LL, LS, SL, SS, false);
  
  return mat;
}    

template void Matrix<double>::componentScatter(Matrix<double> & LL,
    Matrix<double> & LS, Matrix<double> & SL, Matrix<double> & SS, bool increment) const;
template void Matrix<double>::componentScatter(Matrix<dcomplex> & LL,
    Matrix<dcomplex> & LS, Matrix<dcomplex> & SL, Matrix<dcomplex> & SS, bool increment) const;
template void Matrix<dcomplex>::componentScatter(Matrix<dcomplex> & LL,
    Matrix<dcomplex> & LS, Matrix<dcomplex> & SL, Matrix<dcomplex> & SS, bool increment) const;

template void Matrix<double>::componentGather(const Matrix<double> & LL,
    const Matrix<double> & LS, const Matrix<double> & SL, const Matrix<double> & SS, bool increment);
template void Matrix<dcomplex>::componentGather(const Matrix<double> & LL,
    const Matrix<double> & LS, const Matrix<double> & SL, const Matrix<double> & SS, bool increment);
template void Matrix<dcomplex>::componentGather(const Matrix<dcomplex> & LL,
    const Matrix<dcomplex> & LS, const Matrix<dcomplex> & SL, const Matrix<dcomplex> & SS, bool increment);

template Matrix<double> Matrix<double>::componentGatherBuild(const Matrix<double> & LL,
    const Matrix<double> & LS, const Matrix<double> & SL, const Matrix<double> & SS);
template Matrix<dcomplex> Matrix<dcomplex>::componentGatherBuild(const Matrix<double> & LL,
    const Matrix<double> & LS, const Matrix<double> & SL, const Matrix<double> & SS);
template Matrix<dcomplex> Matrix<dcomplex>::componentGatherBuild(const Matrix<dcomplex> & LL,
    const Matrix<dcomplex> & LS, const Matrix<dcomplex> & SL, const Matrix<dcomplex> & SS);

template <typename MatsT>
template <typename MatsU>
void PauliSpinorMatrices<MatsT>::componentScatter(
    PauliSpinorMatrices<MatsU> & LL, 
    PauliSpinorMatrices<MatsU> & LS, 
    PauliSpinorMatrices<MatsU> & SL, 
    PauliSpinorMatrices<MatsU> & SS, 
    bool increment) const {
    
  this->S().componentScatter(LL.S(), LS.S(), SL.S(), SS.S(), increment); 
  
  Matrix<MatsU> dummy(0);
  if (this->hasZ()) {
    Matrix<MatsU> & LLZ = LL.hasZ() ? LL.Z(): dummy;
    Matrix<MatsU> & LSZ = LS.hasZ() ? LS.Z(): dummy;
    Matrix<MatsU> & SLZ = SL.hasZ() ? SL.Z(): dummy;
    Matrix<MatsU> & SSZ = SS.hasZ() ? SS.Z(): dummy;
    this->Z().componentScatter(LLZ, LSZ, SLZ, SSZ, increment);
  }

  if (this->hasXY()) {
    Matrix<MatsU> & LLX = LL.hasXY() ? LL.X(): dummy;
    Matrix<MatsU> & LLY = LL.hasXY() ? LL.Y(): dummy;
    Matrix<MatsU> & LSX = LS.hasXY() ? LS.X(): dummy;
    Matrix<MatsU> & LSY = LS.hasXY() ? LS.Y(): dummy;
    Matrix<MatsU> & SLX = SL.hasXY() ? SL.X(): dummy;
    Matrix<MatsU> & SLY = SL.hasXY() ? SL.Y(): dummy;
    Matrix<MatsU> & SSX = SS.hasXY() ? SS.X(): dummy;
    Matrix<MatsU> & SSY = SS.hasXY() ? SS.Y(): dummy;
    this->X().componentScatter(LLX, LSX, SLX, SSX, increment);
    this->Y().componentScatter(LLY, LSY, SLY, SSY, increment);
  }
}

template <typename MatsT>
template <typename MatsU>
void PauliSpinorMatrices<MatsT>::componentGather(
    const PauliSpinorMatrices<MatsU> & LL, 
    const PauliSpinorMatrices<MatsU> & LS, 
    const PauliSpinorMatrices<MatsU> & SL, 
    const PauliSpinorMatrices<MatsU> & SS, 
    bool increment) {
  
  this->S().componentGather(LL.S(), LS.S(), SL.S(), SS.S(), increment); 
  
  Matrix<MatsU> dummy(0);
  if (this->hasZ()) {
    const Matrix<MatsU> & LLZ = LL.hasZ() ? LL.Z(): dummy;
    const Matrix<MatsU> & LSZ = LS.hasZ() ? LS.Z(): dummy;
    const Matrix<MatsU> & SLZ = SL.hasZ() ? SL.Z(): dummy;
    const Matrix<MatsU> & SSZ = SS.hasZ() ? SS.Z(): dummy;
    this->Z().componentGather(LLZ, LSZ, SLZ, SSZ, increment);
  }
  
  if (this->hasXY()) {
    const Matrix<MatsU> & LLX = LL.hasXY() ? LL.X(): dummy;
    const Matrix<MatsU> & LLY = LL.hasXY() ? LL.Y(): dummy;
    const Matrix<MatsU> & LSX = LS.hasXY() ? LS.X(): dummy;
    const Matrix<MatsU> & LSY = LS.hasXY() ? LS.Y(): dummy;
    const Matrix<MatsU> & SLX = SL.hasXY() ? SL.X(): dummy;
    const Matrix<MatsU> & SLY = SL.hasXY() ? SL.Y(): dummy;
    const Matrix<MatsU> & SSX = SS.hasXY() ? SS.X(): dummy;
    const Matrix<MatsU> & SSY = SS.hasXY() ? SS.Y(): dummy;
    this->X().componentGather(LLX, LSX, SLX, SSX, increment);
    this->Y().componentGather(LLY, LSY, SLY, SSY, increment);
  }
}

template <typename MatsT>
template <typename MatsU>
PauliSpinorMatrices<MatsT> 
PauliSpinorMatrices<MatsT>::componentGatherBuild(
    const PauliSpinorMatrices<MatsU> & LL, 
    const PauliSpinorMatrices<MatsU> & LS, 
    const PauliSpinorMatrices<MatsU> & SL, 
    const PauliSpinorMatrices<MatsU> & SS) {
    
  // get dimension
  size_t nRow, nCol;
  if (LL.S().pointer()) {
    nRow = LL.nRows();
    nCol = LL.nColumns();
  } else if (LS.S().pointer()) {
    nRow = LS.nRows();
    nCol = LS.nColumns();
  } else if (SL.S().pointer()) {
    nRow = SL.nRows();
    nCol = SL.nColumns();
  } else if (SS.S().pointer()) {
    nRow = SS.nRows();
    nCol = SS.nColumns();
  } else {
    CErr("Nothing to build in componentGatherBuild");
  }

  bool any_hasZ  = LL.hasZ()  or LS.hasZ()  or SL.hasZ()  or SS.hasZ();
  bool any_hasXY = LL.hasXY() or LS.hasXY() or SL.hasXY() or SS.hasXY();
  
  PauliSpinorMatrices<MatsT> pauli(2 * nRow, 2 * nCol, any_hasXY, any_hasZ);

  pauli.S() = Matrix<MatsT>::componentGatherBuild(LL.S(), LS.S(), SL.S(), SS.S());
  
  if(any_hasZ)
    pauli.Z() = Matrix<MatsT>::componentGatherBuild(LL.Z(), LS.Z(), SL.Z(), SS.Z());
  
  if(any_hasXY) {
    pauli.X() = Matrix<MatsT>::componentGatherBuild(LL.X(), LS.X(), SL.X(), SS.X());
    pauli.Y() = Matrix<MatsT>::componentGatherBuild(LL.Y(), LS.Y(), SL.Y(), SS.Y());
  }

  return pauli;
}


template void PauliSpinorMatrices<double>::componentScatter(
    PauliSpinorMatrices<double> & LL, PauliSpinorMatrices<double> & LS, 
    PauliSpinorMatrices<double> & SL, PauliSpinorMatrices<double> & SS, bool increment) const;
template void PauliSpinorMatrices<double>::componentScatter(
    PauliSpinorMatrices<dcomplex> & LL, PauliSpinorMatrices<dcomplex> & LS, 
    PauliSpinorMatrices<dcomplex> & SL, PauliSpinorMatrices<dcomplex> & SS, bool increment) const;
template void PauliSpinorMatrices<dcomplex>::componentScatter(
    PauliSpinorMatrices<dcomplex> & LL, PauliSpinorMatrices<dcomplex> & LS, 
    PauliSpinorMatrices<dcomplex> & SL, PauliSpinorMatrices<dcomplex> & SS, bool increment) const;

template void PauliSpinorMatrices<double>::componentGather(
    const PauliSpinorMatrices<double> & LL, const PauliSpinorMatrices<double> & LS, 
    const PauliSpinorMatrices<double> & SL, const PauliSpinorMatrices<double> & SS, bool increment);
template void PauliSpinorMatrices<dcomplex>::componentGather(
    const PauliSpinorMatrices<double> & LL, const PauliSpinorMatrices<double> & LS, 
    const PauliSpinorMatrices<double> & SL, const PauliSpinorMatrices<double> & SS, bool increment);
template void PauliSpinorMatrices<dcomplex>::componentGather(
    const PauliSpinorMatrices<dcomplex> & LL, const PauliSpinorMatrices<dcomplex> & LS, 
    const PauliSpinorMatrices<dcomplex> & SL, const PauliSpinorMatrices<dcomplex> & SS, bool increment);

template PauliSpinorMatrices<double> PauliSpinorMatrices<double>::componentGatherBuild(
    const PauliSpinorMatrices<double> & LL, const PauliSpinorMatrices<double> & LS, 
    const PauliSpinorMatrices<double> & SL, const PauliSpinorMatrices<double> & SS);
template PauliSpinorMatrices<dcomplex> PauliSpinorMatrices<dcomplex>::componentGatherBuild(
    const PauliSpinorMatrices<double> & LL, const PauliSpinorMatrices<double> & LS, 
    const PauliSpinorMatrices<double> & SL, const PauliSpinorMatrices<double> & SS);
template PauliSpinorMatrices<dcomplex> PauliSpinorMatrices<dcomplex>::componentGatherBuild(
    const PauliSpinorMatrices<dcomplex> & LL, const PauliSpinorMatrices<dcomplex> & LS, 
    const PauliSpinorMatrices<dcomplex> & SL, const PauliSpinorMatrices<dcomplex> & SS);
  
template <typename MatsT>
template <typename MatsU>
void PauliSpinorMatrices<MatsT>::componentAdd(
    const char TRANS, MatsU scale, const std::string & comp,
    const PauliSpinorMatrices<MatsU> & pauli) {
   
  if (not this->isSquareMatrix() or not pauli.isSquareMatrix()) {
    CErr("componentAdd only supported for square matrix");
  }

  if (this->nRows() != pauli.nRows()*2) 
    CErr("Dimension mismatch in componentAdd");
  
  size_t N = pauli.nRows();
  size_t twoN = 2*N;
  size_t offset;
  if (not comp.compare("LL")) offset = 0;
  else if (not comp.compare("LS")) offset = N*twoN;
  else if (not comp.compare("SL")) offset = N;
  else if (not comp.compare("SS")) offset = N + N*twoN;
  else CErr("Unsupported component types in componentAdd"); 
  
  MatAdd(TRANS,'N',N,N,scale,pauli.S().pointer(),N,
    MatsT(1.),this->S().pointer()+offset,twoN,this->S().pointer()+offset,twoN); 
  
  if (this->hasZ() and pauli.hasZ()) {
    MatAdd(TRANS,'N',N,N,scale,pauli.Z().pointer(),N,
      MatsT(1.),this->Z().pointer()+offset,twoN,this->Z().pointer()+offset,twoN); 
  }
  
  if (this->hasXY() and pauli.hasXY()) {
    MatAdd(TRANS,'N',N,N,scale,pauli.X().pointer(),N,
      MatsT(1.),this->X().pointer()+offset,twoN,this->X().pointer()+offset,twoN); 
    MatAdd(TRANS,'N',N,N,scale,pauli.Y().pointer(),N,
      MatsT(1.),this->Y().pointer()+offset,twoN,this->Y().pointer()+offset,twoN); 
  }
}

template void PauliSpinorMatrices<double>::componentAdd(const char TRANS, 
    double scale, const std::string & comp, const PauliSpinorMatrices<double> & pauli);
template void PauliSpinorMatrices<dcomplex>::componentAdd(const char TRANS,
    double scale, const std::string & comp, const PauliSpinorMatrices<double> & pauli);
template void PauliSpinorMatrices<dcomplex>::componentAdd(const char TRANS,
    dcomplex scale, const std::string & comp, const PauliSpinorMatrices<dcomplex> & pauli);

template <typename MatsT>
void PauliSpinorMatrices<MatsT>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS) {
    
  if (not this->isSquareMatrix()) {
    CErr("symmetrizeLSSL only supported for square matrix");
  }
  
  size_t twoN = this->nRows();
  if (twoN % 2 == 1) {
   CErr("symmetrizeLSSL is only supported for SqaureMatrix with even size"); 
  }
  size_t N = twoN / 2;
  size_t P1, P2;

  if (get_SL_from_LS) {
    P1 = N*twoN; // LS
    P2 = N;      // SL
  } else { 
    P1 = N;      // SL
    P2 = N*twoN; // LS
  }
  
  SetMat(TRANS,N,N,MatsT(1.),this->S().pointer()+P1,twoN,this->S().pointer()+P2,twoN); 
  if(hasZ()) SetMat(TRANS,N,N,MatsT(1.),this->Z().pointer()+P1,twoN,this->Z().pointer()+P2,twoN); 
  if(hasXY()) {
    SetMat(TRANS,N,N,MatsT(1.),this->Y().pointer()+P1,twoN,this->Y().pointer()+P2,twoN);
    SetMat(TRANS,N,N,MatsT(1.),this->X().pointer()+P1,twoN,this->X().pointer()+P2,twoN);
  }
}

template void PauliSpinorMatrices<double>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS);
template void PauliSpinorMatrices<dcomplex>::symmetrizeLSSL(char TRANS, bool get_SL_from_LS);

} // namespace cqmatrix 
} // namespace ChronusQ
