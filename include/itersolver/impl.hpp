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

#include <itersolver.hpp>
#include <itersolver/iterlinearsolver.hpp>
#include <itersolver/iterdiagonalizer.hpp>
#include <itersolver/gmres.hpp>
#include <itersolver/gplhr.hpp>
#include <itersolver/davidson.hpp>

namespace ChronusQ {

  template <typename _F>
  void RawVectors<_F>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                       _F alpha, _F const *B, int64_t ldb,
                                       _F beta, SolverVectors<_F> &C, size_t shiftC) const {
    ROOT_ONLY(comm_);

    if (k + shiftA > size() or n + shiftC > C.size())
      CErr("RawVectors does not contain enough number of vectors"
           " to match requested matrix multiplication.", std::cout);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,transB,
               length(), n, k, alpha, getPtr(shiftA), length(),
               B, ldb, beta, C.getPtr(shiftC), C.length());

  }

  template <typename _F>
  void SolverVectorsView<_F>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                              _F alpha, _F const *B, int64_t ldb,
                                              _F beta, SolverVectors<_F> &C, size_t shiftC) const {

    vecs_.multiply_matrix(shift() + shiftA, transB, n, k, alpha, B, ldb, beta, C, shiftC);

  }

  template <typename _F>
  void RawVectors<_F>::dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                                   int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA) const {

    if (MPIRank(comm_) == 0){
      if (m + shiftA > size() or n + shiftB > B.size())
        CErr("RawVectors does not contain enough number of vectors"
             " to match requested dot product.", std::cout);

      if (length() != B.length())
        CErr("Lengths of vectors does not match for dot product.");

      blas::gemm(blas::Layout::ColMajor,
                 conjA ? blas::Op::ConjTrans : blas::Op::Trans,
                 blas::Op::NoTrans,
                 m,n,length(),_F(1.),getPtr(shiftA),length(),
                 B.getPtr(shiftB),length(),_F(0.),C,ldc);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1) {
      bool localC = C == nullptr;
      if (localC)
        C = memManager_.template malloc<_F>(ldc * n);
      if (ldc == m)
        MPIBCast(C, m*n, 0, comm_);
      else
        for (size_t i = 0; i < n; i++)
          MPIBCast(C + i * ldc, m, 0, comm_);
      if (localC)
        memManager_.free(C);
    }
#endif

  }

  template <typename _F>
  void SolverVectorsView<_F>::dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                                          int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA) const {

    vecs_.dot_product(shift() + shiftA, B, shiftB, m, n, C, ldc, conjA);

  }

  template <typename _F>
  void RawVectors<_F>::set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB) {
    ROOT_ONLY(comm_);

    if (length() != B.length())
      CErr("Lengths of vectors does not match for set_data.");

    getPtr(shiftA + nVec - 1);

    std::copy_n(B.getPtr(shiftB), length() * nVec, getPtr(shiftA));

  }

  template <typename _F>
  void SolverVectorsView<_F>::set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB) {

    vecs_.set_data(shift() + shiftA, nVec, B, shiftB);

  }

  template <typename _F>
  void RawVectors<_F>::scale(_F scalar, size_t shiftA, size_t nVec) {
    ROOT_ONLY(comm_);

    getPtr(shiftA + nVec - 1);

    blas::scal(length()*nVec,scalar,getPtr(shiftA),1);

  }

  template <typename _F>
  void SolverVectorsView<_F>::scale(_F scalar, size_t shiftA, size_t nVec) {

    vecs_.scale(scalar, shift() + shiftA, nVec);

  }

  template <typename _F>
  void RawVectors<_F>::conjugate(size_t shiftA, size_t nVec) {

    ROOT_ONLY(comm_);

    if (std::is_same<_F, double>::value) return;

    getPtr(shiftA + nVec - 1);

    std::for_each(getPtr(shiftA), getPtr(shiftA) + length()*nVec, [](_F &v) {
      v = SmartConj(v);
    });

  }

  template <typename _F>
  void SolverVectorsView<_F>::conjugate(size_t shiftA, size_t nVec) {

    vecs_.conjugate(shift() + shiftA, nVec);

  }

  template <typename _F>
  void RawVectors<_F>::axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) {
    ROOT_ONLY(comm_);

    if (length() != X.length())
      CErr("Lengths of vectors does not match for axpy.");

    getPtr(shiftY + nVec - 1);

    blas::axpy(length() * nVec, alpha, X.getPtr(shiftX), 1, getPtr(shiftY), 1);

  }

  template <typename _F>
  void SolverVectorsView<_F>::axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) {

    vecs_.axpy(shift() + shiftY, nVec, alpha, X, shiftX);

  }

  template <typename _F>
  size_t SolverVectors<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                                     size_t NRe, double eps) {

    _F * SCR = mem.template malloc<_F>(Mold + Mnew);

    if( Mold == 0 ) {
      // Normalize the first vector
      double inner = norm2F(shift, 1);
      if(std::abs(inner) < eps) CErr("Zero inner product incurred!");
      scale(1.0/inner, shift, 1);
    }

    // Orthonormalize the rest of the matrix using GS
    size_t iOrtho = (Mold == 0) ? Mold + 1: Mold;
    for(auto k = iOrtho; k < (Mold + Mnew); k++) {

      if( k != iOrtho ) set_data(iOrtho + shift, 1, *this, k + shift);

      // Project out the inner products
      for(auto iRe = 0; iRe < (NRe+1); iRe++) {
        dot_product(shift, *this, shift + iOrtho, iOrtho,1, SCR, iOrtho);
        multiply_matrix(shift, blas::Op::NoTrans, 1, iOrtho, _F(-1.0), SCR, iOrtho, _F(1.0), *this, shift + iOrtho);
      }

      // Normalize the new vector
      double inner = norm2F(iOrtho + shift, 1);
      std::cout << k << " " << inner << std::endl;
      if(std::abs(inner) < length()*eps) {
        std::cout << "Zero inner product incurred! " << k << "\n";
        scale(0.0, iOrtho + shift, 1);
      } else {
        scale(1.0 / inner, iOrtho + shift, 1);
        iOrtho++;
      }
    }

    mem.free(SCR);

    return iOrtho;

  }

  template <typename _F>
  size_t RawVectors<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                                     size_t NRe, double eps) {
    size_t n = 0;

    if (MPIRank(comm_) == 0) {
      getPtr(shift + Mold + Mnew - 1);

      n = ChronusQ::GramSchmidt(length(), Mold, Mnew, getPtr(shift), length(), mem, NRe, eps);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(n, 0, comm_);
#endif

    return n;

  }

  template <typename _F>
  size_t SolverVectorsView<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                                            size_t NRe, double eps) {

    return vecs_.GramSchmidt(this->shift() + shift, Mold, Mnew, mem, NRe, eps);

  }

  template <typename _F>
  void RawVectors<_F>::trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) {
    ROOT_ONLY(comm_);

    getPtr(shift + n - 1);

    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
               length(), n, alpha, A, lda, getPtr(shift), length());

  }

  template <typename _F>
  void SolverVectorsView<_F>::trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) {

    vecs_.trsm(this->shift() + shift, n, alpha, A, lda);

  }

  template <typename _F>
  int RawVectors<_F>::QR(size_t shift, size_t nVec, CQMemManager &mem, _F *R, int LDR) {
    int n = 0;

    if (MPIRank(comm_) == 0) {

      getPtr(shift + nVec - 1);

      if (R)
        n = ChronusQ::QR(length(), nVec, getPtr(shift), length(), R, LDR, mem);
      else
        n = ChronusQ::QR(length(), nVec, getPtr(shift), length(), mem);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(n, 0, comm_);
#endif

    return n;

  }

  template <typename _F>
  int SolverVectorsView<_F>::QR(size_t shift, size_t nVec, CQMemManager &mem, _F *R, int LDR) {

    return vecs_.QR(this->shift() + shift, nVec, mem, R, LDR);

  }

  template <typename _F>
  double RawVectors<_F>::norm2F(size_t shift, size_t nVec) const {
    double v = 0.0;

    if (MPIRank(comm_) == 0) {

      getPtr(shift + nVec - 1);

      v = blas::nrm2(length() * nVec,getPtr(shift),1);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(v, 0, comm_);
#endif

    return v;

  }

  template <typename _F>
  double SolverVectorsView<_F>::norm2F(size_t shift, size_t nVec) const {

    return vecs_.norm2F(this->shift() + shift, nVec);

  }

  template <typename _F>
  double RawVectors<_F>::maxNormElement(size_t shift, size_t nVec) const {
    double v = 0.0;

    if (MPIRank(comm_) == 0) {

      getPtr(shift + nVec - 1);

      v = std::abs(*std::max_element(getPtr(shift), getPtr(shift) + nVec * length(),
                                     [&] (_F A, _F B) { return std::norm(A) < std::norm(B); }));
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(v, 0, comm_);
#endif

    return v;

  }

  template <typename _F>
  double SolverVectorsView<_F>::maxNormElement(size_t shift, size_t nVec) const {

    return vecs_.maxNormElement(this->shift() + shift, nVec);

  }

};


