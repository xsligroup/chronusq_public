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

#include <itersolver/solvervectors.hpp>
#include <cqlinalg/ortho.hpp>

namespace ChronusQ {

  template <typename _F>
  void RawVectors<_F>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                       _F alpha, _F const *B, int64_t ldb,
                                       _F beta, SolverVectors<_F> &C, size_t shiftC) const {
    ROOT_ONLY(comm_);

    this->sizeCheck(shiftA + n, "A during RawVectors<_F>::multiply_matrix");
    C.sizeCheck(shiftC + n, "C during RawVectors<_F>::multiply_matrix");
    
    tryDowncastReferenceTo<RawVectors<_F>>(C,
        [&] (auto& CRef, size_t extraShiftC) {
          shiftC += extraShiftC;
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,transB,
              length(), n, k, alpha, getPtr(shiftA), length(),
              B, ldb, beta, CRef.getPtr(shiftC), CRef.length());
        }
    );
  }

  template <typename _F>
  void DistributedVectors<_F>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                               _F alpha, _F const *B, int64_t ldb,
                                               _F beta, SolverVectors<_F> &C, size_t shiftC) const {
    
    this->sizeCheck(shiftA + n, "A during DistributedVectors<_F>::multiply_matrix");
    C.sizeCheck(shiftC + n, "C during DistributedVectors<_F>::multiply_matrix");

    // downcasting C
    tryDowncastReferenceTo<DistributedVectors<_F>>(C,
                                                   [&](auto& CRef, size_t extraShiftC) {
          if (localLength() != CRef.localLength()) {
            CErr("C Doesn't have same local length as A in multiply_matrix"); 
          }
          shiftC += extraShiftC;
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,transB,
                     localLength(), n, k, alpha, getLocalPtr(shiftA), localLength(),
                     B, ldb, beta, CRef.getLocalPtr(shiftC), CRef.localLength());
        }
    );

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
    if (m * n == 0) return;

    if (MPIRank(comm_) == 0) {
      this->sizeCheck(shiftA + m, "A during RawVectors<_F>::dot_product");
      B.sizeCheck(shiftB + n, "B during RawVectors<_F>::dot_product");

      if (length() != B.length())
        CErr("Lengths of vectors does not match for dot product.");

      tryDowncastReferenceTo<RawVectors<_F>>(B,
          [&] (auto& BRef, size_t extraShiftB) {
            shiftB += extraShiftB;
            blas::gemm(blas::Layout::ColMajor,
                       conjA ? blas::Op::ConjTrans : blas::Op::Trans,
                       blas::Op::NoTrans,
                       m,n,length(),_F(1.),getPtr(shiftA),length(),
                       BRef.getPtr(shiftB),length(),_F(0.),C,ldc);
          }
       );
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1) {
      bool localC = C == nullptr;
      if (localC)
        C = CQMemManager::get().malloc<_F>(ldc * n);
      if (ldc == m)
        MPIBCast(C, m*n, 0, comm_);
      else
        for (size_t i = 0; i < n; i++)
          MPIBCast(C + i * ldc, m, 0, comm_);
      if (localC)
        CQMemManager::get().free(C);
    }
#endif

  }

  template <typename _F>
  void DistributedVectors<_F>::dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                                           int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA) const {
    if (m * n == 0) return;

    this->sizeCheck(shiftA + m, "A during DistributedVectors<_F>::dot_product");
    B.sizeCheck(shiftB + n, "B during DistributedVectors<_F>::dot_product");

    if (length() != B.length())
      CErr("Lengths of vectors does not match for dot product.");

#ifdef CQ_ENABLE_MPI
    _F* CRes = CQMemManager::get().malloc<_F>(m * n);
#endif

    tryDowncastReferenceTo<DistributedVectors<_F>>(B,
                                                   [&](auto& BRef, size_t extraShiftB) {
          if (localLength() != BRef.localLength()) {
            CErr("B Doesn't have same local length as A in dot product"); 
          }
          shiftB += extraShiftB;
          blas::gemm(blas::Layout::ColMajor,
                     conjA ? blas::Op::ConjTrans : blas::Op::Trans,
                     blas::Op::NoTrans,
                     m,n,localLength(),_F(1.),getLocalPtr(shiftA),localLength(),
                     BRef.getLocalPtr(shiftB),localLength(),_F(0.),
#ifdef CQ_ENABLE_MPI
                     CRes, m);                  
#else                     
                     C,ldc);
#endif
        }
    );
    
#ifdef CQ_ENABLE_MPI
    bool localC = C == nullptr;
    if (localC) C = CQMemManager::get().malloc<_F>(ldc * n);
    
    // AllReduce 
    if (m == ldc) {
      MPIAllReduce(CRes, m * n, C, comm_); 
    } else {
      for (size_t i = 0; i < n; i++) {
        MPIAllReduce(CRes + i * m, m, C + i * ldc, comm_);
      }
    }

    CQMemManager::get().free(CRes);
    if (localC) CQMemManager::get().free(C);
#endif
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                                          int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA) const {

    vecs_.dot_product(shift() + shiftA, B, shiftB, m, n, C, ldc, conjA);

  }

  template <typename _F>
  void RawVectors<_F>::set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable) {
    if (nVec == 0) return;
    ROOT_ONLY(comm_);

    if (length() != B.length())
      CErr("Lengths of vectors does not match for set_data.");
    
    this->sizeCheck(shiftA + nVec, "A during RawVectors<_F>::set_data");
    B.sizeCheck(shiftB + nVec, "B during RawVectors<_F>::set_data");

    tryDowncastReferenceTo<RawVectors<_F>>(B,
        [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          if (getPtr(shiftA) == BRef.getPtr(shiftB)) return; 
          
          std::copy_n(BRef.getPtr(shiftB), length() * nVec, getPtr(shiftA));
        }
    );
  }

  template <typename _F>
  void DistributedVectors<_F>::set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable) {
    
    // check more length
    if (length() != B.length())
      CErr("Lengths of vectors does not match for set_data.");
    
    this->sizeCheck(shiftA + nVec, "A during DistributedVectors<_F>::set_data");
    B.sizeCheck(shiftB + nVec, "B during DistributedVectors<_F>::set_data");
    
    tryDowncastReferenceTo<DistributedVectors<_F>>(B,
                                                   [&](auto& BRef, size_t extraShiftB) {
          if (BRef.localLength() != localLength() ) { 
            CErr("A and B don't have matching size locally to perform set data");
          }
          shiftB += extraShiftB;
          std::copy_n(BRef.getLocalPtr(shiftB), localLength() * nVec, getLocalPtr(shiftA));
        }
    );
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable) {

    vecs_.set_data(shift() + shiftA, nVec, B, shiftB, moveable);

  }

  template <typename _F>
  void RawVectors<_F>::swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) {
    if (nVec == 0) return;
    ROOT_ONLY(comm_);

    if (length() != B.length())
      CErr("Lengths of vectors does not match for swap_data.");

    this->sizeCheck(shiftA + nVec, "A during RawVectors<_F>::swap_data");
    B.sizeCheck(shiftB + nVec, "B during RawVectors<_F>::swap_data");

    tryDowncastReferenceTo<RawVectors<_F>>(B,
        [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          blas::swap(length() * nVec, getPtr(shiftA), 1, BRef.getPtr(shiftB), 1);
        }
    );
  }

  template <typename _F>
  void DistributedVectors<_F>::swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) {
     
    // check more length
    if (length() != B.length())
      CErr("Lengths of vectors does not match for set_data.");
    
    this->sizeCheck(shiftA + nVec, "A during DistributedVectors<_F>::swap_data");
    B.sizeCheck(shiftB + nVec, "B during DistributedVectors<_F>::swap_data");
    
    tryDowncastReferenceTo<DistributedVectors<_F>>(B,
                                                   [&](auto& BRef, size_t extraShiftB) {
          if (BRef.localLength() != localLength() ) { 
            CErr("A and B don't have matching size locally to perform swapping data");
          }
          shiftB += extraShiftB;
          blas::swap(localLength() * nVec, getLocalPtr(shiftA), 1, BRef.getLocalPtr(shiftB), 1);
        }
    );
  } 
  
  template <typename _F>
  void SolverVectorsView<_F>::swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) {

    vecs_.swap_data(shift() + shiftA, nVec, B, shiftB);

  }

  template <typename _F>
  void RawVectors<_F>::scale(_F scalar, size_t shiftA, size_t nVec) {
    if (nVec == 0) return;
    ROOT_ONLY(comm_);

    this->sizeCheck(shiftA + nVec, "A during RawVectors<_F>::scale");

    blas::scal(length()*nVec,scalar,getPtr(shiftA),1);

  }

  template <typename _F>
  void DistributedVectors<_F>::scale(_F scalar, size_t shiftA, size_t nVec) {
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "A during DistributedVectors<_F>::scale");
    blas::scal(localLength() * nVec, scalar ,getLocalPtr(shiftA), 1);
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::scale(_F scalar, size_t shiftA, size_t nVec) {

    vecs_.scale(scalar, shift() + shiftA, nVec);

  }

  template <typename _F>
  void RawVectors<_F>::conjugate(size_t shiftA, size_t nVec) {

    ROOT_ONLY(comm_);

    if (std::is_same<_F, double>::value) return;
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "A during RawVectors<_F>::conjugate");

    std::for_each(getPtr(shiftA), getPtr(shiftA) + length()*nVec, [](_F &v) {
      v = SmartConj(v);
    });

  }

  template <typename _F>
  void DistributedVectors<_F>::conjugate(size_t shiftA, size_t nVec) {
    
    if (std::is_same<_F, double>::value) return;
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "A during DistributedVectors<_F>::conjugate");

    std::for_each(getLocalPtr(shiftA), getLocalPtr(shiftA) + localLength() * nVec,
        [] (_F &v) { v = SmartConj(v);} );
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::conjugate(size_t shiftA, size_t nVec) {

    vecs_.conjugate(shift() + shiftA, nVec);

  }

  template <typename _F>
  void RawVectors<_F>::axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) {
    if (nVec == 0) return;
    ROOT_ONLY(comm_);

    if (length() != X.length())
      CErr("Lengths of vectors does not match for axpy.");

    this->sizeCheck(shiftY + nVec, "Y during RawVectors<_F>::axpy");
    X.sizeCheck(shiftX + nVec, "X during RawVectors<_F>::axpy");

    tryDowncastReferenceTo<RawVectors<_F>>(X,
        [&] (auto& XRef, size_t extraShiftX) {
          shiftX += extraShiftX;
          blas::axpy(length() * nVec, alpha, XRef.getPtr(shiftX), 1, getPtr(shiftY), 1);
        }
    );
  }

  template <typename _F>
  void DistributedVectors<_F>::axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) {
    
    if (length() != X.length())
      CErr("Lengths of vectors does not match for axpy.");
    
    this->sizeCheck(shiftY + nVec, "Y during DistributedVectors<_F>::axpy");
    X.sizeCheck(shiftX + nVec, "X during DistributedVectors<_F>::axpy");
    
    tryDowncastReferenceTo<DistributedVectors<_F>>(X,
                                                   [&](auto& XRef, size_t extraShiftX) {
          if (XRef.localLength() != localLength() ) { 
            CErr("X and Y don't have matching size locally to perform axpy");
          }
          shiftX += extraShiftX;
          blas::axpy(localLength() * nVec, alpha, XRef.getLocalPtr(shiftX), 1, getLocalPtr(shiftY), 1);
        }
    );
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) {

    vecs_.axpy(shift() + shiftY, nVec, alpha, X, shiftX);

  }

  template <typename _F>
  size_t SolverVectors<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                                     size_t NRe, double eps) {
    if (Mnew == 0) return Mold;

    _F * SCR = CQMemManager::get().malloc<_F>(Mold + Mnew);

    if( Mold == 0 ) {
      // Normalize the first vector
      double inner = norm2F(shift, 1);
      if(std::abs(inner) < std::sqrt(length())*eps) CErr("Zero inner product incurred!");
      scale(1.0/inner, shift, 1);
    }

    // Orthonormalize the rest of the matrix using GS
    size_t iOrtho = (Mold == 0) ? 1: Mold;
    for(auto k = iOrtho; k < (Mold + Mnew); k++) {

      if( k != iOrtho ) swap_data(iOrtho + shift, 1, *this, k + shift);

      // Project out the inner products
      for(auto iRe = 0; iRe < (NRe+1); iRe++) {
        dot_product(shift, *this, shift + iOrtho, iOrtho,1, SCR, iOrtho);
        multiply_matrix(shift, blas::Op::NoTrans, 1, iOrtho, _F(-1.0), SCR, iOrtho, _F(1.0), *this, shift + iOrtho);
      }

      // Normalize the new vector
      double inner = norm2F(iOrtho + shift, 1);
      std::cout << k << " " << inner << std::endl;
      if(std::abs(inner) < std::sqrt(length())*eps) {
        std::cout << "Zero inner product incurred! " << k << "\n";
        scale(0.0, iOrtho + shift, 1);
      } else {
        scale(1.0 / inner, iOrtho + shift, 1);
        iOrtho++;
      }
    }

    CQMemManager::get().free(SCR);

    return iOrtho;

  }

  template <typename _F>
  size_t RawVectors<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                                     size_t NRe, double eps) {
    size_t n = Mold;

    if (MPIRank(comm_) == 0 and Mnew > 0) {
      this->sizeCheck(shift + Mold + Mnew, "in RawVectors<_F>::GramSchimit");

      n = ChronusQ::GramSchmidt(length(), Mold, Mnew, getPtr(shift), length(), NRe, eps);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(n, 0, comm_);
#endif

    return n;

  }

  template <typename _F>
  size_t SolverVectorsView<_F>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                                            size_t NRe, double eps) {

    return vecs_.GramSchmidt(this->shift() + shift, Mold, Mnew, NRe, eps);

  }

  template <typename _F>
  void RawVectors<_F>::trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) {
    ROOT_ONLY(comm_);
    if (n > 0) this->sizeCheck(shift + n, "in RawVectors<_F>::trsm");

    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
               length(), n, alpha, A, lda, getPtr(shift), length());

  }

  template <typename _F>
  void DistributedVectors<_F>::trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) {
    CErr("trsm NYI for DistributedVectors");
  }
  
  template <typename _F>
  void SolverVectorsView<_F>::trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) {

    vecs_.trsm(this->shift() + shift, n, alpha, A, lda);

  }

  template <typename _F>
  int RawVectors<_F>::QR(size_t shift, size_t nVec, _F *R, int LDR) {
    int n = 0;

    if (MPIRank(comm_) == 0) {
      if (n > 0) this->sizeCheck(shift + nVec, "in RawVectors<_F>::QR"); 

      if (R)
        n = ChronusQ::QR(length(), nVec, getPtr(shift), length(), R, LDR);
      else
        n = ChronusQ::QR(length(), nVec, getPtr(shift), length());
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(n, 0, comm_);
#endif

    return n;

  }

  template <typename _F>
  int DistributedVectors<_F>::QR(size_t shift, size_t nVec, _F *R, int LDR) {
    CErr("QR NYI for DistributedVectors");
    return 0;
  }
  
  template <typename _F>
  int SolverVectorsView<_F>::QR(size_t shift, size_t nVec, _F *R, int LDR) {

    return vecs_.QR(this->shift() + shift, nVec, R, LDR);

  }

  template <typename _F>
  double RawVectors<_F>::norm2F(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;
    double v = 0.0;

    if (MPIRank(comm_) == 0) {

      this->sizeCheck(shift + nVec, "in RawVectors<_F>::norm2F");

      v = blas::nrm2(length() * nVec,getPtr(shift),1);
    }

#ifdef CQ_ENABLE_MPI
    if (MPISize(comm_) > 1)
      MPIBCast(v, 0, comm_);
#endif

    return v;

  }

  template <typename _F>
  double DistributedVectors<_F>::norm2F(size_t shift, size_t nVec) const {
    
    this->sizeCheck(shift + nVec, "in DistributedVectors<_F>::norm2F");
    
    double sqrtv = blas::nrm2(localLength() * nVec, getLocalPtr(shift), 1);
    double v = sqrtv * sqrtv;
    
    return std::sqrt(MPIAllReduce(v, comm_)); 
  }
  
  template <typename _F>
  double SolverVectorsView<_F>::norm2F(size_t shift, size_t nVec) const {

    return vecs_.norm2F(this->shift() + shift, nVec);

  }

  template <typename _F>
  double RawVectors<_F>::maxNormElement(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;
    double v = 0.0;

    if (MPIRank(comm_) == 0) {

      this->sizeCheck(shift + nVec, "in RawVectors<_F>::maxNormElement");

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
  double DistributedVectors<_F>::maxNormElement(size_t shift, size_t nVec) const {
    // TODO: need to sync the vectors across different nodes before calling this
    this->sizeCheck(shift + nVec, "in DistributedVectors<_F>::maxNormElement");
    double v = 0.0;
    if (nVec == 0) return v;
    
    v = std::max(v, std::abs(*std::max_element(getLocalPtr(shift), 
        getLocalPtr(shift) + nVec * localLength(),
        [] (_F A, _F B) { return std::norm(A) < std::norm(B); })));
    
    // all reduce
#ifdef CQ_ENABLE_MPI
    //return MPIAllReduce(v, [](double x, double y) { return std::max(x, y);}, comm_);
    double m;
    MPI_Allreduce(&v, &m, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return m;
#else
    return v;
#endif
  }
  
  template <typename _F>
  double SolverVectorsView<_F>::maxNormElement(size_t shift, size_t nVec) const {

    return vecs_.maxNormElement(this->shift() + shift, nVec);

  }

} // namespace ChrounsQ
