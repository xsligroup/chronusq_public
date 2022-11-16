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
#include <cqlinalg.hpp>
#include <cerr.hpp>

#include <util/mpi.hpp>
#include <util/matout.hpp>

namespace ChronusQ {


  template <typename _F>
  class SolverVectors {

  public:
    // Length of each vector
    virtual size_t length() const = 0;
    // Number of vectors in container
    virtual size_t size() const = 0;
    // Get raw pointer of vectors, only works for RawVectors
    virtual _F* getPtr(size_t i = 0) = 0;
    const _F* getPtr(size_t i = 0) const {
      return const_cast<SolverVectors<_F>*>(this)->getPtr(i);
    }
    // Get element
    virtual _F get(size_t i, size_t j) const = 0;
    // Set element
    virtual void set(size_t i, size_t j, _F value) = 0;
    // Clear elements
    void clear(size_t shift = 0) {
      clear(shift, size() - shift);
    }
    virtual void clear(size_t shift, size_t nVec) = 0;
    // Print elements
    void print(std::ostream& out, std::string str, size_t shift = 0) const {
      print(out, str, shift, size() - shift);
    }
    virtual void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const = 0;
    // Get underlying type
    virtual const std::type_info& underlyingType() const {
      return typeid(*this);
    }

    /**
     * C = alpha * this * op(B) + beta * C
     * A wrapper for
     * blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,transB,
     *            length(), n, k, alpha, getPtr(), length(),
     *            B, ldb, beta, C_ptr, C.length());
     * Use to linear combine the vectors in this
     * @param shiftA The beginning vector in this
     * @param transB specifies op(B)
     * @param n      Number of vector after linear combination
     * @param k      Number of vector for linear combination in this
     * @param alpha  Scalar factor for this * op(B)
     * @param B      Linear transformation matrix
     * @param ldb    Leading dimension of B
     * @param beta   Scalar factor for C
     * @param C      Result vectors
     * @param shiftC The beginning vector in C
     */
    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 _F alpha, _F const *B, int64_t ldb,
                                 _F beta, SolverVectors<_F> &C, size_t shiftC) const = 0;

    /**
     * C = conj(this) * B
     * A wrapper for
     * blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
     *            m,n,length(),_F(1.),getPtr(),length(),
     *            B_ptr,length(),_F(0.),C,ldc);
     * @param shiftA The beginning vector in this
     * @param B      Another set of vectors
     * @param shiftB The beginning vector in B
     * @param m      Number of vectors for dot product in this
     * @param n      Number of vectors for dot product in B
     * @param C      Result matrix
     * @param ldc    Leading dimension of C
     */
    virtual void dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                             int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA = true) const = 0;

    /**
     * this[shiftA : shiftA+nVec] = B[shiftB : shiftB+nVec]
     * Copy nVec number of vectors beginning from the shiftB-th vector in B to
     * this beginning from shiftA
     * @param shiftA Beginning vector index to write in this
     * @param nVec   Number of vectors to copy
     * @param B      Source vectors
     * @param shiftB Beginning vector index to copy from in B
     */
    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable = false) = 0;

    /**
     * this[shiftA : shiftA+nVec] <-> B[shiftB : shiftB+nVec]
     * Swap nVec number of vectors beginning from the shiftB-th vector in B with
     * this beginning from shiftA
     * @param shiftA Beginning vector index to write in this
     * @param nVec   Number of vectors to swap
     * @param B      The other vectors
     * @param shiftB Beginning vector index to swap in B
     */
    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) = 0;

    /**
     * this[shiftA : shiftA+nVec] *= scalar
     * A wrapper for
     * blas::scal(length()*nVec,scalar,getPtr(shiftA),1);
     * @param shiftA Beginning vector index to scale in this
     * @param nVec   Number of vectors to scale
     * @param scalar Scalar factor
     */
    void scale(_F scalar, size_t shift = 0) {
      scale(scalar, shift, size() - shift);
    }
    virtual void scale(_F scalar, size_t shift, size_t nVec) = 0;

    /**
     * this[shiftA : shiftA+nVec] = conj(this[shiftA : shiftA+nVec])
     * @param shiftA Beginning vector index to conjugate in this
     * @param nVec   Number of vectors to conjugate
     */
    void conjugate(size_t shift = 0) {
      conjugate(shift, size() - shift);
    }
    virtual void conjugate(size_t shift, size_t nVec) = 0;

    /**
     * this[shiftY : shiftY+nVec] += alpha * X[shiftX : shiftX+nVec]
     * A wrapper for
     * blas::axpy(length() * nVec, alpha, X_ptr, 1, getPtr(shiftY), 1);
     * @param shiftY Beginning vector index to add in this
     * @param nVec   Number of vectors to add
     * @param alpha  Scalar factor for X[shiftX : shiftX+nVec]
     * @param X      Vectors to add
     * @param shiftX Beginning vector index to add in X
     */
    virtual void axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) = 0;

    /**
     * GramSchmidt orthogonalization of vectors
     * @param shift Beginning vector index
     * @param Mold  Number of vectors that are already orthonormal
     * @param Mnew  Number of vectors to be orthonormalize
     * @param mem   MemManager referece
     * @param NRe   Number of repeats of projection
     * @param eps   Threshold for linear dependency
     * @return
     */
    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                               size_t NRe = 0, double eps = 1e-12);

    /**
     * Solve trangular linear system  X * A = alpha * B
     * A wrapper for
     * blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
     *            length(), n, alpha, A, lda, getPtr(), length());
     * @param shift Beginning vector index
     * @param n     Number of vectors
     * @param alpha Scalar for rhs
     * @param A     Triangular matrix
     * @param lda   Leading dimension of A
     */
    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) = 0;

    /**
     * QR factorization
     * A wrapper for
     * ChronusQ::QR(length(), nVec, getPtr(), length(), R, LDR, mem);
     * @param shift Beginning vector index
     * @param nVec  Number of vectors
     * @param mem   Reference to CQMemManager
     * @param R     returns R
     * @param LDR   Leading dimension of R
     * @return      Lapack information
     */
    virtual int QR(size_t shift, size_t nVec, CQMemManager &mem, _F *R = nullptr, int LDR = 0) = 0;

    /**
     * 2-norm of vector this[shift] or F-norm of matrix this[shift : shift+nVec]
     * A wrapper for
     * blas::nrm2(length() * nVec,getPtr(shift),1);
     * @param shift Beginning vector index to compute norm in this
     * @param nVec  Number of vectors to compute norm
     * @return      2(F)-Norm of vector(s)
     */
    double norm2F(size_t shift = 0) const {
      return norm2F(shift, size() - shift);
    }
    virtual double norm2F(size_t shift, size_t nVec) const = 0;

    /**
     * The norm of the element with the greatest norm in vector(s)
     * For vector, this is the inf-norm
     * @param shift Beginning vector index to find
     * @param nVec  Number of vectors to find
     * @return      The norm of the element with the greatest norm
     */
    double maxNormElement(size_t shift = 0) const {
      return maxNormElement(shift, size() - shift);
    }
    virtual double maxNormElement(size_t shift, size_t nVec) const = 0;

    virtual ~SolverVectors() {}

  }; // class SolverVectors


  template <typename _F>
  class RawVectors : public SolverVectors<_F> {

  protected:

    MPI_Comm      comm_;
    CQMemManager &memManager_;
    _F* data_ = nullptr;
    size_t len_;
    size_t size_ = 0;

  public:
    RawVectors(MPI_Comm c, CQMemManager &mem, size_t len, size_t size)
    : comm_(c), memManager_(mem), len_(len), size_(size) {
      if (MPIRank(comm_) == 0 and size > 0)
        data_ = memManager_.malloc<_F>(len_ * size_);
    }
    RawVectors(const RawVectors<_F> &other)
    : comm_(other.comm_), memManager_(other.memManager_),
    len_(other.len_), size_(other.size_) {
      if (MPIRank(comm_) == 0 and size_ > 0 and other.data_ != nullptr) {
        data_ = memManager_.malloc<_F>(len_ * size_);
        std::copy_n(other.data_, len_ * size_, data_);
      }
    }
    RawVectors(RawVectors<_F> &&other)
    : comm_(other.comm_), memManager_(other.memManager_),
    data_(other.data_), len_(other.len_), size_(other.size_) {
      other.data_ = nullptr;
    }

    MPI_Comm getMPIcomm() const { return comm_; }
    CQMemManager& getMem() const { return memManager_; }

    virtual size_t length() const override { return len_; }
    virtual size_t size() const override { return size_; }

    using SolverVectors<_F>::getPtr;
    virtual _F* getPtr(size_t i = 0) override {
#ifdef CQ_ENABLE_MPI
      if (MPIRank(comm_) != 0 or size() == 0)
        return nullptr;
#endif
      if (i >= size())
        CErr("Requesting invalid pointer in RawVectors object.");
      return data_ + i * len_;
    }
    // Get element
    virtual _F get(size_t i, size_t j) const override {
#ifdef CQ_ENABLE_MPI
      _F v;
      if (MPIRank(comm_) == 0)
        v = getPtr(j)[i];
      if (MPISize(comm_) > 1)
        MPIBCast(v, 0, comm_);
      return v;
#else
      return getPtr(j)[i];
#endif
    }
    // Set element
    virtual void set(size_t i, size_t j, _F value) override {
      ROOT_ONLY(comm_);
      getPtr(j)[i] = value;
    }

    using SolverVectors<_F>::clear;
    void clear(size_t shift, size_t nVec) override {
      if (nVec == 0) return;
      ROOT_ONLY(comm_);
      getPtr(shift + nVec - 1);
      std::fill_n(getPtr(shift), nVec*length(), 0.);
    }

    using SolverVectors<_F>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      ROOT_ONLY(comm_);
      if (nVec == 0) {
        out << std::endl << str + ": " << std::endl;
        return;
      }
      getPtr(shift + nVec - 1);
      prettyPrintSmart(out, str, getPtr(shift), length(), nVec, length());
    }

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 _F alpha, _F const *B, int64_t ldb,
                                 _F beta, SolverVectors<_F> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                             int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) override;

    using SolverVectors<_F>::scale;
    virtual void scale(_F scalar, size_t shift, size_t nVec) override;

    using SolverVectors<_F>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) override;

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, CQMemManager &mem, _F *R = nullptr, int LDR = 0) override;

    using SolverVectors<_F>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<_F>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual ~RawVectors() {
      if (data_)
        memManager_.free(data_);
    }

  }; // class RawVectors


  template <typename _F>
  class SolverVectorsView : public SolverVectors<_F> {

  protected:

    SolverVectors<_F> &vecs_;
    size_t shift_;

  public:
    SolverVectorsView(SolverVectors<_F> &vecs, size_t shift = 0)
    : vecs_(typeid(vecs) == typeid(SolverVectorsView<_F>) ? dynamic_cast<SolverVectorsView<_F>&>(vecs).vecs_ : vecs),
    shift_(shift) {
      if (typeid(vecs) == typeid(SolverVectorsView<_F>))
        shift_ += dynamic_cast<SolverVectorsView<_F>&>(vecs).shift_;
      if (shift_ >= vecs_.size())
        CErr("Creating a view out of the vectors' capacity.");
    }

    size_t shift() const { return shift_; }
    virtual size_t length() const override { return vecs_.length(); }
    virtual size_t size() const override { return vecs_.size() - shift(); }

    SolverVectors<_F>& getVecs() {
      return vecs_;
    }
    const SolverVectors<_F>& getVecs() const {
      return vecs_;
    }

    using SolverVectors<_F>::getPtr;
    virtual _F* getPtr(size_t i = 0) override {
      return vecs_.getPtr(shift() + i);
    }
    // Get element
    virtual _F get(size_t i, size_t j) const override {
      return vecs_.get(i, shift() + j);
    }
    // Set element
    virtual void set(size_t i, size_t j, _F value) override {
      vecs_.set(i, shift() + j, value);
    }

    using SolverVectors<_F>::clear;
    void clear(size_t shift, size_t nVec) override {
      vecs_.clear(this->shift() + shift, nVec);
    }

    using SolverVectors<_F>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      out << "Printing a SolverVectorsView object with shift: " << this->shift() << std::endl;
      vecs_.print(out, str, shift + this->shift(), nVec);
    }

    // Get underlying type
    virtual const std::type_info& underlyingType() const override{
      return typeid(vecs_);
    }

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 _F alpha, _F const *B, int64_t ldb,
                                 _F beta, SolverVectors<_F> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                             int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) override;

    using SolverVectors<_F>::scale;
    virtual void scale(_F scalar, size_t shift, size_t nVec) override;

    using SolverVectors<_F>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) override;

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew, CQMemManager &mem,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, CQMemManager &mem, _F *R = nullptr, int LDR = 0) override;

    using SolverVectors<_F>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<_F>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual ~SolverVectorsView() {}

  }; // class SolverVectorsView

}; // namespace
