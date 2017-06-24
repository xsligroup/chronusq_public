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
    
    void sizeCheck(size_t i, const std::string& str = "SolverVectors") const {
      if (i > size()) CErr("Doesn't have enough vectors in" + str);
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

    virtual bool isView() const { return false; }

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
     * @param NRe   Number of repeats of projection
     * @param eps   Threshold for linear dependency
     * @return
     */
    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
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
     * ChronusQ::QR(length(), nVec, getPtr(), length(), R, LDR);
     * @param shift Beginning vector index
     * @param nVec  Number of vectors
     * @param R     returns R
     * @param LDR   Leading dimension of R
     * @return      Lapack information
     */
    virtual int QR(size_t shift, size_t nVec, _F *R = nullptr, int LDR = 0) = 0;

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
    _F* data_ = nullptr;
    size_t len_;
    size_t size_ = 0;

  public:
    RawVectors(MPI_Comm c, size_t len, size_t size)
    : comm_(c), len_(len), size_(size) {
      if (MPIRank(comm_) == 0 and size > 0) {
        data_ = CQMemManager::get().malloc<_F>(len_ * size_);
        clear();
      }
    }
    RawVectors(const RawVectors<_F> &other)
    : comm_(other.comm_),
    len_(other.len_), size_(other.size_) {
      if (MPIRank(comm_) == 0 and size_ > 0 and other.data_ != nullptr) {
        data_ = CQMemManager::get().malloc<_F>(len_ * size_);
        std::copy_n(other.data_, len_ * size_, data_);
      }
    }
    RawVectors(RawVectors<_F> &&other)
    : comm_(other.comm_),
    data_(other.data_), len_(other.len_), size_(other.size_) {
      other.data_ = nullptr;
    }

    MPI_Comm getMPIcomm() const { return comm_; }

    virtual size_t length() const override { return len_; }
    virtual size_t size() const override { return size_; }

    _F* getPtr(size_t i = 0) {
#ifdef CQ_ENABLE_MPI
      if (MPIRank(comm_) != 0 or size() == 0)
        return nullptr;
#endif
      this->sizeCheck(i, "RawVectors<_F>::getPtr");
      return data_ + i * len_;
    }
    
    const _F* getPtr(size_t i = 0) const {
      return const_cast<RawVectors*>(this)->getPtr(i);
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
      this->sizeCheck(shift + nVec, "RawVectors<_F>::clear");
      std::fill_n(getPtr(shift), nVec*length(), 0.);
    }

    using SolverVectors<_F>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      ROOT_ONLY(comm_);
      if (nVec == 0) {
        out << std::endl << str + ": " << std::endl;
        return;
      }
      this->sizeCheck(shift + nVec, "RawVectors<_F>::print");
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

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, _F *R = nullptr, int LDR = 0) override;

    using SolverVectors<_F>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<_F>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual ~RawVectors() {
      if (data_)
        CQMemManager::get().free(data_);
    }

  }; // class RawVectors

   /*
    * \brief DistributedVectors
    * 
    * The actual storage of the data is divided into blocks
    *   across different nodes.
    *
    *   Default is split evenly, but can be initialized with input
    *
    */
  template <typename _F>
  class DistributedVectors : public SolverVectors<_F> {

  protected:

    MPI_Comm comm_;
    size_t len_;
    size_t size_ = 0;
    
    // lengths across all nodes 
    std::vector<size_t> lens_; // lengths at each node
    std::vector<size_t> accLens_; // accumulated lengths at each node
    
    // local data storage 
    _F* data_ = nullptr;
    size_t localLen_;
    size_t localOffset_;

  public:
    
    // dividing the vector evenly 
    explicit DistributedVectors(MPI_Comm c,
                                size_t len, size_t size) :
        comm_(c), len_(len), size_(size) {
      
      size_t nNodes = MPISize(comm_);
      
      size_t blockLen = std::ceil( double(len_) / nNodes);  
      lens_ = {blockLen};
      accLens_ = {blockLen};
      
      for (auto i = 1ul; i < nNodes; ++i) {
        lens_.push_back(std::min(blockLen + accLens_.back(), len_) - accLens_.back()); 
        accLens_.push_back(accLens_.back() + lens_.back());
      }
      localLen_ = lens_[MPIRank(comm_)];
      localOffset_ = MPIRank(comm_) == 0 ? 0ul : accLens_[MPIRank(comm_) - 1];

      alloc();
    }
     
    // dividing the vector with inputs
    explicit DistributedVectors(MPI_Comm c,
                                const std::vector<size_t>& lens, size_t size):
        comm_(c), lens_(lens), size_(size) {
      
      size_t nNodes = MPISize(comm_);

      while (lens_.size() < nNodes) {
        lens_.push_back(0ul);
      }

      if (lens_.size() > nNodes ) {
        CErr("Can't distribute data to blocks more than number of nodes");
      }
      
      accLens_.resize(nNodes);
      accLens_[0] = lens[0];
      for (auto i = 1ul; i < nNodes; ++i) {
        accLens_[i] = accLens_[i - 1] + lens_[i];
      }
      len_ = accLens_.back(); 
      localLen_ = lens_[MPIRank(comm_)];
      localOffset_ = MPIRank(comm_) == 0 ? 0ul : accLens_[MPIRank(comm_) - 1];
      
      alloc();
    }
    
    // scatter data from root
    void scatter(size_t shift, size_t nVec, const _F* A, size_t LDA, int root) {
      this->sizeCheck(nVec + shift, "DistributedVectors<_F>::scatter");
      if (LDA < length()) {
        CErr("Can't scatter from a matrix with leading dimension smaller than the BlockVector length");
      }

      for (auto iVec = 0ul; iVec < nVec; ++iVec) {
        MPIScatterV(A + iVec * LDA, lens_, getLocalPtr(iVec + shift), localLen_, root, comm_);
      }
    }

    // gather data to root
    void gather(size_t shift, size_t nVec, _F* A, size_t LDA, int root) const {
      this->sizeCheck(nVec + shift, "DistributedVectors<_F>::gather");
      if (LDA < length()) {
        CErr("Can't gather to a matrix with leading dimension smaller than the BlockVector length");
      }
      for (auto iVec = 0ul; iVec < nVec; ++iVec) {
        MPIGatherV(getLocalPtr(iVec + shift), localLen_, A + iVec * LDA, lens_, root, comm_);
      }
    }
    
    // gather data to all process
    void allGather(size_t shift, size_t nVec, _F* A, size_t LDA) const {
      this->sizeCheck(nVec + shift, "DistributedVectors<_F>::allGather");
      if (LDA < length()) {
        CErr("Can't gather to a matrix with leading dimension smaller than the BlockVector length");
      }
      for (auto iVec = 0ul; iVec < nVec; ++iVec) {
        MPIAllGatherV(getLocalPtr(iVec + shift), localLen_, A + iVec * LDA,lens_, comm_);
      }
    }

    // set data from RawVectors
    void fromRawVectors(size_t shiftA, const RawVectors<_F>& vecs, size_t shiftB, size_t nVec) {
      vecs.sizeCheck(nVec + shiftB, "B during DistributedVectors<_F>::fromRawVectors");
      this->sizeCheck(nVec + shiftA, "A during DistributedVectors<_F>::fromRawVectors");
      scatter(shiftA, nVec, vecs.getPtr(shiftB), vecs.length(), 0);
    }

    explicit DistributedVectors(const RawVectors<_F>& vecs, size_t shift, size_t nVec):
        DistributedVectors(vecs.getMPIcomm(), vecs.length(), nVec) {
      fromRawVectors(0ul, vecs, shift, nVec); 
    }
   
    explicit DistributedVectors(const RawVectors<_F>& vecs):
        DistributedVectors(vecs, 0ul, vecs.size()) { }
     
    RawVectors<_F> toRawVectors(size_t shift, size_t nVec) const {
      RawVectors<_F> vecs(comm_, length(), nVec);
      setRawVectors(0ul, vecs, shift, nVec);
      return vecs;
    }
    
    void setRawVectors(size_t shiftA, RawVectors<_F>& vecs, size_t shiftB, size_t nVec) const {
      vecs.sizeCheck(nVec + shiftB, "B during DistributedVectors<_F>::setRawVectors");
      this->sizeCheck(nVec + shiftA, "A during DistributedVectors<_F>::setRawVectors");
      gather(shiftA, nVec, vecs.getPtr(shiftB), vecs.length(), 0);
    }
    
    RawVectors<_F> toRawVectors() const {
      return toRawVectors(0ul, size());
    }

    DistributedVectors(const DistributedVectors<_F> &other):
        DistributedVectors(other.comm_, other.other.lens_, other.size_) {
      set_data(0, size_, other, 0ul, false);
    }

    DistributedVectors(DistributedVectors<_F> &&other):
        comm_(other.comm_),
        data_(std::move(other.data_)), len_(other.len_), lens_(other.lens_), 
        accLens_(other.accLens_), size_(other.size_) { }

    virtual ~DistributedVectors() { dealloc(); }

    MPI_Comm getMPIcomm() const { return comm_; }

    void dealloc() {
      if (data_) CQMemManager::get().free(data_);
    }
    
    void alloc() {
      dealloc();
      //std::cout << "Allocating local data, with localLen_ = " << localLen_ << std::endl;
      // data_ = CQMemManager::get().malloc<_F>(localLength() * size_);
      if (not data_) {
        try {
          data_ = CQMemManager::get().malloc<_F>(localLength() * size_);;
        } catch (...) {
          std::cout << std::fixed;
          std::cout << "Insufficient memory for DistributedVectors object, "
                    <<  " (" << (localLength() * size_ / 1e9) * sizeof(_F) << " GB)"
                    << std::endl;
          std::cout << CQMemManager::get() << std::endl;
          CErr();
        }
      }

    }
    
    virtual size_t length() const override { return len_; }
    virtual size_t size() const override { return size_; }
    
    size_t lengthAtNode(size_t i) const { return lens_[i]; }
    size_t localLength() const { return localLen_; }
    size_t localOffset() const { return localOffset_; }
    
    _F* getLocalPtr(size_t i = 0ul) {
      this->sizeCheck(i, "DistributedVectors<_F>::getPtr");
      return data_ + i * localLength();
    }

    const _F* getLocalPtr(size_t i = 0ul) const {
      return const_cast<DistributedVectors*>(this)->getLocalPtr(i);
    }
    
    // Get element
    virtual _F get(size_t i, size_t j) const override {
      if (i >= length() or j >= size()) {
        CErr("Geting invalid place in DistributedVectors object.");
      }

      // find where it's stored
      size_t nodeId = std::distance(accLens_.begin(),
          std::upper_bound(accLens_.begin(), accLens_.end(), i));
      
      _F result = _F(0);

      if (MPIRank(comm_) == nodeId) {
        result =  getLocalPtr(j)[i - localOffset_];
      }
      MPIBCast(result, nodeId, comm_);
      
      return result;
    }

    // Set element
    virtual void set(size_t i, size_t j, _F value) override {
      if (i >= length() or j >= size()) {
        CErr("Seting invalid place in DistributedVectors object.");
      }

      if (i >= localOffset_ and i < accLens_[MPIRank(comm_)]) {
        getLocalPtr(j)[i - localOffset_] = value;
      }
    }

    using SolverVectors<_F>::clear;
    void clear(size_t shift, size_t nVec) override {
      if (nVec == 0) return;
      
      this->sizeCheck(shift + nVec, "DistributedVectors<_F>::clear");
      
      std::fill_n(getLocalPtr(shift), nVec * localLength(), 0.);
    }
     
    DistributedVectors<_F> copy(size_t shift, size_t nVec) const {
      
      this->sizeCheck(shift + nVec, "DistributedVectors<_F>::copy");
      
      DistributedVectors<_F> vecs(comm_, lens_, nVec);
      std::copy_n(getLocalPtr(shift), localLength() * nVec, vecs.getLocalPtr(0ul));
      
      return vecs;
    }

    // printfull on root 
    // but only blocks allocated in place for other nodes
    using SolverVectors<_F>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
     
      this->sizeCheck(shift + nVec, "DistributedVectors<_F>::print");
      std::string output_str = "";
      if (str != "") {
        output_str = "[" + str + "]";
      } else {
        output_str = "[DistributedVectors]";
      }
      
      output_str += " Block " + std::to_string(MPIRank(comm_)) + ", with offset = " + 
          std::to_string(localOffset_);
      prettyPrintSmart(out, output_str, getLocalPtr(shift), localLength(), nVec, localLength());
    }
    
    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 _F alpha, _F const *B, int64_t ldb,
                                 _F beta, SolverVectors<_F> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                             int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA = true) const override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<_F> &B, size_t shiftB) override;
    
    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB, bool moveable = false) override;

    using SolverVectors<_F>::scale;
    virtual void scale(_F scalar, size_t shift, size_t nVec) override;

    using SolverVectors<_F>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, _F alpha, const SolverVectors<_F> &X, size_t shiftX) override;

    // virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
    //                           size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, _F *R = nullptr, int LDR = 0) override;

    using SolverVectors<_F>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<_F>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;
  
    template <typename Compare>
    void getKIndicesAndValues(size_t K, size_t iVec, 
        std::vector<size_t>& kIndices, std::vector<_F>& kValues,
        Compare comp) const {
       
       std::vector<size_t> localIndices(localLength());
       std::iota(localIndices.begin(), localIndices.end(), 0ul);
       auto vecPointer = getLocalPtr(iVec);

       std::stable_sort(localIndices.begin(), localIndices.end(),
         [&] (size_t i, size_t j) {
           return comp(vecPointer[i], vecPointer[j]);
         });
       
       size_t nLocalK = std::min(K, localLength());
       
       if (MPISize() == 1) {
         for (auto i = 0ul; i < nLocalK; ++i) {
           kIndices.push_back(localIndices[i]);
           kValues.push_back(vecPointer[kIndices.back()]);
         }
         return;
       }
       
       /*
        * MPI case
        */
       std::vector<size_t> kLocalIndices;
       std::vector<_F> kLocalValues;
       for (auto i = 0ul; i < nLocalK; ++i) {
         kLocalIndices.push_back(localIndices[i] + localOffset_);
         kLocalValues.push_back(vecPointer[localIndices[i]]);
       }

       // reduction
       size_t nResultK = std::min(K, length());

       // gather sizes
       std::vector<size_t> recv_sizes = MPIGather(nLocalK, 0, this->comm_);
       size_t totalGatheredSize = (MPIRank(this->comm_) == 0) ? 
           std::accumulate(recv_sizes.begin(), recv_sizes.end(), size_t(0ul)) : 1ul;
       
       std::vector<size_t> gatheredIndices(totalGatheredSize);
       std::vector<_F> gatheredValues(totalGatheredSize);
       MPIGatherV(&kLocalIndices[0], nLocalK, &gatheredIndices[0], recv_sizes, 0, this->comm_);
       MPIGatherV(&kLocalValues[0], nLocalK, &gatheredValues[0], recv_sizes, 0, this->comm_);
       
       kIndices.resize(nResultK);
       kValues.resize(nResultK);
       if (MPIRank(this->comm_) == 0) {
         localIndices.resize(totalGatheredSize);
         std::iota(localIndices.begin(), localIndices.end(), 0ul);
         std::stable_sort(localIndices.begin(), localIndices.end(),
           [&] (size_t i, size_t j) {
             return comp(gatheredValues[i], gatheredValues[j]);
           });
         
         for (auto i = 0ul; i < nResultK; ++i) {
           kIndices[i] = gatheredIndices[localIndices[i]];
           kValues[i] = gatheredValues[localIndices[i]];
         }
       }
       MPIBCast(&kIndices[0], nResultK, 0, this->comm_);
       MPIBCast(&kValues[0], nResultK, 0, this->comm_);
    } // getKIndicesAndValues
  
  }; // class DistributedVectors

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

    virtual bool isView() const override { return true; }

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

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, _F alpha, _F const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, _F *R = nullptr, int LDR = 0) override;

    using SolverVectors<_F>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<_F>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual ~SolverVectorsView() {}

  }; // class SolverVectorsView

  template <typename T, typename _F, typename Operation>
  void tryDowncastReferenceTo(SolverVectors<_F>& vecs, Operation op) {
    try {
      if (not vecs.isView()) {
        op(dynamic_cast<T&>(vecs), 0ul);
      } else {
        SolverVectorsView<_F>& vecs_view = dynamic_cast<SolverVectorsView<_F>&>(vecs);
        op(dynamic_cast<T&>(vecs_view.getVecs()), vecs_view.shift());
      }
    } catch (const std::bad_cast& e) {
      CErr("SolverVectors Downcast failed!"); 
    }
  }    
  
  template <typename T, typename _F, typename Operation>
  void tryDowncastReferenceTo(const SolverVectors<_F> & vecs, Operation op) {
    try {
      if (not vecs.isView()) {
        op(dynamic_cast<const T&>(vecs), 0ul);
      } else {
        const SolverVectorsView<_F>& vecs_view = dynamic_cast<const SolverVectorsView<_F>&>(vecs);
        op(dynamic_cast<const T&>(vecs_view.getVecs()), vecs_view.shift());
      }
    } catch (const std::bad_cast& e) {
      CErr("const SolverVectors Downcast failed!"); 
    }
  }
  
  template <typename _F>
  _F* tryGetRawVectorsPointer(SolverVectors<_F>& vecs, size_t shift = 0) {
    _F* pointer = nullptr;
    tryDowncastReferenceTo<RawVectors<_F>>(vecs,
        [&] (auto& vecsRef, size_t extraShift) {
          shift += extraShift;
          pointer = vecsRef.getPtr(shift);
        }
    );
    return pointer;
  }

  template <typename _F>
  const _F* tryGetRawVectorsPointer(const SolverVectors<_F>& vecs, size_t shift = 0) {
    const _F* pointer = nullptr;
    tryDowncastReferenceTo<RawVectors<_F>>(vecs,
        [&] (auto& vecsRef, size_t extraShift) {
          shift += extraShift;
          pointer = vecsRef.getPtr(shift);
        }
    );
    return pointer;
  }

  template <typename _F>
  _F* tryGetDistributedVectorsLocalPointer(SolverVectors<_F>& vecs, size_t shift = 0) {
    _F* pointer = nullptr;
    tryDowncastReferenceTo<DistributedVectors<_F>>(vecs,
        [&] (auto& vecsRef, size_t extraShift) {
          shift += extraShift;
          pointer = vecsRef.getLocalPtr(shift);
        }
    );
    return pointer;
  }

  template <typename _F>
  const _F* tryGetDistributedVectorsLocalPointer(const SolverVectors<_F>& vecs, size_t shift = 0) {
    const _F* pointer = nullptr;
    tryDowncastReferenceTo<DistributedVectors<_F>>(vecs,
        [&] (auto& vecsRef, size_t extraShift) {
          shift += extraShift;
          pointer = vecsRef.getLocalPtr(shift);
        }
    );
    return pointer;
  }

} // namespace ChronusQ
