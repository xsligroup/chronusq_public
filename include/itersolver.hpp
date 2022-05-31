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

#include <chronusq_sys.hpp>
#include <memmanager.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/solve.hpp>
#include <cerr.hpp>

#include <util/mpi.hpp>
#include <util/math.hpp>
#include <util/matout.hpp>

#include <iostream>

namespace ChronusQ {


  template <typename _F>
  class SolverVectors {

  public:
    // Length of each vector
    virtual size_t length() const = 0;
    // Number of vectors in container
    virtual size_t size() const = 0;
    // Get raw pointer of vectors, only works for RawVectors
    virtual _F* getPtr(size_t i = 0) const = 0;
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
    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB) = 0;

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
    double norm2F(size_t shift) const {
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
    double maxNormElement(size_t shift) const {
      return maxNormElement(shift, size() - shift);
    }
    virtual double maxNormElement(size_t shift, size_t nVec) const = 0;

    virtual ~SolverVectors() {}

  };


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
    virtual _F* getPtr(size_t i = 0) const override {
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
      ROOT_ONLY(comm_);
      if (nVec == std::numeric_limits<size_t>::max()) {
        getPtr(shift);
        nVec = size() - shift;
      } else
        getPtr(shift + nVec - 1);
      std::fill_n(getPtr(shift), nVec*length(), 0.);
    }

    using SolverVectors<_F>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      ROOT_ONLY(comm_);
      if (nVec == std::numeric_limits<size_t>::max()) {
        getPtr(shift);
        nVec = size() - shift;
      } else
        getPtr(shift + nVec - 1);
      prettyPrintSmart(out, str, getPtr(shift), length(), nVec, length());
    }

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 _F alpha, _F const *B, int64_t ldb,
                                 _F beta, SolverVectors<_F> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<_F> &B, size_t shiftB,
                             int64_t m, int64_t n, _F *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB) override;

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

  };


  template <typename _F>
  class SolverVectorsView : public SolverVectors<_F> {

  protected:

    SolverVectors<_F> &vecs_;
    size_t shift_;

  public:
    SolverVectorsView(SolverVectors<_F> &vecs, size_t shift = 0)
        : vecs_(vecs), shift_(shift) {
      if (shift_ >= vecs_.size())
        CErr("Creating a view out of the vectors' capacity.");
    }
    SolverVectorsView(SolverVectorsView<_F> &view, size_t shift = 0)
        : vecs_(view.vecs_), shift_(shift + view.shift_) {
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

    virtual _F* getPtr(size_t i = 0) const override {
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

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<_F> &B, size_t shiftB) override;

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

  };


  template <typename _F>
  class IterSolver {

  protected:

    typedef std::function< std::shared_ptr<SolverVectors<_F>>(size_t) >            VecsGen_t;
    typedef std::function< void(size_t,SolverVectors<_F>&,SolverVectors<_F>&) >    LinearTrans_t;
    typedef std::function< void(size_t,_F,SolverVectors<_F>&,SolverVectors<_F>&) > Shift_t;

    CQMemManager &memManager_;
    MPI_Comm     comm_;

    size_t N_;       ///< Problem dimension
    size_t mSS_;     ///< Max subspace dimension
    size_t maxMacroIter_; ///< Maximum number of Macro iterations
    size_t maxMicroIter_; ///< Maximum number of Micro iterations


    double convCrit_; ///< Convergence criteria

    VecsGen_t vecGen_; ///< Vector generator

    LinearTrans_t linearTrans_;    ///< AX Product
    LinearTrans_t preCondNoShift_; ///< Unshifted preconditioner

    Shift_t shiftVec_;       ///< (A - sB)X given AX
    Shift_t preCondWShift_;  ///< Shifted preconditioner

    inline std::shared_ptr<SolverVectors<_F>> rawVecsGenerator_(size_t nVec) {

      return std::make_shared<RawVectors<_F>>(comm_, this->memManager_, this->N_, nVec);

    };

    inline void defaultShiftVec_(size_t nVec, _F shift, SolverVectors<_F> &V, SolverVectors<_F> &AV) {

      AV.axpy(0, nVec, shift, V, 0);

    };

    inline void defaultPreCondShift_(size_t nVec, _F shift, SolverVectors<_F> &V, SolverVectors<_F> &AV) {

      preCondNoShift_(nVec,V,AV);
      if(std::abs(std::abs(shift)) > 1e-15) shiftVec_(nVec,-1./shift,V,AV);

    };

  public:

    // Constructor (unshifted preconditioner)
    IterSolver(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
        comm_(c), memManager_(mem), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), vecGen_(vecGen),
        linearTrans_(linearTrans), shiftVec_(shiftVec),
        preCondNoShift_(preNoShift) {

      using namespace std::placeholders;
      if( not vecGen_ ) {
        vecGen_ =
            std::bind(&IterSolver<_F>::rawVecsGenerator_,this,_1);
      }

      if( not shiftVec_ ) {
        shiftVec_ =
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4);
      }

      preCondWShift_ =
        std::bind(&IterSolver<_F>::defaultPreCondShift_,this,_1,_2,_3,_4);

    }

    // Constructor (shifted preconditioner)
    IterSolver(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
        comm_(c), memManager_(mem), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), vecGen_(vecGen),
        linearTrans_(linearTrans), shiftVec_(shiftVec),
        preCondWShift_(preShift) {

      using namespace std::placeholders;
      if( not vecGen_ ) {
        vecGen_ =
            std::bind(&IterSolver<_F>::rawVecsGenerator_,this,_1);
      }

      if( not shiftVec_ ) {
        shiftVec_ =
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4);
      }

    }



  };




  template <typename _F>
  class IterLinearSolver : public IterSolver<_F> {

  protected:

    size_t nRHS_;
    std::shared_ptr<SolverVectors<_F>> RHS_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> SOL_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> V_   = nullptr;
    std::shared_ptr<SolverVectors<_F>> AV_  = nullptr;
    std::shared_ptr<SolverVectors<_F>> RES_ = nullptr;

    std::vector<_F> shifts_;

    std::vector<double>              rhsNorm_;
    std::vector<std::vector<double>> resNorm_;

  public:


    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    size_t shiftBS = 1;
    size_t rhsBS   = 1;

    IterLinearSolver(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preNoShift,
                    vecGen,shiftVec){ }

    IterLinearSolver(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preShift,
                    vecGen,shiftVec){ }


    ~IterLinearSolver() {};




    template <typename T>
    void setRHS(size_t nRHS, T*RHS, size_t LDRHS);

    template <typename T>
    void setShifts(size_t nShift, T* shifts);


    std::shared_ptr<SolverVectors<_F>> getSol() const {
      return SOL_;
    };




    virtual void alloc() {

//      size_t nRHSnSN       = nRHS_ * shifts_.size() * this->N_;
//      size_t nRHSnSN_batch = rhsBS * shiftBS * this->N_;

      SOL_ = this->vecGen_(nRHS_ * shifts_.size());
      V_   = this->vecGen_(rhsBS * shiftBS * this->mSS_);
      AV_  = this->vecGen_(rhsBS * shiftBS * this->mSS_);
      RES_  = this->vecGen_(rhsBS * shiftBS);

    }



    void run();

    virtual void runBatch(size_t nRHS, size_t nShift,
                          std::shared_ptr<SolverVectors<_F>> RHS, _F *shifts,
                          std::shared_ptr<SolverVectors<_F>> SOL, double *RHSNorm);

  };




  template <typename _F>
  class IterDiagonalizer : public IterSolver<_F> {


  protected:

    size_t nRoots_;
    size_t nGuess_;
    std::shared_ptr<SolverVectors<_F>> VR_   = nullptr;
    _F     *VL_   = nullptr;
    std::shared_ptr<SolverVectors<_F>> AVR_  = nullptr;
    _F     *AVL_  = nullptr;
    _F     *RESR_ = nullptr;
    _F     *RESL_ = nullptr;


    dcomplex * eigVal_ = nullptr;

    void RayleighRitz();



  public:

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    IterDiagonalizer(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      size_t nR,
      size_t nG,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
      IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                     preNoShift,vecGen,shiftVec),
        nRoots_(nR), nGuess_(nG){ }

    IterDiagonalizer(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      size_t nR,
      size_t nG,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,mem,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                    preShift,vecGen,shiftVec),
       nRoots_(nR), nGuess_(nG){ }


    ~IterDiagonalizer() {

      if(VL_)   this->memManager_.free(VL_);
      if(AVL_)  this->memManager_.free(AVL_);
      if(RESR_) this->memManager_.free(RESR_);
      if(RESL_) this->memManager_.free(RESL_);
      if(eigVal_) this->memManager_.free(eigVal_);

    };




    virtual void alloc() {

      VR_     = this->vecGen_(this->mSS_);

      // NO MPI
      // ROOT_ONLY(this->comm_);

//      size_t NNR = this->N_ * this->mSS_;

      eigVal_ = this->memManager_.template malloc<dcomplex>(this->nRoots_);
      //AVR_    = this->memManager_.template malloc<_F>(NNR);
      //RESR_   = this->memManager_.template malloc<_F>(NNR);

    }


    const dcomplex* eigVal() const { return eigVal_; }
    std::shared_ptr<SolverVectors<_F>> VR() const { return VR_; }


    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)>) = 0;

    void run();

    virtual bool runMicro() = 0;
    virtual void restart() = 0;

  };








  // GPLHR
  template <typename _F>
  class GPLHR : public IterDiagonalizer<_F> {

    double *RelRes = nullptr;
    std::shared_ptr<SolverVectors<_F>> Guess = nullptr;


  protected:


    bool useAdaptiveSigma() const { return false; }


    void getResidualNorms(size_t N, size_t nR, const SolverVectors<_F> &WMAT,
                          double *RelRes, dcomplex *LAMBDA, double nrmA);

    void getTriU(size_t N, _F *TRIUA, size_t LDTRIUA, _F *TRIUB,
      size_t LDTRIUB, _F *MA, size_t LDMA, _F *MB, size_t LDMB);

    bool checkConv(size_t nR, double *RelRes);




    void halfProj(size_t nV, size_t nS, const SolverVectors<_F> &V,
                  SolverVectors<_F> &S, _F *smSCR, size_t LDSCR);

    /**
     *  \brief Overload of halfProj where nV == nS
     */
    inline void halfProj(size_t nV, const SolverVectors<_F> &V,
                         SolverVectors<_F> &S, _F *smSCR, size_t LDSCR) {

      halfProj(nV,nV,V,S,smSCR,LDSCR);

    }

    void halfProj2(size_t nV, size_t nS, const SolverVectors<_F> &V, const SolverVectors<_F> &AV,
                   SolverVectors<_F> &S, SolverVectors<_F> &AS, _F *smSCR, size_t LDSCR);

    /**
     *  \brief Overload of halfProj2 where nV == nS
     */
    inline void halfProj2(size_t nV, const SolverVectors<_F> &V, const SolverVectors<_F> &AV,
                          SolverVectors<_F> &S, SolverVectors<_F> &AS, _F *smSCR, size_t LDSCR) {

      halfProj2(nV,nV,V,AV,S,AS,smSCR,LDSCR);

    }


    void newSMatrix(size_t nR, const SolverVectors<_F> &V, const SolverVectors<_F> &Q,
                    SolverVectors<_F> &S, _F *smSCR, size_t LDSCR);


  public:


    size_t m = 1;
    _F sigma = 0.;
    _F hardLim = -std::numeric_limits<double>::infinity();

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    GPLHR(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MAXITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
      IterDiagonalizer<_F>(c,mem,N,2 + m*nR,1,MAXITER,conv,nR,nR,
                           linearTrans,preNoShift,vecGen,shiftVec){ }

    GPLHR(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MAXITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
      IterDiagonalizer<_F>(c,mem,N,2 + m*nR,1,MAXITER,conv,nR,nR,
                           linearTrans,preShift,vecGen,shiftVec){ }

    ~GPLHR() {

      if( RelRes ) this->memManager_.free(RelRes);

    }

    void setM(size_t _m) {

      m = _m;
      this->mSS_ = (3 + m)*this->nRoots_;

    }

    void alloc() {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      // ROOT_ONLY(this->comm_);


      // Allocate GPLHR specific Memory

      this->RelRes = this->memManager_.template malloc<double>(this->nRoots_);

    }

    bool runMicro();

    void restart();

    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)> func) override {

      if( nGuess != this->nRoots_ )
        CErr("GPLHR Requires nGuess = nRoots",std::cout);

      Guess = this->vecGen_(nGuess);

      // NO MPI
      // ROOT_ONLY(this->comm_);

      func(nGuess, *Guess, this->N_);

    }

  }; // GPLHR





  // Davidson
  template <typename _F>
  class Davidson : public IterDiagonalizer<_F> {

    bool DoLeftEigVec   = false;
    bool sortByDistance = false;
    size_t GramSchmidt_NRe = 1;
    double GramSchmidt_eps = 1e-12;

    // Energy specific options
    bool EnergySpecific = false;
    size_t nHighERoots  = 0;
    size_t  nLowERoots  = 0;
    double EnergyRef    = 0.;
    bool adaptiveERef   = false; // use loweset eigenvalues at current iteration

    double   *RelRes  = nullptr;
    std::shared_ptr<SolverVectors<_F>> Guess = nullptr;
    dcomplex *EigForT = nullptr; // Eigenvalues can be used for preconditioner

    // creating double(dcomplex) operator to accomodate the
    // complex eigenvalue to real arithmetric
    _F dcomplexTo_F(dcomplex &a) { return * reinterpret_cast<_F*>(&a);}

  protected:

  public:

    size_t m = 50;
    size_t whenSc = 2;
    size_t kG = 3;

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    Davidson(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MAXMACROITER,
      const size_t MAXMICROITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()):
      IterDiagonalizer<_F>(c,mem,N,m*nR,MAXMACROITER,MAXMICROITER,conv,nR,nR*kG,
                           linearTrans,preNoShift,vecGen,shiftVec){ }

    ~Davidson() {

      if( RelRes  ) this->memManager_.free(RelRes);

    }

    void setM(size_t _m) {
      m = _m;
      this->mSS_ = m*this->nRoots_;
    }

    void setkG(size_t _kG) {
      kG = _kG;
      this->nGuess_ = kG*this->nRoots_;
    }

    void setWhenSc(size_t _WhenSc) { whenSc = _WhenSc;}

    void setGramSchmidtRepeat(size_t nRe) { GramSchmidt_NRe = nRe;}

    void setGramSchmidtEps(double eps) { GramSchmidt_eps = eps;}

    void setEigForT(dcomplex * _Eig) {

      if( this->memManager_.getSize(_Eig) < this->nGuess_ )
        CErr("Davison EigForT requires a memory block with size at least nGuess ",std::cout);

      EigForT = _Eig;
    }

    void alloc() {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      // ROOT_ONLY(this->comm_);

      // Allocate Davidson specific Memory
      this->RelRes  = this->memManager_.template malloc<double>(this->nGuess_);
    }

    bool runMicro();

    void restart();

    void useEnergySpecific(size_t nHER, double energyThd, double ERef) {
      assert (this->nRoots_ >= nHER);
//      assert (energyThd > 0);

      this->EnergySpecific = true;
      this->nHighERoots    = nHER;
      this->nLowERoots     = this->nRoots_ - nHER;
      this->EnergyRef      = ERef;
    }

    void useEnergySpecific(size_t nHER, double energyThd) {
      this->adaptiveERef = true;
      this->useEnergySpecific(nHER, energyThd, 0.);
    }

    void doLeftEigenvector() {this->DoLeftEigVec = true; }

    void setSortByDistance(){
      this->sortByDistance = true;
    }

    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)> func) override {

      if( nGuess != this->nGuess_ )
        CErr("Davison Requires nGuess = nGuess_",std::cout);

      Guess = this->vecGen_(nGuess);

      // NO MPI
      // ROOT_ONLY(this->comm_);

      func(nGuess, *Guess, this->N_);

    }


  }; // Davidson


  template <typename _F>
  class GMRES : public IterLinearSolver<_F> {


    std::shared_ptr<SolverVectors<_F>> HHR_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> U_   = nullptr;
    _F * W_   = nullptr;
    _F * J_   = nullptr;
    _F * R_   = nullptr;

  public:

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterLinearSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterLinearSolver<_F>::Shift_t;


    GMRES(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterLinearSolver<_F>(c,mem,N,MSS,1,MSS,conv,linearTrans,preNoShift,
                          vecGen,shiftVec){ }

    GMRES(
      MPI_Comm c,
      CQMemManager &mem,
      const size_t N,
      const size_t MSS,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterLinearSolver<_F>(c,mem,N,MSS,1,MSS,conv,linearTrans,preShift,
                          vecGen,shiftVec){ }


    ~GMRES() {

      if(W_  )  this->memManager_.free(W_);
      if(J_  )  this->memManager_.free(J_);
      if(R_  )  this->memManager_.free(R_);

    }


    void alloc() {

      // Standard linear solver allocations
      IterLinearSolver<_F>::alloc();

      size_t nBatch = this->rhsBS * this->shiftBS;
      size_t MSSnBatch = this->mSS_ * nBatch;

      HHR_ = this->vecGen_(nBatch);
      U_   = this->vecGen_(MSSnBatch);

      // No MPI for GMRES
      // ROOT_ONLY(this->comm_);

      J_   = this->memManager_.template malloc<_F>(2 * MSSnBatch);
      R_   = this->memManager_.template malloc<_F>(MSSnBatch * this->mSS_);
      W_   = this->memManager_.template malloc<_F>(MSSnBatch + nBatch);

    }

    void runBatch(size_t nRHS, size_t nShift,
                  std::shared_ptr<SolverVectors<_F>> RHS, _F *shifts,
                  std::shared_ptr<SolverVectors<_F>> SOL, double *RHSNorm );

  };

  

};

