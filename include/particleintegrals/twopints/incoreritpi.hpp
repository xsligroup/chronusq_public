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

#include <cmath>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cxxapi/output.hpp>
#include <particleintegrals/twopints/eri3j.hpp>
namespace ChronusQ {

  enum class CHOLESKY_ALG {
    TRADITIONAL,      // The traditional one-step algorithm
    DYNAMIC_ALL,      // Individual ERI and CD vectors are all dynamically controlled
    SPAN_FACTOR,      // Koch's span-factor algorithm
    SPAN_FACTOR_REUSE,// Span-factor algorithm that reuses residual matrix from previous steps
    DYNAMIC_ERI,      // Span-factor algorithm that reuses ERI vectors from previous steps
    READ_PIVOTS       // Read pivots from file
  };


  struct SpanFactorShellPair {
    std::pair<size_t,size_t> PQ;
    std::vector<size_t> indices;
    bool evaluated;

    bool operator<(const SpanFactorShellPair &rhs) const {
      return indices[0] < rhs.indices[0];
    }
  };


  template <typename IntsT>
  class InCoreRITPI : public InCoreTPI<IntsT> {

    template <typename IntsU>
    friend class InCoreRITPI;

  protected:
    size_t NBRI, NBNBRI;
    std::shared_ptr<ERI3JBase<IntsT>> eri3j_ = nullptr;
    MPI_Comm comm_ = MPI_COMM_NULL;
    bool distributed_ = false;
    bool redistribute_ = false;
    std::shared_ptr<cqmatrix::Matrix<IntsT>> twocenterERI_ = nullptr;// L=(P|Q)^-1/2
    bool saveRawERI_ = false; ///< Save raw ERI before contraction with (P|Q)^-1/2
    std::shared_ptr<ERI3JBase<IntsT>> rawERI3J_ = nullptr; ///< raw Electron-Electron repulsion integrals (P|rs), compound rs
    std::shared_ptr<cqmatrix::Matrix<IntsT>> rawERI2C_ = nullptr;// (P|Q)
    void saveRawERI3J();

  public:

    // Constructor
    InCoreRITPI() = delete;
    InCoreRITPI(size_t nb, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false,
                TPI_KERNEL kernel = TPI_KERNEL::Coulomb, double omega = 0.):
        InCoreTPI<IntsT>(nb, 0, kernel, omega), NBRI(0), NBNBRI(0), comm_(comm), distributed_(distributed), redistribute_(redistribute) {
      malloc();
    }
    InCoreRITPI(size_t nb, size_t nbri, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false,
                TPI_KERNEL kernel = TPI_KERNEL::Coulomb, double omega = 0.):
        InCoreTPI<IntsT>(nb, 0, kernel, omega), NBRI(nbri), NBNBRI(nb*nbri), comm_(comm), distributed_(distributed), redistribute_(redistribute) {
      malloc();
    }
    InCoreRITPI( const InCoreRITPI &other ):
        InCoreRITPI(other.nBasis(), other.nRIBasis(), MPI_COMM_NULL, false, false,
                    other.kernel(), other.rangeSeparationParameter()) {
      std::copy_n(other.eri3j_->data(), this->nBasis()*NBNBRI, eri3j_->data());
    }
    template <typename IntsU>
    InCoreRITPI( const InCoreRITPI<IntsU> &other, int = 0 ):
        InCoreRITPI(other.nBasis(), other.nRIBasis(), MPI_COMM_NULL, false, false,
                    other.kernel(), other.rangeSeparationParameter()) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real InCoreRITPI from a Complex one.");
      std::copy_n(other.eri3j_->data(), this->nBasis()*NBNBRI, eri3j_->data());
    }
    InCoreRITPI( InCoreRITPI &&other ) noexcept: InCoreTPI<IntsT>(std::move(other)),
        NBRI(other.NBRI), NBNBRI(other.NBNBRI), comm_(other.comm_), distributed_(other.distributed_), 
        eri3j_(std::move(other.eri3j_)), redistribute_(other.redistribute_) {
    }

    InCoreRITPI& operator=( const InCoreRITPI &other ) {
      if (this != &other) { // self-assignment check expected
        if (this->nBasis() != other.nBasis()) {
          this->NB = other.NB;
          NBRI = other.NBRI;
          NBNBRI = other.NBNBRI;
          comm_ = other.comm_;
          distributed_ = other.distributed_;
          redistribute_ = other.redistribute_;
          malloc(); // reallocate memory
        }
        std::copy_n(other.eri3j_->data(), this->nBasis()*NBNBRI, eri3j_->data());
      }
      return *this;
    }
    InCoreRITPI& operator=( InCoreRITPI &&other ) noexcept {
      if (this != &other) { // self-assignment check expected
        this->NB = other.NB;
        NBRI = other.NBRI;
        NBNBRI = other.NBNBRI;
        comm_ = other.comm_;
        distributed_ = other.distributed_;
        redistribute_ = other.redistribute_;
        eri3j_ = std::move(other.eri3j_);
      }
      return *this;
    }

    size_t nRIBasis() const { return NBRI; }
    void setNRIBasis(size_t nbri) {
      if(NBRI != nbri) {
        NBRI = nbri;
        NBNBRI = this->NB * NBRI;
        this->malloc();
      }
    }

    // Single element interfaces
    virtual IntsT operator()(size_t p, size_t q, size_t r, size_t s) const {
      return operator()(p+q*this->nBasis(), r+s*this->nBasis());
    }
    virtual IntsT operator()(size_t pq, size_t rs) const {
      return blas::dot(NBRI, &eri3j_->data()[pq*NBRI], 1, &eri3j_->data()[rs*NBRI], 1);
    }
    IntsT& operator()(size_t L, size_t p, size_t q) {
      return eri3j_->data()[L + p*NBRI + q*NBNBRI];
    }
    IntsT operator()(size_t L, size_t p, size_t q) const {
      return eri3j_->data()[L + p*NBRI + q*NBNBRI];
    }
    std::shared_ptr<ERI3JBase<IntsT>> eri3j() const { return eri3j_; }

    // Tensor direct access
    virtual IntsT* pointer() override { return eri3j_->data(); }
    virtual const IntsT* pointer() const override { return eri3j_->data(); }
    MPI_Comm comm() const { return comm_; }
    bool isDistributed() const { return distributed_; }
    bool redistribute() const { return redistribute_; }

    // Raw ERI access
    void setSaveRawERI(bool save) { saveRawERI_ = save; }
    void clearRawERI() { rawERI3J_ = nullptr; rawERI2C_ = nullptr; }
    IntsT* rawERI3J() { return rawERI3J_->data(); }
    const IntsT* rawERI3J() const { return rawERI3J_->data(); }
    std::shared_ptr<cqmatrix::Matrix<IntsT>> rawERI2C() const { return rawERI2C_; }

    // 2-index ERI
    std::shared_ptr<cqmatrix::Matrix<IntsT>> twoIndexERI() { return twocenterERI_; } // returns L=(P|Q)
    const std::shared_ptr<cqmatrix::Matrix<IntsT>> twoIndexERI() const { return twocenterERI_; } // returns L=(P|Q)

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) {
      CErr("AO integral evaluation is NOT implemented in super class InCoreRITPI.");
    }

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) {
      CErr("AO integral evaluation with two basis sets is NOT implemented in super class InCoreRITPI.");
    }

    static void halfInverse2CenterERI(cqmatrix::Matrix<IntsT> &S); ///< forms S^{-1/2}, destroys S

    void contract2CenterERI(bool useCompoundIndex = true); ///< forms S^{-1/2}(Q|ij)

    virtual void clear() { eri3j_->clear(); }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (s == "")
        out << "  Electron repulsion integral:" << std::endl;
      else
        out << "  ERI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "INCORE RI (Gemm)";
      out << std::endl;
      if (printFull) {
        out << bannerTop << std::endl;
        size_t NB = this->nBasis(), NBRI = nRIBasis();
        out << std::scientific << std::left << std::setprecision(8);
        for(auto L = 0ul; L < NBRI; L++)
        for(auto i = 0ul; i < NB; i++)
        for(auto j = 0ul; j < NB; j++){
          out << "    (" << L << "|" << i << "," << j << ")  ";
          out << (*this)(L,i,j) << std::endl;
        };
        out << bannerEnd << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
      InCoreTPI<IntsT>::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t NBRI_bcast = NBRI;
        size_t NBNBRI_bcast = NBNBRI;
        MPIBCast(NBRI_bcast,root,comm);
        MPIBCast(NBNBRI_bcast,root,comm);

        if (NBRI_bcast != NBRI or NBNBRI_bcast != NBNBRI) {
          NBRI = NBRI_bcast;
          NBNBRI = NBNBRI_bcast;
          malloc();
        }

        MPIBCast(eri3j_->data(),this->nBasis()*NBNBRI,root,comm);
      }
#endif

    }

    InCore4indexTPI<IntsT> to4indexERI() {

      InCore4indexTPI<IntsT> eri4i(this->nBasis());
      
      size_t NB2 = this->nBasis() * this->nBasis();
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB2,NB2,NBRI,IntsT(1.),pointer(),NBRI,
           pointer(),NBRI,IntsT(0.),eri4i.pointer(),NB2);
      return eri4i;
    }

    template <typename IntsU>
    InCoreRITPI<IntsU> spatialToSpinBlock() const;

    template <typename TransT>
    InCoreRITPI<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const;

    template <typename TransT, typename OutT>
    void subsetTransform(
        char TRANS, const TransT* T, int LDT,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        OutT* out, bool increment = false) const;

    void malloc() {
      if (distributed_)
        eri3j_ = std::make_shared<DistributedERI3J<IntsT>>(comm(),this->nBasis(), NBRI);
      else
        eri3j_ = std::make_shared<IncoreERI3J<IntsT>>(this->nBasis(), NBRI);
    }

    virtual ~InCoreRITPI() {}

  }; // class InCoreRITPI

  template <typename IntsT>
  class InCoreAuxBasisRIERI : public InCoreRITPI<IntsT> {

    template <typename IntsU>
    friend class InCoreAuxBasisRIERI;

  protected:
    std::shared_ptr<BasisSet> auxBasisSet_ = nullptr; ///< BasisSet for the GTO basis defintion

  public:

    // Constructor
    InCoreAuxBasisRIERI() = delete;
    InCoreAuxBasisRIERI(size_t nb, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false):
        InCoreRITPI<IntsT>(nb, comm, distributed, redistribute) {}
    InCoreAuxBasisRIERI(size_t nb, size_t nbri, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false):
        InCoreRITPI<IntsT>(nb, nbri, comm, distributed, redistribute) {}
    InCoreAuxBasisRIERI(size_t nb,
        std::shared_ptr<BasisSet> auxBasisSet, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false):
        InCoreRITPI<IntsT>(nb, auxBasisSet->nBasis, comm, distributed, redistribute),
        auxBasisSet_(auxBasisSet) {}
    InCoreAuxBasisRIERI( const InCoreAuxBasisRIERI& ) = default;
    template <typename IntsU>
    InCoreAuxBasisRIERI( const InCoreAuxBasisRIERI<IntsU> &other, int = 0 ):
        InCoreRITPI<IntsT>(other) {
      auxBasisSet_ = other.auxBasisSet_;
    }
    InCoreAuxBasisRIERI( InCoreAuxBasisRIERI &&other ) = default;

    InCoreAuxBasisRIERI& operator=( const InCoreAuxBasisRIERI& ) = default;
    InCoreAuxBasisRIERI& operator=( InCoreAuxBasisRIERI&& ) = default;

    void setAuxBasisSet(std::shared_ptr<BasisSet> auxbasisSet) {
      auxBasisSet_ = auxbasisSet;
      this->setNRIBasis(auxBasisSet_->nBasis);
    }
    std::shared_ptr<BasisSet> auxbasisSet() const { return auxBasisSet_; }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);

    void compute3CenterERI(BasisSet &basisSet, BasisSet &auxBasisSet);

    void compute2CenterERI(BasisSet &auxBasisSet, IntsT *S) const;

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (s == "")
        out << "  Electron repulsion integral:" << std::endl;
      else
        out << "  ERI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "INCORE auxiliary basis RI (Gemm)";
      out << std::endl;
      if (printFull) {
        out << bannerTop << std::endl;
        size_t NB = this->nBasis(), NBRI = this->nRIBasis();
        out << std::scientific << std::left << std::setprecision(8);
        for(auto L = 0ul; L < NBRI; L++)
        for(auto i = 0ul; i < NB; i++)
        for(auto j = 0ul; j < NB; j++){
          out << "    (" << L << "|" << i << "," << j << ")  ";
          out << (*this)(L,i,j) << std::endl;
        };
        out << bannerEnd << std::endl;
      }
    }

    virtual ~InCoreAuxBasisRIERI() {}

  }; // class InCoreAuxBasisRIERI

  template <typename IntsT>
  class InCoreCholeskyRIERI : public InCoreRITPI<IntsT> {

    template <typename IntsU>
    friend class InCoreCholeskyRIERI;

  public:

    using Kernel = TPI_KERNEL;

  protected:
    CHOLESKY_ALG alg_; // Algorithm for building 3-index Cholesky vectors
    double tau_; // Maximum error allowed in the decomposition
    double sigma_; // sigma in the span factor algorithm
    size_t maxQual_; // Max qualified candidates per span
    size_t minShrinkCycle_; // Mininum # of iterations between shrinks
    bool generalContraction_; // Uncontract basis then run
    bool build4I_ = false; // Explicitly build four-index TPI
    std::shared_ptr<InCore4indexTPI<IntsT>> eri4I_ = nullptr; // Four-index TPI
    std::vector<size_t> pivots_; // List of selected pivots
    std::vector<std::vector<libint2::Shell>> shellPrims_; // Mappings from primitives to CGTOs
    std::vector<IntsT*> coefBlocks_; // Mappings from primitives to CGTOs
    bool updatePivots_ = true; // Update pivots

    // Libint
    size_t maxNcontrAMSize_ = 1, maxNprimAMSize_ = 1, maxAMSize_ = 1;
    std::vector<libint2::Engine> engines;
    std::vector<double*> workBlocks;

    // Libcint
    int *atm, *bas;
    double *env, *buffAll, *cacheAll;
    size_t buffN4, cache_size;
    int nAtoms, nShells;
    bool libcint_ = false; // Using libcint

    // ERI stat
    size_t c1ERI = 0, c2ERI = 0;
    double t1ERI = 0.0, t2ERI = 0.0;

  public:

    // Constructor
    InCoreCholeskyRIERI() = delete;
    InCoreCholeskyRIERI(size_t nb, double tau,
        CHOLESKY_ALG alg = CHOLESKY_ALG::DYNAMIC_ERI, bool genContr = true,
        double sigma = 0.01, size_t maxQual = 1000, size_t minShrink = 10,
        bool build4I = false, MPI_Comm comm = MPI_COMM_NULL, bool distributed = false, bool redistribute = false,
        Kernel kernel = Kernel::Coulomb, double omega = 0.):
        InCoreRITPI<IntsT>(nb, comm, distributed, redistribute, kernel, omega), alg_(alg), tau_(tau),
        sigma_(sigma), maxQual_(maxQual), minShrinkCycle_(minShrink),
        generalContraction_(genContr), build4I_(build4I) {}

    InCoreCholeskyRIERI( const InCoreCholeskyRIERI& ) = default;
    template <typename IntsU>
    InCoreCholeskyRIERI( const InCoreCholeskyRIERI<IntsU> &other, int = 0 ):
        InCoreRITPI<IntsT>(other), alg_(other.alg_), tau_(other.tau_),
        sigma_(other.sigma_), maxQual_(other.maxQual_),
        minShrinkCycle_(other.minShrinkCycle_),
        generalContraction_(other.generalContraction_),
        build4I_(other.build4I_),
        eri4I_(std::make_shared<InCore4indexTPI<IntsT>>(*other.eri4I_)),
        pivots_(other.pivots_){}
    InCoreCholeskyRIERI( InCoreCholeskyRIERI &&other ) = default;

    InCoreCholeskyRIERI& operator=( const InCoreCholeskyRIERI& ) = default;
    InCoreCholeskyRIERI& operator=( InCoreCholeskyRIERI&& ) = default;

    void setTau( double tau ) { tau_ = tau; }
    double tau() const { return tau_; }
    //std::vector<size_t> pivots() { return pivots_; }// List of selected pivots
    void setPivots(const std::vector<size_t>& pivots) { pivots_ = pivots; }
    const std::vector<size_t>& pivots() const { return pivots_; }// List of selected pivots
    void setUpdatePivots(bool updatePivots) { updatePivots_ = updatePivots; }
    bool updatePivots() const { return updatePivots_; }

    libint2::Operator libintOperator() const {
      return ChronusQ::libintOperator(this->kernel());
    }

    // Rebuild the Cholesky RI TPI over a different interaction kernel
    std::shared_ptr<TwoPInts<IntsT>> createWithKernel(TPI_KERNEL kernel, double omega) const override {
      return std::make_shared<InCoreCholeskyRIERI<IntsT>>(
          this->nBasis(), tau_, alg_, generalContraction_,
          sigma_, maxQual_, minShrinkCycle_, build4I_, this->comm(),
          this->isDistributed(), this->redistribute(), kernel, omega);
    }

    // 4-index ERI direct access
    void setFourIndexERI(std::shared_ptr<InCore4indexTPI<IntsT>> eri4I) {
      eri4I_ = eri4I;
    }
    std::shared_ptr<InCore4indexTPI<IntsT>> fourIndexERI() { return eri4I_; }
    const std::shared_ptr<InCore4indexTPI<IntsT>> fourIndexERI() const { return eri4I_; }


    // Computation interfaces

    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
                               OPERATOR, const HamiltonianOptions&);

    void computeCD_Traditional(BasisSet&);

    void computeCDPivots_DynamicAll(BasisSet&);

    void computeCDPivots_SpanFactor(BasisSet&);

    void computeCDPivots_DynamicERI(BasisSet&);

    void computeCDPivots_SpanFactorReuse(BasisSet&);

    std::map<std::pair<size_t, size_t>, IntsT*>
    computeDiagonalLibint(BasisSet&, IntsT *diag, bool saveDiagBlocks = false);

    std::map<std::pair<size_t, size_t>, IntsT*>
    computeDiagonalLibcint(BasisSet&, IntsT *diag, bool saveDiagBlocks = false);

    void computeTraditionalERIVector(BasisSet&, size_t pivot, IntsT *L);

    void computeTraditionalERIVectorByShellLibint(BasisSet&,
        size_t r, size_t s, size_t R, size_t S, IntsT *L);

    void computeTraditionalERIVectorByShellLibcint(BasisSet&,
        size_t r, size_t s, size_t R, size_t S, IntsT *L);

    void computeSpanFactorERI(BasisSet&,
        const std::vector<size_t> &D,
        std::vector<SpanFactorShellPair> &shellPair_Dindices,
        size_t lenQ, IntsT* ERIvecAlloc,
        std::map<std::pair<size_t, size_t>, IntsT*> &diagBlocks,
        bool Mpq_only = false);

    void computeSpanFactorERIVectorLibint(BasisSet&,
        const std::vector<size_t> &D,
        std::vector<SpanFactorShellPair> &shellPair_Dindices,
        size_t lenQ, size_t qShellPairIndex, IntsT* ERIvecAlloc,
        std::map<std::pair<size_t, size_t>, IntsT*> &diagBlocks,
        bool Mpq_only = false);

    void computeSpanFactorERIVectorLibcint(BasisSet&,
        const std::vector<size_t> &D,
        std::vector<SpanFactorShellPair> &shellPair_Dindices,
        size_t lenQ, size_t qShellPairIndex, IntsT* ERIvecAlloc,
        std::map<std::pair<size_t, size_t>, IntsT*> &diagBlocks,
        bool Mpq_only = false);

    void computePivotRI(BasisSet &basisSet);

    void computePivotRI3indexERILibint(BasisSet &basisSet, bool useCompoundIndex = true);
    void computePivotRI3indexERILibintMPI(BasisSet &basisSet, bool useCompoundIndex = true);

    void computePivotRI3indexERILibcint(BasisSet &basisSet, bool useCompoundIndex = true);
    void computePivotRI3indexERILibcintMPI(BasisSet &basisSet, bool useCompoundIndex = true);

    static void extractTwoCenterSubsetFrom3indexERI(
        const std::vector<size_t> &pivots, size_t NBRI, size_t NB,
        const double* eri3J, size_t LD3J, double* S, size_t LDS,
        bool upperTriOnly = true, bool useCompoundIndex = true);
    static void extractTwoCenterSubsetFrom3indexERIMPI(
        const std::vector<size_t> &pivots, size_t NBRI, size_t NB,
        const DistributedERI3J<double>& eri3J, double* S, size_t LDS,
        bool upperTriOnly = true);

    static std::map<std::pair<size_t,size_t>, std::vector<size_t>>
    groupPivotsByShell(BasisSet&, const std::vector<size_t> &pivots);

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (s == "")
        out << "  Electron repulsion integral:" << std::endl;
      else
        out << "  ERI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "INCORE Cholesky decomposition RI (Gemm)";
      out << std::endl;
      if (this->kernel() == Kernel::ShortRangeErfc)
        out << "    * Short-range erfc kernel, omega = " << this->rangeSeparationParameter() << std::endl;
      if (printFull) {
        out << bannerTop << std::endl;
        size_t NB = this->nBasis(), NBRI = this->nRIBasis();
        out << std::scientific << std::left << std::setprecision(8);
        for(auto L = 0ul; L < NBRI; L++)
        for(auto i = 0ul; i < NB; i++)
        for(auto j = 0ul; j < NB; j++){
          out << "    (" << L << "|" << i << "," << j << ")  ";
          out << (*this)(L,i,j) << std::endl;
        };
        out << bannerEnd << std::endl;
      }
    }

    virtual ~InCoreCholeskyRIERI() {}

  }; // class InCoreCholeskyRIERI



  template <typename MatsT, typename IntsT>
  class RITPIContraction : public InCoreTPIContraction<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class RITPIContraction;

  public:

    // Constructors

    RITPIContraction() = delete;
    RITPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      InCoreTPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    RITPIContraction(
      const RITPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      RITPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
      this->isCross = other.isCross;
    }
    template <typename MatsU>
    RITPIContraction(
      RITPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      RITPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
      this->isCross = other.isCross;
    }

    RITPIContraction( const RITPIContraction &other ):
      RITPIContraction(other, 0) {
      this->contractSecond = other.contractSecond;
      this->isCross = other.isCross;
    }
    RITPIContraction( RITPIContraction &&other ):
      RITPIContraction(std::move(other), 0) {
      this->contractSecond = other.contractSecond;
      this->isCross = other.isCross;
    }

    virtual ~RITPIContraction() {}

  }; // class RITPIContraction



  template <typename MatsT, typename IntsT>
  class InCoreRITPIContraction : public RITPIContraction<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCoreRITPIContraction;

  public:

    // Constructors

    InCoreRITPIContraction() = delete;
    InCoreRITPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      RITPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCoreRITPIContraction(
        const InCoreRITPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCoreRITPIContraction(other.ints_) {}
    template <typename MatsU>
    InCoreRITPIContraction(
        InCoreRITPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCoreRITPIContraction(other.ints_) {}

    InCoreRITPIContraction( const InCoreRITPIContraction &other ):
      InCoreRITPIContraction(other, 0) {}
    InCoreRITPIContraction( InCoreRITPIContraction &&other ):
      InCoreRITPIContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    void KCoefContract(
        MPI_Comm, size_t nO, MatsT *X, MatsT *AX) const;

    virtual ~InCoreRITPIContraction() {}

  
  private:

    void KContract_real_impl(
        MPI_Comm, const double *Xr, double *AXr, double *Ktemp) const;

  }; // class InCoreRITPIContraction

}; // namespace ChronusQ
