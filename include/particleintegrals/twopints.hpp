/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <particleintegrals.hpp>
#include <matrix.hpp>
#include <util/mpi.hpp>

namespace ChronusQ {

  /// Two-particle interaction kernel: full Coulomb or short-range erfc-attenuated
  enum class TPI_KERNEL {
    Coulomb,
    ShortRangeErfc
  };

  enum TWOBODY_CONTRACTION_TYPE {
    COULOMB, ///< (mn | kl) X(lk)
    EXCHANGE,///< (mn | kl) X(nk)
    PAIR,    ///< (mn | kl) X(nl)
    BARE_COULOMB,
    LLLL,
    LLSS,
    SSSS,
    GAUNT,
    GAUGE,
    DC_COULOMB,
    DC_EXCHANGE,
    SSSS_COULOMB,
    SSSS_EXCHANGE,
    GAUNT_COULOMB,
    GAUNT_EXCHANGE,
    GAUGE_COULOMB,
    GAUGE_EXCHANGE,
  }; ///< 2-Body Tensor Contraction Specification


  // ERI transpose type
  enum INTEGRAL_TRANSPOSE {
    TRANS_NONE,
    TRANS_MN_TRANS_KL,
    TRANS_MNKL,
    TRANS_KL,
    TRANS_MN
  };


  /**
   *  The TwoBodyContraction struct. Stores information
   *  pertinant for a two body operator contraction with
   *  a one body (2 index) operator. z.B. The density matrix.
   */
  template <typename T>
  struct TwoBodyContraction {


    T*  X;  ///< 1-Body (2 index) operator to contraction
    T*  AX; ///< 1-Body (2 index) storage for the contraction

    bool HER; ///< Whether or not X is hermetian

    TWOBODY_CONTRACTION_TYPE contType;

    int ERI4Ind = -1;

    INTEGRAL_TRANSPOSE intTrans;

    double* ERI4 = nullptr;


  }; // struct TwoBodyContraction
  
  /**
   *  The RelTwoBodyContraction struct. Stores information
   *  pertinant for a relativisitc two body operator contraction with
   *  a one body (2 index) operator. z.B. The density matrix.
   */
  template <typename T>
  struct TwoBodyRelContraction {

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<T>>  X;  ///< 1-Body (2 index) operator to contraction
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<T>>  AX; ///< 1-Body (2 index) storage for the contraction

    bool HER; ///< Whether or not X is hermetian

    TWOBODY_CONTRACTION_TYPE contType;

    double* ERI4 = nullptr;

    INTEGRAL_TRANSPOSE intTrans;

  }; // struct RelTwoBodyContraction

  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  electron-electron repulsion integral tensors in a finite basis
   *  set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class TwoPInts : public ParticleIntegrals {

  protected:

    // second basis numbers
    size_t sNB;

    TPI_KERNEL kernel_ = TPI_KERNEL::Coulomb; ///< Interaction kernel (Coulomb or short-range erfc)
    double omega_ = 0.; ///< Range-separation parameter for short-range integrals

  public:

    // Constructors

    TwoPInts() = delete;
    TwoPInts( const TwoPInts & ) = default;
    TwoPInts( TwoPInts && ) = default;

    TwoPInts(size_t nb, size_t snb = 0,
             TPI_KERNEL kernel = TPI_KERNEL::Coulomb, double omega = 0.):
        ParticleIntegrals(nb), sNB(snb), kernel_(kernel), omega_(omega) {
        // if the second basis does not exist, set it to be the same as the first one
        if (snb == 0) 
          sNB = nb; 

        // Validate the range-separation parameter for attenuated kernels
        if (kernel_ == TPI_KERNEL::ShortRangeErfc and
            (not std::isfinite(omega_) or omega_ <= 0.))
          CErr("Short-range erfc integrals require finite omega > 0.");
        if (kernel_ == TPI_KERNEL::Coulomb) omega_ = 0.;
     }

    template <typename IntsU>
    TwoPInts( const TwoPInts<IntsU> &other, int = 0 ):
        TwoPInts(other.nBasis(), other.snBasis(),
                 other.kernel(), other.rangeSeparationParameter()) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real TwoPInts from a Complex one.");
    }

    // return the basis numbers of the second particle
    size_t snBasis() const { return sNB; }

    // Interaction kernel accessors
    TPI_KERNEL kernel() const { return kernel_; }
    double rangeSeparationParameter() const { return omega_; }

    // Single element interfaces
    virtual IntsT operator()(size_t, size_t, size_t, size_t) const = 0;
    virtual IntsT operator()(size_t, size_t) const = 0;

    /**
     *  \brief Construct a fresh integral object of the same settings
     *  but evaluated with a different interaction kernel
     *
     *  The base implementation errors out. 
     *  Supported derived classes (Direct, Incore, DynamicERI) implement this in their own classes.
     */
    virtual std::shared_ptr<TwoPInts<IntsT>> createWithKernel(TPI_KERNEL kernel, double omega) const {
      CErr("Short-range/attenuated kernels are not implemented for this "
           "integral representation; use ALG=DIRECT, ALG=INCORE, or RI=DYNAMICERI.");
      return nullptr;
    }

    void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        MPIBCast(sNB,root,comm);
        int kernel = static_cast<int>(kernel_);
        MPIBCast(kernel, root, comm);
        kernel_ = static_cast<TPI_KERNEL>(kernel);
        MPIBCast(omega_, root, comm);
      }
#endif
    }

    //virtual TensorContraction ERITensor();

    virtual ~TwoPInts() {}

  }; // class ERInts

  /**
   *  \brief Templated class to define the interface to perform
   *  transformations and contractions of ERInts. Handles the
   *  contraction of 2-body (3,4 index) integrals with
   *  1-body (2 index) operators.
   *
   *  Templated over matrix type (MatsT) to allow for a seamless
   *  interface to both real- and complex-valued coefficients
   *  and density.
   */
  template <typename MatsT, typename IntsT>
  class TPIContractions {

    template <typename MatsU, typename IntsU>
    friend class TPIContractions;

  protected:
    std::shared_ptr<TwoPInts<IntsT>> ints_;
    bool oneCenterK_ = false;
    std::vector<size_t> mapCen2BfSt_;

  public:

    // Constructors

    TPIContractions() = delete;
    TPIContractions(std::shared_ptr<TwoPInts<IntsT>> tpi): ints_(tpi) {}
    template <typename MatsU>
    TPIContractions( const TPIContractions<MatsU,IntsT> &other, int dummy = 0 ):
      TPIContractions(other.ints_) {
      contractSecond = other.contractSecond;
      isCross = other.isCross;
    }
    template <typename MatsU>
    TPIContractions( TPIContractions<MatsU,IntsT> &&other, int dummy = 0 ):
      TPIContractions(other.ints_) {
      contractSecond = other.contractSecond;
      isCross = other.isCross;
    }

    TPIContractions( const TPIContractions &other ):
      TPIContractions(other, 0) {
      contractSecond = other.contractSecond;
      isCross = other.isCross;
    }
    TPIContractions( TPIContractions &&other ):
      TPIContractions(std::move(other), 0) {
      contractSecond = other.contractSecond;
      isCross = other.isCross;
    }

    std::shared_ptr<TwoPInts<IntsT>> ints() const { return ints_; }

    bool oneCenterK() const { return oneCenterK_; }
    const std::vector<size_t>& mapCen2BfSt() const { return mapCen2BfSt_; }
    void setOneCenterK(bool oneCenterK) { oneCenterK_ = oneCenterK; }
    void setMapCen2BfSt(const std::vector<size_t>& mapCen2BfSt) { mapCen2BfSt_ = mapCen2BfSt; }

    // Computation interfaces

    /**
     *  Contract the two body potential with one body (2 index) operators.
     *
     *  Smartly determines whether to do the contraction directly, incore
     *  or using density fitting depending on context
     *
     *  \param [in/ont] contList List of one body operators for contraction.
     */
    virtual void twoBodyContract(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        EMPerturbation&) const = 0;

    inline void twoBodyContract(
        MPI_Comm comm,
        const bool screen,
        std::vector<TwoBodyContraction<MatsT>> &contList) const {
      EMPerturbation pert;
      twoBodyContract(comm,screen,contList,pert);
    }

    inline void twoBodyContract(
        MPI_Comm comm,
        std::vector<TwoBodyContraction<MatsT>> &contList,
        EMPerturbation &pert) const {
      twoBodyContract(comm,true,contList,pert);
    }

    inline void twoBodyContract(
        MPI_Comm comm,
        std::vector<TwoBodyContraction<MatsT>> &contList) const {
      twoBodyContract(comm,true,contList);
    }

    void printTiming(const char* label, double elapsed, MPI_Comm comm, int rank, int size) const {
      #ifdef CQ_ENABLE_MPI
        double minT, maxT, avgT;
        MPI_Reduce(&elapsed, &minT, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&elapsed, &maxT, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        MPI_Reduce(&elapsed, &avgT, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avgT /= size;
        if (rank == 0) {
          std::cout << "        " << std::left << std::setw(38) << label
                    << "min = " << minT
                    << ", max = " << maxT
                    << ", avg = " << avgT << std::endl;
        }
      #else
        std::cout << "        " << std::left << std::setw(38) << label << elapsed << " s " << std::endl;
      #endif
    }

    // Destructor
    virtual ~TPIContractions() {}

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<TPIContractions<MatsU,IntsT>>
    convert(const std::shared_ptr<TPIContractions<MatsT,IntsT>>&);

    // Whether the contraction is done in the first or second basis
    bool contractSecond = false;
    // Whether the interaction integral is a cross inter-particle integral
    bool isCross = false;
    // Whether to time contractions when building Fock matrices and print to output
    bool printContractionTiming = false;

  }; // class TPIContractions
  
  template <typename IntsT>
  class InCoreTPI : public TwoPInts<IntsT> {

    template <typename IntsU>
    friend class InCoreTPI;

  public:

    // Constructor
    InCoreTPI() = delete;
    InCoreTPI(size_t nb, size_t snb = 0,
              TPI_KERNEL kernel = TPI_KERNEL::Coulomb, double omega = 0.):
      TwoPInts<IntsT>(nb, snb, kernel, omega) {}
    InCoreTPI( const InCoreTPI &other ) = default;
    InCoreTPI( InCoreTPI &&other ) = default;
    template <typename IntsU>
    InCoreTPI( const InCoreTPI<IntsU> &other, int = 0 ):
    InCoreTPI(other.nBasis(), other.snBasis(),
              other.kernel(), other.rangeSeparationParameter()) {}

    InCoreTPI& operator=( const InCoreTPI &other ) = default;
    InCoreTPI& operator=( InCoreTPI &&other ) = default;

    // Tensor direct access
    virtual IntsT* pointer() = 0;
    virtual const IntsT* pointer() const = 0;

    // Computation interfaces
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) = 0;


    virtual void computeAOInts(BasisSet &basisSet, BasisSet &basisSet2, 
                               Molecule &mol, EMPerturbation &emPert, OPERATOR op, 
                               const HamiltonianOptions &hamiltonianOptions) = 0;

    virtual void clear() = 0;

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const = 0;

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
      TwoPInts<IntsT>::broadcast(comm, root);
    }

    virtual ~InCoreTPI() {}

  }; // class InCoreTPI
  

  template <typename MatsT, typename IntsT>
  class InCoreTPIContraction : public TPIContractions<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCoreTPIContraction;

  public:

    // Constructors

    InCoreTPIContraction() = delete;
    InCoreTPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      TPIContractions<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCoreTPIContraction(
        const InCoreTPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
        InCoreTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }
    template <typename MatsU>
    InCoreTPIContraction(
        InCoreTPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
        InCoreTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }

    InCoreTPIContraction( const InCoreTPIContraction &other ):
    InCoreTPIContraction(other, 0) {
      this->contractSecond = other.contractSecond;
    }
    InCoreTPIContraction( InCoreTPIContraction &&other ):
    InCoreTPIContraction(std::move(other), 0) {
      this->contractSecond = other.contractSecond;
    }

    // Computation interfaces
    virtual void twoBodyContract(
        MPI_Comm comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        EMPerturbation&) const;

    // Computation interfaces
    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const = 0;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const = 0;

    virtual ~InCoreTPIContraction() {}

  }; // class InCoreTPIContraction

}; // namespace ChronusQ
