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
#include <particleintegrals.hpp>
#include <matrix.hpp>
#include <util/mpi.hpp>

namespace ChronusQ {

  enum TWOBODY_CONTRACTION_TYPE {
    COULOMB, ///< (mn | kl) X(lk)
    EXCHANGE,///< (mn | kl) X(nk)
    PAIR,    ///< (mn | kl) X(nl)
    BARE_COULOMB,
    LLLL,
    LLSS,
    SSSS,
    GAUNT,
    GAUGE
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

    double* ERI4 = nullptr;

    INTEGRAL_TRANSPOSE intTrans;


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

  public:

    // Constructors

    TwoPInts() = delete;
    TwoPInts( const TwoPInts & ) = default;
    TwoPInts( TwoPInts && ) = default;

    TwoPInts(size_t nb, size_t snb = 0):
        ParticleIntegrals(nb), sNB(snb) {
        
        // if the second basis does not exist, set it to be the same as the first one
        if (snb == 0) 
          sNB = nb; 
     }

    template <typename IntsU>
    TwoPInts( const TwoPInts<IntsU> &other, int = 0 ):
        TwoPInts(other.nBasis(), other.snBasis()) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real TwoPInts from a Complex one.");
    }

    // return the basis numbers of the second particle
    size_t snBasis() const { return sNB; }

    // Single element interfaces
    virtual IntsT operator()(size_t, size_t, size_t, size_t) const = 0;
    virtual IntsT operator()(size_t, size_t) const = 0;

    void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        MPIBCast(sNB,root,comm);
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

  public:

    // Constructors

    TPIContractions() = delete;
    TPIContractions(std::shared_ptr<TwoPInts<IntsT>> tpi): ints_(tpi) {}
    template <typename MatsU>
    TPIContractions( const TPIContractions<MatsU,IntsT> &other, int dummy = 0 ):
      TPIContractions(other.ints_) {
      contractSecond = other.contractSecond;
    }
    template <typename MatsU>
    TPIContractions( TPIContractions<MatsU,IntsT> &&other, int dummy = 0 ):
      TPIContractions(other.ints_) {
      contractSecond = other.contractSecond;
    }

    TPIContractions( const TPIContractions &other ):
      TPIContractions(other, 0) {
      contractSecond = other.contractSecond;
    }
    TPIContractions( TPIContractions &&other ):
      TPIContractions(std::move(other), 0) {
      contractSecond = other.contractSecond;
    }

    std::shared_ptr<TwoPInts<IntsT>> ints() const { return ints_; }

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

    // Destructor
    virtual ~TPIContractions() {}

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<TPIContractions<MatsU,IntsT>>
    convert(const std::shared_ptr<TPIContractions<MatsT,IntsT>>&);

    // Whether the contraction is done in the first or second basis
    bool contractSecond = false;
    
    // Whether to time contractions when building Fock matrices and print to output
    bool printContractionTiming = false;

  }; // class TPIContractions

}; // namespace ChronusQ
