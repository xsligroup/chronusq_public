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
//#define REAL_SPACE_X2C_ALGORITHM
//#define UDU_ATOMIC_X2C_ALGORITHM

#include <molecule.hpp>
#include <basisset.hpp>
#include <fields.hpp>
#include <integrals.hpp>
#include <singleslater.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  class AtomicX2C;

  /**
   *  \brief The X2C class. A class to compute X2C Core Hamiltonian.
   *  Stores intermediate matrices.
   */
  template <typename MatsT, typename IntsT>
  class X2C {

    template <typename MatsU, typename IntsU>
    friend class X2C;
    template <typename MatsU, typename IntsU>
    friend class AtomicX2C;

  protected:

    Integrals<IntsT>   &aoints_;            ///< AOIntegrals for contracted basis
    SingleSlaterOptions ssOptions_;         ///< Options to build 4C singleslater object
    Molecule            molecule_;          ///< Molecule object for nuclear potential
    BasisSet            basisSet_;          ///< BasisSet for original basis defintion
    BasisSet            uncontractedBasis_; ///< BasisSet for uncontracted basis defintion
    Integrals<IntsT>    uncontractedInts_;  ///< AOIntegrals for uncontracted basis
    size_t              nPrimUse_;          ///< Number of primitives used in p space

  public:

    // Operator storage
    IntsT*  mapPrim2Cont = nullptr;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> W  = nullptr; ///< W = (\sigma p) V (\sigma p)

    // Transformations for momentum space
    IntsT*  UK = nullptr; ///< K transformation between p- and R-space
    double* p  = nullptr; ///< p momentum eigens

    // X and Y matrices, means differently in real- and momentum-spaces.
    std::shared_ptr<cqmatrix::Matrix<MatsT>> X = nullptr; ///< X = S * L^-1
    std::shared_ptr<cqmatrix::Matrix<MatsT>> Y = nullptr; ///< Y = 1/sqrt(1 + X**H * X)
    ///< In real-space non-orthogonal basis, Y = ( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2

    // Picture-change U matrics from primitives to contracted basis
    MatsT*  UL = nullptr; ///< Picture change matrix of large component
    MatsT*  US = nullptr; ///< Picture change matrix of small component


    // Constructors

    // Disable default constructor
    X2C() = delete;

    // Default copy and move constructors
    X2C(const X2C<MatsT,IntsT> &);
    X2C(X2C<MatsT,IntsT> &&);

    /**
     * \brief Constructor
     *
     *  \param [in] aoints             Reference to the global AOIntegrals
     *  \param [in] mol                Molecule object for molecular specification
     *  \param [in] basis              The GTO basis for integral evaluation
     *  \param [in] hamiltonianOptions Flags for AO integrals evaluation
     */
    X2C(Integrals<IntsT> &aoints,
        const Molecule &mol, const BasisSet &basis, SingleSlaterOptions ssOptions) :
      aoints_(aoints), ssOptions_(ssOptions),
      molecule_(mol), basisSet_(basis),
      uncontractedBasis_(basisSet_.uncontractBasis()) {}

    // Different type
    template <typename MatsU>
    X2C(const X2C<MatsU,IntsT> &other, int dummy = 0);
    template <typename MatsU>
    X2C(X2C<MatsU,IntsT> &&     other, int dummy = 0);

    /**
     *  Destructor.
     *
     *  Destructs a X2C object
     */
    virtual ~X2C() { dealloc(); }


    // Public Member functions

    // Deallocation (see include/x2c/impl.hpp for docs)
    virtual void dealloc();

    // Compute core Hamitlonian
    virtual void computeOneEX2C(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void computeOneEX2C_Umatrix();
    virtual void computeOneEX2C_UDU(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void computeOneEX2C_corr(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void saveX2C(std::shared_ptr<SingleSlaterBase>);
    void SNSOScale(std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>, SNSO_TYPE);
    void RowDepDCB_SNSO(std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);

    // Compute Fock X2C
    virtual void computeFockX2C(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>,
        std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>>,
        bool incore = true, double threshSchwarz = 1e-12);
    void computeFockX2C_Umatrix(const cqmatrix::Matrix<MatsT> &fourCompMOSpin);

    static void compute_CoreH_Fock(Molecule &mol,
        BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
        EMPerturbation &emPert,
        std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);

    // Compute the gradient
    virtual std::vector<double> getGrad(EMPerturbation&, SingleSlater<MatsT,IntsT>&) {
      CErr("X2C CoreH gradient NYI",std::cout);
      abort();
    }

  };

  void compute_X2C_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);

}; // namespace ChronusQ
