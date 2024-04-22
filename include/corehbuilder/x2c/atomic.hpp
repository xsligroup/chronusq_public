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

#include <corehbuilder/x2c.hpp>

namespace ChronusQ {

  /**
   *  \brief The AtomicX2C class. A class to compute X2C Core Hamiltonian.
   *  Stores intermediate matrices.
   */
  template <typename MatsT, typename IntsT>
  class AtomicX2C : public X2C<MatsT, IntsT> {

    template <typename MatsU, typename IntsU>
    friend class AtomicX2C;

  protected:
    ATOMIC_X2C_TYPE                type_;    ///< Type of atomic approximation
    std::vector<X2C<MatsT, IntsT>> atoms_;   ///< X2C of unique atoms in the molecule
    std::vector<size_t>            atomIdx_; ///< Index of unique atom in atoms_

  public:

    // Constructors

    // Disable default constructor
    AtomicX2C() = delete;

    // Default copy and move constructors
    AtomicX2C(const AtomicX2C<MatsT,IntsT> &);
    AtomicX2C(AtomicX2C<MatsT,IntsT> &&);

    /**
     * \brief Constructor
     *
     *  \param [in] aoints     Reference to the global AOIntegrals
     *  \param [in] mol        Molecule object for molecular specification
     *  \param [in] basis      The GTO basis for integral evaluation
     *  \param [in] type       The type of atomic X2C
     */
    AtomicX2C(Integrals<IntsT> &aoints,
        const Molecule &mol, const BasisSet &basis,
        SingleSlaterOptions ssOptions) :
      X2C<MatsT,IntsT>(aoints, mol, basis, ssOptions),
      type_(ssOptions.hamiltonianOptions.AtomicX2CType) {}

    // Different type
    template <typename MatsU>
    AtomicX2C(const AtomicX2C<MatsU,IntsT> &other, int dummy = 0);
    template <typename MatsU>
    AtomicX2C(AtomicX2C<MatsU,IntsT> &&     other, int dummy = 0);


    // Public Member functions

    // Deallocation (see include/x2c/impl.hpp for docs)
    virtual void dealloc();

    // Compute core Hamitlonian
    virtual void computeOneEX2C(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void computeOneEX2C_Umatrix();

    // Compute the gradient
    virtual void getGrad() {
      CErr("X2C CoreH gradient NYI",std::cout);
    }

  };

}; // namespace ChronusQ
