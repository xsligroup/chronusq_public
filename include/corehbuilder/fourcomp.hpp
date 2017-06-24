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

#include <corehbuilder.hpp>

namespace ChronusQ {

  /**
   * \brief The FourComponent class
   */
  template <typename MatsT, typename IntsT>
  class FourComponent : public CoreHBuilder<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class FourComponent;

  public:

    // Constructors

    // Disable default constructor
    FourComponent() = delete;

    // Default copy and move constructors
    FourComponent(const FourComponent<MatsT,IntsT> &);
    FourComponent(FourComponent<MatsT,IntsT> &&);

    /**
     * \brief Constructor
     *
     *  \param [in] aoints             Reference to the global AOIntegrals
     *  \param [in] memManager         Memory manager for matrix allocation
     *  \param [in] mol                Molecule object for molecular specification
     *  \param [in] basis              The GTO basis for integral evaluation
     *  \param [in] hamiltonianOptions Flags for AO integrals evaluation
     */
    FourComponent(Integrals<IntsT> &aoints, HamiltonianOptions hamiltonianOptions) :
      CoreHBuilder<MatsT,IntsT>(aoints, hamiltonianOptions) {}

    // Different type
    template <typename MatsU>
    FourComponent(const FourComponent<MatsU,IntsT> &other, int dummy = 0);
    template <typename MatsU>
    FourComponent(FourComponent<MatsU,IntsT> &&     other, int dummy = 0);

    /**
     *  Destructor.
     *
     *  Destructs a FourComponent object
     */
    virtual ~FourComponent() {}


    // Public Member functions

    // Compute core Hamitlonian
    virtual void computeCoreH(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void compute4CCH(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);

    // Compute the gradient
    virtual std::vector<double> getGrad(EMPerturbation&, SingleSlater<MatsT,IntsT>&) {
      CErr("4C CoreH gradient NYI",std::cout);
      abort();
    }

  };

}
