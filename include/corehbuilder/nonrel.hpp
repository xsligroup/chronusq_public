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
   * \brief The NRCoreH class
   */
  template <typename MatsT, typename IntsT>
  class NRCoreH : public CoreHBuilder<MatsT,IntsT> {

  public:

    // Constructors

    // Disable default constructor
    NRCoreH() = delete;
    NRCoreH(Integrals<IntsT> &aoints, HamiltonianOptions hamiltonianOptions):
      CoreHBuilder<MatsT,IntsT>(aoints, hamiltonianOptions) {}

    // Same or Different type
    template <typename MatsU>
    NRCoreH(const NRCoreH<MatsU,IntsT> &other):
      CoreHBuilder<MatsT,IntsT>(other) {}
    template <typename MatsU>
    NRCoreH(NRCoreH<MatsU,IntsT> &&other):
      CoreHBuilder<MatsT,IntsT>(other) {}

    // Virtual destructor
    virtual ~NRCoreH() {}

    // Public member functions

    // Compute core Hamitlonian
    void addMagPert(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void computeCoreH(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);
    virtual void computeNRCH(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);

    // Compute the gradient
    virtual std::vector<double> getGrad(EMPerturbation&,
      SingleSlater<MatsT,IntsT>&);

  };

}
