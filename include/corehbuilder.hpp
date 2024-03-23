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

#include <fields.hpp>
#include <singleslater.hpp>
#include <integrals.hpp>

namespace ChronusQ {

  /**
   * \brief The CoreHBuilder class
   */
  template <typename MatsT, typename IntsT>
  class CoreHBuilder {

    template <typename MatsU, typename IntsU>
    friend class CoreHBuilder;

  protected:
    Integrals<IntsT>  &aoints_;
    HamiltonianOptions hamiltonianOptions_; ///< One electron terms to be computed

  public:

    // Constructors

    // Disable default constructor
    CoreHBuilder() = delete;
    CoreHBuilder(Integrals<IntsT> &aoints, HamiltonianOptions hamiltonianOptions):
      aoints_(aoints), hamiltonianOptions_(hamiltonianOptions) {}

    // Same or Different type
    template <typename MatsU>
    CoreHBuilder(const CoreHBuilder<MatsU,IntsT> &other):
      aoints_(other.aoints_), hamiltonianOptions_(other.hamiltonianOptions_) {}
    template <typename MatsU>
    CoreHBuilder(CoreHBuilder<MatsU,IntsT> &&other):
      aoints_(other.aoints_), hamiltonianOptions_(other.hamiltonianOptions_) {}

    // Virtual destructor
    virtual ~CoreHBuilder() {}


    // Public member functions
    const HamiltonianOptions& getHamiltonianOptions() const {
      return hamiltonianOptions_;
    }

    // Compute various core Hamitlonian
    virtual void computeCoreH(EMPerturbation&,
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>) = 0;

    // Compute the gradient
    virtual std::vector<double> getGrad(EMPerturbation&,
      SingleSlater<MatsT,IntsT>&) = 0;

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<CoreHBuilder<MatsU,IntsT>>
    convert(const std::shared_ptr<CoreHBuilder<MatsT,IntsT>>&);

  };

}
