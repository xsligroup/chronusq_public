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
   * \brief The MatrixCoreH class
   */
  template <typename MatsT, typename IntsT>
  class MatrixCoreH : public CoreHBuilder<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class MatrixCoreH;

  protected:

    cqmatrix::PauliSpinorMatrices<MatsT> coreHMatrix;

  public:

    // Constructors

    // Disable default constructor
    MatrixCoreH() = delete;
    MatrixCoreH(Integrals<IntsT> &aoints, HamiltonianOptions hamiltonianOptions,
                const cqmatrix::PauliSpinorMatrices<MatsT> &matrix):
        CoreHBuilder<MatsT,IntsT>(aoints, hamiltonianOptions),
        coreHMatrix(matrix) {}
    MatrixCoreH(Integrals<IntsT> &aoints, HamiltonianOptions hamiltonianOptions,
                cqmatrix::PauliSpinorMatrices<MatsT> &&matrix):
        CoreHBuilder<MatsT,IntsT>(aoints, hamiltonianOptions),
        coreHMatrix(matrix) {}

    // Same or Different type
    template <typename MatsU>
    MatrixCoreH(const MatrixCoreH<MatsU,IntsT> &other):
    CoreHBuilder<MatsT,IntsT>(other), coreHMatrix(other.coreHMatrix) {}
    template <typename MatsU>
    MatrixCoreH(MatrixCoreH<MatsU,IntsT> &&other):
    CoreHBuilder<MatsT,IntsT>(other), coreHMatrix(std::move(other.coreHMatrix)) {}

    // Virtual destructor
    virtual ~MatrixCoreH() {}

    // Public member functions

    // Compute core Hamitlonian
    virtual void computeCoreH(EMPerturbation&,
                              std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>);

    // Compute the gradient
    virtual std::vector<double> getGrad(EMPerturbation&,
      SingleSlater<MatsT,IntsT>&) {
      CErr("Matrix CoreH gradient NYI",std::cout);
      abort();
    }

  };

}
