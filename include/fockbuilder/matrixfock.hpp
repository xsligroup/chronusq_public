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

#include <fockbuilder.hpp>


namespace ChronusQ {

  /**
   * \brief The FourCompFock class
   */
  template <typename MatsT, typename IntsT>
  class MatrixFock : public FockBuilder<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class MatrixFock;

  protected:

    cqmatrix::PauliSpinorMatrices<MatsT> fockMatrix;

  public:

    // Constructors
    MatrixFock() = delete;
    MatrixFock(HamiltonianOptions hamiltonianOptions,
               const cqmatrix::PauliSpinorMatrices<MatsT> &matrix):
        FockBuilder<MatsT,IntsT>(hamiltonianOptions),
        fockMatrix(matrix) {}
    MatrixFock(HamiltonianOptions hamiltonianOptions,
               cqmatrix::PauliSpinorMatrices<MatsT> &&matrix):
        FockBuilder<MatsT,IntsT>(hamiltonianOptions),
        fockMatrix(matrix) {}

    // Different type
    template <typename MatsU>
    MatrixFock(const MatrixFock<MatsU,IntsT> &other):
    FockBuilder<MatsT,IntsT>(other), fockMatrix(other.fockMatrix){}
    template <typename MatsU>
    MatrixFock(MatrixFock<MatsU,IntsT> &&other):
    FockBuilder<MatsT,IntsT>(other), fockMatrix(std::move(other.fockMatrix)){}

    // Virtual destructor
    virtual ~MatrixFock() {}

    // Form a fock matrix (see include/fockbuilder/impl.hpp for docs)
    virtual void formFock(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1.);

  };

}