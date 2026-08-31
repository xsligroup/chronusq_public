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
#include <basisset.hpp>
#include <matrix.hpp>
#include <hamiltonianoptions.hpp>
#include <cxxapi/input.hpp>

namespace ChronusQ {

  void read_option_and_remove_linear_dependency(
      const CQInputFile &input, BasisSet &basis, Molecule &mol, std::ostream &output);

  template <typename IntsT>
  void remove_linear_dependency(
      BasisSet &basis, Molecule &mol, const HamiltonianOptions &hamiltonianOptions,
      double linearDependencyThreshold, bool is4C = false, bool atomicOnly = false);

  template<typename IntsT>
  std::vector<size_t> two_step_remove_linearly_dependent_shells(
      BasisSet &originalBasis, Molecule &mol,
      std::shared_ptr<cqmatrix::Matrix<IntsT>> overlapMatrix,
      std::shared_ptr<cqmatrix::Matrix<IntsT>> kineticMatrix,
      double linearDependencyThreshold,
      bool is4C = false, bool atomicOnly = false);

  template<typename IntsT>
  void remove_linearly_dependent_shells(
      const BasisSet &originalBasis,
      cqmatrix::Matrix<IntsT>& overlapMatrix,
      std::vector<size_t>& keptShells,
      double linearDependencyThreshold);

}; // namespace ChronusQ