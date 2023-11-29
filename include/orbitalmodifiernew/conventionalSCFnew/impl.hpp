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

#include <orbitalmodifiernew/conventionalSCFnew/extrap.hpp>

namespace ChronusQ {

/**
 *  \brief Obtain a new set of orbitals given a Fock matrix.
 *
 *  Currently implements the fixed-point SCF procedure.
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void ConventionalSCFNew<singleSlaterT,MatsT,IntsT>::getNewOrbitals(EMPerturbation& pert) {

  // Transform AO fock into the orthonormal basis (on root MPI process)
  this->ao2orthoFock();

  // Modify fock matrix if requested (on root MPI process)
  if( this->scfControls.doExtrap ) modifyFock(pert);

  // Diagonalize the orthonormal fock Matrix (on root MPI process)
  this->diagOrthoFock();
  this->ortho2aoMOs();

};    //ConventionalSCF<MatsT>::getNewOrbitals


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
  void ConventionalSCFNew<singleSlaterT,MatsT,IntsT>::printRunHeader(std::ostream& out, EMPerturbation& pert) const {
  OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::printRunHeader(out, pert);

  // Print DIIS Algorithm info
  if( this->scfControls.doExtrap ) {
    if( this->scfControls.doDamp ) {
      out << std::setw(38) << std::left << "  Static Damping Factor:" << this->scfControls.dampParam << std::endl;
      out << std::setw(38) << std::left << "  Damping Error:" << this->scfControls.dampError << std::endl;
    }

    if( this->scfControls.diisAlg != NONE ) {
      out << std::setw(38) << std::left << "  DIIS Extrapolation Algorithm:";
      if( this->scfControls.diisAlg == CDIIS ) out << "CDIIS";
      if( this->scfControls.diisAlg == EDIIS ) out << "EDIIS";
      if( this->scfControls.diisAlg == CEDIIS ) out << "CEDIIS";
      out << std::endl;

      if( this->scfControls.diisAlg == CEDIIS ) {
        out << std::left << "    * CDIIS will track up to " << this->scfControls.nKeep << " previous iterations" << std::endl;
        out << std::left << "    * EDIIS will track up to " << this->scfControls.nKeep << " previous iterations" << std::endl;
        out << std::left << "    * will switch at " << std::fixed << std::setprecision(8) << this->scfControls.cediisSwitch << " for max([F,D])"
            << std::endl;
      } else {
        out << std::left << "    * DIIS will track up to " << this->scfControls.nKeep << " previous iterations" << std::endl;
      }
    }
  }
};

};   // namespace ChronusQ
