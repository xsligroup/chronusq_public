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

#include <modifyorbitals/conventionalSCF/extrap.hpp>

namespace ChronusQ {

/**
 *  \brief Obtain a new set of orbitals given a Fock matrix.
 *
 *  Currently implements the fixed-point SCF procedure.
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::getNewOrbitals(EMPerturbation& pert, VecMORef<MatsT>& mo,
                                            VecEPtr& eps) {

  // Form the Fock matrix D(k) -> F(k)
  ProgramTimer::timeOp("Form Fock", [&]() { this->modOrbOpt.formFock(pert); });

  // Transform AO fock into the orthonormal basis (on root MPI process)
  ao2orthoFock();

  // Compute the orbital gradient for convergence tests and DIIS
  computeOrbGradient(orbGrad);

  // Modify fock matrix if requested (on root MPI process)
  if( this->scfControls.doExtrap ) modifyFock(pert);

  // Diagonalize the orthonormal fock Matrix (on root MPI process)
  diagOrthoFock(mo, eps);

  ortho2aoMOs(mo);

  this->modOrbOpt.formDensity();

};   // NewtonRaphsonSCF<MatsT,IntsT>::getNewOrbitals

template<typename MatsT>
void ConventionalSCF<MatsT>::ao2orthoFock(VecShrdPtrMat<MatsT> fock) {

  ROOT_ONLY(this->comm);

  if( fock.empty() ) fock = this->modOrbOpt.getFock();

  VecShrdPtrOrtho<MatsT> ortho = this->modOrbOpt.getOrtho();
  for( size_t i = 0; i < fock.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) {
      fockMatrixOrtho[i] = ortho[i]->nonortho2ortho(*fock[i]);
    } else {
      fockMatrixOrtho[i] = *(fock[i]);
    }
  }
};

template<typename MatsT>
void ConventionalSCF<MatsT>::ao2orthoDen(VecShrdPtrMat<MatsT> den) {

  ROOT_ONLY(this->comm);

  if( den.empty() ) den = this->modOrbOpt.getOnePDM();

  VecShrdPtrOrtho<MatsT> ortho = this->modOrbOpt.getOrtho();
  for( size_t i = 0; i < den.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) {
      onePDMOrtho[i] = ortho[i]->ortho2nonortho(*den[i]);
    } else {
      onePDMOrtho[i] = *(den[i]);
    }
  }
};

template<typename MatsT>
void ConventionalSCF<MatsT>::diagOrthoFock(VecMORef<MatsT>& mo, VecEPtr& eps) {

  ROOT_ONLY(this->comm);

  for( size_t i = 0; i < mo.size(); i++ ) {
    size_t NB = mo[i].get().dimension();
    if( NB != fockMatrixOrtho[i].dimension() ) CErr("Ortho Fock and MO dimensions do not match");
    std::copy_n(fockMatrixOrtho[i].pointer(),NB*NB,mo[i].get().pointer());
    int INFO  = HermetianEigen('V', 'L', NB, mo[i].get().pointer(), NB, eps[i], this->memManager);
    if( INFO != 0 ) {
      std::cout << "Attempted to diagonalize " << i << "th Fock Matrix" << std::endl;
      CErr("HermetianEigen failed in Fock", std::cout);
    }
  }
}

template<typename MatsT>
void ConventionalSCF<MatsT>::ortho2aoMOs(VecMORef<MatsT>& mo) {

  ROOT_ONLY(this->comm);

  VecShrdPtrOrtho<MatsT> ortho = this->modOrbOpt.getOrtho();
  for( size_t i = 0; i < mo.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) { ortho[i]->ortho2nonorthoCoeffs(mo[i]); }
  }
}

template<typename MatsT>
void ConventionalSCF<MatsT>::printRunHeader(std::ostream& out, EMPerturbation& pert) const {
  OptimizeOrbitals<MatsT>::printRunHeader(out, pert);

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
