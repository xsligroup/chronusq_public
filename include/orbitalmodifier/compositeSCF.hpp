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

/*
 *    This is a proposed object for carrying out a composite SCF calculation
 *    (where you can use a combination of several methods). However, this
 *    object has not been implemented or tested in anyway. It is included
 *    here for someone in the future to use it if they want to. 
 *
 */
#include <orbitalmodifier.hpp>

namespace ChronusQ {

template<typename MatsT>
class CompositeSCF : public OrbitalModifier<MatsT> {
protected:
  std::vector<std::shared_ptr<OrbitalModifier<MatsT>>> orbitalModifier;
  std::vector<SCFControls> scfControls;

public:
  // Constructor
  CompositeSCF(std::vector<SCFControls> sC, MPI_Comm comm, OrbitalModifierDrivers<MatsT> modOpt):
    ModifyOrbitals(comm, modOpt) {

    for( auto s : scfControls ) {
      if( s.scfAlg == _CONVENTIONAL_SCF ) {
        orbitalModifier.emplace_back(std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(std::make_shared<ConventionalSCF<MatsT>>(s, comm, modOpt)));
      } else if( s.scfAlg == _NEWTON_RAPHSON_SCF ) {
        orbitalModifier.emplace_back(
          std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(std::make_shared<NewtonRaphsonSCF<MatsT>>(s, comm, modOpt)));
      } else {
        orbitalModifier.emplace_back(std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(std::make_shared<SkipSCF<MatsT>>(s, comm, modOpt)));
      }
    }
  };

  // Member Functions

  void runModifyOrbitals(EMPerturbation& pert, vecMORef& mo, vecEPtr eps) {
    for( auto mO : orbitalModifier )
      mO->runModifyOrbitals(pert, mo, eps);
  };

  void getNewOrbitals(EMPerturbation& pert, vecMORef& mo, vecEPtr eps, bool frmFock = true);

  // Printing functions
  virtual void printRunHeader(std::ostream& out, EMPerturbation&);
  virtual void printIteration(std::ostream& out = std::cout, bool printDiff = true) {};
}

template<typename MatsT>
void CompositeSCF<MatsT>::getNewOrbitals(EMPerturbation& pert, vecMORef& mo, vecEPtr eps, frmFock) {
  for( auto mO : modifyOrbitals )
    mO->runModifyOrbitals(pert, mo, eps, frmFock);
};

template<typename MatsT>
void CompositeSCF<MatsT> :: printRunHeader( std::ostream& out, EMPerturbation& pert){
    std::cout << "Performing a Composite SCF calculation with the following Algorithms:" << std::endl;
    // FIXME: complete this function
};
};   // NameSpace ChronusQ
