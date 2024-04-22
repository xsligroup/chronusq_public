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

#include <cqlinalg/blas3.hpp>
#include <singleslater/kohnsham.hpp>

#include <response.hpp>

namespace ChronusQ {

template<typename MatsT, typename IntsT>
void KohnSham<MatsT, IntsT>::buildOrbitalModifierOptions() {
  // Modify SCFControls
  this->scfControls.printLevel    = this->printLevel;
  this->scfControls.refLongName_  = this->refLongName_;
  this->scfControls.refShortName_ = this->refShortName_;

  // Initialize ModifyOrbitalOptions
  OrbitalModifierDrivers<MatsT> modOrbOpt;

  // Bind Lambdas to std::functions
  modOrbOpt.printProperties   = [this]() { this->printProperties(); };
  modOrbOpt.saveCurrentState  = [this]() { this->saveCurrentState(); };
  modOrbOpt.formFock          = [this](EMPerturbation& pert) { this->formFock(pert); };
  modOrbOpt.computeProperties = [this](EMPerturbation& pert) { this->computeProperties(pert); };
  modOrbOpt.computeEnergy     = [this](EMPerturbation& pert) { this->computeEnergy(); };
  modOrbOpt.formDensity       = [this]() { this->formDensity(); };
  modOrbOpt.getFock           = [this]() { return this->getFock(); };
  modOrbOpt.getOnePDM         = [this]() { return this->getOnePDM(); };
  modOrbOpt.getOrtho          = [this]() { return this->getOrtho(); };
  modOrbOpt.setDenEqCoeff     = [this](bool val) { this->setDenEqCoeff(val); };
  modOrbOpt.getTotalEnergy    = [this]() { return this->getTotalEnergy(); };

  // Make OrbitalModifier based on scfControls
  if( this->scfControls.scfAlg == _CONVENTIONAL_SCF ) {
    this->orbitalModifier = std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(
        std::make_shared<ConventionalSCF<MatsT>>(this->scfControls, this->comm, modOrbOpt));
  } else if( this->scfControls.scfAlg == _NEWTON_RAPHSON_SCF ) {

    bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT,IntsT>>(this->fockBuilder) != nullptr);
    // The orbital gradient is wrong for ROHF. Also, it neglects the mixing
    // between the doubly occupied orbitals and the singly occupied orbitals
    if( iRO ) CErr("NewtonRaphson SCF is not yet implemented for Restricted Open Shell");

    // Generate NRRotationOptions
    std::vector<NRRotOptions> rotOpt = this->buildRotOpt();

    this->orbitalModifier = std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(
        std::make_shared<NewtonRaphsonSCF<MatsT>>(rotOpt, this->scfControls, this->comm, modOrbOpt));
  } else {
    //this->scfControls.doExtrap = false;
    //this->orbitalModifier       = std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(
    //    std::make_shared<SkipSCF<MatsT>>(this->scfControls, this->comm, modOrbOpt));
  }
};   // KohnSham<MatsT,IntsT> :: buildOrbitalModifierOptions

template<typename MatsT, typename IntsT>
void KohnSham<MatsT, IntsT>::computeFullNRStep(MatsT* orbRot) {

  CErr("Full Newton Raphson SCF is not yet implemented for Kohn-Sham DFT");
}

template<typename MatsT, typename IntsT>
std::pair<double, MatsT*> KohnSham<MatsT, IntsT>::getStab() {
  CErr();
  return {0., nullptr};
};
};   // namespace ChronusQ
