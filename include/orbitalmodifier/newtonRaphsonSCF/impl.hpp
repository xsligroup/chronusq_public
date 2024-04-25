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

#include <orbitalmodifier/newtonRaphsonSCF/quasi-newton.hpp>
#include <orbitalmodifier/newtonRaphsonSCF/rotate.hpp>

namespace ChronusQ {

/**
 *  \Brief: Computes the fock Matrix and then computes a new set of orbitals
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::getNewOrbitals(EMPerturbation& pert, vecMORef<MatsT>& mo,
                                             vecEPtr& eps) {

  // Form the Fock matrix D(k) -> F(k)
  ProgramTimer::timeOp("Form Fock", [&]() { this->orbitalModifierDrivers.formFock(pert); });

  if( not storedRef ) saveRefMOs(mo);

  computeGradient(mo);
  computeDiagHess(eps);

  // Compute NR step and update orbRot
  NewtonRaphsonIteration();

  // Rotate reference MO's to new MO's
  rotateMOs(mo);

  // Compute new Eigenvalues
  this->computeEigenvalues(pert, mo, eps);

  this->orbitalModifierDrivers.formDensity();
}

/*
 *   Brief: Function to choose which algorithm/approximation is used
 *          for the Newton step.
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::NewtonRaphsonIteration() {
  if( this->scfControls.nrAlg == FULL_NR ) {

    fullNRStep();

  } else if( this->scfControls.nrAlg == QUASI_BFGS ) {

    qnBFGSStep();

  } else if( this->scfControls.nrAlg == QUASI_SR1 ) {

    qnSR1Step();

  } else if( this->scfControls.nrAlg == GRAD_DESCENT ) {

    gradDescentStep();

  } else {

    CErr("Requested Newton-Raphson Step not yet implemented");

  }
}

template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::gradDescentStep(){

    MatsT* dx = CQMemManager::get().malloc<MatsT>(nParam);
    for( size_t i=0; i<nParam; ++i )
        dx[i] = orbGrad[i] / orbDiagHess[i];
    takeStep(dx);
    CQMemManager::get().free(dx);

}


/*
 *  Brief: Computes the gradient convergence criteria for optimization
 *         In this case it is the maximum element of the orbital gradient
 */
template<typename MatsT>
double NewtonRaphsonSCF<MatsT>::computeFDCConv() {
  // Compute the Maximum of the gradient or norm of gradient
  double maxGrad = 0.;
  for( size_t i = 0; i < nParam; ++i )
    if( maxGrad < std::abs(orbGrad[i]) ) maxGrad = std::abs(orbGrad[i]);
  return maxGrad;
};// NewtonRaphsonSCF<MatsT> :: computeFDCConv

/*
 *  Brief: Print Algorithm specific information at the end of the run header
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::printRunHeader(std::ostream& out, EMPerturbation& pert) const {
  OrbitalOptimizer<MatsT>::printRunHeader(out, pert);

  // Print Damping and Level Shift Info
  out << std::setw(38) << std::left << "  Static Trust Region:" << this->scfControls.nrTrust << std::endl;
  if( this->scfControls.doDamp ) {
    out << std::setw(38) << std::left << "  Static Damping Factor:" << this->scfControls.dampParam << std::endl;
    out << std::setw(38) << std::left << "  Damping Error:" << this->scfControls.dampError << std::endl;
  }
  if( this->scfControls.nrLevelShift != 0. ){
    out << std::setw(38) << std::left << "  Level Shift:" << this->scfControls.nrLevelShift << std::endl;
  }

  // Print QN Algorithm info
  if( this->scfControls.nrAlg == GRAD_DESCENT ) {
    out << std::setw(38) << std::left << "  Newton-Raphson Step Approximation: Gradient-Descent" << std::endl;
  } else if( this->scfControls.nrAlg == QUASI_BFGS ) {
    out << std::setw(38) << std::left << "  Newton-Raphson Step Approximation: Quasi-Newton Broyden-Fletcher-Goldfarb-Shanno" << std::endl;
    out << std::left << "      Saving " << this->scfControls.nKeep << " previous iterations for Quasi-Newton" << std::endl;
  } else if( this->scfControls.nrAlg == QUASI_SR1 ) {
    out << std::setw(38) << std::left << "  Newton-Raphson Step Approximation: Quasi-Newton Symmetric Rank One" << std::endl;
    out << std::left << "      Saving " << this->scfControls.nKeep << " previous iterations for Quasi-Newton" << std::endl;
  } else if( this->scfControls.nrAlg == FULL_NR ) {
    out << std::setw(38) << std::left << "  Computing Full Hessian Matrix" << std::endl;
  }
};

/*
 *  Brief: Perform sanity checks 
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::sanityChecks(){

  // Check self-rotations
  for( auto& rO : rotOpt )
    for( auto& iO : rO.rotIndices.first )
      for( auto& iV : rO.rotIndices.second )
        if( iO == iV ) CErr("Cannot rotate an orbital into itself");
  
  // Check that there are no repeated rotations
  std::unordered_map<std::string,bool> checkMap;
  for( auto& rO : rotOpt )
    for( auto& iO : rO.rotIndices.first )
      for( auto& iV : rO.rotIndices.second ){
        std::string indStr = std::to_string(rO.spaceIndex) + "_" + std::to_string(iO)
            + "_" + std::to_string(iV);
        if( checkMap.find(indStr) != checkMap.end() ){
          CErr("Rotation parameters were included twice");
        } else {
          checkMap.insert(std::pair<std::string,bool>(indStr,true));
        }
      }

  // Check that the spaceIndex is correct
  vecShrdPtrMat<MatsT> den = this->orbitalModifierDrivers.getOnePDM();
  size_t nMat = den.size();
  for( auto& rO : rotOpt )
    if( rO.spaceIndex >= nMat ) CErr("Space Index in NRRotOptions is larger than the input vectors");
}


std::pair<std::set<size_t>,std::set<size_t>> ssNRRotIndices( size_t nOcc, size_t NB, size_t shift ){

  std::pair<std::set<size_t>, std::set<size_t>> ssRot;
  for( size_t i=0; i<nOcc; ++i )
    ssRot.first.insert(i+shift);
  for( size_t i=nOcc; i<NB; ++i )
    ssRot.second.insert(i+shift);
  return ssRot;
};

/*
 *  Brief: Allocates the parameter vectors
 */
template<typename MatsT>
void NewtonRaphsonSCF<MatsT>::alloc() {
  orbRot  = CQMemManager::get().malloc<MatsT>(nParam);
  orbGrad = CQMemManager::get().malloc<MatsT>(nParam);
  orbDiagHess = CQMemManager::get().malloc<MatsT>(nParam);
  std::fill_n(orbRot, nParam, MatsT(0.));
  std::fill_n(orbGrad, nParam, MatsT(0.));
  std::fill_n(orbDiagHess, nParam, MatsT(0.));
  if( this->scfControls.nrAlg == QUASI_BFGS or this->scfControls.nrAlg == QUASI_SR1 ) {
    for( size_t i = 0; i < this->scfControls.nKeep; i++ ) {
      qnOrbRot.emplace_back(CQMemManager::get().malloc<MatsT>(nParam));
      qnOrbGrad.emplace_back(CQMemManager::get().malloc<MatsT>(nParam));
      std::fill_n(qnOrbRot[i], nParam, MatsT(0.));
      std::fill_n(qnOrbGrad[i], nParam, MatsT(0.));
    }
  }
};

};   // namespace ChronusQ
// namespace ChronusQ
