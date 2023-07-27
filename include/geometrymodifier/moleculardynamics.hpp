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
#include <molecule.hpp>
#include <geometrymodifier.hpp>
#include <geometrymodifier/moleculardynamics/enums.hpp>

namespace ChronusQ {

  struct NEODynamicsOptions {
    QuantumProtonMoveAlg QProtMoveAlg = FIXED; // Default NEO-Ehrenfest alg is to fix quantum proton basis center
    bool includeQProtKE = false;               // Whether to include the translational KE associated with changes in protonic basis centers         
  }; // struct NEODynamicsOptions

  /**
   * \brief The MolecularDynamics class
   */
  class MolecularDynamics : public GeometryModifier {

    IntegrationProgress curState;
    double tMax;

  public:

    std::function<std::vector<double>()> gradientGetter;
    std::function<double(double)> finalMidpointFock;
    std::vector<double> gradientCurrent;
    std::vector<double> velocityHalfTime;  ///< nuclear velocity at half time, (t-1/2) upon entry and (t+1/2) upon exist
    std::vector<double> velocityCurrent;   ///< nuclear velocity at the current time (t)
    std::vector<double> acceleration;      ///< acceleration at the current time (t)

    double   nuclearKineticEnergy;      ///< nuclear kinetic energy
    double   totalEnergy0;              ///< total energy at step 0
    double   previousTotalEnergy;       ///< total energy at the previous step

    NEODynamicsOptions NEODynamicsOpts;


    // Constructors
    MolecularDynamics() = delete;
    MolecularDynamics(MolecularOptions molecularOptions, Molecule& molecule) :
      GeometryModifier(molecularOptions),
      gradientCurrent(3*molecule.nAtoms, 0.),
      velocityHalfTime(3*molecule.nAtoms, 0.),
      velocityCurrent(3*molecule.nAtoms, 0.),
      acceleration(3*molecule.nAtoms, 0.)
    {
      curState.stepSize = molecularOptions.timeStepAU;
      tMax = molecularOptions.timeStepAU * molecularOptions.nNuclearSteps;
    };

    // Different type
    MolecularDynamics(const MolecularDynamics &other):
        GeometryModifier(other){}
    MolecularDynamics(MolecularDynamics &&other):
        GeometryModifier(other){}

    // Virtual destructor
    virtual ~MolecularDynamics() {}

    virtual bool hasNext() override;

    virtual void update(bool print, Molecule &molecule, bool firstStep) override;

    void initializeMD(Molecule& molecule);

    void velocityVV(Molecule &molecule, std::vector<double> gradientCurrent, double timeStep, bool firstStep);

    void geometryVV(Molecule &molecule, std::vector<double> gradientCurrent, double timeStep);

    void computeKineticEnergy(Molecule& molecule);
  };
}
