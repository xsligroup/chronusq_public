/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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
#include <realtime.hpp>
#include <geometrymodifier/moleculardynamics.hpp>

namespace ChronusQ {

  void MolecularDynamics::initializeMD(Molecule& molecule, std::shared_ptr<SingleSlaterBase> ss){

    if (mdOptions.restoreFromNuclearStep != 0) {
      restoreState(molecule, ss);
    } else {
      // TODO: set arbitrary velocity
      std::vector<double> velocityHalfTime(3*molecule.nAtoms, 0.);
      std::vector<double> velocityCurrent(3*molecule.nAtoms, 0.);
      std::vector<double> acceleration(3*molecule.nAtoms, 0.);
      nuclearKineticEnergy = 0.0;

      createMDDataSets(molecule, ss);
    }

  }


  /**
   *  \brief Updates the positions of the classical nuclei. Recompute
   *  all member data to fit with the new positions
   *
   *  \param[in] pos  An array holding the new positions
   */
  void MolecularDynamics::update(bool print,
              Molecule &molecule,
              bool firstStep, TDSCFOptions& tdSCFOptions, std::shared_ptr<SingleSlaterBase> ss,
              EMPerturbation& emPert, std::vector<std::shared_ptr<CubeGen>> cubes)
  {

    // Update step and time
    if ( firstStep ) {
      curState.iStep = 0;
      initializeMD(molecule, ss);
      if(mdOptions.pertFirstAtom) pertFirstAtom();
      electronicPotentialEnergy = finalMidpointFock();
    } else {
      curState.iStep++;
    }

    // Determine if we re-calculate gradient at this step
    bool doGrad = mdOptions.nMidpointFockSteps == 0 || 
                  curState.iStep % mdOptions.nMidpointFockSteps == 0;

    // Update gradient whenever we restart midpoint fock
    if( doGrad ) {

      // If we have midpoint fock steps, we need to take the final fock step to make sure geometry is at full \Delta t_N step
      if( mdOptions.nMidpointFockSteps != 0 && !firstStep ) {
        // Move geometry by half \Delta t_{N_q} step
        double half_fock_dt = (mdOptions.timeStepAU/mdOptions.nMidpointFockSteps) / 2 ;
        geometryVV(molecule, gradientCurrent, half_fock_dt);
        double totalTimeCur = curState.time;
        curState.time += half_fock_dt;
        molecule.update();
        electronicPotentialEnergy = finalMidpointFock();
        std::cout << "  *** Moving Nuclei from x( t = " << totalTimeCur << " au) to x( t = " << curState.time << " au) ***"<< std::endl;
      }

      gradientCurrent = gradientGetter();
      std::cout << "  *** Calculating Gradient at g( t = " << curState.time << " au) ***" << std::endl;

    }

    bool moveGeometry = true;
    bool moveVelocity = doGrad;

    if(print and doGrad) printCurrentGeometry(molecule);


    // if velocity Verlet
    if(moveVelocity) {
      // compute p( t + 1/2 \Delta t_N) 
      velocityVV(molecule, gradientCurrent, mdOptions.timeStepAU, firstStep);
      std::cout << "  *** Calculating Velocity at p( t = " << curState.time + 0.5*mdOptions.timeStepAU << " au ) ***" << std::endl;
      // compute kinetic energy
      computeKineticEnergy(molecule);
    }
    
    // Compute total energies
    currentTotalEnergy = electronicPotentialEnergy + nuclearKineticEnergy;
    if(firstStep and mdOptions.restoreFromNuclearStep == 0) totalEnergy0 = currentTotalEnergy;
    if(firstStep) previousTotalEnergy = currentTotalEnergy;



    // At this point we have full-step geom, time, g, v, and half-step v. Save these to bin file
    if(doGrad or mdOptions.saveAllGeometry) saveState(molecule, ss);

    if (std::any_of(cubes.begin(), cubes.end(), [](const std::shared_ptr<CubeGen>& ptr) { return ptr != nullptr; })) 
      ss->runCube(cubes, emPert, "_MDStep"+std::to_string(curState.iStep), std::make_shared<Molecule>(molecule));



    // Determine the stepsize at which we move the geometry:
    // If JobType is BOMD, geometry move in full \Delta t_N step
    // If JobType is Ehrenfest, then geometry move in \Delta t_{N_q} step ( \Delta t_{N_q} = \dfrac{ \Delta t_N }  { m } )
    curState.currentStepSize = mdOptions.nMidpointFockSteps == 0 ? 
        mdOptions.timeStepAU : mdOptions.timeStepAU/ mdOptions.nMidpointFockSteps;
    // For Ehrenfest job, if this step evaluates gradient, that means we are at full \Delta t_N step 
    // If so, we need take half \Delta t_{N_q} to start mid-point fock algorithm
    if(mdOptions.nMidpointFockSteps != 0 and doGrad) curState.currentStepSize /= 2;


    if(moveGeometry) {
      // compute x( t + \Delta t_Nq)
      geometryVV(molecule, gradientCurrent, curState.currentStepSize); 
      // update molecular properties after changing the geometry
      molecule.update();
      double totalTimeCur = curState.time;
      curState.time += curState.currentStepSize;
      std::cout << "  *** Moving Nuclei from x( t = " << totalTimeCur << " au) to x( t = " << curState.time << " au) ***"<< std::endl;
    }

    // Update proton velocity
    double timestep = curState.time - (curState.iStep / mdOptions.nMidpointFockSteps) * mdOptions.timeStepAU;
    updateProtonVelocity(molecule, gradientCurrent, timestep);

    



    // output important dynamic information
    if(print and doGrad) printMDInfo(molecule, currentTotalEnergy);

    previousTotalEnergy = currentTotalEnergy;

    // Set next RT simulation length by increasing maxSteps
    if(firstStep and mdOptions.restoreFromNuclearStep != 0) {
      tdSCFOptions.restoreFromStep = (curState.lastSavePoint-1);
      if (not mdOptions.saveAllGeometry)
        tdSCFOptions.restoreFromStep *= mdOptions.nMidpointFockSteps;
    }
  
    tdSCFOptions.maxSteps = (curState.iStep+1) * mdOptions.nElectronicSteps;
    //tdSCFOptions.maxSteps += mdOptions.nElectronicSteps;
    tdSCFOptions.tMax      = tdSCFOptions.maxSteps * tdSCFOptions.deltaT;

  }


  /**
  *  \brief Compute the nuclear-nuclear repulsion energy for classical
  *  point nuclei using the Atoms contained in the atoms array
  *
  *  \f[
  *    V_{NN} = \sum_{A < B} \frac{Z_A Z_B}{R_{AB}}
  *  \f]
  */ 
  void MolecularDynamics::computeKineticEnergy(Molecule& molecule) {

    nuclearKineticEnergy = 0.;

    size_t i = -3;
    for( Atom& atom : molecule.atoms ) {
      
      i += 3;

      if(atom.quantum && !NEODynamicsOpts.includeQProtKE) continue;

      nuclearKineticEnergy += 0.5*velocityCurrent[i  ]*velocityCurrent[i  ]*atom.atomicMass*AUPerAMU;
      nuclearKineticEnergy += 0.5*velocityCurrent[i+1]*velocityCurrent[i+1]*atom.atomicMass*AUPerAMU;
      nuclearKineticEnergy += 0.5*velocityCurrent[i+2]*velocityCurrent[i+2]*atom.atomicMass*AUPerAMU;
      
    }

    molecule.nucKinEnergy = nuclearKineticEnergy;
  }

  bool MolecularDynamics::hasNext() {
    // Because this is at the beginning of the previous iteration, do i-1
    double nextStep = mdOptions.timeStepAU;
    if( mdOptions.nMidpointFockSteps != 0 )
      nextStep /= mdOptions.nMidpointFockSteps;
    return (curState.time + nextStep) <= tMax;
  }

}