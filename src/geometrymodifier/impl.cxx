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

    // Update step
    if ( firstStep ) {
      curState.iStep = 0;
      initializeMD(molecule, ss);
      if(mdOptions.pertFirstAtom) pertFirstAtom();
      if(mdOptions.projectOrthoDen or mdOptions.pertFirstAtom) electronicPotentialEnergy = finalMidpointFock();
    } else {
      curState.iStep++;
    }

    // Determine if we re-calculate gradient at this step
    bool isBOMD = mdOptions.nMidpointFockSteps == 0;
    bool isEhrenfest = mdOptions.nMidpointFockSteps != 0;
    bool doGrad = isBOMD || curState.iStep % mdOptions.nMidpointFockSteps == 0;
    bool isHalfNFockStep(false), OddNFockSteps(false);
    if ( isEhrenfest ) {
      // Special case for NMidpointFockSteps=2
      if (mdOptions.nMidpointFockSteps == 2) {
        isHalfNFockStep = curState.iStep%mdOptions.nMidpointFockSteps == 1;
        OddNFockSteps = false;
      } else {
        // Determine if we are at a half-point of all midpoint fock steps (i.e. t = Δt_N/2 + iΔt_N )
        isHalfNFockStep = curState.iStep%mdOptions.nMidpointFockSteps == (mdOptions.nMidpointFockSteps)/2+1;
        // Determine if the number of midpoint fock steps is odd
        OddNFockSteps = mdOptions.nMidpointFockSteps%2 == 1;
      }
    }

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    // =========================================================================================
    // Update gradient (if needed)
    // =========================================================================================
    if ( doGrad ) {

      // If we have midpoint fock steps, we need to take the final fock step before updating gradient
      // This will make sure geometry, velocity, and electronic density are all at full Δt_N step
      if ( isEhrenfest && !firstStep ) {
        
        // Update velocity (at half-step) to be full-step (where geometry is at)
        double dt = curState.currentStepSize; 
        std::cout << "  *** Updating Velocity from p( t = " << curState.ptime << " au) to p( t = " << curState.ptime + dt/2 << " au) ***"<< std::endl;
        //printCurrentVelocity(molecule);
        if ((magAmp[0] != 0.) || (magAmp[1] != 0.) || (magAmp[2] != 0.))
          velocityVV_EXPK1(molecule, velocity, velocity, gradient, dt, emPert);
        else
          velocityVV(molecule, velocity, velocity, gradient, dt);
        curState.ptime += dt/2;
        
        // Set geometry step to be half Δt_{N_q} step
        double half_fock_dt = (mdOptions.timeStepAU/mdOptions.nMidpointFockSteps) / 2 ;
        
        // Update velocity to be half-step
        std::cout << "  *** Updating Velocity from p( t = " << curState.ptime << " au) to p( t = " << curState.ptime + half_fock_dt/2 << " au) ***"<< std::endl;
        //printCurrentVelocity(molecule);
        if ((magAmp[0] != 0.) || (magAmp[1] != 0.) || (magAmp[2] != 0.))
          velocityVV_EXPK1(molecule, velocity, velocity, gradient, half_fock_dt, emPert);
        else
          velocityVV(molecule, velocity, velocity, gradient, half_fock_dt);
        curState.ptime += half_fock_dt/2;

        // Update geometry using velocity at half-step
        geometryVV(molecule, half_fock_dt);
        double totalTimeCur = curState.time;
        curState.time += half_fock_dt;
        molecule.update();
        electronicPotentialEnergy = finalMidpointFock();
        std::cout << "  *** Updating Geometry from x( t = " << totalTimeCur << " au) to x( t = " << curState.time << " au) ***"<< std::endl;
        
      }

      // Obtain new gradient 
      std::cout << "  *** Calculating Gradient at g( t = " << curState.time << " au) ***" << std::endl;
      gradient = gradientGetter();

      // If we have midpoint fock steps, we calculate velocity at full Δt_N step
      // p(t+Δt_N) = p(t+0.5Δt_N) + 0.5 * g(t+Δt_N) / m * Δt_N
      if ( isEhrenfest && !firstStep ) {
        double dt = mdOptions.timeStepAU; 
        std::cout << "  *** Updating Velocity from p( t = " << curState.ptimeHalf << " au) to p( t = " << curState.ptimeHalf + dt/2 << " au) ***"<< std::endl;
        //printCurrentVelocity(molecule);
        if ((magAmp[0] != 0.) || (magAmp[1] != 0.) || (magAmp[2] != 0.))
          velocityVV_EXPK1(molecule, velocityHalfTN, velocity, gradient, dt, emPert);
        else
          velocityVV(molecule, velocityHalfTN, velocity, gradient, dt);
        curState.ptime = curState.ptimeHalf + dt/2 ;
      }

    }



    // =========================================================================================
    // Update half-step velocity to full-step velocity
    // Here we save the velocity at t+0.5Δt_N step, which we later use to compute next full-step velocity
    // For even number of midpoint fock steps, p will arrive to t+0.5Δt_N before the velocity update
    // For odd number of midpoint fock steps,  p will arrive to t+0.5Δt_N after the velocity update
    // =========================================================================================
    // Save velocity at half Δt_N step (for even number of midpoint fock steps NMidpointFockSteps>2)
    if(isHalfNFockStep and not OddNFockSteps and mdOptions.nMidpointFockSteps != 2) {
      velocityHalfTN = velocity;
      curState.ptimeHalf = curState.ptime;
      std::cout << "  *** Saving Half-Step Velocity from p( t = " << curState.ptimeHalf << " au) ***" << std::endl;
    }

    // Update velocity (at half-step) to be full-step (where geometry is at)
    // p(t+Δt) = p(t+0.5Δt) + 0.5 * g(t+Δt) / m * Δt
    if ( (isEhrenfest && !doGrad) || (isBOMD && !firstStep) ) {
      double dt = curState.currentStepSize;
      std::cout << "  *** Updating Velocity from p( t = " << curState.ptime << " au) to p( t = " << curState.ptime + dt/2 << " au) ***"<< std::endl;
      //printCurrentVelocity(molecule);
      if ((magAmp[0] != 0.) || (magAmp[1] != 0.) || (magAmp[2] != 0.))
        velocityVV_EXPK1(molecule, velocity, velocity, gradient, dt, emPert);
      else
        velocityVV(molecule, velocity, velocity, gradient, dt);
      curState.ptime += dt/2;
    }

    // Save velocity at half Δt_N step (for odd number of midpoint fock steps NMidpointFockSteps>=3)
    if(isHalfNFockStep and OddNFockSteps) {
      velocityHalfTN = velocity;
      curState.ptimeHalf = curState.ptime;
      std::cout << "  *** Saving Half-Step Velocity from p( t = " << curState.ptimeHalf << " au) ***" << std::endl;
    }



    // =========================================================================================
    // Compute important quantities and print/save information
    // =========================================================================================
    // compute kinetic energy
    computeKineticEnergy(molecule);

    // Compute total energies
    currentTotalEnergy = electronicPotentialEnergy + nuclearKineticEnergy;
    if(firstStep and mdOptions.restoreFromNuclearStep == 0) totalEnergy0 = currentTotalEnergy;
    if(firstStep) previousTotalEnergy = currentTotalEnergy;
    
    // output important dynamic information
    if(print and doGrad) {
      printMDInfo(molecule, currentTotalEnergy);
      if (mdOptions.printProperty) {
        ss->computeProperties(emPert);
        ss->printProperties();
      }
    }
    previousTotalEnergy = currentTotalEnergy;

    // At this point we have full-step geom, time, g, v. Save these to bin file
    if(doGrad or mdOptions.saveAllGeometry) saveState(molecule, ss);

    if (std::any_of(cubes.begin(), cubes.end(), [](const std::shared_ptr<CubeGen>& ptr) { return ptr != nullptr; })) 
      ss->runCube(cubes, "_MDStep"+std::to_string(curState.iStep), std::make_shared<Molecule>(molecule));



    // =========================================================================================
    // Compute half-step velocity &&
    // Update geometry
    // =========================================================================================
    // Determine the geometry update step:
    // If JobType is BOMD, geometry move in full Δt_N step
    // If JobType is Ehrenfest, then geometry move in Δt_{N_q} step ( Δt_{N_q} = \dfrac{ Δt_N }  { m } )
    curState.currentStepSize = mdOptions.nMidpointFockSteps == 0 ? 
        mdOptions.timeStepAU : mdOptions.timeStepAU/ mdOptions.nMidpointFockSteps;
    // For Ehrenfest job, if this step evaluates gradient, that means we are at full Δt_N step 
    // If so, we need take half Δt_{N_q} to start mid-point fock algorithm
    if(mdOptions.nMidpointFockSteps != 0 and doGrad) curState.currentStepSize /= 2;

    // Update velocity to be half-step
    // p(t+0.5Δt) = p(t) + 0.5 * g(t) / m * Δt
    double dt = curState.currentStepSize;
    std::cout << "  *** Updating & saving Velocity from p( t = " << curState.ptime << " au) to p( t = " << curState.ptime + dt/2 << " au) ***"<< std::endl;
    //printCurrentVelocity(molecule);
    if ((magAmp[0] != 0.) || (magAmp[1] != 0.) || (magAmp[2] != 0.))
      velocityVV_EXPK1(molecule, velocity, velocity, gradient, dt, emPert);
    else
      velocityVV(molecule, velocity, velocity, gradient, dt);
    curState.ptime += dt/2;

    // Save velocity at half Δt_N step (special case for NMidpointFockSteps=2)
    // For NMidpointFockSteps=2, p will arrive to t+0.5Δt_N after this half-step velocity update
    if(isHalfNFockStep and mdOptions.nMidpointFockSteps == 2) {
      velocityHalfTN = velocity;
      curState.ptimeHalf = curState.ptime;
      std::cout << "  *** Saving Half-Step Velocity from p( t = " << curState.ptimeHalf << " au) ***" << std::endl;
    }

    // Update geometry using velocity at half-step
    // x(t+Δt) = x(t) + 0.5 * p(t+0.5Δt) * Δt
    std::cout << "  *** Updating Geometry from x( t = " << curState.time << " au) to x( t = " << curState.time+curState.currentStepSize << " au) ***"<< std::endl;
    geometryVV(molecule, curState.currentStepSize); 
    molecule.update();
    curState.time += curState.currentStepSize;



    // =========================================================================================
    // Set up for next RT simulation
    // =========================================================================================
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

      nuclearKineticEnergy += 0.5*velocity[i  ]*velocity[i  ]*atom.atomicMass*AUPerAMU;
      nuclearKineticEnergy += 0.5*velocity[i+1]*velocity[i+1]*atom.atomicMass*AUPerAMU;
      nuclearKineticEnergy += 0.5*velocity[i+2]*velocity[i+2]*atom.atomicMass*AUPerAMU;
      
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