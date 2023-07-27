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

    void MolecularDynamics::initializeMD(Molecule& molecule){

      std::vector<double> velocityHalfTime(3*molecule.nAtoms, 0.);
      std::vector<double> velocityCurrent(3*molecule.nAtoms, 0.);
      std::vector<double> acceleration(3*molecule.nAtoms, 0.);
      nuclearKineticEnergy = 0.0;

    }

    
    /**
     *  \brief Updates the positions of the classical nuclei. Recompute
     *  all member data to fit with the new positions
     *
     *  \param[in] pos  An array holding the new positions
     */
    void MolecularDynamics::update(bool print,
                Molecule &molecule,
                bool firstStep)
    {

      // Update time data
      if( firstStep )
        curState.iStep = 0;
      else
        curState.iStep++;

      curState.xTime = curState.iStep * curState.stepSize;
      if( molecularOptions_.nMidpointFockSteps != 0 )
        curState.xTime /= molecularOptions_.nMidpointFockSteps;

      auto doGrad = molecularOptions_.nMidpointFockSteps == 0 || 
                    curState.iStep % molecularOptions_.nMidpointFockSteps == 0;

      // Update gradient whenever we restart midpoint fock
      if( doGrad ) {

        // If we have midpoint fock steps, we need to take the final fock step
        if( molecularOptions_.nMidpointFockSteps != 0 && !firstStep ) {
          double fock_dt = molecularOptions_.timeStepAU/(molecularOptions_.nMidpointFockSteps + 1);
          geometryVV(molecule, gradientCurrent, fock_dt);
          molecule.update();
          electronicPotentialEnergy = finalMidpointFock(curState.iStep*fock_dt);
        }

        gradientCurrent = gradientGetter();

      }

      bool moveGeometry = true;
      bool moveVelocity = doGrad;

      size_t i = 0;

      if(print and doGrad) {
        std::cout << std::endl;
        std::cout << "MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD"<<std::endl;
        std::cout << "Molecular Dynamics Information for Step "<<std::setw(8)<<curState.iStep<<std::endl;
        std::cout << std::scientific << std::setprecision(12);

        std::cout << std::endl<<"Molecular Geometry: (Bohr)"<<std::endl;
        for( Atom& atom : molecule.atoms ) {

          std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                    << std::right <<"  X= "<< std::setw(19) << atom.coord[0]
                    << std::right <<"  Y= "<< std::setw(19) << atom.coord[1]
                    << std::right <<"  Z= "<< std::setw(19) << atom.coord[2] <<std::endl;
          i += 3;
 
        }
      }


      // if velocity Verlet
      if(moveVelocity) {
        velocityVV(molecule, gradientCurrent, molecularOptions_.timeStepAU, firstStep);
        // compute kinetic energy
        computeKineticEnergy(molecule);
      }
      if(moveGeometry) {
        geometryVV(molecule, gradientCurrent, molecularOptions_.timeStepAU/(molecularOptions_.nMidpointFockSteps + 1)); 
      }

      // update other quantities
      molecule.update();


      // Compute total energies
      double totalEnergy = electronicPotentialEnergy + nuclearKineticEnergy;
      if(firstStep) totalEnergy0 = totalEnergy;
      if(firstStep) previousTotalEnergy = totalEnergy;

      // output important dynamic information
      if(print and doGrad) {

        std::cout << std::setprecision(12);
        std::cout << "Velocity:"<<std::endl;
        i = 0;
        for( Atom& atom : molecule.atoms ) {

          std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                    << std::right <<"  X= "<< std::setw(19) <<  velocityCurrent[i  ]
                    << std::right <<"  Y= "<< std::setw(19) <<  velocityCurrent[i+1]
                    << std::right <<"  Z= "<< std::setw(19) <<  velocityCurrent[i+2]<<std::endl;
          i += 3;
 
        }

        std::cout << "Forces: (Hartrees/Bohr)"<<std::endl;
        i = 0;
        for( Atom& atom : molecule.atoms ) {

          std::cout <<"AtomicNumber = " << std::setw(4) <<  atom.atomicNumber 
                    << std::right <<"  X= "<< std::setw(19) <<  -gradientCurrent[i  ]
                    << std::right <<"  Y= "<< std::setw(19) <<  -gradientCurrent[i+1]
                    << std::right <<"  Z= "<< std::setw(19) <<  -gradientCurrent[i+2]<<std::endl;

          i += 3;
 
        }

        std::cout << std::setprecision(8);
        std::cout << std::endl<<"Time (fs): "<< std::right<<std::setw(16)<<curState.xTime*FSPerAUTime
                  << "  Time (au): "<<std::right<< std::setw(16)<<curState.xTime<<std::endl;

        std::cout <<  "EKin= " << std::right << std::setw(16) << nuclearKineticEnergy
                  << " EPot= " << std::right << std::setw(16) << electronicPotentialEnergy
                  << " ETot= " << std::right << std::setw(16) << totalEnergy << " a.u."<<std::endl;

        std::cout << "ΔETot (current-previous)= " << std::right << std::setw(16)<< totalEnergy-previousTotalEnergy
                  << " ΔETot (cumulative)= " << std::right << std::setw(16)<< totalEnergy-totalEnergy0<< " a.u."<<std::endl;


        std::cout << std::setprecision(12);
        std::cout << std::endl<<"Predicted Molecular Geometry: (Bohr)"<<std::endl;
        i = 0;
        for( Atom& atom : molecule.atoms ) {

          std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                    << std::right <<"  X= "<< std::setw(19) << atom.coord[0]
                    << std::right <<"  Y= "<< std::setw(19) << atom.coord[1]
                    << std::right <<"  Z= "<< std::setw(19) << atom.coord[2] <<std::endl;
          i += 3;
 
        }
      }

      previousTotalEnergy = totalEnergy;

      if( print && moveGeometry )
        std::cout << "  *** Moving Nuclei ***" << std::endl;

    }

    // Advance the velocity using the Verlet
    void MolecularDynamics::velocityVV(Molecule &molecule, std::vector<double> gradientCurrent, double timeStep, bool firstStep){

      size_t i = 0;

      // loop over atoms
      for( Atom& atom : molecule.atoms ) {
        //compute acceleration = -g/m
        acceleration[i  ] = -gradientCurrent[i  ]/(AUPerAMU*atom.atomicMass);
        acceleration[i+1] = -gradientCurrent[i+1]/(AUPerAMU*atom.atomicMass);
        acceleration[i+2] = -gradientCurrent[i+2]/(AUPerAMU*atom.atomicMass);

        //advance the half-time velocity to the current step 
        //v(t+1) = v(t+1/2) + 1/2dT∙a(t+1)
        if(not firstStep) {
          velocityCurrent[i  ] = velocityHalfTime[i  ] + 0.5*timeStep*acceleration[i  ];
          velocityCurrent[i+1] = velocityHalfTime[i+1] + 0.5*timeStep*acceleration[i+1];
          velocityCurrent[i+2] = velocityHalfTime[i+2] + 0.5*timeStep*acceleration[i+2];
        }

        //prepare the next velocity at half-time
        //v(t+1/2) = v(t) + 1/2dT∙a(t)
        velocityHalfTime[i  ] = velocityCurrent[i  ] + 0.5*timeStep*acceleration[i  ];
        velocityHalfTime[i+1] = velocityCurrent[i+1] + 0.5*timeStep*acceleration[i+1];
        velocityHalfTime[i+2] = velocityCurrent[i+2] + 0.5*timeStep*acceleration[i+2];

        i+=3;
      }

    }

    // Advance the geometry using the velocity
    void MolecularDynamics::geometryVV(Molecule &molecule, std::vector<double> gradientCurrent, double timeStep){

      size_t i = -3;

      // loop over atoms
      for( Atom& atom : molecule.atoms ) {
        i+=3;

        if (atom.quantum && NEODynamicsOpts.QProtMoveAlg == FIXED ) continue;
        //advance the geometry to the next time 
        //r(t+1) = r(t) + dT∙v(t+1/2)
        atom.coord[0] += timeStep*velocityHalfTime[i  ]; // x
        atom.coord[1] += timeStep*velocityHalfTime[i+1]; // y
        atom.coord[2] += timeStep*velocityHalfTime[i+2]; // z

        
      }
 
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
      double nextStep = curState.stepSize;
      if( molecularOptions_.nMidpointFockSteps != 0 )
        nextStep /= molecularOptions_.nMidpointFockSteps;
      return (curState.xTime + nextStep) < tMax;
    }

}