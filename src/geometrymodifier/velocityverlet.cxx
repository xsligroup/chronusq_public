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
  // Advance the velocity using the Verlet
  void MolecularDynamics::velocityVV(Molecule &molecule, std::vector<double>& vIn, std::vector<double>& vOut, std::vector<double> gradient, double timeStep) {
  
    size_t i = 0;
  
    // loop over atoms
    for( Atom& atom : molecule.atoms ) {
      //compute acceleration = -g/m
      acceleration[i  ] = -gradient[i  ]/(AUPerAMU*atom.atomicMass);
      acceleration[i+1] = -gradient[i+1]/(AUPerAMU*atom.atomicMass);
      acceleration[i+2] = -gradient[i+2]/(AUPerAMU*atom.atomicMass);
  
      vOut[i  ] = vIn[i  ] + 0.5*timeStep*acceleration[i  ];
      vOut[i+1] = vIn[i+1] + 0.5*timeStep*acceleration[i+1];
      vOut[i+2] = vIn[i+2] + 0.5*timeStep*acceleration[i+2];
  
      // copy values to each atom
      atom.velocity[0] = vOut[i  ];
      atom.velocity[1] = vOut[i+1];
      atom.velocity[2] = vOut[i+2];
      
      i+=3;
    }
  
  }
  
  // Advance the geometry using the velocity
  void MolecularDynamics::geometryVV(Molecule &molecule, double timeStep){
  
    size_t i = -3;
  
    // loop over atoms
    for( Atom& atom : molecule.atoms ) {
      i+=3;
  
      if (atom.quantum && !NEODynamicsOpts.tpb) continue;
      //advance the geometry to the next time 
      //r(t+1) = r(t) + dT∙v(t+1/2)
      atom.coord[0] += timeStep*velocity[i  ]; // x
      atom.coord[1] += timeStep*velocity[i+1]; // y
      atom.coord[2] += timeStep*velocity[i+2]; // z
      
    }
  
  }

} 