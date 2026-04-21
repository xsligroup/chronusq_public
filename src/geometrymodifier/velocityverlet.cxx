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

//#define EXP_VV_DEBUG

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

  // Exponential (EXP) method Order 1
  // v(t+0.5Δt) = u(t)^1/2 * v(t) + u(t)^1/4 * f(t)
  // See Helgaker (2024) Mol. Phys. 122(5) e2259008: Algorithm 4
  void MolecularDynamics::velocityVV_EXPK1(Molecule &molecule, std::vector<double>& vIn, std::vector<double>& vOut, 
    std::vector<double> gradient, double timeStep, EMPerturbation& emPert) {
  
    size_t i = 0;

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    double w_matrix[9];
    double u_half_matrix[9]; // u(t)^1/2
    double u_quar_matrix[9]; // u(t)^1/4

#ifdef EXP_VV_DEBUG
    std::vector<double> vIn2(vIn);
    std::vector<double> vOut2(vOut);
#endif
    
    // loop over atoms
    // std::cout << "  *** Calculating Velocity-Dependent Forces ***" << std::endl;
    for( Atom& atom : molecule.atoms ) {

      if (atom.quantum && !NEODynamicsOpts.tpb) continue;

      // build w matrix
      w_matrix[0] =  0.0;
      w_matrix[1] =  magAmp[2] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass); 
      w_matrix[2] = -magAmp[1] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass); 
      w_matrix[3] = -magAmp[2] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass);
      w_matrix[4] =  0.0;
      w_matrix[5] =  magAmp[0] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass);
      w_matrix[6] =  magAmp[1] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass);
      w_matrix[7] = -magAmp[0] * (-1.0) * atom.nucCharge / (AUPerAMU*atom.atomicMass);
      w_matrix[8] =  0.0;

      // U = exp(wΔt)
      for (size_t j = 0; j < 9; ++j) w_matrix[j] *= (0.5*timeStep);
      MatExp(3, w_matrix, 3, u_half_matrix, 3);
      for (size_t j = 0; j < 9; ++j) w_matrix[j] *= 0.5;
      MatExp(3, w_matrix, 3, u_quar_matrix, 3);

#ifdef EXP_VV_DEBUG
      std::cout << "Gradient "  << i << " " << gradient[i  ] << " " << gradient[i+1] << " " << gradient[i+2] << std::endl;
      std::cout << "Gradient "  << i << " " << atom.nucCharge * (vIn[i+1] * magAmp[2] - vIn[i+2] * magAmp[1])
                                     << " " << atom.nucCharge * (vIn[i+2] * magAmp[0] - vIn[i+0] * magAmp[2]) 
                                     << " " << atom.nucCharge * (vIn[i+0] * magAmp[1] - vIn[i+1] * magAmp[0]) << std::endl;
#endif

      //compute acceleration = -g/m
      acceleration[i  ] = -gradient[i  ]/(AUPerAMU*atom.atomicMass);
      acceleration[i+1] = -gradient[i+1]/(AUPerAMU*atom.atomicMass);
      acceleration[i+2] = -gradient[i+2]/(AUPerAMU*atom.atomicMass);

      // v(t+0.5Δt) = u(t)^1/2 * v(t) + u(t)^1/4 * f(t)
      vOut[i  ] = u_half_matrix[0] * vIn[i  ] 
	              + u_half_matrix[3] * vIn[i+1] 
	              + u_half_matrix[6] * vIn[i+2]
                + 0.5*timeStep*u_quar_matrix[0] * acceleration[i  ] 
				        + 0.5*timeStep*u_quar_matrix[3] * acceleration[i+1] 
				        + 0.5*timeStep*u_quar_matrix[6] * acceleration[i+2];
      vOut[i+1] = u_half_matrix[1] * vIn[i  ] 
	              + u_half_matrix[4] * vIn[i+1]
	              + u_half_matrix[7] * vIn[i+2]
                + 0.5*timeStep*u_quar_matrix[1] * acceleration[i  ] 
				        + 0.5*timeStep*u_quar_matrix[4] * acceleration[i+1] 
				        + 0.5*timeStep*u_quar_matrix[7] * acceleration[i+2];
      vOut[i+2] = u_half_matrix[2] * vIn[i  ] 
	              + u_half_matrix[5] * vIn[i+1] 
	              + u_half_matrix[8] * vIn[i+2]
                + 0.5*timeStep*u_quar_matrix[2] * acceleration[i  ] 
			         	+ 0.5*timeStep*u_quar_matrix[5] * acceleration[i+1] 
			        	+ 0.5*timeStep*u_quar_matrix[8] * acceleration[i+2];

      // copy values to each atom
      atom.velocity[0] = vOut[i  ];
      atom.velocity[1] = vOut[i+1];
      atom.velocity[2] = vOut[i+2];
      
      i+=3;
    }

#ifdef EXP_VV_DEBUG

std::cout << "EXP_VV_DEBUG: vIn:" << std::endl;
i=0;
for( Atom& atom : molecule.atoms ) {
  std::cout << vIn[i  ] << " " << vIn[i+1] << " " << vIn[i+2] << std::endl;
  i+=3;
}

std::cout << "EXP_VV_DEBUG: dT:" << std::endl;
std::cout << timeStep << std::endl;

    std::cout << "EXP_VV_DEBUG: VOut:" << std::endl;
    i=0;
    for( Atom& atom : molecule.atoms ) {
      std::cout << vOut[i  ] << " " << vOut[i+1] << " " << vOut[i+2] << std::endl;
      i+=3;
    }

    std::cout << "EXP_VV_DEBUG: VOut (Pure VV):" << std::endl;
    i=0;
    for( Atom& atom : molecule.atoms ) {

      if (atom.quantum && !NEODynamicsOpts.tpb) continue;

//      std::cout << "Gradient "  << i << " " << gradient[i  ] << " " << gradient[i+1] << " " << gradient[i+2] << std::endl;
//      std::cout << "Gradient "  << i << " " << atom.nucCharge * (vIn2[i+1] * magAmp[2] - vIn2[i+2] * magAmp[1])
//                                     << " " << atom.nucCharge * (vIn2[i+2] * magAmp[0] - vIn2[i+0] * magAmp[2]) 
//                                     << " " << atom.nucCharge * (vIn2[i+0] * magAmp[1] - vIn2[i+1] * magAmp[0]) << std::endl;

      gradient[i  ] -= atom.nucCharge * (vIn2[i+1] * magAmp[2] - vIn2[i+2] * magAmp[1]);
      gradient[i+1] -= atom.nucCharge * (vIn2[i+2] * magAmp[0] - vIn2[i+0] * magAmp[2]);
      gradient[i+2] -= atom.nucCharge * (vIn2[i+0] * magAmp[1] - vIn2[i+1] * magAmp[0]);

      acceleration[i  ] = -gradient[i  ]/(AUPerAMU*atom.atomicMass);
      acceleration[i+1] = -gradient[i+1]/(AUPerAMU*atom.atomicMass);
      acceleration[i+2] = -gradient[i+2]/(AUPerAMU*atom.atomicMass);

      vOut2[i  ] = vIn2[i  ] + 0.5*timeStep*acceleration[i  ];
      vOut2[i+1] = vIn2[i+1] + 0.5*timeStep*acceleration[i+1];
      vOut2[i+2] = vIn2[i+2] + 0.5*timeStep*acceleration[i+2];

      std::cout << vOut2[i  ] << " " << vOut2[i+1] << " " << vOut2[i+2] << std::endl;
      i+=3;
    }

#endif
  
  } //void MolecularDynamics::velocityVV_EXPK1
  
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