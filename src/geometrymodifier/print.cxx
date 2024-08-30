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
  // Print current geometry 
  void MolecularDynamics::printCurrentGeometry(Molecule &molecule){

    size_t i = 0;

    std::cout << std::endl;
    std::cout << "MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD-MD"<<std::endl;
    std::cout << "Molecular Dynamics Information for Step "<<std::setw(8)<<curState.iStep<<std::endl;
    std::cout << std::scientific << std::setprecision(16);

    std::cout << std::endl<<"Molecular Geometry: (Bohr)"<<std::endl;
    for( Atom& atom : molecule.atoms ) {

      std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                << std::right <<"  X= "<< std::setw(24) << atom.coord[0]
                << std::right <<"  Y= "<< std::setw(24) << atom.coord[1]
                << std::right <<"  Z= "<< std::setw(24) << atom.coord[2] <<std::endl;
      i += 3;
 
    }
  }

  // Print velocities, forces, energies, and predicted geometry 
  void MolecularDynamics::printMDInfo(Molecule &molecule, double totalEnergy){
    std::cout << std::setprecision(16);
    std::cout << "Velocity:"<<std::endl;
    size_t i = 0;
    for( Atom& atom : molecule.atoms ) {

      std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                << std::right <<"  X= "<< std::setw(24) <<  velocityCurrent[i  ]
                << std::right <<"  Y= "<< std::setw(24) <<  velocityCurrent[i+1]
                << std::right <<"  Z= "<< std::setw(24) <<  velocityCurrent[i+2]<<std::endl;
      i += 3;
 
    }

    std::cout << "Forces: (Hartrees/Bohr)"<<std::endl;
    i = 0;
    for( Atom& atom : molecule.atoms ) {

      std::cout <<"AtomicNumber = " << std::setw(4) <<  atom.atomicNumber 
                << std::right <<"  X= "<< std::setw(24) <<  -gradientCurrent[i  ]
                << std::right <<"  Y= "<< std::setw(24) <<  -gradientCurrent[i+1]
                << std::right <<"  Z= "<< std::setw(24) <<  -gradientCurrent[i+2]<<std::endl;

      i += 3;
 
    }

    std::cout << std::setprecision(8);
    std::cout << std::endl<<"Time (fs): "<< std::right<<std::setw(16)<<curState.time*FSPerAUTime
              << "  Time (au): "<<std::right<< std::setw(16)<<curState.time<<std::endl;

    std::cout << std::setprecision(16);
    std::cout <<  "EKin= " << std::right << std::setw(24) << nuclearKineticEnergy
              << " EPot= " << std::right << std::setw(24) << electronicPotentialEnergy
              << " ETot= " << std::right << std::setw(24) << totalEnergy << " a.u."<<std::endl;

    std::cout << "ΔETot (current-previous)= " << std::right << std::setw(24)<< totalEnergy-previousTotalEnergy
              << " ΔETot (cumulative)= " << std::right << std::setw(24)<< totalEnergy-totalEnergy0<< " a.u."<<std::endl;


    std::cout << std::setprecision(16);
    std::cout << std::endl<<"Predicted Molecular Geometry: (Bohr)"<<std::endl;
    i = 0;
    for( Atom& atom : molecule.atoms ) {

      std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                << std::right <<"  X= "<< std::setw(24) << atom.coord[0]
                << std::right <<"  Y= "<< std::setw(24) << atom.coord[1]
                << std::right <<"  Z= "<< std::setw(24) << atom.coord[2] <<std::endl;
      i += 3;
 
    }
  }



}