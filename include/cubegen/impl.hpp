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

#include <chronusq_sys.hpp>
#include <cubegen.hpp>

namespace ChronusQ {

    /**
     * @brief Writes header of cube file
     * 
     * This creates the header at the beginning of the cube file,
     *  diplaying the dimensions of the surface in the cube format
    */
    void CubeGen::writeSummary(std::string fileSummary) {
      std::vector<double> centerPoint = calcCenter();

      *cubeFile_ << fileSummary << "\n";
      *cubeFile_ << std::fixed;

      *cubeFile_ << "GENERATED BY CHRONUSQ\n";


      size_t nAtoms = mol_->nAtoms;

      double x_edge = 0.0 - voxelUnits_[0]*(double)voxelGrid_[0]/2.0;
      double y_edge = 0.0 - voxelUnits_[1]*(double)voxelGrid_[1]/2.0;
      double z_edge = 0.0 - voxelUnits_[2]*(double)voxelGrid_[2]/2.0;

      *cubeFile_ << std::setprecision(6);

      *cubeFile_ << std::setw(6) << std::right << nAtoms;
      *cubeFile_ << std::setw(15) << x_edge;
      *cubeFile_ << std::setw(15) << y_edge;
      *cubeFile_ << std::setw(15) << z_edge;
      *cubeFile_ << "\n";


      *cubeFile_ << std::setw(6) << std::right << voxelGrid_[0];
      *cubeFile_ << std::setw(15) << voxelUnits_[0];
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << "\n";

      *cubeFile_ << std::setw(6) << std::right << voxelGrid_[1];
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << std::setw(15) << voxelUnits_[1];
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << "\n";

      *cubeFile_ << std::setw(6) << std::right << voxelGrid_[2];
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << std::setw(15) << 0.;
      *cubeFile_ << std::setw(15) << voxelUnits_[2];
      *cubeFile_ << "\n";

      for( auto &atom : mol_->atoms ) {

        *cubeFile_ << std::setw(6) << atom.atomicNumber;
        *cubeFile_ << std::setw(15) << 0.;

        *cubeFile_ << std::setw(15) << atom.coord[0];
        *cubeFile_ << std::setw(15) << atom.coord[1];
        *cubeFile_ << std::setw(15) << atom.coord[2];

        *cubeFile_ << "\n";

      }

    }

    /**
     * @brief create density cubegen files
     * 
    */
    template <typename LocMatsT>
    void CubeGen::evalDenCube(std::string fileNamePrefix,std::shared_ptr<cqmatrix::PauliSpinorMatrices<LocMatsT>> oPDM, double particleCharge, bool skipoutput) {

      if (not skipoutput) {
        std::cout << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "Generating Density Cube files" << std::endl;
        std::cout << "Using " << voxelGrid_[0] << "," << voxelGrid_[1] << "," << voxelGrid_[2] << " Points" << std::endl;
        std::cout << "With steps: " << voxelUnits_[0] << "," << voxelUnits_[1] << "," << voxelUnits_[2] << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << std::endl;
      }
      
      ProgramTimer::tick("Cube Eval");

      std::string denFileName = fileNamePrefix + "_DEN";

      // Scalar density
      std::string scalDenFileName = denFileName + "_S";
      createNewCube(scalDenFileName);
      if (not skipoutput) std::cout << "Writing scalar density to file: " << scalDenFileName << std::endl;
      writeSummary("Scalar Density");
      evalDenCompCube(oPDM->S().pointer(),particleCharge);

      // MZ density
      if( oPDM->hasZ() ){

        std::string mzDenFileName = denFileName + "_MZ";
        createNewCube(mzDenFileName);
        if (not skipoutput) std::cout << "Writing MZ density to file: " << mzDenFileName << std::endl;
        writeSummary("MZ Density");
        evalDenCompCube(oPDM->Z().pointer(),particleCharge);

      }

      if( oPDM->hasXY() ){

          // MX density
          std::string mxDenFileName = denFileName + "_MX";
          createNewCube(mxDenFileName);
          if (not skipoutput) std::cout << "Writing MX density to file: " << mxDenFileName << std::endl;
          writeSummary("MX Density");
          evalDenCompCube(oPDM->X().pointer(),particleCharge);

          // MY density
          std::string myDenFileName = denFileName + "_MY";
          createNewCube(myDenFileName);
          if (not skipoutput) std::cout << "Writing MY density to file: " << myDenFileName << std::endl;
          writeSummary("MY Density");
          evalDenCompCube(oPDM->Y().pointer(),particleCharge);

      }

      ProgramTimer::tock("Cube Eval");

      if (not skipoutput) {
        std::cout << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << std::endl;
      }
      
    }

    template <typename LocMatsT, typename ValManipOp>
    void CubeGen::evalOrbCube(std::string filenamePrefix, LocMatsT * MOBase, size_t LDMO, std::vector<size_t> whichMOs, ValManipOp op)
    {
      std::cout << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << "Generating Orbital Cube files" << std::endl;
      std::cout << "For Component: " << filenamePrefix << std::endl;
      std::cout << "Using " << voxelGrid_[0] << "," << voxelGrid_[1] << "," << voxelGrid_[2] << " Points" << std::endl;
      std::cout << "With steps: " << voxelUnits_[0] << "," << voxelUnits_[1] << "," << voxelUnits_[2] << std::endl;
      std::cout << "----------------------------------------------------" << std::endl;
      std::cout << std::endl;

      ProgramTimer::tick("Cube Eval");

      // Loop over requested orbitals
      for(const size_t & MOIndex : whichMOs)
      {
        createNewCube(filenamePrefix + "_MO_"+std::to_string(MOIndex+1));
        writeSummary(filenamePrefix + " MO " + std::to_string(MOIndex+1));
        evalOrbCompCube(MOBase,LDMO,MOIndex,op);
      }

      ProgramTimer::tock("Cube Eval");

    }

    void CubeGen::updateMolAndBasis(std::shared_ptr<Molecule> newMol){
      mol_ = newMol;
      basis_->updateNuclearCoordinates(*mol_);
      calculateVoxelDimensions();
      ComputeBasis();
    }

}