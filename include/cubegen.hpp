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
#ifndef __INCLUDED_CUBEGEN_HPP__
#define __INCLUDED_CUBEGEN_HPP__
#pragma once

#include <chronusq_sys.hpp>
#include <quantum.hpp>
#include <singleslater.hpp>

#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <ctime> 

namespace ChronusQ {

  // Different types of cube outputs
  enum class CUBE_TYPE {
    _CHARGE_DENSITY,
    _SPIN_DENSITY,
    _MAGNETIZATION,
    _ELECTROSTATIC_POTENTIAL
  };
  
  // Resolution types for cube
  enum class RES_TYPE {
    COARSE,
    MEDIUM,
    FINE
  };


  template <typename MatsT, typename IntsT>
  /**
   * @brief Generates a cubefile of specfied surface.
   * 
   * The CubeGen class creates a cubefile of a specified surface
   * for a calculation. It takes the information from the calculation
   * and creates a cubefile, which can be used to visualize surfaces
   * through other scientific visualization programs
  */
  class CubeGen {

    private : 
      std::string cubeTitle;
      std::shared_ptr<std::ofstream> cubeFile;
      CUBE_TYPE cubeType;
      std::array<size_t,3> voxelGrid;
      std::array<double,3> voxelUnits;
      std::shared_ptr<SingleSlaterBase> ref_;
      RES_TYPE res;
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM_;
      double x_edge;
      double y_edge;
      double z_edge;
      static constexpr double CUBE_PADDING = 3.0;
      
    public:
      // default constructor
      CubeGen() {
        std::string cubeTitle = "";
        std::shared_ptr<std::ofstream> cubeFile = nullptr;
        RES_TYPE res = RES_TYPE::COARSE;
        std::array<size_t,3> voxelGrid = {80,80,80};
        std::array<double,3> voxelUnits = {0.1,0.1,0.1};
        std::shared_ptr<SingleSlaterBase> ref_ = nullptr;
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM_ = nullptr;
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses specified
       * grid and step input
       * 
       * @param ref_ reference to calculation type 
       * @param onePDM_ reference to the matrices to visualize
       * @param cubeTitle the name of the cubefile to be outputted
       * @param voxelGrid the dimensions of the grid that holds the information
       *  of the surface
       * @param voxelUnits the increments between datapoints for the voxelGrid
      */
      CubeGen(std::shared_ptr<SingleSlaterBase> ref_,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM_,
      std::string cubeTitle ,
      std::array<size_t, 3> voxelGrid, std::array<double, 3> voxelUnits) {
        this->cubeTitle = cubeTitle;
        this->res = RES_TYPE::COARSE;
        this->voxelGrid = voxelGrid;
        this->voxelUnits = voxelUnits;
        this->cubeFile = std::make_shared<std::ofstream>(cubeTitle);
        this->ref_ = ref_;
        this->onePDM_ = onePDM_;
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses
       * resolution input
       * 
       * @param ref_ reference to calculation type 
       * @param onePDM_ reference to the matrices to visualize
       * @param cubeTitle the name of the cubefile to be outputted
       * @param res resolution of visualization specified by user
      */
      CubeGen(std::shared_ptr<SingleSlaterBase> ref_,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM_,
      std::string cubeTitle , std::string res) {
        this->cubeTitle = cubeTitle;
        this->res = inputToRes(res);
        this->voxelGrid = {80,80,80};
        this->voxelUnits = {0.1,0.1,0.1};
        this->cubeFile = std::make_shared<std::ofstream>(cubeTitle);
        this->ref_ = ref_;
        this->onePDM_ = onePDM_;
        // call this to create grid
        calculateVoxelDimensions();
      }


      RES_TYPE inputToRes(const std::string& input) {

        std::map<std::string, RES_TYPE> strToRes = {
          {"COARSE", RES_TYPE::COARSE},
          {"MEDIUM", RES_TYPE::MEDIUM},
          {"FINE", RES_TYPE::FINE}
        };

        return strToRes[input];
      }


      /**
       * @brief Generates the dimensions of voxel grid
       * 
       * This uses the user input of coarse, medium, or fine
       * in the RES option, and uses that to calculate
       * the dimensions of the voxel grid, creating correct spacing
       * and quality of the output
      */
      void calculateVoxelDimensions() {

        std::array<double,3> maxDimensions = {0.0,0.0,0.0};
        // 1 is coarse
        // 2 is med
        // 3 is fine
        std::map<RES_TYPE,double> resolution = {
            {RES_TYPE::COARSE, 1.0/3},
            {RES_TYPE::MEDIUM, 1.0/6},
            {RES_TYPE::FINE, 1.0/12}
        };
        // set steps based on resolution
        this->voxelUnits = {resolution[res], resolution[res], resolution[res]};
        
        // get array of the greatest dimensions to create
        // total grid lengths
        for( auto &atom : ref_->molecule().atoms ) {
          for(int i = 0; i < std::size(maxDimensions); i++) {
            maxDimensions[i] = std::max(maxDimensions[i], std::abs(atom.coord[i]));
          }
        }

        std::array<double,3> maxGrid = {0,0,0};
        // create the correct amount of space for molecule
        for(int i = 0; i < std::size(maxGrid); i++) {
          // add 1 angstrom of padding
            maxGrid[i] = std::abs(maxDimensions[i]) * 2 + CUBE_PADDING;
        }
        // using largest dimensions, convert space to the voxel spaces
        for(int i = 0; i < std::size(maxDimensions); i++) {
            this->voxelGrid[i] = std::ceil(maxGrid[i] * 2 / voxelUnits[i]);
        }

      }
      
      /**
       * @brief Calculates the 3d coordinates of center of system
       * 
       * This uses the inputs of particles to calculate the center
       * point of the system. This is not currently used, but could
       * have useful future application
       * @returns a vector of the center of system
      */
      std::vector<double> calcCenter() {
        std::vector<double> maxDimensions = {0.0,0.0,0.0};

        for( auto &atom : ref_->molecule().atoms ) {
          for(int i = 0; i < 3; i++) {
            maxDimensions[i] += atom.coord[i];
          }
        }

        size_t num_atoms = ref_->molecule().atoms.size();

        std::vector<double> centerPoint = {maxDimensions[0] / num_atoms,
                                             maxDimensions[1] / num_atoms,
                                             maxDimensions[2] / num_atoms};

        return centerPoint;
      }

      // call create cube file depending on type
      /**
       * @brief Generates the cubefile for specified cube
       * 
       * In procedural, this will be called for each of the surfaces
       *  the user initialized in the input
      */
      void evalCube(CUBE_TYPE cubeCall) {

        // smgargner229
        std::cout << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "Generating Cube file for: " << cubeTitle << std::endl;
        std::cout << "Using " << voxelGrid[0] << "," << voxelGrid[1] << "," << voxelGrid[2] << " Points" << std::endl;
        std::cout << "With steps: " << voxelUnits[0] << "," << voxelUnits[1] << "," << voxelUnits[2] << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << std::endl;

        
        ProgramTimer::tick("Cube Eval");

        // assign to object field
        this->cubeType = cubeCall;

        // write the cubefile header 
        writeSummary();

        if (cubeType == CUBE_TYPE::_CHARGE_DENSITY) {
          evalCDCube();
        }
        else
        {
           CErr("CubeGen for non-Charge Density Cubes NYI!");
         }
        ProgramTimer::tock("Cube Eval");

        std::cout << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << "Finished generating Cube file for: " << cubeTitle << std::endl;
        std::cout << "----------------------------------------------------" << std::endl;
        std::cout << std::endl;

        
      }

      /**
       * @brief Writes header of cube file
       * 
       * This creates the header at the beginning of the cube file,
       *  diplaying the dimensions of the surface in the cube format
      */
      void writeSummary() {
        std::map<CUBE_TYPE, std::string> cubeTypeOut = {
            {CUBE_TYPE::_CHARGE_DENSITY, "CHARGE DENSITY"},
            {CUBE_TYPE::_SPIN_DENSITY, "SPIN DENSITY"},
            {CUBE_TYPE::_ELECTROSTATIC_POTENTIAL, "ELECTROSTATIC POTENTIAL"}
        };
        std::vector<double> centerPoint = calcCenter();

        *cubeFile << cubeTitle << "\n";
        *cubeFile << std::fixed;

        *cubeFile << cubeTypeOut[cubeType];
        *cubeFile << " : GENERATED BY CHRONUSQ\n";


        size_t nAtoms = ref_->molecule().nAtoms;
        Molecule &mol = ref_->molecule();

        double x_edge = 0.0 - voxelUnits[0]*(double)voxelGrid[0]/2.0;
        double y_edge = 0.0 - voxelUnits[1]*(double)voxelGrid[1]/2.0;
        double z_edge = 0.0 - voxelUnits[2]*(double)voxelGrid[2]/2.0;

        *cubeFile << std::setprecision(6);

        *cubeFile << std::setw(6) << std::right << nAtoms;
        *cubeFile << std::setw(15) << x_edge;
        *cubeFile << std::setw(15) << y_edge;
        *cubeFile << std::setw(15) << z_edge;
        *cubeFile << "\n";


        *cubeFile << std::setw(6) << std::right << voxelGrid[0];
        *cubeFile << std::setw(15) << voxelUnits[0];
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << "\n";

        *cubeFile << std::setw(6) << std::right << voxelGrid[1];
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << std::setw(15) << voxelUnits[1];
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << "\n";

        *cubeFile << std::setw(6) << std::right << voxelGrid[2];
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << std::setw(15) << 0.;
        *cubeFile << std::setw(15) << voxelUnits[2];
        *cubeFile << "\n";

        for( auto &atom : ref_->molecule().atoms ) {

          *cubeFile << std::setw(6) << atom.atomicNumber;
          *cubeFile << std::setw(15) << 0.;

          *cubeFile << std::setw(15) << atom.coord[0];
          *cubeFile << std::setw(15) << atom.coord[1];
          *cubeFile << std::setw(15) << atom.coord[2];

          *cubeFile << "\n";

        }

      }

      /**
       * @brief Generates the cubefile for charge density
       * 
       * This is called from evalCube.
      */
      void evalCDCube()
      {
        // gets number of basis sets
        size_t NB = ref_->nAlphaOrbital();

        std::vector<libint2::Shell> &shells = ref_->basisSet().shells;

        for(auto ix = 0l; ix < voxelGrid[0]; ix++) {
          for(auto iy = 0l; iy < voxelGrid[1]; iy++) {
            for(auto iz = 0l; iz < voxelGrid[2]; iz++) {

              // voxel calculation
              std::array<double,3> pt = {
              (ix-(int)voxelGrid[0]/2) * voxelUnits[0],
              (iy-(int)voxelGrid[1]/2) * voxelUnits[1],
              (iz-(int)voxelGrid[2]/2) * voxelUnits[2]
              };

              std::vector<double> BASIS(NB,0.);
              std::vector<MatsT> SCR(NB,0.);

              MatsT val = 0;

              evalShellSet(ref_->memManager,NOGRAD,ref_->basisSet().shells,&pt[0],1,&BASIS[0],false);

              blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, 1, NB, NB,1.,
              &BASIS[0], NB, onePDM_->S().pointer(), NB,0., &SCR[0], 1);

              val = blas::dot(NB,&SCR[0],1,&BASIS[0],1);

              *cubeFile << std::right << std::setw(15) << std::setprecision(5)
              << std::scientific << std::uppercase << std::real(ref_->particle.charge *  val);

              if( iz % 6 == 5 ) *cubeFile << "\n";

            }
            *cubeFile << "\n";
          }
        }
      };

  };

};

#endif //  __INCLUDED_CUBEGEN_HPP__
