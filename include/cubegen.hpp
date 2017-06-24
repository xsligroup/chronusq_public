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
#include <quantum.hpp>
#include <cubegenoptions.hpp>

#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <ctime> 

namespace ChronusQ {
  
  // Resolution types for cube
  enum class RES_TYPE {
    COARSE,
    MEDIUM,
    FINE
  };

  // Particle type for cube
  enum PAR_TYPE {
    ELECTRONIC,
    PROTONIC,
  };

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
      std::shared_ptr<std::ofstream> cubeFile_ = nullptr;
      // Allows for the generation of real and imaginary of an orbital
      // in a single call when generating orbital cubes
      std::array<size_t,3> voxelGrid_;
      std::array<double,3> voxelUnits_;
      RES_TYPE res_;
      double cubePadding_ = 5.0;
      std::shared_ptr<Molecule> mol_;
      std::shared_ptr<BasisSet> basis_;

      // container for general workflow cubegen options
      CubeGenOptions cubeOpts_;

      // Helper functions for precomputing the basis (if memory permits)
      void ComputeBasis();
      std::vector<double> evaluated_basis;
      double * EvalShellSetAtPoint(int,int,int);
      
    public:

      /**
       * @brief Constructor that assumes default dimensions
       * 
      */
      CubeGen(std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> basis){
        res_ = RES_TYPE::COARSE;
        voxelGrid_ = {80,80,80};
        voxelUnits_ = {0.1,0.1,0.1};
        mol_ = mol;
        basis_ = basis;
        ComputeBasis();
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses specified
       * grid and step input
       * 
       * @param mol_ std::shared_ptr<Molecule> 
       * @param basis_ BasisSet 
       * @param voxelGrid the dimensions of the grid that holds the information
       *  of the surface
       * @param voxelUnits the increments between datapoints for the voxelGrid
      */
      CubeGen(std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> basis,
      std::array<size_t, 3> voxelGrid, std::array<double, 3> voxelUnits) {
        voxelGrid_ = voxelGrid;
        voxelUnits_ = voxelUnits;
        mol_ = mol;
        basis_ = basis;
        ComputeBasis();
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses
       * resolution input
       * 
       * @param mol_ std::shared_ptr<Molecule> 
       * @param basis_ BasisSet 
       * @param res resolution of visualization specified by user
      */
      CubeGen(std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> basis,
      std::string resString,
      double cubePadding = 3.0) {
        res_ = inputToRes(resString);
        cubePadding_ = cubePadding;
        mol_ = mol;
        basis_ = basis;
        // call this to create grid
        calculateVoxelDimensions();
        ComputeBasis();
      }


      RES_TYPE inputToRes(const std::string& input) {

        std::map<std::string, RES_TYPE> strToRes = {
          {"COARSE", RES_TYPE::COARSE},
          {"MEDIUM", RES_TYPE::MEDIUM},
          {"FINE", RES_TYPE::FINE}
        };

        return strToRes[input];
      }

      // >>> Some general helper functions

      /**
       * @brief Creates cube file for a given title 
       * 
      */
      void createNewCube(std::string cubeT) {
        cubeFile_ = std::make_shared<std::ofstream>(cubeT+".cube");
      }

      // There is no setter for cube padding because
      // it is only used in the constructor atm.
      /**
       * @brief returns cube padding 
       * 
      */
      double getCubePad() {
        return cubePadding_;
      }

      /**
       * @brief returns cube options 
       * 
      */
      CubeGenOptions &getCubeOptions() {
        return cubeOpts_;
      }

      // >>> Grid-related functions
      // see src/cubegen
      void calculateVoxelDimensions();
      std::vector<double> calcCenter();
      void updateMolAndBasis(std::shared_ptr<Molecule> mol);

      // >>> High-level functions
      void writeSummary(std::string fileSum);
      template <typename LocMatsT>
      void evalDenCube(std::string filePref,std::shared_ptr<cqmatrix::PauliSpinorMatrices<LocMatsT>>, double particleCharge = -1.0, bool skipoutput = false);
      template <typename LocMatsT, typename ValManipOp>
      void evalOrbCube(std::string filePref, LocMatsT* MOBase, size_t LDMO, std::vector<size_t> whichMOs, ValManipOp op);

      // >>> Property evaluation functions
      template <typename LocMatsT>
      void evalDenCompCube(LocMatsT*, double pCharge=-1.0);
      template <typename LocMatsT, typename ValManipOp>
      void evalOrbCompCube(LocMatsT *,size_t,size_t,ValManipOp);
    

  };

};

