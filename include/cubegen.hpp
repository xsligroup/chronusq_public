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
#include <singleslater.hpp>

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
      std::shared_ptr<std::ofstream> cubeFile_;
      std::array<size_t,3> voxelGrid_;
      std::array<double,3> voxelUnits_;
      std::shared_ptr<SingleSlaterBase> ref_;
      RES_TYPE res_;
      bool denCube_ = false;
      std::string cubeFileName_;
      double cubePadding_ = 3.0;
      
    public:

      /**
       * @brief Constructor that assumes default dimensions but has ref 
       * 
      */
      CubeGen(std::shared_ptr<SingleSlaterBase> ref){
        cubeFileName_ = "";
        cubeFile_ = nullptr;
        res_ = RES_TYPE::COARSE;
        voxelGrid_ = {80,80,80};
        voxelUnits_ = {0.1,0.1,0.1};
        ref_ = ref;
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses specified
       * grid and step input
       * 
       * @param ref_ reference to calculation type 
       * @param cubeFileName the name of the cubefile to be outputted
       * @param voxelGrid the dimensions of the grid that holds the information
       *  of the surface
       * @param voxelUnits the increments between datapoints for the voxelGrid
      */
      CubeGen(std::shared_ptr<SingleSlaterBase> ref,
      std::string cubeFileName,
      std::array<size_t, 3> voxelGrid, std::array<double, 3> voxelUnits) {
        cubeFileName_ = cubeFileName;
        voxelGrid_ = voxelGrid;
        voxelUnits_ = voxelUnits;
        cubeFile_ = std::make_shared<std::ofstream>(cubeFileName_);
        ref_ = ref;
      }

      /**
       * @brief Creates a cubefile of a specified surface
       * 
       * This function is called in procedural to output a cube file of a 
       * specified surface to be visualized in a 3D format. Uses
       * resolution input
       * 
       * @param ref_ reference to calculation type 
       * @param cubeFileName the name of the cubefile to be outputted
       * @param res resolution of visualization specified by user
      */
      CubeGen(std::shared_ptr<SingleSlaterBase> ref,
      std::string cubeFileName , std::string resString,
      double cubePadding = 3.0) {
        cubeFileName_ = cubeFileName;
        res_ = inputToRes(resString);
        cubeFile_ = std::make_shared<std::ofstream>(cubeFileName_);
        cubePadding_ = cubePadding;
        ref_ = ref;
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

      // >>> Some general helper functions

      /**
       * @brief updates cube title 
       * 
      */
      void setCubeFileName(std::string cubeT) {
        cubeFileName_ = cubeT;
      }

      /**
       * @brief returns current cube title 
       * 
      */
      std::string getCubeFileName() {
        return cubeFileName_;
      }

      /**
       * @brief Creates cube file for a given title 
       * 
      */
      void createNewCube(std::string cubeT) {
        setCubeFileName(cubeT);
        cubeFile_ = std::make_shared<std::ofstream>(getCubeFileName()+".cube");
      }


      /**
       * @brief set density boolean 
       * 
      */
      void setDenEval(bool denEval) {
        denCube_ = denEval;
      }


      /**
       * @brief returns density boolean 
       * 
      */
      bool getDenEval() {
        return denCube_;
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

      // >>> Grid-related functions
      // see src/cubegen
      void calculateVoxelDimensions();
      std::vector<double> calcCenter();

      // >>> High-level functions
      void writeSummary(std::string fileSum);
      template <typename LocMatsT>
      void evalCube(std::string filePref,std::shared_ptr<cqmatrix::PauliSpinorMatrices<LocMatsT>> );

      // >>> Property evaluation functions
      template <typename LocMatsT>
      void evalDenCube(LocMatsT*);

  };

};

