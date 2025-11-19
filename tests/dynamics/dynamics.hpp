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

#include <ut.hpp>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>
#include <util/mpi.hpp>
#include <iostream>

#include <fstream>
#include <cstdio>


// Directory containing reference files
#define DYNAMICS_TEST_REF TEST_ROOT "/dynamics/reference/"

using namespace ChronusQ;

// HTG Dynamics test
static void CQDYNAMICSTEST(std::string in, std::string ref, 
   double tol = 1e-8,
   bool readBin = false ){ 
  
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",TEST_OUT + in + ".bin","",readBin);
  
  // print the reference file path
  std::cout << "Reference file: " << DYNAMICS_TEST_REF + ref << std::endl;
  std::cout << "Result file: " << TEST_OUT + in + ".bin" << std::endl;
  
  SafeFile refFile(DYNAMICS_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
  \
  std::vector<double> xDummy, yDummy;\
  \
  auto energyDim1 = resFile.getDims("/MD/ETOT");\
  auto energyDim2 = refFile.getDims("/MD/ETOT");\
  ASSERT_EQ( energyDim1.size(), 1 );\
  ASSERT_EQ( energyDim2.size(), 1 );\
  ASSERT_EQ( energyDim1[0], energyDim2[0] );\
  \
  auto forceDim1 = resFile.getDims("/MD/FORCES");\
  auto forceDim2 = refFile.getDims("/MD/FORCES");\
  ASSERT_EQ( forceDim1.size(), 1 );\
  ASSERT_EQ( forceDim2.size(), 1 );\
  ASSERT_EQ( forceDim1[0], forceDim2[0] );\
  \
  auto velocityDim1 = resFile.getDims("/MD/VELOCITY_FULLSTEP");\
  auto velocityDim2 = refFile.getDims("/MD/VELOCITY_FULLSTEP");\
  ASSERT_EQ( velocityDim1.size(), 1 );\
  ASSERT_EQ( velocityDim2.size(), 1 );\
  ASSERT_EQ( velocityDim1[0], velocityDim2[0] );\
  \
  auto trajectoryDim1 = resFile.getDims("/MD/TRAJECTORY");\
  auto trajectoryDim2 = refFile.getDims("/MD/TRAJECTORY");\
  ASSERT_EQ( trajectoryDim1.size(), 1 );\
  ASSERT_EQ( trajectoryDim2.size(), 1 );\
  ASSERT_EQ( trajectoryDim1[0], trajectoryDim2[0] );\
  \
  std::cout << "Checking Dynamics Total Energy" << std::endl;\
  xDummy.resize(energyDim1[0]); yDummy.resize(energyDim1[0]);\
  resFile.readData("/MD/ETOT",&xDummy[0]);\
  refFile.readData("/MD/ETOT",&yDummy[0]);\
  \
  for(auto i = 0; i < energyDim1[0]; i++){\
    EXPECT_NEAR(xDummy[i], yDummy[i], tol);\
  }\
  \
  std::cout << "Checking Dynamics Forces" << std::endl;\
  xDummy.resize(forceDim1[0]); yDummy.resize(forceDim2[0]);\
  resFile.readData("/MD/FORCES",&xDummy[0]);\
  refFile.readData("/MD/FORCES",&yDummy[0]);\
  \
  for(auto i = 0; i < forceDim1[0]; i++) {\
    EXPECT_NEAR(xDummy[i], yDummy[i], tol); \
  }\
  \
  std::cout << "Checking Dynamics Velocity" << std::endl;\
  xDummy.resize(velocityDim1[0]); yDummy.resize(velocityDim2[0]);\
  resFile.readData("/MD/VELOCITY_FULLSTEP",&xDummy[0]);\
  refFile.readData("/MD/VELOCITY_FULLSTEP",&yDummy[0]);\
  \
  for(auto i = 0; i < velocityDim1[0]; i++) {\
    EXPECT_NEAR(xDummy[i], yDummy[i], tol); \
  }\
  \
  std::cout << "Checking Dynamics Trajectory" << std::endl;\
  xDummy.resize(trajectoryDim1[0]); yDummy.resize(trajectoryDim2[0]);\
  resFile.readData("/MD/TRAJECTORY",&xDummy[0]);\
  refFile.readData("/MD/TRAJECTORY",&yDummy[0]);\
  \
  for(auto i = 0; i < trajectoryDim1[0]; i++) {\
    EXPECT_NEAR(xDummy[i], yDummy[i], tol); \
  }\
}

static void CQDYNAMICSRESTARTTEST( std::string midr, std::string in, std::string ref, double tol = 1e-8 ) {

  std::ifstream oldFile( DYNAMICS_TEST_REF + midr, std::ios::binary );\
  std::ofstream newFile( TEST_OUT + in + ".bin" );\
  newFile << oldFile.rdbuf();
  newFile.flush();

  CQDYNAMICSTEST(in, ref, tol, true);
}

