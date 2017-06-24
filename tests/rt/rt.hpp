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
#define RT_TEST_REF TEST_ROOT "/rt/reference/"

using namespace ChronusQ;


#ifdef _CQ_GENERATE_TESTS

// Run CQ job and write reference files
// *** WARNING: This will overwrite existing reference files ***
// ***                       USE AOR                         ***
#define CQRTTEST( in, ref ) \
  RunChronusQ(TEST_ROOT #in ".inp","STDOUT", \
    RT_TEST_REF #ref,TEST_OUT #in ".scr");

#define CQRTRESTARTTEST( midi, midr, in, ref ) \
  CQRTTEST( midi, midr.temp ) \
  CQRTTEST( in ## _ref, ref ) \
  \
  SafeFile tempFile(RT_TEST_REF #midr ".temp", true);\
  std::remove(RT_TEST_REF #midr);\
  SafeFile midFile(RT_TEST_REF #midr, false);\
  midFile.createFile();\
  \
  auto timeDims = tempFile.getDims("/RTNEW/TIME");\
  auto moDims = tempFile.getDims("/SCF/MO1");\
  auto denDims = tempFile.getDims("/SCF/1PDM_SCALAR");\
  \
  bool exists_(true);\
  std::string fName_(RT_TEST_REF #midr ".temp");\
  OpenDataSet(file, obj, "/SCF/MO1");\
  H5::DataType moType = obj.getDataType();\
  \
  if ( moType.getSize() > sizeof(double) ) { \
    std::vector<dcomplex> mo(moDims[0]*moDims[1], dcomplex(0.));\
    tempFile.readData("/SCF/MO1", mo.data());\
    midFile.safeWriteData("/SCF/MO1", mo.data(), moDims);\
    \
    try { \
      tempFile.readData("/SCF/MO2", mo.data());\
      midFile.safeWriteData("/SCF/MO2", mo.data(), moDims);\
    }\
    catch(...) { }\
  }\
  else { \
    std::vector<double> mo(moDims[0]*moDims[1], 0.);\
    tempFile.readData("/SCF/MO1", mo.data());\
    midFile.safeWriteData("/SCF/MO1", mo.data(), moDims);\
    \
    try { \
      tempFile.readData("/SCF/MO2", mo.data());\
      midFile.safeWriteData("/SCF/MO2", mo.data(), moDims);\
    }\
    catch(...) { }\
  }\
  \
  size_t savHash;\
  tempFile.readData("/SCF/FIELD_TYPE", &savHash);\
  midFile.safeWriteData("/SCF/FIELD_TYPE", &savHash, {1});\
  \
  size_t newTimeD = timeDims[0]*2 - 1;\
  std::vector<double> tArr(newTimeD, 0.);\
  std::vector<double> tArr3(newTimeD*3, 0.);\
  std::vector<dcomplex> tdDen(denDims[0]*denDims[1], 0.);\
  \
  tempFile.readData("/RTNEW/TIME", tArr.data());\
  midFile.safeWriteData("/RTNEW/TIME", tArr.data(), {newTimeD});\
  \
  tempFile.readData("/RTNEW/ENERGY", tArr.data());\
  midFile.safeWriteData("/RTNEW/ENERGY", tArr.data(), {newTimeD});\
  \
  tempFile.readData("/RTNEW/LEN_ELEC_DIPOLE", tArr3.data());\
  midFile.safeWriteData("/RTNEW/LEN_ELEC_DIPOLE", tArr3.data(), {newTimeD,3});\
  \
  tempFile.readData("/RTNEW/LEN_ELEC_DIPOLE_FIELD", tArr3.data());\
  midFile.safeWriteData("/RTNEW/LEN_ELEC_DIPOLE_FIELD", tArr3.data(), {newTimeD,3});\
  \
  tempFile.readData("/RTNEW/TD_1PDM_ORTHO0", tArr3.data());\
  midFile.safeWriteData("/RTNEW/TD_1PDM_ORTHO0", tArr3.data(), {newTimeD,denDims,denDims});\
  \
  std::vector<std::string> decomp{"ORTHO0", "ORTHO1"};\
  for ( auto &str : decomp ) { \
    try { \
      tempFile.readData("/RT/TD_1PDM_" + str, tdDen.data());\
      midFile.safeWriteData("/RT/TD_1PDM_" + str, tdDen.data(), moDims);\
      \
    } catch(...) { }\
  }\
  \
  std::remove(RT_TEST_REF #midr ".temp" );

#else

// HTG RT test
static void CQRTTEST(std::string in, std::string ref, 
   double tol = 1e-8,
   bool readBin = false ){ 
  
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",TEST_OUT + in + ".bin","",readBin);
  
  SafeFile refFile(RT_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
  \
  std::vector<double> xDummy, yDummy;\
  \
  auto energyDim1 = resFile.getDims("/RTNEW/ENERGY");\
  auto energyDim2 = refFile.getDims("/RTNEW/ENERGY");\
  ASSERT_EQ( energyDim1.size(), 1 );\
  ASSERT_EQ( energyDim2.size(), 1 );\
  ASSERT_EQ( energyDim1[0], energyDim2[0] );\
  \
  auto dipoleDim1 = resFile.getDims("/RTNEW/LEN_ELEC_DIPOLE");\
  auto dipoleDim2 = refFile.getDims("/RTNEW/LEN_ELEC_DIPOLE");\
  ASSERT_EQ( dipoleDim1.size(), 1 );\
  ASSERT_EQ( dipoleDim2.size(), 1 );\
  ASSERT_EQ( dipoleDim1[0], dipoleDim2[0] );\
  \
  std::cout << "Checking RT energy" << std::endl;\
  xDummy.resize(energyDim1[0]); yDummy.resize(energyDim1[0]);\
  resFile.readData("/RTNEW/ENERGY",&xDummy[0]);\
  refFile.readData("/RTNEW/ENERGY",&yDummy[0]);\
  \
  for(auto i = 0; i < energyDim1[0]; i++){\ 
    EXPECT_NEAR(xDummy[i], yDummy[i], tol);\
  }\
\
  std::cout << "Checking RT Dipole" << std::endl;\
  xDummy.resize(dipoleDim1[0]); yDummy.resize(dipoleDim2[0]);\
  resFile.readData("/RTNEW/LEN_ELEC_DIPOLE",&xDummy[0]);\
  refFile.readData("/RTNEW/LEN_ELEC_DIPOLE",&yDummy[0]);\
  \
  for(auto i = 0; i < 3 * energyDim1[0]; i++) {\
    EXPECT_NEAR(xDummy[i], yDummy[i], tol); \   
  }\
}

static void CQRTRESTARTTEST( std::string midr, std::string in, std::string ref, double tol = 1e-8 ) {

  std::ifstream oldFile( RT_TEST_REF + midr, std::ios::binary );\
  std::ofstream newFile( TEST_OUT + in + ".bin" );\
  newFile << oldFile.rdbuf();
  newFile.flush();

  CQRTTEST(in, ref, tol, true);
}
#endif

