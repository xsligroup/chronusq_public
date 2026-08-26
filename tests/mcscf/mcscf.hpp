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
//#include <util/mpi.hpp>

// Directory containing reference files
#define MCSCF_TEST_REF TEST_ROOT "mcscf/reference/"

using namespace ChronusQ;

inline void CQNORMALMCSCF( std::string in, std::string ref ) {


#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MCSCF_TEST_REF + ref, "", false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    "", false);

#endif

};

inline void CQBINMCSCF( std::string in, std::string ref ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MCSCF_TEST_REF + ref, "", true);

#else

  std::ifstream  src(MCSCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(CQTestOut(in,".bin"), std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    "", true);

#endif

};

inline void CQSCRMCSCF( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MCSCF_TEST_REF + ref, MCSCF_TEST_REF + scr, false);

#else

  // Run off a private copy of the reference scratch file: the job
  // writes to its scratch file, and tests run side by side.
  std::ifstream  scrSrc(MCSCF_TEST_REF + scr, std::ios::binary);
  std::ofstream  scrDst(CQTestOut(in,".scr.bin"), std::ios::binary);
  scrDst << scrSrc.rdbuf();
  scrDst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    CQTestOut(in,".scr.bin"), false);

#endif

};

inline void CQBINSCRMCSCF( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MCSCF_TEST_REF + ref, MCSCF_TEST_REF + scr, true);

#else

  std::ifstream  src(MCSCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(CQTestOut(in,".bin"), std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  // Run off a private copy of the reference scratch file: the job
  // writes to its scratch file, and tests run side by side.
  std::ifstream  scrSrc(MCSCF_TEST_REF + scr, std::ios::binary);
  std::ofstream  scrDst(CQTestOut(in,".scr.bin"), std::ios::binary);
  scrDst << scrSrc.rdbuf();
  scrDst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    CQTestOut(in,".scr.bin"), true);

#endif

};

inline void CQMCSCFTEST( std::string in, std::string ref,
  bool readBin      = false,
  std::string scr   = "",
  double tol        = 1e-8,
  bool checkOctLen  = false,
  bool checkQuadLen = false,
  bool checkDipLen  = false,
  bool checkOsc     = false,
  bool checkEne     = true,
  bool checkSEXP    = false,
  bool checkSSq     = false,
  bool checkLExp    = false,
  bool checkLSq     = false,
  bool checkSL      = false,
  bool checkJExp    = false,
  bool checkJSq     = false,
  bool checkOrbIdx  = false ) {

  if( !readBin and scr=="" ) CQNORMALMCSCF(in,ref);
  else if( readBin and scr=="" ) CQBINMCSCF(in,ref);
  else if( !readBin and scr!="" ) CQSCRMCSCF(in,ref,scr);
  else CQBINSCRMCSCF(in,ref,scr);

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(MCSCF_TEST_REF + ref,  true,  true);
  SafeFile resFile(CQTestOut(in,".bin"),true);
  auto datasetPath = [](SafeFile &file, const std::string &name) {
    const std::string postHFPath = "POSTHF/" + name;
    if (file.exists(postHFPath)) return postHFPath;
    return "MCWFN/" + name;
  };
  
  double xNS, yNS;
  std::cout << " * PERFORMING MCWFN ENERGY CHECK " << std::endl;
  std::cout << "MCSCF_TEST_REF=" << MCSCF_TEST_REF <<std::endl;
  
  refFile.readData(datasetPath(refFile, "NSTATES"), &xNS);
  resFile.readData(datasetPath(resFile, "NSTATES"), &yNS);

  EXPECT_NEAR( xNS, yNS, tol ) << "NUMBER OF STATES TEST FAILED ";

  std::vector<double> xStateEnergy(xNS);
  std::vector<double> yStateEnergy(yNS);
  std::vector<std::array<double,3>> xDummy3(xNS), yDummy3(yNS);
  std::vector<std::array<std::array<double,3>,3>> xDummy33(xNS), yDummy33(yNS);
  std::vector<std::array<std::array<std::array<double,3>,3>,3>> xDummy333(xNS), yDummy333(yNS);
  
  /* Check Energy */ 
  if( checkEne ) {

    refFile.readData(datasetPath(refFile, "STATE_ENERGY"), &xStateEnergy[0]);
    resFile.readData(datasetPath(resFile, "STATE_ENERGY"), &yStateEnergy[0]);

    for(auto i = 0; i < xNS; i++)
      EXPECT_NEAR(xStateEnergy[i], yStateEnergy[i], tol ) << 
        "ENERGY TEST FAILED ISTATE = " << i; 
  }
  
  /* Check Oscillator Strength */
  if( checkOsc ) {

    std::cout << " * PERFORMING MCSCF OSCILLATOR STRENGTH CHECK " << std::endl;

    const auto oscPath = datasetPath(resFile, "OSC_STR");
    auto oscDim     = resFile.getDims(oscPath);
    ASSERT_EQ(oscDim.size(),2);
    std::vector<double> xDummy, yDummy;
    xDummy.resize(oscDim[0] * oscDim[1]);
    yDummy.resize(oscDim[0] * oscDim[1]);
    refFile.readData(datasetPath(refFile, "OSC_STR"), &xDummy[0]);
    resFile.readData(oscPath, &yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++)
    for(auto j = 0; j < oscDim[1]; j++)
      EXPECT_NEAR(yDummy[i*oscDim[1]+j], xDummy[i*oscDim[1]+j], tol ) <<
        "OSCILLATOR STRENGTH TEST FAILED STATE 1 = " << i
                                     << "STATE 2 = " << j;

  }

  /* Check Multipoles */
  if( checkDipLen ) {

    std::cout << " * PERFORMING MCSCF DIPOLE (LEN) CHECK " << std::endl;
    refFile.readData(datasetPath(refFile, "LEN_ELECTRIC_DIPOLE"),&xDummy3[0][0]);
    resFile.readData(datasetPath(resFile, "LEN_ELECTRIC_DIPOLE"),&yDummy3[0][0]);
    for( auto iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
        EXPECT_NEAR(yDummy3[iSt][i], xDummy3[iSt][i], tol ) <<
          "DIPOLE TEST FAILED IXYZ = " << i << " FOR STATE = " << iSt;

  }


  if( checkQuadLen ) {

    std::cout << " * PERFORMING MCSCF QUADRUPOLE (LEN) CHECK " << std::endl;

    refFile.readData(datasetPath(refFile, "LEN_ELECTRIC_QUADRUPOLE"),&xDummy33[0][0][0]);
    resFile.readData(datasetPath(resFile, "LEN_ELECTRIC_QUADRUPOLE"),&yDummy33[0][0][0]);
    for( auto iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
      for(auto j = 0; j < 3; j++)
        EXPECT_NEAR(yDummy33[iSt][i][j], xDummy33[iSt][i][j],  tol) <<
          "QUADRUPOLE TEST FAILED IXYZ = " << i
                             << " JXYZ = " << j
                             << " FOR STATE = " << iSt;

  }

  if( checkOctLen ) {

    std::cout << " * PERFORMING MCSCF OCTUPOLE (LEN) CHECK " << std::endl;

    refFile.readData(datasetPath(refFile, "LEN_ELECTRIC_OCTUPOLE"),&xDummy333[0][0][0][0]);
    resFile.readData(datasetPath(resFile, "LEN_ELECTRIC_OCTUPOLE"),&yDummy333[0][0][0][0]);
    for( auto iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
      for(auto j = 0; j < 3; j++)
      for(auto k = 0; k < 3; k++)
        EXPECT_NEAR(yDummy333[iSt][i][j][k],  xDummy333[iSt][i][j][k],  tol) <<
          "OCTUPOLE TEST FAILED IXYZ = " << i
                                         << " JXYZ = " << j
                                         << " KXYZ = " << k
                                         << " FOR STATE = " << iSt;

  }

  if( checkOrbIdx ) {

    std::cout << " * PERFORMING POSTHF ORBITAL INDEX CHECK " << std::endl;

    const auto orbIndexPath = datasetPath(refFile, "ORB_INDEX");
    const auto orbIndexDims = refFile.getDims(orbIndexPath);
    ASSERT_EQ(orbIndexDims.size(), 1);
    std::vector<int> xOrbIndex(orbIndexDims[0]), yOrbIndex(orbIndexDims[0]);
    refFile.readData(orbIndexPath, xOrbIndex.data());
    resFile.readData(datasetPath(resFile, "ORB_INDEX"), yOrbIndex.data());
    ASSERT_EQ(yOrbIndex.size(), xOrbIndex.size());
    for(size_t iOrb = 0; iOrb < xOrbIndex.size(); iOrb++)
      EXPECT_EQ(yOrbIndex[iOrb], xOrbIndex[iOrb])
        << "ORBITAL INDEX TEST FAILED INDEX = " << iOrb;

  }

  if( checkSEXP ) {

    std::cout << " * PERFORMING POSTHF <S> CHECK " << std::endl;

    std::vector<std::array<double,3>> xSExpect(xNS), ySExpect(yNS);
    refFile.readData(datasetPath(refFile, "S_EXPECT"), &xSExpect[0][0]);
    resFile.readData(datasetPath(resFile, "S_EXPECT"), &ySExpect[0][0]);
    for(size_t iState = 0; iState < xNS; iState++)
      for(auto i = 0; i < 3; i++)
        EXPECT_NEAR(ySExpect[iState][i], xSExpect[iState][i], tol)
          << "<S> TEST FAILED STATE = " << iState << " IXYZ = " << i;

  }

  if( checkSSq ) {

    std::cout << " * PERFORMING POSTHF <S^2> CHECK " << std::endl;

    std::vector<double> xSSq(xNS), ySSq(yNS);
    refFile.readData(datasetPath(refFile, "S_SQUARED"), xSSq.data());
    resFile.readData(datasetPath(resFile, "S_SQUARED"), ySSq.data());
    for(size_t iState = 0; iState < xNS; iState++)
      EXPECT_NEAR(ySSq[iState], xSSq[iState], tol)
        << "<S^2> TEST FAILED STATE = " << iState;

  }

  if( checkLExp ) {

    std::cout << " * PERFORMING POSTHF <L> CHECK " << std::endl;

    std::vector<std::array<double,3>> xLExpect(xNS), yLExpect(yNS);
    refFile.readData(datasetPath(refFile, "L_EXPECT"), &xLExpect[0][0]);
    resFile.readData(datasetPath(resFile, "L_EXPECT"), &yLExpect[0][0]);
    for(size_t iState = 0; iState < xNS; iState++)
      for(auto i = 0; i < 3; i++)
        EXPECT_NEAR(yLExpect[iState][i], xLExpect[iState][i], tol)
          << "<L> TEST FAILED STATE = " << iState << " IXYZ = " << i;

  }

  if( checkLSq ) {

    std::cout << " * PERFORMING POSTHF <L^2> CHECK " << std::endl;

    std::vector<double> xLSq(xNS), yLSq(yNS);
    refFile.readData(datasetPath(refFile, "L_SQUARED"), xLSq.data());
    resFile.readData(datasetPath(resFile, "L_SQUARED"), yLSq.data());
    for(size_t iState = 0; iState < xNS; iState++)
      EXPECT_NEAR(yLSq[iState], xLSq[iState], tol)
        << "<L^2> TEST FAILED STATE = " << iState;

  }

  if( checkSL ) {

    std::cout << " * PERFORMING POSTHF <SL> CHECK " << std::endl;

    std::vector<double> xSL(xNS), ySL(yNS);
    refFile.readData(datasetPath(refFile, "SL"), xSL.data());
    resFile.readData(datasetPath(resFile, "SL"), ySL.data());
    for(size_t iState = 0; iState < xNS; iState++)
      EXPECT_NEAR(ySL[iState], xSL[iState], tol)
        << "<S.L> TEST FAILED STATE = " << iState;

  }

  if( checkJExp ) {

    std::cout << " * PERFORMING POSTHF <J> CHECK " << std::endl;

    std::vector<std::array<double,3>> xJExpect(xNS), yJExpect(yNS);
    refFile.readData(datasetPath(refFile, "J_EXPECT"), &xJExpect[0][0]);
    resFile.readData(datasetPath(resFile, "J_EXPECT"), &yJExpect[0][0]);
    for(size_t iState = 0; iState < xNS; iState++)
      for(auto i = 0; i < 3; i++)
        EXPECT_NEAR(yJExpect[iState][i], xJExpect[iState][i], tol)
          << "<J> TEST FAILED STATE = " << iState << " IXYZ = " << i;

  }

  if( checkJSq ) {

    std::cout << " * PERFORMING POSTHF <J^2> CHECK " << std::endl;

    std::vector<double> xJSq(xNS), yJSq(yNS);
    refFile.readData(datasetPath(refFile, "J_SQUARED"), xJSq.data());
    resFile.readData(datasetPath(resFile, "J_SQUARED"), yJSq.data());
    for(size_t iState = 0; iState < xNS; iState++)
      EXPECT_NEAR(yJSq[iState], xJSq[iState], tol)
        << "<J^2> TEST FAILED STATE = " << iState;

  }

#endif

}




