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
#define DASSCF_TEST_REF TEST_ROOT "/dasscf/reference/"

using namespace ChronusQ;

inline void CQNORMALDASSCF( std::string in, std::string ref ) {


#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    DASSCF_TEST_REF + ref, "", false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    "", false);

#endif

};

inline void CQBINDASSCF( std::string in, std::string ref ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    DASSCF_TEST_REF + ref, "", true);

#else

  std::ifstream  src(DASSCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    "", true);

#endif

};

inline void CQSCRDASSCF( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    DASSCF_TEST_REF + ref, DASSCF_TEST_REF + scr, false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    DASSCF_TEST_REF + scr, false);

#endif

};

inline void CQBINSCRDASSCF( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    DASSCF_TEST_REF + ref, DASSCF_TEST_REF + scr, true);

#else

  std::ifstream  src(DASSCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    DASSCF_TEST_REF + scr, true);

#endif

};

inline void CQDASSCFTEST( std::string in, std::string ref,
  bool readBin      = false,
  std::string scr   = "",
  double tol        = 1e-8,
  bool checkOctLen  = false,
  bool checkQuadLen = false,
  bool checkDipLen  = false,
  bool checkOsc     = false,
  bool checkEne     = true ) {

  if( !readBin and scr=="" ) CQNORMALDASSCF(in,ref);
  else if( readBin and scr=="" ) CQBINDASSCF(in,ref);
  else if( !readBin and scr!="" ) CQSCRDASSCF(in,ref,scr);
  else CQBINSCRDASSCF(in,ref,scr);

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(DASSCF_TEST_REF + ref,  true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
  
  size_t xNS, yNS;
  std::cout << " * PERFORMING MSWFN ENERGY CHECK " << std::endl;
  std::cout << "DASSCF_TEST_REF=" << DASSCF_TEST_REF <<std::endl;

  refFile.readData("POSTHF/NSTATES", &xNS);
  resFile.readData("POSTHF/NSTATES", &yNS);

  EXPECT_EQ( xNS, yNS ) << "NUMBER OF STATES TEST FAILED ";

  std::vector<double> xStateEnergy(xNS);
  std::vector<double> yStateEnergy(yNS);
  std::vector<std::array<double,3>> xDummy3(xNS), yDummy3(yNS);
  std::vector<std::array<std::array<double,3>,3>> xDummy33(xNS), yDummy33(yNS);
  std::vector<std::array<std::array<std::array<double,3>,3>,3>> xDummy333(xNS), yDummy333(yNS);
  
  /* Check Energy */ 
  if( checkEne ) {

    refFile.readData("POSTHF/STATE_ENERGY", &xStateEnergy[0]);
    resFile.readData("POSTHF/STATE_ENERGY", &yStateEnergy[0]);

    for(size_t i = 0; i < xNS; i++)
      EXPECT_NEAR(xStateEnergy[i], yStateEnergy[i], tol ) << 
        "ENERGY TEST FAILED ISTATE = " << i; 
  }
  
  /* Check Oscillator Strength */
  if( checkOsc ) {

    std::cout << " * PERFORMING DASSCF OSCILLATOR STRENGTH CHECK " << std::endl;

    auto oscDim     = resFile.getDims("POSTHF/OSC_STR_ZERO_ORDER");
    ASSERT_EQ(oscDim.size(),2);
    std::vector<double> xDummy, yDummy;
    xDummy.resize(oscDim[0] * oscDim[1]);
    yDummy.resize(oscDim[0] * oscDim[1]);
    refFile.readData("POSTHF/OSC_STR_ZERO_ORDER", &xDummy[0]);
    resFile.readData("POSTHF/OSC_STR_ZERO_ORDER", &yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++)
    for(auto j = 0; j < oscDim[1]; j++)
      EXPECT_NEAR(yDummy[i*oscDim[1]+j], xDummy[i*oscDim[1]+j], tol ) <<
        "OSCILLATOR STRENGTH TEST FAILED STATE 1 = " << i
                                     << "STATE 2 = " << j;

  }

  /* Check Multipoles */
  if( checkDipLen ) {

    std::cout << " * PERFORMING DASSCF DIPOLE (LEN) CHECK " << std::endl;
    refFile.readData("MCWFN/LEN_ELECTRIC_DIPOLE",&xDummy3[0][0]);
    resFile.readData("MCWFN/LEN_ELECTRIC_DIPOLE",&yDummy3[0][0]);
    for( size_t iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
        EXPECT_NEAR(yDummy3[iSt][i], xDummy3[iSt][i], tol ) <<
          "DIPOLE TEST FAILED IXYZ = " << i << " FOR STATE = " << iSt;

  }


  if( checkQuadLen ) {

    std::cout << " * PERFORMING DASSCF QUADRUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("MCWFN/LEN_ELECTRIC_QUADRUPOLE",&xDummy33[0][0][0]);
    resFile.readData("MCWFN/LEN_ELECTRIC_QUADRUPOLE",&yDummy33[0][0][0]);
    for( size_t iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
      for(auto j = 0; j < 3; j++)
        EXPECT_NEAR(yDummy33[iSt][i][j], xDummy33[iSt][i][j],  tol) <<
          "QUADRUPOLE TEST FAILED IXYZ = " << i
                             << " JXYZ = " << j
                             << " FOR STATE = " << iSt;

  }

  if( checkOctLen ) {

    std::cout << " * PERFORMING DASSCF OCTUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("MCWFN/LEN_ELECTRIC_OCTUPOLE",&xDummy333[0][0][0][0]);
    resFile.readData("MCWFN/LEN_ELECTRIC_OCTUPOLE",&yDummy333[0][0][0][0]);
    for( size_t iSt = 0; iSt < xNS; iSt++ )
      for(auto i = 0; i < 3; i++)
      for(auto j = 0; j < 3; j++)
      for(auto k = 0; k < 3; k++)
        EXPECT_NEAR(yDummy333[iSt][i][j][k],  xDummy333[iSt][i][j][k],  tol) <<
          "OCTUPOLE TEST FAILED IXYZ = " << i
                                         << " JXYZ = " << j
                                         << " KXYZ = " << k
                                         << " FOR STATE = " << iSt;

  }

//  /* Check Spin */
//
//  if( checkSEXP ) {
//
//    std::cout << " * PERFORMING SCF <S> CHECK " << std::endl;
//
//    refFile.readData("SCF/S_EXPECT",&xDummy3[0]);
//    resFile.readData("SCF/S_EXPECT",&yDummy3[0]);
//    for(auto i = 0; i < 3; i++)
//      EXPECT_NEAR(yDummy3[i], xDummy3[i], tol ) << 
//        "<S> TEST FAILED IXYZ = " << i; 
//
//  }
//  
//  if( checkSSq ) {
//
//    std::cout << " * PERFORMING SCF <S^2> CHECK " << std::endl;
//
//    refFile.readData("SCF/S_SQUARED",&xDummy);
//    resFile.readData("SCF/S_SQUARED",&yDummy);
//    EXPECT_NEAR(yDummy,  xDummy, tol) << "<S^2> TEST FAILED " ;
//
//  }

#endif

}




