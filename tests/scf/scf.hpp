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

// Directory containing reference files
#define SCF_TEST_REF TEST_ROOT "/scf/reference/"

using namespace ChronusQ;

static void CQNORMALSCF( std::string in, std::string ref ) {


#ifdef _CQ_GENERATE_TESTS

  MPI_Barrier(MPI_COMM_WORLD);

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    SCF_TEST_REF + ref, "", false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    "", false);

#endif

};

static void CQBINSCF( std::string in, std::string ref ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    SCF_TEST_REF + ref, "", true);

#else

  std::ifstream  src(SCF_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    "", true);

#endif

};

static void CQSCRSCF( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    SCF_TEST_REF + ref, SCF_TEST_REF + scr, false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    TEST_OUT + in + ".bin",
    SCF_TEST_REF + scr, false);

#endif

};
 
static void CQSCFTEST( std::string in, std::string ref,
  double tol        = 1e-6,
  bool checkSEXP    = true,
  bool checkSSq     = true,
  bool checkOctLen  = true,
  bool checkQuadLen = true,
  bool checkDipLen  = true,
  bool checkEne     = true,
  bool readBin      = false,
  std::string scr  = "no",
  bool looserPropretyThreshold = false ) {

  MPI_Barrier(MPI_COMM_WORLD);

  if( !readBin and scr=="no" ) CQNORMALSCF(in,ref);
  else if( readBin ) CQBINSCF(in,ref);
  else CQSCRSCF(in,ref,scr);

  MPI_Barrier(MPI_COMM_WORLD);
  if(MPIRank(MPI_COMM_WORLD) != 0) return;


#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(SCF_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  double xDummy, yDummy;
  std::array<double,3> xDummy3, yDummy3;
  std::array<std::array<double,3>,3> xDummy33, yDummy33;
  std::array<std::array<std::array<double,3>,3>,3> xDummy333, yDummy333;

  /* Check Energy */
  if( checkEne ) {

    std::cout << " * PERFORMING SCF ENERGY CHECK " << std::endl;

    refFile.readData("SCF/TOTAL_ENERGY",&xDummy);
    resFile.readData("SCF/TOTAL_ENERGY",&yDummy);

    EXPECT_NEAR( xDummy, yDummy, tol ) << "ENERGY TEST FAILED ";

  }

  double property_tol = looserPropretyThreshold ? tol*100.0 : tol ;

  /* Check Multipoles */

  if( checkDipLen ) {

    std::cout << " * PERFORMING SCF DIPOLE (LEN) CHECK " << std::endl;

    refFile.readData("SCF/LEN_ELECTRIC_DIPOLE",&xDummy3[0]);
    resFile.readData("SCF/LEN_ELECTRIC_DIPOLE",&yDummy3[0]);
    for(auto i = 0; i < 3; i++)
      EXPECT_NEAR(yDummy3[i], xDummy3[i], property_tol ) <<
        "DIPOLE TEST FAILED IXYZ = " << i;


  }


  if( checkQuadLen ) {

    std::cout << " * PERFORMING SCF QUADRUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&xDummy33[0][0]);
    resFile.readData("SCF/LEN_ELECTRIC_QUADRUPOLE",&yDummy33[0][0]);
    for(auto i = 0; i < 3; i++)
    for(auto j = 0; j < 3; j++)
      EXPECT_NEAR(yDummy33[i][j], xDummy33[i][j],  property_tol) <<
        "QUADRUPOLE TEST FAILED IXYZ = " << i
                           << " JXYZ = " << j;

  }

  if( checkOctLen ) {

    std::cout << " * PERFORMING SCF OCTUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&xDummy333[0][0][0]);
    resFile.readData("SCF/LEN_ELECTRIC_OCTUPOLE",&yDummy333[0][0][0]);
    for(auto i = 0; i < 3; i++)
    for(auto j = 0; j < 3; j++)
    for(auto k = 0; k < 3; k++)
      EXPECT_NEAR(yDummy333[i][j][k],  xDummy333[i][j][k],  property_tol) <<
        "OCTUPOLE TEST FAILED IXYZ = " << i
                                       << " JXYZ = " << j
                                       << " KXYZ = " << k;

  }

  /* Check Spin */

  if( checkSEXP ) {

    std::cout << " * PERFORMING SCF <S> CHECK " << std::endl;

    refFile.readData("SCF/S_EXPECT",&xDummy3[0]);
    resFile.readData("SCF/S_EXPECT",&yDummy3[0]);
    for(auto i = 0; i < 3; i++)
      EXPECT_NEAR(yDummy3[i], xDummy3[i], property_tol ) <<
        "<S> TEST FAILED IXYZ = " << i;

  }

  if( checkSSq ) {

    std::cout << " * PERFORMING SCF <S^2> CHECK " << std::endl;

    refFile.readData("SCF/S_SQUARED",&xDummy);
    resFile.readData("SCF/S_SQUARED",&yDummy);
    EXPECT_NEAR(yDummy,  xDummy, property_tol) << "<S^2> TEST FAILED " ;

  }

#endif

}


static void CQNEOSCFTEST( std::string in, std::string ref,
  double tol        = 1e-6,
  bool checkOctLen  = true,
  bool checkQuadLen = true,
  bool checkDipLen  = true,
  bool checkEne     = true,
  bool readBin      = false,
  std::string scr  = "no",
  bool looserPropretyThreshold = false,
  bool checkSubSSEne = true ) {

  MPI_Barrier(MPI_COMM_WORLD);

  if( !readBin and scr=="no" ) CQNORMALSCF(in,ref);
  else if( readBin ) CQBINSCF(in,ref);
  else CQSCRSCF(in,ref,scr);

  MPI_Barrier(MPI_COMM_WORLD);
  if(MPIRank(MPI_COMM_WORLD) != 0) return;

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(SCF_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  double xDummyE, yDummyE, xDummyP, yDummyP, xDummyT, yDummyT;
  std::array<double,3> xDummy3, yDummy3;
  std::array<std::array<double,3>,3> xDummy33, yDummy33;
  std::array<std::array<std::array<double,3>,3>,3> xDummy333, yDummy333;

  /* Check Energy */
  if( checkEne ) {

    std::cout << " * PERFORMING NEO-SCF ENERGY CHECK " << std::endl;

    if(checkSubSSEne){
    
      std::cout << "     * CHECKING ELECTRONIC ENERGY" << std::endl;
      refFile.readData("NEO/ELEC_ENERGY",&xDummyE);
      resFile.readData("NEO/ELEC_ENERGY",&yDummyE);
      EXPECT_NEAR( xDummyE, yDummyE, tol ) << "NEO ELEC-ENERGY TEST FAILED ";

      std::cout << "     * CHECKING PROTONIC ENERGY" << std::endl;
      refFile.readData("NEO/PROT_ENERGY",&xDummyP);
      resFile.readData("NEO/PROT_ENERGY",&yDummyP);
      EXPECT_NEAR( xDummyP, yDummyP, tol ) << "NEO PROT-ENERGY TEST FAILED ";

    }
    
    std::cout << "     * CHECKING TOTAL NEO ENERGY" << std::endl;
    refFile.readData("NEO/TOTAL_ENERGY",&xDummyT);
    resFile.readData("NEO/TOTAL_ENERGY",&yDummyT);
    EXPECT_NEAR( xDummyT, yDummyT, tol ) << "NEO TOTAL-ENERGY TEST FAILED ";

  }

  double property_tol = looserPropretyThreshold ? tol*100.0 : tol ;
  /* Check Multipoles */

  if( checkDipLen ) {

    std::cout << " * PERFORMING NEO-SCF DIPOLE (LEN) CHECK " << std::endl;

    refFile.readData("NEO/LEN_ELECTRIC_DIPOLE",&xDummy3[0]);
    resFile.readData("NEO/LEN_ELECTRIC_DIPOLE",&yDummy3[0]);
    for(auto i = 0; i < 3; i++)
      EXPECT_NEAR(yDummy3[i], xDummy3[i], property_tol) <<
        "DIPOLE TEST FAILED IXYZ = " << i;

  }


  if( checkQuadLen ) {

    std::cout << " * PERFORMING NEO-SCF QUADRUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("NEO/LEN_ELECTRIC_QUADRUPOLE",&xDummy33[0][0]);
    resFile.readData("NEO/LEN_ELECTRIC_QUADRUPOLE",&yDummy33[0][0]);
    for(auto i = 0; i < 3; i++)
    for(auto j = 0; j < 3; j++)
      EXPECT_NEAR(yDummy33[i][j], xDummy33[i][j],  property_tol) <<
        "QUADRUPOLE TEST FAILED IXYZ = " << i
                           << " JXYZ = " << j;

  }

  if( checkOctLen ) {

    std::cout << " * PERFORMING NEO-SCF OCTUPOLE (LEN) CHECK " << std::endl;

    refFile.readData("NEO/LEN_ELECTRIC_OCTUPOLE",&xDummy333[0][0][0]);
    resFile.readData("NEO/LEN_ELECTRIC_OCTUPOLE",&yDummy333[0][0][0]);
    for(auto i = 0; i < 3; i++)
    for(auto j = 0; j < 3; j++)
    for(auto k = 0; k < 3; k++)
      EXPECT_NEAR(yDummy333[i][j][k],  xDummy333[i][j][k],  property_tol) <<
        "OCTUPOLE TEST FAILED IXYZ = " << i
                                       << " JXYZ = " << j
                                       << " KXYZ = " << k;

  }

#endif

}



