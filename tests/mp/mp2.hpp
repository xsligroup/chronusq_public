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

#pragma once

#include <ut.hpp>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>
//#include <util/mpi.hpp>

// Directory containing reference files
#define MP_TEST_REF TEST_ROOT "/mp/reference/"

using namespace ChronusQ;

inline void CQNORMALMP( std::string in, std::string ref ) {


#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MP_TEST_REF + ref, "", false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    "", false);

#endif

};

inline void CQBINMP( std::string in, std::string ref ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MP_TEST_REF + ref, "", true);

#else

  std::ifstream  src(MP_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(CQTestOut(in,".bin"), std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    "", true);

#endif

};

inline void CQSCRMP( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MP_TEST_REF + ref, MP_TEST_REF + scr, false);

#else

  // Run off a private copy of the reference scratch file: the job
  // writes to its scratch file, and tests run side by side.
  std::ifstream  scrSrc(MP_TEST_REF + scr, std::ios::binary);
  std::ofstream  scrDst(CQTestOut(in,".scr.bin"), std::ios::binary);
  scrDst << scrSrc.rdbuf();
  scrDst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    CQTestOut(in,".scr.bin"), false);

#endif

};

inline void CQBINSCRMP( std::string in, std::string ref, std::string scr ) {

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    MP_TEST_REF + ref, MP_TEST_REF + scr, true);

#else

  std::ifstream  src(MP_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(CQTestOut(in,".bin"), std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  // Run off a private copy of the reference scratch file: the job
  // writes to its scratch file, and tests run side by side.
  std::ifstream  scrSrc(MP_TEST_REF + scr, std::ios::binary);
  std::ofstream  scrDst(CQTestOut(in,".scr.bin"), std::ios::binary);
  scrDst << scrSrc.rdbuf();
  scrDst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT",
    CQTestOut(in,".bin"),
    CQTestOut(in,".scr.bin"), true);

#endif

};

inline void CQMPTEST( std::string in, std::string ref,
  bool readBin      = false,
  std::string scr   = "",
  double tol        = 1e-8,
  bool checkEne     = true ) {

  if( !readBin and scr=="" ) CQNORMALMP(in,ref);
  else if( readBin and scr=="" ) CQBINMP(in,ref);
  else if( !readBin and scr!="" ) CQSCRMP(in,ref,scr);
  else CQBINSCRMP(in,ref,scr);

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(MP_TEST_REF + ref,  true,  true);
  SafeFile resFile(CQTestOut(in,".bin"),true);
  
  std::cout << " * PERFORMING MP2 ENERGY CHECK " << std::endl;

  double xDummy, yDummy;
  
  /* Check Energy */ 
  if( checkEne ) {

    refFile.readData("MP2/MP2_ENERGY", &xDummy);
    resFile.readData("MP2/MP2_ENERGY", &yDummy);
    
    EXPECT_NEAR( xDummy, yDummy, tol ) << "ENERGY TEST FAILED ";

  }

#endif

}

