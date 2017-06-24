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

#include <cxxapi/input.hpp>
#include <cxxapi/procedural.hpp>
#include <cxxapi/boilerplate.hpp>

// Workaround to intel compiler bug that breaks libint ERIs
#if !LIBINT2_CONSTEXPR_STATICS
  #include <libint2/statics_definition.h>
#endif

#include <fstream>

#include <cerr.hpp>

using namespace ChronusQ;

int main(int argc, char *argv[]) {

  ChronusQ::initialize();

  int rank = MPIRank();
  int size = MPISize();

  std::string inFileName, outFileName;
  std::string rstFileName, scrFileName;

  std::string oldRstFileName;

  bool rstExist;

  // Parse command line options
  if(argc < 2) { // No Options

    CErr("No Command Line Arguments Given",std::cout);
  //inFileName = "distfromroot_reduced.inp";
  //std::vector<std::string> tokens;
  //split(tokens,inFileName,".");

  //outFileName = tokens[0] + ".out";
  //rstFileName = tokens[0] + ".bin";

  } else if(argc == 2) { // Just Input File

    inFileName = argv[1];

  } else { // Variable Argc

    int c;
    while((c = getopt(argc,argv,"i:o:b:z:s:")) != -1) {
      switch(c) {
        case('i'):
          inFileName = optarg;
          break;
        case('o'):
          outFileName = optarg;
          break;
        case('b'):
          rstFileName = optarg;
          break;
        case('z'):
          oldRstFileName = optarg;
          break;
        case('s'):
          scrFileName = optarg;
          break;
        default:
          abort();
      };
    };

  }

  if( inFileName.empty() ) CErr("No Input File Specified!");

  std::vector<std::string> tokens;
  split(tokens,inFileName,".");

  if( outFileName.empty() ) outFileName = tokens[0] + ".out";
  if( rstFileName.empty() ){
    rstFileName = tokens[0] + ".bin";
    rstExist = false;
  } else {
    rstExist = true;
  }

  if( rstFileName == oldRstFileName ) 
    CErr("Old (-z) and current (-b) rstFile cannot have the same name!");

  if (not oldRstFileName.empty() and rank == 0) {

    rstExist = true;

    std::ifstream oldRstFile( oldRstFileName.c_str(), std::ios::binary );

    if( !oldRstFile ){
      CErr("Cannot find old rstFile!");
    }

    // Remove destination file if it exists
    if(std::ifstream(rstFileName.c_str(), std::ios::binary))
      std::remove(rstFileName.c_str() );

    std::ofstream rstFile(rstFileName.c_str(), std::ios::binary);

    // Copy over "old" rst file into new rst file
    std::cout << "  * Copying " << oldRstFileName << "  -->  "
              << rstFileName << std::endl;
    rstFile << oldRstFile.rdbuf();
  }

  RunChronusQ(inFileName,outFileName,rstFileName,scrFileName,rstExist);

  ChronusQ::finalize();

  return 0;
}

