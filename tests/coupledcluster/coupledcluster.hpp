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

#define CC_TEST_REF TEST_ROOT "/coupledcluster/reference/"

using namespace ChronusQ;



static void CQCCTEST( std::string in, std::string ref,
                      bool checkExcitedEnergy = false,
                      bool checkOsc = false,
                      double etol = 1e-7,
                      double osctol = 1e-5
                    ) {

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    CC_TEST_REF + ref, TEST_OUT + in + ".scr");

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin",TEST_OUT + in + ".scr");

  if(MPIRank(MPI_COMM_WORLD) != 0) return;

          
  SafeFile refFile(CC_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  double testE, refE;
  std::cout << " * PERFORMING CC ENERGY CHECK " << std::endl;
  std::cout << "CC_TEST_REF=" << CC_TEST_REF <<std::endl;

  resFile.readData("/CC/CORRELATION_ENERGY", &testE);
  refFile.readData("/CC/CORRELATION_ENERGY", &refE);

  EXPECT_NEAR( testE, refE, etol ) << "CC CORRELATION ENERGY TEST FAILED";

  std::vector<double> xDummy, yDummy;

  if (checkExcitedEnergy){
    // Check eigenvalues
    std::cerr << "PERFORMING EXCITATION ENERGY TEST\n";

    auto evDim     = resFile.getDims("/CC/EXCITATION_ENERGIES");
    auto evDim_ref = refFile.getDims("/CC/EXCITATION_ENERGIES");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(evDim[0]); yDummy.resize(evDim_ref[0]);


    resFile.readData("/CC/EXCITATION_ENERGIES",&xDummy[0]);
    refFile.readData("/CC/EXCITATION_ENERGIES",&yDummy[0]);

    for(auto i = 0; i < evDim[0]; i++)
      EXPECT_NEAR( xDummy[i], yDummy[i], etol ) << "EXCITATION ENERGY TEST FAILED IN STATE = " << i;
  }

  if (checkOsc){
    // Check Osc Strength
    std::cerr << "PERFORMING OSC STRENGTH TEST\n";

    auto oscDim     = resFile.getDims("/CC/OSCILLATOR_STRENGTHS");
    auto oscDim_ref = refFile.getDims("/CC/OSCILLATOR_STRENGTHS");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(oscDim[0]); yDummy.resize(oscDim_ref[0]);

    resFile.readData("/CC/OSCILLATOR_STRENGTHS",&xDummy[0]);
    refFile.readData("/CC/OSCILLATOR_STRENGTHS",&yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++)
      EXPECT_NEAR( xDummy[i], yDummy[i], osctol ) << "OSC STRENGTH TEST FAILED IN STATE = " << i;
  }


#endif


};
