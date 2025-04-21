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

#define CC_TEST_REF TEST_ROOT "coupledcluster/reference/"

using namespace ChronusQ;



static void CQCCTEST( std::string in, std::string ref, std::string restart_bin = "",
                      bool checkReferenceEnergy = false,
                      bool checkExcitedEnergy = false,
                      bool checkOsc = false,
                      bool checkTriplesCorrection = false,
                      double etol = 1e-7,
                      double osctol = 1e-5
                    ) {

  MPI_Barrier(MPI_COMM_WORLD);

  bool rstExist = false;
  if ( ! restart_bin.empty() ) {
    std::ifstream  src(CC_TEST_REF + restart_bin, std::ios::binary);
    std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
    dst << src.rdbuf();
    dst.flush();
    rstExist = true;
  }

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    CC_TEST_REF + ref, TEST_OUT + in + ".scr", rstExist);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin",TEST_OUT + in + ".scr", rstExist);

  if(MPIRank(MPI_COMM_WORLD) != 0) return;

          
  SafeFile refFile(CC_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);

  dcomplex testE, refE;
  std::cout << " * PERFORMING CC ENERGY CHECK " << std::endl;
  std::cout << "CC_TEST_REF=" << CC_TEST_REF+ref <<std::endl;

  resFile.readData("/CC/CORRELATION_ENERGY", &testE);
  refFile.readData("/CC/CORRELATION_ENERGY", &refE);

  EXPECT_NEAR( testE.real(), refE.real(), etol ) << "CC CORRELATION ENERGY TEST FAILED";
  EXPECT_NEAR( testE.imag(), refE.imag(), etol ) << "CC CORRELATION ENERGY TEST FAILED";

  if (checkReferenceEnergy){
    double testE, refE;
    // Check reference energy values
    std::cout << " * PERFORMING REFERENCE ENERGY TEST\n";

    resFile.readData("/CC/REFERENCE_ENERGY",&testE);
    refFile.readData("/CC/REFERENCE_ENERGY",&refE);

    EXPECT_NEAR(testE, refE, etol) << "CC REFERENCE ENERGY TEST FAILED";
  }
  std::vector<dcomplex> xDummy, yDummy;
  if (checkExcitedEnergy){
    // Check eigenvalues
    std::cout << "PERFORMING EXCITATION ENERGY TEST\n";

    auto evDim     = resFile.getDims("/CC/EXCITATION_ENERGIES");
    auto evDim_ref = refFile.getDims("/CC/EXCITATION_ENERGIES");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(evDim[0]); yDummy.resize(evDim_ref[0]);


    resFile.readData("/CC/EXCITATION_ENERGIES",&xDummy[0]);
    refFile.readData("/CC/EXCITATION_ENERGIES",&yDummy[0]);

    for(auto i = 0; i < evDim[0]; i++) {
      EXPECT_NEAR(xDummy[i].real(), yDummy[i].real(), etol) << "EXCITATION ENERGY TEST FAILED IN STATE = " << i;
      EXPECT_NEAR(xDummy[i].imag(), yDummy[i].imag(), etol) << "EXCITATION ENERGY TEST FAILED IN STATE = " << i;
    }
  }

  if (checkOsc){
    // Check Osc Strength
    std::cout << "PERFORMING OSC STRENGTH TEST\n";

    auto oscDim     = resFile.getDims("/CC/OSCILLATOR_STRENGTHS");
    auto oscDim_ref = refFile.getDims("/CC/OSCILLATOR_STRENGTHS");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(oscDim[0]); yDummy.resize(oscDim_ref[0]);

    resFile.readData("/CC/OSCILLATOR_STRENGTHS",&xDummy[0]);
    refFile.readData("/CC/OSCILLATOR_STRENGTHS",&yDummy[0]);

    for(auto i = 0; i < oscDim[0]; i++) {
      EXPECT_NEAR(xDummy[i].real(), yDummy[i].real(), osctol) << "OSC STRENGTH TEST FAILED IN STATE = " << i;
      EXPECT_NEAR(xDummy[i].imag(), yDummy[i].imag(), osctol) << "OSC STRENGTH TEST FAILED IN STATE = " << i;
    }
  }

  if (checkTriplesCorrection){
    // Check triples correction
    std::cout << "PERFORMING TRIPLES CORRECTION ENERGY TEST\n";

    dcomplex testE, refE;
    std::cout << " * PERFORMING CCSD(T) ENERGY CHECK " << std::endl;

    resFile.readData("/CC/CCSD(T)_CORRECTION", &testE);
    refFile.readData("/CC/CCSD(T)_CORRECTION", &refE);

    EXPECT_NEAR( testE.real(), refE.real(), etol ) << "CC CORRELATION ENERGY TEST FAILED";
    EXPECT_NEAR( testE.imag(), refE.imag(), etol ) << "CC CORRELATION ENERGY TEST FAILED";

    std::cout << " * PERFORMING CR-CC(2,3) (A-D variants) ENERGY CHECK " << std::endl;
    auto evDim     = resFile.getDims("/CC/CR-CC(2,3)_CORRECTION");
    auto evDim_ref = refFile.getDims("/CC/CR-CC(2,3)_CORRECTION");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(evDim[0]); yDummy.resize(evDim_ref[0]);


    resFile.readData("/CC/CR-CC(2,3)_CORRECTION",&xDummy[0]);
    refFile.readData("/CC/CR-CC(2,3)_CORRECTION",&yDummy[0]);

    std::vector<std::string> abcd = {"A", "B", "C", "D"};
    for(auto i = 0; i < evDim[0]; i++) {
      double cretol;
      // CR-C and -D are *NOT* orbital invariant (up to ~50-60 microhartree, maybe more)
      // Meanwhile, X2C orbitals may not be the same with different no. of threads
      i > 1 ? cretol = 6e-5 : cretol = etol;
      EXPECT_NEAR(xDummy[i].real(), yDummy[i].real(), cretol) << "CR-CC(2,3) CORRECTION TEST FAILED IN VARIANT = " << abcd[i];
      EXPECT_NEAR(xDummy[i].imag(), yDummy[i].imag(), cretol) << "CR-CC(2,3) CORRECTION TEST FAILED IN VARIANT = " << abcd[i];
    }
  }


#endif


};
