/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <string>

#include <cxxapi/procedural.hpp>
#include <util/files.hpp>
#include <util/mpi.hpp>

#define CC_TEST_REF TEST_ROOT "coupledcluster/reference/"

using namespace ChronusQ;



static void CQCCTEST( std::string in, std::string ref, std::string restart_bin = "",
                      bool checkReferenceEnergy = false,
                      bool checkExcitedEnergy = false,
                      bool checkOsc = false,
                      bool checkDipoles = false,
                      bool checkTriplesCorrection = false,
                      bool ifRHF = false,
                      double etol = 1e-7,
                      double osctol = 1e-5
                    ) {

  MPI_Barrier(MPI_COMM_WORLD);

  bool rstExist = false;
  if ( ! restart_bin.empty() ) {
    std::ifstream  src(CC_TEST_REF + restart_bin, std::ios::binary);
    std::ofstream  dst(CQTestOut(in,".bin"), std::ios::binary);
    dst << src.rdbuf();
    dst.flush();
    rstExist = true;
  }

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    CC_TEST_REF + ref, CQTestOut(in,".scr"), rstExist);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    CQTestOut(in,".bin"),CQTestOut(in,".scr"), rstExist);

  if(MPIRank(MPI_COMM_WORLD) != 0) return;

          
  SafeFile refFile(CC_TEST_REF + ref,true,true);
  SafeFile resFile(CQTestOut(in,".bin"),true);


  // CC ENERGY CHECK
  dcomplex testE, refE;
  std::cout << " * PERFORMING CC ENERGY CHECK " << std::endl;
  std::cout << "CC_TEST_REF=" << CC_TEST_REF+ref <<std::endl;

  // handle double vs dcomplex values
  double testEReal = 0.0, refEReal = 0.0;
  if (ifRHF) {
    resFile.readData("/CC/CORRELATION_ENERGY", &testEReal);
    refFile.readData("/CC/CORRELATION_ENERGY", &refEReal);
  } else {
    resFile.readData("/CC/CORRELATION_ENERGY", &testE);
    testEReal = testE.real();
    refFile.readData("/CC/CORRELATION_ENERGY", &refE);
    refEReal = refE.real();
  }

  EXPECT_NEAR( testEReal, refEReal, etol ) << "CC CORRELATION ENERGY TEST FAILED";
  if (not ifRHF) {
    EXPECT_NEAR( testE.imag(), refE.imag(), etol ) << "CC CORRELATION ENERGY TEST FAILED";
  }


  // REFERENCE ENERGY CHECK
  if (checkReferenceEnergy){
    // Check reference energy values
    std::cout << " * PERFORMING REFERENCE ENERGY TEST\n";

    resFile.readData("/CC/REFERENCE_ENERGY",&testEReal);
    refFile.readData("/CC/REFERENCE_ENERGY",&refEReal);

    EXPECT_NEAR(testEReal, refEReal, etol) << "CC REFERENCE ENERGY TEST FAILED";
  }

  // EXCITATION ENERGY CHECK
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

  // OSCILLATOR STRENGTH CHECK
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

  // GS, ES-GS/GS-ES, AND ES-ES DIPOLES CHECK
  if (checkDipoles) {
    // Check GS dipole
    std::cout << "PERFORMING GS DIPOLE TEST\n";

    auto dipDim     = resFile.getDims("CC/GROUND_STATE_DIPOLE");
    auto dipDim_ref = refFile.getDims("CC/GROUND_STATE_DIPOLE");

    xDummy.clear(); yDummy.clear();
    xDummy.resize(dipDim[0]); yDummy.resize(dipDim_ref[0]);

    std::vector<double> xDummyReal, yDummyReal;
    xDummyReal.clear(); yDummyReal.clear();
    xDummyReal.resize(dipDim[0]); yDummyReal.resize(dipDim_ref[0]);

    if (ifRHF) {
      resFile.readData("CC/GROUND_STATE_DIPOLE",&xDummyReal[0]);
      refFile.readData("CC/GROUND_STATE_DIPOLE",&yDummyReal[0]);
    } else {
      resFile.readData("CC/GROUND_STATE_DIPOLE",&xDummy[0]);
      refFile.readData("CC/GROUND_STATE_DIPOLE",&yDummy[0]);
      for (auto i = 0; i < dipDim_ref[0]; i++) {
        xDummyReal[i] = xDummy[i].real();
        yDummyReal[i] = yDummy[i].real();
      }
    }

    for(auto i = 0; i < dipDim_ref[0]; i++) {
      EXPECT_NEAR(xDummyReal[i], yDummyReal[i], osctol) << "GS DIPOLE TEST FAILED IN COMPONENT = " << static_cast<char>('X' + i);
      EXPECT_NEAR(xDummy[i].imag(), yDummy[i].imag(), osctol) << "GS DIPOLE TEST FAILED IN COMPONENT = " << static_cast<char>('X' + i);
    }

    // Check GS-ES dipole
    std::cout << "PERFORMING GS-ES DIPOLE TEST\n";

    dipDim.clear(); dipDim_ref.clear();
    dipDim     = resFile.getDims("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE");
    dipDim_ref = refFile.getDims("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE");
    ASSERT_EQ(dipDim.size(), 2);
    auto totDim = dipDim[0] * dipDim[1];
    ASSERT_EQ(totDim, dipDim_ref[0] * dipDim_ref[1]);

    xDummy.clear(); yDummy.clear();
    xDummy.resize(totDim); yDummy.resize(totDim);

    xDummyReal.clear(); yDummyReal.clear();
    xDummyReal.resize(totDim); yDummyReal.resize(totDim);

    if (ifRHF) {
      resFile.readData("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE",&xDummyReal[0]);
      refFile.readData("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE",&yDummyReal[0]);
    } else {
      resFile.readData("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE",&xDummy[0]);
      refFile.readData("CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE",&yDummy[0]);
      for (auto i = 0; i < totDim; i++) {
        xDummyReal[i] = xDummy[i].real();
        yDummyReal[i] = yDummy[i].real();
      }
    }

    for(auto i = 0; i < dipDim_ref[0]; i++) { //cartesian
      for (auto j = 0; j < dipDim_ref[1]; j++) { //state
        auto idx = i * dipDim_ref[1] + j;
        EXPECT_NEAR(abs(xDummyReal[idx]), abs(yDummyReal[idx]), osctol)
        << "GS-ES DIPOLE TEST FAILED IN STATE = " << j << " AND COMPONENT = " << static_cast<char>('X'+i);
        EXPECT_NEAR(abs(xDummy[idx].imag()), abs(yDummy[idx].imag()), osctol)
        << "GS-ES DIPOLE TEST FAILED IN STATE = " << j << " AND COMPONENT = " << static_cast<char>('X'+i);
      }
    }

    // Check ES-GS dipole
    std::cout << "PERFORMING ES-GS DIPOLE TEST\n";

    dipDim.clear(); dipDim_ref.clear();
    dipDim     = resFile.getDims("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE");
    dipDim_ref = refFile.getDims("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE");
    ASSERT_EQ(dipDim.size(), 2);
    totDim = dipDim[0] * dipDim[1];
    ASSERT_EQ(totDim, dipDim_ref[0] * dipDim_ref[1]);

    xDummy.clear(); yDummy.clear();
    xDummy.resize(totDim); yDummy.resize(totDim);

    xDummyReal.clear(); yDummyReal.clear();
    xDummyReal.resize(totDim); yDummyReal.resize(totDim);

    if (ifRHF) {
      resFile.readData("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE",&xDummyReal[0]);
      refFile.readData("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE",&yDummyReal[0]);
    } else {
      resFile.readData("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE",&xDummy[0]);
      refFile.readData("CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE",&yDummy[0]);
      for (auto i = 0; i < totDim; i++) {
        xDummyReal[i] = xDummy[i].real();
        yDummyReal[i] = yDummy[i].real();
      }
    }

    for(auto i = 0; i < dipDim_ref[0]; i++) { //cartesian
      for (auto j = 0; j < dipDim_ref[1]; j++) { //state
        auto idx = i * dipDim_ref[1] + j;
        EXPECT_NEAR(abs(xDummyReal[idx]), abs(yDummyReal[idx]), osctol)
        << "ES-GS DIPOLE TEST FAILED IN STATE = " << j << " AND COMPONENT = " << static_cast<char>('X'+i);
        EXPECT_NEAR(abs(xDummy[idx].imag()), abs(yDummy[idx].imag()), osctol)
        << "ES-GS DIPOLE TEST FAILED IN STATE = " << j << " AND COMPONENT = " << static_cast<char>('X'+i);
      }
    }

    // Check ES-ES dipole
    std::cout << "PERFORMING ES-ES DIPOLE TEST\n";

    dipDim.clear(); dipDim_ref.clear();
    dipDim     = resFile.getDims("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE");
    dipDim_ref = refFile.getDims("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE");
    ASSERT_EQ(dipDim.size(), 3);
    totDim = dipDim[0] * dipDim[1] * dipDim[2];
    ASSERT_EQ(totDim, dipDim_ref[0] * dipDim_ref[1] * dipDim_ref[2]);

    xDummy.clear(); yDummy.clear();
    xDummy.resize(totDim); yDummy.resize(totDim);

    xDummyReal.clear(); yDummyReal.clear();
    xDummyReal.resize(totDim); yDummyReal.resize(totDim);

    if (ifRHF) {
      resFile.readData("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE",&xDummyReal[0]);
      refFile.readData("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE",&yDummyReal[0]);
    } else {
      resFile.readData("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE",&xDummy[0]);
      refFile.readData("CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE",&yDummy[0]);
      for (auto i = 0; i < totDim; i++) {
        xDummyReal[i] = xDummy[i].real();
        yDummyReal[i] = yDummy[i].real();
      }
    }

    for(auto i = 0; i < dipDim_ref[0]; i++) { //cartesian
      for (auto j = 0; j < dipDim_ref[1]; j++) { //bra state
        for (auto k = 0; k < dipDim_ref[2]; k++) { //ket state
          auto idx = (i * dipDim_ref[1] + j) * dipDim_ref[2] + k;
          EXPECT_NEAR(abs(xDummyReal[idx]), abs(yDummyReal[idx]), osctol)
          << "ES-ES DIPOLE TEST FAILED IN STATE PAIR = " << j << "," << k << " AND COMPONENT = " << static_cast<char>('X'+i);
          EXPECT_NEAR(abs(xDummy[idx].imag()), abs(yDummy[idx].imag()), osctol)
          << "ES-ES DIPOLE TEST FAILED IN STATE PAIR = " << j << "," << k << " AND COMPONENT = " << static_cast<char>('X'+i);
        }
      }
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
