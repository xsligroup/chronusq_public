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
#include <cerr.hpp>

// Directory containing reference files
#define RESP_TEST_REF TEST_ROOT "/resp/reference/"

using namespace ChronusQ;



static void CQRESPTEST( std::string in, std::string ref ) {

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    RESP_TEST_REF + ref, "", false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin","", false);

#endif


};


static void CQRESPREFTEST( std::string in, std::string ref ) {

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _CQ_GENERATE_TESTS

  // Assumes user added bin to RESP_TEST_REF
  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    RESP_TEST_REF + ref, "", true);

#else

  std::ifstream  src(RESP_TEST_REF + ref, std::ios::binary);
  std::ofstream  dst(TEST_OUT + in + ".bin", std::ios::binary);
  dst << src.rdbuf();
  dst.flush();

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin","", true);

#endif


};

static void CQRESPSCRTEST( std::string in, std::string ref, std::string scr ) {

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef _CQ_GENERATE_TESTS

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    RESP_TEST_REF + ref, RESP_TEST_REF + scr, false);

#else

  RunChronusQ(TEST_ROOT + in + ".inp","STDOUT", 
    TEST_OUT + in + ".bin", RESP_TEST_REF + scr, false);

#endif


};


static void CQRESTEST( bool checkProp, std::string in, std::string ref, 
    bool runCQ = true, double tol = 1e-06 ) {

  if( runCQ ) CQRESPTEST(in,ref);
  if( MPIRank() != 0 ) return;

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(RESP_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
   
  std::vector<double> xDummy, yDummy; 

  // Get DEMIN from RESFILE
  double DEMIN = 0;
  resFile.readData("/RESP/RESIDUE/DEMIN", &DEMIN);

  // See if RESFILE is full problem
  int isFull;
  resFile.readData("/RESP/DOFULL", &isFull);

  // Check eigenvalues
  std::cerr << "PERFORMING EIGENVALUE TEST\n";

  auto evDim     = resFile.getDims("/RESP/RESIDUE/EIGENVALUES");
  auto evDim_ref = refFile.getDims("/RESP/RESIDUE/EIGENVALUES");

  xDummy.clear(); yDummy.clear();
  xDummy.resize(evDim[0]); yDummy.resize(evDim_ref[0]);
  

  resFile.readData("/RESP/RESIDUE/EIGENVALUES",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/EIGENVALUES",&yDummy[0]);


  // Find where DEMIN cuts off
  size_t eigOff = isFull ? 0 : 
    std::distance( 
      yDummy.begin(),
      std::find_if(yDummy.begin(), yDummy.end(), 
        [&](double x){ return x > DEMIN; })
    );


  if( eigOff )
    std::cerr << "  * EIGENVALUE OFFSET = " << eigOff << std::endl; 
  


  for(auto i = 0; i < evDim[0]; i++)
    EXPECT_NEAR( xDummy[i], yDummy[i + eigOff], tol ) <<
      "EIGENVALUE TEST FAILED IO = " << i;

  if( not checkProp ) return;

  // Check Osc Strength
  std::cerr << "PERFORMING OSC STRENGTH TEST\n";

  auto oscDim     = resFile.getDims("/RESP/RESIDUE/OSC_STRENGTH");
  auto oscDim_ref = refFile.getDims("/RESP/RESIDUE/OSC_STRENGTH");

  xDummy.clear(); yDummy.clear();
  xDummy.resize(oscDim[0]); yDummy.resize(oscDim_ref[0]);
  
  resFile.readData("/RESP/RESIDUE/OSC_STRENGTH",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/OSC_STRENGTH",&yDummy[0]);

  for(auto i = 0; i < oscDim[0]; i++)
    EXPECT_NEAR( xDummy[i], yDummy[i + eigOff], tol ) <<
      "OSC STRENGTH TEST FAILED IO = " << i;


  // Check TDipole (Length)
  std::cerr << "PERFORMING TDIPOLE (LENGTH) TEST\n";

  auto tDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH");
  auto tDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH");
  ASSERT_EQ(tDipoleDim.size(),2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tDipoleDim[0]*tDipoleDim[1]); 
  yDummy.resize(tDipoleDim_ref[0]*tDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH",&yDummy[0]);

  for(auto i = 0; i < tDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[3*i]), std::abs(yDummy[3*(i + eigOff)]), tol ) <<
      "TDIPOLE (LENGTH,X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+1]), std::abs(yDummy[3*(i + eigOff)+1]), tol ) <<
      "TDIPOLE (LENGTH,Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+2]), std::abs(yDummy[3*(i + eigOff)+2]), tol ) <<
      "TDIPOLE (LENGTH,Z) TEST FAILED IO = " << i;

  }

  // Check TDipole (Velocity)
  std::cerr << "PERFORMING TDIPOLE (VELOCITY) TEST\n";

  tDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY");
  tDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY");
  ASSERT_EQ(tDipoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tDipoleDim[0]*tDipoleDim[1]); 
  yDummy.resize(tDipoleDim_ref[0]*tDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY",&yDummy[0]);

  for(auto i = 0; i < tDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[3*i]), std::abs(yDummy[3*(i + eigOff)]), tol ) <<
      "TDIPOLE (VELOCITY,X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+1]), std::abs(yDummy[3*(i + eigOff)+1]), tol ) <<
      "TDIPOLE (VELOCITY,Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+2]), std::abs(yDummy[3*(i + eigOff)+2]), tol ) <<
      "TDIPOLE (VELOCITY,Z) TEST FAILED IO = " << i;

  }

  // Check TQuadrupole (Length)
  std::cerr << "PERFORMING TQUADRUPOLE (LENGTH) TEST\n";

  auto tQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH");
  auto tQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH");
  ASSERT_EQ(tQuadrupoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tQuadrupoleDim[0]*tQuadrupoleDim[1]); 
  yDummy.resize(tQuadrupoleDim_ref[0]*tQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH",&yDummy[0]);

  for(auto i = 0; i < tQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[6*i]), std::abs(yDummy[6*(i + eigOff)]), tol ) <<
      "TQUADRUPOLE (LENGTH,XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+1]), std::abs(yDummy[6*(i + eigOff)+1]), tol ) <<
      "TQUADRUPOLE (LENGTH,XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+2]), std::abs(yDummy[6*(i + eigOff)+2]), tol ) <<
      "TQUADRUPOLE (LENGTH,XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+3]), std::abs(yDummy[6*(i + eigOff)+3]), tol ) <<
      "TQUADRUPOLE (LENGTH,YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+4]), std::abs(yDummy[6*(i + eigOff)+4]), tol ) <<
      "TQUADRUPOLE (LENGTH,YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+5]), std::abs(yDummy[6*(i + eigOff)+5]), tol ) <<
      "TQUADRUPOLE (LENGTH,ZZ) TEST FAILED IO = " << i;

  }


  // Check TQuadrupole (Velocity)
  std::cerr << "PERFORMING TQUADRUPOLE (VELOCITY) TEST\n";

  tQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY");
  tQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY");
  ASSERT_EQ(tQuadrupoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tQuadrupoleDim[0]*tQuadrupoleDim[1]); 
  yDummy.resize(tQuadrupoleDim_ref[0]*tQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY",&yDummy[0]);

  for(auto i = 0; i < tQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[6*i]), std::abs(yDummy[6*(i + eigOff)]), tol ) <<
      "TQUADRUPOLE (VELOCITY,XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+1]), std::abs(yDummy[6*(i + eigOff)+1]), tol ) <<
      "TQUADRUPOLE (VELOCITY,XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+2]), std::abs(yDummy[6*(i + eigOff)+2]), tol ) <<
      "TQUADRUPOLE (VELOCITY,XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+3]), std::abs(yDummy[6*(i + eigOff)+3]), tol ) <<
      "TQUADRUPOLE (VELOCITY,YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+4]), std::abs(yDummy[6*(i + eigOff)+4]), tol ) <<
      "TQUADRUPOLE (VELOCITY,YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+5]), std::abs(yDummy[6*(i + eigOff)+5]), tol ) <<
      "TQUADRUPOLE (VELOCITY,ZZ) TEST FAILED IO = " << i;

  }




  // Check TOctupole (Length)
  std::cerr << "PERFORMING TOCTUPOLE (LENGTH) TEST\n";

  auto tOctupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH");
  auto tOctupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH");
  ASSERT_EQ(tOctupoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tOctupoleDim[0]*tOctupoleDim[1]); 
  yDummy.resize(tOctupoleDim_ref[0]*tOctupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH",&yDummy[0]);

  for(auto i = 0; i < tOctupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[10*i]), std::abs(yDummy[10*(i + eigOff)]), tol ) <<
      "TOCTUPOLE (LENGTH,XXX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+1]), std::abs(yDummy[10*(i + eigOff)+1]), tol ) <<
      "TOCTUPOLE (LENGTH,XXY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+2]), std::abs(yDummy[10*(i + eigOff)+2]), tol ) <<
      "TOCTUPOLE (LENGTH,XXZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+3]), std::abs(yDummy[10*(i + eigOff)+3]), tol ) <<
      "TOCTUPOLE (LENGTH,XYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+4]), std::abs(yDummy[10*(i + eigOff)+4]), tol ) <<
      "TOCTUPOLE (LENGTH,XYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+5]), std::abs(yDummy[10*(i + eigOff)+5]), tol ) <<
      "TOCTUPOLE (LENGTH,XZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+6]), std::abs(yDummy[10*(i + eigOff)+6]), tol ) <<
      "TOCTUPOLE (LENGTH,YYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+7]), std::abs(yDummy[10*(i + eigOff)+7]), tol ) <<
      "TOCTUPOLE (LENGTH,YYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+8]), std::abs(yDummy[10*(i + eigOff)+8]), tol ) <<
      "TOCTUPOLE (LENGTH,YZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+9]), std::abs(yDummy[10*(i + eigOff)+9]), tol ) <<
      "TOCTUPOLE (LENGTH,ZZZ) TEST FAILED IO = " << i;

  }



  // Check TOctupole (Velocity)
  std::cerr << "PERFORMING TOCTUPOLE (VELOCITY) TEST\n";

  tOctupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY");
  tOctupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY");
  ASSERT_EQ(tOctupoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tOctupoleDim[0]*tOctupoleDim[1]); 
  yDummy.resize(tOctupoleDim_ref[0]*tOctupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY",&yDummy[0]);

  for(auto i = 0; i < tOctupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[10*i]), std::abs(yDummy[10*(i + eigOff)]), tol ) <<
      "TOCTUPOLE (VELOCITY,XXX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+1]), std::abs(yDummy[10*(i + eigOff)+1]), tol ) <<
      "TOCTUPOLE (VELOCITY,XXY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+2]), std::abs(yDummy[10*(i + eigOff)+2]), tol ) <<
      "TOCTUPOLE (VELOCITY,XXZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+3]), std::abs(yDummy[10*(i + eigOff)+3]), tol ) <<
      "TOCTUPOLE (VELOCITY,XYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+4]), std::abs(yDummy[10*(i + eigOff)+4]), tol ) <<
      "TOCTUPOLE (VELOCITY,XYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+5]), std::abs(yDummy[10*(i + eigOff)+5]), tol ) <<
      "TOCTUPOLE (VELOCITY,XZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+6]), std::abs(yDummy[10*(i + eigOff)+6]), tol ) <<
      "TOCTUPOLE (VELOCITY,YYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+7]), std::abs(yDummy[10*(i + eigOff)+7]), tol ) <<
      "TOCTUPOLE (VELOCITY,YYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+8]), std::abs(yDummy[10*(i + eigOff)+8]), tol ) <<
      "TOCTUPOLE (VELOCITY,YZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[10*i+9]), std::abs(yDummy[10*(i + eigOff)+9]), tol ) <<
      "TOCTUPOLE (VELOCITY,ZZZ) TEST FAILED IO = " << i;

  }






  // Check TMagDipole
  std::cerr << "PERFORMING T-MAGDIPOLE TEST\n";

  auto tMagDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE");
  auto tMagDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE");
  ASSERT_EQ(tMagDipoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tMagDipoleDim[0]*tMagDipoleDim[1]); 
  yDummy.resize(tMagDipoleDim_ref[0]*tMagDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE",&yDummy[0]);

  for(auto i = 0; i < tMagDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[3*i]), std::abs(yDummy[3*(i + eigOff)]), tol ) <<
      "TMAGDIPOLE (X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+1]), std::abs(yDummy[3*(i + eigOff)+1]), tol ) <<
      "TMAGDIPOLE (Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[3*i+2]), std::abs(yDummy[3*(i + eigOff)+2]), tol ) <<
      "TMAGDIPOLE (Z) TEST FAILED IO = " << i;

  }


  // Check TMagQuadrupole (Length)
  std::cerr << "PERFORMING T-MAGQUADRUPOLE TEST\n";

  auto tMagQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE");
  auto tMagQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE");
  ASSERT_EQ(tMagQuadrupoleDim.size(), 2);

  xDummy.clear(); yDummy.clear();
  xDummy.resize(tMagQuadrupoleDim[0]*tMagQuadrupoleDim[1]); 
  yDummy.resize(tMagQuadrupoleDim_ref[0]*tMagQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE",&yDummy[0]);

  for(auto i = 0; i < tMagQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummy[6*i]), std::abs(yDummy[6*(i + eigOff)]), tol ) <<
      "TMAGQUADRUPOLE (XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+1]), std::abs(yDummy[6*(i + eigOff)+1]), tol ) <<
      "TMAGQUADRUPOLE (XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+2]), std::abs(yDummy[6*(i + eigOff)+2]), tol ) <<
      "TMAGQUADRUPOLE (XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+3]), std::abs(yDummy[6*(i + eigOff)+3]), tol ) <<
      "TMAGQUADRUPOLE (YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+4]), std::abs(yDummy[6*(i + eigOff)+4]), tol ) <<
      "TMAGQUADRUPOLE (YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummy[6*i+5]), std::abs(yDummy[6*(i + eigOff)+5]), tol ) <<
      "TMAGQUADRUPOLE (ZZ) TEST FAILED IO = " << i;

  }
#endif

};



template <typename T>
static void CQFDRTEST( bool checkLHerOps, bool checkLAntiHerOps,
    std::string in, std::string ref, double tol = 1e-05) {

  CQRESPTEST(in,ref);
  if( MPIRank() != 0 ) return;

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(RESP_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
   
  std::vector<double> xDummy, yDummy; 
  std::vector<T> txDummy, tyDummy; 


  double damp;
  resFile.readData("/RESP/FDR/DAMP",&damp);

  // Check cases where LHS op is Hermetian
  if( checkLHerOps ) {

    // Check OPA cross setction (LENGTH)
    std::cerr << "PERFORMING OPA CROSS SECTION TEST\n";

    if( std::abs(damp) > tol) {  
      auto opaDim = resFile.getDims("/RESP/FDR/OPA_CROSS_SECTION_EDA");

      xDummy.clear(); yDummy.clear();
      xDummy.resize(opaDim[0]); yDummy.resize(opaDim[0]);
      
      resFile.readData("/RESP/FDR/OPA_CROSS_SECTION_EDA",&xDummy[0]);
      refFile.readData("/RESP/FDR/OPA_CROSS_SECTION_EDA",&yDummy[0]);

      for(auto i = 0; i < opaDim[0]; i++)
        EXPECT_NEAR(xDummy[i], yDummy[i],  tol) << 
          "OPA CROSS SECTION TEST FAILED IO = " << i; 
    }



    // Check ED_ED_POLAR (Length)
    std::cerr << "PERFORMING ED-ED POLARIZABILITY (LENGTH) TEST\n";

    auto ededPolarDim = resFile.getDims("/RESP/FDR/ED_ED_POLARIZABILITY_LENGTH");
    ASSERT_EQ(ededPolarDim.size(), 3);

    txDummy.clear(); tyDummy.clear();
    txDummy.resize(ededPolarDim[0]*ededPolarDim[1]*ededPolarDim[1]); 
    tyDummy.resize(ededPolarDim[0]*ededPolarDim[1]*ededPolarDim[1]); 
    
    resFile.readData("/RESP/FDR/ED_ED_POLARIZABILITY_LENGTH",&txDummy[0]);
    refFile.readData("/RESP/FDR/ED_ED_POLARIZABILITY_LENGTH",&tyDummy[0]);

    for(auto i = 0; i < ededPolarDim[0]; i++) {

      EXPECT_NEAR( std::abs(txDummy[9*i]), std::abs(tyDummy[9*i]), tol ) <<
        "ED ED POLAR (LENGTH,X;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+1]), std::abs(tyDummy[9*i+1]), tol ) <<
        "ED ED POLAR (LENGTH,Y;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+2]), std::abs(tyDummy[9*i+2]), tol ) <<
        "ED ED POLAR (LENGTH,Z;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+3]), std::abs(tyDummy[9*i+3]), tol ) <<
        "ED ED POLAR (LENGTH,X;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+4]), std::abs(tyDummy[9*i+4]), tol ) <<
        "ED ED POLAR (LENGTH,Y;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+5]), std::abs(tyDummy[9*i+5]), tol ) <<
        "ED ED POLAR (LENGTH,Z;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+6]), std::abs(tyDummy[9*i+6]), tol ) <<
        "ED ED POLAR (LENGTH,X;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+7]), std::abs(tyDummy[9*i+7]), tol ) <<
        "ED ED POLAR (LENGTH,Y;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+8]), std::abs(tyDummy[9*i+8]), tol ) <<
        "ED ED POLAR (LENGTH,Z;Z) TEST FAILED IO = " << i; 

    }



    // Check EQ_ED_POLAR (Length)
    std::cerr << "PERFORMING EQ-ED POLARIZABILITY (LENGTH) TEST\n";

    auto eqedPolarDim = resFile.getDims("/RESP/FDR/EQ_ED_POLARIZABILITY_LENGTH");
    ASSERT_EQ(eqedPolarDim.size(),3);

    txDummy.clear(); tyDummy.clear();
    txDummy.resize(eqedPolarDim[0]*eqedPolarDim[1]*eqedPolarDim[2]); 
    tyDummy.resize(eqedPolarDim[0]*eqedPolarDim[1]*eqedPolarDim[2]); 
    
    resFile.readData("/RESP/FDR/EQ_ED_POLARIZABILITY_LENGTH",&txDummy[0]);
    refFile.readData("/RESP/FDR/EQ_ED_POLARIZABILITY_LENGTH",&tyDummy[0]);

    for(auto i = 0; i < eqedPolarDim[0]; i++) {

      EXPECT_NEAR( std::abs(txDummy[18*i]), std::abs(tyDummy[18*i]), tol ) <<
        "EQ ED POLAR (LENGTH,XX;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+1]), std::abs(tyDummy[18*i+1]), tol ) <<
        "EQ ED POLAR (LENGTH,XY;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+2]), std::abs(tyDummy[18*i+2]), tol ) <<
        "EQ ED POLAR (LENGTH,XZ;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+3]), std::abs(tyDummy[18*i+3]), tol ) <<
        "EQ ED POLAR (LENGTH,YY;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+4]), std::abs(tyDummy[18*i+4]), tol ) <<
        "EQ ED POLAR (LENGTH,YZ;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+5]), std::abs(tyDummy[18*i+5]), tol ) <<
        "EQ ED POLAR (LENGTH,ZZ;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+6]), std::abs(tyDummy[18*i+6]), tol ) <<
        "EQ ED POLAR (LENGTH,XX;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+7]), std::abs(tyDummy[18*i+7]), tol ) <<
        "EQ ED POLAR (LENGTH,XY;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+8]), std::abs(tyDummy[18*i+8]), tol ) <<
        "EQ ED POLAR (LENGTH,XZ;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+9]), std::abs(tyDummy[18*i+9]), tol ) <<
        "EQ ED POLAR (LENGTH,YY;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+10]), std::abs(tyDummy[18*i+10]), tol ) <<
        "EQ ED POLAR (LENGTH,YZ;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+11]), std::abs(tyDummy[18*i+11]), tol ) <<
        "EQ ED POLAR (LENGTH,ZZ;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+12]), std::abs(tyDummy[18*i+12]), tol ) <<
        "EQ ED POLAR (LENGTH,XX;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+13]), std::abs(tyDummy[18*i+13]), tol ) <<
        "EQ ED POLAR (LENGTH,XY;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+14]), std::abs(tyDummy[18*i+14]), tol ) <<
        "EQ ED POLAR (LENGTH,XZ;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+15]), std::abs(tyDummy[18*i+15]), tol ) <<
        "EQ ED POLAR (LENGTH,YY;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+16]), std::abs(tyDummy[18*i+16]), tol ) <<
        "EQ ED POLAR (LENGTH,YZ;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[18*i+17]), std::abs(tyDummy[18*i+17]), tol ) <<
        "EQ ED POLAR (LENGTH,ZZ;Z) TEST FAILED IO = " << i; 

    }
  }



  // Check cases where LHS op is Anti-Hermetian
  if( checkLAntiHerOps ) {
    // Check MD_ED_POLAR (Length)
    std::cerr << "PERFORMING MD-ED POLARIZABILITY (LENGTH) TEST\n";

    auto mdedPolarDim = resFile.getDims("/RESP/FDR/MD_ED_POLARIZABILITY_LENGTH");
    ASSERT_EQ(mdedPolarDim.size(), 3);

    txDummy.clear(); tyDummy.clear();
    txDummy.resize(mdedPolarDim[0]*mdedPolarDim[1]*mdedPolarDim[1]); 
    tyDummy.resize(mdedPolarDim[0]*mdedPolarDim[1]*mdedPolarDim[1]); 
    
    resFile.readData("/RESP/FDR/MD_ED_POLARIZABILITY_LENGTH",&txDummy[0]);
    refFile.readData("/RESP/FDR/MD_ED_POLARIZABILITY_LENGTH",&tyDummy[0]);

    for(auto i = 0; i < mdedPolarDim[0]; i++) {

      EXPECT_NEAR( std::abs(txDummy[9*i]), std::abs(tyDummy[9*i]), tol ) <<
        "MD ED POLAR (LENGTH,X;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+1]), std::abs(tyDummy[9*i+1]), tol ) <<
        "MD ED POLAR (LENGTH,Y;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+2]), std::abs(tyDummy[9*i+2]), tol ) <<
        "MD ED POLAR (LENGTH,Z;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+3]), std::abs(tyDummy[9*i+3]), tol ) <<
        "MD ED POLAR (LENGTH,X;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+4]), std::abs(tyDummy[9*i+4]), tol ) <<
        "MD ED POLAR (LENGTH,Y;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+5]), std::abs(tyDummy[9*i+5]), tol ) <<
        "MD ED POLAR (LENGTH,Z;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+6]), std::abs(tyDummy[9*i+6]), tol ) <<
        "MD ED POLAR (LENGTH,X;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+7]), std::abs(tyDummy[9*i+7]), tol ) <<
        "MD ED POLAR (LENGTH,Y;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+8]), std::abs(tyDummy[9*i+8]), tol ) <<
        "MD ED POLAR (LENGTH,Z;Z) TEST FAILED IO = " << i; 

    }


    // Check MD_MD_POLAR
    std::cerr << "PERFORMING MD-MD POLARIZABILITY TEST\n";

    auto mdmdPolarDim = resFile.getDims("/RESP/FDR/MD_MD_POLARIZABILITY");
    ASSERT_EQ(mdmdPolarDim.size(), 3);

    txDummy.clear(); tyDummy.clear();
    txDummy.resize(mdmdPolarDim[0]*mdmdPolarDim[1]*mdmdPolarDim[1]); 
    tyDummy.resize(mdmdPolarDim[0]*mdmdPolarDim[1]*mdmdPolarDim[1]); 
    
    resFile.readData("/RESP/FDR/MD_MD_POLARIZABILITY",&txDummy[0]);
    refFile.readData("/RESP/FDR/MD_MD_POLARIZABILITY",&tyDummy[0]);

    for(auto i = 0; i < mdmdPolarDim[0]; i++) {

      EXPECT_NEAR( std::abs(txDummy[9*i]), std::abs(tyDummy[9*i]), tol ) <<
        "MD MD POLAR (LENGTH,X;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+1]), std::abs(tyDummy[9*i+1]), tol ) <<
        "MD MD POLAR (LENGTH,Y;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+2]), std::abs(tyDummy[9*i+2]), tol ) <<
        "MD MD POLAR (LENGTH,Z;X) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+3]), std::abs(tyDummy[9*i+3]), tol ) <<
        "MD MD POLAR (LENGTH,X;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+4]), std::abs(tyDummy[9*i+4]), tol ) <<
        "MD MD POLAR (LENGTH,Y;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+5]), std::abs(tyDummy[9*i+5]), tol ) <<
        "MD MD POLAR (LENGTH,Z;Y) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+6]), std::abs(tyDummy[9*i+6]), tol ) <<
        "MD MD POLAR (LENGTH,X;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+7]), std::abs(tyDummy[9*i+7]), tol ) <<
        "MD MD POLAR (LENGTH,Y;Z) TEST FAILED IO = " << i; 

      EXPECT_NEAR( std::abs(txDummy[9*i+8]), std::abs(tyDummy[9*i+8]), tol ) <<
        "MD MD POLAR (LENGTH,Z;Z) TEST FAILED IO = " << i; 

    }
  }
#endif

};


static void CQMORTEST( bool checkLHerOps, bool checkLAntiHerOps,
    std::string in, std::string ref, double tol = 1e-02) {

  CQFDRTEST<dcomplex>(checkLHerOps, checkLAntiHerOps, in, ref, tol);

}


static void CQCRESTEST( bool checkProp, std::string in, std::string ref, 
    bool runCQ = true, double tol = 1e-06, bool readBin = false, std::string scr="no" ) {

  if( runCQ ){ 
    if( readBin ) CQRESPREFTEST(in,ref);
    else if(scr != "no") CQRESPSCRTEST(in,ref,scr);
    else CQRESPTEST(in,ref);
  }
  if( MPIRank() != 0 ) return;

#ifndef _CQ_GENERATE_TESTS

  SafeFile refFile(RESP_TEST_REF + ref,true);
  SafeFile resFile(TEST_OUT + in + ".bin",true);
   
  std::vector<double> xDummy, yDummy; 

  // Get DEMIN from RESFILE
  double DEMIN = 0;
  resFile.readData("/RESP/RESIDUE/DEMIN", &DEMIN);

  // See if RESFILE is full problem
  int isFull;
  resFile.readData("/RESP/DOFULL", &isFull);

  // Check eigenvalues
  std::cerr << "PERFORMING EIGENVALUE TEST\n";

  auto evDim     = resFile.getDims("/RESP/RESIDUE/EIGENVALUES");
  auto evDim_ref = refFile.getDims("/RESP/RESIDUE/EIGENVALUES");

  xDummy.clear(); yDummy.clear();
  xDummy.resize(evDim[0]); yDummy.resize(evDim_ref[0]);
  

  resFile.readData("/RESP/RESIDUE/EIGENVALUES",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/EIGENVALUES",&yDummy[0]);


  // Find where DEMIN cuts off
  size_t eigOff = isFull ? 0 : 
    std::distance( 
      yDummy.begin(),
      std::find_if(yDummy.begin(), yDummy.end(), 
        [&](double x){ return x > DEMIN; })
    );


  if( eigOff )
    std::cerr << "  * EIGENVALUE OFFSET = " << eigOff << std::endl; 
  


  for(auto i = 0; i < evDim[0]; i++)
    EXPECT_NEAR( xDummy[i], yDummy[i + eigOff], tol ) <<
      "EIGENVALUE TEST FAILED IO = " << i;

  if( not checkProp ) return;

  // Check Osc Strength
  std::cerr << "PERFORMING OSC STRENGTH TEST\n";

  auto oscDim     = resFile.getDims("/RESP/RESIDUE/OSC_STRENGTH");
  auto oscDim_ref = refFile.getDims("/RESP/RESIDUE/OSC_STRENGTH");

  xDummy.clear(); yDummy.clear();
  xDummy.resize(oscDim[0]); yDummy.resize(oscDim_ref[0]);
  
  resFile.readData("/RESP/RESIDUE/OSC_STRENGTH",&xDummy[0]);
  refFile.readData("/RESP/RESIDUE/OSC_STRENGTH",&yDummy[0]);

  for(auto i = 0; i < oscDim[0]; i++)
    EXPECT_NEAR( xDummy[i], yDummy[i + eigOff], tol ) <<
      "OSC STRENGTH TEST FAILED IO = " << i;

  // Since properties are complex
  std::vector<dcomplex> xDummyC, yDummyC; 
  dcomplex tolC(tol,tol); 
   

  // Check TDipole (Length)
  std::cerr << "PERFORMING TDIPOLE (LENGTH) TEST\n";

  auto tDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH");
  auto tDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH");
  ASSERT_EQ(tDipoleDim.size(),2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tDipoleDim[0]*tDipoleDim[1]); 
  yDummyC.resize(tDipoleDim_ref[0]*tDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH",&yDummyC[0]);

  for(auto i = 0; i < tDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[3*i]), std::abs(yDummyC[3*(i + eigOff)]), std::abs(tolC) ) <<
      "TDIPOLE (LENGTH,X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+1]), std::abs(yDummyC[3*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TDIPOLE (LENGTH,Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+2]), std::abs(yDummyC[3*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TDIPOLE (LENGTH,Z) TEST FAILED IO = " << i;

  }

  // Check TDipole (Velocity)
  std::cerr << "PERFORMING TDIPOLE (VELOCITY) TEST\n";

  tDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY");
  tDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY");
  ASSERT_EQ(tDipoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tDipoleDim[0]*tDipoleDim[1]); 
  yDummyC.resize(tDipoleDim_ref[0]*tDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY",&yDummyC[0]);

  for(auto i = 0; i < tDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[3*i]), std::abs(yDummyC[3*(i + eigOff)]), std::abs(tolC) ) <<
      "TDIPOLE (VELOCITY,X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+1]), std::abs(yDummyC[3*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TDIPOLE (VELOCITY,Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+2]), std::abs(yDummyC[3*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TDIPOLE (VELOCITY,Z) TEST FAILED IO = " << i;

  }

  // Check TQuadrupole (Length)
  std::cerr << "PERFORMING TQUADRUPOLE (LENGTH) TEST\n";

  auto tQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH");
  auto tQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH");
  ASSERT_EQ(tQuadrupoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tQuadrupoleDim[0]*tQuadrupoleDim[1]); 
  yDummyC.resize(tQuadrupoleDim_ref[0]*tQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH",&yDummyC[0]);

  for(auto i = 0; i < tQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[6*i]), std::abs(yDummyC[6*(i + eigOff)]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+1]), std::abs(yDummyC[6*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+2]), std::abs(yDummyC[6*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+3]), std::abs(yDummyC[6*(i + eigOff)+3]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+4]), std::abs(yDummyC[6*(i + eigOff)+4]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+5]), std::abs(yDummyC[6*(i + eigOff)+5]), std::abs(tolC) ) <<
      "TQUADRUPOLE (LENGTH,ZZ) TEST FAILED IO = " << i;

  }


  // Check TQuadrupole (Velocity)
  std::cerr << "PERFORMING TQUADRUPOLE (VELOCITY) TEST\n";

  tQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY");
  tQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY");
  ASSERT_EQ(tQuadrupoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tQuadrupoleDim[0]*tQuadrupoleDim[1]); 
  yDummyC.resize(tQuadrupoleDim_ref[0]*tQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY",&yDummyC[0]);

  for(auto i = 0; i < tQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[6*i]), std::abs(yDummyC[6*(i + eigOff)]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+1]), std::abs(yDummyC[6*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+2]), std::abs(yDummyC[6*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+3]), std::abs(yDummyC[6*(i + eigOff)+3]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+4]), std::abs(yDummyC[6*(i + eigOff)+4]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+5]), std::abs(yDummyC[6*(i + eigOff)+5]), std::abs(tolC) ) <<
      "TQUADRUPOLE (VELOCITY,ZZ) TEST FAILED IO = " << i;

  }




  // Check TOctupole (Length)
  std::cerr << "PERFORMING TOCTUPOLE (LENGTH) TEST\n";

  auto tOctupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH");
  auto tOctupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH");
  ASSERT_EQ(tOctupoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tOctupoleDim[0]*tOctupoleDim[1]); 
  yDummyC.resize(tOctupoleDim_ref[0]*tOctupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH",&yDummyC[0]);

  for(auto i = 0; i < tOctupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[10*i]), std::abs(yDummyC[10*(i + eigOff)]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XXX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+1]), std::abs(yDummyC[10*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XXY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+2]), std::abs(yDummyC[10*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XXZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+3]), std::abs(yDummyC[10*(i + eigOff)+3]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+4]), std::abs(yDummyC[10*(i + eigOff)+4]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+5]), std::abs(yDummyC[10*(i + eigOff)+5]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,XZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+6]), std::abs(yDummyC[10*(i + eigOff)+6]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,YYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+7]), std::abs(yDummyC[10*(i + eigOff)+7]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,YYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+8]), std::abs(yDummyC[10*(i + eigOff)+8]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,YZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+9]), std::abs(yDummyC[10*(i + eigOff)+9]), std::abs(tolC) ) <<
      "TOCTUPOLE (LENGTH,ZZZ) TEST FAILED IO = " << i;

  }



  // Check TOctupole (Velocity)
  std::cerr << "PERFORMING TOCTUPOLE (VELOCITY) TEST\n";

  tOctupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY");
  tOctupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY");
  ASSERT_EQ(tOctupoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tOctupoleDim[0]*tOctupoleDim[1]); 
  yDummyC.resize(tOctupoleDim_ref[0]*tOctupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY",&yDummyC[0]);

  for(auto i = 0; i < tOctupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[10*i]), std::abs(yDummyC[10*(i + eigOff)]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XXX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+1]), std::abs(yDummyC[10*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XXY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+2]), std::abs(yDummyC[10*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XXZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+3]), std::abs(yDummyC[10*(i + eigOff)+3]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+4]), std::abs(yDummyC[10*(i + eigOff)+4]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+5]), std::abs(yDummyC[10*(i + eigOff)+5]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,XZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+6]), std::abs(yDummyC[10*(i + eigOff)+6]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,YYY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+7]), std::abs(yDummyC[10*(i + eigOff)+7]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,YYZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+8]), std::abs(yDummyC[10*(i + eigOff)+8]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,YZZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[10*i+9]), std::abs(yDummyC[10*(i + eigOff)+9]), std::abs(tolC) ) <<
      "TOCTUPOLE (VELOCITY,ZZZ) TEST FAILED IO = " << i;

  }






  // Check TMagDipole
  std::cerr << "PERFORMING T-MAGDIPOLE TEST\n";

  auto tMagDipoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE");
  auto tMagDipoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE");
  ASSERT_EQ(tMagDipoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tMagDipoleDim[0]*tMagDipoleDim[1]); 
  yDummyC.resize(tMagDipoleDim_ref[0]*tMagDipoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE",&yDummyC[0]);

  for(auto i = 0; i < tMagDipoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[3*i]), std::abs(yDummyC[3*(i + eigOff)]), std::abs(tolC) ) <<
      "TMAGDIPOLE (X) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+1]), std::abs(yDummyC[3*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TMAGDIPOLE (Y) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[3*i+2]), std::abs(yDummyC[3*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TMAGDIPOLE (Z) TEST FAILED IO = " << i;

  }


  // Check TMagQuadrupole (Length)
  std::cerr << "PERFORMING T-MAGQUADRUPOLE TEST\n";

  auto tMagQuadrupoleDim     = resFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE");
  auto tMagQuadrupoleDim_ref = refFile.getDims("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE");
  ASSERT_EQ(tMagQuadrupoleDim.size(), 2);

  xDummyC.clear(); yDummyC.clear();
  xDummyC.resize(tMagQuadrupoleDim[0]*tMagQuadrupoleDim[1]); 
  yDummyC.resize(tMagQuadrupoleDim_ref[0]*tMagQuadrupoleDim_ref[1]); 
  
  resFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE",&xDummyC[0]);
  refFile.readData("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE",&yDummyC[0]);

  for(auto i = 0; i < tMagQuadrupoleDim[0]; i++) {

    EXPECT_NEAR( std::abs(xDummyC[6*i]), std::abs(yDummyC[6*(i + eigOff)]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (XX) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+1]), std::abs(yDummyC[6*(i + eigOff)+1]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (XY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+2]), std::abs(yDummyC[6*(i + eigOff)+2]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (XZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+3]), std::abs(yDummyC[6*(i + eigOff)+3]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (YY) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+4]), std::abs(yDummyC[6*(i + eigOff)+4]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (YZ) TEST FAILED IO = " << i;

    EXPECT_NEAR( std::abs(xDummyC[6*i+5]), std::abs(yDummyC[6*(i + eigOff)+5]), std::abs(tolC) ) <<
      "TMAGQUADRUPOLE (ZZ) TEST FAILED IO = " << i;

  }
#endif

};


