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

#include "scf.hpp"

#include <cstdio>

// Water 6-31G(d) {0., 0.01, 0.} Electric Field test
TEST( MISC_SCF, Water_631Gd_ed_0_0pt01_0 ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0_0.01_0", 
    "water_6-31Gd_ed_0_0.01_0.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true  );
 
};

// Water 6-31G(d) {0.01, 0., 0.} Electric Field test
TEST( MISC_SCF, Water_631Gd_ed_0pt01_0_0 ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0.01_0_0", 
    "water_6-31Gd_ed_0.01_0_0.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// Water 6-31G(d) {0., 0., 0.01} Electric Field test
TEST( MISC_SCF, Water_631Gd_ed_0_0_0pt01 ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_ed_0_0_0.01", 
    "water_6-31Gd_ed_0_0_0.01.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// O2 Minimal basis
TEST( MISC_SCF, O2_STO3G ) {

  CQSCFTEST( "scf/serial/uhf/oxygen_sto-3g", "oxygen_sto-3g.bin.ref" );

};

#ifndef _CQ_GENERATE_TESTS

// Water 6-31G(d) read basis 
TEST( RHF, Water_631Gd_FindBasis ) {

  std::string newFileName = TEST_OUT "scf/serial/rhf/basisFile.gbs";

  remove( newFileName.c_str() );

  std::ifstream oldFile( TEST_ROOT "scf/serial/rhf/basisFile.gbs" );
  std::ofstream newFile( newFileName ); 

  newFile << oldFile.rdbuf();
  newFile.flush();

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_findBasis", "water_6-31Gd.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};


// Water 6-31G(d) input basis 
TEST( RHF, Water_631Gd_InputBasis ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd_inBasis", "water_6-31Gd.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

#endif
