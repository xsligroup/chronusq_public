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

#include "scf.hpp"



// Water RKS/cc-pvdz with RKS/6-31Gd guess (READDEN)
TEST( PROJ_KS, Water_ccpVDZ_RKS_6_31Gd_RKSGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_cc-pvdz_RKS_6-31Gd_RKSGuess",
    "water_cc-pvdz_RKS_6-31Gd_RKSGuess_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_6-31Gd_RKS.scr.bin" );

};

// Water RKS/cc-pvdz with UKS/6-31Gd guess (READDEN)
TEST( PROJ_KS, Water_ccpVDZ_RKS_6_31Gd_UKSGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_cc-pvdz_RKS_6-31Gd_UKSGuess",
    "water_cc-pvdz_RKS_6-31Gd_UKSGuess_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_6-31Gd_UKS.scr.bin" );

};

// Water RKS/cc-pvdz with GKS/6-31Gd guess (READDEN)
TEST( PROJ_KS, Water_ccpVDZ_RKS_6_31Gd_GKSGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_cc-pvdz_RKS_6-31Gd_GKSGuess",
    "water_cc-pvdz_RKS_6-31Gd_GKSGuess_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_6-31Gd_GKS.scr.bin" );

};

// Water GKS/cc-pvdz with RKS/6-31Gd guess (READDEN)
TEST( PROJ_KS, Water_ccpVDZ_GKS_6_31Gd_RKSGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_cc-pvdz_GKS_6-31Gd_RKSGuess",
    "water_cc-pvdz_GKS_6-31Gd_RKSGuess_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_6-31Gd_RKS.scr.bin" );

};

// Water UKS/cc-pvdz with RKS/6-31Gd guess (READDEN)
TEST( PROJ_KS, Water_ccpVDZ_UKS_6_31Gd_RKSGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_cc-pvdz_UKS_6-31Gd_RKSGuess",
    "water_cc-pvdz_UKS_6-31Gd_RKSGuess_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_6-31Gd_RKS.scr.bin" );

};

#ifdef _CQ_DO_PARTESTS


#endif




