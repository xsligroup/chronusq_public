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


// B3LYP / cc-pVTZ
TEST( RKS, Water_ccpVTZ_B3LYP ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref" );

}

// BLYP / cc-pVTZ
TEST( RKS, Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref" );

}

// LSDA / cc-pVTZ
TEST( RKS, Water_ccpVTZ_LSDA ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS

// SMP B3LYP / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_B3LYP ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref" );

}

// SMP BLYP / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref" );

}

// SMP LSDA / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_LSDA ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref", 2e-8 );

}

#endif



