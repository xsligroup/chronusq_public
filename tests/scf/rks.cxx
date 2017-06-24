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

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );

}

// BLYP / cc-pVTZ
TEST( RKS, Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// LSDA / cc-pVTZ
TEST( RKS, Water_ccpVTZ_LSDA ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// PBEXPBEC / cc-pVTZ
TEST( RKS, Water_ccpVTZ_PBEXPBEC ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_PBEXPBEC", "water_cc-pVTZ_PBEXPBEC.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// PBE0 / cc-pVTZ
TEST( RKS, Water_ccpVTZ_PBE0 ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_PBE0", "water_cc-pVTZ_PBE0.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

#ifdef _CQ_DO_PARTESTS

// SMP B3LYP / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_B3LYP ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP BLYP / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP LSDA / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_LSDA ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_LSDA", "water_cc-pVTZ_LSDA.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );

}

// SMP PBEXPBEC / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_PBEXPBEC ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_PBEXPBEC", "water_cc-pVTZ_PBEXPBEC.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP PBE0 / cc-pVTZ
TEST( RKS, PAR_Water_ccpVTZ_PBE0 ) {

  CQSCFTEST( "scf/parallel/rks/water_cc-pVTZ_PBE0", "water_cc-pVTZ_PBE0.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

#endif



