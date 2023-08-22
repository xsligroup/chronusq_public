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


// B3LYP / 6-311pG**
TEST( UKS, Oxygen_6311pGss_B3LYP ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// BLYP / 6-311pG**
TEST( UKS, Oxygen_6311pGss_BLYP ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// LSDA / 6-311pG**
TEST( UKS, Oxygen_6311pGss_LSDA ) {

  CQSCFTEST( "scf/serial/uks/oxygen_6-311pG**_LSDA", "oxygen_6-311pG**_LSDA.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS

// SMP B3LYP / 6-311pG**
TEST( UKS, PAR_Oxygen_6311pGss_B3LYP ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP BLYP / 6-311pG**
TEST( UKS, PAR_Oxygen_6311pGss_BLYP ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref" );

}

// SMP LSDA / 6-311pG**
TEST( UKS, PAR_Oxygen_6311pGss_LSDA ) {

  CQSCFTEST( "scf/parallel/uks/oxygen_6-311pG**_LSDA", "oxygen_6-311pG**_LSDA.bin.ref" );

}

#endif


