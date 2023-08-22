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



// B3LYP
TEST( KS_KEYWORD, KEYWORD_B3LYP ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_B3LYP", "water_sto-3g_B3LYP.bin.ref" );

}

// BHANDH
TEST( KS_KEYWORD, KEYWORD_BHANDH ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_BHANDH", "water_sto-3g_BHANDH.bin.ref" );

}

// BHANDHLYP
TEST( KS_KEYWORD, KEYWORD_BHANDHLYP ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_BHANDHLYP", "water_sto-3g_BHANDHLYP.bin.ref" );

}

// BLYP
TEST( KS_KEYWORD, KEYWORD_BLYP ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_BLYP", "water_sto-3g_BLYP.bin.ref" );

}

// LSDA
TEST( KS_KEYWORD, KEYWORD_LSDA ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_LSDA", "water_sto-3g_LSDA.bin.ref" );

}

// PBE0
TEST( KS_KEYWORD, KEYWORD_PBE0 ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_PBE0", "water_sto-3g_PBE0.bin.ref" );

}

// PBEXPBEC
TEST( KS_KEYWORD, KEYWORD_PBEXPBEC ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_PBEXPBEC", "water_sto-3g_PBEXPBEC.bin.ref" );

}

// SLATER
TEST( KS_KEYWORD, KEYWORD_SLATER ) {

  CQSCFTEST( "scf/serial/rks/water_sto-3g_SLATER", "water_sto-3g_SLATER.bin.ref" );

}


// B3LYP
TEST( KS_FUNC, KS_CART_B3LYP ) {

  CQSCFTEST( "scf/serial/rks/water_cc-pVTZ_cart_B3LYP", "water_cc-pVTZ_cart_B3LYP.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );

}



