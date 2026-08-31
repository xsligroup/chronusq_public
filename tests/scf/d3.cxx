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


// PBE0-D3 / def2-svp
TEST( D3, Water_def2svp_PBE0_D3 ) {

  CQSCFTEST( "scf/serial/d3/water_def2svp_PBE0-D3", "water_def2svp_PBE0-D3.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// B3LYP-D3Zero / def2-svp
TEST( D3, Water_def2svp_B3LYP_D3Zero ) {

  CQSCFTEST( "scf/serial/d3/water_def2svp_B3LYP-D3Zero", "water_def2svp_B3LYP-D3Zero.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// X2CB3LYP-D3 / 6-31G
TEST( D3, O2_631G_X2CB3LYP_D3 ) {

  CQSCFTEST( "scf/serial/d3/o2_631G_X2CB3LYP-D3", "o2_631G_X2CB3LYP-D3.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

// NEO: PBE0-D3, EPC17/ ccpvdz, prot-pb4-d
TEST( D3, Water_ccpvdz_pb4d_PBE0_D3_epc17 ) {

  CQNEOSCFTEST( "scf/serial/d3/water_ccpvdz_pb4d_PBE0-D3_epc17", "water_ccpvdz_pb4d_PBE0-D3_epc17.bin.ref", 1e-6,
    true, true, true, true, false, "no", true, false );

}

#ifdef _CQ_DO_PARTESTS

// PBE0-D3 / def2-svp
TEST( D3, PAR_Water_def2svp_PBE0_D3 ) {

  CQSCFTEST( "scf/parallel/d3/water_def2svp_PBE0-D3", "water_def2svp_PBE0-D3.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// B3LYP-D3Zero / def2-svp
TEST( D3, PAR_Water_def2svp_B3LYP_D3Zero ) {

  CQSCFTEST( "scf/parallel/d3/water_def2svp_B3LYP-D3Zero", "water_def2svp_B3LYP-D3Zero.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// NEO: PBE0-D3, EPC17/ ccpvdz, prot-pb4-d
TEST( D3, PAR_Water_ccpvdz_pb4d_PBE0_D3_epc17 ) {

  CQNEOSCFTEST( "scf/parallel/d3/water_ccpvdz_pb4d_PBE0-D3_epc17", "water_ccpvdz_pb4d_PBE0-D3_epc17.bin.ref", 1e-6,
    true, true, true, true, false, "no", true, false );

}

#endif


