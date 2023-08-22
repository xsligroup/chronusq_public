/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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


// HF NEO-Test with minimal basis set
TEST( NEO_RHF, water_sto3g_protsp ) {

  CQNEOSCFTEST( "scf/serial/neo_rhf/water_sto-3g_prot-sp_rhf", "water_sto-3g_prot-sp_rhf.bin.ref" );
 
};

TEST( NEO_RHF, coh2_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/serial/neo_rhf/coh2_ccpvdz_pb4d", "coh2_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

TEST( NEO_RHF, hcn_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/serial/neo_rhf/hcn_ccpvdz_pb4d", "hcn_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

TEST( NEO_RHF, h2o_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/serial/neo_rhf/h2o_ccpvdz_pb4d", "h2o_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

#ifdef _CQ_DO_PARTESTS

// HF NEO-Test with minimal basis set, parallel job
TEST( NEO_RHF, par_water_sto3g_protsp ) {

  CQNEOSCFTEST( "scf/parallel/neo_rhf/water_sto-3g_prot-sp_rhf", "water_sto-3g_prot-sp_rhf.bin.ref" );
 
};

TEST( NEO_RHF, par_coh2_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/parallel/neo_rhf/coh2_ccpvdz_pb4d", "coh2_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

TEST( NEO_RHF, par_hcn_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/parallel/neo_rhf/hcn_ccpvdz_pb4d", "hcn_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

TEST( NEO_RHF, par_h2o_ccpvdz_pb4d ) {

  CQNEOSCFTEST( "scf/parallel/neo_rhf/h2o_ccpvdz_pb4d", "h2o_ccpvdz_pb4d.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false );
 
};

#endif




