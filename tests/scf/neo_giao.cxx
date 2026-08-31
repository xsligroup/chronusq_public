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


// NEO-GIAO 1AU (Tanner test)
TEST( NEO_GIAO, hcn_ccpvdz_prot4d) {

  CQNEOSCFTEST( "scf/serial/neo_giao/hcn_ccpvdz_prot4d", "hcn_ccpvdz_prot4d.bin.ref" );
 
}

// NEO-GIAO Rotation Test, only test energy
TEST( NEO_GIAO, heh2_ghf_xz) {

  CQNEOSCFTEST( "scf/serial/neo_giao/heh2_ghf_xz", "heh2_ghf_z.bin.ref" , 1e-6, 
    false, false, false, false, false, "no", true, false );
 
}

#ifdef _CQ_DO_PARTESTS

TEST( NEO_GIAO, par_hcn_ccpvdz_prot4d) {

  CQNEOSCFTEST( "scf/parallel/neo_giao/hcn_ccpvdz_prot4d", "hcn_ccpvdz_prot4d.bin.ref" );
 
}

TEST( NEO_GIAO, par_heh2_ghf_xz) {

  CQNEOSCFTEST( "scf/parallel/neo_giao/heh2_ghf_xz", "heh2_ghf_z.bin.ref" , 1e-6, 
    false, false, false, false, false, "no", true, false );
 
}

#endif




