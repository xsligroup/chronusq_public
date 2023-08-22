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

#include "resp.hpp"

#define CQRESTEST_IMPL(TNAME, IN, REF) \
TEST( RKS_RESIDUE, TNAME ) { CQRESTEST( true, IN, REF, true, 1e-5 ); }


// FULL DIMENSIONAL TESTS

// Water 6-31G(d) TDBLYP (RESIDUE)
CQRESTEST_IMPL( Water_631Gd_BLYP_RESIDUE,  
    "resp/serial/rresp_ks/water_6-31Gd_rblyp_residue",
    "water_6-31Gd_rblyp_residue.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) TDBLYP (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( Water_631Gd_BLYP_RESIDUE_GPLHR_DIRECT,  
      "resp/serial/rresp_ks/water_6-31Gd_rblyp_residue_gplhr_direct",
      "water_6-31Gd_rblyp_residue.bin.ref" )

#endif




#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) TDBLYP (RESIDUE)
  CQRESTEST_IMPL( PAR_Water_631Gd_BLYP_RESIDUE,  
      "resp/parallel/rresp_ks/water_6-31Gd_rblyp_residue",
      "water_6-31Gd_rblyp_residue.bin.ref") 

  // SMP Water 6-31G(d) TDBLYP (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( PAR_Water_631Gd_BLYP_RESIDUE_GPLHR_DIRECT,  
      "resp/parallel/rresp_ks/water_6-31Gd_rblyp_residue_gplhr_direct",
      "water_6-31Gd_rblyp_residue.bin.ref" )

#endif


#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
    
  // Water 6-31G(d) TDBLYP (RESIDUE) A+B / A-B
  CQRESTEST_IMPL( Water_631Gd_BLYP_RESIDUE_APB_AMB,  
      "resp/serial/rresp_ks/water_6-31Gd_rblyp_residue_apb_amb",
      "water_6-31Gd_rblyp_residue.bin.ref") 
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) TDBLYP (RESIDUE) A+B / A-B
    CQRESTEST_IMPL( PAR_Water_631Gd_BLYP_RESIDUE_APB_AMB,  
        "resp/parallel/rresp_ks/water_6-31Gd_rblyp_residue_apb_amb",
        "water_6-31Gd_rblyp_residue.bin.ref") 
  
  #endif


  // REDUCED TESTS
    
  // Water 6-31G(d) TDBLYP (RESIDUE) REDUCED
  CQRESTEST_IMPL( Water_631Gd_BLYP_RESIDUE_RED,  
      "resp/serial/rresp_ks/water_6-31Gd_rblyp_residue_reduced",
      "water_6-31Gd_rblyp_residue.bin.ref") 
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) TDBLYP (RESIDUE) REDUCED
    CQRESTEST_IMPL( PAR_Water_631Gd_BLYP_RESIDUE_RED,  
        "resp/parallel/rresp_ks/water_6-31Gd_rblyp_residue_reduced",
        "water_6-31Gd_rblyp_residue.bin.ref") 
  
  #endif


#endif
