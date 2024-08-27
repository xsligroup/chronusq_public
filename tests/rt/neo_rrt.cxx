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

#include "rt.hpp"



// Water sto3g/prot-sp Delta Spike (along Y), restarting with mmut algotithm
TEST( NEO_RHF_RT, water_sto3g_protsp_rhf_mmut_delta_y ) {

  CQRTTEST( "rt/serial/neo_rrt/water_sto3g_protsp_rhf_mmut_delta_y",
    "water_sto3g_protsp_rhf_mmut_delta_y.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with forwardEuler algotithm
TEST( NEO_RHF_RT, water_sto3g_protsp_rhf_forwardeuler_delta_y ) {

  CQRTTEST( "rt/serial/neo_rrt/water_sto3g_protsp_rhf_forwardeuler_delta_y",
    "water_sto3g_protsp_rhf_forwardeuler_delta_y.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm
TEST( NEO_RHF_RT, water_sto3g_protsp_rhf_magnus2_delta_y ) {

  CQRTTEST( "rt/serial/neo_rrt/water_sto3g_protsp_rhf_magnus2_delta_y",
    "water_sto3g_protsp_rhf_magnus2_delta_y.bin.ref" );

}

// Water sto3g/prot-sp non-delta electric field (along Y), restarting with magnus2 algotithm   
TEST( NEO_RHF_RT, water_sto3g_protsp_rhf_magnus2_nondelta_y ) {

  CQRTTEST( "rt/serial/neo_rrt/water_sto3g_protsp_rhf_magnus2_nondelta_y",
    "water_sto3g_protsp_rhf_magnus2_nondelta_y.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS
// Water sto3g/prot-sp Delta Spike (along Y), restarting with mmut algotithm
TEST( NEO_RHF_RT, PAR_water_sto3g_protsp_rhf_mmut_delta_y ) {

  CQRTTEST( "rt/parallel/neo_rrt/water_sto3g_protsp_rhf_mmut_delta_y",
    "water_sto3g_protsp_rhf_mmut_delta_y.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with forwardEuler algotithm
TEST( NEO_RHF_RT, PAR_water_sto3g_protsp_rhf_forwardeuler_delta_y ) {

  CQRTTEST( "rt/parallel/neo_rrt/water_sto3g_protsp_rhf_forwardeuler_delta_y",
    "water_sto3g_protsp_rhf_forwardeuler_delta_y.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm
TEST( NEO_RHF_RT, PAR_water_sto3g_protsp_rhf_magnus2_delta_y ) {

  CQRTTEST( "rt/parallel/neo_rrt/water_sto3g_protsp_rhf_magnus2_delta_y",
    "water_sto3g_protsp_rhf_magnus2_delta_y.bin.ref" );

}

// Water sto3g/prot-sp non-delta electric field (along Y), restarting with magnus2 algotithm   
TEST( NEO_RHF_RT, PAR_water_sto3g_protsp_rhf_magnus2_nondelta_y ) {

  CQRTTEST( "rt/parallel/neo_rrt/water_sto3g_protsp_rhf_magnus2_nondelta_y",
    "water_sto3g_protsp_rhf_magnus2_nondelta_y.bin.ref" );

}

#endif
