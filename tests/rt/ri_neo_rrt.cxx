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




// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with EAUX algorithm
TEST( RI_NEO_RHF_RT, hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux ) {

  CQRTTEST( "rt/serial/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux",
    "hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with PAUX algorithm
TEST( RI_NEO_RHF_RT, hcn_sto3g_protsp_rhf_magnus2_delta_y_paux ) {

  CQRTTEST( "rt/serial/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_paux",
    "hcn_sto3g_protsp_rhf_magnus2_delta_y_paux.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with CONNECTOR algorithm
TEST( RI_NEO_RHF_RT, hcn_sto3g_protsp_rhf_magnus2_delta_y_connector ) {

  CQRTTEST( "rt/serial/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_connector",
    "hcn_sto3g_protsp_rhf_magnus2_delta_y_connector.bin.ref" );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with COMBINEAUXBASIS(two-component RI) algorithm
TEST( RI_NEO_RHF_RT, hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis ) {

  CQRTTEST( "rt/serial/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis",
    "hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis.bin.ref" );

}




/***
#ifdef _CQ_DO_PARTESTS

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with EAUX algorithm
TEST( RI_NEO_RHF_RT, par_hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux ) {

  CQRTTEST( rt/parallel/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux,
    hcn_sto3g_protsp_rhf_magnus2_delta_y_eaux.bin.ref );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with PAUX algorithm
TEST( RI_NEO_RHF_RT, par_hcn_sto3g_protsp_rhf_magnus2_delta_y_paux ) {

  CQRTTEST( rt/parallel/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_paux,
    hcn_sto3g_protsp_rhf_magnus2_delta_y_paux.bin.ref );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with CONNECTOR algorithm
TEST( RI_NEO_RHF_RT, par_hcn_sto3g_protsp_rhf_magnus2_delta_y_connector ) {

  CQRTTEST( rt/parallel/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_connector,
    hcn_sto3g_protsp_rhf_magnus2_delta_y_connector.bin.ref );

}

// Water sto3g/prot-sp Delta Spike (along Y), restarting with magnus2 algotithm, NEO-CD with COMBINEAUXBASIS(two-component RI) algorithm
TEST( RI_NEO_RHF_RT, par_hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis ) {

  CQRTTEST( rt/parallel/ri_neo_rrt/hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis,
    hcn_sto3g_protsp_rhf_magnus2_delta_y_combineauxbasis.bin.ref );

}
#endif
***/


