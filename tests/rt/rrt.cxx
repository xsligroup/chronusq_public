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



// Water RHF
TEST( RHF_RT, water_631Gd_rhf_mmut ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rhf_mmut",
    "water_6-31Gd_rhf_mmut.bin.ref",1e-6 );

}

TEST( RHF_RT, water_631Gd_rhf_magnus2 ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rhf_magnus2",
   "water_6-31Gd_rhf_magnus2.bin.ref",1e-6 );

}

TEST( RHF_RT, water_631Gd_rhf_forwardeuler ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rhf_forwardeuler",
    "water_6-31Gd_rhf_forwardeuler.bin.ref",1e-6 );

}

TEST( RHF_RT, water_631Gd_rhf_mmut_delta ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rhf_mmut_delta",
    "water_6-31Gd_rhf_mmut_delta.bin.ref",1e-6 );

}

#ifdef _CQ_DO_PARTESTS
// Parallel Water RHF
TEST( RHF_RT, PAR_water_631Gd_rhf_mmut ) {

  CQRTTEST( "rt/parallel/rrt/water_6-31Gd_rhf_mmut",
    "water_6-31Gd_rhf_mmut.bin.ref",1e-6 );

}

TEST( RHF_RT, PAR_water_631Gd_rhf_magnus2 ) {

  CQRTTEST( "rt/parallel/rrt/water_6-31Gd_rhf_magnus2",
    "water_6-31Gd_rhf_magnus2.bin.ref",1e-6 );

}

//TEST( RHF_RT, PAR_water_631Gd_rhf_forwardeuler,1e-6 ) {
//
//  CQRTTEST( rt/parallel/rrt/water_6-31Gd_rhf_forwardeuler,
//    water_6-31Gd_rhf_forwardeuler.bin );
//
//}

// // Parallel Water RB3LYP
TEST( RKS_RT, PAR_water_631Gd_rb3lyp_mmut ) {

  CQRTTEST( "rt/parallel/rrt/water_6-31Gd_rb3lyp_mmut",
    "water_6-31Gd_rb3lyp_mmut.bin.ref",1e-6 );

}

TEST( RKS_RT, PAR_water_631Gd_rb3lyp_magnus2 ) {

  CQRTTEST( "rt/parallel/rrt/water_6-31Gd_rb3lyp_magnus2",
    "water_6-31Gd_rb3lyp_magnus2.bin.ref",1e-6 );

}

#endif
//RESTART
TEST( RESTART_RT, Restart_water_631Gd_rhf_mmut ) {
  CQRTRESTARTTEST("water_6-31Gd_rhf_mmut_restart_mid.bin.ref",
    "rt/serial/rrt/water_6-31Gd_rhf_mmut_restart",
    "water_6-31Gd_rhf_mmut_restart.bin.ref" );
}














// // Water RB3LYP
TEST( RKS_RT, water_631Gd_rb3lyp_mmut ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rb3lyp_mmut",
    "water_6-31Gd_rb3lyp_mmut.bin.ref",1e-6 );

}

TEST( RKS_RT, water_631Gd_rb3lyp_magnus2 ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rb3lyp_magnus2",
    "water_6-31Gd_rb3lyp_magnus2.bin.ref",1e-6 );

}

//TEST( RKS_RT, water_631Gd_rb3lyp_forwardeuler ) {
//
//  CQRTTEST( rt/serial/rrt/water_6-31Gd_rb3lyp_forwardeuler,
//    water_6-31Gd_rb3lyp_forwardeuler.bin );
//
//}

TEST( RKS_RT, water_631Gd_rb3lyp_mmut_delta ) {

  CQRTTEST( "rt/serial/rrt/water_6-31Gd_rb3lyp_mmut_delta",
    "water_6-31Gd_rb3lyp_mmut_delta.bin.ref",1e-6 );

}






//TEST( RKS_RT, PAR_water_631Gd_rb3lyp_forwardeuler ) {
//
//  CQRTTEST( rt/parallel/rrt/water_6-31Gd_rb3lyp_forwardeuler,
//    water_6-31Gd_rb3lyp_forwardeuler.bin );
//
//}




