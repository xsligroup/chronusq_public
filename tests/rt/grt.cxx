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



// Oxygen 6-31G(d) MMUT, Magnus2, ForwardEuler
//TEST( X2CHF_RT, oxygen_x2chf_mmut ) {
//
//  CQRTTEST( rt/serial/grt/oxygen_6-31Gd_x2chf_mmut,
//    oxygen_6-31Gd_x2chf_mmut.bin );
//}
//
//TEST( X2CHF_RT, oxygen_x2chf_magnus2 ) {
//
//  CQRTTEST( rt/serial/grt/oxygen_6-31Gd_x2chf_magnus2,
//    oxygen_6-31Gd_x2chf_magnus2.bin );
//}
//
//TEST( X2CHF_RT, oxygen_x2chf_forwardeuler ) {
//
//  CQRTTEST( rt/serial/grt/oxygen_6-31Gd_x2chf_forwardeuler,
//    oxygen_6-31Gd_x2chf_forwardeuler.bin );
//}



//TEST( X2CHF_RT, oxygen_x2chf_magnus2_readmo ) {
//
//  CQRTRESTARTTEST("oxygen_6-31Gd_x2chf.bin",
//    "rt/serial/grt/oxygen_6-31Gd_x2chf_magnus2_readmo",
//    "oxygen_6-31Gd_x2chf_magnus2_readmo.bin" );
//}
//
//TEST( X2CHF_RT, oxygen_x2chf_forwardeuler_readmo ) {
//
//  CQRTRESTARTTEST("oxygen_6-31Gd_x2chf.bin",
//    "rt/serial/grt/oxygen_6-31Gd_x2chf_forwardeuler_readmo",
//    "oxygen_6-31Gd_x2chf_forwardeuler_readmo.bin" );
//}

//// Water 6-31G(d) MMUT, Magnus2, ForwardEuler
//TEST( X2CHF_RT, water_x2chf_mmut ) {
//
//  CQRTTEST( rt/serial/grt/water_6-31Gd_x2chf_mmut,
//    water_6-31Gd_x2chf_mmut.bin );
//}
//
//TEST( X2CHF_RT, water_x2chf_magnus2 ) {
//
//  CQRTTEST( rt/serial/grt/water_6-31Gd_x2chf_magnus2,
//    water_6-31Gd_x2chf_magnus2.bin );
//}
//
//TEST( X2CHF_RT, water_x2chf_forwardeuler ) {
//
//  CQRTTEST( rt/serial/grt/water_6-31Gd_x2chf_forwardeuler,
//    water_6-31Gd_x2chf_forwardeuler.bin );
//}

// Oxygen 6-31G(d) MMUT, readmo SCF
TEST( X2CHF_RT, oxygen_x2chf_mmut_readmo ) {

  CQRTRESTARTTEST("oxygen_6-31Gd_x2chf.bin.ref",
    "rt/serial/grt/oxygen_6-31Gd_x2chf_mmut_readmo",
    "oxygen_6-31Gd_x2chf_mmut_readmo.bin.ref" );
}

TEST( X2CHF_RT, ag2_sto3g_x2chf_mmut ) {

  CQRTTEST( "rt/serial/grt/ag2_sto3g_x2chf_mmut",
    "ag2_sto3g_x2chf_mmut.bin.ref" );
}

TEST( X2CHF_RT, ag2_sto3g_x2chf_mmut_delta ) {

  CQRTTEST( "rt/serial/grt/ag2_sto3g_x2chf_mmut_delta",
    "ag2_sto3g_x2chf_mmut_delta.bin.ref" );
}

TEST( X2CHF_RT, ag2_sto3g_x2chf_magnus2 ) {

  CQRTTEST( "rt/serial/grt/ag2_sto3g_x2chf_magnus2",
    "ag2_sto3g_x2chf_magnus2.bin.ref" );
}

//TEST( X2CHF_RT, ag2_sto3g_x2chf_forwardeuler ) {
//
//  CQRTTEST( rt/serial/grt/ag2_sto3g_x2chf_forwardeuler,
//    ag2_sto3g_x2chf_forwardeuler.bin.ref );
//}

// Oxygen 6-31G(d) MMUT, Magnus2, ForwardEuler
TEST( X2CKS_RT, oxygen_x2cb3lyp_mmut ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gd_x2cb3lyp_mmut",
    "oxygen_6-31Gd_x2cb3lyp_mmut.bin.ref" );
}

TEST( X2CKS_RT, oxygen_x2cb3lyp_magnus2 ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gd_x2cb3lyp_magnus2",
    "oxygen_6-31Gd_x2cb3lyp_magnus2.bin.ref" );
}

TEST( X2CKS_RT, oxygen_x2cb3lyp_mmut_delta ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gd_x2cb3lyp_mmut_delta",
    "oxygen_6-31Gd_x2cb3lyp_mmut_delta.bin.ref" );
}

//TEST( X2CKS_RT, oxygen_x2cb3lyp_forwardeuler ) {
//
//  CQRTTEST( rt/serial/grt/oxygen_6-31Gd_x2cb3lyp_forwardeuler,
//    oxygen_6-31Gd_x2cb3lyp_forwardeuler.bin );
//}

// 4CHF
TEST( FOURCHF_RT, oxygen_4chf_mmut ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gdunc_4chf_mmut",
    "oxygen_6-31Gdunc_4chf_mmut.bin.ref" );
}

TEST( FOURCHF_RT, oxygen_4chf_mmut_gaunt ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gdunc_4chf_mmut_gaunt",
    "oxygen_6-31Gdunc_4chf_mmut_gaunt.bin.ref" );
}

TEST( FOURCHF_RT, oxygen_4chf_mmut_breit ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gdunc_4chf_mmut_breit",
    "oxygen_6-31Gdunc_4chf_mmut_breit.bin.ref" );
}

TEST( FOURCHF_RT, oxygen_4chf_mmut_gaunt_rtgaunt ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gdunc_4chf_mmut_gaunt_rtgaunt",
    "oxygen_6-31Gdunc_4chf_mmut_gaunt_rtgaunt.bin.ref" );
}

TEST( FOURCHF_RT, oxygen_4chf_mmut_breit_rtbreit ) {

  CQRTTEST( "rt/serial/grt/oxygen_6-31Gdunc_4chf_mmut_breit_rtbreit",
    "oxygen_6-31Gdunc_4chf_mmut_breit_rtbreit.bin.ref" );
}

TEST( FOURCHF_RT, ag2_sto2gunc_4chf_mmut ) {

  CQRTTEST( "rt/serial/grt/ag2_sto2gunc_4chf_mmut",
    "ag2_sto2gunc_4chf_mmut.bin.ref",1e-6 );
}

TEST( FOURCHF_RT, ag2_sto2gunc_4chf_mmut_delta ) {

  CQRTTEST( "rt/serial/grt/ag2_sto2gunc_4chf_mmut_delta",
    "ag2_sto2gunc_4chf_mmut_delta.bin.ref",1e-6 );
}


#ifdef _CQ_DO_PARTESTS

//TEST( X2CHF_RT, PAR_oxygen_x2chf_mmut ) {
//
//  CQRTTEST( rt/parallel/grt/oxygen_6-31Gd_x2chf_mmut,
//    oxygen_6-31Gd_x2chf_mmut.bin );
//}

//TEST( X2CHF_RT, PAR_water_x2chf_mmut ) {
//
//  CQRTTEST( rt/parallel/grt/water_6-31Gd_x2chf_mmut,
//    water_6-31Gd_x2chf_mmut.bin );
//}

TEST( X2CHF_RT, PAR_oxygen_x2chf_mmut_readmo ) {

  CQRTRESTARTTEST("oxygen_6-31Gd_x2chf.bin.ref",
    "rt/parallel/grt/oxygen_6-31Gd_x2chf_mmut_readmo",
    "oxygen_6-31Gd_x2chf_mmut_readmo.bin.ref" );
}

TEST( X2CKS_RT, PAR_oxygen_x2cb3lyp_mmut ) {

  CQRTTEST( "rt/parallel/grt/oxygen_6-31Gd_x2cb3lyp_mmut",
    "oxygen_6-31Gd_x2cb3lyp_mmut.bin.ref" );
}

TEST( FOURCHF_RT, PAR_oxygen_4chf_mmut_breit_rtbreit ) {

  CQRTTEST( "rt/parallel/grt/oxygen_6-31Gdunc_4chf_mmut_breit_rtbreit",
    "oxygen_6-31Gdunc_4chf_mmut_breit_rtbreit.bin.ref" );
}
#endif

//RESTART
TEST( RESTART_RT, Restart_oxygen_631Gd_x2chf_mmut ) {

 CQRTRESTARTTEST("oxygen_6-31Gd_x2chf_mmut_restart_mid.bin.ref",
    "rt/serial/grt/oxygen_6-31Gd_x2chf_mmut_restart",
    "oxygen_6-31Gd_x2chf_mmut_restart.bin.ref" );
}
