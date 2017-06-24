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
TEST( UHF_RT, oxygen_uhf_mmut ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_uhf_mmut",
    "oxygen_6-31Gd_uhf_mmut.bin.ref" );

}

TEST( UHF_RT, oxygen_uhf_magnus2 ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_uhf_magnus2",
    "oxygen_6-31Gd_uhf_magnus2.bin.ref" );

}

TEST( UHF_RT, oxygen_uhf_forwardeuler ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_uhf_forwardeuler",
    "oxygen_6-31Gd_uhf_forwardeuler.bin.ref" );

}

TEST( UHF_RT, oxygen_uhf_mmut_delta ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_uhf_mmut_delta",
    "oxygen_6-31Gd_uhf_mmut_delta.bin.ref" );

}



#ifdef _CQ_DO_PARTESTS
// Oxygen 6-31G(d) MMUT, Magnus2, ForwardEuler
TEST( UHF_RT, PAR_oxygen_uhf_mmut ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_uhf_mmut",
    "oxygen_6-31Gd_uhf_mmut.bin.ref" );

}

TEST( UHF_RT, PAR_oxygen_magnus2 ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_uhf_magnus2",
    "oxygen_6-31Gd_uhf_magnus2.bin.ref" );

}

TEST( UHF_RT, PAR_oxygen_uhf_forwardeuler ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_uhf_forwardeuler",
    "oxygen_6-31Gd_uhf_forwardeuler.bin.ref" );

}
#endif





// Oxygen 6-31G(d) MMUT, Magnus2, ForwardEuler
TEST( UKS_RT, oxygen_ub3lyp_mmut ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_ub3lyp_mmut",
    "oxygen_6-31Gd_ub3lyp_mmut.bin.ref" );

}

TEST( UKS_RT, oxygen_ub3lyp_magnus2 ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_ub3lyp_magnus2",
    "oxygen_6-31Gd_ub3lyp_magnus2.bin.ref" );

}

TEST( UKS_RT, oxygen_ub3lyp_forwardeuler ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_ub3lyp_forwardeuler",
    "oxygen_6-31Gd_ub3lyp_forwardeuler.bin.ref" );

}

TEST( UKS_RT, oxygen_ub3lyp_mmut_delta ) {

  CQRTTEST( "rt/serial/urt/oxygen_6-31Gd_ub3lyp_mmut_delta",
    "oxygen_6-31Gd_ub3lyp_mmut_delta.bin.ref" );

}

#ifdef _CQ_DO_PARTESTS
// Oxygen 6-31G(d) MMUT, Magnus2, ForwardEuler
TEST( UKS_RT, PAR_oxygen_ub3lyp_mmut ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_ub3lyp_mmut",
    "oxygen_6-31Gd_ub3lyp_mmut.bin.ref" );

}

TEST( UKS_RT, PAR_oxygen_ub3lyp_magnus2 ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_ub3lyp_magnus2",
    "oxygen_6-31Gd_ub3lyp_magnus2.bin.ref" );

}

TEST( UKS_RT, PAR_oxygen_ub3lyp_forwardeuler ) {

  CQRTTEST( "rt/parallel/urt/oxygen_6-31Gd_ub3lyp_forwardeuler" ,
    "oxygen_6-31Gd_ub3lyp_forwardeuler.bin.ref" );
}

#endif

//RESTART
TEST( RESTART_RT, Restart_oxygen_631Gd_ub3lyp_mmut ) {

  CQRTRESTARTTEST("oxygen_6-31Gd_ub3lyp_mmut_restart_mid.bin.ref",
    "rt/serial/urt/oxygen_6-31Gd_ub3lyp_mmut_restart",
    "oxygen_6-31Gd_ub3lyp_mmut_restart.bin.ref" );

}

TEST( RESTART_RT, Restart_oxygen_631Gd_ub3lyp_magnus2 ) {

  CQRTRESTARTTEST("oxygen_6-31Gd_ub3lyp_magnus2_restart_mid.bin.ref",
    "rt/serial/urt/oxygen_6-31Gd_ub3lyp_magnus2_restart",
    "oxygen_6-31Gd_ub3lyp_magnus2_restart.bin.ref" );

}

