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

#include "dynamics.hpp"

// rk4 propagator
// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, ch2o_hf_hf_rk4_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/ch2o_hf_hf_rk4_tvpb",
    "ch2o_hf_hf_rk4_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics   
TEST( BORT_DYNAMICS, hcn_hf_hf_rk4_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/hcn_hf_hf_rk4_tvpb",
    "hcn_hf_hf_rk4_tvpb.bin.ref" );

}


// mmut propagator
// ch2o b3lyp/epc17 neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, ch2o_b3lyp_epc17_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/ch2o_b3lyp_epc17_mmut_tvpb",
    "ch2o_b3lyp_epc17_mmut_tvpb.bin.ref" );

}

// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, ch2o_hf_hf_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/ch2o_hf_hf_mmut_tvpb",
    "ch2o_hf_hf_mmut_tvpb.bin.ref" );

}

// hcn b3lyp/epc17 neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, hcn_b3lyp_epc17_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/hcn_b3lyp_epc17_mmut_tvpb",
    "hcn_b3lyp_epc17_mmut_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics   
TEST( BORT_DYNAMICS, hcn_hf_hf_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/hcn_hf_hf_mmut_tvpb",
    "hcn_hf_hf_mmut_tvpb.bin.ref" );

}


// magnus2 propagator
// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, ch2o_hf_hf_magnus2_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/ch2o_hf_hf_magnus2_tvpb",
    "ch2o_hf_hf_magnus2_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics   
TEST( BORT_DYNAMICS, hcn_hf_hf_magnus2_tvpb ) {

  CQDYNAMICSTEST( "dynamics/serial/bort-ehrenfest/hcn_hf_hf_magnus2_tvpb",
    "hcn_hf_hf_magnus2_tvpb.bin.ref" );

}



#ifdef _CQ_DO_PARTESTS
// rk4 propagator
// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_ch2o_hf_hf_rk4_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/ch2o_hf_hf_rk4_tvpb",
    "ch2o_hf_hf_rk4_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics   
TEST( BORT_DYNAMICS, PAR_hcn_hf_hf_rk4_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/hcn_hf_hf_rk4_tvpb",
    "hcn_hf_hf_rk4_tvpb.bin.ref" );

}

// mmut propagator
// ch2o b3lyp/epc17 neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_ch2o_b3lyp_epc17_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/ch2o_b3lyp_epc17_mmut_tvpb",
    "ch2o_b3lyp_epc17_mmut_tvpb.bin.ref" );

}

// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_ch2o_hf_hf_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/ch2o_hf_hf_mmut_tvpb",
    "ch2o_hf_hf_mmut_tvpb.bin.ref" );

}

// hcn b3lyp/epc17 neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_hcn_b3lyp_epc17_mmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/hcn_b3lyp_epc17_mmut_tvpb",
    "hcn_b3lyp_epc17_mmut_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_hcn_hf_hf_mmmut_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/hcn_hf_hf_mmut_tvpb",
    "hcn_hf_hf_mmut_tvpb.bin.ref" );

}

// magnus2 propagator
// ch2o hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_ch2o_hf_hf_magnus2_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/ch2o_hf_hf_magnus2_tvpb",
    "ch2o_hf_hf_magnus2_tvpb.bin.ref" );

}

// hcn hf neo bort-ehrenfest rt dynamics
TEST( BORT_DYNAMICS, PAR_hcn_hf_hf_magnus2_tvpb ) {

  CQDYNAMICSTEST( "dynamics/parallel/bort-ehrenfest/hcn_hf_hf_magnus2_tvpb",
    "hcn_hf_hf_magnus2_tvpb.bin.ref" );

}

#endif
