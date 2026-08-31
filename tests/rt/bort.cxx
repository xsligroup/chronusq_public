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

#include "rt.hpp"

// For one proton cases, the new code ignores J/K build which causes roundoff error to accumulate
// For these cases we use larger threshold to compare against legacy references.
static constexpr double oneProtonTol = 1e-6;

//rk4 propagator
// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, ch2o_hf_hf_rk4 ) {

  CQRTTEST( "rt/serial/bort/ch2o_hf_hf_rk4",
    "ch2o_hf_hf_rk4.bin.ref" );

}

// hcn hf neo bort rt dynamics   
TEST( BORT_NEO_RT, hcn_hf_hf_rk4 ) {

  CQRTTEST( "rt/serial/bort/hcn_hf_hf_rk4",
    "hcn_hf_hf_rk4.bin.ref", oneProtonTol );

}



//magnus2 propagator
// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, ch2o_hf_hf_magnus2 ) {

  CQRTTEST( "rt/serial/bort/ch2o_hf_hf_magnus2",
    "ch2o_hf_hf_magnus2.bin.ref" );

}

// hcn hf neo bort rt dynamics   
TEST( BORT_NEO_RT, hcn_hf_hf_magnus2 ) {

  CQRTTEST( "rt/serial/bort/hcn_hf_hf_magnus2",
    "hcn_hf_hf_magnus2.bin.ref", oneProtonTol );

}



//mmut propagator
// ch2o b3lyp/epc17 neo bort rt dynamics
TEST( BORT_NEO_RT, ch2o_b3lyp_epc17_mmut ) {

  CQRTTEST( "rt/serial/bort/ch2o_b3lyp_epc17_mmut",
    "ch2o_b3lyp_epc17_mmut.bin.ref" );

}

// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, ch2o_hf_hf_mmut ) {

  CQRTTEST( "rt/serial/bort/ch2o_hf_hf_mmut",
    "ch2o_hf_hf_mmut.bin.ref" );

}

// hcn b3lyp/epc17 neo bort rt dynamics
TEST( BORT_NEO_RT, hcn_b3lyp_epc17_mmut ) {

  CQRTTEST( "rt/serial/bort/hcn_b3lyp_epc17_mmut",
    "hcn_b3lyp_epc17_mmut.bin.ref", oneProtonTol );

}

// hcn hf neo bort rt dynamics
TEST( BORT_NEO_RT, hcn_hf_hf_mmut ) {

  CQRTTEST( "rt/serial/bort/hcn_hf_hf_mmut",
    "hcn_hf_hf_mmut.bin.ref", oneProtonTol );

}

#ifdef _CQ_DO_PARTESTS
// rk4 propagator
// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, PAR_ch2o_hf_hf_rk4 ) {

  CQRTTEST( "rt/parallel/bort/ch2o_hf_hf_rk4",
    "ch2o_hf_hf_rk4.bin.ref" );

}

// hcn hf neo bort rt dynamics   
TEST( BORT_NEO_RT, PAR_hcn_hf_hf_rk4 ) {

  CQRTTEST( "rt/parallel/bort/hcn_hf_hf_rk4",
    "hcn_hf_hf_rk4.bin.ref", oneProtonTol );

}



// magnus2 propagator
// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, PAR_ch2o_hf_hf_magnus2 ) {

  CQRTTEST( "rt/parallel/bort/ch2o_hf_hf_magnus2",
    "ch2o_hf_hf_magnus2.bin.ref" );

}

// hcn hf neo bort rt dynamics   
TEST( BORT_NEO_RT, PAR_hcn_hf_hf_magnus2 ) {

  CQRTTEST( "rt/parallel/bort/hcn_hf_hf_magnus2",
    "hcn_hf_hf_magnus2.bin.ref", oneProtonTol );

}



// mmut propagator
// ch2o b3lyp/epc17 neo bort rt dynamics
TEST( BORT_NEO_RT, PAR_ch2o_b3lyp_epc17_mmut ) {

  CQRTTEST( "rt/parallel/bort/ch2o_b3lyp_epc17_mmut",
    "ch2o_b3lyp_epc17_mmut.bin.ref" );

}

// ch2o hf neo bort rt dynamics
TEST( BORT_NEO_RT, PAR_ch2o_hf_hf_mmut ) {

  CQRTTEST( "rt/parallel/bort/ch2o_hf_hf_mmut",
    "ch2o_hf_hf_mmut.bin.ref" );

}

// hcn b3lyp/epc17 neo bort rt dynamics
TEST( BORT_NEO_RT, PAR_hcn_b3lyp_epc17_mmut ) {

  CQRTTEST( "rt/parallel/bort/hcn_b3lyp_epc17_mmut",
    "hcn_b3lyp_epc17_mmut.bin.ref", oneProtonTol );

}

// hcn hf neo bort rt dynamics   
TEST( BORT_NEO_RT, PAR_hcn_hf_hf_mmut ) {

  CQRTTEST( "rt/parallel/bort/hcn_hf_hf_mmut",
    "hcn_hf_hf_mmut.bin.ref", oneProtonTol );

}






#endif
