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



// Classical Hydrogen Flouride RHF
TEST( EHRENFEST_DYNAMICS_GIAO, hf_ehrenfest_rhf_mmut_giao ) {

  CQDYNAMICSTEST( "dynamics/serial/ehrenfest_giao/hf_ehrenfest_rhf_mmut_giao",
    "hf_ehrenfest_rhf_mmut_giao.bin.ref");

}
TEST( EHRENFEST_DYNAMICS_GIAO, h2_ehrenfest_ghf_mmut_giao ) {

  CQDYNAMICSTEST( "dynamics/serial/ehrenfest_giao/h2_ehrenfest_ghf_mmut_giao",
    "h2_ehrenfest_ghf_mmut_giao.bin.ref");

}

// NEO Water RHF/UHF Fixed Proton Basis
TEST( EHRENFEST_DYNAMICS_GIAO, h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao ) {

  CQDYNAMICSTEST( "dynamics/serial/ehrenfest_giao/h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao",
    "h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao.bin.ref");

}

// NEO Water RB3LYP/UEPC17 Traveling Proton Basis
TEST( EHRENFEST_DYNAMICS_GIAO, h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao ) {

  CQDYNAMICSTEST( "dynamics/serial/ehrenfest_giao/h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao",
    "h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao.bin.ref");

}



#ifdef _CQ_DO_PARTESTS
TEST( EHRENFEST_DYNAMICS_GIAO, PAR_hf_ehrenfest_rhf_mmut_giao ) {

  CQDYNAMICSTEST( "dynamics/parallel/ehrenfest_giao/hf_ehrenfest_rhf_mmut_giao",
    "hf_ehrenfest_rhf_mmut_giao.bin.ref" );

}
TEST( EHRENFEST_DYNAMICS_GIAO, PAR_h2_ehrenfest_ghf_mmut_giao ) {

  CQDYNAMICSTEST( "dynamics/parallel/ehrenfest_giao/h2_ehrenfest_ghf_mmut_giao",
    "h2_ehrenfest_ghf_mmut_giao.bin.ref");

}

TEST( EHRENFEST_DYNAMICS_GIAO, PAR_h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao ) {

  CQDYNAMICSTEST( "dynamics/parallel/ehrenfest_giao/h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao",
    "h2o_neoehrenfest_rhf_uhf_mmut_fpb_giao.bin.ref" );

}


TEST( EHRENFEST_DYNAMICS_GIAO, PAR_h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao ) {

  CQDYNAMICSTEST( "dynamics/parallel/ehrenfest_giao/h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao",
    "h2o_neoehrenfest_rb3lyp_epc17_mmut_tpb_giao.bin.ref" );

}
#endif
