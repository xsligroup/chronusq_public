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

#include "dynamics.hpp"



// HF PBEXPBEC-D3 BOMD
TEST( D3, hf_bomd_pbexpbec_d3 ) {

  CQDYNAMICSTEST( "dynamics/serial/d3/hf_bomd_pbexpbec-d3",
    "hf_bomd_pbexpbec-d3.bin.ref");

}

// HF PBEXPBEC-D3 Ehrenfest
TEST( D3, hf_ehrenfest_pbexpbec_d3 ) {

  CQDYNAMICSTEST( "dynamics/serial/d3/hf_ehrenfest_pbexpbec-d3",
    "hf_ehrenfest_pbexpbec-d3.bin.ref");

}

// H2O NEO B3LYP-D3/EPC17 Ehrenfest
TEST( D3, h2o_neoehrenfest_b3lyp_d3_epc17 ) {

  CQDYNAMICSTEST( "dynamics/serial/d3/h2o_neoehrenfest_b3lyp-d3_epc17",
    "h2o_neoehrenfest_b3lyp-d3_epc17.bin.ref");

}




#ifdef _CQ_DO_PARTESTS
TEST( D3, PAR_hf_bomd_pbexpbec_d3 ) {

  CQDYNAMICSTEST( "dynamics/parallel/d3/hf_bomd_pbexpbec-d3",
    "hf_bomd_pbexpbec-d3.bin.ref");

}

TEST( D3, PAR_hf_ehrenfest_pbexpbec_d3 ) {

  CQDYNAMICSTEST( "dynamics/parallel/d3/hf_ehrenfest_pbexpbec-d3",
    "hf_ehrenfest_pbexpbec-d3.bin.ref");

}

TEST( D3, PAR_h2o_neoehrenfest_b3lyp_d3_epc17 ) {

  CQDYNAMICSTEST( "dynamics/parallel/d3/h2o_neoehrenfest_b3lyp-d3_epc17",
    "h2o_neoehrenfest_b3lyp-d3_epc17.bin.ref");

}

#endif
