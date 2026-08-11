/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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

#include "scf.hpp"

// Test regenerated after NEO Refactor. Benchmarked against Gaussian:
//h2o_sto3g_protsp_b3lyp_epc17_99590.log: SCF Done:  E(RB3LYP+NEO-LENC) =  -75.2314904855     A.U. after   14 cycles
//h2o_sto3g_protsp_b3lyp_epc19_99590.log: SCF Done:  E(RB3LYP+NEO-LYPENC) =  -75.2370104938     A.U. after   14 cycles

//NEO-DFT with minimal basis set, using epc17 functional
TEST( NEO_RKS, water_sto3g_protsp_rb3lyp_uepc17) {

  CQNEOSCFTEST( "scf/serial/neo_rks/water_sto-3g_prot-sp_rb3lyp_uepc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref" );
 
}

TEST( NEO_RKS, water_stepwise_epc17) {

  CQNEOSCFTEST( "scf/serial/neo_rks/water_stepwise_epc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref",  1e-6,
              true, true, true, true, false, "no", true, false);
 
}

//NEO-DFT with minimal basis set, using epc19 functional
TEST( NEO_RKS, water_sto3g_protsp_rb3lyp_uepc19) {

  CQNEOSCFTEST( "scf/serial/neo_rks/water_sto-3g_prot-sp_rb3lyp_uepc19", "water_sto-3g_prot-sp_rb3lyp_uepc19.bin.ref" );
 
}

// H2O: electronic R-HSE06/cc-pVDZ with U-EPC17/prot-pb4-d.
TEST( NEO_RKS, h2o_ccpvdz_pb4d_rhse06_uepc17_incore) {

  CQNEOSCFTEST( "scf/serial/neo_rks/h2o_cc-pVDZ_pb4d_rhse06_uepc17",
    "h2o_cc-pVDZ_pb4d_rhse06_uepc17.bin.ref" );

}

// HCN: electronic R-CAM-B3LYP/cc-pVDZ with U-EPC17/prot-pb4-d.
TEST( NEO_RKS, hcn_ccpvdz_pb4d_rcamb3lyp_uepc17_direct) {

  CQNEOSCFTEST( "scf/serial/neo_rks/hcn_cc-pVDZ_pb4d_rcamb3lyp_uepc17",
    "hcn_cc-pVDZ_pb4d_rcamb3lyp_uepc17.bin.ref" );

}

TEST( NEO_RKS, hcn_ccpvdz_pb4d_rcamb3lyp_uepc17_stepwise) {

  CQNEOSCFTEST( "scf/serial/neo_rks/hcn_cc-pVDZ_pb4d_rcamb3lyp_uepc17_stepwise",
    "hcn_cc-pVDZ_pb4d_rcamb3lyp_uepc17.bin.ref" );

}

// COH2: electronic R-wB97X/cc-pVDZ with U-EPC17/prot-pb4-d.
TEST( NEO_RKS, coh2_ccpvdz_pb4d_rwb97x_uepc17_incore) {

  CQNEOSCFTEST( "scf/serial/neo_rks/coh2_cc-pVDZ_pb4d_rwb97x_uepc17",
    "coh2_cc-pVDZ_pb4d_rwb97x_uepc17.bin.ref" );

}


#ifdef _CQ_DO_PARTESTS

//NEO-DFT with minimal basis set, using epc17 functional, parallel job
TEST( NEO_RKS, par_water_sto3g_protsp_rb3lyp_uepc17) {

  CQNEOSCFTEST( "scf/parallel/neo_rks/par_water_sto-3g_prot-sp_rb3lyp_uepc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref" );
 
}

TEST( NEO_RKS, par_water_stepwise_epc17) {

  CQNEOSCFTEST( "scf/parallel/neo_rks/water_stepwise_epc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref",  1e-6,
              true, true, true, true, false, "no", true, false);
 
}


// FIXME: for epc19, parallel jobs doesn't pass. Need to debug epc19
//NEO-DFT with minimal basis set, using epc19 functional, parallel job
TEST( NEO_RKS, par_water_sto3g_protsp_rb3lyp_uepc19) {

  CQNEOSCFTEST( "scf/parallel/neo_rks/par_water_sto-3g_prot-sp_rb3lyp_uepc19", "water_sto-3g_prot-sp_rb3lyp_uepc19.bin.ref" );
 
}

#endif

