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

#include "dasscf.hpp"


// Al 6-31G(d) tests
TEST(X2C_DASSCF_FULLMATRIX, Al_631G_direct_n6 ) {

  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_x2c_dasscf_full_direct_n6", "al_6-31G_x2c_dasscf.bin.ref"); //,1e-7);
};

#ifndef _CQ_GENERATE_TESTS
TEST(X2C_DASSCF_FULLMATRIX, Al_631G_incore_n6 ) {

  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_x2c_dasscf_full_incore_n6", "al_6-31G_x2c_dasscf.bin.ref");
};
#endif

// GIAO + DASSCF
// TEST(X2C_DASSCF_GIAO, NO_631G ) {
//   CQDASSCFTEST( "dasscf/serial/cas/nox2chf_sacas", "nox2chf_sacas.bin.ref");
// 
// #ifndef _CQ_GENERATE_TESTS
//   CQDASSCFTEST( "dasscf/serial/cas/nox2chf_sacas", "nox2chf_sacas.bin.ref");
// #endif 
// };


// TEST(FourC_DASSCF_FULLMATRIX, Al_631G ) {
// 
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_bc_dasscf_full_incore_n6", "al_6-31G_4c_bc_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dc_dasscf_full_incore_n6", "al_6-31G_4c_dc_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcssss_dasscf_full_incore_n6", "al_6-31G_4c_dcssss_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcg_dasscf_full_incore_n6", "al_6-31G_4c_dcg_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcb_dasscf_full_incore_n6", "al_6-31G_4c_dcb_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/U91+_DCB_fci", "U91+_DCB_fci.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/U91+_DCB_fci_no_nagetive_rotation", "U91+_DCB_fci_no_nagetive_rotation.bin.ref");
// 
// #ifndef _CQ_GENERATE_TESTS
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_bc_dasscf_full_incore_n5", "al_6-31G_4c_bc_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dc_dasscf_full_incore_n5", "al_6-31G_4c_dc_dasscf.bin.ref");
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_bc_dasscf_full_direct_n6", "al_6-31G_4c_bc_dasscf.bin.ref"); // ,1e-7);
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dc_dasscf_full_direct_n6", "al_6-31G_4c_dc_dasscf.bin.ref"); //,1e-7);
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcssss_dasscf_full_direct_n6", "al_6-31G_4c_dcssss_dasscf.bin.ref"); //, 1e-7);
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcg_dasscf_full_direct_n6", "al_6-31G_4c_dcg_dasscf.bin.ref"); //, 1e-7);
//   CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcb_dasscf_full_direct_n6", "al_6-31G_4c_dcb_dasscf.bin.ref"); //, 1e-7);
// #endif
//  
// };

// #ifndef _CQ_GENERATE_TESTS
// // Al readmo tests - NYI: readci not yet implemented
// TEST(DASCI_READMO_SKIPSCF, Al_631G) {
//   CQDASSCFTEST("dasscf/serial/cas/al_6-31G_1c_dasci_readmo_skipscf", "al_6-31G_1c_dasscf.bin.ref", true);
//   CQDASSCFTEST("dasscf/serial/cas/al_6-31G_2c_dasci_readmo_skipscf", "al_6-31G_x2c_dasscf.bin.ref", true);
//   CQDASSCFTEST("dasscf/serial/cas/al_6-31G_4c_dasci_readmo_skipscf", "al_6-31G_4c_dcb_dasscf.bin.ref", true);
// }
// #endif

// Be swap test
TEST(TwoC_DAS_SWAP, Be_2c_SWAP_sto3G ) {

  CQDASSCFTEST( "dasscf/serial/cas/be_sto-3G_2c_dasci_swap", "be_sto-3G_2c_dasci_swap.bin.ref", true );

};

// oscillator strength test
TEST(GHF_DAS_OSC, Al_GHF_OSC_STR) {
  CQDASSCFTEST( "dasscf/serial/cas/al_ghf_6-31g_dasci_osc_str",
        "al_ghf_6-31g_dasci_osc_str.bin.ref", false, "", 1e-8, false, false, false, true, true, true, true, true, true, true, true );

}

// NYI: ethylene_* inputs have not been ported into serial/cas yet
/*
TEST(MCSCF_FIELD, Ethylene_MCSCF_W_FIELD) {
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_useSCFfield",                          "ethylene_TwoC_x2c_1DetMCSCF_wfield_useSCFfield.bin.ref", true);
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwofieldMCSCFwfield",                "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwofieldMCSCFwfield.bin.ref", true);
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwofield",                "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwofield.bin.ref", true);
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield", "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield.bin.ref", true);

  // 2 electron 2 orbital DASSCF
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_nr_2e2oDASSCF_wfield_xyzfield", "ethylene_TwoC_nr_2e2oDASSCF_wfield_xyzfield.bin.ref", true);
  CQDASSCFTEST( "dasscf/serial/cas/ethylene_TwoC_x2c_2e2oDASSCF_wfield_xyzfield",    "ethylene_TwoC_x2c_2e2oDASSCF_wfield_xyzfield.bin.ref", true);
}
*/

// CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_1c_dasci_davidson",     "al_6-31G_1c_dasci.bin.ref", true);
TEST(DASCI_DAVIDSON, Al_631G_x2c ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_x2c_dasci_davidson",       "al_6-31G_x2c_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, Al_631G_4c_bc ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_bc_dasci_davidson",     "al_6-31G_4c_bc_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, Al_631G_4c_dc ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dc_dasci_davidson",     "al_6-31G_4c_dc_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, Al_631G_4c_dcssss ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcssss_dasci_davidson", "al_6-31G_4c_dcssss_dasci.bin.ref"); };

TEST(DASCI_DAVIDSON, Al_631G_4c_dcg ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcg_dasci_davidson",    "al_6-31G_4c_dcg_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, Al_631G_4c_dcb ) {
  CQDASSCFTEST( "dasscf/serial/cas/al_6-31G_4c_dcb_dasci_davidson",    "al_6-31G_4c_dcb_dasci.bin.ref" ); };

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

// SMP Al 6-31G(d) test

TEST(X2C_DASSCF_FULLMATRIX, PAR_Al_631G_incore_n6 ) {

  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_x2c_dasscf_full_incore_n6", "al_6-31G_x2c_dasscf.bin.ref");

};

TEST(X2C_DASSCF_FULLMATRIX, PAR_Al_631G_direct_n6 ) {

  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_x2c_dasscf_full_direct_n6", "al_6-31G_x2c_dasscf.bin.ref"); //,1e-7);

};

TEST(DASCI_DAVIDSON, PAR_Al_631G_x2c ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_x2c_dasci_davidson",       "al_6-31G_x2c_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, PAR_Al_631G_4c_bc ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_4c_bc_dasci_davidson",     "al_6-31G_4c_bc_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, PAR_Al_631G_4c_dc ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_4c_dc_dasci_davidson",     "al_6-31G_4c_dc_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, PAR_Al_631G_4c_dcssss ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_4c_dcssss_dasci_davidson", "al_6-31G_4c_dcssss_dasci.bin.ref"); };

TEST(DASCI_DAVIDSON, PAR_Al_631G_4c_dcg ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_4c_dcg_dasci_davidson",    "al_6-31G_4c_dcg_dasci.bin.ref" ); };

TEST(DASCI_DAVIDSON, PAR_Al_631G_4c_dcb ) {
  CQDASSCFTEST( "dasscf/parallel/cas/al_6-31G_4c_dcb_dasci_davidson",    "al_6-31G_4c_dcb_dasci.bin.ref" ); };

TEST(TwoC_DAS_SWAP, PAR_Be_2c_SWAP_sto3G ) {

  CQDASSCFTEST( "dasscf/parallel/cas/be_sto-3G_2c_dasci_swap", "be_sto-3G_2c_dasci_swap.bin.ref", true );

};

TEST(GHF_DAS_OSC, PAR_Al_GHF_OSC_STR) {

  CQDASSCFTEST( "dasscf/parallel/cas/al_ghf_6-31g_dasci_osc_str",
        "al_ghf_6-31g_dasci_osc_str.bin.ref", false, "", 1e-8, false, false, false, true );

}

#endif
#endif

// X2C H2O sto-3g for testing READCI (from restart file) - NYI: readci not yet implemented
// TEST(X2C_DASSCF_READCI, Water_STO3G_READCI_RESTART ) {
//
//   CQDASSCFTEST( "dasscf/serial/cas/water_x2cDASSCF_sto-3g_readci", "water_x2cDASSCF_sto-3g.bin.ref",true);
//
// #ifndef _CQ_GENERATE_TESTS
//
//   // CI from restart file and MOs from scratch file
//   CQDASSCFTEST( "dasscf/serial/cas/water_x2cDASSCF_sto-3g_readci", "water_x2cDASSCF_sto-3g.bin.ref",true,"water_x2cDASSCF_sto-3g.scr.ref");
//   // CI and MOs from scratch file
//   CQDASSCFTEST( "dasscf/serial/cas/water_x2cDASSCF_sto-3g_readci", "water_x2cDASSCF_sto-3g.bin.ref",false,"water_x2cDASSCF_sto-3g.scr.ref");
//
// #endif
//
// };

// DC H2O sto-3g for testing READCI (from restart file) - NYI: readci not yet implemented
// TEST(FourC_DASSCF_READCI, Water_STO3G_READCI_RESTART ) {
//
//   CQDASSCFTEST( "dasscf/serial/cas/water_4cDASSCF_sto-3g_readci", "water_4cDASSCF_sto-3g.bin.ref",true);
//
// #ifndef _CQ_GENERATE_TESTS
//
//   // CI from restart file and MOs from scratch file
//   CQDASSCFTEST( "dasscf/serial/cas/water_4cDASSCF_sto-3g_readci", "water_4cDASSCF_sto-3g.bin.ref",true,"water_4cDASSCF_sto-3g.scr.ref");
//   // CI and MOs from scratch file
//   CQDASSCFTEST( "dasscf/serial/cas/water_4cDASSCF_sto-3g_readci", "water_4cDASSCF_sto-3g.bin.ref",false,"water_4cDASSCF_sto-3g.scr.ref");
//
// #endif
//
// };
