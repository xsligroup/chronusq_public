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

#include "mcscf.hpp"


// Al 6-31G(d) test
TEST(OneC_CASSCF_FULLMATRIX, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_1c_casscf_full_incore_n6", "al_6-31G_1c_casscf.bin.ref",true);

#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_1c_casscf_full_direct_n6", "al_6-31G_1c_casscf.bin.ref",true);
#endif

};

TEST(X2C_CASSCF_FULLMATRIX, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_incore_n6", "al_6-31G_x2c_casscf.bin.ref");
 
#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_incore_n5", "al_6-31G_x2c_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_direct_n6", "al_6-31G_x2c_casscf.bin.ref"); //,1e-7);
#endif
};

// GIAO + CASSCF
TEST(X2C_CASSCF_GIAO, NO_631G ) {
  CQMCSCFTEST( "mcscf/serial/cas/nox2chf_sacas", "nox2chf_sacas.bin.ref");

#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/nox2chf_sacas", "nox2chf_sacas.bin.ref");
#endif 
};


TEST(FourC_CASSCF_FULLMATRIX, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_bc_casscf_full_incore_n6", "al_6-31G_4c_bc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dc_casscf_full_incore_n6", "al_6-31G_4c_dc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcssss_casscf_full_incore_n6", "al_6-31G_4c_dcssss_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcg_casscf_full_incore_n6", "al_6-31G_4c_dcg_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcb_casscf_full_incore_n6", "al_6-31G_4c_dcb_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/U91+_DCB_fci", "U91+_DCB_fci.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/U91+_DCB_fci_no_nagetive_rotation", "U91+_DCB_fci_no_nagetive_rotation.bin.ref");

#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_bc_casscf_full_incore_n5", "al_6-31G_4c_bc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dc_casscf_full_incore_n5", "al_6-31G_4c_dc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_bc_casscf_full_direct_n6", "al_6-31G_4c_bc_casscf.bin.ref"); // ,1e-7);
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dc_casscf_full_direct_n6", "al_6-31G_4c_dc_casscf.bin.ref"); //,1e-7);
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcssss_casscf_full_direct_n6", "al_6-31G_4c_dcssss_casscf.bin.ref"); //, 1e-7);
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcg_casscf_full_direct_n6", "al_6-31G_4c_dcg_casscf.bin.ref"); //, 1e-7);
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcb_casscf_full_direct_n6", "al_6-31G_4c_dcb_casscf.bin.ref"); //, 1e-7);
#endif
 
};

#ifndef _CQ_GENERATE_TESTS
// Al readmo tests
TEST(CASCI_READMO_SKIPSCF, Al_631G) {
  CQMCSCFTEST("mcscf/serial/cas/al_6-31G_1c_casci_readmo_skipscf", "al_6-31G_1c_casscf.bin.ref", true);
  CQMCSCFTEST("mcscf/serial/cas/al_6-31G_2c_casci_readmo_skipscf", "al_6-31G_x2c_casscf.bin.ref", true);
  CQMCSCFTEST("mcscf/serial/cas/al_6-31G_4c_casci_readmo_skipscf", "al_6-31G_4c_dcb_casscf.bin.ref", true);
}
#endif

// Be swap test
TEST(OneC_CAS_SWAP, Be_1c_SWAP_sto3G ) {

  CQMCSCFTEST( "mcscf/serial/cas/be_sto-3G_1c_casci_swap", "be_sto-3G_1c_casci_swap.bin.ref", true );

};

TEST(TwoC_CAS_SWAP, Be_2c_SWAP_sto3G ) {

  CQMCSCFTEST( "mcscf/serial/cas/be_sto-3G_2c_casci_swap", "be_sto-3G_2c_casci_swap.bin.ref", true );

};

// oscillator strength test
TEST(GHF_CAS_OSC, Al_GHF_OSC_STR) {

  CQMCSCFTEST( "mcscf/serial/cas/al_ghf_6-31g_casci_osc_str",
        "al_ghf_6-31g_casci_osc_str.bin.ref", false, "", 1e-8, false, false, false, true );

}

TEST(MCSCF_FIELD, Ethylene_MCSCF_W_FIELD) {
  // 1 DET
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_OneC_nr_1DetMCSCF_wfield_useSCFfield",                          "ethylene_OneC_nr_1DetMCSCF_wfield_useSCFfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_OneC_nr_1DetMCSCF_wfield_SCFwofieldMCSCFwfield",                "ethylene_OneC_nr_1DetMCSCF_wfield_SCFwofieldMCSCFwfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_OneC_nr_1DetMCSCF_wfield_SCFwfieldMCSCFwofield",                "ethylene_OneC_nr_1DetMCSCF_wfield_SCFwfieldMCSCFwofield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_OneC_nr_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield", "ethylene_OneC_nr_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield.bin.ref", true);

  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_useSCFfield",                          "ethylene_TwoC_x2c_1DetMCSCF_wfield_useSCFfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwofieldMCSCFwfield",                "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwofieldMCSCFwfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwofield",                "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwofield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield", "ethylene_TwoC_x2c_1DetMCSCF_wfield_SCFwfieldMCSCFwfield_dontusescffield.bin.ref", true);

  // 2 electron 2 orbital CASSCF
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_OneC_nr_2e2oCASSCF_wfield_xyzfield", "ethylene_OneC_nr_2e2oCASSCF_wfield_xyzfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_nr_2e2oCASSCF_wfield_xyzfield", "ethylene_TwoC_nr_2e2oCASSCF_wfield_xyzfield.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/ethylene_TwoC_x2c_2e2oCASSCF_wfield_xyzfield",    "ethylene_TwoC_x2c_2e2oCASSCF_wfield_xyzfield.bin.ref", true);
}

// H2O 6-31G(d) test (for multipole properties)
TEST(OneC_CASSCF_FULLMATRIX, Water_631Gd ) {

  CQMCSCFTEST( "mcscf/serial/cas/water_1cCASSCF_6-31Gd", "water_1cCASSCF_6-31Gd.bin.ref",false,"",1e-6,false,true,true,false,true);

#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/water_1cCASSCF_6-31Gd", "water_1cCASSCF_6-31Gd.bin.ref",false,"",1e-6,false,true,true,false,true);
#endif

};

TEST(CASCI_DAVIDSON, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_1c_casci_davidson",     "al_6-31G_1c_casci.bin.ref", true);
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casci_davidson",       "al_6-31G_x2c_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_bc_casci_davidson",     "al_6-31G_4c_bc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dc_casci_davidson",     "al_6-31G_4c_dc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcssss_casci_davidson", "al_6-31G_4c_dcssss_casci.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcg_casci_davidson",    "al_6-31G_4c_dcg_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_4c_dcb_casci_davidson",    "al_6-31G_4c_dcb_casci.bin.ref" );
 
};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

// SMP Al 6-31G(d) test
TEST(OneC_CASSCF_FULLMATRIX, PAR_Al_631G ) {

  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_incore_n6", "al_6-31G_1c_casscf.bin.ref", true);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_incore_n5", "al_6-31G_1c_casscf.bin.ref", true);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_direct_n6", "al_6-31G_1c_casscf.bin.ref", true);

};

TEST(OneC_CASSCF_FULLMATRIX, PAR_Water_631Gd ) {

  CQMCSCFTEST( "mcscf/parallel/cas/water_1cCASSCF_6-31Gd", "water_1cCASSCF_6-31Gd.bin.ref",false,"",1e-6,false,true,true,false,true);

};

TEST(X2C_CASSCF_FULLMATRIX, PAR_Al_631G ) {

  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_x2c_casscf_full_incore_n6", "al_6-31G_x2c_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_x2c_casscf_full_incore_n5", "al_6-31G_x2c_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_x2c_casscf_full_direct_n6", "al_6-31G_x2c_casscf.bin.ref"); //,1e-7);
 
};

TEST(FourC_CASSCF_FULLMATRIX, PAR_Al_631G ) {

  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_bc_casscf_full_incore_n6", "al_6-31G_4c_bc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dc_casscf_full_incore_n6", "al_6-31G_4c_dc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcssss_casscf_full_incore_n6", "al_6-31G_4c_dcssss_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcg_casscf_full_incore_n6", "al_6-31G_4c_dcg_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcb_casscf_full_incore_n6", "al_6-31G_4c_dcb_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_bc_casscf_full_incore_n5", "al_6-31G_4c_bc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dc_casscf_full_incore_n5", "al_6-31G_4c_dc_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_bc_casscf_full_direct_n6", "al_6-31G_4c_bc_casscf.bin.ref"); //,1e-7);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dc_casscf_full_direct_n6", "al_6-31G_4c_dc_casscf.bin.ref"); //,1e-7);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcssss_casscf_full_direct_n6", "al_6-31G_4c_dcssss_casscf.bin.ref"); //,1e-7);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcg_casscf_full_direct_n6", "al_6-31G_4c_dcg_casscf.bin.ref"); //,1e-7);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcb_casscf_full_direct_n6", "al_6-31G_4c_dcb_casscf.bin.ref"); //,1e-7);
 
};

// Al readmo tests
TEST(CASCI_READMO_SKIPSCF, PAR_Al_631G) {
  CQMCSCFTEST("mcscf/parallel/cas/al_6-31G_1c_casci_readmo_skipscf", "al_6-31G_1c_casscf.bin.ref", true);
  CQMCSCFTEST("mcscf/parallel/cas/al_6-31G_2c_casci_readmo_skipscf", "al_6-31G_x2c_casscf.bin.ref", true);
  CQMCSCFTEST("mcscf/parallel/cas/al_6-31G_4c_casci_readmo_skipscf", "al_6-31G_4c_dcb_casscf.bin.ref", true);
}

TEST(CASCI_DAVIDSON, PAR_Al_631G) {

  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_1c_casci_davidson", "al_6-31G_1c_casci.bin.ref", true);
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_x2c_casci_davidson",  "al_6-31G_x2c_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_bc_casci_davidson",  "al_6-31G_4c_bc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dc_casci_davidson",  "al_6-31G_4c_dc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcssss_casci_davidson",  "al_6-31G_4c_dcssss_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcg_casci_davidson",     "al_6-31G_4c_dcg_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcb_casci_davidson",     "al_6-31G_4c_dcb_casci.bin.ref" );
 
};

#endif
#endif


// 1-component H2O sto-3g for testing READCI (from restart file)
TEST(OneC_CASSCF_READCI, Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/serial/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",true);

#ifndef _CQ_GENERATE_TESTS

  // CI from restart file and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",true,"water_1cCASSCF_sto-3g.scr.ref");
  // CI and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",false,"water_1cCASSCF_sto-3g.scr.ref");

#endif

};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

TEST(OneC_CASSCF_READCI, PAR_Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/parallel/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",true);

  CQMCSCFTEST( "mcscf/parallel/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",true,"water_1cCASSCF_sto-3g.scr.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/water_1cCASSCF_sto-3g_readci", "water_1cCASSCF_sto-3g.bin.ref",false,"water_1cCASSCF_sto-3g.scr.ref");

};

#endif
#endif

// X2C H2O sto-3g for testing READCI (from restart file)
TEST(X2C_CASSCF_READCI, Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/serial/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",true);

#ifndef _CQ_GENERATE_TESTS

  // CI from restart file and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",true,"water_x2cCASSCF_sto-3g.scr.ref");
  // CI and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",false,"water_x2cCASSCF_sto-3g.scr.ref");

#endif

};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

TEST(X2C_CASSCF_READCI, PAR_Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/parallel/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",true);

  CQMCSCFTEST( "mcscf/parallel/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",true,"water_x2cCASSCF_sto-3g.scr.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/water_x2cCASSCF_sto-3g_readci", "water_x2cCASSCF_sto-3g.bin.ref",false,"water_x2cCASSCF_sto-3g.scr.ref");

};

#endif
#endif

// DC H2O sto-3g for testing READCI (from restart file)
TEST(FourC_CASSCF_READCI, Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/serial/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",true);

#ifndef _CQ_GENERATE_TESTS

  // CI from restart file and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",true,"water_4cCASSCF_sto-3g.scr.ref");
  // CI and MOs from scratch file
  CQMCSCFTEST( "mcscf/serial/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",false,"water_4cCASSCF_sto-3g.scr.ref");

#endif

};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

TEST(FourC_CASSCF_READCI, PAR_Water_STO3G_READCI_RESTART ) {

  CQMCSCFTEST( "mcscf/parallel/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",true);

  CQMCSCFTEST( "mcscf/parallel/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",true,"water_4cCASSCF_sto-3g.scr.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/water_4cCASSCF_sto-3g_readci", "water_4cCASSCF_sto-3g.bin.ref",false,"water_4cCASSCF_sto-3g.scr.ref");

};

#endif
#endif
