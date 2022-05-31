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

  CQMCSCFREFTEST( "mcscf/serial/cas/al_6-31G_1c_casscf_full_incore_n6", "al_6-31G_1c_casscf.bin.ref");

#ifndef _CQ_GENERATE_TESTS
  CQMCSCFREFTEST( "mcscf/serial/cas/al_6-31G_1c_casscf_full_direct_n6", "al_6-31G_1c_casscf.bin.ref");
#endif

};

TEST(X2C_CASSCF_FULLMATRIX, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_incore_n6", "al_6-31G_x2c_casscf.bin.ref");
 
#ifndef _CQ_GENERATE_TESTS
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_incore_n5", "al_6-31G_x2c_casscf.bin.ref");
  CQMCSCFTEST( "mcscf/serial/cas/al_6-31G_x2c_casscf_full_direct_n6", "al_6-31G_x2c_casscf.bin.ref"); //,1e-7);
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
  CQMCSCFREFTEST("mcscf/serial/cas/al_6-31G_1c_casci_readmo_skipscf", "al_6-31G_1c_casscf.bin.ref"); 
  CQMCSCFREFTEST("mcscf/serial/cas/al_6-31G_2c_casci_readmo_skipscf", "al_6-31G_x2c_casscf.bin.ref"); 
  CQMCSCFREFTEST("mcscf/serial/cas/al_6-31G_4c_casci_readmo_skipscf", "al_6-31G_4c_dcb_casscf.bin.ref"); 
}
#endif

// Be swap test
TEST(OneC_CAS_SWAP, Be_1c_SWAP_sto3G ) {

  CQMCSCFREFTEST( "mcscf/serial/cas/be_sto-3G_1c_casci_swap", "be_sto-3G_1c_casci_swap.bin.ref" );

};

TEST(TwoC_CAS_SWAP, Be_2c_SWAP_sto3G ) {

  CQMCSCFREFTEST( "mcscf/serial/cas/be_sto-3G_2c_casci_swap", "be_sto-3G_2c_casci_swap.bin.ref" );

};

// oscillator strength test
TEST(GHF_CAS_OSC, Al_GHF_OSC_STR) {

  CQMCSCFTEST( "mcscf/serial/cas/al_ghf_6-31g_casci_osc_str",
        "al_ghf_6-31g_casci_osc_str.bin.ref", 1e-8, true );

}

TEST(CASCI_DAVIDSON, Al_631G ) {

  CQMCSCFREFTEST( "mcscf/serial/cas/al_6-31G_1c_casci_davidson",     "al_6-31G_1c_casci.bin.ref");
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

  CQMCSCFREFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_incore_n6", "al_6-31G_1c_casscf.bin.ref");
  CQMCSCFREFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_incore_n5", "al_6-31G_1c_casscf.bin.ref");
  CQMCSCFREFTEST( "mcscf/parallel/cas/al_6-31G_1c_casscf_full_direct_n6", "al_6-31G_1c_casscf.bin.ref");

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
  CQMCSCFREFTEST("mcscf/parallel/cas/al_6-31G_1c_casci_readmo_skipscf", "al_6-31G_1c_casscf.bin.ref"); 
  CQMCSCFREFTEST("mcscf/parallel/cas/al_6-31G_2c_casci_readmo_skipscf", "al_6-31G_x2c_casscf.bin.ref"); 
  CQMCSCFREFTEST("mcscf/parallel/cas/al_6-31G_4c_casci_readmo_skipscf", "al_6-31G_4c_dcb_casscf.bin.ref"); 
}

TEST(CASCI_DAVIDSON, PAR_Al_631G) {

  CQMCSCFREFTEST( "mcscf/parallel/cas/al_6-31G_1c_casci_davidson", "al_6-31G_1c_casci.bin.ref");
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_x2c_casci_davidson",  "al_6-31G_x2c_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_bc_casci_davidson",  "al_6-31G_4c_bc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dc_casci_davidson",  "al_6-31G_4c_dc_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcssss_casci_davidson",  "al_6-31G_4c_dcssss_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcg_casci_davidson",     "al_6-31G_4c_dcg_casci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/cas/al_6-31G_4c_dcb_casci_davidson",     "al_6-31G_4c_dcb_casci.bin.ref" );
 
};

#endif
#endif


