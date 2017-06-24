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

#include "scf.hpp"

//========================
// BFGS TESTS
//========================
// Water RHF/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_RHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rhf_bfgs", "water_6-31Gd_qn_rhf.bin.ref" );
 
};

// Water UHF/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_UHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_uhf_bfgs", "water_6-31Gd_qn_uhf.bin.ref" );
 
};

// Water GHF/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_GHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ghf_bfgs", "water_6-31Gd_qn_ghf.bin.ref" );
 
};

// Water X2CHF/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_X2C ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_x2chf_bfgs", "water_6-31Gd_qn_x2chf.bin.ref" );
 
};

// Water 4CHF/6-31G(d) test
TEST( FOURC_QN_SCF_BFGS, Water_631Gd_4CHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_4chf_bfgs", "water_6-31Gd_qn_4chf.bin.ref" );
 
};

// Water RBLYP/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_RBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rblyp_bfgs", "water_6-31Gd_qn_rblyp.bin.ref" );
 
};

// Water UBLYP/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_UBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ublyp_bfgs", "water_6-31Gd_qn_ublyp.bin.ref" );
 
};

// Water GBLYP/6-31G(d) test
TEST( QN_SCF_BFGS, Water_631Gd_GBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_gblyp_bfgs", "water_6-31Gd_qn_gblyp.bin.ref" );
 
};

//========================
// SR1 TESTS
//========================
// Water RHF/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_RHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rhf_sr1", "water_6-31Gd_qn_rhf.bin.ref" );
 
};

// Water UHF/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_UHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_uhf_sr1", "water_6-31Gd_qn_uhf.bin.ref" );
 
};

// Water GHF/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_GHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ghf_sr1", "water_6-31Gd_qn_ghf.bin.ref" );
 
};

// Water X2CHF/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_X2C ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_x2chf_sr1", "water_6-31Gd_qn_x2chf.bin.ref" );
 
};

// Water 4CHF/6-31G(d) test
TEST( FOURC_QN_SCF_SR1, Water_631Gd_4CHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_4chf_sr1", "water_6-31Gd_qn_4chf.bin.ref" );
 
};

// Water RBLYP/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_RBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rblyp_sr1", "water_6-31Gd_qn_rblyp.bin.ref" );
 
};

// Water UBLYP/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_UBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ublyp_sr1", "water_6-31Gd_qn_ublyp.bin.ref" );
 
};

// Water GBLYP/6-31G(d) test
TEST( QN_SCF_SR1, Water_631Gd_GBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_gblyp_sr1", "water_6-31Gd_qn_gblyp.bin.ref" );
 
};

//========================
// Gradient Descent TESTS
//========================
// Water RHF/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_RHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rhf_grad", "water_6-31Gd_qn_rhf.bin.ref" );
 
};

// Water UHF/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_UHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_uhf_grad", "water_6-31Gd_qn_uhf.bin.ref" );
 
};

// Water GHF/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_GHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ghf_grad", "water_6-31Gd_qn_ghf.bin.ref" );
 
};

// Water X2CHF/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_X2C ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_x2chf_grad", "water_6-31Gd_qn_x2chf.bin.ref" );
 
};

// Water 4CHF/6-31G(d) test
TEST( FOURC_QN_SCF_GradDescent, Water_631Gd_4CHF ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_4chf_grad", "water_6-31Gd_qn_4chf.bin.ref" );
 
};

// Water RBLYP/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_RBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_rblyp_grad", "water_6-31Gd_qn_rblyp.bin.ref" );
 
};

// Water UBLYP/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_UBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_ublyp_grad", "water_6-31Gd_qn_ublyp.bin.ref" );
 
};

// Water GBLYP/6-31G(d) test
TEST( QN_SCF_GradDescent, Water_631Gd_GBLYP ) {

  CQSCFTEST( "scf/serial/qn/water_6-31Gd_gblyp_grad", "water_6-31Gd_qn_gblyp.bin.ref" );
 
};

#ifdef _CQ_DO_PARTESTS

#endif



