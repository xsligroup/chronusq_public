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

#include "coupledcluster.hpp"


//TEST( EOM_CCSD, H2O_631G_GHF_EOM_CCSD) {
//  CQCCTEST("coupledcluster/serial/eom-ccsd/h2o_631g_ghf_eom_ccsd",
//    "h2o_631g_ghf_eom_ccsd.bin.ref", "", true, true, true, false);
//}
//
//TEST( EOM_CCSD, Na_plus_STO3G_X2C_EOM_CCSD) {
//  CQCCTEST("coupledcluster/serial/eom-ccsd/na_plus_sto3g_x2c_eom_ccsd",
//           "na_plus_sto3g_x2c_eom_ccsd.bin.ref", "", true, true, true, false);
//}

TEST( EOM_CCSD, Na_631G_X2C_EOM_CCSD) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/na_631G_x2c_eom_ccsd",
           "na_631G_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, Na_631G_X2C_EOM_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/na_631G_x2c_eom_ccsd_withdenomshift",
           "na_631G_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, Na_631G_CationRef_X2C_CVS_EOM_CCSD) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/na_631G_cationRef_x2c_eom_ccsd",
           "na_631G_cationRef_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, Na_631G_CationRef_X2C_CVS_EOM_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/na_631G_cationRef_x2c_eom_ccsd_withdenomshift",
           "na_631G_cationRef_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, H2O_STO3G_GHF_CVS) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/h2o_sto3g_ghf_cvs",
           "h2o_sto3g_ghf_cvs.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, H2O_STO3G_GHF_CVS_DENOMSHIFT) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/h2o_sto3g_ghf_cvs_withdenomshift",
           "h2o_sto3g_ghf_cvs.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, BE_631G_GHF_CVS_R_RESTART) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/be_631g_ghf_cvs_r_restart",
           "be_631g_ghf_cvs_r_restart.bin.ref","be_631g_ghf_cvs.bin.ref", true, true, false, false, false, false);
}

TEST( EOM_CCSD, NA_PLUS_631G_DIP_CCSGUESS) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/na_plus_631g_dip_ccsguess",
           "na_plus_631g_dip_ccsguess.bin.ref","", true, true, false, false, false, false);
}

TEST( EOM_CCSD, H2O_STO3G_GHF_CVS_L_RESTART) {
  CQCCTEST("coupledcluster/serial/eom-ccsd/h2o_sto3g_ghf_cvs_l_restart",
           "h2o_sto3g_ghf_cvs_l_restart.bin.ref","h2o_sto3g_ghf_cvs.bin.ref", true, true, true, false, false, false);
}

#ifdef _CQ_DO_PARTESTS

//TEST( EOM_CCSD, PAR_H2O_631G_GHF_EOM_CCSD) {
//  CQCCTEST("coupledcluster/parallel/eom-ccsd/h2o_631g_ghf_eom_ccsd",
//           "h2o_631g_ghf_eom_ccsd.bin.ref", "", true, true, true, false);
//}
//
//TEST( EOM_CCSD, PAR_Na_plus_STO3G_X2C_EOM_CCSD) {
//  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_plus_sto3g_x2c_eom_ccsd",
//           "na_plus_sto3g_x2c_eom_ccsd.bin.ref", "", true, true, true, false);
//}

TEST( EOM_CCSD, PAR_Na_631G_X2C_EOM_CCSD) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_631G_x2c_eom_ccsd",
           "na_631G_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_Na_631G_X2C_EOM_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_631G_x2c_eom_ccsd_withdenomshift",
           "na_631G_x2c_eom_ccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_Na_631G_CationRef_X2C_CVS_EOM_CCSD) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_631G_cationRef_x2c_eom_ccsd",
           "na_631G_cationRef_x2c_eom_ccsd.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_Na_631G_CationRef_X2C_CVS_EOM_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_631G_cationRef_x2c_eom_ccsd_withdenomshift",
           "na_631G_cationRef_x2c_eom_ccsd.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_H2O_STO3G_GHF_CVS) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/h2o_sto3g_ghf_cvs",
           "h2o_sto3g_ghf_cvs.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_H2O_STO3G_GHF_CVS_DENOMSHIFT) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/h2o_sto3g_ghf_cvs_withdenomshift",
           "h2o_sto3g_ghf_cvs.bin.ref","", true, true, true, false, false, false);
}

TEST( EOM_CCSD, PAR_BE_631G_GHF_CVS_R_RESTART) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/be_631g_ghf_cvs_r_restart",
           "be_631g_ghf_cvs_r_restart.bin.ref","be_631g_ghf_cvs.bin.ref", true, true, false, false, false, false);
}

TEST( EOM_CCSD, PAR_NA_PLUS_631G_DIP_CCSGUESS) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/na_plus_631g_dip_ccsguess",
           "na_plus_631g_dip_ccsguess.bin.ref","", true, true, false, false, false, false);
}

TEST( EOM_CCSD, PAR_H2O_STO3G_GHF_CVS_L_RESTART) {
  CQCCTEST("coupledcluster/parallel/eom-ccsd/h2o_sto3g_ghf_cvs_l_restart",
           "h2o_sto3g_ghf_cvs_l_restart.bin.ref","h2o_sto3g_ghf_cvs.bin.ref", true, true, true, false, false, false);
}

#endif

