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

#include "coupledcluster.hpp"


TEST( CCSDT, H2O_631G_GHF_PT_CR_IJK) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_ghf_ccsdpt_crcc23_loopijk",
           "h2o_6-31g_ghf_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, H2O_631G_GHF_PT_CR_ABC) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_ghf_ccsdpt_crcc23_loopabc",
           "h2o_6-31g_ghf_ccsdpt_crcc23_loopabc.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, H2O_631G_X2C_PT_CR_IJK) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopijk",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, H2O_631G_X2C_PT_CR_IJK_DENOMSHIFT) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopijk_withdenomshift",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, H2O_631G_X2C_PT_CR_ABC) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopabc",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopabc.bin.ref", "", false, false, false, false, true, false);
}

// Stephen: EAIP contains DIP-EOMCCSDT for H4/6-31G, which is checked against CCPy.
TEST( CCSDT, H2O_631G_GHF_CCSDT) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_ghf_ccsdt",
           "h2o_6-31g_ghf_ccsdt.bin.ref", "", false, false, false, false, false, false);
}

TEST( CCSDT, H2O_631G_X2C_CCSDT) {
  CQCCTEST("coupledcluster/serial/ccsdt/h2o_6-31g_x2c_ccsdt",
           "h2o_6-31g_x2c_ccsdt.bin.ref", "", false, false, false, false, false, false);
}

#ifdef _CQ_DO_PARTESTS

TEST( CCSDT, PAR_H2O_631G_GHF_PT_CR_IJK) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_ghf_ccsdpt_crcc23_loopijk",
           "h2o_6-31g_ghf_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, PAR_H2O_631G_GHF_PT_CR_ABC) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_ghf_ccsdpt_crcc23_loopabc",
           "h2o_6-31g_ghf_ccsdpt_crcc23_loopabc.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, PAR_H2O_631G_X2C_PT_CR_IJK) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopijk",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, PAR_H2O_631G_X2C_PT_CR_IJK_DENOMSHIFT) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopijk_withdenomshift",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopijk.bin.ref", "", false, false, false, false, true, false);
}

TEST( CCSDT, PAR_H2O_631G_X2C_PT_CR_ABC) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_x2c_ccsdpt_crcc23_loopabc",
           "h2o_6-31g_x2c_ccsdpt_crcc23_loopabc.bin.ref", "", false, false, false, false, true, false);
}

// Stephen: EAIP contains DIP-EOMCCSDT for H4/6-31G, which is checked against CCPy.
TEST( CCSDT, PAR_H2O_631G_GHF_CCSDT) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_ghf_ccsdt",
           "h2o_6-31g_ghf_ccsdt.bin.ref", "", false, false, false, false, false, false);
}

TEST( CCSDT, PAR_H2O_631G_X2C_CCSDT) {
  CQCCTEST("coupledcluster/parallel/ccsdt/h2o_6-31g_x2c_ccsdt",
           "h2o_6-31g_x2c_ccsdt.bin.ref", "", false, false, false, false, false, false);
}

#endif

