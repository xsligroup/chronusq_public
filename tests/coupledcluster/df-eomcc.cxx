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

TEST( DF_EOMCC, H2O_631G_X2C_4IDX_DFCC_DIPEOM3) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_dipeom3",
           "h2o_6-31g_x2c_4index_dfccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_4IDX_DFCC_EAEOM2) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_eaeom2",
           "h2o_6-31g_x2c_4index_dfccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_4IDX_DFCC_IPEOM2) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_ipeom2",
           "h2o_6-31g_x2c_4index_dfccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_4IDX_DFEOMCCSD) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_4index_dfeomccsd",
           "h2o_6-31g_x2c_4index_dfeomccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_3IDX_DFCC_DIPEOM3) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_dipeom3",
           "h2o_6-31g_x2c_4index_dfccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_3IDX_DFCC_EAEOM2) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_eaeom2",
           "h2o_6-31g_x2c_4index_dfccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_3IDX_DFCC_IPEOM2) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_ipeom2",
           "h2o_6-31g_x2c_4index_dfccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, H2O_631G_X2C_3IDX_DFEOMCCSD) {
  CQCCTEST("coupledcluster/serial/df-eomcc/h2o_6-31g_x2c_3index_dfeomccsd",
           "h2o_6-31g_x2c_4index_dfeomccsd.bin.ref", "", true, true, true, false, false, false);
}

#ifdef _CQ_DO_PARTESTS

TEST( DF_EOMCC, PAR_H2O_631G_X2C_4IDX_DFCC_DIPEOM3) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_dipeom3",
           "h2o_6-31g_x2c_4index_dfccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_4IDX_DFCC_EAEOM2) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_eaeom2",
           "h2o_6-31g_x2c_4index_dfccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_4IDX_DFCC_IPEOM2) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_4index_dfccsd_ipeom2",
           "h2o_6-31g_x2c_4index_dfccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_4IDX_DFEOMCCSD) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_4index_dfeomccsd",
           "h2o_6-31g_x2c_4index_dfeomccsd.bin.ref", "", true, true, true, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_3IDX_DFCC_DIPEOM3) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_dipeom3",
           "h2o_6-31g_x2c_4index_dfccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_3IDX_DFCC_EAEOM2) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_eaeom2",
           "h2o_6-31g_x2c_4index_dfccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_3IDX_DFCC_IPEOM2) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_3index_dfccsd_ipeom2",
           "h2o_6-31g_x2c_4index_dfccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( DF_EOMCC, PAR_H2O_631G_X2C_3IDX_DFEOMCCSD) {
  CQCCTEST("coupledcluster/parallel/df-eomcc/h2o_6-31g_x2c_3index_dfeomccsd",
           "h2o_6-31g_x2c_4index_dfeomccsd.bin.ref", "", true, true, true, false, false, false);
}

#endif

