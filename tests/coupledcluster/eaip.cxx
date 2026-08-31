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

TEST( EAIP_EOMCC, H2O_631G_GHF_CCSD_IPEOM2) {
  CQCCTEST("coupledcluster/serial/eaip/h2o_6-31g_ghf_ccsd_ipeom2",
           "h2o_6-31g_ghf_ccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, H2O_631G_GHF_CCSD_EAEOM2) {
  CQCCTEST("coupledcluster/serial/eaip/h2o_6-31g_ghf_ccsd_eaeom2",
           "h2o_6-31g_ghf_ccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, H2O_631G_GHF_CCSD_DIPEOM3) {
  CQCCTEST("coupledcluster/serial/eaip/h2o_6-31g_ghf_ccsd_dipeom3",
           "h2o_6-31g_ghf_ccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, H4_631G_GHF_CCSD_IPEOM3) {
  CQCCTEST("coupledcluster/serial/eaip/h4_6-31g_ghf_ccsd_ipeom3",
           "h4_6-31g_ghf_ccsd_ipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, H4_631G_GHF_CCSDT_IPEOM3) {
  CQCCTEST("coupledcluster/serial/eaip/h4_6-31g_ghf_ccsdt_ipeom3",
           "h4_6-31g_ghf_ccsdt_ipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, H4_631G_GHF_CCSDT_DIPEOM4) {
  CQCCTEST("coupledcluster/serial/eaip/h4_6-31g_ghf_ccsdt_dipeom4",
           "h4_6-31g_ghf_ccsdt_dipeom4.bin.ref", "", true, true, false, false, false, false);
}

#ifdef _CQ_DO_PARTESTS

TEST( EAIP_EOMCC, PAR_H2O_631G_GHF_CCSD_IPEOM2) {
  CQCCTEST("coupledcluster/parallel/eaip/h2o_6-31g_ghf_ccsd_ipeom2",
           "h2o_6-31g_ghf_ccsd_ipeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, PAR_H2O_631G_GHF_CCSD_EAEOM2) {
  CQCCTEST("coupledcluster/parallel/eaip/h2o_6-31g_ghf_ccsd_eaeom2",
           "h2o_6-31g_ghf_ccsd_eaeom2.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, PAR_H2O_631G_GHF_CCSD_DIPEOM3) {
  CQCCTEST("coupledcluster/parallel/eaip/h2o_6-31g_ghf_ccsd_dipeom3",
           "h2o_6-31g_ghf_ccsd_dipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, PAR_H4_631G_GHF_CCSD_IPEOM3) {
  CQCCTEST("coupledcluster/parallel/eaip/h4_6-31g_ghf_ccsd_ipeom3",
           "h4_6-31g_ghf_ccsd_ipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, PAR_H4_631G_GHF_CCSDT_IPEOM3) {
  CQCCTEST("coupledcluster/parallel/eaip/h4_6-31g_ghf_ccsdt_ipeom3",
           "h4_6-31g_ghf_ccsdt_ipeom3.bin.ref", "", true, true, false, false, false, false);
}

TEST( EAIP_EOMCC, PAR_H4_631G_GHF_CCSDT_DIPEOM4) {
  CQCCTEST("coupledcluster/parallel/eaip/h4_6-31g_ghf_ccsdt_dipeom4",
           "h4_6-31g_ghf_ccsdt_dipeom4.bin.ref", "", true, true, false, false, false, false);
}

#endif

