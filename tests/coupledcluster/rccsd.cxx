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


TEST( RCCSD, H2O_STO3G_RHF_CCSD) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_sto3g_rhf_ccsd",
    "h2o_sto3g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_STO3G_RHF_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_sto3g_rhf_ccsd_denomshift",
    "h2o_sto3g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_STO3G_RHF_CCSD_RESTART) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_sto3g_rhf_ccsd_t_restart",
           "h2o_sto3g_rhf_ccsd.bin.ref", "h2o_sto3g_rhf_ccsd.bin.ref", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd",
           "h2o_631g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD_FROZENCORE) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd_frozencore",
           "h2o_631g_rhf_ccsd_frozencore.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD_FROZENVIRTUAL) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd_frozenvirt",
           "h2o_631g_rhf_ccsd_frozenvirt.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD_FROZENCOREVIRTUAL) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd_frozencorevirt",
           "h2o_631g_rhf_ccsd_frozencorevirt.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD_REBUILDFOCK) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd_rebuildfock",
           "h2o_631g_rhf_ccsd_frozencore.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, H2O_631G_RHF_CCSD_RESTART_DIFFERENT_GEOMETRY) {
  CQCCTEST("coupledcluster/serial/rccsd/h2o_631g_rhf_ccsd_restart_different_geometry",
           "h2o_631g_rhf_ccsd_restart_different_geometry.bin.ref", "h2o_631g_rhf_ccsd.bin.ref", true, false, false, false, false, true);
}

#ifdef _CQ_DO_PARTESTS

TEST( RCCSD, PAR_H2O_STO3G_RHF_CCSD) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_sto3g_rhf_ccsd",
           "h2o_sto3g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_STO3G_RHF_CCSD_DENOMSHIFT) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_sto3g_rhf_ccsd_denomshift",
           "h2o_sto3g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_STO3G_RHF_CCSD_RESTART) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_sto3g_rhf_ccsd_t_restart",
           "h2o_sto3g_rhf_ccsd.bin.ref", "h2o_sto3g_rhf_ccsd.bin.ref", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd",
           "h2o_631g_rhf_ccsd.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD_FROZENCORE) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd_frozencore",
           "h2o_631g_rhf_ccsd_frozencore.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD_FROZENVIRTUAL) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd_frozenvirt",
           "h2o_631g_rhf_ccsd_frozenvirt.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD_FROZENCOREVIRTUAL) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd_frozencorevirt",
           "h2o_631g_rhf_ccsd_frozencorevirt.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD_REBUILDFOCK) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd_rebuildfock",
           "h2o_631g_rhf_ccsd_frozencore.bin.ref", "", false, false, false, false, false, true);
}

TEST( RCCSD, PAR_H2O_631G_RHF_CCSD_RESTART_DIFFERENT_GEOMETRY) {
  CQCCTEST("coupledcluster/parallel/rccsd/h2o_631g_rhf_ccsd_restart_different_geometry",
           "h2o_631g_rhf_ccsd_restart_different_geometry.bin.ref", "h2o_631g_rhf_ccsd.bin.ref", true, false, false, false, false, true);
}

#endif

