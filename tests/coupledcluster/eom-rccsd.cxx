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


TEST( EOM_RCCSD, H2O_631G_RHF_EOM_CCSD) {
  CQCCTEST("coupledcluster/serial/eom-rccsd/h2o_631g_rhf_eom_ccsd",
           "h2o_631g_rhf_eom_ccsd.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, H2O_631G_RHF_EOM_CCSD_FROZENCORE) {
  CQCCTEST("coupledcluster/serial/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozencore",
           "h2o_631g_rhf_eom_ccsd_frozencore.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, H2O_631G_RHF_EOM_CCSD_FROZENVIRT) {
  CQCCTEST("coupledcluster/serial/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozenvirt",
           "h2o_631g_rhf_eom_ccsd_frozenvirt.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, H2O_631G_RHF_EOM_CCSD_FROZENCOREVIRT) {
  CQCCTEST("coupledcluster/serial/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozencorevirt",
           "h2o_631g_rhf_eom_ccsd_frozencorevirt.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, H2O_631G_RHF_EOM_CCSD_DIPOLES) {
  CQCCTEST("coupledcluster/serial/eom-rccsd/h2o_631g_rhf_eom_ccsd_dipoles",
           "h2o_631g_rhf_eom_ccsd_dipoles.bin.ref", "", true, true, true, true, false, true);
}

#ifdef _CQ_DO_PARTESTS

TEST( EOM_RCCSD, PAR_H2O_631G_RHF_EOM_CCSD) {
  CQCCTEST("coupledcluster/parallel/eom-rccsd/h2o_631g_rhf_eom_ccsd",
           "h2o_631g_rhf_eom_ccsd.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, PAR_H2O_631G_RHF_EOM_CCSD_FROZENCORE) {
  CQCCTEST("coupledcluster/parallel/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozencore",
           "h2o_631g_rhf_eom_ccsd_frozencore.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, PAR_H2O_631G_RHF_EOM_CCSD_FROZENVIRT) {
  CQCCTEST("coupledcluster/parallel/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozenvirt",
           "h2o_631g_rhf_eom_ccsd_frozenvirt.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, PAR_H2O_631G_RHF_EOM_CCSD_FROZENCOREVIRT) {
  CQCCTEST("coupledcluster/parallel/eom-rccsd/h2o_631g_rhf_eom_ccsd_frozencorevirt",
           "h2o_631g_rhf_eom_ccsd_frozencorevirt.bin.ref", "", true, true, true, false, false, true);
}

TEST( EOM_RCCSD, PAR_H2O_631G_RHF_EOM_CCSD_DIPOLES) {
  CQCCTEST("coupledcluster/parallel/eom-rccsd/h2o_631g_rhf_eom_ccsd_dipoles",
           "h2o_631g_rhf_eom_ccsd_dipoles.bin.ref", "", true, true, true, true, false, true);
}


#endif

