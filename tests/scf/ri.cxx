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



// Water RHF/6-31G(d)/cc-pVDZ-rifit test
TEST( RI_RHF, water_631Gd_ccpvdzrifit ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cc-pvdz-rifit",
             "water_6-31Gd_cc-pvdz-rifit.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) default Cholesky test
TEST( RI_RHF, water_631Gd_cd ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) traditional Cholesky test
TEST( RI_RHF, water_631Gd_cd_traditional ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd_traditional",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) dyanmic Cholesky test
TEST( RI_RHF, water_631Gd_cd_dynamicall ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd_dynamicall",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) original span-factor Cholesky test
TEST( RI_RHF, water_631Gd_cd_spanfactor ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd_spanfactor",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) dynamic-ERI span-factor Cholesky test
TEST( RI_RHF, water_631Gd_cd_dynamiceri ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd_dynamiceri",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) reuse-Mpq span-factor Cholesky test
TEST( RI_RHF, water_631Gd_cd_spanfactorreuse ) {

  CQSCFTEST( "scf/serial/ri_rhf/water_6-31Gd_cd_spanfactorreuse",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ traditional Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_traditional ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_traditional",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ traditional libcint Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_traditional_libcint ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_traditional_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-ERI span-factor Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_dynamiceri ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_dynamiceri",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-ERI span-factor libcint Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_dynamiceri_libcint ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_dynamiceri_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-all Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_dynamicall ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_dynamicall",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-all libcint Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_dynamicall_libcint ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_dynamicall_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_spanfactor ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_spanfactor",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor libcint Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_spanfactor_libcint ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_spanfactor_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor-reuse Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_spanfactorreuse ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_spanfactorreuse",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor-reuse libcint Cholesky test
TEST( RI_RHF, H2S_ccpVDZ_cd_spanfactorreuse_libcint ) {

  CQSCFTEST( "scf/serial/ri_rhf/H2S_cc-pvdz_cd_spanfactorreuse_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RB3LYP/6-31G(d)/cc-pVDZ-rifit test
TEST( RI_RKS, water_631Gd_ccpvdzrifit_B3LYP ) {

  CQSCFTEST( "scf/serial/ri_rks/water_6-31Gd_cc-pvdz-rifit_b3lyp",
             "water_6-31Gd_cc-pvdz-rifit_b3lyp.bin.ref" );

};

// O2 ROHF/cc-pVTZ/cc-pVTZ-jkfit test
TEST( RI_ROHF, Oxygen_ccpvtzjkfit ) {

  CQSCFTEST( "scf/serial/ri_rohf/oxygen_cc-pvtz_cc-pvtz-jkfit",
             "oxygen_cc-pvtz_cc-pvtz-jkfit.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true);

};

// O2 UHF/def2-tzvpd/def2-tzvpd-rifit test
TEST( RI_UHF, Oxygen_def2tzvp_rifit ) {

  CQSCFTEST( "scf/serial/ri_uhf/oxygen_def2-tzvp_def2-tzvp-rifit",
             "oxygen_def2-tzvp_def2-tzvp-rifit.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true);

};

// O2 UPBE0/6-31++g(d)/aug-cc-pvdz-rifit test
TEST( RI_UKS, Oxygen_631ppGd_augccpvdzrifit_PBE0 ) {

  CQSCFTEST( "scf/serial/ri_uks/oxygen_631++Gd_aug-cc-pvdz-rifit_pbe0",
             "oxygen_631++Gd_aug-cc-pvdz-rifit_pbe0.bin.ref" );

};

// Cd Scalar X2C UHF/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CHF, cd_sap_dz_dkh3_2012_sp_x2c_jfit ) {

  CQSCFTEST( "scf/serial/ri_x2c/cd_sap-dz-dkh3-2012-sp_x2c-jfit",
             "cd_sap-dz-dkh3-2012-sp_x2c-jfit.bin.ref",
             1e-6 );

};

// Hg X2C UHF/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CHF, hg_sap_dz_dkh3_2012_sp_x2c_jfit ) {

  CQSCFTEST( "scf/serial/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c-jfit",
             "hg_sap-dz-dkh3-2012-sp_x2c-jfit.bin.ref", 1e-6, 
            true, true, true, true, true, true, false, "no", true );

};

// Hg X2C HF/sapporo_dz_dkh3_2012_sp default Cholesky test
TEST( RI_X2CHF, hg_sap_dz_dkh3_2012_sp_x2c_cd ) {

  CQSCFTEST( "scf/serial/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c_cd",
             "hg_sap-dz-dkh3-2012-sp_x2c_cd.bin.ref",
             1e-6 );

};

// Hg X2C B3LYP/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CKS, hg_sap_dz_dkh3_2012_sp_x2c_jfit_b3lyp ) {

  CQSCFTEST( "scf/serial/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c-jfit_b3lyp",
             "hg_sap-dz-dkh3-2012-sp_x2c-jfit_b3lyp.bin.ref", 1e-6, 
            true, true, true, true, true, true, false, "no", true );

};

#ifdef _CQ_DO_PARTESTS
// Water RHF/6-31G(d)/cc-pVDZ-rifit test
TEST( RI_RHF, PAR_water_631Gd_ccpvdzrifit ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cc-pvdz-rifit",
             "water_6-31Gd_cc-pvdz-rifit.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) default Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true);

};

// Water RHF/6-31G(d) traditional Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd_traditional ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd_traditional",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) dynamic Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd_dynamicall ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd_dynamicall",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) original span-factor Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd_spanfactor ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd_spanfactor",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) dynamic-ERI span-factor Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd_dynamiceri ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd_dynamiceri",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true );

};

// Water RHF/6-31G(d) reuse-Mpq span-factor Cholesky test
TEST( RI_RHF, PAR_water_631Gd_cd_spanfactorreuse ) {

  CQSCFTEST( "scf/parallel/ri_rhf/water_6-31Gd_cd_spanfactorreuse",
             "water_6-31Gd_cd.bin.ref", 1e-6, 
             true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ traditional Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_traditional ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_traditional",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ traditional libcint Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_traditional_libcint ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_traditional_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-ERI span-factor Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_dynamiceri ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_dynamiceri",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-ERI span-factor libcint Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_dynamiceri_libcint ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_dynamiceri_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-all Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_dynamicall ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_dynamicall",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ dynamic-all libcint Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_dynamicall_libcint ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_dynamicall_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_spanfactor ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_spanfactor",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor libcint Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_spanfactor_libcint ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_spanfactor_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor-reuse Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_spanfactorreuse ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_spanfactorreuse",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// H2S RHF/cc-pVDZ span-factor-reuse libcint Cholesky test
TEST( RI_RHF, PAR_H2S_ccpVDZ_cd_spanfactorreuse_libcint ) {

  CQSCFTEST( "scf/parallel/ri_rhf/H2S_cc-pvdz_cd_spanfactorreuse_libcint",
             "H2S_cc-pvdz_cd.bin.ref", 1e-6, 
              true, true, true, true, true, true, false, "no", true );

};

// Water RB3LYP/6-31G(d)/cc-pVDZ-rifit test
TEST( RI_RKS, PAR_water_631Gd_ccpvdzrifit_B3LYP ) {

  CQSCFTEST( "scf/parallel/ri_rks/water_6-31Gd_cc-pvdz-rifit_b3lyp",
             "water_6-31Gd_cc-pvdz-rifit_b3lyp.bin.ref" );

};

// O2 ROHF/cc-pVTZ/cc-pVTZ-jkfit test
TEST( RI_ROHF, PAR_Oxygen_ccpvtzjkfit ) {

  CQSCFTEST( "scf/parallel/ri_rohf/oxygen_cc-pvtz_cc-pvtz-jkfit",
             "oxygen_cc-pvtz_cc-pvtz-jkfit.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );

};

// O2 UHF/def2-tzvpd/def2-tzvpd-rifit test
TEST( RI_UHF, PAR_Oxygen_def2tzvp_rifit ) {

  CQSCFTEST( "scf/parallel/ri_uhf/oxygen_def2-tzvp_def2-tzvp-rifit",
             "oxygen_def2-tzvp_def2-tzvp-rifit.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );

};

// O2 UPBE0/6-31++g(d)/aug-cc-pvdz-rifit test
TEST( RI_UKS, PAR_Oxygen_631ppGd_augccpvdzrifit_PBE0 ) {

  CQSCFTEST( "scf/parallel/ri_uks/oxygen_631++Gd_aug-cc-pvdz-rifit_pbe0",
             "oxygen_631++Gd_aug-cc-pvdz-rifit_pbe0.bin.ref" );

};

// Cd Scalar X2C UHF/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CHF, PAR_cd_sap_dz_dkh3_2012_sp_x2c_jfit ) {

  CQSCFTEST( "scf/parallel/ri_x2c/cd_sap-dz-dkh3-2012-sp_x2c-jfit",
             "cd_sap-dz-dkh3-2012-sp_x2c-jfit.bin.ref",
             1e-6 );

};

// Hg X2C UHF/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CHF, PAR_hg_sap_dz_dkh3_2012_sp_x2c_jfit ) {

  CQSCFTEST( "scf/parallel/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c-jfit",
             "hg_sap-dz-dkh3-2012-sp_x2c-jfit.bin.ref", 1e-6, 
            true, true, true, true, true, true, false, "no", true  );

};

// Hg X2C HF/sapporo_dz_dkh3_2012_sp default Cholesky test
TEST( RI_X2CHF, PAR_hg_sap_dz_dkh3_2012_sp_x2c_cd ) {

  CQSCFTEST( "scf/parallel/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c_cd",
             "hg_sap-dz-dkh3-2012-sp_x2c_cd.bin.ref",
             1e-6 );

};

// Hg X2C B3LYP/sapporo_dz_dkh3_2012_sp/x2c-jfit test
TEST( RI_X2CKS, PAR_hg_sap_dz_dkh3_2012_sp_x2c_jfit_b3lyp ) {

  CQSCFTEST( "scf/parallel/ri_x2c/hg_sap-dz-dkh3-2012-sp_x2c-jfit_b3lyp",
             "hg_sap-dz-dkh3-2012-sp_x2c-jfit_b3lyp.bin.ref", 1e-6, 
            true, true, true, true, true, true, false, "no", true  );

};
#endif




/**********************************/
// RI_NEO_RHF TESTS, SERIAL
/**********************************/

// EAUX TEST:
// do CD for (ee|ee), do incore 4-index for (pp|pp), use elec aux for (ee|pp)
TEST( RI_NEO_RHF, coh2_ccpvdz_pb4d_cde_4I_eaux ) {

  CQNEOSCFTEST( "scf/serial/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_4I_eaux",
             "coh2_ccpvdz_pb4d_cde_4I_eaux.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};


// PAUX TEST:
// do incore 4-index for (ee|ee), do CD for (pp|pp), use prot aux for (ee|pp)
TEST( RI_NEO_RHF, coh2_ccpvdz_pb4d_4I_cdp_paux ) {

  CQNEOSCFTEST( "scf/serial/ri_neo_rhf/coh2_ccpvdz_pb4d_4I_cdp_paux",
             "coh2_ccpvdz_pb4d_4I_cdp_paux.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};


// EAUX + PAUX TEST, CONNECTOR:
// do CD for (ee|ee), do CD for (pp|pp), use elec + prot aux CONNECTOR for (ee|pp)
TEST( RI_NEO_RHF, coh2_ccpvdz_pb4d_cde_cdp_connector ) {

  CQNEOSCFTEST( "scf/serial/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_cdp_connector",
             "coh2_ccpvdz_pb4d_cde_cdp_connector.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};

// EAUX + PAUX TEST, COMBINEAUXBASIS:
// do CD for (ee|ee), do CD for (pp|pp), use elec + prot aux COMBINEAUXBASIS for (ee|pp)
TEST( RI_NEO_RHF, coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis ) {

  CQNEOSCFTEST( "scf/serial/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis",
             "coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};

#ifdef _CQ_DO_PARTESTS
/**********************************/
// RI_NEO_RHF TESTS, PARALLEL
/**********************************/

// EAUX TEST:
// do CD for (ee|ee), do incore 4-index for (pp|pp), use elec aux for (ee|pp)
TEST( RI_NEO_RHF, par_coh2_ccpvdz_pb4d_cde_4I_eaux ) {

  CQNEOSCFTEST( "scf/parallel/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_4I_eaux",
             "coh2_ccpvdz_pb4d_cde_4I_eaux.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};


// PAUX TEST:
// do incore 4-index for (ee|ee), do CD for (pp|pp), use prot aux for (ee|pp)
TEST( RI_NEO_RHF, par_coh2_ccpvdz_pb4d_4I_cdp_paux ) {

  CQNEOSCFTEST( "scf/parallel/ri_neo_rhf/coh2_ccpvdz_pb4d_4I_cdp_paux",
             "coh2_ccpvdz_pb4d_4I_cdp_paux.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};


// EAUX + PAUX TEST, CONNECTOR:
// do CD for (ee|ee), do CD for (pp|pp), use elec + prot aux CONNECTOR for (ee|pp)
TEST( RI_NEO_RHF, par_coh2_ccpvdz_pb4d_cde_cdp_connector ) {

  CQNEOSCFTEST( "scf/parallel/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_cdp_connector",
             "coh2_ccpvdz_pb4d_cde_cdp_connector.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};

// EAUX + PAUX TEST, COMBINEAUXBASIS:
// do CD for (ee|ee), do CD for (pp|pp), use elec + prot aux COMBINEAUXBASIS for (ee|pp)
TEST( RI_NEO_RHF, par_coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis ) {

  CQNEOSCFTEST( "scf/parallel/ri_neo_rhf/coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis",
             "coh2_ccpvdz_pb4d_cde_cdp_combineauxbasis.bin.ref", 1e-6, 
              true, true, true, true, false, "no", true, false ) ;

};
#endif




