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

// X2C Water 6-311+G(d,p) B3LYP (Spherical) test
TEST( X2CKS, Water_6311pGdp_x2c_b3lyp_sph ) {

  CQSCFTEST( "scf/serial/x2c/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// Water 6-311+G(d,p) B3LYP (Cartesian) test
TEST( X2CKS, Water_6311pGdp_x2c_b3lyp_cart ) {

  CQSCFTEST( "scf/serial/x2c/water_6-311+Gdp_b3lyp_cart", 
    "water_6-311+Gdp_cart_x2c_b3lyp.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// X2C O2 6-31G B3LYP test
TEST( X2CKS, O2_triplet_x2cks_631G ) {

  CQSCFTEST( "scf/serial/x2c/o2_triplet_x2cks_631G", 
    "o2_triplet_x2cks_631G.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// Hg SAPPORO DZP DKH_2012 SP SLATER
TEST( X2CKS, Hg_SAP_DZP_DKH3_2012_SP_SLATER  ) {

  CQSCFTEST( "scf/serial/x2c/hg_sap_dz_dkh3_2012_sp_slater", 
    "hg_sap_dz_dkh3_2012_sp_slater.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

// Cd SAPPORO DZP DKH_2012 SP SLATER
TEST( X2CKS, Cd_SAP_DZP_DKH3_2012_SP_SLATER  ) {

  CQSCFTEST( "scf/serial/x2c/cd_sap_dz_dkh3_2012_sp_slater", 
    "cd_sap_dz_dkh3_2012_sp_slater.bin.ref",1e-6 );
 
};



// Hg SAPPORO DZP DKH_2012 SP B3LYP
// (Turned off Quad and Oct check)
TEST( X2CKS, Hg_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/serial/x2c/hg_sap_dz_dkh3_2012_sp_b3lyp", 
    "hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6,true,true,
               false,false,true,true);

};

// Cd SAPPORO DZP DKH_2012 SP B3LYP
TEST( X2CKS, Cd_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/serial/x2c/cd_sap_dz_dkh3_2012_sp_b3lyp", 
    "cd_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};

#ifdef _CQ_DO_PARTESTS

// SMP X2C Water 6-311+G(d,p) B3LYP (Spherical) test
TEST( X2CKS, PAR_Water_6311pGdp_x2c_b3lyp_sph ) {

  CQSCFTEST( "scf/parallel/x2c/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};


/*
// SMP Hg SAPPORO DZP DKH_2012 SP SLATER
TEST( X2CKS, PAR_Hg_SAP_DZP_DKH3_2012_SP_SLATER  ) {

  CQSCFTEST( scf/parallel/x2c/hg_sap_dz_dkh3_2012_sp_slater, 
    hg_sap_dz_dkh3_2012_sp_slater.bin.ref );
 
};

// SMP Hg SAPPORO DZP DKH_2012 SP B3LYP
TEST( X2CKS, PAR_Hg_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( scf/parallel/x2c/hg_sap_dz_dkh3_2012_sp_b3lyp, 
    hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref );
 
};
*/

#endif



