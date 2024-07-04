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

// -------------------START RKS TESTS-------------------
// B3LYP / cc-pVTZ
TEST( GAUXC_KS, Water_ccpVTZ_B3LYP ) {

  CQSCFTEST( "scf/serial/gauxc_ks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );

}

// BLYP / cc-pVTZ
TEST( GAUXC_KS, Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/serial/gauxc_ks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// PBEXPBEC / cc-pVTZ
TEST( GAUXC_KS, Water_ccpVTZ_PBEXPBEC ) {

  CQSCFTEST( "scf/serial/gauxc_ks/water_cc-pVTZ_PBEXPBEC", "water_cc-pVTZ_PBEXPBEC.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// PBE0 / cc-pVTZ
TEST( GAUXC_KS, Water_ccpVTZ_PBE0 ) {

  CQSCFTEST( "scf/serial/gauxc_ks/water_cc-pVTZ_PBE0", "water_cc-pVTZ_PBE0.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}
// -------------------END RKS TESTS-------------------



// -------------------START UKS TESTS-------------------
// B3LYP / 6-311pG**
TEST( GAUXC_KS, Oxygen_6311pGss_B3LYP ) {

  CQSCFTEST( "scf/serial/gauxc_ks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// BLYP / 6-311pG**
TEST( GAUXC_KS, Oxygen_6311pGss_BLYP ) {

  CQSCFTEST( "scf/serial/gauxc_ks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}
// -------------------END UKS TESTS-------------------



// -------------------START X2CKS TESTS-------------------
// Water 6-311+G(d,p) B3LYP (Spherical) test
TEST( GAUXC_KS, Water_6311pGdp_x2c_b3lyp_sph ) {

  CQSCFTEST( "scf/serial/gauxc_ks/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// X2C O2 6-31G B3LYP test
TEST( GAUXC_KS, O2_triplet_x2cks_631G ) {

  CQSCFTEST( "scf/serial/gauxc_ks/o2_triplet_x2cks_631G", 
    "o2_triplet_x2cks_631G.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// Hg SAPPORO DZP DKH_2012 SP B3LYP
// (Turned off Quad and Oct check)
TEST( GAUXC_KS, Hg_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/serial/gauxc_ks/hg_sap_dz_dkh3_2012_sp_b3lyp", 
    "hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6,true,true,
               false,false,true,true);

};

// Cd SAPPORO DZP DKH_2012 SP B3LYP
TEST( GAUXC_KS, Cd_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/serial/gauxc_ks/cd_sap_dz_dkh3_2012_sp_b3lyp", 
    "cd_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};
// -------------------END X2CKS TESTS-------------------







#ifdef _CQ_DO_PARTESTS
// -------------------START RKS TESTS-------------------
// SMP B3LYP / cc-pVTZ
TEST( GAUXC_KS, PAR_Water_ccpVTZ_B3LYP ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/water_cc-pVTZ_B3LYP", "water_cc-pVTZ_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP BLYP / cc-pVTZ
TEST( GAUXC_KS, PAR_Water_ccpVTZ_BLYP ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/water_cc-pVTZ_BLYP", "water_cc-pVTZ_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP PBEXPBEC / cc-pVTZ
TEST( GAUXC_KS, PAR_Water_ccpVTZ_PBEXPBEC ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/water_cc-pVTZ_PBEXPBEC", "water_cc-pVTZ_PBEXPBEC.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// SMP PBE0 / cc-pVTZ
TEST( GAUXC_KS, PAR_Water_ccpVTZ_PBE0 ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/water_cc-pVTZ_PBE0", "water_cc-pVTZ_PBE0.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}
// -------------------END RKS TESTS-------------------



// -------------------START UKS TESTS-------------------
// B3LYP / 6-311pG**
TEST( GAUXC_KS, PAR_Oxygen_6311pGss_B3LYP ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/oxygen_6-311pG**_B3LYP", "oxygen_6-311pG**_B3LYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}

// BLYP / 6-311pG**
TEST( GAUXC_KS, PAR_Oxygen_6311pGss_BLYP ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/oxygen_6-311pG**_BLYP", "oxygen_6-311pG**_BLYP.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true  );

}
// -------------------END UKS TESTS-------------------



// -------------------START X2CKS TESTS-------------------
// Water 6-311+G(d,p) B3LYP (Spherical) test
TEST( GAUXC_KS, PAR_Water_6311pGdp_x2c_b3lyp_sph ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/water_6-311+Gdp_b3lyp_sph", 
    "water_6-311+Gdp_sph_x2c_b3lyp.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// X2C O2 6-31G B3LYP test
TEST( GAUXC_KS, PAR_O2_triplet_x2cks_631G ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/o2_triplet_x2cks_631G", 
    "o2_triplet_x2cks_631G.bin.ref", 1e-6, 
    true, true, true, true, true, true, false, "no", true );
 
};

// Hg SAPPORO DZP DKH_2012 SP B3LYP
// (Turned off Quad and Oct check)
TEST( GAUXC_KS, PAR_Hg_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/hg_sap_dz_dkh3_2012_sp_b3lyp", 
    "hg_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6,true,true,
               false,false,true,true);

};

// Cd SAPPORO DZP DKH_2012 SP B3LYP
TEST( GAUXC_KS, PAR_Cd_SAP_DZP_DKH3_2012_SP_B3LYP  ) {

  CQSCFTEST( "scf/parallel/gauxc_ks/cd_sap_dz_dkh3_2012_sp_b3lyp", 
    "cd_sap_dz_dkh3_2012_sp_b3lyp.bin.ref",1e-6 );
 
};
// -------------------END X2CKS TESTS-------------------
#endif



