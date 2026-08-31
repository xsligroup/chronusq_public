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

#include "scf.hpp"

/* 
 * Note: Currently only SCF energies are tested. Properties need
 *       to be implemented.
 */

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, B3LYP Direct
TEST( DKS, UPu_184_plus_P_DC_B3LYP_Direct ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_B3LYP_Direct",
    "UPu_184+_P_DC_B3LYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, BLYP Direct
TEST( DKS, UPu_184_plus_P_DC_BLYP_Direct ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_BLYP_Direct",
    "UPu_184+_P_DC_BLYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, B3LYP Incore
TEST( DKS, UPu_184_plus_P_DC_B3LYP_Incore ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_B3LYP_Incore",
    "UPu_184+_P_DC_B3LYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, BLYP Incore
TEST( DKS, UPu_184_plus_P_DC_BLYP_Incore ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_BLYP_Incore",
    "UPu_184+_P_DC_BLYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};
////////
// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, B3LYP Direct
TEST( DKS, UPu_184_plus_P_DC_B3LYP_Direct_VLL ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_B3LYP_Direct_VLL",
    "UPu_184+_P_DC_B3LYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, BLYP Direct
TEST( DKS, UPu_184_plus_P_DC_BLYP_Direct_VLL ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_BLYP_Direct_VLL",
    "UPu_184+_P_DC_BLYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, B3LYP Incore
TEST( DKS, UPu_184_plus_P_DC_B3LYP_Incore_VLL ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_B3LYP_Incore_VLL",
    "UPu_184+_P_DC_B3LYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, BLYP Incore
TEST( DKS, UPu_184_plus_P_DC_BLYP_Incore_VLL ) {

  CQSCFTEST( "scf/serial/dks/UPu_184+_P_DC_BLYP_Incore_VLL",
    "UPu_184+_P_DC_BLYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#ifndef _CQ_GENERATE_TESTS

#endif

#ifdef _CQ_DO_PARTESTS

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, B3LYP Direct
TEST( DKS, PAR_UPu_184_plus_P_DC_B3LYP_Direct ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_B3LYP_Direct",
    "UPu_184+_P_DC_B3LYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, BLYP Direct
TEST( DKS, PAR_UPu_184_plus_P_DC_BLYP_Direct ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_BLYP_Direct",
    "UPu_184+_P_DC_BLYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, B3LYP Incore
TEST( DKS, PAR_UPu_184_plus_P_DC_B3LYP_Incore ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_B3LYP_Incore",
    "UPu_184+_P_DC_B3LYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham, BLYP Incore
TEST( DKS, PAR_UPu_184_plus_P_DC_BLYP_Incore ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_BLYP_Incore",
    "UPu_184+_P_DC_BLYP_Direct.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, B3LYP Direct
TEST( DKS, PAR_UPu_184_plus_P_DC_B3LYP_Direct_VLL ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_B3LYP_Direct_VLL",
    "UPu_184+_P_DC_B3LYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, BLYP Direct
TEST( DKS, PAR_UPu_184_plus_P_DC_BLYP_Direct_VLL ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_BLYP_Direct_VLL",
    "UPu_184+_P_DC_BLYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, B3LYP Incore
TEST( DKS, PAR_UPu_184_plus_P_DC_B3LYP_Incore_VLL ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_B3LYP_Incore_VLL",
    "UPu_184+_P_DC_B3LYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Kohn-Sham VXCLL Approx, BLYP Incore
TEST( DKS, PAR_UPu_184_plus_P_DC_BLYP_Incore_VLL ) {

  CQSCFTEST( "scf/parallel/dks/UPu_184+_P_DC_BLYP_Incore_VLL",
    "UPu_184+_P_DC_BLYP_Direct_VLL.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#endif


