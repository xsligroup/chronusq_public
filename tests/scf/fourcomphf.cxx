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

/* 
 * Note: Currently only SCF energies are tested. Properties need
 *       to be implemented.
 */

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, UPu_184_plus_P_NR_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_pointnuc",
    "UPu_184+_P_NR_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, UPu_184_plus_P_DC_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_pointnuc",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#ifndef _CQ_GENERATE_TESTS
// Two electron U-Pu 184+ test Dirac-Couloumb with readmo
TEST( FOURCHF, UPu_184_plus_P_DC_pointnuc_readmo ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_pointnuc_readmo",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true, true);

};
#endif

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, UPu_184_plus_P_DCG_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCG_pointnuc",
    "UPu_184+_P_DCG_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, UPu_184_plus_P_NR_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_finitenuc",
    "UPu_184+_P_NR_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, UPu_184_plus_P_DC_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DC_finitenuc",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, UPu_184_plus_P_DCG_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCG_finitenuc",
    "UPu_184+_P_DCG_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS
TEST( FOURCHF, UPu_184_plus_P_DCGGS_finitenuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_DCGGS_finitenuc",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Distorted CuH3 to test C1 symmetry Dirac-Coulomb-Gaunt
TEST( FOURCHF, CuH3_321g_DCG ) {

  CQSCFTEST( "scf/serial/fourcomp/CuH3_321g_DCG",
    "CuH3_321g_DCG.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Ag Neutral doublet Dirac-Coulomb-Gaunt
TEST( FOURCHF, Ag_sapporoDZ_DCG ) {

  CQSCFTEST( "scf/serial/fourcomp/Ag_sapporoDZ_DCG",
    "Ag_sapporoDZ_DCG.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Water DC-HF/cc-pVDZ with RHF guess (READDEN)
TEST( FOURCHF, Water_ccpVDZ_DCHF_RHFGuess_READDEN ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_RHFGuess_READDEN.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_RHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with RHF guess (READMO)
TEST( FOURCHF, Water_ccpVDZ_DCHF_RHFGuess_READMO ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_RHFGuess_READMO.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_RHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with ROHF guess (READDEN)
TEST( FOURCHF, Water_ccpVDZ_DCHF_ROHFGuess_READDEN ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_ROHFGuess_READDEN.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_ROHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with ROHF guess (READMO)
TEST( FOURCHF, Water_ccpVDZ_DCHF_ROHFGuess_READMO ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_ROHFGuess_READMO.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_ROHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with UHF guess (READDEN)
TEST( FOURCHF, Water_ccpVDZ_DCHF_UHFGuess_READDEN ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_UHFGuess_READDEN.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_UHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with UHF guess (READMO)
TEST( FOURCHF, Water_ccpVDZ_DCHF_UHFGuess_READMO ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_UHFGuess_READMO.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_UHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with GHF guess (READMO)
TEST( FOURCHF, Water_ccpVDZ_DCHF_GHFGuess_READMO ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_GHFGuess_READMO.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_GHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with X2C guess (READDEN)
TEST( FOURCHF, Water_ccpVDZ_DCHF_X2CGuess_READDEN ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_X2CGuess_READDEN.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_X2C.scr.bin" );

};

// Water DC-HF/cc-pVDZ with X2C guess (READMO)
TEST( FOURCHF, Water_ccpVDZ_DCHF_X2CGuess_READMO ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_X2CGuess_READMO.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_X2C.scr.bin" );

};

// Water DC-HF/cc-pVDZ with DC-HF guess (READDEN from scratch)
TEST( FOURCHF, Water_ccpVDZ_DCHF_DCHFGuess_READDEN_scratch ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_DCHFGuess_READDEN_scratch.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_DCHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with DC-HF guess (READMO from scratch)
TEST( FOURCHF, Water_ccpVDZ_DCHF_DCHFGuess_READMO_scratch ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_DCHFGuess_READMO_scratch.bin.ref",1e-6,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_DCHF.scr.bin" );

};

// Water DC-HF/cc-pVDZ with DC-HF guess (READDEN from restart)
TEST( FOURCHF, Water_ccpVDZ_DCHF_DCHFGuess_READDEN_restart ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READDEN",
    "water_cc-pVDZ_DC-HF_DCHFGuess_READDEN_restart.bin.ref",1e-6,
    false, false, false, false, false, true, true );

};

// Water DC-HF/cc-pVDZ with DC-HF guess (READMO from restart)
TEST( FOURCHF, Water_ccpVDZ_DCHF_DCHFGuess_READMO_restart ) {

  CQSCFTEST( "scf/serial/fourcomp/water_cc-pVDZ_DC-HF_READMO",
    "water_cc-pVDZ_DC-HF_DCHFGuess_READMO_restart.bin.ref",1e-6,
    false, false, false, false, false, true, true );

};

#ifdef _CQ_DO_PARTESTS

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, PAR_UPu_184_plus_P_NR_pointnuc ) {

  CQSCFTEST( "scf/serial/fourcomp/UPu_184+_P_NR_pointnuc",
    "UPu_184+_P_NR_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_pointnuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_pointnuc",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, PAR_UPu_184_plus_P_DCG_pointnuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCG_pointnuc",
    "UPu_184+_P_DCG_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#ifndef _CQ_GENERATE_TESTS
// Two electron U-Pu 184+ test Dirac-Couloumb with readmo
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_pointnuc_readmo ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_pointnuc_readmo",
    "UPu_184+_P_DC_pointnuc.bin.ref",1e-8,
    false, false, false, false, false, true, true);

};
#endif

// Two electron U-Pu 184+ test Dirac-HF no 2ERI relativitic, LLLL only
TEST( FOURCHF, PAR_UPu_184_plus_P_NR_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_NR_finitenuc",
    "UPu_184+_P_NR_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb
TEST( FOURCHF, PAR_UPu_184_plus_P_DC_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DC_finitenuc",
    "UPu_184+_P_DC_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt
TEST( FOURCHF, PAR_UPu_184_plus_P_DCG_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCG_finitenuc",
    "UPu_184+_P_DCG_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

// Two electron U-Pu 184+ test Dirac-Couloumb-Gaunt-Gauge-SSSS
TEST( FOURCHF, PAR_UPu_184_plus_P_DCGGS_finitenuc ) {

  CQSCFTEST( "scf/parallel/fourcomp/UPu_184+_P_DCGGS_finitenuc",
    "UPu_184+_P_DCGGS_finitenuc.bin.ref",1e-8,
    false, false, false, false, false, true);

};

#endif


