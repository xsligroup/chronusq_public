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


// Water RHF/cc-pvdz(unc.) with RHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, Water_RHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_RHF_cc-pvdz_sto-3g_READDEN",
    "water_RHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_RHF.scr.bin" );

};

// Water RHF/cc-pvdz(unc.) with RHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, Water_RHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/serial/proj/water_RHF_cc-pvdz_sto-3g_READMO",
    "water_RHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_RHF.scr.bin" );

};

// Water ROHF/cc-pvdz(unc.) with ROHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, Water_ROHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_ROHF_cc-pvdz_sto-3g_READDEN",
    "water_ROHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_ROHF.scr.bin" );

};

// Water ROHF/cc-pvdz(unc.) with ROHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, Water_ROHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/serial/proj/water_ROHF_cc-pvdz_sto-3g_READMO",
    "water_ROHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_ROHF.scr.bin" );

};

// Water UHF/cc-pvdz(unc.) with UHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, Water_UHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_UHF_cc-pvdz_sto-3g_READDEN",
    "water_UHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_UHF.scr.bin" );

};

// Water UHF/cc-pvdz(unc.) with UHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, Water_UHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/serial/proj/water_UHF_cc-pvdz_sto-3g_READMO",
    "water_UHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_UHF.scr.bin" );

};

// Water X2CHF/cc-pvdz(unc.) with X2CHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, Water_X2CHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_X2CHF_cc-pvdz_sto-3g_READDEN",
    "water_X2CHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_X2CHF.scr.bin" );

};

// Water X2CHF/cc-pvdz(unc.) with X2CHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, Water_X2CHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/serial/proj/water_X2CHF_cc-pvdz_sto-3g_READMO",
    "water_X2CHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_X2CHF.scr.bin" );

};

// Water DCHF/cc-pvdz(unc.) with DCHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, Water_DCHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/serial/proj/water_DCHF_cc-pvdz_sto-3g_READDEN",
    "water_DCHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_DCHF.scr.bin" );

};

// Water DCHF/cc-pvdz(unc.) with DCHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, Water_DCHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/serial/proj/water_DCHF_cc-pvdz_sto-3g_READMO",
    "water_DCHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_DCHF.scr.bin" );

};

#ifdef _CQ_DO_PARTESTS

// Water RHF/cc-pvdz(unc.) with RHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, PAR_Water_RHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/parallel/proj/water_RHF_cc-pvdz_sto-3g_READDEN",
    "water_RHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_RHF.scr.bin" );

};

// Water RHF/cc-pvdz(unc.) with RHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, PAR_Water_RHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/parallel/proj/water_RHF_cc-pvdz_sto-3g_READMO",
    "water_RHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_RHF.scr.bin" );

};

// Water ROHF/cc-pvdz(unc.) with ROHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, PAR_Water_ROHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/parallel/proj/water_ROHF_cc-pvdz_sto-3g_READDEN",
    "water_ROHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_ROHF.scr.bin" );

};

// Water ROHF/cc-pvdz(unc.) with ROHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, PAR_Water_ROHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/parallel/proj/water_ROHF_cc-pvdz_sto-3g_READMO",
    "water_ROHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_ROHF.scr.bin" );

};

// Water UHF/cc-pvdz(unc.) with UHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, PAR_Water_UHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/parallel/proj/water_UHF_cc-pvdz_sto-3g_READDEN",
    "water_UHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_UHF.scr.bin" );

};

// Water UHF/cc-pvdz(unc.) with UHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, PAR_Water_UHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/parallel/proj/water_UHF_cc-pvdz_sto-3g_READMO",
    "water_UHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_UHF.scr.bin" );

};

// Water X2CHF/cc-pvdz(unc.) with X2CHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, PAR_Water_X2CHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/parallel/proj/water_X2CHF_cc-pvdz_sto-3g_READDEN",
    "water_X2CHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_X2CHF.scr.bin" );

};

// Water X2CHF/cc-pvdz(unc.) with X2CHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, PAR_Water_X2CHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/parallel/proj/water_X2CHF_cc-pvdz_sto-3g_READMO",
    "water_X2CHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_X2CHF.scr.bin" );

};

// Water DCHF/cc-pvdz(unc.) with DCHF/sto-3g(unc.) guess (READDEN)
TEST( PROJ_HF, PAR_Water_DCHF_ccpVDZ_sto3gGuess_READDEN ) {

  CQSCFTEST( "scf/parallel/proj/water_DCHF_cc-pvdz_sto-3g_READDEN",
    "water_DCHF_cc-pvdz_sto-3g_READDEN.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_DCHF.scr.bin" );

};

// Water DCHF/cc-pvdz(unc.) with DCHF/sto-3g(unc.) guess (READMO)
TEST( PROJ_HF, PAR_Water_DCHF_ccpVDZ_sto3gGuess_READMO ) {

  CQSCFTEST( "scf/parallel/proj/water_DCHF_cc-pvdz_sto-3g_READMO",
    "water_DCHF_cc-pvdz_sto-3g_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_sto-3g_DCHF.scr.bin" );

};

#endif




