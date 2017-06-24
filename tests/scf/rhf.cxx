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


// Water 6-31G(d) test
TEST( RHF, Water_631Gd ) {

  CQSCFTEST( "scf/serial/rhf/water_6-31Gd", "water_6-31Gd.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

// NeHe sto-3g (specifying atoms by atomic number) test
TEST( RHF, NeHe_RHF_STO3G_ATNUM ) {

  CQSCFTEST( "scf/serial/rhf/NeHe_rhf_sto-3g_atNum", "NeHe_rhf_sto-3g_atNum.bin.ref" );
 
};

// NeHe sto-3g (specifying atoms by default isotope) test
TEST( RHF, NeHe_RHF_STO3G_ISO ) {

  CQSCFTEST( "scf/serial/rhf/NeHe_rhf_sto-3g_iso", "NeHe_rhf_sto-3g_iso.bin.ref" );
 
};

// HeKr sto-3g test for fchk parsing
TEST( RHF, HeKr_sto3G ) {

  CQSCFTEST( "scf/serial/rhf/HeKr_sto-3G", "HeKr_sto-3G.bin.ref",
      1e-8, true, true, true, true, true, true,
      false, "HeKr_sto-3G.fchk" );

};

// H2S cc-pVDZ incore ERI test
TEST( RHF, H2S_ccpvdz ) {

  CQSCFTEST( "scf/serial/rhf/H2S_cc-pvdz", "H2S_cc-pvdz.bin.ref" );

};

// H2S cc-pVDZ incore Libcint ERI test
TEST( RHF, H2S_ccpvdz_libcint ) {

  CQSCFTEST( "scf/serial/rhf/H2S_cc-pvdz_libcint", "H2S_cc-pvdz.bin.ref" );

};

// Water RHF/cc-pVDZ with RHF guess (READDEN from scratch)
TEST( RHF, Water_ccpVDZ_RHF_RHFGuess_READDEN_scratch ) {

  CQSCFTEST( "scf/serial/rhf/water_cc-pVDZ_RHF_READDEN",
    "water_cc-pVDZ_RHF_RHFGuess_READDEN_scratch.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_RHF.scr.bin" );

};

// Water RHF/cc-pVDZ with RHF guess (READMO from scratch)
TEST( RHF, Water_ccpVDZ_RHF_RHFGuess_READMO_scratch ) {

  CQSCFTEST( "scf/serial/rhf/water_cc-pVDZ_RHF_READMO",
    "water_cc-pVDZ_RHF_RHFGuess_READMO_scratch.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "water_cc-pVDZ_RHF.scr.bin" );

};

// Water RHF/cc-pVDZ with RHF guess (READDEN from restart)
TEST( RHF, Water_ccpVDZ_RHF_RHFGuess_READDEN_restart ) {

  CQSCFTEST( "scf/serial/rhf/water_cc-pVDZ_RHF_READDEN",
    "water_cc-pVDZ_RHF_RHFGuess_READDEN_restart.bin.ref",1e-8,
    false, false, false, false, false, true, true );

};

// Water RHF/cc-pVDZ with RHF guess (READMO from restart)
TEST( RHF, Water_ccpVDZ_RHF_RHFGuess_READMO_restart ) {

  CQSCFTEST( "scf/serial/rhf/water_cc-pVDZ_RHF_READMO",
    "water_cc-pVDZ_RHF_RHFGuess_READMO_restart.bin.ref",1e-8,
    false, false, false, false, false, true, true );

};

#ifdef _CQ_DO_PARTESTS

// SMP Water 6-31G(d) test
TEST( RHF, PAR_Water_631Gd ) {

  CQSCFTEST( "scf/parallel/rhf/water_6-31Gd", "water_6-31Gd.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

// H2S cc-pVDZ incore ERI test
TEST( RHF, PAR_H2S_ccpvdz ) {

  CQSCFTEST( "scf/parallel/rhf/H2S_cc-pvdz", "H2S_cc-pvdz.bin.ref" );

};

// H2S cc-pVDZ incore Libcint ERI test
TEST( RHF, PAR_H2S_ccpvdz_libcint ) {

  CQSCFTEST( "scf/parallel/rhf/H2S_cc-pvdz_libcint", "H2S_cc-pvdz.bin.ref" );

};

#endif



