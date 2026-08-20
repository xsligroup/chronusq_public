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

// KCaKrK sto-3g test for fchk parsing
TEST( GHF, KCaKr_sto3G ) {

  CQSCFTEST( "scf/serial/ghf/KCaKr_sto-3G", "KCaKr_sto-3G.bin.ref", 2e-8,
      true, true, true, true, true, true,
      false, "KCaKr_sto-3G.fchk" );

};

// KCaKr sto-3g spin and angular momentum property test
TEST( GHF, KCaKr_sto3G_SpinAngular ) {

  CQSCFTEST( "scf/serial/ghf/KCaKr_sto-3G", "KCaKr_sto-3G.bin.ref", 2e-8,
      true, true, true, true, true, true,
      false, "KCaKr_sto-3G.fchk", false,
      true, true, true, true, true );

};

// B sto-3g test for MO swapping
TEST( GHF, B_swap_GHF_sto3G ) {

  CQSCFTEST( "scf/serial/ghf/B_swap_GHF_sto-3g", "B_swap_GHF_sto-3g.bin.ref",
    2e-8, true, true, true, true, true, true,
    true );

};

// N GHF/sto-3g with UHF guess (READMO)
TEST( GHF, N_STO3G_GHF_UHFGuess_READMO ) {

  CQSCFTEST( "scf/serial/ghf/N_STO-3G_GHF_READMO",
    "N_STO-3G_GHF_UHFGuess_READMO.bin.ref",1e-8,
    false, false, false, false, false, true, false,
    "N_STO-3G_UHF.scr.bin" );

};

TEST( GHF, Water_ccpVDZ_GHF_HOMO2_LUMO2_MOM ) {

  CQSCFTEST( "scf/serial/ghf/water_cc-pVDZ_GHF_HOMO2_LUMO2_MOM",
    "water_cc-pVDZ_GHF_HOMO2_LUMO2_MOM.bin.ref",1e-6,
    false, false, false, false, false, true, false, "water_cc-pVDZ_GHF.scr.bin");

};


#ifdef _CQ_DO_PARTESTS


#endif
