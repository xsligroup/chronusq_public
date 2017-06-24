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

// B sto-3g test for MO swapping
TEST( GHF, B_swap_GHF_sto3G ) {

  CQSCFTEST( "scf/serial/ghf/B_swap_GHF_sto-3g", "B_swap_GHF_sto-3g.bin.ref",
    2e-8, true, true, true, true, true, true,
    true );

};

#ifdef _CQ_DO_PARTESTS


#endif
