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
 *    E-Mail: xsli@@uw.edu
 *  
 */

#include "resp.hpp"

#define CQRESTEST_IMPL(TNAME, IN, REF) \
TEST( UHF_PP_RESIDUE, TNAME ) { CQRESTEST( false, IN, REF ); }


// Li cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( Li_ccpVDZ_ppRPA_AA_RESIDUE,  
    "resp/serial/upp/Li_ccpVDZ_AA_ppRPA",
    "Li_ccpVDZ_AA_ppRPA.bin.ref" )

// Li cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( Li_ccpVDZ_ppRPA_AB_RESIDUE,  
    "resp/serial/upp/Li_ccpVDZ_AB_ppRPA",
    "Li_ccpVDZ_AB_ppRPA.bin.ref" )

// Li cc-pVDZ PP-RPA-HF (BB) (RESIDUE)
CQRESTEST_IMPL( Li_ccpVDZ_ppRPA_BB_RESIDUE,  
    "resp/serial/upp/Li_ccpVDZ_BB_ppRPA",
    "Li_ccpVDZ_BB_ppRPA.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Li cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_Li_ccpVDZ_ppRPA_AA_RESIDUE,  
      "resp/parallel/upp/Li_ccpVDZ_AA_ppRPA",
      "Li_ccpVDZ_AA_ppRPA.bin.ref" )
  
  // SMP Li cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_Li_ccpVDZ_ppRPA_AB_RESIDUE,  
      "resp/parallel/upp/Li_ccpVDZ_AB_ppRPA",
      "Li_ccpVDZ_AB_ppRPA.bin.ref" )

  // SMP Li cc-pVDZ PP-RPA-HF (BB) (RESIDUE)
  CQRESTEST_IMPL( PAR_Li_ccpVDZ_ppRPA_BB_RESIDUE,  
      "resp/parallel/upp/Li_ccpVDZ_BB_ppRPA",
      "Li_ccpVDZ_BB_ppRPA.bin.ref" )

#endif
