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
TEST( RHF_PP_RESIDUE, TNAME ) { CQRESTEST( false, IN, REF, true, 1e-5); }


// AA PP-RPA

// BH cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AA_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AA_ppRPA",
    "BH_ccpVDZ_AA_ppRPA.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // BH cc-pVDZ PP-RPA-HF (AA) (RESIDUE, GPLHR)
  CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AA_RESIDUE_GPLHR,  
      "resp/serial/rpp/BH_ccpVDZ_AA_ppRPA_gplhr",
      "BH_ccpVDZ_AA_ppRPA.bin.ref" )

  // BH cc-pVDZ PP-RPA-HF (AA) (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AA_RESIDUE_GPLHR_DIRECT,  
      "resp/serial/rpp/BH_ccpVDZ_AA_ppRPA_gplhr_direct",
      "BH_ccpVDZ_AA_ppRPA.bin.ref" )

#endif





// AB PP-RPA

// BH cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AB_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AB_ppRPA",
    "BH_ccpVDZ_AB_ppRPA.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // BH cc-pVDZ PP-RPA-HF (AB) (RESIDUE, GPLHR)
  CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AB_RESIDUE_GPLHR,  
      "resp/serial/rpp/BH_ccpVDZ_AB_ppRPA_gplhr",
      "BH_ccpVDZ_AB_ppRPA.bin.ref" )

  // BH cc-pVDZ PP-RPA-HF (AB) (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( BH_ccpVDZ_ppRPA_AB_RESIDUE_GPLHR_DIRECT,  
      "resp/serial/rpp/BH_ccpVDZ_AB_ppRPA_gplhr_direct",
      "BH_ccpVDZ_AB_ppRPA.bin.ref" )

#endif





// AA PP-TDA


// BH cc-pVDZ PP-TDA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AA_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AA_ppTDA",
    "BH_ccpVDZ_AA_ppTDA.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // BH cc-pVDZ PP-TDA-HF (AA) (RESIDUE, GPLHR)
  CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AA_RESIDUE_GPLHR,  
      "resp/serial/rpp/BH_ccpVDZ_AA_ppTDA_gplhr",
      "BH_ccpVDZ_AA_ppTDA.bin.ref" )

  // BH cc-pVDZ PP-TDA-HF (AA) (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AA_RESIDUE_GPLHR_DIRECT,  
      "resp/serial/rpp/BH_ccpVDZ_AA_ppTDA_gplhr_direct",
      "BH_ccpVDZ_AA_ppTDA.bin.ref" )

#endif











// AB PP-TDA
  
  
// BH cc-pVDZ PP-TDA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AB_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AB_ppTDA",
    "BH_ccpVDZ_AB_ppTDA.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // BH cc-pVDZ PP-TDA-HF (AB) (RESIDUE, GPLHR)
  CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AB_RESIDUE_GPLHR,  
      "resp/serial/rpp/BH_ccpVDZ_AB_ppTDA_gplhr",
      "BH_ccpVDZ_AB_ppTDA.bin.ref" )

  // BH cc-pVDZ PP-TDA-HF (AB) (RESIDUE, GPLHR + DIRECT)
  CQRESTEST_IMPL( BH_ccpVDZ_ppTDA_AB_RESIDUE_GPLHR_DIRECT,  
      "resp/serial/rpp/BH_ccpVDZ_AB_ppTDA_gplhr_direct",
      "BH_ccpVDZ_AB_ppTDA.bin.ref" )

#endif






// HH-TDA



// BH cc-pVDZ HH-TDA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_hhTDA_AA_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AA_hhTDA",
    "BH_ccpVDZ_AA_hhTDA.bin.ref" )

// BH cc-pVDZ HH-TDA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( BH_ccpVDZ_hhTDA_AB_RESIDUE,  
    "resp/serial/rpp/BH_ccpVDZ_AB_hhTDA",
    "BH_ccpVDZ_AB_hhTDA.bin.ref" )






// O2 cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_ppRPA_AA_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AA_ppRPA",
    "O2_ccpVDZ_AA_ppRPA.bin.ref" )

// O2 cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_ppRPA_AB_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AB_ppRPA",
    "O2_ccpVDZ_AB_ppRPA.bin.ref" )

// O2 cc-pVDZ PP-TDA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_ppTDA_AA_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AA_ppTDA",
    "O2_ccpVDZ_AA_ppTDA.bin.ref" )

// O2 cc-pVDZ PP-TDA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_ppTDA_AB_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AB_ppTDA",
    "O2_ccpVDZ_AB_ppTDA.bin.ref" )

// O2 cc-pVDZ HH-TDA-HF (AA) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_hhTDA_AA_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AA_hhTDA",
    "O2_ccpVDZ_AA_hhTDA.bin.ref" )

// O2 cc-pVDZ HH-TDA-HF (AB) (RESIDUE)
CQRESTEST_IMPL( O2_ccpVDZ_hhTDA_AB_RESIDUE,  
    "resp/serial/rpp/O2_ccpVDZ_AB_hhTDA",
    "O2_ccpVDZ_AB_hhTDA.bin.ref" )














#ifdef _CQ_DO_PARTESTS

  // SMP BH cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_ppRPA_AA_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AA_ppRPA",
      "BH_ccpVDZ_AA_ppRPA.bin.ref" )
  
  // SMP BH cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_ppRPA_AB_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AB_ppRPA",
      "BH_ccpVDZ_AB_ppRPA.bin.ref" )
  
  // SMP BH cc-pVDZ PP-TDA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_ppTDA_AA_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AA_ppTDA",
      "BH_ccpVDZ_AA_ppTDA.bin.ref" )
  
  // SMP BH cc-pVDZ PP-TDA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_ppTDA_AB_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AB_ppTDA",
      "BH_ccpVDZ_AB_ppTDA.bin.ref" )
  
  
  // SMP BH cc-pVDZ HH-TDA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_hhTDA_AA_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AA_hhTDA",
      "BH_ccpVDZ_AA_hhTDA.bin.ref" )
  
  // SMP BH cc-pVDZ HH-TDA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_BH_ccpVDZ_hhTDA_AB_RESIDUE,  
      "resp/parallel/rpp/BH_ccpVDZ_AB_hhTDA",
      "BH_ccpVDZ_AB_hhTDA.bin.ref" )






  // SMP O2 cc-pVDZ PP-RPA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_ppRPA_AA_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AA_ppRPA",
      "O2_ccpVDZ_AA_ppRPA.bin.ref" )
  
  // SMP O2 cc-pVDZ PP-RPA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_ppRPA_AB_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AB_ppRPA",
      "O2_ccpVDZ_AB_ppRPA.bin.ref" )

  // SMP O2 cc-pVDZ PP-TDA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_ppTDA_AA_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AA_ppTDA",
      "O2_ccpVDZ_AA_ppTDA.bin.ref" )
  
  // SMP O2 cc-pVDZ PP-TDA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_ppTDA_AB_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AB_ppTDA",
      "O2_ccpVDZ_AB_ppTDA.bin.ref" )

  // SMP O2 cc-pVDZ HH-TDA-HF (AA) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_hhTDA_AA_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AA_hhTDA",
      "O2_ccpVDZ_AA_hhTDA.bin.ref" )
  
  // SMP O2 cc-pVDZ HH-TDA-HF (AB) (RESIDUE)
  CQRESTEST_IMPL( PAR_O2_ccpVDZ_hhTDA_AB_RESIDUE,  
      "resp/parallel/rpp/O2_ccpVDZ_AB_hhTDA",
      "O2_ccpVDZ_AB_hhTDA.bin.ref" )


#endif

