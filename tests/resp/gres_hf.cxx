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

#include "resp.hpp"

#define CQCRESTEST_IMPL(TNAME, IN, REF) \
TEST( GHF_RESIDUE, TNAME ) { CQCRESTEST( true, IN, REF ); }

#define CQCRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GHF_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF ); }

#define CQCRESREFTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GHF_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF, true, 1e-06, true ); }

#define CQCRESFCHKTEST_IMPL_NP(TNAME, IN, REF, FCHK) \
TEST( GHF_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF, true, 1e-06, false, FCHK ); }

// FULL DIMENSIONAL TESTS

// H2O 3-21G TD-GHF (RESIDUE)
CQCRESTEST_IMPL_NP( H2O_321G_GHF_GTO_RESIDUE,
    "resp/serial/gresp_hf/H2O_3-21G_ghf_gto_residue",
    "H2O_3-21G_ghf_gto_residue.bin.ref" )

#ifdef _CQ_DO_PARTESTS

// H2O 3-21G TD-GHF (RESIDUE)
CQCRESTEST_IMPL_NP( PAR_H2O_321G_GHF_GTO_RESIDUE,
    "resp/parallel/gresp_hf/H2O_3-21G_ghf_gto_residue",
    "H2O_3-21G_ghf_gto_residue.bin.ref" )

#endif
