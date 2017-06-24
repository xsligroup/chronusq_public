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
TEST( GLSDA_RESIDUE, TNAME ) { CQCRESTEST( true, IN, REF ); }

#define CQCRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GLSDA_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF ); }

#define CQCRESREFTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GLSDA_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF, true, 1e-06, true ); }

#define CQCRESFCHKTEST_IMPL_NP(TNAME, IN, REF, FCHK) \
TEST( GLSDA_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF, true, 1e-06, false, FCHK ); }

// FULL DIMENSIONAL TESTS

// H2O 3-21G TD-GLSDA (RESIDUE)
CQCRESREFTEST_IMPL_NP( H2O_321G_GLSDA_GTO_RESIDUE,
    "resp/serial/gresp_ks/H2O_3-21G_glsda_gto_residue",
    "H2O_3-21G_glsda_gto_residue.bin.ref" )

// H3 3-21g TD-GLSDA (RESIDUE)
CQCRESFCHKTEST_IMPL_NP( H3_321g_GLSDA_GTO_RESIDUE,
    "resp/serial/gresp_ks/H3_3-21g_glsda_gto_residue",
    "H3_3-21g_glsda_gto_residue.bin.ref", 
    "H3_3-21g_glsda_gto_residue.fchk" )

#ifdef _CQ_DO_PARTESTS

// H2O 3-21G TD-GLSDA (RESIDUE)
CQCRESREFTEST_IMPL_NP( PAR_H2O_321G_GLSDA_GTO_RESIDUE,
    "resp/parallel/gresp_ks/H2O_3-21G_glsda_gto_residue",
    "H2O_3-21G_glsda_gto_residue.bin.ref" )

// H3 3-21g TD-GLSDA (RESIDUE)
CQCRESFCHKTEST_IMPL_NP( PAR_H3_321g_GLSDA_GTO_RESIDUE,
    "resp/parallel/gresp_ks/H3_3-21g_glsda_gto_residue",
    "H3_3-21g_glsda_gto_residue.bin.ref", 
    "H3_3-21g_glsda_gto_residue.fchk" )


#endif
