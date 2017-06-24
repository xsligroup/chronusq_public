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
TEST( GHGGA_RESIDUE, TNAME ) { CQCRESTEST( true, IN, REF ); }

#define CQCRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GHGGA_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF ); }

#define CQCRESTEST_IMPL_NP_MPI(TNAME, IN, REF) \
TEST( GHGGA_RESIDUE_MPI, TNAME ) { CQCRESTEST( false, IN, REF ); }

#define CQCRESREFTEST_IMPL_NP(TNAME, IN, REF) \
TEST( GHGGA_RESIDUE, TNAME ) { CQCRESTEST( false, IN, REF, true, 1e-06, true ); }

// FULL DIMENSIONAL TESTS

// H2O 3-21G TD-GB3LYP (RESIDUE)
CQCRESTEST_IMPL_NP( H2O_321G_GB3LYP_GTO_RESIDUE,
    "resp/serial/gresp_ks/H2O_3-21G_gb3lyp_gto_residue",
    "H2O_3-21G_gb3lyp_gto_residue.bin.ref" )

// O2 3-21G TD-GB3LYP (RESIDUE)
CQCRESREFTEST_IMPL_NP( O2_321G_GB3LYP_GTO_RESIDUE,
    "resp/serial/gresp_ks/O2_3-21G_gb3lyp_gto_residue",
    "O2_3-21G_gb3lyp_gto_residue.bin.ref" )

// GPLHR 
 
// H2O 6-31G TD-GB3LYP (RESIDUE)
CQCRESTEST_IMPL_NP( H2O_631G_GB3LYP_GPLHR_RESIDUE,
    "resp/serial/gresp_ks/H2O_6-31G_gb3lyp_gplhr_residue",
    "H2O_6-31G_gb3lyp_gplhr_residue.bin.ref" )


#ifdef _CQ_DO_PARTESTS

// H2O 3-21G TD-GB3LYP (RESIDUE)
CQCRESTEST_IMPL_NP( PAR_H2O_321G_GB3LYP_GTO_RESIDUE,
    "resp/parallel/gresp_ks/H2O_3-21G_gb3lyp_gto_residue",
    "H2O_3-21G_gb3lyp_gto_residue.bin.ref" )

// O2 3-21G TD-GB3LYP (RESIDUE)
CQCRESREFTEST_IMPL_NP( PAR_O2_321G_GB3LYP_GTO_RESIDUE,
    "resp/parallel/gresp_ks/O2_3-21G_gb3lyp_gto_residue",
    "O2_3-21G_gb3lyp_gto_residue.bin.ref" )

// H2O 6-31G TD-GB3LYP (RESIDUE)
CQCRESTEST_IMPL_NP( PAR_H2O_631G_GB3LYP_GPLHR_RESIDUE,
    "resp/parallel/gresp_ks/H2O_6-31G_gb3lyp_gplhr_residue",
    "H2O_6-31G_gb3lyp_gplhr_residue.bin.ref" )

#endif

#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

// H2O 6-31G TD-GB3LYP (RESIDUE)
CQCRESTEST_IMPL_NP_MPI( H2O_631G_GB3LYP_DISTMATFROMROOT,
    "resp/parallel/gresp_ks/H2O_6-31G_gb3lyp_distmatfromroot",
    "H2O_6-31G_gb3lyp_gplhr_residue.bin.ref" )

#endif
