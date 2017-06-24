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

#define CQRESTEST_IMPL(TNAME, IN, REF) \
TEST( UHF_RESIDUE, TNAME ) { CQRESTEST( true, IN, REF ); }

#define CQRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( UHF_RESIDUE, TNAME ) { CQRESTEST( false, IN, REF ); }

// FULL DIMENSIONAL TESTS

// Na STO-3G TD-UHF (RESIDUE)
CQRESTEST_IMPL_NP( Na_STO3G_UHF_GTO_RESIDUE,
    "resp/serial/uresp_hf/Na_STO-3G_uhf_gto_residue",
    "Na_STO-3G_uhf_gto_residue.bin.ref" )

// H2O 3-21G TD-UHF (RESIDUE)
CQRESTEST_IMPL_NP( H2O_321G_UHF_GTO_RESIDUE,
    "resp/serial/uresp_hf/H2O_3-21G_uhf_gto_residue",
    "H2O_3-21G_uhf_gto_residue.bin.ref" )

#ifdef _CQ_DO_PARTESTS

// Na STO-3G TD-UHF (RESIDUE)
CQRESTEST_IMPL_NP( PAR_Na_STO3G_UHF_GTO_RESIDUE,
    "resp/parallel/uresp_hf/Na_STO-3G_uhf_gto_residue",
    "Na_STO-3G_uhf_gto_residue.bin.ref" )

// H2O 3-21G TD-UHF (RESIDUE)
CQRESTEST_IMPL_NP( PAR_H2O_321G_UHF_GTO_RESIDUE,
    "resp/parallel/uresp_hf/H2O_3-21G_uhf_gto_residue",
    "H2O_3-21G_uhf_gto_residue.bin.ref" )

#endif
