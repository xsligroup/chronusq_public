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
TEST( RHF_GIAO_RESIDUE, TNAME ) { CQRESTEST( true, IN, REF ); }

#define CQRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( RHF_GIAO_RESIDUE, TNAME ) { CQRESTEST( false, IN, REF ); }

// FULL DIMENSIONAL TESTS

// Water STO3G+GIAO TD-NEOHF (RESIDUE)
CQRESTEST_IMPL_NP( H2O_STO3G_GIAO_RESIDUE,
    "resp/serial/rresp_giao/water_sto3g_giao_residue",
    "water_sto3g_giao_residue.bin.ref" )

#ifdef _CQ_DO_PARTESTS

// HCN STO-3G/SP TD-NEOHF (RESIDUE)
CQRESTEST_IMPL_NP( PAR_H2O_STO3G_GIAO_RESIDUE,
    "resp/parallel/rresp_giao/water_sto3g_giao_residue",
    "water_sto3g_giao_residue.bin.ref" )

#endif

