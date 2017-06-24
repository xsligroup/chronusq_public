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

#define CQRESTEST_IMPL_NP(TNAME, IN, REF) \
TEST( X2CHF_GIAO_RESIDUE, TNAME ) { CQRESTEST( false, IN, REF ); }

// RESIDUE TESTS

// O2/631G GIAO TD-NEOHF (RESIDUE)
CQRESTEST_IMPL_NP( O2_X2C_GIAO_RESIDUE,
    "resp/serial/x2cresp_giao_hf/o2_triplet_631g_x2c_giao_residue",
    "o2_triplet_631g_x2c_giao_residue.bin.ref" )

#ifdef _CQ_DO_PARTESTS

// O2/631G GIAO TD-NEOHF (RESIDUE)
CQRESTEST_IMPL_NP( PAR_O2_X2C_GIAO_RESIDUE,
    "resp/parallel/x2cresp_giao_hf/o2_triplet_631g_x2c_giao_residue",
    "o2_triplet_631g_x2c_giao_residue.bin.ref" )

#endif

