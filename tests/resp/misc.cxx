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
TEST( MISC_RESP, TNAME ) { CQRESTEST( false, IN, REF, true, 1e-5 ); }


// Water 6-31G(d) TDHF (RESIDUE, GPLHR + DIRECT) Oxygen K-Edge (540 eV)
CQRESTEST_IMPL( Water_631Gd_RESIDUE_XRAY_540eV,  
    "resp/serial/misc/water_6-31Gd_rhf_residue_xray_540eV",
    "water_6-31Gd_rhf_residue.bin.ref" )

