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

// BH GIAO STO-3G test  B = 0, 0, 0.02
TEST( RHF_GIAO, BH_GIAO_STO3G ) {

  CQSCFTEST( "scf/serial/rhf_giao/BH_giao_sto3g", "BH_giao_sto3g.bin.ref", 1e-6 );
 
};

#ifndef _CQ_GENERATE_TESTS

// BH GIAO STO-3G test  B = 0, 0, 0.02 (DIRECT)
TEST( RHF_GIAO, BH_GIAO_STO3G_DIRECT ) {

  CQSCFTEST( "scf/serial/rhf_giao/BH_giao_sto3g_direct", "BH_giao_sto3g.bin.ref", 1e-6 );
 
};

#endif



#ifdef _CQ_DO_PARTESTS

// SMP BH GIAO STO-3G test  B = 0, 0, 0.02
TEST( RHF_GIAO, Par_BH_GIAO_STO3G ) {

  CQSCFTEST( "scf/parallel/rhf_giao/BH_giao_sto3g", "BH_giao_sto3g.bin.ref", 1e-6 );
 
};


// SMP BH GIAO STO-3G test  B = 0, 0, 0.02 (DIRECT)
TEST( RHF_GIAO, Par_BH_GIAO_STO3G_DIRECT ) {

  CQSCFTEST( "scf/parallel/rhf_giao/BH_giao_sto3g_direct", "BH_giao_sto3g.bin.ref", 1e-6 );
 
};

#endif
