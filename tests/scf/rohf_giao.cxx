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

// H2 TRIPLET GIAO 6-31G test  B = 0, 0, -0.001
TEST( ROHF_GIAO, H2_TRIPLET_ROHF_GIAO_631G ) {

  CQSCFTEST( "scf/serial/rohf_giao/h2_triplet_rohf_giao_631G", "h2_triplet_rohf_giao_631G.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true);
 
};


#ifndef _CQ_GENERATE_TESTS

// H2 TRIPLET GIAO 6-31G test  B = 0, 0, -0.001 (DIRECT)
TEST( ROHF_GIAO, H2_TRIPLET_ROHF_GIAO_631G_DIRECT ) {

  CQSCFTEST( "scf/serial/rohf_giao/h2_triplet_rohf_giao_631G_direct", "h2_triplet_rohf_giao_631G.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

#endif

#ifdef _CQ_DO_PARTESTS

// SMP H2 TRIPLET GIAO 6-31G test  B = 0, 0, -0.001
TEST( ROHF_GIAO, Par_H2_TRIPLET_ROHF_GIAO_631G ) {

  CQSCFTEST( "scf/parallel/rohf_giao/h2_triplet_rohf_giao_631G", "h2_triplet_rohf_giao_631G.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true);
 
};


// SMP H2 TRIPLET GIAO 6-31G test  B = 0, 0, -0.001 (DIRECT)
TEST( ROHF_GIAO, Par_H2_TRIPLET_ROHF_GIAO_631G_DIRECT ) {

  CQSCFTEST( "scf/parallel/rohf_giao/h2_triplet_rohf_giao_631G_direct", "h2_triplet_rohf_giao_631G.bin.ref", 1e-6, 
      true, true, true, true, true, true, false, "no", true );
 
};

#endif
