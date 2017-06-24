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
TEST( GHF_GIAO, H2_TRIPLET_GHF_GIAO_631G ) {

  CQSCFTEST( 
      "scf/serial/ghf_giao/h2_triplet_ghf_giao_631G", 
      "h2_triplet_ghf_giao_631G.bin.ref", 1e-6,
      false, false, false, false, false);
 
};

#ifdef _CQ_DO_PARTESTS

// SMP H2 TRIPLET GIAO 6-31G test  B = 0, 0, -0.001
//TEST( GHF_GIAO, Par_H2_TRIPLET_GHF_GIAO_631G ) {
//
//  CQSCFTEST( 
//      "scf/parallel/ghf_giao/h2_triplet_ghf_giao_631G", 
//      "h2_triplet_ghf_giao_631G.bin.ref", 9e-8,
//      false, false, false, false, false);
// 
//};

#endif
