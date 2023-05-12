/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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


//NEO-DFT with minimal basis set, using epc17 functional
TEST( NEO_RKS, water_sto3g_protsp_rb3lyp_uepc17) {

  CQNEOSCFTEST( "scf/serial/neo_rks/water_sto-3g_prot-sp_rb3lyp_uepc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref" );
 
}


//NEO-DFT with minimal basis set, using epc19 functional
TEST( NEO_RKS, water_sto3g_protsp_rb3lyp_uepc19) {

  CQNEOSCFTEST( "scf/serial/neo_rks/water_sto-3g_prot-sp_rb3lyp_uepc19", "water_sto-3g_prot-sp_rb3lyp_uepc19.bin.ref" );
 
}


#ifdef _CQ_DO_PARTESTS

//NEO-DFT with minimal basis set, using epc17 functional, parallel job
TEST( NEO_RKS, par_water_sto3g_protsp_rb3lyp_uepc17) {

  CQNEOSCFTEST( "scf/parallel/neo_rks/par_water_sto-3g_prot-sp_rb3lyp_uepc17", "water_sto-3g_prot-sp_rb3lyp_uepc17.bin.ref" );
 
}


// FIXME: for epc19, parallel jobs doesn't pass. Need to debug epc19
//
// //NEO-DFT with minimal basis set, using epc19 functional, parallel job
// TEST( NEO_RKS, par_water_sto3g_protsp_rb3lyp_uepc19) {
// 
//   CQNEOSCFTEST( "scf/parallel/neo_rks/par_water_sto-3g_prot-sp_rb3lyp_uepc19", "water_sto-3g_prot-sp_rb3lyp_uepc19.bin.ref" );
//  
// }

#endif




