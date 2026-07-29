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

#include "dasscf.hpp"


// Al 6-31G(d) test

TEST(DAS_RASCI_DAVIDSON, Al_631G ) {

  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_x2c_dasci_excrestrictions_davidson", "al_6-31G_x2c_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_4c_bc_dasci_excrestrictions_davidson", "al_6-31G_4c_bc_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_4c_dc_dasci_excrestrictions_davidson", "al_6-31G_4c_dc_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_4c_dcssss_dasci_excrestrictions_davidson", "al_6-31G_4c_dcssss_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_4c_dcg_dasci_excrestrictions_davidson", "al_6-31G_4c_dcg_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/serial/ras/al_6-31G_4c_dcb_dasci_excrestrictions_davidson", "al_6-31G_4c_dcb_dasci_excrestrictions.bin.ref" );
 
};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

// SMP Al 6-31G(d) test

TEST(DAS_RASCI_DAVIDSON, PAR_Al_631G ) {

  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_x2c_dasci_excrestrictions_davidson", "al_6-31G_x2c_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_4c_bc_dasci_excrestrictions_davidson", "al_6-31G_4c_bc_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_4c_dc_dasci_excrestrictions_davidson", "al_6-31G_4c_dc_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_4c_dcssss_dasci_excrestrictions_davidson", "al_6-31G_4c_dcssss_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_4c_dcg_dasci_excrestrictions_davidson", "al_6-31G_4c_dcg_dasci_excrestrictions.bin.ref" );
  CQDASSCFTEST( "dasscf/parallel/ras/al_6-31G_4c_dcb_dasci_excrestrictions_davidson", "al_6-31G_4c_dcb_dasci_excrestrictions.bin.ref" );

};
#endif
#endif



