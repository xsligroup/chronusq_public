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

#include "mcscf.hpp"


// Al 6-31G(d) test

TEST(RASCI_DAVIDSON, Al_631G ) {

  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_x2c_rasci_davidson", "al_6-31G_x2c_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_4c_bc_rasci_davidson", "al_6-31G_4c_bc_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_4c_dc_rasci_davidson", "al_6-31G_4c_dc_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_4c_dcssss_rasci_davidson", "al_6-31G_4c_dcssss_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_4c_dcg_rasci_davidson", "al_6-31G_4c_dcg_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/serial/ras/al_6-31G_4c_dcb_rasci_davidson", "al_6-31G_4c_dcb_rasci.bin.ref" );
 
};

#ifndef _CQ_GENERATE_TESTS
#ifdef _CQ_DO_PARTESTS

// SMP Al 6-31G(d) test

TEST(RASCI_DAVIDSON, PAR_Al_631G ) {

  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_x2c_rasci_davidson", "al_6-31G_x2c_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_4c_bc_rasci_davidson", "al_6-31G_4c_bc_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_4c_dc_rasci_davidson", "al_6-31G_4c_dc_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_4c_dcssss_rasci_davidson", "al_6-31G_4c_dcssss_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_4c_dcg_rasci_davidson", "al_6-31G_4c_dcg_rasci.bin.ref" );
  CQMCSCFTEST( "mcscf/parallel/ras/al_6-31G_4c_dcb_rasci_davidson", "al_6-31G_4c_dcb_rasci.bin.ref" );
 
};
#endif
#endif



