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

#ifndef _CQ_GENERATE_TESTS

#include "resp.hpp"

#define CQMORTEST_IMPL(TNAME, HER, AHER, IN, REF) \
TEST( RHF_MOR, TNAME ){ CQMORTEST(HER,AHER,IN,REF); }


// Full Dimensional Tests

// Water 6-31G(d) MOR
CQMORTEST_IMPL( Water_631Gd_MOR,  true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor",
    "water_6-31Gd_rhf_interp.bin.ref" )


// Water 6-31G(d) MOR (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_GMRES,  true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR (GMRES + DIRECT)
CQMORTEST_IMPL( Water_631Gd_MOR_GMRES_DIRECT,  true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres_direct",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR,  true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES,  true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES + DIRECT)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES_DIRECT,  true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres_direct",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif





// A+B / A-B TESTS
   
// Water 6-31G(d) MOR A+B / A-B
CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB,  true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb",
    "water_6-31Gd_rhf_interp.bin.ref" )


// Water 6-31G(d) MOR A+B / A-B (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_GMRES,  true, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR A+B /A-B
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB,  true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B /A-B (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_GMRES,  
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif



// REDUCED HERMETIAN TESTS 
   
// Water 6-31G(d) MOR REDUCED HERMETIAN 
CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER,  true, false, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_GMRES,  true, false, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER,  true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_GMRES,  
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif


#if 0

// REDUCED ANTI-HERMETIAN TESTS 
   
// Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN 
CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER,  false, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher",
    "water_6-31Gd_rhf_interp.bin.ref" )

// Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN (GMRES)
CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_GMRES,  false, true, 
    "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres",
    "water_6-31Gd_rhf_interp.bin.ref" )

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER,  
      false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTI-HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_GMRES,  
      false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres",
      "water_6-31Gd_rhf_interp.bin.ref" )

#endif

#endif





#if defined(CQ_ENABLE_MPI)

  // Distributed Full Dimensional Tests

  // Water 6-31G(d) MOR
  CQMORTEST_IMPL( Water_631Gd_MOR_DISTMATFROMROOT,  true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )




  // Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) MOR A+B / A-B
  CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B / A-B
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR A+B / A-B (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_APB_AMB_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR A+B / A-B (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_APB_AMB_GMRES_DISTMATFROMROOT, 
      true, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )



  // Distributed Reduced Hermetian Tests

  // Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_HER_GMRES_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED HERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_HER_GMRES_DISTMATFROMROOT, 
      true, false, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

#if 0


  // Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) MOR REDUCED ANTIHERMETIAN
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_DISTMATFROMROOT,  
      false, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTIHERMETIAN
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // Water 6-31G(d) MOR REDUCED ANTIHERMETIAN (GMRES)
  CQMORTEST_IMPL( Water_631Gd_MOR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      false, true, 
      "resp/serial/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )

  // SMP Water 6-31G(d) MOR REDUCED ANTIHERMETIAN (GMRES)
  CQMORTEST_IMPL( PAR_Water_631Gd_MOR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rmor/water_6-31Gd_rhf_mor_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_interp.bin.ref" )


#endif

#endif


#endif
