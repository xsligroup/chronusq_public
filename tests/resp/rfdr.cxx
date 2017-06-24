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


#define CQFDRTEST_IMPL(TNAME, HER, AHER, IN, REF) \
TEST( RHF_FDR, TNAME ){ CQFDRTEST<double>(HER,AHER,IN,REF); }

#define CQDFDRTEST_IMPL(TNAME, HER, AHER, IN, REF) \
TEST( RHF_DFDR, TNAME ){ CQFDRTEST<dcomplex>(HER,AHER,IN,REF); }


// FULL DIMENSIONAL TESTS
  
// Water 6-31G(d) FDR
CQFDRTEST_IMPL( Water_631Gd_FDR,  true, true, 
    "resp/serial/rresp/water_6-31Gd_rhf_fdr",
    "water_6-31Gd_rhf_fdr.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR (GMRES + DIRECT)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES_DIRECT,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres_direct",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) FDR
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES + DIRECT)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES_DIRECT,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres_direct",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif


#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
  
  // Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb",
      "water_6-31Gd_rhf_fdr.bin.ref" )


  // Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_GMRES,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR A+B / A-B
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB,  true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb",
        "water_6-31Gd_rhf_fdr.bin.ref" )


    // SMP Water 6-31G(d) FDR A+B / A-B (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_GMRES,  true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif



  // REDUCED HERMETIAN TESTS
  
  // Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER,  true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her",
      "water_6-31Gd_rhf_fdr.bin.ref" )


  // Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_GMRES,  true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR REDUCED HERMETIAN
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER,  true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her",
        "water_6-31Gd_rhf_fdr.bin.ref" )

    // SMP Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_GMRES,  
        true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif




  // REDUCED ANTI-HERMETIAN TESTS
  
  // Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER,  false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_GMRES,  false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres",
      "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER,  false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher",
        "water_6-31Gd_rhf_fdr.bin.ref" )

    // SMP Water 6-31G(d) FDR REDUCED ANTI-HERMETIAN (GMRES)
    CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_GMRES,  
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres",
        "water_6-31Gd_rhf_fdr.bin.ref" )
  
  #endif

#endif


#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

  // Distributed Full Dimensional Tests

  // Water 6-31G(d) FDR
  CQFDRTEST_IMPL( Water_631Gd_FDR_DISTMATFROMROOT,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_DISTMATFROMROOT,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_GMRES_DISTMATFROMROOT,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )




  // Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR A+B / A-B
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_APB_AMB_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR A+B / A-B (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_APB_AMB_GMRES_DISTMATFROMROOT, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )



  // Distributed Reduced Hermetian Tests

  // Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED HERMETIAN
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_HER_GMRES_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED HERMETIAN (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_HER_GMRES_DISTMATFROMROOT, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )



  // Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) FDR REDUCED ANTIHERMETIAN
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_DISTMATFROMROOT,  
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED ANTIHERMETIAN
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // Water 6-31G(d) FDR REDUCED ANTIHERMETIAN (GMRES)
  CQFDRTEST_IMPL( Water_631Gd_FDR_RED_ANTIHER_GMRES_DISTMATFROMROOT,  
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

  // SMP Water 6-31G(d) FDR REDUCED ANTIHERMETIAN (GMRES)
  CQFDRTEST_IMPL( PAR_Water_631Gd_FDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_fdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_fdr.bin.ref" )

#endif




























// Full Dimensional Tests

// Water 6-31G(d) DFDR
CQDFDRTEST_IMPL( Water_631Gd_DFDR,  true, true, 
    "resp/serial/rresp/water_6-31Gd_rhf_dfdr",
    "water_6-31Gd_rhf_dfdr.bin.ref" )

#ifndef _CQ_GENERATE_TESTS

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES + DIRECT)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_DIRECT,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_direct",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif

#ifdef _CQ_DO_PARTESTS

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES + DIRECT)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_DIRECT,  true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_direct",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif




#ifndef _CQ_GENERATE_TESTS

  // A+B / A-B TESTS
     
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR A+B /A-B
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB,  true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR A+B /A-B (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES,  
        true, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif


  // REDUCED HERMETIAN TESTS 
     
  // Water 6-31G(d) DFDR REDUCED HERMETIAN 
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER,  true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES,  true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER,  true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES,  
        true, false, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif



  // REDUCED ANTI-HERMETIAN TESTS 
     
  // Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN 
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER,  false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES,  false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres",
      "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #ifdef _CQ_DO_PARTESTS
  
    // SMP Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER,  
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher",
        "water_6-31Gd_rhf_dfdr.bin.ref" )

    // SMP Water 6-31G(d) DFDR REDUCED ANTI-HERMETIAN (GMRES)
    CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES,  
        false, true, 
        "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres",
        "water_6-31Gd_rhf_dfdr.bin.ref" )
  
  #endif


#endif



#if defined(CQ_ENABLE_MPI) && !defined(_CQ_GENERATE_TESTS)

  // Distributed-From-Root Full Dimensional Tests

  // Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_DISTMATFROMROOT,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )




  // Distributed-From-Root A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_DISTMATFROMROOT,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES_DISTMATFROMROOT,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES_DISTMATFROMROOT, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Distributed-From-Root Reduced Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_DISTMATFROMROOT,  
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES_DISTMATFROMROOT,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES_DISTMATFROMROOT, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Distributed-From-Root Reduced Anti-Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_DISTMATFROMROOT,  
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES_DISTMATFROMROOT, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_distfromroot",
      "water_6-31Gd_rhf_dfdr.bin.ref" )












  // Build-Distributed Full Dimensional Tests

  // Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_BUILDDIST,  true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_BUILDDIST,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_GMRES_BUILDDIST,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_GMRES_BUILDDIST,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )




  // Build-Distributed A+B / A-B Dimensional Tests
    
  // Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_BUILDDIST,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_BUILDDIST,  
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_APB_AMB_GMRES_BUILDDIST,  
      true, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR A+B / A-B (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_APB_AMB_GMRES_BUILDDIST, 
      true, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_apb_amb_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Build-Distributed Reduced Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_BUILDDIST,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_BUILDDIST,  
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_HER_GMRES_BUILDDIST,  
      true, false, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED HERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_HER_GMRES_BUILDDIST, 
      true, false, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_her_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )



  // Build-Distributed Reduced Anti-Hermetian Tests

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_BUILDDIST,  
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_BUILDDIST, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( Water_631Gd_DFDR_RED_ANTIHER_GMRES_BUILDDIST, 
      false, true, 
      "resp/serial/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

  // SMP Water 6-31G(d) DFDR REDUCED ANTIHERMETIAN (GMRES)
  CQDFDRTEST_IMPL( PAR_Water_631Gd_DFDR_RED_ANTIHER_GMRES_BUILDDIST, 
       false, true, 
      "resp/parallel/rresp/water_6-31Gd_rhf_dfdr_reduced_antiher_gmres_builddist",
      "water_6-31Gd_rhf_dfdr.bin.ref" )

#endif

