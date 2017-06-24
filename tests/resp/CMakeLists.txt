#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2022 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
#

# Directories for RESP tests
set(RESP_TEST_SOURCE_ROOT "${TEST_ROOT}/resp" )
set(RESP_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/resp" )

# Set up compilation of RESP test exe
add_executable(resptest ../ut.cxx rres.cxx rfdr.cxx rmor.cxx
  rres_lda.cxx rres_gga.cxx rres_hgga.cxx rpp.cxx upp.cxx gpp.cxx
  ures_hf.cxx ures_lsda.cxx ures_gga.cxx ures_hgga.cxx 
  gres_hf.cxx gres_lsda.cxx gres_gga.cxx gres_hgga.cxx 
  x2cres_hf.cxx x2cres_lsda.cxx x2cres_gga.cxx x2cres_hgga.cxx
  neohf_res.cxx rres_giao.cxx x2cres_giao.cxx 
  misc.cxx)

target_include_directories(resptest PUBLIC ${RESP_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(resptest PUBLIC cq)
target_link_libraries(resptest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rresp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/uresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/uresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/gresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/gresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/x2cresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/x2cresp_giao_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/x2cresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/neoresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rresp_giao)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rresp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/uresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/uresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/gresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/gresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/x2cresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/x2cresp_giao_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/x2cresp_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/neoresp_hf)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rresp_giao)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rmor)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rmor_ks)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rmor)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rmor_ks)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/rpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/upp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/gpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/rpp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/upp)
file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/parallel/gpp)

file(MAKE_DIRECTORY ${RESP_TEST_BINARY_ROOT}/serial/misc)

# Add the Tests

include( CQTestGeneration )

add_cq_test( RHF_RESIDUE       resptest "RHF_RESIDUE.*"      )
add_cq_test( RKS_RESIDUE       resptest "RKS_RESIDUE.*"      )
add_cq_test( UHF_RESIDUE       resptest "UHF_RESIDUE.*"      )
add_cq_test( ULSDA_RESIDUE     resptest "ULSDA_RESIDUE.*"    )
add_cq_test( UGGA_RESIDUE      resptest "UGGA_RESIDUE.*"     )
add_cq_test( UHGGA_RESIDUE     resptest "UHGGA_RESIDUE.*"    )
add_cq_test( GHF_RESIDUE       resptest "GHF_RESIDUE.*"      )
add_cq_test( GLSDA_RESIDUE     resptest "GLSDA_RESIDUE.*"    )
add_cq_test( GGGA_RESIDUE      resptest "GGGA_RESIDUE.*"     )
add_cq_test( GHGGA_RESIDUE     resptest "GHGGA_RESIDUE.*"    )
add_cq_test( X2CHF_RESIDUE     resptest "X2CHF_RESIDUE.*"    )
add_cq_test( X2CLSDA_RESIDUE   resptest "X2CLSDA_RESIDUE.*"  )
add_cq_test( X2CGGA_RESIDUE    resptest "X2CGGA_RESIDUE.*"   )
add_cq_test( X2CHGGA_RESIDUE   resptest "X2CHGGA_RESIDUE.*"  )
add_cq_test( RHF_PP_RESIDUE    resptest "RHF_PP_RESIDUE.*"   )
add_cq_test( UHF_PP_RESIDUE    resptest "UHF_PP_RESIDUE.*"   )
add_cq_test( RHF_FDR           resptest "RHF_FDR.*"          )
add_cq_test( RHF_DFDR          resptest "RHF_DFDR.*"         )
add_cq_test( RHF_MOR           resptest "RHF_MOR.*"          )
add_cq_test( MISC_RESP         resptest "MISC_RESP.*"        )
#add_cq_test( NEOHF_RESIDUE     resptest "NEOHF_RESIDUE.*"    )
add_cq_test( RHF_GIAO_RESIDUE  resptest "RHF_GIAO_RESIDUE.*" ) 
add_cq_test( X2CHF_GIAO_RESIDUE     resptest "X2CHF_GIAO_RESIDUE.*"    )

if( CQ_ENABLE_MPI )

  #add_cq_mpi_test( RHF_RESIDUE_MPI  2  resptest "RHF_RESIDUE.*" )
  #add_cq_mpi_test( RKS_RESIDUE_MPI  2  resptest "RKS_RESIDUE.*" )
  #add_cq_mpi_test( GHGGA_RESIDUE_MPI  2  resptest "GHGGA_RESIDUE_MPI.*" )

  #add_cq_mpi_test( RHF_FDR_MPI  2  resptest "RHF_FDR.*"  )
  #add_cq_mpi_test( RHF_DFDR_MPI 2  resptest "RHF_DFDR.*" )
  #add_cq_mpi_test( RHF_MOR_MPI  2  resptest "RHF_MOR.*"  )

endif()

#if( CQ_ENABLE_MPI )
#
#  add_test( NAME RHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_RESIDUE" )
# #add_test( NAME UHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_RESIDUE" )
# #add_test( NAME GHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_RESIDUE" )
# #add_test( NAME X2CHF_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_RESIDUE" )
#  add_test( NAME RKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_RESIDUE" )
# #add_test( NAME UKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_RESIDUE" )
# #add_test( NAME GKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_RESIDUE" )
# #add_test( NAME X2CKS_RESIDUE_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_RESIDUE" )
#
#
#  add_test( NAME RHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_FDR" )
# #add_test( NAME UHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_FDR" )
# #add_test( NAME GHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_FDR" )
# #add_test( NAME X2CHF_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_FDR" )
# #add_test( NAME RKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_FDR" )
# #add_test( NAME UKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_FDR" )
# #add_test( NAME GKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_FDR" )
# #add_test( NAME X2CKS_FDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_FDR" )
#
#  add_test( NAME RHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_RESP/RHF_DFDR" )
# #add_test( NAME UHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_RESP/UHF_DFDR" )
# #add_test( NAME GHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_RESP/GHF_DFDR" )
# #add_test( NAME X2CHF_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_RESP/X2CHF_DFDR" )
# #add_test( NAME RKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RKS_RESP/RKS_DFDR" )
# #add_test( NAME UKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UKS_RESP/UKS_DFDR" )
# #add_test( NAME GKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GKS_RESP/GKS_DFDR" )
# #add_test( NAME X2CKS_DFDR_MPI COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CKS_RESP/X2CKS_DFDR" )
#
#
#
#  add_test( NAME RHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=RHF_MOR"  )
# #add_test( NAME UHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=UHF_MOR"  )
# #add_test( NAME GHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=GHF_MOR"  )
# #add_test( NAME X2CHF_MOR_MPI      COMMAND ${MPIEXEC} -np 2 "./resptest" "--report_level=detailed" "--run_test=X2CHF_MOR"  )
#
#endif()
