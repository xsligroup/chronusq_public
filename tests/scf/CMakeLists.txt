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

# Directories for SCF tests
set(SCF_TEST_SOURCE_ROOT "${TEST_ROOT}/scf" )
set(SCF_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/scf" )


# Set up compilation of SCF test exe
add_executable(scftest ../ut.cxx
	rhf.cxx uhf.cxx ghf.cxx rohf.cxx x2chf.cxx fourcomphf.cxx
        rhf_giao.cxx uhf_giao.cxx rohf_giao.cxx ghf_giao.cxx x2c_giao.cxx
	rks.cxx uks.cxx x2cks.cxx ks.cxx gks_giao.cxx x2cks_giao.cxx 
        neo_rhf.cxx neo_rks.cxx
        gauxc_ks.cxx
        misc.cxx ri.cxx cediis.cxx qn.cxx
)

target_include_directories(scftest PUBLIC ${SCF_TEST_SOURCE_ROOT}
  ${TEST_BINARY_ROOT})
target_link_libraries(scftest PUBLIC cq)
target_link_libraries(scftest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rohf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ghf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/fourcomp)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/uhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/rohf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ghf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/x2c_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/gks_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/x2cks_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_rohf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/ri_neo_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/cediis)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/qn)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/neo_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/neo_rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/serial/gauxc_ks)

file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rohf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/fourcomp)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/uhf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/rohf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ghf_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/x2c_giao)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_uhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_rohf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_uks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_x2c)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/ri_neo_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/neo_rhf)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/neo_rks)
file(MAKE_DIRECTORY ${SCF_TEST_BINARY_ROOT}/parallel/gauxc_ks)


include( CQTestGeneration )

# Add tests
add_cq_test( RHF_SCF      scftest "RHF.*"      )
add_cq_test( UHF_SCF      scftest "UHF.*"      )
add_cq_test( GHF_SCF      scftest "GHF.*"      )
add_cq_test( ROHF_SCF     scftest "ROHF.*"     )
add_cq_test( X2CHF_SCF    scftest "X2CHF.*"    )
add_cq_test( FOURCHF_SCF  scftest "FOURCHF.*"    )
add_cq_test( RHF_GIAO_SCF scftest "RHF_GIAO.*" )
add_cq_test( UHF_GIAO_SCF scftest "UHF_GIAO.*" )
add_cq_test( ROHF_GIAO_SCF scftest "ROHF_GIAO.*")
add_cq_test( GHF_GIAO_SCF scftest "GHF_GIAO.*" )
add_cq_test( GKS_GIAO_SCF scftest "GKS_GIAO.*" )
add_cq_test( X2C_GIAO_SCF scftest "X2C_GIAO.*" )
add_cq_test( X2CKS_GIAO_SCF scftest "X2CKS_GIAO.*" )
add_cq_test( RI_RHF_SCF   scftest "RI_RHF.*"   )
add_cq_test( RI_UHF_SCF   scftest "RI_UHF.*"   )
add_cq_test( RI_ROHF_SCF  scftest "RI_ROHF.*"  )
add_cq_test( RI_X2CHF_SCF scftest "RI_X2CHF.*" )
add_cq_test( RI_NEO_RHF_SCF scftest "RI_NEO_RHF.*" )

add_cq_test( RKS_SCF      scftest "RKS.*"        )
add_cq_test( UKS_SCF      scftest "UKS.*"        )
add_cq_test( X2CKS_SCF    scftest "X2CKS.*"      )
add_cq_test( KS_KEYWORD   scftest "KS_KEYWORD.*" )
add_cq_test( KS_FUNC      scftest "KS_FUNC.*"    )
add_cq_test( RI_RKS_SCF   scftest "RI_RKS.*"     )
add_cq_test( RI_UKS_SCF   scftest "RI_UKS.*"     )
add_cq_test( RI_X2CKS_SCF scftest "RI_X2CKS.*"   )
#add_cq_test( NEO_RKS_SCF  scftest "NEO_RKS.*"  )
add_cq_test( NEO_RHF_SCF  scftest "NEO_RHF.*"  )
add_cq_test( GAUXC_KS_SCF  scftest "GAUXC_KS.*"  )
#add_cq_test( CEDIIS_SCF   scftest "CEDIIS.*"  )
#add_cq_test( QN_SCF       scftest "QN_SCF*"  )
#add_cq_test( FOURC_QN_SCF scftest "FOURC_QN_SCF*"  )

add_cq_test( MISC_SCF     scftest "MISC_SCF.*" )


if( CQ_ENABLE_MPI )
  add_cq_mpi_test( RHF_SCF_MPI      2   scftest "RHF.*"   )
  add_cq_mpi_test( UHF_SCF_MPI      2   scftest "UHF.*"   )
  add_cq_mpi_test( ROHF_SCF_MPI      2   scftest "ROHF.*" )
  add_cq_mpi_test( X2CHF_SCF_MPI    2   scftest "X2CHF.*" )
  add_cq_mpi_test( RHF_GIAO_SCF_MPI 2   scftest "RHF_GIAO.*"   )
  #add_cq_mpi_test( UHF_GIAO_SCF_MPI 2   scftest "UHF_GIAO.*"   )
  #add_cq_mpi_test( ROHF_GIAO_SCF_MPI 2  scftest "ROHF_GIAO.*"  )
  add_cq_mpi_test( X2C_GIAO_SCF_MPI 2  scftest "X2C_GIAO.*"  )
  add_cq_mpi_test( RI_RHF_SCF_MPI   2   scftest "RI_RHF.*"   )
  add_cq_mpi_test( RI_UHF_SCF_MPI   2   scftest "RI_UHF.*"   )
  add_cq_mpi_test( RI_ROHF_SCF_MPI  2   scftest "RI_ROHF.*"  )
  add_cq_mpi_test( RI_X2CHF_SCF_MPI 2   scftest "RI_X2CHF.*" )

  add_cq_mpi_test( RKS_SCF_MPI 2   scftest "RKS.*"   )
  add_cq_mpi_test( UKS_SCF_MPI 2   scftest "UKS.*"   )
  #add_cq_mpi_test( X2CKS_SCF_MPI 2 scftest "X2CKS.*" )
  add_cq_mpi_test( RI_RKS_SCF_MPI   2   scftest "RI_RKS.*"   )
  add_cq_mpi_test( RI_UKS_SCF_MPI   2   scftest "RI_UKS.*"   )
  add_cq_mpi_test( RI_X2CKS_SCF_MPI 2   scftest "RI_X2CKS.*" )
  #add_cq_mpi_test( QN_SCF_MPI       2   scftest "QN_SCF*"  )
  #add_cq_mpi_test( NEO_RKS_MPI          2    scftest "NEO_RKS.*"  )
  add_cq_mpi_test( NEO_RHF_MPI          2    scftest "NEO_RHF.*"  )
  add_cq_mpi_test( RI_NEO_RHF_SCF_MPI   2   scftest "RI_NEO_RHF.*"   )
  add_cq_mpi_test( GAUXC_KS_SCF_MPI     2   scftest "GAUXC_KS.*"  )
endif()
