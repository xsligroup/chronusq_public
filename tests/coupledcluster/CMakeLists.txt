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
if(CQ_HAS_TA)
# Directories for Coupled-Cluster tests
set(CC_TEST_SOURCE_ROOT "${TEST_ROOT}/coupledcluster" )
set(CC_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/coupledcluster" )


# Set up compilation of Coupled-Cluster test exe
add_executable(cctest ../ut.cxx ccsd.cxx eom-ccsd.cxx)

target_include_directories(cctest PUBLIC ${CC_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(cctest PUBLIC cq)
target_link_libraries(cctest PUBLIC gtest)

file(MAKE_DIRECTORY ${CC_TEST_BINARY_ROOT}/serial/ccsd)
file(MAKE_DIRECTORY ${CC_TEST_BINARY_ROOT}/serial/eom-ccsd)
file(MAKE_DIRECTORY ${CC_TEST_BINARY_ROOT}/parallel/ccsd)
file(MAKE_DIRECTORY ${CC_TEST_BINARY_ROOT}/parallel/eom-ccsd)
include( CQTestGeneration )

add_cq_test( CCSD cctest "CCSD.*" )
add_cq_test( EOM_CCSD cctest "EOM_CCSD.*" )

  if( CQ_ENABLE_MPI )
#    add_cq_mpi_test( CCSD_MPI          2   cctest "CCSD.*" )
#    add_cq_mpi_test( EOM_CCSD_MPI      2   cctest "EOM_CCSD.*" )
  endif()
endif()

