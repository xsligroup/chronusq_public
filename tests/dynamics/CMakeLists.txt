#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

# Directories for DYNAMICS tests
set(DYNAMICS_TEST_SOURCE_ROOT "${TEST_ROOT}/dynamics" )
set(DYNAMICS_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/dynamics" )
  
# Set up compilation of DYNAMICS test exe
add_library(dynamicstest_objs OBJECT ../ut.cxx bomd.cxx ehrenfest.cxx bort-ehrenfest.cxx d3.cxx bomd_giao.cxx ehrenfest_giao.cxx)

target_include_directories(dynamicstest_objs PUBLIC ${DYNAMICS_TEST_SOURCE_ROOT}
  ${TEST_BINARY_ROOT})
target_link_libraries(dynamicstest_objs PUBLIC cq)
target_link_libraries(dynamicstest_objs PUBLIC gtest)

add_executable(dynamicstest $<TARGET_OBJECTS:dynamicstest_objs>)
target_link_libraries(dynamicstest PUBLIC dynamicstest_objs)

# Generate directories
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/ehrenfest)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/ehrenfest_giao)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/bomd)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/bomd_giao)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/bort-ehrenfest)

file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/ehrenfest)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/ehrenfest_giao)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/bomd)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/bomd_giao)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/bort-ehrenfest)

include( CQTestGeneration )

add_cq_test( BOMD   dynamicstest "BOMD.*"   )
add_cq_test( EHRENFEST_DYNAMICS   dynamicstest "EHRENFEST_DYNAMICS.*"   )
add_cq_test( BOMD_GIAO   dynamicstest "BOMD_GIAO.*"   )
add_cq_test( EHRENFEST_DYNAMICS_GIAO   dynamicstest "EHRENFEST_DYNAMICS_GIAO.*"   )
add_cq_test( BORT_DYNAMICS   dynamicstest "BORT_DYNAMICS.*"   )

if(CQ_HAS_D3)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/d3)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/d3)
add_cq_test(D3_DYNAMICS  dynamicstest "D3.*")
endif()



if( CQ_ENABLE_MPI )
  add_executable(dynamicstest_mpi $<TARGET_OBJECTS:dynamicstest_objs>)
  target_link_libraries(dynamicstest_mpi PUBLIC dynamicstest_objs)

  add_cq_mpi_test( BOMD_MPI      2   dynamicstest "BOMD.*"   )
  add_cq_mpi_test( EHRENFEST_DYNAMICS_MPI      2   dynamicstest "EHRENFEST_DYNAMICS.*"   )
  add_cq_mpi_test( BORT_DYNAMICS_MPI      2   dynamicstest "BORT_DYNAMICS.*"   )
endif()