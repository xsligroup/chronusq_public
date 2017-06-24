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

# Directories for basic functionality tests
set(FUNC_TEST_SOURCE_ROOT "${TEST_ROOT}/func" )
set(FUNC_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/func" )

# Set up compilation of Functionality test exe
add_executable(functest ../ut.cxx contract.cxx ordqz.cxx gplhr.cxx davidson.cxx)

target_compile_definitions(functest PUBLIC CQ_FUNC_TEST)
target_include_directories(functest PUBLIC ${FUNC_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(functest PUBLIC cq)
target_link_libraries(functest PUBLIC gtest)

include( CQTestGeneration )

add_cq_test( DIRECT_CONTRACTION functest "DIRECT_CONTRACTION*" )
add_cq_test( ORDQZ              functest "ORDQZ.*" )
add_cq_test( GPLHR              functest "GPLHR.*" )
add_cq_test( DAVIDSON           functest "DAVIDSON.*" )
add_cq_test( CQMEMMANAGER       functest "CQMEM.*" )

if( CQ_ENABLE_MPI )
  add_cq_mpi_test( GPLHR_MPI 2 functest "GPLHR.*" )
  add_cq_mpi_test( DAVIDSON_MPI 2 functest "DAVIDSON.*" )
endif()


