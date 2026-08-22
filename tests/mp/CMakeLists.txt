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

# Directories for RT tests
set(MP_TEST_SOURCE_ROOT "${TEST_ROOT}/mp" )
set(MP_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/mp" )
  

# Set up compilation of MCSCF test exe
add_executable(mptest ../ut.cxx mp2.cxx)

target_include_directories(mptest PUBLIC ${MP_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(mptest PUBLIC cq)
target_link_libraries(mptest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${MP_TEST_BINARY_ROOT}/serial/mp2)

file(MAKE_DIRECTORY ${MP_TEST_BINARY_ROOT}/parallel/mp2)

include( CQTestGeneration )

# MP2 test
add_cq_test(MP2 mptest "MP2.*")

if( CQ_ENABLE_MPI )

endif()

