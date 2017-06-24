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

# Directories for DYNAMICS tests
set(DYNAMICS_TEST_SOURCE_ROOT "${TEST_ROOT}/dynamics" )
set(DYNAMICS_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/dynamics" )
  
# Set up compilation of DYNAMICS test exe
add_executable(dynamicstest ../ut.cxx ehrenfest.cxx)

target_include_directories(dynamicstest PUBLIC ${DYNAMICS_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(dynamicstest PUBLIC cq)
target_link_libraries(dynamicstest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/ehrenfest)
file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/serial/bomd)

#file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/ehrenfest)
#file(MAKE_DIRECTORY ${DYNAMICS_TEST_BINARY_ROOT}/parallel/bomd)

include( CQTestGeneration )

add_cq_test( EHRENFEST_DYNAMICS   dynamicstest "EHRENFEST_DYNAMICS.*"   )
#add_cq_test( BOMD_DYNAMICS   dynamicstest "BOMD_DYNAMICS.*"   )
