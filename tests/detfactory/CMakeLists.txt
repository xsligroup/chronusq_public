#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2020 Li Research Group (University of Washington)
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
set(DETFACTORY_TEST_SOURCE_ROOT "${TEST_ROOT}/detfactory" )
set(DETFACTORY_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/detfactory" )


# Set up compilation of Functionality test exe
add_executable(detfactorytest ../ut.cxx detfactory.cxx)

target_include_directories(detfactorytest PUBLIC ${DETFACTORY_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(detfactorytest PUBLIC cq)
target_link_libraries(detfactorytest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${DETFACTORY_TEST_BINARY_ROOT}/serial/cas)
file(MAKE_DIRECTORY ${DETFACTORY_TEST_BINARY_ROOT}/serial/ras)

file(MAKE_DIRECTORY ${DETFACTORY_TEST_BINARY_ROOT}/parallel/cas)
file(MAKE_DIRECTORY ${DETFACTORY_TEST_BINARY_ROOT}/parallel/ras)

if(CQEX_DEP)
  add_dependencies(detfactorytest ${CQEX_DEP})
endif()



include( CQTestGeneration )


add_cq_test(DETFACTORY detfactorytest "DETFACTORY.*")



