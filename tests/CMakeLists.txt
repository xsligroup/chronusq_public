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

# Define global test variables
set( TEST_ROOT ${PROJECT_SOURCE_DIR}/tests )
set( TEST_BINARY_ROOT ${PROJECT_BINARY_DIR}/tests )

# Common UT configuration header
configure_file(
  "${TEST_ROOT}/ut.hpp.in"
  "${TEST_BINARY_ROOT}/ut.hpp"
)

#find_package( GTest REQUIRED )
include( HandleGTest )

add_subdirectory(scf)
add_subdirectory(rt)
add_subdirectory(resp)
add_subdirectory(func)
add_subdirectory(mcscf)
add_subdirectory(coupledcluster)
add_subdirectory(pt)
add_subdirectory(dynamics)