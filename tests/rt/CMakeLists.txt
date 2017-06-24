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
set(RT_TEST_SOURCE_ROOT "${TEST_ROOT}/rt" )
set(RT_TEST_BINARY_ROOT "${TEST_BINARY_ROOT}/rt" )
  
# Set up compilation of RT test exe
add_executable(rttest ../ut.cxx rrt.cxx urt.cxx grt.cxx neo_rrt.cxx ri_neo_rrt.cxx)

target_include_directories(rttest PUBLIC ${RT_TEST_SOURCE_ROOT} 
  ${TEST_BINARY_ROOT})
target_link_libraries(rttest PUBLIC cq)
target_link_libraries(rttest PUBLIC gtest)

# Generate directories
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/urt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/grt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/neo_rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/serial/ri_neo_rrt)

file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/urt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/grt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/neo_rrt)
file(MAKE_DIRECTORY ${RT_TEST_BINARY_ROOT}/parallel/rineo_rrt)

include( CQTestGeneration )

add_cq_test( RHF_RT   rttest "RHF_RT.*"   )
add_cq_test( UHF_RT   rttest "UHF_RT.*"   )
add_cq_test( X2CHF_RT rttest "X2CHF_RT.*" )
#
add_cq_test( RKS_RT   rttest "RKS_RT.*"   )
add_cq_test( UKS_RT   rttest "UKS_RT.*"   )
add_cq_test( X2CKS_RT rttest "X2CKS_RT.*" )
#
add_cq_test( NEO_RHF_RT   rttest "NEO_RHF_RT.*"   )
##add_cq_test( NEO_UHF_RT   rttest "NEO_UHF_RT.*"   )
#add_cq_test( RI_NEO_RHF_RT   rttest "RI_NEO_RHF_RT.*"   )
#
add_cq_test( RESTART_RT rttest "RESTART_RT.*" )
