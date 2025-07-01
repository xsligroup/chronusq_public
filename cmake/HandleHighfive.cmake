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

message ( "\n == Highfive ==" )

 
# Method A(default): Pull from github branch
include(FetchContent)
# HDF5
set(HIGHFIVE_UNIT_TESTS OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_EXAMPLES OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_BUILD_DOCS OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_SPAN OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_BOOST OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_BOOST_SPAN OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_EIGEN OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_OPENCV OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_XTENSOR OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_TEST_HALF_FLOAT OFF CACHE INTERNAL "")  # Forces the value
set(HIGHFIVE_HAS_CONCEPTS ON CACHE INTERNAL "")  # Forces the value
set(HDF5_USE_STATIC_LIBRARIES ON CACHE INTERNAL "")  # Forces the value
set(HDF5_PREFER_PARALLEL ON CACHE INTERNAL "")  # Forces the value

FetchContent_Declare(
  HighFive
  GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
  GIT_TAG v3.0.0-beta2)
FetchContent_MakeAvailable(HighFive)

add_library(HighFive::HighFive ALIAS HighFive) # Highfive doesn't export the namespaced target?

target_link_libraries(cq PUBLIC HighFive::HighFive)

# Method B: Local Highfive install discovery (NOT SUPPORTED)
#message("Linking a local version of Highfive")

message ( " == End HighFive ==\n" )
