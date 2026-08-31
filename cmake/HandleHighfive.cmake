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

message("\n == HighFive ==")

include(FetchContent)

# If GauXC or something else already added HighFive, just reuse it
if(TARGET HighFive OR TARGET HighFive::HighFive)
  message(STATUS "HighFive target already available (reusing)")
else()
  FetchContent_GetProperties(HighFive)
  if(NOT HighFive_POPULATED)
    # Only set options if *we* are going to fetch/build HighFive
    set(HIGHFIVE_UNIT_TESTS        OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_EXAMPLES          OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_BUILD_DOCS        OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_SPAN         OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_BOOST        OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_BOOST_SPAN   OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_EIGEN        OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_OPENCV       OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_XTENSOR      OFF CACHE BOOL "" FORCE)
    set(HIGHFIVE_TEST_HALF_FLOAT   OFF CACHE BOOL "" FORCE)

    # If you want parallel HDF5 when CQ uses MPI:
    if(CQ_ENABLE_MPI)
      set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" FORCE)
    endif()

    FetchContent_Declare(
      HighFive
      GIT_REPOSITORY https://github.com/highfive-devs/highfive.git
      GIT_TAG v3.0.0-beta2
    )
    FetchContent_MakeAvailable(HighFive)
  endif()
endif()

if(TARGET HighFive AND NOT TARGET HighFive::HighFive)
  add_library(HighFive::HighFive ALIAS HighFive)
elseif(TARGET HighFive::HighFive AND NOT TARGET HighFive)
  add_library(HighFive ALIAS HighFive::HighFive)
endif()

target_link_libraries(cq PUBLIC HighFive::HighFive)

message(" == End HighFive ==\n")

