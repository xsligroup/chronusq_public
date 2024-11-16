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

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(ChronusQ
    VERSION 0.8.0
    DESCRIPTION
    "Chronus Quantum is an ab initio computational chemistry \
software package that focuses on explicitly time-dependent and unconventional \
post-SCF quantum mechanical methods."
    LANGUAGES CXX C Fortran)

# Append local cmake directory to find CMAKE Modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

option(CQ_ENABLE_MPI    "Enable MPI parallelism"          OFF)
option(CQ_ENABLE_TA     "Enable TA"  OFF)
option(ENABLE_COVERAGE  "Enable coverage and profiling"   OFF)

# -----------------------------------------------------------------------------
# Build the CQ Libaray
# -----------------------------------------------------------------------------

add_library( cq STATIC )
target_compile_features(cq PUBLIC cxx_std_20)

# Parallelism
message( "\n\n" )
message( "ChronusQ Parallelism Settings" )

# OpenMP
message( "" )
find_package(OpenMP REQUIRED)
target_link_libraries( cq INTERFACE OpenMP::OpenMP_CXX )

# MPI
include(HandleMPI)
message( "\n\n" )

# add includers
target_include_directories( cq PUBLIC 
  "${PROJECT_SOURCE_DIR}/include" 
  "${PROJECT_BINARY_DIR}/include" )

# Non-trivial dependencies
message( "ChronusQ Dependencies" )

# TiledArray
include(HandleTiledArray)

# Linear Algebra
include(HandleLA)

# integral and DFT engine
include(HandleLibint)
include(HandleLibcint)
# GauXC
include(HandleGauXC)
include(FindLibXC)

# HDF5
# HDF5 targets aren't until CMake 3.19, so mock them here
find_package(HDF5 REQUIRED COMPONENTS CXX)
target_link_libraries(cq INTERFACE HDF5::HDF5 )
target_include_directories(cq PRIVATE ${HDF5_INCLUDE_DIRS})

find_package(ZLIB)
target_link_libraries(cq INTERFACE ZLIB::ZLIB)


message( "\n\n" )

# Compiling with coverage report
if(ENABLE_COVERAGE)

  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "ENABLE_COVERAGE Only Compatible with GNU Compilers!")
  else()
    message(STATUS "Enabling Coverage and Profiling -- Warning! May impede performance")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  endif()
endif()

set(BASIS_PATH     "${PROJECT_SOURCE_DIR}/basis"     )
set(AUTHOR_PATH    "${PROJECT_SOURCE_DIR}/AUTHORS"   )

configure_file(
  "${PROJECT_SOURCE_DIR}/include/chronusq_config.hpp.in"
  "${PROJECT_BINARY_DIR}/include/chronusq_config.hpp"
)

# Recurse into subdirectories
add_subdirectory(src)

enable_testing()
add_subdirectory(tests)

set(CQEX_LINK "${CQEX_LINK}" CACHE STRING "List of CQ targets for linking")
