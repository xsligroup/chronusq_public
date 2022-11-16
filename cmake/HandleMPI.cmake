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

include(FetchContent)

if(CQ_ENABLE_MPI)

  message( "" )
  
  # FindMPI
  find_package(MPI REQUIRED)
  target_link_libraries( cq PUBLIC MPI::MPI_CXX )
  
  message( "" )
  
  # Print out extraneous information
  message( STATUS "MPIEXEC found to be: ${MPIEXEC}" )
  message( STATUS "MPIEXEC_NUMPROC_FLAG found to be: ${MPIEXEC_NUMPROC_FLAG}" )
  message( STATUS "MPI_INCLUDE_PATH found to be: ${MPI_INCLUDE_PATH}" )
  
  message( "" )
  
  # MXX
  message( STATUS "Adding CMake Target for MXX" )
  FetchContent_Declare(
    mxx
    GIT_REPOSITORY https://github.com/patflick/mxx.git 
    GIT_TAG e1f4acd8f5dc91da4945b5f6b6e7828991afeb0a
    PATCH_COMMAND git apply "${PROJECT_SOURCE_DIR}/cmake/mxx_complex_datatype.patch"
  )  

  FetchContent_MakeAvailable ( mxx )
  target_include_directories(cq PUBLIC "${mxx_SOURCE_DIR}/include")
  
  message( "" )
endif()
