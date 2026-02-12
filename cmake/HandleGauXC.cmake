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

message ( "\n == GauXC ==" )

 
# Method A(default): Pull from github branch
include(FetchContent)

# Propagate ENABLE_MPI flag to GAUXC
if(CQ_ENABLE_MPI)
  set(GAUXC_ENABLE_MPI ON CACHE BOOL "Enable MPI Bindings" FORCE)
else()
  set(GAUXC_ENABLE_MPI OFF CACHE BOOL "Enable MPI Bindings" FORCE)
endif()
#
# Propagate ENABLE_CUDA flag to GAUXC
if(CQ_ENABLE_CUDA)
  set(GAUXC_ENABLE_CUDA ON CACHE BOOL "Enable CUDA Bindings" FORCE)
else()
  set(GAUXC_ENABLE_CUDA OFF CACHE BOOL "Enable CUDA Bindings" FORCE)
endif()
if(CQ_ENABLE_MAGMA)
  set(GAUXC_ENABLE_MAGMA ON CACHE BOOL "Enable MAGMA Bindings" FORCE)
else()
  set(GAUXC_ENABLE_MAGMA OFF CACHE BOOL "Enable MAGMA Bindings" FORCE)
endif()

FetchContent_Declare(
  gauxc
  # Using Aodong's merge_neo branch (master + NEO CPU)
  GIT_REPOSITORY https://github.com/aodongliu/GauXC.git
  GIT_TAG e144e5c73f09f559da1bacab55b3ae8b480adbbb
)

#CMAKE_ARGS "-DGAUXC_ENABLE_MPI=${CQ_ENABLE_MPI} -DGAUXC_ENABLE_CUDA=${CQ_ENABLE_CUDA}"

FetchContent_MakeAvailable( gauxc )

# Link to target
target_link_libraries( cq PUBLIC gauxc::gauxc )




# Method B: Local Gauxc install discovery
#message("Linking a local version of GAUXC")

#include_directories( path/to/GauXC/installation/usr/local/include )
#link_directories( path/to/GauXC/installation/usr/local/lib )
#target_link_libraries( cq PUBLIC -lgauxc -lexchcxx)

message ( " == End GauXC ==\n" )
