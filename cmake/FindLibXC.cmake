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

message ( "\n == LibXC ==" )


## Find LibXC
include(FetchContent)

# Check if libxc has already been fetched and is available
FetchContent_GetProperties(libxc)
# Check if libxc content has already been populated
if(libxc_POPULATED) 

  message(STATUS "Libxc has already been populated")

else()

  FetchContent_Declare(
    libxc
    GIT_REPOSITORY https://gitlab.com/libxc/libxc.git
    GIT_TAG 6.2.0  # v6.2.0
  )
  FetchContent_MakeAvailable(libxc)

endif()

if(NOT TARGET xc) # Check if the target already exists to avoid redefinition
  add_library(Libxc::xc ALIAS xc)
endif()

# Configure target properties only if they haven't been set before
if(TARGET xc)
  target_include_directories(xc
    PUBLIC 
      $<BUILD_INTERFACE:${libxc_SOURCE_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/src>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}>
      $<BUILD_INTERFACE:${libxc_BINARY_DIR}/gen_funcidx>
  )

  # disable unity builds for libxc
  if (CMAKE_UNITY_BUILD)
    set_target_properties(xc PROPERTIES UNITY_BUILD OFF)
    message(STATUS "Will disable unity-build for Libxc::xc")
  endif()
endif()

set( BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" FORCE )

target_link_libraries( cq PUBLIC Libxc::xc )


message ( " == End LibXC ==\n" )