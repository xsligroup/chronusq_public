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

include(FetchContent)

FetchContent_Declare(
  libxc
  GIT_REPOSITORY https://gitlab.com/eduard1/libxc.git
  GIT_TAG        v5.0.0-plus-prs-324-351
)
set( Libxc_VERSION 5.0.0 )

set( OLD_BUILD_TESTING ${BUILD_TESTING} )
set( BUILD_TESTING OFF CACHE BOOL "" FORCE )

FetchContent_MakeAvailable( libxc )
if( TARGET xc ) 
  message( "XC IS A TARGET!" )
endif()
add_library( Libxc::xc ALIAS xc )
target_include_directories( xc 
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

set( BUILD_TESTING ${OLD_BUILD_TESTING} CACHE BOOL "" FORCE )

target_link_libraries( cq PUBLIC Libxc::xc )
