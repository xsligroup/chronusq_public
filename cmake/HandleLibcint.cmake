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

message ( "\n == Libcint ==" )

FetchContent_Declare (
  Libcint
  GIT_REPOSITORY "https://github.com/sunqm/libcint"
  GIT_TAG "v5.1.2"
  PATCH_COMMAND git reset --hard && git clean -f -d && git apply "${PROJECT_SOURCE_DIR}/cmake/libcint_breit_sf.patch"
)

FetchContent_MakeAvailable ( Libcint )

FetchContent_GetProperties(Libcint
  SOURCE_DIR Libcint_SOURCE_DIR
  BINARY_DIR Libcint_BINARY_DIR
)

target_include_directories(cint PUBLIC $<BUILD_INTERFACE:${Libcint_BINARY_DIR}/include>)

target_link_libraries( cq PUBLIC cint )

message ( " == End Libcint ==\n" )
