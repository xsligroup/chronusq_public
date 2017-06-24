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

include(CheckCXXSourceCompiles)
include(FetchContent)

function ( check_libint_am )

   set ( LIBINT_MAX_AM_SOURCE
   "
   #include <libint2/config.h>

   #if LIBINT_MAX_AM < 5
   #error Libint Angular Momentum must be >= 5
   #endif

   int main() { return 0; }
   " )

   set ( CMAKE_REQUIRED_INCLUDES ${LIBINT2_INCLUDE_DIRS} )
   check_cxx_source_compiles("${LIBINT_MAX_AM_SOURCE}" LIBINT_MAX_AM_GOOD )
   set ( LIBINT_MAX_AM_GOOD PARENT_SCOPE )

endfunction()

message ( "\n == Libint ==" )

#
#  Find preinstalled Libint unless turned off
#

find_package ( Libint2 CONFIG QUIET )


# Prefer CMake installed
if ( TARGET Libint2::cxx )

  check_libint_am()

  if ( LIBINT_MAX_AM_GOOD )
    message(STATUS "found Libint2::cxx target that passed Angular Momentum Test")
  endif()

# Otherwise create dummy target
else()

  FetchContent_Declare (
    Libint2
    #PREFIX ${CUSTOM_LIBINT_PREFIX}
    GIT_REPOSITORY "https://github.com/xsligroup/libint-cq.git"
    GIT_TAG "2.7.0-beta.6"
  )

  FetchContent_MakeAvailable ( Libint2 )

  # ALWAYS build libint using unity build by setting UNITY_BUILD for libint2_obj
  # this saves time and works around issues with Libint library being too big for Ninja on MacOS
  if (NOT DEFINED CMAKE_UNITY_BUILD)
      set_target_properties(libint2_obj PROPERTIES UNITY_BUILD ON)
      message(STATUS "Will unity-build Libint2")
  endif()

  # Libint2::cxx is available in the install tree only, make an alias to the target in the build tree
  if (TARGET libint2_cxx AND NOT TARGET Libint2::cxx)
      add_library(Libint2::cxx ALIAS libint2_cxx)
  endif()

endif()

target_link_libraries( cq PUBLIC Libint2::cxx )

# Workaround to intel compiler bug breaking libint2 ERI evaluation
if( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
  set_property( TARGET cq APPEND PROPERTY
    INTERFACE_COMPILE_DEFINITIONS "LIBINT2_CONSTEXPR_STATICS=0"
  )
endif()

message ( " == End Libint ==\n" )
