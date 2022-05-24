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

include(CheckCSourceCompiles)
include(FetchContent)


# Checks if quadmath can be used and sets the result to QUADMATH_FOUND
function ( check_libcint_quadmath )

  set ( QUADMATH_TEST_SOURCE
  "
  #include <quadmath.h>
  int main() { fabsq(1); return 0; }
  ")
  
  set ( CMAKE_REQUIRED_LIBRARIES quadmath )
  check_c_source_compiles ( "${QUADMATH_TEST_SOURCE}" QUADMATH_FOUND )

  set ( QUADMATH_FOUND ${LIBCINT_QUADMATH_REQUIRED} PARENT_SCOPE )

endfunction()


message ( "\n == Libcint ==" )

#
#  Find preinstalled Libcint unless turned off
#

if ( NOT CQ_BUILD_LIBCINT_TYPE STREQUAL "FORCE" )

  # Otherwise create dummy target
  if ( DEFINED Libcint_ROOT )

    find_library ( LIBCINT_LIBRARIES
                   NAMES cint
                   PATHS "${Libcint_ROOT}/lib64"
                 )

    if ( LIBCINT_LIBRARIES AND
         EXISTS "${Libcint_ROOT}/include/cint.h" )
      
      set ( LIBCINT_INCLUDE_DIRS
          ${Libcint_ROOT}/include
      )

      set ( libcint_location ${LIBCINT_LIBRARIES} )

      check_libcint_quadmath()
      if ( QUADMATH_FOUND )
        list( APPEND LIBCINT_LIBRARIES quadmath )
      endif()

      add_library ( ChronusQ::Libcint INTERFACE IMPORTED )
      set_target_properties ( ChronusQ::Libcint PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LIBCINT_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES      "${LIBCINT_LIBRARIES}"
      )

      set ( LIBCINT_FOUND TRUE )

    endif()

  endif()

  if ( LIBCINT_FOUND )
    message ( STATUS "Found External Libcint Installation: ${libcint_location}" )
  endif()

endif()


#
#  Build Libcint if a suitable libcint hasn't been found
#

if ( NOT TARGET ChronusQ::Libcint )

  if ( NOT CQ_BUILD_LIBCINT_TYPE STREQUAL "NONE" )

    if ( NOT CQ_BUILD_LIBCINT_TYPE STREQUAL "FORCE" )

      message( STATUS "Checking for a previously built Libcint for CQ" )
 
      find_library ( LIBCINT_LIBRARIES
                     NAMES cint
                     PATHS "${FETCHCONTENT_BASE_DIR}/libcint-build"
                     NO_DEFAULT_PATH )
 
      if ( LIBCINT_LIBRARIES AND
           EXISTS "${FETCHCONTENT_BASE_DIR}/libcint-build/include/cint.h" AND
           EXISTS "${FETCHCONTENT_BASE_DIR}/libcint-src/include/cint_funcs.h" )
 
 
        set ( LIBCINT_INCLUDE_DIRS
              "${FETCHCONTENT_BASE_DIR}/libcint-build/include"
              "${FETCHCONTENT_BASE_DIR}/libcint-build/src"
              "${FETCHCONTENT_BASE_DIR}/libcint-src/include" 
              "${FETCHCONTENT_BASE_DIR}/libcint-src/src"
        )


        set(libcint_location ${LIBCINT_LIBRARIES})

        check_libcint_quadmath()
        if ( QUADMATH_FOUND )
          list( APPEND LIBCINT_LIBRARIES quadmath )
        endif()
 
        add_library ( ChronusQ::Libcint INTERFACE IMPORTED )
        set_target_properties ( ChronusQ::Libcint PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${LIBCINT_INCLUDE_DIRS}"
          INTERFACE_LINK_LIBRARIES      "${LIBCINT_LIBRARIES}"
        )
        
        message ( STATUS "Found CQ Libcint Installation: ${libcint_location}" )
      
      endif()

    endif()


    # If we're forcing building or failed finding a previously built version of
    #   libcint for CQ
    if ( NOT TARGET ChronusQ::Libcint )

      message( STATUS "Opting to build a copy of Libcint" )

      # Update project policy for Libcint
      set ( CMAKE_POLICY_DEFAULT_CMP0048 NEW )

      FetchContent_Declare (
        Libcint
        GIT_REPOSITORY "https://github.com/sunqm/libcint"
        GIT_TAG "v5.1.2"
        UPDATE_COMMAND cd ${FETCHCONTENT_BASE_DIR}/libcint-src && patch -N < ${PROJECT_SOURCE_DIR}/external/libcint/patch/CMakeLists.txt.patch || patch -N < ${PROJECT_SOURCE_DIR}/external/libcint/patch/CMakeLists.txt.patch | grep "Skipping patch" -q
      )
  
      FetchContent_GetProperties ( Libcint )

      if ( NOT Libcint_POPULATED )

        message ( STATUS "Downloading Libcint..." )
        FetchContent_Populate ( Libcint )
        message ( STATUS "Downloading Libcint - Done" )

        add_subdirectory ( ${libcint_SOURCE_DIR} ${libcint_BINARY_DIR} )

      endif()

      add_library ( ChronusQ::Libcint ALIAS cint )

    endif()
  
  else()
  
    message ( FATAL_ERROR "Suitable Libcint installation could not be found! \
Set Libcint_ROOT to the prefix of the Libcint installation or set \
CQ_BUILD_LIBCINT_TYPE to ALLOW or FORCE."
    )
  
  endif()

endif()


target_link_libraries( ChronusQ::Dependencies INTERFACE ChronusQ::Libcint )
copy_header_properties( ChronusQ::Libcint ChronusQ::DepHeaders )
list(APPEND CQ_EXT_LINK ChronusQ::Libcint)

message ( " == End Libcint ==\n" )

