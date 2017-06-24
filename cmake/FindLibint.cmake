#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2017 Li Research Group (University of Washington)
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

include(ExternalProject)
set( LIBINT2_PREFIX ${PROJECT_SOURCE_DIR}/external/libint2 )
set( LIBINT2_INCLUDEDIR ${LIBINT2_PREFIX}/include 
  ${LIBINT2_PREFIX}/include/libint2 )
set( LIBINT2_LIBDIR ${LIBINT2_PREFIX}/lib )

include_directories("${LIBINT2_INCLUDEDIR}")
link_directories("${LIBINT2_LIBDIR}")

if( NOT EXISTS "${LIBINT2_PREFIX}/include/libint2.hpp" )

  ExternalProject_Add(libint
    PREFIX ${LIBINT2_PREFIX}
    URL "${LIBINT2_PREFIX}/libint-2.3.0-beta.3.tgz"
    CONFIGURE_COMMAND ./configure 
      --prefix=${LIBINT2_PREFIX} 
      CXX=${CMAKE_CXX_COMPILER} 
      CXXFLAGS=${CMAKE_CXX_FLAGS} 
    BUILD_COMMAND make -j2
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )

  list(APPEND CQEX_DEP libint)
  
  message( STATUS "Opting to build a copy of Libint2" )
else()
  message( STATUS "Found Libint2 installation!" )
endif()

list(APPEND CQ_EXT_LINK ${LIBINT2_LIBDIR}/libint2.a)
