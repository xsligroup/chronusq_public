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

message( "\n\n" )
message( "ChronusQ Linear Algebra Settings:\n" )

# Eigen3
find_package(Eigen3 CONFIG REQUIRED)
set_property( TARGET cq APPEND PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR}
)

if( NOT CQ_ENABLE_TA )
if( CQ_ENABLE_MPI )

  # If MPI, we need ScaLAPACK
  if(NOT TARGET scalapackpp::scalapackpp )
    FetchContent_Declare( scalapackpp
      GIT_REPOSITORY      https://github.com/wavefunction91/scalapackpp.git
      GIT_TAG             6397f52cf11c0dfd82a79698ee198a2fce515d81
    )
    FetchContent_MakeAvailable( scalapackpp )
    
    # propagate MPI_CXX_SKIP_MPICXX=ON
    if (DEFINED MPI_CXX_COMPILE_DEFINITIONS)
      target_compile_definitions( blacspp     PRIVATE ${MPI_CXX_COMPILE_DEFINITIONS} )
      target_compile_definitions( scalapackpp PRIVATE ${MPI_CXX_COMPILE_DEFINITIONS} )
    endif()
    
    # set {blacspp,scalapackpp}_CONFIG to the install location so that we know where to find it
    set(blacspp_CONFIG ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/blacspp/blacspp-config.cmake)
    set(scalapackpp_CONFIG ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/cmake/scalapackpp/scalapackpp-config.cmake)
  endif()

  target_link_libraries( cq PUBLIC scalapackpp::scalapackpp )

else()

  ## Better BLAS discovery
  FetchContent_Declare( la_cmake
    GIT_REPOSITORY https://github.com/wavefunction91/linalg-cmake-modules.git
    GIT_TAG 290e2e080418e15fc35c80daad23e78c220fd6b4
  )
  FetchContent_MakeAvailable( la_cmake )
  list( APPEND CMAKE_MODULE_PATH ${la_cmake_SOURCE_DIR} )

  find_package( LAPACK REQUIRED )
  if( NOT BLAS_LIBRARIES )
    set(BLAS_LIBRARIES "${LAPACK_LIBRARIES}" CACHE STRING "BLAS Libraries" FORCE)
  endif()

endif()

if( NOT TARGET blaspp )
  
     FetchContent_Declare( blaspp 
       GIT_REPOSITORY https://bitbucket.org/icl/blaspp.git
       GIT_TAG ed392fe
     )
     FetchContent_MakeAvailable( blaspp )

endif()

if( NOT TARGET lapackpp )
  
     FetchContent_Declare( lapackpp 
       GIT_REPOSITORY https://bitbucket.org/icl/lapackpp.git
       GIT_TAG dbcf60f
     )
     FetchContent_MakeAvailable( lapackpp )

endif()

target_link_libraries( cq PUBLIC lapackpp )

endif()

message( "\n\n\n" )
