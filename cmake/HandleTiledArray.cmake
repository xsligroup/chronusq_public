# 
#  This file is part of the Chronus Quantum (ChronusQ) software package
#  
#  Copyright (C) 2014-2022 Li Research Group (University of Washington)
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#  
#  Contact the Developers:
#    E-Mail: xsli@uw.edu
#  
#

include( FetchContent )

# Do nothing if TA is disabled
if( CQ_ENABLE_TA )

  if( NOT CQ_ENABLE_MPI )
    message(FATAL_ERROR "TA requires MPI. Please rerun with CQ_ENABLE_MPI=On" )
  endif()

  set( ENABLE_SCALAPACK ON CACHE BOOL "Enable ScaLAPACK" FORCE )
  FetchContent_Declare( tiledarray
    GIT_REPOSITORY https://github.com/ValeevGroup/tiledarray.git
    GIT_TAG d72357a06384c673f82ba12d8a2c94746a0ca379
  )
  
  FetchContent_MakeAvailable( tiledarray )
  target_link_libraries( cq PUBLIC tiledarray)
  set(CQ_HAS_TA ON CACHE BOOL "" FORCE)
      
else( CQ_ENABLE_TA )
  message( STATUS "TiledArray not enabled; Coupled-cluster functionality in CQ is disabled." )
endif( CQ_ENABLE_TA )
