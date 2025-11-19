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

# Do nothing if D3 is disabled
  if( CQ_ENABLE_D3 )
  message("\n == Simple-DFTD3 ==")
  
  FetchContent_Declare(
    simple_dftd3
    GIT_REPOSITORY https://github.com/dftd3/simple-dftd3.git
    GIT_TAG        v1.2.1
  )
  
  FetchContent_MakeAvailable(simple_dftd3)
  
  if(TARGET s-dftd3)
    add_library(SimpleDFTD3::dftd3 ALIAS s-dftd3)
  elseif(TARGET dftd3)
    add_library(SimpleDFTD3::dftd3 ALIAS dftd3)
  else()
    message(FATAL_ERROR "simple-dftd3: expected target 's-dftd3' (or 'dftd3') not found.")
  endif()
  
  target_link_libraries(cq PUBLIC SimpleDFTD3::dftd3)
  
  set(CQ_HAS_D3 ON CACHE BOOL "" FORCE)
  
  message(" == End Simple-DFTD3 ==\n")

else()
  message(STATUS "simple-dftd3 not enabled; dispersion corrections (D3) are disabled.")
endif()


