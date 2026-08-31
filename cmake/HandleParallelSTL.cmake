#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

message( "\n\n" )
message( "ChronusQ Parallel STL Settings Check:\n" )
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX20_STANDARD_COMPILE_OPTION}")
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
check_cxx_source_compiles("
  #include <algorithm>
  #include <execution>
  #include <vector>
  void f() {
    std::vector<int> v{3,1,2};
    std::sort(std::execution::par, v.begin(), v.end());
  }
" CQ_HAS_PARALLEL_STL)
unset(CMAKE_TRY_COMPILE_TARGET_TYPE)
unset(CMAKE_REQUIRED_FLAGS)

if(CQ_HAS_PARALLEL_STL)
  target_compile_definitions( cq PUBLIC
    CQ_HAS_PARALLEL_STL
    CQ_EXEC_PAR=std::execution::par, )
else()
  message(STATUS "std::execution not available -- parallel STL algorithms will run serially")
  target_compile_definitions( cq PUBLIC CQ_EXEC_PAR= )
endif()
