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

if(CQ_ENABLE_MPI)

  message( "" )
  
  # FindMPI
  find_package(MPI REQUIRED)
  target_link_libraries( cq PUBLIC MPI::MPI_CXX )
  
  message( "" )
  
  # Print out extraneous information
  message( STATUS "MPIEXEC found to be: ${MPIEXEC}" )
  message( STATUS "MPIEXEC_NUMPROC_FLAG found to be: ${MPIEXEC_NUMPROC_FLAG}" )
  message( STATUS "MPI_INCLUDE_PATH found to be: ${MPI_INCLUDE_PATH}" )

  # Define the source code of the test program
  set(TEST_PROGRAM_SOURCE_CODE "
#include <cstdint>
#include <iostream>
#include <type_traits>

int main() {
    if (std::is_same<long int, int64_t>::value) {
        return 1; // Indicates long int is typedef'd as int64_t
    } else {
        return 0; // Indicates long int is not typedef'd as int64_t
    }
}
")

  # Specify the file path for the test program source
  set(TEST_PROGRAM_SOURCE_FILE "${CMAKE_BINARY_DIR}/check_long_int_typedef.cpp")

  # Write the source code to the file
  file(WRITE ${TEST_PROGRAM_SOURCE_FILE} "${TEST_PROGRAM_SOURCE_CODE}")

  # Try to compile and run the test program
  try_run(RUN_RESULT_VAR COMPILE_RESULT_VAR
          ${CMAKE_BINARY_DIR}
          ${TEST_PROGRAM_SOURCE_FILE}
          COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check if compilation was successful
  if(NOT COMPILE_RESULT_VAR)
    message(FATAL_ERROR "Failed to compile the test program: ${COMPILE_OUTPUT}")
  endif()

  # Interpret the result based on return value
  if("${RUN_RESULT_VAR}" STREQUAL "1")
    set(LONG_INT_IS_INT64_T ON CACHE BOOL "" FORCE)
  endif()

  # Output the result for confirmation
  if(LONG_INT_IS_INT64_T)
    message(STATUS "long int is typedef'd as int64_t.")
  else()
    message(STATUS "long int is not typedef'd as int64_t.")
  endif()
  
  message( "" )
endif()
