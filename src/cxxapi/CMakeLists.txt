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

set(INPUT_SRC input/parse.cxx input/freeparsers.cxx)

set(OPT_SRC input/molopts.cxx input/basisopts.cxx 
    input/singleslateropts.cxx input/scfopts.cxx input/rtopts.cxx
    input/intsopts.cxx input/miscopts.cxx input/respopts.cxx 
    input/coupledclusteropts.cxx input/geomopts.cxx input/mcscfopts.cxx
    input/cubeopts.cxx input/scfcubeopts.cxx
    input/saverefs.cxx input/perturbopts.cxx
    input/gauxcopts.cxx procedural.cxx output.cxx)

target_sources(cq PRIVATE ${INPUT_SRC} ${OPT_SRC})

add_executable(chronusq chronusq.cxx)
target_link_libraries(chronusq PUBLIC cq)

add_custom_target(linkexe ALL
  ${CMAKE_COMMAND} -E create_symlink ${PROJECT_BINARY_DIR}/src/cxxapi/chronusq ${PROJECT_BINARY_DIR}/chronusq )

install(TARGETS chronusq RUNTIME DESTINATION "bin")
