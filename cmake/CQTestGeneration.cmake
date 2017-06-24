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

# Add a regular test
function( add_cq_test _test_name _test_exe _filter )

  add_test( NAME ${_test_name} COMMAND ${_test_exe} --gtest_filter=${_filter} )

endfunction()


# Add an MPI test
function( add_cq_mpi_test _test_name _np _test_exe _filter )

  add_test( NAME ${_test_name} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_np}
    ${MPIEXEC_PREFLAGS} "./${_test_exe}" ${MPIEXEC_POSTFLAGS} "--gtest_filter=${_filter}" )

endfunction()
