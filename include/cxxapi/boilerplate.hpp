/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#pragma once


#include <chronusq_sys.hpp>
#include <libint2/cxxapi.h>
#include <basisset/basisset_def.hpp>
#include <util/threads.hpp>
#include <H5Cpp.h>

#include <util/mpi.hpp>

#ifdef CQ_HAS_TA
  #include <tiledarray.h>
#endif

namespace ChronusQ {

  namespace detail {
    inline bool& initialized_mpi_accessor() { static bool value = false; return value; }
    inline bool& initialized_tiledarray_accessor() { static bool value = false; return value; }
    inline bool& initialized_libint_accessor() { static bool value = false; return value; }
  }
  inline bool initialized_mpi() { return detail::initialized_mpi_accessor(); }
  inline void initialized_mpi(bool im) { detail::initialized_mpi_accessor() = im; }
  inline bool initialized_tiledarray() { return detail::initialized_tiledarray_accessor(); }
  inline void initialized_tiledarray(bool it) { detail::initialized_tiledarray_accessor() = it; }
  inline bool initialized_libint() { return detail::initialized_libint_accessor(); }
  inline void initialized_libint(bool it) { detail::initialized_libint_accessor() = it; }

  void generateFmTTable();

  /**
   *  \brief Initialize the ChronusQ environment
   *
   *  Sets up to the environment for a ChronusQ calculation.
   *  Currently initialized the libint2 environment and sets
   *  the default thread pool to a single thread (serial 
   *  calculation)
   */ 
  inline void initialize() {

    // Bootstrap libint2 env
    if (!libint2::initialized()) {
      libint2::initialize();
      initialized_libint(true);
    }

    // SS start
    pop_cart_ang_list();  // populate cartesian angular momentum list  
    pop_car2sph_matrix(); // populate cartesian to spherical transform matrix
    generateFmTTable();
    // SS end

    // Parallelism

    // SMP
    SetNumThreads(1);

    // MPI
#if defined(CQ_ENABLE_MPI)
    int mpi_already_initialized;
    MPI_Initialized(&mpi_already_initialized);
    if (!mpi_already_initialized) {
#ifdef _OPENMP
      int mpi_th_support;
      MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &mpi_th_support);
      assert(mpi_th_support == MPI_THREAD_MULTIPLE);
#else
      MPI_Init(NULL, NULL);
#endif
      initialized_mpi(true);
    }
#endif

    

    H5::Exception::dontPrint();

  }; // initialize

  /**
   *  \brief Finalized the ChronusQ environment
   *
   *  Cleans up anything that was initialized by ChronusQ::initialize.
   *  Currently finalized the libint2 environment.
   */ 
  inline void finalize() {

    // Finalize libint2 env
    libint2::finalize();

    // MPI
#if defined(CQ_HAS_TA)
    if (initialized_tiledarray()) {
      // Finalize TA
      TA::finalize();
      madness::finalize();
      initialized_tiledarray(false);
    }
#endif
#if defined(CQ_ENABLE_MPI)
    if (initialized_mpi()) {
      MPI_Finalize();
      initialized_mpi(false);
    }
#endif

    if (initialized_libint()) {
      libint2::finalize();
      initialized_libint(false);
    }

  }; // finalize

}; // namespace ChronusQ

