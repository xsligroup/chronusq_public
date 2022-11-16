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
#include <cxxapi/boilerplate.hpp>
// Workaround to intel compiler bug that breaks libint ERIs
#if !LIBINT2_CONSTEXPR_STATICS
  #include <libint2/statics_definition.h>
#endif

// GTest header
#include <gtest/gtest.h>

// UT Headers
#ifdef CQ_FUNC_TEST
  #include <func.hpp>
#endif

int main(int argc, char **argv) {

  // Call CQ::initialize only once
  ChronusQ::initialize();

#if defined(_CQ_GENERATE_TESTS) && defined(CQ_FUNC_TEST)
  // Register contraction environment
  ::testing::Environment * const contr_env = 
     ::testing::AddGlobalTestEnvironment(new ChronusQ::ContractEnvironment);
#endif

#ifdef CQ_ENABLE_MPI
  // Only get print from root
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();
  if(ChronusQ::MPIRank(MPI_COMM_WORLD) != 0) {
      delete listeners.Release(listeners.default_result_printer());
  }
#endif

  // Init GT and run tests
  ::testing::InitGoogleTest(&argc, argv);
  auto gt_return = RUN_ALL_TESTS();

  // Call CQ::finalize only once
  ChronusQ::finalize();

  return gt_return; // return GT result
}
