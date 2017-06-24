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

#include <func.hpp>
#include <cerr.hpp>
#include <memmanager.hpp>

using namespace ChronusQ;

template <typename T>
void CQMemManager_TEST(size_t size = 1) {


  size_t mem     = 256e6; // Default 256 MB allocation
  size_t blkSize = 2048;  // Default 2KB block size

  CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,mem,blkSize);

  size_t  trial_mem = 112e6/sizeof(T);
  CQMemManager::get().malloc<T>(trial_mem);

  size_t find_max = CQMemManager::get().max_avail_allocatable<T>(size);
  size_t max_mem = (mem - 112e6 - 1024)/(sizeof(T) * size);

  EXPECT_TRUE(find_max == max_mem);

}

TEST(CQMEMMANAGER, CQMEMMANAGER_MAX_MEM) {

  CQMemManager_TEST<char>();
  CQMemManager_TEST<dcomplex>();
  CQMemManager_TEST<double>(20);

}





