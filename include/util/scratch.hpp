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
#include <util/threads.hpp>

namespace ChronusQ {

/*
 * Helper class for managing scratch of type T 
 * in multi-threading (OpenMP) environment
 */ 
template <typename T>
class SharedMemoryScratch {
 private:
  CQMemManager& memManager_;
  const size_t nThreads_ = GetNumThreads();
  
  // actual storage of scratch
  size_t len_ = 0ul; 
  std::vector<T*> scratch_;
  
 public:
  SharedMemoryScratch() = delete;
  SharedMemoryScratch(const SharedMemoryScratch&) = delete;
  SharedMemoryScratch(SharedMemoryScratch&&) = delete;

  SharedMemoryScratch(CQMemManager& mem): memManager_(mem) {
    scratch_ = std::vector<T*>(nThreads_, nullptr); 
  }
  ~SharedMemoryScratch() { dealloc(); }

  void dealloc() {
    // root only
    if (GetThreadID() == 0) {
      for (auto& scr : scratch_) {
        if (scr) memManager_.free(scr);
      }
    }
  }

  void alloc(size_t n, std::string badAllocStr = "Scratch") {
    dealloc();
    if (GetThreadID() == 0) {
      try {
        for (auto& scr : scratch_) {
          scr = memManager_.template malloc<T>(n);
        }
      } catch (...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for " << badAllocStr 
                  <<  " (" << (n / 1e9) * GetNumThreads() * sizeof(T) << " GB)" 
                  << std::endl;
        std::cout << memManager_ << std::endl;
        throw std::bad_alloc();
      }
    }
  }

  void reserve(size_t length, std::string badAllocStr = "Scratch") {
    if (length < len_) return;
    resize(length, badAllocStr);
  }
  
  void resize(size_t length, std::string badAllocStr) {
    len_ = length;
    alloc(length, badAllocStr);
  }
  
  T* getPtr() {
    return scratch_[GetThreadID()];
  }

  T* getPtr(size_t i) {
    return scratch_[i];
  }
  
}; // class SharedMemoryScratch

} // namespace ChronusQ
