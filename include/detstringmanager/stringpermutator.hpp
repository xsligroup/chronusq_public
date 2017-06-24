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
#include <memmanager.hpp>
#include <util/math.hpp>

namespace ChronusQ {
  

  /* 
   * Brief class to generate all permutation of strings
   * with n permuted bits within m total bits
   */ 
  class StringPermutator {
  
    size_t m_;      // number of total bits
    size_t n_;      // number of occupied bits  
    size_t nStr_;   // number of string in total
    std::vector<size_t> pos_; // permuted bits locations
    
    size_t str_index_;  // starting from 0
    size_t cur_ptr_;
  
  public:
    
    // disable default constructors
    StringPermutator() = delete;
    StringPermutator(const StringPermutator &) = default;
    StringPermutator(StringPermutator && )     = default;
    
    StringPermutator(size_t m, size_t n): m_(m), n_(n) {
      nStr_ = Comb(m,n);
      pos_.resize(n);
      reset_seed();
    }
    
    void reset_seed() {
      for(auto i = 0ul; i < n_; i++) pos_[i] = i; 
      cur_ptr_ = n_ - 1;
      str_index_ = 0;
    }
    
    size_t stringIndex() const { return str_index_; }
    size_t maxNString()  const { return nStr_; }
    std::vector<size_t> permutedPositions() const { return pos_; }
    std::vector<size_t> & permutedPositions() { return pos_; }

    size_t next() {
      
      if (str_index_ == nStr_ - 1) return nStr_;
      
      if (cur_ptr_ == n_ - 1) {
        pos_[cur_ptr_] += 1;
      } else {
        auto k = pos_[cur_ptr_];
        for (auto i = cur_ptr_, j = 1ul; i < n_; i++, j++) 
          pos_[i] = k + j;
        cur_ptr_ = n_ - 1;
      }
      
      while (pos_[cur_ptr_] + (n_ - cur_ptr_) == m_) cur_ptr_--;
      
      str_index_++;
       
      return str_index_;
    }

    void fullStringList(size_t * strList) {
      
      reset_seed();
      auto strList_ptr = strList;
      std::copy_n(pos_.begin(), n_, strList_ptr);
      
      while (next() < nStr_) { 
        strList_ptr += n_;
        std::copy_n(pos_.begin(), n_, strList_ptr);
      }
    }

  }; // ChronusQ::DetString
  
}; // namespace ChronusQ
