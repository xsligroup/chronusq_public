/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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

namespace ChronusQ {

/*
 * \brief ExcitationListStorage
 *   A data container class handles excitation list storage of various type
 *  
 *   As a column major three-dimentional tensor (n1, n2, n3)
 */
template <typename DetsT,
          typename std::enable_if<std::is_fundamental_v<DetsT>, int>::type = 0>
class ExcitationListStorage {
 private:
  DetsT * data_ = nullptr;
   
  size_t n1_ = 0;
  size_t n2_ = 0;
  size_t n3_ = 0;
  size_t n12_ = 0;
  size_t n123_ = 0;
 
 public:
  ExcitationListStorage() = default; 

  ExcitationListStorage(size_t n1, size_t n2, size_t n3) {
    resize(n1, n2, n3);
  }

  ExcitationListStorage(const ExcitationListStorage& other):
      ExcitationListStorage(other.n1_, other.n2_, other.n3_) {
    std::copy_n(other.data_, n123_, data_);
  }
     
  ExcitationListStorage(ExcitationListStorage&& other):
      n1_(other.n1_), n2_(other.n2_), 
      n3_(other.n3_), n12_(other.n12_), n123_(other.n123_) {
    data_ = other.data_;
    other.data_ = nullptr;
  }

  ~ExcitationListStorage() {
     dealloc();
  } 
  
  void tryAlloc() {
    if (not data_) {
      try {
        data_ = CQMemManager::get().malloc<DetsT>(n123_);
      } catch (...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for ExcitationListStorage" 
                  <<  " (" << (n123_ / 1e9) * sizeof(DetsT) << " GB)" 
                  << std::endl;
        std::cout << CQMemManager::get() << std::endl;
        CErr();
      }
    } 
  }

  void dealloc() {
    if (data_) CQMemManager::get().free(data_);
  }
  
  void resize(size_t n1, size_t n2, size_t n3) {
    n1_ = n1;
    n2_ = n2;
    n3_ = n3;
    n12_ = n1 * n2;
    n123_ = n12_ * n3;
  }

  size_t totalDimension() const { return n123_; }
  size_t dimension1() const { return n1_; }
  size_t dimension2() const { return n2_; }
  size_t dimension3() const { return n3_; }
  
  void clear() { if (data_) std::fill_n(data_, n123_, DetsT(0)); }

  const DetsT* pointer(size_t i, size_t j, size_t k) const {
    return data_ + i + j * n1_ + k * n12_;
  }
  
  DetsT* pointer(size_t i, size_t j, size_t k) {
    return data_ + i + j * n1_ + k * n12_;
  }

  const DetsT& operator() (size_t i, size_t j, size_t k) const {
    return *pointer(i, j, k);
  }

  DetsT& operator() (size_t i, size_t j, size_t k) {
    return *pointer(i, j, k);
  }
}; // class ExcitationListStorage

} // namespace ChronusQ
