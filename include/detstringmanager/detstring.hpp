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

namespace ChronusQ {
  

  /* 
   * Class of binary representations for electronic 
   * configurations in determinants form
   */
  class DetString {

  protected:
    
    bool * ptr_ = nullptr;     ///< Raw string storage
    size_t N_;
    size_t n1s_;

  public:
    
    // Constructors
    DetString() = delete;
    DetString(size_t n, bool b = 0): N_(n) {
      alloc();
      set(b);
    }
    DetString(size_t n, const std::vector<size_t> & electronLocs):
        DetString(n) {
      set(electronLocs);
    }
    DetString(const DetString & other):
        DetString(other.N_) {
      std::copy_n(other.ptr_, N_, ptr_);
    }
    DetString(DetString && other): N_(other.N_),
      ptr_(other.ptr_) { other.ptr_ = nullptr; }
    
    ~DetString() { dealloc(); }

    // helper functions
    void set(bool b) {
      std::fill_n(ptr_, N_, b);
      if(b) n1s_ = N_;
      else n1s_ = 0;
    }
    
    void set(const std::vector<size_t> & Locs) {
      for (const size_t & L: Locs) {
//        if(L >= N_) 
//          CErr("The bit to assign is exceeding the DetString range");
        ptr_[L] = 1;
      }
      n1s_ = Locs.size();
    }
    
    void flip(size_t Loc) { 
//      if(Loc >= N_) 
//        CErr("The bit to flip is exceeding the DetString range");
      ptr_[Loc] = not ptr_[Loc];
      if(ptr_[Loc]) n1s_++;
      else n1s_--;
    } 
    
    size_t size() const { return N_; }
    size_t n1s()  const { return n1s_; }
    size_t n0s()  const { return N_ - n1s_; }
    bool * pointer(size_t i = 0) { return ptr_ + i; }
    
    bool & operator[](size_t Loc) { return ptr_[Loc]; }   
    bool operator[] (size_t Loc) const { return ptr_[Loc]; }
    DetString & operator=( const DetString & other) {
      if (this != &other) {
        N_ = other.N_;
        n1s_ = other.n1s_;
        alloc();
        std::copy_n(other.ptr_, N_, ptr_);
      }
      return *this;
    }
    
    DetString & operator=(DetString && other) {
      if (this != &other) {
        N_ = other.N_;
        n1s_ = other.n1s_;
        dealloc();
        ptr_ = other.ptr_;
        other.ptr_ = nullptr;
      }
      return *this;
    }

    DetString operator+(const DetString & other) const { 
      DetString joint_str(N_ + other.N_); 
      joint_str.n1s_ = n1s_ + other.n1s_; 
      std::copy_n(ptr_, N_, joint_str.ptr_); 
      std::copy_n(other.ptr_, other.N_, joint_str.ptr_ + N_);
      return joint_str;
    }
    
    std::string to_string() {
      std::string s(N_, '0');
      for(auto i =0ul; i < N_; i++)
        if (ptr_[i]) s[i] = '1';

      return s;
    }
    /*
     * Determine the sign for nonzero matrix elements.
     * WARNING: it's assumed that p and q will generate 
     * a nonzero matrix element
     */ 
    int sign1e(const size_t p, const size_t q) const {
      
      bool sign = true;
      size_t l = std::min(p,q);
      size_t r = std::max(p,q);
      for (auto k = l+1; k < r; k++){
        sign ^= ptr_[k]; // using exclusive or 
      }

      return  (sign)? 1: -1; 
    } // sign1e 
    
    /* 
     * in-place 1e Excitation operation
     * p: 0->1, q: 1->0
     */
    void excitation(const size_t p, const size_t q) {
//      if (not ptr_.[q] or (ptr_[p] and p != q)) return;
      ptr_[q] = not ptr_[q];
      ptr_[p] = not ptr_[p];
    } // excitation
    
    template< class InputIt >
    void bitInfo(const bool bitType, const std::pair<size_t, size_t> & occRange, 
      InputIt bitList) const {
      auto bitList_ptr = bitList;
      for(auto i = occRange.first; i < occRange.second; i++) {
        if(bitType == ptr_[i]) {
          *bitList_ptr = i;
          bitList_ptr++;
        }
      }
    }
   
    std::vector<size_t> occupationInfo() const {
      std::vector<size_t> occList(n1s(), 0);
      bitInfo(true, {0,N_}, occList.begin());
      return occList; 
    }
    
    std::vector<size_t> virtualInfo() const {
      std::vector<size_t> virList(n0s(), 0);
      bitInfo(false, {0,N_}, virList.begin());
      return virList;
    }

    // memory management
    void alloc() {
      dealloc();
      try { ptr_ =  CQMemManager::get().malloc<bool>(N_);}
      catch(...) {
        CErr("need more memory to allocate detString" );
      }
    }
    
    void dealloc() { 
      if(ptr_) CQMemManager::get().free(ptr_); 
    }
  
  }; // ChronusQ::DetString
  
}; // namespace ChronusQ
  
