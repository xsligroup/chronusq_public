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
  
  // use macro for now, TODO: change to structured binding in the future
  #define UNPACK_EXCITATIONLIST_4(PTR, _0, _1, _2, _3) \
    _0 = PTR[0]; \
    _1 = PTR[1]; \
    _2 = PTR[2]; \
    _3 = PTR[3]; \
 
  #define UNPACK_EXCITATIONLIST_5(PTR, _0, _1, _2, _3, _4) \
    _0 = PTR[0]; \
    _1 = PTR[1]; \
    _2 = PTR[2]; \
    _3 = PTR[3]; \
    _4 = PTR[4];
 
  /* 
   * A data container class for excitation list to 
   * handle the indexing of the non-zero matrix 
   * elements ... 
   * This could be a three- or four-dimensional integer tensor container 
   */
  class ExcitationList {
      
  protected:
    
    int * ptr_ = nullptr;
    size_t nElement_;
    size_t nNonZero_;
    size_t nStr1_;
    size_t nStr2_;
    
    size_t N1_;
    size_t N2_;
    size_t N3_;
    size_t N_;
  
  public:

    // default Constructors
    ExcitationList() = delete;
    ExcitationList(const ExcitationList & other):
        ExcitationList(other.nElement_,
        other.nNonZero_, other.nStr1_, other.nStr2_) {
      std::copy_n(other.ptr_, N_, ptr_);
    }

    ExcitationList(ExcitationList && other):
      nElement_(other.nElement_),
      nNonZero_(other.nNonZero_), nStr1_(other.nStr1_),
      nStr2_(other.nStr2_), N1_(other.N1_), N2_(other.N2_),
      N3_(other.N3_), ptr_(other.ptr_) { other.ptr_ = nullptr; }
    
    ExcitationList(size_t nElement, 
      size_t nNZ, size_t nStr):
      nElement_(nElement), nNonZero_(nNZ), nStr1_(nStr), 
      nStr2_(1) { alloc(); }
    
    ExcitationList(size_t nElement, 
      size_t nNZ, size_t nStr1, size_t nStr2):
      nElement_(nElement), nNonZero_(nNZ), nStr1_(nStr1), 
      nStr2_(nStr2) { alloc(); }
    
    ~ExcitationList() { dealloc(); };
 
    int * pointer() { return ptr_; }
    int * pointerAtDet(size_t iDet) { return ptr_ + iDet*N2_; }
    int * pointerAtDet(size_t iDet, size_t jDet) { return ptr_ + iDet*N2_ + jDet*N3_;}
    
    const int * pointer() const { return ptr_; }
    const int * pointerAtDet(size_t iDet) const { return ptr_ + iDet*N2_; }
    const int * pointerAtDet(size_t iDet, size_t jDet) const { return ptr_ + iDet*N2_ + jDet*N3_;}
    
    int & operator [](size_t Loc) { return ptr_[Loc]; }
    int operator [] (size_t Loc) const { return ptr_[Loc]; }
    
    // method used by 3-dimensional tensor
    int & operator ()(size_t i, size_t iNZ, size_t iDet) {
        return ptr_[i + iNZ*N1_ + iDet*N2_];
    }
    int operator ()(size_t i, size_t iNZ, size_t iDet) const {
        return ptr_[i + iNZ*N1_ + iDet*N2_];
    }
    // method used by 4-dimensional tensor
    int & operator ()(size_t i, size_t iNZ, size_t iDet, size_t jDet) {
        return ptr_[i + iNZ*N1_ + iDet*N2_ + jDet*N3_];
    }
    int operator ()(size_t i, size_t iNZ, size_t iDet, size_t jDet) const {
        return ptr_[i + iNZ*N1_ + iDet*N2_ + jDet*N3_];
    }
    
    size_t nElement() const { return nElement_;}
    size_t nNonZero() const { return nNonZero_;}
    size_t nString()  const { return nStr1_*nStr2_; }
 
    std::pair<size_t,size_t> nStrings() const { return {nStr1_, nStr2_};}

    void clear()   { std::fill_n(ptr_, N_, 0); }
    
    void output(std::ostream & out) const {
      out << "Excitation List,  ";
      // three dimension tensor output
      if(nElement_ == 4) {
         std::cout << " for total number of string = " << nStr1_ << std::endl;
         for (auto i = 0ul; i < nStr1_; i++)
           prettyPrintSmart(out, "iStr = " + std::to_string(i), 
             pointerAtDet(i), nElement_, nNonZero_, nElement_);
      } else {
         std::cout << " for total number of string1 and string2 = (" << nStr1_ 
                   << "," << nStr2_ << ")" << std::endl;
         for (auto i = 0ul; i < nStr1_; i++)
         for (auto j = 0ul; j < nStr2_; j++)
           prettyPrintSmart(out, 
             "(iStr, jStr) = (" + std::to_string(i) + "," + std::to_string(j) + ")", 
             pointerAtDet(i, j), nElement_, nNonZero_, nElement_);
      } 
    }

    void alloc() {
      dealloc();
      N1_ = nElement_;
      N2_ = N1_ * nNonZero_;
      N3_ = N2_ * nStr1_;
      N_  = N3_ * nStr2_;
      this->ptr_ = CQMemManager::get().malloc<int>(N_);
    }
    
    void dealloc() { if(ptr_) CQMemManager::get().free(ptr_); }
  
  }; // ExcitationList

}; // namespace::ChronusQ
