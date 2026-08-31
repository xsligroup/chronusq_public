/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
 * \brief the TensorLooper Class
 * 
 * Interface of Loop thru a tensor with N dimensions
 * 
 * IMPORTANT! This is only done with one thread
 * TODO: expand this to multiple threads
 *
 */
class TensorLooper {
 protected:
  
  size_t nDim_   = 0ul; 
  size_t index_  = 0ul; 
  size_t nTotal_ = 1ul;

  std::vector<size_t> loopIndices_;

  // stores the size of each dimension (e.g., dim1, dim2)
  // In DAS, the size of each dimension is the total number of determinants in each space
  std::vector<size_t> dimensions_;

  // stores offset of each dimension, i.e., offset1 = dim1, offset2 = dim1*dim2, etc...
  // Here, it assumes that all dimensions are continuous. It is used for looping over all
  // dimensions
  std::vector<size_t> idxOffsets_;

  // stores address offset of each dimension. This can be different from
  // those in idxOffsets_ when the dimensions are not continuous
  // It is used for figuring out the absolute address
  std::vector<size_t> addrOffsets_;
  std::vector<size_t> auxAddrOffsets_;
  
  // upper bound by total address of 2^63
  int64_t address_    = 0ul;
  int64_t auxAddress_ = 0ul;

 public:
  
  TensorLooper() = default;
  TensorLooper(const std::vector<size_t> & dims):
    dimensions_(dims), nDim_(dims.size()) {
    loopIndices_.resize(nDim_);
    // idxOffsets_.resize(nDim + 1); 
    idxOffsets_ = {1ul};
    for (auto & i: dims) {
      idxOffsets_.push_back(idxOffsets_.back() * i); 
    } 
    nTotal_ = idxOffsets_.back();
    idxOffsets_.pop_back();
  }

  TensorLooper(const TensorLooper &) = default;
  TensorLooper(TensorLooper &&)      = default;
  ~TensorLooper()                    = default;

  // virtual increments
  virtual void updateAddress() { };
  virtual void incrementAddress(size_t i) { };

  virtual void increment() {
    size_t i;
    for (i = 0ul; i < nDim_ - 1; i++) {
      if (this->loopIndices_[i] != (this->dimensions_[i] - 1)) break;
      this->loopIndices_[i] = 0ul;
    }
    ++index_;
    ++loopIndices_[i];
    incrementAddress(i);
  } 
  
  // initialization
  // or find the loop indices given a vectorized index
  void setIndex(size_t index = 0ul) {
    index_ = index;
    for (int i = nDim_ - 1; i >= 0; --i) {
      loopIndices_[i] = index / idxOffsets_[i];
      index %= idxOffsets_[i];
    }
    this->updateAddress(); 
  } 
  
  // getters
  const int64_t& address()    const { return address_; }
  const int64_t& auxAddress() const { return auxAddress_; } 
  
  const size_t& dimension(size_t i)  const { return dimensions_[i];} 

  const std::vector<size_t>& indexOffsets() const { return idxOffsets_; }
  const std::vector<size_t>& addressOffsets() const { return addrOffsets_; }
  const std::vector<size_t>& auxilaryAddressOffsets() const { return auxAddrOffsets_; }
  
  const size_t& index()  const { return index_; }
  const size_t& nTotal() const { return nTotal_; }
  const size_t& operator[](size_t i) const { return loopIndices_[i]; }  

  bool isEnd() { return index_ >= nTotal_; }
  
}; //TensorLooper

class ZeroDimTensorLooper: public TensorLooper {

 public:
  
  ZeroDimTensorLooper()                            = default;
  ZeroDimTensorLooper(const ZeroDimTensorLooper &) = default; 
  ZeroDimTensorLooper(ZeroDimTensorLooper &&)      = default; 
  ~ZeroDimTensorLooper()                           = default; 
  
  void increment() override { 
    ++this->index_; 
  }
}; //ZeroDimTensorLooper

class SingleAddrTensorLooper: public TensorLooper {
  
 protected:
   
  std::vector<int64_t> dimIncrements_;
 
 public:
  SingleAddrTensorLooper() = delete;
  SingleAddrTensorLooper(const std::vector<size_t> & dims,
    const std::vector<size_t> & offs):
    TensorLooper(dims) {
      this->addrOffsets_ = offs; 
      this->updateAddress();
      
      // use recursion relations
      // This is the true dimension offset that takes into account
      // dimensions that are not in the loop
      dimIncrements_.push_back(offs[0]);
      for (auto i = 1ul; i < nDim_; ++i) {
        dimIncrements_.push_back(
          dimIncrements_.back() + offs[i] - offs[i-1] * dims[i-1]);
      }      
  }
  SingleAddrTensorLooper(const SingleAddrTensorLooper &) = default; 
  SingleAddrTensorLooper(SingleAddrTensorLooper &&)      = default; 
  
  void updateAddress() override {
    size_t tmp = 0ul;
    for (auto i = 0ul; i < nDim_; ++i) { 
      tmp  += this->loopIndices_[i] * this->addrOffsets_[i];
    }
    this->address_ = tmp;
  }
  
  void incrementAddress(size_t i) override {
    // find current increment index and reset everything previous to curId
    this->address_+= dimIncrements_[i]; 
  }

}; // SingleAddrTensorLooper

class DoubleAddrTensorLooper: public TensorLooper {
  
protected:
  
 std::vector<int64_t> dimIncrements_;
 std::vector<int64_t> auxIncrements_;

public:
 DoubleAddrTensorLooper() = delete;
 DoubleAddrTensorLooper(const std::vector<size_t> & dims,
   const std::vector<size_t> & offs, 
   const std::vector<size_t> & auxOffs):
   TensorLooper(dims) {
     
     this->addrOffsets_ = offs; 
     this->auxAddrOffsets_ = auxOffs; 
     this->updateAddress();
 
     // use recursion relations
     dimIncrements_.push_back(offs[0]);
     auxIncrements_.push_back(auxOffs[0]);
     for (auto i = 1ul; i < nDim_; ++i) {
       dimIncrements_.push_back(
         dimIncrements_.back() + offs[i] - offs[i-1] * dims[i-1]);
       auxIncrements_.push_back(
         auxIncrements_.back() + auxOffs[i] - auxOffs[i-1] * dims[i-1]);
     }      
 }
 DoubleAddrTensorLooper(const DoubleAddrTensorLooper &) = default; 
 DoubleAddrTensorLooper(DoubleAddrTensorLooper &&)      = default; 
 
 void updateAddress() override {
   size_t tmp = 0ul, tmp2 = 0ul;
   for (auto i = 0ul; i < nDim_; ++i) { 
     tmp   += this->loopIndices_[i] * this->addrOffsets_[i];
     tmp2  += this->loopIndices_[i] * this->auxAddrOffsets_[i];
   }
   this->address_    = tmp;
   this->auxAddress_ = tmp2;
 }

 void incrementAddress(size_t i) override {
   this->address_   += dimIncrements_[i]; 
   this->auxAddress_+= auxIncrements_[i]; 
 }

}; // DoubleAddrTensorLooper

// smart constructors
inline std::shared_ptr<TensorLooper> constructTensorLooper(
  const std::vector<size_t> & dims, const std::vector<size_t> & firstOffset = {},
  const std::vector<size_t> & secondOffset = {}) {
  
  if (dims.size() == 0) {
    return std::dynamic_pointer_cast<TensorLooper>(
      std::make_shared<ZeroDimTensorLooper>());
  } else if (firstOffset.size() == 0) {
    return std::make_shared<TensorLooper>(dims);
  } else if (secondOffset.size() == 0) {
    return std::dynamic_pointer_cast<TensorLooper>(
      std::make_shared<SingleAddrTensorLooper>(dims, firstOffset));
  } else {
    return std::dynamic_pointer_cast<TensorLooper>(
      std::make_shared<DoubleAddrTensorLooper>(dims, firstOffset, secondOffset));
  }
} // constructTensorLooper

}; // namespace ChronusQ
