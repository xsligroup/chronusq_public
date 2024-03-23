/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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
 * \brief LocalCIVectorsView class  
 *
 * A class for easier access categorical raw pointers
 *
 */
template <typename MatsT>
class LocalCIVectorsView {

 private:

  MatsT* data_ = nullptr;
  size_t localLen_;
  size_t size_;
  
  size_t localCategoryBegin_;
  size_t localCategoryEnd_;
  std::vector<size_t> localCategoryOffsets_;
 
 public:
  LocalCIVectorsView() = default;
  LocalCIVectorsView(const LocalCIVectorsView&) = default;
  LocalCIVectorsView(LocalCIVectorsView&&) = default;
  
  LocalCIVectorsView(MatsT* data, size_t localLength, size_t size,
      size_t localCategoryBegin, size_t localCategoryEnd, 
      const std::vector<size_t>& localCategoryOffsets):
      localLen_(localLength), size_(size), localCategoryBegin_(localCategoryBegin), 
      localCategoryEnd_(localCategoryEnd), 
      localCategoryOffsets_(localCategoryOffsets), data_(data) {
  }
  
  size_t size() const { return size_; }
  size_t localLength() const { return localLen_; }
  size_t localCategoryBegin() const { return localCategoryBegin_; }
  size_t localCategoryEnd() const { return localCategoryEnd_; }

  bool containsLocalCategory(size_t i) const {
    return i >= localCategoryBegin_ and i < localCategoryEnd_;  
  }
  
  MatsT* getCategoryPointer() const {
    return data_;
  }

  MatsT* getCategoryPointer(size_t iCat, size_t iVec = 0) const {
    if (not containsLocalCategory(iCat) or iVec >= size_) {
       CErr("Can't access pointer that exceeing the definition in LocalCIVectorsView");
    }
    return data_ + iVec * localLen_ + localCategoryOffsets_[iCat - localCategoryBegin_];
  }
  
  void syncViewSchemeToNode(int i, MPI_Comm comm) {
    // patch all data for broadcast
    localCategoryOffsets_.push_back(localLen_);
    localCategoryOffsets_.push_back(localCategoryBegin_);
    localCategoryOffsets_.push_back(localCategoryEnd_);
    size_t nBCast = localCategoryOffsets_.size();
    MPIBCast(nBCast, i, comm);
    localCategoryOffsets_.resize(nBCast);
    MPIBCast(localCategoryOffsets_.data(), nBCast, i, comm);
    
    // unpatch all the data
    localCategoryEnd_ = localCategoryOffsets_.back();
    localCategoryOffsets_.pop_back();
    localCategoryBegin_ = localCategoryOffsets_.back(); 
    localCategoryOffsets_.pop_back();
    localLen_ = localCategoryOffsets_.back();
    localCategoryOffsets_.pop_back();
  }

}; // class DASCIVectorViewer

} // namespace ChronusQ
