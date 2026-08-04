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
#ifdef CQ_ENABLE_SPARSE
  #include <itersolver/cqSparseMatrix.hpp>
#endif

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


/*
 * \brief LocalCISparseVectorsView class
 *
 * A class for easier access categorical raw pointers
 *
 */
#ifdef CQ_ENABLE_SPARSE
template <typename MatsT>
class LocalCISparseVectorsView {

 private:

  size_t shift_;//This is the shift of the DistributedVector, this replaces getLocalPtr(shift)

  size_t localLen_;
  size_t size_;

  std::vector<size_t> nonZerosArray;
  LLSparseMatrix<MatsT>* vecs_ = nullptr;
  std::vector<LLSparseMatrix<MatsT>> vecsByCat_;

  size_t localCategoryBegin_;
  size_t localCategoryEnd_;
  std::vector<size_t> localCategoryOffsets_;

 public:

  LocalCISparseVectorsView() = default;
  LocalCISparseVectorsView(const LocalCISparseVectorsView&) = default;
  LocalCISparseVectorsView(LocalCISparseVectorsView&&) = default;

  LocalCISparseVectorsView(LLSparseMatrix<MatsT>& vecs, size_t localLength, size_t size,
      size_t localCategoryBegin, size_t localCategoryEnd,
      const std::vector<size_t>& localCategoryOffsets, size_t shift = 0):
      localLen_(localLength), size_(size), localCategoryBegin_(localCategoryBegin),
      localCategoryEnd_(localCategoryEnd),
      localCategoryOffsets_(localCategoryOffsets), vecs_(&vecs), shift_(shift) {

  }

  LocalCISparseVectorsView(char* vecsBuffer, size_t dataNonZeros, size_t localLength, size_t size,
      size_t localCategoryBegin, size_t localCategoryEnd,
      const std::vector<size_t>& localCategoryOffsets, size_t shift = 0):
      localLen_(localLength), size_(size), localCategoryBegin_(localCategoryBegin),
      localCategoryEnd_(localCategoryEnd),
      localCategoryOffsets_(localCategoryOffsets), shift_(shift) {

      toVecsByCat(vecsBuffer, dataNonZeros);
  }
  
  size_t size() const { return size_; }
  size_t localLength() const { return localLen_; }
  size_t localCategoryBegin() const { return localCategoryBegin_; }
  size_t localCategoryEnd() const { return localCategoryEnd_; }

  bool containsLocalCategory(size_t i) const {
    return i >= localCategoryBegin_ and i < localCategoryEnd_;
  }

  size_t getCatOffest(size_t iCat) const {
    return localCategoryOffsets_[iCat - localCategoryBegin_];
  }

  size_t getShift() const{
    return shift_;
  }
  
  auto getVecsPtr() {
    return vecs_;
  }

  const auto getVecsPtr() const {
    return vecs_;
  }

  const auto& getVecsByCat() const {
    return vecsByCat_;
  }

  void toVecsByCat(size_t nVec) {
	
    vecsByCat_.resize(localCategoryEnd_ - localCategoryBegin_);
    #pragma omp parallel for
    for(auto& mat : vecsByCat_)
      mat.resize(nVec);

    #pragma omp parallel for schedule(static) default(shared)
    for (size_t col = 0; col < nVec; ++col) {
      size_t iCat = localCategoryBegin_;
      for(auto it = vecs_->cbegin(shift_ + col); it != vecs_->cend(shift_ + col); it++) {
        size_t adjustedOffset = iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1);
        while(it->first >= adjustedOffset) {
          iCat++;
          adjustedOffset = iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1);
        }
        size_t iCatSize = (iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1)) - getCatOffest(iCat);
	vecsByCat_[iCat - localCategoryBegin_].sortedInsert(it->first - getCatOffest(iCat), col, it->second);
      }
    } 
  }
  
  void toVecsByCat(char* vecsBuffer, size_t nonZeros) {

    nonZerosArray.resize(size_);

    vecsByCat_.resize(localCategoryEnd_ - localCategoryBegin_);
    #pragma omp parallel for
    for(auto& mat : vecsByCat_)
      mat.resize(size_);

    size_t* prevBuffRowCols = (size_t*)vecsBuffer;
    MatsT* prevBuffVals = (MatsT*)(vecsBuffer + nonZeros * sizeof(size_t) * 2);

    for(size_t element = 0; element < nonZeros; element++) {

      size_t row = *prevBuffRowCols;
      size_t col = *(prevBuffRowCols + 1);
      nonZerosArray[col]++; //update num of non-zeros at column col
      MatsT value = *prevBuffVals;

      size_t iCat = localCategoryBegin_;
      size_t adjustedOffset = iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1);


      while(row >= adjustedOffset) {
        iCat++;
        adjustedOffset = iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1);
      }

      size_t iCatSize = (iCat == localCategoryEnd_ - 1? localLen_ : getCatOffest(iCat + 1)) - getCatOffest(iCat);
      vecsByCat_[iCat - localCategoryBegin_].sortedInsert(row - getCatOffest(iCat), col, value);

      prevBuffRowCols += 2;
      prevBuffVals++;
    }
  }

  size_t nonZeros(size_t iVec) {
    if(iVec >= nonZerosArray.size()) {
      CErr("localView trying to get non-zeros of non-existing column");
    }
    return nonZerosArray[iVec];
  }

  void reduceSigmaIntermediate(HashSparseMatrix<MatsT>& B) {
    vecs_->setFromHashSparseMatrix(B, shift_, 0, size_);
  }

  void reduceSigmaIntermediate(HashSparseMatrix<MatsT>& A, HashSparseMatrix<MatsT>& B) {
    vecs_->setFromTwoHashSparseMatrices(A, B, shift_, 0, size_);
  }

}; // class DASCISparseVectorViewer
#endif

} // namespace ChronusQ
