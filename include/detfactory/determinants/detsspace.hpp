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

#ifndef DETFACTORY_DETERMINANTS
#error This file may only be included from detfactory/determinants.hpp
#endif

#include <posthartreefock/base/space.hpp>
#include <detfactory/util.hpp>
#include <detfactory/determinants/detscategory.hpp>
#include <detfactory/localcivectorsview.hpp>
#include <itersolver/solvervectors.hpp>
#include <itersolver/solvervectorsimpl.hpp>

namespace ChronusQ {
  
/*
 * \brief CategoricalSpace Class
 *
 * defines determinant basis based on a list of categories
 * with functions to construct non-zero element between active spaces
 *
 */ 
class CategoricalSpace {

 protected:
  
  const std::vector<ActiveSpaceParameters>& activeSpaces_;

  // Unique category ID is a string consisting of space occupation numbers
  // categoricalIDToIndexMap_ is a hash table that maps categorical string id to its index in storage
  std::unordered_map<std::string, size_t> categoricalIDToIndexMap_;
  std::vector<std::shared_ptr<DeterminantCategory>> categories_;

  size_t nDeterminants_ = 0ul;

public:
  
  CategoricalSpace() = delete;
  CategoricalSpace(const std::vector<ActiveSpaceParameters>& spaces):
      activeSpaces_(spaces) { }
  
  CategoricalSpace(const CategoricalSpace &) = default;
  CategoricalSpace(CategoricalSpace &&)      = default;
  
  // getters
  size_t nDeterminants() const { return nDeterminants_; }
  size_t nCategories()   const { return categories_.size(); }
  const std::vector<ActiveSpaceParameters>& activeSpaces() const { return activeSpaces_; }

  std::shared_ptr<DeterminantCategory> getCategory(size_t i) {
    if (i >= categories_.size() ) return nullptr;
    else return categories_[i];
  }

  size_t getCategoryIdx(size_t i) {
    size_t idx = 0;
    if ( i > nDeterminants_)
      CErr("CategoricalSpace::getCategoryIdx error. Searching for i > nDets!");
    while (idx < categories_.size()) {
      if (i >= categories_[idx]->offset() and i < categories_[idx]->offset() + categories_[idx]->nDeterminants()){
        return idx;
      } else {
        idx++;
      }
    }
    CErr("CategoricalSpace::getCategoryIdx error. Could not locate Category");
  }
  
  std::shared_ptr<const DeterminantCategory> getCategory(size_t i) const {
    return const_cast<CategoricalSpace*>(this)->getCategory(i); 
  }

  size_t getCategoryIdx(size_t i) const {
    return const_cast<CategoricalSpace*>(this)->getCategoryIdx(i); 
  }
  
  size_t getCategoricalIndex(const std::string& id) const {
    if (categoricalIDToIndexMap_.count(id) == 0) return categories_.size();
    return categoricalIDToIndexMap_.at(id);
  }
  
  std::string generateCategoryID(const DeterminantCategory& cat) {
    std::string id = "";
    for (const auto& g: cat.detGroups()) {
      id += std::to_string(g.nElectrons()) + "-";
    }
    id.pop_back();
    return id;
  };

  // create a new category ID based on the space occupations
  // Category ID is a string consisting of space occupation number (number of electrons
  // in each space)
  std::string generateCategoryID(const std::vector<size_t>& spaceOcc) {
    std::string id = "";
    for (const auto & iOcc: spaceOcc) {
      id += std::to_string(iOcc) + "-";
    }
    id.pop_back();
    return id;
  };
  
  void addCategory(std::shared_ptr<DeterminantCategory> newCategory) {
    if (newCategory == nullptr) return;
    assert(newCategory->detGroups().size() == activeSpaces_.size());
    auto newCategoryID = generateCategoryID(*newCategory);
    if (categoricalIDToIndexMap_.count(newCategoryID) == 0) {
      categoricalIDToIndexMap_.emplace(newCategoryID, categories_.size());
      categories_.emplace_back(newCategory);
      newCategory->setOffset(nDeterminants_);
      nDeterminants_ += newCategory->nDeterminants();
    }
  };
  
  // Build and add user-defined reference categories.
  // Each category is a collection of distributed active space with initially
  // defined space occupations
  void addReferenceCategory(const std::vector<size_t>& spaceOcc) {
    assert(spaceOcc.size() == activeSpaces_.size());
    auto newCategoryID = generateCategoryID(spaceOcc);
    // search to see if this category already exists
    // only new category is added
    if (categoricalIDToIndexMap_.count(newCategoryID) == 0) {
      auto newCategory = buildFullDeterminantCategory(activeSpaces_, spaceOcc);
      if (newCategory != nullptr) {
        categoricalIDToIndexMap_.emplace(newCategoryID, categories_.size());
        categories_.emplace_back(newCategory);
        newCategory->setOffset(nDeterminants_);
        nDeterminants_ += newCategory->nDeterminants();
      }
    }
  }

  void expandCategory(int exLevel) {
    std::vector<std::shared_ptr<DeterminantCategory>> existingCats(categories_.begin(), categories_.end());
    std::vector<size_t> spaceOcc;
    for (const auto &cat: existingCats) {
      for (auto i = 0ul; i < activeSpaces_.size(); ++i)
        for (auto j = 0ul; j < activeSpaces_.size(); ++j) {
          // Generate complete excitation between space within a same iDASGroup
          // i.e., a DAS group is a complete active space
          if ((activeSpaces_[i].iDASGroup == activeSpaces_[j].iDASGroup) and (i != j)) {
            spaceOcc = cat->SpaceOccupations();
            while (spaceOcc[i] < activeSpaces_[i].nOrbitals and spaceOcc[j] > 0) {
              spaceOcc[i]++;
              spaceOcc[j]--;
              addReferenceCategory(spaceOcc);
            }
            //spaceOcc = cat->nElectronsInEachSpace();
            //while (spaceOcc[j] < activeSpaces_[j].nOrbitals and spaceOcc[i] > 0) {
            //  spaceOcc[j]++;
            //  spaceOcc[i]--;
            //  addFullDetsCategory(spaceOcc);
            //}
          }
        }
    }
    expandCategoryWithInterGroupExcitation(exLevel);
  }

  // Build categories based on excitation operators arising from the base categories.
  // Note that new categories are generated with single excitation at a time from the previous set
  // until reaching the target excitation. Redundant categories are removed.
  void expandCategoryWithInterGroupExcitation(int exLevel = 1) {

    if (exLevel == 1) {
      std::vector<std::shared_ptr<DeterminantCategory>> existingCats(categories_.begin(), categories_.end());
      for (const auto& cat: existingCats) {
        auto spaceOcc = cat->SpaceOccupations();
        for(auto i = 0ul; i < spaceOcc.size(); ++i)
        for(auto j = 0ul; j < spaceOcc.size(); ++j) {
          // ignore same space, same group, or empty space
          if (activeSpaces_[i].iDASGroup == activeSpaces_[j].iDASGroup) continue;

          if (i == j or spaceOcc[j] == 0ul) continue;
          spaceOcc[i]++;
          spaceOcc[j]--;
          //TODO: add capabilities to add non FullDetsCategory in the future
          addReferenceCategory(spaceOcc);
          spaceOcc[j]++;
          spaceOcc[i]--;
        }
      }
    } else if (exLevel < 0 ) {
       // Generate all possible categories under the condition defined by minOcc and maxOcc in ActiveSpaceParameters.
       size_t nCat = 0ul;
       // Keep expanding categories until no more new category is added
       while (nCat != categories_.size()) {
         nCat = categories_.size();
         expandCategoryWithInterGroupExcitation();
       }
    } else {
      // If the excitation level is greater than 1, build categories using multiple single excitations until reach
      //   the desired excitation level.
      for (auto i = 0ul; i < exLevel; ++i)
        expandCategoryWithInterGroupExcitation();
    }

  } // expandCategoryWithInterSpaceExcitation
  
  void output(std::ostream& os, const std::string& s = "", 
    bool printMapping = false) const { 
    std::string outputStr;
    
    if (s == "") {
      outputStr = "  * Determinant Space";
    } else {
      outputStr = "  * Determinant Space Of " + s;
    }

    os << outputStr << std::endl;
    size_t counter = 0ul;
    for (const auto& cat: categories_) {
      os << "    - Category "   << std::setw(5)  << counter 
         << ": nDets = "  << std::setw(10) << cat->nDeterminants()  
         << ", offset = " << std::setw(10) << cat->offset()
         << ", Occupations:";
      for (const auto& g: cat->detGroups()) {
        os << " " << std::setw(3) << g.nElectrons(); 
      }
      os << std::endl;
      if (printMapping) os << (*cat) << std::endl;
      counter++;
    }
    os << "    >>>>Total number of Determinants = " << nDeterminants() << std::endl;
    
    // This should almost never be hit, but we guard against it just in case.
    if (GetNumThreads() > nDeterminants())
      CErr("The number of threads is greater than the number of determinants. Please reduce nsmp for this calculation.");

    os << std::endl;
    if (distributedAccumulatedNCategories_.size() > 0) {
      os << "    - Distributed Map:" << std::endl;
      size_t accNCats = 0ul;
      for (auto i = 0ul; i < distributedAccumulatedNCategories_.size(); ++i) {
        os << "      $ Node " << std::setw(3) << i
           << " handles categories " << std::setw(5) << accNCats
           << " ~ " << std::setw(5) << distributedAccumulatedNCategories_[i] - 1
           << ", contains nDets = " << distributedCategoryLengths_[i]
           << std::endl;
        accNCats = distributedAccumulatedNCategories_[i];
      }
    }
    
    return;
  }
 
 protected:
  
  // Distributed Map and lends for MPI
  // used to reference in DASCIVectors 
  std::vector<size_t> distributedAccumulatedNCategories_;
  std::vector<size_t> distributedCategoryLengths_;
    
  // local variables
  size_t localCategoryBegin_;
  size_t localCategoryEnd_;
  std::vector<size_t> localCategoryOffsets_;
  size_t localNDeterminants_;
 
 public:
  
  // getters
  const std::vector<size_t> & distributedAccumulatedNCategories() const { 
      return distributedAccumulatedNCategories_; 
  }
  const std::vector<size_t> & distributedCategoryLengths() const { 
      return distributedCategoryLengths_; 
  }
  size_t localCategoryBegin() const { return localCategoryBegin_; }
  size_t localCategoryEnd() const { return localCategoryEnd_; }
  size_t localNDeterminants() const { return localNDeterminants_; }
  const std::vector<size_t>& localCategoryOffsets() const { return localCategoryOffsets_; }
  
  bool containsLocalCategory(size_t i) const {
    return i >= localCategoryBegin_ and i < localCategoryEnd_;  
  }
  
  size_t getLocalCategoryOffset(size_t i) const {
    if (not containsLocalCategory(i)) {
      CErr("Category " + std::to_string(i) + " is not on this Node");
    }
    return localCategoryOffsets_[i - localCategoryBegin_];
  }

  // generate CIVector function
  template <typename MatsT>
  std::shared_ptr<DistributedVectors<MatsT>> constructDistributedCIVectors(
      MPI_Comm comm, size_t size) const {
    return std::make_shared<DistributedVectors<MatsT>>(comm, distributedCategoryLengths_, size);
  }
  
  template <typename MatsT>
  LocalCIVectorsView<MatsT> createLocalCIVectorsView(MatsT* data, size_t nVec, MPI_Comm comm, int iNode) const {
    if (MPIRank(comm) == iNode) {
      return createLocalCIVectorsView(data, nVec);
    } else {
      // create alternative local variables
      size_t localCategoryBegin = (iNode == 0ul) ? 0ul : distributedAccumulatedNCategories_[iNode - 1];
      size_t localCategoryEnd = distributedAccumulatedNCategories_[iNode];
      size_t localNDeterminants = distributedCategoryLengths_[iNode];
      std::vector<size_t> localCategoryOffsets = {0ul};
      for (auto i = localCategoryBegin; i < localCategoryEnd - 1; ++i) {
         localCategoryOffsets.push_back(localCategoryOffsets.back() + categories_[i]->nDeterminants());
      }
      return LocalCIVectorsView<MatsT>(data, localNDeterminants, nVec, 
          localCategoryBegin, localCategoryEnd, localCategoryOffsets);
    }
  }
  
  template <typename MatsT>
  LocalCIVectorsView<MatsT> createLocalCIVectorsView(MatsT* data, size_t nVec) const {
    return LocalCIVectorsView<MatsT>(data, localNDeterminants_, nVec, 
        localCategoryBegin_, localCategoryEnd_, localCategoryOffsets_);
  }
      
  template <typename MatsT>
  LocalCIVectorsView<MatsT> createLocalCIVectorsView(DistributedVectors<MatsT>& ciVectors, 
      size_t shift, size_t nVec) const {
    return createLocalCIVectorsView(ciVectors.getLocalPtr(shift), nVec);
  }
  
  template <typename MatsT>
  LocalCIVectorsView<const MatsT> createLocalCIVectorsView(const DistributedVectors<MatsT>& ciVectors, 
      size_t shift, size_t nVec) const {
    return createLocalCIVectorsView(ciVectors.getLocalPtr(shift), nVec);
  }
  
  // initalizer
  void initializeDistributedCatMap(MPI_Comm comm) {
    
    size_t nNodes = MPISize(comm);
     
    // try best to evenly distribute the blocks across different nodes
    std::vector<size_t> categoryAssignment(categories_.size(), 0ul);
    
    if (MPIRank(comm) == 0) {
      std::vector<size_t> nodeIdHeap(nNodes);
      std::vector<size_t> nodeNDets(nNodes);
      std::iota(nodeIdHeap.begin(), nodeIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) { 
          return nodeNDets[i] > nodeNDets[j]; 
      };

      for (auto i = 0ul; i < categoryAssignment.size(); ++i) {
        categoryAssignment[i] = nodeIdHeap[0];
        nodeNDets[nodeIdHeap[0]] += categories_[i]->nDeterminants();
        std::make_heap(nodeIdHeap.begin(), nodeIdHeap.end(), comp);
      }
    }
    
    MPIBCast(categoryAssignment.data(), categories_.size(), 0, comm); 

    // populate distributed accumulated NCats 
    std::vector<size_t> nodeNCats(nNodes, 0ul);
    for (const auto& iAssignment : categoryAssignment) {
      nodeNCats[iAssignment]++;
    }
    std::vector<size_t> iNodeId({0ul});
    for (const auto& nCats : nodeNCats) {
      iNodeId.push_back(iNodeId.back() + nCats); 
      distributedAccumulatedNCategories_.push_back(iNodeId.back());
    }
    
    // bucket sort all categories based on category Assignment
    // and make sure that the categories are continuous in each nodes
    std::vector<std::shared_ptr<DeterminantCategory>> 
        sortedCategories(categories_.size(), nullptr);

    for (auto i = 0ul; i < categories_.size(); ++i) {
      auto iNode = categoryAssignment[i];
      auto j = iNodeId[iNode];
      sortedCategories[j] = categories_[i];
      iNodeId[iNode]++;
    }
    
    // rebuild category offsets and categoricalIDToIndexMap_
    size_t NDets = 0ul;
    size_t iNode = 0ul;
    size_t iNodeNDets = 0ul;
    for (auto i = 0ul; i < sortedCategories.size(); ++i) {
      sortedCategories[i]->setOffset(NDets);
      NDets += sortedCategories[i]->nDeterminants();
      iNodeNDets += sortedCategories[i]->nDeterminants();
      auto categoryId = generateCategoryID(*sortedCategories[i]);
        categoricalIDToIndexMap_[categoryId] = i;
      if ((i + 1) == distributedAccumulatedNCategories_[iNode]) {
        distributedCategoryLengths_.push_back(iNodeNDets);
        iNodeNDets = 0ul;
        iNode++;
      }
    }
    categories_ = std::move(sortedCategories);
    
    // populate local variales
    iNode =  MPIRank(comm);
    localCategoryBegin_ = (iNode == 0ul) ? 0ul : distributedAccumulatedNCategories_[iNode - 1];
    localCategoryEnd_ =  distributedAccumulatedNCategories_[iNode];
    localNDeterminants_ = distributedCategoryLengths_[iNode];
  
    localCategoryOffsets_ = {0ul};
    for (auto i = localCategoryBegin_; i < localCategoryEnd_ - 1; ++i) {
      localCategoryOffsets_.push_back(localCategoryOffsets_.back() + categories_[i]->nDeterminants());
    }
  } // initializeDistributedCatMap

}; // CategoricalSpace

inline std::ostream& operator<<(std::ostream& os, const CategoricalSpace& space) {
  space.output(os, "", true);
  return os;
}

} // namespace ChronusQ
