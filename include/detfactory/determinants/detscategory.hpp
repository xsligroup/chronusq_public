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

#ifndef DETFACTORY_DETERMINANTS
#error This file may only be included from detfactory/determinants.hpp
#endif

#include <detfactory/determinants/detsgroup.hpp>
#include <posthartreefock/base/space.hpp>
#include <detfactory/util.hpp>

namespace ChronusQ {

/*
 * \brief the DeterminantCategory Class
 *
 * Base class of FullDetsCategory and SelectedDetsCategory
 *
 * represents determinants in multiple complete orbital space
 * with nOrbitals and nElectrons
 *
 */ 

class DeterminantCategory {
 protected:
  std::vector<DeterminantGroup> detGroups_;
  // number of determinants in this category
  size_t nDeterminants_ = 0ul;
  size_t offset_ = 0ul; // offset in a space
 
 public:
  DeterminantCategory() = delete;
  DeterminantCategory(const DeterminantCategory &) = default;
  DeterminantCategory(DeterminantCategory &&) = default;
  ~DeterminantCategory() = default;
  
  explicit DeterminantCategory(const std::vector<DeterminantGroup>& detGroups):
  detGroups_(detGroups) {
    assert(detGroups.size() > 0);
  }

  explicit DeterminantCategory(
      const std::vector<std::array<size_t, 2>> & groupsInfo) {
    assert(groupsInfo.size() > 0);
    for (const auto& gInfo: groupsInfo) {
      detGroups_.emplace_back(gInfo);
    }
  }
  
  DeterminantCategory(const std::vector<size_t> & nSpaceOcc,
                      const std::vector<size_t> & nOrbitals) {
    assert(nOrbitals.size() == nSpaceOcc.size());
    assert(nOrbitals.size() > 0);
    for(auto i = 0ul; i < nSpaceOcc.size(); i++) {
      detGroups_.emplace_back(nSpaceOcc[i], nOrbitals[i]);
    }
  }

  const std::vector<DeterminantGroup>& detGroups() const {
    return detGroups_;
  }
  
  // offset control
  size_t offset() const { return offset_; }
  void setOffset(size_t off) { offset_ = off; }
  
  // getters
  size_t nDeterminants() const { return nDeterminants_; }
  
  std::vector<size_t> SpaceOccupations() const {
    std::vector<size_t> nEs;
    for (const auto& g : detGroups_) nEs.push_back(g.nElectrons());
    return nEs; 
  } 
  
  //std::vector<size_t> nOrbitalsInEachSpace() const {
  //  std::vector<size_t> nOs;
  //  for (const auto& g : detGroups_) nOs.push_back(g.nOrbitals());
  //  return nOs;
  //}

  double signOffsetBetween(size_t i, size_t j) const {
    size_t nEs = 0;
    for (auto k = std::min(i, j) + 1; k < std::max(i, j); ++k) {
      nEs += detGroups_[k].nElectrons();
    }
    if (nEs % 2 == 1) return -1.;
    return 1.;
  }    
  
  virtual std::ostream& printAllDeterminants(std::ostream& os) const = 0; 

};// class DeterminantCategory

// interface delacration
template <typename DetsT>
class FullDetsCatAddresser;

template <typename DetsT>
class FullDetsCatGenerator;

class FullDeterminantCategory: public DeterminantCategory {
 
 protected:
  std::vector<size_t> groupOffsets_;

 public:
  FullDeterminantCategory() = delete;
  FullDeterminantCategory(const FullDeterminantCategory &) = default;
  FullDeterminantCategory(FullDeterminantCategory &&)      = default;
  ~FullDeterminantCategory()                        = default;
  
  explicit FullDeterminantCategory(const std::vector<DeterminantGroup>& groups):
          DeterminantCategory(groups) {
    initialization();
  }
  
  explicit FullDeterminantCategory(
      const std::vector<std::array<size_t, 2>>& groupsInfo):
          DeterminantCategory(groupsInfo) {
    initialization();
  }
  
  FullDeterminantCategory(const std::vector<size_t>& nSpaceOcc,
                          const std::vector<size_t>& nOrbitals):
          DeterminantCategory(nSpaceOcc, nOrbitals) {
    initialization();
  }

  void initialization() {
    groupOffsets_ = {1ul};
    for (const auto & g: this->detGroups_) {
      groupOffsets_.push_back(groupOffsets_.back() * g.nDeterminants());
    }
    this->nDeterminants_ = groupOffsets_.back();
  }

  const std::vector<size_t>& groupOffsets() const { return groupOffsets_; }

  // separate and group non-excitation spaces.
  // continuous non-excitation spaces are grouped together and the total dimension
  // is equal to the product of individual dimensions.
  // e.g. [nDet1][nDet2][nDet3][nDet4][nDet5][nDet6][nDet7] is consolidated
  // into [nDet12][nDet3][nDet4][nDet5][nDet67] where nDet12=nDet1*nDet2 and nDet67=nDet6*nDet7
  void separateExAndNonExDimensions(
      const std::unordered_set<size_t>& exSpaces,
      std::vector<size_t>& exOffs, std::vector<size_t>& exDims, 
      std::vector<size_t>& nonExOffs, std::vector<size_t>& nonExDims) const {
    
    exOffs.clear();
    exDims.clear();
    nonExOffs.clear();
    nonExDims.clear();

    bool continuousNonExcitation = false;
    for (auto i = 0ul; i < detGroups_.size(); ++i) {
      size_t nDets_i = this->detGroups_[i].nDeterminants();

      if (exSpaces.count(i) == 1){
        // if space i is one of the excitation spaces
        exDims.push_back(nDets_i);
        exOffs.push_back(groupOffsets_[i]);
        continuousNonExcitation = false;
      } else {
        // if space i is not one of the excitation spaces
        if (nDets_i == 1) continue;
        if (continuousNonExcitation) {
          // if the current space is a continuous block of non-excitation spaces
          // the total dimension is a product of individual dimensions
          nonExDims.back() *= nDets_i;
        } else {
          nonExOffs.push_back(groupOffsets_[i]);
          nonExDims.push_back(nDets_i);
        }
        // if the current space is a continuous block of non-excitation spaces
        continuousNonExcitation = true;
      }
    }
  }

  std::ostream& printAllDeterminants(std::ostream& os) const override;

  template <typename DetsT>
  FullDetsCatAddresser<DetsT> addresser() const {
     return FullDetsCatAddresser<DetsT>(*this);
  } 

  template <typename DetsT>
  FullDetsCatGenerator<DetsT> generator() const {
     return FullDetsCatGenerator<DetsT>(*this);
  } 

}; // class FullDetsCategory

template <typename DetsT = uint64_t>
class FullDetsCatAddresser {
 
 private:
  std::vector<DetGroupAddresser<DetsT>> groupAddressers_;
  std::vector<size_t> groupOffsets_;

 public:
  FullDetsCatAddresser() = default;
  FullDetsCatAddresser(const FullDetsCatAddresser&) = default;
  FullDetsCatAddresser(FullDetsCatAddresser&&) = default;
  
  explicit FullDetsCatAddresser(const FullDeterminantCategory& fullDetCat):
      groupOffsets_(fullDetCat.groupOffsets()) {
    
    for (const auto& g : fullDetCat.detGroups()) {
      groupAddressers_.push_back(g.template addresser<DetsT>());
    }
  }

  size_t detsToAddress(const std::vector<DetsT>& dets) const {
    size_t addr = 0ul;
    for (auto i = 0ul; i < dets.size(); i++) {
      addr += groupOffsets_[i] * 
          size_t(groupAddressers_[i].bitStringToAddress(dets[i]));
    }
    return addr;
  }

  void addressToBitStrings(size_t addr, std::vector<DetsT>& detStrings) const {
    for (int i = detStrings.size() - 1; i >= 0; --i) {
      detStrings[i] = groupAddressers_[i].addressToBitString(DetsT(addr / groupOffsets_[i]));
      addr %= groupOffsets_[i]; 
    }
  }
  
  std::string detsToString(const std::vector<DetsT>& dets) const { 
    std::string s = "";
    for (auto i = 0ul; i < dets.size(); ++i) 
       s += groupAddressers_[i].detToString(dets[i]); 
    return s;
  }
  
  std::vector<size_t> addressToOccInfo(size_t addr) const { 
    std::vector<size_t> occInfo;
    std::vector<size_t> detOccInfo;
    std::vector<DetsT> dets(groupAddressers_.size(), DetsT(0));
    addressToBitStrings(addr, dets);
    size_t orbOffset = 0ul;
    for (auto i = 0ul; i < dets.size(); ++i) {
       detOccInfo.clear();
       detOccInfo = groupAddressers_[i].detOccInfo(dets[i]); 
       for (auto occOrb : detOccInfo) {
         occInfo.push_back(occOrb + orbOffset);
       }
       orbOffset += groupAddressers_[i].nOrbitals();
    }
    return occInfo;
  }

  void stringToDets(const std::string& detString, std::vector<DetsT>& dets) const {
    const auto sIt = detString.begin();
    for (int i = 0ul; i < dets.size(); ++i) {
      dets[i] = groupAddressers_[i].stringToDet(sIt);
    }
  }
}; // class FullDetsCatAddresser

template <typename DetsT = uint64_t> 
class FullDetsCatGenerator {

 private:
  FullDetsCatAddresser<DetsT> addresser_;
  size_t nSpaces_;
  size_t nDeterminants_;
  std::vector<DetsT> firstDetOfGroups_;
  std::vector<DetsT> lastDetOfGroups_;

 public:
  FullDetsCatGenerator() = delete;
  FullDetsCatGenerator(const FullDetsCatGenerator&) = default;
  FullDetsCatGenerator(FullDetsCatGenerator&&) = default;
  
  explicit FullDetsCatGenerator(const FullDeterminantCategory& fullDetCat):
      addresser_(fullDetCat.template addresser<DetsT>()) {
    for (const auto& g : fullDetCat.detGroups()) {
      firstDetOfGroups_.push_back(firstDeterminant<DetsT>(g.nElectrons())); 
      lastDetOfGroups_.push_back(lastDeterminant<DetsT>(g.nElectrons(), g.nOrbitals())); 
    }
    nSpaces_ = firstDetOfGroups_.size();
    nDeterminants_ = fullDetCat.nDeterminants();
  }
  
  const FullDetsCatAddresser<DetsT>& addresser() const { return addresser_; }
  const size_t nSpaces() const { return nSpaces_; }

  template <typename Visitor>
  void visitAllDeterminants(Visitor visitor) {
    visitDeterminants(0ul, nDeterminants_, std::forward<Visitor>(visitor));
  }

  template <typename Visitor>
  void visitDeterminants(size_t braAddrStart, size_t braAddrEnd,
                         Visitor visitor) const {
    std::vector<DetsT> braDetStrings(nSpaces_, 0ul);

    addresser_.addressToBitStrings(braAddrStart, braDetStrings);
    
    for (size_t addr = braAddrStart; addr < braAddrEnd; ++addr) {
      
      visitor(addr, braDetStrings);

      // increment
      // find current iterator index
      size_t i = 0ul;
      for (; i < nSpaces_ - 1; i++) {
        if (braDetStrings[i] != lastDetOfGroups_[i]) break;
      }

      braDetStrings[i] = nextLexicographicBitsPermutation(braDetStrings[i]);

      // reset every determinant previous to i
      for (auto j = 0ul; j < i; j++) {
          braDetStrings[j] = firstDetOfGroups_[j];
      }
    }
  }
  

}; // class FullDetsCatGenerator

/*
 * Helper Functions
 */ 
inline std::ostream& FullDeterminantCategory::printAllDeterminants(std::ostream& os) const {
  os << "    * Mapping of the ";
  for (const auto& g : this->detGroups_) {
    os << "CAS(" << g.nElectrons() << "," << g.nOrbitals() << ")-";
  }
  os << ":" << std::endl;
  
  const auto detGen = this->template generator<uint64_t>();
  const auto& addresser = detGen.addresser();
  detGen.visitDeterminants(0ul, nDeterminants_, 
      [&] (size_t addr, const std::vector<uint64_t>& dets) {
        os << "      - Addr " << std::setw(5) << addr
           << ": " << addresser.detsToString(dets)
           << std::endl;
      }
  );
  os << "--------------------" << std::endl;
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const DeterminantCategory& category) {
  return category.printAllDeterminants(os);
}

/*
 * Build a category with user defined space partitioning and occupation numbers
 * DeterminantCategory is a virtual class
 */
inline std::shared_ptr<DeterminantCategory> buildFullDeterminantCategory(
  const std::vector<ActiveSpaceParameters>& activeSpace,
  // refOccupation is the reference orbital occupation
  const std::vector<size_t>& refOccupation) {

  //size_t accOcc = 0ul;

  std::vector<size_t> nOrbs;
  std::vector<int> iDASGroupOcc(refOccupation.size(),0);
  for (auto i = 0ul; i < refOccupation.size(); i++) {
    // If this is part of a combined excitation, the sum of occupations of relevant spaces must be restricted
    // TODO: there can be different restrictions, e.g. up to triplets in certain virtual spaces and up to doubles in higher virtual spaces
    iDASGroupOcc[activeSpace[i].iDASGroup] += refOccupation[i] - activeSpace[i].nElectrons;
    if (refOccupation[i] > activeSpace[i].nOrbitals) return nullptr;
    nOrbs.push_back(activeSpace[i].nOrbitals);
  }

  for (auto i = 0ul; i < refOccupation.size(); i++)
    if ( (activeSpace[i].eLimit != 0 and (iDASGroupOcc[activeSpace[i].iDASGroup] > activeSpace[i].eLimit) )
      or (activeSpace[i].hLimit != 0 and (iDASGroupOcc[activeSpace[i].iDASGroup] < -activeSpace[i].hLimit)) ) {
        return nullptr;
    }

  // FullDetsCategory includes full expansion in each active space, i.e., multiple complete active spaces
  return std::make_shared<FullDeterminantCategory>(refOccupation, nOrbs);
} // buildFullDetsCategory

} // namespace ChronusQ
