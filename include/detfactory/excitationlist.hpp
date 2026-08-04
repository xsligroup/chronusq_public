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

#define DETFACTORY_EXCITATIONLIST

#include <detfactory/determinants/detscategory.hpp>
#include <detfactory/excitationlist/storage.hpp>

namespace ChronusQ {
  
class NewExcitationList { 
 protected:
  // excitation: Cat L -> Cat K, <K|E|L> 
  std::shared_ptr<DeterminantCategory> ketCategory_ = nullptr;
  std::shared_ptr<DeterminantCategory> braCategory_ = nullptr;
  
 public:
  NewExcitationList() = default;
  NewExcitationList(const NewExcitationList&) = default;
  NewExcitationList(NewExcitationList&&) = default;
  ~NewExcitationList() = default;
  
  // getters
  std::shared_ptr<const DeterminantCategory> ketCategory() const { return ketCategory_; }
  std::shared_ptr<const DeterminantCategory> braCategory() const { return braCategory_; }
  
  // setters
  void setKetCategory(const std::shared_ptr<DeterminantCategory>& catL) { ketCategory_ = catL; }
  void setBraCategory(const std::shared_ptr<DeterminantCategory>& catK) { braCategory_ = catK; }

  // virtual computing and other interface
  virtual void computeExcitationList() = 0;
  virtual void output(std::ostream& os, 
      const std::string& name = "Excitation List") const {
    os << name << ": (";
    for (const auto & j : braCategory_->SpaceOccupations()) {
      os << std::setw(3) << j;
    }
    os << ") <- (";
    for (const auto & j : ketCategory_->SpaceOccupations()) {
      os << std::setw(3) << j;
    }
    os << ")" << std::endl;
    return;
  }
  
  virtual double storageSize() const = 0;

}; // NewExcitationList 

// forward declaration
class FullCD1eExListGenerator;

/*
 * \brief the FullCD1eExList class
 * 
 *  for configuration driven full 1e excitation list
 *
 */  
class FullCD1eExList: public NewExcitationList {
 protected:
  size_t nBraDets_;     // number of determinants in the Bra space
  size_t nNonZeroKetDets_;  // number of non-zero Ket determinant per Bra determinant

  // exList_pgSign_ in here stores the p, q addresses and the sign for a given non-Zero BraKet determinant pair.
  // [nNonZeroKetDets, nBraDets]={p,q,sign}
  // {p, q, sign} are stored as uint8
  ExcitationListStorage<uint8_t> exList_pqSign_;   
 
 public:
  FullCD1eExList() = delete;
  FullCD1eExList(size_t nBraDets, size_t nNonZeroKetDets):
    nBraDets_(nBraDets), nNonZeroKetDets_(nNonZeroKetDets),
    exList_pqSign_(3ul, nNonZeroKetDets, nBraDets) { };
  FullCD1eExList(const FullCD1eExList &) = default;
  FullCD1eExList(FullCD1eExList &&) = default;
  ~FullCD1eExList() = default; 

  virtual void tryAlloc() {
    exList_pqSign_.tryAlloc();
  }

  virtual void dealloc() {
    exList_pqSign_.dealloc();
  }
  
  size_t nNonZeroExcitations() const { return nNonZeroKetDets_; }
  size_t nBraDeterminants() const { return nBraDets_; }

  const ExcitationListStorage<uint8_t>& exList_pqSign() const { 
    return exList_pqSign_; 
  }
  
  virtual double storageSize() const override {
    return sizeof(uint8_t) * exList_pqSign_.totalDimension(); 
  }
  
  virtual void computeExcitationList() override = 0;

  // virtual generator interface
  virtual std::shared_ptr<FullCD1eExListGenerator> generator() const = 0;
  
  virtual std::shared_ptr<FullCD1eExListGenerator> generator(
      const std::pair<size_t, size_t>& LExOffs,
      const std::pair<size_t, size_t>& KExOffs) const = 0;


}; // FullCD1eExList

inline std::ostream& operator<<(std::ostream& os, const FullCD1eExList& exList) {
  exList.output(os, "1e Excitation List");
  return os;
}

class FullCD1eExListGenerator {
 protected:
  const ExcitationListStorage<uint8_t>& exList_pqSign_;   
  
  // as {p, q, ketEx, sign} * nNonZeroKetDets_;
  // this is used as size_t instead of unit8_t storage type
  std::vector<std::array<size_t, 4>> excitations_; 
  
  size_t curK_;
  size_t curExtraLOff_ = 0ul;
   
  size_t nBraDets_ = 0ul;

  size_t currentBraEx_;
  std::pair<size_t, size_t> braExOffs_;
  std::pair<size_t, size_t> ketExOffs_;
  
  virtual void buildExcitations_(size_t K, size_t extraLOff) {
    const uint8_t* exList_pqSign_ptr = exList_pqSign_.pointer(0ul, 0ul, K);
    for (auto& ex : excitations_) {
      ex[0] = exList_pqSign_ptr[0]; // p
      ex[1] = exList_pqSign_ptr[1]; // q
      ex[2] = extraLOff;
      ex[3] = exList_pqSign_ptr[2]; // sign
      exList_pqSign_ptr += 3;
    }
  };
 
 public:
  FullCD1eExListGenerator() = delete;
  FullCD1eExListGenerator(const FullCD1eExListGenerator&) = delete;
  FullCD1eExListGenerator(FullCD1eExListGenerator&&) = delete;
  
  explicit FullCD1eExListGenerator(size_t nNonZeroKetDets, size_t nBraDets,
      const std::pair<size_t, size_t>& ketExOffs,
      const std::pair<size_t, size_t>& braExOffs,
      const ExcitationListStorage<uint8_t>& exList_pqSign):
    nBraDets_(nBraDets), ketExOffs_(ketExOffs), braExOffs_(braExOffs),
    exList_pqSign_(exList_pqSign) {
    curK_ = nBraDets;
    excitations_.resize(nNonZeroKetDets);
  }

  const std::vector<std::array<size_t, 4>>& excitations() const { 
      return excitations_; 
  }
  
  const size_t& braExAddress() const { return currentBraEx_; }

  void updateExcitations(size_t K, size_t extraLOff) {
    if ( K != curK_) {
      curK_ = K;
      buildExcitations_(K, extraLOff);
    } else {
      for (auto& ex : excitations_) {
        ex[2] += extraLOff;
        ex[2] -= curExtraLOff_;
      }
    }
    curExtraLOff_ = extraLOff;
    return;
  }
  
  template <typename Visitor>
  void visitAllExcitations(Visitor visitor) {
    visitExcitations(0ul, nBraDets_, std::forward<Visitor>(visitor));
  }

  template <typename Visitor>
  void visitExcitations(size_t braBegin, size_t braEnd, Visitor visitor) {
    for (size_t K = braBegin; K < braEnd; ++K) {
      buildExcitations_(K, 0ul);
      visitor(K, excitations_);
    }
  }

#ifdef CQ_ENABLE_SPARSE
  template <typename Visitor>
  void visitAllSparseExcitations(Visitor visitor) {
    visitSparseExcitations(0ul, nBraDets_, std::forward<Visitor>(visitor));
  }

 template <typename Visitor>
  void visitSparseExcitations(size_t braBegin, size_t braEnd, Visitor visitor) {
    for (size_t K = braBegin; K < braEnd; ++K) {
      buildExcitations_(K, 0ul);
      visitor(K, excitations_);
    }
  }
#endif

}; // class FullCD1eExListGenerator

} // namespace ChronusQ

// include implementations
#include <detfactory/excitationlist/intraspacefullcd1eexlist.hpp>
#include <detfactory/excitationlist/interspacefullcd1eexlist.hpp>
#include <detfactory/excitationlist/doublefullcd1eexlistgenerator.hpp>
