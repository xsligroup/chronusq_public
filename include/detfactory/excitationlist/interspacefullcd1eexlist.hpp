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

#ifndef DETFACTORY_EXCITATIONLIST
#error This file may only be included from detfactory/excitationlist.hpp
#endif

#include <util/threads.hpp>

// #define DEBUG_INTER1EEXLIST

namespace ChronusQ {
  
template <typename StorageT>
class InterSpaceFullCD1eExListGenerator;

/* 
 * \brief the InterSpaceFullCD1eExList Class 
 * 
 * as a Kp, Kq excitation as [p, q, <K|E_pq|L>], [Lp, Lq]
 *
 */
template <typename StringT, typename StorageT>
class InterSpaceFullCD1eExList: public FullCD1eExList {
 private:

  // There are two connected excitation lists.
  // exList_Ket_ here stores lexicographical addresses of paired non-Zero Ket determinant for a given Bra determinant.
  // [nNonZeroKetDets, nBraDets]={pNonZeroKetDetAddress, qNonZeroKetDetAddress}
  // Note that the nBraDets = nBraDetsInX + nBraDetsInX' and pNonZeroKetDetAddress, qNonZeroKetDetAddress are
  // paired Ket determinant addresses relative to their own active space X and X' where p and q reside, respectively.
  //
  // exList_pgSign_ in the parent class FullCD1eExList stores the p, q addresses and the sign
  // for a given BraKet determinant pair.
  // [nNonZeroKetDets, nBraDets]={p,q,sign}
  ExcitationListStorage<StorageT> exList_Ket_;
  size_t ntBraDets_;
  size_t ntKetDets_;
  bool reversedSpaceOrder_; 

 public:
  // Constructors
  InterSpaceFullCD1eExList() = delete;
  InterSpaceFullCD1eExList(const InterSpaceFullCD1eExList & other) = default;
  InterSpaceFullCD1eExList(InterSpaceFullCD1eExList && other)      = default;
  ~InterSpaceFullCD1eExList() = default;  
    
  InterSpaceFullCD1eExList(const DeterminantGroup& tBraGroup, const DeterminantGroup& uBraGroup,
                           bool reversedSpaceOrder): reversedSpaceOrder_(reversedSpaceOrder),
                                                     FullCD1eExList(tBraGroup.nDeterminants() * uBraGroup.nDeterminants(),
                                                                    tBraGroup.nElectrons() * uBraGroup.nHoles()), exList_Ket_() {
    
    exList_Ket_.resize(2, this->nNonZeroKetDets_, this->nBraDets_);

    // These are local categories that include active spaces of interest with respect
    // to the excitation list.
    auto localBraCategory = std::make_shared<FullDeterminantCategory>(std::vector<DeterminantGroup>{tBraGroup, uBraGroup});
    this->setBraCategory(localBraCategory);

    // These are local categories that include active spaces of interest with respect
    // to the excitation list.
    DeterminantGroup tKetGroup(tBraGroup), uKetGroup(uBraGroup);
    tKetGroup.removeNElectrons(1);
    uKetGroup.addNElectrons(1);
    auto localKetCategory = std::make_shared<FullDeterminantCategory>(
      std::vector<DeterminantGroup>{tKetGroup, uKetGroup});
    this->setKetCategory(localKetCategory);

    ntBraDets_ = tBraGroup.nDeterminants();
    ntKetDets_ = tKetGroup.nDeterminants();
  }
    
  const ExcitationListStorage<StorageT>& exList_L() const { 
    return exList_Ket_;
  }
  
  void tryAlloc() override {
    FullCD1eExList::tryAlloc();
    exList_Ket_.tryAlloc();
  }

  void dealloc() override {
    FullCD1eExList::dealloc();
    exList_Ket_.dealloc();
  }
  
  double storageSize() const override {
    return sizeof(StorageT) * exList_Ket_.totalDimension() +
           FullCD1eExList::storageSize();
  }
  
  // override virtual methods to populate excitations
  void computeExcitationList() override { 
    
    tryAlloc();

    const auto braCategory = std::dynamic_pointer_cast<const FullDeterminantCategory>(this->braCategory());
    const auto& braDetGroups = braCategory->detGroups();
    const auto& pBraGroup = braDetGroups[0];
    const auto& qBraGroup = braDetGroups[1];

    const auto ketCategory = this->ketCategory();
    const auto& ketDetGroups = ketCategory->detGroups();
    const auto& pKetGroup = ketDetGroups[0];
    const auto& qKetGroup = ketDetGroups[1];

    const auto checkSign = reversedSpaceOrder_ ?
        [](StringT Dp, StringT Dq, size_t p, size_t q) {
          return negativeSign1eExcitation(Dq, Dp, q, p); 
        } : 
        [](StringT Dp, StringT Dq, size_t p, size_t q) {
          return negativeSign1eExcitation(Dp, Dq, p, q); 
        }; 

    const size_t pSpaceNO = pBraGroup.nOrbitals();
    const size_t qSpaceNO = qBraGroup.nOrbitals();
    
    const size_t nBraDets = this->nBraDets_;
    const size_t nBraDetsPerThread = std::ceil(double(nBraDets) / GetNumThreads());
    
    #pragma omp parallel default(shared)
    {
      std::vector<size_t> npBraOcc(pBraGroup.nElectrons());
      std::vector<size_t> nqBraVir(qBraGroup.nHoles());
      const auto pKetAddresser = pKetGroup.template addresser<StringT>();
      const auto qKetAddresser = qKetGroup.template addresser<StringT>();
      auto detCatGen = braCategory->template generator<StringT>();
      
      size_t iBraDetBegin = nBraDetsPerThread * GetThreadID();
      size_t iBraDetEnd   = std::min(nBraDets, iBraDetBegin + nBraDetsPerThread);
      
      detCatGen.visitDeterminants(iBraDetBegin, iBraDetEnd,
                                  [&] (size_t iBraDetAddr, const std::vector<StringT>& iBraDetStrings) {
            auto exList_Ket_ptr = exList_Ket_.pointer(0ul, 0ul, iBraDetAddr);
            auto exList_pqSign_ptr = this->exList_pqSign_.pointer(0ul, 0ul, iBraDetAddr);
            StringT pBraDetString = iBraDetStrings[0], qBraDetString = iBraDetStrings[1], pKetDetString, qKetDetString;
  
            unsafeOccupationInfo<bool(1)>(pBraDetString, pSpaceNO, npBraOcc.data());
            unsafeOccupationInfo<bool(0)>(qBraDetString, qSpaceNO, nqBraVir.data());
      
            // p -> q, p=occ, q=virtual in K
            for(const auto& p: npBraOcc)
            for(const auto& q: nqBraVir) {
              pKetDetString = flipBit(pBraDetString, p);
              qKetDetString = flipBit(qBraDetString, q);
              exList_pqSign_ptr[0] = p;
              exList_pqSign_ptr[1] = q;
              exList_pqSign_ptr[2] = checkSign(pKetDetString, qKetDetString, p, q);
              exList_pqSign_ptr += 3;
              exList_Ket_ptr[0] = static_cast<StorageT>(pKetAddresser.bitStringToAddress(pKetDetString));
              exList_Ket_ptr[1] = static_cast<StorageT>(qKetAddresser.bitStringToAddress(qKetDetString));
              exList_Ket_ptr += 2;
            }
          }
      );
    }
  } // computeExcitationList

  using FullCD1eExList::generator;
  std::shared_ptr<FullCD1eExListGenerator> generator() const override {
    return generator({1ul, ntKetDets_}, {1ul, ntBraDets_});
  }

  std::shared_ptr<FullCD1eExListGenerator> generator(
      const std::pair<size_t, size_t>& LExOffs,
      const std::pair<size_t, size_t>& KExOffs) const override {
    return std::make_shared<InterSpaceFullCD1eExListGenerator<StorageT>>(
      this->nNonZeroKetDets_, this->nBraDets_, ntBraDets_, LExOffs, KExOffs,
      exList_L(), this->exList_pqSign());
  }

}; // class InterSpaceFullCD1eExList

template <typename StorageT>
class InterSpaceFullCD1eExListGenerator: public FullCD1eExListGenerator {
 private:
  size_t KpDim_;

 protected:
  const ExcitationListStorage<StorageT>& exList_L_;   
   
  void buildExcitations_(size_t K, size_t extraLOff) override {
    FullCD1eExListGenerator::buildExcitations_(K, extraLOff);
    this->currentBraEx_ = (K % KpDim_) * this->braExOffs_.first
                          + (K / KpDim_) * this->braExOffs_.second;
    const StorageT* exList_L_ptr = exList_L_.pointer(0ul, 0ul, K);
    const size_t LpOff = this->ketExOffs_.first;
    const size_t LqOff = this->ketExOffs_.second;
    for (auto& ex : this->excitations_) {
      ex[2] += exList_L_ptr[0] * LpOff + exList_L_ptr[1] * LqOff;
      exList_L_ptr += 2;
    }
  }

 public:
  explicit InterSpaceFullCD1eExListGenerator(
      size_t nNonZero, size_t nKDets, size_t nKpDets,
      const std::pair<size_t, size_t>& LExOffs,
      const std::pair<size_t, size_t>& KExOffs,
      const ExcitationListStorage<StorageT>& exList_L,
      const ExcitationListStorage<uint8_t>& exList_pqSign):
      FullCD1eExListGenerator(nNonZero, nKDets, LExOffs, KExOffs, exList_pqSign),
      exList_L_(exList_L), KpDim_(nKpDets) { }

}; // InterSpaceFullCD1eExListGenerator

/*
 *  Utility function to build InterSpaceFullCD1eExList
 */
template <typename DetsT>
std::shared_ptr<NewExcitationList>
constructInterSpaceFullCD1eExList(const DeterminantGroup& tBraGroup,
                                  const DeterminantGroup& uBraGroup, const bool reversedOrder) {
 
  size_t nDets = std::max(tBraGroup.nDetsAfterExcitation(-1),
                          uBraGroup.nDetsAfterExcitation(1));

  // The second data type is the integer address converted from bit-string.
  if (nDets <= std::numeric_limits<uint8_t>::max()) {
    return std::make_shared<InterSpaceFullCD1eExList<DetsT, uint8_t>>(
      tBraGroup, uBraGroup, reversedOrder);
  } else if (nDets <= std::numeric_limits<uint16_t>::max()) {
    return std::make_shared<InterSpaceFullCD1eExList<DetsT, uint16_t>>(
      tBraGroup, uBraGroup, reversedOrder);
  } else if (nDets <= std::numeric_limits<uint32_t>::max()) {
    return std::make_shared<InterSpaceFullCD1eExList<DetsT, uint32_t>>(
      tBraGroup, uBraGroup, reversedOrder);
  }
  
  return std::make_shared<InterSpaceFullCD1eExList<DetsT, uint64_t>>(
    tBraGroup, uBraGroup, reversedOrder);
}

inline std::shared_ptr<NewExcitationList> 
 constructInterSpaceFullCD1eExList(const DeterminantGroup& tBraGroup,
                                  const DeterminantGroup& uBraGroup, const bool reversedOrder) {
  
   size_t nOrbs = std::max(tBraGroup.nOrbitals(), uBraGroup.nOrbitals());

  // Figure out the minimal number of bits needed for describing determinant bit-strings
  // In the new framework, this may not be necessary since bit-strings are never saved.
  // It will be converted and saved as an integer address.
  if (nOrbs <= 8) {
    return constructInterSpaceFullCD1eExList<uint8_t>(tBraGroup, uBraGroup, reversedOrder);
  } else if (nOrbs <= 16) {
    return constructInterSpaceFullCD1eExList<uint16_t>(tBraGroup, uBraGroup, reversedOrder);
  } else if (nOrbs <= 32) {
    return constructInterSpaceFullCD1eExList<uint32_t>(tBraGroup, uBraGroup, reversedOrder);
  } else if (nOrbs <= 64) {
    return constructInterSpaceFullCD1eExList<uint64_t>(tBraGroup, uBraGroup, reversedOrder);
  }
  CErr("number of orbitals in DeterminantGroup is too large to be represented in fundamental data types");
  return nullptr;
}

} // namespace ChronusQ
