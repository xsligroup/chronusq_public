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

// #define DEBUG_INTRA1EEXLIST

namespace ChronusQ {
  
template <typename StorageT>
class IntraSpaceFullCD1eExListGenerator;

/* 
 * Implementation of intra-space full configuration driven
 *   1e excitation list:
 *
 *   as a K list with excitation elements stored as [p, q, <K|E_pq|L>], [L]
 *
 *   The first date type is the bit-string for each determinant.
 *   The second data type is the integer address converted from bit-string.
 */ 
template <typename StringT, typename StorageT>
class IntraSpaceFullCD1eExList: public FullCD1eExList {

 private:
  // There are two connected excitation lists.
  // exList_Ket_ here stores lexicographical addresses of non-Zero Ket determinant for a given Bra determinant.
  // [nNonZeroKetDets, nBraDets]=iNonZeroKetDetAddress
  //
  // exList_pgSign_ in the parent class FullCD1eExList stores the p, q addresses and the sign
  // for a given BraKet determinant pair.
  // [nNonZeroKetDets, nBraDets]={p,q,sign}
  ExcitationListStorage<StorageT> exList_Ket_;

 public:
  
  // Constructors
  IntraSpaceFullCD1eExList() = delete;
  IntraSpaceFullCD1eExList(const IntraSpaceFullCD1eExList & other) = default;
  IntraSpaceFullCD1eExList(IntraSpaceFullCD1eExList && other)      = default;
  ~IntraSpaceFullCD1eExList() = default;  
  
  IntraSpaceFullCD1eExList(const DeterminantGroup& detGroup):
    FullCD1eExList(detGroup.nDeterminants(),
                   detGroup.nElectrons() * (detGroup.nHoles() + 1)), exList_Ket_() {
    
    exList_Ket_.resize(1, this->nNonZeroKetDets_, this->nBraDets_);
    
    auto localCategory = std::make_shared<FullDeterminantCategory>(std::vector<DeterminantGroup>{detGroup});

    // These are local categories that include active spaces of interest with respect
    // to the excitation list.
    this->setKetCategory(localCategory);
    this->setBraCategory(localCategory);
  }
   
  const ExcitationListStorage<StorageT>& exList_Ket() const {
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
    
    const auto& braDetGroup = this->braCategory()->detGroups()[0];
  
    const size_t nOrb = braDetGroup.nOrbitals();
    const size_t nBraDets = this->nBraDets_;
    const size_t nBraDetsPerThread = std::ceil(double(nBraDets) / GetNumThreads());

    // CLion suggests default(none)
    #pragma omp parallel default(shared)
    {
      auto detGenerator = braDetGroup.template generator<StringT>();
      const auto& addresser = detGenerator.addresser(); 
      std::vector<size_t> nBraOcc(braDetGroup.nElectrons());
      std::vector<size_t> nBraVir(braDetGroup.nHoles());
      
      size_t iBraDetBegin = nBraDetsPerThread * GetThreadID();
      size_t iBraDetEnd   = std::min(nBraDets, iBraDetBegin + nBraDetsPerThread);

      detGenerator.visitDeterminants(iBraDetBegin, iBraDetEnd,
                                     [&] (size_t iBraDetAddr, StringT iBraDetString) {
            auto exList_Ket_ptr = exList_Ket_.pointer(0ul, 0ul, iBraDetAddr);
            auto exList_pqSign_ptr = this->exList_pqSign_.pointer(0ul, 0ul, iBraDetAddr);

            // get the occupied bits in q space and empty bits in p space for iBraDetString
            unsafeOccupationInfo<bool(1)>(iBraDetString, nOrb, nBraOcc.data());
            unsafeOccupationInfo<bool(0)>(iBraDetString, nOrb, nBraVir.data());
      
            /*
             * Configuration Driven 1e Excitation List: <iBraDetString|a_p^\dagger a_q |iKetDetString>
             *   (1) nonzero excitation p,
             *   (2) nonzero excitation q,
             *   (3) address of excited string |iKetDetString> in ketCategory
             *   (4) whether it's negative sign
             */
      
             // Case 1: self-excitation
             for(const auto& p: nBraOcc)  {
               exList_pqSign_ptr[0] = p;
               exList_pqSign_ptr[1] = p;
               exList_pqSign_ptr[2] = 0;
               exList_pqSign_ptr += 3;
               *exList_Ket_ptr = static_cast<StorageT>(iBraDetAddr);
               exList_Ket_ptr++;
             }

             // Case 2: p -> q, p=occ, q=virtual in iBraDetString
             StringT ketDetString;
             for(const auto& p: nBraOcc)
             for(const auto& q: nBraVir) {
               ketDetString = flipBit(iBraDetString, p);
               ketDetString = flipBit(ketDetString, q);
               exList_pqSign_ptr[0] = p;
               exList_pqSign_ptr[1] = q;
               exList_pqSign_ptr[2] = negativeSign1eExcitation(ketDetString, p, q);
               exList_pqSign_ptr += 3;
               *exList_Ket_ptr = static_cast<StorageT>(addresser.bitStringToAddress(ketDetString));
               exList_Ket_ptr++;
             }
           }
      ); // visitDeterminants
    }

  } // computeExcitationList

  using FullCD1eExList::generator;
  std::shared_ptr<FullCD1eExListGenerator> generator() const override {
    return generator({1ul, 0ul}, {1ul, 0ul}); 
  }

  std::shared_ptr<FullCD1eExListGenerator> generator(
      const std::pair<size_t, size_t>& ketExOffs,
      const std::pair<size_t, size_t>& braExOffs) const override {
    return std::make_shared<IntraSpaceFullCD1eExListGenerator<StorageT>>(
      this->nNonZeroKetDets_, this->nBraDets_, ketExOffs, braExOffs,
      exList_Ket(), this->exList_pqSign());
  }

}; // class IntraSpaceFullCD1eExList 

template <typename StorageT>
class IntraSpaceFullCD1eExListGenerator: public FullCD1eExListGenerator {
 protected:
  const ExcitationListStorage<StorageT>& exList_Ket_;
   
  void buildExcitations_(size_t K, size_t extraLOff) override {
    FullCD1eExListGenerator::buildExcitations_(K, extraLOff);
    // advance the address of Bra determinant in the category
    this->currentBraEx_ = K * braExOffs_.first;
    const StorageT* exList_Ket_ptr = exList_Ket_.pointer(0ul, 0ul, K);
    const size_t LOff = this->ketExOffs_.first;
    for (auto& ex : this->excitations_) {
      ex[2] += (*exList_Ket_ptr) * LOff;
      exList_Ket_ptr++;
    }
  }

 public:
  explicit IntraSpaceFullCD1eExListGenerator(
    size_t nNonZeroKetDets, size_t nBraDets,
    const std::pair<size_t, size_t>& ketExOffs,
    const std::pair<size_t, size_t>& braExOffs,
    const ExcitationListStorage<StorageT>& exList_Ket,
    const ExcitationListStorage<uint8_t>& exList_pqSign):
    FullCD1eExListGenerator(nNonZeroKetDets, nBraDets, ketExOffs, braExOffs, exList_pqSign),
    exList_Ket_(exList_Ket) { }

}; // class IntraSpaceFullCD1eExListGenerator

/*
 *  Utility function to build IntraSpaceFullCD1eExList
 */
template <typename DetStringT>
std::shared_ptr<NewExcitationList>
constructIntraSpaceFullCD1eExList(const DeterminantGroup& detGroup) {
 
  size_t nDets = detGroup.nDeterminants();

  // The second data type is the integer address converted from bit-string.
  if (nDets <= std::numeric_limits<uint8_t>::max()) {
    return std::make_shared<IntraSpaceFullCD1eExList<DetStringT, uint8_t>>(detGroup);
  } else if (nDets <= std::numeric_limits<uint16_t>::max()) {
    return std::make_shared<IntraSpaceFullCD1eExList<DetStringT, uint16_t>>(detGroup);
  } else if (nDets <= std::numeric_limits<uint32_t>::max()) {
    return std::make_shared<IntraSpaceFullCD1eExList<DetStringT, uint32_t>>(detGroup);
  }
  
  return std::make_shared<IntraSpaceFullCD1eExList<DetStringT, uint64_t>>(detGroup);
}

inline std::shared_ptr<NewExcitationList> 
constructIntraSpaceFullCD1eExList(const DeterminantGroup& detGroup) {
  
  size_t nOrbs = detGroup.nOrbitals();

  // Figure out the minimal number of bits needed for describing determinant bit-strings
  // In the new framework, this may not be necessary since bit-strings are never saved.
  // It will be converted and saved as an integer address.
  if (nOrbs <= 8) {
    return constructIntraSpaceFullCD1eExList<uint8_t>(detGroup);
  } else if (nOrbs <= 16) {
    return constructIntraSpaceFullCD1eExList<uint16_t>(detGroup);
  } else if (nOrbs <= 32) {
    return constructIntraSpaceFullCD1eExList<uint32_t>(detGroup);
  } else if (nOrbs <= 64) {
    return constructIntraSpaceFullCD1eExList<uint64_t>(detGroup);
  }
  CErr("number of orbitals in DeterminantGroup is too large to be represented in fundamental data types");
  return nullptr;
}

} // namespace ChronusQ
