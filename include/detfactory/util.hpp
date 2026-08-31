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

#include <bit>
#include <detfactory/bitmasks.hpp>
#include <detfactory/addressingarray.hpp>
#include <cxxapi/input.hpp>

namespace ChronusQ {
  
template <typename DetsT>
bool checkBit(DetsT d, size_t i) { 
  return d & BITMASKS<DetsT>::bitAt[i]; 
}

template <typename DetsT>
DetsT flipBit(DetsT d, size_t i) { 
  return d ^ BITMASKS<DetsT>::bitAt[i]; 
}

template <typename DetsT>
DetsT createDeterminant(const std::vector<DetsT>& elecPos) {
  DetsT d = 0;
  for (const auto & i: elecPos) d = flipBit(d, i);
  return d; 
}

template <typename DetsT>
DetsT firstDeterminant(size_t nE) { 
  return BITMASKS<DetsT>::bitsPrevTo[nE]; 
}

template <typename DetsT>
DetsT lastDeterminant(size_t nE, size_t nO) {
  return BITMASKS<DetsT>::bitsAfter[nO - nE] & 
    BITMASKS<DetsT>::bitsPrevTo[nO];
} 
  
// helper functions for debug and analysis determinants
template <typename DetsT, std::size_t N = sizeof(DetsT) * 8>
std::string determinantToString(DetsT d, size_t nOrb) {
  assert(nOrb <= N);
  std::bitset<N> bits(d);
  
  std::string s = bits.to_string();
  s = s.substr(N-nOrb, nOrb); 
  std::reverse(s.begin(), s.end());

  return s;
};

template <typename DetsT>
std::string determinantToString(const std::vector<DetsT>& ds, const std::vector<size_t>& nOrbs) {
  assert(ds.size() == nOrbs.size());
  std::string s = ""; 
  for (auto i = 0; i < ds.size(); i++)
    s += determinantToString(ds[i], nOrbs[i]) + " ";
  return s;
} // determinantToString(ds, nOrbs)
  
template <typename DetsT, typename StringIt>
DetsT stringToDeterminant(StringIt sIt, size_t nOrb) {
  DetsT d = 0ul;
  for (auto i = 0ul; i < nOrb; ++sIt, ++i) {
    if (*sIt != '0')  d = flipBit(d, i);
  }
  return d;
} // stringToDeterminant

template <typename DetsT>
DetsT stringToDeterminant(const std::string& s) {
  return stringToDeterminant<DetsT>(s.begin(), std::min(8ul * sizeof(DetsT), s.size())); 
} // stringToDeterminant

template <typename DetsT>
std::vector<DetsT> stringToDeterminant(const std::string& s, const std::vector<size_t>& nOrbs) {
  
  assert(s.size() == std::accumulate(nOrbs.begin(), nOrbs.end(), 0));
  std::vector<DetsT> ds;

  for (auto i = 0, j = 0; i < nOrbs.size(); i++) {
    ds.push_back(stringToDeterminant<DetsT>(s.substr(j, j + nOrbs[i])));
    j += nOrbs[i];
  }
    
  return ds; 
} // stringToDeterminant

template <bool occType, typename DetsT> 
void unsafeOccupationInfo(DetsT d, size_t nOrb, size_t * pos){
  
  size_t iE = 0;
  for (auto iO = 0ul; iO < nOrb; iO++) {
    if (checkBit(d, iO) == occType) {
      pos[iE] = iO;
      iE++;
    }
  }
  return; 
}

template <typename DetsT>
std::vector<size_t> occupiedInfo(DetsT d, size_t nOrb, size_t nElec) {
  assert(nOrb <= sizeof(DetsT) * 8);
  std::vector<size_t> occupiedPos(nElec, 0); 
  unsafeOccupationInfo<bool(1)>(d, nOrb, occupiedPos.data());
  
  return occupiedPos;
}

template <typename DetsT>
std::vector<size_t> virtualInfo(DetsT d, size_t nOrb, size_t nElec) {
  
  assert(nOrb <= sizeof(DetsT) * 8);
  std::vector<size_t> virtualPos(nElec, 0); 
  unsafeOccupationInfo<bool(0)>(d, nOrb, virtualPos.data());
  
  return virtualPos;
}

template <typename DetsT>
DetsT negativeSign1eExcitation(DetsT d, size_t p, size_t q) {
  if (p > q) return negativeSign1eExcitation(d, q, p);

  // count number of "1" from p+1 to q-1 bit
  return std::popcount(static_cast<DetsT>((d & BITMASKS<DetsT>::bitsAfter[p+1]) & BITMASKS<DetsT>::bitsPrevTo[q])) % 2;
}

// assuming dp dq order
template <typename DetsT>
DetsT negativeSign1eExcitation(DetsT dp, DetsT dq, size_t p, size_t q) {
  return (std::popcount(static_cast<DetsT>(dp & BITMASKS<DetsT>::bitsAfter[p+1])) +
                std::popcount(static_cast<DetsT>(dq & BITMASKS<DetsT>::bitsPrevTo[q]))) % 2 ;
}
   
template <typename DetsT>
DetsT nextLexicographicBitsPermutation(DetsT v) {

  // https://graphics.stanford.edu/~seander/bithacks.html
  // t gets v's least significant 0 bits set to 1
  // Next set to 1 the most significant bit to change, 
  // set to 0 the least significant ones, and add the necessary 1 bits.
  
  DetsT t = v | (v - 1);
  return (t + 1) | (((~t & -~t) - 1) >> (std::countr_zero(v) + 1));
}
  
template <typename DetsT>
DetsT detToLexicographicAddr(const std::vector<size_t> & elecPos) {
  DetsT addr = 0;
  for (auto iE = 0ul; iE < elecPos.size(); iE++)
    addr += BITADDRESSING<DetsT>::array[iE][elecPos[iE]];
  return addr;
}

template <typename DetsT>
DetsT detToLexicographicAddr(DetsT d, size_t nE, size_t nO) {
  if (d == 0 or d == std::numeric_limits<DetsT>::max()) return DetsT(0);
  DetsT addr = 0;
  for (auto iE = 0ul, iO = 0ul; iO < nO and iE < nE; iO++) {
    if (checkBit(d, iO)) {
      addr += BITADDRESSING<DetsT>::array[iE][iO];
      iE++;
    }   
  }
  return addr; 
}

template <typename DetsT>
DetsT lexicographicAddrToBitString(DetsT addr, size_t nE, size_t nO) {
  // address 0 is associated with the bit-string where the lowest nE orbitals are occupied.
  if (addr == 0 or nE == nO) return firstDeterminant<DetsT>(nE);
  DetsT d = 0;
  size_t orb_upper = nO - 1;
  for (int iE = nE - 1; iE >= 0; iE--) 
  for (int iO = orb_upper; iO >= iE; iO--) {
    if (addr >= BITADDRESSING<DetsT>::array[iE][iO]) {
      d = flipBit(d, iO);
      addr -= BITADDRESSING<DetsT>::array[iE][iO];
      orb_upper = iO - 1;
      break;
    }
  }
  return d;
}

/*
 * Helper function to find occ  
 */ 
template <typename DetsT>
void determinantsToOccs(const std::vector<DetsT>& dets, 
    const std::vector<size_t>& nEs, const std::vector<size_t>& nOs,
    std::vector<size_t>& occ) {
  size_t * occ_it = occ.data(); 
  size_t nOOff = 0ul;
  for (auto i = 0ul; i < dets.size(); ++i) {
    unsafeOccupationInfo<bool(1)>(dets[i], nOs[i], occ_it);
    for (auto j = 0ul; j < nEs[i]; ++j, ++occ_it) {
      *occ_it += nOOff; 
    }    
    nOOff += nOs[i];
  }
} // determinantsToOccs

/*
 * parse g2e(1,0,2,4)-X
 *       h1e(1,0)-...
 *
 * to an array with numbers and extraTerms
 */ 
inline void parseTermSpan(const std::string& term,
  std::vector<std::string>& individualTerms,
  std::vector<size_t>& firstTermSpans, 
  const std::string& delimiters = "-") {
  
  split(individualTerms, term, delimiters);
  
  // get first term span 
  std::vector<std::string> span_tokens;
  split(span_tokens, individualTerms[0].substr(3), "(, )");
  
  firstTermSpans.clear();
  for (auto const & t: span_tokens) firstTermSpans.push_back(std::stoul(t)); 
  
} // parseTermSpan
  
}; // namespace ChronusQ

