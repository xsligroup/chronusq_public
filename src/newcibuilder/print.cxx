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

#include <detfactory/print.hpp>

/*
 *  Diagnostic Printing Utility Functions 
 */

namespace ChronusQ {

#ifdef EXCITATIONLIST_DIAGNOSIS
bool FullCD1eExList::runDiagnosis(std::ostream& os) const {
  const auto failedTests = examEntries(); 
  if (failedTests.empty()) {
    os << "==========> This 1e EX list has passed ALL tests" << std::endl;
    return true;
  } else {
    os << "==========> Wrong Entries in this 1e EX list:" << std::endl;
    auto count = 1ul;
    for (const auto& T: failedTests) {  
      os  << "   #" << std::setw(4) << count << ":";
      for (const auto& i: T) os << i;
      os  << std::endl;  
    }
  }
  os << "---- catgeory K in 1e ex list" << std::endl << *this->categoryK() << std::endl;
  os << "---- catgeory L in 1e ex list" << std::endl << *this->categoryL() << std::endl;
  return false;
}

template <typename DetsT> 
std::vector<std::array<size_t, 5>> IntraSpaceFullCD1eExList<DetsT>::examEntries() const {
  const auto& group  = this->categoryK()->groups()[0];
  auto addresser = group.template addresser<DetsT>();
  
  auto it = begin();
  auto itEnd = end();
  auto const & K = it.curBatch();
  auto const & [p, q, L, sign] = *it;
  
  std::vector<std::array<size_t, 5>> failedTests{}; 
  auto exGen = this->generator()   

  for (; it != itEnd; ++it) {
    auto DetK = addresser.addressToDet(K);
    auto DetL = addresser.addressToDet(L);
    DetK = flipBit(DetK, p);
    DetK = flipBit(DetK, q);
    if (DetK != DetL or  
        sign != negativeSign1eExcitation(DetL, p, q) ) {
      failedTests.push_back({K, p, q, L, sign});
    }
  }
  return std::move(failedTests);
}

template std::vector<std::array<size_t, 5>> IntraSpaceFullCD1eExList<uint8_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> IntraSpaceFullCD1eExList<uint16_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> IntraSpaceFullCD1eExList<uint32_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> IntraSpaceFullCD1eExList<uint64_t>::examEntries() const;

template <typename DetsT> 
std::vector<std::array<size_t, 5>> InterSpaceFullCD1eExList<DetsT>::examEntries() const {

  auto catK = this->categoryK();
  auto catL = this->categoryL();
  const auto& groupKp = catK->groups()[0];
  const auto& groupKq = catK->groups()[1];
  const auto& groupLp = catL->groups()[0];
  const auto& groupLq = catL->groups()[1];
  
  auto KpAddresser = groupKp.template addresser<DetsT>();
  auto KqAddresser = groupKq.template addresser<DetsT>();
  auto LpAddresser = groupLp.template addresser<DetsT>();
  auto LqAddresser = groupLq.template addresser<DetsT>();


  const auto& KpDim = groupKp.nDeterminants(); 
  const auto& LpDim = groupLp.nDeterminants(); 
  
  auto it = begin();
  auto itEnd = end();
  const auto& K = it.curBatch();
  const auto& [p, q, Lp, Lq, sign] = *it;
  
  const auto checkSign = (groupLp.order() < groupLq.order()) ?
      [](DetsT Dp, DetsT Dq, size_t p, size_t q) {
        return negativeSign1eExcitation(Dp, Dq, p, q); 
      } : 
      [](DetsT Dp, DetsT Dq, size_t p, size_t q) {
        return negativeSign1eExcitation(Dq, Dp, q, p); 
      }; 

  std::vector<std::array<size_t, 5>> failedTests{}; 
  for (; it != itEnd; ++it) {
    auto Kp = K % KpDim;
    auto Kq = K / KpDim;
    auto DetKp = KpAddresser.addressToDet(Kp);
    auto DetKq = KqAddresser.addressToDet(Kq);
    auto DetLp = LpAddresser.addressToDet(Lp);
    auto DetLq = LqAddresser.addressToDet(Lq);
    DetKp = flipBit(DetKp, p);
    DetKq = flipBit(DetKq, q);
    if (DetKp != DetLp or DetKq != DetLq or 
        sign != checkSign(DetLp, DetLq, p, q) ) {
      failedTests.push_back({K, p, q, Lp + Lq * LpDim, sign});
    }
  }
  return std::move(failedTests);
}

template std::vector<std::array<size_t, 5>> InterSpaceFullCD1eExList<uint8_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> InterSpaceFullCD1eExList<uint16_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> InterSpaceFullCD1eExList<uint32_t>::examEntries() const;
template std::vector<std::array<size_t, 5>> InterSpaceFullCD1eExList<uint64_t>::examEntries() const;
#endif

} // namespace ChronusQ
