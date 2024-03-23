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

#include <detfactory/util.hpp>
#include <util/math.hpp>

namespace ChronusQ {
  
// declare addresser and generator 

template <typename DetsT>
class DetGroupAddresser;

template <typename DetsT>
class DetsGroupGenerator;

/*
 * \brief DeterminantGroup Class
 *
 * represents determinants in a Complete orbital space
 * with nOrbitals and nElectrons
 *
 */ 
class DeterminantGroup {

protected:

  size_t nOrbitals_     = 0;
  size_t nElectrons_    = 0;
  size_t nDeterminants_ = 1;

public:
  
  DeterminantGroup() = delete;
  DeterminantGroup(const DeterminantGroup &) = default;
  DeterminantGroup(DeterminantGroup &&) = default;
  ~DeterminantGroup() = default;
  
  DeterminantGroup(size_t nOccupation, size_t nOrbitals):
          nOrbitals_(nOrbitals), nElectrons_(nOccupation) {
    initialization();
  }
  
  explicit DeterminantGroup(const std::array<size_t, 2>& groupInfo):
    DeterminantGroup(groupInfo[0], groupInfo[1]) { }
  
  void initialization() {
    if (nOrbitals_ < nElectrons_) {
      std::string output = "Can't initialize group with nO < nE: nO ="; 
      output += std::to_string(nOrbitals_);
      output += ", nE = ";
      output += std::to_string(nElectrons_);
      throw std::invalid_argument(output);
    }   
    nDeterminants_ = Comb(nOrbitals_, nElectrons_);
  }

  void resetNElectrons(size_t nE) {
    nElectrons_ = nE; 
    initialization();
  }
  
  void addNElectrons(size_t nE) {
    resetNElectrons(nElectrons_ + nE);
  }

  void removeNElectrons(size_t nE) {
    resetNElectrons(nElectrons_ - nE);
  }

  // getters
  size_t nOrbitals() const { return nOrbitals_; }
  size_t nElectrons() const { return nElectrons_; }
  size_t nHoles() const { return nOrbitals_ - nElectrons_; }
  size_t nDeterminants() const { return nDeterminants_; }
  
  size_t nDetsAfterExcitation(int nEx) const { 
    size_t nElec = int(nElectrons_) + nEx;
    return Comb(nOrbitals_, nElec);
  }

  bool operator==(const DeterminantGroup& other) const {
    return nOrbitals()  == other.nOrbitals() and
      nElectrons() == other.nElectrons();
  }
  bool operator!=(const DeterminantGroup& other) const {
    return not operator==(other);
  }
  
  std::ostream& printAllDeterminants(std::ostream& os) const;

  // major interface for addresssing and generating 
  template <typename DetsT>
  DetGroupAddresser<DetsT> addresser() const {
    return DetGroupAddresser<DetsT>(nElectrons_, nOrbitals_);
  }
  
  template <typename DetsT>
  DetsGroupGenerator<DetsT> generator() const {
    return DetsGroupGenerator<DetsT>(nElectrons_, nOrbitals_);
  }

}; // DeterminantGroup

// DetsT = StringT
template <typename DetsT = uint64_t>
class DetGroupAddresser {
 private:
  size_t nE_ = 0ul;
  size_t nO_ = 0ul;
 
 public:
  DetGroupAddresser() = default;
  DetGroupAddresser(const DetGroupAddresser&) = default;
  DetGroupAddresser(DetGroupAddresser&&) = default;

  DetGroupAddresser(size_t nE, size_t nO): nE_(nE), nO_(nO) { }
  
  DetsT bitStringToAddress(DetsT d) const {
    return detToLexicographicAddr(d, nE_, nO_); 
  }

  DetsT addressToBitString(DetsT addr) const {
    return lexicographicAddrToBitString(addr, nE_, nO_);
  }

  std::string detToString(DetsT d) const {
    return determinantToString(d, nO_); 
  }
  
  DetsT stringToDet(const std::string& detString) const { 
    return stringToDet(detString.begin());
  }

  DetsT stringToDet(const std::string::iterator sIt) const {
    return stringToDeterminant<DetsT>(sIt, nO_); 
  }
}; // class DetGroupAddresser

// DetsT = StringT
template <typename DetsT = uint64_t> 
class DetsGroupGenerator {
 
 private:
  DetGroupAddresser<DetsT> addresser_;
  
 public:
  DetsGroupGenerator() = delete;
  DetsGroupGenerator(const DetsGroupGenerator&) = default;
  DetsGroupGenerator(DetsGroupGenerator&&) = default;
  
  DetsGroupGenerator(size_t nE, size_t nO): addresser_(nE, nO) { }
  
  const DetGroupAddresser<DetsT> addresser() const { return addresser_; }

  template <typename Visitor>
  void visitDeterminants(size_t addrStart, size_t addrEnd, 
      Visitor visitor) const { 
    // assert(addrStart <= addrEnd);

    // generate a bit-string unique to the determinant based on the determinant index
    // maybe it should be called indexToBitString
    DetsT detString = addresser_.addressToBitString(addrStart);
    for (size_t addr = addrStart; addr < addrEnd; ++addr) {
      //std::cout << "addr = " << addr << ", string = " << addresser_.detToString(det) << std::endl;
      visitor(addr, detString);
      detString = nextLexicographicBitsPermutation(detString);
    }
  }

}; // class DetsGroupGenerator

/*
 * Helper Functions
 */ 
inline std::ostream& DeterminantGroup::printAllDeterminants(std::ostream& os) const {
  os << "* Mapping of the CAS(" << nElectrons_ << ", " << nOrbitals_ << "):" << std::endl;
  
  const auto detGen = this->template generator<uint64_t>();
  const auto& addresser = detGen.addresser();
  detGen.visitDeterminants(0ul, nDeterminants_, 
      [&] (size_t addr, uint64_t det) {
        os << "      - Addr " << std::setw(5) << addr
           << ": " << addresser.detToString(det)
           << std::endl;
      }
  );
  os << "--------------------" << std::endl;
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const DeterminantGroup& group) {
  return group.printAllDeterminants(os);
}

} // namespace ChronusQ
