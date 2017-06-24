/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
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

#include <mointstransformer.hpp>

namespace ChronusQ {

  /**
   *  \brief set up MO ranges 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::addMORanges(const std::set<char> & symSet,
    const std::pair<size_t, size_t> & range) {
      
      symbol_sets_.push_back(symSet);
      mo_ranges_.push_back(range);

  }; // MOIntsTransformer::addMORanges  
  
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::setMORanges(size_t nFrozenCore, size_t nFrozenVirt) {
    
      // set 4C for no-pair approximation
      size_t fourCompOffset = (ss_.nC == 4) ? ss_.nAlphaOrbital() * 2: 0;
      size_t oneCompFactor  = (ss_.nC == 1) ? 2: 1;

      size_t nO = ss_.nO/oneCompFactor - nFrozenCore;
      size_t nV = ss_.nV/oneCompFactor - nFrozenVirt; 
      size_t nT = nO + nV;
      
      resetMORanges();

      // general indices
      addMORanges({'p','q','r','s'}, {fourCompOffset + nFrozenCore, nT});
      
      // hole indices
      addMORanges({'i','j','k','l'}, {fourCompOffset + nFrozenCore, nO});
      
      // particle indices
      addMORanges({'a','b','c','d'}, {fourCompOffset + ss_.nO, nV}); 
      
  }; // MOIntsTransformer::setMORanges for SingleSlater
  
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::setMORanges(const MCWaveFunctionBase & mcwfn) {
  
      // set 4C for no-pair approximation
      size_t offset = (ss_.nC == 4 and mcwfn.FourCompNoPair) ? ss_.nAlphaOrbital() * 2: 0;
      
      const auto &  mopart = mcwfn.MOPartition;
      offset += mopart.nFCore;
      size_t nInact = mopart.nInact; 
      size_t nCorrO = mopart.nCorrO;
      size_t nFVirt = mopart.nFVirt; 
      size_t nT = mopart.nMO; 
      size_t corrOffset = offset + nInact;

      resetMORanges();
      
      // negative energy MO indices
      if (ss_.nC == 4) addMORanges({'m', 'n'}, {0ul, mopart.nNegMO});

      // general electronic indices
      addMORanges({'p','q','r','s'}, {0ul, nT});
      
      // inactive core indices
      addMORanges({'i','j','k','l'}, {offset, nInact});
      
      // correlated space
      addMORanges({'t','u','v','w'}, {corrOffset, nCorrO});

      // virtual indices
      addMORanges({'a','b','c','d'}, {corrOffset + nCorrO, nFVirt});
      
      if (mopart.scheme == DetScheme::RAS) {
        // '1', '2', '3' for RAS1 RAS2 RAS3 repectively 
        char nRAS = '1';
        size_t nRASOff = corrOffset;
        for (const auto &nActO: mopart.nActOs) {
          addMORanges({nRAS}, {nRASOff, nActO}); 
          nRAS++;
          nRASOff += nActO;
        }
      }

  }; // MOIntsTransformer::setMORanges for MCWaveFunction
  
  /**
   *  \brief parsing the mo ints types to offsizes 
   */
  template <typename MatsT, typename IntsT>
  std::vector<std::pair<size_t,size_t>> 
  MOIntsTransformer<MatsT,IntsT>::parseMOType(const std::string & moType) {
      
      // std::cout << " * parsing moType " << moType << std::endl;

      std::vector<std::pair<size_t,size_t>> off_sizes;
      
      for (const auto& ch: moType) off_sizes.push_back(parseMOType(ch));        
      
      return off_sizes;
  }; // MOIntsTransformer::parseMOIntsType(string)
  
  /**
   *  \brief parsing the mo ints type to offsize
   */
  template <typename MatsT, typename IntsT>
  std::pair<size_t,size_t> 
  MOIntsTransformer<MatsT,IntsT>::parseMOType(const char type) {

      std::pair<size_t,size_t> off_size;
      bool foundMOType = false;
      for (auto i = 0ul; i < symbol_sets_.size(); i++) {
        if (symbol_sets_[i].count(type)) {
          off_size = mo_ranges_[i];
          foundMOType = true;
          break;
        }
      }
      
      if (not foundMOType) CErr("Wrong MO Type in parseMOIntsType");
      
      return off_size;
  }; // MOIntsTransformer::parseMOIntsType(char)

  /*
   * \brief get unique symbols
   */ 
  template <typename MatsT, typename IntsT>
  char MOIntsTransformer<MatsT,IntsT>::getUniqueSymbol(char type) {
  
      for (auto i = 0ul; i < symbol_sets_.size(); i++) {
        if (symbol_sets_[i].count(type)) {
          return *symbol_sets_[i].begin();
        }
      }
      CErr("Wrong MO Type in parseMOIntsType");
      return ' ';
  }; // MOIntsTransformer::getUniqueSymbol

  /**
   *  \brief print offsizes 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::printOffSizes(
    const std::vector<std::pair<size_t,size_t>> & off_size) { 
    std::cout << "----MO Ranges: " << std::endl;
    
    size_t counter = 0;
    for (auto & size: off_size) {
      if (size.second == 0) continue; 
      counter ++;
      std::cout << "Index " << counter << ": from " << std::setw(5)
                << size.first + 1 << " to " << std::setw(5)
                << size.first + size.second << std::endl;
    }
    std::cout << std::endl;
  } // MOIntsTransformer::printOffSizes
  
  /**
   *  \brief print offsizes 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::printMORangesSummary() {
  
    std::cout << "  * Summary of MO ranges and symbols in MOIntsTransformer: " << std::endl;
    
    for (auto i = 0ul; i < mo_ranges_.size(); i++) {
      if (mo_ranges_[i].second == 0) continue; 
      std::cout << "    - MO Ranges " << std::setw(5) << mo_ranges_[i].first + 1 << " ~ "  
                << std::setw(5)  << mo_ranges_[i].first + mo_ranges_[i].second 
                << ", labeled as: ";
      
      for (const auto & c: symbol_sets_[i])
        std::cout << c << " ";
      
      std::cout << std::endl;
    }
  } // MOIntsTransformer::printMORangesSummary

}; // namespace ChronusQ
