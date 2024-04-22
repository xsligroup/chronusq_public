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


#include <chronusq_sys.hpp>
#include <detstringmanager/detstring.hpp>
#include <detstringmanager/excitationlist.hpp>
#include <util/math.hpp>

namespace ChronusQ {
  
  enum ExcitationScheme {
    PRECOMPUTED_CONFIGURATION_DRIVEN_LIST,
    PRECOMPUTED_INTEGRAL_DRIVEN_LIST,  // TODO
    COMPUTING_EXCITATION_ON_THE_FLY    // TODO
  };
  
  class DetStringManager {
  
  protected:  
    
    size_t nE_;
    size_t nOrb_;
    size_t nStr_;
   
    ExcitationScheme scheme_;
  
  public:
    
    typedef std::vector<std::vector<int>> int_matrix; 
    // differentiate with NDet_ in CI class 
    
    // Disable default constructor, 
    // and use default copy and move constructors
    DetStringManager() = delete;
    DetStringManager(const DetStringManager & other): 
        DetStringManager(other.nOrb_, 
        other.nE_, other.scheme_) { };

    DetStringManager(DetStringManager && other):
        DetStringManager(other.nOrb_, 
        other.nE_, other.scheme_) { };

    DetStringManager(size_t nOrb, size_t nE,
        ExcitationScheme scheme):
        nOrb_(nOrb), nE_(nE), scheme_(scheme) {
    };

    virtual ~DetStringManager() { };
    
    // virtual function to do precomputations
    virtual void computeList() = 0;
    
    size_t nElectron() const { return nE_; } 
    size_t nOrbital()  const { return nOrb_; } 
    size_t nString()   const { return nStr_; }
    ExcitationScheme scheme() const { return scheme_; }

    // helper functions to handle string addressing for single graph
    int_matrix buildAddressingArray(size_t, size_t) const;
    int_matrix buildDeAddressingArray(size_t, size_t) const ;
    int_matrix buildDeAddressingArray(int_matrix &) const;
    size_t detString2Address(std::vector<size_t> &, int_matrix &) const;
    size_t detString2Address(bool *, int_matrix &) const;
    std::vector<size_t> address2DetString(size_t, int_matrix &, int_matrix &) const;
    DetString address2DetString(size_t, size_t, int_matrix &, int_matrix &) const;
    
    // helper functions for computing excitation list
    void computeSingleString1eExcitationList(DetString &, int *, size_t, size_t, int_matrix &) const; 
    ExcitationList computeIntraCASExcitationList(size_t, size_t, int_matrix &) const;
    
    void computeTwoString1eExcitationList(DetString &, DetString &, int *,
                    int, bool, int_matrix &, int_matrix &) const;
    ExcitationList computeInterCASExcitationList(size_t, size_t, size_t, size_t,
                    int, bool, int_matrix &, int_matrix &) const;

    // helper functions to initialize integer vectors
    int_matrix init_int_matrix(size_t d1, size_t d2, int val) const {
      return int_matrix(d1, std::vector<int>(d2, val)); };

  }; // class DetStringManager

  class CASStringManager: public DetStringManager {
  
  protected:

    std::shared_ptr<ExcitationList> exList_; 
    int_matrix addrArray_;  

  public:
    CASStringManager()                         = delete;
    CASStringManager(const CASStringManager &) = default;
    CASStringManager(CASStringManager &&)      = default;
    
    CASStringManager(size_t nOrb, size_t nE,
        ExcitationScheme scheme = PRECOMPUTED_CONFIGURATION_DRIVEN_LIST):
        DetStringManager(nOrb, nE, scheme) {
      addrArray_ = this->buildAddressingArray(nE, nOrb);
    }   
    
    ~CASStringManager() { dealloc(); }
    
    std::shared_ptr<const ExcitationList> excitationList() const { return exList_; }
     
    void computeList();
    
    void dealloc() { exList_ = nullptr; }

  }; // class CASStringManager

  class RASStringManager: public DetStringManager {

  protected:

    std::vector<size_t> nActO_; // number of orbitals in RAS 1/2/3 spaces
    std::unordered_map<std::string,
                       std::shared_ptr<ExcitationList>> exList_;
    size_t nCat_; // total number of string categories for RAS
    size_t mxHole_; // maximum number of holes in RAS 1
    size_t mxElec_; // maximum number of electrons in RAS 3
    std::vector<std::vector<size_t>> LCat_; // number of strings in each space for each category, NCat by 4;
    
  public:
    RASStringManager()                         = delete;
    RASStringManager(const RASStringManager &) = default;
    RASStringManager(RASStringManager &&)      = default;

    RASStringManager(std::vector<size_t> nActO, size_t nE,
	                 size_t mxHole, size_t mxElec,
          ExcitationScheme scheme = PRECOMPUTED_CONFIGURATION_DRIVEN_LIST):
	  nActO_(nActO), mxHole_(mxHole), mxElec_(mxElec),
      DetStringManager(std::accumulate(nActO.begin(), nActO.end(), 0),
      nE, scheme) {

      nCat_ = (mxHole_+1) * (mxElec_+1);
      this->genLCat();
    }
 
    ~RASStringManager() { };

    std::vector<std::vector<size_t>> LCategory() const { return LCat_; }
    size_t nCategory() const { return nCat_; }
    size_t maxHole() const { return mxHole_; }
    size_t maxElectron() const { return mxElec_; }

    std::string mapKey(size_t, size_t, size_t, size_t) const;

    std::shared_ptr<const ExcitationList> find_exlist(size_t, size_t, size_t, size_t) const;

    void genLCat();

    void computeList();

    //void dealloc() { exList_ = nullptr; }

  }; // class RASStringManager


}; // namespace ChronusQ
