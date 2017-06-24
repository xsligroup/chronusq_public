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

//#define _DEBUG_DETSTRING

#include <cxxapi/output.hpp>
#include <util/matout.hpp>
#include <util/math.hpp>
#include <detstringmanager.hpp>
#include <detstringmanager/excitationlist.hpp>
#include <detstringmanager/stringpermutator.hpp>

namespace ChronusQ {
  
  typedef DetStringManager::int_matrix int_matrix; 
  
  /*Generate and return the addrssing array for given number of
   * orbitals nOrb and number of electron nE
   */
  
  int_matrix DetStringManager::buildAddressingArray(size_t nE, size_t nOrb) const {
    
    assert (nE <= nOrb and nOrb != 0 ); 
    int_matrix addr_array;
    
    if (nE == 0) addr_array = {{}};
    else {
      addr_array = this->init_int_matrix(nE, nOrb, 0);
      
      for ( auto k = 1ul; k < nE; k++)
      for ( auto l = k  ; l < nOrb-nE+k+1; l++)
      for ( auto m = nOrb-l+1; m < nOrb-k+1; m++) {
        addr_array[k-1][l-1] += Comb(m-1, nE-k);
      }
      
      for ( auto l = nE; l < nOrb+1; l++) {
        addr_array[nE-1][l-1] = l - nE;
      }
    }

    return  addr_array;
  }; // DetStringManager::buildAddressingArray
  
  int_matrix DetStringManager::buildDeAddressingArray(int_matrix & addr_array) const {
    
    size_t nE   = addr_array.size(); 
    size_t nOrb = addr_array[0].size();
    if (nE == 0) return {{}};
    
    int_matrix de_addr_array = this->init_int_matrix(nE, nOrb, 0);
    
    std::copy_n(addr_array[nE - 1].begin(), nOrb, de_addr_array[nE - 1].begin());
     
    for (int iE = nE - 2; iE >= 0; iE--)
    for (int iOrb = 1; iOrb < nOrb - 1; iOrb++) 
      de_addr_array[iE][iOrb] = addr_array[iE][iOrb] + de_addr_array[iE + 1][iOrb + 1];

    return de_addr_array;

  }; // DetStringManager::buildDeAddressingArray(int_matrix)
  
  int_matrix DetStringManager::buildDeAddressingArray(size_t nE, size_t nOrb) const {
    
    auto addr_array = buildAddressingArray(nE, nOrb);
    
    return buildDeAddressingArray(addr_array);
   
  }; //DetStringManager::buildDeAddressingArray(size_t, size_t)
    
  
  // Determine the address of K, CPL 111 p137 
  // !!address start from 0 in C++
  // CAS type addressing
  size_t DetStringManager::detString2Address(std::vector<size_t> & elecPos, 
    int_matrix & addr_array) const {
      
      size_t addr = 0;
      for (auto iE = 0ul; iE < elecPos.size(); iE++)
        addr += addr_array[iE][elecPos[iE]];

      return addr;
  };
  
  size_t DetStringManager::detString2Address(bool * detStr, 
      int_matrix & addr_array) const {
    
    size_t nE   = addr_array.size(); 
    size_t nOrb = addr_array[0].size();
    if (nOrb == 0) return 0; 
    
    std::vector<size_t> elecPos;
    
    for (auto iE = 0ul, iOrb = 0ul; (iOrb < nOrb) && (iE < nE); iOrb++) {
      if (detStr[iOrb]) {
        elecPos.push_back(iOrb);
        iE ++;
      }   
    }
    
    return detString2Address(elecPos, addr_array);

  }; // DetStringManager::detString2Address
  
  std::vector<size_t> DetStringManager::address2DetString(size_t addr, 
      int_matrix & addr_array, int_matrix & de_addr_array) const {

    size_t nE   = addr_array.size(); 
    if (nE == 0) return {};
    
    size_t nOrb = addr_array[0].size();
    std::vector<size_t> elecPos;
    std::pair<int, int> iOrb_search = {nOrb - nE, 0};

    for (int iE = 0; iE < nE; iE++)
    for (int iOrb = iOrb_search.first; iOrb >= iOrb_search.second; iOrb--) {
      if (addr >= de_addr_array[iE][iOrb]) {
        elecPos.push_back(iOrb);
        addr -= addr_array[iE][iOrb];
        iOrb_search.first ++;
        iOrb_search.second = iOrb + 1;
        break;
      }
    }

    return elecPos;
  }; // DetStringManager::address2DetString

  DetString DetStringManager::address2DetString(size_t addr, size_t nOrb, 
      int_matrix & addr_array, int_matrix & de_addr_array) const {
   
    return DetString(nOrb, address2DetString(addr, addr_array, de_addr_array));
     
  }; // DetStringManager::address2DetString
  
  /* 
   * 
   * compute all <K|Epq|L> for a given L 
   * Configuration Driven 1e Excitation List: <K|a_i^\dagger a_j |L> 
   *               (1) nonzero excitation i, 
   *               (2) nonzero excitation j, 
   *               (3) address of excited string |K>
   *               (4) sign (+/- 1)
   */ 
  void DetStringManager::computeSingleString1eExcitationList(
      DetString & detStr, int * exList, size_t nOrb, 
      size_t nE, int_matrix & addrArray) const { 
    
    std::vector<size_t> LOcc = detStr.occupationInfo();
    std::vector<size_t> LVir = detStr.virtualInfo();
    /*
    std::cout << "LOcc.size() = " << LOcc.size() << std::endl;
    for(auto & i: LOcc) std::cout << " " << i;
    std::cout << std::endl;
    std::cout << "LVir.size() = " << LVir.size() << std::endl;
    for(auto & i: LVir) std::cout << " " << i;
    std::cout << std::endl;
    */

    int * exList_ptr = exList;
    size_t p, q;
    size_t LAddr = detString2Address(detStr.pointer(), addrArray); 
    
    // case 1: self-excitation, p=q == occ
    for (auto pOcc = 0ul ; pOcc < nE; pOcc++) {
	  p = LOcc[pOcc];
      exList_ptr[0] = p;
      exList_ptr[1] = p;
      exList_ptr[2] = LAddr;
      exList_ptr[3] = 1;
      exList_ptr += 4;
    }
    
    // Case 2: q -> p, p=virtual, q=occupied in L
    for (auto pVir = 0ul; pVir < nOrb - nE; pVir++)
    for (auto qOcc = 0ul; qOcc < nE; qOcc++) {
      p = LVir[pVir];
	  q = LOcc[qOcc];
      exList_ptr[0] = p;
      exList_ptr[1] = q;
      detStr.excitation(q,p);
      exList_ptr[2] = detString2Address(detStr.pointer(), addrArray);
      exList_ptr[3] = detStr.sign1e(p,q); 
      detStr.excitation(p,q);
      exList_ptr += 4;  
    }
    
    //prettyPrintSmart(std::cout, "iStr = " + std::to_string(LAddr), 
    //  exList, 4, (nOrb - nE + 1)*nE, 4);//, 1, 10, 8);
    
  }; // DetStringManager::computeSingleString1eExcitationList


  /*
   * Compute the excitation list <K|Eij|L> 
   *   L is formed from two DetStrings in difference spaces,
   *   i.e. detsp and detsq
   * Configuration Driven 1e Excitation List: 
   *   <K|a_i^\dagger a_j |L>
   *        (1) nonzero excitation i,
   *        (2) nonzero excitation j,
   *        (3) address of excited string |K> in p block,
   *        (4) address of excited string |K> in q block,
   *        (5) sign (+/- 1)
   */
  void DetStringManager::computeTwoString1eExcitationList(
    DetString & detsp, DetString & detsq, int * exList, 
    int Sgnoff, bool pqseq,
    int_matrix & addrArrayp, int_matrix & addrArrayq) const {

    std::vector<size_t> qOcc = detsq.occupationInfo();
    std::vector<size_t> pVir = detsp.virtualInfo();

    int * exList_ptr = exList;
    size_t i, j;

    for (auto jq = 0ul; jq < detsq.n1s(); jq++)
    for (auto ip = 0ul; ip < detsp.n0s(); ip++) {
      j = qOcc[jq];
      i = pVir[ip];
      detsq.flip(j);
      detsp.flip(i);

      int Sign1e;
      if (pqseq) {
        DetString tmp = detsp + detsq;
        Sign1e = Sgnoff * tmp.sign1e(i, j + detsp.size());
      } else {
        DetString tmp = detsq + detsp;
        Sign1e = Sgnoff * tmp.sign1e(j, i + detsq.size());
      }

      exList_ptr[0] = i;
      exList_ptr[1] = j;
      exList_ptr[2] = detString2Address(detsp.pointer(), addrArrayp);
      exList_ptr[3] = detString2Address(detsq.pointer(), addrArrayq);
      exList_ptr[4] = Sign1e;
      exList_ptr += 5;

      detsq.flip(j);
      detsp.flip(i);
    }
  } // DetStringManager::computeTwoString1eExcitationList

  
  ExcitationList DetStringManager::computeIntraCASExcitationList(
      size_t nOrb, size_t nE, int_matrix & addrArray) const {
     
    size_t nNZ = nE * (nOrb - nE + 1);
    size_t nStr = Comb(nOrb, nE);
    ExcitationList exList(4, nNZ, nStr);
    
    StringPermutator strList(nOrb, nE);
    
    size_t iStr = 0;
    
    // TODO: Parallel this loop
    while (iStr != nStr) { 
      
      const std::vector<size_t> & eLocs = strList.permutedPositions();  
      DetString L(nOrb, eLocs); 
      
      computeSingleString1eExcitationList(L, 
        exList.pointerAtDet(iStr), nOrb, nE, addrArray);
      
      iStr = strList.next();
    }

#ifdef _DEBUG_DETSTRING
    exList.output(std::cout);  
#endif  
    return exList;
  };// DetStringManager::computeIntraCASExcitationList


  /*
   * Compute excitation list between spaces
   *   e.g off-diagonal block in RAS string
   *   exList dimensions: nStr2_: q block |L> address,
   *                      nStr1_: p block |L> address.
   */
  ExcitationList DetStringManager::computeInterCASExcitationList(
                 size_t nOp, size_t nEp, size_t nOq, size_t nEq,
                 int Sgnoff, bool pqseq, int_matrix & addr_array_pK,
                 int_matrix & addr_array_qK) const {

    size_t nVp = nOp - nEp;
    size_t nNZ = nEq * nVp;
    size_t nStrp = Comb(nOp, nEp);
    size_t nStrq = Comb(nOq, nEq);

    ExcitationList exList(5, nNZ, nStrp, nStrq);
    
    StringPermutator strList_qL(nOq, nEq); // SOqL
    StringPermutator strList_pL(nOp, nEp); // SOpL

    size_t iStr_qL = 0;
    size_t iStr_pL = 0;

    // TODO: Parallel this loop
    while (iStr_qL != nStrq) {
      const std::vector<size_t> & eLocs_qL = strList_qL.permutedPositions();
      DetString tmpQ(nOq, eLocs_qL);

      strList_pL.reset_seed();
      iStr_pL = 0;
      
      while (iStr_pL != nStrp) {
        const std::vector<size_t> & eLocs_pL = strList_pL.permutedPositions();
        DetString tmpP(nOp, eLocs_pL);

        computeTwoString1eExcitationList(tmpP, tmpQ,
          exList.pointerAtDet(iStr_pL, iStr_qL),
          Sgnoff, pqseq, addr_array_pK, addr_array_qK);

        iStr_pL = strList_pL.next();
      }
      iStr_qL = strList_qL.next();
    }

#ifdef _DEBUG_DETSTRING
    exList.output(std::cout);
#endif
    return exList;

  } // DetStringManager::computeInterCASExcitationList
  
}; // namespace ChronusQ

// Other headers
#include <detstringmanager/casstringmanager.hpp>
#include <detstringmanager/rasstringmanager.hpp>


