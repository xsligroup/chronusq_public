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

//#define _DEBUG_RASSTRING

#include <detstringmanager.hpp> 

namespace ChronusQ {

  typedef DetStringManager::int_matrix int_matrix; 

  /*
   * Generate unordered_map Key for RAS excitation lists
   *   p, q ...... p, q blocks
   *   L, K ...... Categories for strings L and K
   *   return a string "p_q_L_K"
   */
  std::string RASStringManager::mapKey(size_t p, size_t q,
                                size_t L, size_t K) const {

    return std::to_string(p) + "_" + std::to_string(q)
      + "_" + std::to_string(L) + "_" + std::to_string(K);

  } // RASStringManager::mapKey

  /*
   * find the Excitation list corresponding to a set of input:
   *    p, q ...... p, q blocks
   *    L, K ...... Categories for strings L and K
   *    return a pointer to the ExcitationList
   */
  std::shared_ptr<const ExcitationList> RASStringManager::find_exlist(
                size_t p, size_t q, size_t L, size_t K) const {

    std::string map_key = this->mapKey(p, q, L, K);
    auto search = this->exList_.find(map_key);
    if (search != this->exList_.end() ) return search->second;
    else return nullptr;

  } // RASStringManager::find_exlist

  /* 
   * Generate the offset of each space for each category
   *   nStr: total number of strings for RAS case.
   *   LCat: nCat_ x 4
   *     RAS 1 offset ... 1
   *     RAS 2 offset ... N(RAS1)
   *     RAS 3 offset ... N(RAS1)*N(RAS2)
   *     total number of string for that category iCat ... N(RAS1)*NRAS2)*N(RAS3)
   */ 
  void RASStringManager::genLCat() {

    LCat_ = std::vector<std::vector<size_t>>(nCat_, std::vector<size_t>(4, 0));
    
    this->nStr_ = 0ul;
    for (auto iE = 0ul, iCat = 0ul; iE <= mxElec_; iE++)
    for (auto iH = 0ul; iH <= mxHole_; iH++, iCat++) {
      if((nE_+iH) < (nActO_[0]+iE) or (nE_+iH) > (nActO_[0]+iE+nActO_[1])) continue;
      size_t nEras2 = nE_ - nActO_[0] + iH - iE;
      LCat_[iCat][0] = 1;
      LCat_[iCat][1] = Comb(nActO_[0], iH);
      LCat_[iCat][2] = LCat_[iCat][1] * Comb(nActO_[1], nEras2);
      LCat_[iCat][3] = LCat_[iCat][2] * Comb(nActO_[2], iE);
      this->nStr_ += LCat_[iCat][3];
    }
    
  } // RASStringManager::genLCat

  void RASStringManager::computeList() {

//    auto RASStringSt = tick();

    if(this->scheme_ == PRECOMPUTED_CONFIGURATION_DRIVEN_LIST){
    for (auto pBlk = 0ul; pBlk < 3; pBlk++)
    for (auto qBlk = 0ul; qBlk < 3; qBlk++) {
      int iCatL = 0;
      std::unordered_map<size_t, std::shared_ptr<ExcitationList>> LCache;

      for (auto iE = 0ul; iE <= mxElec_; iE++)
      for (auto iH = 0ul; iH <= mxHole_; iH++, iCatL++) {
        
        if (LCat_[iCatL][3] ==0) continue;
        
        int_matrix nEras = this->init_int_matrix(3, 2, 0);
        int_matrix nDras = this->init_int_matrix(3, 2, 0);
        std::shared_ptr<ExcitationList> ras_exl;

        nEras[0][0] = nActO_[0] - iH;
        nEras[2][0] = iE;
        nEras[1][0] = nE_ - nEras[0][0] - nEras[2][0];

        nDras[0][0] = LCat_[iCatL][1];
        nDras[1][0] = LCat_[iCatL][2]/LCat_[iCatL][1];
        nDras[2][0] = LCat_[iCatL][3]/LCat_[iCatL][2];

        if(pBlk == qBlk) { // diagonal block
          int nVp = nActO_[pBlk] - nEras[pBlk][0];
          int nNZ = nEras[pBlk][0] * (nVp + 1);

          if(nNZ == 0) continue;
          else {
            auto search = LCache.find(nEras[pBlk][0]);
            if (search != LCache.end() ) ras_exl = search->second;
            else {
              int_matrix ras_addr_array =
                  this->buildAddressingArray(nEras[pBlk][0], nActO_[pBlk]);
              ras_exl = std::make_shared<ExcitationList>(
                this->computeIntraCASExcitationList(nActO_[pBlk],
                                nEras[pBlk][0], ras_addr_array));
              LCache.insert({nEras[pBlk][0], ras_exl});
            }
            exList_.insert({mapKey(pBlk, qBlk, iCatL, iCatL), ras_exl});
#ifdef _DEBUG_RASSTRING
            std::cout << "p: " << pBlk+1 << " and q: " << qBlk+1 << std::endl;
            std::cout << "LCat: " << iCatL+1 << " and KCat: " << iCatL + 1 << std::endl;
#endif
          }

        } else { // off-diagonal block
          for (auto i = 0ul; i < 3; i++)
            nEras[i][1] = nEras[i][0];
          nEras[pBlk][1] += 1;
          nEras[qBlk][1] -= 1;

          bool checko = false;
          for (auto i = 0ul; i < 3; i++) {
            if(nEras[i][1] < 0 or nEras[i][1] > nActO_[i]) {
              checko = true;
              break;
            }
          }
          if(not checko) {
            int iHK = nActO_[0] - nEras[0][1];
            int iEK = nEras[2][1];
            if(iHK > mxHole_ or iEK > mxElec_) continue;
            int iCatK = iHK + (mxHole_ + 1) * iEK;
            if(LCat_[iCatK][3] == 0) continue;

            nDras[0][1] = LCat_[iCatK][1];
            nDras[1][1] = LCat_[iCatK][2]/LCat_[iCatK][1];
            nDras[2][1] = LCat_[iCatK][3]/LCat_[iCatK][2];
            size_t nVpL = nActO_[pBlk] - nEras[pBlk][0];
            size_t nNZL = nEras[qBlk][0] * nVpL;
            size_t nVqK = nActO_[qBlk] - nEras[qBlk][1];
            size_t nNZK = nEras[pBlk][1] * nVqK;

            if(nNZL == 0 or nNZK == 0) continue;

#ifdef _DEBUG_RASSTRING
            std::cout << "p: " << pBlk+1 << " and q: " << qBlk+1 << std::endl;
            std::cout << "LCat: " << iCatL+1 << " and KCat: " << iCatK + 1 << std::endl;
#endif

            int Sgnoff = 1;
            if((pBlk == 2 and qBlk == 0) or (pBlk == 0 and qBlk == 2))
              if(nEras[1][1] % 2 == 1) Sgnoff = -1;
            bool pqseq = (pBlk < qBlk) ? true : false;

            // form addressing array
            int_matrix addr_array_pK = buildAddressingArray(nEras[pBlk][1], nActO_[pBlk]);
            int_matrix addr_array_qK = buildAddressingArray(nEras[qBlk][1], nActO_[qBlk]);

            ras_exl = std::make_shared<ExcitationList>(
              this->computeInterCASExcitationList(nActO_[pBlk], nEras[pBlk][0],
                        nActO_[qBlk], nEras[qBlk][0], Sgnoff, pqseq,
                        addr_array_pK, addr_array_qK));
            exList_.insert({mapKey(pBlk, qBlk, iCatL, iCatK),
                            ras_exl});
          }
        }
      }
    }

//    double RASStringdur = tock(RASStringSt);
//    std::cout << "\nRAS String - DURATION = " << std::setprecision(8)
//                << RASStringdur << " s." << std::endl;

#ifdef _DEBUG_RASSTRING
    std::cout << "RAS string all finished." << std::endl;
#endif

  } else CErr("Not Implemented yet");

  } // RASString::buildConfDrivenExcitationList

}; // namespace ChronusQ
