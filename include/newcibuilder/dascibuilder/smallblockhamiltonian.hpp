/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <newcibuilder/dascibuilder.hpp>
#include <detfactory/excitationlist.hpp>
#include <tensor/tensorlooper.hpp>

namespace ChronusQ {
   
/*
 * Loop to build Sigma 2e Part for every CatK <- CatJ <- CatL
 * 
 * exListJK: CatJ <- Ex(q, p) <- CatK
 * exListJL: CatJ <- Ex(r, s) <- CatL
 *
 * Steps: 
 *   1. form X(NZ_LEx, NZ_KEx) <- ∑_rs ∑_pq <KEx|E_pq|JEx><JEx|Ers|LEx>S2e(pq,rs)  
 *   2. form cSCR(NZ_LEx, NonEx, iVec) <- C(L, iVec)  
 *   3. form sSCR(NZ_KEx, NonEx, iVec) <- ∑_L X(NZ_LEx, NZ_KEx) cSCR(NZ_LEx, NonEx, iVec)
 *   4. update Sigma2e(K, iVec) <- sSCR(NZ_KEx, NonEx, iVec)
 */ 
template <typename MatsT>
void DASCISigma2eBuilder::buildSmallBlockHamiltonian(
    size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
    const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
    const std::vector<size_t>& KExOffs, std::shared_ptr<TensorLooper>& KNonExLooper, 
    const std::vector<size_t>& LExOffs, std::shared_ptr<TensorLooper>& LNonExLooper,
    MatsT* cSCR, MatsT* sSCR, MatsT* X, const double symmFactor) {
  
  const size_t nJExDets = double1eExListsGen.totalDimension(); 
  const size_t nJExPerThread = std::ceil(double(nJExDets) / GetNumThreads());
  size_t JBegin = nJExPerThread * GetThreadID();
  size_t JEnd   = std::min(nJExDets, JBegin + nJExPerThread);
  if (JBegin >= JEnd) return;
  
  const size_t nNonExDets = LNonExLooper->nTotal();
  
  // binding addresses
  const auto& KNonEx = KNonExLooper->address();
  const auto& LNonEx = LNonExLooper->address();

  double1eExListsGen.visitExcitations(JBegin, JEnd, LExOffs, KExOffs,
      [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {
        
        // 0. get nonzero LEx and KEx map
        std::map<size_t, size_t> LExToNZLExMap, KExToNZKExMap;
        for (const auto& rsEx : rsExcitations) {
          LExToNZLExMap.try_emplace(rsEx[2], 0ul);
        }
        const size_t nNZLExDets = LExToNZLExMap.size();
        size_t count = 0ul;
        for (auto& [LEx, iNZLEx] : LExToNZLExMap) {
          iNZLEx = count;
          count++;
        }
        for (const auto& qpEx : qpExcitations) {
          KExToNZKExMap.try_emplace(qpEx[2], 0ul);
        }
        const size_t nNZKExDets = KExToNZKExMap.size();
        count = 0ul;
        for (auto& [KEx, iNZKEx] : KExToNZKExMap) {
          iNZKEx = count;
          count++;
        }
        
        // 1. form X(NZ_LEx, NZ_KEx) <- ∑_rs ∑_pq <KEx|E_pq|JEx><JEx|Ers|LEx>S2e(pq,rs)  
        MatsT* XCur = X;
        std::fill_n(X, nNZLExDets * nNZKExDets, MatsT(0.));
        for (const auto& [q, p, KEx, pqSign] : qpExcitations) {
          XCur = X + KExToNZKExMap[KEx] * nNZLExDets;
          for (const auto& [r, s, LEx, rsSign] : rsExcitations) {
            const auto iNZLEx = LExToNZLExMap[LEx];
            XCur[iNZLEx] += (rsSign == pqSign) ? s2e(p, q, r, s) : - s2e(p, q, r, s);      
          }
        }
       
        // 2. form cSCR(nNZ_LEx, NonEx, iVec) <- C(L, iVec)  
        const MatsT* VC = C;
        MatsT* cCur = cSCR;
        for (auto iVec = 0ul; iVec < nVec; ++iVec, VC+=LDC) 
        for (LNonExLooper->setIndex(0ul);
             not LNonExLooper->isEnd();
             LNonExLooper->increment()) 
        for (const auto& [LEx, iNZLEx] : LExToNZLExMap) {
          *cCur = VC[LEx + LNonEx];
          ++cCur;
        }

        // Matrix-Multiplication Step:
        // 3. form sSCR(NZ_KEx, NonEx, iVec) <- ∑_L X(NZ_LEx, NZ_KEx) cSCR(NZ_LEx, NonEx, iVec)
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZKExDets, nNonExDets * nVec, nNZLExDets, MatsT(symmFactor), X, nNZLExDets, cSCR, nNZLExDets, 
            MatsT(0.), sSCR, nNZKExDets);  

        // 4. update Sigma2e(K, iVec) <- sSCR(NZ_KEx, NonEx, iVec)
        MatsT* VS = Sigma;
        const MatsT* sCur = sSCR;
        for (auto iVec = 0ul; iVec < nVec; ++iVec, VS+=LDS) 
        for (KNonExLooper->setIndex(0ul);
             not KNonExLooper->isEnd();
             KNonExLooper->increment()) 
        for (const auto& [KEx, iNZKEx] : KExToNZKExMap) {
          VS[KEx + KNonEx] += *sCur;
          ++sCur;
        }
      } // inner loop 
  ); // visitAllExcitations 
} // DASCISigma2eBuilder<MatsT>::buildFrischLi
  
} // namespace ChronusQ
