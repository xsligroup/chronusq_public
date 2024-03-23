/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
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
 *   1. form X(LEx, KEx) <- ∑_rs ∑_pq ∑_J_Ex <KEx|E_pq|JEx><JEx|Ers|LEx>S2e(pq,rs)  
 *   2. form cSCR(NZ_LEx, NonEx) <- C(L)  
 *   3. form sSCR(NZ_KEx, NonEx) <- ∑_L X(NZ_LEx, NZ_KEx) cSCR(NZ_LEx, NonEx)
 *   4. update Sigma2e(K) <- sSCR(NZ_KEx, NonEx)
 */ 
template <typename MatsT>
void DASCISigma2eBuilder::buildLargeBlockHamiltonian(
    size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
    const GASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
    std::shared_ptr<TensorLooper>& KExLooper, std::shared_ptr<TensorLooper>& KNonExLooper, 
    std::shared_ptr<TensorLooper>& LExLooper, std::shared_ptr<TensorLooper>& LNonExLooper, 
    MatsT* cSCR, MatsT* sSCR, MatsT* X, const double symmFactor) {
  
  const size_t nJExDets = double1eExListsGen.totalDimension(); 
  const size_t nJExPerThread = std::ceil(double(nJExDets) / GetNumThreads());
  size_t JBegin = nJExPerThread * GetThreadID();
  size_t JEnd   = std::min(nJExDets, JBegin + nJExPerThread);
  if (JBegin >= JEnd) return;
  
  const size_t nNonExDets = LNonExLooper->nTotal();
  const size_t nKExDets = KExLooper->nTotal();
  const size_t nLExDets = LExLooper->nTotal();
  
  // 1. form X(LEx, KEx) <- ∑_rs ∑_pq ∑_J_Ex <KEx|E_pq|JEx><JEx|Ers|LEx>S2e(pq,rs)  
  std::fill_n(X, nLExDets * nKExDets, MatsT(0.));
  double1eExListsGen.visitExcitations(JBegin, JEnd,
      LExLooper->indexOffsets(), KExLooper->indexOffsets(),
      [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {
        MatsT* XCur = X;
        for (const auto& [q, p, KEx, pqSign] : qpExcitations) {
          XCur = X + KEx * nLExDets;
          for (const auto& [r, s, LEx, rsSign] : rsExcitations) {
            XCur[LEx] += (rsSign == pqSign) ? s2e(p, q, r, s) : - s2e(p, q, r, s);      
          }
        }
      }
  ); // visitAllExcitations 
  
  // binding addresses
  const auto& KEx = KExLooper->address();
  const auto& LEx = LExLooper->address();
  const auto& KNonEx = KNonExLooper->address();
  const auto& LNonEx = LNonExLooper->address();

  const MatsT* VC = C;
  MatsT* VS = Sigma;
  for (auto iVec = 0ul; iVec < nVec; ++iVec, VC+=LDC, VS+=LDS) {
    // 2. transpose cSCR(LEx, NonEx) <- C(L)  
    MatsT* cCur = cSCR;
    for (LNonExLooper->setIndex(0ul);
         not LNonExLooper->isEnd();
         LNonExLooper->increment()) {
      for (LExLooper->setIndex(0ul);
           not LExLooper->isEnd();
           LExLooper->increment()) {
        *cCur = VC[LEx + LNonEx];
        ++cCur;
      }
    }

    // Matrix-Multiplication Step:
    // 3. form sSCR(KEx, NonEx) <- ∑_L X(LEx, KEx) cSCR(LEx, NonEx)
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
        nKExDets, nNonExDets, nLExDets, MatsT(symmFactor), X, nLExDets, cSCR, nLExDets, 
        MatsT(0.), sSCR, nKExDets);  

    // 4. update Sigma2e(K) <- sSCR(KEx, NonEx)
    MatsT* sCur = sSCR;
    for (KNonExLooper->setIndex(0ul);
         not KNonExLooper->isEnd();
         KNonExLooper->increment()) {
      for (KExLooper->setIndex(0ul);
           not KExLooper->isEnd();
           KExLooper->increment()) {
        VS[KEx + KNonEx] += *sCur;
        ++sCur;
      }
    }
  } // iVec
} // DASCISigma2eBuilder<MatsT>::buildFrischLi
  
} // namespace ChronusQ
