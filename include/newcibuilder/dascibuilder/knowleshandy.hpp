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
 *   1. form X(NZrs, NZqp) <- <KEx|E_pq|JEx><JEx|Ers|LEx>S2e(pq,rs)  
 *   2. form Omega(NZrs, NonEx, iVec) <- ∑_L C(L, iVec) based on <JEx|Ers|LEx>  
 *   3. form Lambda(NZqp, NonEx, iVec) <- ∑_rs X(NZrs, NZqp) Omega(NZrs, NonEx, iVec)
 *   4. update Sigma2e(K, iVec) <- ∑pq Lambda(NZqp, NonEx, iVec) based on <KEx|E_pq|JEx>
 */ 
template <typename MatsT>
void DASCISigma2eBuilder::buildKnowlesHandy(
        size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
        const GASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
        const std::vector<size_t>& braExOffs, std::shared_ptr<TensorLooper>& braNonExLooper,
        const std::vector<size_t>& ketExOffs, std::shared_ptr<TensorLooper>& ketNonExLooper,
        MatsT* omega, MatsT* lambda, MatsT* X, const double symmFactor) {
  
  const size_t nJExDets = double1eExListsGen.totalDimension(); 
  const size_t nJExPerThread = std::ceil(double(nJExDets) / GetNumThreads());
  size_t JBegin = nJExPerThread * GetThreadID();
  size_t JEnd   = std::min(nJExDets, JBegin + nJExPerThread);
  if (JBegin >= JEnd) return;
  
  const size_t nNonExDets = ketNonExLooper->nTotal();
  const size_t nNZrs = double1eExListsGen.rsNNZ(); 
  const size_t nNZqp = double1eExListsGen.qpNNZ();

  // binding addresses
  const auto& braNonEx = braNonExLooper->address();
  const auto& ketNonEx = ketNonExLooper->address();

  double1eExListsGen.visitExcitations(JBegin, JEnd, ketExOffs, braExOffs,
                                      [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {
        
        // 1. form X(NZrs, NZqp) <- S2e(pq,rs) <KEx|E_pq|JEx><JEx|Ers|LEx>
        MatsT* XCur = X;
        for (const auto& qpEx : qpExcitations) {
          const auto& p = qpEx[1];
          const auto& q = qpEx[0]; 
          for (const auto& rsEx : rsExcitations) {
            const auto& r = rsEx[0]; 
            const auto& s = rsEx[1];
            *XCur = qpEx[3] == rsEx[3] ? s2e(p, q, r, s): -s2e(p, q, r, s);
            ++XCur; 
          }
        }
        
        // 2. form Omega(NZrs, NonEx, iVec) <- ∑_L C(L, iVec) based on <Jrs|Ers|Lrs>
        const MatsT* VC = C;
        MatsT* omg = omega;
        for (auto iVec = 0ul; iVec < nVec; ++iVec, VC+=LDC) 
        for (ketNonExLooper->setIndex(0ul);
             not ketNonExLooper->isEnd();
             ketNonExLooper->increment())
        for (const auto& rsEx : rsExcitations) {
            *omg = VC[rsEx[2] + ketNonEx];
            ++omg;
        }    

        // Matrix-Multiplication Step:
        // 3. form Lambda(NZqp, NonEx, iVec) <- ∑_rs X(NZrs, NZqp) * Omega(NZrs, NonEx, iVec)
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZqp, nNonExDets * nVec, nNZrs, MatsT(symmFactor), X, nNZrs, omega, nNZrs,
            MatsT(0.), lambda, nNZqp);  

        // 4. update Sigma2e(K, iVec) <- ∑_pq Lambda(NZqp, NonEx, iVec) based on <Kpq|Epq|Jpq>
        MatsT* VS = Sigma;
        const MatsT* lmb = lambda;
        for (auto iVec = 0ul; iVec < nVec; ++iVec, VS+=LDS) 
        for (braNonExLooper->setIndex(0ul);
             not braNonExLooper->isEnd();
             braNonExLooper->increment())
        for (const auto& qpEx : qpExcitations) {
          VS[qpEx[2] + braNonEx] += *lmb;
          ++lmb;
        } 
      } // inner loop 
  ); // visitAllExcitations 
} // DASCISigma2eBuilder<MatsT>::buildKnowlesHandy
  
} // namespace ChronusQ
