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

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
// #define _DEBUG_CIBuilder_IMPL

namespace ChronusQ {
  
/*
 *  Default Build Sigma one-electron part for DASCI 
 */ 
template <typename MatsT>
void DASCIBuilder<MatsT>::buildSigma1e(
  const LocalCIVectorsView<const MatsT>& C,
  const LocalCIVectorsView<MatsT>& Sigma) const {
  
  assert(C.size() == Sigma.size());
  size_t nVec = C.size();

  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();
   
  for(const auto& oneEEx: this->detFactory_.oneEExcitations()) {
    
    // std::cout << " contraction on term - " << oneEEx.term << std::endl; 
    
    // MPI Parallelism
    if (not C.containsLocalCategory(oneEEx.categoricalIndices.second)) continue;
    
    const auto& h1e = *(this->moints_->template getIntegral<DASOnePInts, MatsT>(oneEEx.term));
    // h1e.output(std::cout, oneEEx.term, true);
    const auto& braCategory = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(oneEEx.categoricalIndices.first));
    const auto& ketCategory = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(oneEEx.categoricalIndices.second));
    // This is the symmetry factor between spaces
    const auto& symmFact = oneEEx.symmetryFactor; 
    const auto& exList = dynamic_cast<const FullCD1eExList&>(*oneEEx.exLists[0]); 
    
    const auto pqExSpaces = oneEEx.exSpaces; 
    std::unordered_set<size_t> exSpaces({pqExSpaces[0], pqExSpaces[1]});

    // separate and group non-excitation spaces.
    // continuous non-excitation spaces are grouped together and the total dimension
    // is equal to the product of individual dimensions.
    std::vector<size_t> KExOffsBuffer, KNonExOffs, LExOffsBuffer, LNonExOffs, dummy, nonExDims;
    braCategory.separateExAndNonExDimensions(exSpaces, KExOffsBuffer, dummy, KNonExOffs, nonExDims);
    ketCategory.separateExAndNonExDimensions(exSpaces, LExOffsBuffer, dummy, LNonExOffs, dummy);
    
    std::pair<size_t, size_t> ketExOffs, braExOffs;
    if (exSpaces.size() == 1) {
        ketExOffs = {LExOffsBuffer[0], 0ul};
        braExOffs = {KExOffsBuffer[0], 0ul};
    } else if (pqExSpaces[0] < pqExSpaces[1]) {
        ketExOffs = {LExOffsBuffer[0], LExOffsBuffer[1]};
        braExOffs = {KExOffsBuffer[0], KExOffsBuffer[1]};
    } else {
        ketExOffs = {LExOffsBuffer[1], LExOffsBuffer[0]};
        braExOffs = {KExOffsBuffer[1], KExOffsBuffer[0]};
    }

    const size_t nExBraDets = exList.nBraDeterminants();
    const size_t nExBraDetsPerThread = std::ceil(double(nExBraDets) / GetNumThreads());

    #pragma omp parallel default(shared)
    {
      // multidimensional index advancer
      auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);
      // bind non-excitation part addresses
      const auto& ketNonEx = nonExLooper->address();
      const auto& braNonEx = nonExLooper->auxAddress();
      
      auto exListGen = exList.generator(ketExOffs, braExOffs);
      // bind braEx with advancer called inside visitExcitations
      const auto& braEx = exListGen->braExAddress();
      size_t iExBraDetBegin = nExBraDetsPerThread * GetThreadID();
      size_t iExBraDetEnd   = std::min(nExBraDets, iExBraDetBegin + nExBraDetsPerThread);
     
      exListGen->visitExcitations(iExBraDetBegin, iExBraDetEnd,
                                  [&] (const auto& braExIter, const auto& pqExcitations) {
            for (const auto& [p, q, ketEx, pqSign]: pqExcitations) {
              // between-space symmetry factor multiplied by the between-orbital symmetry factor
              MatsT h1e_pq = pqSign ? -symmFact * h1e(p, q) : symmFact * h1e(p, q);

              for (auto iVec = 0ul; iVec < nVec; ++iVec) {
                MatsT* sigmaVecAtBraCat = Sigma.getCategoryPointer(oneEEx.categoricalIndices.first, iVec);
                const MatsT* coeffVecAtKetCat = C.getCategoryPointer(oneEEx.categoricalIndices.second, iVec);

                // For a given BraKet pair of determinants in the excitation spaces,
                // TensorLooper loops over all possible determinants in the non-excitation spaces.
                for (nonExLooper->setIndex(0ul);
                     not nonExLooper->isEnd();
                     nonExLooper->increment()) {
                   auto K = braEx + braNonEx;
                   auto L = ketEx + ketNonEx;
                   sigmaVecAtBraCat[K] += h1e_pq * coeffVecAtKetCat[L];
                }
              }
            }
          }
      );
    } // parallel region
  } // oneEExcitations
  // Sigma.print(std::cout, "HH full H after Sigma loop");
} // DASCIBuilder::buildSigma1e


/*
 * Loop to build Sigma 2e Excitations
 * 
 * Naive Loop: useful for debugging!
 *
 */ 
template <typename MatsT>
void DASCISigma2eBuilder::buildNaive(
    size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
    const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
    const std::vector<size_t>& KExOffs, const std::vector<size_t>& LExOffs, 
    std::shared_ptr<TensorLooper>& nonExLooper, const double symmFact) {
    
  const size_t nJExDets = double1eExListsGen.totalDimension(); 
  const size_t nJExPerThread = std::ceil(double(nJExDets) / GetNumThreads());
  size_t JBegin = nJExPerThread * GetThreadID();
  size_t JEnd   = std::min(nJExDets, JBegin + nJExPerThread);
  if (JBegin >= JEnd) return;
  
  // binding non excitation part address
  const auto& LNonEx = nonExLooper->address();
  const auto& KNonEx = nonExLooper->auxAddress();

  double1eExListsGen.visitExcitations(JBegin, JEnd, LExOffs, KExOffs,
      [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {
        
        for (const auto& [q, p, KEx, pqSign]: qpExcitations) 
        for (const auto& [r, s, LEx, rsSign]: rsExcitations) {
          double fc = (pqSign == rsSign) ? symmFact : -symmFact;
          const MatsT* VC = C;
          MatsT* VS = Sigma;
          for (auto iVec = 0ul; iVec < nVec; ++iVec, VC+=LDC, VS+=LDS) {
           
            for (nonExLooper->setIndex(0ul);
                 not nonExLooper->isEnd();
                 nonExLooper->increment()) {
               auto K = KEx + KNonEx;
               auto L = LEx + LNonEx;
            
               VS[K] += fc * s2e(p, q, r, s) * VC[L]; 
           
            }
          }
        }
      }
  ); // visitAllExcitations
} // DASCISigma2eBuilder::buildSigma2eExcitation

} // namespace ChronusQ
