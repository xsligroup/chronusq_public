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

#include <cqlinalg/blasext.hpp>
#include <newcibuilder/dascibuilder.hpp>
#include <detfactory/excitationlist.hpp>

namespace ChronusQ {
   
template <typename MatsT>
void DASCIBuilder<MatsT>::build1TDM(
    const LocalCIVectorsView<const MatsT>& CBra, 
    const LocalCIVectorsView<const MatsT>& CKet, 
    cqmatrix::Matrix<MatsT>& oneTDM,
    const double scale) const {
  
  const auto& orbOffs = this->detFactory_.orbOffsInEachSpace();

  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();
  
  std::vector<cqmatrix::Matrix<MatsT>> SCR;
  for (auto i = 0ul; i < GetNumThreads() - 1; i++) {
    SCR.emplace_back(oneTDM.dimension());
    SCR.back().clear();
  }

  for(const auto& oneEEx: this->detFactory_.oneEExcitations()) {
    // std::cout << " contraction on term - " << oneEEx.term << std::endl; 
    
    // MPI Parallelism
    if (not CKet.containsLocalCategory(oneEEx.categoricalIndices.second)) continue;
    
    const auto& catK = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(oneEEx.categoricalIndices.first));
    const auto& catL = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(oneEEx.categoricalIndices.second));
    const auto& symmFact = oneEEx.symmetryFactor * scale; 
    const auto& exList = dynamic_cast<const FullCD1eExList&>(*oneEEx.exLists[0]); 
    
    const auto pqExSpaces = oneEEx.exSpaces;
    const auto pOrbOff = orbOffs[pqExSpaces[0]];
    const auto qOrbOff = orbOffs[pqExSpaces[1]];
    std::unordered_set<size_t> exSpaces({pqExSpaces[0], pqExSpaces[1]});
    std::vector<size_t> KExOffsBuffer, KNonExOffs, LExOffsBuffer, LNonExOffs, dummy, nonExDims;
    catK.separateExAndNonExDimensions(exSpaces, KExOffsBuffer, dummy, KNonExOffs, nonExDims);
    catL.separateExAndNonExDimensions(exSpaces, LExOffsBuffer, dummy, LNonExOffs, dummy);
    
    std::pair<size_t, size_t> LExOffs, KExOffs;
    if (exSpaces.size() == 1) {
      LExOffs = {LExOffsBuffer[0], 0ul};
      KExOffs = {KExOffsBuffer[0], 0ul};
    } else if (pqExSpaces[0] < pqExSpaces[1]) {
      LExOffs = {LExOffsBuffer[0], LExOffsBuffer[1]};
      KExOffs = {KExOffsBuffer[0], KExOffsBuffer[1]};
    } else {
      LExOffs = {LExOffsBuffer[1], LExOffsBuffer[0]};
      KExOffs = {KExOffsBuffer[1], KExOffsBuffer[0]};
    }

    const MatsT* CBra_ptr = CBra.getCategoryPointer(oneEEx.categoricalIndices.first); 
    const MatsT* CKet_ptr = CKet.getCategoryPointer(oneEEx.categoricalIndices.second);  
    const size_t nKExDets = exList.nBraDeterminants();
    const size_t nKExPerThread = std::ceil(double(nKExDets) / GetNumThreads());
  
    #pragma omp parallel default(shared)
    {
      auto exListGen = exList.generator(LExOffs, KExOffs);
      const auto& KEx = exListGen->braExAddress();
      
      // binding non excitation part address
      auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);
      const auto& LNonEx = nonExLooper->address();
      const auto& KNonEx = nonExLooper->auxAddress();
      size_t iThread = GetThreadID();
      size_t KBegin = nKExPerThread * iThread;
      size_t KEnd   = std::min(nKExDets, KBegin + nKExPerThread);
      auto& oneTDMSCR = (iThread == 0) ? oneTDM : SCR[iThread - 1];

      exListGen->visitExcitations(KBegin, KEnd,
          [&] (const auto& KExIter, const auto& pqExcitations) {
            for (const auto& [pp, qq, LEx, pqSign]: pqExcitations) {
              double fc = pqSign ? -1.: 1.;
              auto p = pp + pOrbOff;
              auto q = qq + qOrbOff;

              for (nonExLooper->setIndex(0ul);
                   not nonExLooper->isEnd();
                   nonExLooper->increment()) {
                 auto K = KEx + KNonEx;
                 auto L = LEx + LNonEx;
                 oneTDMSCR(p, q) += fc * SmartConj(CBra_ptr[K]) * CKet_ptr[L];   
              }
            }
          }
      );
    } // parallel region
  } // oneEExcitations
  
  // data reduction
  for (const auto& oneTDMSCR : SCR) oneTDM += oneTDMSCR;

} // DASCIBuilder::compute1TDMs

/*
 * compute twoTDM_IJ(t, w, u, v) = 
 *   \sum_{KL} conjugate(CI(K)) <K|EtuEwv|L>CJ(L) 
 */  
template <typename MatsT>
void DASCIBuilder<MatsT>::build2TDM(
    const LocalCIVectorsView<const MatsT>& CBra, 
    const LocalCIVectorsView<const MatsT>& CKet, 
    InCore4indexTPI<MatsT>& twoTDM,
    const double scale) const {
  
    // TODO: FIX this for taking 2e permutational symmetries
#ifdef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
  CErr("Need to fix 2tdm build for DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY");
#endif
  
  const auto& orbOffs = this->detFactory_.orbOffsInEachSpace();
   
  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();
  
  // TODO: break down the SCR storage in the future
  std::vector<InCore4indexTPI<MatsT>> SCR;
  for (auto i = 0ul; i < GetNumThreads() - 1; i++) {
    SCR.emplace_back(twoTDM.nBasis());
    SCR.back().clear();
  }
  
  for (const auto& twoEEx: this->detFactory_.twoEExcitations()) {
    // std::cout << " contraction on term - " << twoEEx.term << std::endl; 
    
    // MPI Parallelism
    if (not CKet.containsLocalCategory(twoEEx.categoricalIndices.second)) continue;
    
    const auto& catK = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(twoEEx.categoricalIndices.first));
    const auto& catL = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(twoEEx.categoricalIndices.second));
    const auto symmFact = twoEEx.symmetryFactor * scale; 
    const auto& exList_qp = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[0]); 
    const auto& exList_rs = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[1]); 
    
    const MatsT* CBra_ptr = CBra.getCategoryPointer(twoEEx.categoricalIndices.first); 
    const MatsT* CKet_ptr = CKet.getCategoryPointer(twoEEx.categoricalIndices.second);  

    #pragma omp parallel default(shared)
    {
      DoubleFullCD1eExListGenerator double1eExListsGen(exList_qp, exList_rs, twoEEx.exSpaces);
      const auto& exSpaces = double1eExListsGen.excitationSpaces();
      std::vector<size_t> KExOffs, KNonExOffs, LExOffs, LNonExOffs, dummy, nonExDims;
      catK.separateExAndNonExDimensions(exSpaces, KExOffs, dummy, KNonExOffs, nonExDims);
      catL.separateExAndNonExDimensions(exSpaces, LExOffs, dummy, LNonExOffs, dummy);

      auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);
      // binding non excitation part address
      const auto& LNonEx = nonExLooper->address();
      const auto& KNonEx = nonExLooper->auxAddress();
      
      const auto& pqrsSpaces = double1eExListsGen.pqrsSpaces();
      const auto pOrbOff = orbOffs[pqrsSpaces[0]];
      const auto qOrbOff = orbOffs[pqrsSpaces[1]];
      const auto rOrbOff = orbOffs[pqrsSpaces[2]];
      const auto sOrbOff = orbOffs[pqrsSpaces[3]];
      
      size_t iThread = GetThreadID();
      const size_t nJExDets = double1eExListsGen.totalDimension(); 
      const size_t nJExPerThread = std::ceil(double(nJExDets) / GetNumThreads());
      size_t JBegin = nJExPerThread * iThread;
      size_t JEnd   = std::min(nJExDets, JBegin + nJExPerThread);
      auto& twoTDMSCR = (iThread == 0) ? twoTDM : SCR[iThread - 1];
      
      double1eExListsGen.visitExcitations(JBegin, JEnd, LExOffs, KExOffs,
          [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {
            
            for (const auto& [qq, pp, KEx, pqSign]: qpExcitations) 
            for (const auto& [rr, ss, LEx, rsSign]: rsExcitations) {
              double fc = (pqSign == rsSign) ? 1. : -1.;
              auto p = pp + pOrbOff;
              auto q = qq + qOrbOff;
              auto r = rr + rOrbOff;
              auto s = ss + sOrbOff;
              
              for (nonExLooper->setIndex(0ul);
                   not nonExLooper->isEnd();
                   nonExLooper->increment()) {
                 auto K = KEx + KNonEx;
                 auto L = LEx + LNonEx;
                 twoTDMSCR(p, q, r, s) += fc * SmartConj(CBra_ptr[K]) * CKet_ptr[L];
              }
              
            }
          }
      ); // visitAllExcitations
    } // parallel region
  } // twoEExcitations
  
  // data reduction
  size_t N = twoTDM.nBasis();
  size_t N4 = N * N * N * N; 
  for (const auto& twoTDMSCR : SCR) {
   blas::axpy(N4, MatsT(1.), twoTDMSCR.pointer(), 1, twoTDM.pointer(), 1);
  }

} // DASCIBuilder::compute2TDMs

} // namespace ChronusQ

