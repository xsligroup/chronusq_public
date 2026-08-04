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
#ifdef CQ_ENABLE_SPARSE
#include <newcibuilder/dascibuilder.hpp>
#include <detfactory/excitationlist.hpp>

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <util/scratch.hpp>

#include <itersolver/cqSparseMatrix.hpp>

// #define _DEBUG_CIBuilder_IMPL

namespace ChronusQ {

/*
 *  Default Sparse Build Sigma one-electron part for DASCI
 */
template <typename MatsT>
void DASCIBuilder<MatsT>::formSubA2e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC,
      std::vector<MatsT>& VAV, size_t iVec, size_t jVec) const {

  //number of new vectors
  size_t nVec = C.size();

  //total number of vectors, nVec + nPrev
  size_t nTot = myC.size();

  const auto& twoEExcitations = this->detFactory_.twoEExcitations();

  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();


  //std::unordered_set<size_t> mySet;

    #pragma omp parallel for schedule(dynamic) default(shared)
    for (const auto& twoEEx : twoEExcitations) {

      // skip invalid categories
      if ((not C.containsLocalCategory(twoEEx.categoricalIndices.second)) or C.getVecsByCat()[twoEEx.categoricalIndices.second - C.localCategoryBegin()].nonZeros(iVec) == 0 or myC.getVecsByCat()[twoEEx.categoricalIndices.first - myC.localCategoryBegin()].nonZeros(jVec) == 0) {
      continue;
      }

      const auto& s2e = *this->moints_->template getIntegral<DASTwoPInts, MatsT>(twoEEx.term);

      const auto& braCategory = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(twoEEx.categoricalIndices.first));
      const auto& ketCategory = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(twoEEx.categoricalIndices.second));
      const auto& exList_qp = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[0]);
      const auto& exList_rs = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[1]);

      std::vector<size_t> KExOffs, KNonExOffs, KExDims, LExOffs, LNonExOffs, LExDims, dummy, nonExDims;
      DoubleFullCD1eExListGenerator double1eExListsGen(exList_qp, exList_rs, twoEEx.exSpaces);
      const auto& exSpaces = double1eExListsGen.excitationSpaces();
      braCategory.separateExAndNonExDimensions(exSpaces, KExOffs, KExDims, KNonExOffs, nonExDims);
      ketCategory.separateExAndNonExDimensions(exSpaces, LExOffs, LExDims, LNonExOffs, dummy);

      auto ketNonExLooper = constructTensorLooper(nonExDims, LNonExOffs);
      auto braNonExLooper = constructTensorLooper(nonExDims, KNonExOffs);
      auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);

      DASCISigma2eBuilder::formSubANaive(VAV,
        nVec, C.getVecsByCat()[twoEEx.categoricalIndices.second - C.localCategoryBegin()], C.getShift(),
	nTot, myC.getVecsByCat()[twoEEx.categoricalIndices.first - myC.localCategoryBegin()], myC.getShift(),
        s2e, double1eExListsGen,
        KExOffs, LExOffs, nonExLooper, twoEEx.symmetryFactor * 0.5, iVec, jVec
      );

    } // twoEExcitations



}







/*
 *  Default Sparse Build Sigma two-electron part for DASCI
 */
template <typename MatsT>
void DASCISigma2eBuilder::formSubANaive(std::vector<MatsT>& VAV,
    size_t nVec, const LLSparseMatrix<MatsT>& C, size_t shiftC, 
    size_t nTot, const LLSparseMatrix<MatsT>& myC, size_t shiftMyC,
    const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
    const std::vector<size_t>& KExOffs, const std::vector<size_t>& LExOffs,
    std::shared_ptr<TensorLooper>& nonExLooper, const double symmFact, size_t iVec, size_t jVec) {

  // binding non excitation part address
  const auto& LNonEx = nonExLooper->address();
  const auto& KNonEx = nonExLooper->auxAddress();

  //Excitation visitor with dynamic load balancing
  double1eExListsGen.visitAllSparseExcitations(LExOffs, KExOffs,
      [&] (const auto& JExIter, const auto& qpExcitations, const auto& rsExcitations) {

        //for (auto iVec = 0ul; iVec < nVec; ++iVec) {

          size_t nonZeroCs = C.nonZeros(shiftC + iVec);
          //if(nonZeroCs == 0) {
            //return;//continue;
          //}

          for (const auto& [r, s, LEx, rsSign]: rsExcitations) {

            //Need to find intersection of C[iVec] and LEx + LNonEx
            size_t nonExIdx = 0;
            nonExLooper->setIndex(nonExIdx);
            size_t CIdx = 0;
            auto itC = C.cbegin(shiftC + iVec);
            size_t nonExDets = nonExLooper->nTotal();

            while(nonExIdx < nonExDets and CIdx < nonZeroCs) {
              auto L = LEx + LNonEx;
              size_t deltaC = 1, deltaNonEx = 1;

              while(CIdx < nonZeroCs and (itC + CIdx)->first < L) {
                CIdx += deltaC;
                deltaC *= 2;
              }
              CIdx -= deltaC / 2;

              while(nonExIdx < nonExDets  and L < (itC + CIdx)->first) {
                nonExIdx += deltaNonEx;
                deltaNonEx *= 2;
                if(nonExIdx < nonExDets) {
                  nonExLooper->setIndex(nonExIdx);
                  L = LEx + LNonEx;
                }
              }
              nonExIdx -= deltaNonEx / 2;
              if(deltaNonEx >= 2) {
                nonExLooper->setIndex(nonExIdx);
                L = LEx + LNonEx;
              }

              if((itC + CIdx)->first < L) {
                //binary search for the right index
                //auto lower_bound_it = std::lower_bound((itC + CIdx), (itC + std::min(CIdx + deltaC, nonZeroCs)), {L, MatsT(0.)}, [](const auto& a, const auto& b){return a.first < b.first;});

                size_t mid;
                size_t low = CIdx;
                const size_t max = std::min(CIdx + deltaC / 2, nonZeroCs);
                size_t high = max;

                while (low < high) {
                  mid = low + (high - low) / 2;

                  if (L <= (itC + mid)->first) {
                    high = mid;
                  }
                  else {
                    low = mid + 1;
                  }
                }

                if(low < max and (itC + low)->first < L) {
                  low++;
                }

                CIdx = low;

              }
              else {
                if((itC + CIdx)->first > L) {
                  //binary search for the right index
                  size_t mid;
                  size_t low = nonExIdx;
                  const size_t max = std::min(nonExIdx + deltaNonEx / 2, nonExDets);
                  size_t high = max;

                  while (low < high) {
                    mid = low + (high - low) / 2;
                    nonExLooper->setIndex(mid);

                    if (LEx + LNonEx >= (itC + CIdx)->first) {
                      high = mid;
                    }
                    else {
                      low = mid + 1;
                    }
                  }

                  nonExLooper->setIndex(low);

                  if(low < max and (itC + CIdx)->first > LEx + LNonEx) {
                    low++;
                  }

                  nonExIdx = low;
                  nonExLooper->setIndex(nonExIdx);
                }
                else {
                  for (const auto& [q, p, KEx, pqSign]: qpExcitations) {
                    MatsT h2e_pqrs = (pqSign == rsSign) ? symmFact * s2e(p, q, r, s) : -symmFact * s2e(p, q, r, s);
                    if(h2e_pqrs == MatsT(0.)) {continue;}
                    auto K = KEx + KNonEx;
		    //for(size_t jVec = 0; jVec < nTot; jVec++) {
		      //(VAV.getPtr())[iVec * nTot + jVec] += SmartConj(myC.get(K, jVec)) * ((itC + CIdx)->second) * h2e_pqrs;
		      VAV[iVec * nTot + jVec + GetThreadID() * nTot * nVec] += SmartConj(myC.get(K, jVec)) * ((itC + CIdx)->second) * h2e_pqrs;
		    //}
                    //Sigma.update(K + Sigma_iCatOffset, iVec, ((itC + CIdx)->second) * h2e_pqrs);

                  }
                  CIdx++;
                  nonExIdx++;
                  nonExLooper->increment();
                }
              }
            }
          }
        //}
     }
  ); // visitAllExcitations

} // DASCISigma2eBuilder::formSubANaive







/*
 *  Default Sparse Build Sigma one-electron part for DASCI
 */
template <typename MatsT>
void DASCIBuilder<MatsT>::formSubA1e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC,
      std::vector<MatsT>& VAV, size_t iVec, size_t jVec) const {
      
  
  //number of new vectors
  size_t nVec = C.size();

  //total number of vectors, nVec + nPrev
  size_t nTot = myC.size();


  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();

  const auto& oneEExcitations = this->detFactory_.oneEExcitations();

  #pragma omp parallel for schedule(dynamic) default(shared)
  for(const auto& oneEEx : oneEExcitations) {

    if ((not C.containsLocalCategory(oneEEx.categoricalIndices.second)) or C.getVecsByCat()[oneEEx.categoricalIndices.second - C.localCategoryBegin()].nonZeros(iVec) == 0 or myC.getVecsByCat()[oneEEx.categoricalIndices.first - myC.localCategoryBegin()].nonZeros(jVec) == 0) {
      continue;
    }


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


    // multidimensional index advancer
    auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);
    size_t nNonExDets = nonExLooper->nTotal();

    // bind non-excitation part addresses
    const auto& ketNonEx = nonExLooper->address();
    const auto& braNonEx = nonExLooper->auxAddress();

    auto exListGen = exList.generator(ketExOffs, braExOffs);
    // bind braEx with advancer called inside visitExcitations
    const auto& braEx = exListGen->braExAddress();

    exListGen->visitAllSparseExcitations(
      [&] (const auto& braExIter, const auto& pqExcitations) {

        //for (auto iVec = 0ul; iVec < nVec; ++iVec) {

          size_t nonZeroCs = C.getVecsByCat()[oneEEx.categoricalIndices.second - C.localCategoryBegin()].nonZeros(C.getShift() + iVec);
          //if(nonZeroCs == 0) {
            //continue;
          //}

          for (const auto& [p, q, ketEx, pqSign]: pqExcitations) {

            // between-space symmetry factor multiplied by the between-orbital symmetry factor
            MatsT h1e_pq = pqSign ? -symmFact * h1e(p, q) : symmFact * h1e(p, q);
            if(h1e_pq == MatsT(0.)) {continue;}

            //Need to find intersection of C[iVec] and LEx + LNonEx
            size_t nonExIdx = 0;
            nonExLooper->setIndex(nonExIdx);
            size_t CIdx = 0;
            auto itC = C.getVecsByCat()[oneEEx.categoricalIndices.second - C.localCategoryBegin()].cbegin(C.getShift() + iVec);
            size_t nonExDets = nonExLooper->nTotal();

            while(nonExIdx < nonExDets and CIdx < nonZeroCs) {
              auto L = ketEx + ketNonEx;
              size_t deltaC = 1, deltaNonEx = 1;

              while(CIdx < nonZeroCs and (itC + CIdx)->first < L) {
                CIdx += deltaC;
                deltaC *= 2;
              }
              CIdx -= deltaC / 2;

              while(nonExIdx < nonExDets and L < (itC + CIdx)->first) {
                nonExIdx += deltaNonEx;
                deltaNonEx *= 2;
                if(nonExIdx < nonExDets) {
                  nonExLooper->setIndex(nonExIdx);
                  L = ketEx + ketNonEx;
                }
              }
              nonExIdx -= deltaNonEx / 2;
              if(deltaNonEx >= 2) {
                nonExLooper->setIndex(nonExIdx);
                L = ketEx + ketNonEx;
              }

              if((itC + CIdx)->first < L) {
                //binary search for the right index
                //auto lower_bound_it = std::lower_bound((itC + CIdx), (itC + std::min(CIdx + deltaC, nonZeroCs)), {L, MatsT(0.)}, [](const auto& a, const auto& b){return a.first < b.first;});

                size_t mid;
                size_t low = CIdx;
                const size_t max = std::min(CIdx + deltaC / 2, nonZeroCs);
                size_t high = max;

                while (low < high) {
                  mid = low + (high - low) / 2;

                  if (L <= (itC + mid)->first) {
                    high = mid;
                  }
                  else {
                    low = mid + 1;
                  }
                }

                if(low < max and (itC + low)->first < L) {
                  low++;
                }

                CIdx = low;

              }
              else {
                if((itC + CIdx)->first > L) {
                  //binary search for the right index
                  size_t mid;
                  size_t low = nonExIdx;
                  const size_t max = std::min(nonExIdx + deltaNonEx / 2, nonExDets);
                  size_t high = max;

                  while (low < high) {
                    mid = low + (high - low) / 2;
                    nonExLooper->setIndex(mid);

                    if (ketEx + ketNonEx >= (itC + CIdx)->first) {
                      high = mid;
                    }
                    else {
                      low = mid + 1;
                    }
                  }

                  nonExLooper->setIndex(low);

                  if(low < max and (itC + CIdx)->first > ketEx + ketNonEx) {
                    low++;
                  }

                  nonExIdx = low;
                  nonExLooper->setIndex(nonExIdx);
                }
                else {
                  auto K = braEx + braNonEx;

		  //for(size_t jVec = 0; jVec < nTot; jVec++) {
			  
		    auto& myCForThisCat = myC.getVecsByCat()[oneEEx.categoricalIndices.first - myC.localCategoryBegin()];
		    VAV[iVec * nTot + jVec + GetThreadID() * nTot * nVec] += SmartConj(myCForThisCat.get(K, jVec)) * ((itC + CIdx)->second) * h1e_pq;
                    //(VAV.getPtr())[iVec * nTot + jVec] += SmartConj(myCForThisCat.get(K, jVec)) * ((itC + CIdx)->second) * h1e_pq;
                  //}

                  //SCRSigma_hash.update(K + Sigma_iCatOffset, iVec, ((itC + CIdx)->second) * h1e_pq);
                  CIdx++;
                  nonExIdx++;
                  nonExLooper->increment();
                }
              }
            }
          }
        //}
      }
    );
  } // oneEExcitations
} // DASCIBuilder::formSubA1e





} //namespace 
#endif
