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
#include <newcibuilder/dascibuilder/defaultbuildtdms.hpp>
#include <newcibuilder/dascibuilder/defaultbuildsigma.hpp>
#include <newcibuilder/dascibuilder/knowleshandy.hpp>
#include <newcibuilder/dascibuilder/olsenroos.hpp>
#include <newcibuilder/dascibuilder/frischli.hpp>
#include <newcibuilder/dascibuilder/smallblockhamiltonian.hpp>
#include <newcibuilder/dascibuilder/largeblockhamiltonian.hpp>
#include <util/scratch.hpp>

namespace ChronusQ {

template <typename MatsT>
void DASCIBuilder<MatsT>::setSigma2eContractionAlgorithm(std::string alg) {
  if (alg == "NAIVE" or alg == "NAIVELOOP" or alg == "NL") {
    twoEcontAlg_ = DASCISigma2eContAlg::NAIVE;
  } else if (alg == "KNOWLESHANDY" or alg == "KH" ) {
    twoEcontAlg_ = DASCISigma2eContAlg::KNOWLESHANDY;  
  } else if (alg == "OLSENROOS" or alg == "OR") {
    twoEcontAlg_ = DASCISigma2eContAlg::OLSENROOS;  
  } else if (alg == "FRISCHLI" or alg == "FL") {
    twoEcontAlg_ = DASCISigma2eContAlg::FRISCHLI;  
  } else if (alg == "SMALLBLOCKH" or alg == "SBH" 
      or alg == "SMALLBLOCKHAMILTONIAN") {
    twoEcontAlg_ = DASCISigma2eContAlg::SMALLBLOCKHAMILTONIAN;
  } else if (alg == "LARGEBLOCKH" or alg == "LBH" 
      or alg == "LARGEBLOCKHAMILTONIAN") {
    twoEcontAlg_ = DASCISigma2eContAlg::LARGEBLOCKHAMILTONIAN;
  } else {
    CErr("Unknown Sigma 2e Contraction Algorithm in DASCI");
  }
  this->estimateLocalMemoryForMPIComm();  
  estimateMemoryForIntermediates();
}

/*
 * Determine the memory upper bound of scratching space
 */ 
template <typename MatsT>
void DASCIBuilder<MatsT>::estimateMemoryForIntermediates() {
  
  std::vector<size_t> KExDims, LExDims, dummy, nonExDims;
  size_t nSCR1 = 0ul, nSCR2 = 0ul, nSCR3 = 0ul, nSCRSigma = 0ul; 
  
  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();
  
  for (const auto& twoEEx: this->detFactory_.twoEExcitations()) {
    
    const auto& catK = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(twoEEx.categoricalIndices.first));

    nSCRSigma = std::max(nSCRSigma, catK.nDeterminants());
    
    if (twoEcontAlg_ == DASCISigma2eContAlg::NAIVE) continue; 
    
    const auto& catL = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(twoEEx.categoricalIndices.second));
    const auto& exList_qp = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[0]);
    const auto& exList_rs = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[1]);
    
    DoubleFullCD1eExListGenerator double1eExListsGen(exList_qp, exList_rs, twoEEx.exSpaces);
    const auto& exSpaces = double1eExListsGen.excitationSpaces();
    catK.separateExAndNonExDimensions(exSpaces, dummy, KExDims, dummy, nonExDims);
    catL.separateExAndNonExDimensions(exSpaces, dummy, LExDims, dummy, dummy);
    
    size_t nNonExDets = 1ul;
    for (const auto& n : nonExDims) nNonExDets *= n;
    
    size_t nNZqp = exList_qp.nNonZeroExcitations();
    size_t nNZrs = exList_rs.nNonZeroExcitations();
     
    switch (twoEcontAlg_) {
      case DASCISigma2eContAlg::KNOWLESHANDY:
      case DASCISigma2eContAlg::OLSENROOS:
      case DASCISigma2eContAlg::FRISCHLI:
      case DASCISigma2eContAlg::SMALLBLOCKHAMILTONIAN:
        // for DASCISigma2eContAlg::KNOWLESHANDY, Omega->SCR1, Lambda->SCR2, X->SCR3
        // for DASCISigma2eContAlg::OLSENROOS, Omega->SCR1, sSCR->SCR2, X->SCR3 
        // for DASCISigma2eContAlg::FRISCHLI,  cSCR->SCR1, Lambda->SCR2, X->SCR3
        // for DASCISigma2eContAlg::SMALLBLOCKHAMILTONIAN,  cSCR->SCR1, sSCR->SCR2, X->SCR3
        nSCR1 = std::max(nSCR1, nNZrs * nNonExDets);
        nSCR2 = std::max(nSCR2, nNZqp * nNonExDets);
        nSCR3 = std::max(nSCR3, nNZrs * nNZqp);
        break;
      case DASCISigma2eContAlg::LARGEBLOCKHAMILTONIAN:
        // for DASCISigma2eContAlg::LARGEBLOCKHAMILTONIAN,cSCR->SCR1, sSCR->SCR2, X->SCR3
        {
          size_t nKEx = catK.nDeterminants() / nNonExDets;
          size_t nLEx = catL.nDeterminants() / nNonExDets;
          nSCR1 = std::max(nSCR1, catL.nDeterminants());
          nSCR2 = std::max(nSCR2, catK.nDeterminants());
          nSCR3 = std::max(nSCR3, nKEx * nLEx);
        }
        break;
      default:
        continue;

    } // switch for contration orders

  } // twoEExcitations  
  
  
  this->nSCR_.insert_or_assign("SCR1", nSCR1);
  this->nSCR_.insert_or_assign("SCR2", nSCR2);
  this->nSCR_.insert_or_assign("SCR3", nSCR3);
  this->nSCR_.insert_or_assign("Sigma", nSCRSigma);
  

  std::cout << " - Intermediates in DASCIBuilder Needs "
            << ((nSCR1 + nSCR2 + nSCR3 + nSCRSigma) / 1e9 ) * sizeof(MatsT) * GetNumThreads() 
            << " GB";
  if (static_cast<size_t>(twoEcontAlg_) < 5) {
    std::cout << " per vector";
  }

}

/*
 * Loop to build Sigma 2e Part for every CatK <- CatJ <- CatL
 * 
 * exList_qp: CatJ <- Ex(q, p) <- CatK
 * exList_rs: CatJ <- Ex(r, s) <- CatL
 */
template <typename MatsT>
void DASCIBuilder<MatsT>::buildSigma2e(
    const LocalCIVectorsView<const MatsT>& C,
    const LocalCIVectorsView<MatsT>& Sigma) const {
  
  // std::cout << " twoEcontAlg_ = " << twoEcontAlg_ << std::endl;
  
  assert(C.size() == Sigma.size());
  size_t nVec = C.size();
  
  // allocate intermediates memory
  SharedMemoryScratch<MatsT> SCR1;
  SharedMemoryScratch<MatsT> SCR2;
  SharedMemoryScratch<MatsT> SCR3;
  SharedMemoryScratch<MatsT> SCRSigma;
  auto nSCR1 = this->nSCR_.at("SCR1");
  auto nSCR2 = this->nSCR_.at("SCR2");
  auto nSCR3 = this->nSCR_.at("SCR3");
  auto nSCRSigma = this->nSCR_.at("Sigma");
  if (static_cast<size_t>(twoEcontAlg_) < 5) {
    nSCR1 *= nVec; 
    nSCR2 *= nVec;
    nSCR3 *= nVec;
  }
  nSCRSigma *= nVec;

  SCR1.resize(nSCR1, "SCR1");
  SCR2.resize(nSCR2, "SCR2");
  SCR3.resize(nSCR3, "SCR3");
  SCRSigma.resize(nSCRSigma, "SCRSigma");

  auto nLAThreads = GetLAThreads();
  bool initializeSCRSigma = true;  
  const auto& twoEExcitations = this->detFactory_.twoEExcitations();

  const auto braCategoricalSpace = this->detFactory_.braCategoricalSpace();
  const auto ketCategoricalSpace = this->detFactory_.ketCategoricalSpace();
  
  // pull out the i incrementation to later
  size_t i = 0ul;
  
  // find first valid one
  for(; i < twoEExcitations.size(); ++i) {
    if (C.containsLocalCategory(
        twoEExcitations[i].categoricalIndices.second)) break;
  }

  while (i < twoEExcitations.size()) {
    const auto& twoEEx = twoEExcitations[i];

    // std::cout << " contraction on term - " << twoEEx.term << std::endl; 
    // twoEEx.output(std::cout);
    // find the interaction terms
    const auto& s2e = *this->moints_->template getIntegral<DASTwoPInts, MatsT>(twoEEx.term);
    // s2e.output(std::cout, twoEEx.term, true);
     
    const auto& braCategory = dynamic_cast<const FullDeterminantCategory&>(
        *braCategoricalSpace->getCategory(twoEEx.categoricalIndices.first));
    const auto& ketCategory = dynamic_cast<const FullDeterminantCategory&>(
        *ketCategoricalSpace->getCategory(twoEEx.categoricalIndices.second));
    const auto& exList_qp = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[0]);
    const auto& exList_rs = dynamic_cast<const FullCD1eExList&>(*twoEEx.exLists[1]);
    
    const MatsT* C_ptr = C.getCategoryPointer(twoEEx.categoricalIndices.second); 
    const auto nC = C.localLength();
    const auto nBraCatDets = braCategory.nDeterminants();

    SetLAThreads(1);

    #pragma omp parallel default(shared)
    {      
      std::vector<size_t> KExOffs, KNonExOffs, KExDims, LExOffs, LNonExOffs, LExDims, dummy, nonExDims;
      DoubleFullCD1eExListGenerator double1eExListsGen(exList_qp, exList_rs, twoEEx.exSpaces);
      const auto& exSpaces = double1eExListsGen.excitationSpaces();
      braCategory.separateExAndNonExDimensions(exSpaces, KExOffs, KExDims, KNonExOffs, nonExDims);
      ketCategory.separateExAndNonExDimensions(exSpaces, LExOffs, LExDims, LNonExOffs, dummy);
      
      auto ketNonExLooper = constructTensorLooper(nonExDims, LNonExOffs);
      auto braNonExLooper = constructTensorLooper(nonExDims, KNonExOffs);
      
      MatsT* Sigma_ptr = SCRSigma.getPtr(); 
      if (initializeSCRSigma) std::fill_n(Sigma_ptr, nBraCatDets * nVec, MatsT(0.));

      switch (twoEcontAlg_) {
        case DASCISigma2eContAlg::KNOWLESHANDY:
          DASCISigma2eBuilder::buildKnowlesHandy(
                  nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                  KExOffs, braNonExLooper, LExOffs, ketNonExLooper,
                  SCR1.getPtr(), SCR2.getPtr(), SCR3.getPtr(),
              twoEEx.symmetryFactor * 0.5);
          break;
        case DASCISigma2eContAlg::OLSENROOS:
          DASCISigma2eBuilder::buildOlsenRoos(
                  nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                  KExOffs, braNonExLooper, LExOffs, ketNonExLooper,
                  SCR1.getPtr(), SCR2.getPtr(), SCR3.getPtr(),
              twoEEx.symmetryFactor * 0.5);
          break;
        case DASCISigma2eContAlg::FRISCHLI:
          DASCISigma2eBuilder::buildFrischLi(
                  nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                  KExOffs, braNonExLooper, LExOffs, ketNonExLooper,
                  SCR1.getPtr(), SCR2.getPtr(), SCR3.getPtr(),
              twoEEx.symmetryFactor * 0.5);
          break;
        case DASCISigma2eContAlg::SMALLBLOCKHAMILTONIAN:
          DASCISigma2eBuilder::buildSmallBlockHamiltonian(
                  nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                  KExOffs, braNonExLooper, LExOffs, ketNonExLooper,
                  SCR1.getPtr(), SCR2.getPtr(), SCR3.getPtr(),
              twoEEx.symmetryFactor * 0.5);
          break;
        case DASCISigma2eContAlg::LARGEBLOCKHAMILTONIAN:
          {
            auto KExLooper = constructTensorLooper(KExDims, KExOffs);
            auto LExLooper = constructTensorLooper(LExDims, LExOffs);
            DASCISigma2eBuilder::buildLargeBlockHamiltonian(
                    nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                    KExLooper, braNonExLooper, LExLooper, ketNonExLooper,
                    SCR1.getPtr(), SCR2.getPtr(), SCR3.getPtr(),
                twoEEx.symmetryFactor * 0.5);
          }
          break;
        default:
          {
            auto nonExLooper = constructTensorLooper(nonExDims, LNonExOffs, KNonExOffs);
            DASCISigma2eBuilder::buildNaive(
                    nVec, C_ptr, nC, Sigma_ptr, nBraCatDets, s2e, double1eExListsGen,
                    KExOffs, LExOffs, nonExLooper,
                twoEEx.symmetryFactor * 0.5);
          }
      } // switch for contraction orders
    } // parallel region 
    
    SetLAThreads(nLAThreads);
    
    // increment on i
    for (++i; i < twoEExcitations.size(); ++i) {
      if (C.containsLocalCategory(
          twoEExcitations[i].categoricalIndices.second)) break;

    }
   
    if (i >= twoEExcitations.size() or 
        twoEExcitations[i].categoricalIndices.first != twoEEx.categoricalIndices.first) {
      for (auto j = 0ul; j < nVec; ++j) { 
        MatsT* Sigma_ptr = Sigma.getCategoryPointer(twoEEx.categoricalIndices.first, j);
        for (auto k = 0ul; k < GetNumThreads(); ++k) {
          blas::axpy(nBraCatDets, MatsT(1.), SCRSigma.getPtr(k) + nBraCatDets * j, 1, Sigma_ptr, 1);
        }
      }
      initializeSCRSigma = true;
    } else {
      initializeSCRSigma = false;
    } // Data Reduction when needed
  
  } // twoEExcitations

} // buildSigma2e

} // namespace ChronusQ
