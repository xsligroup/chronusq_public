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

#include <newcibuilder.hpp>
#include <newcibuilder/dascibuilder/impl.hpp>

namespace ChronusQ {

template <typename MatsT>
void NewCIBuilder<MatsT>::estimateLocalMemoryForMPIComm() {

  const auto& localCLens = detFactory_.ketCategoricalSpace()
      ->distributedCategoryLengths();
  nSCR_.insert_or_assign("MaxLocalKetNDim", *std::max_element(localCLens.begin(), localCLens.end()));
  std::cout << " - To BroadCast Ket CI Vectors,  Needs "
            << (nSCR_.at("MaxLocalKetNDim") / 1e9 ) * sizeof(MatsT) 
            << " GB per vector" << std::endl;
  std::cout << std::endl;
}

/*
 *  Default Build Diagonal H, can be overridden
 */ 
template <typename MatsT>
void NewCIBuilder<MatsT>::buildDiagH(DistributedVectors<MatsT>& diagH, 
    const CategoricalSpace& categoricalSpace) const {
  
  auto diagHSt = tick();

  const auto& hCore_tt = *(moints_.template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
  const auto& antiSymmetricERI_ttuu = *(this->moints_.template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));
  
  // std::cout << " hCore_tt_Correlated_Space dim = " << hCore_tt.nBasis1() << std::endl;
  // std::cout << " antiSymmetricERI_ttuu dim = " << antiSymmetricERI_ttuu.nBasis() <<  std::endl;

  // main loop
  const auto& nOrbs = detFactory_.nOrbitalsInEachSpace();
  const size_t nCorrE = detFactory_.nTotalCorrElectrons();
  
  auto diagHLocalView = categoricalSpace.createLocalCIVectorsView(diagH, 0ul, 1ul); 

  for (auto i = diagHLocalView.localCategoryBegin(); i < diagHLocalView.localCategoryEnd(); ++i) {
    const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*categoricalSpace.getCategory(i));
    MatsT *catDiagH = diagHLocalView.getCategoryPointer(i); 
    const auto nEs = cat.SpaceOccupations();
    const size_t nDets = cat.nDeterminants();
    const size_t nDetsPerThread = std::ceil(double(nDets) / GetNumThreads());
    
    #pragma omp parallel default(shared)
    {
      std::vector<size_t> occ(nCorrE);
      auto detCatGen = cat.template generator<uint64_t>();
      size_t iBegin = nDetsPerThread * GetThreadID();
      size_t iEnd   = std::min(nDets, iBegin + nDetsPerThread);

      detCatGen.visitDeterminants(iBegin, iEnd,
          [&] (size_t addr, const auto& dets)  { 
            determinantsToOccs(dets, nEs, nOrbs, occ);
            MatsT tmp = MatsT(0.);
            for(const auto& t : occ) {
              tmp += hCore_tt(t, 0);
              for(const auto& u : occ) {
                tmp += MatsT(0.5) * antiSymmetricERI_ttuu(t, u);
              }
            }
            catDiagH[addr] = tmp;
          }
      );
    } // parallel region
  } // localCategory 

  double diagHdur = tock(diagHSt);
  std::cout << "\nCIBuilder::DiagH - DURATION = " << std::setprecision(8)
            << diagHdur << " s." << std::endl;

} // NewCIBuilder::buildDiagH

template <typename MatsT>
void NewCIBuilder<MatsT>::buildSigma(size_t nVec, 
    const DistributedVectors<MatsT> & C, size_t CShift,
    DistributedVectors<MatsT>& Sigma, size_t SigmaShift) const {
  
  ProgramTimer::tick("NEWSigma");

  // Clear Sigma Vectors
   
  // auto clear_start = tick();
  MatsT *sigma_ptr = Sigma.getLocalPtr(SigmaShift);
  size_t clear_len = Sigma.localLength() * nVec;
  size_t clear_len_per_thread = std::ceil(double(clear_len) / GetNumThreads());
  #pragma omp parallel default(shared)
  {
    size_t clear_begin = clear_len_per_thread * GetThreadID();
    size_t clear_end = std::min(clear_len, clear_begin + clear_len_per_thread);
    std::fill(sigma_ptr + clear_begin, sigma_ptr + clear_end, MatsT(0.));
  }
  //auto dur_clear = tock(clear_start);
  //std::cout << std::endl << "      * Initialized " << nVec << " Sigma Vectors, used " 
  //          << dur_clear << " s" << std::endl;
  
  auto SigmaView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(Sigma, SigmaShift, nVec);
  
  // MPI parallelization here
#ifdef CQ_ENABLE_MPI 
  
  MatsT* localC_MPIBuff1 = nullptr;
  MatsT* localC_MPIBuff2 = nullptr;
  const MatsT* localC = nullptr; 
  MatsT* bcast_ptr = nullptr;
  size_t bcast_len = 0ul;
  std::vector<MPI_Request> bcast_req;
  std::vector<MPI_Status> bcast_status;
  
  if (MPISize(comm_) > 1) {
    localC_MPIBuff1 = memManager_.malloc<MatsT>(nSCR_.at("MaxLocalKetNDim") * nVec);
    localC_MPIBuff2 = memManager_.malloc<MatsT>(nSCR_.at("MaxLocalKetNDim") * nVec);
  }
  
  size_t nTotalDet = detFactory_.braCategoricalSpace()->nDeterminants();
  size_t nContractionDone = 0;
  std::cout << std::endl << "      * Sigma Contraction Progress:" 
            << std::fixed << std::setprecision(2) << std::right << std::endl; 
  auto contraction_start = tick();
  
  double idle_time = 0.;

  // only the first broad cast is blocking
  for (int iBCast = 0; iBCast < MPISize(comm_) + 1; ++iBCast) {
    if (iBCast == 0) { 
      if (MPISize(comm_) > 1) {
        ProgramTimer::tick("Sigma MPI COMM First BCast");
        bcast_ptr = MPIRank(comm_) == 0 ? const_cast<MatsT*>(C.getLocalPtr(CShift)) : localC_MPIBuff1;  
        bcast_len = C.lengthAtNode(iBCast) * nVec;
        MPIBCast(bcast_ptr, bcast_len, 0, comm_);
        ProgramTimer::tock("Sigma MPI COMM First BCast");
      } else {
        bcast_ptr = const_cast<MatsT*>(C.getLocalPtr(CShift));
      }
      continue;
    } 
    
    ProgramTimer::tick("Sigma MPI COMM Create CView");
    localC = bcast_ptr;
    auto CView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView(localC, nVec, comm_, iBCast - 1);
    ProgramTimer::tock("Sigma MPI COMM Create CView");
    
    // send broadcast requests for next
    if (iBCast < MPISize(comm_)) {
      ProgramTimer::tick("Sigma MPI COMM Init IBCast");
      if (MPIRank(comm_) == iBCast) {
        bcast_ptr = const_cast<MatsT*>(C.getLocalPtr(CShift));
      } else {
        bcast_ptr = iBCast % 2 == 0 ? localC_MPIBuff1 : localC_MPIBuff2;
      }
      bcast_len = C.lengthAtNode(iBCast) * nVec;
      // std::cout << "HHDebug Start MPI IBCast, bcast_len = " << bcast_len
      //           << ", mem is " << sizeof(MatsT) * bcast_len / 1e9 << " GB" << std::endl;
      bcast_req = MPIIBCast(bcast_ptr, bcast_len, iBCast, comm_);
      ProgramTimer::tock("Sigma MPI COMM Init IBCast");
    }

#else
  auto CView = detFactory_.ketCategoricalSpace()
      ->createLocalCIVectorsView(C, CShift, nVec);
#endif

    ProgramTimer::tick("Sigma Contraction");
    buildSigma1e(CView, SigmaView);
    buildSigma2e(CView, SigmaView); 
    ProgramTimer::tock("Sigma Contraction");
  
#ifdef CQ_ENABLE_MPI 

    nContractionDone += C.lengthAtNode(iBCast - 1);
    auto ratio = double(nContractionDone) / double(nTotalDet);
    auto dur_contraction = tock(contraction_start);

    std::cout << "        - completed " 
              << std::setw(6) << 100 * ratio 
              << "%, time elapsed " << std::setw(10) <<  dur_contraction << " s" 
              << ", approximatedly needs " << std::setw(10) 
              << dur_contraction / ratio - dur_contraction  <<  " s more" 
              << std::endl;

    if (iBCast < MPISize(comm_)) {
      ProgramTimer::tick("Sigma MPI COMM Root Wait");
      auto idle_start = tick();
      bcast_status = MPIWait(bcast_req);
      idle_time += tock(idle_start);      
      ProgramTimer::tock("Sigma MPI COMM Root Wait");
    }
  } // Broadcast C 
  
  if (localC_MPIBuff1) memManager_.free(localC_MPIBuff1);
  if (localC_MPIBuff2) memManager_.free(localC_MPIBuff2); 
  
  if (MPISize(comm_) > 1) { 
    std::cout << "      * Rank " << MPIRank(comm_) 
              << " idled (and waited for broadcast) " 
              << idle_time << " s" << std::endl; 
  }
#endif

  ProgramTimer::tock("NEWSigma");
}

/*
 *  Default Build Full H hack through buildSigma function, can be overridden
 */ 
template <typename MatsT>
void NewCIBuilder<MatsT>::buildFullH(
  DistributedVectors<MatsT>& fullH) const { 
  
  const auto ketCategoricalSpace = detFactory_.ketCategoricalSpace();
  size_t nKetDets = ketCategoricalSpace->nDeterminants();
  
  auto IdenM = ketCategoricalSpace->template constructDistributedCIVectors<MatsT>(comm_, memManager_, nKetDets);

  IdenM->clear();
 
  for(auto i = 0ul; i < nKetDets; ++i) IdenM->set(i, i, MatsT(1.));
  
  buildSigma(nKetDets, *IdenM, 0ul, fullH, 0ul);
  
} // NewCIBuilder::buildFullH

/*
 * Build TDM function wrappers, with MPI parallelization
 */ 
template <typename MatsT>
void NewCIBuilder<MatsT>::buildTDM(const DistributedVectors<MatsT>& CBra, 
    const DistributedVectors<MatsT>& CKet, 
    const std::pair<size_t, size_t>& stateIndices,  
    std::shared_ptr<cqmatrix::Matrix<MatsT>> oneTDM,
    bool oneTDMIncrement, double oneTDMScale, 
    std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM,
    bool twoTDMIncrement, double twoTDMScale
#ifdef CQ_ENABLE_MPI 
    , bool reduceOneTDM, bool reduceTwoTDM 
#endif
    ) const {

  if (not oneTDM and not twoTDM) return;

  if (oneTDM and not oneTDMIncrement) oneTDM->clear();
  if (twoTDM and not twoTDMIncrement) twoTDM->clear();

#ifdef CQ_ENABLE_MPI 
  if (MPISize(comm_) <= 1) {      
    reduceOneTDM = false;
    reduceTwoTDM = false;
  }

  std::shared_ptr<cqmatrix::Matrix<MatsT>> reducedOneTDM = nullptr;
  if (oneTDM and reduceOneTDM) {
    // put the results to reduceOneTDM instead
    reducedOneTDM = std::make_shared<cqmatrix::Matrix<MatsT>>(
        memManager_, oneTDM->dimension());
    *reducedOneTDM = *oneTDM;
    oneTDM.swap(reducedOneTDM);
  }
  
  std::shared_ptr<InCore4indexTPI<MatsT>> reducedTwoTDM = nullptr;
  if (twoTDM and reduceTwoTDM) {
    // put the results to reduceTwoTDM instead
    reducedTwoTDM = std::make_shared<InCore4indexTPI<MatsT>>(
        memManager_, twoTDM->nBasis());
    *reducedTwoTDM = *twoTDM;
    twoTDM.swap(reducedTwoTDM);
  }
#endif
  
  auto CBraIView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(CBra, stateIndices.first, 1ul);
  
  // MPI parallelization here
#ifdef CQ_ENABLE_MPI 
  MatsT* localCKetJ_MPIBuff1 = nullptr;
  MatsT* localCKetJ_MPIBuff2 = nullptr;
  const MatsT* localCKetJ = nullptr; 
  MatsT* bcast_ptr = nullptr;
  size_t bcast_len = 0ul;
  std::vector<MPI_Request> bcast_req;
  std::vector<MPI_Status> bcast_status;
  
  if (MPISize(comm_) > 1) {
    localCKetJ_MPIBuff1 = memManager_.malloc<MatsT>(nSCR_.at("MaxLocalKetNDim"));
    localCKetJ_MPIBuff2 = memManager_.malloc<MatsT>(nSCR_.at("MaxLocalKetNDim"));
  }
  
  // only the first broad cast is blocking
  for (int iBCast = 0; iBCast < MPISize(comm_) + 1; ++iBCast) {
  
    if (iBCast == 0) { 
      if (MPISize(comm_) > 1) {
        if (MPIRank(comm_) == 0) {
          bcast_ptr = const_cast<MatsT*>(CKet.getLocalPtr(stateIndices.second));
        } else {
          bcast_ptr = localCKetJ_MPIBuff1;
        }
        bcast_len = CKet.lengthAtNode(iBCast);
        MPIBCast(bcast_ptr, bcast_len, 0, comm_);
      } else {
        bcast_ptr = const_cast<MatsT*>(CKet.getLocalPtr(stateIndices.second));
      }
      continue;
    }  
    
    localCKetJ = bcast_ptr;
    auto CKetJView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView(localCKetJ, 1ul, comm_, iBCast - 1);
      
    // send broadcast requests for next
    if (iBCast < MPISize(comm_)) {
      if (MPIRank(comm_) == iBCast) {
        bcast_ptr = const_cast<MatsT*>(CKet.getLocalPtr(stateIndices.second));
      } else {
        bcast_ptr = iBCast % 2 == 0 ? localCKetJ_MPIBuff1 : localCKetJ_MPIBuff2;
      }
      bcast_len = CKet.lengthAtNode(iBCast);
      bcast_req = MPIIBCast(bcast_ptr, bcast_len, iBCast, comm_);
    }
      
#else
    auto CKetJView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView(CKet, stateIndices.second, 1ul);
#endif
  
    if (oneTDM) {
      build1TDM(CBraIView, CKetJView, *oneTDM, oneTDMScale);
    }

    if (twoTDM) {
      build2TDM(CBraIView, CKetJView, *twoTDM, twoTDMScale);
    }

#ifdef CQ_ENABLE_MPI 
    if (iBCast < MPISize(comm_)) {
      bcast_status = MPIWait(bcast_req);
    }
  } // Broadcast C 
  if (localCKetJ_MPIBuff1) memManager_.free(localCKetJ_MPIBuff1);
  if (localCKetJ_MPIBuff2) memManager_.free(localCKetJ_MPIBuff2); 
  
  // data reductions
  if (oneTDM and reduceOneTDM) {
    size_t reductionDim = oneTDM->dimension();
    reductionDim *= reductionDim;
    MPIAllReduce(oneTDM->pointer(), reductionDim, reducedOneTDM->pointer(), comm_);
  }

  if (twoTDM and reduceTwoTDM) {
    size_t reductionDim = twoTDM->nBasis();
    reductionDim *= reductionDim;
    reductionDim *= reductionDim;
    MPIAllReduce(twoTDM->pointer(), reductionDim, reducedTwoTDM->pointer(), comm_);
  }
#endif

} // NewCIBuilder::buildTDM


} // namespace ChronusQ
