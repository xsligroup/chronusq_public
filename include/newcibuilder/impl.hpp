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

  const auto& hCore_tt = *(moints_->template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
  const auto& antiSymmetricERI_ttuu = *(this->moints_->template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));
  
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


/*
 *  Find K largest diagonal elements on the fly
 *  TODO: make this stable_sort ordered. Easy to do by dividing the threads over CATEGORIES (rather than dets) and using stable_sort instead of nth element. Also use stable heaps
 */
#ifdef CQ_ENABLE_SPARSE
template <typename MatsT>
void NewCIBuilder<MatsT>::largestDiagH(DistributedSparseVectors<MatsT>& C,
    const CategoricalSpace& categoricalSpace, size_t K, std::vector<size_t>& kIndices, std::vector<MatsT>& kValues, bool (*comp)(const MatsT&,const MatsT&)) const {

  auto diagHSt = tick();

  std::vector<std::vector<std::pair<size_t, MatsT>>> localElements(GetNumThreads());
  for(size_t thrd = 0; thrd < GetNumThreads(); thrd++)
    localElements[thrd].reserve(K);
  
  const auto& hCore_tt = *(moints_->template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
  const auto& antiSymmetricERI_ttuu = *(this->moints_->template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));

  // main loop
  const auto& nOrbs = detFactory_.nOrbitalsInEachSpace();
  const size_t nCorrE = detFactory_.nTotalCorrElectrons();

  auto CLocalView = categoricalSpace.createLocalCIVectorsView(C, 0ul, 1ul);

  if(K == 1) {
    if(MPIRank(this->comm_) == 0) {
      std::vector<size_t> occ(nCorrE);
      const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*categoricalSpace.getCategory(0));
      const auto nEs = cat.SpaceOccupations();
      auto detCatGen = cat.template generator<uint64_t>();

      detCatGen.visitDeterminants(0, 1,
        [&] (size_t addr, const auto& dets)  {
          determinantsToOccs(dets, nEs, nOrbs, occ);
          MatsT tmp = MatsT(0.);
          for(const auto& t : occ) {
            tmp += hCore_tt(t, 0);
            for(const auto& u : occ) {
              tmp += MatsT(0.5) * antiSymmetricERI_ttuu(t, u);
          }
        }
        kIndices.push_back(0);
        kValues.push_back(tmp);
      });
    }
    else {
      kIndices.resize(1);
      kValues.resize(1);
    }
    MPI_Barrier(this->comm_);

    MPIBCast(&kIndices[0], 1, 0, this->comm_);
    MPIBCast(&kValues[0], 1, 0, this->comm_);

    return;
  }


  for (auto i = CLocalView.localCategoryBegin(); i < CLocalView.localCategoryEnd(); ++i) {
    const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*categoricalSpace.getCategory(i));
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

            //we want vector 'localElements' to store the K largest H values
            size_t index = addr + CLocalView.getCatOffest(i);
            MatsT diagHElement = tmp;

            if(localElements[GetThreadID()].size() < K) {
              localElements[GetThreadID()].push_back(std::make_pair(index, diagHElement));

	      if(localElements[GetThreadID()].size() == K) {
                std::make_heap(localElements[GetThreadID()].begin(), localElements[GetThreadID()].end(),
                    [&](const auto& i, const auto& j) { return comp(i.second, j.second); }
                );
              }
            }
            else
            {
              auto currMax = localElements[GetThreadID()].front();
	      //remove 'minimum' and insert the current diagonal elements, since it's 'larger'
              if(comp(diagHElement, currMax.second)) {
		std::pop_heap(localElements[GetThreadID()].begin(), localElements[GetThreadID()].end(),
                    [&](const auto& i, const auto& j) { return comp(i.second, j.second); }
                );
		localElements[GetThreadID()].pop_back();

	        localElements[GetThreadID()].push_back(std::make_pair(index, diagHElement));
	        std::push_heap(localElements[GetThreadID()].begin(), localElements[GetThreadID()].end(),
                    [&](const auto& i, const auto& j) { return comp(i.second, j.second); }
                );
              }
            }
          }
      );
    } // parallel region
  } // localCategory 
 
  size_t nLocalK = std::min(K, C.localLength());
 
  std::vector<std::pair<size_t, MatsT>> combinedLocalElements;
  for(size_t thrd = 0; thrd < GetNumThreads(); thrd++)
    combinedLocalElements.insert(combinedLocalElements.end(), localElements[thrd].begin(), localElements[thrd].end());

  std::nth_element(combinedLocalElements.begin(), combinedLocalElements.begin() + nLocalK, combinedLocalElements.end(),
    [&] (const auto& i, const auto& j) {
      return comp(i.second, j.second);
    });

  std::stable_sort(combinedLocalElements.begin(), combinedLocalElements.begin() + nLocalK,
    [&] (const auto& i, const auto& j) {
      return comp(i.second, j.second);
    });


  if (MPISize() == 1) {

    for (auto i = 0ul; i < nLocalK; ++i) {
      kIndices.push_back(combinedLocalElements[i].first);
      kValues.push_back(combinedLocalElements[i].second);
    } 
    return;
  }




  /*
   * MPI case
   */
       std::vector<size_t> kLocalIndices;
       std::vector<MatsT> kLocalValues;
       for (auto i = 0ul; i < nLocalK; ++i) {
         kLocalIndices.push_back(combinedLocalElements[i].first + C.localOffset());
         kLocalValues.push_back(combinedLocalElements[i].second);
       }

       // reduction
       size_t nResultK = std::min(K, C.length());

       // gather sizes
       std::vector<size_t> recv_sizes = MPIGather(nLocalK, 0, this->comm_);
       size_t totalGatheredSize = (MPIRank(this->comm_) == 0) ?
           std::accumulate(recv_sizes.begin(), recv_sizes.end(), size_t(0ul)) : 1ul;

       std::vector<size_t> gatheredIndices(totalGatheredSize);
       std::vector<MatsT> gatheredValues(totalGatheredSize);

       MPIGatherV(&kLocalIndices[0], nLocalK, &gatheredIndices[0], recv_sizes, 0, this->comm_);
       MPIGatherV(&kLocalValues[0], nLocalK, &gatheredValues[0], recv_sizes, 0, this->comm_);

       kIndices.resize(nResultK);
       kValues.resize(nResultK);

       if (MPIRank(this->comm_) == 0) {
         combinedLocalElements.resize(totalGatheredSize);
	 for(size_t i = 0; i < combinedLocalElements.size(); i++)
	   combinedLocalElements[i] = std::make_pair(i, MatsT(0.));

         std::nth_element(combinedLocalElements.begin(), combinedLocalElements.begin() + nResultK, combinedLocalElements.end(),
           [&] (const auto& i, const auto& j) {
             return comp(gatheredValues[i.first], gatheredValues[j.first]);
           });

	   std::stable_sort(combinedLocalElements.begin(), combinedLocalElements.begin() + nResultK,
           [&] (const auto& i, const auto& j) {
             return comp(gatheredValues[i.first], gatheredValues[j.first]);
           });

         for (auto i = 0ul; i < nResultK; ++i) {
           kIndices[i] = gatheredIndices[combinedLocalElements[i].first];
           kValues[i] = gatheredValues[combinedLocalElements[i].first];
         }
       }
       MPIBCast(&kIndices[0], nResultK, 0, this->comm_);
       MPIBCast(&kValues[0], nResultK, 0, this->comm_);

  if (MPIRank(this->comm_) == 0) {
    double diagHdur = tock(diagHSt);
    std::cout << "\nCIBuilder::DiagH - DURATION = " << std::setprecision(8)
              << diagHdur << " s." << std::endl;
  }

} // largestDiagH
#endif



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
    localC_MPIBuff1 = CQMemManager::get().malloc<MatsT>(nSCR_.at("MaxLocalKetNDim") * nVec);
    localC_MPIBuff2 = CQMemManager::get().malloc<MatsT>(nSCR_.at("MaxLocalKetNDim") * nVec);
  }
  
  size_t nTotalDet = detFactory_.braCategoricalSpace()->nDeterminants();
  size_t nContractionDone = 0;
  std::cout << std::endl << "      * Sigma Contraction Progress:" 
            << std::fixed << std::setprecision(2) << std::right << std::endl; 
  auto contraction_start = tick();
  
  double idle_time = 0.;

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
  
  if (localC_MPIBuff1) CQMemManager::get().free(localC_MPIBuff1);
  if (localC_MPIBuff2) CQMemManager::get().free(localC_MPIBuff2); 
  
  if (MPISize(comm_) > 1) { 
    std::cout << "      * Rank " << MPIRank(comm_) 
              << " idled (and waited for broadcast) " 
              << idle_time << " s" << std::endl; 
  }
#endif

  ProgramTimer::tock("NEWSigma");
}

#ifdef CQ_ENABLE_SPARSE
template <typename MatsT>
void NewCIBuilder<MatsT>::buildSigma(size_t nVec,
    const DistributedSparseVectors<MatsT> & C, size_t CShift,
    DistributedSparseVectors<MatsT>& Sigma, size_t SigmaShift, dcomplex* curEigenvalues, double eps) const {

  #ifndef CQ_ENABLE_MPI
  CErr("Can't have sparse dist vectors without MPI");
  #endif


  //std::cout << "nVec: " << nVec << " CShift: " << CShift << "sigmaShift: " << SigmaShift << std::endl;
  ProgramTimer::tick("NEWSigma");

  // Clear Sigma Vectors
  Sigma.clear(SigmaShift, nVec);

  auto SigmaView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(Sigma, SigmaShift, nVec);

  auto myCView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(C, CShift, nVec);

  myCView.toVecsByCat(nVec);


  size_t nTotalDet = detFactory_.braCategoricalSpace()->nDeterminants();
  size_t nContractionDone = 0;
  std::cout << std::endl << "      * Sigma Contraction Progress:"
            << std::fixed << std::setprecision(2) << std::right << std::endl;
  auto contraction_start = tick();

  double idle_time = 0.;

  // MPI parallelization here
  std::vector<size_t> nonZerosRanks(MPISize(comm_));
  size_t myNonZeros = C.getVecs().nonZeros(CShift, nVec);

  if (MPISize(comm_) > 1) {
    MPI_Allgather(&myNonZeros, 1, MPI_UINT64_T, nonZerosRanks.data(), 1, MPI_UINT64_T,  comm_);
  }
  else {
    nonZerosRanks[0] = myNonZeros;
  }

  size_t nonZerosRanksMax = *std::max_element(nonZerosRanks.begin(), nonZerosRanks.end());

  //Each ranks puts their own C local view here
  char* localC_MPIBuff0 = CQMemManager::get().malloc<char>(nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)));

  //convert sparse C into a flat buffer
  size_t* localC_MPIBuffRowCols = (size_t*)localC_MPIBuff0;
  MatsT* localC_MPIBuffVals = (MatsT*)(localC_MPIBuff0 + nonZerosRanks[MPIRank(comm_)] * sizeof(size_t) * 2);

  for (size_t col = 0; col < nVec; ++col) {
    #pragma omp parallel for schedule(static) default(shared)
    for(auto it = C.getVecs().cbegin(CShift + col); it != C.getVecs().cend(CShift + col); it++) {
      *(localC_MPIBuffRowCols + 2 * std::distance(C.getVecs().cbegin(CShift + col), it)) = it->first;
      *(localC_MPIBuffRowCols + 2 * std::distance(C.getVecs().cbegin(CShift + col), it) + 1) = col;
      *(localC_MPIBuffVals + std::distance(C.getVecs().cbegin(CShift + col), it)) = it->second;
    }
    localC_MPIBuffRowCols += 2 * C.getVecs().nonZeros(CShift + col);
    localC_MPIBuffVals += C.getVecs().nonZeros(CShift + col);
  }

  //used for one-sided communication
  char* localC_MPIBuff1 = nullptr;
  char* localC_MPIBuff2 = nullptr;
  char* bcast_ptr = nullptr;

  if (MPISize(comm_) > 1) {
    localC_MPIBuff1 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
    localC_MPIBuff2 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
  }

  MPIWin localC_Win(localC_MPIBuff0, nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)), comm_);
  localC_Win.lock_all(0);

  //for reduction
  HashSparseMatrix<MatsT> SCRSigma_hash(nVec);
  
  for(size_t j = 0; j < nVec; j++) {
    for(auto it = C.getVecs().cbegin(CShift + j); it != C.getVecs().cend(CShift + j); it++) {
      SCRSigma_hash.update(it->first, j, - it->second * std::real(curEigenvalues[j]));
    }
  }
 

  // only the first broad cast is blocking
  for (int iBCast = 0; iBCast < MPISize(comm_) + 1; ++iBCast) {

    bcast_ptr = iBCast % 2 == 0 ? localC_MPIBuff1 : localC_MPIBuff2;
    char* bcast_ptr_prev = iBCast % 2 == 0 ? localC_MPIBuff2 : localC_MPIBuff1;

    if (iBCast == 0) {
      if (MPISize(comm_) > 1 and MPIRank(comm_) != 0) {
	localC_Win.get(bcast_ptr, nonZerosRanks[0] * (sizeof(MatsT) + sizeof(size_t) * 2), 0, 0);
      }
      continue;
    }

    // send broadcast requests for next
    if (iBCast < MPISize(comm_) and MPIRank(comm_) != iBCast) {
      ProgramTimer::tick("Sigma MPI COMM Init IBCast");
      localC_Win.get(bcast_ptr, nonZerosRanks[iBCast] * (sizeof(MatsT) + sizeof(size_t) * 2), iBCast, 0);
      ProgramTimer::tock("Sigma MPI COMM Init IBCast");
    }

    
    if (MPIRank(comm_) != iBCast - 1) {
      //This means we are waiting for the first MPI communication
      if(iBCast == 1) {
        ProgramTimer::tick("Sigma MPI COMM First BCast");
        localC_Win.flush(iBCast - 1);
        ProgramTimer::tock("Sigma MPI COMM First BCast");
      }
      else {
	ProgramTimer::tick("Sigma MPI COMM Root Wait");
        auto idle_start = tick();
        localC_Win.flush(iBCast - 1);
        idle_time += tock(idle_start);
        ProgramTimer::tock("Sigma MPI COMM Root Wait");
      }
    }
    

    ProgramTimer::tick("Sigma MPI COMM Create CView");
    //convert flattened buffer to SparseMatrix
    
    char* vecsBuffer = (MPISize(comm_) > 1 and MPIRank(comm_) != iBCast - 1) ? bcast_ptr_prev : localC_MPIBuff0;
    auto CView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView<MatsT>(vecsBuffer, nonZerosRanks[iBCast - 1], nVec, comm_, iBCast - 1);

    ProgramTimer::tock("Sigma MPI COMM Create CView");


    ProgramTimer::tick("Sigma Contraction");
    if ( nonZerosRanks[iBCast - 1] > 0) {
      for(size_t j = 0; j < nVec; j++) 
        SCRSigma_hash.reserve(CView.nonZeros(j), j);

      buildSigma1e(CView, myCView, SigmaView, SCRSigma_hash, curEigenvalues, eps);
      //buildSigma1e(CView, SigmaView, SCRSigma_hash);
      buildSigma2e(CView, myCView, SigmaView, SCRSigma_hash, curEigenvalues, eps);

    }
    ProgramTimer::tock("Sigma Contraction");


    nContractionDone += C.lengthAtNode(iBCast - 1);
    auto ratio = double(nContractionDone) / double(nTotalDet);
    auto dur_contraction = tock(contraction_start);

    std::cout << "        - completed "
              << std::setw(6) << 100 * ratio
              << "%, time elapsed " << std::setw(10) <<  dur_contraction << " s"
              << ", approximatedly needs " << std::setw(10)
              << dur_contraction / ratio - dur_contraction  <<  " s more"
              << std::endl;

  } // Broadcast C


  if (localC_MPIBuff1) CQMemManager::get().free(localC_MPIBuff1);
  if (localC_MPIBuff2) CQMemManager::get().free(localC_MPIBuff2);

  SigmaView.reduceSigmaIntermediate(SCRSigma_hash); 

  if (MPISize(comm_) > 1) {
    std::cout << "      * Rank " << MPIRank(comm_)
              << " waiting for other ranks to finish... "
              << std::endl;
  }

  ProgramTimer::tick("Sigma MPI COMM Root Wait");
  auto idle_start = tick();
  localC_Win.unlock_all();
  localC_Win.free();
  idle_time += tock(idle_start);
  ProgramTimer::tock("Sigma MPI COMM Root Wait");

  if (MPISize(comm_) > 1) {
    std::cout << "      * Rank " << MPIRank(comm_)
              << " idled (and waited for broadcast) "
              << idle_time << " s" << std::endl;
  }
  
  CQMemManager::get().free(localC_MPIBuff0);

  ProgramTimer::tock("NEWSigma");

}// NewCIBuilder::buildSigma
















template <typename MatsT>
void NewCIBuilder<MatsT>::formSubA(size_t nVec, size_t nTot,
    const DistributedSparseVectors<MatsT> & C, size_t CShift,
    std::vector<MatsT>& VAV) const {

#ifndef CQ_ENABLE_MPI
  CErr("Can't have sparse dist vectors without MPI");
#endif


  ProgramTimer::tick("NEWSigma");

  size_t nTotalDet = detFactory_.braCategoricalSpace()->nDeterminants();
  size_t nContractionDone = 0;
  std::cout << std::endl << "      * Sigma Contraction Progress:"
            << std::fixed << std::setprecision(2) << std::right << std::endl;
  auto contraction_start = tick();

  double idle_time = 0.;

  auto myCView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(C, 0, nTot);

  myCView.toVecsByCat(nTot);

  // MPI parallelization here
  std::vector<size_t> nonZerosRanks(MPISize(comm_));
  size_t myNonZeros = C.getVecs().nonZeros(CShift, nVec);

  if (MPISize(comm_) > 1) {
    MPI_Allgather(&myNonZeros, 1, MPI_UINT64_T, nonZerosRanks.data(), 1, MPI_UINT64_T,  comm_);
  }
  else {
    nonZerosRanks[0] = myNonZeros;
  }

  size_t nonZerosRanksMax = *std::max_element(nonZerosRanks.begin(), nonZerosRanks.end());

  //Each ranks puts their own C local view here
  char* localC_MPIBuff0 = CQMemManager::get().malloc<char>(nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)));

  //convert sparse C into a flat buffer
  size_t* localC_MPIBuffRowCols = (size_t*)localC_MPIBuff0;
  MatsT* localC_MPIBuffVals = (MatsT*)(localC_MPIBuff0 + nonZerosRanks[MPIRank(comm_)] * sizeof(size_t) * 2);

  for (size_t col = 0; col < nVec; ++col) {
    #pragma omp parallel for schedule(static) default(shared)
    for(auto it = C.getVecs().cbegin(CShift + col); it != C.getVecs().cend(CShift + col); it++) {
      *(localC_MPIBuffRowCols + 2 * std::distance(C.getVecs().cbegin(CShift + col), it)) = it->first;
      *(localC_MPIBuffRowCols + 2 * std::distance(C.getVecs().cbegin(CShift + col), it) + 1) = col;
      *(localC_MPIBuffVals + std::distance(C.getVecs().cbegin(CShift + col), it)) = it->second;
    }
    localC_MPIBuffRowCols += 2 * C.getVecs().nonZeros(CShift + col);
    localC_MPIBuffVals += C.getVecs().nonZeros(CShift + col);
  }

  //used for one-sided communication
  char* localC_MPIBuff1 = nullptr;
  char* localC_MPIBuff2 = nullptr;
  char* bcast_ptr = nullptr;

  if (MPISize(comm_) > 1) {
    localC_MPIBuff1 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
    localC_MPIBuff2 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
  }

  MPIWin localC_Win(localC_MPIBuff0, nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)), comm_);
  localC_Win.lock_all(0);


  // only the first broad cast is blocking
  for (int iBCast = 0; iBCast < MPISize(comm_) + 1; ++iBCast) {

    bcast_ptr = iBCast % 2 == 0 ? localC_MPIBuff1 : localC_MPIBuff2;
    char* bcast_ptr_prev = iBCast % 2 == 0 ? localC_MPIBuff2 : localC_MPIBuff1;

    if (iBCast == 0) {
      if (MPISize(comm_) > 1 and MPIRank(comm_) != 0) {
	localC_Win.get(bcast_ptr, nonZerosRanks[0] * (sizeof(MatsT) + sizeof(size_t) * 2), 0, 0);
      }
      continue;
    }

    // send broadcast requests for next
    if (iBCast < MPISize(comm_) and MPIRank(comm_) != iBCast) {
      ProgramTimer::tick("Sigma MPI COMM Init IBCast");
      localC_Win.get(bcast_ptr, nonZerosRanks[iBCast] * (sizeof(MatsT) + sizeof(size_t) * 2), iBCast, 0);
      ProgramTimer::tock("Sigma MPI COMM Init IBCast");
    }


    if (MPIRank(comm_) != iBCast - 1) {
      //This means we are waiting for the first MPI communication
      if(iBCast == 1) {
        ProgramTimer::tick("Sigma MPI COMM First BCast");
        localC_Win.flush(iBCast - 1);
        ProgramTimer::tock("Sigma MPI COMM First BCast");
      }
      else {
	ProgramTimer::tick("Sigma MPI COMM Root Wait");
        auto idle_start = tick();
        localC_Win.flush(iBCast - 1);
        idle_time += tock(idle_start);
        ProgramTimer::tock("Sigma MPI COMM Root Wait");
      }
    }


    ProgramTimer::tick("Sigma MPI COMM Create CView");
    //convert flattened buffer to SparseMatrix

    char* vecsBuffer = (MPISize(comm_) > 1 and MPIRank(comm_) != iBCast - 1) ? bcast_ptr_prev : localC_MPIBuff0;
    auto CView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView<MatsT>(vecsBuffer, nonZerosRanks[iBCast - 1], nVec, comm_, iBCast - 1);

    ProgramTimer::tock("Sigma MPI COMM Create CView");


    ProgramTimer::tick("Sigma Contraction");
    //if ( nonZerosRanks[iBCast - 1] > 0) {
      //formSubA1e(CView, myCView, VAV);
      //formSubA2e(CView, myCView, VAV);

    //}
    for(size_t jVec = 0; jVec < nTot; jVec++) {
      for(size_t iVec = 0; iVec < nVec; iVec++) {
      if (C.getVecs().nonZeros(jVec) > 0 and nonZerosRanks[iBCast - 1] > 0) {
	formSubA1e(CView, myCView, VAV, iVec, jVec);
	formSubA2e(CView, myCView, VAV, iVec, jVec);
      }
      }
    }
    ProgramTimer::tock("Sigma Contraction");


    nContractionDone += C.lengthAtNode(iBCast - 1);
    auto ratio = double(nContractionDone) / double(nTotalDet);
    auto dur_contraction = tock(contraction_start);

    std::cout << "        - completed "
              << std::setw(6) << 100 * ratio
              << "%, time elapsed " << std::setw(10) <<  dur_contraction << " s"
              << ", approximatedly needs " << std::setw(10)
              << dur_contraction / ratio - dur_contraction  <<  " s more"
              << std::endl;

  } // Broadcast C


  if (localC_MPIBuff1) CQMemManager::get().free(localC_MPIBuff1);
  if (localC_MPIBuff2) CQMemManager::get().free(localC_MPIBuff2);


  if (MPISize(comm_) > 1) {
    std::cout << "      * Rank " << MPIRank(comm_)
              << " waiting for other ranks to finish... "
              << std::endl;
  }

  ProgramTimer::tick("Sigma MPI COMM Root Wait");
  auto idle_start = tick();
  localC_Win.unlock_all();
  localC_Win.free();
  idle_time += tock(idle_start);
  ProgramTimer::tock("Sigma MPI COMM Root Wait");

  if (MPISize(comm_) > 1) {
    std::cout << "      * Rank " << MPIRank(comm_)
              << " idled (and waited for broadcast) "
              << idle_time << " s" << std::endl;
  }

  CQMemManager::get().free(localC_MPIBuff0);

  ProgramTimer::tock("NEWSigma");

}
#endif









/*
 *  Default Build Full H hack through buildSigma function, can be overridden
 */ 
template <typename MatsT>
void NewCIBuilder<MatsT>::buildFullH(
  DistributedVectors<MatsT>& fullH) const { 
  
  const auto ketCategoricalSpace = detFactory_.ketCategoricalSpace();
  size_t nKetDets = ketCategoricalSpace->nDeterminants();
  
  auto IdenM = ketCategoricalSpace->template constructDistributedCIVectors<MatsT>(comm_, nKetDets);

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
        oneTDM->nRows());
    *reducedOneTDM = *oneTDM;
    oneTDM.swap(reducedOneTDM);
  }
  
  std::shared_ptr<InCore4indexTPI<MatsT>> reducedTwoTDM = nullptr;
  if (twoTDM and reduceTwoTDM) {
    // put the results to reduceTwoTDM instead
    reducedTwoTDM = std::make_shared<InCore4indexTPI<MatsT>>(
        twoTDM->nBasis());
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
    localCKetJ_MPIBuff1 = CQMemManager::get().malloc<MatsT>(nSCR_.at("MaxLocalKetNDim"));
    localCKetJ_MPIBuff2 = CQMemManager::get().malloc<MatsT>(nSCR_.at("MaxLocalKetNDim"));
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
  if (localCKetJ_MPIBuff1) CQMemManager::get().free(localCKetJ_MPIBuff1);
  if (localCKetJ_MPIBuff2) CQMemManager::get().free(localCKetJ_MPIBuff2); 
  
  // data reductions
  if (oneTDM and reduceOneTDM) {
    size_t reductionDim = oneTDM->nRows();
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





/*
 * Build TDM function wrappers, with MPI parallelization
 */
#ifdef CQ_ENABLE_SPARSE
template <typename MatsT>
void NewCIBuilder<MatsT>::buildTDM(const DistributedSparseVectors<MatsT>& CBra,
    const DistributedSparseVectors<MatsT>& CKet,
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
        oneTDM->dimension());
    *reducedOneTDM = *oneTDM;
    oneTDM.swap(reducedOneTDM);
  }

  std::shared_ptr<InCore4indexTPI<MatsT>> reducedTwoTDM = nullptr;
  if (twoTDM and reduceTwoTDM) {
    // put the results to reduceTwoTDM instead
    reducedTwoTDM = std::make_shared<InCore4indexTPI<MatsT>>(
        twoTDM->nBasis());
    *reducedTwoTDM = *twoTDM;
    twoTDM.swap(reducedTwoTDM);
  }
#endif

  auto CBraIView = detFactory_.braCategoricalSpace()
      ->createLocalCIVectorsView(CBra, stateIndices.first, 1ul);
  //splits CI vector to CSC storage per category
  CBraIView.toVecsByCat(1ul);

  // MPI parallelization here
#ifdef CQ_ENABLE_MPI
  std::vector<size_t> nonZerosRanks(MPISize(comm_));
  size_t myNonZeros = CKet.getVecs().nonZeros(stateIndices.second, 1);

  if (MPISize(comm_) > 1) {
    MPI_Allgather(&myNonZeros, 1, MPI_UINT64_T, nonZerosRanks.data(), 1, MPI_UINT64_T,  comm_);
  }
  else {
    nonZerosRanks[0] = myNonZeros;
  }

  size_t nonZerosRanksMax = *std::max_element(nonZerosRanks.begin(), nonZerosRanks.end());

  //Each ranks puts their own CKet local view here
  char* localC_MPIBuff0 = CQMemManager::get().malloc<char>(nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)));

  //convert sparse C into a flat buffer
  size_t* localC_MPIBuffRowCols = (size_t*)localC_MPIBuff0;
  MatsT* localC_MPIBuffVals = (MatsT*)(localC_MPIBuff0 + nonZerosRanks[MPIRank(comm_)] * sizeof(size_t) * 2);

  #pragma omp parallel for schedule(static) default(shared)
  for(auto it = CKet.getVecs().cbegin(stateIndices.second); it != CKet.getVecs().cend(stateIndices.second); it++) {
    *(localC_MPIBuffRowCols + 2 * std::distance(CKet.getVecs().cbegin(stateIndices.second), it)) = it->first;
    *(localC_MPIBuffRowCols + 2 * std::distance(CKet.getVecs().cbegin(stateIndices.second), it) + 1) = 0;
    *(localC_MPIBuffVals + std::distance(CKet.getVecs().cbegin(stateIndices.second), it)) = it->second;
  }

  //used for one-sided communication
  char* localC_MPIBuff1 = nullptr;
  char* localC_MPIBuff2 = nullptr;
  char* bcast_ptr = nullptr;

  if (MPISize(comm_) > 1) {
    localC_MPIBuff1 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
    localC_MPIBuff2 = CQMemManager::get().malloc<char>(nonZerosRanksMax * (sizeof(MatsT) + 2 * sizeof(size_t)));
  }

  MPIWin localC_Win(localC_MPIBuff0, nonZerosRanks[MPIRank(comm_)] * (sizeof(MatsT) + 2 * sizeof(size_t)), comm_);
  localC_Win.lock_all(0);

  for (int iBCast = 0; iBCast < MPISize(comm_) + 1; ++iBCast) {

    bcast_ptr = iBCast % 2 == 0 ? localC_MPIBuff1 : localC_MPIBuff2;
    char* bcast_ptr_prev = iBCast % 2 == 0 ? localC_MPIBuff2 : localC_MPIBuff1;

    if (iBCast == 0) {
      if (MPISize(comm_) > 1 and MPIRank(comm_) != 0) {
	localC_Win.get(bcast_ptr, nonZerosRanks[0] * (sizeof(MatsT) + sizeof(size_t) * 2), 0, 0);
      }
      continue;
    }

    // send broadcast requests for next
    if (iBCast < MPISize(comm_) and MPIRank(comm_) != iBCast) {
      localC_Win.get(bcast_ptr, nonZerosRanks[iBCast] * (sizeof(MatsT) + sizeof(size_t) * 2), iBCast, 0);
    }

    if (MPIRank(comm_) != iBCast - 1) {
      //This means we are waiting for the first MPI communication
      if(iBCast == 1) {
        localC_Win.flush(iBCast - 1);
      }
      else {
        localC_Win.flush(iBCast - 1);
      }
    }
     
    char* vecsBuffer = (MPISize(comm_) > 1 and MPIRank(comm_) != iBCast - 1) ? bcast_ptr_prev : localC_MPIBuff0;
    auto CKetJView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView<MatsT>(vecsBuffer, nonZerosRanks[iBCast - 1], 1, comm_, iBCast - 1);

    if (oneTDM and nonZerosRanks[iBCast - 1] > 0) {
      build1TDM(CBraIView, CKetJView, *oneTDM, oneTDMScale);
    }

    if (twoTDM and nonZerosRanks[iBCast - 1] > 0) {
      build2TDM(CBraIView, CKetJView, *twoTDM, twoTDMScale);
    }

#else
    auto CKetJView = detFactory_.ketCategoricalSpace()
        ->createLocalCIVectorsView(CKet, stateIndices.second, 1ul);

    if (oneTDM) {
      build1TDM(CBraIView, CKetJView, *oneTDM, oneTDMScale);
    }

    if (twoTDM) {
      build2TDM(CBraIView, CKetJView, *twoTDM, twoTDMScale);
    }

#endif

#ifdef CQ_ENABLE_MPI

  } // Broadcast C

  if (localC_MPIBuff1) CQMemManager::get().free(localC_MPIBuff1);
  if (localC_MPIBuff2) CQMemManager::get().free(localC_MPIBuff2);

  if (MPISize(comm_) > 1) {
    std::cout << "      * Rank " << MPIRank(comm_)
              << " waiting for other ranks to finish... " 
	      << std::endl;
  }

  localC_Win.unlock_all();
  localC_Win.free();
  CQMemManager::get().free(localC_MPIBuff0);


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
#endif

} // namespace ChronusQ
