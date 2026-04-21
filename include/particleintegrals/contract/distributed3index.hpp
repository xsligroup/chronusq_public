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


#include <integrals.hpp>
#include <util/matout.hpp>
#include <util/timer.hpp>
#include <util/threads.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/gradints/incore.hpp>

//#define _REPORT_COMM_TIMINGS

#ifdef _REPORT_COMM_TIMINGS
#define TIMED_COMM(label, call)                                            \
  {                                                                        \
    auto __start = tick();                                                 \
    call;                                                                  \
    double __elapsed = tock(__start);                                      \
    std::cout << "Rank " << rank << " :: " << label                        \
              << " = " << __elapsed << " s" << std::endl;                  \
  }
#else
#define TIMED_COMM(label, call) call
#endif

constexpr size_t KCONTRACTSPLITNB_NBUFFERS_ = 2;

namespace ChronusQ {


  /**
   *  \brief Perform a Coulomb-type (34,12) RI-ERI contraction with
   *  a one-body operator.
   *  J_{pq} = ∑_α ∑_{rs} L^α_{pq} · L^α_{rs} · P_{rs}  
   */   
  template<typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT,IntsT>::JContract
        (MPI_Comm comm, TwoBodyContraction<MatsT>& C) const {

    InCoreRITPI<IntsT> &ritpi = *std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_);
    auto eri3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ritpi.eri3j());
    if (eri3j == nullptr)
      CErr("DistributedRITPIContraction::JContract expect a DistributedERI3J integral object");

    if (eri3j->distributionLayout() == DistributionLayout::SplitNB)
      return JContractSplitNB(comm, C, eri3j);
    else
      return JContractSplitNBRI(comm, C, eri3j);

  }; // DistributedRITPIContraction::JContract



  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::JContractSplitNB(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const {
    auto contractStart = tick();
    int rank = 0, size = 1;
#ifdef CQ_ENABLE_MPI
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);
    size_t NB = eri3j->nBasis();
    size_t NB2 = NB*NB;
    size_t NBRI = eri3j->nRIBasis();
    size_t localNB = eri3j->localSize();
    size_t localNBStart = eri3j->localStart(); 

    // Prepare density and result matrices
    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    IntsT *AX = reinterpret_cast<IntsT*>(C.AX);
    // Extract the real part of X if X is Hermetian and if the ints are real
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;
    // Allocate scratch if IntsT and MatsT are different
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;
    if( extractRealPartX ) {
      X = CQMemManager::get().malloc<IntsT>(NB2);
      for(auto k = 0ul; k < NB2; k++) X[k] = std::real(C.X[k]);
    }
    if( allocAXScratch )
      AX = CQMemManager::get().malloc<IntsT>(NB2);
    
    // Step 0: Obtain density for local rows D[r_local][s]
    auto X_local = CQMemManager::get().malloc<IntsT>(localNB*NB);
    size_t idx = 0;
    for (size_t i = 0; i < localNB; i++) {
      size_t globalNB = localNBStart + i;
      for (size_t s = 0; s < NB; s++) {
        X_local[idx++] = X[globalNB + s*NB];
      }
    }

    // Set up scratch space
    auto T_local = CQMemManager::get().malloc<IntsT>(NBRI);
    auto T_global = CQMemManager::get().malloc<IntsT>(NBRI);
    auto AX_local = CQMemManager::get().malloc<IntsT>(localNB*NB);
    std::fill_n(T_local, NBRI, IntsT(0.));
    std::fill_n(T_global, NBRI, IntsT(0.));
    std::fill_n(AX_local, localNB*NB, IntsT(0.));
    std::fill_n(AX,NB2,IntsT(0.));

    // Step 1: Compute T[α] = ∑_{r_local,s} L[α][s][r_local] * D[r_local][s]
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
               NBRI, 1, NB*localNB,
               IntsT(1.), eri3j->data(), NBRI,
               X_local, NB*localNB,
               IntsT(0.), T_local, NBRI);

    // Step 2: Sum T[α] over all ranks 
    TIMED_COMM("JContraction MPIAllReduce",
      MPIAllReduce(T_local, static_cast<int>(NBRI), T_global, comm)
    );

    // Step 3: Compute J_local[p_local][q] = ∑_α L[α][q][p_local] * T[α]
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
               NB * localNB, 1, NBRI,
               IntsT(1.), eri3j->data(), NBRI,
               T_global, NBRI,
               IntsT(0.), AX_local, NB * localNB);

    // Step 4: Assemble global J matrix
    std::vector<int> recvcounts(size), displs(size);
    for (int r = 0; r < size; ++r) {
      size_t r_localNB = eri3j->splitSize(NB, r, size);
      displs[r] = eri3j->splitStart(NB, r, size) * NB;
      recvcounts[r] = static_cast<int>(r_localNB * NB);
    } 
    IntsT *bufPtr  = (rank == 0) ? AX : nullptr;
    int *recvPtr   = (rank == 0) ? recvcounts.data() : nullptr;
    int *displsPtr = (rank == 0) ? displs.data() : nullptr;
    TIMED_COMM("JContraction MPI_Gatherv",
      MPI_Gatherv(AX_local, localNB * NB, mpi_data_type<IntsT>(),
                bufPtr, recvPtr, displsPtr, mpi_data_type<IntsT>(),
                0, comm)
    );

    // Cleanup temporary storage
    CQMemManager::get().free(X_local, T_local, T_global, AX_local);
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {
      std::copy_n(AX,NB2,C.AX);
      CQMemManager::get().free(AX);
    }

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: JContraction (SplitNB) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedRITPIContraction::JContractSplitNB called without MPI support");
#endif

  }; // DistributedRITPIContraction::JContractSplitNB



  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::JContractSplitNBRI(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const {
    auto contractStart = tick();
    int rank = 0, size = 1;
#ifdef CQ_ENABLE_MPI
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);
    size_t NB = eri3j->nBasis();
    size_t NB2 = NB*NB;
    size_t NBRI = eri3j->nRIBasis();
    size_t nQ = eri3j->localSize();

    // Prepare density and result matrices
    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    IntsT *AX = reinterpret_cast<IntsT*>(C.AX);
    // Extract the real part of X if X is Hermetian and if the ints are real
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;
    // Allocate scratch if IntsT and MatsT are different
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;
    if( extractRealPartX ) {
      X = CQMemManager::get().malloc<IntsT>(NB2);
      for(auto k = 0ul; k < NB2; k++) X[k] = std::real(C.X[k]);
    }
    if( allocAXScratch )
      AX = CQMemManager::get().malloc<IntsT>(NB2);

    // Set up scratch space
    auto T_local = CQMemManager::get().malloc<IntsT>(nQ);
    auto AX_local = CQMemManager::get().malloc<IntsT>(NB*NB);
    std::fill_n(T_local, nQ, IntsT(0.));
    std::fill_n(AX_local, NB2, IntsT(0.));
    std::fill_n(AX,NB2,IntsT(0.));

    // Step 1: Compute T[α_local] = ∑_{r,s} L[α_local][s][r] * D[r][s]
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nQ,1,NB2,IntsT(1.),eri3j->data(),nQ,X,NB2,IntsT(0.),T_local,nQ);

    // Step 2: Compute J_local[p][q] = ∑_α L[α][q][p] * T[α]
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB2,1,nQ,IntsT(1.),eri3j->data(),nQ,T_local,nQ,IntsT(0.),AX_local,NB2);
    
    // Step 3: Sum up J_local[p][q] over all ranks 
    TIMED_COMM("JContraction MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB2), AX, 0, comm)
    );
    
    // Cleanup temporary storage
    CQMemManager::get().free(T_local, AX_local);
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {
      std::copy_n(AX,NB2,C.AX);
      CQMemManager::get().free(AX);
    }

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: JContraction (SplitNBRI) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedRITPIContraction::JContractSplitNBRI called without MPI support");
#endif

  }; // DistributedRITPIContraction::JContractSplitNBRI



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  a one-body operator.
   */   
  template<typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT,IntsT>::KContract
        (MPI_Comm comm, TwoBodyContraction<MatsT>& C) const {

    InCoreRITPI<IntsT> &ritpi = *std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_);
    auto eri3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ritpi.eri3j());
    if (eri3j == nullptr)
      CErr("DistributedRITPIContraction::KContract expect a DistributedERI3J integral object");

    if (eri3j->distributionLayout() == DistributionLayout::SplitNB)
      return KContractSplitNB(comm, C, eri3j);
    else
      return KContractSplitNBRI(comm, C, eri3j);

  }; // DistributedRITPIContraction::KContract



  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::KContractSplitNB(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const {
    auto contractStart = tick();
    int rank = 0, size = 1;
#ifdef CQ_ENABLE_MPI
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size);

    size_t NB = eri3j->nBasis();
    size_t NB2 = NB*NB;
    size_t NBRI = eri3j->nRIBasis();
    size_t NTT = NB*(NB+1)/2;
    size_t localNB = eri3j->localSize();
    size_t localNBStart = eri3j->localStart();
    size_t Ktemp_size = localNB*NB*NBRI;

    MatsT *X  = C.X;
    MatsT *AX = C.AX;

    auto Ktemp_local = CQMemManager::get().malloc<MatsT>(Ktemp_size);
    auto AX_local_NTT = CQMemManager::get().malloc<MatsT>(NTT);
    auto AX_global_NTT = CQMemManager::get().malloc<MatsT>(NTT);
    std::fill_n(Ktemp_local, Ktemp_size, MatsT(0.));
    std::fill_n(AX, NB*NB, MatsT(0.));
    std::fill_n(AX_local_NTT, NTT, MatsT(0.));
    std::fill_n(AX_global_NTT, NTT, MatsT(0.));

    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);

    #pragma omp parallel for
    // Step 1: For each local q, compute K^q_[α][r] = ∑_{s} L^q[α][s] · (P[r][s])^T
    for(auto q = 0ul; q < localNB; q++)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBRI,NB,NB,MatsT(1.),eri3j->data()+q*NB*NBRI,NBRI,X,NB,MatsT(0.),Ktemp_local+q*NB*NBRI,NBRI);
    
    SetLAThreads(LAThreads);

    // Initialize timers
    double sendElapsedTotal = 0.0, recvElapsedTotal = 0.0, waitElapsedTotal = 0.0;
    auto& plan = eri3j->getCommPlan();

    // Step 2: Post receives for all ranks to form off-diagonal blocks
    // 2.1 Set up receive buffers (in round-robin fashion)
    size_t nBuffers = KCONTRACTSPLITNB_NBUFFERS_;
    std::vector<MatsT*>      recvBuffers    (nBuffers, nullptr);
    std::vector<size_t>      recvBufferSizes(nBuffers, 0);
    std::vector<std::vector<MPI_Request>> recvRequests(nBuffers); // Extra vector for chunking data
    std::vector<int>         srcRank        (nBuffers, -1);
    int tag = 777; 

    // 2.2 Define a helper function to (re)post a receive
    auto post_recv = [&](size_t slot, int peer) {
      size_t bufferSize = eri3j->splitSize(NB, peer, size) * NB * NBRI;
      // Initialize buffer or resize if needed 
      if (recvBuffers[slot] == nullptr || bufferSize > recvBufferSizes[slot]) {
        if (recvBuffers[slot]) CQMemManager::get().free(recvBuffers[slot]);
        recvBuffers[slot] = CQMemManager::get().malloc<MatsT>(bufferSize);
        recvBufferSizes[slot] = bufferSize;
      }
      recvRequests[slot] = MPIIrecv(recvBuffers[slot], bufferSize, peer, tag * 1000, comm);
      srcRank[slot] = peer;
    };

    // 2.3 Initialize non-blocking receives for the first set of ranks
    size_t nRecvTotal = plan.recvFromRanks.size();
    size_t posted  = 0;
    for (; posted < std::min(nRecvTotal, nBuffers); ++posted) {
      post_recv(posted, plan.recvFromRanks[posted]);
    }
    size_t nRecvDone = 0;
    size_t cursor = 0;

    // Step 3: Send Ktemp_local to the appropriate rank
    std::vector<MPI_Request> sendRequests;
    for (const auto& dst : plan.sendToRanks) {
      auto sendStart = tick();
      auto reqs = MPIIsend(Ktemp_local, Ktemp_size, dst, tag * 1000, comm);
      sendRequests.insert(sendRequests.end(), reqs.begin(), reqs.end());
#ifdef _REPORT_COMM_TIMINGS
      sendElapsedTotal += tock(sendStart);
#endif
    }

    // Step 4: Compute Diagonal blocks of K_{p_local q_local} = ∑_{αr} (L[α][r][q_local])^T K[α][r][p_local]
    std::map<std::pair<int, int>, MatsT*> ownedKBlocks;
    std::pair<int,int> diagBlock = {rank, rank}; 
    ownedKBlocks[diagBlock] = CQMemManager::get().malloc<MatsT>(localNB * localNB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,localNB,localNB,NB*NBRI,MatsT(1.),eri3j->data(),NB*NBRI,Ktemp_local,NB*NBRI,MatsT(0.),ownedKBlocks[diagBlock],localNB);

    // Step 5: Wait for data to arrive and process received data to form off-diagonal blocks
    while (nRecvDone < nRecvTotal) {

      auto waitStart = tick();
      if (not recvRequests[cursor].empty())
        MPIWait(recvRequests[cursor]);
#ifdef _REPORT_COMM_TIMINGS
      waitElapsedTotal += tock(waitStart);
#endif

      int peer = srcRank[cursor];
      // Step 5.1: Determine how to form the off-diagonal blocks (directly or via a transposition)
      std::pair<int, int> fwdKey = {rank, peer};
      std::pair<int, int> revKey = {peer, rank};
      auto owns = [&](const std::pair<int,int>& key) {
        return std::find(plan.ownedBlocks.begin(),plan.ownedBlocks.end(),key) != plan.ownedBlocks.end();};
      // Form K_{p_local, q_remote} upper-triangular block
      bool ownsFwd = owns(fwdKey);
      // Form K_{p_local, q_remote} lower-triangular block, then transpose it
      bool ownsRev = owns(revKey);
      if (!ownsFwd and !ownsRev) 
        CErr("Unexpected: no matching owned block for srcRank in CommPlan.");
      auto outKey = ownsFwd ? fwdKey : revKey;
      int p_size = eri3j->splitSize(NB, rank, size);
      int q_size = eri3j->splitSize(NB, peer, size);
      // Allocate a temporary buffer for the transposition if needed
      MatsT* dest = nullptr;
      MatsT* finalDest = CQMemManager::get().malloc<MatsT>(p_size * q_size);
      if (ownsRev) {
        dest = CQMemManager::get().malloc<MatsT>(p_size * q_size);
      } else {
        dest = finalDest;
      }
      ownedKBlocks[outKey] = finalDest;
      // Step 5.2: K_{p_local, q_remote} = Lᵀ × Ktemp_remote
      blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,p_size, q_size, NB*NBRI,MatsT(1.0),eri3j->data(),NB*NBRI,
                    recvBuffers[cursor],NB*NBRI,MatsT(0.0),dest,p_size);
      // Step 5.3: Transpose the block if the block is in the lower-triangular form
      if (ownsRev) {
        SetMat('C',p_size, q_size, MatsT(1.0), dest, p_size, finalDest, q_size);
        CQMemManager::get().free(dest);
      }

      // 4.5 Update the cursor and post a receive for the next buffer
      ++nRecvDone;
      recvRequests[cursor].clear();
      if (posted < nRecvTotal) post_recv(cursor, plan.recvFromRanks[posted++]);

      // Move to the next buffer (round-robin)
      cursor = (cursor + 1) % nBuffers;
    }

    // Cleanup receive buffers
    for (auto* p : recvBuffers)
      if (p) CQMemManager::get().free(p);

    //printOwnedKBlocks(ownedKBlocks, *eri3j, NB, size);
    int tagBase = 888;
    std::vector<MPI_Request> sendKBlockRequests, recvKBlockRequests;
    std::map<std::pair<int,int>, std::vector<MatsT>> recvKBlockBuffers;
    auto const& ownerMap = eri3j->blockOwnerMap();

    // Rank 0: Insert rank-0 blocks and post receives for all other blocks
    if (rank == 0) {
      for (auto const& [block, owner] : ownerMap) {
        int i = block.first, j = block.second;
        int p_size = eri3j->splitSize(NB, i, size);
        int q_size = eri3j->splitSize(NB, j, size);
        int count = p_size * q_size;
        int tag = tagBase + i * size + j;
        if (owner == rank) {
          // Insert rank-0 blocks into the full K matrix
          auto it = ownedKBlocks.find(block);
          insertBlockIntoFullK(AX, i, j, it->second, *eri3j, NB, size);
        } else {
          // Post receives for all other blocks not owned by rank 0
          recvKBlockBuffers[block].resize(count);
          auto rqs = MPIIrecv(recvKBlockBuffers[block].data(), count, owner, tag, comm);
          recvKBlockRequests.insert(recvKBlockRequests.end(), rqs.begin(), rqs.end());
        }
      }
    } 
    // Other ranks: send their blocks to rank 0
    else {
      for (auto const& [block, data] : ownedKBlocks) {
        int i = block.first, j = block.second;
        int p_size = eri3j->splitSize(NB, i, size);
        int q_size = eri3j->splitSize(NB, j, size);
        int count = p_size * q_size;
        int tag = tagBase + i * size + j;
        auto sqs = MPIIsend(data, count, 0, tag, comm);
        sendKBlockRequests.insert(sendKBlockRequests.end(), sqs.begin(), sqs.end());
      }
    }

    // Wait for all sends and receives to complete
    TIMED_COMM("KContraction MPIWait", {
      if (!sendKBlockRequests.empty())  MPIWait(sendKBlockRequests);
      if (!recvKBlockRequests.empty())  MPIWait(recvKBlockRequests);
    });

    // On Rank 0, insert received blocks into the full K matrix
    if (rank == 0) {
      for (auto const& [block, data] : recvKBlockBuffers) {
        insertBlockIntoFullK(AX, block.first, block.second, data.data(), *eri3j, NB, size);
      }
    }

    //prettyPrintSmart(std::cout,"AX Global",AX,NB,NB,NB);

    TIMED_COMM("KContraction MPI_Waitall", {
      MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
    });

    CQMemManager::get().free(Ktemp_local, AX_local_NTT, AX_global_NTT);
    for (auto& [key, ptr] : ownedKBlocks)
      CQMemManager::get().free(ptr);

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: KContraction (SplitNB) MPI_Isend Total = " << sendElapsedTotal << " s" << std::endl;
    std::cout << "Rank " << rank << " :: KContraction (SplitNB) MPI_Irecv Total = " << recvElapsedTotal << " s" << std::endl;
    std::cout << "Rank " << rank << " :: KContraction (SplitNB) MPI_Wait  Total = " << waitElapsedTotal << " s" << std::endl;
    std::cout << "Rank " << rank << " :: KContraction (SplitNB) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedRITPIContraction::KContractSplitNB called without MPI support");
#endif

  }; // DistributedRITPIContraction::KContractSplitNB

  

  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::KContractSplitNBRI_real_impl(
      MPI_Comm comm, const double *Xr, double *AX_local, double *Ktemp, const DistributedERI3J<double>& eri3j) const {
    int rank, size;
#ifdef CQ_ENABLE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    size_t NB = eri3j.nBasis();
    size_t localNBRI = eri3j.localSize();
    auto mapCen2BfSt = this->mapCen2BfSt();
    size_t nAtoms = mapCen2BfSt.size();

    std::fill_n(Ktemp, localNBRI*NB*NB, double(0.));
    std::fill_n(AX_local, NB*NB, double(0.));
    auto XrT = CQMemManager::get().malloc<double>(NB*NB);
    SetMat('T',NB,NB,double(1.),Xr,NB,XrT,NB);
    

    // Step 1: For each q, compute T^q_[α_local][r] = ∑_{s} L^q[α_local][s] · (P[r][s])^T
    auto step1Start = tick();
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      #pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < nAtoms; ++i) {
        size_t bfStart = mapCen2BfSt[i];
        size_t bfEnd   = (i+1 < nAtoms) ? mapCen2BfSt[i+1] : NB;  
        size_t nBf = bfEnd - bfStart;
        for(auto q = bfStart; q < bfEnd; q++) 
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,localNBRI,nBf,nBf,double(1.),
                     eri3j.data()+q*NB*localNBRI+bfStart*localNBRI,localNBRI,
                     XrT+bfStart*NB+bfStart,NB,double(0.),Ktemp+q*NB*localNBRI+bfStart*localNBRI,localNBRI);
      }
    } else {
      #pragma omp parallel for
      for(auto q = 0ul; q < NB; q++)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,localNBRI,NB,NB,double(1.),eri3j.data()+q*NB*localNBRI,localNBRI,XrT,NB,double(0.),Ktemp+q*NB*localNBRI,localNBRI);
    }
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KContractSplitNBRI_real_impl Step 1 Contraction = " << tock(step1Start) << " s" << std::endl;
    CQMemManager::get().free(XrT);

    // Step 2: Compute K_local[p][q] = ∑_{α_local, r} K^q_[α_local][r] · T^q_[α_local][r]
    auto step2Start = tick();
    SetLAThreads(LAThreads);
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      for (size_t i = 0; i < nAtoms; ++i) {
        size_t bfStart = mapCen2BfSt[i];
        size_t bfEnd   = (i+1 < nAtoms) ? mapCen2BfSt[i+1] : NB;
        size_t nBf = bfEnd - bfStart;
        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,nBf,nBf,NB*localNBRI,double(1.),
                   eri3j.data()+bfStart*NB*localNBRI,localNBRI*NB,
                   Ktemp+bfStart*NB*localNBRI,localNBRI*NB,
                   double(0.),AX_local+bfStart*NB+bfStart,NB);
      }
    } else {
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB,NB,NB*localNBRI,double(1.),eri3j.data(),localNBRI*NB,Ktemp,localNBRI*NB,double(0.),AX_local,NB);
    }     
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KContractSplitNBRI_real_impl Step 2 = " << tock(step2Start) << " s" << std::endl;


#else
    CErr("DistributedRITPIContraction::KContractSplitNBRI_real_impl called without MPI support");
#endif
  }; // DistributedRITPIContraction::KContractSplitNBRI_real_impl



  template <>
  void DistributedRITPIContraction<double, double>::KContractSplitNBRI(
      MPI_Comm comm, TwoBodyContraction<double> &C, const std::shared_ptr<DistributedERI3J<double>>& eri3j) const {
    auto contractStart = tick();
    int rank, size;
#ifdef CQ_ENABLE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    size_t NB = eri3j->nBasis();
    size_t localNBRI = eri3j->localSize();
    size_t Ktemp_size = localNBRI*NB*NB;

    auto Ktemp = CQMemManager::get().malloc<double>(Ktemp_size);
    auto AX_local = CQMemManager::get().malloc<double>(NB*NB);
    std::fill_n(C.AX, NB*NB, double(0.));

    KContractSplitNBRI_real_impl(comm, C.X, AX_local, Ktemp, *eri3j);

    // Step 3: Sum up K_local[p][q] over all ranks
    TIMED_COMM("KContraction MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB*NB), C.AX, 0, comm)
    );

    CQMemManager::get().free(Ktemp, AX_local);

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: KContraction (SplitNBRI) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedRITPIContraction::KContractSplitNBRI called without MPI support");
#endif

  }; // DistributedRITPIContraction<double, double>::KContractSplitNBRI



  template <>
  void DistributedRITPIContraction<dcomplex, double>::KContractSplitNBRI(
      MPI_Comm comm, TwoBodyContraction<dcomplex> &C, const std::shared_ptr<DistributedERI3J<double>>& eri3j) const {
    auto contractStart = tick();
    int rank, size;
#ifdef CQ_ENABLE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    size_t NB = eri3j->nBasis();
    size_t NB2 = NB*NB;
    size_t localNBRI = eri3j->localSize();
    size_t Ktemp_size = localNBRI*NB*NB;

    // Separate real and imaginary parts of density matrix
    double *Xr = CQMemManager::get().malloc<double>(NB*NB);
    double *Xi = CQMemManager::get().malloc<double>(NB*NB);
    #pragma omp parallel for
    for(auto k = 0ul; k < NB2; k++) {
      Xr[k] = std::real(C.X[k]);
      Xi[k] = std::imag(C.X[k]);
    }

    // Allocate result buffer and scratch
    double *AX_localr = CQMemManager::get().malloc<double>(NB*NB);
    double *AX_locali = CQMemManager::get().malloc<double>(NB*NB);
    double *Ktemp = CQMemManager::get().malloc<double>(Ktemp_size);
    std::fill_n(AX_localr, NB*NB, double(0.));
    std::fill_n(AX_locali, NB*NB, double(0.));

    KContractSplitNBRI_real_impl(comm, Xr, AX_localr, Ktemp, *eri3j);
    KContractSplitNBRI_real_impl(comm, Xi, AX_locali, Ktemp, *eri3j);

    // Combine real and imaginary parts
    dcomplex *AX_local = CQMemManager::get().malloc<dcomplex>(NB*NB);
    #pragma omp parallel for
    for(auto k = 0ul; k < NB2; k++) {
      AX_local[k] = dcomplex(AX_localr[k], AX_locali[k]);
    }
    
    // Step 3: Sum up K_local[p][q] over all ranks
    TIMED_COMM("KContraction MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB*NB), C.AX, 0, comm)
    );

    CQMemManager::get().free(Xr, Xi, Ktemp, AX_localr, AX_locali, AX_local);

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: KContraction (SplitNBRI) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedRITPIContraction::KContractSplitNBRI called without MPI support");
#endif

  }; // DistributedRITPIContraction<dcomplex, double>::KContractSplitNBRI
  


  template <>
  void DistributedRITPIContraction<dcomplex, dcomplex>::KContractSplitNBRI(
      MPI_Comm comm, TwoBodyContraction<dcomplex> &C, const std::shared_ptr<DistributedERI3J<dcomplex>>& eri3j) const {
    CErr("DistributedRITPIContraction<dcomplex, dcomplex>::KContractSplitNBRI not implemented");
  }; // DistributedRITPIContraction<dcomplex, dcomplex>::KContractSplitNBRI



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  orbital coefficients.
   */
  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::KCoefContract(
       MPI_Comm comm, size_t NO, MatsT *C, MatsT *AX) const {

    InCoreRITPI<IntsT> &ritpi = *std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_);
    auto eri3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ritpi.eri3j());

    if (eri3j == nullptr)
      CErr("DistributedRITPIContraction::KCoefContract expects a DistributedERI3J integral object");
    if (eri3j->distributionLayout() == DistributionLayout::SplitNB)
      CErr("KCoefContract not implemented for SplitNB layout");
    else
      return KCoefContractSplitNBRI(comm, NO, C, AX, eri3j);

   }; // DistributedRITPIContraction::KCoefContract



    /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  orbital coefficients.
   */
  template <>
  void DistributedRITPIContraction<double, double>::KCoefContractSplitNBRI(
       MPI_Comm comm, size_t NO, double *C, double *AX, const std::shared_ptr<DistributedERI3J<double>>& eri3j) const {
    auto contractStart = tick();
    int rank, size;
    #ifdef CQ_ENABLE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    size_t NB = eri3j->nBasis();
    size_t localNBRI = eri3j->localSize();
    size_t NONBRI_local = localNBRI*NO;

    auto *T  = CQMemManager::get().malloc<double>(NONBRI_local*NB);
    auto *AX_local = CQMemManager::get().malloc<double>(NB*NB);
    std::fill_n(AX_local, NB*NB, double(0.));

    if (this->oneCenterK())  std::cout << "One-center KCoefContractSplitNBRI" << std::endl;

    // 1. T(i, J_local | ν) = C(λ, i)^T @ L(J_local, λ | ν)^T
    //    One-center: restrict λ to atom(ν), i.e. T^J_{νi} = Σ_{λ∈atom(ν)} L^J_{νλ} C_{λi}
    auto step1Start = tick();
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      // Zero T — one-center only writes partial columns per ν
      std::fill_n(T, NONBRI_local * NB, double(0.));
      #pragma omp parallel for schedule(dynamic)
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        for (size_t q = bfStart; q < bfEnd; q++)
          // C_A^T: (NO × nBf) at C + bfStart, lda=NB
          // L_q_A: (localNBRI × nBf) at data + q*NB*localNBRI + bfStart*localNBRI
          blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                     NO, localNBRI, nBf, double(1.),
                     C + bfStart, NB,
                     eri3j->data() + q * NB * localNBRI + bfStart * localNBRI, localNBRI,
                     double(0.), T + q * NONBRI_local, NO);
      }
    } else {
      #pragma omp parallel for
      for (auto nu = 0ul; nu < NB; nu++)
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                   NO, localNBRI, NB, double(1.),
                   C, NB, eri3j->data() + nu * localNBRI * NB, localNBRI,
                   double(0.), T + nu * NONBRI_local, NO);
    }
    SetLAThreads(LAThreads);
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Step 1 = " << tock(step1Start) << " s" << std::endl;

    // 2: K_local(mu, nu) = T(iJ, mu)^T @ T(iJ, nu) 
    //    One-center: K is block-diagonal, per-atom syrk of size nBf_A
    auto step2Start = tick();
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        // T_A: (NONBRI_local × nBf) at T + bfStart * NONBRI_local
        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI_local, double(1.),
                   T + bfStart * NONBRI_local, NONBRI_local,
                   double(0.), AX_local + bfStart * NB + bfStart, NB);
        // Symmetrize within block
        for (size_t q = 0; q < nBf; q++)
          for (size_t p = q + 1; p < nBf; p++)
            AX_local[(bfStart + q) + (bfStart + p) * NB] =
                AX_local[(bfStart + p) + (bfStart + q) * NB];
      }
    } else {
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI_local, double(1.), T, NONBRI_local,
                 double(0.), AX_local, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = q + 1; p < NB; p++)
          AX_local[q + p * NB] = AX_local[p + q * NB];
    }
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Step 2 = " << tock(step2Start) << " s" << std::endl;

    // 3. Sum up K_local over all ranks
    TIMED_COMM("KCoefContractSplitNBRI MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB*NB), AX, 0, comm)
    );
    
    CQMemManager::get().free(T, AX_local);
    
#ifdef _REPORT_COMM_TIMINGS
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Total = " << tock(contractStart) << " s" << std::endl;
#endif
    
#else
    CErr("DistributedRITPIContraction::KCoefContractSplitNBRI called without MPI support");
#endif
 
  }; // DistributedRITPIContraction<double, double>::KCoefContractSplitNBRI
 
 
 
  template <>
  void DistributedRITPIContraction<dcomplex, double>::KCoefContractSplitNBRI(
      MPI_Comm comm, size_t NO, dcomplex *C, dcomplex *AX, const std::shared_ptr<DistributedERI3J<double>>& eri3j) const {
    auto contractStart = tick();
    int rank, size;
#ifdef CQ_ENABLE_MPI
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    size_t NB = eri3j->nBasis();
    size_t localNBRI = eri3j->localSize();
    size_t NONBRI_local = localNBRI*NO;

    auto *Cr  = CQMemManager::get().malloc<double>(NO*NB);
    auto *Ci  = CQMemManager::get().malloc<double>(NO*NB);
    #pragma omp parallel for
    for(auto k = 0ul; k < NO*NB; k++) {
      Cr[k] = std::real(C[k]);
      Ci[k] = std::imag(C[k]);
    }
     
    // Half-transformed buffers (real)
    double *T_re = CQMemManager::get().malloc<double>(NONBRI_local * NB);
    double *T_im = CQMemManager::get().malloc<double>(NONBRI_local * NB);

    // 1. T_re(i, J_local | ν) = Re(C)(λ, i)^T @ L(J_local, λ | ν)^T
    //    T_im(i, J_local | ν) = Im(C)(λ, i)^T @ L(J_local, λ | ν)^T
    //    One-center: restrict λ to atom(ν)
    auto step1Start = tick();
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      std::fill_n(T_re, NONBRI_local * NB, double(0.));
      std::fill_n(T_im, NONBRI_local * NB, double(0.));
      #pragma omp parallel for schedule(dynamic)
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        for (size_t q = bfStart; q < bfEnd; q++) {
          blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                     NO, localNBRI, nBf, double(1.),
                     Cr + bfStart, NB,
                     eri3j->data() + q * NB * localNBRI + bfStart * localNBRI, localNBRI,
                     double(0.), T_re + q * NONBRI_local, NO);
          blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                     NO, localNBRI, nBf, double(1.),
                     Ci + bfStart, NB,
                     eri3j->data() + q * NB * localNBRI + bfStart * localNBRI, localNBRI,
                     double(0.), T_im + q * NONBRI_local, NO);
        }
      }
    } else {
      #pragma omp parallel for
      for (size_t q = 0; q < NB; q++) {
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                   NO, localNBRI, NB, double(1.),
                   Cr, NB, eri3j->data() + q * NB * localNBRI, localNBRI,
                   double(0.), T_re + q * NONBRI_local, NO);
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                   NO, localNBRI, NB, double(1.),
                   Ci, NB, eri3j->data() + q * NB * localNBRI, localNBRI,
                   double(0.), T_im + q * NONBRI_local, NO);
      }
    }
    SetLAThreads(LAThreads);
    CQMemManager::get().free(Cr, Ci);
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Step 1 = " << tock(step1Start) << " s" << std::endl;

    double *K_temp = CQMemManager::get().malloc<double>(NB*NB);
    dcomplex *AX_local = CQMemManager::get().malloc<dcomplex>(NB*NB);
    std::fill_n(K_temp, NB*NB, double(0.));
    std::fill_n(AX_local, NB*NB, dcomplex(0.));
  
    // 2a. Re(K) = T_re^T · T_re + T_im^T · T_im
    auto step2Start = tick();
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
 
      // 2a. Per-atom: Re(K_A) = T_re_A^T · T_re_A + T_im_A^T · T_im_A
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        size_t off = bfStart * NONBRI_local;
        size_t koff = bfStart * NB + bfStart;
 
        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI_local, double(1.),
                   T_re + off, NONBRI_local,
                   double(0.), K_temp + koff, NB);
        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI_local, double(1.),
                   T_im + off, NONBRI_local,
                   double(1.), K_temp + koff, NB);
 
        // Symmetrize Re(K_A) into AX_local
        for (size_t q = 0; q < nBf; q++)
          for (size_t p = 0; p <= q; p++) {
            double val = K_temp[(bfStart + q) + (bfStart + p) * NB];
            AX_local[(bfStart + p) + (bfStart + q) * NB] += val;
            if (p != q)
              AX_local[(bfStart + q) + (bfStart + p) * NB] += val;
          }
      }
 
      // 2b. Per-atom: Im(K_A) = T_re_A^T · T_im_A - (T_re_A^T · T_im_A)^T
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        size_t off = bfStart * NONBRI_local;
        size_t koff = bfStart * NB + bfStart;
 
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                   nBf, nBf, NONBRI_local, double(1.),
                   T_re + off, NONBRI_local,
                   T_im + off, NONBRI_local,
                   double(0.), K_temp + koff, NB);
 
        for (size_t q = 0; q < nBf; q++)
          for (size_t p = 0; p <= q; p++) {
            double imval = K_temp[(bfStart + p) + (bfStart + q) * NB]
                         - K_temp[(bfStart + q) + (bfStart + p) * NB];
            AX_local[(bfStart + p) + (bfStart + q) * NB] += dcomplex(0., imval);
            if (p != q)
              AX_local[(bfStart + q) + (bfStart + p) * NB] += dcomplex(0., -imval);
          }
      }
 
    } else {
 
      // 2a. Re(K) = T_re^T · T_re + T_im^T · T_im
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI_local, double(1.),
                 T_re, NONBRI_local,
                 double(0.), K_temp, NB);
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI_local, double(1.),
                 T_im, NONBRI_local,
                 double(1.), K_temp, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = 0; p <= q; p++) {
          double val = K_temp[q + p * NB];
          AX_local[p + q * NB] += val;
          if (p != q)
            AX_local[q + p * NB] += val;
        }
 
      // 2b. Im(K) = T_re^T · T_im - (T_re^T · T_im)^T
      blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                 NB, NB, NONBRI_local, double(1.),
                 T_re, NONBRI_local, T_im, NONBRI_local,
                 double(0.), K_temp, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = 0; p <= q; p++) {
          double imval = K_temp[p + q * NB] - K_temp[q + p * NB];
          AX_local[p + q * NB] += dcomplex(0., imval);
          if (p != q)
            AX_local[q + p * NB] += dcomplex(0., -imval);
        }
    }
 
    CQMemManager::get().free(T_re, T_im, K_temp);
    if (this->printContractionTiming)
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Step 2 = " << tock(step2Start) << " s" << std::endl;

    // 3. Sum up K_local over all ranks
    TIMED_COMM("KCoefContractSplitNBRI MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB*NB), AX, 0, comm)
    );
  
    CQMemManager::get().free(AX_local);
  
  #ifdef _REPORT_COMM_TIMINGS
      std::cout << "Rank " << rank << " :: KCoefContractSplitNBRI Total = " << tock(contractStart) << " s" << std::endl;
  #endif
  
  #else
      CErr("DistributedRITPIContraction::KCoefContractSplitNBRI called without MPI support");
  #endif
  
  }; // DistributedRITPIContraction<dcomplex, double>::KCoefContractSplitNBRI


  template <>
  void DistributedRITPIContraction<dcomplex, dcomplex>::KCoefContract(
      MPI_Comm comm, size_t NO, dcomplex *C, dcomplex *AX) const {
    CErr("DistributedRITPIContraction<dcomplex, dcomplex>::KCoefContract not implemented");
  }; // DistributedRITPIContraction<dcomplex, dcomplex>::KCoefContract


  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::printOwnedKBlocks(
      const std::map<std::pair<int, int>, MatsT*>& blocks,
      const DistributedERI3J<IntsT>& eri3j, int NB, int size) const {
    for (const auto& [blk, data] : blocks) {
      int i = blk.first;
      int j = blk.second;
      size_t mu_start = eri3j.splitStart(NB, i, size);
      size_t nu_start = eri3j.splitStart(NB, j, size);
      size_t mu_size  = eri3j.splitSize(NB, i, size);
      size_t nu_size  = eri3j.splitSize(NB, j, size);

      std::cout << "Block (" << i << "," << j << ") → [" << mu_size << " x " << nu_size
                << "] at (" << mu_start << "," << nu_start << "), LD=" << mu_size << ":\n";

      for (size_t mu = 0; mu < mu_size; ++mu) {
        for (size_t nu = 0; nu < nu_size; ++nu) {
          std::cout << std::setw(10) << data[nu * mu_size + mu] << " ";
        }
        std::cout << "\n";
      }
      std::cout << std::endl;
    }
  }; // DistributedRITPIContraction::printOwnedKBlocks

  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::placeKBlocksNTT(
      const std::map<std::pair<int, int>, MatsT*>& blocks,
      MatsT* AX_local_NTT, const DistributedERI3J<IntsT>& eri3j, int NB, int size) const {

    MatsT* AX_local_full = CQMemManager::get().malloc<MatsT>(NB*NB);
    std::fill_n(AX_local_full, NB*NB, MatsT(0.));
    for (const auto& [blk, data] : blocks) {
      int i = blk.first;
      int j = blk.second;

      size_t p_start = eri3j.splitStart(NB, i, size);
      size_t q_start = eri3j.splitStart(NB, j, size);
      size_t p_size  = eri3j.splitSize(NB, i, size);
      size_t q_size  = eri3j.splitSize(NB, j, size);

      for (size_t p = 0; p < p_size; ++p) {
        for (size_t q = 0; q < q_size; ++q) {
          size_t global_row = p_start + p;
          size_t global_col = q_start + q;
          AX_local_full[global_col * NB + global_row] = data[q * p_size + p];
        } 
      }
    }

    for (int j=0, idx=0; j<NB; ++j)
      for (int i=0; i<=j; ++i)
          AX_local_NTT[idx++] = AX_local_full[j*NB + i];
    CQMemManager::get().free(AX_local_full);
  }; // DistributedRITPIContraction::placeKBlocksNTT



  template <typename MatsT, typename IntsT>
  void DistributedRITPIContraction<MatsT, IntsT>::insertBlockIntoFullK(MatsT* AX, int i, int j, const MatsT* blockData, 
      const DistributedERI3J<IntsT>& eri3j, int NB, int size) const {
    
    //std::cout << "Inserting block (" << i << "," << j << ") into full K matrix" << std::endl;
    int p_size = eri3j.splitSize(NB, i, size);
    int q_size = eri3j.splitSize(NB, j, size);
    int p_start = eri3j.splitStart(NB, i, size);
    int q_start = eri3j.splitStart(NB, j, size);

    for (int p = 0; p < p_size; ++p) {
      for (int q = 0; q < q_size; ++q) {
        int row = p_start + p;
        int col = q_start + q;
        AX[col*NB + row] = blockData[ q * p_size + p ];
        if (i != j) {
          // also fill symmetric entry
          AX[row*NB + col] = SmartConj(blockData[ q * p_size + p ]);
        }
      }
    }
  }; // DistributedRITPIContraction::insertBlockIntoFullK



  /**
   *  \brief Perform a Coulomb-type (34,12) RI-ERI contraction with
   *  a one-body operator.
   *  J_{pq} = ∑_α ∑_{RS} L^α_{pq} · L^α_{RS} · P_{RS}  
   */   
  template<typename MatsT, typename IntsT>
  void DistributedAsymmRITPIContraction<MatsT,IntsT>::JContract
        (MPI_Comm comm, TwoBodyContraction<MatsT>& C) const {

    auto asymmRITPI = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(this->ints_);

    auto asymmCDalg = asymmRITPI->asymmCDalg();
    std::shared_ptr<DistributedERI3J<IntsT>> eri3j1, eri3j2;
    if (asymmCDalg == ASYMM_CD_ALG::INT1_AUX) {
      eri3j1 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->getAux1()->eri3j());
      eri3j2 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->partialTPI());
    } else if (asymmCDalg == ASYMM_CD_ALG::INT2_AUX) {
      eri3j1 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->partialTPI());
      eri3j2 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->getAux2()->eri3j());
    } else if (asymmCDalg == ASYMM_CD_ALG::COMBINEAUXBASIS) {
      eri3j1 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->getAux1()->eri3j());
      eri3j2 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(asymmRITPI->getAux2()->eri3j());
    } else {
      CErr("DistributedAsymmRITPIContraction::JContract called with unsupported ASYMM_CD_ALG");
    }

    // Sanity check:
    if (eri3j1->nRIBasis() != eri3j2->nRIBasis())
      CErr("DistributedAsymmRITPIContraction::JContract: RI basis functions of aux1 and aux2 do not match");  
    if (eri3j1->distributionLayout() != eri3j2->distributionLayout())
      CErr("DistributedAsymmRITPIContraction::JContract: distribution layout of aux1 and aux2 do not match");  

    if(eri3j1->distributionLayout() == DistributionLayout::SplitNB)
      JContractSplitNB(comm, C, eri3j1, eri3j2);
    else
      JContractSplitNBRI(comm, C, eri3j1, eri3j2);
  }; // DistributedAsymmRITPIContraction::JContract



  template<typename MatsT, typename IntsT>
  void DistributedAsymmRITPIContraction<MatsT,IntsT>::JContractSplitNB
        (MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j1, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j2) const {
    auto contractStart = tick();
    int rank = 0, size = 1;
#ifdef CQ_ENABLE_MPI
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size); 
    auto L3J = eri3j1;
    auto R3J = eri3j2;
    // If contractSecond set to true, need to modify order to create (pp|ee) ints
    if (this->contractSecond) std::swap(L3J, R3J);
    size_t NB = L3J->nBasis();
    size_t snNB = R3J->nBasis();
    size_t localNB = L3J->localSize();
    size_t localsnNB = R3J->localSize();
    size_t NBRI = L3J->nRIBasis();

    // Prepare input and output matrices
    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    IntsT* AX = reinterpret_cast<IntsT*>(C.AX);
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;
    if( extractRealPartX ) {
      X = CQMemManager::get().malloc<IntsT>(snNB*snNB);
      for(auto k = 0ul; k < snNB*snNB; k++) X[k] = std::real(C.X[k]);
    }
    if( allocAXScratch ) AX = CQMemManager::get().malloc<IntsT>(NB*NB);
    std::fill_n(AX,NB*NB,IntsT(0.));
      
    // Load local density matrix P[R_local][S]
    auto X_local = CQMemManager::get().malloc<IntsT>(localsnNB*snNB);
    auto AX_local = CQMemManager::get().malloc<IntsT>(localNB*NB);
    size_t idx = 0; 
      for (size_t i = 0; i < localsnNB; i++) 
      { size_t R_global = R3J->localStart() + i;
        for (size_t S = 0; S < snNB; S++) X_local[idx++] = X[R_global + S*snNB]; }

    // Allocate temporary storage for T[α]
    auto Ttemp = CQMemManager::get().malloc<IntsT>(NBRI);
    auto Ttemp_global = CQMemManager::get().malloc<IntsT>(NBRI);
    std::fill_n(Ttemp, NBRI, IntsT(0.));
    std::fill_n(Ttemp_global, NBRI, IntsT(0.));
      
    // Step 1: T[α] = ∑_{R_local,S} L[α][S][R_local] * P[R_local][S]
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,1,snNB*localsnNB,IntsT(1.),R3J->data(),NBRI,X_local,snNB*localsnNB,IntsT(0.),Ttemp,NBRI);
    // Step 2: Sum T[α] over all ranks 
    TIMED_COMM("JContraction MPIAllReduce",
      MPIAllReduce(Ttemp, static_cast<int>(NBRI), Ttemp_global, comm)
    );
      
    // Step 3: J[p_local][q] = ∑_{α} L[α][p_local][q] * T[α] 
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,localNB*NB,1,NBRI,IntsT(1.),L3J->data(),NBRI,Ttemp_global,NBRI,IntsT(0.),AX_local,localNB*NB);
      
    // Step 4: Assemble global J matrix
    std::vector<int> recvcounts(size), displs(size);
    for (int r = 0; r < size; ++r) {
      size_t r_localNB = L3J->splitSize(NB, r, size);
      displs[r] = L3J->splitStart(NB, r, size) * NB;
      recvcounts[r] = static_cast<int>(r_localNB * NB);
    } 
    IntsT *bufPtr  = (rank == 0) ? AX : nullptr;
    int *recvPtr   = (rank == 0) ? recvcounts.data() : nullptr;
    int *displsPtr = (rank == 0) ? displs.data() : nullptr;
    TIMED_COMM("JContraction MPI_Gatherv",
      MPI_Gatherv(AX_local, localNB * NB, mpi_data_type<IntsT>(),
                bufPtr, recvPtr, displsPtr, mpi_data_type<IntsT>(),
                0, comm)
    );
    CQMemManager::get().free(Ttemp_global, Ttemp, X_local, AX_local);
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {
      std::copy_n(AX,NB*NB,C.AX);
      CQMemManager::get().free(AX);
    }

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: Asymm JContraction (SplitNB) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedAsymmRITPIContraction::JContractSplitNB called without MPI support");
#endif
  }; // DistributedAsymmRITPIContraction::JContractSplitNBRI



  template<typename MatsT, typename IntsT>
  void DistributedAsymmRITPIContraction<MatsT,IntsT>::JContractSplitNBRI
        (MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j1, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j2) const {
    auto contractStart = tick();
    int rank = 0, size = 1;
#ifdef CQ_ENABLE_MPI
MPI_Comm_rank(comm, &rank);
MPI_Comm_size(comm, &size); 
    auto L3J = eri3j1;
    auto R3J = eri3j2;
    // If contractSecond set to true, need to modify order to create (pp|ee) ints
    if (this->contractSecond) std::swap(L3J, R3J);
    size_t NB = L3J->nBasis();
    size_t snNB = R3J->nBasis();
    if (L3J->localSize() != R3J->localSize() or L3J->localStart() != R3J->localStart())
      CErr("DistributedAsymmRITPIContraction::JContractSplitNBRI: RI basis functions distributed differently in aux1 and aux2");
    size_t localNBRI = L3J->localSize();
    
    // Prepare input and output matrices
    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    IntsT* AX = reinterpret_cast<IntsT*>(C.AX);
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;
    if( extractRealPartX ) {
      X = CQMemManager::get().malloc<IntsT>(snNB*snNB);
      for(auto k = 0ul; k < snNB*snNB; k++) X[k] = std::real(C.X[k]);
    }
    if( allocAXScratch ) AX = CQMemManager::get().malloc<IntsT>(NB*NB);
    std::fill_n(AX,NB*NB,IntsT(0.));

    // Allocate temporary storage for T[α]
    auto Ttemp = CQMemManager::get().malloc<IntsT>(localNBRI);
    auto AX_local = CQMemManager::get().malloc<IntsT>(NB*NB);
    std::fill_n(Ttemp, localNBRI, IntsT(0.));
    std::fill_n(AX_local, NB*NB, IntsT(0.));
      
    // Step 1: T_[α_local] = ∑_{R,S} L[α_local][S][R] * P[R][S]
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,localNBRI,1,snNB*snNB,IntsT(1.),R3J->data(),localNBRI,X,snNB*snNB,IntsT(0.),Ttemp,localNBRI);
      
    // Step 2: J_local[p][q] = ∑_{α_local} L[α_local][p][q] * T[α_local] 
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB*NB,1,localNBRI,IntsT(1.),L3J->data(),localNBRI,Ttemp,localNBRI,IntsT(0.),AX_local,NB*NB);
      
    // Step 3: Sum up J_local[p][q] over all ranks
    TIMED_COMM("JContraction MPIReduce",
      MPIReduce(AX_local, static_cast<int>(NB*NB), AX, 0, comm)
    );

    // Clean up temporary storage
    CQMemManager::get().free(Ttemp, AX_local);
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {
      std::copy_n(AX,NB*NB,C.AX);
      CQMemManager::get().free(AX);
    }

#ifdef _REPORT_COMM_TIMINGS
    std::cout << "Rank " << rank << " :: Asymm JContraction (SplitNBRI) Total = " << tock(contractStart) << " s" << std::endl;
#endif

#else
    CErr("DistributedAsymmRITPIContraction::JContractSplitNB called without MPI support");
#endif
  }; // DistributedAsymmRITPIContraction::JContractSplitNBRI


      

}; // namespace ChronusQ