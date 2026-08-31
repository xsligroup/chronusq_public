/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#ifdef CQ_ENABLE_MPI
#include <mpi.h>
#endif
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <cxxapi/output.hpp>
#include <matrix/ndarray.hpp>
#include <vector>
#include <utility>

// Distribution layout for ERI3J
enum class DistributionLayout {
  SplitNB,      ///< split over basis functions
  SplitNBRI     ///< split over auxiliary basis functions
};
namespace ChronusQ {

template <typename IntsT> class InCoreRITPI; // Forward declaration

struct CommPlan {
  std::vector<std::pair<int, int>> ownedBlocks;
  std::vector<int> recvFromRanks; 
  std::vector<int> sendToRanks;    
};

template <typename IntsT>
class ERI3JBase {
  public:
    virtual ~ERI3JBase() = default;

	virtual IntsT* data() = 0;
	virtual const IntsT* data() const = 0;

	virtual size_t nBasis() const = 0;
	virtual size_t nRIBasis() const = 0;
	virtual void clear() = 0;

};


template <typename IntsT>
class IncoreERI3J : public ERI3JBase<IntsT> {
  private:
	std::unique_ptr<cqmatrix::NDArray<IntsT>> data_;         // Flat array (μν, α)
	size_t NB_, NBRI_;

  public:
	IncoreERI3J(size_t NB, size_t NBRI)
		: data_(std::make_unique<cqmatrix::NDArray<IntsT>>(std::vector<size_t>{NB*NB, NBRI})), 
		  NB_(NB), NBRI_(NBRI) {}

	virtual ~IncoreERI3J() {}

	IntsT* data() override { return data_->pointer(); }
	const IntsT* data() const override { return data_->pointer(); }

	size_t nBasis() const override { return NB_; }
	size_t nRIBasis() const override { return NBRI_; }

	virtual void clear() override { data_->clear(); }

};



template <typename IntsT>
class DistributedERI3J : public ERI3JBase<IntsT> {
  private:
	std::unique_ptr<cqmatrix::NDArray<IntsT>> localData_; 
	size_t NB_, NBRI_, localStart_, localSize_;
	DistributionLayout distLayout_ = DistributionLayout::SplitNB; 
	MPI_Comm comm_;
	std::map<std::pair<int, int>, int> blockOwner_;
	CommPlan commPlan_;

  public:
	DistributedERI3J(MPI_Comm comm, size_t nb, size_t nbri)
	    : comm_(comm),  NB_(nb), NBRI_(nbri) {
#ifdef CQ_ENABLE_MPI
      int rank, size;
      MPI_Comm_rank(comm_, &rank);
      MPI_Comm_size(comm_, &size);

	  // === Initially split by basis functions ===
      localStart_ = splitStart(NB_, rank, size);
      localSize_ = splitSize(NB_, rank, size);
      localData_ = std::make_unique<cqmatrix::NDArray<IntsT>>(std::vector<size_t>{NBRI_, NB_, localSize_});
	  
	  // === For KContractSplitNB, build block to rank map to assign K blocks to each MPI process ===
	  blockOwner_ = buildBlockToRankMap(size);
	  //printBlockToRankMap(blockOwner_, size);
	  // === For KContractSplitNB, generate communication plan to send/recv blocks ===
	  commPlan_ = generateCommPlan(rank, size);
	  //printCommPlan(commPlan_, rank);

#else
      CErr("DistributedERI3J called without MPI support");
#endif
	}

	size_t nBasis() const override { return NB_; }
	size_t nRIBasis() const override { return NBRI_; }
	virtual void clear() override { localData_->clear(); }

	DistributionLayout distributionLayout() const { return distLayout_; }
  MPI_Comm comm() const { return comm_; }
	const IntsT* data() const override { return localData_->pointer(); }
	IntsT* data() override { return localData_->pointer(); }
	const CommPlan& getCommPlan() const { return commPlan_; }
	const std::map<std::pair<int, int>, int>& blockOwnerMap() const { return blockOwner_; }
	size_t localStart() const { return localStart_; }
	size_t localSize() const { return localSize_; }

	inline constexpr size_t splitSize(size_t N, int rank, int size) const {
	  size_t chunk = N / size;
	  size_t remainder = N % size;
	  return chunk + (rank < remainder ? 1 : 0);
	}

	inline constexpr size_t splitStart(size_t N, int rank, int size) const {
	  size_t chunk = N / size;
	  size_t remainder = N % size;
	  return rank * chunk + std::min<size_t>(rank, remainder);
	}

	// Helper function to build a map of block ownership
	std::map<std::pair<int, int>, int> buildBlockToRankMap(int nProcs, bool diagOnly=false) const {
	  std::map<std::pair<int, int>, int> ownership;

	  // If diagonal only, set all diagonal blocks to be owned by the rank
	  if (diagOnly) {
	  	for (int i = 0; i < nProcs; ++i) ownership[{i, i}] = i;
		return ownership;
	  } 

	  // Distribute all upper triangular blocks
	  for (int offset = 0; offset < nProcs; ++offset) {
	  	for (int i = 0; i < nProcs - offset; ++i) {
	  	  int j = i + offset;
	  	  if (i == j) 
		  	ownership[{i, j}] = i;  // Diagonal blocks are owned by the rank
	  	  else if (offset < (nProcs + 1) / 2) 
		  	ownership[{i, j}] = i;  // Lower triangular blocks are owned by the lower rank
	  	  else 
		  	ownership[{i, j}] = j;  // Upper triangular blocks are owned by the upper rank
	  	}
	  }
	  return ownership;
	}

	CommPlan generateCommPlan(int myRank, int nProcs, bool diagOnly=false) const {
	  CommPlan plan;

	  for (const auto& [blk, owner] : blockOwner_) {
	  	int i = blk.first;
	  	int j = blk.second;

		if (diagOnly && i != j) continue;

	  	if (owner == myRank) {
	  	  // 1. Record ownership
	  	  plan.ownedBlocks.emplace_back(i, j);
	  	  // 2. Record which blocks to receive from other ranks
	  	  if (i != myRank)
	  	  	plan.recvFromRanks.push_back(i);
	  	  if (j != myRank && j != i)
	  	  	plan.recvFromRanks.push_back(j);
	  	}
	  	// 3. Record which blocks to send to other ranks
		if (i == myRank && owner != myRank)
		  plan.sendToRanks.push_back(owner);
		if (j == myRank && owner != myRank && j != i)
    	  plan.sendToRanks.push_back(owner);
	  }

	  // 4: Remove duplicates
	  std::sort(plan.recvFromRanks.begin(), plan.recvFromRanks.end());
	  plan.recvFromRanks.erase(std::unique(plan.recvFromRanks.begin(), plan.recvFromRanks.end()), plan.recvFromRanks.end());
	  std::sort(plan.sendToRanks.begin(), plan.sendToRanks.end());
	  plan.sendToRanks.erase(std::unique(plan.sendToRanks.begin(), plan.sendToRanks.end()), plan.sendToRanks.end());

	  return plan;
	}

	void printBlockToRankMap(const std::map<std::pair<int, int>, int>& blockOwner, int size) const {
	  std::cout << size << " processors:\n";
	  for (int i = 0; i < size; ++i, std::cout << "\n") {
		for (int j = 0; j < size; ++j) {
		  if (i <= j) {
			auto it = blockOwner.find({i, j});
			if (it != blockOwner.end()) {
			  std::cout << ' ' << std::setw(2) << it->second;
			} else {
			  std::cout << " X ";
			}
		  } else {
			std::cout << "   ";
		  }
	    }
	  }
	}

	void printCommPlan(const CommPlan& plan, int myRank) const{
	  std::cout << "Rank " << myRank << " CommPlan:\n";

	  std::cout << "  Owned Blocks:\n";
	  for (auto [i, j] : plan.ownedBlocks)
	  	std::cout << "    (" << i << "," << j << ")\n";

	  std::cout << "  Send to:\n";
	  for (auto& dst : plan.sendToRanks) {
	  	std::cout << "    → Rank " << dst << "\n";
	  }

	  std::cout << "  Recv from:\n";
	  for (auto& src : plan.recvFromRanks) {
	  	std::cout << "    ← Rank " << src << "\n";
	  }
	}

	void redistributeToSplitNBRI() {
#ifdef CQ_ENABLE_MPI
	  if (distLayout_ == DistributionLayout::SplitNBRI) {
		std::cout << "    * ERI3J already split over auxiliary basis functions" << std::endl;
		return;
	  }

	  int rank, size;
	  MPI_Comm_rank(comm_, &rank);
	  MPI_Comm_size(comm_, &size);

	  // 1. Compute current (NB-split) and target (NBRI-split) slice info
	  size_t myNBStart = splitStart(NB_, rank, size);
	  size_t myNBSize  = splitSize(NB_, rank, size);
	  size_t newLocalStart = splitStart(NBRI_, rank, size);
	  size_t newLocalSize  = splitSize(NBRI_, rank, size);

	  // 2. Prepare the counts and displacements for an MPI_Alltoall
	  std::vector<size_t> sendCounts(size,0); // elements to send to rank p
	  std::vector<size_t> sendDispls(size,0); // offset (in elements) in the packed send buffer
	  std::vector<size_t> recvCounts(size,0); // elements to receive from rank p
	  std::vector<size_t> recvDispls(size,0); // offset (in elements) in the packed receive buffer

	  // 2. Compute the send and receive counts
	  size_t totalSendElements = 0;
	  size_t totalRecvElements = 0;
	  for (int p = 0; p < size; ++p) {
		// how many aux functions rank p owns
		size_t pNBRISize  = splitSize(NBRI_, p, size);
		size_t pNBSize    = splitSize(NB_, p, size);
		sendCounts[p] = pNBRISize * NB_ * myNBSize;
		recvCounts[p] = newLocalSize * NB_ * pNBSize;
		// Set next block starts right after the previous block
		sendDispls[p] = (p == 0) ? 0 : sendDispls[p-1] + sendCounts[p-1];
		recvDispls[p] = (p == 0) ? 0 : recvDispls[p-1] + recvCounts[p-1];
		totalSendElements += sendCounts[p];
		totalRecvElements += recvCounts[p];
	  }

	  // 3. Pack the send buffer
	  IntsT* sendBuffer = CQMemManager::get().malloc<IntsT>(totalSendElements);
      std::vector<size_t> cur(size); // keep track of current position in sendBuffer
      for (int p = 0; p < size; ++p) cur[p] = sendDispls[p];
    
	  for (size_t mu_local = 0; mu_local < myNBSize; ++mu_local) {
	  	for (size_t nu = 0; nu < NB_; ++nu) {
	  	  for (int p = 0; p < size; ++p) {
	  	  	size_t pNBRIStart = splitStart(NBRI_, p, size);
	  	  	size_t pNBRISize  = splitSize(NBRI_, p, size);
	  	  	const IntsT* src = &(*localData_)(pNBRIStart, nu, mu_local);
	  	  	std::memcpy(&sendBuffer[cur[p]], src, pNBRISize * sizeof(IntsT));
	  	  	cur[p] += pNBRISize;
	  	  }
	  	}
	  }
	  // check if the send buffer is packed correctly
	  for (int p = 0; p < size; ++p)
		assert(cur[p] == sendDispls[p] + sendCounts[p]);
	  
	  // 4.Free the old local data (split NB layout)
	  localData_ = nullptr;

	  // 5. Perform the alltoall
	  IntsT* recvBuffer = CQMemManager::get().malloc<IntsT>(totalRecvElements);
	  MPIAlltoallv(sendBuffer, sendCounts, sendDispls,
				   recvBuffer, recvCounts, recvDispls, 
				   comm_);
	  CQMemManager::get().free(sendBuffer);
	  
	  // 6. Allocate new local data (split NBRI layout)
	  localData_ = std::make_unique<cqmatrix::NDArray<IntsT>>(std::vector<size_t>{newLocalSize, NB_, NB_});
	  localData_->clear();

	  // 7. Unpack the receive buffer
	  size_t pos = 0;
	  for (int p = 0; p < size; ++p) {
	    size_t pNBStart = splitStart(NB_, p, size);
	    size_t pNBSize  = splitSize(NB_, p, size);
	    for (size_t pMu = 0; pMu < pNBSize; ++pMu) {
	      size_t globalMu = pNBStart + pMu;
	      for (size_t nu = 0; nu < NB_; ++nu) {
	        IntsT* dst = &(*localData_)(0, nu, globalMu);
	        std::memcpy(dst, &recvBuffer[pos], newLocalSize * sizeof(IntsT));
	        pos += newLocalSize;
	      }
	    }
	  }
	  assert(pos == totalRecvElements);
	  CQMemManager::get().free(recvBuffer);

	  // 8. Update stored variables
	  localStart_ = newLocalStart;
	  localSize_ = newLocalSize;
	  distLayout_ = DistributionLayout::SplitNBRI;
#else
	  CErr("DistributedERI3J::redistributeToAux requires MPI");
#endif
	}

};

} // namespace ChronusQ
