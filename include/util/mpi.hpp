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

#include <chronusq_sys.hpp>
#ifdef CQ_ENABLE_MPI 
#include <mpi.h>
#endif

namespace ChronusQ {

#ifndef CQ_ENABLE_MPI

  struct MPI_Comm { 

    int internal = 0;


    MPI_Comm(int c) : internal(c){ }
    MPI_Comm() : MPI_Comm(0) { } ;

    static inline int size() { return 1; }
    static inline int rank() { return 0; }

    operator int() const { return internal; }

  };
  
  struct MPI_Request {};
  struct MPI_Status  {};

  static inline bool operator==(MPI_Comm c, MPI_Comm d){ 
    return c.internal == d.internal; }
  static inline bool operator!=(MPI_Comm c, MPI_Comm d){ 
    return c.internal != d.internal; }

  static inline bool operator==(MPI_Comm c, int x){ 
    return c.internal == x; }
  static inline bool operator!=(MPI_Comm c, int x){ 
    return not operator==(c,x); }

  static inline bool operator==(int x, MPI_Comm c){ 
    return operator==(c,x); }
  static inline bool operator!=(int x, MPI_Comm c){ 
    return operator!=(c,x); }

  static MPI_Comm MPI_COMM_WORLD{0 };
  static MPI_Comm MPI_COMM_NULL {-1};

#define MPI_UNDEFINED 1

  static inline void MPI_Barrier(MPI_Comm c) { };
  static inline int MPI_Wait(MPI_Request *request, MPI_Status *status) { return 0; }

#else // MPI is enabled
  template <typename T> MPI_Datatype mpi_data_type();
  #define REGISTER_MPI_TYPE(CXXTYPE, MPITYPE) \
  template <> inline MPI_Datatype mpi_data_type<CXXTYPE>() { return MPITYPE; }

  REGISTER_MPI_TYPE(double, MPI_DOUBLE)
  REGISTER_MPI_TYPE(int,    MPI_INT   )
  REGISTER_MPI_TYPE(char,    MPI_CHAR   )
  REGISTER_MPI_TYPE(int64_t,MPI_INT64_T)
  REGISTER_MPI_TYPE(size_t, MPI_UINT64_T)

  // For compilers that are known to have `long int` equivalent to `int64_t`, avoid redefinition
  #if !defined(LONG_INT_IS_INT64_T)
  REGISTER_MPI_TYPE(long int, MPI_LONG)
  #endif

  REGISTER_MPI_TYPE(std::complex<double>, MPI_C_DOUBLE_COMPLEX)

  #undef REGISTER_MPI_TYPE
#endif

#define MPI_MAX_INT std::numeric_limits<int32_t>::max()


  static inline int MPIRank(MPI_Comm comm = MPI_COMM_WORLD) {

#ifdef CQ_ENABLE_MPI
    int rank;
    MPI_Comm_rank(comm,&rank);
    return rank;
#else
    return comm.rank();
#endif

  }



  static inline int MPISize(MPI_Comm comm = MPI_COMM_WORLD) {

#ifdef CQ_ENABLE_MPI
    int size;
    MPI_Comm_size(comm,&size);
    return size;
#else
    return comm.size();
#endif

  }


  static inline MPI_Comm MPICommSplit(MPI_Comm comm, int color, int key) {

    MPI_Comm c;
#ifdef CQ_ENABLE_MPI
    MPI_Comm_split(comm,color,key,&c);
#endif
    return c;

  }

  static inline void MPICommFree(MPI_Comm &comm) {

#ifdef CQ_ENABLE_MPI
    if( comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
#endif

  }

#ifdef ENABLE_BCAST_COUNTER
  extern int bcastCounter;
#endif

  template <typename T>
  static inline void MPIBCast(T* msg, size_t count, int root, MPI_Comm c) {

#ifdef CQ_ENABLE_MPI

#ifdef ENABLE_BCAST_COUNTER
    bcastCounter++;
#endif
    int int_count = count;
    int bcast_count = (count > MPI_MAX_INT) ? MPI_MAX_INT : std::min(int_count, MPI_MAX_INT);
    MPI_Bcast(msg, bcast_count, mpi_data_type<T>(), root, c);
    if (bcast_count < count) MPIBCast(msg + bcast_count, count - bcast_count, root, c);
#endif

  }

  template <typename T>
  static inline void MPIBCast(T& msg, int root, MPI_Comm c) {
    MPIBCast(&msg,1,root,c);
  }

  template <>
  inline void MPIBCast(bool& msg, int root, MPI_Comm c) {
    int i = msg;
    MPIBCast(i, root, c);
    msg = bool(i);
  }
  
  template <typename T>
  static inline void MPIIBCast(T* msg, int count, int root, MPI_Comm c, MPI_Request* r) {
#ifdef CQ_ENABLE_MPI
    // might be wrong if count > std::numeric_limit<int32_t>::max() for mpich
    MPI_Ibcast(msg, count, mpi_data_type<T>(), root, c, r);
#endif
  }

  template <typename T>
  static inline std::vector<MPI_Request> MPIIBCast(T* msg, size_t count, int root, MPI_Comm c) {
    constexpr size_t max_T = MPI_MAX_INT / sizeof(T);
    size_t n_requests = (count + max_T - 1) / max_T;
    std::vector<MPI_Request> requests(n_requests);
    for (auto i = 0ul; i < n_requests; i++) {
      int count_i = (i == n_requests - 1) ? count - (max_T * (n_requests - 1)) : max_T;
      MPIIBCast(msg + i * max_T, count_i, root, c, &requests[i]);
    }
    return requests;
  }

  template <typename T>
  static inline void MPIIsend(const T* buf, int count, int dst, int tag, MPI_Comm comm, MPI_Request* req) {
  #ifdef CQ_ENABLE_MPI
    MPI_Isend(buf, count, mpi_data_type<T>(), dst, tag, comm, req);
  #endif
  }

  template <typename T>
  static inline std::vector<MPI_Request> MPIIsend(const T* buf, size_t count, int dst, int tag, MPI_Comm comm) {
    constexpr size_t max_T = MPI_MAX_INT / sizeof(T);
    size_t n_chunks = (count + max_T - 1) / max_T;
    std::vector<MPI_Request> requests(n_chunks);
    for (size_t i = 0; i < n_chunks; ++i) {
      size_t offset = i * max_T;
      int count_i = static_cast<int>(std::min(max_T, count - offset));
      MPIIsend(buf + offset, count_i, dst, tag + static_cast<int>(i), comm, &requests[i]);
    }
    return requests;
  }

  template <typename T>
  static inline void MPIIrecv(T* buf, int count, int src, int tag, MPI_Comm comm, MPI_Request* req) {
  #ifdef CQ_ENABLE_MPI
    MPI_Irecv(buf, count, mpi_data_type<T>(), src, tag, comm, req);
  #endif
  }

  template <typename T>
  static inline std::vector<MPI_Request> MPIIrecv(T* buf, size_t count, int src, int tag, MPI_Comm comm) {
    constexpr size_t max_T = MPI_MAX_INT / sizeof(T);
    size_t n_chunks = (count + max_T - 1) / max_T;
    std::vector<MPI_Request> requests(n_chunks);
    for (size_t i = 0; i < n_chunks; ++i) {
      size_t offset = i * max_T;
      int count_i = static_cast<int>(std::min(max_T, count - offset));
      MPIIrecv(buf + offset, count_i, src, tag + static_cast<int>(i), comm, &requests[i]);
    }
    return requests;
  }



  static inline std::vector<MPI_Status> MPIWait(std::vector<MPI_Request>& requests) {
    std::vector<MPI_Status> statuses(requests.size());
    for (auto i = 0ul; i < requests.size(); i++) {
      MPI_Wait(&requests[i], &statuses[i]);
    }
    return statuses;
  }
  
  template <typename T>
  static inline void MPIReduce(const T* in, int n, T* out, int root, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    int nReduce = std::min(n, MPI_MAX_INT);  
    MPI_Reduce(in, out, nReduce, mpi_data_type<T>(), MPI_SUM, root, c);
    if (nReduce < n) MPIReduce(in + nReduce, n - nReduce, out + nReduce, root, c);
#else
    std::copy_n(in, n, out);
#endif
  }
  
  template <typename T>
  static inline T MPIReduce(const T& x, int root, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    T out; MPIReduce(&x, 1, &out, root, c); 
    return out;
#else
    return x;
#endif
  }
  
  template <typename T>
  static inline void MPIAllReduce(const T* in, int n, T* out, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    int nReduce = std::min(n, MPI_MAX_INT);
    if (in == out) {
      MPI_Allreduce(MPI_IN_PLACE, out, nReduce, mpi_data_type<T>(), MPI_SUM, c);
    } else {
      MPI_Allreduce(in, out, nReduce, mpi_data_type<T>(), MPI_SUM, c);
    }
    if (nReduce < n) MPIAllReduce(in + nReduce, n - nReduce, out + nReduce, c);
#else
    std::copy_n(in, n, out);
#endif
  }
  
  
  template <typename T>
  static inline T MPIAllReduce(const T& x, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    T out;
    MPIAllReduce(&x, 1, &out, c);
    return out;
#else
    return x;
#endif
  }

  template <typename T>
  static inline void MPIScatterV(const T* x, const std::vector<size_t>& sizes, T* out, size_t recv_size, int root, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    int rankId = 0, numRanks = 1;
    MPI_Comm_rank(c, &rankId);
    MPI_Comm_size(c, &numRanks);
    constexpr size_t MAX_ELEMENTS_PER_ROUND = static_cast<size_t>(MPI_MAX_INT) / sizeof(T);
    // Build 64‑bit displacement array, to track where each rank starts in the global array
    std::vector<size_t> disp64(numRanks, 0);
    for (int r = 1; r < numRanks; ++r) disp64[r] = disp64[r - 1] + sizes[r - 1];
    // Keep track of how many elements each rank still needs
    std::vector<size_t> remainingToSend = sizes;
    // Keep track of where the next unsent element is for each rank in the global array
    std::vector<size_t> nextReadPos     = disp64;
    // Keep track of the offset into the receiving buffer
    size_t elementsReceivedLocally = 0;

    // Loop over the global array, in a sliding window fashion 
    while (true) {
      // If every rank has received all the elements, break
      size_t largestRemainder = 0;
      for (size_t left : remainingToSend) largestRemainder = std::max(largestRemainder, left);
      if (largestRemainder == 0) break;

      // Choose window size ≤ MAX_ELEMENTS_PER_ROUND
      const size_t windowSize = std::min(MAX_ELEMENTS_PER_ROUND, largestRemainder);

      // Set windowStart = smallest nextReadPos among ranks that still need data.
      // This is the first element of the window that will be scattered in *this* round.
      size_t windowStart = std::numeric_limits<size_t>::max();
      for (int r = 0; r < numRanks; ++r)
        if (remainingToSend[r] > 0) windowStart = std::min(windowStart, nextReadPos[r]);
      // Slide base pointer
      const T* sendBase = x + windowStart;

      // Build 32‑bit arrays for counts and displacements, relative to sendBase.
      std::vector<int> sendCounts32(numRanks, 0);
      std::vector<int> displs32    (numRanks, 0);
      const size_t windowEnd = windowStart + windowSize;

      // Loop over ranks to calculate send counts and displacements
      for (int r = 0; r < numRanks; ++r) {
        if (remainingToSend[r] == 0) continue;           // rank has no more data to send
        if (nextReadPos[r] >= windowEnd) continue;       // rank is not in this current window

        size_t sendNow  = std::min(windowEnd - nextReadPos[r], remainingToSend[r]);
        sendCounts32[r] = static_cast<int>(sendNow);                
        displs32[r]     = static_cast<int>(nextReadPos[r] - windowStart); //
      }

      MPI_Scatterv(sendBase,
                    sendCounts32.data(), displs32.data(), mpi_data_type<T>(),
                    out + elementsReceivedLocally, sendCounts32[rankId], mpi_data_type<T>(),
                    root, c);

      // Update trackers
      for (int r = 0; r < numRanks; ++r) {
        nextReadPos[r]     += static_cast<size_t>(sendCounts32[r]);
        remainingToSend[r] -= static_cast<size_t>(sendCounts32[r]);
      }
      elementsReceivedLocally += static_cast<size_t>(sendCounts32[rankId]);
    }

    // Sanity check
    if (elementsReceivedLocally != recv_size)
      throw std::runtime_error{
          "Error in MPIScatterV: Received " + std::to_string(elementsReceivedLocally)
          + " elements, Expected " + std::to_string(recv_size)
          + " on rank " + std::to_string(rankId)
      };
#else
    std::copy_n(x, recv_size, out);
#endif
}
   
  template <typename T>
  static inline std::vector<T> MPIGather(const T& x, 
      int root, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
   std::vector<T> out;
   if (MPIRank(c) == root) out.resize(MPISize(c));
   MPI_Gather(&x, 1, mpi_data_type<T>(), out.data(), 1, mpi_data_type<T>(), root, c);
   return out;
#else
   return {x}; 
#endif
  }

  // might be wrong if count > std::numeric_limit<int32_t>::max() for mpich
  template <typename T>
  static inline void MPIGatherV(const T* x,
      size_t size,
      T* out,
      const std::vector<size_t>& recv_sizes,
      int root,
      MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
   std::vector<int> _sizes(recv_sizes.begin(), recv_sizes.end());
   std::vector<int> _displs(_sizes.size());
   std::exclusive_scan(_sizes.begin(), _sizes.end(), _displs.begin(), 0);
   MPI_Gatherv(x, size, mpi_data_type<T>(), out, _sizes.data(), 
     _displs.data(), mpi_data_type<T>(), root, c);
#else
   std::copy_n(x, size, out);
#endif
  }

  template <typename T>
  static inline void MPIAllGatherV(const T* x, size_t local_size, T* out,
                                  const std::vector<size_t>& recv_sizes, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    int rankId  = 0, numRanks = 1;
    MPI_Comm_rank(c, &rankId);
    MPI_Comm_size(c, &numRanks);

    constexpr size_t MAX_CHUNK = static_cast<size_t>(MPI_MAX_INT) / sizeof(T);

    // Build 64‑bit displacement array into the global array
    std::vector<size_t> disp64(numRanks, 0);
    for (int r = 1; r < numRanks; ++r) disp64[r] = disp64[r - 1] + recv_sizes[r - 1];
    const size_t max_recv = *std::max_element(recv_sizes.begin(), recv_sizes.end());
    size_t offset = 0;           
    
    while (offset < max_recv) {
      // Choose chunk so that *every* displacement & count fits in int
      const size_t chunk = std::min(MAX_CHUNK, max_recv - offset);
      size_t windowEnd = offset + chunk;

      // Build 32‑bit arrays for this chunk
      std::vector<int> chunkCounts32(numRanks, 0);
      std::vector<int> chunkDispls32(numRanks, 0);
      for (int r = 0; r < numRanks; ++r) {
        // skip ranks that begins after this current window
        if (disp64[r] >= windowEnd) continue;
        size_t remaining = (recv_sizes[r] > offset) ? recv_sizes[r] - offset : 0;
        size_t sendNow   = std::min(chunk, remaining);
        chunkCounts32[r] = static_cast<int>(sendNow);
        chunkDispls32[r] = static_cast<int>(disp64[r] > offset ? disp64[r] - offset : 0);
      }

      // Choose a safe send buffer for *this* rank
      const T* sendbuf = (offset < local_size) ? x + offset: x;
      
      MPI_Allgatherv(sendbuf, chunkCounts32[rankId], mpi_data_type<T>(),
                      out + offset,
                      chunkCounts32.data(), chunkDispls32.data(),
                      mpi_data_type<T>(), c);
      offset += chunk;
    }
#else
    std::copy_n(x, local_size, out);
#endif
}

template <typename T>
static inline void MPIAlltoallv(const T                         *sendbuf,
                                const std::vector<size_t>       &sendCounts64,
                                const std::vector<size_t>       &sendDispls64,
                                T                               *recvbuf,
                                const std::vector<size_t>       &recvCounts64,
                                const std::vector<size_t>       &recvDispls64,
                                MPI_Comm                        comm) {
#ifdef CQ_ENABLE_MPI
  int rankId = 0, numRanks = 1;
  MPI_Comm_rank(comm, &rankId);
  MPI_Comm_size(comm, &numRanks);

  // Sanity checks
  assert(sendCounts64.size() == numRanks);
  assert(recvCounts64.size() == numRanks);
  assert(sendDispls64.size() == numRanks);
  assert(recvDispls64.size() == numRanks);

  // Hard per‑peer cap to limit temporary buffers and RDMA registration size. 128 MiB ÷ sizeof(T)
  const size_t PER_PEER_CAP = 128ULL * 1024 * 1024 / sizeof(T);  // 128 MiB

  // 32‑bit vectors required by MPI‑3 signature
  std::vector<int> sendCounts32(numRanks, 0), recvCounts32(numRanks, 0);
  std::vector<int> sendDispls32(numRanks, 0), recvDispls32(numRanks, 0);

  // For single-shot fast path, all local counts/displacements must fit the 32‑bit and be less than PER_PEER_CAP
  bool canUseFastPath = true;
  for (int p = 0; p < numRanks && canUseFastPath; ++p) {
    canUseFastPath &= (sendCounts64[p] <= PER_PEER_CAP) &&
                      (recvCounts64[p] <= PER_PEER_CAP) &&
                      (sendDispls64[p] <= static_cast<size_t>(MPI_MAX_INT)) &&
                      (recvDispls64[p] <= static_cast<size_t>(MPI_MAX_INT));
  }

  // Collect over all ranks to make the decision
  int local_ok = canUseFastPath ? 1 : 0, global_ok = 0;
  MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_LAND, comm);

  // Fast path: one direct MPI_Alltoallv call
  if (global_ok) {
    for (int p = 0; p < numRanks; ++p) {
      sendCounts32[p] = static_cast<int>(sendCounts64[p]);
      recvCounts32[p] = static_cast<int>(recvCounts64[p]);
      sendDispls32[p] = static_cast<int>(sendDispls64[p]);
      recvDispls32[p] = static_cast<int>(recvDispls64[p]);
    }
    MPI_Alltoallv(sendbuf,  sendCounts32.data(), sendDispls32.data(), mpi_data_type<T>(),
                  recvbuf,  recvCounts32.data(), recvDispls32.data(), mpi_data_type<T>(),
                  comm);
    return;
  }

  // Chunked path: multiple MPI_Alltoallv calls after packing/unpacking
  std::vector<size_t> sendRemaining = sendCounts64;  // Remaining elements to send
  std::vector<size_t> recvRemaining = recvCounts64;  // Remaining elements to receive
  std::vector<size_t> sendOffset(numRanks, 0);       // Offset into original user send buffer
  std::vector<size_t> recvOffset(numRanks, 0);       // Offset into original user recv buffer

  // Progress bookkeeping
  size_t totalSendAmount = std::accumulate(sendCounts64.begin(), sendCounts64.end(), size_t(0));
  size_t totalRecvAmount = std::accumulate(recvCounts64.begin(), recvCounts64.end(), size_t(0));
  size_t totalSendComplete = 0, totalRecvComplete = 0;

  // Print progress if rank 0
  bool print_progress = (rankId == 0);

  for (int round = 0; ; ++round) {
    if (print_progress) {
      // show both percent and MiB to give meaningful resolution
      double sentMiB     = totalSendComplete * sizeof(T) / (1024.0*1024.0);
      double needMiB     = totalSendAmount   * sizeof(T) / (1024.0*1024.0);
      double recvMiB     = totalRecvComplete * sizeof(T) / (1024.0*1024.0);
      double needRecvMiB = totalRecvAmount   * sizeof(T) / (1024.0*1024.0);

      double sPct = needMiB     ? 100.0 * sentMiB / needMiB     : 100.0;
      double rPct = needRecvMiB ? 100.0 * recvMiB / needRecvMiB : 100.0;

      std::cout << std::fixed << std::setprecision(2)
                << "MPIAlltoallv: "
                << std::left  << std::setw(9) << ("Round " + std::to_string(round))
                << " - Send: " << std::right << std::setw(6) << sPct << "% ("
                << sentMiB << "/" << needMiB << " MiB), "
                << "Receive: " << std::setw(6) << rPct << "% ("
                << recvMiB << "/" << needRecvMiB << " MiB)" << std::endl;
    }

    bool sendDone = std::all_of(sendRemaining.begin(), sendRemaining.end(),
                                [](size_t v){ return v == 0; });
    bool recvDone = std::all_of(recvRemaining.begin(), recvRemaining.end(),
                                [](size_t v){ return v == 0; });
    if (sendDone && recvDone) break;

    // Calculate chunk sizes for this round
    size_t roundSend = 0, roundRecv = 0;
    for (int p = 0; p < numRanks; ++p) {
      size_t s = std::min(sendRemaining[p], PER_PEER_CAP);
      size_t r = std::min(recvRemaining[p], PER_PEER_CAP);
      sendCounts32[p] = static_cast<int>(s);
      recvCounts32[p] = static_cast<int>(r);
      roundSend += s;
      roundRecv += r;
    }

    // Safeguard: if no data to send/receive this round, exit
    if (roundSend == 0 && roundRecv == 0) break;

    // Calculate 32-bit displacements
    sendDispls32[0] = recvDispls32[0] = 0;
    for (int p = 1; p < numRanks; ++p) {
      sendDispls32[p] = sendDispls32[p-1] + sendCounts32[p-1];
      recvDispls32[p] = recvDispls32[p-1] + recvCounts32[p-1];
    }

    // Sanity check
    assert(sendDispls32.back() + sendCounts32.back() <= MPI_MAX_INT);
    assert(recvDispls32.back() + recvCounts32.back() <= MPI_MAX_INT);

    std::vector<T> tempSendBuf(roundSend);
    std::vector<T> tempRecvBuf(roundRecv);

    // Pack user data into tempSendBuf
    if (roundSend) {
      size_t offs = 0;
      for (int p = 0; p < numRanks; ++p) {
        if (sendCounts32[p]) {
          const T* src = sendbuf + sendDispls64[p] + sendOffset[p];
          std::copy_n(src, sendCounts32[p], tempSendBuf.data() + offs);
          offs += sendCounts32[p];
        }
      }
    }

    // MPI_Alltoallv: exchange data with other ranks
    MPI_Alltoallv( roundSend ? tempSendBuf.data() : nullptr,
                   sendCounts32.data(), sendDispls32.data(), mpi_data_type<T>(),
                   roundRecv ? tempRecvBuf.data() : nullptr,
                   recvCounts32.data(), recvDispls32.data(), mpi_data_type<T>(),
                   comm );

    // Unpack user data from tempRecvBuf
    if (roundRecv) {
      for (int p = 0; p < numRanks; ++p) {
        if (recvCounts32[p]) {
          T* dst = recvbuf + recvDispls64[p] + recvOffset[p];
          std::copy_n(tempRecvBuf.data() + recvDispls32[p], recvCounts32[p], dst);
        }
      }
    }

    // Bookkeeping for next round
    for (int p = 0; p < numRanks; ++p) {
      sendRemaining[p] -= static_cast<size_t>(sendCounts32[p]);
      recvRemaining[p] -= static_cast<size_t>(recvCounts32[p]);
      sendOffset[p]    += static_cast<size_t>(sendCounts32[p]);
      recvOffset[p]    += static_cast<size_t>(recvCounts32[p]);
      totalSendComplete += static_cast<size_t>(sendCounts32[p]);
      totalRecvComplete += static_cast<size_t>(recvCounts32[p]);
    }
  } // end round loop

  if (print_progress) std::cout << "MPIAlltoallv: Transfer complete" << std::endl;
#else
  /* serial fallback */
  if (!sendCounts64.empty()) std::copy_n(sendbuf, sendCounts64[0], recvbuf);
#endif
}







  static inline bool MPIAnyOf(bool x, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    int i = x ? 1 : 0;
    i = MPIAllReduce(i, c);
    return bool(i);
#else
    return x;
#endif
  }

#define ROOT_ONLY(comm) if(MPIRank(comm) != 0) return;

  static inline MPI_Comm CreateRootComm(MPI_Comm c) {

#ifdef CQ_ENABLE_MPI
    return MPICommSplit(c, (MPIRank(c) == 0) ? 1 : MPI_UNDEFINED, 0);
#else
    return MPI_COMM_WORLD;
#endif

  }

#ifdef CQ_ENABLE_SPARSE
  template <typename T>
  class MPIWin {

    std::vector<MPI_Win> windows_;
    // std::vector<size_t> all_sizes_;
    const size_t max_window_size_ = MPI_MAX_INT / sizeof(T);   
    size_t size_;
    std::vector<MPI_Request> requests_;

   public:
    MPIWin() = delete;
    MPIWin(T* msg, size_t size, MPI_Comm c): size_(size) {
      size_t n_windows = (size_ + max_window_size_ - 1) / max_window_size_;
      size_t n_windows_local = n_windows;
      MPI_Allreduce(MPI_IN_PLACE, &n_windows, 1, mpi_data_type<size_t>(), MPI_MAX, c);   
      windows_.resize(n_windows);
      // std::cout << "n_windows = " << n_windows << std::endl;
      // std::cout << "size_ = " << size_ << std::endl;
      for (auto i = 0ul; i < n_windows; i++) {
        T* msg_i = nullptr;
        size_t s_i = 0ul;
        if (i < n_windows_local) {
          s_i = (i == n_windows_local - 1) ? size_ - (max_window_size_ * (n_windows_local - 1)) : max_window_size_;
          msg_i = msg + i * max_window_size_; 
        }
        MPI_Win_create(msg_i, s_i * sizeof(T), sizeof(T), MPI_INFO_NULL, c, &windows_[i]);
        MPI_Win_fence(0, windows_[i]);
      }
    }
    ~MPIWin() {}

    size_t size() { return size_; }

    // locks
    void lock(int lock_type, int rank, int assert) {
      for (auto& win: windows_) MPI_Win_lock(lock_type, rank, assert, win);
    }
    void unlock(int rank) {
      for (auto& win: windows_) MPI_Win_unlock(rank, win);
    }
    void lock_all(int assert) {
      for (auto& win : windows_) MPI_Win_lock_all(assert, win);
    }
    void unlock_all() {
      for (auto& win : windows_) MPI_Win_unlock_all(win);
    }
    void flush(size_t rank) {
      for (auto& win : windows_) MPI_Win_flush(rank, win);
    } 
    void free() {
      for (auto& win : windows_)  MPI_Win_free(&win);
    }
    void get(T* buffer, size_t len, size_t target, size_t displacement) {
      size_t w = displacement / max_window_size_;
      size_t d = displacement % max_window_size_;
      size_t l = std::min(max_window_size_ - d, len);
      MPI_Get(buffer, l, mpi_data_type<T>(), target, d, l, mpi_data_type<T>(), windows_[w]);
      if (l != len) get(buffer + l, len - l, target, displacement + l); 
    } // get
    void put(const T* buffer, size_t len, size_t target, size_t displacement) {
      size_t w = displacement / max_window_size_;
      size_t d = displacement % max_window_size_;
      size_t l = std::min(max_window_size_ - d, len);
      MPI_Put(buffer, l, mpi_data_type<T>(), target, d, l, mpi_data_type<T>(), windows_[w]);
      if (l != len) put(buffer + l, len - l, target, displacement + l); 
    } // put

    void rget(T* buffer, size_t len, size_t target, size_t displacement) {
      size_t w = displacement / max_window_size_;
      size_t d = displacement % max_window_size_;
      size_t l = std::min(max_window_size_ - d, len);
      requests_.push_back(MPI_REQUEST_NULL);
      MPI_Rget(buffer, l, mpi_data_type<T>(), target, d, l, mpi_data_type<T>(), windows_[w], &requests_.back());
      if (l != len) rget(buffer + l, len - l, target, displacement + l);
    } // rget
    void rput(T* buffer, size_t len, size_t target, size_t displacement) {
      size_t w = displacement / max_window_size_;
      size_t d = displacement % max_window_size_;
      size_t l = std::min(max_window_size_ - d, len);
      requests_.push_back(MPI_REQUEST_NULL);
      MPI_Rput(buffer, l, mpi_data_type<T>(), target, d, l, mpi_data_type<T>(), windows_[w], &requests_.back());
      if (l != len) rput(buffer + l, len - l, target, displacement + l);
    } // rput
    void wait() {
      MPIWait(requests_);
      requests_.clear();
    }

  }; // class MPIWin
#endif

}; // namespace ChronusQ

