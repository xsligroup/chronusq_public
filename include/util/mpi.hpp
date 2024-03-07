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


#else // MPI is enabled
  template <typename T> MPI_Datatype mpi_data_type();
  #define REGISTER_MPI_TYPE(CXXTYPE, MPITYPE) \
  template <> inline MPI_Datatype mpi_data_type<CXXTYPE>() { return MPITYPE; }

  REGISTER_MPI_TYPE(double, MPI_DOUBLE)
  REGISTER_MPI_TYPE(int,    MPI_INT   )
  REGISTER_MPI_TYPE(int64_t,MPI_INT64_T)
  REGISTER_MPI_TYPE(size_t, MPI_UINT64_T)

  // For compilers that are known to have `long int` equivalent to `int64_t`, avoid redefinition
  #if !defined(LONG_INT_IS_INT64_T)
  REGISTER_MPI_TYPE(long int, MPI_LONG)
  #endif

  REGISTER_MPI_TYPE(std::complex<double>, MPI_C_DOUBLE_COMPLEX)

  #undef REGISTER_MPI_TYPE
#endif


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
  static inline void MPIBCast(T* msg, int count, int root, MPI_Comm c) {

#ifdef CQ_ENABLE_MPI

#ifdef ENABLE_BCAST_COUNTER
    bcastCounter++;
#endif
    MPI_Bcast(msg, count, mpi_data_type<T>(), root, c);
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
  static inline void MPIReduce(const T* in, size_t n, T* out, int root, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    MPI_Reduce(in, out, n, mpi_data_type<T>(), MPI_SUM, root, c);
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
  static inline void MPIAllReduce(const T* in, size_t n, T* out, MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
    MPI_Allreduce(in, out, n, mpi_data_type<T>(), MPI_SUM, c);
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
  static inline void MPIScatterV(const T* x, 
      const std::vector<size_t>& sizes,
      T* out,
      size_t recv_size,
      int root,
      MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
   std::vector<int> _sizes(sizes.begin(), sizes.end());
   std::vector<int> _displs(sizes.size());
   std::exclusive_scan(_sizes.begin(), _sizes.end(), _displs.begin(), 0);
   MPI_Scatterv(x, _sizes.data(), _displs.data(), mpi_data_type<T>(),
     out, recv_size, mpi_data_type<T>(), root, c);
#else
   std::copy_n(x, recv_size, out);
#endif
  }

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
  static inline void MPIAllGatherV(const T* x,
      size_t size,
      T* out,
      const std::vector<size_t>& recv_sizes,
      MPI_Comm c) {
#ifdef CQ_ENABLE_MPI
   std::vector<int> _sizes(recv_sizes.begin(), recv_sizes.end());
   std::vector<int> _displs(_sizes.size());
   std::exclusive_scan(_sizes.begin(), _sizes.end(), _displs.begin(), 0);
   MPI_Allgatherv(x, size, mpi_data_type<T>(), _sizes.data(),
     _displs.data(), out, mpi_data_type<T>(), c);
#else
   std::copy_n(x, size, out);
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
}; // namespace ChronusQ

