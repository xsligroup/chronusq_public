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

#include <func.hpp>
#include <cerr.hpp>
#include <memmanager.hpp>
#include <itersolver.hpp>
#include <itersolver/solvervectorsimpl.hpp>
#include <util/files.hpp>
#include <util/timer.hpp>

#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/eig.hpp>

using namespace ChronusQ;


template <typename ReadT, typename EigT>
void DAVIDSON_RAWVECTORS_TEST(size_t nRoots, size_t m, size_t kG,
  std::string fname, bool doPre = false, double conver = 1e-10, 
  double etol = 8e-8, size_t block_size = 128) {
  
  MPI_Barrier(MPI_COMM_WORLD);

  std::string refName( FUNC_REFERENCE + fname );

  bool isRoot = MPIRank(MPI_COMM_WORLD) == 0;
  bool isMPI  = MPISize(MPI_COMM_WORLD) > 1;

  
  if( isRoot ) {
    std::cout << "  * Will use " << MPISize(MPI_COMM_WORLD) << 
      " MPI Processes" << std::endl;
    std::cout << "DAVIDSON_RAWVECTORS_TEST with file:" << refName << std::endl;
  }

  SafeFile matFile(refName,true);

  CQMemManager mem(2e9,256);
  
  
  size_t N = 0;
  if (isRoot) { 
    // Dummy checks on
    auto dims = matFile.getDims("/matrix");
    EXPECT_TRUE( dims.size() == 2 ) << "dims.size() = " << dims.size();
    EXPECT_TRUE( dims[0] == dims[1] ) << " dims[0] = " << dims[0] <<", dims[1] =" << dims[1];
    N = dims[0];
  }
  MPIBCast(N, 0, MPI_COMM_WORLD);

  int64_t MLoc = N, NLoc = N;

#ifdef CQ_ENABLE_MPI
  std::shared_ptr<scalapackpp::BlockCyclicDist2D> grid = nullptr;
  scalapackpp::scalapack_desc descA;

  if( isMPI ) {
    grid = std::make_shared<scalapackpp::BlockCyclicDist2D>(
        blacspp::Grid::square_grid(MPI_COMM_WORLD),
        block_size, block_size);

    std::tie(MLoc,NLoc) = grid->get_local_dims(N,N);
    descA = grid->descinit_noerror(N,N,MLoc);
  }
#endif

  // Allocate and read the matrix
  ReadT * A = (MLoc and NLoc) ? mem.malloc<ReadT>(MLoc * NLoc) : nullptr;

  ReadT *AREAD = nullptr;
  ReadT *DIAG  = nullptr;
  if( isRoot ) {

    AREAD = isMPI ? mem.malloc<ReadT>(N*N) : A;
    matFile.readData("/matrix",AREAD); 

    DIAG  = mem.malloc<ReadT>(N);
    for( size_t k = 0; k < N; k++ ) DIAG[k] = AREAD[k*(N+1)];

  }

  if( isRoot ) std::cout << "  * Load Data Successfully " << std::endl;

#ifdef CQ_ENABLE_MPI
  if( isMPI ) grid->scatter(N,N,AREAD,N,A,MLoc,0,0);
#endif


  EigT *ALOC = reinterpret_cast<EigT*>(A);

#ifdef CQ_ENABLE_MPI
  if( isMPI ) {

    if( AREAD ) mem.free(AREAD);


    if( not std::is_same<EigT,ReadT>::value ) {
      ALOC = mem.malloc<EigT>(MLoc * NLoc);
      std::copy_n(A,MLoc*NLoc,ALOC);
      mem.free(A);
    }

  }
#endif


typename Davidson<EigT>::LinearTrans_t func = [&]( size_t nVec, SolverVectors<EigT> &V,
    SolverVectors<EigT> &AV) {

    auto V_ptr = tryGetRawVectorsPointer(V);
    auto AV_ptr = tryGetRawVectorsPointer(AV);

#ifdef CQ_ENABLE_MPI
    if( isMPI ) {

      int64_t MLoc_V = N, NLoc_V = nVec;
      std::tie(MLoc_V, NLoc_V) = grid->get_local_dims(N,nVec);
      auto descV = grid->descinit_noerror(N,nVec,MLoc_V);

      bool alloc = MLoc_V and NLoc_V;

      EigT *VLOC = nullptr, *AVLOC = nullptr;

      if( alloc ) {

        VLOC  = mem.malloc<EigT>(MLoc_V * NLoc_V);
        AVLOC = mem.malloc<EigT>(MLoc_V * NLoc_V);

      }

      grid->scatter(N,nVec,V_ptr,N,VLOC,MLoc_V,0,0);

      Gemm_MPI('N','N',N,nVec,N,EigT(1.),ALOC,1,1,descA,VLOC,1,1,descV,
          EigT(0.),AVLOC,1,1,descV);

      grid->gather(N,nVec,AV_ptr,N,AVLOC,MLoc_V,0,0);

    } else 
#endif
  blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
      N,nVec,N,EigT(1.),A,N,V_ptr,N,EigT(0.),AV_ptr,N);

  };

  typename Davidson<EigT>::LinearTrans_t PC = [&]( size_t nVec, SolverVectors<EigT> &V,
      SolverVectors<EigT> &AV) {

    ROOT_ONLY(MPI_COMM_WORLD);

    auto V_ptr = tryGetRawVectorsPointer(V);
    auto AV_ptr = tryGetRawVectorsPointer(AV);
    
    if( V_ptr != AV_ptr )
      AV.set_data(0, nVec, V, 0);

    if( doPre )
    // Scale by inverse diagonals
    for( auto k = 0ul; k < N; k++ ) 
      blas::scal(nVec, EigT(1.) / DIAG[k], AV_ptr + k, N);
  };

#ifndef _CQ_GENERATE_TESTS

  size_t nThreads = omp_get_num_threads();
  ProgramTimer::initialize("Davidson test", nThreads);
  Davidson<EigT> davidson(MPI_COMM_WORLD,mem,N,5,128,conver,nRoots,
    func,PC);

  davidson.setM(m);
  davidson.setkG(kG);
  
  davidson.run();

  ROOT_ONLY(MPI_COMM_WORLD);

  dcomplex *refW = mem.malloc<dcomplex>(N);
  matFile.readData("/W",refW);

  const dcomplex *W = davidson.eigVal();

  for(auto k = 0ul; k < nRoots; k++) {
    double diff1 = std::abs((W[k] - refW[k])/refW[k]);
    double diff2 = std::abs((W[k] - std::conj(refW[k]))/refW[k]);
    EXPECT_TRUE( diff1 < etol or diff2 < etol) <<
      "DIFF1 = " << diff1 << ", DIFF2 = " << diff2;
  }

#else

  dcomplex *ACMPLX = mem.malloc<dcomplex>(N*N);
  std::copy_n(A,N*N,ACMPLX);

  dcomplex *W = mem.malloc<dcomplex>(N);
  dcomplex *VR = mem.malloc<dcomplex>(N*N);
  dcomplex *VL = mem.malloc<dcomplex>(N*N);

  GeneralEigen('V','V',N,ACMPLX,N,W,VL,N,VR,N);

  matFile.safeWriteData("/W",W,{N});

#endif

}

template <typename ReadT, typename EigT>
void DAVIDSON_DISTRIBUTEDVECTORS_TEST(size_t nRoots, size_t m, size_t kG,
  std::string fname, bool doPre = false, double conver = 1e-10, 
  double etol = 8e-8) {

  MPI_Barrier(MPI_COMM_WORLD);
  std::string refName( FUNC_REFERENCE + fname );

  bool isRoot = MPIRank(MPI_COMM_WORLD) == 0;
  bool isMPI  = MPISize(MPI_COMM_WORLD) > 1;


  if( isRoot ) {
    std::cout << "DAVIDSON_DISTRIBUTEDVECTORS_TEST with file:" << refName 
             << ", doPre = " << doPre << std::endl;
    std::cout << "  * Will use " << MPISize(MPI_COMM_WORLD) << 
    " MPI Processes" << std::endl;
  }

  SafeFile matFile(refName,true);

  CQMemManager mem(2e9,256);
  MPI_Comm comm(MPI_COMM_WORLD);

  size_t N = 0;
  if (isRoot) {
    // Dummy checks on
    auto dims = matFile.getDims("/matrix");
    EXPECT_TRUE( dims.size() == 2 ) << "dims.size() = " << dims.size();
    EXPECT_TRUE( dims[0] == dims[1] ) << " dims[0] = " << dims[0] <<", dims[1] =" << dims[1];
    N = dims[0];
  }
  MPIBCast(N, 0, comm);

  size_t N2 = N * N;

  ReadT* AREAD = nullptr;
  EigT* ADIAG = mem.malloc<EigT>(N);
  
  if (isRoot) {
    AREAD = mem.malloc<ReadT>(N2);
    matFile.readData("/matrix", AREAD);
    
    for( size_t k = 0; k < N; k++ ) ADIAG[k] = AREAD[k*(N+1)];
  }

  MPIBCast(ADIAG, N, 0, comm);

  EigT* A = reinterpret_cast<EigT*>(AREAD); 
  
  typename Davidson<EigT>::LinearTrans_t func = 
      [&]( size_t nVec, SolverVectors<EigT> &V,
        SolverVectors<EigT> &AV) {

        // copy the V out
        EigT* VRaw = mem.malloc<EigT>(N*nVec);
        EigT* AVRaw = mem.malloc<EigT>(N*nVec);
        
        tryDowncastReferenceTo<DistributedVectors<EigT>>(V,
                                                         [&] (auto& VRef, size_t shiftV) {
              VRef.gather(shiftV, nVec, VRaw, N, 0);
            }
        );
        
        if (isRoot) {
          // prettyPrintSmart(std::cout, "VRaw", VRaw, N, nVec, N);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
              N,nVec,N,EigT(1.),A,N,VRaw,N,EigT(0.),AVRaw,N);
          // prettyPrintSmart(std::cout, "AVRaw", AVRaw, N, nVec, N);
        }
        
        MPI_Barrier(comm);
        
        // scatter AVRaw Back 
        tryDowncastReferenceTo<DistributedVectors<EigT>>(AV,
                                                         [&] (auto& AVRef, size_t shiftAV) {
              AVRef.scatter(shiftAV, nVec, AVRaw, N, 0); 
            }
         );

      };

  typename Davidson<EigT>::LinearTrans_t PC = 
      [&](size_t nVec, SolverVectors<EigT> &V,
        SolverVectors<EigT> &AV) {

        AV.set_data(0, nVec, V, 0);
        
        if( doPre ) {
          tryDowncastReferenceTo<DistributedVectors<EigT>>(AV,
                                                           [&] (auto& AVRef, size_t shiftAV) {
                const EigT* ADIAG_ptr = ADIAG + AVRef.localOffset();
                auto AV_ptr = AVRef.getLocalPtr(shiftAV);
                for (auto i = 0ul; i < AVRef.localLength(); ++i, ++ADIAG_ptr, ++AV_ptr) {
                  blas::scal(nVec, EigT(1.) / (*ADIAG_ptr), AV_ptr, AVRef.localLength()); 
                }
              }
          );              
        }

      };
   
  std::function<std::shared_ptr<SolverVectors<EigT>>(size_t)> distributedVecsGenerator = 
      [&] (size_t nVec) {
         return std::make_shared<DistributedVectors<EigT>>(comm, mem, N, nVec);
      };
  
  size_t nThreads = omp_get_num_threads();
  ProgramTimer::initialize("Davidson test", nThreads);
  Davidson<EigT> davidson(MPI_COMM_WORLD,mem,N,5,128,conver,nRoots,
    func,PC, distributedVecsGenerator);

  davidson.setM(m);
  davidson.setkG(kG);
  
  davidson.run();
  
  ROOT_ONLY(MPI_COMM_WORLD);

  dcomplex *refW = mem.malloc<dcomplex>(N);
  matFile.readData("/W",refW);

  const dcomplex *W = davidson.eigVal();

  for(auto k = 0ul; k < nRoots; k++) {
    double diff1 = std::abs((W[k] - refW[k])/refW[k]);
    double diff2 = std::abs((W[k] - std::conj(refW[k]))/refW[k]);
    EXPECT_TRUE( diff1 < etol or diff2 < etol) <<
      "DIFF1 = " << diff1 << ", DIFF2 = " << diff2;
  }

  if (AREAD) mem.free(AREAD);
  if (ADIAG) mem.free(ADIAG);
}

// 
// Hermittian Davidson
//
TEST(DAVIDSON, DAVIDSON_REAL_HERMITIAN) {
  
  DAVIDSON_RAWVECTORS_TEST<double,double>(3,50,3,"real_Hermitian.hdf5");
  DAVIDSON_RAWVECTORS_TEST<double,double>(3,50,3,"real_Hermitian.hdf5",true);
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<double,double>(3,50,3,"real_Hermitian.hdf5");
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<double,double>(3,50,3,"real_Hermitian.hdf5",true);

}

TEST(DAVIDSON, DAVIDSON_COMPLEX_HERMITIAN) {
  
  DAVIDSON_RAWVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_Hermitian.hdf5");
  DAVIDSON_RAWVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_Hermitian.hdf5",true);
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_Hermitian.hdf5");
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_Hermitian.hdf5",true);

}

// 
// Non-Hermittian
//
TEST(DAVIDSON, DAVIDSON_REAL_NONHERMITIAN) {
  
  DAVIDSON_RAWVECTORS_TEST<double,double>(3,50,3,"real_nonHermitian.hdf5");
  DAVIDSON_RAWVECTORS_TEST<double,double>(3,50,3,"real_nonHermitian.hdf5",true);
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<double,double>(3,50,3,"real_nonHermitian.hdf5");
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<double,double>(3,50,3,"real_nonHermitian.hdf5",true);

}

TEST(DAVIDSON, DAVIDSON_COMPLEX_NONHERMITIAN) {
  
  DAVIDSON_RAWVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_nonHermitian.hdf5");
  DAVIDSON_RAWVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_nonHermitian.hdf5",true);
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_nonHermitian.hdf5");
  DAVIDSON_DISTRIBUTEDVECTORS_TEST<dcomplex,dcomplex>(3,50,3,"complex_nonHermitian.hdf5",true);

}



