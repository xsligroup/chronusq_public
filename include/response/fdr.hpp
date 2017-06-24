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

#include <response/tbase.hpp>


#include <cqlinalg/solve.hpp>
#include <cqlinalg/blas1.hpp>
#include <itersolver.hpp>

#include <lapack.hh>

namespace ChronusQ {

  template <typename T>
  template <typename U>
  void ResponseTBase<T>::runFullFDR(FDResponseResults<T,U> &results) {

    ProgramTimer::tick("Full FDR");

    bool isRoot = MPIRank(comm_) == 0;
    bool isDist = this->genSettings.isDist();

    size_t N = nSingleDim_;
    size_t nOmega = fdrSettings.bFreq.size();
    size_t nRHS   = fdrSettings.nRHS;

    if( this->genSettings.printLevel >= 0 and isRoot ) {
      std::cout << "  ** PERFORMING INCORE LINEAR SOLVES FOR " << nOmega 
                << " FREQUENCIES AND " << nRHS << " RHS\n\n";
      std::cout << "    * PROBLEM DIMENSION = " << N << "\n";
      std::cout << "    * SOLVER            = ";
      if( isDist )
        std::cout << "PXGESV (ScaLAPACK)";
      else
        std::cout << "XGESV (LAPACK)";

      std::cout << "\n\n";
    }

    size_t localDimM = N, localDimN = N;
#ifdef CQ_ENABLE_MPI
    if( isDist ) 
      std::tie(localDimM,localDimN) = fullMatGrid_->get_local_dims(N,N);
#endif
    


    // Allocate space for the shifted matrix
    U* shiftedMat = nullptr;

    if( localDimM and localDimN )
      shiftedMat = CQMemManager::get().malloc<U>(localDimM*localDimN);


    U* distSOL = nullptr;
    int64_t localRHS;
#ifdef CQ_ENABLE_MPI
    scalapackpp::scalapack_desc DescA, DescB;
    if( isDist ) {

      // Allocate space for distributed RHS / Solution
      std::tie(localDimM, localRHS) = fullMatGrid_->get_local_dims(N,nRHS);
      if( localDimM and localRHS)
        distSOL = CQMemManager::get().malloc<U>(localDimM*localRHS);


      DescB = fullMatGrid_->descinit_noerror(N,nRHS,localDimM);
      DescA = fullMatGrid_->descinit_noerror(N,N,localDimM);

    }
#endif

    T* fMatUse = formFullFromMemory();

    // Loop over omegas
    for(size_t iOmega = 0; iOmega < nOmega; iOmega++) {

      ProgramTimer::tick("Omega");

      U omega = results.shifts[iOmega];

      if( this->genSettings.printLevel >= 0 and isRoot )
      std::cout << "  * SOLVING IOMEGA = " << std::setw(5) << iOmega 
                << ": W = " << std::setprecision(8) << std::scientific
                << omega << " (AU) \n";

      U* SOL = results.SOL + iOmega*nRHS*N;


      // Copy over the RHS into the solution storage
      if( isRoot ) std::copy_n(results.RHS,N*nRHS,SOL);

#ifdef CQ_ENABLE_MPI
      if( isDist ) 
        fullMatGrid_->scatter(N,nRHS,SOL,N,distSOL,localDimM,0,0);
#endif

      // Copy over the full matrix and shift
      if( shiftedMat ) std::copy_n(fMatUse,localDimM*localDimN,shiftedMat);

#ifdef CQ_ENABLE_MPI
      if( isDist and shiftedMat )
        for(auto j = 0ul; j < localDimN; j++) 
        for(auto i = 0ul; i < localDimM; i++) {

          int64_t I,J;
          std::tie(I,J) = fullMatGrid_->global_idx(i,j);
          if(I == J)
            shiftedMat[i + j*localDimM] -= omega;

        } 
      else
#endif
        for(auto k = 0ul; k < N; k++) shiftedMat[k*(N+1)] -= omega;

      /*
      if( isDist )
        prettyPrintSmart(std::cout,"RHS " + std::to_string(iOmega),distSOL,
          N,nRHS,N);
      else
        prettyPrintSmart(std::cout,"RHS " + std::to_string(iOmega),SOL,
          N,nRHS,N);
     // prettyPrintSmart(std::cout,"FM " + std::to_string(iOmega),shiftedMat,
     //   N,N,N);
       */



      // Solve the linear System
      ProgramTimer::tick("Solve Linear System");
#ifdef CQ_ENABLE_MPI
      if( isDist )
        LinSolve(N,nRHS,shiftedMat,1,1,DescA,distSOL,1,1,DescB);
      else
#endif
      { 
        int64_t* IPIV = CQMemManager::get().malloc<int64_t>(N);
        lapack::gesv(N,nRHS,shiftedMat,N,IPIV,SOL,N);
        CQMemManager::get().free(IPIV);
      }
      
      ProgramTimer::tock("Solve Linear System");


      // Gather solution to root process
#ifdef CQ_ENABLE_MPI
      if( isDist )
        fullMatGrid_->gather(N,nRHS,SOL,N,distSOL,localDimM,0,0);
#endif

    //prettyPrintSmart(std::cout,"SOL " + std::to_string(iOmega),SOL,
    //  N,nRHS,N);
      
      ProgramTimer::tock("Omega");

    }


    if( shiftedMat ) CQMemManager::get().free(shiftedMat);
    if( fMatUse and (fMatUse != fullMatrix_) ) CQMemManager::get().free(fMatUse);
    if( isDist and distSOL ) CQMemManager::get().free(distSOL);


    ProgramTimer::tock("Full FDR");

  }

  template <typename T>
  template <typename U>
  void ResponseTBase<T>::runIterFDR(FDResponseResults<T,U> &results,
      std::function< void(size_t,U,SolverVectors<U>&,SolverVectors<U>&) > &preCond ) {

    bool isRoot = MPIRank(comm_) == 0;
    bool isDist = this->genSettings.isDist();

    ProgramTimer::tick("Iter FDR");

    typename GMRES<U>::LinearTrans_t lt = [&](size_t nVec, SolverVectors<U> &V, SolverVectors<U> &AV) {

      iterLinearTrans(nVec,V,AV);

    };


    typename GMRES<U>::Shift_t pc = 
      bool(preCond) ? preCond :
      [&](size_t nVec, U shift, SolverVectors<U> &V, SolverVectors<U> &AV) {

      if( not this->fullMatrix_ ) CErr();

      AV.set_data(0, nVec, V, 0);

    };

    MPI_Comm gmresComm = (isDist or not genSettings.formFullMat) 
      ? comm_ : rcomm_;
    
    GMRES<U> gmres(gmresComm,nSingleDim_,
      genSettings.maxIter,genSettings.convCrit,lt,pc);


    // Set the RHS and shifts
    gmres.setRHS(fdrSettings.nRHS,results.RHS,this->nSingleDim_);
    gmres.setShifts(results.shifts.size(),&results.shifts[0]);

    gmres.rhsBS   = fdrSettings.nRHS;
    gmres.shiftBS = results.shifts.size();

    gmres.run();

    if( isRoot )
      std::copy_n(tryGetRawVectorsPointer(*gmres.getSol()),
                  fdrSettings.nRHS * results.shifts.size() * nSingleDim_,
                  results.SOL);

    ProgramTimer::tock("Iter FDR");

  }

};

