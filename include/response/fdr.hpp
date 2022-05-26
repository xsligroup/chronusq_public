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
      std::tie(localDimM,localDimN) = fullMatGrid_->getLocalDims(N,N);
#endif
    


    // Allocate space for the shifted matrix
    U* shiftedMat = nullptr;

    if( localDimM and localDimN )
      shiftedMat = memManager_.template malloc<U>(localDimM*localDimN);


    U* distSOL = nullptr;
    CB_INT localRHS;
#ifdef CQ_ENABLE_MPI
    CXXBLACS::ScaLAPACK_Desc_t DescA, DescB;
    if( isDist ) {

      // Allocate space for distributed RHS / Solution
      std::tie(localDimM, localRHS) = fullMatGrid_->getLocalDims(N,nRHS);
      if( localDimM and localRHS)
        distSOL = memManager_.template malloc<U>(localDimM*localRHS);


      DescB = fullMatGrid_->descInit(N,nRHS,0,0,localDimM);
      DescA = fullMatGrid_->descInit(N,N,0,0,localDimM);

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
        fullMatGrid_->Scatter(N,nRHS,SOL,N,distSOL,localDimM,0,0);
#endif

      // Copy over the full matrix and shift
      if( shiftedMat ) std::copy_n(fMatUse,localDimM*localDimN,shiftedMat);

#ifdef CQ_ENABLE_MPI
      if( isDist and shiftedMat )
        for(auto j = 0ul; j < localDimN; j++) 
        for(auto i = 0ul; i < localDimM; i++) {

          CB_INT I,J;
          std::tie(I,J) = fullMatGrid_->globalFromLocal(i,j);
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
        int64_t* IPIV = memManager_.malloc<int64_t>(N);
        lapack::gesv(N,nRHS,shiftedMat,N,IPIV,SOL,N);
        memManager_.free(IPIV);
      }
      
      ProgramTimer::tock("Solve Linear System");


      // Gather solution to root process
#ifdef CQ_ENABLE_MPI
      if( isDist )
        fullMatGrid_->Gather(N,nRHS,SOL,N,distSOL,localDimM,0,0);
#endif

    //prettyPrintSmart(std::cout,"SOL " + std::to_string(iOmega),SOL,
    //  N,nRHS,N);
      
      ProgramTimer::tock("Omega");

    }


    if( shiftedMat ) memManager_.free(shiftedMat);
    if( fMatUse and (fMatUse != fullMatrix_) ) memManager_.free(fMatUse);
    if( isDist and distSOL ) memManager_.free(distSOL);


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
    
    GMRES<U> gmres(gmresComm,this->memManager_,nSingleDim_,
      genSettings.maxIter,genSettings.convCrit,lt,pc);


    // Set the RHS and shifts
    gmres.setRHS(fdrSettings.nRHS,results.RHS,this->nSingleDim_);
    gmres.setShifts(results.shifts.size(),&results.shifts[0]);

    gmres.rhsBS   = fdrSettings.nRHS;
    gmres.shiftBS = results.shifts.size();

    gmres.run();

    if( isRoot )
      std::copy_n(gmres.getSol()->getPtr(),
                  fdrSettings.nRHS * results.shifts.size() * nSingleDim_,
                  results.SOL);

    ProgramTimer::tock("Iter FDR");

  }

};

