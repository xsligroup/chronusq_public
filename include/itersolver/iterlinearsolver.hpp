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

#include <itersolver.hpp>

#include <util/timer.hpp>

namespace ChronusQ {






  template <typename _F>
  template <typename T >
  void IterLinearSolver<_F>::setRHS(size_t nRHS, T* RHS, size_t LDRHS) {

    // Copy RHS
    nRHS_ = nRHS;
    RHS_ = this->vecGen_(nRHS);

    // NO MPI
//    ROOT_ONLY(this->comm_);
    
    // TODO: try to make it work for other solver vectors
    if( MPIRank(this->comm_) == 0 ) {
      tryDowncastReferenceTo<RawVectors<_F>>(*RHS_,
          [&](auto& RHS_Ref, size_t RHS_shift) {
            SetMat('N', this->N_, nRHS, _F(1.), RHS, LDRHS,
                RHS_Ref.getPtr(RHS_shift), this->N_);
          }
      );
    }

    // Compute norms of RHS
    for(auto iRHS = 0; iRHS < nRHS; iRHS++)
      rhsNorm_.emplace_back(
          RHS_->norm2F(iRHS, 1)
      );

    std::cout << "\n  * IterLinearSolver has recieved " << nRHS 
              << " Right Hand Sides with norms:\n";

    for(auto iRHS = 0; iRHS < nRHS; iRHS++)
      std::cout << "    | RHS(" << std::setw(3) << iRHS << ") | = "
                << std::scientific << std::setprecision(8)
                << rhsNorm_[iRHS] << "\n";

  };

  template <typename _F>
  template <typename T >
  void IterLinearSolver<_F>::setShifts(size_t nShift, T* shifts) {

    // NO MPI
//    ROOT_ONLY(this->comm_);

    std::cout << "\n  * IterLinearSolver has recieved " << nShift
              << "  shifts:\n";

    std::copy_n(shifts, nShift, std::back_inserter(shifts_) );
    for(auto iShift = 0; iShift < nShift; iShift++)
      std::cout << "    Shift(" << std::setw(4) << iShift << ") = "
                << std::scientific << std::setprecision(8)
                << shifts_[iShift] << "\n";

  };


  template <typename _F>
  void IterLinearSolver<_F>::run() {

    bool isRoot = MPIRank(this->comm_) == 0;

    size_t nOmega = shifts_.size();
    rhsBS   = std::min(nRHS_ ,rhsBS);
    shiftBS = std::min(nOmega,shiftBS);
    this->mSS_ = std::min(this->mSS_,this->N_);

    alloc(); // Allocate Scratch space

    if( isRoot ) {
      std::cout << "\n  * IterLinearSolver will solve " << nRHS_*nOmega
                << " linear systems consisting of " 
                << nOmega << " linear shifts and " << nRHS_ << " RHS\n";

      std::cout << std::left;
      std::cout << std::setw(30) << "    * RHS Batch Size   = " << rhsBS 
        << "\n";
      std::cout << std::setw(30) << "    * Shift Batch Size = " << shiftBS 
        << "\n";
      std::cout << std::setw(30) << "    * Maximum Subspace = " << this->mSS_ 
        << "\n";

      std::cout << "\n\n\n";
    }
    
    if( MPISize(this->comm_) > 1 ) {

      MPIBCast(nOmega,0,this->comm_);
      MPIBCast(shiftBS,0,this->comm_);
      MPIBCast(nRHS_,0,this->comm_);
      MPIBCast(rhsBS,0,this->comm_);
      MPIBCast(this->N_,0,this->comm_);
      MPIBCast(this->mSS_,0,this->comm_);

    }

    // Shift batch loop
    for(auto iOmega = 0, iBatch = 0; iOmega < nOmega; iOmega += shiftBS) {

      ProgramTimer::tick("Omega");

      // Shift batch size information 
      size_t nOmegaDo = std::min(shiftBS, nOmega - iOmega);
      std::vector<_F> shiftBatch(nOmegaDo);

//      if( isRoot )
        std::copy_n(&shifts_[iOmega], nOmegaDo, shiftBatch.begin());


    // RHS batch loop
    for(auto iRHS = 0; iRHS < nRHS_; iRHS += rhsBS, iBatch++) {

      // RHS batch size information
      size_t nRHSDo   = std::min(rhsBS, nRHS_ - iRHS);
      std::shared_ptr<SolverVectors<_F>> RHSBatch =
          RHS_ == nullptr ? nullptr : std::make_shared<SolverVectorsView<_F>>(*RHS_, iRHS);

      // Pointer to solution
      std::shared_ptr<SolverVectors<_F>> SOLBatch =
          SOL_ == nullptr ? nullptr : std::make_shared<SolverVectorsView<_F>>(*SOL_, iRHS + iOmega * nRHS_);


      if( isRoot ) {
        std::cout << "  * IterLinearSolver Starting Batch " 
          << std::setw(6) << std::right << iBatch << "\n";
        std::cout << "    --------------------------------" << "------\n";

        std::cout << "    * IShift = ";
        if(nOmegaDo > 1)
          std::cout << "[" << iOmega << ", " << iOmega + nOmegaDo - 1 << "]";
        else
          std::cout << iOmega;
        std::cout << "\n";

        std::cout << "    * IRHS   = ";
        if(nRHSDo > 1)
          std::cout << "[" << iRHS << ", " << iRHS + nRHSDo - 1 << "]";
        else
          std::cout << iRHS;
        std::cout << "\n";

        std::cout << "\n";
      }


      runBatch(nRHSDo,nOmegaDo,RHSBatch,&shiftBatch[0],SOLBatch,
          &rhsNorm_[iRHS]);


      if( isRoot ) std::cout << "\n\n\n\n";

    } // end RHS batch loop

    ProgramTimer::tock("Omega");

    } // end Shift batch loop

  };





  template <typename _F>
  void IterLinearSolver<_F>::runBatch(size_t nRHS, size_t nShift,
                                      std::shared_ptr<SolverVectors<_F>> RHS, _F *shifts,
                                      std::shared_ptr<SolverVectors<_F>> SOL, double *RHSNorm) {

    // No MPI
//    ROOT_ONLY(this->comm_);

    // Zero out guess
    V_->clear();

    // Initial residual is RHS for zero guess 
    // FIXME: not for general guess
    for(auto iOmega = 0; iOmega < nShift; iOmega++)
      RES_->set_data(iOmega*nRHS, nRHS, *RHS, 0);

    // Precondition the residuals
    for(auto iOmega = 0; iOmega < nShift; iOmega++) {
      SolverVectorsView<_F> RES_view(*RES_, iOmega * nRHS);
      this->preCondWShift_(nRHS,shifts[iOmega], RES_view, RES_view);
    }


    // Clear out and compute the preconditioned residual norms
    resNorm_.clear();
    resNorm_.emplace_back(nShift * nRHS,0.);
    for(auto iDo = 0; iDo < nShift * nRHS; iDo++)
      resNorm_.back()[iDo] = RES_->norm2F(iDo, 1);

    //prettyPrintSmart(std::cout,"R0",RES_,this->N_,1,this->N_);

    // Output the residual norms
    std::cout << "    * Initial Residual Norms:\n";
    for(auto iDo = 0; iDo < nShift * nRHS; iDo++)
      std::cout << "      | RES( " << std::setw(4) << iDo << ") | = " 
        << std::scientific << std::setprecision(8) 
        << resNorm_.back()[iDo] << "\n";
    std::cout << "\n\n";
    
  }










};


