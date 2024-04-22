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
  void GMRES<_F>::runBatch(size_t nRHS, size_t nShift,
      std::shared_ptr<SolverVectors < _F>> RHS, _F *shifts,
      std::shared_ptr<SolverVectors < _F>> SOL, double *RHSNorm ) {

    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;

    auto topGMRES = tick();

    // Do the standard stuff...
    IterLinearSolver<_F>::runBatch(nRHS,nShift,RHS,shifts,SOL,RHSNorm);


    // Zero out the scratch allocations
//    if( isRoot ) {

      std::fill_n(J_, 2 * this->mSS_ * nRHS * nShift         , 0.);
      U_->clear();
      std::fill_n(R_, this->mSS_ * this->mSS_ * nRHS * nShift, 0.);
      std::fill_n(W_, (this->mSS_ + 1) * nRHS * nShift       , 0.);

//    }






    // Construct the initial Householder reflectors in place
      
//    if( isRoot ) {

      // Copy over residuals to HHR
      HHR_->set_data(0, nRHS * nShift, *this->RES_, 0);

      for(auto iDo = 0ul; iDo < nShift * nRHS; iDo++) {

        std::shared_ptr<SolverVectors < _F>> curHHR = std::make_shared<SolverVectorsView < _F>>(*HHR_, iDo);

        _F beta = curHHR->get(0, 0);
        if(std::abs(beta) < std::numeric_limits<double>::epsilon())
          beta = this->resNorm_.back()[iDo];
        else
          beta *= this->resNorm_.back()[iDo] / std::abs(beta);

        curHHR->set(0, 0, curHHR->get(0, 0) + beta);

        // Normalize the HHR column
        double norm = curHHR->norm2F(0, 1);
        curHHR->scale(1/norm, 0, 1);

        // Copy the HHR to the first column of U
        U_->set_data(iDo*this->mSS_, 1, *curHHR, 0);

        // Apply HHR projection to residual
        W_[iDo * (this->mSS_ + 1)] = -beta;

        //std::cout << "W0 " << W_[iDo * (this->mSS_ + 1)] << "\n";

      }

//    }



    // Update tracking / counting
    std::vector<size_t> mDim( nRHS * nShift, 0 );
    std::vector<bool>   solConv( nRHS * nShift, false );

    // AX Scratch
    std::shared_ptr<SolverVectors<_F>> VContract  = nullptr;
    std::shared_ptr<SolverVectors<_F>> AVContract = nullptr;
    if( nRHS * nShift > 1 ) {

      VContract = this->vecGen_(nRHS*nShift);
      AVContract = this->vecGen_(nRHS*nShift);

    }


    // Start the micro iterations
    size_t maxMicroIter = this->mSS_;

    if( isRoot ) std::cout << "    * Starting GMRES iterations\n\n";
    size_t iMicro;


    for(iMicro = 0; iMicro < maxMicroIter; iMicro++) {

      ProgramTimer::tick("Lin Solve Iter");

      auto topMicro = tick();

      for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) mDim[iDo]++;

      if( isConverged ) break;

      size_t nConv(0), nNotConv(0);
//      if( isRoot ) {
        nNotConv = std::count(solConv.begin(),solConv.end(),false);
        nConv = solConv.size() - nNotConv;
//      }


//      if( isRoot ) // Only root process
      for(auto iDo = 0ul; iDo < nRHS * nShift; iDo++) {

        if( solConv[iDo] ) continue;

        std::shared_ptr<SolverVectors < _F>> curV = std::make_shared<SolverVectorsView < _F>>(*this->V_, iDo);
        std::shared_ptr<SolverVectors < _F>> curHHR = std::make_shared<SolverVectorsView < _F>>(*HHR_, iDo);

        curV->set_data(0, 1, *curHHR, 0);
        curV->scale(-2. * SmartConj(curHHR->get(iMicro, 0)), 0, 1);
        curV->set(iMicro, 0, curV->get(iMicro, 0) + 1.);

        if( iMicro > 0 )
        for(int k = iMicro-1; k >= 0; k--) {

          std::shared_ptr<SolverVectors < _F>> curU = std::make_shared<SolverVectorsView < _F>>(*U_, k+iDo*this->mSS_);

          _F inner = 0;
          curU->dot_product(0, *curV, 0, 1, 1, &inner, 1);
          curV->axpy(0, 1, -2.*inner, *curU, 0);

        }

        // Explicitly normalize V column
        curV->scale(1.0 / curV->norm2F(0, 1), 0, 1);

      }

      bool CopyBuffer = VContract and (nConv != 0);// isRoot and VContract and (nConv != 0);

      size_t nContract = nShift * nRHS;
      std::shared_ptr<SolverVectors<_F>> VSend  = this->V_;
      std::shared_ptr<SolverVectors<_F>> AVRecv = this->AV_;
      if( CopyBuffer ) {

        VSend  = VContract;
        AVRecv = AVContract;

        nContract = 0;
        for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) {

          VContract->set_data(nContract++, 1, *this->V_, iDo);

        }


      }

      if( MPISize(this->comm_) > 1 ) MPIBCast(nContract,0,this->comm_);

      // Form (A-sB)X product and precondition product
      auto topLT = tick();

      this->linearTrans_(nContract, *VSend, *AVRecv);

      double durLT = tock(topLT);
      MPI_Barrier(this->comm_);

      if( CopyBuffer ) {

        size_t iContract = 0;
        for(auto iDo = 0; iDo < nRHS * nShift; iDo++)
        if( not solConv[iDo] ) {

          this->AV_->set_data(iDo, 1, *AVRecv, iContract++);

        }

      }

//      if( isRoot )
      for(auto iS = 0; iS < nShift; iS++) {
        std::shared_ptr<SolverVectors < _F>> curV  = std::make_shared<SolverVectorsView < _F>>(*this->V_, iS*nRHS);
        std::shared_ptr<SolverVectors < _F>> curAV = std::make_shared<SolverVectorsView < _F>>(*this->AV_, iS*nRHS);
        this->shiftVec_     (nRHS, -shifts[iS], *curV, *curAV);
        this->preCondWShift_(nRHS, shifts[iS], *curAV, *curAV);
      }


      // Copy AV -> V
//      if( isRoot )
        this->V_->set_data(0, nRHS * nShift, *this->AV_, 0);

      
//      if( isRoot )
        this->resNorm_.emplace_back(this->resNorm_.back());

//      if( isRoot )
      for(auto iDo = 0ul; iDo < nRHS * nShift; iDo++) {

        if( solConv[iDo] ) continue;

        std::shared_ptr<SolverVectors < _F>> curV   = std::make_shared<SolverVectorsView < _F>>(*this->V_, iDo);
        std::shared_ptr<SolverVectors < _F>> curHHR = std::make_shared<SolverVectorsView < _F>>(*HHR_, iDo);
        _F * curJ   = J_       + iDo*2*this->mSS_;
        _F * curW   = W_       + iDo*(this->mSS_ + 1);
        _F * curR   = R_       + iDo*this->mSS_*this->mSS_;

        for(auto k = 0; k <= iMicro; k++) {

          std::shared_ptr<SolverVectors < _F>> curU = std::make_shared<SolverVectorsView < _F>>(*U_, k+iDo*this->mSS_);

          _F inner = 0;
          curU->dot_product(0, *curV, 0, 1, 1, &inner, 1);
          curV->axpy(0, 1, -2.*inner, *curU, 0);

        }


        std::shared_ptr<SolverVectors < _F>> curU = std::make_shared<SolverVectorsView < _F>>(*U_, iMicro + 1 + iDo*this->mSS_);

        // Determine next projector
        if( iMicro < maxMicroIter - 1 ) {

          curHHR->set_data(0, 1, *curV, 0);
          for (size_t i = 0; i < iMicro+1; i++)
            curHHR->set(i, 0, 0.0);

          _F alpha = curHHR->norm2F(0, 1);
          if( std::abs(alpha) > 1e-10 ) {

            if (std::abs(curV->get(iMicro+1, 0)) > 0)
            alpha *= curV->get(iMicro+1, 0) / std::abs(curV->get(iMicro+1, 0));

            curHHR->set(iMicro+1, 0, curHHR->get(iMicro+1, 0) + alpha);
            double norm = curHHR->norm2F(0, 1);
            curHHR->scale(1/norm, 0, 1);

            curU->set_data(0, 1, *curHHR, 0);

            for (size_t i = 0; i < this->N_ - iMicro - 2; i++)
              curV->set(i + iMicro + 2, 0, 0.0);
            curV->set(iMicro+1, 0, -alpha);

          }

        }


        // Apply Given's rotation
        if( iMicro > 0 ) 
        for(auto j = 0ul; j < iMicro; j++) {

          auto tmp = curV->get(j, 0);
          
          curV->set(j, 0,
            SmartConj(curJ[2*j    ]) * tmp + 
            SmartConj(curJ[2*j + 1]) * curV->get(j+1, 0));

          curV->set(j + 1, 0, curJ[2*j] * curV->get(j+1, 0) - curJ[2*j + 1] * tmp);

        }

        if( iMicro < maxMicroIter - 1 ) {

          double rho = std::sqrt(std::real(SmartConj(curV->get(iMicro, 0))*curV->get(iMicro, 0)
              + SmartConj(curV->get(iMicro+1, 0))*curV->get(iMicro+1, 0)));

          curJ[2*iMicro]   = curV->get(iMicro, 0)   / rho;
          curJ[2*iMicro+1] = curV->get(iMicro+1, 0) / rho;

          curW[iMicro+1]   = -curJ[2*iMicro+1]         * curW[iMicro];
          curW[iMicro]     = SmartConj(curJ[2*iMicro]) * curW[iMicro];

          curV->set(iMicro, 0, rho);
          curV->set(iMicro + 1, 0, 0.);

        }

        for (size_t i = 0; i < this->mSS_; i++)
          curR[i + iMicro*this->mSS_] = curV->get(i, 0);

        double normR = std::abs(curW[iMicro+1]);
        this->resNorm_.back()[iDo] = normR;


      } // loop over vectors


      double durMicro = tock(topMicro);



      if( isRoot ) {
        std::cout << "      GMRESIter " << std::setw(5) << iMicro + 1;
        if( nRHS * nShift > 1 ) {

          std::cout << "  :  ";
          std::cout << "  IHAVE   = " << nConv ;
          std::cout << "  NUPDATE = " << nNotConv;
          std::cout << "  DURATION = " << std::scientific 
            << durMicro << " s ( " << std::fixed 
            << durLT * 100. / durMicro << "% LT )\n";

        }
      }


      if( isRoot )
      for(auto iDo = 0; iDo < nShift*nRHS; iDo++) {

        double normR = this->resNorm_.back()[iDo];
        if(nRHS * nShift > 1) 
          std::cout << "        iDo = " << std::setw(6) << iDo << "  ";

        std::cout << std::scientific << std::setprecision(8);
        std::cout << "  ResNorm = " << normR;
        std::cout << "  RelResNorm = " << normR / RHSNorm[iDo % nRHS]; 
        std::cout << "\n";

      }
      if( isRoot and (nRHS * nShift > 1) ) std::cout << "\n";


      // Evaluate convergence
//      if( isRoot ) {

        auto OldSolConv(solConv);
        size_t nNewConv(0);
        for(auto iDo = 0; iDo < nRHS*nShift; iDo++) {
          solConv[iDo] = 
            (this->resNorm_.back()[iDo] / RHSNorm[iDo % nRHS]) < 
            this->convCrit_;

          if( nRHS * nShift > 1 and solConv[iDo] and not OldSolConv[iDo] ) {

            nNewConv++;
            std::cout << "          "
              <<  "*** CEASING TO UPDATE IDO = " << std::setw(4) << iDo 
              << " ***\n";

          }
        }

        if(nNewConv) std::cout << "\n";



        isConverged = std::all_of(solConv.begin(),solConv.end(),
            [&](bool x){ return x; });

//      }

      // Broadcast the convergence result to all the mpi processes
      if(MPISize(this->comm_) > 1) MPIBCast(isConverged,0,this->comm_);

      ProgramTimer::tock("Lin Solve Iter");

    } // Micro iterations


    double durGMRES = tock(topGMRES);

    // Cleanup memory
//    if( VContract )  CQMemManager::get().free(VContract);
//    if( AVContract ) CQMemManager::get().free(AVContract);

    if( isRoot )
      std::cout << "\n    * GMRES Converged in " << iMicro  
        << " Iterations (" << durGMRES << " s) \n\n";


    // Reconstruct the solution for the batch
          
//    if( isRoot )
    for(auto iDo = 0; iDo < nShift*nRHS; iDo++) {

      _F *curR   = R_   + iDo * this->mSS_ * this->mSS_;
      _F *curW   = W_   + iDo * (this->mSS_ + 1);
      std::shared_ptr<SolverVectors < _F>> curU = std::make_shared<SolverVectorsView < _F>>(*U_, iDo * this->mSS_);
      std::shared_ptr<SolverVectors < _F>> curHHR = std::make_shared<SolverVectorsView < _F>>(*HHR_, iDo);

      int nMicro = mDim[iDo];

      // Linear Solve
      int64_t* IPIV = CQMemManager::get().malloc<int64_t>(nMicro);
      lapack::gesv(nMicro,1,curR,this->mSS_,IPIV,curW,this->mSS_+1);
      CQMemManager::get().free(IPIV);

      curHHR->set_data(0, 1, *curU, nMicro-1);

      _F fact = curW[nMicro-1] * SmartConj(curU->get(nMicro-1, nMicro-1));
      curHHR->scale(-2.*fact, 0, 1);
      curHHR->set(nMicro-1, 0, curHHR->get(nMicro-1, 0) + curW[nMicro-1]);

      if( nMicro > 1 )
      for(int k = nMicro - 2; k >= 0; k--) {

        curHHR->set(k, 0, curHHR->get(k, 0) + curW[k]);
        _F inner = 0;
        curU->dot_product(k, *curHHR, 0, 1, 1, &inner, 1);
        curHHR->axpy(0, 1, -2.*inner, *curU, k);

      }

      // FIXME: assumes 0 guess
      std::shared_ptr<SolverVectors < _F>> curSOL = std::make_shared<SolverVectorsView < _F>>(*SOL, iDo);
      curSOL->set_data(0, 1, *curHHR, 0);

    }


  }




};

