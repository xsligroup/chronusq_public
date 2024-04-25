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
#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/eig.hpp>
#include <cqlinalg.hpp>
#include <cerr.hpp>

namespace ChronusQ {


  template <typename _F>
  bool GPLHR<_F>::runMicro() {


    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;


    if( isRoot ) {
      std::cout << "\n\n\n";
      std::cout << "  * GPLHR Settings\n"
                << "    * Sigma                        = " << std::real(sigma) << "\n"
                << "    * Min                          = " << std::real(hardLim) << "\n"
                << "    * Kyrlov-Arnoldi Parameter (M) = " << this->m
                << std::endl;


      std::cout << "\n\n" << std::endl;
    }
    
    const size_t N = this->N_;
    const size_t MSS = this->mSS_;
    const size_t nR = this->nRoots_;


    const size_t MSS2 = MSS * MSS;
    const size_t NMSS = N * MSS;
    const size_t NNR  = N * nR;
    const size_t M_NR = this->m * nR;


    // FIXME Reassign sigma (TEMP)
    if( std::imag(sigma) > 1e-10 )
      CErr("GPLHR Only Supports Real Sigma for the Time Being",std::cout);
    double sigmaD   = std::real(sigma);
    double hardLimD = std::real(hardLim);
    //sigmaD = hardLimD;


    // Initialize pointers
    std::shared_ptr<SolverVectors<_F>> V = nullptr;
    std::shared_ptr<SolverVectors<_F>> W = nullptr , S = nullptr , P = nullptr ;
    std::shared_ptr<SolverVectors<_F>> AV = nullptr;
    std::shared_ptr<SolverVectors<_F>> AW = nullptr, AS = nullptr, AP = nullptr;
    std::shared_ptr<SolverVectors<_F>> Q = nullptr;
    std::shared_ptr<SolverVectors<_F>> Qp = nullptr, Q1 = nullptr, Q2 = nullptr,
       Q3_p = nullptr;


    _F *BETA = nullptr; dcomplex *ALPHA = nullptr;
    _F *RMAT = nullptr, *PHI = nullptr, *PSI = nullptr;
    _F *VSR  = nullptr, *VSL = nullptr, *MA  = nullptr, *MB = nullptr;
    std::shared_ptr<SolverVectors<_F>> VSCR = nullptr, VSCR2 = nullptr;

      // V  = [V W S1...Sm P]
//      V = CQMemManager::get().malloc<_F>(NMSS);
      V = this->vecGen_(MSS);
      W = std::make_shared<SolverVectorsView<_F>>(*V, nR);
      S = std::make_shared<SolverVectorsView<_F>>(*W, nR);
      P = std::make_shared<SolverVectorsView<_F>>(*S, this->m * nR);

      // AV = [AV AW AS1...ASm AP]
      AV = this->vecGen_(MSS);
      AW = std::make_shared<SolverVectorsView<_F>>(*AV, nR);
      AS = std::make_shared<SolverVectorsView<_F>>(*AW, nR);
      AP = std::make_shared<SolverVectorsView<_F>>(*AS, this->m * nR);


      // Q = [Q Q1 Q2 Q3]
      Q = this->vecGen_(MSS);
      Qp = std::make_shared<SolverVectorsView<_F>>(*Q, nR);

      Q1   = std::make_shared<SolverVectorsView<_F>>(*Qp);
      Q2   = std::make_shared<SolverVectorsView<_F>>(*Q1, nR);
      Q3_p = std::make_shared<SolverVectorsView<_F>>(*Q2, this->m * nR);





//      if( isRoot ) {

      // EVAL(I) = ALPHA(I) / BETA(I)
      ALPHA = CQMemManager::get().malloc<dcomplex>(MSS);
      BETA  = CQMemManager::get().malloc<_F>(MSS);
      //double *RITZ  = CQMemManager::get().malloc<double>(MSS);
      //double *nRITZ = CQMemManager::get().malloc<double>(nR);

      RMAT = CQMemManager::get().malloc<_F>(MSS2); 
      PHI  = CQMemManager::get().malloc<_F>(MSS2); 
      PSI  = CQMemManager::get().malloc<_F>(MSS2); 

      VSR = CQMemManager::get().malloc<_F>(MSS2);
      VSL = CQMemManager::get().malloc<_F>(MSS2);
      MA  = CQMemManager::get().malloc<_F>(MSS2);
      MB  = CQMemManager::get().malloc<_F>(MSS2);

//    } // ROOT only






    // Need a NNR scratch on all processes
    VSCR  = this->vecGen_(nR);
    VSCR2 = this->vecGen_(nR);








//    if( isRoot ) {

      // Initailize V as Guess, if no Guess set, init to identity
      if( Guess ) V->set_data(0, nR, *Guess, 0);
      else {
        V->clear();
        for(auto i = 0ul; i < nR; i++) V->set(i, i, 1.0);
      }


      // V <- QR(V)
      V->QR(0, nR);

//    } // ROOT only





    // Sync processes
    MPI_Barrier(this->comm_);


    std::shared_ptr<SolverVectors<_F>> VSend  = V;
    std::shared_ptr<SolverVectors<_F>> AVRecv = AV;


    // AV <- A * V
    this->linearTrans_(nR, *VSend, *AVRecv);

    // Sync processes
    MPI_Barrier(this->comm_);





    double nrmA = 0.;

//    if( isRoot ) {


      // Q <- AV - sig * V
      Q->set_data(0, nR, *AV, 0);
      Q->axpy(0, nR, -sigma, *V, 0);

      // Q <- QR(Q)
      Q->QR(0, nR);
        



      // PHI <- Q**H * AV
      // PSI <- Q**H * V
      Q->dot_product(0, *AV, 0, nR,nR,PHI,nR);
      Q->dot_product(0, *V, 0, nR,nR,PSI,nR);



      // VSR, VSL, ALPHA, BETA <- ORDQZ(PHI,PSI,sigma)
      OrdQZ2('V','V',nR,PHI,nR,PSI,nR,ALPHA,BETA,hardLimD,sigmaD,
             VSL,nR,VSR,nR);


      // Update eigs with ALPHA/BETA
      for (auto i = 0ul; i < nR; i++)
        this->eigVal_[i] = ALPHA[i] / BETA[i];

      // Print warinings if needed
      if( std::is_same<double, _F>::value ) {

        for( auto i = 0ul; i < nR; i++ )
          if( std::abs(std::imag(this->eigVal_[i])) > 1e-10 )
            std::cout << "  *** WARNING: INTERMEDIATE COMPLEX EIGENVALUE "
                      << "INCURRED ***" << std::endl;

      }

      /*
      // Adaptive sigma
      double limDiff = std::numeric_limits<double>::infinity();
      double omega0 = 0.;
      double midPt = 0.;
      bool belowTh = false;
      double maxW = -std::numeric_limits<double>::infinity();
      */


      // Get MA and MB ('Q-free' Schur form)
      getTriU(nR,PHI,nR,PSI,nR,MA,nR,MB,nR);


      // Right and Left Schur Vectors
      // VR  <- VR  * VSR
      // VL  <- VL  * VSL
      // AVR <- AVR * VSR
      V->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSR,nR,_F(0.),*VSCR, 0);
      V->set_data(0, nR, *VSCR, 0);
      
      Q->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSL,nR,_F(0.),*VSCR, 0);
      Q->set_data(0, nR, *VSCR, 0);
      
      AV->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSR,nR,_F(0.),*VSCR, 0);
      AV->set_data(0, nR, *VSCR, 0);

      // Form initial residuals in W
      // W = AV * MB - V * MA
      AV->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),MB,nR,_F(0.),*W, 0);
      V->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(-1.),MA,nR,_F(1.),*W, 0);


      // Get Residual Norms
      nrmA = AV->norm2F(0, nR);
      getResidualNorms(N,nR,*W,RelRes,this->eigVal_,nrmA);




      std::cout << "  * Initial Matrix Norm Estimate = " << nrmA << std::endl;

      std::cout << "\n\n  * Initial Eigenvalues\n\n";
      std::cout << std::setprecision(10) << std::scientific;
      for(auto iRt = 0ul; iRt < nR; iRt++) {

        std::cout << "    IRt = " << std::setw(5) << std::left << iRt;

        std::cout << std::right;
        if( std::is_same<double, _F>::value ) std::cout << "EigVal = "; 
        else                                  std::cout << "Re(EigVal) = ";


        std::cout << std::setw(20) << std::real(this->eigVal_[iRt]);

        if( std::is_same<dcomplex,_F>::value ) {
          std::cout << "    Im(EigVal) = ";
          std::cout << std::setw(20) << std::imag(this->eigVal_[iRt]);
        }


        std::cout << "    RelResNorm = ";
        std::cout << std::setw(20) << RelRes[iRt];


        std::cout << std::endl;
      }


      std::cout << "\n\n\n" << std::endl;

      std::cout << "  * Starting GPLHR Iterations\n" << std::endl;



//    } // ROOT only

    size_t iter = 0;

    // ****************************
    // ** Begin GPLHR iterations **
    // ****************************
    for( iter = 0; iter < this->maxMicroIter_; iter++ ) {

      ProgramTimer::tick("Diagonalize Iter");

//      if( isRoot ) {

        std::cout << "    GPLHRIter " << std::setw(5) << iter+1;
    
        // V, RMAT <- QR(V)
        V->QR(0, nR, RMAT, nR);

        // AV <- X : [X * RMAT = AV]
        AV->trsm(0, nR,_F(1.),RMAT,nR);

        // Q <- QR(Q)
        Q->QR(0, nR);


        // W = (I - V * V**H) * T * (I - V * V**H) * W 
        newSMatrix(nR,*V,*V,*W,RMAT,nR);

        // W <- QR(W)
        W->QR(0, nR);

//      } // ROOT only


      // Sync processes
      MPI_Barrier(this->comm_);

      // AW <- A * W
      auto LTst = tick();

      std::shared_ptr<SolverVectors<_F>> WSend  = W;
      std::shared_ptr<SolverVectors<_F>> AWRecv = AW;

      this->linearTrans_(nR, *WSend, *AWRecv);


      double LTdur = tock(LTst);

      // Sync processes
      MPI_Barrier(this->comm_);
    

      // Form S-blocks
        
      // S(0) = W, AS(0) = AW
      for( auto k = 1; k <= this->m; k++ ) {

        std::shared_ptr<SolverVectors<_F>> SSend = nullptr, ASRecv = nullptr;

        // S(k-1) / AS(k-1)
        std::shared_ptr<SolverVectors<_F>> Sprev  = std::make_shared<SolverVectorsView<_F>>(*W, (k - 1) * nR);
        std::shared_ptr<SolverVectors<_F>> ASprev = std::make_shared<SolverVectorsView<_F>>(*AW, (k - 1) * nR);

        // S(k) / AS(k)
        std::shared_ptr<SolverVectors<_F>> Scur  = std::make_shared<SolverVectorsView<_F>>(*W, k * nR);
        std::shared_ptr<SolverVectors<_F>> AScur = std::make_shared<SolverVectorsView<_F>>(*AW, k * nR);

        SSend  = Scur;
        ASRecv = AScur;



//        if( isRoot ) {
          // S(k) = AS(k-1) * MB - S(k-1) * MA
          ASprev->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.) ,MB,nR,_F(0.),*Scur, 0);
          Sprev->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(-1.),MA,nR,_F(1.),*Scur, 0);

          // S(k) = (I - V * V**H) * T * (I - V * V**H) * S(k-1) 
          newSMatrix(nR,*V,*V,*Scur,RMAT,nR);


          // Project out previous S's
          // S(k) = (I - S(l) * S(l)**H) * S(k)  for l = [0,k)
          for(auto l = 0; l < k; l++) 
            halfProj(nR, SolverVectorsView<_F>(*W, l * nR), *Scur, RMAT, nR);
       

          // S(k) = QR(S(k))
          Scur->QR(0, nR);

//        }



        // Sync processes
        MPI_Barrier(this->comm_);

        // AS(k) = A * S(k) 
        LTst = tick();

        this->linearTrans_(nR, *SSend, *ASRecv);

        LTdur += tock(LTst);

        // Sync processes
        MPI_Barrier(this->comm_);

      }




      const size_t nQp = iter ? 2 + this->m : 1 + this->m;
      const size_t nQ  = nQp + 1;

//      if( isRoot ) {

        // Conjugate direction orthogonalization
        if (iter) {
          
          // P = (I - V * V**H) * P
          // P = (I - W * W**H) * P
          // P = (I - S * S**H) * P
          halfProj2(     nR,*V,*AV,*P,*AP,RMAT,nR  ); // Project out V
          halfProj2(     nR,*W,*AW,*P,*AP,RMAT,nR  ); // Project out W
          halfProj2(M_NR,nR,*S,*AS,*P,*AP,RMAT,M_NR); // Project out S

          // P, RMAT <- QR(P)
          P->QR(0, nR, RMAT, nR);

          // AP <- X : [X * RMAT = AP]
          AP->trsm(0, nR,_F(1.),RMAT,nR);

        }

        std::shared_ptr<SolverVectors<_F>> Q3 = iter ? Q3_p : nullptr;

        // Q' = (A - sigma * I) [W, S, [P]]
        Qp->set_data(0, nQp*nR, *AW, 0);
        Qp->axpy(0, nQp * nR, -sigma, *W, 0);




        // Q1 = (I - Q * Q**H) * Q1
        // Q1 = QR(Q1)
        halfProj(nR,*Q,*Q1,RMAT,nR);
        Q1->QR(0, nR);

        // Q2 = (I - Q  * Q**H ) * Q2
        // Q2 = (I - Q1 * Q1**H) * Q2
        // Q2 = QR(Q2)
        halfProj(nR,M_NR,*Q,*Q2,RMAT,nR);
        halfProj(nR,M_NR,*Q1,*Q2,RMAT,nR);
        Q2->QR(0, M_NR);


        if( Q3 ) {
          // Q3 = (I - Q  * Q**H ) * Q3
          // Q3 = (I - Q1 * Q1**H) * Q3
          // Q3 = (I - Q2 * Q2**H) * Q3
          // Q3 = QR(Q3)
          halfProj(nR,*Q,*Q3,RMAT,nR);
          halfProj(nR,*Q1,*Q3,RMAT,nR);
          halfProj(M_NR,nR,*Q2,*Q3,RMAT,M_NR);
          Q3->QR(0, nR);
        }



        // PHI <- Q**H * AV
        // PSI <- Q**H * V
        Q->dot_product(0, *AV, 0, nQ*nR, nQ*nR, PHI, nQ*nR);
        Q->dot_product(0, *V, 0, nQ*nR, nQ*nR, PSI, nQ*nR);

        // VSR, VSL, ALPHA, BETA <- ORDQZ(PHI,PSI,sigma)
        OrdQZ2('V','V',nQ*nR,PHI,nQ*nR,PSI,nQ*nR,ALPHA,BETA,hardLimD,sigmaD,
          VSL,nQ*nR,VSR,nQ*nR);


        // VSRt is thick-restart
        _F * VSRt = VSR + nR * (nQ * nR);

        _F * VSR_V = VSR;
        _F * VSR_W = VSR_V + nR;
        _F * VSR_S = VSR_W + nR;
        _F * VSR_P = VSR_S + M_NR;

        _F * VSL_V = VSL;
        _F * VSL_W = VSL_V + nR;
        _F * VSL_S = VSL_W + nR;
        _F * VSL_P = VSL_S + M_NR;

        _F * VSRt_V = VSRt;
        _F * VSRt_W = VSRt_V + nR;
        _F * VSRt_S = VSRt_W + nR;
        _F * VSRt_P = VSRt_S + M_NR;

        // VSCR = V * VSR_V + W * VSR_W + S * VSR_S + P * VSR_P
        V->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSR_V,nQ*nR,_F(0.),*VSCR, 0);
        W->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSR_W,nQ*nR,_F(1.),*VSCR, 0);
        S->multiply_matrix(0, blas::Op::NoTrans,nR,M_NR,_F(1.),VSR_S,nQ*nR,_F(1.),*VSCR, 0);
        if( iter )
          P->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSR_P,nQ*nR,_F(1.),*VSCR, 0);


        // VSCR2 = V * VSRt_V + W * VSRt_W + S * VSRt_S + P * VSRt_P
        V->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSRt_V,nQ*nR,_F(0.),*VSCR2, 0);
        W->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSRt_W,nQ*nR,_F(1.),*VSCR2, 0);
        S->multiply_matrix(0, blas::Op::NoTrans,nR,M_NR,_F(1.),VSRt_S,nQ*nR,_F(1.),*VSCR2, 0);
        if( iter )
          P->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSRt_P,nQ*nR,_F(1.),*VSCR2, 0);

        // V = VSCR
        // P = VSCR2
        V->set_data(0, nR, *VSCR, 0);
        P->set_data(0, nR, *VSCR2, 0);


        // VSCR = AV * VSR_V + AW * VSR_W + AS * VSR_S + AP * VSR_P
        AV->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSR_V,nQ*nR,_F(0.),*VSCR, 0);
        AW->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSR_W,nQ*nR,_F(1.),*VSCR, 0);
        AS->multiply_matrix(0, blas::Op::NoTrans,nR,M_NR,_F(1.),VSR_S,nQ*nR,_F(1.),*VSCR, 0);
        if( iter )
          AP->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSR_P,nQ*nR,_F(1.),*VSCR, 0);

        // VSCR2 = AV * VSRt_V + AW * VSRt_W + AS * VSRt_S + AP * VSRt_P
        AV->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSRt_V,nQ*nR,_F(0.),*VSCR2, 0);
        AW->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSRt_W,nQ*nR,_F(1.),*VSCR2, 0);
        AS->multiply_matrix(0, blas::Op::NoTrans,nR,M_NR,_F(1.),VSRt_S,nQ*nR,_F(1.),*VSCR2, 0);
        if( iter )
          AP->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSRt_P,nQ*nR,_F(1.),*VSCR2, 0);


        // AV = VSCR
        // AP = VSCR2
        AV->set_data(0, nR, *VSCR, 0);
        AP->set_data(0, nR, *VSCR2, 0);

        // VSCR = Q * VSR_V + Q1 * VSR_W + Q2 * VSR_S + Q3 * VSR_P
        Q->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSL_V,nQ*nR,_F(0.),*VSCR, 0);
        Q1->multiply_matrix(0, blas::Op::NoTrans,nR,nR  ,_F(1.),VSL_W,nQ*nR,_F(1.),*VSCR, 0);
        Q2->multiply_matrix(0, blas::Op::NoTrans,nR,M_NR,_F(1.),VSL_S,nQ*nR,_F(1.),*VSCR, 0);
        if( Q3 )
          Q3->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSL_P,nQ*nR,_F(1.),*VSCR, 0);

        // Q = VSCR
        Q->set_data(0, nR, *VSCR, 0);

//      } // ROOT only

      // Refresh AV
      if( (iter+1) % 100 == 0 ) {

        if( isRoot ) 
          std::cout << "  * Refreshing AV at iteration " << iter+1 
            << std::endl;

        // Sync processes
        MPI_Barrier(this->comm_);

        LTst = tick();


        this->linearTrans_(nR,*VSend,*AVRecv);

        LTdur += tock(LTst);

        // Sync processes
        MPI_Barrier(this->comm_);
      }

//      if( isRoot ) {

        // Update MA, MB
        getTriU(nR,PHI,nQ*nR,PSI,nQ*nR,MA,nR,MB,nR);
        
        // Update eigenvalues
        for(auto i = 0; i < nR; i++)
          this->eigVal_[i] = ALPHA[i] / BETA[i];
        

        // W = AV * MB - V * MA
        AV->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.) ,MB,nR,_F(0.),*W, 0);
        V->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(-1.),MA,nR,_F(1.),*W, 0);
        
        /*
        // Adaptive sigma
        limDiff = std::numeric_limits<double>::infinity();
        maxW = -std::numeric_limits<double>::infinity();

        
        for(auto i = 0ul; i < nQ*nR; i++)
          RITZ[i] = std::real(ALPHA[i] / BETA[i]);

        std::sort(RITZ,RITZ + nQ*nR);
        //prettyPrintSmart(std::cout,"RITZ SORT",RITZ,MSS,1,MSS);

        for( auto i = 0ul; i < MSS; i++ ) {    
          if( RITZ[i] > hardLimD) {
            std::copy(RITZ + i, RITZ + i + nR, nRITZ);
            omega0 = RITZ[i-1];
            maxW = nRITZ[nR-1];
            break;
          }
        }
        
        belowTh = true;
        

        
        midPt = (maxW + omega0) / 2.;
        if (midPt > hardLimD and belowTh) {
          sigmaD = midPt;
        } else {
          sigmaD = hardLimD;
        }
        std::cout << "midPt: " << midPt << std::endl;
        std::cout << "New Sigma: " << sigmaD << std::endl;

        */


        // Get Residual norms
        getResidualNorms(N,nR,*W,RelRes,this->eigVal_,nrmA);
        

        // Check convergence
        isConverged = checkConv(nR,RelRes);

//      } // ROOT only


      // Bcast converged
      MPIBCast(isConverged,0,this->comm_);





      ProgramTimer::tock("Diagonalize Iter");


      if( isRoot ) {

        auto GPLHRdur = ProgramTimer::getDurationTotal<CQSecond>(
          "Diagonalize Iter").count();

        double perLT = LTdur * 100 / GPLHRdur;

        std::cout << "  DURATION = " << std::setprecision(8) << GPLHRdur 
          << " s  ( " << perLT << " % LT )" << std::endl;

        std::cout << std::setprecision(10) << std::scientific;
        for(auto iRt = 0ul; iRt < nR; iRt++) {

          std::cout << "      IRt = " << std::setw(5) << std::left << iRt;

          std::cout << std::right;
          if( std::is_same<double, _F>::value ) std::cout << "EigVal = "; 
          else                                  std::cout << "Re(EigVal) = ";


          std::cout << std::setw(20) << std::real(this->eigVal_[iRt]);

          if( std::is_same<dcomplex,_F>::value ) {
            std::cout << "    Im(EigVal) = ";
            std::cout << std::setw(20) << std::imag(this->eigVal_[iRt]);
          }


          std::cout << "    RelResNorm = ";
          std::cout << std::setw(20) << RelRes[iRt];


          std::cout << std::endl;
        }

        std::cout << "\n" << std::endl;

      } // ROOR only



      if( isConverged ) break; // Break loop on convergence


      

    } // end for


    MPI_Barrier(this->comm_); // Sync processes

    


//    if( isRoot ) {

      // Reconstruct Eigen vectors
      V->dot_product(0, *AV, 0, nR, nR, PSI, nR);
      GeneralEigen('N', 'V', nR, PSI, nR, ALPHA, VSL, nR, VSR, nR);
      V->multiply_matrix(0, blas::Op::NoTrans,nR,nR,_F(1.),VSR,nR,_F(0.),*this->VR_, 0);

//    }


    // Free Scratch space

    if(ALPHA)  CQMemManager::get().free(ALPHA);
    if(BETA)   CQMemManager::get().free(BETA);
    if(RMAT)   CQMemManager::get().free(RMAT);
    if(PHI)    CQMemManager::get().free(PHI);
    if(PSI)    CQMemManager::get().free(PSI);
    if(VSL)    CQMemManager::get().free(VSL);
    if(VSR)    CQMemManager::get().free(VSR);
    if(MA)     CQMemManager::get().free(MA);
    if(MB)     CQMemManager::get().free(MB);

    if( isRoot ) {

      std::cout << "  * ";
      if( isConverged )
        std::cout << "GPLHR Converged in " << iter+1 << " Iterations" 
                  << std::endl;
      else
        std::cout << "GPLHR Failed to Converged in " << iter+1 << " Iterations" 
                  << std::endl;

    }


    return isConverged;

  }



  /*
    Calculate relative norm for each residual vector 
    Assuming B = Identity.
    Assuming eigenvalues are finite.
  */
  template <typename _F>
  void GPLHR<_F>::getResidualNorms(size_t N, size_t nR, const SolverVectors<_F> &WMAT, double *RelRes, dcomplex *LAMBDA, double nrmA) {

    // ROOT_ONLY(this->comm_);


    for (auto i = 0; i < nR; i++) {

      double nrmI = WMAT.norm2F(i, 1);
      RelRes[i] = nrmI / (nrmA + std::abs(LAMBDA[i]));

    }

  }

  template <typename _F>
  bool GPLHR<_F>::checkConv(size_t nR, double *RelRes) { 

    const bool conv = std::none_of(RelRes,RelRes + nR,
        [&]( double x ) -> bool { 
          return x > this->convCrit_ or std::isnan(x);
        });

    return conv;
  }



  template <typename _F>
  void GPLHR<_F>::getTriU(size_t N, _F *TRIUA, size_t LDTRIUA, _F *TRIUB, 
    size_t LDTRIUB, _F *MA, size_t LDMA, _F *MB, size_t LDMB){

    // ROOT_ONLY(this->comm_);
#if 0

    // G(TRIUB) = TRIUA * CA + TRIUB * CB
    for(auto i = 0ul; i < N; i++){ 
      blas::scal(N,CB[i*(N+1)],                TRIUB + i*ldm,1);
      blas::axpy( N,CA[i*(N+1)],TRIUA + i*ldm,1,TRIUB + i*ldm,1);
    }

    // TRIUA = inv(TRIUB) * TRIUA
    blas::trsm(blas::Layout::ColMajor,blas::Side::Left,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),TRIUB,ldm,TRIUA,ldm);

    //for(auto i = 0ul; i < N; i++){ 

    //  std::copy_n(TRIUA + i*ldm, N, MA + i*N);
    //  std::copy_n(TRIUA + i*ldm, N, MB + i*N);

    //  blas::scal(N, CB[i*(N+1)],MA + i*N,1);
    //  blas::scal(N,-CA[i*(N+1)],MB + i*N,1);

    //  MB[ i*(N + 1) ] += 1.;

    //}

#else

    // Initialize CA and CB as identity FIXME: memory
    _F *CA    = CQMemManager::get().malloc<_F>(N*N);
    _F *CB    = CQMemManager::get().malloc<_F>(N*N);
    _F *Ident = CQMemManager::get().malloc<_F>(N*N);
    _F *G     = CQMemManager::get().malloc<_F>(N*N);
    std::fill_n(CA,N*N,_F(0.)); 
    std::fill_n(CB,N*N,_F(0.)); 
    std::fill_n(Ident,N*N,_F(0.)); 
    for (auto i = 0; i < N; i++) {
      CA[i + i*N]    = _F(1.0);
      CB[i + i*N]    = _F(1.0);
      Ident[i + i*N] = _F(1.0);
    }

    // Create CA and CB
    for (auto i = 0; i < N; i++) {
      if (std::abs(TRIUA[i + i*LDTRIUA]) >= std::abs(TRIUB[i + i*LDTRIUB])) {
        CA[i + i*N] = (_F(1.) - TRIUB[i + i*LDTRIUB]) / TRIUA[i + i*LDTRIUA];
      } else {
        CA[i + i*N] = _F(0.);
        CB[i + i*N] = _F(1.) / TRIUB[i + i*LDTRIUA];
      }
    }

    // Form G
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),TRIUA,LDTRIUA,CA,N,_F(0.),G,N);  
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),TRIUB,LDTRIUB,CB,N,_F(1.),G,N);
    
    // Form MA
    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),G,N,CB,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),CB,N,TRIUA,LDTRIUA,_F(0.),MA,LDMA);

    // Form MB   
    blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Upper,blas::Op::NoTrans,blas::Diag::NonUnit,
      N,N,_F(1.),G,N,CA,N);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,_F(1.),CA,N,TRIUA,LDTRIUA,_F(0.),MB,LDMB);
    blas::axpy(N*N,_F(-1.),MB,1,Ident,1);
    std::copy_n(Ident, N*N, MB);


    // Free SCR Mem
    CQMemManager::get().free(CA);
    CQMemManager::get().free(CB);
    CQMemManager::get().free(Ident);
    CQMemManager::get().free(G);
#endif


  }








  /**
   *  \brief Projects out set of vectors, V, from another set of vectors, S
   *
   *  S = (I - V * V**H) * S
   */
  template <typename _F>
  void GPLHR<_F>::halfProj(size_t nV, size_t nS, const SolverVectors<_F> &V,
                           SolverVectors<_F> &S, _F *SCR, size_t LDSCR) {

    if( LDSCR < nV ) CErr("nV MUST be >= LDSCR");

    // ROOT_ONLY(this->comm_);

    // SCR = V**H * S
    V.dot_product(0, S, 0, nV, nS, SCR, LDSCR);

    // S = S - V * SCR
    V.multiply_matrix(0, blas::Op::NoTrans, nS, nV, _F(-1.), SCR, LDSCR, _F(1.), S, 0);

  }

  /**
   *  \brief Projects out set of vectors, V, from another set of vectors, S
   *  and updates the linear transformed AS with the same projection on AV
   *
   *  S  = (I - V * V**H) * S
   *  AS = AS - V * V**H  * S
   */
  template <typename _F>
  void GPLHR<_F>::halfProj2(size_t nV, size_t nS, const SolverVectors<_F> &V, const SolverVectors<_F> &AV,
                            SolverVectors<_F> &S, SolverVectors<_F> &AS, _F *SCR, size_t LDSCR) {

    if( LDSCR < nV ) CErr("nV MUST be >= LDSCR");

    // ROOT_ONLY(this->comm_);

    // SCR = V**H * S
    V.dot_product(0, S, 0, nV, nS, SCR, LDSCR);

    // S = S - V * SCR
    V.multiply_matrix(0, blas::Op::NoTrans, nS, nV, _F(-1.), SCR, LDSCR, _F(1.), S, 0);

    // AS = AS - AV * SCR
    AV.multiply_matrix(0, blas::Op::NoTrans, nS, nV, _F(-1.), SCR, LDSCR, _F(1.), AS, 0);

  }



  /**
   *  \brief Forms a new S matrix in the Krylov-Arnoldi space
   *
   *  S = (I - V * V**H) * T * (I - U * U**H ) * S
   *
   *  where T is the preconditioner
   */
  template <typename _F>
  void GPLHR<_F>::newSMatrix(size_t nR, const SolverVectors<_F> &V, const SolverVectors<_F> &Q,
                             SolverVectors<_F> &S, _F *SCR, size_t LDSCR) {


    // ROOT_ONLY(this->comm_);

    halfProj(nR,V,S,SCR,LDSCR);
    this->preCondWShift_(nR,sigma,S,S);
    halfProj(nR,Q,S,SCR,LDSCR);

  }



  template <typename _F>
  void GPLHR<_F>::restart() { }

}; // namespace ChronusQ

