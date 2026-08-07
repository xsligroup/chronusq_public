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
#include <particleintegrals/twopints/distributedritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/gradints/incore.hpp>

// Use stupid but bullet proof incore contraction for debug
//#define _BULLET_PROOF_INCORE

//#define _REPORT_JCon
//#define _REPORT_KCon

namespace ChronusQ {

  /**
   *  \brief Perform various tensor contractions of the full ERI
   *  tensor in core. Wraps other helper functions and provides
   *  loop structure
   *
   *  Currently supports
   *    - Coulomb-type (34,12) contractions
   *    - Exchange-type (23,14) contractions
   *
   *  Works with both real and complex matricies
   *
   *  \param [in/out] list Contains information pertinent to the
   *    matricies to be contracted with. See TwoBodyContraction
   *    for details
   */
  template <typename MatsT, typename IntsT>
  void InCoreTPIContraction<MatsT, IntsT>::twoBodyContract(
      MPI_Comm comm,
      const bool,
      std::vector<TwoBodyContraction<MatsT>> &list,
      EMPerturbation&) const {
    //ROOT_ONLY(comm);
    MPI_Comm workComm = comm;
    int workRank = 0, workSize = 1;
#ifdef CQ_ENABLE_MPI
    // Handle MPI with non-distributed ERI, where non-root ranks exits prematurely
    bool requiresAllRanks = false;
    if (auto ritpi = std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_))
      requiresAllRanks = ritpi->isDistributed();
    else if (auto ritpi = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(this->ints_))
      requiresAllRanks = ritpi->isDistributed();
    workComm = requiresAllRanks ? comm : MPI_COMM_SELF;
    MPI_Comm_rank(workComm, &workRank);
    MPI_Comm_size(workComm, &workSize);
    if (!requiresAllRanks && workRank != 0) return;
#endif

    if (typeid(*this) == typeid(InCoreRITPIContraction<MatsT,IntsT>)
        or typeid(*this) == typeid(DistributedRITPIContraction<MatsT,IntsT>))
      if ( std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_) == nullptr){
        CErr("RITPIContraction expect a InCoreRITPI reference.");
      }

    if (typeid(*this) == typeid(InCoreAsymmRITPIContraction<MatsT,IntsT>)
        or typeid(*this) == typeid(DistributedAsymmRITPIContraction<MatsT,IntsT>))
      if ( std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(this->ints_) == nullptr){
        CErr("RITPIContraction expect a InCoreAsymmRITPI reference.");
      }

    if (typeid(*this) == typeid(InCoreRelERIContraction<MatsT,IntsT>)
        and typeid( *(this->ints_)) != typeid(InCoreRelERI<IntsT>))
      CErr("InCore4indexRelTPIContraction expect a InCore4indexRelTPI reference.");

    if (typeid(*this) == typeid(InCore4indexTPIContraction<MatsT,IntsT>))
      if ( std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->ints_) == nullptr)
        CErr("InCore4indexTPIContraction expect a InCore4indexTPI reference.");

    ProgramTimer::timeOp("Contraction Total", [&](){

      // Loop over matricies to contract with
      for(auto &C : list) {

        auto beginContract = tick();
        switch (C.contType) {
          case COULOMB:
          case DC_COULOMB:
          case SSSS_COULOMB:
          case GAUNT_COULOMB:
          case GAUGE_COULOMB:
            // Coulomb-type (34,12) ERI contraction
            // AX(mn) = (mn | kl) X(kl)
            JContract(workComm,C);
            if(this->printContractionTiming)
              this->printTiming("J-Contraction duration (s): ", tock(beginContract), workComm, workRank, workSize);
            break;
          case EXCHANGE:
          case DC_EXCHANGE:
          case SSSS_EXCHANGE:
          case GAUNT_EXCHANGE:
          case GAUGE_EXCHANGE:
            // Exchange-type (23,12) ERI contraction
            // AX(mn) = (mk |ln) X(kl)
            KContract(workComm,C);
            if(this->printContractionTiming)
              this->printTiming("K-Contraction duration (s): ", tock(beginContract), workComm, workRank, workSize);
            break;
          default:
            CErr("Unsupported two-body contraction type for InCoreTPIContraction.");
            break;
        }

      } // loop over matricies

    });

  } // InCore4indexTPIContraction::twoBodyContract


  /**
   *  \brief Perform a Coulomb-type (34,12) ERI contraction with
   *  a one-body operator.
   */   
  template <typename MatsT, typename IntsT>
  void InCore4indexTPIContraction<MatsT, IntsT>::JContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    ProgramTimer::tick("J Contract");

    InCore4indexTPI<IntsT> &tpi4I = *std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->ints_);
    size_t NB   = tpi4I.nBasis();
    size_t NB2  = NB*NB;
    size_t sNB  = tpi4I.snBasis();
    size_t sNB2 = sNB*sNB;

    // need to swap NB and sNB if contraction is done in aux
    if (this->contractSecond) {
      std::swap(NB,sNB);
      std::swap(NB2,sNB2);
    }
    
    const bool sameIntsTMatsT = std::is_same<IntsT,MatsT>::value;

    // for Hermitian densities, output Coulomb Matrix is same type as IntsT
    if (C.HER or sameIntsTMatsT) {

      IntsT *X  = reinterpret_cast<IntsT*>(C.X);
      IntsT *AX = reinterpret_cast<IntsT*>(C.AX);

      // Extract the real part of X if X is Hermetian and if the ints are
      // real
      const bool extractRealPartX = 
        C.HER and std::is_same<IntsT,double>::value and 
        std::is_same<MatsT,dcomplex>::value;

      // Allocate scratch if IntsT and MatsT are different
      const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;

      if( extractRealPartX ) {

      X = CQMemManager::get().malloc<IntsT>(sNB2);
      for(auto k = 0ul; k < sNB2; k++) X[k] = std::real(C.X[k]);
    }

      if( allocAXScratch ) {

        AX = CQMemManager::get().malloc<IntsT>(NB2);
        std::fill_n(AX,NB2,0.);

      }


      #ifdef _BULLET_PROOF_INCORE

    size_t NB3 = NB * NB2;
    #pragma omp parallel for
    for(auto i = 0; i < NB; ++i)
    for(auto j = 0; j < NB; ++j)
    for(auto k = 0; k < sNB; ++k)
    for(auto l = 0; l < sNB; ++l)

      AX[i + j*NB] += tpi4I.pointer()[i+j*NB+k*NB2+l*NB3] * X[l + k*NB];

      #else

    //if( std::is_same<IntsT,dcomplex>::value )
    //  blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,sNB2,IntsT(1.),tpi4I.pointer(),NB2,X,sNB2,IntsT(0.),AX,NB2);
    //else 

    // GIAO
    // J(μν) = ((μν|κλ).H * D(λκ)).H = ((νμ|λκ) * D(λκ)).H = J(νμ).H = J(μν) 
    if( std::is_same<IntsT,dcomplex>::value ){
      if(not this->isCross){
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,sNB2,IntsT(1.),tpi4I.pointer(),NB2,X,sNB2,IntsT(0.),AX,NB2);
      }else{
        if (not this->contractSecond)
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,1,NB2,sNB2,IntsT(1.),X,1,tpi4I.pointer(),NB2,IntsT(0.),AX,1);
        else
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,sNB2,IntsT(1.),tpi4I.pointer(),sNB2,X,sNB2,IntsT(0.),AX,NB2);
      }    
    }else{
      if (not this->contractSecond)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB2,1,sNB2,IntsT(1.),tpi4I.pointer(),NB2,X,sNB2,IntsT(0.),AX,NB2);
      else
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,sNB2,IntsT(1.),tpi4I.pointer(),sNB2,X,sNB2,IntsT(0.),AX,NB2);
    }
    
      // if Complex ints + Hermitian, conjugate
      if( std::is_same<IntsT,dcomplex>::value and C.HER )
        IMatCopy('R',NB,NB,IntsT(1.),AX,NB,NB);

      // If non-hermetian, transpose
      if( not C.HER )  {

        IMatCopy('T',NB,NB,IntsT(1.),AX,NB,NB);

      }

      #endif

      // Cleanup temporaries
      if( extractRealPartX ) CQMemManager::get().free(X);
      if( allocAXScratch ) {

        std::copy_n(AX,NB2,C.AX);
        CQMemManager::get().free(AX);

      }
    
    // for non-hermiatin and MatsT = dcomplex, IntsT = double
    } else {
      if (not this->contractSecond)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB2,1,sNB2,MatsT(1.),tpi4I.pointer(),NB2,C.X,sNB2,MatsT(0.),C.AX,NB2);
      else
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB2,1,sNB2,MatsT(1.),tpi4I.pointer(),sNB2,C.X,sNB2,MatsT(0.),C.AX,NB2);
    }

    ProgramTimer::tock("J Contract");

  }; // InCore4indexTPIContraction::JContract

  template <typename MatsT, typename IntsT>
  void InCore4indexTPIContraction<MatsT, IntsT>::KContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    ProgramTimer::tick("K Contract");
 
    InCore4indexTPI<IntsT> &tpi4I = *std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->ints_);
    // check to see whether the two basis sets are different
    size_t NB = tpi4I.nBasis();
    size_t NB2 = NB*NB;
    size_t NB3 = NB * NB2;

    //#ifdef _BULLET_PROOF_INCORE
    #if 0

    std::fill_n(C.AX,NB2,0.);

    #pragma omp parallel for
    for(auto i = 0; i < NB; ++i)
    for(auto j = 0; j < NB; ++j)
    for(auto k = 0; k < NB; ++k)
    for(auto l = 0; l < NB; ++l) {
      C.AX[i + j*NB] += tpi4I.pointer()[i+l*NB+k*NB2+j*NB3] * C.X[l + k*NB];
    }

    #else

    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);

    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      for (size_t i = 0; i < nAtoms; ++i) {
        size_t bfStart = mapCen2BfSt[i];
        size_t bfEnd   = (i+1 < nAtoms) ? mapCen2BfSt[i+1] : NB;
        size_t nBf = bfEnd - bfStart;
        for (size_t nu = bfStart; nu < bfEnd; ++nu) 
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nBf,1,NB2,MatsT(1.),tpi4I.pointer()+nu*NB3+bfStart,NB,C.X,NB2,MatsT(0.),C.AX+nu*NB+bfStart,NB);
      }
    } else {
      #pragma omp parallel for
      for(auto nu = 0; nu < NB; nu++) 
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,1,NB2,MatsT(1.),tpi4I.pointer()+nu*NB3,NB,C.X,NB2,MatsT(0.),C.AX+nu*NB,NB);
    }

    SetLAThreads(LAThreads);

    #endif

    ProgramTimer::tock("K Contract");

  }; // InCore4indexTPIContraction::KContract


  /**
   *  \brief Perform a Coulomb-type (34,12) ERI contraction with
   *  a one-body operator.
   */
  template <typename MatsT, typename IntsT>
  void InCoreRelERIContraction<MatsT, IntsT>::JContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    ProgramTimer::tick("J Contract");

    InCoreRelERI<IntsT> &tpi4I = *std::dynamic_pointer_cast<InCoreRelERI<IntsT>>(this->ints_);
    size_t NB = tpi4I.nBasis();
    size_t NB2 = NB * NB;
    size_t NB3 = NB * NB2;

    TPIContractionPointers<IntsT> ERI4s = tpi4I.getPointers(C.contType, C.ERI4Ind);

    memset(C.AX,0,NB2*sizeof(MatsT));

    if (not ERI4s.isRI()) {
      double *ERI4 = reinterpret_cast<double*>(ERI4s.pointers[0]);
    if( C.intTrans == TRANS_KL ) {

  #ifdef _REPORT_JCon
        auto topJ5 = tick();
  #endif

      // D(μν) = D(λκ)(μν|[κλ]^T) = D(λκ)(μν|λκ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[m + n*NB + l*NB2 + k*NB3] * C.X[l + k*NB];

      }
  #ifdef _REPORT_JCon
        auto durJ5 = tock(topJ5);
        std::cout << "J5 Contract  duration   = " << durJ5 << std::endl;
  #endif

    } else if (C.intTrans == TRANS_MN_TRANS_KL) {

  #ifdef _REPORT_JCon
        auto topJ4 = tick();
  #endif
      
      // D(μν) = D(λκ)([μν]^T|[κλ]^T) = D(λκ)(νμ|λκ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[n + m*NB + l*NB2 + k*NB3]*C.X[l + k*NB];

      }
  #ifdef _REPORT_JCon
        auto durJ4 = tock(topJ4);
        std::cout << "J4 Contract  duration   = " << durJ4 << std::endl;
  #endif
    
    } else if (C.intTrans == TRANS_MN) {

  #ifdef _REPORT_JCon
        auto topJ3 = tick();
  #endif
      
      // D(μν) = D(λκ)([μν]^T|κλ) = D(λκ)(νμ|κλ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[n + m*NB + k*NB2 + l*NB3]*C.X[l + k*NB];

      }
  #ifdef _REPORT_JCon
        auto durJ3 = tock(topJ3);
        std::cout << "J3 Contract  duration   = " << durJ3 << std::endl;
  #endif
    
    } else if( C.intTrans == TRANS_MNKL ) {

  #ifdef _REPORT_JCon
        auto topJ2 = tick();
  #endif

      // D(μν) = D(λκ)(μν|κλ)^T = D(λκ)(κλ|μν)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[k + l*NB + m*NB2 + n*NB3]*C.X[l + k*NB];

      }
  #ifdef _REPORT_JCon
        auto durJ2 = tock(topJ2);
        std::cout << "J2 Contract  duration   = " << durJ2 << std::endl;
  #endif

    } else if( C.intTrans == TRANS_NONE ) {

  #ifdef _REPORT_JCon
        auto topJ1 = tick();
  #endif

      // D(μν) = D(λκ)(μν|κλ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[m + n*NB + k*NB2 + l*NB3] * C.X[l + k*NB];

      }
  #ifdef _REPORT_JCon
        auto durJ1 = tock(topJ1);
        std::cout << "J1 Contract  duration   = " << durJ1 << std::endl;
  #endif
    }
    } // NOT RI
    else if (ERI4s.pointers.size() == 2) {
      size_t NBRI = ERI4s.aux_dims[0];
      size_t NBNBRI = NB * NBRI;
      IntsT *ERI3_1 = ERI4s.pointers[0];
      IntsT *ERI3_2 = ERI4s.pointers[1];

      auto temp = CQMemManager::get().template malloc<MatsT>(NBRI);
      memset(temp,0,NBRI*sizeof(MatsT));


      if( C.intTrans == TRANS_KL ) {

  #ifdef _REPORT_JCon
        auto topJ5 = tick();
  #endif

        // D(μν) = D(λκ)(μν|[κλ]^T) = D(λκ)(μν|λκ)
        // temp(P) = B2(P|lk)D(lk)
        // D(mn) = temp(P)B1(P|mn)
        #pragma omp parallel for
        for(auto P = 0; P < NBRI; ++P) 
        for(auto k = 0; k < NB; ++k)
        for(auto l = 0; l < NB; ++l) {

          temp[P] += ERI3_2[P + l*NBRI + k*NBNBRI] * C.X[l + k*NB];

        }

        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P] * ERI3_1[P + m*NBRI + n*NBNBRI];

        }

  #ifdef _REPORT_JCon
        auto durJ5 = tock(topJ5);
        std::cout << "J5 Contract  duration   = " << durJ5 << std::endl;
  #endif

      } else if (C.intTrans == TRANS_MN_TRANS_KL) {

  #ifdef _REPORT_JCon
        auto topJ4 = tick();
  #endif

        // D(μν) = D(λκ)([μν]^T|[κλ]^T) = D(λκ)(νμ|λκ)
        // temp(P) = B2(P|lk)D(lk)
        // D(mn) = temp(P)B1(P|nm)
        #pragma omp parallel for
        for(auto P = 0; P < NBRI; ++P) 
        for(auto k = 0; k < NB; ++k)
        for(auto l = 0; l < NB; ++l) {

          temp[P] += ERI3_2[P + l*NBRI + k*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P] * ERI3_1[P + n*NBRI + m*NBNBRI];

        }
  #ifdef _REPORT_JCon
        auto durJ4 = tock(topJ4);
        std::cout << "J4 Contract  duration   = " << durJ4 << std::endl;
  #endif

      } else if (C.intTrans == TRANS_MN) {

  #ifdef _REPORT_JCon
        auto topJ3 = tick();
  #endif

        // D(μν) = D(λκ)([μν]^T|κλ) = D(λκ)(νμ|κλ)
        // temp(P) = B2(P|kl)D(lk)
        // D(mn) = temp(P)B1(P|nm)
        #pragma omp parallel for
        for(auto P = 0; P < NBRI; ++P) 
        for(auto l = 0; l < NB; ++l)
        for(auto k = 0; k < NB; ++k) {

          temp[P] += ERI3_2[P + k*NBRI + l*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P] * ERI3_1[P + n*NBRI + m*NBNBRI];

        }
  #ifdef _REPORT_JCon
        auto durJ3 = tock(topJ3);
        std::cout << "J3 Contract  duration   = " << durJ3 << std::endl;
  #endif

      } else if( C.intTrans == TRANS_MNKL ) {

  #ifdef _REPORT_JCon
        auto topJ2 = tick();
  #endif

        // D(μν) = D(λκ)(μν|κλ)^T = D(λκ)(κλ|μν)
        // temp(P) = B1(P|kl)D(lk)
        // D(mn) = temp(P)B2(P|mn)
        #pragma omp parallel for
        for(auto P = 0; P < NBRI; ++P) 
        for(auto l = 0; l < NB; ++l)
        for(auto k = 0; k < NB; ++k) {

          temp[P] += ERI3_1[P + k*NBRI + l*NBNBRI] * C.X[l + k*NB];

        }
//        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,1,NB2,IntsT(1.),ERI3_1,NBRI,C.X,NB2,IntsT(0.),temp,NBRI);

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P] * ERI3_2[P + m*NBRI + n*NBNBRI];

        }
  #ifdef _REPORT_JCon
        auto durJ2 = tock(topJ2);
        std::cout << "J2 Contract  duration   = " << durJ2 << std::endl;
  #endif

//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB2,1,NBRI,IntsT(1.),ERI3_2,NBRI,temp,NBRI,IntsT(0.),C.AX,NB2);

      } else if( C.intTrans == TRANS_NONE ) {

  #ifdef _REPORT_JCon
        auto topJ1 = tick();
  #endif

        // D(μν) = D(λκ)(μν|κλ)
        // temp(P) = B2(P|kl)D(lk)
        // D(mn) = temp(P)B1(P|mn)
        #pragma omp parallel for
        for(auto P = 0; P < NBRI; ++P)
        for(auto l = 0; l < NB; ++l)
        for(auto k = 0; k < NB; ++k) {

          temp[P] += ERI3_2[P + k*NBRI + l*NBNBRI] * C.X[l + k*NB];

        }
//        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,1,NB2,IntsT(1.),ERI3_2,NBRI,C.X,NB2,IntsT(0.),temp,NBRI);

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P] * ERI3_1[P + m*NBRI + n*NBNBRI];

        }
  #ifdef _REPORT_JCon
        auto durJ1 = tock(topJ1);
        std::cout << "J1 Contract  duration   = " << durJ1 << std::endl;
  #endif

//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB2,1,NBRI,IntsT(1.),ERI3_1,NBRI,temp,NBRI,IntsT(0.),C.AX,NB2);

      }
      CQMemManager::get().free(temp);
     }
     ProgramTimer::tock("J Contract");

  }; // InCoreRelERIContraction::JContract


  template <typename MatsT, typename IntsT>
  void InCoreRelERIContraction<MatsT, IntsT>::KContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    ProgramTimer::tick("K Contract");

    InCoreRelERI<IntsT> &tpi4I = *std::dynamic_pointer_cast<InCoreRelERI<IntsT>>(this->ints_);
    size_t NB = tpi4I.nBasis();
    size_t NB2 = NB * NB;
    size_t NB3 = NB * NB2;

    TPIContractionPointers<IntsT> ERI4s = tpi4I.getPointers(C.contType, C.ERI4Ind);

    memset(C.AX,0,NB2*sizeof(MatsT));

    if (not ERI4s.isRI()) {
      double *ERI4 = reinterpret_cast<double*>(ERI4s.pointers[0]);
    if ( C.intTrans == TRANS_MN_TRANS_KL ) {

  #ifdef _REPORT_KCon
        auto topK4 = tick();
  #endif

      // D(μν) = D(λκ)([μλ]^T|[κν]^T) = D(λκ)(λμ|νκ)

      #pragma omp parallel for
      for (auto m = 0; m < NB; ++m)
      for (auto n = 0; n < NB; ++n)
      for (auto k = 0; k < NB; ++k)
      for (auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[l + m * NB + n * NB2 + k * NB3] * C.X[l + k * NB];

      }

  #ifdef _REPORT_KCon
        auto durK4 = tock(topK4);
        std::cout << "K4 Contract  duration   = " << durK4 << std::endl;
  #endif

    } else if( C.intTrans == TRANS_KL ) {

  #ifdef _REPORT_KCon
        auto topK3 = tick();
  #endif

      // D(μν) = D(λκ)(μλ|[κν]^T) = D(λκ)(μλ|νκ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[m + l*NB + n*NB2 + k*NB3] * C.X[l + k*NB];

      }

  #ifdef _REPORT_KCon
        auto durK3 = tock(topK3);
        std::cout << "K3 Contract  duration   = " << durK3 << std::endl;
  #endif

    } else if( C.intTrans == TRANS_MNKL ) {

  #ifdef _REPORT_KCon
        auto topK2 = tick();
  #endif

      // D(μν) = D(λκ)(μλ|κν)^T = D(λκ)(κν|μλ)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + l*NB] += ERI4[k + l*NB + m*NB2 + n*NB3] * C.X[n + k*NB];

      }

  #ifdef _REPORT_KCon
        auto durK2 = tock(topK2);
        std::cout << "K2 Contract  duration   = " << durK2 << std::endl;
  #endif

    } else if( C.intTrans == TRANS_NONE ) {

  #ifdef _REPORT_KCon
        auto topK1 = tick();
  #endif

      // D(μν) = D(λκ)(μλ|κν)
      #pragma omp parallel for
      for(auto m = 0; m < NB; ++m)
      for(auto n = 0; n < NB; ++n)
      for(auto k = 0; k < NB; ++k)
      for(auto l = 0; l < NB; ++l) {

        C.AX[m + n*NB] += ERI4[m + l*NB + k*NB2 + n*NB3] * C.X[l + k*NB];

      }

  #ifdef _REPORT_KCon
        auto durK1 = tock(topK1);
        std::cout << "K1 Contract  duration   = " << durK1 << std::endl;
  #endif
    }
    } // NOT RI
    else if (ERI4s.pointers.size() == 2) {
      size_t NBRI = ERI4s.aux_dims[0];
      size_t NBNBRI = NB * NBRI;
      IntsT *ERI3_1 = ERI4s.pointers[0];
      IntsT *ERI3_2 = ERI4s.pointers[1];

      auto temp = CQMemManager::get().template malloc<MatsT>(NBRI*NB2);
      memset(temp,0,NBRI*NB2*sizeof(MatsT));

      if ( C.intTrans == TRANS_MN_TRANS_KL ) {

  #ifdef _REPORT_KCon
        auto topK4 = tick();
  #endif

        // D(μν) = D(λκ)([μλ]^T|[κν]^T) = D(λκ)(λμ|νκ)
        // temp(P|nl) = B2(P|nk)D(lk)
        // D(mn) = temp(P|nl)B1(P|lm)
        #pragma omp parallel for
        for(auto l = 0; l < NB; ++l)
        for(auto n = 0; n < NB; ++n)
        for(auto P = 0; P < NBRI; ++P) 
        for(auto k = 0; k < NB; ++k) {

          temp[P + n*NBRI + l*NBNBRI] += ERI3_2[P + n*NBRI + k*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto l = 0; l < NB; ++l)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P + n*NBRI + l*NBNBRI] * ERI3_2[P + l*NBRI + m*NBNBRI];

        }
  #ifdef _REPORT_KCon
        auto durK4 = tock(topK4);
        std::cout << "K4 Contract  duration   = " << durK4 << std::endl;
  #endif

      } else if( C.intTrans == TRANS_KL ) {

  #ifdef _REPORT_KCon
        auto topK3 = tick();
  #endif

        // D(μν) = D(λκ)(μλ|[κν]^T) = D(λκ)(μλ|νκ)
        // temp(P|nl) = B2(P|nk)D(lk)
        // D(mn) = temp(P|nl)B1(P|ml)
        #pragma omp parallel for
        for(auto l = 0; l < NB; ++l)
        for(auto n = 0; n < NB; ++n)
        for(auto P = 0; P < NBRI; ++P) 
        for(auto k = 0; k < NB; ++k) {

          temp[P + n*NBRI + l*NBNBRI] += ERI3_2[P + n*NBRI + k*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto l = 0; l < NB; ++l)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P + n*NBRI + l*NBNBRI] * ERI3_2[P + m*NBRI + l*NBNBRI];

  #ifdef _REPORT_KCon
        auto durK3 = tock(topK3);
        std::cout << "K3 Contract  duration   = " << durK3 << std::endl;
  #endif

        }
      } else if( C.intTrans == TRANS_MNKL ) {

  #ifdef _REPORT_KCon
        auto topK2 = tick();
  #endif

        // D(μν) = D(λκ)(μλ|κν)^T = D(λκ)(κν|μλ)
        // temp(P|ln) = B1(P|kn)D(lk)
        // D(mn) = temp(P|ln)B2(P|ml)
        #pragma omp parallel for
        for(auto l = 0; l < NB; ++l)
        for(auto n = 0; n < NB; ++n)
        for(auto P = 0; P < NBRI; ++P) 
        for(auto k = 0; k < NB; ++k) {

          temp[P + l*NBRI + n*NBNBRI] += ERI3_1[P + k*NBRI + n*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto l = 0; l < NB; ++l)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P + l*NBRI + n*NBNBRI] * ERI3_2[P + m*NBRI + l*NBNBRI];

        }
  #ifdef _REPORT_KCon
        auto durK2 = tock(topK2);
        std::cout << "K2 Contract  duration   = " << durK2 << std::endl;
  #endif

//        size_t LAThreads = GetLAThreads();
//        SetLAThreads(1);
//
//        #pragma omp parallel for
//        for(auto n = 0; n < NB; ++n)
//          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBRI,NB,NB,MatsT(1.),ERI3_1+n*NBNBRI,NBRI,C.X,NB,MatsT(0.),temp+n*NBNBRI,NBRI);
//
//        SetLAThreads(LAThreads);
//
//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB,NB,NBNBRI,MatsT(1.),ERI3_2,NBNBRI,temp,NBNBRI,MatsT(0.),C.AX,NB);

      } else if( C.intTrans == TRANS_NONE ) {

  #ifdef _REPORT_KCon
        auto topK1 = tick();
  #endif

        // D(μν) = D(λκ)(μλ|κν)
        // temp(P|ln) = B2(P|kn)D(lk)
        // D(mn) = temp(P|ln)B1(P|ml)
        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto l = 0; l < NB; ++l)
        for(auto P = 0; P < NBRI; ++P)
        for(auto k = 0; k < NB; ++k) {

          temp[P + l*NBRI + n*NBNBRI] += ERI3_2[P + k*NBRI + n*NBNBRI] * C.X[l + k*NB];

        }

        #pragma omp parallel for
        for(auto n = 0; n < NB; ++n)
        for(auto m = 0; m < NB; ++m)
        for(auto l = 0; l < NB; ++l)
        for(auto P = 0; P < NBRI; ++P) {

          C.AX[m + n*NB] += temp[P + l*NBRI + n*NBNBRI] * ERI3_1[P + m*NBRI + l*NBNBRI];

        }
  #ifdef _REPORT_KCon
        auto durK1 = tock(topK1);
        std::cout << "K1 Contract  duration   = " << durK1 << std::endl;
  #endif

//        size_t LAThreads = GetLAThreads();
//        SetLAThreads(1);
//
//        #pragma omp parallel for
//        for(auto n = 0; n < NB; ++n) {
//          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBRI,NB,NB,MatsT(1.),ERI3_2+n*NBNBRI,NBRI,C.X,NB,MatsT(0.),temp+n*NBNBRI,NBRI);
//        } 
//
//	SetLAThreads(LAThreads);
//
//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB,NB,NBNBRI,MatsT(1.),ERI3_1,NBNBRI,temp,NBNBRI,MatsT(0.),C.AX,NB);

//        // (ij|Q)S^{-1/2} -> ERI3J
//        size_t LAThreads = GetLAThreads();
//        SetLAThreads(1);
//
//        #pragma omp parallel for
//        for(auto nu = 0ul; nu < NB; nu++)
//          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBRI,NB,NB,MatsT(1.),ERI3_2+nu*NBNBRI,NBRI,C.X,NB,MatsT(0.),temp+nu*NBNBRI,NBRI);
//
//        SetLAThreads(LAThreads);
//
//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB,NB,NBNBRI,MatsT(1.),ERI3_1,NBNBRI,temp,NBNBRI,MatsT(0.),C.AX,NB);
      }
      CQMemManager::get().free(temp);
    }
    ProgramTimer::tock("K Contract");

  }; // InCoreRelERIContraction::KContract
  


  // Contraction into separate storages
  template <typename MatsT, typename IntsT>
  void InCore4indexGradContraction<MatsT,IntsT>::gradTwoBodyContract(
    MPI_Comm comm,
    const bool screen,
    std::vector<std::vector<TwoBodyContraction<MatsT>>>& list,
    EMPerturbation& pert) const {

    // Contract over each 3N gradient component
    size_t nGrad = this->grad_.size();
    assert(nGrad == list.size());

    for (auto i = 0; i < nGrad; i++) {
      auto casted = std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(this->grad_[i]);
      InCore4indexTPIContraction<MatsT,IntsT> contraction(casted);
      contraction.contractSecond = this->contractSecond;
      contraction.isCross = this->isCross;
      contraction.twoBodyContract(comm, screen, list[i], pert);
    }

  }; // InCore4indexGradContraction::gradTwoBodyContract

}; // namespace ChronusQ
