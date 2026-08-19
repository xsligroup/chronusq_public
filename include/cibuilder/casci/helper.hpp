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

#include <mcwavefunction.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <detstringmanager.hpp>
#include <cibuilder/casci.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildFullHOneParticle(MCWaveFunction<MatsT,IntsT> & mcwfn, 
                                                     MatsT * HBlock,
                                                     const std::shared_ptr<const ExcitationList> exList,
                                                     const std::string OneBodyIntsString,
                                                     const std::string TwoBodyIntsString)
  {

    size_t nStr = exList->nString();
    size_t nNZ = exList->nNonZero();

    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts, MatsT>(OneBodyIntsString));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    // Allocate SCR
    size_t nThreads = GetNumThreads();
    MatsT * SCR  = CQMemManager::get().malloc<MatsT>(nStr * nThreads);

    size_t i, j, k, l, La, Ka, Ja;
    double signij; 
    double signkl;
    double small_number = std::numeric_limits<double>::epsilon();

    size_t nStr_nThread = nStr * nThreads;
    MatsT *CIHCol = nullptr; 
    MatsT *SCR_ith = nullptr;

    // fullH as (Ka, La), CIHCol as a column
#pragma omp parallel default(shared) private(CIHCol, SCR_ith, La, k, l, Ka, signkl, i, j, Ja, signij)  
    { 
      auto iThread = GetThreadID();
      CIHCol  = HBlock + nStr * iThread;
      SCR_ith = SCR   + nStr * iThread;
      for (La = iThread; La < nStr; La+=nThreads, CIHCol+=nStr_nThread) {
        
        std::fill_n(SCR_ith, nStr, MatsT(0.));
        const int * exList_La = exList->pointerAtDet(La);
        for (auto Ekl = 0ul; Ekl < nNZ; Ekl++, exList_La+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);
	      SCR_ith[Ka] += signkl * hCoreP(k, l);
         
          const int * exList_Ka = exList->pointerAtDet(Ka);
	      for (auto Eij = 0ul; Eij < nNZ; Eij++, exList_Ka+=4) {
             
            UNPACK_EXCITATIONLIST_4(exList_Ka, i, j, Ja, signij);
	        SCR_ith[Ja] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	      }
        }
        
        // passive screening and updating CHCol
        for (Ka = 0 ; Ka < nStr; Ka++) {
          if (std::abs(SCR_ith[Ka]) > small_number) CIHCol[Ka] = SCR_ith[Ka];
        }
      
      }  // La
    } 

    CQMemManager::get().free(SCR);
 
  }; // CASCIHelper::BuildFullHOneParticle

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildFullHTwoParticle(MCWaveFunction<MatsT,IntsT> & mcwfn, 
                                                     MatsT * HBlock,
                                                     const std::shared_ptr<const ExcitationList> exList_a,
                                                     const std::shared_ptr<const ExcitationList> exList_b,
                                                     const std::string TwoBodyIntsString)
  {

    size_t nStr_a = exList_a->nString();
    size_t nNZa = exList_a->nNonZero();
    size_t nStr_b = exList_b->nString();
    size_t nNZb = exList_b->nNonZero();
    size_t NDet = nStr_a * nStr_b;

    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    // Allocate SCR
    size_t nSCR = std::max(nStr_a, nStr_b);

    size_t i, j, k, l, La, Ka, Lb, Kb;
    double signij; 
    double signkl;
    double small_number = std::numeric_limits<double>::epsilon();

    MatsT *CIHCol = nullptr; 

    // 1C Continued: Alpha-Beta and Beta-Aphla Part 
#pragma omp parallel for schedule(static) default(shared) \
  private(CIHCol, Lb, k, l, Kb, signkl, La, i, j, Ka, signij)  
    for(Lb = 0; Lb < nStr_b; Lb++) {
      
      CIHCol = HBlock + NDet*(Lb*nStr_a); // La = 0 
      const int * exList_Lb_head = exList_b->pointerAtDet(Lb);
      
      for(La = 0; La < nStr_a; La++, CIHCol+=NDet) {
         
        const int * exList_La_head = exList_a->pointerAtDet(La);
        
        const int * exList_Lb = exList_Lb_head;
        for(auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Lb+=4) {
	
          UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
	    
          const int * exList_La = exList_La_head;
          for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_La+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_La, i, j, Ka, signij);
            CIHCol[Ka + Kb*nStr_a] += signij * signkl * moERI(i, j, k, l); 
          }
        }
      }
    }

  }; // CASHelper::BuildFullHTwoParticle

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildDiagHOneParticle(MCWaveFunction<MatsT,IntsT> & mcwfn, 
                                                     MatsT * HDiag,
                                                     const std::shared_ptr<const ExcitationList> exList,
                                                     const std::string OneBodyIntsString,
                                                     const std::string TwoBodyIntsString)
  {
    size_t nStr = exList->nString();
    size_t nNZ = exList->nNonZero();

    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts, MatsT>(OneBodyIntsString));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    int i, j, k, l, La, Ka, Kb, Ja;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();

    MatsT SCR;

#pragma omp parallel for schedule(static) default(shared) private(SCR, La, k, l, Ka, signkl, i, j, Ja, signij)  
    for (La = 0; La < nStr; La++) {
      
      const int * exList_La = exList->pointerAtDet(La);
      SCR = MatsT(0.); 
      for (auto Ekl = 0ul; Ekl < nNZ; Ekl++, exList_La+=4) {
        
        UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);
	    if(Ka == La) SCR += signkl * hCoreP(k, l);

        const int * exList_Ka = exList->pointerAtDet(Ka);
	    for (auto Eij = 0ul; Eij < nNZ; Eij++, exList_Ka+=4) {
            
          UNPACK_EXCITATIONLIST_4(exList_Ka, i, j, Ja, signij);
	      if(Ja == La) SCR += 0.5 * signij * signkl * moERI(i, j, k, l); 
	    }
       }
      
      // update diagH
      HDiag[La] = SCR;
    
    }  // La

  }; // CASHelper::BuildDiagHOneParticle

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildDiagHTwoParticle(MCWaveFunction<MatsT,IntsT> & mcwfn, 
                                                     MatsT * HDiag,
                                                     const size_t nActA,
                                                     const size_t nActB,
                                                     const std::shared_ptr<const ExcitationList> exList_a,
                                                     const std::shared_ptr<const ExcitationList> exList_b,
                                                     const std::string TwoBodyIntsString)
  {

    size_t nStr_a = exList_a->nString();
    size_t nNZa = exList_a->nNonZero();
    size_t nStr_b = exList_b->nString();
    size_t nNZb = exList_b->nNonZero();
    size_t NDet = nStr_a * nStr_b;

    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    MatsT * dH;

#pragma omp parallel for schedule(static) default(shared) \
  private(dH, Lb, k, l, Kb, signkl, La, i, j, Ka, signij)  
    for(Lb = 0; Lb < nStr_b; Lb++) {
      
      dH = HDiag + Lb*nStr_a;
      const int * exList_Lb_head = exList_b->pointerAtDet(Lb);

      for(La = 0; La < nStr_a;  La++, dH++) {
        
        const int * exList_La_head = exList_a->pointerAtDet(La);
      
        auto exList_Lb = exList_Lb_head; 
        for(auto Ekl = 0ul; Ekl < nActB; Ekl++, exList_Lb+=4) {
	
          UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
          const int * exList_La = exList_La_head; 
        
	      for (auto Eij = 0ul; Eij < nActA; Eij++, exList_La+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_La, i, j, Ka, signij);
            *dH += signij * signkl * moERI(i, j, k, l); 
          }
        }
      }
    }

  }; // CASHelper::BuildDiagHTwoParticle

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildSigmaOneParticle(MCWaveFunction<MatsT,IntsT> & mcwfn,
                                                     MatsT * C,
                                                     MatsT * Sigma,
                                                     size_t nVec,
                                                     size_t nAuxDet,
                                                     std::shared_ptr<const ExcitationList> exList,
                                                     const std::string OneBodyIntsString,
                                                     const std::string TwoBodyIntsString)
  {
    size_t nStr = exList->nString();
    size_t nNZ = exList->nNonZero();

    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts, MatsT>(OneBodyIntsString));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    int i, j, k, l, La, Ka, Ja, Kb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();

    size_t nThreads = GetNumThreads();
    MatsT * SCR  = CQMemManager::get().malloc<MatsT>(nStr* nThreads);
    
    MatsT *HC, *Ci, *SCR_ith;

#pragma omp parallel default(shared) private(HC, Ci, SCR_ith, La, k, l, Ka, signkl, \
  i, j, Ja, signij, Kb)
    {
      auto iThread = GetThreadID();
      SCR_ith = SCR + nStr * iThread;
      for (Ka = iThread; Ka < nStr; Ka+=nThreads) {
      
        std::fill_n(SCR_ith, nStr, MatsT(0.));
        const int * exList_Ka = exList->pointerAtDet(Ka);
        for (auto Ekl = 0ul; Ekl < nNZ; Ekl++, exList_Ka+=4) {
        
          UNPACK_EXCITATIONLIST_4(exList_Ka, l, k, La, signkl);
	      SCR_ith[La] += signkl * hCoreP(k, l);
       
          const int * exList_La = exList->pointerAtDet(La);
	      for (auto Eij = 0ul; Eij < nNZ; Eij++, exList_La+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_La, j, i, Ja, signij);
	        SCR_ith[Ja] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	      }
        }
        
        // passive screening 
        std::vector<int> SCR_nonZero_ith;
        for (La = 0 ; La < nStr;  La++) {
          if (std::abs(SCR_ith[La]) > small_number) 
	        SCR_nonZero_ith.push_back(La);
        }
        
        // update Sigma
        HC = Sigma;
        Ci = C;
        for(auto iVec = 0ul; iVec < nVec; iVec++)
        for(auto KAux = 0ul; KAux < nAuxDet; KAux++, HC+=nStr, Ci+=nStr)
        for(auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++)
        {
          La = SCR_nonZero_ith[iSCR];
          HC[Ka] += SCR_ith[La] * Ci[La];
        }
      }  // La
    }


    CQMemManager::get().free(SCR);
  }

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::buildSigmaTwoParticle(MCWaveFunction<MatsT,IntsT> & mcwfn,
                                                     MatsT * C,
                                                     MatsT * Sigma,
                                                     size_t nVec,
                                                     size_t nAuxDet,
                                                     std::shared_ptr<const ExcitationList> exList_a,
                                                     std::shared_ptr<const ExcitationList> exList_b,
                                                     const std::string TwoBodyIntsString,
                                                     const double chargeproduct)
  {
    size_t nStr_a = exList_a->nString();
    size_t nNZa = exList_a->nNonZero();
    size_t nStr_b = exList_b->nString();
    size_t nNZb = exList_b->nNonZero();
    size_t NDet = nStr_a*nStr_b;

    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>(TwoBodyIntsString));

    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    size_t nThreads = GetNumThreads();

    MatsT * dH;
    MatsT *HC, *Ci, *SCR_ith;

#pragma omp parallel default(shared) private(HC, Ci, SCR_ith, Lb, k, l, Kb, signkl, \
  i, j, La, signij, Ka)
    {  
      MatsT tmpH;
      int KAddr, LAddr;
      auto iThread = GetThreadID();
      for(Kb = iThread; Kb < nStr_b; Kb+=nThreads) { 
        
        const int * exList_Kb_head = exList_b->pointerAtDet(Kb);
       
        for(Ka = 0, KAddr = Kb * nStr_a; Ka < nStr_a; Ka++, KAddr++) {
          
          const int * exList_Ka_head = exList_a->pointerAtDet(Ka);
          
          auto exList_Kb = exList_Kb_head;
          for(auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Kb+=4) {
         
            UNPACK_EXCITATIONLIST_4(exList_Kb, l, k, Lb, signkl);
            
            auto exList_Ka = exList_Ka_head;
	        for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_Ka+=4) {

              UNPACK_EXCITATIONLIST_4(exList_Ka, j, i, La, signij);

              tmpH  = chargeproduct * signij * signkl * moERI(i, j, k, l);
	  
              if (std::abs(tmpH) > small_number) {
                HC    = Sigma;
                Ci    = C;
                LAddr = La + Lb*nStr_a;
                for (auto iVec = 0ul; iVec < nVec; iVec++)
                for (auto KAux = 0ul; KAux < nAuxDet; KAux++, HC+=NDet, Ci+=NDet)
                  HC[KAddr] += tmpH*Ci[LAddr]; 
              }
            }
          }
        }  
      }
    }

  }

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::computeTDM(MCWaveFunction<MatsT,IntsT> & mcwfn,
                                          MatsT * Cm,
                                          MatsT * Cn,
                                          size_t nAuxDet,
                                          std::shared_ptr<const ExcitationList> exList,
                                          cqmatrix::Matrix<MatsT> & TDM)
  {
    size_t nThreads = GetNumThreads();
    std::vector<cqmatrix::Matrix<MatsT>> SCR;
    for(size_t i = 0; i < nThreads; i++)
      SCR.emplace_back(TDM.nRows());

    size_t nStr = exList->nString();
    size_t nNZ = exList->nNonZero();

    int k, l, La, Ka, Lb;
    double signkl;

#pragma omp parallel default(shared) private(La, Lb, k, l, Ka, signkl)
  {
    auto iThread = GetThreadID();
    auto & tmpRDM = SCR[iThread];
    tmpRDM.clear();
    for(La = iThread; La < nStr; La+=nThreads)
    {

      const int * exList_La = exList->pointerAtDet(La);

      for(auto Ekl = 0ul; Ekl < nNZ; Ekl++, exList_La+=4)
      {
        UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);

        auto tmp = MatsT(0.);
        for (Lb = 0; Lb < nAuxDet; Lb++)
          tmp += SmartConj(Cm[Ka + Lb*nStr]) * Cn[La + Lb*nStr];

        tmpRDM(k, l) += tmp * signkl;

      }
    }
  }

    for (auto i = 0ul; i < nThreads; i++) TDM += SCR[i];

  }


  template <typename MatsT, typename IntsT>
  template <typename... MatsArgs>
  void CASHelper<MatsT,IntsT>::transposeVectors(size_t nVec,
                                                size_t dimA,
                                                size_t dimB,
                                                MatsArgs ... Vecs)
  {
    ([&](MatsT* vec){
      for(size_t iVec = 0; iVec < nVec; iVec++, vec += dimA*dimB)
        IMatCopy('T',dimA,dimB,MatsT(1.0),vec,dimA,dimB);
    }(Vecs), ...);
  } // CASHelper::transposeVectors

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::addBlockToMatrixBlockDiagonal(size_t NDetToAdd,
                                                             size_t NDetFull,
                                                             MatsT scale,
                                                             MatsT * SubBlock,
                                                             MatsT * FullMat)
  {
    size_t NAuxDet = NDetFull / NDetToAdd;
    MatsT * tmp = FullMat;
    for(size_t nAux = 0; nAux < NAuxDet; nAux++, tmp += NDetFull*NDetToAdd + NDetToAdd)
    {
      MatAdd('N','N',NDetToAdd,NDetToAdd,scale,SubBlock,NDetToAdd,MatsT(1.0),tmp,NDetFull,tmp,NDetFull);
    }
  }

  template <typename MatsT, typename IntsT>
  void CASHelper<MatsT,IntsT>::addVecToBlockedVector(size_t NDetToAdd,
                                                     size_t NDetFull,
                                                     MatsT scale,
                                                     MatsT * VecToAdd,
                                                     MatsT * FullVec)
  {
    size_t NAuxDet = NDetFull / NDetToAdd;
    MatsT * tmp = FullVec;
    for(size_t nAux = 0; nAux < NAuxDet; nAux++, tmp += NDetToAdd)
    {
      blas::axpy(NDetToAdd,scale,VecToAdd,1,tmp,1);
    }

  }

}; // namespace ChronusQ

