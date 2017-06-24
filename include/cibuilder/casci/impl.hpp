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

// #define DEBUG_CI_MU
// #define DEBUG_CI_SIGMA

#define CASCI_LOOP_INIT() \
  size_t nC = mcwfn.reference().nC; \
  size_t NDet = mcwfn.NDet; \
  std::shared_ptr<const ExcitationList> exList_a = \
    std::dynamic_pointer_cast<CASStringManager>( \
      mcwfn.detStr)->excitationList(); \
  std::shared_ptr<const ExcitationList> exList_b = \
    (nC == 1) ?  std::dynamic_pointer_cast<CASStringManager>( \
      mcwfn.detStrBeta)->excitationList() : nullptr; \
  size_t nStr_a  = exList_a->nString(); \
  size_t nStr_b = (exList_b) ? exList_b->nString(): 1; \
  size_t nNZa = exList_a->nNonZero(); \
  size_t nNZb = (exList_b) ? exList_b->nNonZero(): 0; \


namespace ChronusQ {
  /*  
   *  Build full CASCI Hamiltonian  as H(K, L)
   */ 
  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::buildFullH(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * fullH) {  
    
    CASCI_LOOP_INIT(); // check top for variable definitions
    
    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("hCoreP_Correlated_Space"));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));

    // Allocate SCR
    size_t nSCR = std::max(nStr_a, nStr_b);
    size_t nThreads = GetNumThreads();
    MatsT * SCR  = CQMemManager::get().malloc<MatsT>(nSCR * nThreads);
     
    // empty CI Hamiltonian
    std::fill_n(fullH, nStr_a*nStr_a, MatsT(0.));
    // Alpha Part for 1C or the whole build for 2C and 4C
    // NOTE: ex_list_ are row-majored c++ objects 
    
    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    
    size_t nStr_a_nThread = nStr_a * nThreads;
    
    MatsT *CIHCol = nullptr, *SCR_ith = nullptr;
    
    // alpha part
    
    // fullH as (Ka, La), CIHCol as a column
#pragma omp parallel default(shared) private(CIHCol, SCR_ith, La, k, l, Ka, signkl, i, j, Ja, signij)  
    { 
      auto iThread = GetThreadID();
      CIHCol  = fullH + nStr_a * iThread;
      SCR_ith = SCR   + nSCR * iThread;
      for (La = iThread; La < nStr_a; La+=nThreads, CIHCol+=nStr_a_nThread) {
        
        std::fill_n(SCR_ith, nStr_a, MatsT(0.));
        const int * exList_La = exList_a->pointerAtDet(La);
        for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_La+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);
	      SCR_ith[Ka] += signkl * hCoreP(k, l);
         
          const int * exList_Ka = exList_a->pointerAtDet(Ka);
	      for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_Ka+=4) {
             
            UNPACK_EXCITATIONLIST_4(exList_Ka, i, j, Ja, signij);
	        SCR_ith[Ja] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	      }
        }
        
        // passive screening and updating CHCol
        for (Ka = 0 ; Ka < nStr_a; Ka++) {
          if (std::abs(SCR_ith[Ka]) > small_number) CIHCol[Ka] = SCR_ith[Ka];
        }
      
      }  // La
    } 
    
    if (nC != 1) {
       
#ifdef _DEBUG_CIBUILDER_CASCI_IMPL
    prettyPrintSmart(std::cout,"HH full CASCI Hamiltonian",fullH, NDet, NDet, NDet);
#endif
       CQMemManager::get().free(SCR);
       return;
    } 
    
    // 1C Continued: Expand Alpha Part 
    //    (Ka, La) -> (Ka, Kb, La, Kb) 
    size_t nStr_a2 = nStr_a * nStr_a;
    size_t nStr_b2 = nStr_b * nStr_b;
    MatsT * tmpH  = CQMemManager::get().malloc<MatsT>(nSCR*nSCR);
    std::copy_n(fullH, nStr_a2, tmpH);
    std::fill_n(fullH, NDet*NDet, MatsT(0.));
    
    // transpose tmpH as (La, Ka)
    IMatCopy('T', nStr_a, nStr_a, MatsT(1.), tmpH, nStr_a, nStr_a);  
     
    // update fullH as (La, Ka, Kb, Kb)
    // then transpose back to (Ka, Kb, La, Kb) in each last Kb
    size_t lastDimOff = nStr_a2*nStr_b;
    CIHCol = fullH;
    for (Kb = 0; Kb < nStr_b; Kb++, CIHCol+=lastDimOff) {
      std::copy_n(tmpH, nStr_a2, CIHCol+Kb*nStr_a2); 
      IMatCopy('T', nStr_a, NDet, MatsT(1.), CIHCol, nStr_a, NDet);  
    }

    
    // 1C Continued: build Beta part as (Kb, Lb)
    std::fill_n(tmpH, nStr_b2, MatsT(0.));
    size_t nStr_b_nThread = nStr_b * nThreads;
#pragma omp parallel default(shared) private(CIHCol, SCR_ith, Lb, k, l, Kb, signkl, i, j, Jb, signij)  
    {
      auto iThread = GetThreadID();
      CIHCol  = tmpH + nStr_b * iThread;
      SCR_ith = SCR + nSCR * iThread;
      for (Lb = iThread; Lb < nStr_b; Lb+=nThreads, CIHCol+=nStr_b_nThread) {
        
        std::fill_n(SCR_ith, nStr_b, MatsT(0.));
        const int * exList_Lb = exList_b->pointerAtDet(Lb);
        for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Lb+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
          SCR_ith[Kb] += signkl * hCoreP(k, l);
         
          const int * exList_Kb = exList_b->pointerAtDet(Kb);
	      for (auto Eij = 0ul; Eij < nNZb; Eij++, exList_Kb+=4) {
               
            UNPACK_EXCITATIONLIST_4(exList_Kb, i, j, Jb, signij);
	        SCR_ith[Jb] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	  
          }
        }
       
        // passive screening and updating CHCol
        for (Kb = 0 ; Kb < nStr_b; Kb++) {
          if (std::abs(SCR_ith[Kb]) > small_number) CIHCol[Kb] = SCR_ith[Kb];
        }
      }  // Lb
    }
    // 1C Continued: Expand Beta Part 
    //    (Kb, Lb) -> (Ka, Kb, Ka, Lb) 
    
    // TODO: try to used non-transposed algorithm
    // transpose fullH: (Ka, Kb, La, Lb) -> (Kb, La, Lb, Ka) 
    IMatCopy('T', nStr_a, nStr_b*NDet, MatsT(1.), fullH, nStr_a, nStr_b*NDet);  
    
    // transpose tmpH as (Lb, Kb)
    IMatCopy('T', nStr_b, nStr_b, MatsT(1.), tmpH, nStr_b, nStr_b);  
   
    // update fullH as (Lb, Kb, Ka, Ka)
    // then transpose back to (Kb, Ka, Lb, Ka) in each last Ka
    lastDimOff = nStr_b2*nStr_a;
    CIHCol = fullH;
    for (Ka = 0; Ka < nStr_a; Ka++, CIHCol+=lastDimOff) {
      IMatCopy('T', NDet, nStr_b, MatsT(1.), CIHCol, NDet, nStr_b);  
      MatAdd('N', 'N', nStr_b, nStr_b, MatsT(1.), CIHCol+Ka*nStr_b2, nStr_b,
        MatsT(1.), tmpH, nStr_b, CIHCol+Ka*nStr_b2, nStr_b); 
      IMatCopy('T', nStr_b, NDet, MatsT(1.), CIHCol, nStr_b, NDet);  
    }
     
    // transpose fullH: ((Kb, La, Lb, Ka) -> (Ka, Kb, La, Lb) 
    IMatCopy('T', nStr_b*NDet, nStr_a, MatsT(1.), fullH, nStr_b*NDet, nStr_a);  

    CQMemManager::get().free(SCR, tmpH);
    
    // 1C Continued: Alpha-Beta and Beta-Aphla Part 
#pragma omp parallel for schedule(static) default(shared) \
  private(CIHCol, Lb, k, l, Kb, signkl, La, i, j, Ka, signij)  
    for(Lb = 0; Lb < nStr_b; Lb++) {
      
      CIHCol = fullH + NDet*(Lb*nStr_a); // La = 0 
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
    
#ifdef _DEBUG_CIENGINE_CASCI_IMPL
    prettyPrintSmart(std::cout,"HH full CASCI Hamiltonian", fullH, NDet, NDet, NDet);
#endif
  
  }; // CASCI::buildFullH

  /*  
   *  Build Diagonal CASCI Hamiltonian 
   */ 
  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::buildDiagH(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * diagH) {  
    
    CASCI_LOOP_INIT(); // check top for variable definitions
    
    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts, MatsT>("hCoreP_Correlated_Space"));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI, MatsT>("ERI_Correlated_Space"));

	// Allocate SCR
    MatsT SCR;

    // empty CI Diagonal Hamiltonian
    std::fill_n(diagH, nStr_a, MatsT(0.));
    // Alpha Part for 1C or the whole build for 2C and 4C
    
    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    
#pragma omp parallel for schedule(static) default(shared) private(SCR, La, k, l, Ka, signkl, i, j, Ja, signij)  
    for (La = 0; La < nStr_a; La++) {
      
      const int * exList_La = exList_a->pointerAtDet(La);
      SCR = MatsT(0.); 
      for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_La+=4) {
        
        UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);
	    if(Ka == La) SCR += signkl * hCoreP(k, l);

        const int * exList_Ka = exList_a->pointerAtDet(Ka);
	    for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_Ka+=4) {
            
          UNPACK_EXCITATIONLIST_4(exList_Ka, i, j, Ja, signij);
	      if(Ja == La) SCR += 0.5 * signij * signkl * moERI(i, j, k, l); 
	    }
       }
      
      // update diagH
      diagH[La] = SCR;
    
    }  // La
    
    if (nC != 1) { 
#ifdef _DEBUG_CIENGINE_CASCI_IMPL
      prettyPrintSmart(std::cout,"HH full CASCI Hamiltonian",diagH, NDet, 1, NDet);
#endif
      return;
    }
    // 1C Continued: Expand Alpha Part 
    //    (Ka) -> (Ka, Kb)
    MatsT * dH = diagH + nStr_a;
    for (Kb = 1; Kb < nStr_b; Kb++, dH+=nStr_a) std::copy_n(diagH, nStr_a, dH);
     
    MatsT * tmpdH  = CQMemManager::get().malloc<MatsT>(nStr_b);

    // 1C Continued: build Beta part
    std::fill_n(tmpdH, nStr_b, MatsT(0.));
#pragma omp parallel for schedule(static) default(shared) private(SCR, Lb, k, l, Kb, signkl, i, j, Jb, signij)  
    for (Lb = 0; Lb < nStr_b; Lb++) {
      
      SCR = MatsT(0.); 
      const int * exList_Lb = exList_b->pointerAtDet(Lb);
      for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Lb+=4) {
        
        UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
	    if(Kb == Lb) SCR += signkl * hCoreP(k, l);
       
        const int * exList_Kb = exList_b->pointerAtDet(Kb);
	    for (auto Eij = 0ul; Eij < nNZb; Eij++, exList_Kb+=4) {
            
          UNPACK_EXCITATIONLIST_4(exList_Kb, i, j, Jb, signij);
	      if(Jb == Lb) SCR += 0.5 * signij * signkl * moERI(i, j, k, l); 
	    }
      }
     
      // update diagH
      tmpdH[Lb] = SCR;
    
    }  // Lb
    
    // 1C Continued: Expand Beta Part 
    // transpose diagH: (Ka, Kb) -> (Kb, Ka) 
    IMatCopy('T', nStr_a, nStr_b, MatsT(1.), diagH, nStr_a, nStr_b);  
    
    // Expand
    dH = diagH;
    for (Ka = 0; Ka < nStr_a; Ka++)  
    for (Kb = 0; Kb < nStr_b; Kb++, dH++)
      *dH += tmpdH[Kb];
    
    // transpose diagH: (Kb, Ka) -> (Ka, Kb) 
    IMatCopy('T', nStr_b, nStr_a, MatsT(1.), diagH, nStr_b, nStr_a);  
    
    CQMemManager::get().free(tmpdH);
    
    // 1C Continued: Alpha-Beta and Beta-Aphla Part 
    size_t nActEa = mcwfn.MOPartition.nCorrEA;
    size_t nActEb = mcwfn.MOPartition.nCorrEB;
#pragma omp parallel for schedule(static) default(shared) \
  private(dH, Lb, k, l, Kb, signkl, La, i, j, Ka, signij)  
    for(Lb = 0; Lb < nStr_b; Lb++) {
      
      dH = diagH + Lb*nStr_a;
      const int * exList_Lb_head = exList_b->pointerAtDet(Lb);

      for(La = 0; La < nStr_a;  La++, dH++) {
        
        const int * exList_La_head = exList_a->pointerAtDet(La);
      
        auto exList_Lb = exList_Lb_head; 
        for(auto Ekl = 0ul; Ekl < nActEa; Ekl++, exList_Lb+=4) {
	
          UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
          const int * exList_La = exList_La_head; 
        
	      for (auto Eij = 0ul; Eij < nActEb; Eij++, exList_La+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_La, i, j, Ka, signij);
            *dH += signij * signkl * moERI(i, j, k, l); 
          }
        }
      }
    }
  }; // CASCI::buildDiagH
  
  /*
   *  Sigma, Matrix-vector product
   */

  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::buildSigma(MCWaveFunction<MatsT, IntsT> & mcwfn, 
    size_t nVec, MatsT * C, MatsT * Sigma) {
    
    CASCI_LOOP_INIT(); // check top for variable definitions
    
    auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("hCoreP_Correlated_Space"));
    auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));

#ifdef DEBUG_CI_SIGMA
    prettyPrintSmart(std::cout,"HH CASCI Sigma Build -- C", C, NDet, nVec, NDet);
#endif
    
    // Allocate SCR
    size_t nSCR = std::max(nStr_a, nStr_b);
    size_t nThreads = GetNumThreads();
    MatsT * SCR  = CQMemManager::get().malloc<MatsT>(nSCR * nThreads);

    // empty Sigma
    std::fill_n(Sigma, NDet*nVec, MatsT(0.));
    // Alpha Part for 1C or the whole build for 2C and 4C
    // NOTE: exList_ are row-majored c++ objects 
    
    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    
    MatsT *HC, *Ci, *SCR_ith;
    
    // SCR_ith is a row of H instead of a column now
#pragma omp parallel default(shared) private(HC, Ci, SCR_ith, La, k, l, Ka, signkl, \
  i, j, Ja, signij, Kb)
    {
      auto iThread = GetThreadID();
      SCR_ith = SCR + nSCR * iThread;
      for (Ka = iThread; Ka < nStr_a; Ka+=nThreads) {
      
        std::fill_n(SCR_ith, nStr_a, MatsT(0.));
        const int * exList_Ka = exList_a->pointerAtDet(Ka);
        for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_Ka+=4) {
        
          UNPACK_EXCITATIONLIST_4(exList_Ka, l, k, La, signkl);
	      SCR_ith[La] += signkl * hCoreP(k, l);
       
          const int * exList_La = exList_a->pointerAtDet(La);
	      for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_La+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_La, j, i, Ja, signij);
	        SCR_ith[Ja] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	      }
        }
        
        // passive screening 
        std::vector<int> SCR_nonZero_ith;
        for (La = 0 ; La < nStr_a;  La++) {
          if (std::abs(SCR_ith[La]) > small_number) 
	        SCR_nonZero_ith.push_back(La);
        }
        
        // update Sigma
        HC = Sigma;
        Ci = C;
        if(nC != 1) {  // for more than 1C
          for (auto iVec = 0ul; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet) 
          for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
	        La = SCR_nonZero_ith[iSCR]; 
	        HC[Ka] += SCR_ith[La] * Ci[La]; 
          }	
        } else { // 1C
          for (auto iVec = 0ul; iVec < nVec; iVec++) 
          for (Kb = 0; Kb < nStr_b; Kb++, HC+=nStr_a, Ci+=nStr_a)
	      for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
	        La = SCR_nonZero_ith[iSCR]; 
	        HC[Ka] += SCR_ith[La] * Ci[La]; 
          }	
        } 
      }  // La
    }

    if (nC != 1) {
       
#ifdef DEBUG_CI_SIGMA
    prettyPrintSmart(std::cout,"HH CASCI Sigma Build -- Sigma", Sigma, NDet, nVec, NDet);
#endif
       CQMemManager::get().free(SCR);
       return;
    } 
    
    // 1C Continued: Build Beta part
    
    // transpose sigma and C to make update continuous in inner loop
    // (a, b) -> (b, a) 
    HC = Sigma;
    Ci = C;
    for (auto iVec = 0ul; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet) {
      IMatCopy('T', nStr_a, nStr_b, MatsT(1.), HC, nStr_a, nStr_b);  
      IMatCopy('T', nStr_a, nStr_b, MatsT(1.), Ci, nStr_a, nStr_b);  
    }

#pragma omp parallel default(shared) private(HC, Ci, SCR_ith, Lb, k, l, Kb, signkl, \
  i, j, Jb, signij, Ka)
    {
      auto iThread = GetThreadID();
      SCR_ith = SCR + nSCR * iThread;
      for (Kb = iThread; Kb < nStr_b; Kb+=nThreads) {
        
        std::fill_n(SCR_ith, nStr_b, MatsT(0.));
        const int * exList_Kb = exList_b->pointerAtDet(Kb);
        for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Kb+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_Kb, l, k, Lb, signkl);
	      SCR_ith[Lb] += signkl * hCoreP(k,l);
         
          const int * exList_Lb = exList_b->pointerAtDet(Lb);
	      for (auto Eij = 0ul; Eij < nNZb; Eij++, exList_Lb+=4) {
              
            UNPACK_EXCITATIONLIST_4(exList_Lb, j, i, Jb, signij);
	        SCR_ith[Jb] += 0.5 * signij * signkl * moERI(i, j, k, l); 
	      }
        }
       
        // passive screening 
        std::vector<int> SCR_nonZero_ith;
        for (Lb = 0 ; Lb < nStr_b;  Lb++) {
          if (std::abs(SCR_ith[Lb]) > small_number) 
	        SCR_nonZero_ith.push_back(Lb);
        }
        
        // update Sigma
        HC = Sigma;
        Ci = C;
        for (auto iVec = 0ul; iVec < nVec; iVec++) 
        for (Ka = 0; Ka < nStr_a; Ka++, HC+=nStr_b, Ci+=nStr_b)
        for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
          Lb = SCR_nonZero_ith[iSCR]; 
          HC[Kb] += SCR_ith[Lb] * Ci[Lb]; 
        }	
      }  // Lb
    }

    // transpose sigma and C back
    HC = Sigma;
    Ci = C;
    // (b, a) -> (a, b)
    for (auto iVec = 0ul; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet) { 
      IMatCopy('T', nStr_b, nStr_a, MatsT(1.), HC, nStr_b, nStr_a);  
      IMatCopy('T', nStr_b, nStr_a, MatsT(1.), Ci, nStr_b, nStr_a);  
    }

    CQMemManager::get().free(SCR);

    // 1C Continued: Alpha-Beta and Beta-Aphla Part 
    // TODO: try to vectorized the loop
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

              tmpH  = signij * signkl * moERI(i, j, k, l);
	  
              if (std::abs(tmpH) > small_number) {
                HC    = Sigma;
                Ci    = C;
                LAddr = La + Lb*nStr_a;
                for (auto iVec = 0ul; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet) 
                  HC[KAddr] += tmpH*Ci[LAddr]; 
              }
            }
          }
        }  
      }
    }

#ifdef DEBUG_CI_SIGMA
    prettyPrintSmart(std::cout,"HH Sigma Hamiltonian -- Sigma", Sigma, NDet, nVec, NDet);
#endif
  
  } // CASCI::buildSigma

  /*
   *  Mu, Matrix-vector product
   */

  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::buildMu(MCWaveFunction<MatsT, IntsT> & mcwfn, 
    size_t nVec, MatsT * C, MatsT * Mu, EMPerturbation & pert) {
    
    auto dipAmp = pert.getDipoleAmp(Electric);
    CASCI_LOOP_INIT(); // check top for variable definitions

    // dont! empty Mu. CONSUMERS MUST EMPTY ON THEIR END
    //std::fill_n(Mu, NDet*nVec, MatsT(0.));
    bool return_early = true;
    for(auto i = 0;    i < 3;     i++){
      if (dipAmp[i] != 0.0) {
	return_early = false;
      }
    }
    if (return_early)
      return;

    // dipole AO -> MO transformation
    auto MOdipole = mcwfn.moints->template getIntegral<VectorInts,MatsT>("MOdipole");
    size_t nAO = mcwfn.reference().nAlphaOrbital() * mcwfn.reference().nC;
    size_t nCorrO = mcwfn.MOPartition.nCorrO;
    size_t nInact = mcwfn.MOPartition.nInact;
    if (not MOdipole) {
      std::shared_ptr<VectorInts<IntsT>> AOdipole =
                std::make_shared<VectorInts<IntsT>>( nAO, 1, true);
      std::shared_ptr<VectorInts<MatsT>> MOdipole_scr =
                std::make_shared<VectorInts<MatsT>>( nCorrO, 1, true);

      std::vector<std::pair<size_t, size_t>> active(2, {mcwfn.MOPartition.nFCore+nInact, nCorrO});

      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (mcwfn.referenceWaveFunction().nC == 1)
            (*AOdipole)[iXYZ] = std::make_shared<OnePInts<IntsT>>( *((*mcwfn.reference().aoints_->lenElectric)[iXYZ]) );
        else if (mcwfn.referenceWaveFunction().nC == 2)
            (*AOdipole)[iXYZ] = std::make_shared<OnePInts<IntsT>>( (*mcwfn.reference().aoints_->lenElectric)[iXYZ]->template spatialToSpinBlock<IntsT>() ) ;
        
        (*AOdipole)[iXYZ]->subsetTransform('N',mcwfn.reference().mo[0].pointer(),
            nAO, active, (*MOdipole_scr)[iXYZ]->pointer(), false);
      }

      mcwfn.moints->addIntegral("MOdipole", MOdipole_scr);
    }

    MOdipole = mcwfn.moints->template getIntegral<VectorInts,MatsT>("MOdipole");

#ifdef DEBUG_CI_SIGMA
    prettyPrintSmart(std::cout,"SU CASCI Mu Build -- C", C, NDet, nVec, NDet);
#endif
    
    // Allocate SCR
    size_t nSCR = std::max(nStr_a, nStr_b);
    size_t nThreads = GetNumThreads();
    MatsT * SCR  = CQMemManager::get().malloc<MatsT>(nSCR * nThreads);

    // Alpha Part for 1C or the whole build for 2C and 4C
    // NOTE: exList_ are row-majored c++ objects 
    
    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
    double small_number = std::numeric_limits<double>::epsilon();
    
    MatsT *MuC, *Ci, *SCR_ith;
    
    // SCR_ith is a row of H instead of a column now
#pragma omp parallel default(shared) private(MuC, Ci, SCR_ith, La, k, l, Ka, signkl, \
  i, j, Ja, signij, Kb)
    {
      auto iThread = GetThreadID();
      SCR_ith = SCR + nSCR * iThread;
      for (Ka = iThread; Ka < nStr_a; Ka+=nThreads) {
      
        std::fill_n(SCR_ith, nStr_a, MatsT(0.));
        const int * exList_Ka = exList_a->pointerAtDet(Ka);
        for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_Ka+=4) {
        
          UNPACK_EXCITATIONLIST_4(exList_Ka, l, k, La, signkl);
          for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
              if (dipAmp[iXYZ] != 0.0) {
	      SCR_ith[La] += signkl * (-dipAmp[iXYZ]) *(*MOdipole)[iXYZ]->pointer()[k *nCorrO +  l];
              }
	  }
        }

        // passive screening 
        std::vector<int> SCR_nonZero_ith;
        for (La = 0 ; La < nStr_a;  La++) {
          if (std::abs(SCR_ith[La]) > small_number) 
	        SCR_nonZero_ith.push_back(La);
        }

        // update Mu
        MuC = Mu;
        Ci = C;
        if(nC != 1) {  // for more than 1C
          for (auto iVec = 0ul; iVec < nVec; iVec++, MuC+=NDet, Ci+=NDet) 
          for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
	        La = SCR_nonZero_ith[iSCR]; 
	        MuC[Ka] += SCR_ith[La] * Ci[La]; 
          }	
        } else { // 1C
          for (auto iVec = 0ul; iVec < nVec; iVec++) 
          for (Kb = 0; Kb < nStr_b; Kb++, MuC+=nStr_a, Ci+=nStr_a)
	      for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
	        La = SCR_nonZero_ith[iSCR]; 
	        MuC[Ka] += SCR_ith[La] * Ci[La]; 
          }	
        } 
      }  // La
    }

    if (nC != 1) {
       
#ifdef DEBUG_CI_SIGMA
    prettyPrintSmart(std::cout,"SU Mu Build -- Mu", Mu, NDet, nVec, NDet);
#endif
       CQMemManager::get().free(SCR);
       return;
    } 
    
    // 1C Continued: Build Beta part
    
    // transpose sigma and C to make update continuous in inner loop
    // (a, b) -> (b, a) 
    MuC = Mu;
    Ci = C;
    for (auto iVec = 0ul; iVec < nVec; iVec++, MuC+=NDet, Ci+=NDet) {
      IMatCopy('T', nStr_a, nStr_b, MatsT(1.), MuC, nStr_a, nStr_b);  
      IMatCopy('T', nStr_a, nStr_b, MatsT(1.), Ci, nStr_a, nStr_b);  
    }

#pragma omp parallel default(shared) private(MuC, Ci, SCR_ith, Lb, k, l, Kb, signkl, \
  i, j, Jb, signij, Ka)
    {
      auto iThread = GetThreadID();
      SCR_ith = SCR + nSCR * iThread;
      for (Kb = iThread; Kb < nStr_b; Kb+=nThreads) {
        
        std::fill_n(SCR_ith, nStr_b, MatsT(0.));
        const int * exList_Kb = exList_b->pointerAtDet(Kb);
        for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Kb+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_Kb, l, k, Lb, signkl);
	   for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
              if (dipAmp[iXYZ] != 0.0) {
	      SCR_ith[Lb] += signkl * (-dipAmp[iXYZ]) *(*MOdipole)[iXYZ]->pointer()[k *nCorrO +  l];
              }
	  }
        }
        // passive screening 
        std::vector<int> SCR_nonZero_ith;
        for (Lb = 0 ; Lb < nStr_b;  Lb++) {
          if (std::abs(SCR_ith[Lb]) > small_number) 
	        SCR_nonZero_ith.push_back(Lb);
        }
        
        // update Mu
        MuC = Mu;
        Ci = C;
        for (auto iVec = 0ul; iVec < nVec; iVec++) 
        for (Ka = 0; Ka < nStr_a; Ka++, MuC+=nStr_b, Ci+=nStr_b)
        for (auto iSCR = 0ul; iSCR < SCR_nonZero_ith.size(); iSCR++) {
          Lb = SCR_nonZero_ith[iSCR]; 
          MuC[Kb] += SCR_ith[Lb] * Ci[Lb]; 
        }	
      }  // Lb
    }

    // transpose sigma and C back
    MuC = Mu;
    Ci = C;
    // (b, a) -> (a, b)
    for (auto iVec = 0ul; iVec < nVec; iVec++, MuC+=NDet, Ci+=NDet) { 
      IMatCopy('T', nStr_b, nStr_a, MatsT(1.), MuC, nStr_b, nStr_a);  
      IMatCopy('T', nStr_b, nStr_a, MatsT(1.), Ci, nStr_b, nStr_a);  
    }

    CQMemManager::get().free(SCR);


#ifdef DEBUG_CI_MU
    prettyPrintSmart(std::cout,"SU Mu Hamiltonian -- Mu", Mu, NDet, nVec, NDet);
#endif
  
  } // CASCI::buildMu
  

  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::computeOneRDM(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * C, 
    cqmatrix::Matrix<MatsT> & oneRDM) {
  
    computeTDM(mcwfn, C, C, oneRDM);
  
  } // CASCI::computeOneRDM 

  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::computeTwoRDM(MCWaveFunction<MatsT, IntsT> & mcwfn, 
    MatsT * C, InCore4indexTPI<MatsT> & twoRDM) {
       
    CASCI_LOOP_INIT(); // check top for variable definitions
    
    
    size_t nThreads = GetNumThreads();
    std::vector<InCore4indexTPI<MatsT>> SCR;
    auto nDim  = twoRDM.nBasis();
    for (auto i = 0ul; i < nThreads; i++)
      SCR.emplace_back(nDim);

    // alpha-alpha part
    int i, j, k, l, La, Lb, Ka, Kb, Ja, Jb;
    double signij, signkl;
#pragma omp parallel default(shared) private(i, j, k, l, La, Ka, Ja, Lb,\
  signij, signkl)
    {
      auto iThread = GetThreadID();
      auto & tmpRDM = SCR[iThread];
      tmpRDM.clear();
      for (La = iThread; La < nStr_a; La+=nThreads) {
        
        const int * exList_La = exList_a->pointerAtDet(La);
        for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_La+=4) {
          
          UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);
          const int * exList_Ka = exList_a->pointerAtDet(Ka);
   	  
          for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_Ka+=4) {
              
            UNPACK_EXCITATIONLIST_4(exList_Ka, i, j, Ja, signij);
          
            auto tmp = MatsT(0.); 
            for (auto Lb = 0; Lb < nStr_b; Lb++) 
              tmp += SmartConj(C[Ja + Lb*nStr_a]) * C[La + Lb*nStr_a];
            
            tmpRDM(i, j, k, l) += tmp * signkl * signij; 
	      }
        }
      }
    }


    if (nC == 1) {
    
#pragma omp parallel default(shared) private(i, j, k, l, Lb, Kb, Jb, La, Ka, \
  signij, signkl)
      {
        auto iThread = GetThreadID();
        auto & tmpRDM = SCR[iThread];
        for (Lb = iThread; Lb < nStr_b; Lb+=nThreads) {
          
          const int * exList_Lb = exList_b->pointerAtDet(Lb);
          for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Lb+=4) {
            
            UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);
            
            // beta-beta part
            const int * exList_Kb = exList_b->pointerAtDet(Kb);
            for (auto Eij = 0ul; Eij < nNZb; Eij++, exList_Kb+=4) {
              
              UNPACK_EXCITATIONLIST_4(exList_Kb, i, j, Jb, signij);
              
              auto tmp = MatsT(0.);
              for (La = 0; La < nStr_a; La++) 
                tmp += SmartConj(C[La + Jb*nStr_a]) * C[La + Lb*nStr_a];

              tmpRDM(i, j, k, l) += tmp * signkl * signij; 
            }
          
            // alpha-beta and beta-alpha part
            for(La = 0; La < nStr_a; La++) {
              const int * exList_La = exList_a->pointerAtDet(La);
	          for (auto Eij = 0ul; Eij < nNZa; Eij++, exList_La+=4) {
                
                UNPACK_EXCITATIONLIST_4(exList_La, i, j, Ka, signij);
                tmpRDM(i, j, k, l) += 2.0 * signkl * signij *  SmartConj(C[Ka + Kb*nStr_a]) * C[La + Lb*nStr_a];
              }
            }
          
          }
        }
      } 
    } 
    
    twoRDM.clear();
    auto nDim2 = nDim * nDim;
    auto nDim4 = nDim2 * nDim2;
    for (auto i = 0ul; i < nThreads; i++) 
      blas::axpy(nDim4, 1.0, SCR[i].pointer(), 1, twoRDM.pointer(), 1);
    
    return;
  } // CASCI::computeTwoRDM 

  template <typename MatsT, typename IntsT>
  void CASCI<MatsT,IntsT>::computeTDM(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * Cm,
        MatsT * Cn, cqmatrix::Matrix<MatsT> & TDM) {

    CASCI_LOOP_INIT(); // check top for variable definitions

    size_t nThreads = GetNumThreads();
    std::vector<cqmatrix::Matrix<MatsT>> SCR;
    for (auto i = 0ul; i < nThreads; i++)
      SCR.emplace_back(TDM.dimension());

    // alpha part 
    int k, l, La, Lb, Ka, Kb;
    double signkl;
#pragma omp parallel default(shared) private(La, Lb, k, l, Ka, signkl)
    {
      auto iThread = GetThreadID();
      auto & tmpRDM = SCR[iThread];
      tmpRDM.clear();
      for (La = iThread; La < nStr_a; La+=nThreads) {

        const int * exList_La = exList_a->pointerAtDet(La);

        for (auto Ekl = 0ul; Ekl < nNZa; Ekl++, exList_La+=4) {

          UNPACK_EXCITATIONLIST_4(exList_La, k, l, Ka, signkl);

          auto tmp = MatsT(0.);
          for (Lb = 0; Lb < nStr_b; Lb++)
            tmp += SmartConj(Cm[Ka + Lb*nStr_a]) * Cn[La + Lb*nStr_a];

          tmpRDM(k, l) += tmp * signkl;
        }
      }
    }

    if (nC == 1) {
    // beta part
#pragma omp parallel default(shared) private(La, Lb, k, l, Kb, signkl)
      {
        auto iThread = GetThreadID();
        auto & tmpRDM = SCR[iThread];
        for (Lb = iThread; Lb < nStr_b; Lb+=nThreads) {

          const int * exList_Lb = exList_b->pointerAtDet(Lb);
          for (auto Ekl = 0ul; Ekl < nNZb; Ekl++, exList_Lb+=4) {

            UNPACK_EXCITATIONLIST_4(exList_Lb, k, l, Kb, signkl);

            auto tmp = MatsT(0.);
            for (La = 0; La < nStr_a; La++)
              tmp += SmartConj(Cm[La + Kb*nStr_a]) * Cn[La + Lb*nStr_a];

            tmpRDM(k, l) += tmp * signkl;
          }
        }
      }
    }

    TDM.clear();
    for (auto i = 0ul; i < nThreads; i++) TDM += SCR[i];

    return;

  } // CASCI::computeTDM


}; // namespace ChronusQ

