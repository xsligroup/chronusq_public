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

#include <mointstransformer.hpp>
#include <util/timer.hpp>
#include <fockbuilder/fourcompfock/batchgd.hpp>
#include <cqlinalg.hpp>
#include <matrix.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasutil.hpp>

namespace ChronusQ {

  #define ALLOCATE_AND_CLEAR_CACHE_IF_NECESSARY(ALLOCATION) \
    try { ALLOCATION; } \
    catch (...) { \
      ints_cache_.clear(); \
      try { ALLOCATION; } \
      catch (...) { CErr("Not Enough Memory in subsetTransformTPISSFockN6");} \
    }
 
  /**
   *  \brief transform AO TPI to form MO TPI
   *  thru SingleSlater formfock 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::subsetTransformTPISSFockN6(EMPerturbation & pert, 
    const std::vector<std::pair<size_t,size_t>> &off_sizes, MatsT* MOTPI, 
    const std::string & moType, bool cacheIntermediates,
    TPI_TRANS_DELTA_TYPE delta) {

    size_t poff = off_sizes[0].first;
    size_t qoff = off_sizes[1].first;
    size_t roff = off_sizes[2].first;
    size_t soff = off_sizes[3].first;
    size_t np = off_sizes[0].second;
    size_t nq = off_sizes[1].second;
    size_t nr = off_sizes[2].second;
    size_t ns = off_sizes[3].second;

    size_t nAO  = ss_.nAlphaOrbital() * ss_.nC;
    
    // Initial Batching if can not hold all the HalfTransTPI
    
    size_t npq = (delta == KRONECKER_DELTA_PQ or 
                  delta == KRONECKER_DELTA_PQ_RS ) ? np: np*nq;

    size_t nrs = (delta == KRONECKER_DELTA_RS or 
                  delta == KRONECKER_DELTA_PQ_RS ) ? nr: nr*ns;
    
    size_t halfTMOTPISize = nAO * nAO * std::min(npq, nrs); 
    
    auto nBatch = CQMemManager::get().max_avail_allocatable<MatsT>(halfTMOTPISize, 2);
    
    // clear cache for more memory 
    if ( nBatch <= 1) {
      ints_cache_.clear();
      nBatch = CQMemManager::get().max_avail_allocatable<MatsT>(halfTMOTPISize, 2);
    }

    // No batch at this level if remaining memory can allocate at 
    // least two halfTransTPI 
    if (nBatch > 1) {
      performSubsetTransformTPISSFockN6(pert, off_sizes, MOTPI, 
        moType, cacheIntermediates, delta);
      return;
    } 
    
    // TODO: preserve the rs symmetry if presented
    // simple batch over s here 
    size_t nsMax = ns;
    halfTMOTPISize = nAO * nAO * nr * nsMax;
    nBatch = CQMemManager::get().max_avail_allocatable<MatsT>(halfTMOTPISize, 1);
    
    while (nBatch < 1) {
      
      nsMax /= 2;
      
      std::cout << "nsMax = " << nsMax << std::endl;
      
      if (nsMax < 1) {
        std::cout << "Memory not enough to batch over last index" << std::endl;
        double mem_avail  = CQMemManager::get().template max_avail_allocatable<double>(1,halfTMOTPISize)
                           * sizeof(double) / 1e9; 
        double mem_needed = nAO * nAO * nr * sizeof(MatsT) / 1e9;
            
        std::cout << "  - Memory available before batch integral transformation: " 
                  << mem_avail << " GB" << std::endl; 
        std::cout << "  - Memory needed at least for nsMax = 1: " 
                  << mem_needed << " GB" << std::endl; 
        
        throw std::bad_alloc();
      }
      
      halfTMOTPISize = nAO * nAO * nr * nsMax;
      nBatch = CQMemManager::get().max_avail_allocatable<MatsT>(halfTMOTPISize, 1);
    }
    
    // balance batching inside the work loop
    if (nsMax > 1) {
      size_t nMatsTAvail = CQMemManager::get().max_avail_allocatable<MatsT>(1,halfTMOTPISize*2);
      // with extra memory as 20%
      size_t fockGDSCRSize = ss_.fockBuilder->formRawGDSCRSizePerBatch(ss_, false, false) * 1.2; 
      
      // try different nsMax and use the smallest total batch
      size_t nrnsMax = nr * nsMax; 
      size_t maxNBatch = (nMatsTAvail - halfTMOTPISize)/fockGDSCRSize; 
      size_t minTotalBatch = ((nrnsMax - 1)/maxNBatch + 1) * nBatch;
      size_t bestnsMax = nsMax;
      size_t curTotalBatch = minTotalBatch;

      for (auto i = 1ul; i < nsMax; i++) {
        nrnsMax = nr * i;
        nBatch  = (ns - 1) / i + 1;
        halfTMOTPISize = nAO * nAO * nr * i;
        maxNBatch = (nMatsTAvail - halfTMOTPISize)/fockGDSCRSize;
        curTotalBatch = ((nrnsMax - 1)/maxNBatch + 1) * nBatch;
        if (curTotalBatch < minTotalBatch) {
          minTotalBatch = curTotalBatch;
          bestnsMax = i;
        }
      }
      
      // update nsMax
      nsMax = bestnsMax;
      halfTMOTPISize = nAO * nAO * nr * nsMax;
      nBatch  = (ns - 1) / nsMax + 1;
    }

    std::cout << "* Batch over last indice s: with maxNs =  " << nsMax << std::endl;

    size_t npqr = np * nq * nr; 

    for (auto si = 0ul, sioff = soff, iBatch = 1ul; si < ns; iBatch++) {
      
      size_t nsi = std::min(nsMax, ns - si); 
      
      std::cout << "   - Batch " << iBatch << ": s range from " 
                << std::setw(5) << sioff + 1 << " ~ " 
                << std::setw(5) << sioff + nsi <<std::endl;  

      std::vector<std::pair<size_t,size_t>> batch_off_sizes = off_sizes; 
      batch_off_sizes[3].first = sioff;
      batch_off_sizes[3].second = nsi;

      performSubsetTransformTPISSFockN6(pert, batch_off_sizes, MOTPI + npqr * si,
        moType, false, delta);
      
      si += nsi;
      sioff += nsi;
    }

    return; 
  
  } // subsetTransformTPISSFockN6

  /**
   *  \brief transform AO TPI to form MO TPI
   *  thru SingleSlater formfock 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::performSubsetTransformTPISSFockN6(EMPerturbation & pert, 
    const std::vector<std::pair<size_t,size_t>> &off_sizes, MatsT* MOTPI, 
    const std::string & moType, bool cacheIntermediates,
    TPI_TRANS_DELTA_TYPE delta) {

    size_t poff = off_sizes[0].first;
    size_t qoff = off_sizes[1].first;
    size_t roff = off_sizes[2].first;
    size_t soff = off_sizes[3].first;
    size_t np = off_sizes[0].second;
    size_t nq = off_sizes[1].second;
    size_t nr = off_sizes[2].second;
    size_t ns = off_sizes[3].second;

    bool pqSymm = (poff == qoff) and (np == nq); 
    std::string moType_cache = "HalfTMOTPI-";
    moType_cache += getUniqueSymbol(moType[0]);
    moType_cache += getUniqueSymbol(moType[1]);
    std::string rs_moType_cache = "HalfTMOTPI-";
    rs_moType_cache += getUniqueSymbol(moType[2]);
    rs_moType_cache += getUniqueSymbol(moType[3]);
    
    if (delta == KRONECKER_DELTA_PQ or delta == KRONECKER_DELTA_PQ_RS) 
      moType_cache += "-delta";
    if (delta == KRONECKER_DELTA_RS or delta == KRONECKER_DELTA_PQ_RS) 
      rs_moType_cache += "-delta";

    // check pq and rs can be swapped to accelarate the computation 
    bool rsSymm = (roff == soff) and (nr == ns); 
    bool swap_pq_rs = false;
    auto pqCache = ints_cache_.getIntegral<InCoreRITPI, MatsT>(moType_cache); 
    auto rsCache = ints_cache_.getIntegral<InCoreRITPI, MatsT>(rs_moType_cache); 
    
    if (not pqCache and rsCache)     swap_pq_rs = true;
    if (delta == KRONECKER_DELTA_RS) swap_pq_rs = true;
    
    //if (rsSymm and not pqSymm)   swap_pq_rs = true; 
    
    // swap pq rs when nrs_batch < npq_batch because first half is N6 and second is N5
    if (not swap_pq_rs) {
      
      size_t npq_batch = pqSymm ? np * (np + 1) / 2: np*nq; 
      size_t nrs_batch = rsSymm ? nr * (nr + 1) / 2: nr*ns; 
      
      if (delta == KRONECKER_DELTA_PQ or delta == KRONECKER_DELTA_PQ_RS) 
        npq_batch = np;
      if (delta == KRONECKER_DELTA_RS or delta == KRONECKER_DELTA_PQ_RS) 
        nrs_batch = nr;

      if (npq_batch > nrs_batch) swap_pq_rs = true;  
    }
    
    if (swap_pq_rs) {
      std::swap(poff, roff); 
      std::swap(qoff, soff); 
      std::swap(np, nr); 
      std::swap(nq, ns); 
      std::swap(pqSymm, rsSymm);
      moType_cache = rs_moType_cache; 
      
      if      (delta == KRONECKER_DELTA_PQ) delta = KRONECKER_DELTA_RS;
      else if (delta == KRONECKER_DELTA_RS) delta = KRONECKER_DELTA_PQ;
      else if (delta == KRONECKER_DELTA_PS) delta = KRONECKER_DELTA_RQ;
      else if (delta == KRONECKER_DELTA_RQ) delta = KRONECKER_DELTA_PS;
    }
    
#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
    std::cout << "moType = " << moType   << ", kronecker_delta = " << delta 
              << ", swap_pq_rs = " << swap_pq_rs << std::endl;
    std::cout << "moType_cache = " << std::setw(20) << moType_cache;
#endif

    size_t npr  = np * nr;
    size_t npq  = np * nq;
    size_t npqr = npq * nr;
    size_t npqDim = (delta == KRONECKER_DELTA_PQ or 
                     delta == KRONECKER_DELTA_PQ_RS) ? np: npq;

    size_t nAO  = ss_.nAlphaOrbital() * ss_.nC;
    size_t NB   = ss_.nAlphaOrbital();
    size_t nAO2 = nAO * nAO;
    
    MatsT * dummy_ptr = nullptr;
    std::shared_ptr<InCoreRITPI<MatsT>> halfTMOTPI = nullptr; 
    
    if (cacheIntermediates) {
      halfTMOTPI = ints_cache_.getIntegral<InCoreRITPI, MatsT>(moType_cache);
#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
      if (halfTMOTPI) std::cout << "----Find cache!!!" << std::endl;
#endif    
    } 
    
    //
    // 1/2 transformation to obtain halfTMOTPI(mu, nu, p, q)
    //
    if (not halfTMOTPI) {
      
#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
      std::cout << "----Not find cache, do transformation" << std::endl;
#endif

      ALLOCATE_AND_CLEAR_CACHE_IF_NECESSARY(
        halfTMOTPI = std::make_shared<InCoreRITPI<MatsT>>(nAO, npqDim);
      )

      cqmatrix::Matrix<MatsT> SCR(nAO);
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> pq1PDMs, pqMOTPIs, dummy;
      MatsT *MOTPIpq_ptr = nullptr, *density_ptr = nullptr; 
      bool is4C = ss_.nC == 4;
      bool is2C = ss_.nC == 2;
      bool is1C = ss_.nC == 1;
      size_t pqSCRSize = is4C ? 2*NB: NB;
      
      // find out maximum batch size base on current memory limit
      
      // SCR size form fockbuilder
      size_t fockGDSCRSize = ss_.fockBuilder->formRawGDSCRSizePerBatch(ss_, false, false); 
      
      // SCR size for spinor density and half-transformed integrals
      fockGDSCRSize += is4C ? pqSCRSize*pqSCRSize: 2*pqSCRSize*pqSCRSize;  
      size_t maxNBatch = 0;
      double allow_extra = 0.2;
      
      ALLOCATE_AND_CLEAR_CACHE_IF_NECESSARY(
        maxNBatch = CQMemManager::get().max_avail_allocatable<MatsT>(size_t(fockGDSCRSize*(1.+allow_extra)), 1);
        if(maxNBatch == 0) 
          CErr(" Memory is not enough for 1 density in subsetTransformTPISSFockN6");
      )

      std::vector<std::pair<size_t,size_t>> pqJobs;
      
      if (delta == KRONECKER_DELTA_PQ or delta == KRONECKER_DELTA_PQ_RS) {
        for (auto p = 0ul; p <  np; p++) pqJobs.push_back({p,p});
      } else {
        for (auto q = 0ul; q <  nq; q++) 
        for (auto p = 0ul; p <  np; p++) {
          pqJobs.push_back({p,q});
          if (pqSymm and p == q) break; 
        }
      }

      size_t NJob = pqJobs.size();
      size_t NJobComplete  = 0ul;
      size_t NJobToDo = 0ul;

#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
      std::cout << "First Half transfromation with maxNBatch = " << maxNBatch << std::endl;
#endif
      
      while (NJobComplete < NJob) { 
        
        size_t p, q;
        NJobToDo = std::min(maxNBatch, NJob-NJobComplete);

#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
        std::cout << "----" << std::endl;
        std::cout << "NJob         = " << NJob << std::endl;
        std::cout << "NJobComplete = " << NJobComplete << std::endl;
        std::cout << "NJobToDo     = " << NJobToDo << std::endl;
#endif
        
        // build fake densities in batch
        for (auto i = 0ul; i < NJobToDo; i++) { 
          
          p = pqJobs[i + NJobComplete].first; 
          q = pqJobs[i + NJobComplete].second; 
          
          pq1PDMs.push_back(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(pqSCRSize, is4C, is4C));  
          
          // 4C will reuse pq1PDM as densities are component scattered anyways
          if (is4C) pqMOTPIs.push_back(pq1PDMs.back());
          else pqMOTPIs.push_back(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(pqSCRSize, is4C, is4C));  
          
          if (is1C) {
            density_ptr = pq1PDMs.back()->S().pointer();
          } else {
            density_ptr = SCR.pointer();
          }
        
          blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
            nAO, nAO, 1, MatsT(1.), ss_.mo[0].pointer() + (q+qoff)*nAO, nAO,
            ss_.mo[0].pointer() + (p+poff)*nAO, nAO, MatsT(0.), density_ptr, nAO);
          
          if (not is1C) *pq1PDMs.back() = SCR.template spinScatter<MatsT>(is4C, is4C); 
        
// DEBUG DENSITY***************************        
//          cqmatrix::Matrix<MatsT> denSCR = pq1PDMs.back()->S();
//          denSCR += denSCR.T();
//          denSCR.output(std::cout, "p = " + std::to_string(p) + ", q = " + std::to_string(q), true);  
// DEBUG DENSITY***************************    
        }
        
        // do contraction to get half-transformed integrals
        if (is1C) ss_.fockBuilder->formRawGDInBatches(ss_, pert, false, 0., false, pq1PDMs, pqMOTPIs, dummy, dummy);
        else      ss_.fockBuilder->formRawGDInBatches(ss_, pert, false, 0., false, pq1PDMs, dummy, dummy, pqMOTPIs);
        
        // copy over to SCR
        for (auto i = 0ul; i < NJobToDo; i++) { 
           
          if (not is1C) {
            SCR = pqMOTPIs[i]->template spinGather<MatsT>();
            MOTPIpq_ptr = SCR.pointer(); 
          } else MOTPIpq_ptr = pqMOTPIs[i]->S().pointer();
          
          p = pqJobs[i + NJobComplete].first; 
          
          if (delta == KRONECKER_DELTA_PQ or delta == KRONECKER_DELTA_PQ_RS) {
            SetMat('N', nAO, nAO, MatsT(1.), MOTPIpq_ptr, nAO, halfTMOTPI->pointer() + p*nAO2, nAO);
          } else {  
            q = pqJobs[i + NJobComplete].second; 
            SetMat('N', nAO, nAO, MatsT(1.), MOTPIpq_ptr, nAO, halfTMOTPI->pointer() + (p + q*np)*nAO2, nAO);
            if (pqSymm and p < q)  
              SetMat('C', nAO, nAO, MatsT(1.), MOTPIpq_ptr, nAO, halfTMOTPI->pointer() + (q + p*np)*nAO2, nAO); 
          }
        }

        // increment and clear SCR after job done
        NJobComplete += NJobToDo; 
        pq1PDMs.clear();
        pqMOTPIs.clear();
      
      } // main loop
      
      if (cacheIntermediates) ints_cache_.addIntegral(moType_cache, halfTMOTPI);
     
    } //  if there is no cache
    
    // allocate intermediate SCR
    MatsT * SCR = nullptr;
    size_t SCRSize;
    if      (delta == NO_KRONECKER_DELTA)    SCRSize = nAO*npqr;  
    else if (delta == KRONECKER_DELTA_PQ)    SCRSize = nAO*np*nr;  
    else if (delta == KRONECKER_DELTA_RS)    SCRSize = nAO*npq;  
    else if (delta == KRONECKER_DELTA_PS)    SCRSize = nAO*npqr;  
    else if (delta == KRONECKER_DELTA_RQ)    SCRSize = nAO*npr;  
    else if (delta == KRONECKER_DELTA_PQ_RS) SCRSize = nAO*np;  
    else if (delta == KRONECKER_DELTA_PS_RQ) SCRSize = nAO*npr;  
    
    SCR = CQMemManager::get().malloc<MatsT>(SCRSize);
    
    if (delta == NO_KRONECKER_DELTA or delta == KRONECKER_DELTA_PQ) {
      
      char TransMOTPI = swap_pq_rs ? 'N': 'T';
      size_t ndim = delta == NO_KRONECKER_DELTA ? npq: np;     

      PairTransformation('N', ss_.mo[0].pointer(), nAO, roff, soff,
        'N', halfTMOTPI->pointer(), nAO, nAO, ndim, 
        TransMOTPI, MOTPI, nr, ns, dummy_ptr, SCR, false); 
    
    } else if (delta == KRONECKER_DELTA_RS or delta == KRONECKER_DELTA_PQ_RS) {

      for (auto r = 0ul; r < nr; r++) {

        auto rMO = ss_.mo[0].pointer() + (r + roff) * nAO;      
        
        // SCR(nu p q, r) = (mu, nu p q)^H MO(mu, r) 
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
          nAO*npqDim, 1, nAO, MatsT(1.), halfTMOTPI->pointer(), nAO, 
          rMO, nAO, MatsT(0.), SCR, nAO*npqDim);
      
        // MOTPI(p q, r) = SCR(nu, p q)^H MO(nu, r) 
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
          npqDim, 1, nAO, MatsT(1.), SCR, nAO, 
          rMO, nAO, MatsT(0.), MOTPI + r*npqDim, npqDim);
      
      }

      if(swap_pq_rs) IMatCopy('T', npqDim, nr, MatsT(1.), MOTPI, npqDim, nr); 

    } else if (delta == KRONECKER_DELTA_PS) {
        
      // SCR(nu p q, r) = (mu, nu p q)^H MO(mu, r) 
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
        nAO*npq, nr, nAO, MatsT(1.), halfTMOTPI->pointer(), nAO, 
        ss_.mo[0].pointer() + roff*nAO, nAO, MatsT(0.), SCR, nAO*npq);
    
      MatsT * SCR2 = nullptr;
      SCR2 = CQMemManager::get().malloc<MatsT>(nAO*nq*nr);
      
      size_t nqr = nq*nr;

      for (auto p = 0ul; p < np; p++) {
        
        // Gather SCR2(nu, q, r) from SCR(nu p q, r)
        #pragma omp parallel for
        for (auto r = 0ul; r < nr; r++) 
        for (auto q = 0ul; q < nq; q++) {
          std::copy_n(SCR + p*nAO + q*nAO*np + r*nAO*npq,
            nAO, SCR2 + q*nAO + r*nAO*nq);
        }

        // MOTPI(q r, p) = SCR2(nu, q r)^H MO(nu, p) 
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
          nqr, 1, nAO, MatsT(1.), SCR2, nAO,
          ss_.mo[0].pointer() + (p+poff)*nAO, nAO, 
          MatsT(0.), MOTPI + p*nqr, nqr);  
      }
      
      // MOTPI(q r, p) -> MOTPI(p,q,r)
      if (not swap_pq_rs) IMatCopy('T', nqr, np, MatsT(1.), MOTPI, nqr, np);
      // MOTPI(q r, p) -> MOTPI(r,p,q)
      else                IMatCopy('T', nq, nr*np, MatsT(1.), MOTPI, nq, nr*np);  

      if (SCR2) CQMemManager::get().free(SCR2);
    
    } else if (delta == KRONECKER_DELTA_RQ or delta == KRONECKER_DELTA_PS_RQ) {
      
      size_t nAOp = nAO * np;
      for (auto r = 0ul; r < nr; r++) {
        // SCR(nu p, r) = (mu, nu p)^H MO(mu, r) 
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
          nAOp, 1, nAO, MatsT(1.), halfTMOTPI->pointer() + r*nAO2*np, nAO, 
          ss_.mo[0].pointer() + (r+roff)*nAO, nAO, MatsT(0.), SCR+r*nAOp, nAOp);
      }
      
      if (delta == KRONECKER_DELTA_RQ) {
        // MOTPI(p,r,s) = SCR(nu p, r)^H MO(nu, s)
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
          npr, ns, nAO, MatsT(1.), SCR, nAO, 
          ss_.mo[0].pointer() + soff*nAO, nAO, MatsT(0.), MOTPI, npr);
        
        // MOTPI(p,r,s) -> MOTPI(r,s,p)
        if (swap_pq_rs) IMatCopy('T', np, ns*nr, MatsT(1.), MOTPI, np, ns*nr);
       
       } else if (delta == KRONECKER_DELTA_PS_RQ) {
       
         MatsT * SCR2 = nullptr;
         SCR2 = CQMemManager::get().malloc<MatsT>(nAO*nr);
         
         for (auto p = 0ul; p < np; p++) {
           
           // Gather SCR2(nu, r) from SCR(nu, p, r)
           #pragma omp parallel for
           for (auto r = 0ul; r < nr; r++) { 
             std::copy_n(SCR + p*nAO + r*nAO*np, nAO, SCR2 + r*nAO);
           }
           
           // MOTPI(r,p) = SCR2(nu, r)^H MO(nu, p)
           blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
             nr, 1, nAO, MatsT(1.), SCR2, nAO, 
             ss_.mo[0].pointer() + (p+poff)*nAO, nAO, MatsT(0.), MOTPI + p*nr, nr);
           
         }
         
         // MOTPI(r,p) -> MOTPI(p,r)
         if (not swap_pq_rs) IMatCopy('T', nr, np, MatsT(1.), MOTPI, nr, np);
         
         if (SCR2) CQMemManager::get().free(SCR2);
       }

    } // delta 
    
    // free memories
    halfTMOTPI = nullptr;
    if (SCR) CQMemManager::get().free(SCR);
 
  }; // MOIntsTransformer::perfromSubsetTransformTPISSFockN6
  

}; // namespace ChronusQ
