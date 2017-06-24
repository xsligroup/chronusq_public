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
#include <cibuilder/rasci.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cibuilder/rasci/sigma.hpp>

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

//#define _DEBUG_CIBUILDER_RASCI_IMPL
//#define _DEBUG_CIBUILDER_RASCI_BATCH

#define RASCI_LOOP_INIT() \
  size_t nC = mcwfn.reference().nC; \
  auto RASString = std::dynamic_pointer_cast< \
    RASStringManager>(mcwfn.detStr); \
  auto RASStringBeta = (nC == 1) ? std::dynamic_pointer_cast< \
    RASStringManager>(mcwfn.detStrBeta) : nullptr; \
  std::vector<size_t> nActO = mcwfn.MOPartition.nActOs; \
  std::vector<int> fCat = mcwfn.MOPartition.fCat; \
  size_t nCorrE = RASString->nElectron(); \
  size_t nStr = RASString->nString(); \
  size_t nStr_b = (nC == 1) ? RASStringBeta->nString() : 1; \
  size_t nCat = RASString->nCategory(); \
  size_t mxHole = RASString->maxHole(); \
  size_t mxElec = RASString->maxElectron(); \
  std::vector<std::vector<size_t>> LCat = RASString->LCategory(); \
  auto & hCoreP = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("hCoreP_Correlated_Space")); \
  auto & moERI  = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));


namespace ChronusQ {

  typedef std::vector<std::vector<int>> int_matrix; 

  /*
   * Build full RASCI Hamiltonian 
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::buildFullH(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * fullH) {

    CErr("Full Hamiltonian build for RASCI is not implemented.");

  } // RASCI::buildFullH


  /*
   * Build Diagonal RASCI Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::buildDiagH(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * diagH) {

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    std::cout << "LL RAS diagonal Hamilton build." << std::endl;
#endif
    auto RASdiagHSt = tick();

    RASCI_LOOP_INIT(); // check top for variable definitions
    auto & hCore = *(mcwfn.moints->template getIntegral<OnePInts, MatsT>("hCore_Correlated_Space"));

    // empty CI Diagonal Hamiltonian
    std::fill_n(diagH, nStr, MatsT(0.));
    std::vector<size_t> j(3, 0);
    std::vector<size_t> nEras(3, 0);
    std::vector<size_t> nDras(3, 0);
    std::vector<size_t> elecPos(nCorrE, 0);
    int_matrix empty = {{}};

    // Alpha Part for 1C or the whole build for 2C and 4C
    size_t iCatJ = 0;
    for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
    for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {

      nEras[0] = nActO[0] - iHJ;
      nEras[2] = iEJ;
      if (nEras[0] + nEras[2] > nCorrE or nActO[1] + nEras[0] + nEras[2] < nCorrE)
        continue;

      nEras[1] = nCorrE - nEras[0] - nEras[2];
      nDras[0] = LCat[iCatJ][1];
      nDras[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
      nDras[2] = LCat[iCatJ][3]/LCat[iCatJ][2];

      MatsT* diagH_J = diagH + fCat[iCatJ];

      std::vector<int_matrix> addr_array;
      std::vector<int_matrix> de_addr_array;
      for (auto i = 0ul; i < 3; i++){
        if (nActO[i] == 0) {
          addr_array.push_back(empty);
          de_addr_array.push_back(empty);
        } else {
        addr_array.push_back(RASString->buildAddressingArray(nEras[i], nActO[i]));
        de_addr_array.push_back(RASString->buildDeAddressingArray(addr_array[i]));
        }
      }

      for (j[2] = 0; j[2] < nDras[2]; j[2]++)
      for (j[1] = 0; j[1] < nDras[1]; j[1]++)
      for (j[0] = 0; j[0] < nDras[0]; j[0]++) {
        size_t j_addr = j[0] * LCat[iCatJ][0] + j[1] * LCat[iCatJ][1] + j[2] * LCat[iCatJ][2];
        size_t orb_offset = 0;
        size_t e_offset = 0;
        for (auto i = 0ul; i < 3; i++){
          if (nActO[i] == 0) continue;
          std::vector<size_t> eP = RASString->address2DetString(j[i], addr_array[i],
                                            de_addr_array[i]);
          for (auto pos = 0ul; pos < nEras[i]; pos++)
            elecPos[pos+e_offset] = eP[pos] + orb_offset;
          orb_offset += nActO[i];
          e_offset += nEras[i];
        }
        // Update 1e part
        // hCore_ii
        for (auto pos_i = 0ul; pos_i < nCorrE; pos_i++) {
          diagH_J[j_addr] += hCore(elecPos[pos_i], elecPos[pos_i]);
          // Update 2e part
          // moERI: (ii|jj) or -(ij|ji)
          for (auto pos_j = 0ul; pos_j < nCorrE; pos_j++) {
            diagH_J[j_addr] += MatsT(0.5) * moERI(elecPos[pos_i], elecPos[pos_i],
                                        elecPos[pos_j], elecPos[pos_j]);
            diagH_J[j_addr] -= MatsT(0.5) * moERI(elecPos[pos_i], elecPos[pos_j],
                                            elecPos[pos_j], elecPos[pos_i]);
          }
        }
      }
    }

    if (nC != 1) {
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
      prettyPrintSmart(std::cout,"LL RASCI diag Hamiltonian", diagH, nStr, 1, nStr);
#endif

      double RASdiagHdur = tock(RASdiagHSt);
      std::cout << "\nRASDiagH - DURATION = " << std::setprecision(8)
                << RASdiagHdur << " s." << std::endl;

    } else CErr("Only 2C RASCI diag Hamiltonian is implemented");

  } // RASCI::buildDiagH

  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::buildMu(MCWaveFunction<MatsT, IntsT> & mcwfn,
    size_t nVec, MatsT * C, MatsT * Mu, EMPerturbation & pert){
     CErr("RAS Mu operator NYI");
   }


  /*
   * Build RASCI Sigma, Matrix-vector product
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::buildSigma(MCWaveFunction<MatsT, IntsT> & mcwfn,
    size_t nVec, MatsT * C, MatsT * Sigma) {

//    auto RASSigmaTSt = tick();

    RASCI_LOOP_INIT(); // check top for variable definitions

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    std::cout<<"Sigma build starts"<<std::endl;
    prettyPrintSmart(std::cout,"LL RASCI Sigma Build -- C", C, nStr, nVec, nStr);
#endif

    // Allocate SCR for D and G subblocks    
    size_t maxdim = 0; // max dimension of D and G subblocks
    size_t nNZatMD = 0; // number of nonzero at maxdim
    size_t nDetatMD = 0; // number of determinants at maxdim
    size_t maxnNZ = 0; // max number of nonzero elements
    size_t maxdSCR = 0; // max dimension of D and G scratch array
    size_t mindSCR = 0; // min dimension of D and G scratch array (when nNonzero = 1)
    size_t maxDetCat = 0; // max number of determinants in one catergory

    std::vector<size_t> nDetL(3, 0);
    std::vector<size_t> nDetJ(3, 0);
    std::vector<size_t> nDetK(3, 0);

    for (auto p = 0ul; p < 3; p++)
    for (auto q = 0ul; q < 3; q++)
    for (auto iCatL = 0ul; iCatL < nCat; iCatL++)
    for (auto iCatK = 0ul; iCatK < nCat; iCatK++) {
      if (auto exList = RASString->find_exlist(p, q, iCatL, iCatK)) {
        size_t nNonzero = exList->nNonZero();
        if (nNonzero * LCat[iCatL][3] > maxdim) {
          maxdim = nNonzero * LCat[iCatL][3];
          nNZatMD = nNonzero;
          nDetatMD = LCat[iCatL][3];
        }
        if (nNonzero > maxnNZ) maxnNZ = nNonzero;

        nDetL[0] = LCat[iCatL][1];
        nDetL[1] = LCat[iCatL][2]/LCat[iCatL][1];
        nDetL[2] = LCat[iCatL][3]/LCat[iCatL][2];

        size_t nDetSub = (p == q) ? LCat[iCatL][3]/nDetL[p]
              : LCat[iCatL][3]/(nDetL[p]*nDetL[q]);
        if (nNonzero * nDetSub > maxdSCR) maxdSCR = nNonzero * nDetSub;
        if (nDetSub > mindSCR) mindSCR = nDetSub;
      }
    }    
//    std::cout<<"maxdim: "<<maxdim<<std::endl;
//    std::cout<<"nNZatMD: "<<nNZatMD<<", nDetatMD: "<<nDetatMD<<std::endl;
//    std::cout<<"maxdSCR: "<<maxdSCR<<", mindSCR: "<<mindSCR<<std::endl;

    for (auto iCatL = 0ul; iCatL < nCat; iCatL++)
    if (LCat[iCatL][3] > maxDetCat) maxDetCat = LCat[iCatL][3];
//    std::cout<<"maxDetCat: "<<maxDetCat<<std::endl;

    bool DoBatch = false;
    MatsT* Ci = C;
    MatsT* HC = Sigma;
    size_t maxdimD = maxdim;
    size_t maxdimG = maxdim;
    size_t ndSCR = maxdSCR;
    size_t nrSCR= int((maxdim-1)/ndSCR) + 1; // ratio of maxdim over maxdSCR
    size_t maxAvlMem, nNZAvl, nBatch;
    int nNZG = 0;
    int nNZD = 0;

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
//    size_t totalSCR = (maxnNZ*maxnNZ + maxdSCR + maxdSCR) * nthreads;

    /*-------------------------------------*/
    /*  Determine whether there is enough  */
    /*  memory to store the maxdim block D */
    /*  and G. If not, do batching.        */
    /*-------------------------------------*/

//    int check_mem = 0;
    MatsT* DBlk = nullptr;
    MatsT* GBlk = nullptr; 
    MatsT* ERIscr = nullptr;
    MatsT* Dscr = nullptr;
    MatsT* Gscr = nullptr;

//    std::cout << mem << std::endl;
    try{
      ERIscr = CQMemManager::get().malloc<MatsT>(maxnNZ*maxnNZ*nthreads);
      Dscr = CQMemManager::get().malloc<MatsT>(ndSCR*nthreads);
      Gscr = CQMemManager::get().malloc<MatsT>(ndSCR*nthreads);
      DBlk = CQMemManager::get().malloc<MatsT>(maxdimD);
      GBlk = CQMemManager::get().malloc<MatsT>(maxdimG);
    } catch (std::bad_alloc& ba) {
      if (not ERIscr) CErr ("Not enough memory for RAS sigma.");
      std::cout<<"---- modifying batching ----"<<std::endl;
      if (Dscr) CQMemManager::get().free(Dscr);
      if (Gscr) CQMemManager::get().free(Gscr);
      if (DBlk) CQMemManager::get().free(DBlk);
      if (GBlk) CQMemManager::get().free(GBlk);

//      CQMemManager::get().print_free();
      size_t totalSCR = ndSCR * nthreads + maxdimD + maxdimG;
      maxAvlMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, totalSCR);
      // to determine ndSCR
      if ((maxAvlMem / (ndSCR*2)) < (nrSCR + nthreads))
        ndSCR = std::max(maxAvlMem / ((nrSCR + nthreads)*2), mindSCR);
      Dscr = CQMemManager::get().malloc<MatsT>(ndSCR*nthreads);
      Gscr = CQMemManager::get().malloc<MatsT>(ndSCR*nthreads);
      // to determine maxdimD and maxdimG
      totalSCR = maxdimD + maxdimG;
      maxAvlMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, totalSCR);
//      std::cout<<"maxAvlMem: "<<maxAvlMem<<std::endl;
      if (maxAvlMem < 2 * maxDetCat) CErr("Not enough memory for RAS sigma.");
      nNZAvl = maxAvlMem / nDetatMD;
      nBatch = nNZatMD / nNZAvl;
      while (nNZG < 1 or (2*nNZG) < nNZD) {
        nBatch += 1;
        nNZD = nNZatMD / nBatch + 1;
        nNZG = nNZAvl - nNZD;
      }
//      std::cout<<"nNZD: "<<nNZD<<", nNZG: "<<nNZG<<std::endl;
      maxdimG = nDetatMD * nNZG;
      GBlk = CQMemManager::get().malloc<MatsT>(maxdimG);
      maxdimD = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxdimD);
      if (maxdimD < maxDetCat) CErr("Not enough memory for RAS sigma.");
      DBlk = CQMemManager::get().malloc<MatsT>(maxdimD);;
//      std::cout<<"maxdimG: "<<maxdimG<<", maxdimD: "<<maxdimD<<std::endl;
      DoBatch = true;
//      std::cout<<"Batching done!"<<std::endl;
    }

    if (DoBatch) std::cout<<"Doing batch"<<std::endl;

    // Compute nBatch
    std::vector<size_t> nBatchrs(9*nCat, 1);
    std::vector<size_t> nBatchpq(9*nCat, 1);
    std::shared_ptr<const ExcitationList> rasList;
    size_t nNZL, ndCat, ndim;
    if (DoBatch) {
      for (auto I = 0ul; I < 3; I++)
      for (auto J = 0ul; J < 3; J++)
      for (auto L = 0ul; L < nCat; L++) {
        for (auto iCat = 0ul; iCat < nCat; iCat++) {
          rasList = RASString->find_exlist(I, J, L, iCat);
          if (rasList) break;
        }
        if (rasList == nullptr) continue;
        nNZL = rasList->nNonZero();
        ndCat = LCat[L][3];
        ndim = nNZL * ndCat;
        nDetL[0] = LCat[L][1];
        nDetL[1] = LCat[L][2]/LCat[L][1];
        nDetL[2] = LCat[L][3]/LCat[L][2];
        size_t nDetSub = (I == J) ? ndCat/nDetL[I] : ndCat/(nDetL[I] * nDetL[J]);

        nNZAvl = ndSCR / nDetSub;
        nBatchrs[L+J*nCat+I*3*nCat] = 1 + int((nNZL - 1) / nNZAvl);
        nBatchpq[L+J*nCat+I*3*nCat] = 1 + int((nNZL - 1) / nNZAvl);
        if (ndim > maxdimD) {
          nNZAvl = maxdimD / ndCat;
          size_t prenBrs = nBatchrs[L+J*nCat+I*3*nCat];
          nBatchrs[L+J*nCat+I*3*nCat] = std::max(1 + ((nNZL - 1) / nNZAvl), prenBrs);
        }
        if (ndim > maxdimG) {
          nNZAvl = maxdimG / ndCat;;
          size_t prenBpq = nBatchpq[L+J*nCat+I*3*nCat];
          nBatchpq[L+J*nCat+I*3*nCat] = std::max(1 + ((nNZL - 1) / nNZAvl), prenBpq);
        }  
      }
    }


    // Initalize Sigma
    std::fill_n(Sigma, nStr*nVec, MatsT(0.));

    size_t nNZrs, iCatJ, iCatL, nrsBatch, nNZrsof, nNZrsB, nNZpq, npqBatch, nNZpqB, nNZpqof;
    double symmFc;
    std::shared_ptr<const ExcitationList> LrsList, LpqList;

    auto nLAThreads = GetLAThreads();
    SetLAThreads(1);
    for (auto iVec = 0ul; iVec < nVec; iVec++, HC+=nStr, Ci+=nStr) {

//    auto RASSigmaVSt = tick();

    // Alpha part for 1C or the whole build for 2C and 4C
    for (auto rBlk = 0ul; rBlk < 3; rBlk++)
    for (auto sBlk = 0ul; sBlk < 3; sBlk++) {
      if (nActO[rBlk] == 0 or nActO[sBlk] == 0) continue;
      iCatJ = 0;
      for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
      for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {
        if (LCat[iCatJ][3] == 0) continue;
        for (iCatL = 0; iCatL < nCat; iCatL++) {
          LrsList = RASString->find_exlist(sBlk, rBlk, iCatJ, iCatL);
          if (LrsList) break;
        }
        if (fCat[iCatL] < 0 or LrsList == nullptr) continue;
        nDetJ[0] = LCat[iCatJ][1];
        nDetJ[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
        nDetJ[2] = LCat[iCatJ][3]/LCat[iCatJ][2];
        nNZrs = LrsList->nNonZero();
        MatsT * Ci_L = Ci + fCat[iCatL];

        // Batch rs block
        nrsBatch = nBatchrs[iCatJ + rBlk*nCat + sBlk*3*nCat];
        nNZrsB = 1 + ((nNZrs-1)/nrsBatch);
        for (auto irsB = 0ul; irsB < nrsBatch; irsB++) {
          nNZrsof = irsB * nNZrsB;
          if (nNZrsof >= nNZrs) break;
          else {
            nNZrsB = std::min(nNZrsB, nNZrs-nNZrsof);
//            std::cout<<std::endl;
//            std::cout<<"irsBatch: "<<irsB<<", nNZrsBatch: "<<nNZrsB<<std::endl;


            /*-------------------------------------*/
            /* Build Block D as [nNZrs, J]         */
            /*   Use excitation list as<J|E_sr|L>  */
            /*-------------------------------------*/
            //auto RASSigma1St = tick();
            std::fill_n(DBlk, nNZrsB*LCat[iCatJ][3], MatsT(0.));
            if (rBlk == sBlk) // diagonal
              SigmaD_diag(Ci_L, DBlk, nNZrsof, nNZrsB, LrsList, rBlk, LCat[iCatJ], nDetJ);
            else // off-diagonal
              SigmaD_offdiag(Ci_L, DBlk, nNZrsof, nNZrsB, LrsList, rBlk, sBlk,
                            LCat[iCatJ], LCat[iCatL], nDetJ);
            //double RASSigma1dur = tock(RASSigma1St);
            //std::cout << "\nRAS Sigma Step 1 - DURATION = " << std::setprecision(8)
            //            << RASSigma1dur << " s." << std::endl;

            /*-------------------------------------*/
            /* Update Sigma_1e[K] =                */
            /*       \sum_pq hCorP_pq  * D[pq, K]; */
            /*-------------------------------------*/
            //auto RASSigma2St = tick();
            MatsT * HC_J = HC + fCat[iCatJ];
            int thread_id;
            MatsT *ERISCR, *DSCR, *GSCR;
            #pragma omp parallel default(shared) firstprivate(ERISCR, DSCR, GSCR) private(thread_id)
            {
            thread_id = GetThreadID();
            ERISCR = ERIscr + maxnNZ*maxnNZ*thread_id;
            DSCR = Dscr + ndSCR*thread_id;
            GSCR = Gscr + ndSCR*thread_id;
            
            if (rBlk == sBlk) { // diagonal
              int rsoff = 0;
              for (auto i = 0ul; i < rBlk; i++) rsoff += nActO[i];
              Sigma1e_diag(HC_J, GSCR, DBlk, DSCR, hCoreP, ERISCR, nNZrsof,
                            nNZrsB, LrsList, rBlk, LCat[iCatJ], nDetJ, rsoff);
            }
            else { // off-diagonal
              int roff = 0;
              int soff = 0;
              for (auto i = 0ul; i < rBlk; i++) roff += nActO[i];
              for (auto i = 0ul; i < sBlk; i++) soff += nActO[i];
              Sigma1e_offdiag(HC_J, GSCR, DBlk, DSCR, hCoreP, ERISCR, nNZrsof,
                            nNZrsB, LrsList, rBlk, sBlk, LCat[iCatJ], nDetJ, roff, soff);
            }
            } // End omp parallel
            //double RASSigma2dur = tock(RASSigma2St);
            //std::cout << "\nRAS Sigma Step 2 - DURATION = " << std::setprecision(8)
            //            << RASSigma2dur << " s." << std::endl;

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
            prettyPrintSmart(std::cout,"LL Sigma Hamiltonian -- Sigma1e", Sigma, nStr, nVec, nStr);
#endif

            /*-------------------------------------*/
            /* Update Sigma 2e parts               */
            /* 2 steps:                            */
            /*      a. Gemm: for nonzero pq, rs    */
            /*         G(pq, J) = (pq|rs) * D(rs,J)*/
            /*      b. Update: Sigma_2e(K) =       */
            /*         \sum_pq <K|E_pq|J> * G(pq,J)*/
            /*-------------------------------------*/

            for (auto pBlk = 0ul; pBlk < 3; pBlk++)
            for (auto qBlk = 0ul; qBlk < 3; qBlk++) {
              if (nActO[pBlk] == 0 or nActO[qBlk] == 0) continue;
              size_t iCatK;
              for (iCatK = 0; iCatK < nCat; iCatK++) {
                LpqList = RASString->find_exlist(pBlk, qBlk, iCatJ, iCatK);
                if (LpqList) break;
              }
              if (LpqList == nullptr) continue;
              if ((pBlk+3*qBlk) > (rBlk+3*sBlk)) continue;
              if (pBlk != rBlk and qBlk != sBlk and pBlk > rBlk) continue;
              if ((pBlk+3*qBlk) == (rBlk+3*sBlk)) symmFc = 0.5;
              else symmFc = 1.;

              nDetK[0] = LCat[iCatK][1];
              nDetK[1] = LCat[iCatK][2]/LCat[iCatK][1];
              nDetK[2] = LCat[iCatK][3]/LCat[iCatK][2];
              nNZpq = LpqList->nNonZero();
              MatsT * HC_K = HC + fCat[iCatK];

              // Batch pq block
              npqBatch = nBatchpq[iCatJ+qBlk*nCat+pBlk*3*nCat];
              nNZpqB = 1 + (nNZpq - 1) / npqBatch;
              for (auto ipqB = 0ul; ipqB < npqBatch; ipqB++) {
                nNZpqof = ipqB * nNZpqB;
                if (nNZpqof >= nNZpq) break;
                else {
                  nNZpqB = std::min(nNZpqB, nNZpq - nNZpqof);
//                  std::cout<<"r: "<<rBlk<<", s: "<<sBlk<<", p: "<<pBlk<<", q: "<<qBlk<<std::endl;
//                  std::cout<<"ipqBatch: "<<ipqB<<", nNZpqBatch: "<<nNZpqB<<std::endl;


                  // a. G(pq, J) = (pq|rs) * D(rs,J)
                  //auto RASSigma3St = tick();
                  std::fill_n(GBlk, nNZpqB*LCat[iCatJ][3], MatsT(0.));

                  int poff = 0;
                  int qoff = 0;
                  int roff = 0;
                  int soff = 0;
                  for (auto i = 0ul; i < rBlk; i++) roff += nActO[i];
                  for (auto i = 0ul; i < sBlk; i++) soff += nActO[i];
                  for (auto i = 0ul; i < pBlk; i++) poff += nActO[i];
                  for (auto i = 0ul; i < qBlk; i++) qoff += nActO[i];

                  #pragma omp parallel default(shared) firstprivate(ERISCR, DSCR, GSCR) private(thread_id)
                  {
                  thread_id = GetThreadID();
                  ERISCR = ERIscr + maxnNZ*maxnNZ*thread_id;
                  DSCR = Dscr + ndSCR*thread_id;
                  GSCR = Gscr + ndSCR*thread_id;                

                  if (rBlk == sBlk and pBlk == qBlk) // diagonal x diagonal block
                    SigmaG_dd(GBlk, GSCR, DBlk, DSCR, ERISCR, moERI, nNZrsof, nNZrsB,
                        nNZpqof, nNZpqB, LpqList, LrsList, rBlk, pBlk, LCat[iCatJ],
                        nDetJ, roff, poff);
                  else if (rBlk == sBlk) // diagonal x off-diagonal block
                    SigmaG_do(GBlk, GSCR, DBlk, DSCR, ERISCR, moERI, nNZrsof, nNZrsB,
                        nNZpqof, nNZpqB, LpqList, LrsList, rBlk, pBlk, qBlk, LCat[iCatJ],
                        nDetJ, roff, poff, qoff);
                  else if (pBlk == qBlk) // off-diagonal x diagonal block
                    SigmaG_od(GBlk, GSCR, DBlk, DSCR, ERISCR, moERI, nNZrsof, nNZrsB,
                        nNZpqof, nNZpqB, LpqList, LrsList, rBlk, sBlk, pBlk, LCat[iCatJ],
                        nDetJ, roff, soff, poff);
                  else // off-diagonal x off-diagonal block
                    SigmaG_oo(GBlk, GSCR, DBlk, DSCR, ERISCR, moERI, nNZrsof, nNZrsB,
                        nNZpqof, nNZpqB, LpqList, LrsList, rBlk, sBlk, pBlk, qBlk,
                        LCat[iCatJ], nDetJ, roff, soff, poff, qoff);
                  } // End omp parallel
                  //double RASSigma3dur = tock(RASSigma3St);
                  //std::cout << "\nRAS Sigma step 3 - DURATION = " << std::setprecision(8)
                  //              << RASSigma3dur << " s." << std::endl;

                  //auto RASSigma4St = tick();
                  // b. Update: Sigma_2e(K)  = \sum_pq <K|E_pq|J> * G(pq,J)
                  if (pBlk == qBlk) // diagonal block
                    Sigma2e_diag(HC_K, GBlk, nNZpqof, nNZpqB, LpqList, pBlk,
                        LCat[iCatJ], nDetJ, symmFc);
                  else // off-diagonal block
                    Sigma2e_offdiag(HC_K, GBlk, nNZpqof, nNZpqB, LpqList, pBlk,
                        qBlk, LCat[iCatJ], LCat[iCatK], nDetJ, symmFc);
                  //double RASSigma4dur = tock(RASSigma4St);
                  //std::cout << "\nRAS Sigma step 4 - DURATION = " << std::setprecision(8)
                  //   << RASSigma4dur << " s." << std::endl;
                }
              }
            }
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
            prettyPrintSmart(std::cout,"LL Sigma Hamiltonian -- Sigma2e", Sigma, nStr, nVec, nStr);
#endif
          }  
        }
      }
    }
//    double RASSigmaVdur = tock(RASSigmaVSt);
//    std::cout << "\nRAS Sigma single vector - DURATION = " << std::setprecision(8)
//                << RASSigmaVdur << " s." << std::endl;
    }
    CQMemManager::get().free(DBlk, GBlk, ERIscr, Dscr, Gscr);
    SetLAThreads(nLAThreads);
    
//    double RASSigmaTdur = tock(RASSigmaTSt);
//    std::cout << "\nRAS Sigma total - DURATION = " << std::setprecision(8)
//                << RASSigmaTdur << " s." << std::endl;



#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL Sigma Hamiltonian -- Sigma", Sigma, nStr, nVec, nStr); 
#endif
 
  } // void RASCI::buildSigma


  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::computeOneRDM(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * C,
    cqmatrix::Matrix<MatsT> & oneRDM) {

    computeTDM(mcwfn, C, C, oneRDM);

  } // RASCI::computeoneRDM


  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::computeTDM(MCWaveFunction<MatsT, IntsT> & mcwfn, MatsT * Cm,
        MatsT * Cn, cqmatrix::Matrix<MatsT> & TDM) {

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    std::cout << "LL RAS compute TDM." << std::endl;
#endif
    auto RAS1rdmSt = tick();

    RASCI_LOOP_INIT(); // check top for variable definitions
    TDM.clear();

    size_t nThreads = GetNumThreads();
    std::vector<cqmatrix::Matrix<MatsT>> SCR;
    for (auto i = 0ul; i < nThreads; i++) {
      SCR.emplace_back(TDM.dimension());
      SCR.back().clear();
    }

    // 2C or 1C alpha part
    std::shared_ptr<const ExcitationList> LrsList;
    size_t iCatJ, iCatL, r, s, l, l1, l2, j, j1, j2, j3, j23off, iNZ;
    double sign;
    std::vector<size_t> nDetJ(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> LCatAlt2(3, 0);
    std::vector<size_t> nDetAlt(3, 0);

    for (auto rBlk = 0ul; rBlk < 3; rBlk++)
    for (auto sBlk = 0ul; sBlk < 3; sBlk++) {
      if (nActO[rBlk] == 0 or nActO[sBlk] == 0) continue;
      iCatJ = 0;
      for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
      for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {
        if (LCat[iCatJ][3] == 0) continue;
        for (iCatL = 0; iCatL < nCat; iCatL++) {
          LrsList = RASString->find_exlist(rBlk, sBlk, iCatJ, iCatL);
          if (LrsList) break;
        }
        if (fCat[iCatL] < 0 or LrsList == nullptr) continue;
        nDetJ[0] = LCat[iCatJ][1];
        nDetJ[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
        nDetJ[2] = LCat[iCatJ][3]/LCat[iCatJ][2];

        if (rBlk == sBlk) {// diagonal
          reassembleCatDet1RAS(LCat[iCatJ], nDetJ, rBlk, LCatAlt, nDetAlt);
          int rsoff = 0;
          for (auto i = 0ul; i < rBlk; i++) rsoff += nActO[i];
          #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, l1, l, j23off, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++) {
            auto iThread = GetThreadID();
            j23off = j3 * LCatAlt[2] + j2 * LCatAlt[1];
            j = j1 * LCatAlt[0] + j23off;
            const int * exList_j1 = LrsList->pointerAtDet(j1);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j1+=4) {
              UNPACK_EXCITATIONLIST_4(exList_j1, r, s, l1, sign);
              l = l1 * LCatAlt[0] + j23off; // <L|rs|J>
              SCR[iThread](r+rsoff, s+rsoff) += SmartConj(Cm[l+fCat[iCatL]]) * sign * Cn[j+fCat[iCatJ]];
            }
          } // End openMP
        } else { // off-diagonal
          reassembleCatDet2RAS(LCat[iCatJ], nDetJ, rBlk, sBlk, LCatAlt, nDetAlt, LCatAlt2, LCat[iCatL]);

          int roff = 0;
          int soff = 0;
          for (auto i = 0ul; i < rBlk; i++) roff += nActO[i];
          for (auto i = 0ul; i < sBlk; i++) soff += nActO[i];

          #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, l1, l, l2, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++) {
            auto iThread = GetThreadID();
            j = j1*LCatAlt[0] + j2*LCatAlt[1] + j3*LCatAlt[2];
            const int * exList_j = LrsList->pointerAtDet(j1, j2);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j+=5) {
              UNPACK_EXCITATIONLIST_5(exList_j, r, s, l1, l2, sign);
              l = l1*LCatAlt2[0] + l2*LCatAlt2[1] + j3*LCatAlt2[2];
              SCR[iThread](r+roff, s+soff) += SmartConj(Cm[l+fCat[iCatL]]) * sign * Cn[j+fCat[iCatJ]];
            }
          } // End openMP
        }
      }
    }

    for (auto i = 0ul; i < nThreads; i++) TDM += SCR[i];
    size_t nO = mcwfn.MOPartition.nCorrO;
//    prettyPrintSmart(std::cout,"LL RAS TDM -- ", TDM.pointer(), nO, nO, nO);

    return;

  } // RASCI::computeTDM

  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::computeTwoRDM(MCWaveFunction<MatsT, IntsT> & mcwfn,
    MatsT * C, InCore4indexTPI<MatsT> & twoRDM) {

    CErr("not implemented");

  } // RASCI::computeTwoRDM



}; // namespace ChronusQ
