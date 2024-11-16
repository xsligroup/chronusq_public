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
#include <particleintegrals/inhouseaointegral.hpp>
#include <util/matout.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>

#include <util/threads.hpp>
#include <util/mpi.hpp>
#include <util/timer.hpp>
#include <util/math.hpp>

#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <particleintegrals/contract/direct.hpp>

#include <chrono>

#define _CONTRACTION_

#define _SEPARATED_SHZ_SCREEN_4C

#include <libcint.hpp>

namespace ChronusQ {


  /**************************/
  /* Libcint 4C-direct Implementation With Coulomb only-type of contraction*/
  /**************************/

  // For DCB Hamiltonian,
  // 12 density matrices upon input stored as
  // LL(MS,MX,MY,MZ), SS(MS,MX,MY,MZ), LS(MS,MX,MY,MZ)
  //
  // 12 contrated matrices upon output stored as
  // LL(MS,MX,MY,MZ), SS(MS,MX,MY,MZ), LS(MS,MX,MY,MZ)
  //
  //
  // Density Matrices are NOT assumed as Hermitian here
  //

  /*******************************************************************************/
  /*                                                                             */
  /* Libcint Batch 4C-direct Implementation With Coulomb only-type of contraction*/
  /*                                                                             */
  /*******************************************************************************/

  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintCoulombOnly(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyRelContraction<MatsT>> &matList,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    
    const size_t mMat = matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL ? 2: 1;
    const size_t nMat = matList.size();
    const size_t nBatch = nMat / mMat;

    DirectTPI<IntsT> &originalERI = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& originalBasisSet_ = originalERI.basisSet();
    Molecule& molecule_ = originalERI.molecule();

    if (originalBasisSet_.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet_.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    DirectTPI<IntsT> &eri = originalERI;
    
    // Determine the number of OpenMP threads
    size_t nThreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    // no need to do that as LA functions are not used in parallel region
    // SetLAThreads(1); // Turn off parallelism in LA functions
  
    bool HerDen = matList[0].HER;
  
    /****************************/
    /* Format Basis for Libcint */
    /****************************/

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;
    int iAtom, iShell, off;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);
    size_t cache_size = libcintCacheSize(matList[0].contType, atm, nAtoms, bas, nShells, env);

    enum ERI_2ND_DERIV {
      AxBx,
      AxBy,
      AxBz,
      AyBx,
      AyBy,
      AyBz,
      AzBx,
      AzBy,
      AzBz
    };

    /*-------------------------------------*/
    /* End of Basis Formatting for Libcint */
    /*-------------------------------------*/


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;

    // 81 is for fourth-derivative; 9 for second derivative
    int nERI;
    double *buffAll = nullptr;
    double *cacheAll = nullptr;
    double *ERIBuffer = nullptr;

    size_t maxShellSize = buffSize;
    size_t NB  = maxShellSize*4;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;
    size_t NB4_2 = 2*NB4;
    size_t NB4_3 = 3*NB4;

    const size_t nBasis   = basisSet_.nBasis;
    const size_t nShell   = basisSet_.nShell;


    
    // initialize AX and allocate scratch space

    for(auto iMat = 0; iMat < nMat; iMat++) {
      matList[iMat].AX->clear(); 
    }
   
    // make copies for each thread
    std::vector<std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>>> AXthreads;
    for(auto iTh = 0; iTh < nThreads; iTh++) {
      AXthreads.emplace_back();
      for(auto iMat = 0; iMat < nMat; iMat++) {
        AXthreads.back().push_back(
          std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(nBasis, 
            matList[iMat].AX->hasZ(), matList[iMat].AX->hasXY())); 
      }
    }

    SetLAThreads(1); // Turn off parallelism in LA functions
    
    #pragma omp parallel for
    for(auto iTh = 0; iTh < nThreads; iTh++) 
    for(auto iMat = 0; iMat < nMat; iMat++)
      AXthreads[iTh][iMat]->clear();
    
    SetLAThreads(LAThreads);// Turn threads for LA back on
    
    
    
    /************************************/
    /*                                  */
    /* Preparation of Schwarz ERI       */
    /*                                  */
    /************************************/


    double *SchwarzERI = nullptr, *SchwarzSSSS = nullptr, 
      *SchwarzGaunt = nullptr, *SchwarzGauge = nullptr;

    if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::BARE_COULOMB or
       matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL or
       matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
//      if(SchwarzERI == nullptr) eri.computeSchwarz();
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topERIchwarz = tick();
#endif
  
      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      SchwarzERI = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzERI,0,nShell*nShell*sizeof(double));
  
      #pragma omp parallel
      {
  
        size_t thread_id = GetThreadID();
        size_t n1,n2;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        // set up Schwarz screening
        for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s12 % nThreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s1);
          shls[3] = int(s2);
  
          auto nQuad = n1*n2*n1*n2;
  
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

          for(size_t iERI(0); iERI < nERI*nQuad; iERI++) 
            SchwarzERI[s1+s2*nShell] = std::max(SchwarzERI[s1+s2*nShell], std::abs(buff[iERI]));
  
          SchwarzERI[s1+s2*nShell] = std::sqrt(SchwarzERI[s1+s2*nShell]);
          SchwarzERI[s2+s1*nShell] = SchwarzERI[s1+s2*nShell];
  
        }
        }
  
      };

      CQMemManager::get().free(buffAll, cacheAll);
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durERIchwarz = tock(topERIchwarz);
      // std::cout << "ERI Schwarz took " <<  durERIchwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
  
    }; // SchwarzERI

    if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::SSSS or 
       matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL or
       matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topSSSSSchwarz = tick();
#endif
  
      nERI = 81;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      SchwarzSSSS = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzSSSS,0,nShell*nShell*sizeof(double));
  
      #pragma omp parallel
      {
  
        double C2 = 1./(4*SpeedOfLight*SpeedOfLight);
        size_t thread_id = GetThreadID();
        size_t n1,n2;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        // set up Schwarz screening
        for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s12 % nThreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s1);
          shls[3] = int(s2);
  
          auto nQuad = n1*n2*n1*n2;
  
          if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

          for(size_t iERI(0); iERI < nERI*nQuad; iERI++) 
            SchwarzSSSS[s1+s2*nShell] = std::max(SchwarzSSSS[s1+s2*nShell], std::abs(buff[iERI]));
  
          SchwarzSSSS[s1+s2*nShell] = C2*std::sqrt(SchwarzSSSS[s1+s2*nShell]);
          SchwarzSSSS[s2+s1*nShell] = SchwarzSSSS[s1+s2*nShell];
  
        }
        }
  
      };

      CQMemManager::get().free(buffAll, cacheAll);
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durSSSSSchwarz = tock(topSSSSSchwarz);
//      std::cout << "∇A∇B∇C∇D Schwarz took " <<  durSSSSSchwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
  
    }; //ShwarzSSSS
    
    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUNT ) {
    
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topGauntSchwarz = tick();
#endif
      
      nERI = 9;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      SchwarzGaunt = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzGaunt,0,nShell*nShell*sizeof(double));
    
      #pragma omp parallel
      {
  
        double C1 = 1./(2*SpeedOfLight);
        size_t thread_id = GetThreadID();
        size_t n1,n2;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        // set up Schwarz screening
        for(size_t s1(0), bf1_s(0), s12(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++, s12++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s12 % nThreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s1);
          shls[3] = int(s2);
  
          auto nQuad = n1*n2*n1*n2;
  
          if(int2e_ip1ip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

          for(size_t iERI(0); iERI < nERI*nQuad; iERI++) 
              SchwarzGaunt[s1+s2*nShell] = std::max(SchwarzGaunt[s1+s2*nShell], std::abs(buff[iERI]));
  
          SchwarzGaunt[s1+s2*nShell] = C1*std::sqrt(SchwarzGaunt[s1+s2*nShell]);
  
        }
        }
  
      };
      
      CQMemManager::get().free(buffAll, cacheAll);

#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durGauntSchwarz = tock(topGauntSchwarz);
//      std::cout << "Gaunt Schwarz took " <<  durGauntchwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
    
    }; // ShwarzGaunt 

    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUGE ) {
    
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topGaugeSchwarz = tick();
#endif
    
      nERI = 16;
      buffAll = CQMemManager::get().malloc<double>(2*nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      std::cout << " cache_size = " << cache_size << std::endl;
      
      SchwarzGauge = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzGauge,0,nShell*nShell*sizeof(double));
      
      #pragma omp parallel
      {
  
        double C1 = 1./(2*SpeedOfLight);
        size_t thread_id = GetThreadID();
        size_t n1,n2;
        int skiperi1,skiperi2;
  
        int shls[4];
        double *buffr1 = buffAll+nERI*buffN4*thread_id;
        double *buffr2 = buffAll+nERI*nThreads*buffN4+nERI*buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
  
        // set up Schwarz screening
        for(size_t s1(0), bf1_s(0), s12(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++, s12++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s12 % nThreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s1);
          shls[3] = int(s2);
  
          auto nQuad = n1*n2*n1*n2;

          //∇B∇D
          skiperi1 = int2e_gauge_r1_ssp1ssp2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi2 = int2e_gauge_r2_ssp1ssp2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);

          if(skiperi1==0 and skiperi2==0) continue;

          for(size_t iERI(0); iERI < nERI*nQuad; iERI++) 
              SchwarzGauge[s1+s2*nShell] = std::max(SchwarzGauge[s1+s2*nShell], std::abs(buffr1[iERI]-buffr2[iERI]));
  
          SchwarzGauge[s1+s2*nShell] = C1*std::sqrt(SchwarzGauge[s1+s2*nShell]);
  
        }
        }
  
      };
      
      CQMemManager::get().free(buffAll, cacheAll);

#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durGaugeSchwarz = tock(topGaugeSchwarz);
//      std::cout << "Gauge Schwarz took " <<  durGaugechwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
    
    }; // ShwarzGauge 
    
    /****************************************/
    /*                                      */
    /* End Preparation of Schwarz ERIs      */
    /*                                      */
    /****************************************/
    

#ifdef _SHZ_SCREEN_4C
    /******************************************************/
    /*                                                    */
    /*Compute shell block norms (∞-norm) of all matList.X */
    /*                                                    */
    /******************************************************/
    
  
      // check Densitry matrices
      double * ShBlkNorms_raw = CQMemManager::get().malloc<double>((nBatch+1)*nShell*nShell);
      double * ShBlkNorms = ShBlkNorms_raw;
      std::vector<double*> ShBlkNorms_batch;
      { 
        double * ShBlkNormsSCR = CQMemManager::get().malloc<double>(nShell*nShell*nThreads);
        std::vector<double*> ShBlkNormsSCR_batch;
        for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) 
          ShBlkNorms_batch.push_back(ShBlkNorms_raw+(iBatch+1)*nShell*nShell);
        
        for( auto iTh = 0ul; iTh < nThreads; iTh++)
          ShBlkNormsSCR_batch.push_back(ShBlkNormsSCR+iTh*nShell*nShell); 
        
        memset(ShBlkNorms_raw,0.,(nBatch+1)*nShell*nShell*sizeof(double));
      
        #define POPMAX_SHELLBLOCKNORMS(DEN) \
          ShellBlockNorm(basisSet_.shells, DEN, nBasis, ShBlkNormsSCR_i); \
          for (auto i = 0ul; i < nShell*nShell; i++) { \
            ShBlkNorms_i[i] = std::max(ShBlkNorms_i[i],std::abs(ShBlkNormsSCR_i[i])); \
          }
        
        SetLAThreads(1); // Turn off parallelism in LA functions
        
        #pragma omp parallel for 
        for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) { 
          auto matOff = iBatch*mMat;
          auto ShBlkNorms_i    = ShBlkNorms_batch[iBatch]; 
          auto ShBlkNormsSCR_i = ShBlkNormsSCR_batch[GetThreadID()]; 
          for (auto iMat = 0ul; iMat < mMat; iMat++) {
            auto matX = matList[iMat + matOff].X;  
            POPMAX_SHELLBLOCKNORMS(matX->S().pointer()); 
            if (matX->hasZ()) POPMAX_SHELLBLOCKNORMS(matX->Z().pointer());
            if (matX->hasXY()) { 
              POPMAX_SHELLBLOCKNORMS(matX->Y().pointer());
              POPMAX_SHELLBLOCKNORMS(matX->X().pointer());
            } 
          }
          
          // symmetrize nonHermitian ShBlkNorms  ?? 
          if (not HerDen) {
            for (auto k = 0ul; k < nShell; k++)
            for (auto l = 0ul; l < k     ; l++) {
              double mx = std::max(ShBlkNorms_i[k + l*nShell], ShBlkNorms_i[l + k*nShell]);
              ShBlkNorms_i[k + l*nShell] = mx;
              ShBlkNorms_i[l + k*nShell] = mx;
            }
          }
          // prettyPrintSmart(std::cout, "shBlkNorm[" + std::to_string(iBatch) + "]", ShBlkNorms_i, nShell, nShell, nShell);
        }
        
        SetLAThreads(LAThreads);// Turn threads for LA back on
        CQMemManager::get().free(ShBlkNormsSCR);
        
        // get the maximum from all batches 
        #pragma omp parallel for 
        for (auto i = 0ul; i < nShell*nShell; i++) 
        for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) 
          ShBlkNorms[i] = std::max(ShBlkNorms_batch[iBatch][i], ShBlkNorms[i]); 
        
      } 
      // prettyPrintSmart(std::cout, "maxShBlkNormsSymmDenLLMS", ShBlkNorms, nShell, nShell, nShell);
    /******************************************************/
    /*                                                    */
    /*End of Compute shell block norms (∞-norm) of all matList.X */
    /*                                                    */
    /******************************************************/
#endif
    
    /***********************************/
    /*                                 */
    /* Start of Bare-Coulomb           */
    /*                                 */
    /***********************************/
  
    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::BARE_COULOMB ) {
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirect = tick();
      std::vector<double> durDirectInt(nThreads, 0.);
      std::vector<double> durDirectCon(nThreads, 0.);
#endif
  
      // Keeping track of number of integrals and contraction skipped
      std::vector<size_t> nIntSkip(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C      
      std::vector<size_t> nConSkip(nThreads,0);
#endif    
      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      
      #pragma omp parallel
      {
        int thread_id = GetThreadID();
        
        auto &AX_loc = AXthreads[thread_id];
        
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t bf43,bf34,bf21,bf12;
        size_t s4_max;
        int shls[4];
        double *buff = &buffAll[buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
        
        MatsT *ADCLLMS, *symmDCLLMS; 
        std::vector<size_t> contract_batch(nBatch);
        size_t iCon, nCon, matOff;

#if defined(_SHZ_SCREEN_4C) && defined(_SEPARATED_SHZ_SCREEN_4C)
        std::vector<double> shMax123_batch(nBatch);
#else
        std::iota(contract_batch.begin(), contract_batch.end(), 0);
        nCon = nBatch;
#endif

        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; 
            bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
          // Deneneracy factor for s1,s2 pair
          double s12_deg = (s1 == s2) ? 1.0 : 2.0;

#ifdef _SHZ_SCREEN_4C
          double shMax12 = ShBlkNorms[s1 + s2*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

#ifdef _SHZ_SCREEN_4C
          double shMax123 = std::max(ShBlkNorms[s1 + s3*nShell], 
                                     ShBlkNorms[s2 + s3*nShell]);

          shMax123 = std::max(shMax123,shMax12);
#ifdef _SEPARATED_SHZ_SCREEN_4C
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++)
             shMax123_batch[iBatch] = 
               std::max(ShBlkNorms_batch[iBatch][s1 + s2*nShell],
               std::max(ShBlkNorms_batch[iBatch][s1 + s3*nShell],
                        ShBlkNorms_batch[iBatch][s2 + s3*nShell]));
#endif
#endif 

        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nThreads != thread_id ) continue;
          #endif
  
#ifdef _SHZ_SCREEN_4C
          double shMax = std::max(ShBlkNorms[s1 + s4*nShell],
                         std::max(ShBlkNorms[s2 + s4*nShell],
                                  ShBlkNorms[s3 + s4*nShell]));

          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzERI[s1+s2*nShell]*SchwarzERI[s3+s4*nShell]) <
              eri.threshSchwarz()) {  
            nIntSkip[thread_id]++;
#ifdef _SEPARATED_SHZ_SCREEN_4C
            nConSkip[thread_id] += nBatch;
#endif            
            continue; 
          }
          
        // std::cout << "Shell Quads = [" << s1 << ", " << s2 << ", " << s3 << ", " << s4 << "]" << std::endl;
#ifdef _SEPARATED_SHZ_SCREEN_4C
          nCon = 0;
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) {
            shMax = std::max(ShBlkNorms_batch[iBatch][s1 + s4*nShell],
                    std::max(ShBlkNorms_batch[iBatch][s2 + s4*nShell],
                             ShBlkNorms_batch[iBatch][s3 + s4*nShell]));
            shMax = std::max(shMax, shMax123_batch[iBatch]); 
             
            if((shMax*SchwarzERI[s1+s2*nShell]*SchwarzERI[s3+s4*nShell]) <
              eri.threshSchwarz()) {
              nConSkip[thread_id] ++;
            } else {
              contract_batch[nCon] = iBatch; 
              nCon++;  
            }
          }
#endif
#endif 

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectInt = tick();
#endif
          // Degeneracy factor for s3,s4 pair
          double s34_deg = (s3 == s4) ? 1.0 : 2.0;
          // Degeneracy factor for s1, s2, s3, s4 quartet
          double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          // Total degeneracy factor
          double s1234_deg = s12_deg * s34_deg * s12_34_deg;
   
          auto nQuad = n1*n2*n3*n4;

          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);

          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

          // Scale the buffer by the degeneracy factor and store
          for(auto i = 0ul; i < nQuad; i++) buff[i] *= 0.5*s1234_deg;

#ifdef _REPORT_INTEGRAL_TIMINGS
          durDirectInt[thread_id] += tock(topDirectInt); 
#endif

#ifdef _CONTRACTION_ // Contraction
          
#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectCon = tick();
#endif
          for (auto iCon = 0ul; iCon < nCon; iCon++) {  
            
            matOff = contract_batch[iCon] * mMat; 
            ADCLLMS = AX_loc[matOff]->S().pointer();
            symmDCLLMS  = matList[matOff].X->S().pointer();

          for(l = 0ul, bf4 = bf4_s, mnkl=0ul ; l < n4; ++l, bf4++) {
  
          for(k = 0ul, bf3 = bf3_s           ; k < n3; ++k, bf3++) {
            bf43 = bf4 + bf3*nBasis;
            bf34 = bf3 + bf4*nBasis;
  
          for(n = 0ul, bf2 = bf2_s           ; n < n2; ++n, bf2++) {
  
          for(m = 0ul, bf1 = bf1_s           ; m < n1; ++m, bf1++, ++mnkl) {
            bf21 = bf2 + bf1*nBasis;
            bf12 = bf1 + bf2*nBasis;
 
            ADCLLMS[bf12] += buff[mnkl]*symmDCLLMS[bf43];
            ADCLLMS[bf34] += buff[mnkl]*symmDCLLMS[bf21];
            
          }; // bf1
          }; // bf2
          }; // bf3
          }; // bf4  
          }; // iCon

#endif // Contraction
        
#ifdef _REPORT_INTEGRAL_TIMINGS
        durDirectCon[thread_id] += tock(topDirectCon);
#endif
        }; // s4
        }; // s3
        }; // s2
        }; // s1
          
      
      } // OpenMP context 
      
      // Take care of the symmetry for CLLMS
      
      MatsT scale = MatsT(0.25);

      for( auto iMat = 0; iMat < nMat;  iMat++ ) 
      for( auto iTh  = 0; iTh < nThreads; iTh++) {
        matList[iMat].AX->S() += scale * (AXthreads[iTh][iMat]->S() + AXthreads[iTh][iMat]->S().T());   
      }
     
      CQMemManager::get().free(buffAll, cacheAll);
    
#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipAcc = std::accumulate(nIntSkip.begin(),nIntSkip.end(),0);
      double durDirectIntAcc = std::accumulate(durDirectInt.begin(), durDirectInt.end(),0.);
      double durDirectConAcc = std::accumulate(durDirectCon.begin(), durDirectCon.end(),0.);
      
      std::cout << std::endl;
      std::cout << "Bare-Coulomb-Exchange Skipped Integral   : " << nIntSkipAcc << std::endl;
#ifdef _SEPARATED_SHZ_SCREEN_4C
      size_t nConSkipAcc = std::accumulate(nConSkip.begin(),nConSkip.end(),0);
      std::cout << "                      Skipped Contraction: " << nConSkipAcc << std::endl;
#endif

      auto durDirect = tock(topDirect);
      std::cout << "Bare-Coulomb-Exchange AO Direct Contraction took " <<  durDirect 
                << " s for " << nBatch << " Batch" << std::endl; 
      std::cout << "                            Build Integrals took " <<  durDirectIntAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << "                            Contractions    took " <<  durDirectConAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << std::endl;
#endif
    
    };

    /*********************************/
    /*                               */
    /* End of Bare-Coulomb           */
    /*                               */
    /*********************************/
  
    /******************************************/
    /*                                        */
    /* Start of Dirac-Coulomb C(2)            */
    /* includes DC-LLLL and DC-LLSs/SSLL      */      
    /*                                        */
    /******************************************/
    
    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL ) {

 
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectLL = tick();
      std::vector<double> durDirectLLInt(nThreads, 0.);
      std::vector<double> durDirectLLCon(nThreads, 0.);
#endif
      enum X_COMP {
        sDSS,
        sDLLMS,
      };

      enum AX_COMP {
        CLLMS,
        CSS,
      };

      // 81 is for fourth-derivative; 9 for second derivative
      nERI = 9;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(2*4*NB4*nThreads);

      // Keeping track of number of integrals skipped
      std::vector<size_t> nIntSkipLL(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C
      std::vector<size_t> nConSkipLL(nThreads,0);
#endif

      #pragma omp parallel
      {

        double C2 = 1./(4*SpeedOfLight*SpeedOfLight);
  
        size_t thread_id = GetThreadID();
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *ERIBuffAB   = &ERIBuffer[thread_id*4*NB4];
        double *ERIBuffCD   = &ERIBuffer[nThreads*4*NB4 + thread_id*4*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t bf3nB,bf34,bf43,bf1nB,bf21,bf12;
        size_t mNB2,mnNBkNB2,mnNB,mNB2nNB3,kmNB2nNB,kmNB2nNB3; 
        size_t KLMN, DotPrdKLMN, CrossZKLMN, CrossXKLMN, CrossYKLMN,
          MNKL, DotPrdMNKL, CrossXMNKL, CrossYMNKL, CrossZMNKL;

        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
        
        MatsT *ADCLLMS, *ADCSSMS, *ADCSSMX, *ADCSSMY, *ADCSSMZ,
          *symmDCLLMS, *symmDCSSMS, *symmDCSSMX, *symmDCSSMY, *symmDCSSMZ;

        std::vector<size_t> contract_batch(nBatch);
        size_t iCon, nCon, matOff;

#if !defined(_SHZ_SCREEN_4C) || !defined(_SEPARATED_SHZ_SCREEN_4C)
        std::iota(contract_batch.begin(), contract_batch.end(), 0);
        nCon = nBatch;
#endif
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
          // Deneneracy factor for s1,s2 pair
          double s12_deg = (s1 == s2) ? 1.0 : 2.0;
  
        for(size_t s3(0), bf3_s(0); s3 < nShell; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4
  
        for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4

          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nThreads != thread_id ) continue;
          #endif
          
          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
            if (not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) or bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) {
              nIntSkipLL[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipLL[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) { 
            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) { 
              nIntSkipLL[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipLL[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) { 
            if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)  
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) ) {
              nIntSkipLL[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipLL[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

#ifdef _SHZ_SCREEN_4C

          double shMax = std::max(ShBlkNorms[s1 + s2*nShell],ShBlkNorms[s3 + s4*nShell]);

          if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzERI[s3 + s4*nShell]) <
             eri.threshSchwarz()) { 
            nIntSkipLL[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
            nConSkipLL[thread_id] += nBatch; 
#endif            
            continue; 
          }

#ifdef _SEPARATED_SHZ_SCREEN_4C
          nCon = 0;
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) {
             
            shMax = std::max(ShBlkNorms_batch[iBatch][s1 + s2*nShell],
                             ShBlkNorms_batch[iBatch][s3 + s4*nShell]);
            
            if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzERI[s3+s4*nShell]) <
              eri.threshSchwarz()) {
              nConSkipLL[thread_id]++;
            } else {
              contract_batch[nCon] = iBatch;
              nCon++;
            }
          }
#endif

#endif

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectLLInt = tick();
#endif
          
          // Degeneracy factor for s3,s4 pair
          double s34_deg = (s3 == s4) ? 1.0 : 2.0;

          auto nQuad = n1*n2*n3*n4;
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);
  
          if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

 
          for(l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l) {
            auto lNB3 = l*NB3;
            auto lNB  = l*NB;

          for(k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k) {
            auto kNB2lNB3 = k*NB2 + lNB3;
            auto klNB  = k + lNB;

          for(n =   maxShellSize            ; n <   maxShellSize + n2; ++n) {
            auto nNBkNB2lNB3 = n*NB + kNB2lNB3;
            auto klNBnNB3    = klNB + n*NB3;

          for(m = 0                         ; m <                  n1; ++m, ++mnkl) {
   
            /* Dirac-Coulomb */
            // ∇A∙∇B(mn|kl)
            auto dAdotdB = buff[AxBx*nQuad+mnkl] + buff[AyBy*nQuad+mnkl] + buff[AzBz*nQuad+mnkl];
            // ∇Ax∇B(mn|kl)
            auto dAcrossdB_x =  buff[AyBz*nQuad+mnkl] - buff[AzBy*nQuad+mnkl];
            auto dAcrossdB_y = -buff[AxBz*nQuad+mnkl] + buff[AzBx*nQuad+mnkl];
            auto dAcrossdB_z =  buff[AxBy*nQuad+mnkl] - buff[AyBx*nQuad+mnkl];
  
            //auto MNKL = m + n*NB + k*NB2 + l*NB3;
            MNKL = m + nNBkNB2lNB3;
            //auto KLMN = k + l*NB + m*NB2 + n*NB3;
            KLMN = m*NB2 + klNBnNB3;
  
  
            // ∇A∙∇B(mn|kl) followed by ∇Ax∇B(mn|kl) X, Y, and Z
            // (mn|kl)
            ERIBuffAB[       MNKL] =  (double)s34_deg*dAdotdB;
            ERIBuffAB[   NB4+MNKL] =  (double)s34_deg*dAcrossdB_x;
            ERIBuffAB[ NB4_2+MNKL] =  (double)s34_deg*dAcrossdB_y;
            ERIBuffAB[ NB4_3+MNKL] =  (double)s34_deg*dAcrossdB_z;
  
            // ∇C∙∇D(kl|nm) followed by ∇Cx∇D(kl|nm) X, Y, and Z
            // (kl|mn)
            ERIBuffCD[       KLMN] =  (double)s12_deg*dAdotdB;
            ERIBuffCD[   NB4+KLMN] =  (double)s12_deg*dAcrossdB_x;
            ERIBuffCD[ NB4_2+KLMN] =  (double)s12_deg*dAcrossdB_y;
            ERIBuffCD[ NB4_3+KLMN] =  (double)s12_deg*dAcrossdB_z;
 
          }
          }
          } 
          } // ∇A∇B integral preparation loop

#ifdef _REPORT_INTEGRAL_TIMINGS
         durDirectLLInt[thread_id] += tock(topDirectLLInt);
#endif

#ifdef _CONTRACTION_ // Contraction

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectLLCon = tick();
#endif

          for (iCon = 0ul; iCon < nCon; iCon++) {  
            
            matOff = contract_batch[iCon] * mMat;
            
            ADCLLMS = AX_loc[CLLMS + matOff]->S().pointer();
            ADCSSMS = AX_loc[CSS + matOff]->S().pointer();
            ADCSSMX = AX_loc[CSS + matOff]->X().pointer();
            ADCSSMY = AX_loc[CSS + matOff]->Y().pointer();
            ADCSSMZ = AX_loc[CSS + matOff]->Z().pointer();

            symmDCLLMS  = matList[sDLLMS + matOff].X->S().pointer();
            symmDCSSMS  = matList[sDSS + matOff].X->S().pointer();
            symmDCSSMX  = matList[sDSS + matOff].X->X().pointer();
            symmDCSSMY  = matList[sDSS + matOff].X->Y().pointer();
            symmDCSSMZ  = matList[sDSS + matOff].X->Z().pointer();

          for(m = 0ul,            bf1 = bf1_s; m <                  n1; ++m, bf1++) {
            mNB2 = m*NB2;
            bf1nB = bf1*nBasis;

          for(n =   maxShellSize, bf2 = bf2_s; n <   maxShellSize + n2; ++n, bf2++) {
            mNB2nNB3 = mNB2 + n*NB3;
            mnNB     = m + n*NB;
            bf21 = bf2 + bf1nB;
            bf12 = bf1 + bf2*nBasis;

          for(k = 2*maxShellSize, bf3 = bf3_s; k < 2*maxShellSize + n3; ++k, bf3++) {
            mnNBkNB2  = mnNB + k*NB2;
            kmNB2nNB3 = k + mNB2nNB3;
            bf3nB = bf3*nBasis;

          for(l = 3*maxShellSize, bf4 = bf4_s; l < 3*maxShellSize + n4; ++l, bf4++) {
  
            MNKL = mnNBkNB2 + l*NB3;
            KLMN = kmNB2nNB3 + l*NB;
  
            bf43 = bf4 + bf3nB;
            bf34 = bf3 + bf4*nBasis;
  
            DotPrdMNKL = MNKL;
            CrossXMNKL = MNKL+NB4;
            CrossYMNKL = MNKL+NB4_2;
            CrossZMNKL = MNKL+NB4_3;
  
            DotPrdKLMN = KLMN;
            CrossXKLMN = KLMN+NB4;
            CrossYKLMN = KLMN+NB4_2;
            CrossZKLMN = KLMN+NB4_3;
  
            /*++++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb (LL|LL) Contraction */
            /*++++++++++++++++++++++++++++++++++++++++++++*/
  
            //KLMN
            if(bf3 >= bf4 ) {
              ADCLLMS[bf34] +=  ERIBuffCD[DotPrdKLMN] * symmDCSSMS[bf21]
                               + (  ERIBuffCD[CrossZKLMN] * symmDCSSMZ[bf21]
                                  + ERIBuffCD[CrossXKLMN] * symmDCSSMX[bf21]
                                  + ERIBuffCD[CrossYKLMN] * symmDCSSMY[bf21]) * dcomplex(0.,1.);
            }
           
            /*------------------------------------------*/
            /* End of Dirac-Coulomb (LL|LL) Contraction */
            /*------------------------------------------*/
  
  
            /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb C(2)-(SS|SS) Contraction */
            /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
  
            /* MNKL */
            if (bf1 >= bf2) {
              ADCSSMS[bf12] += ERIBuffAB[DotPrdMNKL] * symmDCLLMS[bf43];
              ADCSSMX[bf12] += ERIBuffAB[CrossXMNKL] * symmDCLLMS[bf43];
              ADCSSMY[bf12] += ERIBuffAB[CrossYMNKL] * symmDCLLMS[bf43];
              ADCSSMZ[bf12] += ERIBuffAB[CrossZMNKL] * symmDCLLMS[bf43];
            }
  
            /*-----------------------------------------------*/
            /* End of Dirac-Coulomb C(2)-(SS|SS) Contraction */
            /*-----------------------------------------------*/
          };
          };
          };  
          }; 
          }; // iBatch

#endif // Contraction

#ifdef _REPORT_INTEGRAL_TIMINGS
        durDirectLLCon[thread_id] += tock(topDirectLLCon); 
#endif
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1

      } // OpenMP context

   
      MatsT  scale = MatsT(0.5);
      MatsT iscale = scale * dcomplex(0.0, 1.0);
      
      for( auto iTh  = 0; iTh < nThreads; iTh++) {
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        matList[CLLMS + matOff].AX->S() +=  scale * AXthreads[iTh][CLLMS + matOff]->S();  
        matList[CSS   + matOff].AX->S() +=  scale * AXthreads[iTh][CSS   + matOff]->S();  
        matList[CSS   + matOff].AX->X() += iscale * AXthreads[iTh][CSS   + matOff]->X();  
        matList[CSS   + matOff].AX->Y() += iscale * AXthreads[iTh][CSS   + matOff]->Y();  
        matList[CSS   + matOff].AX->Z() += iscale * AXthreads[iTh][CSS   + matOff]->Z();  
      }}
    
      // Take care of the symmetry in the LL and SS blocks
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        auto ADCLLMS = matList[CLLMS + matOff].AX->S().pointer();
        auto ADCSSMS = matList[CSS   + matOff].AX->S().pointer();
        auto ADCSSMX = matList[CSS   + matOff].AX->X().pointer();
        auto ADCSSMY = matList[CSS   + matOff].AX->Y().pointer();
        auto ADCSSMZ = matList[CSS   + matOff].AX->Z().pointer();
        for( auto i = 0; i < nBasis; i++ ) 
        for( auto j = 0; j < i; j++ ) {
          ADCLLMS[j + i*nBasis] =  ADCLLMS[i + j*nBasis];
          ADCSSMS[j + i*nBasis] =  ADCSSMS[i + j*nBasis];
          ADCSSMX[j + i*nBasis] = -ADCSSMX[i + j*nBasis];
          ADCSSMY[j + i*nBasis] = -ADCSSMY[i + j*nBasis];
          ADCSSMZ[j + i*nBasis] = -ADCSSMZ[i + j*nBasis];
        }
      } 
      
      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
    
#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipLLAcc = std::accumulate(nIntSkipLL.begin(),nIntSkipLL.end(),0);
      double durDirectLLIntAcc = std::accumulate(durDirectLLInt.begin(), durDirectLLInt.end(),0.);
      double durDirectLLConAcc = std::accumulate(durDirectLLCon.begin(), durDirectLLCon.end(),0.);
      
      std::cout << std::endl;
      std::cout << "Dirac-Coulomb-LL/C(2)-SS Skipped Integral   : " << nIntSkipLLAcc << std::endl;
#ifdef _SEPARATED_SHZ_SCREEN_4C
      size_t nConSkipLLAcc = std::accumulate(nConSkipLL.begin(),nConSkipLL.end(),0);
      std::cout << "                         Skipped Contraction: " << nConSkipLLAcc << std::endl;
#endif
      auto durDirectLL = tock(topDirectLL);
      std::cout << "Dirac-Coulomb-LL/C(2)-SS AO Direct Contraction took " <<  durDirectLL 
                << " s for " << nBatch << " Batch " << std::endl; 
      std::cout << "                               Build Integrals took " <<  durDirectLLIntAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << "                               Contractions    took " <<  durDirectLLConAcc << " s"
                << " from all " << nThreads << " thread(s)" << std::endl;

      std::cout << std::endl;
#endif
    
    } // if(LLLL)
    
    /******************************************/
    /*                                        */
    /*   End of Dirac-Coulomb LL and C(2)-SS  */
    /*                                        */
    /******************************************/
    
    
    
    /****************************************/
    /*                                      */
    /* Start of Dirac-Coulomb C(4)-SSSS     */
    /*                                      */
    /****************************************/

    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::SSSS ) {
  
      enum ERI_4TH_DERIV {
        AxBxCxDx,
        AxBxCxDy,
        AxBxCxDz,
        AxBxCyDx,
        AxBxCyDy,
        AxBxCyDz,
        AxBxCzDx,
        AxBxCzDy,
        AxBxCzDz,
        AxByCxDx,
        AxByCxDy,
        AxByCxDz,
        AxByCyDx,
        AxByCyDy,
        AxByCyDz,
        AxByCzDx,
        AxByCzDy,
        AxByCzDz,
        AxBzCxDx,
        AxBzCxDy,
        AxBzCxDz,
        AxBzCyDx,
        AxBzCyDy,
        AxBzCyDz,
        AxBzCzDx,
        AxBzCzDy,
        AxBzCzDz,
        AyBxCxDx,
        AyBxCxDy,
        AyBxCxDz,
        AyBxCyDx,
        AyBxCyDy,
        AyBxCyDz,
        AyBxCzDx,
        AyBxCzDy,
        AyBxCzDz,
        AyByCxDx,
        AyByCxDy,
        AyByCxDz,
        AyByCyDx,
        AyByCyDy,
        AyByCyDz,
        AyByCzDx,
        AyByCzDy,
        AyByCzDz,
        AyBzCxDx,
        AyBzCxDy,
        AyBzCxDz,
        AyBzCyDx,
        AyBzCyDy,
        AyBzCyDz,
        AyBzCzDx,
        AyBzCzDy,
        AyBzCzDz,
        AzBxCxDx,
        AzBxCxDy,
        AzBxCxDz,
        AzBxCyDx,
        AzBxCyDy,
        AzBxCyDz,
        AzBxCzDx,
        AzBxCzDy,
        AzBxCzDz,
        AzByCxDx,
        AzByCxDy,
        AzByCxDz,
        AzByCyDx,
        AzByCyDy,
        AzByCyDz,
        AzByCzDx,
        AzByCzDy,
        AzByCzDz,
        AzBzCxDx,
        AzBzCxDy,
        AzBzCxDz,
        AzBzCyDx,
        AzBzCyDy,
        AzBzCyDz,
        AzBzCzDx,
        AzBzCzDy,
        AzBzCzDz,
      };
  
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectSSSS = tick();
      std::vector<double> durDirectSSSSInt(nThreads, 0.);
      std::vector<double> durDirectSSSSCon(nThreads, 0.);
#endif

      nERI = 81;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(16*NB4*nThreads);

      // Keeping track of number of integrals skipped
      std::vector<size_t> nIntSkipSSSS(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C
      std::vector<size_t> nConSkipSSSS(nThreads,0);
#endif
      
      auto NB4_4  = 4*NB4;  
      auto NB4_5  = 5*NB4;  
      auto NB4_6  = 6*NB4;  
      auto NB4_7  = 7*NB4;  
      auto NB4_8  = 8*NB4;  
      auto NB4_9  = 9*NB4;  
      auto NB4_10 =10*NB4;  
      auto NB4_11 =11*NB4;  
      auto NB4_12 =12*NB4;  
      auto NB4_13 =13*NB4;  
      auto NB4_14 =14*NB4;  
      auto NB4_15 =15*NB4;  

      #pragma omp parallel
      {
  
        dcomplex iscale = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *ERIBuffABCD = &ERIBuffer[thread_id*16*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        MatsT *ADCSSMS, *ADCSSMX, *ADCSSMY, *ADCSSMZ,
          *symmDCSSMS, *symmDCSSMX, *symmDCSSMY, *symmDCSSMZ; 
        
        std::vector<size_t> contract_batch(nBatch);
        size_t iCon, nCon, matOff;

#if defined(_SHZ_SCREEN_4C) && defined(_SEPARATED_SHZ_SCREEN_4C)
        std::vector<double> shMax123_batch(nBatch);
#else
        std::iota(contract_batch.begin(), contract_batch.end(), 0);
        nCon = nBatch;
#endif
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet_.nShell; 
            bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C
          double shMax12 = ShBlkNorms[s1 + s2*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C
          double shMax123 = std::max(ShBlkNorms[s1 + s3*nShell], 
                                     ShBlkNorms[s2 + s3*nShell]);

          shMax123 = std::max(shMax123,shMax12);
#ifdef _SEPARATED_SHZ_SCREEN_4C
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++)
             shMax123_batch[iBatch] = 
               std::max(ShBlkNorms_batch[iBatch][s1 + s2*nShell],
               std::max(ShBlkNorms_batch[iBatch][s1 + s3*nShell],
                        ShBlkNorms_batch[iBatch][s2 + s3*nShell]));
#endif
#endif

        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nThreads != thread_id ) continue;
          #endif

          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
            if (not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) or bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) {
              nIntSkipSSSS[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipSSSS[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) { 
            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) { 
              nIntSkipSSSS[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipSSSS[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) { 
            if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) ) {
              nIntSkipSSSS[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipSSSS[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

#ifdef _SHZ_SCREEN_4C
          double shMax = std::max(ShBlkNorms[s1 + s4*nShell],
                         std::max(ShBlkNorms[s2 + s4*nShell],
                                  ShBlkNorms[s3 + s4*nShell]));

          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzSSSS[s3+s4*nShell]) <
              eri.threshSchwarz()) {  
            nIntSkipSSSS[thread_id]++;
#ifdef _SEPARATED_SHZ_SCREEN_4C
            nConSkipSSSS[thread_id] += nBatch;
#endif            
            continue; 
          }

#ifdef _SEPARATED_SHZ_SCREEN_4C
          nCon = 0;
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) {
            shMax = std::max(ShBlkNorms_batch[iBatch][s1 + s4*nShell],
                    std::max(ShBlkNorms_batch[iBatch][s2 + s4*nShell],
                             ShBlkNorms_batch[iBatch][s3 + s4*nShell]));
            shMax = std::max(shMax, shMax123_batch[iBatch]); 
             
            if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzSSSS[s3+s4*nShell]) <
              eri.threshSchwarz()) {
              nConSkipSSSS[thread_id] ++;
            } else {
              contract_batch[nCon] = iBatch; 
              nCon++;  
            }
          }
#endif
#endif

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectSSSSInt = tick();
#endif

          auto nQuad = n1*n2*n3*n4;
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);
  
          if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
  
          for(l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
          for(k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k) 
          for(n =   maxShellSize            ; n <   maxShellSize + n2; ++n) 
          for(m = 0                         ; m <                  n1; ++m, ++mnkl) {
  
            // (∇A∙∇B)(∇C∙∇D)(mnkl)
            auto dAdotdBdCdotdD =  buff[AxBxCxDx*nQuad+mnkl] + buff[AxBxCyDy*nQuad+mnkl] + buff[AxBxCzDz*nQuad+mnkl]
                                 + buff[AyByCxDx*nQuad+mnkl] + buff[AyByCyDy*nQuad+mnkl] + buff[AyByCzDz*nQuad+mnkl]
                                 + buff[AzBzCxDx*nQuad+mnkl] + buff[AzBzCyDy*nQuad+mnkl] + buff[AzBzCzDz*nQuad+mnkl];
  
            // (∇Ax∇B)(∇C∙∇D)(mnkl)
            auto dAcrossdB_xdCdotdD =  buff[AyBzCxDx*nQuad+mnkl] - buff[AzByCxDx*nQuad+mnkl]
                                     + buff[AyBzCyDy*nQuad+mnkl] - buff[AzByCyDy*nQuad+mnkl]
                                     + buff[AyBzCzDz*nQuad+mnkl] - buff[AzByCzDz*nQuad+mnkl];
  
            auto dAcrossdB_ydCdotdD = -buff[AxBzCxDx*nQuad+mnkl] + buff[AzBxCxDx*nQuad+mnkl]
                                      -buff[AxBzCyDy*nQuad+mnkl] + buff[AzBxCyDy*nQuad+mnkl]
                                      -buff[AxBzCzDz*nQuad+mnkl] + buff[AzBxCzDz*nQuad+mnkl];
                                      
            auto dAcrossdB_zdCdotdD =  buff[AxByCxDx*nQuad+mnkl] - buff[AyBxCxDx*nQuad+mnkl]
                                     + buff[AxByCyDy*nQuad+mnkl] - buff[AyBxCyDy*nQuad+mnkl]
                                     + buff[AxByCzDz*nQuad+mnkl] - buff[AyBxCzDz*nQuad+mnkl];
    
            // (∇A∙∇B)(∇Cx∇D)(mnkl)
            auto dAdotdBdCcrossdD_x =  buff[AxBxCyDz*nQuad+mnkl] - buff[AxBxCzDy*nQuad+mnkl]
                                     + buff[AyByCyDz*nQuad+mnkl] - buff[AyByCzDy*nQuad+mnkl]
                                     + buff[AzBzCyDz*nQuad+mnkl] - buff[AzBzCzDy*nQuad+mnkl];
  
            auto dAdotdBdCcrossdD_y = -buff[AxBxCxDz*nQuad+mnkl] + buff[AxBxCzDx*nQuad+mnkl]
                                      -buff[AyByCxDz*nQuad+mnkl] + buff[AyByCzDx*nQuad+mnkl]
                                      -buff[AzBzCxDz*nQuad+mnkl] + buff[AzBzCzDx*nQuad+mnkl];
                                      
            auto dAdotdBdCcrossdD_z =  buff[AxBxCxDy*nQuad+mnkl] - buff[AxBxCyDx*nQuad+mnkl]
                                     + buff[AyByCxDy*nQuad+mnkl] - buff[AyByCyDx*nQuad+mnkl]
                                     + buff[AzBzCxDy*nQuad+mnkl] - buff[AzBzCyDx*nQuad+mnkl];
  
            // (∇Ax∇B)(∇Cx∇D)(mnkl)
            auto dAcrossdB_xdCcrossdD_x =  buff[AyBzCyDz*nQuad+mnkl] - buff[AzByCyDz*nQuad+mnkl]
                                         - buff[AyBzCzDy*nQuad+mnkl] + buff[AzByCzDy*nQuad+mnkl];
  
            auto dAcrossdB_xdCcrossdD_y =  buff[AyBzCzDx*nQuad+mnkl] - buff[AzByCzDx*nQuad+mnkl]
                                         - buff[AyBzCxDz*nQuad+mnkl] + buff[AzByCxDz*nQuad+mnkl];
  
            auto dAcrossdB_xdCcrossdD_z =  buff[AyBzCxDy*nQuad+mnkl] - buff[AzByCxDy*nQuad+mnkl]
                                         - buff[AyBzCyDx*nQuad+mnkl] + buff[AzByCyDx*nQuad+mnkl];
  
            auto dAcrossdB_ydCcrossdD_x =  buff[AzBxCyDz*nQuad+mnkl] - buff[AxBzCyDz*nQuad+mnkl]
                                         - buff[AzBxCzDy*nQuad+mnkl] + buff[AxBzCzDy*nQuad+mnkl];
  
            auto dAcrossdB_ydCcrossdD_y =  buff[AzBxCzDx*nQuad+mnkl] - buff[AxBzCzDx*nQuad+mnkl]
                                         - buff[AzBxCxDz*nQuad+mnkl] + buff[AxBzCxDz*nQuad+mnkl];
  
            auto dAcrossdB_ydCcrossdD_z =  buff[AzBxCxDy*nQuad+mnkl] - buff[AxBzCxDy*nQuad+mnkl]
                                         - buff[AzBxCyDx*nQuad+mnkl] + buff[AxBzCyDx*nQuad+mnkl];
  
            auto dAcrossdB_zdCcrossdD_x =  buff[AxByCyDz*nQuad+mnkl] - buff[AyBxCyDz*nQuad+mnkl]
                                         - buff[AxByCzDy*nQuad+mnkl] + buff[AyBxCzDy*nQuad+mnkl];
  
            auto dAcrossdB_zdCcrossdD_y =  buff[AxByCzDx*nQuad+mnkl] - buff[AyBxCzDx*nQuad+mnkl]
                                         - buff[AxByCxDz*nQuad+mnkl] + buff[AyBxCxDz*nQuad+mnkl];
  
            auto dAcrossdB_zdCcrossdD_z =  buff[AxByCxDy*nQuad+mnkl] - buff[AyBxCxDy*nQuad+mnkl]
                                         - buff[AxByCyDx*nQuad+mnkl] + buff[AyBxCyDx*nQuad+mnkl];
  
  
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto KLMN = k + l*NB + m*NB2 + n*NB3;

            // (mn|kl)
            ERIBuffABCD[         MNKL] =  dAdotdBdCdotdD;
            ERIBuffABCD[   NB4 + MNKL] =  dAcrossdB_xdCdotdD;
            ERIBuffABCD[ NB4_2 + MNKL] =  dAcrossdB_ydCdotdD;
            ERIBuffABCD[ NB4_3 + MNKL] =  dAcrossdB_zdCdotdD;
            ERIBuffABCD[ NB4_4 + MNKL] =  dAdotdBdCcrossdD_x;
            ERIBuffABCD[ NB4_5 + MNKL] =  dAdotdBdCcrossdD_y;
            ERIBuffABCD[ NB4_6 + MNKL] =  dAdotdBdCcrossdD_z;
            ERIBuffABCD[ NB4_7 + MNKL] =  dAcrossdB_xdCcrossdD_x;
            ERIBuffABCD[ NB4_8 + MNKL] =  dAcrossdB_xdCcrossdD_y;
            ERIBuffABCD[ NB4_9 + MNKL] =  dAcrossdB_xdCcrossdD_z;
            ERIBuffABCD[NB4_10 + MNKL] =  dAcrossdB_ydCcrossdD_x;
            ERIBuffABCD[NB4_11 + MNKL] =  dAcrossdB_ydCcrossdD_y;
            ERIBuffABCD[NB4_12 + MNKL] =  dAcrossdB_ydCcrossdD_z;
            ERIBuffABCD[NB4_13 + MNKL] =  dAcrossdB_zdCcrossdD_x; 
            ERIBuffABCD[NB4_14 + MNKL] =  dAcrossdB_zdCcrossdD_y;
            ERIBuffABCD[NB4_15 + MNKL] =  dAcrossdB_zdCcrossdD_z;
  
            // (kl|mn)
            ERIBuffABCD[         KLMN] =  dAdotdBdCdotdD;
            ERIBuffABCD[   NB4 + KLMN] =  dAdotdBdCcrossdD_x;
            ERIBuffABCD[ NB4_2 + KLMN] =  dAdotdBdCcrossdD_y;
            ERIBuffABCD[ NB4_3 + KLMN] =  dAdotdBdCcrossdD_z;
            ERIBuffABCD[ NB4_4 + KLMN] =  dAcrossdB_xdCdotdD;
            ERIBuffABCD[ NB4_5 + KLMN] =  dAcrossdB_ydCdotdD;
            ERIBuffABCD[ NB4_6 + KLMN] =  dAcrossdB_zdCdotdD;
            ERIBuffABCD[ NB4_7 + KLMN] =  dAcrossdB_xdCcrossdD_x;
            ERIBuffABCD[ NB4_8 + KLMN] =  dAcrossdB_ydCcrossdD_x;
            ERIBuffABCD[ NB4_9 + KLMN] =  dAcrossdB_zdCcrossdD_x;
            ERIBuffABCD[NB4_10 + KLMN] =  dAcrossdB_xdCcrossdD_y;
            ERIBuffABCD[NB4_11 + KLMN] =  dAcrossdB_ydCcrossdD_y;
            ERIBuffABCD[NB4_12 + KLMN] =  dAcrossdB_zdCcrossdD_y;
            ERIBuffABCD[NB4_13 + KLMN] =  dAcrossdB_xdCcrossdD_z; 
            ERIBuffABCD[NB4_14 + KLMN] =  dAcrossdB_ydCcrossdD_z;
            ERIBuffABCD[NB4_15 + KLMN] =  dAcrossdB_zdCcrossdD_z;
  
          } // ∇A∇B∇C∇D (SSSS) integrals 
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          durDirectSSSSInt[thread_id] += tock(topDirectSSSSInt); 
#endif
  
  
#ifdef _CONTRACTION_ // Contraction
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectSSSSCon = tick();
#endif
          for (auto iCon = 0ul; iCon < nCon; iCon++) {  
          
            matOff = contract_batch[iCon] * mMat; 
            symmDCSSMS = matList[matOff].X->S().pointer();
            symmDCSSMX = matList[matOff].X->X().pointer();
            symmDCSSMY = matList[matOff].X->Y().pointer();
            symmDCSSMZ = matList[matOff].X->Z().pointer();
  
            ADCSSMS  = AX_loc[matOff]->S().pointer();
            ADCSSMX  = AX_loc[matOff]->X().pointer();
            ADCSSMY  = AX_loc[matOff]->Y().pointer();
            ADCSSMZ  = AX_loc[matOff]->Z().pointer();
  
          for(m = 0ul,            bf1 = bf1_s; m <                  n1; ++m, bf1++) {

          for(n =   maxShellSize, bf2 = bf2_s; n <   maxShellSize + n2; ++n, bf2++) {
            auto bf1nB = bf1*nBasis;
            auto bf2nB = bf2*nBasis;
            auto bf21 = bf2 + bf1nB;
            auto bf12 = bf1 + bf2nB;
            auto mn1 = m + n*NB;
            auto mn2 = m*NB2 + n*NB3;

          for(k = 2*maxShellSize, bf3 = bf3_s; k < 2*maxShellSize + n3; ++k, bf3++) {
            auto bf3nB = bf3*nBasis;
            auto kmn1 = mn1 + k*NB2;
            auto kmn2 = k + mn2;

          for(l = 3*maxShellSize, bf4 = bf4_s; l < 3*maxShellSize + n4; ++l, bf4++) {
            auto MNKL = kmn1 + l*NB3;
            auto KLMN = kmn2 + l*NB;

            auto bf4nB = bf4*nBasis;
  
            auto bf14 = bf1 + bf4nB;
            auto bf24 = bf2 + bf4nB;
            auto bf23 = bf2 + bf3nB;
            auto bf13 = bf1 + bf3nB;

            auto bf41 = bf4 + bf1nB;
            auto bf31 = bf3 + bf1nB;
            auto bf32 = bf3 + bf2nB;
            auto bf42 = bf4 + bf2nB;
 

            auto bf43 = bf4 + bf3nB;
            auto bf34 = bf3 + bf4nB;
  
 

            /*********************************/
            /* Dirac-Coulomb (SSSS) Integral */
            /*********************************/
      
  
            auto MNKLdAdotdBdCdotdD          = ERIBuffABCD[         MNKL];
            auto MNKLdAcrossdB_xdCdotdD      = ERIBuffABCD[   NB4 + MNKL];
            auto MNKLdAcrossdB_ydCdotdD      = ERIBuffABCD[ NB4_2 + MNKL];
            auto MNKLdAcrossdB_zdCdotdD      = ERIBuffABCD[ NB4_3 + MNKL];
            auto MNKLdAdotdBdCcrossdD_x      = ERIBuffABCD[ NB4_4 + MNKL];
            auto MNKLdAdotdBdCcrossdD_y      = ERIBuffABCD[ NB4_5 + MNKL];
            auto MNKLdAdotdBdCcrossdD_z      = ERIBuffABCD[ NB4_6 + MNKL];
            auto MNKLdAcrossdB_xdCcrossdD_x  = ERIBuffABCD[ NB4_7 + MNKL];
            auto MNKLdAcrossdB_xdCcrossdD_y  = ERIBuffABCD[ NB4_8 + MNKL];
            auto MNKLdAcrossdB_xdCcrossdD_z  = ERIBuffABCD[ NB4_9 + MNKL];
            auto MNKLdAcrossdB_ydCcrossdD_x  = ERIBuffABCD[NB4_10 + MNKL];
            auto MNKLdAcrossdB_ydCcrossdD_y  = ERIBuffABCD[NB4_11 + MNKL];
            auto MNKLdAcrossdB_ydCcrossdD_z  = ERIBuffABCD[NB4_12 + MNKL];
            auto MNKLdAcrossdB_zdCcrossdD_x  = ERIBuffABCD[NB4_13 + MNKL];
            auto MNKLdAcrossdB_zdCcrossdD_y  = ERIBuffABCD[NB4_14 + MNKL];
            auto MNKLdAcrossdB_zdCcrossdD_z  = ERIBuffABCD[NB4_15 + MNKL];

            auto KLMNdAdotdBdCdotdD          = ERIBuffABCD[         KLMN];
            auto KLMNdAcrossdB_xdCdotdD      = ERIBuffABCD[   NB4 + KLMN];
            auto KLMNdAcrossdB_ydCdotdD      = ERIBuffABCD[ NB4_2 + KLMN];
            auto KLMNdAcrossdB_zdCdotdD      = ERIBuffABCD[ NB4_3 + KLMN];
            auto KLMNdAdotdBdCcrossdD_x      = ERIBuffABCD[ NB4_4 + KLMN];
            auto KLMNdAdotdBdCcrossdD_y      = ERIBuffABCD[ NB4_5 + KLMN];
            auto KLMNdAdotdBdCcrossdD_z      = ERIBuffABCD[ NB4_6 + KLMN];
            auto KLMNdAcrossdB_xdCcrossdD_x  = ERIBuffABCD[ NB4_7 + KLMN];
            auto KLMNdAcrossdB_xdCcrossdD_y  = ERIBuffABCD[ NB4_8 + KLMN];
            auto KLMNdAcrossdB_xdCcrossdD_z  = ERIBuffABCD[ NB4_9 + KLMN];
            auto KLMNdAcrossdB_ydCcrossdD_x  = ERIBuffABCD[NB4_10 + KLMN];
            auto KLMNdAcrossdB_ydCcrossdD_y  = ERIBuffABCD[NB4_11 + KLMN];
            auto KLMNdAcrossdB_ydCcrossdD_z  = ERIBuffABCD[NB4_12 + KLMN];
            auto KLMNdAcrossdB_zdCcrossdD_x  = ERIBuffABCD[NB4_13 + KLMN];
            auto KLMNdAcrossdB_zdCcrossdD_y  = ERIBuffABCD[NB4_14 + KLMN];
            auto KLMNdAcrossdB_zdCcrossdD_z  = ERIBuffABCD[NB4_15 + KLMN];
 
            // COULOMB
            // MNKL
            auto ScaleF = bf3_s==bf4_s ? 0.5 : 1.0;
  
            if(bf1 >= bf2) {

              /* Equation 70 in the paper */
              ADCSSMS[bf12] += ScaleF * (
                     symmDCSSMS[bf43] * MNKLdAdotdBdCdotdD
                 + ( symmDCSSMZ[bf43] * MNKLdAdotdBdCcrossdD_z
                    +symmDCSSMX[bf43] * MNKLdAdotdBdCcrossdD_x
                    +symmDCSSMY[bf43] * MNKLdAdotdBdCcrossdD_y) * dcomplex(0.,1.) );
               
              /* Equation 71 in the paper */
              ADCSSMZ[bf12] += ScaleF * (
                  symmDCSSMS[bf43] * MNKLdAcrossdB_zdCdotdD * dcomplex(0.,1.)
                 -symmDCSSMZ[bf43] * MNKLdAcrossdB_zdCcrossdD_z
                 -symmDCSSMX[bf43] * MNKLdAcrossdB_zdCcrossdD_x
                 -symmDCSSMY[bf43] * MNKLdAcrossdB_zdCcrossdD_y );
       
              /* Equation 72 in the paper */
              ADCSSMX[bf12] += ScaleF * (
                  symmDCSSMS[bf43] * MNKLdAcrossdB_xdCdotdD * dcomplex(0.,1.)
                 -symmDCSSMZ[bf43] * MNKLdAcrossdB_xdCcrossdD_z
                 -symmDCSSMX[bf43] * MNKLdAcrossdB_xdCcrossdD_x
                 -symmDCSSMY[bf43] * MNKLdAcrossdB_xdCcrossdD_y );
        
              /* Equation 73 in the paper */
              ADCSSMY[bf12] += ScaleF * (
                  symmDCSSMS[bf43] * MNKLdAcrossdB_ydCdotdD * dcomplex(0.,1.)
                 -symmDCSSMZ[bf43] * MNKLdAcrossdB_ydCcrossdD_z
                 -symmDCSSMX[bf43] * MNKLdAcrossdB_ydCcrossdD_x
                 -symmDCSSMY[bf43] * MNKLdAcrossdB_ydCcrossdD_y );
            }
  
            /* KLMN */
            if((bf1_s!=bf3_s or bf2_s!=bf4_s) and bf3 >= bf4){
 
              ScaleF = bf1_s==bf2_s ? 0.5 : 1.0;
              
              /* Equation 70 in the paper */
              ADCSSMS[bf34] += ScaleF * (
                     symmDCSSMS[bf21] * KLMNdAdotdBdCdotdD
                 + ( symmDCSSMZ[bf21] * KLMNdAdotdBdCcrossdD_z
                    +symmDCSSMX[bf21] * KLMNdAdotdBdCcrossdD_x
                    +symmDCSSMY[bf21] * KLMNdAdotdBdCcrossdD_y) * dcomplex(0.,1.) );
  
              /* Equation 71 in the paper */
              ADCSSMZ[bf34] += ScaleF * (
                  symmDCSSMS[bf21] * KLMNdAcrossdB_zdCdotdD * dcomplex(0.,1.) 
                 -symmDCSSMZ[bf21] * KLMNdAcrossdB_zdCcrossdD_z
                 -symmDCSSMX[bf21] * KLMNdAcrossdB_zdCcrossdD_x
                 -symmDCSSMY[bf21] * KLMNdAcrossdB_zdCcrossdD_y );
      
              /* Equation 72 in the paper */
              ADCSSMX[bf34] += ScaleF * (
                  symmDCSSMS[bf21] * KLMNdAcrossdB_xdCdotdD * dcomplex(0.,1.)
                 -symmDCSSMZ[bf21] * KLMNdAcrossdB_xdCcrossdD_z
                 -symmDCSSMX[bf21] * KLMNdAcrossdB_xdCcrossdD_x
                 -symmDCSSMY[bf21] * KLMNdAcrossdB_xdCcrossdD_y );
      
              /* Equation 73 in the paper */
              ADCSSMY[bf34] += ScaleF * (
                  symmDCSSMS[bf21] * KLMNdAcrossdB_ydCdotdD * dcomplex(0.,1.)
                 -symmDCSSMZ[bf21] * KLMNdAcrossdB_ydCcrossdD_z
                 -symmDCSSMX[bf21] * KLMNdAcrossdB_ydCcrossdD_x
                 -symmDCSSMY[bf21] * KLMNdAcrossdB_ydCcrossdD_y);
            }
  
	  }
	  }
	  }
      } 
      } // contraction loop
  
#endif // Contraction
  
#ifdef _REPORT_INTEGRAL_TIMINGS
        durDirectSSSSCon[thread_id] += tock(topDirectSSSSCon); 
#endif
  
      }; // loop s4
      }; // loop s3
      }; // loop s2
      }; // loop s1
  
  
      }; // OpenMP context
  
      for( auto iTh  = 0; iTh < nThreads; iTh++) {
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        matList[matOff].AX->S() +=  AXthreads[iTh][matOff]->S();  
        matList[matOff].AX->X() +=  AXthreads[iTh][matOff]->X();  
        matList[matOff].AX->Y() +=  AXthreads[iTh][matOff]->Y();  
        matList[matOff].AX->Z() +=  AXthreads[iTh][matOff]->Z();  
      }}
      
      // take care of symmetries
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        auto ADCSSMS = matList[matOff].AX->S().pointer();
        auto ADCSSMX = matList[matOff].AX->X().pointer();
        auto ADCSSMY = matList[matOff].AX->Y().pointer();
        auto ADCSSMZ = matList[matOff].AX->Z().pointer();
        for( auto i = 0; i < nBasis; i++ ) 
        for( auto j = 0; j < i; j++ ) {
          ADCSSMS[j + i*nBasis] =  ADCSSMS[i + j*nBasis];
          ADCSSMZ[j + i*nBasis] = -ADCSSMZ[i + j*nBasis];
          ADCSSMX[j + i*nBasis] = -ADCSSMX[i + j*nBasis];
          ADCSSMY[j + i*nBasis] = -ADCSSMY[i + j*nBasis];
        }
      } 

      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);

#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipSSSSAcc = std::accumulate(nIntSkipSSSS.begin(),nIntSkipSSSS.end(),0);
      double durDirectSSSSIntAcc = std::accumulate(durDirectSSSSInt.begin(), durDirectSSSSInt.end(),0.);
      double durDirectSSSSConAcc = std::accumulate(durDirectSSSSCon.begin(), durDirectSSSSCon.end(),0.);
      
      std::cout << std::endl;
      std::cout << "Dirac-Coulomb-SSSS Skipped Integral   : " << nIntSkipSSSSAcc << std::endl;
  
#ifdef _SEPARATED_SHZ_SCREEN_4C
      size_t nConSkipSSSSAcc = std::accumulate(nConSkipSSSS.begin(),nConSkipSSSS.end(),0);
      std::cout << "                         Skipped Contraction: " << nConSkipSSSSAcc << std::endl;
#endif
      auto durDirectSSSS = tock(topDirectSSSS);
      std::cout << "Dirac-Coulomb-SSSS AO Direct Contraction took " <<  durDirectSSSS
                << " s for " << nBatch << " Batch " << std::endl; 
      std::cout << "                          Build Integrals took " <<  durDirectSSSSIntAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << "                          Contractions    took " <<  durDirectSSSSConAcc << " s"
                << " from all " << nThreads << " thread(s)" << std::endl;

      std::cout << std::endl;
#endif
  
    } // if(SSSS)

    /*************************************/
    /*                                   */
    /* End of Dirac-Coulomb C(4)-SSSS    */
    /*                                   */
    /*************************************/


    /*******************/
    /*                 */
    /* Start of Gaunt  */
    /*                 */
    /*******************/

    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUNT ) {

      enum Gaunt_2ND_DERIV_AC {
        AxCx,
        AxCy,
        AxCz,
        AyCx,
        AyCy,
        AyCz,
        AzCx,
        AzCy,
        AzCz
      };
  
      enum Gaunt_2ND_DERIV_BC {
        BxCx,
        BxCy,
        BxCz,
        ByCx,
        ByCy,
        ByCz,
        BzCx,
        BzCy,
        BzCz
      };
  
      auto NB4_4  = 4*NB4;  
      auto NB4_5  = 5*NB4;  
      auto NB4_6  = 6*NB4;  
      auto NB4_7  = 7*NB4;  
      auto NB4_8  = 8*NB4;  
      auto NB4_9  = 9*NB4;  
      auto NB4_10 =10*NB4;  
      auto NB4_11 =11*NB4;  
      auto NB4_12 =12*NB4;  

#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectGaunt = tick();
      std::vector<double> durDirectGauntInt(nThreads, 0.);
      std::vector<double> durDirectGauntCon(nThreads, 0.);
#endif

      nERI = 9;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);

      int nSave = 13;
      ERIBuffer = CQMemManager::get().malloc<double>(nSave*NB4*nThreads);
      
      std::vector<size_t> nIntSkipGaunt(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C
      std::vector<size_t> nConSkipGaunt(nThreads,0);
#endif

      #pragma omp parallel
      {
  
        dcomplex iS = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *BC   = &ERIBuffer[thread_id*nSave*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        std::vector<size_t> contract_batch(nBatch);
        size_t iCon, nCon, matOff;
        
#if defined(_SHZ_SCREEN_4C) && defined(_SEPARATED_SHZ_SCREEN_4C)
        std::vector<double> shMax123_batch(nBatch);
#else
        std::iota(contract_batch.begin(), contract_batch.end(), 0);
        nCon = nBatch;
#endif
        
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C
          double shMax12 = ShBlkNorms[s2 + s1*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s2 ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s2 == s3) ? s1 : nShell-1; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C
          double shMax123 = std::max(ShBlkNorms[s1 + s3*nShell], 
                                     ShBlkNorms[s2 + s3*nShell]);
          
          shMax123 = std::max(shMax123,shMax12);
#ifdef _SEPARATED_SHZ_SCREEN_4C
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++)
             shMax123_batch[iBatch] = 
               std::max(ShBlkNorms_batch[iBatch][s1 + s2*nShell],
               std::max(ShBlkNorms_batch[iBatch][s1 + s3*nShell],
                        ShBlkNorms_batch[iBatch][s2 + s3*nShell]));
#endif
#endif
   
  
        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nThreads != thread_id ) continue;
          #endif
  
          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
            if (not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) or bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) {
              nIntSkipGaunt[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGaunt[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) { 
            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) { 
              nIntSkipGaunt[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGaunt[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) { 
            if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) ) {
              nIntSkipGaunt[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGaunt[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

#ifdef _SHZ_SCREEN_4C
  
          double shMax = std::max(ShBlkNorms[s1 + s4*nShell],
                         std::max(ShBlkNorms[s2 + s4*nShell],
                                  ShBlkNorms[s3 + s4*nShell]));
  
          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzGaunt[s3+s4*nShell]*SchwarzGaunt[s2 + s1*nShell]) <
             eri.threshSchwarz()) { 
            nIntSkipGaunt[thread_id]++;
#ifdef _SEPARATED_SHZ_SCREEN_4C
            nConSkipGaunt[thread_id] += nBatch;
#endif            
            continue; 
          }

#ifdef _SEPARATED_SHZ_SCREEN_4C
          nCon = 0;
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) {
            shMax = std::max(ShBlkNorms_batch[iBatch][s1 + s4*nShell],
                    std::max(ShBlkNorms_batch[iBatch][s2 + s4*nShell],
                             ShBlkNorms_batch[iBatch][s3 + s4*nShell]));
            shMax = std::max(shMax, shMax123_batch[iBatch]); 
             
            if((shMax*SchwarzGaunt[s2+s1*nShell]*SchwarzGaunt[s3+s4*nShell]) <
              eri.threshSchwarz()) {
              nConSkipGaunt[thread_id] ++;
            } else {
              contract_batch[nCon] = iBatch; 
              nCon++;  
            }
          }
#endif
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectGauntInt = tick();
#endif
  
          shls[0] = int(s2);
          shls[1] = int(s1);
          shls[2] = int(s3);
          shls[3] = int(s4);

          if(int2e_ip1ip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
 
          auto nQuad = n1*n2*n3*n4;

          for(l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l) 
          for(k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k) 
          for(m = 0                         ; m <                  n1; ++m)
          for(n =   maxShellSize            ; n <   maxShellSize + n2; ++n, ++mnkl) {

  
            /* Gaunt */
            // ∇A∙∇C(mn|kl)
            auto dAdotdC = buff[AxCx*nQuad+mnkl] + buff[AyCy*nQuad+mnkl] + buff[AzCz*nQuad+mnkl];
            // ∇Ax∇C(mn|kl)
            auto dAcrossdC_x =  buff[AyCz*nQuad+mnkl] - buff[AzCy*nQuad+mnkl];
            auto dAcrossdC_y = -buff[AxCz*nQuad+mnkl] + buff[AzCx*nQuad+mnkl];
            auto dAcrossdC_z =  buff[AxCy*nQuad+mnkl] - buff[AyCx*nQuad+mnkl];


            // Change the index so that we do ∇B∙∇C(ij|kl) using the ∇B∇C engine
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;
  
            // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
            // (kl|mn)
            BC[      MNKL] = dAdotdC;    
            BC[  NB4+MNKL] = dAcrossdC_x;
            BC[NB4_2+MNKL] = dAcrossdC_y;
            BC[NB4_3+MNKL] = dAcrossdC_z;

            // (lk|nm)
            BC[      LKNM] = dAdotdC;
            BC[  NB4+LKNM] = -dAcrossdC_x;
            BC[NB4_2+LKNM] = -dAcrossdC_y;
            BC[NB4_3+LKNM] = -dAcrossdC_z;
  
            // ∇B_x∇C_x(mn|kl) - ∇B∙∇C(mn|kl)
            // (mn|kl)
            BC[NB4_4+MNKL] = buff[BxCx*nQuad+mnkl] - dAdotdC;
            // (lk|nm)
            BC[NB4_4+LKNM] = buff[BxCx*nQuad+mnkl] - dAdotdC;
  
  
            // ∇B_y∇C_x(mn|kl)
            // (mn|kl)
            BC[NB4_5+MNKL] = buff[ByCx*nQuad+mnkl];
            // (lk|nm)
            BC[NB4_5+LKNM] = buff[BxCy*nQuad+mnkl];
  
  
            // ∇B_z∇C_x(mn|kl)
            // (mn|kl)
            BC[NB4_6+MNKL] = buff[BzCx*nQuad+mnkl];
            // (lk|nm)
            BC[NB4_6+LKNM] = buff[BxCz*nQuad+mnkl];
  
  
            // ∇B_x∇C_y(mn|kl)
            // (mn|kl)
            BC[NB4_7+MNKL] = buff[BxCy*nQuad+mnkl]; 
            // (lk|nm)
            BC[NB4_7+LKNM] = buff[ByCx*nQuad+mnkl];
  
  
            // ∇B_y∇C_y(mn|kl) - ∇B∙∇C(mn|kl)
            // (mn|kl)
            BC[NB4_8+MNKL] = buff[ByCy*nQuad+mnkl] - dAdotdC;
            // (lk|nm)
            BC[NB4_8+LKNM] = buff[ByCy*nQuad+mnkl] - dAdotdC;
  
  
            // ∇B_z∇C_y(mn|kl)
            // (mn|kl)
            BC[NB4_9+MNKL] = buff[BzCy*nQuad+mnkl];
            // (lk|nm)
            BC[NB4_9+LKNM] = buff[ByCz*nQuad+mnkl];
  
            // ∇B_x∇C_z(mn|kl)
            // (mn|kl)
            BC[NB4_10+MNKL] = buff[BxCz*nQuad+mnkl];
            // (lk|nm)
            BC[NB4_10+LKNM] = buff[BzCx*nQuad+mnkl];
  
            // ∇B_y∇C_z(mn|kl)
            // (mn|kl)
            BC[NB4_11+MNKL] = buff[ByCz*nQuad+mnkl];
            // (lk|nm)
            BC[NB4_11+LKNM] = buff[BzCy*nQuad+mnkl];
  
            // ∇B_z∇C_z(mn|kl) - ∇B∙∇C(mn|kl) 
            // (mn|kl)
            BC[NB4_12+MNKL] = buff[BzCz*nQuad+mnkl] - dAdotdC;
            // (lk|nm)
            BC[NB4_12+LKNM] = buff[BzCz*nQuad+mnkl] - dAdotdC;
  
          } // ∇C∙∇D integral preparation loop
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          durDirectGauntInt[thread_id] += tock(topDirectGauntInt); 
#endif
  
#ifdef _CONTRACTION_ // Contraction

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectGauntCon = tick();
#endif
          for (auto iCon = 0ul; iCon < nCon; iCon++) {  
            
            matOff = contract_batch[iCon] * mMat; 
            
            auto DCLSpmSLMS = matList[matOff].X->S().pointer();
            auto DCLSpmSLMX = matList[matOff].X->X().pointer();
            auto DCLSpmSLMY = matList[matOff].X->Y().pointer();
            auto DCLSpmSLMZ = matList[matOff].X->Z().pointer();
            
            auto ADCLSMS  = AX_loc[matOff]->S().pointer();
            auto ADCLSMX  = AX_loc[matOff]->X().pointer();
            auto ADCLSMY  = AX_loc[matOff]->Y().pointer();
            auto ADCLSMZ  = AX_loc[matOff]->Z().pointer();
            
          for(m = 0ul,            bf1 = bf1_s; m <                  n1; ++m, bf1++) {
          for(n =   maxShellSize, bf2 = bf2_s; n <   maxShellSize + n2; ++n, bf2++) {
            auto bf1nB = bf1*nBasis;
            auto bf2nB = bf2*nBasis;

          for(k = 2*maxShellSize, bf3 = bf3_s; k < 2*maxShellSize + n3; ++k, bf3++) {
            auto bf3nB = bf3*nBasis;

          for(l = 3*maxShellSize, bf4 = bf4_s; l < 3*maxShellSize + n4; ++l, bf4++) {
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;
            auto NMLK = n + m*NB + l*NB2 + k*NB3;
            auto KLMN = k + l*NB + m*NB2 + n*NB3;

            auto bf4nB= bf4*nBasis;
            auto bf14 = bf1 + bf4nB;
            auto bf24 = bf2 + bf4nB;
            auto bf23 = bf2 + bf3nB;
            auto bf13 = bf1 + bf3nB;
            auto bf41 = bf4 + bf1nB;
            auto bf31 = bf3 + bf1nB;
            auto bf32 = bf3 + bf2nB;
            auto bf42 = bf4 + bf2nB;
            auto bf34 = bf3 + bf4nB;
            auto bf43 = bf4 + bf3nB;
            auto bf12 = bf1 + bf2nB;
            auto bf21 = bf2 + bf1nB;

            //ERI0 :   ∇B∙∇C(mn|kl)
            //ERI1 :   ∇Bx∇C(mn|kl)-X
            //ERI2 :   ∇Bx∇C(mn|kl)-Y
            //ERI3 :   ∇Bx∇C(mn|kl)-Z
            //ERI4 :   ∇B_x∇C_x(mn|kl) - ∇B∙∇C(mn|kl)
            //ERI5 :   ∇B_y∇C_x(mn|kl)
            //ERI6 :   ∇B_z∇C_x(mn|kl) 
            //ERI7 :   ∇B_x∇C_y(mn|kl)
            //ERI8 :   ∇B_y∇C_y(mn|kl) - ∇B∙∇C(mn|kl)
            //ERI9 :   ∇B_z∇C_y(mn|kl) 
            //ERI10:   ∇B_x∇C_z(mn|kl)
            //ERI11:   ∇B_y∇C_z(mn|kl)
            //ERI12:   ∇B_z∇C_z(mn|kl) - ∇B∙∇C(mn|kl)
            
            /*++++++++++++++++++++++++*/
            /* Start of Gaunt (LL|SS) */
            /*++++++++++++++++++++++++*/
            
            // MNKL
            // Coulomb //TRANSKL

            ADCLSMS[bf12]+=    - DCLSpmSLMS[bf43] * BC[MNKL]
                            + (  DCLSpmSLMX[bf43] * BC[MNKL+NB4] 
                               + DCLSpmSLMY[bf43] * BC[MNKL+NB4_2] 
                               + DCLSpmSLMZ[bf43] * BC[MNKL+NB4_3]) * iS;
            
            ADCLSMX[bf12]+=  DCLSpmSLMS[bf43] * BC[MNKL+NB4] * iS 
                           + DCLSpmSLMX[bf43] * BC[MNKL+NB4_4]
                           + DCLSpmSLMY[bf43] * BC[MNKL+NB4_5] 
                           + DCLSpmSLMZ[bf43] * BC[MNKL+NB4_6];
            
            ADCLSMY[bf12]+=  DCLSpmSLMS[bf43] * BC[MNKL+NB4_2] * iS 
                           + DCLSpmSLMX[bf43] * BC[MNKL+NB4_7] 
                           + DCLSpmSLMY[bf43] * BC[MNKL+NB4_8]
                           + DCLSpmSLMZ[bf43] * BC[MNKL+NB4_9];
            
            ADCLSMZ[bf12]+=  DCLSpmSLMS[bf43] * BC[MNKL+NB4_3] * iS 
                           + DCLSpmSLMX[bf43] * BC[MNKL+NB4_10]
                           + DCLSpmSLMY[bf43] * BC[MNKL+NB4_11] 
                           + DCLSpmSLMZ[bf43] * BC[MNKL+NB4_12];
            // LKNM
            // Coulomb //TRANSKL
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
            
              ADCLSMS[bf43]+=   - DCLSpmSLMS[bf12] * BC[LKNM] 
                              +(  DCLSpmSLMX[bf12] * BC[LKNM+NB4] 
                                + DCLSpmSLMY[bf12] * BC[LKNM+NB4_2] 
                                + DCLSpmSLMZ[bf12] * BC[LKNM+NB4_3]) * iS;
            
              ADCLSMX[bf43]+=  DCLSpmSLMS[bf12] * BC[LKNM+NB4] * iS 
                             + DCLSpmSLMX[bf12] * BC[LKNM+NB4_4]
                             + DCLSpmSLMY[bf12] * BC[LKNM+NB4_5] 
                             + DCLSpmSLMZ[bf12] * BC[LKNM+NB4_6];
              
              ADCLSMY[bf43]+=  DCLSpmSLMS[bf12] * BC[LKNM+NB4_2] * iS 
                             + DCLSpmSLMX[bf12] * BC[LKNM+NB4_7] 
                             + DCLSpmSLMY[bf12] * BC[LKNM+NB4_8]
                             + DCLSpmSLMZ[bf12] * BC[LKNM+NB4_9];
            
              ADCLSMZ[bf43]+=  DCLSpmSLMS[bf12] * BC[LKNM+NB4_3] * iS 
                             + DCLSpmSLMX[bf12] * BC[LKNM+NB4_10]
                             + DCLSpmSLMY[bf12] * BC[LKNM+NB4_11]
                             + DCLSpmSLMZ[bf12] * BC[LKNM+NB4_12];
            
            }

          }
          }
          }
          }
          } // contraction loop
#endif // Contraction
  
#ifdef _REPORT_INTEGRAL_TIMINGS
        durDirectGauntCon[thread_id] += tock(topDirectGauntCon); 
#endif
  
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1
  
  
      } // OpenMP context
  
  
      for( auto iTh  = 0; iTh < nThreads; iTh++) 
      for (auto iMat = 0ul; iMat < nMat; iMat++) {  
        matList[iMat].AX->S() +=  AXthreads[iTh][iMat]->S();  
        matList[iMat].AX->X() +=  AXthreads[iTh][iMat]->X();  
        matList[iMat].AX->Y() +=  AXthreads[iTh][iMat]->Y();  
        matList[iMat].AX->Z() +=  AXthreads[iTh][iMat]->Z();  
      }
      
      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);

#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipGauntAcc = std::accumulate(nIntSkipGaunt.begin(),nIntSkipGaunt.end(),0);
      double durDirectGauntIntAcc = std::accumulate(durDirectGauntInt.begin(), durDirectGauntInt.end(),0.);
      double durDirectGauntConAcc = std::accumulate(durDirectGauntCon.begin(), durDirectGauntCon.end(),0.);
      
      std::cout << std::endl;
      std::cout << "Gaunt Skipped Integral   : " << nIntSkipGauntAcc << std::endl;
  
#ifdef _SEPARATED_SHZ_SCREEN_4C
      size_t nConSkipGauntAcc = std::accumulate(nConSkipGaunt.begin(),nConSkipGaunt.end(),0);
      std::cout << "       Skipped Contraction: " << nConSkipGauntAcc << std::endl;
#endif
      auto durDirectGaunt = tock(topDirectGaunt);
      std::cout << "Gaunt AO Direct Contraction took " <<  durDirectGaunt
                << " s for " << nBatch << " Batch " << std::endl; 
      std::cout << "            Build Integrals took " <<  durDirectGauntIntAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << "            Contractions    took " <<  durDirectGauntConAcc << " s"
                << " from all " << nThreads << " thread(s)" << std::endl;

      std::cout << std::endl;
#endif

    } // if(Gaunt)

    /*******************/
    /*                 */
    /*   End of Gaunt  */
    /*                 */
    /*******************/


    /*******************/
    /*                 */
    /* Start of Gauge  */
    /*                 */
    /*******************/

    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUGE ) {

      enum Gauge_ERI {
        SxSx, // σ_x * σ_x     0
        SySx, // σ_y * σ_x     1
        SzSx, // σ_z * σ_x     2
        ISx,  // I   * σ_x     3
        SxSy, // σ_x * σ_y     4
        SySy, // σ_y * σ_y     5
        SzSy, // σ_z * σ_y     6
        ISy,  // I   * σ_y     7
        SxSz, // σ_x * σ_z     8
        SySz, // σ_y * σ_z     9
        SzSz, // σ_z * σ_z     10
        ISz,  // I   * σ_z     11
        SxI,  // σ_x * I       12
        SyI,  // σ_y * I       13
        SzI,  // σ_z * I       14
        II    // I   * I       15
      };
  
      auto NB4_4  = 4*NB4;  
      auto NB4_5  = 5*NB4;  
      auto NB4_6  = 6*NB4;  
      auto NB4_7  = 7*NB4;  
      auto NB4_8  = 8*NB4;  
      auto NB4_9  = 9*NB4;  
      auto NB4_10 =10*NB4;  
      auto NB4_11 =11*NB4;  
      auto NB4_12 =12*NB4;  
      auto NB4_13 =13*NB4;  
      auto NB4_14 =14*NB4;  
      auto NB4_15 =15*NB4;  

      auto ss = 0     ;    // ERI00 (ss)(ij|kl)
      auto sx = NB4   ;    // ERI01 (sσ)_x(ijkl)
      auto sy = NB4_2 ;    // ERI02 (sσ)_y(ijkl)
      auto sz = NB4_3 ;    // ERI03 (sσ)_z(ijkl)
      auto xs = NB4_4 ;    // ERI04 (σs)_x(ijkl)
      auto ys = NB4_5 ;    // ERI05 (σs)_y(ijkl)
      auto zs = NB4_6 ;    // ERI06 (σs)_z(ijkl)
      auto xx = NB4_7;     // ERI07 (σ_x σ_x)(ijkl)
      auto xy = NB4_8;     // ERI08 (σ_x σ_y)(ijkl)
      auto xz = NB4_9;     // ERI09 (σ_x σ_z)(ijkl)
      auto yx = NB4_10;    // ERI10 (σ_y σ_x)(ijkl)
      auto yy = NB4_11;    // ERI11 (σ_y σ_y)(ijkl)
      auto yz = NB4_12;    // ERI12 (σ_y σ_z)(ijkl)
      auto zx = NB4_13;    // ERI13 (σ_z σ_x)(ijkl)
      auto zy = NB4_14;    // ERI14 (σ_z σ_y)(ijkl)
      auto zz = NB4_15;    // ERI15 (σ_z σ_z)(ijkl)
   
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectGauge = tick();
      std::vector<double> durDirectGaugeInt(nThreads, 0.);
      std::vector<double> durDirectGaugeCon(nThreads, 0.);
#endif

      enum X_COMP {
        DLS,
        DSL,
      };

      enum AX_COMP {
        CLS,
        CSL,
      };
      
      nERI = 16;
      buffAll = CQMemManager::get().malloc<double>(2*nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);

      int nSave = 16;
      ERIBuffer = CQMemManager::get().malloc<double>(nSave*NB4*nThreads);
      
      std::vector<size_t> nIntSkipGauge(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C
      std::vector<size_t> nConSkipGauge(nThreads,0);
#endif
  
      #pragma omp parallel
      {
  
        dcomplex iS = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *BC   = &ERIBuffer[thread_id*nSave*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;

        int skiperi1,skiperi2;
  
        int shls[4];
        double *buffr1 = buffAll+nERI*buffN4*thread_id;
        double *buffr2 = buffAll+nERI*nThreads*buffN4+nERI*buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
  
        std::vector<size_t> contract_batch(nBatch);
        size_t iCon, nCon, matOff;

#if defined(_SHZ_SCREEN_4C) && defined(_SEPARATED_SHZ_SCREEN_4C)
        std::vector<double> shMax123_batch(nBatch);
#else
        std::iota(contract_batch.begin(), contract_batch.end(), 0);
        nCon = nBatch;
#endif
        
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C
          double shMax12 = ShBlkNorms[s2 + s1*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s2 ; bf3_s+=n3, s3++) {
        //for(size_t s3(0), bf3_s(0); s3 < nShell ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s2 == s3) ? s1 : nShell-1; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C
          double shMax123 = std::max(ShBlkNorms[s1 + s3*nShell], 
                                     ShBlkNorms[s2 + s3*nShell]);
          
          shMax123 = std::max(shMax123,shMax12);
#ifdef _SEPARATED_SHZ_SCREEN_4C
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++)
             shMax123_batch[iBatch] = 
               std::max(ShBlkNorms_batch[iBatch][s1 + s2*nShell],
               std::max(ShBlkNorms_batch[iBatch][s1 + s3*nShell],
                        ShBlkNorms_batch[iBatch][s2 + s3*nShell]));
#endif
#endif
   
  
        //for(size_t s4(0), bf4_s(0); s4 < nShell; bf4_s+=n4, s4++, s1234++) {
        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nThreads != thread_id ) continue;
          #endif
  
          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
            if (not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) or bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) {
              nIntSkipGauge[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGauge[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) { 
            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4))) { 
              nIntSkipGauge[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGauge[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) { 
            if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) ) {
              nIntSkipGauge[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN_4C
              nConSkipGauge[thread_id] += nBatch; 
#endif            
              continue;
            }
          }

#ifdef _SHZ_SCREEN_4C
  
          double shMax = std::max(ShBlkNorms[s1 + s4*nShell],
                         std::max(ShBlkNorms[s2 + s4*nShell],
                                  ShBlkNorms[s3 + s4*nShell]));
  
          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzGauge[s1+s2*nShell]*SchwarzGauge[s4+s3*nShell]) <
             eri.threshSchwarz()) { 
            nIntSkipGauge[thread_id]++;
#ifdef _SEPARATED_SHZ_SCREEN_4C
            nConSkipGauge[thread_id] += nBatch;
#endif            
            continue; 
          }

#ifdef _SEPARATED_SHZ_SCREEN_4C
          nCon = 0;
          for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) {
            shMax = std::max(ShBlkNorms_batch[iBatch][s1 + s4*nShell],
                    std::max(ShBlkNorms_batch[iBatch][s2 + s4*nShell],
                             ShBlkNorms_batch[iBatch][s3 + s4*nShell]));
            shMax = std::max(shMax, shMax123_batch[iBatch]); 
             
            if((shMax*SchwarzGauge[s1+s2*nShell]*SchwarzGauge[s4+s3*nShell]) <
              eri.threshSchwarz()) {
              nConSkipGauge[thread_id] ++;
            } else {
              contract_batch[nCon] = iBatch; 
              nCon++;  
            }
          }
#endif
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectGaugeInt = tick();
#endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);

          //∇B∇C
          skiperi1 = int2e_gauge_r1_ssp1sps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi2 = int2e_gauge_r2_ssp1sps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);

          if(skiperi1==0 and skiperi2==0) continue;

          auto nQuad = n1*n2*n3*n4;

	  mnkl = 0ul;
          for(l = 3*maxShellSize; l < 3*maxShellSize + n4; ++l) 
          for(k = 2*maxShellSize; k < 2*maxShellSize + n3; ++k) 
          for(n =   maxShellSize; n <   maxShellSize + n2; ++n)
          for(m = 0             ; m <                  n1; ++m) {

            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;

            //enum Gauge_ERI {
            //  SxSx,  // σ_x * σ_x     0
            //  SySx,  // σ_y * σ_x     1
            //  SzSx,  // σ_z * σ_x     2
            //  ISx,   // I   * σ_x     3
            //  SxSy,  // σ_x * σ_y     4
            //  SySy,  // σ_y * σ_y     5
            //  SzSy,  // σ_z * σ_y     6
            //  ISy,   // I   * σ_y     7
            //  SxSz,  // σ_x * σ_z     8
            //  SySz,  // σ_y * σ_z     9
            //  SzSz,  // σ_z * σ_z     10
            //  ISz,   // I   * σ_z     11
            //  SxI,   // σ_x * I       12
            //  SyI,   // σ_y * I       13
            //  SzI,   // σ_z * I       14
            //  II     // I   * I       15
            //};
            
            // (ss)
            BC[MNKL] = buffr1[15*nQuad+mnkl] - buffr2[15*nQuad+mnkl];
            BC[LKNM] = buffr1[15*nQuad+mnkl] - buffr2[15*nQuad+mnkl];
            
            // (sσ)_x
            BC[NB4+MNKL] = buffr1[3*nQuad+mnkl] - buffr2[3*nQuad+mnkl];
            BC[NB4+LKNM] = buffr2[12*nQuad+mnkl] - buffr1[12*nQuad+mnkl];
            
            // (sσ)_y
            BC[NB4_2+MNKL] = buffr1[7*nQuad+mnkl] - buffr2[7*nQuad+mnkl];
            BC[NB4_2+LKNM] = buffr2[13*nQuad+mnkl] - buffr1[13*nQuad+mnkl];
            
            // (sσ)_z
            BC[NB4_3+MNKL] = buffr1[11*nQuad+mnkl] - buffr2[11*nQuad+mnkl];
            BC[NB4_3+LKNM] = buffr2[14*nQuad+mnkl] - buffr1[14*nQuad+mnkl];
            
            // (σs)_x
            BC[NB4_4+MNKL] = buffr1[12*nQuad+mnkl] - buffr2[12*nQuad+mnkl];
            BC[NB4_4+LKNM] = buffr2[3*nQuad+mnkl] - buffr1[3*nQuad+mnkl];
            
            // (σs)_y
            BC[NB4_5+MNKL] = buffr1[13*nQuad+mnkl] - buffr2[13*nQuad+mnkl];
            BC[NB4_5+LKNM] = buffr2[7*nQuad+mnkl] - buffr1[7*nQuad+mnkl];
            
            // (σs)_z
            BC[NB4_6+MNKL] = buffr1[14*nQuad+mnkl] - buffr2[14*nQuad+mnkl];
            BC[NB4_6+LKNM] = buffr2[11*nQuad+mnkl] - buffr1[11*nQuad+mnkl];
            
            // σ_x*σ_x
            BC[NB4_7+MNKL] = -(buffr1[mnkl] - buffr2[mnkl]);
            BC[NB4_7+LKNM] = -(buffr1[mnkl] - buffr2[mnkl]);
            
            // σ_x*σ_y
            BC[NB4_8+MNKL] = -(buffr1[4*nQuad+mnkl] - buffr2[4*nQuad+mnkl]);
            BC[NB4_8+LKNM] = -(buffr1[1*nQuad+mnkl] - buffr2[1*nQuad+mnkl]);
            
            // σ_x*σ_z
            BC[NB4_9+MNKL] = -(buffr1[8*nQuad+mnkl] - buffr2[8*nQuad+mnkl]);
            BC[NB4_9+LKNM] = -(buffr1[2*nQuad+mnkl] - buffr2[2*nQuad+mnkl]);
            
            // σ_y*σ_x
            BC[NB4_10+MNKL] = -(buffr1[1*nQuad+mnkl] - buffr2[1*nQuad+mnkl]);
            BC[NB4_10+LKNM] = -(buffr1[4*nQuad+mnkl] - buffr2[4*nQuad+mnkl]);
            
            // σ_y*σ_y
            BC[NB4_11+MNKL] = -(buffr1[5*nQuad+mnkl] - buffr2[5*nQuad+mnkl]);
            BC[NB4_11+LKNM] = -(buffr1[5*nQuad+mnkl] - buffr2[5*nQuad+mnkl]);
            
            // σ_y*σ_z
            BC[NB4_12+MNKL] = -(buffr1[9*nQuad+mnkl] - buffr2[9*nQuad+mnkl]);
            BC[NB4_12+LKNM] = -(buffr1[6*nQuad+mnkl] - buffr2[6*nQuad+mnkl]);
            
            // σ_z*σ_x
            BC[NB4_13+MNKL] = -(buffr1[2*nQuad+mnkl] - buffr2[2*nQuad+mnkl]);
            BC[NB4_13+LKNM] = -(buffr1[8*nQuad+mnkl] - buffr2[8*nQuad+mnkl]);
            
            // σ_z*σ_y
            BC[NB4_14+MNKL] = -(buffr1[6*nQuad+mnkl] - buffr2[6*nQuad+mnkl]);
            BC[NB4_14+LKNM] = -(buffr1[9*nQuad+mnkl] - buffr2[9*nQuad+mnkl]);
            
            // σ_z*σ_z
            BC[NB4_15+MNKL] = -(buffr1[10*nQuad+mnkl] - buffr2[10*nQuad+mnkl]);
            BC[NB4_15+LKNM] = -(buffr1[10*nQuad+mnkl] - buffr2[10*nQuad+mnkl]);
            
            mnkl++;

          } // ∇B∇C integral preparation loop
  
#ifdef _REPORT_INTEGRAL_TIMINGS
          durDirectGaugeInt[thread_id] += tock(topDirectGaugeInt); 
#endif
  
#ifdef _CONTRACTION_ // Contraction

#ifdef _REPORT_INTEGRAL_TIMINGS
          auto topDirectGaugeCon = tick();
#endif
          for (auto iCon = 0ul; iCon < nCon; iCon++) {  
            
            matOff = contract_batch[iCon] * mMat; 
            
            auto DCLSpmSLMS = matList[matOff].X->S().pointer();
            auto DCLSpmSLMX = matList[matOff].X->X().pointer();
            auto DCLSpmSLMY = matList[matOff].X->Y().pointer();
            auto DCLSpmSLMZ = matList[matOff].X->Z().pointer();
            
            auto ADCLSMS  = AX_loc[matOff]->S().pointer();
            auto ADCLSMX  = AX_loc[matOff]->X().pointer();
            auto ADCLSMY  = AX_loc[matOff]->Y().pointer();
            auto ADCLSMZ  = AX_loc[matOff]->Z().pointer();
 
          for(m = 0ul,            bf1 = bf1_s; m <                  n1; ++m, bf1++) {
          for(n =   maxShellSize, bf2 = bf2_s; n <   maxShellSize + n2; ++n, bf2++) {
            auto bf1nB = bf1*nBasis;
            auto bf2nB = bf2*nBasis;

          for(k = 2*maxShellSize, bf3 = bf3_s; k < 2*maxShellSize + n3; ++k, bf3++) {
            auto bf3nB = bf3*nBasis;

          for(l = 3*maxShellSize, bf4 = bf4_s; l < 3*maxShellSize + n4; ++l, bf4++) {
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;

            auto bf4nB= bf4*nBasis;
            auto bf14 = bf1 + bf4nB;
            auto bf24 = bf2 + bf4nB;
            auto bf23 = bf2 + bf3nB;
            auto bf13 = bf1 + bf3nB;
            auto bf41 = bf4 + bf1nB;
            auto bf31 = bf3 + bf1nB;
            auto bf32 = bf3 + bf2nB;
            auto bf42 = bf4 + bf2nB;
            auto bf34 = bf3 + bf4nB;
            auto bf43 = bf4 + bf3nB;
            auto bf12 = bf1 + bf2nB;
            auto bf21 = bf2 + bf1nB;

            /*++++++++++++++++++++++++*/
            /* Start of Gauge (LL|SS) */
            /*++++++++++++++++++++++++*/
             
            
            // Coulomb
            // MNKL
            
            /* Equations (232) and (233) */
            ADCLSMS[bf12]+=    - DCLSpmSLMS[bf43] * BC[MNKL]
                            - (  DCLSpmSLMX[bf43] * BC[MNKL+sx] 
                               + DCLSpmSLMY[bf43] * BC[MNKL+sy] 
                               + DCLSpmSLMZ[bf43] * BC[MNKL+sz]) * iS;
            
            ADCLSMX[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+xs] * iS 
                           - DCLSpmSLMX[bf43] * BC[MNKL+xx]
                           - DCLSpmSLMY[bf43] * BC[MNKL+xy]
                           - DCLSpmSLMZ[bf43] * BC[MNKL+xz]; 
            
            ADCLSMY[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+ys] * iS 
                           - DCLSpmSLMX[bf43] * BC[MNKL+yx]
                           - DCLSpmSLMY[bf43] * BC[MNKL+yy]
                           - DCLSpmSLMZ[bf43] * BC[MNKL+yz]; 
            
            ADCLSMZ[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+zs] * iS 
                           - DCLSpmSLMX[bf43] * BC[MNKL+zx]
                           - DCLSpmSLMY[bf43] * BC[MNKL+zy]
                           - DCLSpmSLMZ[bf43] * BC[MNKL+zz]; 
            
            // LKNM
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              
              ADCLSMS[bf43]+=    - DCLSpmSLMS[bf12] * BC[LKNM]
                              - (  DCLSpmSLMX[bf12] * BC[LKNM+sx] 
                                 + DCLSpmSLMY[bf12] * BC[LKNM+sy] 
                                 + DCLSpmSLMZ[bf12] * BC[LKNM+sz]) * iS;
              
              ADCLSMX[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+xs] * iS 
                             - DCLSpmSLMX[bf12] * BC[LKNM+xx]
                             - DCLSpmSLMY[bf12] * BC[LKNM+xy]
                             - DCLSpmSLMZ[bf12] * BC[LKNM+xz]; 
              
              ADCLSMY[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+ys] * iS 
                             - DCLSpmSLMX[bf12] * BC[LKNM+yx]
                             - DCLSpmSLMY[bf12] * BC[LKNM+yy]
                             - DCLSpmSLMZ[bf12] * BC[LKNM+yz]; 
              
              ADCLSMZ[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+zs] * iS 
                             - DCLSpmSLMX[bf12] * BC[LKNM+zx]
                             - DCLSpmSLMY[bf12] * BC[LKNM+zy]
                             - DCLSpmSLMZ[bf12] * BC[LKNM+zz]; 
              
            }
             
          }
          }
          }
          }
          } // contraction loop
  
#endif // Contraction
  
#ifdef _REPORT_INTEGRAL_TIMINGS
        durDirectGaugeCon[thread_id] += tock(topDirectGaugeCon); 
#endif
  
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1
  
  
      } // OpenMP context
  
  
      for( auto iTh  = 0; iTh < nThreads; iTh++) 
      for (auto iMat = 0ul; iMat < nMat; iMat++) {  
        matList[iMat].AX->S() +=  AXthreads[iTh][iMat]->S();  
        matList[iMat].AX->X() +=  AXthreads[iTh][iMat]->X();  
        matList[iMat].AX->Y() +=  AXthreads[iTh][iMat]->Y();  
        matList[iMat].AX->Z() +=  AXthreads[iTh][iMat]->Z();  
      }

      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipGaugeAcc = std::accumulate(nIntSkipGauge.begin(),nIntSkipGauge.end(),0);
      double durDirectGaugeIntAcc = std::accumulate(durDirectGaugeInt.begin(), durDirectGaugeInt.end(),0.);
      double durDirectGaugeConAcc = std::accumulate(durDirectGaugeCon.begin(), durDirectGaugeCon.end(),0.);
      
      std::cout << std::endl;
      std::cout << "Gauge Skipped Integral   : " << nIntSkipGaugeAcc << std::endl;
  
#ifdef _SEPARATED_SHZ_SCREEN_4C
      size_t nConSkipGaugeAcc = std::accumulate(nConSkipGauge.begin(),nConSkipGauge.end(),0);
      std::cout << "      Skipped Contraction: " << nConSkipGaugeAcc << std::endl;
#endif
      auto durDirectGauge = tock(topDirectGauge);
      std::cout << "Gauge AO Direct Contraction took " <<  durDirectGauge
                << " s for " << nBatch << " Batch " << std::endl; 
      std::cout << "            Build Integrals took " <<  durDirectGaugeIntAcc << " s" 
                << " from all " << nThreads << " thread(s)" << std::endl;
      std::cout << "            Contractions    took " <<  durDirectGaugeConAcc << " s"
                << " from all " << nThreads << " thread(s)" << std::endl;

      std::cout << std::endl;
#endif

    } // if(Gauge)

    /*******************/
    /*                 */
    /*   End of Gaunt  */
    /*                 */
    /*******************/


#ifdef _SHZ_SCREEN_4C
    if(ShBlkNorms_raw) CQMemManager::get().free(ShBlkNorms_raw);
    if(SchwarzSSSS)    CQMemManager::get().free(SchwarzSSSS);
    if(SchwarzERI)     CQMemManager::get().free(SchwarzERI);
    if(SchwarzGaunt)   CQMemManager::get().free(SchwarzGaunt);
    if(SchwarzGauge)   CQMemManager::get().free(SchwarzGauge);
#endif
  
    CQMemManager::get().free(atm);
    CQMemManager::get().free(bas);
    CQMemManager::get().free(env);
  
    // Turn threads for LA back on
    // SetLAThreads(LAThreads);
  
  }; // GTODirectRelERIContraction::directRelScaffoldLibcintCoulombOnly
  
  template <>
  void GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintCoulombOnly(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyRelContraction<double>> &matList,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);  
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintCoulombOnly(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyRelContraction<dcomplex>> &matList,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);  
  }
  
  template <typename MatsT, typename IntsT>
  size_t GTODirectRelERIContraction<MatsT,IntsT>::libcintCacheSize(
    const TWOBODY_CONTRACTION_TYPE & contType, int * atm, const int nAtoms, 
    int * bas, const int nShells, double * env) const {
    
    size_t cache_size = 0ul;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      if(contType == TWOBODY_CONTRACTION_TYPE::BARE_COULOMB or
         contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(contType == TWOBODY_CONTRACTION_TYPE::SSSS or
         contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_ipvip1ipvip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_ipvip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(contType == TWOBODY_CONTRACTION_TYPE::GAUNT) {
        n = int2e_ip1ip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(contType == TWOBODY_CONTRACTION_TYPE::GAUGE) {
        n = int2e_gauge_r1_ssp1sps2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        n = int2e_gauge_r2_ssp1sps2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        n = int2e_gauge_r1_ssp1ssp2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        n = int2e_gauge_r2_ssp1ssp2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
    }

    return cache_size; 
  }

  /*******************************************************************************/
  /*                                                                             */
  /* Compute memory requirement for Libcint Batch 4C-direct                      */
  /* Returns:                                                                    */
  /*   size_t SCR size needed for one batch                                      */
  /*   IMPORTANT HERE: size are all in MatsT(dcomplex)                           */
  /*******************************************************************************/
  
  template <typename MatsT, typename IntsT>
  size_t GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintSCRSize(
    const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {
    
    // disable NYI part
    if (computeExchange or contType == TWOBODY_CONTRACTION_TYPE::LLSS)
      CErr("NYI contType in directRelScaffoldSCRSize"); 
    
    size_t threadSCRSize  = 0ul;
    size_t generalSCRSize = 0ul;
  
    DirectTPI<IntsT> &originalERI = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& originalBasisSet_ = originalERI.basisSet();
    Molecule& molecule_ = originalERI.molecule();
    
    if (originalBasisSet_.forceCart)
      CErr("Libcint + cartesian GTO NYI.");
    
    BasisSet  basisSet_ = originalBasisSet_.groupGeneralContractionBasis();
    
  
    int nAtoms      = molecule_.nAtoms;
    int nShells     = basisSet_.nShell;
    size_t nBasis   = basisSet_.nBasis;
  
    size_t basisSCRSize = nAtoms * ATM_SLOTS
                         + nShells * BAS_SLOTS
                         + basisSet_.getLibcintEnvLength(molecule_);
  
    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();
    
    int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().malloc<double>(basisSet_.getLibcintEnvLength(molecule_));
    
    basisSet_.setLibcintEnv(molecule_, atm, bas, env);
    
    size_t cacheSize = libcintCacheSize(contType, atm, nAtoms, bas, nShells, env);
    
    CQMemManager::get().free(atm, bas, env);
    
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    
    int nERI, nERIbuff, nSHZ, mRawMat;
    size_t maxShellSize = buffSize;
    size_t NB  = maxShellSize*4;
    size_t NB4 = NB*NB*NB*NB;
    
    // Coulomb Only
    if (not computeExchange) {
      if (contType == TWOBODY_CONTRACTION_TYPE::BARE_COULOMB) {
        mRawMat = 1;
        nERI = 1;
        nERIbuff = 0;
      } else if (contType == TWOBODY_CONTRACTION_TYPE::LLLL) {
        mRawMat = 5;
        nERI = 9;
        nERIbuff = 2*4;
      } else if (contType == TWOBODY_CONTRACTION_TYPE::SSSS) {
        mRawMat = 4;
        nERI = 81;
        nERIbuff = 16;
      } else if (contType == TWOBODY_CONTRACTION_TYPE::GAUNT) {
        mRawMat = 4;
        nERI = 9;
        nERIbuff = 13;
      } else if (contType == TWOBODY_CONTRACTION_TYPE::GAUGE) {
        mRawMat = 4;
        nERI = 16*2;
        nERIbuff = 16;
      } else {
        std::cout << "contType = " << contType << std::endl;
        CErr("NYI contType in directRelScaffoldSCRSize"); 
      }
    } else{
      CErr("computeExchange NYI in directRelScaffoldSCRSize"); 
    }
    // SCR needed for integrals  
    threadSCRSize += nERI*buffN4 + cacheSize + nERIbuff*NB4; 
    // SCR needed for AX, scale 2 for complex 
    threadSCRSize += size_t(sizeof(MatsT)/sizeof(double))*mRawMat*nBasis*nBasis; 

    // SCR needed for Schwarz Screening

#ifdef _SHZ_SCREEN_4C
    nSHZ = 3; // 1 for SchwarZERI, 1 for ShBlkNorm, 1 for ShBlkNorm_raw 
    if (contType == TWOBODY_CONTRACTION_TYPE::LLLL or 
        contType == TWOBODY_CONTRACTION_TYPE::LLSS) nSHZ += 1;
    
    size_t SHZThreadSCRSize  = std::max(nERI*buffN4 + cacheSize, size_t(nShells*nShells));  
    size_t SHZGeneralSCRSize = nShells*nShells*nSHZ; 
    
    threadSCRSize   = std::max(threadSCRSize, SHZThreadSCRSize);
    generalSCRSize += SHZGeneralSCRSize; 
#endif
    
    // scale to dcomplex
    // size_t scale = size_t(sizeof(dcomplex)/sizeof(double));

    return (threadSCRSize * GetNumThreads() + generalSCRSize) / 2;
  }

  template <>
  size_t GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintSCRSize(
    const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);  
    abort();
  }

  template <>
  size_t GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintSCRSize(
    const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {
    CErr("Complex  is an invalid option",std::cout);  
    abort();
  }

}; // namespace ChronusQ
