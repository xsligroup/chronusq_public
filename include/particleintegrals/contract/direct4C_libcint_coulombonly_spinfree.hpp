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
  void GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintCoulombOnlySpinFree(
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
  
      nERI = 1;
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
  
          //if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          if(int2e_pp1pp2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

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
      
      nERI = 1;
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
  
          //if(int2e_ip1ip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          if(int2e_gaunt_ps1ps2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

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
    
      nERI = 4;
      buffAll = CQMemManager::get().malloc<double>(2*nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      
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
          // skiperi1 = int2e_gauge_r1_ssp1ssp2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          // skiperi2 = int2e_gauge_r2_ssp1ssp2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi1 = int2e_gauge_r1_sp1sp2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi2 = int2e_gauge_r2_sp1sp2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);

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
        }
        
        SetLAThreads(LAThreads);// Turn threads for LA back on
        CQMemManager::get().free(ShBlkNormsSCR);
        
        // get the maximum from all batches 
        #pragma omp parallel for 
        for (auto i = 0ul; i < nShell*nShell; i++) 
        for (auto iBatch = 0ul; iBatch < nBatch; iBatch++) 
          ShBlkNorms[i] = std::max(ShBlkNorms_batch[iBatch][i], ShBlkNorms[i]); 
        
      } 
    /******************************************************/
    /*                                                    */
    /*End of Compute shell block norms (∞-norm) of all matList.X */
    /*                                                    */
    /******************************************************/
#endif
    
    /******************************************/
    /*                                        */
    /* Start of Dirac-Coulomb C(2)            */
    /* includes DC-LLLL and DC-LLSs/SSLL      */      
    /*                                        */
    /******************************************/
    
    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS ) {

 
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
      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(2*NB4*nThreads);

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
  
        double *ERIBuffAB   = &ERIBuffer[thread_id*NB4];
        double *ERIBuffCD   = &ERIBuffer[nThreads*NB4 + thread_id*NB4];
  
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
  
          // if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          if(int2e_pp1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

 
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
            auto dAdotdB = buff[mnkl];
            
            //auto MNKL = m + n*NB + k*NB2 + l*NB3;
            MNKL = m + nNBkNB2lNB3;
            //auto KLMN = k + l*NB + m*NB2 + n*NB3;
            KLMN = m*NB2 + klNBnNB3;
  
            // ∇A∙∇B(mn|kl) followed by ∇Ax∇B(mn|kl) X, Y, and Z
            // (mn|kl)
            ERIBuffAB[       MNKL] =  (double)s34_deg*dAdotdB;
  
            // ∇C∙∇D(kl|nm) followed by ∇Cx∇D(kl|nm) X, Y, and Z
            // (kl|mn)
            ERIBuffCD[       KLMN] =  (double)s12_deg*dAdotdB;
 
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

            symmDCLLMS  = matList[sDLLMS + matOff].X->S().pointer();
            symmDCSSMS  = matList[sDSS + matOff].X->S().pointer();

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
            DotPrdKLMN = KLMN;
  
            /*++++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb (LL|LL) Contraction */
            /*++++++++++++++++++++++++++++++++++++++++++++*/
  
            //KLMN
            if(bf3 >= bf4 ) {
              ADCLLMS[bf34] +=  ERIBuffCD[DotPrdKLMN] * symmDCSSMS[bf21];
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
      }}
    
      // Take care of the symmetry in the LL and SS blocks
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        auto ADCLLMS = matList[CLLMS + matOff].AX->S().pointer();
        auto ADCSSMS = matList[CSS   + matOff].AX->S().pointer();
        for( auto i = 0; i < nBasis; i++ ) 
        for( auto j = 0; j < i; j++ ) {
          ADCLLMS[j + i*nBasis] =  ADCLLMS[i + j*nBasis];
          ADCSSMS[j + i*nBasis] =  ADCSSMS[i + j*nBasis];
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

      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(NB4*nThreads);

      // Keeping track of number of integrals skipped
      std::vector<size_t> nIntSkipSSSS(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN_4C
      std::vector<size_t> nConSkipSSSS(nThreads,0);
#endif
      
      #pragma omp parallel
      {
  
        dcomplex iscale = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *ERIBuffABCD = &ERIBuffer[thread_id*NB4];
  
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
  
          //if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          if(int2e_pp1pp2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
  
          for(l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
          for(k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k) 
          for(n =   maxShellSize            ; n <   maxShellSize + n2; ++n) 
          for(m = 0                         ; m <                  n1; ++m, ++mnkl) {
  
            // (∇A∙∇B)(∇C∙∇D)(mnkl)
            auto dAdotdBdCdotdD =  buff[mnkl];
  
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto KLMN = k + l*NB + m*NB2 + n*NB3;

            // (mn|kl)
            ERIBuffABCD[         MNKL] =  dAdotdBdCdotdD;
  
            // (kl|mn)
            ERIBuffABCD[         KLMN] =  dAdotdBdCdotdD;
  
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
            auto KLMNdAdotdBdCdotdD          = ERIBuffABCD[         KLMN];
 
            // COULOMB
            // MNKL
            /* Equation (C19-C22) in Xiaosong's note */
            if(bf1 >= bf2) {
              auto ScaleF = bf3_s==bf4_s ? 0.5 : 1.0;
              ADCSSMS[bf12] += ScaleF * (symmDCSSMS[bf43] * MNKLdAdotdBdCdotdD);
            }
  
            /* KLMN */
            if((bf1_s!=bf3_s or bf2_s!=bf4_s) and bf3 >= bf4){
              auto ScaleF = bf1_s==bf2_s ? 0.5 : 1.0;
              ADCSSMS[bf34] += ScaleF * (symmDCSSMS[bf21] * KLMNdAdotdBdCdotdD);
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
      }}
      
      // take care of symmetries
      for (auto iBatch = 0ul, matOff = 0ul; iBatch < nBatch; iBatch++, matOff+=mMat) {  
        auto ADCSSMS = matList[matOff].AX->S().pointer();
        for( auto i = 0; i < nBasis; i++ ) 
        for( auto j = 0; j < i; j++ ) {
          ADCSSMS[j + i*nBasis] =  ADCSSMS[i + j*nBasis];
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
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectGaunt = tick();
      std::vector<double> durDirectGauntInt(nThreads, 0.);
      std::vector<double> durDirectGauntCon(nThreads, 0.);
#endif

      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);

      int nSave = 1;
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

          //if(int2e_ip1ip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          if(int2e_gaunt_ps1ps2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
 
          auto nQuad = n1*n2*n3*n4;

          for(l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l) 
          for(k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k) 
          for(m = 0                         ; m <                  n1; ++m)
          for(n =   maxShellSize            ; n <   maxShellSize + n2; ++n, ++mnkl) {

  
            /* Gaunt */
            // ∇A∙∇C(mn|kl)
            auto dAdotdC = -buff[mnkl];
            //auto dAdotdC = buff[AxCx*nQuad+mnkl] + buff[AyCy*nQuad+mnkl] + buff[AzCz*nQuad+mnkl];

            // Change the index so that we do ∇B∙∇C(ij|kl) using the ∇B∇C engine
            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;
  
            // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
            // (kl|mn)
            BC[      MNKL] = dAdotdC;    
            // (lk|nm)
            BC[      LKNM] = dAdotdC;
  
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
            ADCLSMS[bf12]+= - DCLSpmSLMS[bf43] * BC[MNKL];
            ADCLSMX[bf12]+= - DCLSpmSLMX[bf43] * BC[MNKL];
            ADCLSMY[bf12]+= - DCLSpmSLMY[bf43] * BC[MNKL];
            ADCLSMZ[bf12]+= - DCLSpmSLMZ[bf43] * BC[MNKL];
            // LKNM
            // Coulomb //TRANSKL
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              ADCLSMS[bf43]+= - DCLSpmSLMS[bf12] * BC[LKNM]; 
              ADCLSMX[bf43]+= - DCLSpmSLMX[bf12] * BC[LKNM];
              ADCLSMY[bf43]+= - DCLSpmSLMY[bf12] * BC[LKNM];
              ADCLSMZ[bf43]+= - DCLSpmSLMZ[bf12] * BC[LKNM];
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
  
      auto ss = 0     ;    // ERI00 (ss)(ij|kl)
      auto dot = NB4 ;         // ERI01 (σ∙σ)(ijkl)
   
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
      
      nERI = 4;
      buffAll = CQMemManager::get().malloc<double>(2*nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);

      int nSave = 2;
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
          //skiperi1 = int2e_gauge_r1_ssp1sps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          //skiperi2 = int2e_gauge_r2_ssp1sps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi1 = int2e_gauge_r1_sp1ps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi2 = int2e_gauge_r2_sp1ps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);

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
            BC[MNKL] = buffr1[3*nQuad+mnkl] - buffr2[3*nQuad+mnkl];
            BC[LKNM] = buffr1[3*nQuad+mnkl] - buffr2[3*nQuad+mnkl];
            
            // (σ∙σ)
            BC[NB4+MNKL] = -(buffr1[mnkl] - buffr2[mnkl])
                           -(buffr1[nQuad+mnkl] - buffr2[nQuad+mnkl])
                           -(buffr1[2*nQuad+mnkl] - buffr2[2*nQuad+mnkl]);
            BC[NB4+LKNM] = -(buffr1[mnkl] - buffr2[mnkl])
                           -(buffr1[nQuad+mnkl] - buffr2[nQuad+mnkl])
                           -(buffr1[2*nQuad+mnkl] - buffr2[2*nQuad+mnkl]);
            
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
            ADCLSMS[bf12]+= - DCLSpmSLMS[bf43] * BC[MNKL];
            
            // ADCLSMX[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+xs] * iS 
            //                - DCLSpmSLMX[bf43] * BC[MNKL+xx]
            //                - DCLSpmSLMY[bf43] * BC[MNKL+xy]
            //                - DCLSpmSLMZ[bf43] * BC[MNKL+xz]; 
            // 
            // ADCLSMY[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+ys] * iS 
            //                - DCLSpmSLMX[bf43] * BC[MNKL+yx]
            //                - DCLSpmSLMY[bf43] * BC[MNKL+yy]
            //                - DCLSpmSLMZ[bf43] * BC[MNKL+yz]; 
            // 
            // ADCLSMZ[bf12]+=- DCLSpmSLMS[bf43] * BC[MNKL+zs] * iS 
            //                - DCLSpmSLMX[bf43] * BC[MNKL+zx]
            //                - DCLSpmSLMY[bf43] * BC[MNKL+zy]
            //                - DCLSpmSLMZ[bf43] * BC[MNKL+zz]; 
            
            // LKNM
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              
              ADCLSMS[bf43]+= - DCLSpmSLMS[bf12] * BC[LKNM];
              
              // ADCLSMX[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+xs] * iS 
              //                - DCLSpmSLMX[bf12] * BC[LKNM+xx]
              //                - DCLSpmSLMY[bf12] * BC[LKNM+xy]
              //                - DCLSpmSLMZ[bf12] * BC[LKNM+xz]; 
              // 
              // ADCLSMY[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+ys] * iS 
              //                - DCLSpmSLMX[bf12] * BC[LKNM+yx]
              //                - DCLSpmSLMY[bf12] * BC[LKNM+yy]
              //                - DCLSpmSLMZ[bf12] * BC[LKNM+yz]; 
              // 
              // ADCLSMZ[bf43]+=- DCLSpmSLMS[bf12] * BC[LKNM+zs] * iS 
              //                - DCLSpmSLMX[bf12] * BC[LKNM+zx]
              //                - DCLSpmSLMY[bf12] * BC[LKNM+zy]
              //                - DCLSpmSLMZ[bf12] * BC[LKNM+zz]; 
              
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
        // matList[iMat].AX->X() +=  AXthreads[iTh][iMat]->X();  
        // matList[iMat].AX->Y() +=  AXthreads[iTh][iMat]->Y();  
        // matList[iMat].AX->Z() +=  AXthreads[iTh][iMat]->Z();  
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
    /*   End of Gauge  */
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
  void GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintCoulombOnlySpinFree(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyRelContraction<double>> &matList,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);  
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintCoulombOnlySpinFree(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyRelContraction<dcomplex>> &matList,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);  
  }

}; // namespace ChronusQ
