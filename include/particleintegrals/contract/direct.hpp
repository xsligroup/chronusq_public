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

#include <chrono>
#include <array>
#include <numeric>
#include <iomanip>

#include <integrals.hpp>
#include <particleintegrals/inhouseaointegral.hpp>
#include <util/matout.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>

#include <util/threads.hpp>
#include <util/mpi.hpp>
#include <util/timer.hpp>
#include <util/math.hpp>

#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/giaodirecteri.hpp>
#include <particleintegrals/gradints/direct.hpp>

#define _FULL_DIRECT
//#define _SUB_TIMINGS
//#define _REPORT_INTEGRAL_TIMINGS

//#define _PRECOMPUTE_SHELL_PAIRS

#define _SHZ_SCREEN
#define _SEPARATED_SHZ_SCREEN


#ifndef _FULL_DIRECT
  #define _BATCH_DIRECT
#endif

#if defined(_FULL_DIRECT) 
  #define _USE_EIGHT_FOLD
#else
  #define _USE_FOUR_FOLD
#endif

#ifndef _FULL_DIRECT
  #warning "Batch Direct ERI contraction is broken for complex"
#endif

#define GetRealPtr(X,I,J,N) reinterpret_cast<double*>(X + I + J*N)

#define bottomupGIAO //SS

//#define _PROFILE_DIRECT_GRAD

namespace ChronusQ {

  inline double ReProd(const double &a, const double &b) { return a*b; }
  inline double ReProd(const dcomplex &a, const dcomplex &b) {
    return a.real()*b.real() - a.imag()*b.imag();
  }

  template <typename MatsT>
  void ShellBlockNorm(std::vector<libint2::Shell> &shSet, MatsT *MAT, 
    size_t LDM, double *ShBlk) {

    size_t nShell = shSet.size();

    size_t n1,n2;
    for(auto s1(0ul), bf1(0ul); s1 < nShell; s1++, bf1 += n1) {
      n1 = shSet[s1].size();
    for(auto s2(0ul), bf2(0ul); s2 < nShell; s2++, bf2 += n2) {
      n2 = shSet[s2].size();

      MatsT *block = MAT + bf1 + bf2*LDM;
      ShBlk[s1 + s2*nShell] = lapack::lange(lapack::Norm::Inf,n1,n2,block,LDM);

    }
    }

  };


  template <typename T>
  double * ShellBlockNorm(std::vector<libint2::Shell> &shSet, T *MAT, 
    size_t LDM) {

    size_t nShell = shSet.size();
    double *ShBlk = CQMemManager::get().malloc<double>(nShell*nShell);

    ShellBlockNorm(shSet,MAT,LDM,ShBlk);

    return ShBlk;

  };


  template <typename MatsT, typename IntsT>
  void GTODirectTPIContraction<MatsT,IntsT>::directScaffold(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<MatsT>> &list,EMPerturbation&) const {

    DirectTPI<IntsT> &eri = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& basisSet_ = eri.basisSet();

    size_t nthreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(comm);
    size_t mpiSize   = MPISize(comm);

    SetLAThreads(1); // Turn off parallelism in LA functions

    const size_t NB   = basisSet_.nBasis;
    const size_t NMat = list.size();
    const size_t NS   = basisSet_.nShell;


#ifdef _SHZ_SCREEN
    // Check whether any of the contractions are non-hermetian
    const bool AnyNonHer = std::any_of(list.begin(),list.end(),
      []( TwoBodyContraction<MatsT> & x ) -> bool { return not x.HER; });

    // Compute schwarz bounds if we haven't already
    if(eri.schwarz() == nullptr) eri.computeSchwarz();
#endif


/*
    if( mpiSize > 1 )
      for(auto &C : list )
        prettyPrintSmart(std::cerr,"X in Direct",C.X,NB,NB,NB);
*/


    // Create thread-safe libint2::Engine's
      
    std::vector<libint2::Engine> engines(nthreads);

    // Construct engine for master thread
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      basisSet_.maxPrim, basisSet_.maxL, 0);





    // Allocate scratch for raw integral batches
    size_t maxShellSize = 
      std::max_element(basisSet_.shells.begin(),basisSet_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    size_t lenIntBuffer = 
      maxShellSize * maxShellSize * maxShellSize * maxShellSize; 

    lenIntBuffer *= sizeof(MatsT) / sizeof(double);

    size_t nBuffer = 2;

    size_t nAlloc = nBuffer*lenIntBuffer*nthreads*sizeof(double) + 
      nthreads*NMat*NB*NB*sizeof(MatsT) +
      list.size()*NS*NS*sizeof(double);
//  std::cerr << "DIRECT CONTRACTION " << nAlloc / 1e9 << std::endl;


    double * intBuffer = 
      CQMemManager::get().malloc<double>(nBuffer*lenIntBuffer*nthreads);
   
    double *intBuffer2 = intBuffer + nthreads*lenIntBuffer;


    // Allocate thread local storage to store integral contractions
    // XXX: Don't allocate anything if serial
    std::vector<std::vector<MatsT*>> AXthreads;
    MatsT *AXRaw = nullptr;
    if(nthreads != 1) {
      AXRaw = CQMemManager::get().malloc<MatsT>(nthreads*NMat*NB*NB);
      memset(AXRaw,0,nthreads*NMat*NB*NB*sizeof(MatsT));
    }

    for(auto ithread = 0, iMat = 0; ithread < nthreads; ithread++) {
      AXthreads.emplace_back();
      for(auto jMat = 0; jMat < NMat; jMat++, iMat++) {
        if(nthreads == 1) {
          AXthreads.back().push_back(list[jMat].AX);
        } else {
          AXthreads.back().push_back(AXRaw + iMat*NB*NB);
        }
      }
    }


#ifdef _SHZ_SCREEN
    // Compute shell block norms
    double *ShBlkNorms_raw = 
      CQMemManager::get().malloc<double>(list.size()*NS*NS);

    std::vector<double*> ShBlkNorms;
    for(auto iMat = 0, iOff = 0; iMat < NMat; iMat++, 
      iOff += NS*NS ) {

      ShellBlockNorm(basisSet_.shells,list[iMat].X,NB,
        ShBlkNorms_raw + iOff);

      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

    }

    double maxShBlk = 0.;
    for(auto iMat = 0; iMat < NMat; iMat++)
      maxShBlk = std::max(maxShBlk,
        *std::max_element(ShBlkNorms[iMat],ShBlkNorms[iMat] + NS*NS) ); 


    size_t NP4 = 
      basisSet_.maxPrim * basisSet_.maxPrim * basisSet_.maxPrim * 
      basisSet_.maxPrim;

    engines[0].set_precision(
      std::min(
        std::numeric_limits<double>::epsilon(),
        eri.threshSchwarz() / maxShBlk
      ) / NP4
    );



    // Get the max over all the matricies for
    // the shell block norms
    // OVERWRITES ShBlkNorms[0]
    for(auto k = 0; k < NS*NS; k++) {

      double mx = std::abs(ShBlkNorms[0][k]);
      for(auto iMat = 1; iMat < NMat; iMat++)
        mx = std::max(mx,std::abs(ShBlkNorms[iMat][k]));
      ShBlkNorms[0][k] = mx;

    }


    if( AnyNonHer )
    for(auto i = 0; i < NS; i++)
    for(auto j = 0; j <= i; j++) {
      double mx = 
        std::max(std::abs(ShBlkNorms[0][i + j*NS]),
                 std::abs(ShBlkNorms[0][j + i*NS]));

      for(auto iMat = 1; iMat < NMat; iMat++)
        mx = std::max(mx,
          std::max(std::abs(ShBlkNorms[iMat][i + j*NS]),
                   std::abs(ShBlkNorms[iMat][j + i*NS])));

      ShBlkNorms[0][i + j*NS] = mx;
      ShBlkNorms[0][j + i*NS] = mx;

    }


#else
    // Set precision
    engines[0].set_precision(std::numeric_limits<double>::epsilon());
#endif


    // Copy master thread engine to other threads
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

#ifdef _SUB_TIMINGS
    std::chrono::duration<double> durInner(0.), durCont(0.), durSymm(0.),
      durZero(0.);
#endif

    // Keeping track of number of integrals skipped
    std::vector<size_t> nSkip(nthreads,0);


    // MPI info
    size_t mpiChunks = (NS * (NS + 1) / 2) / mpiSize;
    size_t mpiS12St  = mpiRank * mpiChunks;
    size_t mpiS12End = (mpiRank + 1) * mpiChunks;
    if( mpiRank == (mpiSize - 1) ) mpiS12End = (NS * (NS + 1) / 2);

/*
    double t1 = MPI_Wtime();
    auto topDirect = std::chrono::high_resolution_clock::now();
*/
    auto topDirect = tick();
    #pragma omp parallel
    {

    // Set up thread local storage

    // SMP info
    size_t thread_id = GetThreadID();

    auto &engine = engines[thread_id];
    const auto& buf_vec = engine.results();
    
    auto &AX_loc = AXthreads[thread_id];


    double * intBuffer_loc  = intBuffer  + thread_id*lenIntBuffer;
    double * intBuffer2_loc = intBuffer2 + thread_id*lenIntBuffer;


    size_t n1,n2;

    // Always Loop over s2 <= s1
    for(size_t s1(0ul), bf1_s(0ul), s12(0ul); s1 < NS; bf1_s+=n1, s1++) { 
      n1 = basisSet_.shells[s1].size(); // Size of Shell 1

    auto sigPair12_it = basisSet_.shellData.shData.at(s1).begin();
    for( const size_t& s2 : basisSet_.shellData.sigShellPair[s1] ) {
      size_t bf2_s = basisSet_.mapSh2Bf[s2];
      n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      const auto * sigPair12 = sigPair12_it->get();
      sigPair12_it++;

#ifdef CQ_ENABLE_MPI
      // MPI partition s12 blocks
      if( (s12 < mpiS12St) or (s12 >= mpiS12End) ) { s12++; continue; }
#endif

      // Round-Robbin work distribution
      if( (s12++) % nthreads != thread_id ) continue;


      // Cache variables for shells 1 and 2
        
#ifdef _FULL_DIRECT
      // Deneneracy factor for s1,s2 pair
      double s12_deg = (s1 == s2) ? 1.0 : 2.0;
#endif

#ifdef _SHZ_SCREEN
      double shz12 = 0, shMax12 = 0;
      if( screen ) {
        shz12 = eri.schwarz()[s1 + s2*NS];
        shMax12 = ShBlkNorms[0][s1 + s2*NS];
      }
#endif



#ifdef _BATCH_DIRECT

#ifdef _SUB_TIMINGS
      auto topZero = std::chrono::high_resolution_clock::now();
#endif

      // Zero out the integral buffer (hot spot)
      memset(intBuffer_loc,0,lenIntBuffer);

#ifdef _SUB_TIMINGS
      auto botZero = std::chrono::high_resolution_clock::now();
      durZero += botZero - topZero;
#endif

      double *intBuffCur = intBuffer_loc;

#endif


#ifdef _SUB_TIMINGS
      auto topInner = std::chrono::high_resolution_clock::now();
#endif


// The upper bound of s3 is s1 for the 8-fold symmetry and
// nShell for 4-fold.
#ifdef _USE_EIGHT_FOLD
  #define S3_MAX s1
#elif defined(_USE_FOUR_FOLD)
  // the "-" is for the <= in the loop
  #define S3_MAX NS - 1
#endif

      size_t n3,n4;

      for(size_t s3(0ul), bf3_s(0ul), s34(0ul); s3 <= S3_MAX; s3++, bf3_s += n3) { 
        n3 = basisSet_.shells[s3].size(); // Size of Shell 3

#ifdef _SHZ_SCREEN

        double shMax123 = 0;
        if( screen ) {
          // Pre-calculate shell-block norm max's that only
          // depend on shells 1,2 and 3
          shMax123 = 
            std::max(ShBlkNorms[0][s1 + s3*NS], 
                     ShBlkNorms[0][s2 + s3*NS]);

          shMax123 = std::max(shMax123,shMax12);
        }

#endif
        
// The upper bound of s4 is either s2 or s3 based on s1 and s3 for
// the 8-fold symmetry and s3 for the 4-fold symmetry
#ifdef _USE_EIGHT_FOLD
        size_t s4_max = (s1 == s3) ? s2 : s3;
#elif defined(_USE_FOUR_FOLD)
        size_t s4_max =  s3;
#endif

      auto sigPair34_it = basisSet_.shellData.shData.at(s3).begin();
      for( const size_t& s4 : basisSet_.shellData.sigShellPair[s3] ) {

        if (s4 > s4_max)
          break;  // for each s3, s4 are stored in monotonically increasing
                  // order

        const auto * sigPair34 = sigPair34_it->get();
        sigPair34_it++;
                    
        size_t bf4_s = basisSet_.mapSh2Bf[s4];
        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

#ifdef _SHZ_SCREEN

        double shMax = 0;

        if( screen ) {
          // Compute Shell norm max
          shMax = 
            std::max(ShBlkNorms[0][s1 + s4*NS],
            std::max(ShBlkNorms[0][s2 + s4*NS],
                     ShBlkNorms[0][s3 + s4*NS]));

          shMax = std::max(shMax,shMax123);

          if((shMax * shz12 * eri.schwarz()[s3 + s4*NS]) <
             eri.threshSchwarz()) { nSkip[thread_id]++; continue; }
        }
#endif
      

#ifdef _FULL_DIRECT

        // Degeneracy factor for s3,s4 pair
        double s34_deg = (s3 == s4) ? 1.0 : 2.0;

        // Degeneracy factor for s1, s2, s3, s4 quartet
        double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;

        // Total degeneracy factor
        double s1234_deg = s12_deg * s34_deg * s12_34_deg;

#endif

        // Evaluate ERI for shell quartet (s1 s2 | s3 s4)
        engine.compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
          basisSet_.shells[s1],
          basisSet_.shells[s2],
          basisSet_.shells[s3],
          basisSet_.shells[s4]
#ifdef _PRECOMPUTE_SHELL_PAIRS
          ,sigPair12,sigPair34
#endif
        );

        // Libint2 internal screening
        const double *buff = buf_vec[0];

        if(buff == nullptr) { nSkip[thread_id]++; continue; }

#ifdef _BATCH_DIRECT

        // Copy over buffer
        //std::copy_n(buff,n1*n2*n3*n4,intBuffCur);
        memcpy(intBuffCur,buff,n1*n2*n3*n4*sizeof(double));
        intBuffCur += n1*n2*n3*n4;

#elif defined(_FULL_DIRECT)

// Flag to turn contraction on and off
#if 1
        // Scale the buffer by the degeneracy factor and store
        // in infBuffer
        std::transform(buff,buff + n1*n2*n3*n4,intBuffer_loc,
          [&](auto& x){ return x*0.5*s1234_deg; });

        size_t b1,b2,b3,b4;
        double *Xp1, *Xp2;
        double X1,X2;
        MatsT      T1,T2,T3,T4;
        MatsT      *Tp1,*Tp2;

        for(auto iMat = 0; iMat < NMat; iMat++) {
          
          // Hermetian contraction
          if( list[iMat].HER ) { 
            if( list[iMat].contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // Cache i,j variables
              b1 = bf1 + NB*bf2; 
              X1 = *reinterpret_cast<double*>(list[iMat].X  + b1);
              Xp1 = reinterpret_cast<double*>(AX_loc[iMat] + b1);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // J(1,2) += I * X(4,3)
              *Xp1 += *GetRealPtr(list[iMat].X,bf4,bf3,NB) * intBuffer_loc[ijkl];

              // J(4,3) += I * X(1,2)
              *GetRealPtr(AX_loc[iMat],bf4,bf3,NB) +=  X1 * intBuffer_loc[ijkl];

              // J(2,1) and J(3,4) are handled on symmetrization after
              // contraction
            } // kl loop
            } // ij loop

            else if( list[iMat].contType == EXCHANGE )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Cache i,j,k variables
              b1 = bf1 + bf3*NB;
              b2 = bf2 + bf3*NB;

              T1 = 0.5 * SmartConj(list[iMat].X[b1]);
              T2 = 0.5 * SmartConj(list[iMat].X[b2]);

            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // Indicies are swapped here to loop over contiguous memory
                
              // K(1,3) += 0.5 * I * X(2,4) = 0.5 * I * CONJ(X(4,2)) (**HER**)
              AX_loc[iMat][b1]           += 0.5 * SmartConj(list[iMat].X[bf4+NB*bf2]) * intBuffer_loc[ijkl];

              // K(4,2) += 0.5 * I * X(3,1) = 0.5 * I * CONJ(X(1,3)) (**HER**)
              AX_loc[iMat][bf4 + bf2*NB] += T1 * intBuffer_loc[ijkl];

              // K(4,1) += 0.5 * I * X(3,2) = 0.5 * I * CONJ(X(2,3)) (**HER**)
              AX_loc[iMat][bf4 + bf1*NB] += T2 * intBuffer_loc[ijkl];

              // K(2,3) += 0.5 * I * X(1,4) = 0.5 * I * CONJ(X(4,1)) (**HER**)
              AX_loc[iMat][b2]           += 0.5 * SmartConj(list[iMat].X[bf4+NB*bf1]) * intBuffer_loc[ijkl];

            } // l loop
            } // ijk

          // Nonhermetian contraction
          } else {

            if( list[iMat].contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // Cache i,j variables
              b1 = bf1 + NB*bf2; 
              T1 = *(list[iMat].X  + b1);
              Tp1 = (AX_loc[iMat] + b1);

              b2 = bf2 + NB*bf1; 
              T2 = *(list[iMat].X  + b2);
              Tp2 = (AX_loc[iMat] + b2);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // J(1,2) += I * X(4,3)
              *Tp1 += 0.5*( list[iMat].X[bf4 + bf3*NB] + list[iMat].X[bf3 + bf4*NB]) * intBuffer_loc[ijkl];

              // J(3,4) += I * X(2,1)
              AX_loc[iMat][bf3 + bf4*NB] +=  0.5*(T2+T1) * intBuffer_loc[ijkl];

              // J(2,1) += I * X(3,4)
              *Tp2 += 0.5*( list[iMat].X[bf4 + bf3*NB] + list[iMat].X[bf3 + bf4*NB]) * intBuffer_loc[ijkl];

              // J(4,3) += I * X(1,2)
              AX_loc[iMat][bf4 + bf3*NB] +=  0.5*(T2+T1) * intBuffer_loc[ijkl];

            } // kl loop
            } // ij loop

            else if( list[iMat].contType == EXCHANGE )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Cache i,j,k variables
              b1 = bf1 + bf3*NB;
              b2 = bf2 + bf3*NB;

              T1 = 0.5 * list[iMat].X[b1];
              T2 = 0.5 * list[iMat].X[b2];

              b3 = bf3 + bf1*NB;
              b4 = bf3 + bf2*NB;

              T3 = 0.5 * list[iMat].X[b3];
              T4 = 0.5 * list[iMat].X[b4];
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // K(3,1) += 0.5 * I * X(4,2)
              AX_loc[iMat][b3]           += 0.5 * list[iMat].X[bf4+NB*bf2] * intBuffer_loc[ijkl];

              // K(4,2) += 0.5 * I * X(3,1)
              AX_loc[iMat][bf4 + bf2*NB] += T3 * intBuffer_loc[ijkl];
 
              // K(4,1) += 0.5 * I * X(3,2)
              AX_loc[iMat][bf4 + bf1*NB] += T4 * intBuffer_loc[ijkl];

              // K(3,2) += 0.5 * I * X(4,1)
              AX_loc[iMat][b4]           += 0.5 * list[iMat].X[bf4+NB*bf1] * intBuffer_loc[ijkl];

              // K(1,3) += 0.5 * I * X(2,4)
              AX_loc[iMat][b1]           += 0.5 * list[iMat].X[bf2+NB*bf4] * intBuffer_loc[ijkl];

              // K(2,4) += 0.5 * I * X(1,3)
              AX_loc[iMat][bf2 + bf4*NB] += T1 * intBuffer_loc[ijkl];
 
              // K(1,4) += 0.5 * I * X(2,3)
              AX_loc[iMat][bf1 + bf4*NB] += T2 * intBuffer_loc[ijkl];

              // K(2,3) += 0.5 * I * X(1,4)
              AX_loc[iMat][b2]           += 0.5 * list[iMat].X[bf1+NB*bf4] * intBuffer_loc[ijkl];

            } // l loop
            } // ijk

          } // Symmetry check

        } // iMat loop

#endif

#endif

      } // loop s4
      } // loop s3

#ifdef _SUB_TIMINGS
      auto botInner = std::chrono::high_resolution_clock::now();

      durInner += botInner - topInner;
#endif

#ifdef _BATCH_DIRECT
#if 0
      assert(nthreads == 1);

#ifdef _SUB_TIMINGS
      auto topSymm = std::chrono::high_resolution_clock::now();
#endif

      // Reorder and expand integrals into square matricies
      for(auto s3 = 0ul, bf3_s = 0ul, ijkl = 0ul; s3 < NS; s3++, 
        bf3_s += n3) { 
        n3 = basisSet_.shells[s3].size();

      for(auto s4 = 0ul, bf4_s = 0ul; s4 <= s3; s4++, bf4_s += n4) { 
        n4 = basisSet_.shells[s4].size();

        for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)      
        for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
        for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
        for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

          intBuffer2_loc[bf4 + bf3*NB + j*nSQ_ + i*nSQ_*n2] = intBuffer_loc[ijkl];
          intBuffer2_loc[bf3 + bf4*NB + j*nSQ_ + i*nSQ_*n2] = intBuffer_loc[ijkl];

        }

      }
      }

#ifdef _SUB_TIMINGS
      auto botSymm = std::chrono::high_resolution_clock::now();
      durSymm += botSymm - topSymm;
#endif
      
     
#ifdef _SUB_TIMINGS
      auto topCont = std::chrono::high_resolution_clock::now();
#endif



      // Perform batched contractions
      for(auto &C : list ) {

        // J Contraction
        if( C.contType == COULOMB ) {

          blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,n1*n2,1,nSQ_,T(1.),intBuffer2_loc,nSQ_,C.X,nSQ_,
            T(0.),reinterpret_cast<G*>(intBuffer_loc),n1*n2);

          // Populate the lower triangle of J contraction storage
          for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)      
          for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)
            C.AX[bf1 + bf2*NB] = intBuffer_loc[j + i*n2]; 

        // K Contraction
        } else if( C.contType == EXCHANGE ) {

          for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)      
          for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) {

/*
            // T(m,n) = I(m,k) * X(n,k)
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NB,NB,NB,T(1.),intBuffer2 +j*nSQ_ + i*n2*nSQ_,NB,
              C.X,NB,T(0.),reinterpret_cast<G*>(intBuffer),NB);

            for(auto nu = 0; nu < NB; nu++) {
              C.AX[nu + bf1*NB] += intBuffer[nu + bf2*NB];
              if(s1 != s2) { 
                C.AX[nu + bf2*NB] += intBuffer[nu + bf1*NB];
              }
            }
*/

            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NB,1,NB,T(1.),intBuffer2_loc +j*nSQ_ + i*n2*nSQ_,NB,
              C.X + bf2,NB,T(0.),reinterpret_cast<G*>(intBuffer_loc),NB);
            for(auto nu = 0; nu < NB; nu++) 
              C.AX[nu + bf1*NB] += intBuffer_loc[nu];

            if( s1 != s2 ) {

              blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NB,1,NB,T(1.),intBuffer2_loc +j*nSQ_ + i*n2*nSQ_,NB,
                C.X + bf1,NB,T(0.),reinterpret_cast<G*>(intBuffer_loc),NB);
              for(auto nu = 0; nu < NB; nu++) 
                C.AX[nu + bf2*NB] += intBuffer_loc[nu];

            }

          } // ij loop

        } // Exchange check
      } // Loop over contractions

#ifdef _SUB_TIMINGS
      auto botCont = std::chrono::high_resolution_clock::now();
      durCont += botCont - topCont;
#endif

#endif
#endif

    }; // s2
    }; // s1


    }; // OpenMP context

/*
    auto botDirect = std::chrono::high_resolution_clock::now();

    double t2 = MPI_Wtime();
*/

    auto durDirect = tock(topDirect);

#ifdef _REPORT_INTEGRAL_TIMINGS
    size_t nIntSkip = std::accumulate(nSkip.begin(),nSkip.end(),0);
    std::cerr << "Screened " << nIntSkip << std::endl;

//  std::chrono::duration<double> durDirect = botDirect - topDirect;
    //std::cerr << "Direct Contraction took " << durDirect.count() << " s\n"; 
    std::cerr << "Direct Contraction took " <<  durDirect << " s\n"; 

#ifdef _SUB_TIMINGS
    std::cerr << "  " << durInner.count() << " (" << durInner.count() / durDirect.count() * 100 
              << "%) Inner loop" << std::endl;
#ifndef _FULL_DIRECT
    std::cerr << "  " << durZero.count() << " (" << durZero.count() / durDirect.count() * 100 
              << "%) Zeroing buffer" << std::endl;
    std::cerr << "  " << durSymm.count() << " (" << durSymm.count() / durDirect.count() * 100 
              << "%) Symmetrization loop" << std::endl;
    std::cerr << "  " << durCont.count() << " (" << durCont.count() / durDirect.count() * 100 
              << "%) Contraction loop" << std::endl;
#endif
#endif
    std::cerr << std::endl;
#endif


#ifdef _FULL_DIRECT

    MatsT* SCR = CQMemManager::get().malloc<MatsT>(NB*NB);
    for( auto iMat = 0; iMat < NMat;  iMat++ ) 
    for( auto iTh  = 0; iTh < nthreads; iTh++) {
  
    //prettyPrintSmart(std::cerr,"AX " + std::to_string(iMat) + " " + std::to_string(iTh),
    //  AXthreads[iTh][iMat],NB,NB,NB);

      if( list[iMat].HER ) {

        MatAdd('N','C',NB,NB,MatsT(0.5),AXthreads[iTh][iMat],NB,MatsT(0.5),
          AXthreads[iTh][iMat],NB,SCR,NB);

        if( nthreads != 1 )
          MatAdd('N','N',NB,NB,MatsT(1.),SCR,NB,MatsT(1.), list[iMat].AX,NB,list[iMat].AX,NB);
        else
          SetMat('N',NB,NB,MatsT(1.),SCR,NB,list[iMat].AX,NB);

      } else {

        if( nthreads != 1 )
          MatAdd('N','N',NB,NB,MatsT(0.5),AXthreads[iTh][iMat],NB,
            MatsT(1.), list[iMat].AX,NB,list[iMat].AX,NB);
        else 
          blas::scal(NB*NB,MatsT(0.5),list[iMat].AX,1);


      //std::transform(AXthreads[iTh][iMat], AXthreads[iTh][iMat] + NB*NB, 
      //  list[iMat].AX, []( G x ) -> G { return x / 4.; } );
      }

    };
    CQMemManager::get().free(SCR);
    
#else

    for( auto &C : list ) {

      // Symmetrize J contraction
      if( C.contType == COULOMB ) 
        HerMat('L',NB,C.AX,NB);
  
      // Inplace transpose of K contraction
      if( C.contType == EXCHANGE ) 
        IMatCopy('C',NB,NB,MatsT(1.),C.AX,NB,NB);

    } // Loop over contractions

#endif


#ifdef CQ_ENABLE_MPI
    // Combine all G[X] contributions onto Root process
    if( mpiSize > 1 ) {

      // FIXME: This should be able to be done with MPI_IN_PLACE for
      // the root process
        
      MatsT* mpiScr;
      if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<MatsT>(NB*NB);

      for( auto &C : list ) {
//      prettyPrintSmart(std::cerr,"AX in Direct",C.AX,NB,NB,NB);

        MPIReduce( C.AX, NB*NB, mpiScr, 0, comm );

        // Copy over the output buffer on root
        if( mpiRank == 0 ) std::copy_n(mpiScr,NB*NB,C.AX);

      }

      if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

    }

#endif




#ifdef _SUB_TIMINGS
    auto topFree = std::chrono::high_resolution_clock::now();
#endif

    // Free scratch space
    CQMemManager::get().free(intBuffer);
#ifdef _SHZ_SCREEN
    CQMemManager::get().free(ShBlkNorms_raw);
#endif
    if(AXRaw != nullptr) CQMemManager::get().free(AXRaw);

#ifdef _SUB_TIMINGS
    auto botFree = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> durFree = botFree - topFree;

    std::cerr << "Free took " << durFree.count() << "s" << std::endl;
#endif


    // Turn threads for LA back on
    SetLAThreads(LAThreads);

  };



  template <>
  void GTODirectTPIContraction<dcomplex,dcomplex>::directScaffold(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<dcomplex>> &list, EMPerturbation &pert) const {

    DirectTPI<dcomplex> &eri = *std::dynamic_pointer_cast<DirectTPI<dcomplex>>(this->ints_);
    BasisSet& basisSet_ = eri.basisSet();
     

    size_t nthreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(comm);
    size_t mpiSize   = MPISize(comm);

    const size_t NB   = basisSet_.nBasis;
    const size_t NMat = list.size();
    const size_t NS   = basisSet_.nShell;

    auto magAmp = pert.getDipoleAmp(Magnetic);  

    // Allocate scratch for raw integral batches
    size_t maxShellSize = 
      std::max_element(basisSet_.shells.begin(),basisSet_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    size_t lenIntBuffer = 
      maxShellSize * maxShellSize * maxShellSize * maxShellSize; 

    size_t nBuffer = 2;

    dcomplex * intBuffer = 
      CQMemManager::get().malloc<dcomplex>(nBuffer*lenIntBuffer*nthreads);
   
    // double *intBuffer2 = intBuffer + nthreads*lenIntBuffer;

    dcomplex * alterintBuffer = 
      CQMemManager::get().malloc<dcomplex>(nBuffer*lenIntBuffer*nthreads);

    // Allocate thread local storage to store integral contractions
    // XXX: Don't allocate anything if serial
    std::vector<std::vector<dcomplex*>> AXthreads;
    dcomplex *AXRaw = nullptr;
    if(nthreads != 1) {
      AXRaw = CQMemManager::get().malloc<dcomplex>(nthreads*NMat*NB*NB);    
      memset(AXRaw,0,nthreads*NMat*NB*NB*sizeof(dcomplex));
    }

    for(auto ithread = 0, iMat = 0; ithread < nthreads; ithread++) {
      AXthreads.emplace_back();
      for(auto jMat = 0; jMat < NMat; jMat++, iMat++) {
        if(nthreads == 1) {
          AXthreads.back().push_back(list[jMat].AX);
        } else {
          AXthreads.back().push_back(AXRaw + iMat*NB*NB);
        }
      }
    }


    // MPI info
    size_t mpiChunks = (NS * (NS + 1) / 2) / mpiSize;
    size_t mpiS12St  = mpiRank * mpiChunks;
    size_t mpiS12End = (mpiRank + 1) * mpiChunks;
    if( mpiRank == (mpiSize - 1) ) mpiS12End = (NS * (NS + 1) / 2);






    // start parallel
    #pragma omp parallel
    {

    // Set up thread local storage

    // SMP info
    size_t thread_id = GetThreadID();

    auto &AX_loc = AXthreads[thread_id];


    dcomplex * intBuffer_loc  = intBuffer  + thread_id*lenIntBuffer;
    dcomplex * alterintBuffer_loc  = alterintBuffer  + thread_id*lenIntBuffer;


    size_t n1,n2;

    // Always Loop over s2 <= s1
    for(size_t s1(0ul), bf1_s(0ul), s12(0ul); s1 < NS; bf1_s+=n1, s1++) { 
      n1 = basisSet_.shells[s1].size(); // Size of Shell 1

    for ( int s2 = 0 ; s2 <= s1 ; s2++ ) {
      size_t bf2_s = basisSet_.mapSh2Bf[s2];
      n2 = basisSet_.shells[s2].size(); // Size of Shell 2

#ifdef CQ_ENABLE_MPI
      // MPI partition s12 blocks
      if( (s12 < mpiS12St) or (s12 >= mpiS12End) ) { s12++; continue; }
#endif

      // Round-Robbin work distribution
      if( (s12++) % nthreads != thread_id ) continue;

#ifdef _FULL_DIRECT
      // Deneneracy factor for s1,s2 pair
      double s12_deg = (s1 == s2) ? 1.0 : 2.0;
#endif


      //SS Start generate shellpair1 

      libint2::ShellPair pair1_to_use;
      pair1_to_use.init( basisSet_.shells[s1],basisSet_.shells[s2],-1000);

      libint2::ShellPair pair1_to_use_switch;
      pair1_to_use_switch.init( basisSet_.shells[s2],basisSet_.shells[s1],-1000); 

// The upper bound of s3 is s1 for the 8-fold symmetry and
// nShell for 4-fold.
#ifdef _USE_EIGHT_FOLD
  #define S3_MAX s1
#elif defined(_USE_FOUR_FOLD)
  // the "-" is for the <= in the loop
  #define S3_MAX NS - 1
#endif

      size_t n3,n4;

      for(size_t s3(0), bf3_s(0), s34(0); s3 <= S3_MAX; s3++, bf3_s += n3) {
        n3 = basisSet_.shells[s3].size(); // Size of Shell 3





// The upper bound of s4 is either s2 or s3 based on s1 and s3 for
// the 8-fold symmetry and s3 for the 4-fold symmetry
#ifdef _USE_EIGHT_FOLD
        size_t s4_max = (s1 == s3) ? s2 : s3;
#elif defined(_USE_FOUR_FOLD)
        size_t s4_max =  s3;
#endif

      for ( int s4 = 0 ; s4 <=s4_max ; s4++ ){     
        if (s4 > s4_max)
          break;  // for each s3, s4 are stored in monotonically increasing
                  // order

#if 0
        const auto * sigPair34 = sigPair34_it->get();
        sigPair34_it++;
#endif 
                    
        size_t bf4_s = basisSet_.mapSh2Bf[s4];
        n4 = basisSet_.shells[s4].size(); // Size of Shell 4


        //SS start generate shellpair2 and calculate GIAO ERI

        libint2::ShellPair pair2_to_use;
        
        pair2_to_use.init( basisSet_.shells[s3],basisSet_.shells[s4],-1000);

/*
        libint2::ShellPair pair2_to_use_switch;
        // switch s3 and s4
        pair2_to_use_switch.init( basisSet_.shells[s4],basisSet_.shells[s3],-1000); 
*/

#ifdef _FULL_DIRECT

        // Degeneracy factor for s3,s4 pair
        double s34_deg = (s3 == s4) ? 1.0 : 2.0;

        // Degeneracy factor for s1, s2, s3, s4 quartet
        double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;

        // Total degeneracy factor
        double s1234_deg = s12_deg * s34_deg * s12_34_deg;

#endif

        // Evaluate ERI for shell quartet (s1 s2 | s3 s4)  

// std::cout<<"s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;
#ifdef bottomupGIAO

        // calculate integral (s1,s2|s3,s4)
        auto two2buff = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use,pair2_to_use,
          basisSet_.shells[s1],basisSet_.shells[s2],
          basisSet_.shells[s3],basisSet_.shells[s4],&magAmp[0],-1.0,-1.0);


        // calculate integral (s1,s2|s4,s3)
        auto two2buff_switch = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use_switch,pair2_to_use,
          basisSet_.shells[s2],basisSet_.shells[s1],
          basisSet_.shells[s3],basisSet_.shells[s4],&magAmp[0],-1.0,-1.0);
#else 

        // calculate integral (s1,s2|s3,s4)
        auto two2buff = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use,pair2_to_use,
          basisSet_.shells[s1],basisSet_.shells[s2],
          basisSet_.shells[s3],basisSet_.shells[s4],&magAmp[0],-1.0,-1.0);


        // calculate integral (s1,s2|s4,s3)
        auto two2buff_switch = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use_switch,pair2_to_use,
          basisSet_.shells[s2],basisSet_.shells[s1],
          basisSet_.shells[s3],basisSet_.shells[s4],&magAmp[0],-1.0,-1.0);

#endif 

        const dcomplex *buff = &(two2buff[0]); 
        const dcomplex *buffswitch = &(two2buff_switch[0]); 

#ifdef _FULL_DIRECT

// Flag to turn contraction on and off
#if 1
        // Scale the buffer by the degeneracy factor and store
        // in infBuffer

        std::transform(buff,buff + n1*n2*n3*n4 , intBuffer_loc,
          [&](auto& x){ return x*0.5*s1234_deg; });

        std::transform(buffswitch,buffswitch+n1*n2*n3*n4 ,alterintBuffer_loc,
          [&](auto& x){ return x*0.5*s1234_deg; });

        size_t b1,b2,b3,b4;

        for(auto iMat = 0; iMat < NMat; iMat++) {
          auto& C = list[iMat];
          
          // Hermetian contraction
          if( C.HER ) { 

            if( C.contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // bf1, bf2 are the index in the matrix
              // b1 is the index in the whole trunk of number 
              // Cache i,j variables
              b1 = bf1 + NB*bf2; 
              
              // in GIAO, J is complex. So X1 and Xp1 are not required 
              
              // X1 = *reinterpret_cast<double*>(list[iMat].X  + b1);
              // Xp1 = reinterpret_cast<double*>(AX_loc[iMat] + b1);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) {

              int jikl;
              jikl = j*n1*n3*n4 + i*n3*n4 + k*n4 + l;  


              // J(1,2) += 1/2[I(1,2|3,4) * X(4,3) + I(1,2|4,3) * X(3,4)]
              AX_loc[iMat][b1] += 0.5 * (C.X[bf4+bf3*NB] * intBuffer_loc[ijkl]
                            + C.X[bf3+bf4*NB] * std::conj(alterintBuffer_loc[jikl]));     

              // J(4,3) += 1/2[I(4,3|2,1) * X(1,2)+I(4,3|1,2) * X(2,1)
              //         = 1/2[I(1,2|3,4)* *X(1,2)+I(1,2|4,3) * X(2,1)
              AX_loc[iMat][bf4+bf3*NB] +=  0.5 *( C.X[b1] 
                                * std::conj(intBuffer_loc[ijkl])
                            + C.X[bf2+bf1*NB] * std::conj(alterintBuffer_loc[jikl]));
                                                                             
              // J(2,1) and J(3,4) are handled on symmetrization after
              // contraction
                
            } // kl loop
            } // ij loop

            else if( C.contType == EXCHANGE )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Cache i,j,k variables
              b1 = bf1 + bf3*NB;
              b2 = bf2 + bf3*NB;

              dcomplex T1 = 0.5 * SmartConj(C.X[b1]);
              dcomplex T2 = 0.5 * SmartConj(C.X[b2]);

            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              int jikl;
              jikl = j*n1*n4*n3 + i*n4*n3 + k*n4 + l;  
              // Indicies are swapped here to loop over contiguous memory
                
              // K(1,3) += 0.5 * I(1,2|4,3) * X(2,4) = 0.5 * I * CONJ(X(4,2)) (**HER**)
              AX_loc[iMat][b1]           += 0.5 * SmartConj(C.X[bf4+NB*bf2]) * std::conj(alterintBuffer_loc[jikl]);

              // K(4,2) += 0.5 * I(4,3|1,2) * X(3,1) = 0.5 * I * CONJ(X(1,3)) (**HER**)
              AX_loc[iMat][bf4 + bf2*NB] += T1 * std::conj(alterintBuffer_loc[jikl]);

              // K(4,1) += 0.5 * I(4,3|2,1) * X(3,2) = 0.5 * I * CONJ(X(2,3)) (**HER**)
              AX_loc[iMat][bf4 + bf1*NB] += T2 * std::conj( intBuffer_loc[ijkl] );

              // K(2,3) += 0.5 * I(2,1|4,3) * X(1,4) = 0.5 * I * CONJ(X(4,1)) (**HER**)
              AX_loc[iMat][b2]           += 0.5 * SmartConj(C.X[bf4+NB*bf1]) * std::conj( intBuffer_loc[ijkl] );

            } // l loop
            } // ijk

          } else {    // here is non Hermitian contraction   


            if( C.contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // Cache i,j variables
              b1 = bf1 + NB*bf2; 
              // T1 = *(C.X  + b1);
              // Tp1 = (AX_loc[iMat] + b1);

              b2 = bf2 + NB*bf1; 
              // T2 = *(list[iMat].X  + b2);
              // Tp2 = (AX_loc[iMat] + b2);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 


              int jikl;
              jikl = j*n1*n3*n4 + i*n3*n4 + k*n4 + l;  


              // J(1,2) += 1/2*(I(1,2|3,4) * X(4,3) + I(1,2|4,3) * X(3,4) 
              AX_loc[iMat][b1] += 0.5*( C.X[bf4 + bf3*NB]*intBuffer_loc[ijkl] 
                         +C.X[bf3 + bf4*NB] * std::conj(alterintBuffer_loc[jikl]));

              // J(3,4) += 1/2(I(3,4|1,2) * X(2,1) + I(3,4|2,1) * X(1,2))
              AX_loc[iMat][bf3 + bf4*NB] += 0.5*( C.X[b2] * intBuffer_loc[ijkl]
                                    +C.X[b1] * alterintBuffer_loc[jikl]) ;

              // J(2,1) += 1/2(I(2,1|3,4) * X(4,3) + I(2,1|4,3)* X(3,4))
              AX_loc[iMat][b2] += 0.5*( C.X[bf4 + bf3*NB] * alterintBuffer_loc[jikl] 
                                + C.X[bf3 + bf4*NB] * std::conj(intBuffer_loc[ijkl]));

              // J(4,3) += 1/2[I(4,3|2,1) * X(1,2)+I(4,3|1,2) * X(2,1)
              //         = 1/2[I(1,2|3,4)* *X(1,2)+I(1,2|4,3) * X(2,1)
              AX_loc[iMat][bf4+bf3*NB] +=  0.5 * (C.X[b1] 
                                * std::conj(intBuffer_loc[ijkl])
                            + C.X[bf2+bf1*NB] * std::conj(alterintBuffer_loc[jikl]));

            } // kl loop
            } // ij loop

            else if( C.contType == EXCHANGE )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Cache i,j,k variables
              b1 = bf1 + bf3*NB;
              b2 = bf2 + bf3*NB;

              // T1 = 0.5 * list[iMat].X[b1];
              // T2 = 0.5 * list[iMat].X[b2];

              b3 = bf3 + bf1*NB;
              b4 = bf3 + bf2*NB;

              // T3 = 0.5 * list[iMat].X[b3];
              // T4 = 0.5 * list[iMat].X[b4];
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              int jikl;
              jikl = j*n1*n3*n4 + i*n3*n4 + k*n4 + l;  

              // K(3,1) += 0.5 * I(3,4|2,1) * X(4,2)
              AX_loc[iMat][b3]           += 0.5 * C.X[bf4+NB*bf2] * alterintBuffer_loc[jikl];

              // K(4,2) += 0.5 * I(4,3|1,2) * X(3,1)
              AX_loc[iMat][bf4 + bf2*NB] += 0.5 * C.X[b3] * std::conj(alterintBuffer_loc[jikl]);
 
              // K(4,1) += 0.5 * I(4,3|2,1) * X(3,2)
              AX_loc[iMat][bf4 + bf1*NB] += 0.5 * C.X[b4] * std::conj(intBuffer_loc[ijkl]);

              // K(3,2) += 0.5 * I(3,4|1,2) * X(4,1)
              AX_loc[iMat][b4]           += 0.5 * C.X[bf4+NB*bf1] * intBuffer_loc[ijkl];

              // K(1,3) += 0.5 * I(1,2|4,3) * X(2,4)
              AX_loc[iMat][b1]           += 0.5 * C.X[bf2+NB*bf4] * std::conj(alterintBuffer_loc[jikl]);

              // K(2,4) += 0.5 * I(2,1|3,4) * X(1,3)
              AX_loc[iMat][bf2 + bf4*NB] += 0.5 * C.X[b1] * alterintBuffer_loc[jikl];
 
              // K(1,4) += 0.5 * I(1,2|3,4) * X(2,3)
              AX_loc[iMat][bf1 + bf4*NB] += 0.5 * C.X[b2] * intBuffer_loc[ijkl];

              // K(2,3) += 0.5 * I(2,1|4,3) * X(1,4)
              AX_loc[iMat][b2]           += 0.5 * C.X[bf1+NB*bf4] * std::conj(intBuffer_loc[ijkl]);

            } // l loop
            } // ijk

          } // non Hermitian finished 


        } // iMat loop

#endif
// this is for if 1

#endif
// this is for ifdef defined(_FULL_DIRECT) 


      } // loop s4
      } // loop s3



    }; // s2
    }; // s1

    } // end parallel

#ifdef _FULL_DIRECT

    dcomplex* SCR = CQMemManager::get().malloc<dcomplex>(NB*NB);
    for( auto iMat = 0; iMat < NMat;  iMat++ ) 
    for( auto iTh  = 0; iTh < nthreads; iTh++) {
  
    //prettyPrintSmart(std::cerr,"AX " + std::to_string(iMat) + " " + std::to_string(iTh),
    //  AXthreads[iTh][iMat],NB,NB,NB);

      if( list[iMat].HER ) {

        MatAdd('N','C',NB,NB,dcomplex(0.5),AXthreads[iTh][iMat],NB,dcomplex(0.5),
          AXthreads[iTh][iMat],NB,SCR,NB);

        if( nthreads != 1 )
          MatAdd('N','N',NB,NB,dcomplex(1.),SCR,NB,dcomplex(1.), list[iMat].AX,NB,list[iMat].AX,NB);
        else
          SetMat('N',NB,NB,dcomplex(1.),SCR,NB,list[iMat].AX,NB);

      } else {

        if( nthreads != 1 )
          MatAdd('N','N',NB,NB,dcomplex(0.5),AXthreads[iTh][iMat],NB,
            dcomplex(1.), list[iMat].AX,NB,list[iMat].AX,NB);
        else 
          blas::scal(NB*NB,dcomplex(0.5),list[iMat].AX,1);


      }

    };
    CQMemManager::get().free(SCR);
    

#endif

#ifdef CQ_ENABLE_MPI
    // Combine all G[X] contributions onto Root process
    if( mpiSize > 1 ) {

      // FIXME: This should be able to be done with MPI_IN_PLACE for
      // the root process
        
      dcomplex* mpiScr;
      if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<dcomplex>(NB*NB);

      for( auto &C : list ) {
//      prettyPrintSmart(std::cerr,"AX in Direct",C.AX,NB,NB,NB);

        MPIReduce( C.AX, NB*NB, mpiScr, 0, comm );

        // Copy over the output buffer on root
        if( mpiRank == 0 ) std::copy_n(mpiScr,NB*NB,C.AX);

      }

      if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

    }

#endif


    // Free scratch space
    CQMemManager::get().free(intBuffer);
    CQMemManager::get().free(alterintBuffer);

    if(AXRaw != nullptr) CQMemManager::get().free(AXRaw);





  }

  template <>
  void GTODirectTPIContraction<double,double>::directScaffold(
    MPI_Comm c, const bool b, 
    std::vector<TwoBodyContraction<double>> &list, EMPerturbation &pert) const {
    CErr("GIAO + Real is an invalid option",std::cout);  
  }

  template <>
  void GTODirectTPIContraction<dcomplex,double>::directScaffold(
    MPI_Comm c, const bool b, 
    std::vector<TwoBodyContraction<dcomplex>> &list, EMPerturbation &pert) const {
    CErr("GIAO + Real is an invalid option",std::cout);  
  }


  // New Direct Code for 2-Particle contraction
  template <typename MatsT, typename IntsT>
  void GTODirectTPIContraction<MatsT,IntsT>::directScaffoldNew(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<MatsT>> &matList) const {

    size_t parentId(0);
    size_t callLevel(0);
//    parentId = ProgramTimer::tick("Contract Total");
//    callLevel = ProgramTimer::getCallLevel();

    DirectTPI<IntsT> &tpi = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& basisSet_  = this->contractSecond ? tpi.basisSet2() : tpi.basisSet();
    BasisSet& basisSet2_ = this->contractSecond ? tpi.basisSet()  : tpi.basisSet2();

    size_t nThreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(comm);
    size_t mpiSize   = MPISize(comm);

    SetLAThreads(1); // Turn off parallelism in LA functions

    const size_t nBasis   = basisSet_.nBasis;
    const size_t snBasis  = basisSet2_.nBasis;
    const size_t nMat     = matList.size();
    const size_t nShell   = basisSet_.nShell;
    const size_t snShell  = basisSet2_.nShell;
    const typename DirectTPI<IntsT>::Kernel eriKernel = tpi.kernel();


    // Check whether any of the contractions are non-hermetian
    const bool NonHermitian = std::any_of(matList.begin(),matList.end(),
      []( TwoBodyContraction<MatsT> & x ) -> bool { return not x.HER; });


    bool sameBasisSet12 = &basisSet_ == &basisSet2_;
#ifdef _SHZ_SCREEN
    // Compute schwarz bounds if we haven't already
    if(tpi.schwarz() == nullptr or tpi.schwarz2() == nullptr) 
      tpi.computeSchwarz();

    double * schwarz1 = this->contractSecond ? tpi.schwarz2() : tpi.schwarz();
    double * schwarz2 = this->contractSecond ? tpi.schwarz()  : tpi.schwarz2();

    if (sameBasisSet12) schwarz2 = schwarz1;
#endif


    // Create thread-safe libint2::Engines
    std::vector<libint2::Engine> engines(nThreads);

    // Construct engine for master thread
    engines[0] = libint2::Engine(tpi.libintOperator(),
      std::max(basisSet_.maxPrim, basisSet2_.maxPrim), 
      std::max(basisSet_.maxL, basisSet2_.maxL),0);

    if (eriKernel == DirectTPI<IntsT>::Kernel::ShortRangeErfc)
      engines[0].set_params(tpi.rangeSeparationParameter());


    // Allocate scratch for raw integral batches
    size_t maxShellSize = 
      std::max_element(basisSet_.shells.begin(),basisSet_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    size_t maxShellSize2 = 
      std::max_element(basisSet2_.shells.begin(),basisSet2_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    // lenIntBuffer is allocated to be able to store EPAI's of the 
    // shell with the highest angular momentum
    size_t lenIntBuffer = 
      maxShellSize * maxShellSize * maxShellSize2 * maxShellSize2; 

    lenIntBuffer *= sizeof(MatsT) / sizeof(double);

    size_t nBuffer = 2;

    double * intBuffer = 
      CQMemManager::get().malloc<double>(nBuffer*lenIntBuffer*nThreads);
   
    double *intBuffer2 = intBuffer + nThreads*lenIntBuffer;


    // Allocate thread local storage to store integral contractions
    std::vector<std::vector<MatsT*>> AXthreads;
    MatsT *AXRaw = nullptr;
    if(nThreads != 1) {
      AXRaw = CQMemManager::get().malloc<MatsT>(nThreads*nMat*nBasis*nBasis);    
      memset(AXRaw,0,nThreads*nMat*nBasis*nBasis*sizeof(MatsT));
    }

    if(nThreads == 1) {
      AXthreads.emplace_back();
      for(auto iMat = 0; iMat < nMat; iMat++)
        AXthreads.back().push_back(matList[iMat].AX);
    } else {
      for(auto iThread = 0; iThread < nThreads; iThread++) {
        AXthreads.emplace_back();
        for(auto iMat = 0; iMat < nMat; iMat++) 
          AXthreads.back().push_back(AXRaw + iThread*nMat*nBasis*nBasis + iMat*nBasis*nBasis);
      }
    }

#ifdef _SHZ_SCREEN
    // Compute shell block norms (∞-norm) of matList.X
    // for all matrix
    size_t nShBlkNormsMat = nMat + (nMat == 1 ? 0: 1);
    double *ShBlkNorms_raw = CQMemManager::get().malloc<double>(nShBlkNormsMat*snShell*snShell);
    double *ShBlkNorms = ShBlkNorms_raw; 
    std::vector<double*> ShBlkNorms_Mat(nMat, nullptr);
    
    if (NonHermitian and not sameBasisSet12)
      CErr("EPAI Contraction does not support non-Hermitian type.");
    
    size_t ShBlkNorms_Mat_Off = nMat == 1 ? 0: 1;
    
    // #pragma omp parallel for 
    for(auto iMat = 0; iMat < nMat; iMat++) {
      
      ShBlkNorms_Mat[iMat] = ShBlkNorms_raw + (iMat+ShBlkNorms_Mat_Off)*snShell*snShell;
      
      double * ShBlkNorms_i = ShBlkNorms_Mat[iMat];
      
      ShellBlockNorm(basisSet2_.shells,matList[iMat].X,snBasis,ShBlkNorms_i);
      for(auto j = 0; j < snShell*snShell; j++)
        ShBlkNorms_i[j] = std::abs(ShBlkNorms_i[j]);

      // symmetrize nonHermitian ShBlkNorms
      if (not matList[iMat].HER) {
        for(auto k = 0; k < nShell; k++)
        for(auto l = 0; l < k;      l++) {
          double mx = std::max(ShBlkNorms_i[k + l*nShell],
                               ShBlkNorms_i[l + k*nShell]);
          ShBlkNorms_i[k + l*nShell] = mx;
          ShBlkNorms_i[l + k*nShell] = mx;
        }
      }
    }
    
    // Get the max over all the matricies for the shell block ∞-norms
    if (nMat != 1) { 
      memset(ShBlkNorms,0.,snShell*snShell*sizeof(double));
      #pragma omp parallel for 
      for(auto i = 0; i < snShell*snShell; i++) 
      for(auto iMat = 0; iMat < nMat; iMat++)
        ShBlkNorms[i] = std::max(ShBlkNorms[i], ShBlkNorms_Mat[iMat][i]);
    }
    
    // Find the max value of shell block ∞-norms of all matList.X
    double maxShBlkNorm = *std::max_element(ShBlkNorms, ShBlkNorms + snShell*snShell);

    size_t maxnPrim4 = 
      basisSet2_.maxPrim * basisSet2_.maxPrim * basisSet2_.maxPrim * 
      basisSet2_.maxPrim;

    // Set Libint precision
#if 0
    engines[0].set_precision(
      std::min(
        std::numeric_limits<double>::epsilon(),
        threshSchwarz/maxShBlkNorm
      )/maxnPrim4
    );
#else
    engines[0].set_precision(
      std::max(
        std::numeric_limits<double>::epsilon(),
        tpi.threshSchwarz()/(maxShBlkNorm*maxnPrim4))
      );
#endif

#else
    // Set Linbint precision
    engines[0].set_precision(std::numeric_limits<double>::epsilon());
#endif

    // Copy master thread engine to other threads
    for(size_t i = 1; i < nThreads; i++) engines[i] = engines[0];

#ifdef _SUB_TIMINGS
    std::chrono::duration<double> durInner(0.), durCont(0.), durSymm(0.), durZero(0.);
#endif

    // Keeping track of number of integrals and contration skipped
    std::vector<size_t> nIntSkip(nThreads,0);
#ifdef _SEPARATED_SHZ_SCREEN    
    std::vector<size_t> nConSkip(nThreads,0);
#endif
    // MPI info
    size_t mpiChunks = (nShell * (nShell + 1) / 2) / mpiSize;
    size_t mpiS12St  = mpiRank * mpiChunks;
    size_t mpiS12End = (mpiRank + 1) * mpiChunks;
    if( mpiRank == (mpiSize - 1) ) mpiS12End = (nShell * (nShell + 1) / 2);

    auto topDirect = tick();
    
    #pragma omp parallel
    {

//    ProgramTimer::setContext(parentId, callLevel);

    // Set up thread local storage

    // SMP info
    size_t thread_id = GetThreadID();

    auto &engine = engines[thread_id];
    const auto& buf_vec = engine.results();
    
    auto &AX_loc = AXthreads[thread_id];

    double * intBuffer_loc  = intBuffer  + thread_id*lenIntBuffer;
    double * intBuffer2_loc = intBuffer2 + thread_id*lenIntBuffer;

    size_t n1,n2;
    
    std::vector<size_t> contract_Mat(nMat);
    size_t iCon, nCon;

#if defined(_SHZ_SCREEN) && defined(_SEPARATED_SHZ_SCREEN)
    std::vector<double> shMax123_Mat(nMat);
    if (nMat == 1) {
      contract_Mat[0] = 0;
      nCon = 1;
    }
#else
    std::iota(contract_Mat.begin(), contract_Mat.end(), 0);
    nCon = nMat;
#endif

    // Always Loop over s2 <= s1
    for(size_t s1(0ul), bf1_s(0ul), s12(0ul); s1 < nShell; bf1_s+=n1, s1++) { 
      n1 = basisSet_.shells[s1].size(); // Size of Shell 1

    auto sigPair12_it = basisSet_.shellData.shData.at(s1).begin();
    for( const size_t& s2 : basisSet_.shellData.sigShellPair[s1] ) {
      size_t bf2_s = basisSet_.mapSh2Bf[s2];
      n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      const auto * sigPair12 = sigPair12_it->get();
      sigPair12_it++;

#ifdef CQ_ENABLE_MPI
      // MPI partition s12 blocks
      // if( (s12 < mpiS12St) or (s12 >= mpiS12End) ) { s12++; continue; }
#endif

      // Round-Robin work distribution
      if( (s12++) % nThreads != thread_id ) continue;

      // Cache variables for shells 1 and 2
        
#ifdef _FULL_DIRECT
      // Deneneracy factor for s1,s2 pair
      double s12_deg = (s1 == s2) ? 1.0 : 2.0;
#endif

#ifdef _SHZ_SCREEN
      double shz12 = 0, shMax12 = 0;
      if( screen ) {
        shz12 = schwarz1[s1 + s2*nShell];
        shMax12 = ShBlkNorms[s1 + s2*nShell];
      }
#endif

// The upper bound of s3 is s1 for the 8-fold symmetry and
// nShell for 4-fold.
#ifdef _USE_EIGHT_FOLD
  #define S3_MAX s1
#elif defined(_USE_FOUR_FOLD)
  // the "-" is for the <= in the loop
  #define S3_MAX nShell - 1
#endif

      size_t n3,n4;
      size_t s3_max = (&basisSet_ == &basisSet2_) ? S3_MAX : snShell - 1;

      for(size_t s3(0ul), bf3_s(0ul), s34(0ul); s3 <= s3_max; s3++, bf3_s += n3) { 
        n3 = basisSet2_.shells[s3].size(); // Size of Shell 3

#ifdef _SHZ_SCREEN

        double shMax123 = 0;
        if( screen and sameBasisSet12 ) {
          // Pre-calculate shell-block norm max's that only
          // depend on shells 1,2 and 3
          shMax123 = 
            std::max(ShBlkNorms[s1 + s3*nShell], 
                     ShBlkNorms[s2 + s3*nShell]);

          shMax123 = std::max(shMax123,shMax12);
#ifdef _SEPARATED_SHZ_SCREEN          
          if (nMat != 1)
          for (auto iMat = 0; iMat < nMat; iMat++)
            shMax123_Mat[iMat] = 
              std::max(ShBlkNorms_Mat[iMat][s1 + s2*nShell],
              std::max(ShBlkNorms_Mat[iMat][s1 + s3*nShell],
                       ShBlkNorms_Mat[iMat][s2 + s3*nShell]));
#endif
        }

#endif

// The upper bound of s4 is either s2 or s3 based on s1 and s3 for
// the 8-fold symmetry and s3 for the 4-fold symmetry
#ifdef _USE_EIGHT_FOLD
        size_t s4_max = (s1 == s3) ? s2 : s3;
#elif defined(_USE_FOUR_FOLD)
        size_t s4_max =  s3;
#endif

      if (&basisSet_ != &basisSet2_)
        s4_max =  s3;

      auto sigPair34_it = basisSet2_.shellData.shData.at(s3).begin();
      for( const size_t& s4 : basisSet2_.shellData.sigShellPair[s3] ) {

        if (s4 > s4_max)
          break;  // for each s3, s4 are stored in monotonically increasing
                  // order

        const auto * sigPair34 = sigPair34_it->get();
        sigPair34_it++;
                    
        size_t bf4_s = basisSet2_.mapSh2Bf[s4];
        n4 = basisSet2_.shells[s4].size(); // Size of Shell 4

        const bool sameCenter1234 =
          basisSet_.shells[s1].O == basisSet_.shells[s2].O &&
          basisSet_.shells[s1].O == basisSet_.shells[s3].O &&
          basisSet_.shells[s1].O == basisSet_.shells[s4].O;

#ifdef _SHZ_SCREEN

        double shMax = 0;

        if( screen ) {
          // Compute Shell norm max
          shMax = ShBlkNorms[s3 + s4*snShell];
          
          if (sameBasisSet12) {
            shMax = std::max(shMax,
                      std::max(ShBlkNorms[s1 + s4*nShell],
                               ShBlkNorms[s2 + s4*nShell]));
            shMax = std::max(shMax,shMax123);
          }
          
          // for same basissets, schwarz2 has been changed to schwarz1
          if((shMax * shz12 * schwarz2[s3 + s4*snShell]) <
             tpi.threshSchwarz()) { 
            nIntSkip[thread_id]++; 
#ifdef _SEPARATED_SHZ_SCREEN
            nConSkip[thread_id] += nMat;
#endif
            continue; 
          }

#ifdef _SEPARATED_SHZ_SCREEN
          if (nMat != 1) {
            nCon = 0;
            for (auto iMat = 0ul; iMat < nMat; iMat++) {
              shMax = ShBlkNorms_Mat[iMat][s3 + s4*snShell]; 
              if (sameBasisSet12) {
                shMax = std::max(shMax, 
                          std::max(ShBlkNorms_Mat[iMat][s1 + s4*nShell],
                                   ShBlkNorms_Mat[iMat][s2 + s4*nShell]));
                shMax = std::max(shMax,shMax123_Mat[iMat]);
              } 
              
              // for same basissets, schwarz2 has been changed to schwarz1
              if((shMax * shz12 * schwarz2[s3 + s4*snShell]) <
                 tpi.threshSchwarz()) { 
                 nConSkip[thread_id]++; 
              } else {   
                 contract_Mat[nCon] = iMat;
                 nCon++; 
              }
            }
          }
#endif
        }
      
#endif

#ifdef _FULL_DIRECT

        // Degeneracy factor for s3,s4 pair
        double s34_deg = (s3 == s4) ? 1.0 : 2.0;

        // Degeneracy factor for s1, s2, s3, s4 quartet
        double s12_34_deg = 2.0;
        if (&basisSet_ == &basisSet2_)
          s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;

        // Total degeneracy factor
        double s1234_deg = s12_deg * s34_deg * s12_34_deg;

#endif
#ifdef _REPORT_INTEGRAL_TIMINGS
//        ProgramTimer::tick("Direct Int Form");
#endif

#if 1
        // Evaluate ERI for shell quartet (s1 s2 | s3 s4)
        switch (eriKernel) {
          case DirectTPI<IntsT>::Kernel::Coulomb:
            engine.compute2<
              libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              basisSet_.shells[s1],
              basisSet_.shells[s2],
              basisSet2_.shells[s3],
              basisSet2_.shells[s4]
#ifdef _PRECOMPUTE_SHELL_PAIRS
              ,sigPair12,sigPair34
#endif
            );
            break;
          case DirectTPI<IntsT>::Kernel::ShortRangeErfc:
            engine.compute2<
              libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xx, 0>(
              basisSet_.shells[s1],
              basisSet_.shells[s2],
              basisSet2_.shells[s3],
              basisSet2_.shells[s4]
#ifdef _PRECOMPUTE_SHELL_PAIRS
              ,sigPair12,sigPair34
#endif
            );
            break;
          default:
            CErr("Unrecognized two-electron kernel in direct contraction.");
        }
#endif
#ifdef _REPORT_INTEGRAL_TIMINGS
//        ProgramTimer::tock("Direct Int Form");
#endif

        // Libint2 internal screening
        const double *buff = buf_vec[0];

        if(buff == nullptr) { 
          nIntSkip[thread_id]++; 
          
#ifdef _SEPARATED_SHZ_SCREEN
          nConSkip[thread_id] += nCon;
#endif          
          continue; 
        }

#ifdef _FULL_DIRECT

// Flag to turn contraction on and off
#if 1

#ifdef _REPORT_INTEGRAL_TIMINGS
//        ProgramTimer::tick("Direct Den Contract");
#endif

        // Scale the buffer by the degeneracy factor and store
        // in infBuffer
        std::transform(buff,buff + n1*n2*n3*n4,intBuffer_loc,
          [&](auto& x) { return x*0.5*s1234_deg; });

        size_t b1,b2,b3,b4;
        double *Xp1, *Xp2;
        double X1,X2;
        MatsT      T1,T2,T3,T4;
        MatsT      *Tp1,*Tp2;

        for(iCon = 0; iCon < nCon; iCon++) {
          
          auto iMat = contract_Mat[iCon]; 
          
          // Hermetian contraction
          if ( matList[iMat].HER ) {
            if( matList[iMat].contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // Cache i,j variables
              b1 = bf1 + nBasis*bf2; 
              // X is stored in basisSet2_, so b1 is valid only for the same-basis reverse digest below.
              X1 = 0.;
              if(&basisSet_ == &basisSet2_) X1 = *reinterpret_cast<double*>(matList[iMat].X + b1);
              Xp1 = reinterpret_cast<double*>(AX_loc[iMat] + b1);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // J(1,2) += I * X(4,3)
              *Xp1 += *GetRealPtr(matList[iMat].X,bf4,bf3,snBasis) * intBuffer_loc[ijkl];

              if (&basisSet_ == &basisSet2_)
                *GetRealPtr(AX_loc[iMat],bf4,bf3,nBasis) += X1 * intBuffer_loc[ijkl];

              // J(2,1) and J(3,4) are handled on symmetrization after
              // contraction
            } // kl loop
            } // ij loop

            else if( matList[iMat].contType == EXCHANGE ) {
              if (&basisSet_ != &basisSet2_)
                CErr("No exchange contraction between two different basis!", std::cout);

              if (this->oneCenterK() and not sameCenter1234) continue;

              for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
              for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
              for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

                // Cache i,j,k variables
                b1 = bf1 + bf3*nBasis;
                b2 = bf2 + bf3*nBasis;

                T1 = 0.5 * SmartConj(matList[iMat].X[b1]);
                T2 = 0.5 * SmartConj(matList[iMat].X[b2]);

              for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

                // Indicies are swapped here to loop over contiguous memory
                  
                // K(1,3) += 0.5 * I * X(2,4) = 0.5 * I * CONJ(X(4,2)) (**HER**)
                AX_loc[iMat][b1]           += 0.5 * SmartConj(matList[iMat].X[bf4+nBasis*bf2]) * intBuffer_loc[ijkl];

                // K(4,2) += 0.5 * I * X(3,1) = 0.5 * I * CONJ(X(1,3)) (**HER**)
                AX_loc[iMat][bf4 + bf2*nBasis] += T1 * intBuffer_loc[ijkl];

                // K(4,1) += 0.5 * I * X(3,2) = 0.5 * I * CONJ(X(2,3)) (**HER**)
                AX_loc[iMat][bf4 + bf1*nBasis] += T2 * intBuffer_loc[ijkl];

                // K(2,3) += 0.5 * I * X(1,4) = 0.5 * I * CONJ(X(4,1)) (**HER**)
                AX_loc[iMat][b2]           += 0.5 * SmartConj(matList[iMat].X[bf4+nBasis*bf1]) * intBuffer_loc[ijkl];

              } // l loop
              } // ijk
            }
          // Nonhermetian contraction
          } else {

            if (&basisSet_ != &basisSet2_)
              CErr("No non-Hermitian contraction between two different basis!", std::cout);

            if( matList[iMat].contType == COULOMB )
            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) { 
              // Cache i,j variables
              b1 = bf1 + nBasis*bf2; 
              T1 = *(matList[iMat].X  + b1);
              Tp1 = (AX_loc[iMat] + b1);

              b2 = bf2 + nBasis*bf1; 
              T2 = *(matList[iMat].X  + b2);
              Tp2 = (AX_loc[iMat] + b2);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) 
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // J(1,2) += I * X(4,3)
              *Tp1 += 0.5*( matList[iMat].X[bf4 + bf3*nBasis] + matList[iMat].X[bf3 + bf4*nBasis]) * intBuffer_loc[ijkl];

              // J(3,4) += I * X(2,1)
              AX_loc[iMat][bf3 + bf4*nBasis] +=  0.5*(T2+T1) * intBuffer_loc[ijkl];

              // J(2,1) += I * X(3,4)
              *Tp2 += 0.5*( matList[iMat].X[bf4 + bf3*nBasis] + matList[iMat].X[bf3 + bf4*nBasis]) * intBuffer_loc[ijkl];

              // J(4,3) += I * X(1,2)
              AX_loc[iMat][bf4 + bf3*nBasis] +=  0.5*(T2+T1) * intBuffer_loc[ijkl];

            } // kl loop
            } // ij loop

            else if( matList[iMat].contType == EXCHANGE ) {

            if (this->oneCenterK() and not sameCenter1234) continue;

            for(auto i = 0ul, bf1 = bf1_s, ijkl(0ul); i < n1; i++, bf1++)      
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Cache i,j,k variables
              b1 = bf1 + bf3*nBasis;
              b2 = bf2 + bf3*nBasis;

              T1 = 0.5 * matList[iMat].X[b1];
              T2 = 0.5 * matList[iMat].X[b2];

              b3 = bf3 + bf1*nBasis;
              b4 = bf3 + bf2*nBasis;

              T3 = 0.5 * matList[iMat].X[b3];
              T4 = 0.5 * matList[iMat].X[b4];
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

              // K(3,1) += 0.5 * I * X(4,2)
              AX_loc[iMat][b3]           += 0.5 * matList[iMat].X[bf4+nBasis*bf2] * intBuffer_loc[ijkl];

              // K(4,2) += 0.5 * I * X(3,1)
              AX_loc[iMat][bf4 + bf2*nBasis] += T3 * intBuffer_loc[ijkl];
 
              // K(4,1) += 0.5 * I * X(3,2)
              AX_loc[iMat][bf4 + bf1*nBasis] += T4 * intBuffer_loc[ijkl];

              // K(3,2) += 0.5 * I * X(4,1)
              AX_loc[iMat][b4]           += 0.5 * matList[iMat].X[bf4+nBasis*bf1] * intBuffer_loc[ijkl];

              // K(1,3) += 0.5 * I * X(2,4)
              AX_loc[iMat][b1]           += 0.5 * matList[iMat].X[bf2+nBasis*bf4] * intBuffer_loc[ijkl];

              // K(2,4) += 0.5 * I * X(1,3)
              AX_loc[iMat][bf2 + bf4*nBasis] += T1 * intBuffer_loc[ijkl];
 
              // K(1,4) += 0.5 * I * X(2,3)
              AX_loc[iMat][bf1 + bf4*nBasis] += T2 * intBuffer_loc[ijkl];

              // K(2,3) += 0.5 * I * X(1,4)
              AX_loc[iMat][b2]           += 0.5 * matList[iMat].X[bf1+nBasis*bf4] * intBuffer_loc[ijkl];

            } // l loop
            } // ijk
          }
          } // Symmetry check

        } // iMat loop

#ifdef _REPORT_INTEGRAL_TIMINGS
//        ProgramTimer::tock("Direct Den Contract");
#endif

#endif

#endif
      } // loop s4
      } // loop s3

    }; // s2
    }; // s1


    }; // OpenMP context


#ifdef _REPORT_INTEGRAL_TIMINGS
    size_t nIntSkipAcc = std::accumulate(nIntSkip.begin(),nIntSkip.end(),0);
    std::cout << "Skipped Intgral:     " << nIntSkipAcc << std::endl;
#ifdef _SEPARATED_SHZ_SCREEN    
    size_t nConSkipAcc = std::accumulate(nConSkip.begin(),nConSkip.end(),0);
    std::cout << "Skipped Contraction: " << nConSkipAcc << std::endl;
#endif

    auto durDirect = tock(topDirect);
    std::cout << "Coulomb-Exchange AO Direct Contraction took " <<  durDirect << " s\n"; 

    std::cout << std::endl;
#endif


#ifdef _FULL_DIRECT

    MatsT* SCR = CQMemManager::get().malloc<MatsT>(nBasis * nBasis);
    for( auto iMat = 0; iMat < nMat;  iMat++ ) 
    for( auto iTh  = 0; iTh < nThreads; iTh++) {
  
      if( matList[iMat].HER ) {

        MatAdd('N','C',nBasis,nBasis,MatsT(0.5),AXthreads[iTh][iMat],nBasis,MatsT(0.5),
          AXthreads[iTh][iMat],nBasis,SCR,nBasis);

        if( nThreads != 1 )
          MatAdd('N','N',nBasis,nBasis,MatsT(1.),SCR,nBasis,MatsT(1.), matList[iMat].AX,nBasis,matList[iMat].AX,nBasis);
        else
          SetMat('N',nBasis,nBasis,MatsT(1.),SCR,nBasis,matList[iMat].AX,nBasis);

      } else {
        
        if (&basisSet_ != &basisSet2_)
          CErr("No non-Hermitian contraction between two different basis!", std::cout);

        if( nThreads != 1 )
          MatAdd('N','N',nBasis,nBasis,MatsT(0.5),AXthreads[iTh][iMat],nBasis,
            MatsT(1.), matList[iMat].AX,nBasis,matList[iMat].AX,nBasis);
        else 
          blas::scal(nBasis*nBasis,MatsT(0.5),matList[iMat].AX,1);

      }

    };
    CQMemManager::get().free(SCR);
    
#else

    for( auto &C : matList ) {

      // Symmetrize J contraction
      if( C.contType == COULOMB ) 
        HerMat('L',nBasis,C.AX,nBasis);

      // Inplace transpose of K contraction
      if( C.contType == EXCHANGE ) 
        IMatCopy('C',nBasis,nBasis,MatsT(1.),C.AX,nBasis,nBasis);
  
    } // Loop over contractions

#endif


#ifdef CQ_ENABLE_MPI
    // Combine all G[X] contributions onto Root process
    if( mpiSize > 1 ) {

    //   // FIXME: This should be able to be done with MPI_IN_PLACE for
    //   // the root process
    //   MatsT* mpiScr;
    //   if( mpiRank == 0 ) mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);

    //   for( auto &C : matList ) {
//  //     prettyPrintSmart(std::cerr,"AX in Direct",C.AX,nBasis,nBasis,nBasis);

    //     MPIReduce( C.AX, nBasis*nBasis, mpiScr, 0, comm );

    //     // Copy over the output buffer on root
    //     if( mpiRank == 0 ) std::copy_n(mpiScr,nBasis*nBasis,C.AX);

    //   }

    //   if( mpiRank == 0 ) CQMemManager::get().free(mpiScr);

    }

#endif

    // Free scratch space
    CQMemManager::get().free(intBuffer);
#ifdef _SHZ_SCREEN
    CQMemManager::get().free(ShBlkNorms_raw);
#endif
    if(AXRaw) CQMemManager::get().free(AXRaw);

    // Turn threads for LA back on
    SetLAThreads(LAThreads);

//    ProgramTimer::tock("Contract Total");

  };

  template <typename MatsT, typename IntsT>
  size_t GTODirectTPIContraction<MatsT,IntsT>::directScaffoldNewSCRSize() const {

    size_t threadSCRSize  = 0ul;
    size_t generalSCRSize = 0ul; 
    
    // SCR needed for integrals
    DirectTPI<IntsT> &tpi = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& basisSet_  = this->contractSecond ? tpi.basisSet2() : tpi.basisSet();
    BasisSet& basisSet2_ = this->contractSecond ? tpi.basisSet()  : tpi.basisSet2();
    
    const size_t nBasis   = basisSet_.nBasis;
    const size_t snBasis  = basisSet2_.nBasis;
    const size_t nShell   = basisSet_.nShell;
    const size_t snShell  = basisSet2_.nShell;
    size_t nThreads  = GetNumThreads();
    
    // create a dummy engine to figure out sizes
    libint2::Engine engine(libint2::Operator::coulomb, 
      std::max(basisSet_.maxPrim, basisSet2_.maxPrim),
      std::max(basisSet_.maxL, basisSet2_.maxL),0);

    // Allocate scratch for raw integral batches
    size_t maxShellSize = 
      std::max_element(basisSet_.shells.begin(),basisSet_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    size_t maxShellSize2 = 
      std::max_element(basisSet2_.shells.begin(),basisSet2_.shells.end(),
        [](libint2::Shell &sh1, libint2::Shell &sh2) {
          return sh1.size() < sh2.size();
        })->size();

    // lenIntBuffer is allocated to be able to store EPAI's of the 
    // shell with the highest angular momentum
    // seems that lenIntBuffer is already in MatsT
    size_t lenIntBuffer = 
      maxShellSize * maxShellSize * maxShellSize2 * maxShellSize2; 

    size_t nBuffer = 2;
    
    threadSCRSize += nBuffer*lenIntBuffer; 
    
    // SCR needed for contraction storage in each thread
    if (nThreads != 1) threadSCRSize += nBasis*nBasis;

#ifdef _SHZ_SCREEN
    // 1 for general shell block ∞-norms and 1 for each matrix
    generalSCRSize += snShell*snShell*2;  

#endif

    return threadSCRSize * nThreads + generalSCRSize;
  }; // GTODirectTPIContraction::directScaffoldNewSCRSize()
  
  void GIAODirectERIContraction::twoBodyContract(
      MPI_Comm c,
      const bool screen,
      std::vector<TwoBodyContraction<dcomplex>> &list,
      EMPerturbation &pert) const {
    // Only use GIAOs if GIAOs are selected and if
    // a Magnetic field is in the EMPerturbation

    if( pert_has_type(pert,Magnetic) ) {
      directScaffold(c, screen, list, pert);
    } else {
      directScaffoldNew(c, screen, list);
    }
  }

#ifdef _PROFILE_DIRECT_GRAD
  namespace {
    struct DirectGradProfile {
      // Counters
      size_t nQuartetsVisited  = 0;
      size_t nQuartetsSkipped  = 0;   // wired in for screening; 0 in this version
      size_t nDerivBufsNonNull = 0;   // # of non-null buffers (out of 12 per quartet) summed
      size_t nDerivIntegrals   = 0;   // sum of n1*n2*n3*n4 over non-null buffers
      size_t nJOps             = 0;   // FMA count for Coulomb contractions
      size_t nKOps             = 0;   // FMA count for Exchange contractions
      // Timers (seconds)
      double tCompute2         = 0.0; // engine.compute2 only
      double tScale            = 0.0; // scale-and-copy of derivative buffer
      double tContraction      = 0.0; // J + K kernels
      char _pad[56];                  // pad 9*8=72 to 128 (2 cache lines, no false sharing)
    };
    static_assert(sizeof(DirectGradProfile) == 128, "DirectGradProfile size should be 128 bytes for cache alignment");
  }
#endif

  template <typename MatsT, typename IntsT>
  void DirectGradContraction<MatsT,IntsT>::directScaffoldGrad(
      MPI_Comm comm,
      const bool screen,
      std::vector<std::vector<TwoBodyContraction<MatsT>>>& cList) const {

    directScaffoldGradImpl(comm, screen, cList, {}, {}, nullptr);

  }

  template <typename MatsT, typename IntsT>
  void DirectGradContraction<MatsT,IntsT>::directScaffoldGradImpl(
      MPI_Comm comm,
      const bool screen,
      std::vector<std::vector<TwoBodyContraction<MatsT>>>& cList,
      const std::vector<const MatsT*>& traceDensities,
      const std::vector<double>& traceCoeffs,
      std::vector<double>* gradientOut) const {

    // Determine mode/output here.
    // If gradientOut is not nullPtr, directly trace with density for each dF/dI and output gradient
    // Otherwise, output each dF/dI and trace with density outside
    const bool TRACE = (gradientOut != nullptr);

    // -----------------------------------------------------------------
    // SCREENING DESIGN NOTE (Horn, Weiss, Haeser, Ehrig, Ahlrichs,
    // J. Comput. Chem. 12 (1991) 1058, DOI 10.1002/jcc.540120903)
    //
    // Quartet skip test (Horn eq. 21a):
    //   |E^i_{νμχλ}| <= (Q_NM R_KL + Q_KL R_NM) D_{NM,KL}
    //   D_{NM,KL}   =  4 D_NM D_KL + D_NK D_ML + D_NL D_MK   (Horn eq. 14)
    //
    // This is the chart-1 (energy-gradient) density weight, not the
    // chart-2 gradient-Fock weight D~ of eq. 15. That is deliberate: this
    // routine builds gradient Fock matrices F^I, but every caller only
    // traces them against the same density used in the contraction, so
    // screening the traced (energy-gradient) contribution is both tighter
    // and sufficient. Cross-basis (NEO ep) callers must supply that trace
    // density via GradContractions::traceDensity.
    //
    // Rigor of the integral bound (Horn eqs. 8-11): for every libint
    // derivative buffer the element-wise bound holds,
    //     |buf_vec[0..5][ijkl]|  <= R_NM * Q_KL   (derivs on shells 1,2)
    //     |buf_vec[6..11][ijkl]| <= Q_NM * R_KL   (derivs on shells 3,4)
    // where R is the FD-based derivative Schwarz bound assembled in
    // DirectTPI::computeSchwarzGrad (rigorous up to O(h^2) FD truncation,
    // ~1e-4 relative at h = 1e-4). Both Q and R are geometry dependent and
    // are refreshed each gradient evaluation by
    // GradInts<TwoPInts,double>::computeAOInts.
    //
    // The bound is additionally multiplied by n1*n2*n3*n4 and the quartet
    // degeneracy so it bounds the *summed* block contribution to a single
    // gradient component; the caller-side 0.25 energy prefactor is not
    // included, making the test conservative by ~4x.
    // -----------------------------------------------------------------

    DirectTPI<IntsT> &tpi = dynamic_cast<DirectTPI<IntsT>&>(*this->grad_[0]);
    BasisSet& basisSet_  = this->contractSecond ? tpi.basisSet2() : tpi.basisSet();
    BasisSet& basisSet2_ = this->contractSecond ? tpi.basisSet()  : tpi.basisSet2();
    const typename DirectTPI<IntsT>::Kernel eriKernel = tpi.kernel();

    size_t nThreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank   = MPIRank(comm);
    size_t mpiSize   = MPISize(comm);

    SetLAThreads(1); // Turn off parallelism in LA functions

#ifdef _PROFILE_DIRECT_GRAD
    using _prof_clock = std::chrono::steady_clock;
    auto _prof_t_start         = _prof_clock::now();
    auto _prof_t_setup_end     = _prof_t_start;
    auto _prof_t_parallel_end  = _prof_t_start;
    std::vector<DirectGradProfile> _prof(nThreads);
#endif

    const size_t nBasis   = basisSet_.nBasis;
    const size_t snBasis  = basisSet2_.nBasis;
    const size_t nShell   = basisSet_.nShell;
    const size_t snShell  = basisSet2_.nShell;
    const size_t nTotGrad = cList.size();
    const size_t nMat     = cList[0].size();

    bool NonHermitian = false;
    for (auto& gradComp: cList)
      for (auto& x: gradComp)
        NonHermitian |= not x.HER;

    const bool sameBasisSet12 = (&basisSet_ == &basisSet2_);

    if (TRACE) {
      if (traceDensities.size() != nMat or traceCoeffs.size() != nMat)
        CErr("gradTwoBodyTraceContract needs one trace density and one trace "
             "coefficient per contraction.", std::cout);
      // The hermitization equivalence above only holds for Hermitian A and D
      if (NonHermitian)
        CErr("Trace-mode gradient contraction requires Hermitian densities.",
             std::cout);
      gradientOut->assign(nTotGrad, 0.);
    }

#ifdef _PROFILE_DIRECT_GRAD
    // Precompute J/K op factors per non-null gradient buffer
    // (so per-quartet ops = nNonNull * n1*n2*n3*n4 * factor)
    size_t _prof_jFactor = 0, _prof_kFactor = 0;
    for (auto& _m : cList[0]) {
      if (_m.HER) {
        if      (_m.contType == COULOMB)  _prof_jFactor += sameBasisSet12 ? 2 : 1;
        else if (_m.contType == EXCHANGE) _prof_kFactor += 4;
      }
    }
#endif

    // ----------------- Screening setup -----------------
    double *Q1 = nullptr, *Q2 = nullptr;
    double *R1 = nullptr, *R2 = nullptr;

    // Ket-side (basisSet2_) block norms, from the contraction density X.
    // For same-basis this also serves as the bra-side (D_NM) array.
    double *ShBlkNorms_raw = nullptr;
    double *ShBlkNorms     = nullptr;
    std::vector<double*> ShBlkNorms_Mat;

    // Bra-side (basisSet_) block norms, from the TRACE density the caller
    // dots AX against. Same-basis: aliases ShBlkNorms. Cross-basis (NEO):
    // a distinct single-matrix array (nMat == 1), REQUIRED for a correct
    // energy-gradient bound.
    double *ShBlkNorms1_raw = nullptr;
    double *ShBlkNorms1     = nullptr;

    double maxShBlkNorm = 0.0;

    if (screen) {
      // Make sure both Schwarz arrays are computed
      if (tpi.schwarz()     == nullptr) tpi.computeSchwarz();
      if (tpi.schwarzGrad() == nullptr) tpi.computeSchwarzGrad();

      Q1 = this->contractSecond ? tpi.schwarz2()     : tpi.schwarz();
      Q2 = this->contractSecond ? tpi.schwarz()      : tpi.schwarz2();
      R1 = this->contractSecond ? tpi.schwarzGrad2() : tpi.schwarzGrad();
      R2 = this->contractSecond ? tpi.schwarzGrad()  : tpi.schwarzGrad2();
      if (sameBasisSet12) { Q2 = Q1; R2 = R1; }

      // Shell-block ∞-norms for each density matrix in the contraction.
      // Densities are shared across iGrad (only AX differs), so cList[0]
      // gives all the densities we need.
      const size_t nSlots = nMat + (nMat == 1 ? 0 : 1);
      const size_t off    = (nMat == 1) ? 0 : 1;

      // ---- ket-side (basisSet2_, snShell x snShell) from X ----
      ShBlkNorms_raw = CQMemManager::get().malloc<double>(nSlots * snShell * snShell);
      ShBlkNorms = ShBlkNorms_raw;
      ShBlkNorms_Mat.assign(nMat, nullptr);
      for (size_t iMat = 0; iMat < nMat; iMat++) {
        ShBlkNorms_Mat[iMat] = ShBlkNorms_raw + (iMat + off) * snShell * snShell;
        ShellBlockNorm(basisSet2_.shells, cList[0][iMat].X, snBasis, ShBlkNorms_Mat[iMat]);
        for (size_t j = 0; j < snShell * snShell; j++)
          ShBlkNorms_Mat[iMat][j] = std::abs(ShBlkNorms_Mat[iMat][j]);
      }

      // Max over matrices into the front-slot ShBlkNorms (skip when nMat==1
      // since ShBlkNorms then aliases ShBlkNorms_Mat[0]).
      if (nMat != 1) {
        std::memset(ShBlkNorms, 0, snShell * snShell * sizeof(double));
        #pragma omp parallel for
        for (size_t i = 0; i < snShell * snShell; i++)
        for (size_t iMat = 0; iMat < nMat; iMat++)
          ShBlkNorms[i] = std::max(ShBlkNorms[i], ShBlkNorms_Mat[iMat][i]);
      }

      // ---- bra-side (basisSet_, nShell x nShell) ----
      if (sameBasisSet12) {
        ShBlkNorms1 = ShBlkNorms;          // X is also the traced density
      } else {
        if (this->traceDensity == nullptr)
          CErr("Cross-basis (NEO) gradient screening needs traceDensity "
              "(the density AX is traced against, over basisSet_).", std::cout);
        if (nMat != 1)
          CErr("Cross-basis gradient screening assumes a single Coulomb "
              "contraction (nMat == 1).", std::cout);

        ShBlkNorms1_raw = CQMemManager::get().malloc<double>(nShell * nShell);
        ShBlkNorms1 = ShBlkNorms1_raw;
        ShellBlockNorm(basisSet_.shells, this->traceDensity->S().pointer(),
                      nBasis, ShBlkNorms1);
        for (size_t j = 0; j < nShell * nShell; j++)
          ShBlkNorms1[j] = std::abs(ShBlkNorms1[j]);
      }

      maxShBlkNorm = *std::max_element(ShBlkNorms, ShBlkNorms + snShell * snShell);
      if (!sameBasisSet12)
        maxShBlkNorm = std::max(maxShBlkNorm,
            *std::max_element(ShBlkNorms1, ShBlkNorms1 + nShell * nShell));
    }

    // Global ket-side maxima for the pair-level early-out (Horn chart 1,
    // eqs. 16-19). Each is an upper bound of the corresponding per-quartet
    // quantity, so a pair skipped here would have had every one of its
    // quartets skipped by the quartet-level test -- results are identical.
    double maxQ2 = 0.0, maxR2 = 0.0, maxD2 = 0.0, maxN34 = 0.0;
    if (screen) {
      maxQ2 = *std::max_element(Q2, Q2 + snShell * snShell);
      maxR2 = *std::max_element(R2, R2 + snShell * snShell);
      maxD2 = *std::max_element(ShBlkNorms, ShBlkNorms + snShell * snShell);
      maxN34 = static_cast<double>(
        std::max_element(basisSet2_.shells.begin(), basisSet2_.shells.end(),
          [](libint2::Shell &sh1, libint2::Shell &sh2) {
            return sh1.size() < sh2.size();
          })->size());
    }
    // ----------------- End screening setup -----------------

    std::vector<libint2::Engine> engines(nThreads);

    // Construct engine for master thread
    engines[0] = libint2::Engine(tpi.libintOperator(),
      std::max(basisSet_.maxPrim, basisSet2_.maxPrim), 
      std::max(basisSet_.maxL, basisSet2_.maxL),1);
    if (eriKernel == DirectTPI<IntsT>::Kernel::ShortRangeErfc)
      engines[0].set_params(tpi.rangeSeparationParameter());

    // 12 derivatives per integral (3 xyz * 4 shells)
    size_t nGrad = 12;

    // Allocate thread local storage to store integral contractions
    // Threads, Gradients, Matrices, Basis, Basis
    std::vector<std::vector<std::vector<MatsT*>>> AXthreads;
    MatsT *AXRaw = nullptr;

    std::vector<std::vector<double>> gradThreads;

    if (TRACE) {

    gradThreads.assign(nThreads, {});
    AXthreads.assign(nThreads, {});

    } else {

    if(nThreads != 1) {
      AXRaw = CQMemManager::get().malloc<MatsT>(nTotGrad*nThreads*nMat*nBasis*nBasis);
      const size_t totalSize = nTotGrad*nThreads*nMat*nBasis*nBasis;
      #pragma omp parallel for schedule(static)
      for (size_t k = 0; k < totalSize; k++) {
        AXRaw[k] = MatsT(0);
      }
    }

    if(nThreads == 1) {
      AXthreads.emplace_back();
      for(auto& gradComp: cList) {
        AXthreads.back().emplace_back();
        for(auto& mat: gradComp) AXthreads.back().back().push_back(mat.AX);
      }
    } else {
      for(auto iThread = 0; iThread < nThreads; iThread++) {
        AXthreads.emplace_back();
        for(auto iGrad = 0; iGrad < nTotGrad; iGrad++) {
          AXthreads.back().emplace_back();
          for(auto iMat = 0; iMat < nMat; iMat++)
            AXthreads.back().back().push_back(
              AXRaw +
              iThread*nMat*nBasis*nBasis*nTotGrad +
              iGrad*nMat*nBasis*nBasis +
              iMat*nBasis*nBasis
            );
        }
      }
    }

    } // if (TRACE) ... else

    // Set Linbint precision
    engines[0].set_precision(std::numeric_limits<double>::epsilon());

    // Copy master thread engine to other threads
    for(size_t i = 1; i < nThreads; i++) engines[i] = engines[0];

    // Keeping track of number of integrals skipped
    std::vector<size_t> nSkip(nThreads,0);

#ifdef _PROFILE_DIRECT_GRAD
    _prof_t_setup_end = _prof_clock::now();
#endif

    //
    // Parallel region - start work
    //
    #pragma omp parallel
    {

    // Set up thread local storage

    // SMP info
    size_t thread_id = GetThreadID();

    auto &engine = engines[thread_id];
    const auto& buf_vec = engine.results();
    
    auto &AX_loc = AXthreads[thread_id];

    double *gradLoc = nullptr;
    if (TRACE) {
      gradThreads[thread_id].assign(nTotGrad, 0.);
      gradLoc = gradThreads[thread_id].data();
    }

    size_t n1,n2;
    size_t shell_atoms[4];

    // Always Loop over s2 <= s1
    for(size_t s1(0ul), bf1_s(0ul); s1 < nShell; bf1_s+=n1, s1++) { 

      n1 = basisSet_.shells[s1].size(); // Size of Shell 1
      shell_atoms[0] = basisSet_.mapSh2Cen[s1]; // Atomic center of shell 1

    auto sigPair12_it = basisSet_.shellData.shData.at(s1).begin();
    for( const size_t& s2 : basisSet_.shellData.sigShellPair[s1] ) {
      size_t bf2_s = basisSet_.mapSh2Bf[s2];

      n2 = basisSet_.shells[s2].size(); // Size of Shell 2
      shell_atoms[1] = basisSet_.mapSh2Cen[s2]; // Atomic center of shell 2

      const auto * sigPair12 = sigPair12_it->get();
      sigPair12_it++;

      // Deterministic work distribution over (MPI rank, thread): map each
      // (s1,s2) pair to one of mpiSize*nThreads global workers. globalId maps
      // 1:1 to (rank, thread), so every pair is processed exactly once. For
      // mpiSize == 1 this is identical to the old (s1*nShell+s2) % nThreads.
      const size_t s12id    = s1 * nShell + s2;
      const size_t globalId = s12id % (mpiSize * nThreads);
      if ( globalId / nThreads != mpiRank ||
           globalId % nThreads != thread_id ) continue;

#ifdef _FULL_DIRECT
      // Deneneracy factor for s1,s2 pair
      double s12_deg = (s1 == s2) ? 1.0 : 2.0;
#endif

// The upper bound of s3 is s1 for the 8-fold symmetry and
// nShell for 4-fold.

      // (s1,s2) screening quantities — hoisted out of s3,s4 loops
      double Q12 = 0.0, R12 = 0.0, D12 = 0.0;
      if (screen) {
        Q12 = Q1[s1 + s2*nShell];
        R12 = R1[s1 + s2*nShell];
        D12 = ShBlkNorms1[s1 + s2*nShell];   // D_NM (bra-side density)

        // Pair-level early-out (Horn chart 1): upper-bound the quartet test
        // over all possible (s3,s4). Density weight: same-basis
        // D_{NM,KL} <= 4 D12 maxD + 2 maxD^2, cross-basis D12 * maxD.
        // Block factor: n1*n2*maxN34^2 and degeneracies s12d*s34d*s1234d
        // <= s12d*4.
        const double s12d = (s1 == s2) ? 1.0 : 2.0;
        const double maxDWeight = sameBasisSet12
            ? 4.0 * D12 * maxD2 + 2.0 * maxD2 * maxD2
            : D12 * maxD2;
        const double maxBlockFac =
            static_cast<double>(n1 * n2) * maxN34 * maxN34 * s12d * 4.0;
        if ((R12 * maxQ2 + Q12 * maxR2) * maxDWeight * maxBlockFac < tpi.threshSchwarz())
          continue;
      }

#ifdef _USE_EIGHT_FOLD
  #define S3_MAX s1
#elif defined(_USE_FOUR_FOLD)
  // the "-" is for the <= in the loop
  #define S3_MAX nShell - 1
#endif


      size_t n3,n4;
      size_t s3_max = (&basisSet_ == &basisSet2_) ? S3_MAX : snShell - 1;

      for(size_t s3(0ul), bf3_s(0ul), s34(0ul); s3 <= s3_max; s3++, bf3_s += n3) { 

        n3 = basisSet2_.shells[s3].size(); // Size of Shell 3
        shell_atoms[2] = basisSet2_.mapSh2Cen[s3]; // Atomic center of shell 3

        // (s1,s2,s3) screening — hoist D[N,K] and D[M,K] out of s4 loop.
        // (Only needed when sameBasisSet12; otherwise K-channel terms are absent.)
        double D_NK = 0.0, D_MK = 0.0;
        if (screen && sameBasisSet12) {
          D_NK = ShBlkNorms[s1 + s3*nShell];
          D_MK = ShBlkNorms[s2 + s3*nShell];
        }

// The upper bound of s4 is either s2 or s3 based on s1 and s3 for
// the 8-fold symmetry and s3 for the 4-fold symmetry
#ifdef _USE_EIGHT_FOLD
        size_t s4_max = (s1 == s3) ? s2 : s3;
#elif defined(_USE_FOUR_FOLD)
        size_t s4_max =  s3;
#endif
      if (&basisSet_ != &basisSet2_)
        s4_max =  s3;

      auto sigPair34_it = basisSet2_.shellData.shData.at(s3).begin();
      for( const size_t& s4 : basisSet2_.shellData.sigShellPair[s3] ) {

        if (s4 > s4_max)
          break;  // for each s3, s4 are stored in monotonically increasing
                  // order

        const auto * sigPair34 = sigPair34_it->get();
        sigPair34_it++;
                    
        size_t bf4_s = basisSet2_.mapSh2Bf[s4];

        n4 = basisSet2_.shells[s4].size(); // Size of Shell 4
        shell_atoms[3] = basisSet2_.mapSh2Cen[s4]; // Atomic center of shell 4
 
        // ----------- Quartet-level screening test (Horn eq. 21a) -----------
        if (screen) {
          // Shell-block density norm for shells (K,L)
          const double D_KL = ShBlkNorms[s3 + s4*snShell];   // ket-side (basisSet2_)

          double D_NM_KL = 0.0;
          if (sameBasisSet12) {
            // Same-basis: full gradient density weight, eq. 14
            //   D = 4 D_NM D_KL + D_NK D_ML + D_NL D_MK
            const double D_NL = ShBlkNorms[s1 + s4*nShell];
            const double D_ML = ShBlkNorms[s2 + s4*snShell];
            D_NM_KL = 4.0 * D12 * D_KL + D_NK * D_ML + D_NL * D_MK;
          } else {
            // Cross-basis (NEO): distinguishable particles -> Coulomb only,
            // no exchange, no factor of 4. Energy-gradient weight is
            //   D = D_NM(bra) * D_KL(ket)
            D_NM_KL = D12 * D_KL;
          }

          const double Q34 = Q2[s3 + s4*snShell];
          const double R34 = R2[s3 + s4*snShell];

          // n1n2n3n4 * degeneracy (mirrors contraction weight, incl.
          // s12_34_deg = 2 for cross-basis)
          const double s12d = (s1 == s2) ? 1.0 : 2.0;
          const double s34d = (s3 == s4) ? 1.0 : 2.0;
          double s1234d = 2.0;
          if (sameBasisSet12)
            s1234d = (s1 == s3) ? ((s2 == s4) ? 1.0 : 2.0) : 2.0;
          const double blockFac =
              static_cast<double>(n1*n2*n3*n4) * s12d * s34d * s1234d;

          if ((R12 * Q34 + Q12 * R34) * D_NM_KL * blockFac < tpi.threshSchwarz()) {
        #ifdef _PROFILE_DIRECT_GRAD
            _prof[thread_id].nQuartetsSkipped++;
        #endif
            continue;
          }
        }

#ifdef _FULL_DIRECT

        // Degeneracy factor for s3,s4 pair
        double s34_deg = (s3 == s4) ? 1.0 : 2.0;

        // Degeneracy factor for s1, s2, s3, s4 quartet
        double s12_34_deg = 2.0;
        if (&basisSet_ == &basisSet2_)
          s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;

        // Total degeneracy factor and contraction weight
        double s1234_deg     = s12_deg * s34_deg * s12_34_deg;
        const double w       = 0.5 * s1234_deg;     // J prefactor
        const double w_half  = 0.5 * w;             // K prefactor (was 0.5 * w*I)
#endif

#ifdef _PROFILE_DIRECT_GRAD
        auto _prof_t0 = _prof_clock::now();
#endif

        switch (eriKernel) {
          case DirectTPI<IntsT>::Kernel::Coulomb:
            engine.compute2<
              libint2::Operator::coulomb, libint2::BraKet::xx_xx, 1>(
              basisSet_.shells[s1],
              basisSet_.shells[s2],
              basisSet2_.shells[s3],
              basisSet2_.shells[s4]);
            break;
          case DirectTPI<IntsT>::Kernel::ShortRangeErfc:
            engine.compute2<
              libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xx, 1>(
              basisSet_.shells[s1],
              basisSet_.shells[s2],
              basisSet2_.shells[s3],
              basisSet2_.shells[s4]);
            break;
          default:
            CErr("Unrecognized two-electron kernel in direct gradient contraction.");
        }

        // libint internal screening: buf_vec[0] == nullptr signals the whole
        // derivative set was screened. libint leaves buf_vec[1..11] pointing at
        // a PREVIOUS quartet (stale). 
        // Thus here we should skip the entire thing
        if (buf_vec[0] == nullptr) {
#ifdef _PROFILE_DIRECT_GRAD
          _prof[thread_id].nQuartetsSkipped++;
#endif
          continue;
        }

#ifdef _PROFILE_DIRECT_GRAD
        auto _prof_t1 = _prof_clock::now();
        auto _prof_t2 = _prof_clock::now();
        size_t _prof_nNonNull = 0;
        for (size_t d = 0; d < buf_vec.size(); d++)
          if (buf_vec[d] != nullptr) _prof_nNonNull++;
        _prof[thread_id].tCompute2         += std::chrono::duration<double>(_prof_t1 - _prof_t0).count();
        _prof[thread_id].tScale            += std::chrono::duration<double>(_prof_t2 - _prof_t1).count();
        _prof[thread_id].nQuartetsVisited  += 1;
        _prof[thread_id].nDerivBufsNonNull += _prof_nNonNull;
        _prof[thread_id].nDerivIntegrals   += _prof_nNonNull * n1*n2*n3*n4;
#endif

        size_t b1,b2,b3,b4;
        double *Xp1, *Xp2;
        double X1,X2;
        MatsT  T1,T2,T3,T4;
        MatsT  *Tp1,*Tp2;

        // Loop over gradient components
        for ( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {

        // Libint internal screening (for each gradient component)
        const double* buff = buf_vec[iGrad];
        if ( buff == nullptr ) continue;

        const size_t xyz = iGrad % 3; // Cartesian component of gradient
        const size_t iSh = iGrad / 3; // Shell on which the gradient is taken

        // Gradient component that is relevant for this contraction
        const size_t gComp = shell_atoms[iSh]*3 + xyz;
        std::vector<TwoBodyContraction<MatsT>>& gradList = cList[gComp];

        if (TRACE) {

        for(size_t iMat = 0; iMat < nMat; iMat++) {

          const MatsT* __restrict__ Xmat = gradList[iMat].X;
          const MatsT* __restrict__ Dmat = traceDensities[iMat];
          double acc  = 0.;
          double pref = 0.;

          if ( gradList[iMat].contType == COULOMB ) {

            pref = w;
            size_t ijkl = 0ul;
            for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) {

              // J(1,2) += w * I * Re X(4,3)  ->  weighted by Re D(1,2)
              const double DR12 = std::real(Dmat[bf1 + bf2*nBasis]);

              if (sameBasisSet12) {
                // J(4,3) += w * I * Re X(1,2)  ->  weighted by Re D(4,3)
                const double XR12 = std::real(Xmat[bf1 + bf2*nBasis]);
                for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++)
                for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++)
                  acc += ( DR12 * std::real(Xmat[bf4 + bf3*snBasis])
                         + XR12 * std::real(Dmat[bf4 + bf3*nBasis]) )
                         * buff[ijkl];
              } else {
                for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++)
                for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++)
                  acc += DR12 * std::real(Xmat[bf4 + bf3*snBasis]) * buff[ijkl];
              }

            } // ij loop

          } else if( gradList[iMat].contType == EXCHANGE ) {

            if (not sameBasisSet12)
              CErr("No exchange contraction between two different basis!", std::cout);

            pref = w_half;
            size_t ijkl = 0ul;
            for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

              // Loop invariant in l
              const MatsT D13 = Dmat[bf1 + bf3*nBasis];
              const MatsT X13 = Xmat[bf1 + bf3*nBasis];
              const MatsT D23 = Dmat[bf2 + bf3*nBasis];
              const MatsT X23 = Xmat[bf2 + bf3*nBasis];

            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) {

              // K(1,3) += w_half * I * conj X(4,2)   ->  weighted by conj D(1,3)
              // K(4,2) += w_half * I * conj X(1,3)   ->  weighted by conj D(4,2)
              // K(4,1) += w_half * I * conj X(2,3)   ->  weighted by conj D(4,1)
              // K(2,3) += w_half * I * conj X(4,1)   ->  weighted by conj D(2,3)
              acc += ( ReProd(D13, Xmat[bf4 + bf2*nBasis])
                     + ReProd(Dmat[bf4 + bf2*nBasis], X13)
                     + ReProd(Dmat[bf4 + bf1*nBasis], X23)
                     + ReProd(D23, Xmat[bf4 + bf1*nBasis]) ) * buff[ijkl];

            } // l loop
            } // ijk

          } // EXCHANGE

          gradLoc[gComp] += traceCoeffs[iMat] * pref * acc;

        } // Matrices

        } else {

        // Thread local storage for this contraction
        auto& AX_Grad_loc = AX_loc[gComp];


        // loop over matrices in contraction
        for(auto iMat = 0; iMat < nMat; iMat++) {

          // Hermetian contraction
          if( gradList[iMat].HER ) { 
            if ( gradList[iMat].contType == COULOMB ) {
            // loop over basis functions in the shell quartet
            size_t ijkl = 0ul; // *** fixed: no re‑initialisation in i‑loop ***
            for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)
            for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++) {
              // Cache i,j variables
              b1 = bf1 + nBasis*bf2;
              // X is stored in basisSet2_, so b1 is valid only for the same-basis reverse digest below.
              X1 = 0.;
              if(&basisSet_ == &basisSet2_) X1 = w * (*reinterpret_cast<double*>(gradList[iMat].X + b1));
              Xp1 = reinterpret_cast<double*>(AX_Grad_loc[iMat] + b1);
            for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++)
            for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) {

              // J(1,2) += w * I * X(4,3)
              *Xp1 += w * (*GetRealPtr(gradList[iMat].X,bf4,bf3,snBasis)) * buff[ijkl];

              // J(4,3) += w * I * X(1,2)   (w already baked into X1)
              if (&basisSet_ == &basisSet2_)
                *GetRealPtr(AX_Grad_loc[iMat],bf4,bf3,nBasis) += X1 * buff[ijkl];

              // J(2,1) and J(3,4) are handled on symmetrization after
              // contraction
            } // kl loop
            } // ij loop

            } else if( gradList[iMat].contType == EXCHANGE ) {
              if (&basisSet_ != &basisSet2_)
                CErr("No exchange contraction between two different basis!", std::cout);
              size_t ijkl = 0ul; // *** fixed counter here as well ***
              for(auto i = 0ul, bf1 = bf1_s; i < n1; i++, bf1++)      
              for(auto j = 0ul, bf2 = bf2_s; j < n2; j++, bf2++)       
              for(auto k = 0ul, bf3 = bf3_s; k < n3; k++, bf3++) {

                // Cache i,j,k variables
                b1 = bf1 + bf3*nBasis;
                b2 = bf2 + bf3*nBasis;

                T1 = w_half * SmartConj(gradList[iMat].X[b1]);
                T2 = w_half * SmartConj(gradList[iMat].X[b2]);

              for(auto l = 0ul, bf4 = bf4_s; l < n4; l++, bf4++, ijkl++) { 

                // K(1,3) += w_half * I * X(2,4)  =  w_half * I * CONJ(X(4,2))
                AX_Grad_loc[iMat][b1]               += w_half * SmartConj(gradList[iMat].X[bf4+nBasis*bf2]) * buff[ijkl];

                // K(4,2) += w_half * I * X(3,1)  =  w_half * I * CONJ(X(1,3))   (w_half in T1)
                AX_Grad_loc[iMat][bf4 + bf2*nBasis] += T1 * buff[ijkl];

                // K(4,1) += w_half * I * X(3,2)  =  w_half * I * CONJ(X(2,3))   (w_half in T2)
                AX_Grad_loc[iMat][bf4 + bf1*nBasis] += T2 * buff[ijkl];

                // K(2,3) += w_half * I * X(1,4)  =  w_half * I * CONJ(X(4,1))
                AX_Grad_loc[iMat][b2]               += w_half * SmartConj(gradList[iMat].X[bf4+nBasis*bf1]) * buff[ijkl];

              } // l loop
              } // ijk
            } // EXCHANGE

          // Nonhermitian
          } else {
            CErr("Nonhermetian NYI!");

          } // Symmetry

        } // Matrices

        } // if (TRACE) ... else


        } // Gradient components

#ifdef _PROFILE_DIRECT_GRAD
        auto _prof_t3 = _prof_clock::now();
        _prof[thread_id].tContraction += std::chrono::duration<double>(_prof_t3 - _prof_t2).count();
        _prof[thread_id].nJOps        += _prof_nNonNull * n1*n2*n3*n4 * _prof_jFactor;
        _prof[thread_id].nKOps        += _prof_nNonNull * n1*n2*n3*n4 * _prof_kFactor;
#endif

      } // s4
      } // s3

    } // s2
    } // s1
    
    } // omp parallel

#ifdef _PROFILE_DIRECT_GRAD
    _prof_t_parallel_end = _prof_clock::now();
#endif

    if (TRACE) {

    // Post-parallel reduction
    for(size_t iThread = 0ul; iThread < nThreads; iThread++) {
      if (gradThreads[iThread].empty()) continue;
      for(size_t iGrad = 0ul; iGrad < nTotGrad; iGrad++)
        (*gradientOut)[iGrad] += gradThreads[iThread][iGrad];
    }

    } else {

    // Post-parallel reduction + symmetrization
    #pragma omp parallel for collapse(2) schedule(static)
    for(size_t iGrad = 0ul; iGrad < nTotGrad; iGrad++)
    for(size_t iMat  = 0ul; iMat  < nMat;     iMat++ ) {

      MatsT* __restrict__ out = cList[iGrad][iMat].AX;

      // Sum thread-local contributions into 'out'
      if (nThreads > 1) {
        for (size_t iThread = 0; iThread < nThreads; iThread++) {
          const MatsT* __restrict__ src = AXthreads[iThread][iGrad][iMat];
          for (size_t k = 0; k < nBasis*nBasis; k++) out[k] += src[k];
        }
      }

      // Finalize
      if (cList[iGrad][iMat].HER) {
        // In-place hermitize: out := 0.5 * (out + out^H)
        for (size_t j = 0; j < nBasis; j++) {
          // Diagonal: zero out imaginary part for complex
          MatsT d = out[j + nBasis*j];
          out[j + nBasis*j] = MatsT(0.5) * (d + SmartConj(d));
          // Off-diagonal pair (i,j) and (j,i) updated together
          for (size_t i = j+1; i < nBasis; i++) {
            MatsT a = out[i + nBasis*j];   // (i,j)
            MatsT b = out[j + nBasis*i];   // (j,i)
            out[i + nBasis*j] = MatsT(0.5) * (a + SmartConj(b));
            out[j + nBasis*i] = MatsT(0.5) * (b + SmartConj(a));
          }
        }
      } else {
        // In-place scale by 0.5
        for (size_t k = 0; k < nBasis*nBasis; k++) out[k] *= MatsT(0.5);
      }
    }

    } // if (TRACE) ... else

    if (AXRaw          != nullptr) CQMemManager::get().free(AXRaw);
    if (ShBlkNorms_raw != nullptr) CQMemManager::get().free(ShBlkNorms_raw);
    if (ShBlkNorms1_raw != nullptr) CQMemManager::get().free(ShBlkNorms1_raw);

#ifdef CQ_ENABLE_MPI
    // Combine gradient-Fock contributions across MPI ranks onto root.
    // The hermitization above is linear, so summing the locally-hermitized
    // partials equals hermitizing the global sum.
    if (mpiSize > 1) {
      if (TRACE) {
        std::vector<double> mpiScr(nTotGrad, 0.);
        MPIAllReduce(gradientOut->data(), static_cast<int>(nTotGrad),
                     mpiScr.data(), comm);
        *gradientOut = std::move(mpiScr);
      } else {
      MatsT* mpiScr = nullptr;
      if (mpiRank == 0) mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);
      for (size_t iGrad = 0; iGrad < nTotGrad; iGrad++)
      for (size_t iMat  = 0; iMat  < nMat;     iMat++) {
        MPIReduce(cList[iGrad][iMat].AX, nBasis*nBasis, mpiScr, 0, comm);
        if (mpiRank == 0) std::copy_n(mpiScr, nBasis*nBasis, cList[iGrad][iMat].AX);
      }
      if (mpiRank == 0) CQMemManager::get().free(mpiScr);
      }
    }
#endif
    
    // Turn threads for LA back on
    SetLAThreads(LAThreads);

#ifdef _PROFILE_DIRECT_GRAD
    auto _prof_t_end = _prof_clock::now();

    if (mpiRank == 0) {
      DirectGradProfile total;
      double tMaxC2 = 0.0, tMinC2 =  std::numeric_limits<double>::infinity();
      double tMaxSc = 0.0, tMinSc =  std::numeric_limits<double>::infinity();
      double tMaxCt = 0.0, tMinCt =  std::numeric_limits<double>::infinity();
      for (auto& p : _prof) {
        total.nQuartetsVisited  += p.nQuartetsVisited;
        total.nQuartetsSkipped  += p.nQuartetsSkipped;
        total.nDerivBufsNonNull += p.nDerivBufsNonNull;
        total.nDerivIntegrals   += p.nDerivIntegrals;
        total.nJOps             += p.nJOps;
        total.nKOps             += p.nKOps;
        total.tCompute2         += p.tCompute2;
        total.tScale            += p.tScale;
        total.tContraction      += p.tContraction;
        tMaxC2 = std::max(tMaxC2, p.tCompute2);    tMinC2 = std::min(tMinC2, p.tCompute2);
        tMaxSc = std::max(tMaxSc, p.tScale);       tMinSc = std::min(tMinSc, p.tScale);
        tMaxCt = std::max(tMaxCt, p.tContraction); tMinCt = std::min(tMinCt, p.tContraction);
      }

      using sec = std::chrono::duration<double>;
      double tTotal     = sec(_prof_t_end          - _prof_t_start).count();
      double tSetup     = sec(_prof_t_setup_end    - _prof_t_start).count();
      double tParallel  = sec(_prof_t_parallel_end - _prof_t_setup_end).count();
      double tReduction = sec(_prof_t_end          - _prof_t_parallel_end).count();

      auto pct = [&](double t){ return tTotal > 0.0 ? 100.0*t/tTotal : 0.0; };
      auto rat = [](double a, double b){ return b > 1e-12 ? a/b : 0.0; };
      auto gflops = [](size_t ops, double t){
        return t > 1e-12 ? (2.0 * static_cast<double>(ops) / t) * 1e-9 : 0.0;
      };

      std::cout << "\n========================= directScaffoldGrad profile =========================\n";
      std::cout << std::fixed << std::setprecision(4);
      std::cout << "  Threads: " << nThreads
                << "    MPI ranks: " << mpiSize << " (rank " << mpiRank << ")"
                << "    screen: " << (screen ? "ON" : "OFF");
      std::cout << "\n";
      if (screen)
        std::cout << "    threshSchwarz: " << std::scientific
                  << std::setprecision(2) << tpi.threshSchwarz()
                  << std::fixed << std::setprecision(4);
      std::cout << "\n";
      std::cout << "\n  --- wall time breakdown -----------------------------------------------------\n";
      std::cout << "    total:                " << std::setw(12) << tTotal     << " s\n";
      std::cout << "    setup:                " << std::setw(12) << tSetup     << " s ("
                << std::setw(7) << pct(tSetup)     << "%)\n";
      std::cout << "    parallel region:      " << std::setw(12) << tParallel  << " s ("
                << std::setw(7) << pct(tParallel)  << "%)\n";
      std::cout << "    reduction/symmetrize: " << std::setw(12) << tReduction << " s ("
                << std::setw(7) << pct(tReduction) << "%)\n";

      std::cout << "\n  --- counters (aggregated over threads) --------------------------------------\n";
      std::cout << "    quartets visited:        " << std::setw(18) << total.nQuartetsVisited  << "\n";
      std::cout << "    quartets skipped:        " << std::setw(18) << total.nQuartetsSkipped;
      if (total.nQuartetsVisited + total.nQuartetsSkipped > 0)
        std::cout << "  (skip rate: "
                  << 100.0*total.nQuartetsSkipped/(total.nQuartetsVisited+total.nQuartetsSkipped)
                  << "%)";
      std::cout << "\n";
      std::cout << "    deriv bufs computed:     " << std::setw(18) << total.nDerivBufsNonNull
                << "  (out of " << 12*total.nQuartetsVisited
                << ", " << (total.nQuartetsVisited > 0 ?
                            100.0*total.nDerivBufsNonNull/(12.0*total.nQuartetsVisited) : 0.0)
                << "%)\n";
      std::cout << "    deriv integrals (ERIs):  " << std::setw(18) << total.nDerivIntegrals   << "\n";
      std::cout << "    J ops (FMAs):            " << std::setw(18) << total.nJOps             << "\n";
      std::cout << "    K ops (FMAs):            " << std::setw(18) << total.nKOps             << "\n";
      std::cout << "==============================================================================\n\n";
    }
#endif

  }

}; // namespace ChronusQ

