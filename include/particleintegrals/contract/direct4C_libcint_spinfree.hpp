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


#define _SHZ_SCREEN_4C_LIBCINT

#define _CONTRACTION_

#include <libcint.hpp>

namespace ChronusQ {


  /************************************/
  /* Libcint 4C-direct Implementation */
  /************************************/
  // For DCB Hamiltonian,
  // 12 density matrices upon input stored as
  // LL(MS,MX,MY,MZ), SS(MS,MX,MY,MZ), LS(MS,MX,MY,MZ)
  //
  // 12 contrated matrices upon output stored as
  // LL(MS,MX,MY,MZ), SS(MS,MX,MY,MZ), LS(MS,MX,MY,MZ)
  //
  //
  // Density Matrices are assumed as Hermitian here
  //
 
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directScaffoldLibcintSpinFree(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<MatsT>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {

    if (not matList[0].HER) 
      CErr("Non-Hermitian Density in 4C Contraction (Couloumb + Exchange) is NYI");
    
    DirectTPI<IntsT> &originalERI = *std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->ints_);
    BasisSet& originalBasisSet_ = originalERI.basisSet();
    Molecule& molecule_ = originalERI.molecule();

    if (originalBasisSet_.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet_.groupGeneralContractionBasis();
    const size_t nBasis   = basisSet_.nBasis;

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    DirectTPI<IntsT> &eri = originalERI;
    //DirectTPI<IntsT> eri(basisSet_, molecule_, originalERI.threshSchwarz());
    
    // Determine the number of OpenMP threads
    size_t nThreads  = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1); // Turn off parallelism in LA functions


    // MPI info
    size_t mpiRank   = MPIRank(comm);
    size_t mpiSize   = MPISize(comm);
    size_t mpi_thread_size = mpiSize * nThreads;

#ifdef CQ_ENABLE_MPI
    // Broadcast X to all process
    if( mpiSize > 1 )
      for( auto &C : matList )
        MPIBCast( C.X, nBasis*nBasis, 0, comm );
#endif


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


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::BARE_COULOMB or
         matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::SSSS or
         matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_ipvip1ipvip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLLL or
         matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS) {
        n = int2e_ipvip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUNT) {
        n = int2e_ip1ip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
      if(matList[0].contType == TWOBODY_CONTRACTION_TYPE::GAUGE) {
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

    const size_t nMat     = matList.size();
    const size_t nShell   = basisSet_.nShell;


    enum DIRAC_PAULI_SPINOR_COMP {
      CLLMS,
      XLLMS,
      XLLMX,
      XLLMY,
      XLLMZ,
      CSSMS,
      CSSMX,
      CSSMY,
      CSSMZ,
      XSSMS,
      XSSMX,
      XSSMY,
      XSSMZ,
      XLSMS,
      XLSMX,
      XLSMY,
      XLSMZ,
//    Gaunt
      CLSMS,
      CLSMX,
      CLSMY,
      CLSMZ,
      CSLMS,
      CSLMX,
      CSLMY,
      CSLMZ,
      XSLMS,
      XSLMX,
      XSLMY,
      XSLMZ
    };


    for(auto iMat = 0; iMat < nMat; iMat++)
      memset(matList[iMat].AX,0.,nBasis*nBasis*sizeof(MatsT));

    // Set up scratch space
    std::vector<std::vector<MatsT*>> AXthreads;

    MatsT *AXRaw = nullptr;
    AXRaw = CQMemManager::get().malloc<MatsT>(nThreads*nMat*nBasis*nBasis);    

    for(auto iThread = 0; iThread < nThreads; iThread++) {
      AXthreads.emplace_back();
      for(auto iMat = 0; iMat < nMat; iMat++) 
        AXthreads.back().push_back(AXRaw + iThread*nMat*nBasis*nBasis + iMat*nBasis*nBasis);
    }




    /************************************/
    /*                                  */
    /* Preparation of Schwarz Screening */
    /*                                  */
    /************************************/

    double *SchwarzSSSS = nullptr;

    double *SchwarzERI = nullptr;

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
        size_t mpi_thread_id = mpiRank * nThreads + thread_id;
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
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s12 % mpi_thread_size != mpi_thread_id ) continue;
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

#ifdef CQ_ENABLE_MPI
      MPIAllReduceInPlace(SchwarzERI, nShell*nShell, comm);
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durERIchwarz = tock(topERIchwarz);
//      std::cout << "ERI Schwarz took " <<  durERIchwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
  
    };
 

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
        size_t mpi_thread_id = mpiRank * nThreads + thread_id;
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
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s12 % mpi_thread_size != mpi_thread_id ) continue;
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

#ifdef CQ_ENABLE_MPI
      MPIAllReduceInPlace(SchwarzSSSS, nShell*nShell, comm);
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto durSSSSSchwarz = tock(topSSSSSchwarz);
//      std::cout << "∇A∇B∇C∇D Schwarz took " <<  durSSSSSchwarz << " s\n"; 
  
      //std::cout << std::endl;
#endif
  
    };
    /****************************************/
    /*                                      */
    /* End Preparation of Schwarz Screening */
    /*                                      */
    /****************************************/











    /***********************************/
    /*                                 */
    /* Start of Dirac-Coulomb C(2)     */
    /*                                 */
    /***********************************/

    if( matList[0].contType == TWOBODY_CONTRACTION_TYPE::LLSS ) {

 
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectLL = tick();
#endif

#ifdef _SHZ_SCREEN_4C_LIBCINT
      // Compute shell block norms (∞-norm) of matList.X
      // CLLMS, CSSMS, CSSMX, CSSMY, CSSMZ Densitry matrices
      int mMat = 9;
      double *ShBlkNorms_raw = CQMemManager::get().malloc<double>(mMat*nShell*nShell);
      std::vector<double*> ShBlkNorms;
  
      auto iOff = 0;
      ShellBlockNorm(basisSet_.shells,matList[CLLMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
 
      // Get the max over all the matricies for the shell block ∞-norms
      for(auto k = 0; k < nShell*nShell; k++) {
  
        double mx = std::abs(ShBlkNorms[0][k]);
        for(auto iMat = 0; iMat < mMat; iMat++)
          mx = std::max(mx,std::abs(ShBlkNorms[iMat][k]));
        ShBlkNorms[0][k] = mx;
   
      }
#endif

      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(2*NB4*nThreads);
      memset(AXRaw,0,nThreads*nMat*nBasis*nBasis*sizeof(MatsT));
      // Keeping track of number of integrals skipped
      std::vector<size_t> nSkipLL(nThreads,0);

#ifdef _THREAD_TIMING_
      std::vector<double> durThread(nThreads,0);
#endif

      #pragma omp parallel
      {
#ifdef _THREAD_TIMING_
        auto LLSSBegin = tick();
#endif

        dcomplex iscale = dcomplex(0.0, 1.0);

        size_t thread_id = GetThreadID();
        int mpi_thread_id = mpiRank * nThreads + thread_id;
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *ERIBuffAB   = &ERIBuffer[thread_id*NB4];
        double *ERIBuffCD   = &ERIBuffer[nThreads*NB4 + thread_id*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax12 = ShBlkNorms[0][s1 + s2*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 < nShell; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax123 = std::max(ShBlkNorms[0][s1 + s3*nShell], ShBlkNorms[0][s2 + s3*nShell]);
          shMax123 = std::max(shMax123,shMax12);
#endif
   
        for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4

          // Round Robbin work distribution
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s1234 % mpi_thread_size != mpi_thread_id ) continue;
          #endif

#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax = std::max(ShBlkNorms[0][s1 + s4*nShell],
                                  std::max(ShBlkNorms[0][s2 + s4*nShell],
                                           ShBlkNorms[0][s3 + s4*nShell]));

          shMax = std::max(shMax,shMax123);

//          double shMax = std::max(ShBlkNorms[0][s1 + s2*nShell],ShBlkNorms[0][s3 + s4*nShell]);

          if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzERI[s3+s4*nShell]) <
             eri.threshSchwarz()) { nSkipLL[thread_id]++; continue; }
#endif

          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>3) {nSkipLL[thread_id]++; continue;}
          }
 

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>2) {nSkipLL[thread_id]++; continue;}
          }
 
          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) 
          if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) )
            {nSkipLL[thread_id]++; continue;}
 
          auto nQuad = n1*n2*n3*n4;
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);
  
          //if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
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
            auto MNKL = m + nNBkNB2lNB3;
            //auto KLMN = k + l*NB + m*NB2 + n*NB3;
            auto KLMN = m*NB2 + klNBnNB3;
  
  
            // ∇A∙∇B(mn|kl) followed by ∇Ax∇B(mn|kl) X, Y, and Z
            // (mn|kl)
            ERIBuffAB[       MNKL] =  dAdotdB;
  
            // ∇C∙∇D(kl|mn) followed by ∇Cx∇D(kl|mn) X, Y, and Z
            // (kl|mn)
            ERIBuffCD[       KLMN] =  dAdotdB;

          }
          }
          } 
          } // ∇A∇B integral preparation loop
  


#ifdef _CONTRACTION_ // Contraction

          auto ADCLLMS = AX_loc[CLLMS];
          auto ADCSSMS = AX_loc[CSSMS];
          auto ADCSSMX = AX_loc[CSSMX];
          auto ADCSSMY = AX_loc[CSSMY];
          auto ADCSSMZ = AX_loc[CSSMZ];

          auto DCLLMS  = matList[CLLMS].X;
          auto DCSSMS  = matList[CSSMS].X;
          auto DCSSMX  = matList[CSSMX].X;
          auto DCSSMY  = matList[CSSMY].X;
          auto DCSSMZ  = matList[CSSMZ].X;

          auto ADXLSMS  = AX_loc[XLSMS];
          auto ADXLSMX  = AX_loc[XLSMX];
          auto ADXLSMY  = AX_loc[XLSMY];
          auto ADXLSMZ  = AX_loc[XLSMZ];

          auto DXLSMS = matList[XLSMS].X;
          auto DXLSMX = matList[XLSMX].X;
          auto DXLSMY = matList[XLSMY].X;
          auto DXLSMZ = matList[XLSMZ].X;

          for(m = 0ul,            bf1 = bf1_s; m <                  n1; ++m, bf1++) {
            auto mNB2 = m*NB2;
            auto bf1nB = bf1*nBasis;

          for(n =   maxShellSize, bf2 = bf2_s; n <   maxShellSize + n2; ++n, bf2++) {
            auto mNB2nNB3 = mNB2 + n*NB3;
            auto mnNB     = m + n*NB;
            auto bf21 = bf2 + bf1nB;
            auto bf12 = bf1 + bf2*nBasis;

            auto mn1 = m + n*NB;
            auto mn2 = mNB2 + n*NB3;
            auto bf2nB = bf2*nBasis;

          for(k = 2*maxShellSize, bf3 = bf3_s; k < 2*maxShellSize + n3; ++k, bf3++) {
            auto mnNBkNB2  = mnNB + k*NB2;
            auto kmNB2nNB3 = k + mNB2nNB3;
            auto bf3nB = bf3*nBasis;

            auto kmn1 = mn1 + k*NB2;
            auto kmn2 = k + mn2;

          for(l = 3*maxShellSize, bf4 = bf4_s; l < 3*maxShellSize + n4; ++l, bf4++) {

            //auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto MNKL = mnNBkNB2 + l*NB3;
            //auto KLMN = k + l*NB + m*NB2 + n*NB3;
            auto KLMN = kmNB2nNB3 + l*NB;
  
            auto bf43 = bf4 + bf3nB;
            auto bf34 = bf3 + bf4*nBasis;

            auto bf4nB= bf4*nBasis;  

            auto bf14 = bf1 + bf4nB;
            auto bf41 = bf4 + bf1nB;
            auto bf24 = bf2 + bf4nB;
            auto bf42 = bf4 + bf2nB;
            auto bf23 = bf2 + bf3nB;
            auto bf32 = bf3 + bf2nB;
            auto bf13 = bf1 + bf3nB;
            auto bf31 = bf3 + bf1nB;
  
            auto DotPrdMNKL = MNKL;
            auto DotPrdKLMN = KLMN;
  
            // Deneneracy factor for s1,s2 pair
            double s12_deg = (bf1_s == bf2_s) ? 1.0 : 2.0;
            double s34_deg = (bf3_s == bf4_s) ? 1.0 : 2.0;

            /*++++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb (LL|LL) Contraction */
            /*++++++++++++++++++++++++++++++++++++++++++++*/

           // The LLLL block is all Coulomb type
  
            // KLMN
            // Equation (C6) in Xiaosong's note
            if(bf3 >= bf4) 
              ADCLLMS[bf34] +=  s12_deg*ERIBuffCD[DotPrdKLMN]*DCSSMS[bf21].real();
           
            /*------------------------------------------*/
            /* End of Dirac-Coulomb (LL|LL) Contraction */
            /*------------------------------------------*/
  
  
  
            /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb C(2)-(SS|SS) Contraction */
            /*+++++++++++++++++++++++++++++++++++++++++++++++++*/

            // The SSSS block is all Coulomb type
  
            // MNKL 
            // Equations (C13) and (C14) in Xiaosong's note
            // Complex i for the X, Y, Z components is multiplied at the end
            if (bf1 >= bf2) 
              ADCSSMS[bf12] += s34_deg*ERIBuffAB[DotPrdMNKL]*DCLLMS[bf43].real();
  
            /*-----------------------------------------------*/
            /* End of Dirac-Coulomb C(2)-(SS|SS) Contraction */
            /*-----------------------------------------------*/


            /*++++++++++++++++++++++++++++++++++++++++++*/
            /* Start of Dirac-Coulomb (LL|SS) / (SS|LL) */
            /*++++++++++++++++++++++++++++++++++++++++++*/

            // The LLSS block is all exchange type
#if 1
            //KLMN 3412
            // Equation (C7-C10) in Xiaosong's note
            ADXLSMS[bf32]+= -ERIBuffCD[DotPrdKLMN]*DXLSMS[bf41];
            ADXLSMX[bf32]+= -ERIBuffCD[DotPrdKLMN]*DXLSMX[bf41];
            ADXLSMY[bf32]+= -ERIBuffCD[DotPrdKLMN]*DXLSMY[bf41];
            ADXLSMZ[bf32]+= -ERIBuffCD[DotPrdKLMN]*DXLSMZ[bf41];
  
            // KLMN 3421
            // Since the derivative is on 21, the sign of the cross product will
            // change when the indices of 21 swap.
            if(bf2_s!=bf1_s) {
              ADXLSMS[bf31]+= -ERIBuffCD[DotPrdKLMN]*DXLSMS[bf42];
              ADXLSMX[bf31]+= -ERIBuffCD[DotPrdKLMN]*DXLSMX[bf42];
              ADXLSMY[bf31]+= -ERIBuffCD[DotPrdKLMN]*DXLSMY[bf42];
              ADXLSMZ[bf31]+= -ERIBuffCD[DotPrdKLMN]*DXLSMZ[bf42];
            }
  
            //LKMN 4312
            if(bf3_s!=bf4_s){
              ADXLSMS[bf42]+= -ERIBuffCD[DotPrdKLMN]*DXLSMS[bf31];
              ADXLSMX[bf42]+= -ERIBuffCD[DotPrdKLMN]*DXLSMX[bf31];
              ADXLSMY[bf42]+= -ERIBuffCD[DotPrdKLMN]*DXLSMY[bf31];
              ADXLSMZ[bf42]+= -ERIBuffCD[DotPrdKLMN]*DXLSMZ[bf31];
  
              if(bf1_s!=bf2_s) {
                //LKNM 4321
                ADXLSMS[bf41]+= -ERIBuffCD[DotPrdKLMN]*DXLSMS[bf32];
                ADXLSMX[bf41]+= -ERIBuffCD[DotPrdKLMN]*DXLSMX[bf32];
                ADXLSMY[bf41]+= -ERIBuffCD[DotPrdKLMN]*DXLSMY[bf32];
                ADXLSMZ[bf41]+= -ERIBuffCD[DotPrdKLMN]*DXLSMZ[bf32];
              }
            }
  
            /*------------------------------------------*/
            /*   End of Dirac-Coulomb (LL|SS) / (SS|LL) */
            /*------------------------------------------*/
#endif

          };
          };
          };  
          }; // contraction loop
  
#endif // Contraction
  
  
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1

#ifdef _THREAD_TIMING_
        durThread[thread_id] = tock(LLSSBegin);
#endif
  
      } // OpenMP context
 
      dcomplex iS = dcomplex(0.0, 1.0);

      for( auto iTh  = 1; iTh < nThreads; iTh++) {

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLLMS],nBasis,MatsT(1.0),
            AXthreads[0][CLLMS],nBasis,AXthreads[0][CLLMS],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMS],nBasis,MatsT(1.0),
            AXthreads[0][CSSMS],nBasis,AXthreads[0][CSSMS],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMX],nBasis,MatsT(1.0),
            AXthreads[0][CSSMX],nBasis,AXthreads[0][CSSMX],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMY],nBasis,MatsT(1.0),
            AXthreads[0][CSSMY],nBasis,AXthreads[0][CSSMY],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMZ],nBasis,MatsT(1.0),
            AXthreads[0][CSSMZ],nBasis,AXthreads[0][CSSMZ],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMS],nBasis,MatsT(1.0),
            AXthreads[0][XLSMS],nBasis,AXthreads[0][XLSMS],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMX],nBasis,MatsT(1.0),
            AXthreads[0][XLSMX],nBasis,AXthreads[0][XLSMX],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMY],nBasis,MatsT(1.0),
            AXthreads[0][XLSMY],nBasis,AXthreads[0][XLSMY],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMZ],nBasis,MatsT(1.0),
            AXthreads[0][XLSMZ],nBasis,AXthreads[0][XLSMZ],nBasis);

      };

      #ifdef CQ_ENABLE_MPI
      // Combine all G[X] contributions onto all processes
      if( mpiSize > 1 ) {
        MatsT* mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);

        std::vector<DIRAC_PAULI_SPINOR_COMP> comps{CLLMS, CSSMS, CSSMX, CSSMY, CSSMZ,
                                                   XLSMS, XLSMX, XLSMY, XLSMZ};
        for (auto comp : comps)
          MPIAllReduceInPlace( AXthreads[0][comp], nBasis*nBasis, comm, mpiScr );

        CQMemManager::get().free(mpiScr);

      }
      #endif

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLLMS],nBasis,MatsT(1.0),
             matList[CLLMS].AX,nBasis,matList[CLLMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CSSMS],nBasis,MatsT(1.0),
             matList[CSSMS].AX,nBasis,matList[CSSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis, iS, AXthreads[0][CSSMX],nBasis,MatsT(1.0),
             matList[CSSMX].AX,nBasis,matList[CSSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis, iS, AXthreads[0][CSSMY],nBasis,MatsT(1.0),
             matList[CSSMY].AX,nBasis,matList[CSSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis, iS, AXthreads[0][CSSMZ],nBasis,MatsT(1.0),
             matList[CSSMZ].AX,nBasis,matList[CSSMZ].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMS],nBasis,MatsT(1.0),
             matList[XLSMS].AX,nBasis,matList[XLSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMX],nBasis,MatsT(1.0),
             matList[XLSMX].AX,nBasis,matList[XLSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMY],nBasis,MatsT(1.0),
             matList[XLSMY].AX,nBasis,matList[XLSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMZ],nBasis,MatsT(1.0),
             matList[XLSMZ].AX,nBasis,matList[XLSMZ].AX,nBasis);

      // Take care of the Hermitian symmetry in the LL and SS blocks
      auto ADCLLMS = matList[CLLMS].AX;
      auto ADCSSMS = matList[CSSMS].AX;
      auto ADCSSMX = matList[CSSMX].AX;
      auto ADCSSMY = matList[CSSMY].AX;
      auto ADCSSMZ = matList[CSSMZ].AX;

      for( auto i = 0; i < nBasis; i++ ) 
      for( auto j = 0; j < i; j++ ) {
        ADCLLMS[j + i*nBasis] = std::conj(ADCLLMS[i + j*nBasis]);
        ADCSSMS[j + i*nBasis] = std::conj(ADCSSMS[i + j*nBasis]);
        ADCSSMX[j + i*nBasis] = std::conj(ADCSSMX[i + j*nBasis]);
        ADCSSMY[j + i*nBasis] = std::conj(ADCSSMY[i + j*nBasis]);
        ADCSSMZ[j + i*nBasis] = std::conj(ADCSSMZ[i + j*nBasis]);
      }

      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
#ifdef _SHZ_SCREEN_4C_LIBCINT
      if(ShBlkNorms_raw!=nullptr) CQMemManager::get().free(ShBlkNorms_raw);
#endif

#ifdef _THREAD_TIMING_
      std::cout << "Dirac-Coulomb-C(2) Libcint time on every thread:" << std::endl;
      for (size_t i = 0; i < nThreads; i++)
        std::cout << i << "\t:" << durThread[i] << std::endl;
#endif

#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipLL = std::accumulate(nSkipLL.begin(),nSkipLL.end(),size_t(0));
#ifdef CQ_ENABLE_MPI
      if (mpiSize > 1) {
        std::cout << "Dirac-Coulomb-C(2) Screened "
                  << nIntSkipLL << " on Rank " << mpiRank <<  std::endl;
        nIntSkipLL = MPIAllReduce( nIntSkipLL, comm );
      }
#endif
      std::cout << "Dirac-Coulomb-C(2) Screened " << nIntSkipLL << std::endl;

      auto durDirectLL = tock(topDirectLL);
      std::cout << "Spin-Free Dirac-Coulomb-C(2) AO Direct Contraction took " <<  durDirectLL << " s\n"; 

      std::cout << std::endl;
#endif
    

    } // if(LLSS)
   
    /********************************/
    /*                              */
    /*   End of Dirac-Coulomb C(2)  */
    /*                              */
    /********************************/




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
#endif

#ifdef _SHZ_SCREEN_4C_LIBCINT
      // Compute shell block norms (∞-norm) of matList.X
      // CSSMS, CSSMX, CSSMY, CSSMZ Densitry matrices
      int mMat = 4;
      double *ShBlkNorms_raw = CQMemManager::get().malloc<double>(mMat*nShell*nShell);
      std::vector<double*> ShBlkNorms;
  
      auto iOff = 0;
      ShellBlockNorm(basisSet_.shells,matList[CSSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[CSSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      // Get the max over all the matricies for the shell block ∞-norms
      for(auto k = 0; k < nShell*nShell; k++) {
  
        double mx = std::abs(ShBlkNorms[0][k]);
        for(auto iMat = 0; iMat < mMat; iMat++)
          mx = std::max(mx,std::abs(ShBlkNorms[iMat][k]));
        ShBlkNorms[0][k] = mx;
   
      }
#endif



      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);
      ERIBuffer = CQMemManager::get().malloc<double>(NB4*nThreads);
  
      memset(AXRaw,0,nThreads*nMat*nBasis*nBasis*sizeof(MatsT));
  
      // Keeping track of number of integrals skipped
      std::vector<size_t> nSkipSSSS(nThreads,0);

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

#ifdef _THREAD_TIMING_
      std::vector<double> durThread(nThreads,0);
#endif
  
      #pragma omp parallel
      {
#ifdef _THREAD_TIMING_
        auto SSSSBegin = tick();
#endif
  
        dcomplex iS = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
        int mpi_thread_id = mpiRank * nThreads + thread_id;
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *ERIBuffABCD = &ERIBuffer[thread_id*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet_.nShell; 
            bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax12 = ShBlkNorms[0][s1 + s2*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax123 = std::max(ShBlkNorms[0][s1 + s3*nShell], ShBlkNorms[0][s2 + s3*nShell]);
          shMax123 = std::max(shMax123,shMax12);
#endif

        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s1234 % mpi_thread_size != mpi_thread_id ) continue;
          #endif

#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax = std::max(ShBlkNorms[0][s1 + s4*nShell],
                                  std::max(ShBlkNorms[0][s2 + s4*nShell],
                                           ShBlkNorms[0][s3 + s4*nShell]));

          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzSSSS[s1+s2*nShell]*SchwarzSSSS[s3 + s4*nShell]) <
             eri.threshSchwarz()) { nSkipSSSS[thread_id]++; continue; }
#endif

 
          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>3) {nSkipSSSS[thread_id]++; continue;}
          }
 

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>2) {nSkipSSSS[thread_id]++; continue;}
          }
 
          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) 
          if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) )
            {nSkipSSSS[thread_id]++; continue;}

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
  
  
  
  
  
  
  
#ifdef _CONTRACTION_ // Contraction
  
          auto DCSSMS = matList[CSSMS].X;
          auto DCSSMX = matList[CSSMX].X;
          auto DCSSMY = matList[CSSMY].X;
          auto DCSSMZ = matList[CSSMZ].X;
  
          auto ADCSSMS  = AX_loc[CSSMS];
          auto ADCSSMX  = AX_loc[CSSMX];
          auto ADCSSMY  = AX_loc[CSSMY];
          auto ADCSSMZ  = AX_loc[CSSMZ];
  
          auto DXSSMS = matList[XSSMS].X;
          auto DXSSMX = matList[XSSMX].X;
          auto DXSSMY = matList[XSSMY].X;
          auto DXSSMZ = matList[XSSMZ].X;
  
          auto ADXSSMS  = AX_loc[XSSMS];
          auto ADXSSMX  = AX_loc[XSSMX];
          auto ADXSSMY  = AX_loc[XSSMY];
          auto ADXSSMZ  = AX_loc[XSSMZ];
  
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
            auto MNKLdAdotdBdCdotdD = ERIBuffABCD[MNKL];
            auto KLMNdAdotdBdCdotdD = ERIBuffABCD[KLMN];
  
            //
            // COULOMB
            // MNKL
            if(bf1 >= bf2) {
              /* Hermitian Density Matrices */
              double s34_deg = (bf3_s == bf4_s) ? 1.0 : 2.0;
              /* Equation (C19-C22) in Xiaosong's note */
              ADCSSMS[bf12] += DCSSMS[bf43].real()*MNKLdAdotdBdCdotdD*s34_deg;
            }
  
            /* KLMN */
            if((bf1_s!=bf3_s or bf2_s!=bf4_s) and bf3 >= bf4){
              double s12_deg = (bf1_s == bf2_s) ? 1.0 : 2.0;
              ADCSSMS[bf34] += DCSSMS[bf21].real()*KLMNdAdotdBdCdotdD*s12_deg;
            }
  
   
            /* EXCHANGE */
            // MNKL
            if(bf1 >= bf4) {
              ADXSSMS[bf14] += -DXSSMS[bf23]*MNKLdAdotdBdCdotdD;
              ADXSSMZ[bf14] += -DXSSMZ[bf23]*MNKLdAdotdBdCdotdD;
              ADXSSMX[bf14] += -DXSSMX[bf23]*MNKLdAdotdBdCdotdD;
              ADXSSMY[bf14] += -DXSSMY[bf23]*MNKLdAdotdBdCdotdD;
            }
  
            //MNLK
            if(bf3_s!=bf4_s and bf1 >= bf3) {
              ADXSSMS[bf13] += -DXSSMS[bf24]*MNKLdAdotdBdCdotdD;
              ADXSSMZ[bf13] += -DXSSMZ[bf24]*MNKLdAdotdBdCdotdD;
              ADXSSMX[bf13] += -DXSSMX[bf24]*MNKLdAdotdBdCdotdD;
              ADXSSMY[bf13] += -DXSSMY[bf24]*MNKLdAdotdBdCdotdD;
            }
  
            //NMKL
            if(bf1_s!=bf2_s){
              if(bf2 >= bf4) {
                ADXSSMS[bf24] += -DXSSMS[bf13]*MNKLdAdotdBdCdotdD;
                ADXSSMZ[bf24] += -DXSSMZ[bf13]*MNKLdAdotdBdCdotdD;
                ADXSSMX[bf24] += -DXSSMX[bf13]*MNKLdAdotdBdCdotdD;
                ADXSSMY[bf24] += -DXSSMY[bf13]*MNKLdAdotdBdCdotdD;
              }
    
  
              if(bf3_s!=bf4_s and bf2 >= bf3) {
                ADXSSMS[bf23] += -DXSSMS[bf14]*MNKLdAdotdBdCdotdD;
                ADXSSMZ[bf23] += -DXSSMZ[bf14]*MNKLdAdotdBdCdotdD;
                ADXSSMX[bf23] += -DXSSMX[bf14]*MNKLdAdotdBdCdotdD;
                ADXSSMY[bf23] += -DXSSMY[bf14]*MNKLdAdotdBdCdotdD;
              }
            }
  
  
            if(bf1_s!=bf3_s or bf2_s!=bf4_s){
              if(bf3 >= bf2 ) {
                ADXSSMS[bf32] += -DXSSMS[bf41]*KLMNdAdotdBdCdotdD;
                ADXSSMZ[bf32] += -DXSSMZ[bf41]*KLMNdAdotdBdCdotdD;
                ADXSSMX[bf32] += -DXSSMX[bf41]*KLMNdAdotdBdCdotdD;
                ADXSSMY[bf32] += -DXSSMY[bf41]*KLMNdAdotdBdCdotdD;
              }
   
              
              if(bf1_s!=bf2_s and bf3>=bf1) {
  
                ADXSSMS[bf31] += -DXSSMS[bf42]*KLMNdAdotdBdCdotdD;
                ADXSSMZ[bf31] += -DXSSMZ[bf42]*KLMNdAdotdBdCdotdD;
                ADXSSMX[bf31] += -DXSSMX[bf42]*KLMNdAdotdBdCdotdD;
                ADXSSMY[bf31] += -DXSSMY[bf42]*KLMNdAdotdBdCdotdD;
              }
    
              //NMKL
              if(bf3_s!=bf4_s and bf4>=bf2){
                ADXSSMS[bf42] += -DXSSMS[bf31]*KLMNdAdotdBdCdotdD;
                ADXSSMZ[bf42] += -DXSSMZ[bf31]*KLMNdAdotdBdCdotdD;
                ADXSSMX[bf42] += -DXSSMX[bf31]*KLMNdAdotdBdCdotdD;
                ADXSSMY[bf42] += -DXSSMY[bf31]*KLMNdAdotdBdCdotdD;
    
                if(bf1_s!=bf2_s and bf4>=bf1) {
                  ADXSSMS[bf41] += -DXSSMS[bf32]*KLMNdAdotdBdCdotdD;
                  ADXSSMZ[bf41] += -DXSSMZ[bf32]*KLMNdAdotdBdCdotdD;
                  ADXSSMX[bf41] += -DXSSMX[bf32]*KLMNdAdotdBdCdotdD;
                  ADXSSMY[bf41] += -DXSSMY[bf32]*KLMNdAdotdBdCdotdD;
                }
              }
 
            }

	  }
	  }
	  }
          } // contraction loop
  
#endif // Contraction
  
  
      }; // loop s4
      }; // loop s3
      }; // loop s2
      }; // loop s1

#ifdef _THREAD_TIMING_
durThread[thread_id] = tock(SSSSBegin);
#endif
  
      }; // OpenMP context
  
  
      dcomplex iS = dcomplex(0.0, 1.0);

      for( auto iTh  = 1; iTh < nThreads; iTh++) {

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMS],nBasis,MatsT(1.0),
            AXthreads[0][CSSMS],nBasis,AXthreads[0][CSSMS],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMX],nBasis,MatsT(1.0),
            AXthreads[0][CSSMX],nBasis,AXthreads[0][CSSMX],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMY],nBasis,MatsT(1.0),
            AXthreads[0][CSSMY],nBasis,AXthreads[0][CSSMY],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CSSMZ],nBasis,MatsT(1.0),
            AXthreads[0][CSSMZ],nBasis,AXthreads[0][CSSMZ],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMS],nBasis,MatsT(1.0),
            AXthreads[0][XSSMS],nBasis,AXthreads[0][XSSMS],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMX],nBasis,MatsT(1.0),
            AXthreads[0][XSSMX],nBasis,AXthreads[0][XSSMX],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMY],nBasis,MatsT(1.0),
            AXthreads[0][XSSMY],nBasis,AXthreads[0][XSSMY],nBasis);

        MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMZ],nBasis,MatsT(1.0),
            AXthreads[0][XSSMZ],nBasis,AXthreads[0][XSSMZ],nBasis);

      };

      #ifdef CQ_ENABLE_MPI
      // Combine all G[X] contributions onto all processes
      if( mpiSize > 1 ) {
        MatsT* mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);

        std::vector<DIRAC_PAULI_SPINOR_COMP> comps{CSSMS, CSSMX, CSSMY, CSSMZ,
                                                   XSSMS, XSSMX, XSSMY, XSSMZ};
        for (auto comp : comps)
          MPIAllReduceInPlace( AXthreads[0][comp], nBasis*nBasis, comm, mpiScr );

        CQMemManager::get().free(mpiScr);

      }
      #endif

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CSSMS],nBasis,MatsT(1.0),
             matList[CSSMS].AX,nBasis,matList[CSSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,iS,AXthreads[0][CSSMX],nBasis,MatsT(1.0),
             matList[CSSMX].AX,nBasis,matList[CSSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,iS,AXthreads[0][CSSMY],nBasis,MatsT(1.0),
             matList[CSSMY].AX,nBasis,matList[CSSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,iS,AXthreads[0][CSSMZ],nBasis,MatsT(1.0),
             matList[CSSMZ].AX,nBasis,matList[CSSMZ].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMS],nBasis,MatsT(1.0),
             matList[XSSMS].AX,nBasis,matList[XSSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMX],nBasis,MatsT(1.0),
             matList[XSSMX].AX,nBasis,matList[XSSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMY],nBasis,MatsT(1.0),
             matList[XSSMY].AX,nBasis,matList[XSSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMZ],nBasis,MatsT(1.0),
             matList[XSSMZ].AX,nBasis,matList[XSSMZ].AX,nBasis);

      
      // Take care of the Hermitian symmetry in the LL and SS blocks
      auto ADCSSMS = matList[CSSMS].AX;
      auto ADCSSMX = matList[CSSMX].AX;
      auto ADCSSMY = matList[CSSMY].AX;
      auto ADCSSMZ = matList[CSSMZ].AX;
      auto ADXSSMS = matList[XSSMS].AX;
      auto ADXSSMX = matList[XSSMX].AX;
      auto ADXSSMY = matList[XSSMY].AX;
      auto ADXSSMZ = matList[XSSMZ].AX;

      for( auto i = 0; i < nBasis; i++ )
      for( auto j = 0; j < i; j++ ) {
        ADCSSMS[j + i*nBasis] = std::conj(ADCSSMS[i + j*nBasis]);
        ADCSSMX[j + i*nBasis] = std::conj(ADCSSMX[i + j*nBasis]);
        ADCSSMY[j + i*nBasis] = std::conj(ADCSSMY[i + j*nBasis]);
        ADCSSMZ[j + i*nBasis] = std::conj(ADCSSMZ[i + j*nBasis]);
        ADXSSMS[j + i*nBasis] = std::conj(ADXSSMS[i + j*nBasis]);
        ADXSSMX[j + i*nBasis] = std::conj(ADXSSMX[i + j*nBasis]);
        ADXSSMY[j + i*nBasis] = std::conj(ADXSSMY[i + j*nBasis]);
        ADXSSMZ[j + i*nBasis] = std::conj(ADXSSMZ[i + j*nBasis]);
      }

      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
#ifdef _SHZ_SCREEN_4C_LIBCINT
      if(ShBlkNorms_raw!=nullptr) CQMemManager::get().free(ShBlkNorms_raw);
#endif

#ifdef _THREAD_TIMING_
      std::cout << "Dirac-Coulomb-SSSS Libcint time on every thread:" << std::endl;
      for (size_t i = 0; i < nThreads; i++)
        std::cout << i << "\t:" << durThread[i] << std::endl;
#endif

#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipSSSS = std::accumulate(nSkipSSSS.begin(),nSkipSSSS.end(),size_t(0));
#ifdef CQ_ENABLE_MPI
      if (mpiSize > 1) {
        std::cout << "Dirac-Coulomb-SSSS Screened "
                  << nIntSkipSSSS << " on Rank " << mpiRank <<  std::endl;
        nIntSkipSSSS = MPIAllReduce( nIntSkipSSSS, comm );
      }
#endif
      std::cout << "Dirac-Coulomb-SSSS Screened " << nIntSkipSSSS << std::endl;
  
      auto durDirectSSSS = tock(topDirectSSSS);
      std::cout << "Spin-Free Dirac-Coulomb-SSSS AO Direct Contraction took " <<  durDirectSSSS << " s\n"; 
  
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
      auto NB4_13 =13*NB4;  
      auto NB4_14 =14*NB4;  
      auto NB4_15 =15*NB4;  
      auto NB4_16 =16*NB4;  
      auto NB4_17 =17*NB4;  
      auto NB4_18 =18*NB4;  

                              //ERI0 :   ∇B∙∇C(mn|kl)
      auto crossx    =NB4    ;//ERI1 :   ∇Bx∇C(mn|kl)-X
      auto crossy    =NB4_2  ;//ERI2 :   ∇Bx∇C(mn|kl)-Y
      auto crossz    =NB4_3  ;//ERI3 :   ∇Bx∇C(mn|kl)-Z
      auto xypyx     =NB4_4  ;//ERI4 :   ∇B_x∇C_y(mn|kl) + ∇B_y∇C_x(mn|kl)
      auto yx        =NB4_5  ;//ERI5 :   ∇B_y∇C_x(mn|kl)
      auto xzpzx     =NB4_6  ;//ERI6 :   ∇B_x∇C_z(mn|kl) + ∇B_z∇C_x(mn|kl)
      auto zx        =NB4_7  ;//ERI7 :   ∇B_z∇C_x(mn|kl)
      auto yzpzy     =NB4_8  ;//ERI8 :   ∇B_y∇C_z(mn|kl) + ∇B_z∇C_y(mn|kl)
      auto zy        =NB4_9  ;//ERI9 :   ∇B_z∇C_y(mn|kl)
      auto mxxmyypzz =NB4_10 ;//ERI10: - ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) + ∇B_z∇C_z(mn|kl)
      auto pxxmyymzz =NB4_11 ;//ERI11:   ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
      auto mxxpyymzz =NB4_12 ;//ERI12: - ∇B_x∇C_x(mn|kl) + ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
      auto xx        =NB4_13 ;//ERI13:   ∇B_x∇C_x(mn|kl)
      auto xy        =NB4_14 ;//ERI14:   ∇B_x∇C_y(mn|kl)
      auto xz        =NB4_15 ;//ERI15:   ∇B_x∇C_z(mn|kl)
      auto yy        =NB4_16 ;//ERI16:   ∇B_y∇C_y(mn|kl)
      auto yz        =NB4_17 ;//ERI17:   ∇B_y∇C_z(mn|kl)
      auto zz        =NB4_18 ;//ERI18:   ∇B_z∇C_z(mn|kl)
 
#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectGaunt = tick();
#endif

      nERI = 1;
      buffAll = CQMemManager::get().malloc<double>(nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);


#ifdef _SHZ_SCREEN_4C_LIBCINT
  
      double *SchwarzGaunt = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzGaunt,0,nShell*nShell*sizeof(double));
  
      #pragma omp parallel
      {
  
        double C1 = 1./(2*SpeedOfLight);
        size_t thread_id = GetThreadID();
        size_t mpi_thread_id = mpiRank * nThreads + thread_id;
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
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s12 % mpi_thread_size != mpi_thread_id ) continue;
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

#ifdef CQ_ENABLE_MPI
      MPIAllReduceInPlace(SchwarzGaunt, nShell*nShell, comm);
#endif

  

      // Compute shell block norms (∞-norm) of matList.X
      int mMat = 16;
      double *ShBlkNorms_raw = CQMemManager::get().malloc<double>(mMat*nShell*nShell);
      std::vector<double*> ShBlkNorms;
  
      auto iOff = 0;
      ShellBlockNorm(basisSet_.shells,matList[XLSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
 
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      // Get the max over all the matricies for the shell block ∞-norms
      for(auto k = 0; k < nShell*nShell; k++) {
  
        double mx = std::abs(ShBlkNorms[0][k]);
        for(auto iMat = 0; iMat < mMat; iMat++)
          mx = std::max(mx,std::abs(ShBlkNorms[iMat][k]));
        ShBlkNorms[0][k] = mx;
   
      }
#endif

      int nSave = 1;
      ERIBuffer = CQMemManager::get().malloc<double>(nSave*NB4*nThreads);
      memset(AXRaw,0,nThreads*nMat*nBasis*nBasis*sizeof(MatsT));
      std::vector<size_t> nSkipGaunt(nThreads,0);

#ifdef _THREAD_TIMING_
      std::vector<double> durThread(nThreads,0);
#endif
  
      #pragma omp parallel
      {
#ifdef _THREAD_TIMING_
        auto gauntBegin = tick();
#endif
  
        dcomplex iS = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
        int mpi_thread_id = mpiRank * nThreads + thread_id;
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *BC   = &ERIBuffer[thread_id*nSave*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
  
        int shls[4];
        double *buff = &buffAll[nERI*buffN4*thread_id];
        double *cache = cacheAll+cache_size*thread_id;
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax12 = ShBlkNorms[0][s2 + s1*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s2 ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s2 == s3) ? s1 : nShell-1; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax123 = std::max(ShBlkNorms[0][s1 + s3*nShell], ShBlkNorms[0][s2 + s3*nShell]);
          shMax123 = std::max(shMax123,shMax12);
#endif
   
  
        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s1234 % mpi_thread_size != mpi_thread_id ) continue;
          #endif
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
  
          double shMax = std::max(ShBlkNorms[0][s1 + s4*nShell],
                                  std::max(ShBlkNorms[0][s2 + s4*nShell],
                                           ShBlkNorms[0][s3 + s4*nShell]));
  
          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzGaunt[s3+s4*nShell]*SchwarzGaunt[s2 + s1*nShell]) <
             eri.threshSchwarz()) { nSkipGaunt[thread_id]++; continue; }
#endif

          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>3) {nSkipGaunt[thread_id]++; continue;}
          }
 

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>2) {nSkipGaunt[thread_id]++; continue;}
          }
  
          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) 
          if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) )
            {nSkipGaunt[thread_id]++; continue;}
  
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
            BC[MNKL] = dAdotdC;    
            // (lk|nm)
            BC[LKNM] = dAdotdC;

          } // ∇C∙∇D integral preparation loop
  
  
  
#ifdef _CONTRACTION_ // Contraction

          auto ADXLLMS  = AX_loc[XLLMS];
          auto ADXLLMX  = AX_loc[XLLMX];
          auto ADXLLMY  = AX_loc[XLLMY];
          auto ADXLLMZ  = AX_loc[XLLMZ];
          auto DXLLMS = matList[XLLMS].X;
          auto DXLLMX = matList[XLLMX].X;
          auto DXLLMY = matList[XLLMY].X;
          auto DXLLMZ = matList[XLLMZ].X;

          auto ADXSSMS  = AX_loc[XSSMS];
          auto ADXSSMX  = AX_loc[XSSMX];
          auto ADXSSMY  = AX_loc[XSSMY];
          auto ADXSSMZ  = AX_loc[XSSMZ];
          auto DXSSMS = matList[XSSMS].X;
          auto DXSSMX = matList[XSSMX].X;
          auto DXSSMY = matList[XSSMY].X;
          auto DXSSMZ = matList[XSSMZ].X;

          auto ADCLSMS  = AX_loc[CLSMS];
          auto ADCLSMX  = AX_loc[CLSMX];
          auto ADCLSMY  = AX_loc[CLSMY];
          auto ADCLSMZ  = AX_loc[CLSMZ];
          auto ADXLSMS  = AX_loc[XLSMS];
          auto ADXLSMX  = AX_loc[XLSMX];
          auto ADXLSMY  = AX_loc[XLSMY];
          auto ADXLSMZ  = AX_loc[XLSMZ];

          auto DXLSMS = matList[XLSMS].X;
          auto DXLSMX = matList[XLSMX].X;
          auto DXLSMY = matList[XLSMY].X;
          auto DXLSMZ = matList[XLSMZ].X;
          auto DXSLMS = matList[XSLMS].X;
          auto DXSLMX = matList[XSLMX].X;
          auto DXSLMY = matList[XSLMY].X;
          auto DXSLMZ = matList[XSLMZ].X;

          auto DCLSMS = matList[CLSMS].X;
          auto DCLSMX = matList[CLSMX].X;
          auto DCLSMY = matList[CLSMY].X;
          auto DCLSMZ = matList[CLSMZ].X;
          auto DCSLMS = matList[CSLMS].X;
          auto DCSLMX = matList[CSLMX].X;
          auto DCSLMY = matList[CSLMY].X;
          auto DCSLMZ = matList[CSLMZ].X;
 
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

            //ERI0 :   ∇B∙∇C(mn|kl)
            
            /* Start of Gaunt (LL|LL) */
            /*++++++++++++++++++++++++*/

            // Exchange
            // MNKL

            // Equation (D5-D8) in Xiaosong's note
            if(bf1_s >= bf4_s) 
              ADXLLMS[bf14] += 3.0*DXSSMS[bf23]*BC[MNKL];

            if(bf1_s < bf4_s or (bf1_s==bf4_s and bf2_s!=bf3_s)) 
              ADXLLMS[bf41] += 3.0*DXSSMS[bf32]*BC[LKNM];


            /*++++++++++++++++++++++++*/
            /* Start of Gaunt (SS|SS) */
            /*++++++++++++++++++++++++*/
            // Equation (D9-D12) in Xiaosong's note
            if(bf2_s >= bf3_s) 
              ADXSSMS[bf23]+= 3.0*DXLLMS[bf14]*BC[LKNM];

            if(bf1_s!=bf4_s and bf2_s==bf3_s)
              ADXSSMS[bf32]+= 3.0*DXLLMS[bf41]*BC[MNKL];

           
            /*++++++++++++++++++++++++*/
            /* Start of Gaunt (LL|SS) */
            /*++++++++++++++++++++++++*/

            //Coulomb 
            // MNKL
            // TRANSKL
	    // Equation (D13-D16) in Xiaosong's note
            ADCLSMS[bf12]+= (-DCLSMS[bf43] +DCSLMS[bf34])*BC[MNKL];
            ADCLSMZ[bf12]+= -(DCLSMZ[bf43] +DCSLMZ[bf34])*BC[MNKL];
            ADCLSMX[bf12]+= -(DCLSMX[bf43] +DCSLMX[bf34])*BC[MNKL];
            ADCLSMY[bf12]+= -(DCLSMY[bf43] +DCSLMY[bf34])*BC[MNKL];
            // LKNM
            // Coulomb //TRANSKL
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              ADCLSMS[bf43]+= (-DCLSMS[bf12] +DCSLMS[bf21])*BC[LKNM];
              ADCLSMZ[bf43]+= -(DCLSMZ[bf12] +DCSLMZ[bf21])*BC[LKNM];
              ADCLSMX[bf43]+= -(DCLSMX[bf12] +DCSLMX[bf21])*BC[LKNM];
              ADCLSMY[bf43]+= -(DCLSMY[bf12] +DCSLMY[bf21])*BC[LKNM];
            }
            
            // Exchange            
            // TRANSKL
	    // Equation (D17-D20) in Xiaosong's note
            ADXLSMS[bf13]+=      DXSLMS[bf24]*BC[MNKL];
            ADXLSMZ[bf13]+= -2.0*DXSLMZ[bf24]*BC[MNKL];
            ADXLSMX[bf13]+= -2.0*DXSLMX[bf24]*BC[MNKL];
            ADXLSMY[bf13]+= -2.0*DXSLMY[bf24]*BC[MNKL];
            
            // TRANSKL
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              ADXLSMS[bf42]+=      DXSLMS[bf31]*BC[LKNM];
              ADXLSMZ[bf42]+= -2.0*DXSLMZ[bf31]*BC[LKNM]; 
              ADXLSMX[bf42]+= -2.0*DXSLMX[bf31]*BC[LKNM]; 
              ADXLSMY[bf42]+= -2.0*DXSLMY[bf31]*BC[LKNM]; 
            }

          }
          }
          }
          } // contraction loop
  
#endif // Contraction
  
  
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1

#ifdef _THREAD_TIMING_
        durThread[thread_id] = tock(gauntBegin);
#endif
  
      } // OpenMP context


      for( auto iTh  = 1; iTh < nThreads; iTh++) {

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMS],nBasis,MatsT(1.0),
          AXthreads[0][XLLMS],nBasis,AXthreads[0][XLLMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMX],nBasis,MatsT(1.0),
          AXthreads[0][XLLMX],nBasis,AXthreads[0][XLLMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMY],nBasis,MatsT(1.0),
          AXthreads[0][XLLMY],nBasis,AXthreads[0][XLLMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMZ],nBasis,MatsT(1.0),
          AXthreads[0][XLLMZ],nBasis,AXthreads[0][XLLMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMS],nBasis,MatsT(1.0),
          AXthreads[0][XSSMS],nBasis,AXthreads[0][XSSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMX],nBasis,MatsT(1.0),
          AXthreads[0][XSSMX],nBasis,AXthreads[0][XSSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMY],nBasis,MatsT(1.0),
          AXthreads[0][XSSMY],nBasis,AXthreads[0][XSSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMZ],nBasis,MatsT(1.0),
          AXthreads[0][XSSMZ],nBasis,AXthreads[0][XSSMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMS],nBasis,MatsT(1.0),
          AXthreads[0][CLSMS],nBasis,AXthreads[0][CLSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMX],nBasis,MatsT(1.0),
          AXthreads[0][CLSMX],nBasis,AXthreads[0][CLSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMY],nBasis,MatsT(1.0),
          AXthreads[0][CLSMY],nBasis,AXthreads[0][CLSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMZ],nBasis,MatsT(1.0),
          AXthreads[0][CLSMZ],nBasis,AXthreads[0][CLSMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMS],nBasis,MatsT(1.0),
          AXthreads[0][XLSMS],nBasis,AXthreads[0][XLSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMX],nBasis,MatsT(1.0),
          AXthreads[0][XLSMX],nBasis,AXthreads[0][XLSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMY],nBasis,MatsT(1.0),
          AXthreads[0][XLSMY],nBasis,AXthreads[0][XLSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMZ],nBasis,MatsT(1.0),
          AXthreads[0][XLSMZ],nBasis,AXthreads[0][XLSMZ],nBasis);

      };

      #ifdef CQ_ENABLE_MPI
      // Combine all G[X] contributions onto all processes
      if( mpiSize > 1 ) {
        MatsT* mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);

        std::vector<DIRAC_PAULI_SPINOR_COMP> comps{XLLMS, XLLMX, XLLMY, XLLMZ,
                                                   XSSMS, XSSMX, XSSMY, XSSMZ,
                                                   CLSMS, CLSMX, CLSMY, CLSMZ,
                                                   XLSMS, XLSMX, XLSMY, XLSMZ};
        for (auto comp : comps)
          MPIAllReduceInPlace( AXthreads[0][comp], nBasis*nBasis, comm, mpiScr );

        CQMemManager::get().free(mpiScr);

      }
      #endif

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMS],nBasis,MatsT(1.0),
             matList[XLLMS].AX,nBasis,matList[XLLMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMX],nBasis,MatsT(1.0),
             matList[XLLMX].AX,nBasis,matList[XLLMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMY],nBasis,MatsT(1.0),
             matList[XLLMY].AX,nBasis,matList[XLLMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMZ],nBasis,MatsT(1.0),
             matList[XLLMZ].AX,nBasis,matList[XLLMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMS],nBasis,MatsT(1.0),
             matList[XSSMS].AX,nBasis,matList[XSSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMX],nBasis,MatsT(1.0),
             matList[XSSMX].AX,nBasis,matList[XSSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMY],nBasis,MatsT(1.0),
             matList[XSSMY].AX,nBasis,matList[XSSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMZ],nBasis,MatsT(1.0),
             matList[XSSMZ].AX,nBasis,matList[XSSMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMS],nBasis,MatsT(1.0),
             matList[CLSMS].AX,nBasis,matList[CLSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMX],nBasis,MatsT(1.0),
             matList[CLSMX].AX,nBasis,matList[CLSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMY],nBasis,MatsT(1.0),
             matList[CLSMY].AX,nBasis,matList[CLSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMZ],nBasis,MatsT(1.0),
             matList[CLSMZ].AX,nBasis,matList[CLSMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMS],nBasis,MatsT(1.0),
             matList[XLSMS].AX,nBasis,matList[XLSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMX],nBasis,MatsT(1.0),
             matList[XLSMX].AX,nBasis,matList[XLSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMY],nBasis,MatsT(1.0),
             matList[XLSMY].AX,nBasis,matList[XLSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMZ],nBasis,MatsT(1.0),
             matList[XLSMZ].AX,nBasis,matList[XLSMZ].AX,nBasis);


#if 1

      // Take care of the Hermitian symmetry in the LL and SS blocks
      auto ADXLLMS = matList[XLLMS].AX;
      auto ADXLLMX = matList[XLLMX].AX;
      auto ADXLLMY = matList[XLLMY].AX;
      auto ADXLLMZ = matList[XLLMZ].AX;
      auto ADXSSMS = matList[XSSMS].AX;
      auto ADXSSMX = matList[XSSMX].AX;
      auto ADXSSMY = matList[XSSMY].AX;
      auto ADXSSMZ = matList[XSSMZ].AX;

      for( auto i = 0; i < nBasis; i++ )
      for( auto j = 0; j < i; j++ ) {

        ADXLLMS[j + i*nBasis] = std::conj(ADXLLMS[i + j*nBasis]);
        ADXLLMX[j + i*nBasis] = std::conj(ADXLLMX[i + j*nBasis]);
        ADXLLMY[j + i*nBasis] = std::conj(ADXLLMY[i + j*nBasis]);
        ADXLLMZ[j + i*nBasis] = std::conj(ADXLLMZ[i + j*nBasis]);
        ADXSSMS[j + i*nBasis] = std::conj(ADXSSMS[i + j*nBasis]);
        ADXSSMZ[j + i*nBasis] = std::conj(ADXSSMZ[i + j*nBasis]);
        ADXSSMX[j + i*nBasis] = std::conj(ADXSSMX[i + j*nBasis]);
        ADXSSMY[j + i*nBasis] = std::conj(ADXSSMY[i + j*nBasis]);

      }
#endif
  
      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
#ifdef _SHZ_SCREEN_4C_LIBCINT
      if(ShBlkNorms_raw!=nullptr) CQMemManager::get().free(ShBlkNorms_raw);
      if(SchwarzGaunt!=nullptr) CQMemManager::get().free(SchwarzGaunt);
#endif

#ifdef _THREAD_TIMING_
      std::cout << "Gaunt Libcint time on every thread:" << std::endl;
      for (size_t i = 0; i < nThreads; i++)
        std::cout << i << "\t:" << durThread[i] << std::endl;
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipGaunt = std::accumulate(nSkipGaunt.begin(),nSkipGaunt.end(),size_t(0));
#ifdef CQ_ENABLE_MPI
      if (mpiSize > 1) {
        std::cout << "Gaunt Screened "
                  << nIntSkipGaunt << " on Rank " << mpiRank <<  std::endl;
        nIntSkipGaunt = MPIAllReduce( nIntSkipGaunt, comm );
      }
#endif
      std::cout << "Gaunt Screened " << nIntSkipGaunt << std::endl;
  
      auto durDirectGaunt = tock(topDirectGaunt);
      std::cout << "Spin-Free Gaunt AO Direct Contraction took " <<  durDirectGaunt << " s\n"; 
  
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
  
      auto ss = 0     ;        // ERI00 (ss)+(σ∙σ)(ij|kl)
      auto dot = NB4 ;         // ERI01 (σ∙σ)(ijkl)

#ifdef _REPORT_INTEGRAL_TIMINGS
      auto topDirectGauge = tick();
#endif

      nERI = 4;
      buffAll = CQMemManager::get().malloc<double>(2*nERI*buffN4*nThreads);
      cacheAll = CQMemManager::get().malloc<double>(cache_size*nThreads);


#ifdef _SHZ_SCREEN_4C_LIBCINT
  
      double *SchwarzGauge = CQMemManager::get().malloc<double>(nShell*nShell);
      memset(SchwarzGauge,0,nShell*nShell*sizeof(double));
  
      #pragma omp parallel
      {
  
        double C1 = 1./(2*SpeedOfLight);
        size_t thread_id = GetThreadID();
        size_t mpi_thread_id = mpiRank * nThreads + thread_id;
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
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s12 % mpi_thread_size != mpi_thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s1);
          shls[3] = int(s2);
  
          auto nQuad = n1*n2*n1*n2;

          //∇A∇C
          //int2e_gauge_r1_sps1sps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          //int2e_gauge_r2_sps1sps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
	  
          //∇B∇D
          skiperi1 = int2e_gauge_r1_sp1sp2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
          skiperi2 = int2e_gauge_r2_sp1sp2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);

          if(skiperi1==0 and skiperi2==0) continue;

          for(size_t iERI(0); iERI < nERI*nQuad; iERI++) 
              SchwarzGauge[s1+s2*nShell] = std::max(SchwarzGauge[s1+s2*nShell], std::abs(buffr1[iERI]-buffr2[iERI]));
  
          SchwarzGauge[s1+s2*nShell] = C1*std::sqrt(SchwarzGauge[s1+s2*nShell]);
  
        }
        }
  
      };

#ifdef CQ_ENABLE_MPI
      MPIAllReduceInPlace(SchwarzGauge, nShell*nShell, comm);
#endif

  

      // Compute shell block norms (∞-norm) of matList.X
      int mMat = 16;
      double *ShBlkNorms_raw = CQMemManager::get().malloc<double>(mMat*nShell*nShell);
      std::vector<double*> ShBlkNorms;
  
      auto iOff = 0;
      ShellBlockNorm(basisSet_.shells,matList[XLSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSLMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);

      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XLLMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
 
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMS].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMX].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMY].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      iOff += nShell*nShell;
      ShellBlockNorm(basisSet_.shells,matList[XSSMZ].X,nBasis,ShBlkNorms_raw + iOff);
      ShBlkNorms.emplace_back(ShBlkNorms_raw + iOff);
  
      // Get the max over all the matricies for the shell block ∞-norms
      for(auto k = 0; k < nShell*nShell; k++) {
  
        double mx = std::abs(ShBlkNorms[0][k]);
        for(auto iMat = 0; iMat < mMat; iMat++)
          mx = std::max(mx,std::abs(ShBlkNorms[iMat][k]));
        ShBlkNorms[0][k] = mx;
   
      }

#endif

      int nSave = 2;
      ERIBuffer = CQMemManager::get().malloc<double>(nSave*NB4*nThreads);
      memset(AXRaw,0,nThreads*nMat*nBasis*nBasis*sizeof(MatsT));
      std::vector<size_t> nSkipGauge(nThreads,0);

#ifdef _THREAD_TIMING_
      std::vector<double> durThread(nThreads,0);
#endif
  
      #pragma omp parallel
      {
#ifdef _THREAD_TIMING_
        auto gaugeBegin = tick();
#endif
  
        dcomplex iS = dcomplex(0.0, 1.0);
  
        size_t thread_id = GetThreadID();
        int mpi_thread_id = mpiRank * nThreads + thread_id;
  
        auto &AX_loc = AXthreads[thread_id];
  
        double *BC   = &ERIBuffer[thread_id*nSave*NB4];
  
        size_t n1,n2,n3,n4,m,n,k,l,mnkl,bf1,bf2,bf3,bf4;
        size_t s4_max;

        int skiperi1,skiperi2;
  
        int shls[4];
        double *buffr1 = buffAll+nERI*buffN4*thread_id;
        double *buffr2 = buffAll+nERI*nThreads*buffN4+nERI*buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShell; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShell; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax12 = ShBlkNorms[0][s2 + s1*nShell];
#endif
  
        for(size_t s3(0), bf3_s(0); s3 <= s2 ; bf3_s+=n3, s3++) {
        //for(size_t s3(0), bf3_s(0); s3 < nShell ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s2 == s3) ? s1 : nShell-1; // Determine the unique max of Shell 4
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
          double shMax123 = std::max(ShBlkNorms[0][s1 + s3*nShell], ShBlkNorms[0][s2 + s3*nShell]);
          shMax123 = std::max(shMax123,shMax12);
#endif
   
  
        //for(size_t s4(0), bf4_s(0); s4 < nShell; bf4_s+=n4, s4++, s1234++) {
        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #if defined(_OPENMP) || defined(CQ_ENABLE_MPI)
          if( s1234 % mpi_thread_size != mpi_thread_id ) continue;
          #endif
  
#ifdef _SHZ_SCREEN_4C_LIBCINT
  
          double shMax = std::max(ShBlkNorms[0][s1 + s4*nShell],
                                  std::max(ShBlkNorms[0][s2 + s4*nShell],
                                           ShBlkNorms[0][s3 + s4*nShell]));
  
          shMax = std::max(shMax,shMax123);

          if((shMax*SchwarzGauge[s1+s2*nShell]*SchwarzGauge[s4+s3*nShell]) <
             eri.threshSchwarz()) { nSkipGauge[thread_id]++; continue; }
#endif

  
          if(approximate4C == APPROXIMATION_TYPE_4C::ThreeCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>3) {nSkipGauge[thread_id]++; continue;}
          }
 

          if(approximate4C == APPROXIMATION_TYPE_4C::TwoCenter) {
          //if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4)) ) 
            std::vector<int> atomCenters;
            std::vector<int>::iterator itAtom;

            atomCenters.push_back(bas(ATOM_OF, s1));

            if(not(bas(ATOM_OF, s1)==bas(ATOM_OF, s2))) atomCenters.push_back(bas(ATOM_OF, s2));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s3));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s3));

            itAtom = std::find(atomCenters.begin(), atomCenters.end(), bas(ATOM_OF, s4));
            if (itAtom == atomCenters.end()) atomCenters.push_back(bas(ATOM_OF, s4));

            if(atomCenters.size()>2) {nSkipGauge[thread_id]++; continue;}
          }
  
          if(approximate4C == APPROXIMATION_TYPE_4C::OneCenter) 
          if(not( bas(ATOM_OF, s1)==bas(ATOM_OF, s2) and bas(ATOM_OF, s3)==bas(ATOM_OF, s4) 
                 and bas(ATOM_OF, s1)==bas(ATOM_OF, s3) ) )
            {nSkipGauge[thread_id]++; continue;}

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
          for(m = 0             ; m <                  n1; ++m, ++mnkl) {

            auto MNKL = m + n*NB + k*NB2 + l*NB3;
            auto LKNM = l + k*NB + n*NB2 + m*NB3;

            //enum Gauge_ERI {
            //  SxSx,  // σ_x * σ_x     0
            //  SySy,  // σ_y * σ_y     1
            //  SzSz,  // σ_z * σ_z     2
            //  II     // I   * I       3
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

          } // ∇B∇C integral preparation loop
  
  
  
#ifdef _CONTRACTION_ // Contraction

          auto ADXLLMS  = AX_loc[XLLMS];
          auto ADXLLMX  = AX_loc[XLLMX];
          auto ADXLLMY  = AX_loc[XLLMY];
          auto ADXLLMZ  = AX_loc[XLLMZ];
          auto DXLLMS = matList[XLLMS].X;
          auto DXLLMX = matList[XLLMX].X;
          auto DXLLMY = matList[XLLMY].X;
          auto DXLLMZ = matList[XLLMZ].X;

          auto ADXSSMS  = AX_loc[XSSMS];
          auto ADXSSMX  = AX_loc[XSSMX];
          auto ADXSSMY  = AX_loc[XSSMY];
          auto ADXSSMZ  = AX_loc[XSSMZ];
          auto DXSSMS = matList[XSSMS].X;
          auto DXSSMX = matList[XSSMX].X;
          auto DXSSMY = matList[XSSMY].X;
          auto DXSSMZ = matList[XSSMZ].X;

          auto ADCLSMS  = AX_loc[CLSMS];
          auto ADCLSMX  = AX_loc[CLSMX];
          auto ADCLSMY  = AX_loc[CLSMY];
          auto ADCLSMZ  = AX_loc[CLSMZ];
          auto ADXLSMS  = AX_loc[XLSMS];
          auto ADXLSMX  = AX_loc[XLSMX];
          auto ADXLSMY  = AX_loc[XLSMY];
          auto ADXLSMZ  = AX_loc[XLSMZ];

          auto DXLSMS = matList[XLSMS].X;
          auto DXLSMX = matList[XLSMX].X;
          auto DXLSMY = matList[XLSMY].X;
          auto DXLSMZ = matList[XLSMZ].X;
          auto DXSLMS = matList[XSLMS].X;
          auto DXSLMX = matList[XSLMX].X;
          auto DXSLMY = matList[XSLMY].X;
          auto DXSLMZ = matList[XSLMZ].X;

          auto DCLSMS = matList[CLSMS].X;
          auto DCLSMX = matList[CLSMX].X;
          auto DCLSMY = matList[CLSMY].X;
          auto DCLSMZ = matList[CLSMZ].X;
          auto DCSLMS = matList[CSLMS].X;
          auto DCSLMX = matList[CSLMX].X;
          auto DCSLMY = matList[CSLMY].X;
          auto DCSLMZ = matList[CSLMZ].X;
 
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

            //      auto ss = 0     ;        // ERI00 (ss)(ij|kl)
            //      auto dot = NB4 ;         // ERI01 (σ∙σ)(ijkl)
             
            /*++++++++++++++++++++++++*/
            /* Start of Gauge (LL|LL) */
            /*++++++++++++++++++++++++*/
            
            // Exchange
            // MNKL
	    
            // Follow equations in the gauge publication
	    // For the excahange part, the 1/2 prefactor is taken care of during
	    // the matrix assembly

            if(bf1_s >= bf4_s) {
              /* Equation (E12-E15) in Xiaosong's note */
              ADXLLMS[bf14] += DXSSMS[bf23]*(BC[MNKL] + BC[MNKL+dot]);
              ADXLLMZ[bf14] += DXSSMZ[bf23]*BC[MNKL];
              ADXLLMX[bf14] += DXSSMX[bf23]*BC[MNKL];
              ADXLLMY[bf14] += DXSSMY[bf23]*BC[MNKL];
            }
            
            // LKNM
            if(bf1_s < bf4_s or (bf1_s==bf4_s and bf2_s!=bf3_s)) {
              ADXLLMS[bf41] += DXSSMS[bf32]*(BC[LKNM] + BC[LKNM+dot]);
              ADXLLMZ[bf41] += DXSSMZ[bf32]*BC[LKNM];
              ADXLLMX[bf41] += DXSSMX[bf32]*BC[LKNM];
              ADXLLMY[bf41] += DXSSMY[bf32]*BC[LKNM];
            }
            
            
            /*++++++++++++++++++++++++*/
            /* Start of Gauge (SS|SS) */
            /*++++++++++++++++++++++++*/
            
            // Exchange
            // MNKL: TRANS_MN_TRANS_KL
            // See Equation (A15)/(E10) for integral symmetry
            
            if(bf2_s >= bf3_s) {
              /* Equation (E16-E19) in Xiaosong's note */
              ADXSSMS[bf23] += DXLLMS[bf14]*(BC[MNKL] + BC[MNKL+dot]);
              ADXSSMZ[bf23] += DXLLMZ[bf14]*BC[MNKL];
              ADXSSMX[bf23] += DXLLMX[bf14]*BC[MNKL];
              ADXSSMY[bf23] += DXLLMY[bf14]*BC[MNKL];
            }

            // LKNM, TRANS_MN_TRANS_KL
            if(bf1_s!=bf4_s and bf2_s==bf3_s) {
              ADXSSMS[bf32] += DXLLMS[bf41]*(BC[LKNM] + BC[LKNM+dot]);
              ADXSSMZ[bf32] += DXLLMZ[bf41]*BC[LKNM];
              ADXSSMX[bf32] += DXLLMX[bf41]*BC[LKNM];
              ADXSSMY[bf32] += DXLLMY[bf41]*BC[LKNM];
            }
             
             
            /*++++++++++++++++++++++++*/
            /* Start of Gauge (LL|SS) */
            /*++++++++++++++++++++++++*/
             
            
            // Coulomb
            // MNKL
            // See Equation (A16)/(E11) for integral symmetry
            
            /* Equations (E20-E23) in Xiaosong's note */
            ADCLSMS[bf12] += (DXSLMS[bf34] - DXLSMS[bf43])*BC[MNKL];
           
            // LKNM
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) 
              ADCLSMS[bf43] += (DXSLMS[bf21] - DXLSMS[bf12])*BC[LKNM];
             
            // Exchange
            // MNKL, TRANS_KL
            // See Equation (A16)/(E11) for integral symmetry
            
            /* Equation (E24-E27) in Xiaosong's note */
            ADXLSMS[bf13] +=  DXSLMS[bf24]*(BC[MNKL+dot] - BC[MNKL]);
            ADXLSMZ[bf13] += -DXSLMZ[bf24]*BC[MNKL];
            ADXLSMX[bf13] += -DXSLMX[bf24]*BC[MNKL];
            ADXLSMY[bf13] += -DXSLMY[bf24]*BC[MNKL];
            
            // LKNM, TRANS_KL
            if(bf1_s!=bf4_s or bf2_s!=bf3_s) {
              ADXLSMS[bf42] +=  DXSLMS[bf31]*(BC[LKNM+dot] - BC[LKNM]);
              ADXLSMZ[bf42] += -DXSLMZ[bf31]*BC[LKNM];
              ADXLSMX[bf42] += -DXSLMX[bf31]*BC[LKNM];
              ADXLSMY[bf42] += -DXSLMY[bf31]*BC[LKNM];
            }


          }
          }
          }
          } // contraction loop
  
#endif // Contraction
  
  
        }; // loop s4
        }; // loop s3
        }; // loop s2
        }; // loop s1

#ifdef _THREAD_TIMING_
        durThread[thread_id] = tock(gaugeBegin);
#endif
  
      } // OpenMP context


      for( auto iTh  = 1; iTh < nThreads; iTh++) {

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMS],nBasis,MatsT(1.0),
          AXthreads[0][XLLMS],nBasis,AXthreads[0][XLLMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMX],nBasis,MatsT(1.0),
          AXthreads[0][XLLMX],nBasis,AXthreads[0][XLLMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMY],nBasis,MatsT(1.0),
          AXthreads[0][XLLMY],nBasis,AXthreads[0][XLLMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLLMZ],nBasis,MatsT(1.0),
          AXthreads[0][XLLMZ],nBasis,AXthreads[0][XLLMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMS],nBasis,MatsT(1.0),
          AXthreads[0][XSSMS],nBasis,AXthreads[0][XSSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMX],nBasis,MatsT(1.0),
          AXthreads[0][XSSMX],nBasis,AXthreads[0][XSSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMY],nBasis,MatsT(1.0),
          AXthreads[0][XSSMY],nBasis,AXthreads[0][XSSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XSSMZ],nBasis,MatsT(1.0),
          AXthreads[0][XSSMZ],nBasis,AXthreads[0][XSSMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMS],nBasis,MatsT(1.0),
          AXthreads[0][CLSMS],nBasis,AXthreads[0][CLSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMX],nBasis,MatsT(1.0),
          AXthreads[0][CLSMX],nBasis,AXthreads[0][CLSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMY],nBasis,MatsT(1.0),
          AXthreads[0][CLSMY],nBasis,AXthreads[0][CLSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][CLSMZ],nBasis,MatsT(1.0),
          AXthreads[0][CLSMZ],nBasis,AXthreads[0][CLSMZ],nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMS],nBasis,MatsT(1.0),
          AXthreads[0][XLSMS],nBasis,AXthreads[0][XLSMS],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMX],nBasis,MatsT(1.0),
          AXthreads[0][XLSMX],nBasis,AXthreads[0][XLSMX],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMY],nBasis,MatsT(1.0),
          AXthreads[0][XLSMY],nBasis,AXthreads[0][XLSMY],nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[iTh][XLSMZ],nBasis,MatsT(1.0),
          AXthreads[0][XLSMZ],nBasis,AXthreads[0][XLSMZ],nBasis);

      };

      #ifdef CQ_ENABLE_MPI
      // Combine all G[X] contributions onto all processes
      if( mpiSize > 1 ) {
        MatsT* mpiScr = CQMemManager::get().malloc<MatsT>(nBasis*nBasis);

        std::vector<DIRAC_PAULI_SPINOR_COMP> comps{XLLMS, XLLMX, XLLMY, XLLMZ,
                                                   XSSMS, XSSMX, XSSMY, XSSMZ,
                                                   CLSMS, CLSMX, CLSMY, CLSMZ,
                                                   XLSMS, XLSMX, XLSMY, XLSMZ};
        for (auto comp : comps)
          MPIAllReduceInPlace( AXthreads[0][comp], nBasis*nBasis, comm, mpiScr );

        CQMemManager::get().free(mpiScr);

      }
      #endif

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMS],nBasis,MatsT(1.0),
             matList[XLLMS].AX,nBasis,matList[XLLMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMX],nBasis,MatsT(1.0),
             matList[XLLMX].AX,nBasis,matList[XLLMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMY],nBasis,MatsT(1.0),
             matList[XLLMY].AX,nBasis,matList[XLLMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLLMZ],nBasis,MatsT(1.0),
             matList[XLLMZ].AX,nBasis,matList[XLLMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMS],nBasis,MatsT(1.0),
             matList[XSSMS].AX,nBasis,matList[XSSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMX],nBasis,MatsT(1.0),
             matList[XSSMX].AX,nBasis,matList[XSSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMY],nBasis,MatsT(1.0),
             matList[XSSMY].AX,nBasis,matList[XSSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XSSMZ],nBasis,MatsT(1.0),
             matList[XSSMZ].AX,nBasis,matList[XSSMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMS],nBasis,MatsT(1.0),
             matList[CLSMS].AX,nBasis,matList[CLSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMX],nBasis,MatsT(1.0),
             matList[CLSMX].AX,nBasis,matList[CLSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMY],nBasis,MatsT(1.0),
             matList[CLSMY].AX,nBasis,matList[CLSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][CLSMZ],nBasis,MatsT(1.0),
             matList[CLSMZ].AX,nBasis,matList[CLSMZ].AX,nBasis);



      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMS],nBasis,MatsT(1.0),
             matList[XLSMS].AX,nBasis,matList[XLSMS].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMX],nBasis,MatsT(1.0),
             matList[XLSMX].AX,nBasis,matList[XLSMX].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMY],nBasis,MatsT(1.0),
             matList[XLSMY].AX,nBasis,matList[XLSMY].AX,nBasis);

      MatAdd('N','N',nBasis,nBasis,MatsT(1.0),AXthreads[0][XLSMZ],nBasis,MatsT(1.0),
             matList[XLSMZ].AX,nBasis,matList[XLSMZ].AX,nBasis);

 #if 1

      // Take care of the Hermitian symmetry in the LL and SS blocks
      auto ADXLLMS = matList[XLLMS].AX;
      auto ADXLLMX = matList[XLLMX].AX;
      auto ADXLLMY = matList[XLLMY].AX;
      auto ADXLLMZ = matList[XLLMZ].AX;
      auto ADXSSMS = matList[XSSMS].AX;
      auto ADXSSMX = matList[XSSMX].AX;
      auto ADXSSMY = matList[XSSMY].AX;
      auto ADXSSMZ = matList[XSSMZ].AX;

      for( auto i = 0; i < nBasis; i++ )
      for( auto j = 0; j < i; j++ ) {

        ADXLLMS[j + i*nBasis] = std::conj(ADXLLMS[i + j*nBasis]);
        ADXLLMX[j + i*nBasis] = std::conj(ADXLLMX[i + j*nBasis]);
        ADXLLMY[j + i*nBasis] = std::conj(ADXLLMY[i + j*nBasis]);
        ADXLLMZ[j + i*nBasis] = std::conj(ADXLLMZ[i + j*nBasis]);
        ADXSSMS[j + i*nBasis] = std::conj(ADXSSMS[i + j*nBasis]);
        ADXSSMZ[j + i*nBasis] = std::conj(ADXSSMZ[i + j*nBasis]);
        ADXSSMX[j + i*nBasis] = std::conj(ADXSSMX[i + j*nBasis]);
        ADXSSMY[j + i*nBasis] = std::conj(ADXSSMY[i + j*nBasis]);

      }
#endif
   
      CQMemManager::get().free(ERIBuffer);
      CQMemManager::get().free(buffAll, cacheAll);
#ifdef _SHZ_SCREEN_4C_LIBCINT
      if(ShBlkNorms_raw!=nullptr) CQMemManager::get().free(ShBlkNorms_raw);
      if(SchwarzGauge!=nullptr) CQMemManager::get().free(SchwarzGauge);
#endif

#ifdef _THREAD_TIMING_
      std::cout << "Gauge Libcint time on every thread:" << std::endl;
      for (size_t i = 0; i < nThreads; i++)
        std::cout << i << "\t:" << durThread[i] << std::endl;
#endif
  
#ifdef _REPORT_INTEGRAL_TIMINGS
      size_t nIntSkipGauge = std::accumulate(nSkipGauge.begin(),nSkipGauge.end(),size_t(0));
#ifdef CQ_ENABLE_MPI
      if (mpiSize > 1) {
        std::cout << "Gauge Screened "
                  << nIntSkipGauge << " on Rank " << mpiRank <<  std::endl;
        nIntSkipGauge = MPIAllReduce( nIntSkipGauge, comm );
      }
#endif
      std::cout << "Gauge Screened " << nIntSkipGauge << std::endl;
  
      auto durDirectGauge = tock(topDirectGauge);
      std::cout << "Spin-Free Gauge AO Direct Contraction took " <<  durDirectGauge << " s\n"; 
  
      std::cout << std::endl;
#endif

    } // if(Gauge)

    /*******************/
    /*                 */
    /*   End of Gauge  */
    /*                 */
    /*******************/






#ifdef _SHZ_SCREEN_4C_LIBCINT
    if(SchwarzSSSS!=nullptr) CQMemManager::get().free(SchwarzSSSS);
    if(SchwarzERI!=nullptr) CQMemManager::get().free(SchwarzERI);
#endif

    if(AXRaw!=nullptr) CQMemManager::get().free(AXRaw);

    CQMemManager::get().free(atm);
    CQMemManager::get().free(bas);
    CQMemManager::get().free(env);


    // Turn threads for LA back on
    SetLAThreads(LAThreads);

  }

  template <>
  void GTODirectRelERIContraction<double,double>::directScaffoldLibcintSpinFree(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<double>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);  
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directScaffoldLibcintSpinFree(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<dcomplex>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);  
  }

}; // namespace ChronusQ











