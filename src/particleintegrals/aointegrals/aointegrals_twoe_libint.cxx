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

#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/gradints.hpp>
#include <libint2/engine.h>
#include <cqlinalg.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <particleintegrals/inhouseaointegral.hpp>

#include <util/threads.hpp>
#include <chrono>

#define __IN_HOUSE_INT_GAUGE__
//#define __DEBUGERI__
//#define __DEBUGGAUGE__

namespace ChronusQ {

  typedef std::vector<libint2::Shell> shell_set; 


  /**
   *  \brief Helper function to do the setup and loops over basis functions
   *  for evaluating in core ERIs. Generic to gradient level.
   */
  template <typename InnerFunc>
  void doERIInCore(size_t deriv, std::vector<double*> eris, InnerFunc func,
    BasisSet &basis, BasisSet& basis2, OPERATOR op,
    const HamiltonianOptions &options) {


    if (op != ELECTRON_REPULSION and op != EP_ATTRACTION)
      CErr("Only e-p attraction/e-e/p-p repulsion integrals in incore TPIs",std::cout);
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in InCoreEIRIs",std::cout);

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
 
    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      std::max(basis.maxPrim, basis2.maxPrim),
      std::max(basis.maxL, basis2.maxL), deriv);
    engines[0].set_precision(0.);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

    // Get useful constants
    bool sameBasis = (&basis == &basis2);
    size_t NB = basis.nBasis;
    size_t MB = basis2.nBasis;
    size_t NB2 = NB*NB;
    size_t MB2 = MB*MB;
    size_t NB2MB2 = NB2*MB2;

    // Clear previous ERIs
    std::for_each(eris.begin(), eris.end(),
      [&](double* p){std::fill_n(p,NB2MB2,0.);}
    );

    //
    //  Start parallel scaffold surrounding InnerFunc
    //

    auto libint_start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s3_max, s4_max;
      for(size_t s1(0), bf1_s(0), s1234(0); s1 < basis.nShell;
          bf1_s+=n1, s1++) { 

        n1 = basis.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basis.shells[s2].size(); // Size of Shell 2

#ifdef __IN_HOUSE_INT__
        libint2::ShellPair pair1_to_use;
        pair1_to_use.init( basis.shells[s1],basis.shells[s2],-2000);
#endif

        s3_max = sameBasis ? s1 : basis2.nShell - 1;

      for(size_t s3(0), bf3_s(0); s3 <= s3_max ; bf3_s+=n3, s3++) {

        n3 = basis2.shells[s3].size(); // Size of Shell 3
        s4_max = sameBasis && (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basis2.shells[s4].size(); // Size of Shell 4

#ifdef __IN_HOUSE_INT__
        libint2::ShellPair pair2_to_use;
        pair2_to_use.init( basis.shells[s3],basis.shells[s4],-2000);

#endif

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        //
        // Inner function - this is responsible for calling the appropriate
        //   computation routine and placing the ERI results in the correct
        //   location in memory
        //

        func(s1, s2, s3, s4,
             bf1_s, bf2_s, bf3_s, bf4_s,
             n1, n2, n3, n4,
             engines[thread_id]);

      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    auto libint_stop = std::chrono::high_resolution_clock::now();
    auto libint_duration = std::chrono::duration_cast<std::chrono::milliseconds>(libint_stop - libint_start);
    if ( false )
      std::cout << "Libint duration   = " << libint_duration.count() << std::endl;

    // Debug output of the ERIs

  };

  /**
   *  \brief Compute and store the full rank-4 ERI tensor using
   *  Libint2 over the CGTO basis.
   */ 
  template <>
  void InCore4indexTPI<double>::computeERINR(BasisSet &basisSet, BasisSet &basisSet2, 
      Molecule&, EMPerturbation&, OPERATOR op, const HamiltonianOptions &options) {

    InCore4indexTPI<double> &eri4I = *this;

    bool sameBasis = (&basisSet == &basisSet2);

    // Lambda that does the computation and placing. This gets called by
    //   doERIInCore
    auto computeAndPlace = [&](size_t sh1, size_t sh2, size_t sh3, size_t sh4,
                               size_t b1s, size_t b2s, size_t b3s, size_t b4s,
                               size_t  n1, size_t  n2, size_t  n3, size_t  n4,
                               libint2::Engine& engine) {

#ifndef __IN_HOUSE_INT__
      engine.compute2<
        libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
        basisSet.shells[sh1],
        basisSet.shells[sh2],
        basisSet2.shells[sh3],
        basisSet2.shells[sh4]
      );

      const double* results =  engine.results()[0] ;
      // Libint internal screening
      if(results == nullptr) return;
#else
      auto buff  = RealGTOIntEngine::BottomupHGP(pair1_to_use,pair2_to_use,
        basisSet.shells[sh1],
        basisSet.shells[sh2],
        basisSet2.shells[sh3],
        basisSet2.shells[sh4]
      );
#endif

      // Place shell quartet into persistent storage with
      // permutational symmetry
      for(size_t i = 0ul, bf1 = b1s, ijkl = 0ul ; i < n1; ++i, bf1++) 
      for(size_t j = 0ul, bf2 = b2s             ; j < n2; ++j, bf2++) 
      for(size_t k = 0ul, bf3 = b3s             ; k < n3; ++k, bf3++) 
      for(size_t l = 0ul, bf4 = b4s             ; l < n4; ++l, bf4++, ++ijkl) {

        // (12 | 34)
        eri4I(bf1, bf2, bf3, bf4) = results[ijkl];
        // (12 | 43)
        eri4I(bf1, bf2, bf4, bf3) = results[ijkl];
        // (21 | 34)
        eri4I(bf2, bf1, bf3, bf4) = results[ijkl];
        // (21 | 43)
        eri4I(bf2, bf1, bf4, bf3) = results[ijkl];

        // 8-fold symmetry only if left basis is the same as right basis
        if( sameBasis ) {
          // (34 | 12)
          eri4I(bf3, bf4, bf1, bf2) = results[ijkl];
          // (43 | 12)
          eri4I(bf4, bf3, bf1, bf2) = results[ijkl];
          // (34 | 21)
          eri4I(bf3, bf4, bf2, bf1) = results[ijkl];
          // (43 | 21)
          eri4I(bf4, bf3, bf2, bf1) = results[ijkl];
        }


      }; // ijkl loop

    };

    doERIInCore(0, {TPI}, computeAndPlace, basisSet, basisSet2, op, options);

#ifdef __DEBUGERI__
//#if 1
    std::cout << "Two-Electron Integrals (ERIs)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << TPI[i + j*NB  + k*NB2 + l*NB3] << std::endl;
    };
#endif


  }; // InCore4indexTPI<double>::computeTPI

  /**
   * @brief Compute ERI in a general-contraction shell quartet.
   *        Note layout difference, resPQRS need to be reorganized by
   *        the following code to match the layout of a general
   *        shell quartet output as in libcint:
   *        for (size_t SS = 0, pqrs = 0; SS < sContrSize; SS++)
   *          for (size_t s = 0; s < sAMSize; s++)
   *            for (size_t RR = 0; RR < rContrSize; RR++)
   *              for (size_t r = 0; r < rAMSize; r++)
   *                for (size_t QQ = 0; QQ < qContrSize; QQ++)
   *                  for (size_t q = 0; q < qAMSize; q++)
   *                    for (size_t PP = 0; PP < pContrSize; PP++)
   *                      for (size_t p = 0; p < pAMSize; p++, pqrs++) {
   *                        block[pqrs] = resPQRS[(s + sAMSize * (r + rAMSize * (q + qAMSize * p)))
   *                            + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
   *                      }
   * @return The count of engine.compute2 calls,
   *         and the total time of such calls
   */
  inline std::pair<size_t, double> libintGeneralContractionERI(
      size_t P, size_t Q, size_t R, size_t S,
      BasisSet &basisSet, libint2::Engine &engine,
      std::vector<std::vector<libint2::Shell>> &shellPrims_,
      std::vector<double*> &coefBlocks_, double *workBlock,
      const double *&resPQRS) {

    std::pair<size_t, double> counter_timer(0, 0.0);

    size_t pContrSize = basisSet.shells[P].contr.size();
    size_t qContrSize = basisSet.shells[Q].contr.size();
    size_t rContrSize = basisSet.shells[R].contr.size();
    size_t sContrSize = basisSet.shells[S].contr.size();

    std::vector<libint2::Shell> &primsP = shellPrims_[P];
    std::vector<libint2::Shell> &primsQ = shellPrims_[Q];
    std::vector<libint2::Shell> &primsR = shellPrims_[R];
    std::vector<libint2::Shell> &primsS = shellPrims_[S];

    size_t pAMSize = primsP[0].size();
    size_t qAMSize = primsQ[0].size();
    size_t rAMSize = primsR[0].size();
    size_t sAMSize = primsS[0].size();

    size_t pNprim = pContrSize > 1 ? primsP.size() : 1;
    size_t qNprim = qContrSize > 1 ? primsQ.size() : 1;
    size_t rNprim = rContrSize > 1 ? primsR.size() : 1;
    size_t sNprim = sContrSize > 1 ? primsS.size() : 1;

    size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

    // Compute primitive ERI

    const auto& buf_vec = engine.results();
    const double *buff = buf_vec[0];

    const double *inpP, *inpQ, *inpR, *inpS;
    double *resP, *resQ, *resR, *resS;

    double *sVec, *rVec, *qVec, *pVec;
    sVec = workBlock + pqrsAMSize * pContrSize * qContrSize * rContrSize * sContrSize;
    rVec = sVec + pqrsAMSize * pContrSize * qContrSize * rContrSize * sNprim;
    qVec = rVec + pqrsAMSize * pContrSize * qContrSize * rNprim;
    pVec = qVec + pqrsAMSize * pContrSize * qNprim;

    inpP = pVec;
    inpQ = qVec;
    inpR = rVec;
    inpS = sVec;

    for (size_t SS(0); SS < sNprim; SS++) {
      for (size_t RR(0); RR < rNprim; RR++) {
        for (size_t QQ(0); QQ < qNprim; QQ++) {
          for (size_t PP(0); PP < pNprim; PP++) {

            // Evaluate ERI for shell quartet
            auto beginERI = tick();
            engine.compute2<
              libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                  primsP[PP], primsQ[QQ], primsR[RR], primsS[SS]
            );
            counter_timer.second += tock(beginERI);
            counter_timer.first++;
            buff = buf_vec[0];
            if(buff == nullptr) continue;

            if (pContrSize > 1) {
              std::copy_n(buff, pqrsAMSize, &pVec[PP * pqrsAMSize]);
            } else
              inpP = buff;

          }

          if (qContrSize > 1 or pContrSize > 1) {
            resP = &qVec[QQ * pContrSize * pqrsAMSize];
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, pqrsAMSize, pContrSize, pNprim,
                 1.0, inpP, pqrsAMSize,
                 coefBlocks_[P], pNprim,
                 0.0, resP, pqrsAMSize);
          } else
            inpQ = inpP;

        }

        if (rContrSize > 1 or qContrSize > 1) {
          resQ = &rVec[RR * pContrSize * qContrSize * pqrsAMSize];
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, pqrsAMSize * pContrSize, qContrSize, qNprim,
               1.0, inpQ, pqrsAMSize * pContrSize,
               coefBlocks_[Q], qNprim,
               0.0, resQ, pqrsAMSize * pContrSize);
        } else
          inpR = inpQ;

      }

      if (sContrSize > 1 or rContrSize > 1) {
        resR = &sVec[SS * pContrSize * qContrSize * rContrSize * pqrsAMSize];
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, pqrsAMSize * pContrSize * qContrSize, rContrSize, rNprim,
             1.0, inpR, pqrsAMSize * pContrSize * qContrSize,
             coefBlocks_[R], rNprim,
             0.0, resR, pqrsAMSize * pContrSize * qContrSize);
      } else
        inpS = inpR;

    }

    if (sContrSize > 1) {
      resS = &workBlock[0];
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, pqrsAMSize * pContrSize * qContrSize * rContrSize, sContrSize, sNprim,
           1.0, inpS, pqrsAMSize * pContrSize * qContrSize * rContrSize,
           coefBlocks_[S], sNprim,
           0.0, resS, pqrsAMSize * pContrSize * qContrSize * rContrSize);
      resPQRS = resS;
    } else
      resPQRS = inpS;

    return counter_timer;
  }; // libintGeneralContractionERI


  template <>
  void InCore4indexTPI<dcomplex>::computeERIGCNR(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCore4indexTPI<double>::computeERIGCNR(BasisSet &originalBasisSet, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {

    if (op != ELECTRON_REPULSION)
      CErr("Only Electron repulsion integrals in InCore4indexTPI<double>",std::cout);
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in InCore4indexTPI<double>",std::cout);

    BasisSet basisSet = originalBasisSet.groupGeneralContractionBasis();

    bool segmented = basisSet.nShell == originalBasisSet.nShell;

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      basisSet.maxPrim,basisSet.maxL,0);
    engines[0].set_precision(0.);


    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

    this->clear();
    InCore4indexTPI<double> &eri4I = *this;

    std::vector<std::vector<libint2::Shell>> shellPrims; // Mappings from primitives to CGTOs
    std::vector<double*> coefBlocks; // Mappings from primitives to CGTOs
    std::vector<double*> workBlocks;

    size_t maxNcontrAMSize = 1, maxNprimAMSize = 1, maxAMSize = 1;

    if (not segmented) {

      // Compute the mappings from primitives to CGTOs
      BasisSet primitives(originalBasisSet.uncontractBasis());
      OnePInts<double> overlap(primitives.nBasis);

      overlap.computeAOInts(primitives, mol, emPert, OVERLAP, options);

      double *mapPrim2Cont = CQMemManager::get().malloc<double>(primitives.nBasis*originalBasisSet.nBasis);
      originalBasisSet.makeMapPrim2Cont(overlap.pointer(), mapPrim2Cont);

      coefBlocks.resize(basisSet.nShell, nullptr);
      shellPrims.reserve(basisSet.nShell);

      for (size_t P(0); P < basisSet.nShell; P++) {
        // Gather contraction coefficients
        libint2::Shell &shellP = basisSet.shells[P];
        maxNcontrAMSize = std::max(maxNcontrAMSize, shellP.size());
        size_t pContrSize = shellP.ncontr();
        size_t pAMSize = shellP.contr[0].size();
        if (pContrSize == 1) {
          maxNprimAMSize = std::max(maxNprimAMSize, pAMSize);
          shellPrims.emplace_back(1, shellP);
          coefBlocks[P] = CQMemManager::get().malloc<double>(1);
          coefBlocks[P][0] = 1.0;
          continue;
        }

        size_t pNprim = shellP.nprim();

        std::vector<libint2::Shell> primsP;
        primsP.reserve(pNprim);

        for(auto &a : shellP.alpha) {// Loop over primitives
          libint2::Shell newShell{
            { a },
            { { shellP.contr[0].l, shellP.contr[0].pure, { 1.0 } } },
            { { shellP.O[0], shellP.O[1], shellP.O[2] } }
          };
          primsP.push_back(std::move(newShell));
        }

        maxNprimAMSize = std::max(maxNprimAMSize, pNprim * pAMSize);

        size_t pBegin = basisSet.mapSh2Bf[P];
        coefBlocks[P] = CQMemManager::get().malloc<double>(pContrSize * pNprim);
        for (size_t c = 0; c < pContrSize; c++) {
          for (size_t i = 0; i < pNprim; i++) {
            coefBlocks[P][i + c * pNprim]
                = mapPrim2Cont[pBegin + c * pAMSize +
                                basisSet.primitives[primsP[i]] * basisSet.nBasis];
          }
        }
        shellPrims.push_back(std::move(primsP));
      }
      CQMemManager::get().free(mapPrim2Cont);

      workBlocks.resize(nthreads, nullptr);

      size_t maxL = basisSet.maxL;
      maxAMSize = basisSet.forceCart ?
                      (maxL + 1) * (maxL + 2) / 2 :
                      (2 * maxL + 1);
      size_t primAllocSize = maxNcontrAMSize * maxNcontrAMSize
          * maxNcontrAMSize * maxNcontrAMSize
          + maxNprimAMSize * maxNcontrAMSize
          * maxNcontrAMSize * maxNcontrAMSize
          + maxNprimAMSize * maxNcontrAMSize
          * maxNcontrAMSize * maxAMSize
          + maxNprimAMSize * maxNcontrAMSize
          * maxAMSize * maxAMSize
          + maxNprimAMSize * maxAMSize
          * maxAMSize * maxAMSize;
      for (size_t i = 0; i < nthreads; i++) {
        workBlocks[i] = CQMemManager::get().malloc<double>(primAllocSize);
      }
    }

#ifndef __IN_HOUSE_INT__
    //std::cout<<"  Using Libint "<<std::endl;
#else
    std::cout<<"  Using In-house Integral Engine "<<std::endl;
#endif

    auto topERI4 = tick();
    #pragma omp parallel
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      for(size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet.nShell;
          bf1_s+=n1, s1++) {

        n1 = basisSet.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet.shells[s2].size(); // Size of Shell 2

#ifdef __IN_HOUSE_INT__
        libint2::ShellPair pair1_to_use;
        pair1_to_use.init( basisSet_.shells[s1],basisSet_.shells[s2],-2000);
#endif

      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet.shells[s4].size(); // Size of Shell 4

#ifdef __IN_HOUSE_INT__
        libint2::ShellPair pair2_to_use;
        pair2_to_use.init( basisSet_.shells[s3],basisSet_.shells[s4],-2000);

#endif

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif


        if (segmented or basisSet.shells[s1].ncontr() == 1 and basisSet.shells[s2].ncontr() == 1
            and basisSet.shells[s3].ncontr() == 1 and basisSet.shells[s4].ncontr() == 1) {
          // Evaluate ERI for shell quartet
          engines[thread_id].compute2<
            libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
            basisSet.shells[s1],
            basisSet.shells[s2],
            basisSet.shells[s3],
            basisSet.shells[s4]
          );
          const auto *buff =  buf_vec[0] ;
          if(buff == nullptr) continue;

#ifdef __DEBUGERI__
          prettyPrintSmart(std::cout, "(" + std::to_string(s4) + "," + std::to_string(s3) + "|"
                                          + std::to_string(s2) + "," + std::to_string(s1) + ")",
                           buff, n4*n3, n2*n1, n4*n3);
#endif

          // Libint2 internal screening

          // Place shell quartet into persistent storage with
          // permutational symmetry
          for(i = 0ul, bf1 = bf1_s, ijkl = 0ul ; i < n1; ++i, bf1++)
          for(j = 0ul, bf2 = bf2_s             ; j < n2; ++j, bf2++)
          for(k = 0ul, bf3 = bf3_s             ; k < n3; ++k, bf3++)
          for(l = 0ul, bf4 = bf4_s             ; l < n4; ++l, bf4++, ++ijkl) {


              // (12 | 34)
              eri4I(bf1, bf2, bf3, bf4) = buff[ijkl];
              // (12 | 43)
              eri4I(bf1, bf2, bf4, bf3) = buff[ijkl];
              // (21 | 34)
              eri4I(bf2, bf1, bf3, bf4) = buff[ijkl];
              // (21 | 43)
              eri4I(bf2, bf1, bf4, bf3) = buff[ijkl];
              // (34 | 12)
              eri4I(bf3, bf4, bf1, bf2) = buff[ijkl];
              // (43 | 12)
              eri4I(bf4, bf3, bf1, bf2) = buff[ijkl];
              // (34 | 21)
              eri4I(bf3, bf4, bf2, bf1) = buff[ijkl];
              // (43 | 21)
              eri4I(bf4, bf3, bf2, bf1) = buff[ijkl];


          }; // ijkl loop

        } else {

          size_t P = s1, Q = s2, R = s3, S = s4;

          size_t pContrSize = basisSet.shells[P].contr.size();
          size_t qContrSize = basisSet.shells[Q].contr.size();
          size_t rContrSize = basisSet.shells[R].contr.size();
          size_t sContrSize = basisSet.shells[S].contr.size();

          size_t pAMSize = shellPrims[P][0].size();
          size_t qAMSize = shellPrims[Q][0].size();
          size_t rAMSize = shellPrims[R][0].size();
          size_t sAMSize = shellPrims[S][0].size();

          size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

          const double *resPQRS;
          libintGeneralContractionERI(
              P, Q, R, S,
              basisSet, engines[thread_id],
              shellPrims, coefBlocks, workBlocks[thread_id],
              resPQRS);

  #ifdef __DEBUGERI__
          prettyPrintSmart(std::cout, "resPQRS",
                           resPQRS, pqrsAMSize, pContrSize * qContrSize * rContrSize * sContrSize, pqrsAMSize);
  #endif

          // Reorganize
          for (size_t SS = 0, bf4 = bf4_s, pqrs = 0; SS < sContrSize; SS++)
            for (size_t s = 0; s < sAMSize; s++, bf4++)
              for (size_t RR = 0, bf3 = bf3_s; RR < rContrSize; RR++)
                for (size_t r = 0; r < rAMSize; r++, bf3++)
                  for (size_t QQ = 0, bf2 = bf2_s; QQ < qContrSize; QQ++)
                    for (size_t q = 0; q < qAMSize; q++, bf2++)
                      for (size_t PP = 0, bf1 = bf1_s; PP < pContrSize; PP++)
                        for (size_t p = 0; p < pAMSize; p++, bf1++, pqrs++) {
                          double eriVal = resPQRS[(s + sAMSize * (r + rAMSize * (q + qAMSize * p)))
                              + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                          // (12 | 34)
                          eri4I(bf1, bf2, bf3, bf4) = eriVal;
                          // (12 | 43)
                          eri4I(bf1, bf2, bf4, bf3) = eriVal;
                          // (21 | 34)
                          eri4I(bf2, bf1, bf3, bf4) = eriVal;
                          // (21 | 43)
                          eri4I(bf2, bf1, bf4, bf3) = eriVal;
                          // (34 | 12)
                          eri4I(bf3, bf4, bf1, bf2) = eriVal;
                          // (43 | 12)
                          eri4I(bf4, bf3, bf1, bf2) = eriVal;
                          // (34 | 21)
                          eri4I(bf3, bf4, bf2, bf1) = eriVal;
                          // (43 | 21)
                          eri4I(bf4, bf3, bf2, bf1) = eriVal;
                        }

        }

      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    auto durERI4 = tock(topERI4);
    //std::cout << "  Libint-ERI4 duration   = " << durERI4 << std::endl;

    for (double *p : coefBlocks) {
      if (p) CQMemManager::get().free(p);
    }
    coefBlocks.clear();

    if (not segmented)
      for (size_t i = 0; i < nthreads; i++) {
        CQMemManager::get().free(workBlocks[i]);
      }

    // Debug output of the ERIs
#ifdef __DEBUGERI__
    std::cout << "Two-Electron Integrals (ERIs)" << std::endl;
    std::cout << std::setprecision(12);
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << eri4I(i, j, k, l) << std::endl;
    };
#endif
  }; // InCore4indexTPI<double>::computeERIGCNR


  template <>
  void InCore4indexRelERI<dcomplex>::computeERIDCB(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void InCore4indexRelERI<double>::computeERIDCB(BasisSet &basisSet_, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions &hamiltonianOptions) {

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
 
    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      basisSet_.maxPrim,basisSet_.maxL,2);
    engines[0].set_precision(0.);


    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;

    // There are 78 2nd Derivatives per ERI 
    int AxBx = 3;
    int AxBy = 4;
    int AxBz = 5;
    int AyBx = 14;
    int AyBy = 15;
    int AyBz = 16;
    int AzBx = 24;
    int AzBy = 25;
    int AzBz = 26;

    int CxDx = 60;
    int CxDy = 61;
    int CxDz = 62;
    int CyDx = 65;
    int CyDy = 66;
    int CyDz = 67;
    int CzDx = 69;
    int CzDy = 70;
    int CzDz = 71;

    int AxCx = 6;
    int AxCy = 7;
    int AxCz = 8;
    int AyCx = 17;
    int AyCy = 18;
    int AyCz = 19;
    int AzCx = 27;
    int AzCy = 28;
    int AzCz = 29;

    int AxDx = 9;
    int AxDy = 10;
    int AxDz = 11;
    int AyDx = 20;
    int AyDy = 21;
    int AyDz = 22;
    int AzDx = 30;
    int AzDy = 31;
    int AzDz = 32;

    int BxCx = 36;
    int BxCy = 37;
    int BxCz = 38;
    int ByCx = 44;
    int ByCy = 45;
    int ByCz = 46;
    int BzCx = 51;
    int BzCy = 52;
    int BzCz = 53;

    int BxDx = 39;
    int BxDy = 40;
    int BxDz = 41;
    int ByDx = 47;
    int ByDy = 48;
    int ByDz = 49;
    int BzDx = 54;
    int BzDy = 55;
    int BzDz = 56;

    for (InCore4indexTPI<double>& c : components_)
      c.clear();

    auto topERIDCB = tick();




    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet_.nShell; 
          bf1_s+=n1, s1++) { 

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        // Evaluate ERI for shell quartet
        engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 2>(
          basisSet_.shells[s1],
          basisSet_.shells[s2],
          basisSet_.shells[s3],
          basisSet_.shells[s4]
        );
        const auto& buff = engines[thread_id].results();

        // Place shell quartet into persistent storage with
        // permutational symmetry
        for(i = 0ul, bf1 = bf1_s, ijkl = 0ul ; i < n1; ++i, bf1++) 
        for(j = 0ul, bf2 = bf2_s             ; j < n2; ++j, bf2++) 
        for(k = 0ul, bf3 = bf3_s             ; k < n3; ++k, bf3++) 
        for(l = 0ul, bf4 = bf4_s             ; l < n4; ++l, bf4++, ++ijkl) {

#ifdef __DEBUGERI__

          std::cout << std::scientific << std::setprecision(16); 
#if 0
	  std::cout <<"Libint ∇A∙∇B(ij|kl)"<<std::endl;
	  std::cout<<buff[AxBx][ijkl]<<std::endl; 
	  std::cout<<buff[AxBy][ijkl]<<std::endl; 
	  std::cout<<buff[AxBz][ijkl]<<std::endl;
	  std::cout<<buff[AyBx][ijkl]<<std::endl;
	  std::cout<<buff[AyBy][ijkl]<<std::endl;
	  std::cout<<buff[AyBz][ijkl]<<std::endl;
	  std::cout<<buff[AzBx][ijkl]<<std::endl;
	  std::cout<<buff[AzBy][ijkl]<<std::endl;
	  std::cout<<buff[AzBz][ijkl]<<std::endl;
#endif

#if 0
	  std::cout <<"Libint ∇A∙∇C(ij|kl)"<<std::endl;
	  std::cout<<buff[AxCx][ijkl]<<std::endl; 
	  std::cout<<buff[AxCy][ijkl]<<std::endl; 
	  std::cout<<buff[AxCz][ijkl]<<std::endl;
	  std::cout<<buff[AyCx][ijkl]<<std::endl;
	  std::cout<<buff[AyCy][ijkl]<<std::endl;
	  std::cout<<buff[AyCz][ijkl]<<std::endl;
	  std::cout<<buff[AzCx][ijkl]<<std::endl;
	  std::cout<<buff[AzCy][ijkl]<<std::endl;
	  std::cout<<buff[AzCz][ijkl]<<std::endl;
#endif

#if 0
	  std::cout <<"Libint ∇A∙∇D(ij|kl)"<<std::endl;
	  std::cout<<buff[AxDx][ijkl]<<std::endl; 
	  std::cout<<buff[AxDy][ijkl]<<std::endl; 
	  std::cout<<buff[AxDz][ijkl]<<std::endl;
	  std::cout<<buff[AyDx][ijkl]<<std::endl;
	  std::cout<<buff[AyDy][ijkl]<<std::endl;
	  std::cout<<buff[AyDz][ijkl]<<std::endl;
	  std::cout<<buff[AzDx][ijkl]<<std::endl;
	  std::cout<<buff[AzDy][ijkl]<<std::endl;
	  std::cout<<buff[AzDz][ijkl]<<std::endl;
#endif

#if 0
	  std::cout <<"Libint ∇B∙∇C(ij|kl)"<<std::endl;
	  std::cout<<buff[BxCx][ijkl]<<std::endl; 
	  std::cout<<buff[BxCy][ijkl]<<std::endl; 
	  std::cout<<buff[BxCz][ijkl]<<std::endl;
	  std::cout<<buff[ByCx][ijkl]<<std::endl;
	  std::cout<<buff[ByCy][ijkl]<<std::endl;
	  std::cout<<buff[ByCz][ijkl]<<std::endl;
	  std::cout<<buff[BzCx][ijkl]<<std::endl;
	  std::cout<<buff[BzCy][ijkl]<<std::endl;
	  std::cout<<buff[BzCz][ijkl]<<std::endl;
#endif


#if 0
	  std::cout <<"Libint ∇B∙∇D(ij|kl)"<<std::endl;
	  std::cout<<buff[BxDx][ijkl]<<std::endl; 
	  std::cout<<buff[BxDy][ijkl]<<std::endl; 
	  std::cout<<buff[BxDz][ijkl]<<std::endl;
	  std::cout<<buff[ByDx][ijkl]<<std::endl;
	  std::cout<<buff[ByDy][ijkl]<<std::endl;
	  std::cout<<buff[ByDz][ijkl]<<std::endl;
	  std::cout<<buff[BzDx][ijkl]<<std::endl;
	  std::cout<<buff[BzDy][ijkl]<<std::endl;
	  std::cout<<buff[BzDz][ijkl]<<std::endl;
#endif


#if 0
	  std::cout <<"Libint ∇C∙∇D(ij|kl)"<<std::endl;
	  std::cout<<buff[CxDx][ijkl]<<std::endl; 
	  std::cout<<buff[CxDy][ijkl]<<std::endl; 
	  std::cout<<buff[CxDz][ijkl]<<std::endl;
	  std::cout<<buff[CyDx][ijkl]<<std::endl;
	  std::cout<<buff[CyDy][ijkl]<<std::endl;
	  std::cout<<buff[CyDz][ijkl]<<std::endl;
	  std::cout<<buff[CzDx][ijkl]<<std::endl;
	  std::cout<<buff[CzDy][ijkl]<<std::endl;
	  std::cout<<buff[CzDz][ijkl]<<std::endl;
#endif

#endif //__DEBUGERI__

          auto IJKL = bf1 + bf2*NB + bf3*NB2 + bf4*NB3;
          auto IJLK = bf1 + bf2*NB + bf4*NB2 + bf3*NB3;
          auto JIKL = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
          auto JILK = bf2 + bf1*NB + bf4*NB2 + bf3*NB3;
          auto KLIJ = bf3 + bf4*NB + bf1*NB2 + bf2*NB3;
          auto LKIJ = bf4 + bf3*NB + bf1*NB2 + bf2*NB3;
          auto KLJI = bf3 + bf4*NB + bf2*NB2 + bf1*NB3;
          auto LKJI = bf4 + bf3*NB + bf2*NB2 + bf1*NB3;
#if 1
          /* Coulomb */
          // ∇A∙∇B(ij|kl)
          auto dAdotdB = buff[AxBx][ijkl] + buff[AyBy][ijkl] + buff[AzBz][ijkl];
          // ∇Ax∇B(ijkl)
          auto dAcrossdB_x =  buff[AyBz][ijkl] - buff[AzBy][ijkl];
          auto dAcrossdB_y = -buff[AxBz][ijkl] + buff[AzBx][ijkl];
          auto dAcrossdB_z =  buff[AxBy][ijkl] - buff[AyBx][ijkl];

          // ∇C∙∇D(ij|kl)
          auto dCdotdD = buff[CxDx][ijkl] + buff[CyDy][ijkl] + buff[CzDz][ijkl];
          // ∇Cx∇D(ijkl)
          auto dCcrossdD_x =  buff[CyDz][ijkl] - buff[CzDy][ijkl];
          auto dCcrossdD_y = -buff[CxDz][ijkl] + buff[CzDx][ijkl];
          auto dCcrossdD_z =  buff[CxDy][ijkl] - buff[CyDx][ijkl];

          // ∇A∙∇B(ij|kl) followed by ∇Ax∇B(ij|kl) X, Y, and Z
          // (ij|kl)
          (*this)[0].pointer()[IJKL] = dAdotdB;
          (*this)[1].pointer()[IJKL] = dAcrossdB_x;
          (*this)[2].pointer()[IJKL] = dAcrossdB_y;
          (*this)[3].pointer()[IJKL] = dAcrossdB_z;
          // (ij|lk)
          (*this)[0].pointer()[IJLK] = dAdotdB;
          (*this)[1].pointer()[IJLK] = dAcrossdB_x;
          (*this)[2].pointer()[IJLK] = dAcrossdB_y;
          (*this)[3].pointer()[IJLK] = dAcrossdB_z;
          // (ji|kl)
          (*this)[0].pointer()[JIKL] = dAdotdB;
          (*this)[1].pointer()[JIKL] = -dAcrossdB_x;
          (*this)[2].pointer()[JIKL] = -dAcrossdB_y;
          (*this)[3].pointer()[JIKL] = -dAcrossdB_z;
          // (ji|lk)
          (*this)[0].pointer()[JILK] = dAdotdB;
          (*this)[1].pointer()[JILK] = -dAcrossdB_x;
          (*this)[2].pointer()[JILK] = -dAcrossdB_y;
          (*this)[3].pointer()[JILK] = -dAcrossdB_z;
          // (kl|ij)
          (*this)[0].pointer()[KLIJ] = dCdotdD;
          (*this)[1].pointer()[KLIJ] = dCcrossdD_x;
          (*this)[2].pointer()[KLIJ] = dCcrossdD_y;
          (*this)[3].pointer()[KLIJ] = dCcrossdD_z;
          // (lk|ij)
          (*this)[0].pointer()[LKIJ] = dCdotdD;
          (*this)[1].pointer()[LKIJ] = -dCcrossdD_x;
          (*this)[2].pointer()[LKIJ] = -dCcrossdD_y;
          (*this)[3].pointer()[LKIJ] = -dCcrossdD_z;
          // (kl|ji)
          (*this)[0].pointer()[KLJI] = dCdotdD;
          (*this)[1].pointer()[KLJI] = dCcrossdD_x;
          (*this)[2].pointer()[KLJI] = dCcrossdD_y;
          (*this)[3].pointer()[KLJI] = dCcrossdD_z;
          // (lk|ji)
          (*this)[0].pointer()[LKJI] = dCdotdD;
          (*this)[1].pointer()[LKJI] = -dCcrossdD_x;
          (*this)[2].pointer()[LKJI] = -dCcrossdD_y;
          (*this)[3].pointer()[LKJI] = -dCcrossdD_z;
#endif

          /* Gaunt */
	  if(hamiltonianOptions.Gaunt) {

            // ∇A∙∇C(ij|kl)
            auto dAdotdC = buff[AxCx][ijkl] + buff[AyCy][ijkl] + buff[AzCz][ijkl];
            // ∇Ax∇C(ijkl)
            auto dAcrossdC_x =  buff[AyCz][ijkl] - buff[AzCy][ijkl];
            auto dAcrossdC_y = -buff[AxCz][ijkl] + buff[AzCx][ijkl];
            auto dAcrossdC_z =  buff[AxCy][ijkl] - buff[AyCx][ijkl];
  
            // ∇A∙∇D(ij|kl)
            auto dAdotdD = buff[AxDx][ijkl] + buff[AyDy][ijkl] + buff[AzDz][ijkl];
            // ∇Ax∇D(ijkl)
            auto dAcrossdD_x =  buff[AyDz][ijkl] - buff[AzDy][ijkl];
            auto dAcrossdD_y = -buff[AxDz][ijkl] + buff[AzDx][ijkl];
            auto dAcrossdD_z =  buff[AxDy][ijkl] - buff[AyDx][ijkl];
  
            // ∇B∙∇C(ij|kl)
            auto dBdotdC = buff[BxCx][ijkl]+ buff[ByCy][ijkl]+ buff[BzCz][ijkl];
            // ∇Bx∇C(ijkl)
            auto dBcrossdC_x =  buff[ByCz][ijkl] - buff[BzCy][ijkl];
            auto dBcrossdC_y = -buff[BxCz][ijkl] + buff[BzCx][ijkl];
            auto dBcrossdC_z =  buff[BxCy][ijkl] - buff[ByCx][ijkl];
  
            // ∇B∙∇D(ij|kl)
            auto dBdotdD = buff[BxDx][ijkl] + buff[ByDy][ijkl] + buff[BzDz][ijkl];
            // ∇Bx∇D(ijkl)
            auto dBcrossdD_x =  buff[ByDz][ijkl] - buff[BzDy][ijkl];
            auto dBcrossdD_y = -buff[BxDz][ijkl] + buff[BzDx][ijkl];
            auto dBcrossdD_z =  buff[BxDy][ijkl] - buff[ByDx][ijkl];
  
            // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
            // (ij|kl)
            (*this)[4].pointer()[IJKL] = dBdotdC;
            (*this)[5].pointer()[IJKL] = dBcrossdC_x;
            (*this)[6].pointer()[IJKL] = dBcrossdC_y;
            (*this)[7].pointer()[IJKL] = dBcrossdC_z;
            // (ji|kl)
            (*this)[4].pointer()[JIKL] = dAdotdC;
            (*this)[5].pointer()[JIKL] = dAcrossdC_x;
            (*this)[6].pointer()[JIKL] = dAcrossdC_y;
            (*this)[7].pointer()[JIKL] = dAcrossdC_z;
            // (ij|lk)
            (*this)[4].pointer()[IJLK] = dBdotdD;
            (*this)[5].pointer()[IJLK] = dBcrossdD_x;
            (*this)[6].pointer()[IJLK] = dBcrossdD_y;
            (*this)[7].pointer()[IJLK] = dBcrossdD_z;
            // (ji|lk)
            (*this)[4].pointer()[JILK] = dAdotdD;
            (*this)[5].pointer()[JILK] = dAcrossdD_x;
            (*this)[6].pointer()[JILK] = dAcrossdD_y;
            (*this)[7].pointer()[JILK] = dAcrossdD_z;
            // (kl|ij)
            (*this)[4].pointer()[KLIJ] = dAdotdD;
            (*this)[5].pointer()[KLIJ] = -dAcrossdD_x;
            (*this)[6].pointer()[KLIJ] = -dAcrossdD_y;
            (*this)[7].pointer()[KLIJ] = -dAcrossdD_z;
            // (lk|ij)
            (*this)[4].pointer()[LKIJ] = dAdotdC;
            (*this)[5].pointer()[LKIJ] = -dAcrossdC_x;
            (*this)[6].pointer()[LKIJ] = -dAcrossdC_y;
            (*this)[7].pointer()[LKIJ] = -dAcrossdC_z;
            // (kl|ji)
            (*this)[4].pointer()[KLJI] = dBdotdD;
            (*this)[5].pointer()[KLJI] = -dBcrossdD_x;
            (*this)[6].pointer()[KLJI] = -dBcrossdD_y;
            (*this)[7].pointer()[KLJI] = -dBcrossdD_z;
            // (lk|ji)
            (*this)[4].pointer()[LKJI] = dBdotdC;
            (*this)[5].pointer()[LKJI] = -dBcrossdC_x;
            (*this)[6].pointer()[LKJI] = -dBcrossdC_y;
            (*this)[7].pointer()[LKJI] = -dBcrossdC_z;
  
  
            // ∇B_x∇C_y(ij|kl) + ∇B_y∇C_x(ij|kl)
            // (ij|kl)
            (*this)[8].pointer()[IJKL] = buff[BxCy][ijkl] + buff[ByCx][ijkl];
            // (ji|kl)
            (*this)[8].pointer()[JIKL] = buff[AxCy][ijkl] + buff[AyCx][ijkl];
            // (ij|lk)
            (*this)[8].pointer()[IJLK] = buff[BxDy][ijkl] + buff[ByDx][ijkl];
            // (ji|lk)
            (*this)[8].pointer()[JILK] = buff[AxDy][ijkl] + buff[AyDx][ijkl];
            // (kl|ij)
            (*this)[8].pointer()[KLIJ] = buff[AxDy][ijkl] + buff[AyDx][ijkl];
            // (kl|ji)
            (*this)[8].pointer()[KLJI] = buff[BxDy][ijkl] + buff[ByDx][ijkl];
            // (lk|ij)
            (*this)[8].pointer()[LKIJ] = buff[AxCy][ijkl] + buff[AyCx][ijkl];
            // (lk|ji)
            (*this)[8].pointer()[LKJI] = buff[BxCy][ijkl] + buff[ByCx][ijkl];
  
  
            // ∇B_y∇C_x(ij|kl)
            // (ij|kl)
            (*this)[9].pointer()[IJKL] = buff[ByCx][ijkl];
            // (ji|kl)
            (*this)[9].pointer()[JIKL] = buff[AyCx][ijkl];
            // (ij|lk)
            (*this)[9].pointer()[IJLK] = buff[ByDx][ijkl];
            // (ji|lk)
            (*this)[9].pointer()[JILK] = buff[AyDx][ijkl];
            // (kl|ij)
            (*this)[9].pointer()[KLIJ] = buff[AxDy][ijkl];
            // (kl|ji)
            (*this)[9].pointer()[KLJI] = buff[BxDy][ijkl];
            // (lk|ij)
            (*this)[9].pointer()[LKIJ] = buff[AxCy][ijkl];
            // (lk|ji)
            (*this)[9].pointer()[LKJI] = buff[BxCy][ijkl];
  
  
            // ∇B_x∇C_z(ij|kl) + ∇B_z∇C_x(ij|kl)
            // (ij|kl)
            (*this)[10].pointer()[IJKL] = buff[BxCz][ijkl] + buff[BzCx][ijkl];
            // (ji|kl)
            (*this)[10].pointer()[JIKL] = buff[AxCz][ijkl] + buff[AzCx][ijkl];
            // (ij|lk)
            (*this)[10].pointer()[IJLK] = buff[BxDz][ijkl] + buff[BzDx][ijkl];
            // (ji|lk)
            (*this)[10].pointer()[JILK] = buff[AxDz][ijkl] + buff[AzDx][ijkl];
            // (kl|ij)
            (*this)[10].pointer()[KLIJ] = buff[AxDz][ijkl] + buff[AzDx][ijkl];
            // (kl|ji)
            (*this)[10].pointer()[KLJI] = buff[BxDz][ijkl] + buff[BzDx][ijkl];
            // (lk|ij)
            (*this)[10].pointer()[LKIJ] = buff[AxCz][ijkl] + buff[AzCx][ijkl];
            // (lk|ji)
            (*this)[10].pointer()[LKJI] = buff[BxCz][ijkl] + buff[BzCx][ijkl];
  
  
            // ∇B_z∇C_x(ij|kl)
            // (ij|kl)
            (*this)[11].pointer()[IJKL] = buff[BzCx][ijkl];
            // (ji|kl)
            (*this)[11].pointer()[JIKL] = buff[AzCx][ijkl];
            // (ij|lk)
            (*this)[11].pointer()[IJLK] = buff[BzDx][ijkl];
            // (ji|lk)
            (*this)[11].pointer()[JILK] = buff[AzDx][ijkl];
            // (kl|ij)
            (*this)[11].pointer()[KLIJ] = buff[AxDz][ijkl];
            // (kl|ji)
            (*this)[11].pointer()[KLJI] = buff[BxDz][ijkl];
            // (lk|ij)
            (*this)[11].pointer()[LKIJ] = buff[AxCz][ijkl];
            // (lk|ji)
            (*this)[11].pointer()[LKJI] = buff[BxCz][ijkl];
  
  
            // ∇B_y∇C_z(ij|kl) + ∇B_z∇C_y(ij|kl)
            // (ij|kl)
            (*this)[12].pointer()[IJKL] = buff[ByCz][ijkl] + buff[BzCy][ijkl];
            // (ji|kl)
            (*this)[12].pointer()[JIKL] = buff[AyCz][ijkl] + buff[AzCy][ijkl];
            // (ij|lk)
            (*this)[12].pointer()[IJLK] = buff[ByDz][ijkl] + buff[BzDy][ijkl];
            // (ji|lk)
            (*this)[12].pointer()[JILK] = buff[AyDz][ijkl] + buff[AzDy][ijkl];
            // (kl|ij)
            (*this)[12].pointer()[KLIJ] = buff[AyDz][ijkl] + buff[AzDy][ijkl];
            // (kl|ji)
            (*this)[12].pointer()[KLJI] = buff[ByDz][ijkl] + buff[BzDy][ijkl];
            // (lk|ij)
            (*this)[12].pointer()[LKIJ] = buff[AyCz][ijkl] + buff[AzCy][ijkl];
            // (lk|ji)
            (*this)[12].pointer()[LKJI] = buff[ByCz][ijkl] + buff[BzCy][ijkl];
  
  
            // ∇B_z∇C_y(ij|kl)
            // (ij|kl)
            (*this)[13].pointer()[IJKL] = buff[BzCy][ijkl];
            // (ji|kl)
            (*this)[13].pointer()[JIKL] = buff[AzCy][ijkl];
            // (ij|lk)
            (*this)[13].pointer()[IJLK] = buff[BzDy][ijkl];
            // (ji|lk)
            (*this)[13].pointer()[JILK] = buff[AzDy][ijkl];
            // (kl|ij)
            (*this)[13].pointer()[KLIJ] = buff[AyDz][ijkl];
            // (kl|ji)
            (*this)[13].pointer()[KLJI] = buff[ByDz][ijkl];
            // (lk|ij)
            (*this)[13].pointer()[LKIJ] = buff[AyCz][ijkl];
            // (lk|ji)
            (*this)[13].pointer()[LKJI] = buff[ByCz][ijkl];
  
  
            // - ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) + ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[14].pointer()[IJKL] = - buff[BxCx][ijkl] - buff[ByCy][ijkl] + buff[BzCz][ijkl];
            // (ji|kl)
            (*this)[14].pointer()[JIKL] = - buff[AxCx][ijkl] - buff[AyCy][ijkl] + buff[AzCz][ijkl];
            // (ij|lk)
            (*this)[14].pointer()[IJLK] = - buff[BxDx][ijkl] - buff[ByDy][ijkl] + buff[BzDz][ijkl];
            // (ji|lk)
            (*this)[14].pointer()[JILK] = - buff[AxDx][ijkl] - buff[AyDy][ijkl] + buff[AzDz][ijkl];
            // (kl|ij)
            (*this)[14].pointer()[KLIJ] = - buff[AxDx][ijkl] - buff[AyDy][ijkl] + buff[AzDz][ijkl];
            // (kl|ji)
            (*this)[14].pointer()[KLJI] = - buff[BxDx][ijkl] - buff[ByDy][ijkl] + buff[BzDz][ijkl];
            // (lk|ij)
            (*this)[14].pointer()[LKIJ] = - buff[AxCx][ijkl] - buff[AyCy][ijkl] + buff[AzCz][ijkl];
            // (lk|ji)
            (*this)[14].pointer()[LKJI] = - buff[BxCx][ijkl] - buff[ByCy][ijkl] + buff[BzCz][ijkl];
  
  
            // ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[15].pointer()[IJKL] = buff[BxCx][ijkl] - buff[ByCy][ijkl] - buff[BzCz][ijkl];
            // (ji|kl)
            (*this)[15].pointer()[JIKL] = buff[AxCx][ijkl] - buff[AyCy][ijkl] - buff[AzCz][ijkl];
            // (ij|lk)
            (*this)[15].pointer()[IJLK] = buff[BxDx][ijkl] - buff[ByDy][ijkl] - buff[BzDz][ijkl];
            // (ji|lk)
            (*this)[15].pointer()[JILK] = buff[AxDx][ijkl] - buff[AyDy][ijkl] - buff[AzDz][ijkl];
            // (kl|ij)
            (*this)[15].pointer()[KLIJ] = buff[AxDx][ijkl] - buff[AyDy][ijkl] - buff[AzDz][ijkl];
            // (kl|ji)
            (*this)[15].pointer()[KLJI] = buff[BxDx][ijkl] - buff[ByDy][ijkl] - buff[BzDz][ijkl];
            // (lk|ij)
            (*this)[15].pointer()[LKIJ] = buff[AxCx][ijkl] - buff[AyCy][ijkl] - buff[AzCz][ijkl];
            // (lk|ji)
            (*this)[15].pointer()[LKJI] = buff[BxCx][ijkl] - buff[ByCy][ijkl] - buff[BzCz][ijkl];
  
  
            // - ∇B_x∇C_x(ij|kl) + ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[16].pointer()[IJKL] = - buff[BxCx][ijkl] + buff[ByCy][ijkl] - buff[BzCz][ijkl];
            // (ji|kl)
            (*this)[16].pointer()[JIKL] = - buff[AxCx][ijkl] + buff[AyCy][ijkl] - buff[AzCz][ijkl];
            // (ij|lk)
            (*this)[16].pointer()[IJLK] = - buff[BxDx][ijkl] + buff[ByDy][ijkl] - buff[BzDz][ijkl];
            // (ji|lk)
            (*this)[16].pointer()[JILK] = - buff[AxDx][ijkl] + buff[AyDy][ijkl] - buff[AzDz][ijkl];
            // (kl|ij)
            (*this)[16].pointer()[KLIJ] = - buff[AxDx][ijkl] + buff[AyDy][ijkl] - buff[AzDz][ijkl];
            // (kl|ji)
            (*this)[16].pointer()[KLJI] = - buff[BxDx][ijkl] + buff[ByDy][ijkl] - buff[BzDz][ijkl];
            // (lk|ij)
            (*this)[16].pointer()[LKIJ] = - buff[AxCx][ijkl] + buff[AyCy][ijkl] - buff[AzCz][ijkl];
            // (lk|ji)
            (*this)[16].pointer()[LKJI] = - buff[BxCx][ijkl] + buff[ByCy][ijkl] - buff[BzCz][ijkl];
  
  
            // ∇B_x∇C_x(ij|kl)
            // (ij|kl)
            (*this)[17].pointer()[IJKL] = buff[BxCx][ijkl];
            // (ji|kl)
            (*this)[17].pointer()[JIKL] = buff[AxCx][ijkl];
            // (ij|lk)
            (*this)[17].pointer()[IJLK] = buff[BxDx][ijkl];
            // (ji|lk)
            (*this)[17].pointer()[JILK] = buff[AxDx][ijkl];
            // (kl|ij)
            (*this)[17].pointer()[KLIJ] = buff[AxDx][ijkl];
            // (kl|ji)
            (*this)[17].pointer()[KLJI] = buff[BxDx][ijkl];
            // (lk|ij)
            (*this)[17].pointer()[LKIJ] = buff[AxCx][ijkl];
            // (lk|ji)
            (*this)[17].pointer()[LKJI] = buff[BxCx][ijkl];
  
  
            // ∇B_x∇C_y(ij|kl)
            // (ij|kl)
            (*this)[18].pointer()[IJKL] = buff[BxCy][ijkl];
            // (ji|kl)
            (*this)[18].pointer()[JIKL] = buff[AxCy][ijkl];
            // (ij|lk)
            (*this)[18].pointer()[IJLK] = buff[BxDy][ijkl];
            // (ji|lk)
            (*this)[18].pointer()[JILK] = buff[AxDy][ijkl];
            // (kl|ij)
            (*this)[18].pointer()[KLIJ] = buff[AyDx][ijkl];
            // (kl|ji)
            (*this)[18].pointer()[KLJI] = buff[ByDx][ijkl];
            // (lk|ij)
            (*this)[18].pointer()[LKIJ] = buff[AyCx][ijkl];
            // (lk|ji)
            (*this)[18].pointer()[LKJI] = buff[ByCx][ijkl];
  
  
            // ∇B_x∇C_z(ij|kl)
            // (ij|kl)
            (*this)[19].pointer()[IJKL] = buff[BxCz][ijkl];
            // (ji|kl)
            (*this)[19].pointer()[JIKL] = buff[AxCz][ijkl];
            // (ij|lk)
            (*this)[19].pointer()[IJLK] = buff[BxDz][ijkl];
            // (ji|lk)
            (*this)[19].pointer()[JILK] = buff[AxDz][ijkl];
            // (kl|ij)
            (*this)[19].pointer()[KLIJ] = buff[AzDx][ijkl];
            // (kl|ji)
            (*this)[19].pointer()[KLJI] = buff[BzDx][ijkl];
            // (lk|ij)
            (*this)[19].pointer()[LKIJ] = buff[AzCx][ijkl];
            // (lk|ji)
            (*this)[19].pointer()[LKJI] = buff[BzCx][ijkl];
  
  
  	  // ∇B_y∇C_y(ij|kl)
            // (ij|kl)
            (*this)[20].pointer()[IJKL] = buff[ByCy][ijkl];
            // (ji|kl)
            (*this)[20].pointer()[JIKL] = buff[AyCy][ijkl];
            // (ij|lk)
            (*this)[20].pointer()[IJLK] = buff[ByDy][ijkl];
            // (ji|lk)
            (*this)[20].pointer()[JILK] = buff[AyDy][ijkl];
            // (kl|ij)
            (*this)[20].pointer()[KLIJ] = buff[AyDy][ijkl];
            // (kl|ji)
            (*this)[20].pointer()[KLJI] = buff[ByDy][ijkl];
            // (lk|ij)
            (*this)[20].pointer()[LKIJ] = buff[AyCy][ijkl];
            // (lk|ji)
            (*this)[20].pointer()[LKJI] = buff[ByCy][ijkl];
  
  
  	  // ∇B_y∇C_z(ij|kl)
            // (ij|kl)
            (*this)[21].pointer()[IJKL] = buff[ByCz][ijkl];
            // (ji|kl)
            (*this)[21].pointer()[JIKL] = buff[AyCz][ijkl];
            // (ij|lk)
            (*this)[21].pointer()[IJLK] = buff[ByDz][ijkl];
            // (ji|lk)
            (*this)[21].pointer()[JILK] = buff[AyDz][ijkl];
            // (kl|ij)
            (*this)[21].pointer()[KLIJ] = buff[AzDy][ijkl];
            // (kl|ji)
            (*this)[21].pointer()[KLJI] = buff[BzDy][ijkl];
            // (lk|ij)
            (*this)[21].pointer()[LKIJ] = buff[AzCy][ijkl];
            // (lk|ji)
            (*this)[21].pointer()[LKJI] = buff[BzCy][ijkl];
  
  
            // ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[22].pointer()[IJKL] = buff[BzCz][ijkl];
            // (ji|kl)
            (*this)[22].pointer()[JIKL] = buff[AzCz][ijkl];
            // (ij|lk)
            (*this)[22].pointer()[IJLK] = buff[BzDz][ijkl];
            // (ji|lk)
            (*this)[22].pointer()[JILK] = buff[AzDz][ijkl];
            // (kl|ij)
            (*this)[22].pointer()[KLIJ] = buff[AzDz][ijkl];
            // (kl|ji)
            (*this)[22].pointer()[KLJI] = buff[BzDz][ijkl];
            // (lk|ij)
            (*this)[22].pointer()[LKIJ] = buff[AzCz][ijkl];
            // (lk|ji)
            (*this)[22].pointer()[LKJI] = buff[BzCz][ijkl];
  
  	  } // Gaunt
        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    auto durERIDCB = tock(topERIDCB);
    std::cout << "  Libint-ERI-Dirac-Coulomb-Breit duration   = " << durERIDCB << std::endl;


#ifdef __DEBUGERI__

    std::cout << std::scientific << std::setprecision(16);

    std::cout << "ERI (ab|cd)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)(i, j, k, l) << std::endl;
    };

    std::cout << "ERI00-03: ∇A∙∇B(ab|cd)  ∇Ax∇B(ab|cd)-X  ∇Ax∇B(ab|cd)-Y  ∇Ax∇B(ab|cd)-Z" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[0](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[1](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[2](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[3](i, j, k, l) << std::endl;
    };

    std::cout << "ERI04-07: ∇B∙∇C(ab|cd)  ∇Bx∇C(ab|cd)-X  ∇Bx∇C(ab|cd)-Y  ∇Bx∇C(ab|cd)-Z" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[4](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[5](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[6](i, j, k, l);
      std::cout << "   ";
      std::cout << (*this)[7](i, j, k, l) << std::endl;
    };

    std::cout << "ERI08: ∇B_x∇C_y(ij|kl) + ∇B_y∇C_x(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[8](i, j, k, l) << std::endl;
    };

    std::cout << "ERI09: ∇B_y∇C_x(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[9](i, j, k, l) << std::endl;
    };

    std::cout << "ERI10: ∇B_x∇C_z(ij|kl) + ∇B_z∇C_x(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[10](i, j, k, l) << std::endl;
    };

    std::cout << "ERI11: ∇B_z∇C_x(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[11](i, j, k, l) << std::endl;
    };

    std::cout << "ERI12: ∇B_y∇C_z(ij|kl) + ∇B_z∇C_y(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[12](i, j, k, l) << std::endl;
    };

    std::cout << "ERI13: ∇B_z∇C_y(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[13](i, j, k, l) << std::endl;
    };

    std::cout << "ERI14: - ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) + ∇B_z∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[14](i, j, k, l) << std::endl;
    };

    std::cout << "ERI15: ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[15](i, j, k, l) << std::endl;
    };

    std::cout << "ERI16: - ∇B_x∇C_x(ij|kl) + ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[16](i, j, k, l) << std::endl;
    };

    std::cout << "ERI17: ∇B_x∇C_x(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[17](i, j, k, l) << std::endl;
    };

    std::cout << "ERI18: ∇B_x∇C_y(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[18](i, j, k, l) << std::endl;
    };

    std::cout << "ERI19: ∇B_x∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[19](i, j, k, l) << std::endl;
    };

    std::cout << "ERI20: ∇B_y∇C_y(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[20](i, j, k, l) << std::endl;
    };

    std::cout << "ERI21: ∇B_y∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[21](i, j, k, l) << std::endl;
    };

    std::cout << "ERI22: ∇B_z∇C_z(ij|kl)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)[22](i, j, k, l) << std::endl;
    };

#if 0
    prettyPrintSmart(std::cout,"Rank-2 ERI00 ∇∇(ab|cd)",(*this)[0].pointer(),    NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI01 ∇∇(ab|cd)",(*this)[1].pointer(),       NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI02 ∇∇(ab|cd)",(*this)[2].pointer(),       NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI03 ∇∇(ab|cd)",(*this)[3].pointer(),       NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI04 ∇∇(ab|cd)",(*this)[4].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI05 ∇∇(ab|cd)",(*this)[5].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI06 ∇∇(ab|cd)",(*this)[6].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI07 ∇∇(ab|cd)",(*this)[7].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI08 ∇∇(ab|cd)",(*this)[8].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI09 ∇∇(ab|cd)",(*this)[9].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI10 ∇∇(ab|cd)",(*this)[10].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI11 ∇∇(ab|cd)",(*this)[11].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI12 ∇∇(ab|cd)",(*this)[12].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI13 ∇∇(ab|cd)",(*this)[13].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI14 ∇∇(ab|cd)",(*this)[14].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI15 ∇∇(ab|cd)",(*this)[15].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI16 ∇∇(ab|cd)",(*this)[16].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI17 ∇∇(ab|cd)",(*this)[17].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI18 ∇∇(ab|cd)",(*this)[18].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI19 ∇∇(ab|cd)",(*this)[19].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI20 ∇∇(ab|cd)",(*this)[20].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI21 ∇∇(ab|cd)",(*this)[21].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI22 ∇∇(ab|cd)",(*this)[22].pointer(), NB*NB,NB*NB,NB*NB);
#endif

#endif

  }  // computeERIDCB


  // SS start: compute ERI Gauge integral


  template <>
  void InCore4indexRelERI<dcomplex>::computeERIGauge(BasisSet&, Molecule&,
                                                     EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
      CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void InCore4indexRelERI<double>::computeERIGauge(BasisSet &basisSet_, Molecule &,
                                                   EMPerturbation &, OPERATOR,
                                                   const HamiltonianOptions &hamiltonianOptions) {

    std::cout << "in house gauge integral" << std::endl;

    auto nERIRef = 4; // Dirac-Coulomb
    if (hamiltonianOptions.Gaunt) nERIRef += 19; // Gaunt
    if (hamiltonianOptions.DiracCoulombSSSS) nERIRef += 16; // Dirac-Coulomb-SSSS

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Allocate and zero out ERIs
    size_t NB = basisSet_.nBasis;
    size_t NB2 = NB * NB;
    size_t NB3 = NB2 * NB;
    size_t NB4 = NB2 * NB2;


//SS start
//    double *gaugeSLSLACxx = CQMemManager::get().malloc<double>(pow(NB,4)*16); // 9 elements, xx, xy ... component of derivatives, xx of gauge
// this piece is fake, ignore

    double *gaugeLSSLACsymm = CQMemManager::get().malloc<double>(
      pow(NB, 4) * 16); // 9 elements, xx, xy ... component of derivatives, xx of gauge
    double *gaugeSLSLACxxtrue = CQMemManager::get().malloc<double>(
      pow(NB, 4) * 16); // 9 elements, xx, xy ... component of derivatives, xx of gauge
    double *gaugeSLLSsymm = CQMemManager::get().malloc<double>(
      pow(NB, 4) * 16); // 9 elements, xx, xy .. component of derivatives, xx of gauge
    double *gaugeLSLSsymm = CQMemManager::get().malloc<double>(
      pow(NB, 4) * 16); // 9 elements, xx, xy .. component of derivatives, xx of gauge

//SS end


    auto topERIGauge = tick();

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1, n2, n3, n4, i, j, k, l, ijkl, bf1, bf2, bf3, bf4;
      size_t s4_max;

      int cart_i_size, cart_j_size, cart_k_size, cart_l_size;
      int cart_ip, cart_im, cart_kp, cart_km;
      int cart_ip_size, cart_kp_size;
      int cart_im_size, cart_km_size;

      for (size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet_.nShell;
           bf1_s += n1, s1++) {

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

      for (size_t s2(0), bf2_s(0); s2 <= s1; bf2_s += n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for (size_t s3(0), bf3_s(0); s3 <= s1; bf3_s += n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for (size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s += n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
#ifdef _OPENMP
        if (s1234 % nthreads != thread_id) continue;
#endif


        libint2::ShellPair pair1_to_use;
        pair1_to_use.init(basisSet_.shells[s1], basisSet_.shells[s2], -2000);

        libint2::ShellPair pair2_to_use;
        pair2_to_use.init(basisSet_.shells[s3], basisSet_.shells[s4], -2000);


#ifdef __DEBUGGAUGE__

        std::cout<<"LA "<<basisSet_.shells[s1].contr[0].l
                 <<" LB "<<basisSet_.shells[s2].contr[0].l
                 <<" LC "<<basisSet_.shells[s3].contr[0].l
                 <<" LD "<<basisSet_.shells[s4].contr[0].l<<std::endl;
        std::cout<<"s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;

#endif


        auto gaugeERIgradAC_sph = RealGTOIntEngine::ACgaugederiv(pair1_to_use, pair2_to_use,
                                                                 basisSet_.shells[s1],
                                                                 basisSet_.shells[s2],
                                                                 basisSet_.shells[s3],
                                                                 basisSet_.shells[s4]
        );

        // swap AB
        libint2::ShellPair pair1_to_use_BA;
        pair1_to_use_BA.init(basisSet_.shells[s2], basisSet_.shells[s1], -2000);

        auto gaugeERIgradAC_sphswapAB = RealGTOIntEngine::ACgaugederiv(pair1_to_use_BA, pair2_to_use,
                                                                       basisSet_.shells[s2],
                                                                       basisSet_.shells[s1],
                                                                       basisSet_.shells[s3],
                                                                       basisSet_.shells[s4]
        );

        // swap CD
        libint2::ShellPair pair2_to_use_DC;
        pair2_to_use_DC.init(basisSet_.shells[s4], basisSet_.shells[s3], -2000);

        auto gaugeERIgradAC_sphswapCD = RealGTOIntEngine::ACgaugederiv(pair1_to_use, pair2_to_use_DC,
                                                                       basisSet_.shells[s1],
                                                                       basisSet_.shells[s2],
                                                                       basisSet_.shells[s4],
                                                                       basisSet_.shells[s3]
        );

        // swap AB and CD
        auto gaugeERIgradAC_sphswapABCD = RealGTOIntEngine::ACgaugederiv(pair1_to_use_BA, pair2_to_use_DC,
                                                                         basisSet_.shells[s2],
                                                                         basisSet_.shells[s1],
                                                                         basisSet_.shells[s4],
                                                                         basisSet_.shells[s3]
        );


        for (i = 0ul, bf1 = bf1_s, ijkl = 0ul; i < n1; ++i, bf1++)
        for (j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
        for (k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++)
        for (l = 0ul, bf4 = bf4_s; l < n4; ++l, bf4++, ++ijkl) {

          // AC derivative
          // First dimension: 0-[∇_x∇_x],1-[∇_x∇_y],2-[∇_x∇_z],3-[∇_y∇_x],4-[∇_y∇_y],
    //                  5-[∇_y∇_z],6-[∇_z∇_x],7-[∇_z∇_y],8-[∇_z∇_z]
          // Second dimension: 0-[xx],1-[xy],2-[xz],3-[yy],4-[yz],5-[zz]



// (nabla ij|nabla kl)
// sigma1x sigma2x
          // σ_xσ_x(∇ij|∇kl) = [∇_y∇_y][zz] - [∇_y∇_z][zy] - [∇_z∇_y][yz] + [∇_z∇_z][yy]
          gaugeSLSLACxxtrue[bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[4][5][ijkl] - gaugeERIgradAC_sph[5][4][ijkl]
            - gaugeERIgradAC_sph[7][4][ijkl] + gaugeERIgradAC_sph[8][3][ijkl];

// sigma1x sigma2y
          // σ_xσ_y(∇ij|∇kl) = [∇_y∇_z][zx] - [∇_y∇_x][zz] - [∇_z∇_z][yx] + [∇_z∇_x][yz]
          gaugeSLSLACxxtrue[NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[5][2][ijkl] - gaugeERIgradAC_sph[3][5][ijkl]
            - gaugeERIgradAC_sph[8][1][ijkl] + gaugeERIgradAC_sph[6][4][ijkl];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[3][4][ijkl] - gaugeERIgradAC_sph[4][2][ijkl]
            - gaugeERIgradAC_sph[6][3][ijkl] + gaugeERIgradAC_sph[7][1][ijkl];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[7][2][ijkl] - gaugeERIgradAC_sph[8][1][ijkl]
            - gaugeERIgradAC_sph[1][5][ijkl] + gaugeERIgradAC_sph[2][4][ijkl];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[8][0][ijkl] - gaugeERIgradAC_sph[6][2][ijkl]
            - gaugeERIgradAC_sph[2][2][ijkl] + gaugeERIgradAC_sph[0][5][ijkl];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[6][1][ijkl] - gaugeERIgradAC_sph[7][0][ijkl]
            - gaugeERIgradAC_sph[0][4][ijkl] + gaugeERIgradAC_sph[1][2][ijkl];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[1][4][ijkl] - gaugeERIgradAC_sph[2][3][ijkl]
            - gaugeERIgradAC_sph[4][2][ijkl] + gaugeERIgradAC_sph[5][1][ijkl];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
//gaugeERIgradAC_sph[2][1][ijkl] - gaugeERIgradAC_sph[0][1][ijkl]
            gaugeERIgradAC_sph[2][1][ijkl] - gaugeERIgradAC_sph[0][4][ijkl]
            - gaugeERIgradAC_sph[5][0][ijkl] + gaugeERIgradAC_sph[3][2][ijkl];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sph[0][3][ijkl] - gaugeERIgradAC_sph[1][1][ijkl]
            - gaugeERIgradAC_sph[3][1][ijkl] + gaugeERIgradAC_sph[4][0][ijkl];

// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] = -(
            gaugeERIgradAC_sph[0][0][ijkl] + gaugeERIgradAC_sph[1][1][ijkl]
            + gaugeERIgradAC_sph[2][2][ijkl] + gaugeERIgradAC_sph[3][1][ijkl]
            + gaugeERIgradAC_sph[4][3][ijkl] + gaugeERIgradAC_sph[5][4][ijkl]
            + gaugeERIgradAC_sph[6][2][ijkl] + gaugeERIgradAC_sph[7][4][ijkl]
            + gaugeERIgradAC_sph[8][5][ijkl]);

// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[3][2][ijkl] + gaugeERIgradAC_sph[4][4][ijkl]
              + gaugeERIgradAC_sph[5][5][ijkl] - gaugeERIgradAC_sph[6][1][ijkl]
              - gaugeERIgradAC_sph[7][3][ijkl] - gaugeERIgradAC_sph[8][4][ijkl]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[6][0][ijkl] + gaugeERIgradAC_sph[7][1][ijkl]
              + gaugeERIgradAC_sph[8][2][ijkl] - gaugeERIgradAC_sph[0][2][ijkl]
              - gaugeERIgradAC_sph[1][4][ijkl] - gaugeERIgradAC_sph[2][5][ijkl]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[0][1][ijkl] + gaugeERIgradAC_sph[1][3][ijkl]
              + gaugeERIgradAC_sph[2][4][ijkl] - gaugeERIgradAC_sph[3][0][ijkl]
              - gaugeERIgradAC_sph[4][1][ijkl] - gaugeERIgradAC_sph[5][2][ijkl]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[1][2][ijkl] - gaugeERIgradAC_sph[2][1][ijkl]
              + gaugeERIgradAC_sph[4][4][ijkl] - gaugeERIgradAC_sph[5][3][ijkl]
              + gaugeERIgradAC_sph[7][5][ijkl] - gaugeERIgradAC_sph[8][4][ijkl]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[2][0][ijkl] - gaugeERIgradAC_sph[0][2][ijkl]
              + gaugeERIgradAC_sph[5][1][ijkl] - gaugeERIgradAC_sph[3][4][ijkl]
              + gaugeERIgradAC_sph[8][2][ijkl] - gaugeERIgradAC_sph[6][5][ijkl]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf1 * NB3 + bf2 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sph[0][1][ijkl] - gaugeERIgradAC_sph[1][0][ijkl]
              + gaugeERIgradAC_sph[3][3][ijkl] - gaugeERIgradAC_sph[4][1][ijkl]
              + gaugeERIgradAC_sph[6][4][ijkl] - gaugeERIgradAC_sph[7][2][ijkl]);






// (nabla ji|nabla kl)
          int jikl = j * n4 * n3 * n1 + i * n4 * n3 + k * n4 + l;
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[4][5][jikl] - gaugeERIgradAC_sphswapAB[5][4][jikl]
            - gaugeERIgradAC_sphswapAB[7][4][jikl] + gaugeERIgradAC_sphswapAB[8][3][jikl];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[5][2][jikl] - gaugeERIgradAC_sphswapAB[3][5][jikl]
            - gaugeERIgradAC_sphswapAB[8][1][jikl] + gaugeERIgradAC_sphswapAB[6][4][jikl];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[3][4][jikl] - gaugeERIgradAC_sphswapAB[4][2][jikl]
            - gaugeERIgradAC_sphswapAB[6][3][jikl] + gaugeERIgradAC_sphswapAB[7][1][jikl];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[7][2][jikl] - gaugeERIgradAC_sphswapAB[8][1][jikl]
            - gaugeERIgradAC_sphswapAB[1][5][jikl] + gaugeERIgradAC_sphswapAB[2][4][jikl];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[8][0][jikl] - gaugeERIgradAC_sphswapAB[6][2][jikl]
            - gaugeERIgradAC_sphswapAB[2][2][jikl] + gaugeERIgradAC_sphswapAB[0][5][jikl];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[6][1][jikl] - gaugeERIgradAC_sphswapAB[7][0][jikl]
            - gaugeERIgradAC_sphswapAB[0][4][jikl] + gaugeERIgradAC_sphswapAB[1][2][jikl];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[1][4][jikl] - gaugeERIgradAC_sphswapAB[2][3][jikl]
            - gaugeERIgradAC_sphswapAB[4][2][jikl] + gaugeERIgradAC_sphswapAB[5][1][jikl];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[2][1][jikl] - gaugeERIgradAC_sphswapAB[0][4][jikl]
            - gaugeERIgradAC_sphswapAB[5][0][jikl] + gaugeERIgradAC_sphswapAB[3][2][jikl];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            gaugeERIgradAC_sphswapAB[0][3][jikl] - gaugeERIgradAC_sphswapAB[1][1][jikl]
            - gaugeERIgradAC_sphswapAB[3][1][jikl] + gaugeERIgradAC_sphswapAB[4][0][jikl];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] = -(
            gaugeERIgradAC_sphswapAB[0][0][jikl] + gaugeERIgradAC_sphswapAB[1][1][jikl]
            + gaugeERIgradAC_sphswapAB[2][2][jikl] + gaugeERIgradAC_sphswapAB[3][1][jikl]
            + gaugeERIgradAC_sphswapAB[4][3][jikl] + gaugeERIgradAC_sphswapAB[5][4][jikl]
            + gaugeERIgradAC_sphswapAB[6][2][jikl] + gaugeERIgradAC_sphswapAB[7][4][jikl]
            + gaugeERIgradAC_sphswapAB[8][5][jikl]);


// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[3][2][jikl] + gaugeERIgradAC_sphswapAB[4][4][jikl]
              + gaugeERIgradAC_sphswapAB[5][5][jikl] - gaugeERIgradAC_sphswapAB[6][1][jikl]
              - gaugeERIgradAC_sphswapAB[7][3][jikl] - gaugeERIgradAC_sphswapAB[8][4][jikl]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[6][0][jikl] + gaugeERIgradAC_sphswapAB[7][1][jikl]
              + gaugeERIgradAC_sphswapAB[8][2][jikl] - gaugeERIgradAC_sphswapAB[0][2][jikl]
              - gaugeERIgradAC_sphswapAB[1][4][jikl] - gaugeERIgradAC_sphswapAB[2][5][jikl]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[0][1][jikl] + gaugeERIgradAC_sphswapAB[1][3][jikl]
              + gaugeERIgradAC_sphswapAB[2][4][jikl] - gaugeERIgradAC_sphswapAB[3][0][jikl]
              - gaugeERIgradAC_sphswapAB[4][1][jikl] - gaugeERIgradAC_sphswapAB[5][2][jikl]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[1][2][jikl] - gaugeERIgradAC_sphswapAB[2][1][jikl]
              + gaugeERIgradAC_sphswapAB[4][4][jikl] - gaugeERIgradAC_sphswapAB[5][3][jikl]
              + gaugeERIgradAC_sphswapAB[7][5][jikl] - gaugeERIgradAC_sphswapAB[8][4][jikl]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[2][0][jikl] - gaugeERIgradAC_sphswapAB[0][2][jikl]
              + gaugeERIgradAC_sphswapAB[5][1][jikl] - gaugeERIgradAC_sphswapAB[3][4][jikl]
              + gaugeERIgradAC_sphswapAB[8][2][jikl] - gaugeERIgradAC_sphswapAB[6][5][jikl]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf2 * NB3 + bf1 * NB2 + bf3 * NB + bf4] =
            -(gaugeERIgradAC_sphswapAB[0][1][jikl] - gaugeERIgradAC_sphswapAB[1][0][jikl]
              + gaugeERIgradAC_sphswapAB[3][3][jikl] - gaugeERIgradAC_sphswapAB[4][1][jikl]
              + gaugeERIgradAC_sphswapAB[6][4][jikl] - gaugeERIgradAC_sphswapAB[7][2][jikl]);








// (nabla ij|nabla lk)
          //int ijlk = i*n3*n4*n2 + j*n4*n3 +k*n3+l;
          int ijlk = i * n3 * n4 * n2 + j * n4 * n3 + l * n3 + k;
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[4][5][ijlk] - gaugeERIgradAC_sphswapCD[5][4][ijlk]
            - gaugeERIgradAC_sphswapCD[7][4][ijlk] + gaugeERIgradAC_sphswapCD[8][3][ijlk];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[5][2][ijlk] - gaugeERIgradAC_sphswapCD[3][5][ijlk]
            - gaugeERIgradAC_sphswapCD[8][1][ijlk] + gaugeERIgradAC_sphswapCD[6][4][ijlk];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[3][4][ijlk] - gaugeERIgradAC_sphswapCD[4][2][ijlk]
            - gaugeERIgradAC_sphswapCD[6][3][ijlk] + gaugeERIgradAC_sphswapCD[7][1][ijlk];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[7][2][ijlk] - gaugeERIgradAC_sphswapCD[8][1][ijlk]
            - gaugeERIgradAC_sphswapCD[1][5][ijlk] + gaugeERIgradAC_sphswapCD[2][4][ijlk];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[8][0][ijlk] - gaugeERIgradAC_sphswapCD[6][2][ijlk]
            - gaugeERIgradAC_sphswapCD[2][2][ijlk] + gaugeERIgradAC_sphswapCD[0][5][ijlk];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[6][1][ijlk] - gaugeERIgradAC_sphswapCD[7][0][ijlk]
            - gaugeERIgradAC_sphswapCD[0][4][ijlk] + gaugeERIgradAC_sphswapCD[1][2][ijlk];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[1][4][ijlk] - gaugeERIgradAC_sphswapCD[2][3][ijlk]
            - gaugeERIgradAC_sphswapCD[4][2][ijlk] + gaugeERIgradAC_sphswapCD[5][1][ijlk];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[2][1][ijlk] - gaugeERIgradAC_sphswapCD[0][4][ijlk]
            - gaugeERIgradAC_sphswapCD[5][0][ijlk] + gaugeERIgradAC_sphswapCD[3][2][ijlk];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapCD[0][3][ijlk] - gaugeERIgradAC_sphswapCD[1][1][ijlk]
            - gaugeERIgradAC_sphswapCD[3][1][ijlk] + gaugeERIgradAC_sphswapCD[4][0][ijlk];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] = -(
            gaugeERIgradAC_sphswapCD[0][0][ijlk] + gaugeERIgradAC_sphswapCD[1][1][ijlk]
            + gaugeERIgradAC_sphswapCD[2][2][ijlk] + gaugeERIgradAC_sphswapCD[3][1][ijlk]
            + gaugeERIgradAC_sphswapCD[4][3][ijlk] + gaugeERIgradAC_sphswapCD[5][4][ijlk]
            + gaugeERIgradAC_sphswapCD[6][2][ijlk] + gaugeERIgradAC_sphswapCD[7][4][ijlk]
            + gaugeERIgradAC_sphswapCD[8][5][ijlk]);



// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[3][2][ijlk] + gaugeERIgradAC_sphswapCD[4][4][ijlk]
              + gaugeERIgradAC_sphswapCD[5][5][ijlk] - gaugeERIgradAC_sphswapCD[6][1][ijlk]
              - gaugeERIgradAC_sphswapCD[7][3][ijlk] - gaugeERIgradAC_sphswapCD[8][4][ijlk]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[6][0][ijlk] + gaugeERIgradAC_sphswapCD[7][1][ijlk]
              + gaugeERIgradAC_sphswapCD[8][2][ijlk] - gaugeERIgradAC_sphswapCD[0][2][ijlk]
              - gaugeERIgradAC_sphswapCD[1][4][ijlk] - gaugeERIgradAC_sphswapCD[2][5][ijlk]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[0][1][ijlk] + gaugeERIgradAC_sphswapCD[1][3][ijlk]
              + gaugeERIgradAC_sphswapCD[2][4][ijlk] - gaugeERIgradAC_sphswapCD[3][0][ijlk]
              - gaugeERIgradAC_sphswapCD[4][1][ijlk] - gaugeERIgradAC_sphswapCD[5][2][ijlk]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[1][2][ijlk] - gaugeERIgradAC_sphswapCD[2][1][ijlk]
              + gaugeERIgradAC_sphswapCD[4][4][ijlk] - gaugeERIgradAC_sphswapCD[5][3][ijlk]
              + gaugeERIgradAC_sphswapCD[7][5][ijlk] - gaugeERIgradAC_sphswapCD[8][4][ijlk]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[2][0][ijlk] - gaugeERIgradAC_sphswapCD[0][2][ijlk]
              + gaugeERIgradAC_sphswapCD[5][1][ijlk] - gaugeERIgradAC_sphswapCD[3][4][ijlk]
              + gaugeERIgradAC_sphswapCD[8][2][ijlk] - gaugeERIgradAC_sphswapCD[6][5][ijlk]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf1 * NB3 + bf2 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapCD[0][1][ijlk] - gaugeERIgradAC_sphswapCD[1][0][ijlk]
              + gaugeERIgradAC_sphswapCD[3][3][ijlk] - gaugeERIgradAC_sphswapCD[4][1][ijlk]
              + gaugeERIgradAC_sphswapCD[6][4][ijlk] - gaugeERIgradAC_sphswapCD[7][2][ijlk]);







// (nabla ji|nabla lk)
          int jilk = j * n3 * n4 * n1 + i * n3 * n4 + l * n3 + k;
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[4][5][jilk] - gaugeERIgradAC_sphswapABCD[5][4][jilk]
            - gaugeERIgradAC_sphswapABCD[7][4][jilk] + gaugeERIgradAC_sphswapABCD[8][3][jilk];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[5][2][jilk] - gaugeERIgradAC_sphswapABCD[3][5][jilk]
            - gaugeERIgradAC_sphswapABCD[8][1][jilk] + gaugeERIgradAC_sphswapABCD[6][4][jilk];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[3][4][jilk] - gaugeERIgradAC_sphswapABCD[4][2][jilk]
            - gaugeERIgradAC_sphswapABCD[6][3][jilk] + gaugeERIgradAC_sphswapABCD[7][1][jilk];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[7][2][jilk] - gaugeERIgradAC_sphswapABCD[8][1][jilk]
            - gaugeERIgradAC_sphswapABCD[1][5][jilk] + gaugeERIgradAC_sphswapABCD[2][4][jilk];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[8][0][jilk] - gaugeERIgradAC_sphswapABCD[6][2][jilk]
            - gaugeERIgradAC_sphswapABCD[2][2][jilk] + gaugeERIgradAC_sphswapABCD[0][5][jilk];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[6][1][jilk] - gaugeERIgradAC_sphswapABCD[7][0][jilk]
            - gaugeERIgradAC_sphswapABCD[0][4][jilk] + gaugeERIgradAC_sphswapABCD[1][2][jilk];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[1][4][jilk] - gaugeERIgradAC_sphswapABCD[2][3][jilk]
            - gaugeERIgradAC_sphswapABCD[4][2][jilk] + gaugeERIgradAC_sphswapABCD[5][1][jilk];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[2][1][jilk] - gaugeERIgradAC_sphswapABCD[0][4][jilk]
            - gaugeERIgradAC_sphswapABCD[5][0][jilk] + gaugeERIgradAC_sphswapABCD[3][2][jilk];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            gaugeERIgradAC_sphswapABCD[0][3][jilk] - gaugeERIgradAC_sphswapABCD[1][1][jilk]
            - gaugeERIgradAC_sphswapABCD[3][1][jilk] + gaugeERIgradAC_sphswapABCD[4][0][jilk];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] = -(
            gaugeERIgradAC_sphswapABCD[0][0][jilk] + gaugeERIgradAC_sphswapABCD[1][1][jilk]
            + gaugeERIgradAC_sphswapABCD[2][2][jilk] + gaugeERIgradAC_sphswapABCD[3][1][jilk]
            + gaugeERIgradAC_sphswapABCD[4][3][jilk] + gaugeERIgradAC_sphswapABCD[5][4][jilk]
            + gaugeERIgradAC_sphswapABCD[6][2][jilk] + gaugeERIgradAC_sphswapABCD[7][4][jilk]
            + gaugeERIgradAC_sphswapABCD[8][5][jilk]);


// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[3][2][jilk] + gaugeERIgradAC_sphswapABCD[4][4][jilk]
              + gaugeERIgradAC_sphswapABCD[5][5][jilk] - gaugeERIgradAC_sphswapABCD[6][1][jilk]
              - gaugeERIgradAC_sphswapABCD[7][3][jilk] - gaugeERIgradAC_sphswapABCD[8][4][jilk]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[6][0][jilk] + gaugeERIgradAC_sphswapABCD[7][1][jilk]
              + gaugeERIgradAC_sphswapABCD[8][2][jilk] - gaugeERIgradAC_sphswapABCD[0][2][jilk]
              - gaugeERIgradAC_sphswapABCD[1][4][jilk] - gaugeERIgradAC_sphswapABCD[2][5][jilk]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[0][1][jilk] + gaugeERIgradAC_sphswapABCD[1][3][jilk]
              + gaugeERIgradAC_sphswapABCD[2][4][jilk] - gaugeERIgradAC_sphswapABCD[3][0][jilk]
              - gaugeERIgradAC_sphswapABCD[4][1][jilk] - gaugeERIgradAC_sphswapABCD[5][2][jilk]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[1][2][jilk] - gaugeERIgradAC_sphswapABCD[2][1][jilk]
              + gaugeERIgradAC_sphswapABCD[4][4][jilk] - gaugeERIgradAC_sphswapABCD[5][3][jilk]
              + gaugeERIgradAC_sphswapABCD[7][5][jilk] - gaugeERIgradAC_sphswapABCD[8][4][jilk]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[2][0][jilk] - gaugeERIgradAC_sphswapABCD[0][2][jilk]
              + gaugeERIgradAC_sphswapABCD[5][1][jilk] - gaugeERIgradAC_sphswapABCD[3][4][jilk]
              + gaugeERIgradAC_sphswapABCD[8][2][jilk] - gaugeERIgradAC_sphswapABCD[6][5][jilk]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf2 * NB3 + bf1 * NB2 + bf4 * NB + bf3] =
            -(gaugeERIgradAC_sphswapABCD[0][1][jilk] - gaugeERIgradAC_sphswapABCD[1][0][jilk]
              + gaugeERIgradAC_sphswapABCD[3][3][jilk] - gaugeERIgradAC_sphswapABCD[4][1][jilk]
              + gaugeERIgradAC_sphswapABCD[6][4][jilk] - gaugeERIgradAC_sphswapABCD[7][2][jilk]);








// swap electron 1 and 2
//
// (nabla kl|nabla ij)
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[4][5][ijkl] - gaugeERIgradAC_sph[7][4][ijkl]
            - gaugeERIgradAC_sph[5][4][ijkl] + gaugeERIgradAC_sph[8][3][ijkl];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[7][2][ijkl] - gaugeERIgradAC_sph[1][5][ijkl]
            - gaugeERIgradAC_sph[8][1][ijkl] + gaugeERIgradAC_sph[2][4][ijkl];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[1][4][ijkl] - gaugeERIgradAC_sph[4][2][ijkl]
            - gaugeERIgradAC_sph[2][3][ijkl] + gaugeERIgradAC_sph[5][1][ijkl];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[5][2][ijkl] - gaugeERIgradAC_sph[8][1][ijkl]
            - gaugeERIgradAC_sph[3][5][ijkl] + gaugeERIgradAC_sph[6][4][ijkl];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[8][0][ijkl] - gaugeERIgradAC_sph[2][2][ijkl]
            - gaugeERIgradAC_sph[6][2][ijkl] + gaugeERIgradAC_sph[0][5][ijkl];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[2][1][ijkl] - gaugeERIgradAC_sph[5][0][ijkl]
            - gaugeERIgradAC_sph[0][4][ijkl] + gaugeERIgradAC_sph[3][2][ijkl];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[3][4][ijkl] - gaugeERIgradAC_sph[6][3][ijkl]
            - gaugeERIgradAC_sph[4][2][ijkl] + gaugeERIgradAC_sph[7][1][ijkl];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[6][1][ijkl] - gaugeERIgradAC_sph[0][4][ijkl]
            - gaugeERIgradAC_sph[7][0][ijkl] + gaugeERIgradAC_sph[1][2][ijkl];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sph[0][3][ijkl] - gaugeERIgradAC_sph[3][1][ijkl]
            - gaugeERIgradAC_sph[1][1][ijkl] + gaugeERIgradAC_sph[4][0][ijkl];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] = -(
            gaugeERIgradAC_sph[0][0][ijkl] + gaugeERIgradAC_sph[1][1][ijkl]
            + gaugeERIgradAC_sph[2][2][ijkl] + gaugeERIgradAC_sph[3][1][ijkl]
            + gaugeERIgradAC_sph[4][3][ijkl] + gaugeERIgradAC_sph[5][4][ijkl]
            + gaugeERIgradAC_sph[6][2][ijkl] + gaugeERIgradAC_sph[7][4][ijkl]
            + gaugeERIgradAC_sph[8][5][ijkl]);
// or
// gaugeSLSLACxxtrue[15*NB4 + bf3* NB3+ bf4*NB2 + bf1*NB + bf2]
// =gaugeSLSLACxxtrue[15*NB4 + bf1* NB3+ bf2*NB2 + bf3*NB + bf4] ;



// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[1][2][ijkl] + gaugeERIgradAC_sph[4][4][ijkl]
              + gaugeERIgradAC_sph[7][5][ijkl] - gaugeERIgradAC_sph[2][1][ijkl]
              - gaugeERIgradAC_sph[5][3][ijkl] - gaugeERIgradAC_sph[8][4][ijkl]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[2][0][ijkl] + gaugeERIgradAC_sph[5][1][ijkl]
              + gaugeERIgradAC_sph[8][2][ijkl] - gaugeERIgradAC_sph[0][2][ijkl]
              - gaugeERIgradAC_sph[3][4][ijkl] - gaugeERIgradAC_sph[6][5][ijkl]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[0][1][ijkl] + gaugeERIgradAC_sph[3][3][ijkl]
              + gaugeERIgradAC_sph[6][4][ijkl] - gaugeERIgradAC_sph[1][0][ijkl]
              - gaugeERIgradAC_sph[4][1][ijkl] - gaugeERIgradAC_sph[7][2][ijkl]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[3][2][ijkl] - gaugeERIgradAC_sph[6][1][ijkl]
              + gaugeERIgradAC_sph[4][4][ijkl] - gaugeERIgradAC_sph[7][3][ijkl]
              + gaugeERIgradAC_sph[5][5][ijkl] - gaugeERIgradAC_sph[8][4][ijkl]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[6][0][ijkl] - gaugeERIgradAC_sph[0][2][ijkl]
              + gaugeERIgradAC_sph[7][1][ijkl] - gaugeERIgradAC_sph[1][4][ijkl]
              + gaugeERIgradAC_sph[8][2][ijkl] - gaugeERIgradAC_sph[2][5][ijkl]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf3 * NB3 + bf4 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sph[0][1][ijkl] - gaugeERIgradAC_sph[3][0][ijkl]
              + gaugeERIgradAC_sph[1][3][ijkl] - gaugeERIgradAC_sph[4][1][ijkl]
              + gaugeERIgradAC_sph[2][4][ijkl] - gaugeERIgradAC_sph[5][2][ijkl]);






// (nabla kl|nabla ji)
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[4][5][jikl] - gaugeERIgradAC_sphswapAB[7][4][jikl]
            - gaugeERIgradAC_sphswapAB[5][4][jikl] + gaugeERIgradAC_sphswapAB[8][3][jikl];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[7][2][jikl] - gaugeERIgradAC_sphswapAB[1][5][jikl]
            - gaugeERIgradAC_sphswapAB[8][1][jikl] + gaugeERIgradAC_sphswapAB[2][4][jikl];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[1][4][jikl] - gaugeERIgradAC_sphswapAB[4][2][jikl]
            - gaugeERIgradAC_sphswapAB[2][3][jikl] + gaugeERIgradAC_sphswapAB[5][1][jikl];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[5][2][jikl] - gaugeERIgradAC_sphswapAB[8][1][jikl]
            - gaugeERIgradAC_sphswapAB[3][5][jikl] + gaugeERIgradAC_sphswapAB[6][4][jikl];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[8][0][jikl] - gaugeERIgradAC_sphswapAB[2][2][jikl]
            - gaugeERIgradAC_sphswapAB[6][2][jikl] + gaugeERIgradAC_sphswapAB[0][5][jikl];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[2][1][jikl] - gaugeERIgradAC_sphswapAB[5][0][jikl]
            - gaugeERIgradAC_sphswapAB[0][4][jikl] + gaugeERIgradAC_sphswapAB[3][2][jikl];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[3][4][jikl] - gaugeERIgradAC_sphswapAB[6][3][jikl]
            - gaugeERIgradAC_sphswapAB[4][2][jikl] + gaugeERIgradAC_sphswapAB[7][1][jikl];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[6][1][jikl] - gaugeERIgradAC_sphswapAB[0][4][jikl]
            - gaugeERIgradAC_sphswapAB[7][0][jikl] + gaugeERIgradAC_sphswapAB[1][2][jikl];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapAB[0][3][jikl] - gaugeERIgradAC_sphswapAB[3][1][jikl]
            - gaugeERIgradAC_sphswapAB[1][1][jikl] + gaugeERIgradAC_sphswapAB[4][0][jikl];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] = -(
            gaugeERIgradAC_sphswapAB[0][0][jikl] + gaugeERIgradAC_sphswapAB[1][1][jikl]
            + gaugeERIgradAC_sphswapAB[2][2][jikl] + gaugeERIgradAC_sphswapAB[3][1][jikl]
            + gaugeERIgradAC_sphswapAB[4][3][jikl] + gaugeERIgradAC_sphswapAB[5][4][jikl]
            + gaugeERIgradAC_sphswapAB[6][2][jikl] + gaugeERIgradAC_sphswapAB[7][4][jikl]
            + gaugeERIgradAC_sphswapAB[8][5][jikl]);
//or
//gaugeSLSLACxxtrue[15*NB4 + bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
//gaugeSLSLACxxtrue[15*NB4 + bf2* NB3+ bf1*NB2 + bf3*NB + bf4];



// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[1][2][jikl] + gaugeERIgradAC_sphswapAB[4][4][jikl]
              + gaugeERIgradAC_sphswapAB[7][5][jikl] - gaugeERIgradAC_sphswapAB[2][1][jikl]
              - gaugeERIgradAC_sphswapAB[5][3][jikl] - gaugeERIgradAC_sphswapAB[8][4][jikl]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[2][0][jikl] + gaugeERIgradAC_sphswapAB[5][1][jikl]
              + gaugeERIgradAC_sphswapAB[8][2][jikl] - gaugeERIgradAC_sphswapAB[0][2][jikl]
              - gaugeERIgradAC_sphswapAB[3][4][jikl] - gaugeERIgradAC_sphswapAB[6][5][jikl]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[0][1][jikl] + gaugeERIgradAC_sphswapAB[3][3][jikl]
              + gaugeERIgradAC_sphswapAB[6][4][jikl] - gaugeERIgradAC_sphswapAB[1][0][jikl]
              - gaugeERIgradAC_sphswapAB[4][1][jikl] - gaugeERIgradAC_sphswapAB[7][2][jikl]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[3][2][jikl] - gaugeERIgradAC_sphswapAB[6][1][jikl]
              + gaugeERIgradAC_sphswapAB[4][4][jikl] - gaugeERIgradAC_sphswapAB[7][3][jikl]
              + gaugeERIgradAC_sphswapAB[5][5][jikl] - gaugeERIgradAC_sphswapAB[8][4][jikl]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[6][0][jikl] - gaugeERIgradAC_sphswapAB[0][2][jikl]
              + gaugeERIgradAC_sphswapAB[7][1][jikl] - gaugeERIgradAC_sphswapAB[1][4][jikl]
              + gaugeERIgradAC_sphswapAB[8][2][jikl] - gaugeERIgradAC_sphswapAB[2][5][jikl]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf3 * NB3 + bf4 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapAB[0][1][jikl] - gaugeERIgradAC_sphswapAB[3][0][jikl]
              + gaugeERIgradAC_sphswapAB[1][3][jikl] - gaugeERIgradAC_sphswapAB[4][1][jikl]
              + gaugeERIgradAC_sphswapAB[2][4][jikl] - gaugeERIgradAC_sphswapAB[5][2][jikl]);





// (nabla lk|nabla ij)
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[4][5][ijlk] - gaugeERIgradAC_sphswapCD[7][4][ijlk]
            - gaugeERIgradAC_sphswapCD[5][4][ijlk] + gaugeERIgradAC_sphswapCD[8][3][ijlk];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[7][2][ijlk] - gaugeERIgradAC_sphswapCD[1][5][ijlk]
            - gaugeERIgradAC_sphswapCD[8][1][ijlk] + gaugeERIgradAC_sphswapCD[2][4][ijlk];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[1][4][ijlk] - gaugeERIgradAC_sphswapCD[4][2][ijlk]
            - gaugeERIgradAC_sphswapCD[2][3][ijlk] + gaugeERIgradAC_sphswapCD[5][1][ijlk];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[5][2][ijlk] - gaugeERIgradAC_sphswapCD[8][1][ijlk]
            - gaugeERIgradAC_sphswapCD[3][5][ijlk] + gaugeERIgradAC_sphswapCD[6][4][ijlk];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[8][0][ijlk] - gaugeERIgradAC_sphswapCD[2][2][ijlk]
            - gaugeERIgradAC_sphswapCD[6][2][ijlk] + gaugeERIgradAC_sphswapCD[0][5][ijlk];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[2][1][ijlk] - gaugeERIgradAC_sphswapCD[5][0][ijlk]
            - gaugeERIgradAC_sphswapCD[0][4][ijlk] + gaugeERIgradAC_sphswapCD[3][2][ijlk];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[3][4][ijlk] - gaugeERIgradAC_sphswapCD[6][3][ijlk]
            - gaugeERIgradAC_sphswapCD[4][2][ijlk] + gaugeERIgradAC_sphswapCD[7][1][ijlk];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[6][1][ijlk] - gaugeERIgradAC_sphswapCD[0][4][ijlk]
            - gaugeERIgradAC_sphswapCD[7][0][ijlk] + gaugeERIgradAC_sphswapCD[1][2][ijlk];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            gaugeERIgradAC_sphswapCD[0][3][ijlk] - gaugeERIgradAC_sphswapCD[3][1][ijlk]
            - gaugeERIgradAC_sphswapCD[1][1][ijlk] + gaugeERIgradAC_sphswapCD[4][0][ijlk];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] = -(
            gaugeERIgradAC_sphswapCD[0][0][ijlk] + gaugeERIgradAC_sphswapCD[1][1][ijlk]
            + gaugeERIgradAC_sphswapCD[2][2][ijlk] + gaugeERIgradAC_sphswapCD[3][1][ijlk]
            + gaugeERIgradAC_sphswapCD[4][3][ijlk] + gaugeERIgradAC_sphswapCD[5][4][ijlk]
            + gaugeERIgradAC_sphswapCD[6][2][ijlk] + gaugeERIgradAC_sphswapCD[7][4][ijlk]
            + gaugeERIgradAC_sphswapCD[8][5][ijlk]);
// or
// gaugeSLSLACxxtrue[15*NB4 + bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
// gaugeSLSLACxxtrue[15*NB4 + bf1* NB3+ bf2*NB2 + bf4*NB + bf3];



// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[1][2][ijlk] + gaugeERIgradAC_sphswapCD[4][4][ijlk]
              + gaugeERIgradAC_sphswapCD[7][5][ijlk] - gaugeERIgradAC_sphswapCD[2][1][ijlk]
              - gaugeERIgradAC_sphswapCD[5][3][ijlk] - gaugeERIgradAC_sphswapCD[8][4][ijlk]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[2][0][ijlk] + gaugeERIgradAC_sphswapCD[5][1][ijlk]
              + gaugeERIgradAC_sphswapCD[8][2][ijlk] - gaugeERIgradAC_sphswapCD[0][2][ijlk]
              - gaugeERIgradAC_sphswapCD[3][4][ijlk] - gaugeERIgradAC_sphswapCD[6][5][ijlk]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[0][1][ijlk] + gaugeERIgradAC_sphswapCD[3][3][ijlk]
              + gaugeERIgradAC_sphswapCD[6][4][ijlk] - gaugeERIgradAC_sphswapCD[1][0][ijlk]
              - gaugeERIgradAC_sphswapCD[4][1][ijlk] - gaugeERIgradAC_sphswapCD[7][2][ijlk]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[3][2][ijlk] - gaugeERIgradAC_sphswapCD[6][1][ijlk]
              + gaugeERIgradAC_sphswapCD[4][4][ijlk] - gaugeERIgradAC_sphswapCD[7][3][ijlk]
              + gaugeERIgradAC_sphswapCD[5][5][ijlk] - gaugeERIgradAC_sphswapCD[8][4][ijlk]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[6][0][ijlk] - gaugeERIgradAC_sphswapCD[0][2][ijlk]
              + gaugeERIgradAC_sphswapCD[7][1][ijlk] - gaugeERIgradAC_sphswapCD[1][4][ijlk]
              + gaugeERIgradAC_sphswapCD[8][2][ijlk] - gaugeERIgradAC_sphswapCD[2][5][ijlk]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf4 * NB3 + bf3 * NB2 + bf1 * NB + bf2] =
            -(gaugeERIgradAC_sphswapCD[0][1][ijlk] - gaugeERIgradAC_sphswapCD[3][0][ijlk]
              + gaugeERIgradAC_sphswapCD[1][3][ijlk] - gaugeERIgradAC_sphswapCD[4][1][ijlk]
              + gaugeERIgradAC_sphswapCD[2][4][ijlk] - gaugeERIgradAC_sphswapCD[5][2][ijlk]);










// (nabla lk|nabla ji)
// sigma1x sigma2x
          gaugeSLSLACxxtrue[bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[4][5][jilk] - gaugeERIgradAC_sphswapABCD[7][4][jilk]
            - gaugeERIgradAC_sphswapABCD[5][4][jilk] + gaugeERIgradAC_sphswapABCD[8][3][jilk];

// sigma1x sigma2y
          gaugeSLSLACxxtrue[NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[7][2][jilk] - gaugeERIgradAC_sphswapABCD[1][5][jilk]
            - gaugeERIgradAC_sphswapABCD[8][1][jilk] + gaugeERIgradAC_sphswapABCD[2][4][jilk];

// sigma1x sigma2z
          gaugeSLSLACxxtrue[2 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[1][4][jilk] - gaugeERIgradAC_sphswapABCD[4][2][jilk]
            - gaugeERIgradAC_sphswapABCD[2][3][jilk] + gaugeERIgradAC_sphswapABCD[5][1][jilk];

// sigma1y sigma2x
          gaugeSLSLACxxtrue[3 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[5][2][jilk] - gaugeERIgradAC_sphswapABCD[8][1][jilk]
            - gaugeERIgradAC_sphswapABCD[3][5][jilk] + gaugeERIgradAC_sphswapABCD[6][4][jilk];

// sigma1y sigma2y
          gaugeSLSLACxxtrue[4 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[8][0][jilk] - gaugeERIgradAC_sphswapABCD[2][2][jilk]
            - gaugeERIgradAC_sphswapABCD[6][2][jilk] + gaugeERIgradAC_sphswapABCD[0][5][jilk];

// sigma1y sigma2z
          gaugeSLSLACxxtrue[5 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[2][1][jilk] - gaugeERIgradAC_sphswapABCD[5][0][jilk]
            - gaugeERIgradAC_sphswapABCD[0][4][jilk] + gaugeERIgradAC_sphswapABCD[3][2][jilk];

// sigma1z sigma2x
          gaugeSLSLACxxtrue[6 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[3][4][jilk] - gaugeERIgradAC_sphswapABCD[6][3][jilk]
            - gaugeERIgradAC_sphswapABCD[4][2][jilk] + gaugeERIgradAC_sphswapABCD[7][1][jilk];

// sigma1z sigma2y
          gaugeSLSLACxxtrue[7 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[6][1][jilk] - gaugeERIgradAC_sphswapABCD[0][4][jilk]
            - gaugeERIgradAC_sphswapABCD[7][0][jilk] + gaugeERIgradAC_sphswapABCD[1][2][jilk];

// sigma1z sigma2z
          gaugeSLSLACxxtrue[8 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            gaugeERIgradAC_sphswapABCD[0][3][jilk] - gaugeERIgradAC_sphswapABCD[3][1][jilk]
            - gaugeERIgradAC_sphswapABCD[1][1][jilk] + gaugeERIgradAC_sphswapABCD[4][0][jilk];


// I1 I2
          gaugeSLSLACxxtrue[15 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] = -(
            gaugeERIgradAC_sphswapABCD[0][0][jilk] + gaugeERIgradAC_sphswapABCD[1][1][jilk]
            + gaugeERIgradAC_sphswapABCD[2][2][jilk] + gaugeERIgradAC_sphswapABCD[3][1][jilk]
            + gaugeERIgradAC_sphswapABCD[4][3][jilk] + gaugeERIgradAC_sphswapABCD[5][4][jilk]
            + gaugeERIgradAC_sphswapABCD[6][2][jilk] + gaugeERIgradAC_sphswapABCD[7][4][jilk]
            + gaugeERIgradAC_sphswapABCD[8][5][jilk]);
// or
// gaugeSLSLACxxtrue[15*NB4 + bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
// gaugeSLSLACxxtrue[15*NB4 + bf2* NB3+ bf1*NB2 + bf4*NB + bf3];



// sigma1 x I2
          gaugeSLSLACxxtrue[9 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[1][2][jilk] + gaugeERIgradAC_sphswapABCD[4][4][jilk]
              + gaugeERIgradAC_sphswapABCD[7][5][jilk] - gaugeERIgradAC_sphswapABCD[2][1][jilk]
              - gaugeERIgradAC_sphswapABCD[5][3][jilk] - gaugeERIgradAC_sphswapABCD[8][4][jilk]);

// sigma1 y I2
          gaugeSLSLACxxtrue[10 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[2][0][jilk] + gaugeERIgradAC_sphswapABCD[5][1][jilk]
              + gaugeERIgradAC_sphswapABCD[8][2][jilk] - gaugeERIgradAC_sphswapABCD[0][2][jilk]
              - gaugeERIgradAC_sphswapABCD[3][4][jilk] - gaugeERIgradAC_sphswapABCD[6][5][jilk]);

// sigma1 z I2
          gaugeSLSLACxxtrue[11 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[0][1][jilk] + gaugeERIgradAC_sphswapABCD[3][3][jilk]
              + gaugeERIgradAC_sphswapABCD[6][4][jilk] - gaugeERIgradAC_sphswapABCD[1][0][jilk]
              - gaugeERIgradAC_sphswapABCD[4][1][jilk] - gaugeERIgradAC_sphswapABCD[7][2][jilk]);


// I1 sigma2 x
          gaugeSLSLACxxtrue[12 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[3][2][jilk] - gaugeERIgradAC_sphswapABCD[6][1][jilk]
              + gaugeERIgradAC_sphswapABCD[4][4][jilk] - gaugeERIgradAC_sphswapABCD[7][3][jilk]
              + gaugeERIgradAC_sphswapABCD[5][5][jilk] - gaugeERIgradAC_sphswapABCD[8][4][jilk]);

// I1 sigma2 y
          gaugeSLSLACxxtrue[13 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[6][0][jilk] - gaugeERIgradAC_sphswapABCD[0][2][jilk]
              + gaugeERIgradAC_sphswapABCD[7][1][jilk] - gaugeERIgradAC_sphswapABCD[1][4][jilk]
              + gaugeERIgradAC_sphswapABCD[8][2][jilk] - gaugeERIgradAC_sphswapABCD[2][5][jilk]);

// I1 sigma2 z
          gaugeSLSLACxxtrue[14 * NB4 + bf4 * NB3 + bf3 * NB2 + bf2 * NB + bf1] =
            -(gaugeERIgradAC_sphswapABCD[0][1][jilk] - gaugeERIgradAC_sphswapABCD[3][0][jilk]
              + gaugeERIgradAC_sphswapABCD[1][3][jilk] - gaugeERIgradAC_sphswapABCD[4][1][jilk]
              + gaugeERIgradAC_sphswapABCD[2][4][jilk] - gaugeERIgradAC_sphswapABCD[5][2][jilk]);




          //LSSL

          /*
          //I1I2
          gaugeSLSLACxx[ 15*NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[0][0][jikl]+gaugeERIgradAC_sphswapAB[3][1][jikl]
          +gaugeERIgradAC_sphswapAB[6][2][jikl]+ gaugeERIgradAC_sphswapAB[1][1][jikl]
          +gaugeERIgradAC_sphswapAB[4][3][jikl] + gaugeERIgradAC_sphswapAB[7][4][jikl]
          +gaugeERIgradAC_sphswapAB[2][2][jikl] +gaugeERIgradAC_sphswapAB[5][4][jikl]
          +gaugeERIgradAC_sphswapAB[8][5][jikl] ;


          //swap AB
          gaugeSLSLACxx[ 15*NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[0][0][ijkl]+gaugeERIgradAC_sph[3][1][ijkl]
          +gaugeERIgradAC_sph[6][2][ijkl]+ gaugeERIgradAC_sph[1][1][ijkl]
          +gaugeERIgradAC_sph[4][3][ijkl] + gaugeERIgradAC_sph[7][4][ijkl]
          +gaugeERIgradAC_sph[2][2][ijkl] + gaugeERIgradAC_sph[5][4][ijkl]
          +gaugeERIgradAC_sph[8][5][ijkl] ;

          //swap CD
          gaugeSLSLACxx[ 15*NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[0][0][jilk]+gaugeERIgradAC_sphswapABCD[3][1][jilk]
          +gaugeERIgradAC_sphswapABCD[6][2][jilk]+ gaugeERIgradAC_sphswapABCD[1][1][jilk]
          +gaugeERIgradAC_sphswapABCD[4][3][jilk] + gaugeERIgradAC_sphswapABCD[7][4][jilk]
          +gaugeERIgradAC_sphswapABCD[2][2][jilk] +gaugeERIgradAC_sphswapABCD[5][4][jilk]
          +gaugeERIgradAC_sphswapABCD[8][5][jilk] ;

          //swap ABCD
          gaugeSLSLACxx[ 15*NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[0][0][ijlk]+gaugeERIgradAC_sphswapCD[3][1][ijlk]
          +gaugeERIgradAC_sphswapCD[6][2][ijlk]+ gaugeERIgradAC_sphswapCD[1][1][ijlk]
          +gaugeERIgradAC_sphswapCD[4][3][ijlk] + gaugeERIgradAC_sphswapCD[7][4][ijlk]
          +gaugeERIgradAC_sphswapCD[2][2][ijlk] +gaugeERIgradAC_sphswapCD[5][4][ijlk]
          +gaugeERIgradAC_sphswapCD[8][5][ijlk] ;

          //swap electron 1 2
          gaugeSLSLACxx[ 15*NB4+ bf3* NB3+ bf4*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapAB[0][0][jikl]+gaugeERIgradAC_sphswapAB[3][1][jikl]
          +gaugeERIgradAC_sphswapAB[6][2][jikl]+ gaugeERIgradAC_sphswapAB[1][1][jikl]
          +gaugeERIgradAC_sphswapAB[4][3][jikl] + gaugeERIgradAC_sphswapAB[7][4][jikl]
          +gaugeERIgradAC_sphswapAB[2][2][jikl] +gaugeERIgradAC_sphswapAB[5][4][jikl]
          +gaugeERIgradAC_sphswapAB[8][5][jikl] ;


          //swap AB
          gaugeSLSLACxx[ 15*NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[0][0][ijkl]+gaugeERIgradAC_sph[3][1][ijkl]
          +gaugeERIgradAC_sph[6][2][ijkl]+ gaugeERIgradAC_sph[1][1][ijkl]
          +gaugeERIgradAC_sph[4][3][ijkl] + gaugeERIgradAC_sph[7][4][ijkl]
          +gaugeERIgradAC_sph[2][2][ijkl] +gaugeERIgradAC_sph[5][4][ijkl]
          +gaugeERIgradAC_sph[8][5][ijkl] ;

          //swap CD
          gaugeSLSLACxx[ 15*NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[0][0][jilk]+gaugeERIgradAC_sphswapABCD[3][1][jilk]
          +gaugeERIgradAC_sphswapABCD[6][2][jilk]+ gaugeERIgradAC_sphswapABCD[1][1][jilk]
          +gaugeERIgradAC_sphswapABCD[4][3][jilk] + gaugeERIgradAC_sphswapABCD[7][4][jilk]
          +gaugeERIgradAC_sphswapABCD[2][2][jilk] +gaugeERIgradAC_sphswapABCD[5][4][jilk]
          +gaugeERIgradAC_sphswapABCD[8][5][jilk] ;

          //swap ABCD
          gaugeSLSLACxx[ 15*NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[0][0][ijlk]+gaugeERIgradAC_sphswapCD[3][1][ijlk]
          +gaugeERIgradAC_sphswapCD[6][2][ijlk]+ gaugeERIgradAC_sphswapCD[1][1][ijlk]
          +gaugeERIgradAC_sphswapCD[4][3][ijlk] + gaugeERIgradAC_sphswapCD[7][4][ijlk]
          +gaugeERIgradAC_sphswapCD[2][2][ijlk] +gaugeERIgradAC_sphswapCD[5][4][ijlk]
          +gaugeERIgradAC_sphswapCD[8][5][ijlk] ;



          // sigma1 sigma2


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[7][4][jikl] - gaugeERIgradAC_sphswapAB[8][3][jikl]
          - gaugeERIgradAC_sphswapAB[1][5][jikl] +gaugeERIgradAC_sphswapAB[2][4][jikl] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[1][5][jikl] - gaugeERIgradAC_sphswapAB[2][4][jikl]
          - gaugeERIgradAC_sphswapAB[7][2][jikl] +gaugeERIgradAC_sphswapAB[8][1][jikl] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[4][2][jikl] - gaugeERIgradAC_sphswapAB[1][4][jikl]
          - gaugeERIgradAC_sphswapAB[5][1][jikl] +gaugeERIgradAC_sphswapAB[2][3][jikl] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[8][1][jikl] - gaugeERIgradAC_sphswapAB[5][2][jikl]
          - gaugeERIgradAC_sphswapAB[6][4][jikl] +gaugeERIgradAC_sphswapAB[3][5][jikl] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[2][2][jikl] - gaugeERIgradAC_sphswapAB[8][0][jikl]
          - gaugeERIgradAC_sphswapAB[0][5][jikl] +gaugeERIgradAC_sphswapAB[6][2][jikl] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[5][0][jikl] - gaugeERIgradAC_sphswapAB[2][1][jikl]
          - gaugeERIgradAC_sphswapAB[3][2][jikl] +gaugeERIgradAC_sphswapAB[0][4][jikl] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[6][3][jikl] - gaugeERIgradAC_sphswapAB[3][4][jikl]
          - gaugeERIgradAC_sphswapAB[7][1][jikl] +gaugeERIgradAC_sphswapAB[4][2][jikl] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[0][4][jikl] - gaugeERIgradAC_sphswapAB[6][1][jikl]
          - gaugeERIgradAC_sphswapAB[1][2][jikl] +gaugeERIgradAC_sphswapAB[7][0][jikl] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sphswapAB[3][1][jikl] - gaugeERIgradAC_sphswapAB[0][3][jikl]
          - gaugeERIgradAC_sphswapAB[4][0][jikl] +gaugeERIgradAC_sphswapAB[1][1][jikl] ;


          //swap AB


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[7][4][ijkl] - gaugeERIgradAC_sph[8][3][ijkl]
          - gaugeERIgradAC_sph[1][5][ijkl] +gaugeERIgradAC_sph[2][4][ijkl] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[1][5][ijkl] - gaugeERIgradAC_sph[2][4][ijkl]
          - gaugeERIgradAC_sph[7][2][ijkl] +gaugeERIgradAC_sph[8][1][ijkl] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[4][2][ijkl] - gaugeERIgradAC_sph[1][4][ijkl]
          - gaugeERIgradAC_sph[5][1][ijkl] +gaugeERIgradAC_sph[2][3][ijkl] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[8][1][ijkl] - gaugeERIgradAC_sph[5][2][ijkl]
          - gaugeERIgradAC_sph[6][4][ijkl] +gaugeERIgradAC_sph[3][5][ijkl] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[2][2][ijkl] - gaugeERIgradAC_sph[8][0][ijkl]
          - gaugeERIgradAC_sph[0][5][ijkl] +gaugeERIgradAC_sph[6][2][ijkl] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[5][0][ijkl] - gaugeERIgradAC_sph[2][1][ijkl]
          - gaugeERIgradAC_sph[3][2][ijkl] +gaugeERIgradAC_sph[0][4][ijkl] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[6][3][ijkl] - gaugeERIgradAC_sph[3][4][ijkl]
          - gaugeERIgradAC_sph[7][1][ijkl] +gaugeERIgradAC_sph[4][2][ijkl] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[0][4][ijkl] - gaugeERIgradAC_sph[6][1][ijkl]
          - gaugeERIgradAC_sph[1][2][ijkl] +gaugeERIgradAC_sph[7][0][ijkl] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf2* NB3+ bf1*NB2 + bf3*NB + bf4] =
          gaugeERIgradAC_sph[3][1][ijkl] - gaugeERIgradAC_sph[0][3][ijkl]
          - gaugeERIgradAC_sph[4][0][ijkl] +gaugeERIgradAC_sph[1][1][ijkl] ;


          //swap CD


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[7][4][jilk] - gaugeERIgradAC_sphswapABCD[8][3][jilk]
          - gaugeERIgradAC_sphswapABCD[1][5][jilk] +gaugeERIgradAC_sphswapABCD[2][4][jilk] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[1][5][jilk] - gaugeERIgradAC_sphswapABCD[2][4][jilk]
          - gaugeERIgradAC_sphswapABCD[7][2][jilk] +gaugeERIgradAC_sphswapABCD[8][1][jilk] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[4][2][jilk] - gaugeERIgradAC_sphswapABCD[1][4][jilk]
          - gaugeERIgradAC_sphswapABCD[5][1][jilk] +gaugeERIgradAC_sphswapABCD[2][3][jilk] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[8][1][jilk] - gaugeERIgradAC_sphswapABCD[5][2][jilk]
          - gaugeERIgradAC_sphswapABCD[6][4][jilk] +gaugeERIgradAC_sphswapABCD[3][5][jilk] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[2][2][jilk] - gaugeERIgradAC_sphswapABCD[8][0][jilk]
          - gaugeERIgradAC_sphswapABCD[0][5][jilk] +gaugeERIgradAC_sphswapABCD[6][2][jilk] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[5][0][jilk] - gaugeERIgradAC_sphswapABCD[2][1][jilk]
          - gaugeERIgradAC_sphswapABCD[3][2][jilk] +gaugeERIgradAC_sphswapABCD[0][4][jilk] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[6][3][jilk] - gaugeERIgradAC_sphswapABCD[3][4][jilk]
          - gaugeERIgradAC_sphswapABCD[7][1][jilk] +gaugeERIgradAC_sphswapABCD[4][2][jilk] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[0][4][jilk] - gaugeERIgradAC_sphswapABCD[6][1][jilk]
          - gaugeERIgradAC_sphswapABCD[1][2][jilk] +gaugeERIgradAC_sphswapABCD[7][0][jilk] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf1* NB3+ bf2*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapABCD[3][1][jilk] - gaugeERIgradAC_sphswapABCD[0][3][jilk]
          - gaugeERIgradAC_sphswapABCD[4][0][jilk] +gaugeERIgradAC_sphswapABCD[1][1][jilk] ;


          //swap ABCD


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[7][4][ijlk] - gaugeERIgradAC_sphswapCD[8][3][ijlk]
          - gaugeERIgradAC_sphswapCD[1][5][ijlk] +gaugeERIgradAC_sphswapCD[2][4][ijlk] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[1][5][ijlk] - gaugeERIgradAC_sphswapCD[2][4][ijlk]
          - gaugeERIgradAC_sphswapCD[7][2][ijlk] +gaugeERIgradAC_sphswapCD[8][1][ijlk] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[4][2][ijlk] - gaugeERIgradAC_sphswapCD[1][4][ijlk]
          - gaugeERIgradAC_sphswapCD[5][1][ijlk] +gaugeERIgradAC_sphswapCD[2][3][ijlk] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[8][1][ijlk] - gaugeERIgradAC_sphswapCD[5][2][ijlk]
          - gaugeERIgradAC_sphswapCD[6][4][ijlk] +gaugeERIgradAC_sphswapCD[3][5][ijlk] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[2][2][ijlk] - gaugeERIgradAC_sphswapCD[8][0][ijlk]
          - gaugeERIgradAC_sphswapCD[0][5][ijlk] +gaugeERIgradAC_sphswapCD[6][2][ijlk] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[5][0][ijlk] - gaugeERIgradAC_sphswapCD[2][1][ijlk]
          - gaugeERIgradAC_sphswapCD[3][2][ijlk] +gaugeERIgradAC_sphswapCD[0][4][ijlk] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[6][3][ijlk] - gaugeERIgradAC_sphswapCD[3][4][ijlk]
          - gaugeERIgradAC_sphswapCD[7][1][ijlk] +gaugeERIgradAC_sphswapCD[4][2][ijlk] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[0][4][ijlk] - gaugeERIgradAC_sphswapCD[6][1][ijlk]
          - gaugeERIgradAC_sphswapCD[1][2][ijlk] +gaugeERIgradAC_sphswapCD[7][0][ijlk] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf2* NB3+ bf1*NB2 + bf4*NB + bf3] =
          gaugeERIgradAC_sphswapCD[3][1][ijlk] - gaugeERIgradAC_sphswapCD[0][3][ijlk]
          - gaugeERIgradAC_sphswapCD[4][0][ijlk] +gaugeERIgradAC_sphswapCD[1][1][ijlk] ;


          //swap electron 1 2

          // sigma1x sigma2x
          gaugeSLSLACxx[  bf3*NB3 + bf4*NB2 +bf1* NB + bf2] =
          gaugeERIgradAC_sphswapAB[5][4][jikl] - gaugeERIgradAC_sphswapAB[8][3][jikl]
          - gaugeERIgradAC_sphswapAB[3][5][jikl] +gaugeERIgradAC_sphswapAB[6][4][jikl] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf3*NB3 + bf4*NB2 + bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[3][5][jikl] - gaugeERIgradAC_sphswapAB[6][4][jikl]
          - gaugeERIgradAC_sphswapAB[5][2][jikl] +gaugeERIgradAC_sphswapAB[8][1][jikl] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf3*NB3 + bf4*NB2 + bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[4][2][jikl] - gaugeERIgradAC_sphswapAB[3][4][jikl]
          - gaugeERIgradAC_sphswapAB[7][1][jikl] +gaugeERIgradAC_sphswapAB[6][3][jikl] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[8][1][jikl] - gaugeERIgradAC_sphswapAB[7][2][jikl]
          - gaugeERIgradAC_sphswapAB[2][4][jikl] +gaugeERIgradAC_sphswapAB[1][5][jikl] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[6][2][jikl] - gaugeERIgradAC_sphswapAB[8][0][jikl]
          - gaugeERIgradAC_sphswapAB[0][5][jikl] +gaugeERIgradAC_sphswapAB[2][2][jikl] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[7][0][jikl] - gaugeERIgradAC_sphswapAB[6][1][jikl]
          - gaugeERIgradAC_sphswapAB[1][2][jikl] +gaugeERIgradAC_sphswapAB[0][4][jikl] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[2][3][jikl] - gaugeERIgradAC_sphswapAB[1][4][jikl]
          - gaugeERIgradAC_sphswapAB[5][1][jikl] +gaugeERIgradAC_sphswapAB[4][2][jikl] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[0][4][jikl] - gaugeERIgradAC_sphswapAB[2][1][jikl]
          - gaugeERIgradAC_sphswapAB[3][2][jikl] +gaugeERIgradAC_sphswapAB[5][0][jikl] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf3*NB3 + bf4*NB2+ bf1* NB+ bf2 ] =
          gaugeERIgradAC_sphswapAB[1][1][jikl] - gaugeERIgradAC_sphswapAB[0][3][jikl]
          - gaugeERIgradAC_sphswapAB[4][0][jikl] +gaugeERIgradAC_sphswapAB[3][1][jikl] ;


          //swap AB


          // sigma1x sigma2x
          gaugeSLSLACxx[bf3*NB3 + bf4*NB2 + bf2* NB+ bf1 ] =
          gaugeERIgradAC_sph[5][4][ijkl] - gaugeERIgradAC_sph[8][3][ijkl]
          - gaugeERIgradAC_sph[3][5][ijkl] +gaugeERIgradAC_sph[6][4][ijkl] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[3][5][ijkl] - gaugeERIgradAC_sph[6][4][ijkl]
          - gaugeERIgradAC_sph[5][2][ijkl] +gaugeERIgradAC_sph[8][1][ijkl] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[4][2][ijkl] - gaugeERIgradAC_sph[3][4][ijkl]
          - gaugeERIgradAC_sph[7][1][ijkl] +gaugeERIgradAC_sph[6][3][ijkl] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[8][1][ijkl] - gaugeERIgradAC_sph[7][2][ijkl]
          - gaugeERIgradAC_sph[2][4][ijkl] +gaugeERIgradAC_sph[1][5][ijkl] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[6][2][ijkl] - gaugeERIgradAC_sph[8][0][ijkl]
          - gaugeERIgradAC_sph[0][5][ijkl] +gaugeERIgradAC_sph[2][2][ijkl] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[7][0][ijkl] - gaugeERIgradAC_sph[6][1][ijkl]
          - gaugeERIgradAC_sph[1][2][ijkl] +gaugeERIgradAC_sph[0][4][ijkl] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[2][3][ijkl] - gaugeERIgradAC_sph[1][4][ijkl]
          - gaugeERIgradAC_sph[5][1][ijkl] +gaugeERIgradAC_sph[4][2][ijkl] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[0][4][ijkl] - gaugeERIgradAC_sph[2][1][ijkl]
          - gaugeERIgradAC_sph[3][2][ijkl] +gaugeERIgradAC_sph[5][0][ijkl] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf3* NB3+ bf4*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sph[1][1][ijkl] - gaugeERIgradAC_sph[0][3][ijkl]
          - gaugeERIgradAC_sph[4][0][ijkl] +gaugeERIgradAC_sph[3][1][ijkl] ;


          //swap CD


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[5][4][jilk] - gaugeERIgradAC_sphswapABCD[8][3][jilk]
          - gaugeERIgradAC_sphswapABCD[3][5][jilk] +gaugeERIgradAC_sphswapABCD[6][4][jilk] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[3][5][jilk] - gaugeERIgradAC_sphswapABCD[6][4][jilk]
          - gaugeERIgradAC_sphswapABCD[5][2][jilk] +gaugeERIgradAC_sphswapABCD[8][1][jilk] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[4][2][jilk] - gaugeERIgradAC_sphswapABCD[3][4][jilk]
          - gaugeERIgradAC_sphswapABCD[7][1][jilk] +gaugeERIgradAC_sphswapABCD[6][3][jilk] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[8][1][jilk] - gaugeERIgradAC_sphswapABCD[7][2][jilk]
          - gaugeERIgradAC_sphswapABCD[2][4][jilk] +gaugeERIgradAC_sphswapABCD[1][5][jilk] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[6][2][jilk] - gaugeERIgradAC_sphswapABCD[8][0][jilk]
          - gaugeERIgradAC_sphswapABCD[0][5][jilk] +gaugeERIgradAC_sphswapABCD[2][2][jilk] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[7][0][jilk] - gaugeERIgradAC_sphswapABCD[6][1][jilk]
          - gaugeERIgradAC_sphswapABCD[1][2][jilk] +gaugeERIgradAC_sphswapABCD[0][4][jilk] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[2][3][jilk] - gaugeERIgradAC_sphswapABCD[1][4][jilk]
          - gaugeERIgradAC_sphswapABCD[5][1][jilk] +gaugeERIgradAC_sphswapABCD[4][2][jilk] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[0][4][jilk] - gaugeERIgradAC_sphswapABCD[2][1][jilk]
          - gaugeERIgradAC_sphswapABCD[3][2][jilk] +gaugeERIgradAC_sphswapABCD[5][0][jilk] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf4* NB3+ bf3*NB2 + bf1*NB + bf2] =
          gaugeERIgradAC_sphswapABCD[1][1][jilk] - gaugeERIgradAC_sphswapABCD[0][3][jilk]
          - gaugeERIgradAC_sphswapABCD[4][0][jilk] +gaugeERIgradAC_sphswapABCD[3][1][jilk] ;


          //swap ABCD


          // sigma1x sigma2x
          gaugeSLSLACxx[ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[5][4][ijlk] - gaugeERIgradAC_sphswapCD[8][3][ijlk]
          - gaugeERIgradAC_sphswapCD[3][5][ijlk] +gaugeERIgradAC_sphswapCD[6][4][ijlk] ;

          // sigma1x sigma2y
          gaugeSLSLACxx[1* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[3][5][ijlk] - gaugeERIgradAC_sphswapCD[6][4][ijlk]
          - gaugeERIgradAC_sphswapCD[5][2][ijlk] +gaugeERIgradAC_sphswapCD[8][1][ijlk] ;

          // sigma1x sigma2z
          gaugeSLSLACxx[2* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[4][2][ijlk] - gaugeERIgradAC_sphswapCD[3][4][ijlk]
          - gaugeERIgradAC_sphswapCD[7][1][ijlk] +gaugeERIgradAC_sphswapCD[6][3][ijlk] ;

          // sigma1y sigma2x
          gaugeSLSLACxx[3* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[8][1][ijlk] - gaugeERIgradAC_sphswapCD[7][2][ijlk]
          - gaugeERIgradAC_sphswapCD[2][4][ijlk] +gaugeERIgradAC_sphswapCD[1][5][ijlk] ;

          // sigma1y sigma2y
          gaugeSLSLACxx[4* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[6][2][ijlk] - gaugeERIgradAC_sphswapCD[8][0][ijlk]
          - gaugeERIgradAC_sphswapCD[0][5][ijlk] +gaugeERIgradAC_sphswapCD[2][2][ijlk] ;

          // sigma1y sigma2z
          gaugeSLSLACxx[5* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[7][0][ijlk] - gaugeERIgradAC_sphswapCD[6][1][ijlk]
          - gaugeERIgradAC_sphswapCD[1][2][ijlk] +gaugeERIgradAC_sphswapCD[0][4][ijlk] ;

          // sigma1z sigma2x
          gaugeSLSLACxx[6* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[2][3][ijlk] - gaugeERIgradAC_sphswapCD[1][4][ijlk]
          - gaugeERIgradAC_sphswapCD[5][1][ijlk] +gaugeERIgradAC_sphswapCD[4][2][ijlk] ;

          // sigma1z sigma2y
          gaugeSLSLACxx[7* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[0][4][ijlk] - gaugeERIgradAC_sphswapCD[2][1][ijlk]
          - gaugeERIgradAC_sphswapCD[3][2][ijlk] +gaugeERIgradAC_sphswapCD[5][0][ijlk] ;

          // sigma1z sigma2z
          gaugeSLSLACxx[8* NB4+ bf4* NB3+ bf3*NB2 + bf2*NB + bf1] =
          gaugeERIgradAC_sphswapCD[1][1][ijlk] - gaugeERIgradAC_sphswapCD[0][3][ijlk]
          - gaugeERIgradAC_sphswapCD[4][0][ijlk] +gaugeERIgradAC_sphswapCD[3][1][ijlk] ;
          */






          /*



          // I2 sigma1x
          gaugeSLSLACxx[NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[6][1][ijkl] + gaugeERIgrad_sph[7][3][ijkl] + gaugeERIgrad_sph[8][4][ijkl]
          -gaugeERIgrad_sph[3][2][ijkl] -gaugeERIgrad_sph[4][4][ijkl] - gaugeERIgrad_sph[5][5][ijkl];

          // I2 sigma1y
          gaugeSLSLACxx[2 * NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[0][2][ijkl] + gaugeERIgrad_sph[1][4][ijkl] + gaugeERIgrad_sph[2][5][ijkl]
          -gaugeERIgrad_sph[6][0][ijkl] -gaugeERIgrad_sph[7][1][ijkl] - gaugeERIgrad_sph[8][2][ijkl];

          // I2 sigma1z
          gaugeSLSLACxx[3* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[3][0][ijkl] + gaugeERIgrad_sph[4][1][ijkl] + gaugeERIgrad_sph[5][2][ijkl]
          -gaugeERIgrad_sph[0][1][ijkl] -gaugeERIgrad_sph[1][3][ijkl] - gaugeERIgrad_sph[2][4][ijkl];

          // I1 sigma2x
          gaugeSLSLACxx[4* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[1][2][ijkl] + gaugeERIgrad_sph[4][4][ijkl] + gaugeERIgrad_sph[7][5][ijkl]
          -gaugeERIgrad_sph[2][1][ijkl] -gaugeERIgrad_sph[5][3][ijkl] - gaugeERIgrad_sph[8][4][ijkl];

          // I1 sigma2y
          gaugeSLSLACxx[5* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[2][0][ijkl] + gaugeERIgrad_sph[5][1][ijkl] + gaugeERIgrad_sph[8][2][ijkl]
          -gaugeERIgrad_sph[0][2][ijkl] -gaugeERIgrad_sph[3][4][ijkl] - gaugeERIgrad_sph[6][5][ijkl];

          // I1 sigma2z
          gaugeSLSLACxx[6* NB4+ bf1* NB3+ bf2*NB2 + bf3*NB + bf4] =
          gaugeERIgrad_sph[0][1][ijkl] + gaugeERIgrad_sph[3][3][ijkl] + gaugeERIgrad_sph[6][4][ijkl]
          -gaugeERIgrad_sph[1][0][ijkl] -gaugeERIgrad_sph[4][1][ijkl] - gaugeERIgrad_sph[7][2][ijkl];




          */


          /*
            //std::cout <<"Libint ∇A∙∇C(ij|kl)"<<std::endl
            if (std::abs(gaugeERIgrad_sph[0][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AxCx = "<<gaugeERIgra d_sph[0][0][ijkl]<<std::endl;
                  }

            if (std::abs(gaugeERIgrad_sph[1][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AxCy = "<<gaugeERIgrad_sph[1][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[2][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AxCz = "<<gaugeERIgrad_sph[2][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[3][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AyCx = "<<gaugeERIgrad_sph[3][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[4][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AyCy = "<<gaugeERIgrad_sph[4][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[5][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AyCz = "<<gaugeERIgrad_sph[5][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[6][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AzCx = "<<gaugeERIgrad_sph[6][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[7][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AzCy = "<<gaugeERIgrad_sph[7][0][ijkl]<<std::endl;
                  }
            if (std::abs(gaugeERIgrad_sph[8][0][ijkl])>1.0E-12 ) {
            std::cout<<"i="<<i<<",j="<<j<<",k="<<k<<",l="<<l<<",AzCz = "<<gaugeERIgrad_sph[8][0][ijkl]<<std::endl;
                  }
          */








        } // for ijkl
      } // for s4
      } // for s3
      } // for s2
      } // for s1

    }; // omp region parallel tested ook

    #ifdef __DEBUGGAUGE__
        // print (SL|SL) sigma1 sigma2
    for (int ii = 0 ; ii < 9 ; ii++ ){
    std::cout << "SLSL gauge Integrals AC component sigma1 sigma2 "<<ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    } // for ii

    std::cout << "SLSL gauge Integrals AC component I1 I2 " << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeSLSLACxxtrue[15*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeSLSLACxxtrue[15*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };



    // print sigma1 I2
    for (int ii = 9 ; ii < 12 ; ii++ ){
    std::cout << "SLSL gauge Integrals AC component sigma1 I2 "<<ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    } // for ii

    // print I1 sigma2
    for (int ii = 12 ; ii < 15 ; ii++ ){
    std::cout << "SLSL gauge Integrals AC component I1 sigma2  "<<ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeSLSLACxxtrue[ii*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    } // for ii
    #endif // debug gauge



    // wrong (LS|SL)
    /*
    for ( int ii = 0 ; ii < 9 ; ii++ ) {
      std::cout << "LSSL gauge Integrals BC component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeSLSLACxx[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeSLSLACxx[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    */


// correct (LS|SL)

    // use symmetry
    for (auto i = 0ul, ijkl = 0ul; i < NB; i++)
    for (auto j = 0ul; j < NB; j++)
    for (auto k = 0ul; k < NB; k++)
    for (auto l = 0ul; l < NB; l++, ++ijkl) {

      size_t jikl = j * NB3 + i * NB2 + k * NB + l;

      // sigma1 sigma2
      for (size_t icomp = 0; icomp < 9; icomp++) {
        gaugeLSSLACsymm[icomp * NB4 + jikl] = gaugeSLSLACxxtrue[icomp * NB4 + ijkl];
        //gaugeLSSLACsymm[ icomp * NB4 + jikl ] = -gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
      }

      // sigma1 I2
      for (size_t icomp = 9; icomp < 12; icomp++) {
        gaugeLSSLACsymm[icomp * NB4 + jikl] = gaugeSLSLACxxtrue[icomp * NB4 + ijkl];
        //gaugeLSSLACsymm[ icomp * NB4 + jikl ] = -gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
      }
      // I1 sigma2
      for (size_t icomp = 12; icomp < 15; icomp++) {
        gaugeLSSLACsymm[icomp * NB4 + jikl] = -gaugeSLSLACxxtrue[icomp * NB4 + ijkl];
        //gaugeLSSLACsymm[ icomp * NB4 + jikl ] = gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
      }

      // I1 I2
      gaugeLSSLACsymm[15 * NB4 + jikl] = -gaugeSLSLACxxtrue[15 * NB4 + ijkl];

    };


    #ifdef __DEBUGGAUGE__

        // print LSSL sigma sigma
    for ( int ii = 0 ; ii < 9 ; ii++ ) {
    std::cout << "LSSL gauge Integrals BC use symmetry component" << ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeLSSLACsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeLSSLACsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    }
    // print sigma1 I2
    for (int ii = 9 ; ii < 12 ; ii++ ){
    std::cout << "LSSL gauge Integrals BC component sigma1 I2 "<<ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeLSSLACsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeLSSLACsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    } // for ii

    // print I1 sigma2
    for (int ii = 12 ; ii < 15 ; ii++ ){
    std::cout << "LSSL gauge Integrals BC component I1 sigma2  "<<ii << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeLSSLACsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeLSSLACsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };
    } // for ii

    std::cout << "LSSL gauge Integrals BC component I1 I2  "<< std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
     if (std::abs( gaugeLSSLACsymm[15*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << gaugeLSSLACsymm[15*NB4 + i* NB3+ j*NB2 + k*NB + l] << std::endl;
     }
    };


    #endif





    // the following part validate the symmetry between SLLS, LSLS integrals, not used in actual calculation
    // so comment out

    /*
    // (SL|LS)
      // use symmetry
      for(auto i = 0ul, ijkl = 0ul ; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++, ++ijkl){
       size_t ijlk = i*NB3 + j*NB2 + l*NB + k;

       for ( size_t icomp = 0 ; icomp < 9 ; icomp++ ) {
         gaugeSLLSsymm[ icomp * NB4 + ijlk ] = gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }

       // sigma1 I2
       for ( size_t icomp = 9 ; icomp < 12 ; icomp++ ) {
         gaugeSLLSsymm[ icomp * NB4 + ijlk ] = -gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }
       // I1 sigma2
       for ( size_t icomp = 12 ; icomp < 15 ; icomp++ ) {
         gaugeSLLSsymm[ icomp * NB4 + ijlk ] = gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }

      };

    for ( int ii = 0 ; ii < 9 ; ii++ ) {
      std::cout << "SLLS gauge Integrals use symmetry component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeSLLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeSLLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    for ( int ii = 9 ; ii < 12 ; ii++ ) {
      std::cout << "SLLS gauge Integrals use symmetry sigma1 I2 component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeSLLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeSLLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    for ( int ii = 12 ; ii < 15 ; ii++ ) {
      std::cout << "SLLS gauge Integrals use symmetry I1 sigma2 component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeSLLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeSLLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }


    // (LS|LS)
      // use symmetry
      for(auto i = 0ul, ijkl = 0ul ; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++, ++ijkl){
       size_t jilk = j*NB3 + i*NB2 + l*NB + k;

       for ( size_t icomp = 0 ; icomp < 9 ; icomp++ ) {
         gaugeLSLSsymm[ icomp * NB4 + jilk ] = gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }

       // sigma1 I2
       for ( size_t icomp = 9 ; icomp < 12 ; icomp++ ) {
         gaugeLSLSsymm[ icomp * NB4 + jilk ] = -gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }
       // I1 sigma2
       for ( size_t icomp = 12 ; icomp < 15 ; icomp++ ) {
         gaugeLSLSsymm[ icomp * NB4 + jilk ] = -gaugeSLSLACxxtrue[ icomp*NB4 + ijkl];
       }

      };

    for ( int ii = 0 ; ii < 9 ; ii++ ) {
      std::cout << "LSLS gauge Integrals use symmetry component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeLSLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeLSLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    for ( int ii = 9 ; ii < 12 ; ii++ ) {
      std::cout << "LSLS gauge Integrals use symmetry sigma1 I2 component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeLSLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeLSLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    for ( int ii = 12 ; ii < 15 ; ii++ ) {
      std::cout << "LSLS gauge Integrals use symmetry I1 sigma2 component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( gaugeLSLSsymm[ii*NB4 + i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << gaugeLSLSsymm[ii*NB4+i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    */



// copy LSSL integral to (*this). Gauge integral start from 23, since Gaunt is 4 to 22.
// 23:ss  24:s sigma x 25: s sigma y 26: s sigma z 27: sigma x s 28: sigma y s 29: sigma z s
// 30: sigma X sigma x 31: sigmaXsigma y 32: sigmaXsigma z 33:

    for (auto bf1 = 0ul; bf1 < NB; bf1++)
    for (auto bf2 = 0ul; bf2 < NB; bf2++)
    for (auto bf3 = 0ul; bf3 < NB; bf3++)
    for (auto bf4 = 0ul; bf4 < NB; bf4++) {

      auto IJKL = bf1 + bf2 * NB + bf3 * NB2 + bf4 * NB3;
      auto ijkl = bf4 + bf3 * NB + bf2 * NB2 + bf1 * NB3;
      auto jikl = bf4 + bf3 * NB + bf1 * NB2 + bf2 * NB3;
      // (ss)
    // LSSL obtained from SLSL
    /*
      (*this)[23].pointer()[IJKL] = -gaugeSLSLACxxtrue[15*NB4 + jikl];

      // (s sigma)x
      (*this)[24].pointer()[IJKL] = -gaugeSLSLACxxtrue[12*NB4 + jikl];

      // (s sigma)y
      (*this)[25].pointer()[IJKL] = -gaugeSLSLACxxtrue[13*NB4 + jikl];

      // (s sigma)z
      (*this)[26].pointer()[IJKL] = -gaugeSLSLACxxtrue[14*NB4 + jikl];

      // (sigma s)x
      (*this)[27].pointer()[IJKL] = gaugeSLSLACxxtrue[9*NB4 + jikl];

      // (sigma s)y
      (*this)[28].pointer()[IJKL] = gaugeSLSLACxxtrue[10*NB4 + jikl];

      // (sigma s)z
      (*this)[29].pointer()[IJKL] = gaugeSLSLACxxtrue[11*NB4 + jikl];

      // (sigma dot sigma)
      (*this)[30].pointer()[IJKL] = gaugeSLSLACxxtrue[ jikl] + gaugeSLSLACxxtrue[4*NB4 + jikl]+gaugeSLSLACxxtrue[8*NB4 + jikl];

      // (sigma Vcross sigma)_x : sigma_y sigma_z -sigma_z sigma_y
      (*this)[31].pointer()[IJKL] = gaugeSLSLACxxtrue[5*NB4 + jikl]-gaugeSLSLACxxtrue[7*NB4 + jikl];

      // (sigma cross sigma)_y : sigma_z sigma_x -sigma_x sigma_z
      (*this)[32].pointer()[IJKL] = gaugeSLSLACxxtrue[6*NB4 + jikl]-gaugeSLSLACxxtrue[2*NB4 + jikl];

      // (sigma cross sigma)_z : sigma_x sigma_y -sigma_y sigma_x
      (*this)[33].pointer()[IJKL] = gaugeSLSLACxxtrue[NB4 + jikl]-gaugeSLSLACxxtrue[3*NB4 + jikl];

      // (sigma sigma) (xx - yy - zz)
      (*this)[34].pointer()[IJKL] = gaugeSLSLACxxtrue[ jikl] - gaugeSLSLACxxtrue[4*NB4 + jikl] - gaugeSLSLACxxtrue[8*NB4 + jikl];

      // (sigma sigma) (-xx + yy - zz)
      (*this)[35].pointer()[IJKL] =-gaugeSLSLACxxtrue[ jikl] + gaugeSLSLACxxtrue[4*NB4 + jikl] - gaugeSLSLACxxtrue[8*NB4 + jikl];

      // (sigma sigma) (-xx - yy + zz)
      (*this)[36].pointer()[IJKL] =-gaugeSLSLACxxtrue[ jikl] - gaugeSLSLACxxtrue[4*NB4 + jikl] + gaugeSLSLACxxtrue[8*NB4 + jikl];


      //  sigma_x sigma_y +sigma_y sigma_x
      (*this)[37].pointer()[IJKL] = gaugeSLSLACxxtrue[NB4 + jikl]+gaugeSLSLACxxtrue[3*NB4 + jikl];

      //  sigma_z sigma_x +sigma_x sigma_z
      (*this)[38].pointer()[IJKL] = gaugeSLSLACxxtrue[6*NB4 + jikl]+gaugeSLSLACxxtrue[2*NB4 + jikl];


      //  sigma_y sigma_z +sigma_z sigma_y
      (*this)[39].pointer()[IJKL] = gaugeSLSLACxxtrue[5*NB4 + jikl]+gaugeSLSLACxxtrue[7*NB4 + jikl];

      // sigma_x sigma_x
      (*this)[40].pointer()[IJKL] =gaugeSLSLACxxtrue[jikl];
      // sigma x sigma y
      (*this)[41].pointer()[IJKL] =gaugeSLSLACxxtrue[NB4 + jikl];
      // sigma x sigma z
      (*this)[42].pointer()[IJKL] =gaugeSLSLACxxtrue[2*NB4 + jikl];
      // sigma y sigma x
      (*this)[43].pointer()[IJKL] =gaugeSLSLACxxtrue[3*NB4 + jikl];
      // sigma y sigma y
      (*this)[44].pointer()[IJKL] =gaugeSLSLACxxtrue[4*NB4 + jikl];
      // sigma y sigma z
      (*this)[45].pointer()[IJKL] =gaugeSLSLACxxtrue[5*NB4 + jikl];
      // sigma z sigma x
      (*this)[46].pointer()[IJKL] =gaugeSLSLACxxtrue[6*NB4 + jikl];
      // sigma z sigma y
      (*this)[47].pointer()[IJKL] =gaugeSLSLACxxtrue[7*NB4 + jikl];
      // sigma z sigma z
      (*this)[48].pointer()[IJKL] =gaugeSLSLACxxtrue[8*NB4 + jikl];
    */

      // (ss)
      (*this)[nERIRef + 0].pointer()[IJKL] = gaugeLSSLACsymm[15 * NB4 + ijkl];

      // (sσ)_x
      (*this)[nERIRef + 1].pointer()[IJKL] = gaugeLSSLACsymm[12 * NB4 + ijkl];

      // (sσ)_y
      (*this)[nERIRef + 2].pointer()[IJKL] = gaugeLSSLACsymm[13 * NB4 + ijkl];

      // (sσ)_z
      (*this)[nERIRef + 3].pointer()[IJKL] = gaugeLSSLACsymm[14 * NB4 + ijkl];

      // (σs)_x
      (*this)[nERIRef + 4].pointer()[IJKL] = gaugeLSSLACsymm[9 * NB4 + ijkl];

      // (σs)_y
      (*this)[nERIRef + 5].pointer()[IJKL] = gaugeLSSLACsymm[10 * NB4 + ijkl];

      // (σs)_z
      (*this)[nERIRef + 6].pointer()[IJKL] = gaugeLSSLACsymm[11 * NB4 + ijkl];

      // (σ∙σ)
      (*this)[nERIRef + 7].pointer()[IJKL] =
        gaugeLSSLACsymm[ijkl] + gaugeLSSLACsymm[4 * NB4 + ijkl] + gaugeLSSLACsymm[8 * NB4 + ijkl];

      // (σxσ)_x = σ_y σ_z - σ_z σ_y
      (*this)[nERIRef + 8].pointer()[IJKL] = gaugeLSSLACsymm[5 * NB4 + ijkl] - gaugeLSSLACsymm[7 * NB4 + ijkl];

      // (σxσ)_y = σ_z σ_x - σ_x σ_z
      (*this)[nERIRef + 9].pointer()[IJKL] = gaugeLSSLACsymm[6 * NB4 + ijkl] - gaugeLSSLACsymm[2 * NB4 + ijkl];

      // (σxσ)_z = σ_x σ_y - σ_y σ_x
      (*this)[nERIRef + 10].pointer()[IJKL] = gaugeLSSLACsymm[NB4 + ijkl] - gaugeLSSLACsymm[3 * NB4 + ijkl];

      // (σσ) (xx - yy - zz)
      (*this)[nERIRef + 11].pointer()[IJKL] =
        gaugeLSSLACsymm[ijkl] - gaugeLSSLACsymm[4 * NB4 + ijkl] - gaugeLSSLACsymm[8 * NB4 + ijkl];

      // (σσ) (-xx + yy - zz)
      (*this)[nERIRef + 12].pointer()[IJKL] =
        -gaugeLSSLACsymm[ijkl] + gaugeLSSLACsymm[4 * NB4 + ijkl] - gaugeLSSLACsymm[8 * NB4 + ijkl];

      // (σσ) (-xx - yy + zz)
      (*this)[nERIRef + 13].pointer()[IJKL] =
        -gaugeLSSLACsymm[ijkl] - gaugeLSSLACsymm[4 * NB4 + ijkl] + gaugeLSSLACsymm[8 * NB4 + ijkl];

      //  σ_x σ_y + σ_y σ_x
      (*this)[nERIRef + 14].pointer()[IJKL] = gaugeLSSLACsymm[NB4 + ijkl] + gaugeLSSLACsymm[3 * NB4 + ijkl];

      //  σ_z σ_x + σ_x σ_z
      (*this)[nERIRef + 15].pointer()[IJKL] = gaugeLSSLACsymm[6 * NB4 + ijkl] + gaugeLSSLACsymm[2 * NB4 + ijkl];

      //  σ_y σ_z + σ_z σ_y
      (*this)[nERIRef + 16].pointer()[IJKL] = gaugeLSSLACsymm[5 * NB4 + ijkl] + gaugeLSSLACsymm[7 * NB4 + ijkl];

      // σ_x σ_x
      (*this)[nERIRef + 17].pointer()[IJKL] = gaugeLSSLACsymm[ijkl];

      // σ_x σ_y
      (*this)[nERIRef + 18].pointer()[IJKL] = gaugeLSSLACsymm[NB4 + ijkl];

      // σ_x σ_z
      (*this)[nERIRef + 19].pointer()[IJKL] = gaugeLSSLACsymm[2 * NB4 + ijkl];

      // σ_y σ_x
      (*this)[nERIRef + 20].pointer()[IJKL] = gaugeLSSLACsymm[3 * NB4 + ijkl];

      // σ_y σ_y
      (*this)[nERIRef + 21].pointer()[IJKL] = gaugeLSSLACsymm[4 * NB4 + ijkl];

      // σ_y σ_z
      (*this)[nERIRef + 22].pointer()[IJKL] = gaugeLSSLACsymm[5 * NB4 + ijkl];

      // σ_z σ_x
      (*this)[nERIRef + 23].pointer()[IJKL] = gaugeLSSLACsymm[6 * NB4 + ijkl];

      // σ_z σ_y
      (*this)[nERIRef + 24].pointer()[IJKL] = gaugeLSSLACsymm[7 * NB4 + ijkl];

      // σ_z σ_z
      (*this)[nERIRef + 25].pointer()[IJKL] = gaugeLSSLACsymm[8 * NB4 + ijkl];


    }

    /*
      int icount;
      for (InCore4indexERI<double>& c : components_)
        icount += 1;

    std::cout<<"icount= "<<icount<<std::endl;

    for ( int ii = 23 ; ii < 48 ; ii++ ) {
      std::cout << "LSSL gauge Integrals component" << ii << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
       if (std::abs( (*this)[ii].pointer()[ i* NB3+ j*NB2 + k*NB + l])>1.0e-12 ) {
        std::cout <<" comp "<<ii<< "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout <<(*this)[ii].pointer()[ i* NB3+ j*NB2 + k*NB + l] << std::endl;
       }
      };
    }
    */

    // SLLS
    for(auto index=0; index<26; index++){
      std::cout << "InHouse ("<<index<<")(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+index](i, j, k, l) << std::endl;
      };
    };



    CQMemManager::get().free<double>(gaugeLSSLACsymm);
    CQMemManager::get().free<double>(gaugeSLLSsymm);
    CQMemManager::get().free<double>(gaugeLSLSsymm);

// deallocate SLSL after assigning the integrals
    //CQMemManager::get().free<double>(gaugeSLSLACxx);



    auto durERIGauge = tock(topERIGauge);
    std::cout << "In-House-ERI-Gauge duration   = " << durERIGauge << std::endl;


  } // InCore4indexRelERI<double>::computeERIGauge












  template <>
  void GTODirectRelERIContraction<double,double>::computeERI3Index(size_t s1) {
    CErr("Only complex WFN is allowed",std::cout);
  };

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::computeERI3Index(size_t s1) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void GTODirectRelERIContraction<dcomplex,double>::computeERI3Index(size_t s1){
    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    BasisSet &basisSet_ = std::dynamic_pointer_cast<DirectTPI<double>>(ints())->basisSet();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      basisSet_.maxPrim,basisSet_.maxL,2);
    engines[0].set_precision(0.);


    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;

    // There are 78 2nd Derivatives per ERI
    int AxBx = 3;
    int AxBy = 4;
    int AxBz = 5;
    int AyBx = 14;
    int AyBy = 15;
    int AyBz = 16;
    int AzBx = 24;
    int AzBy = 25;
    int AzBz = 26;

    int CxDx = 60;
    int CxDy = 61;
    int CxDz = 62;
    int CyDx = 65;
    int CyDy = 66;
    int CyDz = 67;
    int CzDx = 69;
    int CzDy = 70;
    int CzDz = 71;

    int AxCx = 6;
    int AxCy = 7;
    int AxCz = 8;
    int AyCx = 17;
    int AyCy = 18;
    int AyCz = 19;
    int AzCx = 27;
    int AzCy = 28;
    int AzCz = 29;

    int AxDx = 9;
    int AxDy = 10;
    int AxDz = 11;
    int AyDx = 20;
    int AyDy = 21;
    int AyDz = 22;
    int AzDx = 30;
    int AzDy = 31;
    int AzDz = 32;

    int BxCx = 36;
    int BxCy = 37;
    int BxCz = 38;
    int ByCx = 44;
    int ByCy = 45;
    int ByCz = 46;
    int BzCx = 51;
    int BzCy = 52;
    int BzCz = 53;

    int BxDx = 39;
    int BxDy = 40;
    int BxDz = 41;
    int ByDx = 47;
    int ByDy = 48;
    int ByDz = 49;
    int BzDx = 54;
    int BzDy = 55;
    int BzDz = 56;

    auto n1 = basisSet_.shells[s1].size();
    size_t nERI3 = 37;

    if(ERI4DCB!=nullptr)
      CQMemManager::get().free(ERI4DCB);

    try { ERI4DCB = CQMemManager::get().malloc<double>(nERI3*NB3*n1); }
    catch(...) {
      std::cout << std::fixed;
      std::cout << "Insufficient memory for the full ERI Dirac-Coulomb-Gaunt tensor ("
                << (nERI3*NB3*n1/1e9) * sizeof(double) << " GB)" << std::endl;
      std::cout << std::endl << CQMemManager::get() << std::endl;
      CErr();
    }
    memset(ERI4DCB, 0.,nERI3*NB3*n1*sizeof(double));

    auto topERIDCB = tick();

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;

      for(size_t s2(0), bf2_s(0), s1234(0); s2 < basisSet_.nShell; bf2_s+=n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for(size_t s3(0), bf3_s(0); s3 < basisSet_.nShell; bf3_s+=n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3

      for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        // Evaluate ERI for shell quartet
        engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 2>(
          basisSet_.shells[s1],
          basisSet_.shells[s2],
          basisSet_.shells[s3],
          basisSet_.shells[s4]
        );
        const auto& buff = engines[thread_id].results();

        // Place shell quartet into persistent storage with
        // permutational symmetry
        for(i = 0ul, bf1 = 0ul, ijkl = 0ul ; i < n1; ++i, bf1++)
        for(j = 0ul, bf2 = bf2_s           ; j < n2; ++j, bf2++)
        for(k = 0ul, bf3 = bf3_s           ; k < n3; ++k, bf3++)
        for(l = 0ul, bf4 = bf4_s           ; l < n4; ++l, bf4++, ++ijkl) {


          auto IJKL =          bf2 + bf3*NB + bf4*NB2    + bf1*nERI3*NB3;
          auto IJLK =          bf2 + bf4*NB + bf3*NB2    + bf1*nERI3*NB3;
          auto JIKL = bf2 +          bf3*NB + bf4*NB2    + bf1*nERI3*NB3;
          auto JILK = bf2 +          bf4*NB + bf3*NB2    + bf1*nERI3*NB3;
          auto KLIJ = bf3 + bf4*NB +          bf2*NB2    + bf1*nERI3*NB3;
          auto LKIJ = bf4 + bf3*NB +          bf2*NB2    + bf1*nERI3*NB3;
          auto KLJI = bf3 + bf4*NB + bf2*NB2             + bf1*nERI3*NB3;
          auto LKJI = bf4 + bf3*NB + bf2*NB2             + bf1*nERI3*NB3;
#if 1
          /* Coulomb */
          // ∇A∙∇B(ij|kl)
          auto dAdotdB = buff[AxBx][ijkl] + buff[AyBy][ijkl] + buff[AzBz][ijkl];
          // ∇Ax∇B(ijkl)
          auto dAcrossdB_x =  buff[AyBz][ijkl] - buff[AzBy][ijkl];
          auto dAcrossdB_y = -buff[AxBz][ijkl] + buff[AzBx][ijkl];
          auto dAcrossdB_z =  buff[AxBy][ijkl] - buff[AyBx][ijkl];

          // ∇C∙∇D(ij|kl)
          auto dCdotdD = buff[CxDx][ijkl] + buff[CyDy][ijkl] + buff[CzDz][ijkl];
          // ∇Cx∇D(ijkl)
          auto dCcrossdD_x =  buff[CyDz][ijkl] - buff[CzDy][ijkl];
          auto dCcrossdD_y = -buff[CxDz][ijkl] + buff[CzDx][ijkl];
          auto dCcrossdD_z =  buff[CxDy][ijkl] - buff[CyDx][ijkl];

          // ∇A∙∇B(ij|kl) followed by ∇Ax∇B(ij|kl) X, Y, and Z
          // (ij|kl)
          ERI4DCB[        IJKL] = dAdotdB;
          ERI4DCB[1*NB3 + IJKL] = dAcrossdB_x;
          ERI4DCB[2*NB3 + IJKL] = dAcrossdB_y;
          ERI4DCB[3*NB3 + IJKL] = dAcrossdB_z;
          // (ij|lk)
          ERI4DCB[        IJLK] = dAdotdB;
          ERI4DCB[1*NB3 + IJLK] = dAcrossdB_x;
          ERI4DCB[2*NB3 + IJLK] = dAcrossdB_y;
          ERI4DCB[3*NB3 + IJLK] = dAcrossdB_z;

          // ∇C∙∇D(ij|kl) followed by ∇Cx∇D(ij|kl) X, Y, and Z
          // (ij|kl)
          ERI4DCB[4*NB3 + IJKL] = dCdotdD;
          ERI4DCB[5*NB3 + IJKL] = dCcrossdD_x;
          ERI4DCB[6*NB3 + IJKL] = dCcrossdD_y;
          ERI4DCB[7*NB3 + IJKL] = dCcrossdD_z;
          // (ij|lk)
          ERI4DCB[4*NB3 + IJLK] = dCdotdD;
          ERI4DCB[5*NB3 + IJLK] = -dCcrossdD_x;
          ERI4DCB[6*NB3 + IJLK] = -dCcrossdD_y;
          ERI4DCB[7*NB3 + IJLK] = -dCcrossdD_z;
#endif

#if 1
          /* Gaunt */
          // ∇A∙∇C(ij|kl)
          auto dAdotdC = buff[AxCx][ijkl] + buff[AyCy][ijkl] + buff[AzCz][ijkl];
          // ∇Ax∇C(ijkl)
          auto dAcrossdC_x =  buff[AyCz][ijkl] - buff[AzCy][ijkl];
          auto dAcrossdC_y = -buff[AxCz][ijkl] + buff[AzCx][ijkl];
          auto dAcrossdC_z =  buff[AxCy][ijkl] - buff[AyCx][ijkl];

          // ∇A∙∇D(ij|kl)
          auto dAdotdD = buff[AxDx][ijkl] + buff[AyDy][ijkl] + buff[AzDz][ijkl];
          // ∇Ax∇D(ijkl)
          auto dAcrossdD_x =  buff[AyDz][ijkl] - buff[AzDy][ijkl];
          auto dAcrossdD_y = -buff[AxDz][ijkl] + buff[AzDx][ijkl];
          auto dAcrossdD_z =  buff[AxDy][ijkl] - buff[AyDx][ijkl];

          // ∇B∙∇C(ij|kl)
          auto dBdotdC = buff[BxCx][ijkl]+ buff[ByCy][ijkl]+ buff[BzCz][ijkl];
          // ∇Bx∇C(ijkl)
          auto dBcrossdC_x =  buff[ByCz][ijkl] - buff[BzCy][ijkl];
          auto dBcrossdC_y = -buff[BxCz][ijkl] + buff[BzCx][ijkl];
          auto dBcrossdC_z =  buff[BxCy][ijkl] - buff[ByCx][ijkl];

          // ∇B∙∇D(ij|kl)
          auto dBdotdD = buff[BxDx][ijkl] + buff[ByDy][ijkl] + buff[BzDz][ijkl];
          // ∇Bx∇D(ijkl)
          auto dBcrossdD_x =  buff[ByDz][ijkl] - buff[BzDy][ijkl];
          auto dBcrossdD_y = -buff[BxDz][ijkl] + buff[BzDx][ijkl];
          auto dBcrossdD_z =  buff[BxDy][ijkl] - buff[ByDx][ijkl];

          // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
          // (ij|kl)
          ERI4DCB[ 8*NB3 + IJKL] = dBdotdC;
          ERI4DCB[ 9*NB3 + IJKL] = dBcrossdC_x;
          ERI4DCB[10*NB3 + IJKL] = dBcrossdC_y;
          ERI4DCB[11*NB3 + IJKL] = dBcrossdC_z;
          // (ij|lk)
          ERI4DCB[ 8*NB3 + IJLK] = dBdotdD;
          ERI4DCB[ 9*NB3 + IJLK] = dBcrossdD_x;
          ERI4DCB[10*NB3 + IJLK] = dBcrossdD_y;
          ERI4DCB[11*NB3 + IJLK] = dBcrossdD_z;

          // ∇A∙∇D(ij|kl) followed by ∇Ax∇D(ij|kl) X, Y, and Z
          // (ij|kl)
          ERI4DCB[27*NB3 + IJKL] = dAdotdD;
          ERI4DCB[28*NB3 + IJKL] = dAcrossdD_x;
          ERI4DCB[29*NB3 + IJKL] = dAcrossdD_y;
          ERI4DCB[30*NB3 + IJKL] = dAcrossdD_z;
          // (ij|lk)
          ERI4DCB[27*NB3 + IJLK] = dAdotdC;
          ERI4DCB[28*NB3 + IJLK] = dAcrossdC_x;
          ERI4DCB[29*NB3 + IJLK] = dAcrossdC_y;
          ERI4DCB[30*NB3 + IJLK] = dAcrossdC_z;


          // ∇B_x∇C_y(ij|kl) + ∇B_y∇C_x(ij|kl)
          // (ij|kl)
          ERI4DCB[12*NB3 + IJKL] = buff[BxCy][ijkl] + buff[ByCx][ijkl];
          // (ij|lk)
          ERI4DCB[12*NB3 + IJLK] = buff[BxDy][ijkl] + buff[ByDx][ijkl];

          // ∇A_x∇D_y(ij|kl) + ∇A_y∇D_x(ij|kl)
          // (ij|kl)
          ERI4DCB[31*NB3 + IJKL] = buff[AxDy][ijkl] + buff[AyDx][ijkl];
          // (ij|lk)
          ERI4DCB[31*NB3 + IJLK] = buff[AxCy][ijkl] + buff[AyCx][ijkl];


          // ∇B_y∇C_x(ij|kl)
          // (ij|kl)
          ERI4DCB[13*NB3 + IJKL] = buff[ByCx][ijkl];
          // (ij|lk)
          ERI4DCB[13*NB3 + IJLK] = buff[ByDx][ijkl];


          // ∇B_x∇C_z(ij|kl) + ∇B_z∇C_x(ij|kl)
          // (ij|kl)
          ERI4DCB[14*NB3 + IJKL] = buff[BxCz][ijkl] + buff[BzCx][ijkl];
          // (ij|lk)
          ERI4DCB[14*NB3 + IJLK] = buff[BxDz][ijkl] + buff[BzDx][ijkl];

          // ∇A_x∇D_z(ij|kl) + ∇A_z∇D_x(ij|kl)
          // (ij|kl)
          ERI4DCB[33*NB3 + IJKL] = buff[AxDz][ijkl] + buff[AzDx][ijkl];
          // (ij|lk)
          ERI4DCB[33*NB3 + IJLK] = buff[AxCz][ijkl] + buff[AzCx][ijkl];


          // ∇B_z∇C_x(ij|kl)
          // (ij|kl)
          ERI4DCB[15*NB3 + IJKL] = buff[BzCx][ijkl];
          // (ij|lk)
          ERI4DCB[15*NB3 + IJLK] = buff[BzDx][ijkl];


          // ∇B_y∇C_z(ij|kl) + ∇B_z∇C_y(ij|kl)
          // (ij|kl)
          ERI4DCB[16*NB3 + IJKL] = buff[ByCz][ijkl] + buff[BzCy][ijkl];
          // (ij|lk)
          ERI4DCB[16*NB3 + IJLK] = buff[ByDz][ijkl] + buff[BzDy][ijkl];

          // ∇A_y∇D_z(ij|kl) + ∇A_z∇D_y(ij|kl)
          // (ij|kl)
          ERI4DCB[32*NB3 + IJKL] = buff[AyDz][ijkl] + buff[AzDy][ijkl];
          // (ij|lk)
          ERI4DCB[32*NB3 + IJLK] = buff[AyCz][ijkl] + buff[AzCy][ijkl];


          // ∇B_z∇C_y(ij|kl)
          // (ij|kl)
          ERI4DCB[17*NB3 + IJKL] = buff[BzCy][ijkl];
          // (ij|lk)
          ERI4DCB[17*NB3 + IJLK] = buff[BzDy][ijkl];


          // - ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) + ∇B_z∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[18*NB3 + IJKL] = - buff[BxCx][ijkl] - buff[ByCy][ijkl] + buff[BzCz][ijkl];
          // (ij|lk)
          ERI4DCB[18*NB3 + IJLK] = - buff[BxDx][ijkl] - buff[ByDy][ijkl] + buff[BzDz][ijkl];

          // - ∇A_x∇D_x(ij|kl) - ∇A_y∇D_y(ij|kl) + ∇A_z∇D_z(ij|kl)
          // (ij|kl)
          ERI4DCB[34*NB3 + IJKL] = - buff[AxDx][ijkl] - buff[AyDy][ijkl] + buff[AzDz][ijkl];
          // (ij|lk)
          ERI4DCB[34*NB3 + IJLK] = - buff[AxCx][ijkl] - buff[AyCy][ijkl] + buff[AzCz][ijkl];



          // ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[19*NB3 + IJKL] = buff[BxCx][ijkl] - buff[ByCy][ijkl] - buff[BzCz][ijkl];
          // (ij|lk)
          ERI4DCB[19*NB3 + IJLK] = buff[BxDx][ijkl] - buff[ByDy][ijkl] - buff[BzDz][ijkl];

          // ∇A_x∇D_x(ij|kl) - ∇A_y∇D_y(ij|kl) - ∇A_z∇D_z(ij|kl)
          // (ij|kl)
          ERI4DCB[35*NB3 + IJKL] = buff[AxDx][ijkl] - buff[AyDy][ijkl] - buff[AzDz][ijkl];
          // (ij|lk)
          ERI4DCB[35*NB3 + IJLK] = buff[AxCx][ijkl] - buff[AyCy][ijkl] - buff[AzCz][ijkl];



          // - ∇B_x∇C_x(ij|kl) + ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[20*NB3 + IJKL] = - buff[BxCx][ijkl] + buff[ByCy][ijkl] - buff[BzCz][ijkl];
          // (ij|lk)
          ERI4DCB[20*NB3 + IJLK] = - buff[BxDx][ijkl] + buff[ByDy][ijkl] - buff[BzDz][ijkl];

          // - ∇A_x∇D_x(ij|kl) + ∇A_y∇D_y(ij|kl) - ∇A_z∇D_z(ij|kl)
          // (ij|kl)
          ERI4DCB[36*NB3 + IJKL] = - buff[AxDx][ijkl] + buff[AyDy][ijkl] - buff[AzDz][ijkl];
          // (ij|lk)
          ERI4DCB[36*NB3 + IJLK] = - buff[AxCx][ijkl] + buff[AyCy][ijkl] - buff[AzCz][ijkl];



          // ∇B_x∇C_x(ij|kl)
          // (ij|kl)
          ERI4DCB[21*NB3 + IJKL] = buff[BxCx][ijkl];
          // (ij|lk)
          ERI4DCB[21*NB3 + IJLK] = buff[BxDx][ijkl];


          // ∇B_x∇C_y(ij|kl)
          // (ij|kl)
          ERI4DCB[22*NB3 + IJKL] = buff[BxCy][ijkl];
          // (ij|lk)
          ERI4DCB[22*NB3 + IJLK] = buff[BxDy][ijkl];

          // ∇B_x∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[23*NB3 + IJKL] = buff[BxCz][ijkl];
          // (ij|lk)
          ERI4DCB[23*NB3 + IJLK] = buff[BxDz][ijkl];

          // ∇B_y∇C_y(ij|kl)
          // (ij|kl)
          ERI4DCB[24*NB3 + IJKL] = buff[ByCy][ijkl];
          // (ij|lk)
          ERI4DCB[24*NB3 + IJLK] = buff[ByDy][ijkl];

          // ∇B_y∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[25*NB3 + IJKL] = buff[ByCz][ijkl];
          // (ij|lk)
          ERI4DCB[25*NB3 + IJLK] = buff[ByDz][ijkl];

          // ∇B_z∇C_z(ij|kl)
          // (ij|kl)
          ERI4DCB[26*NB3 + IJKL] = buff[BzCz][ijkl];
          // (ij|lk)
          ERI4DCB[26*NB3 + IJLK] = buff[BzDz][ijkl];

#endif


        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
    }; // omp region


    auto durERIDCB = tock(topERIDCB);
    //std::cout << "  Libint-ERI4-Dirac-Coulomb-Breit duration   = " << durERIDCB << std::endl;

  }  // computeERI3Index


  // Gradient integrals
  template<>
  void GradInts<TwoPInts,double>::computeAOInts(BasisSet& basisSet,
    BasisSet& basisSet2, Molecule& mol, EMPerturbation& pert, OPERATOR op,
    const HamiltonianOptions &options)
  {

    if (std::dynamic_pointer_cast<DirectTPI<double>>(components_[0]))
      return;

    if (std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(components_[0]) or
        std::dynamic_pointer_cast<InCoreAuxBasisRIERI<double>>(components_[0]))
      CErr("Gradients using RI integrals NYI!");

    // Get vector of internal storages
    std::vector<double*> eris;
    std::transform(components_.begin(), components_.end(),
      std::back_inserter(eris),
      [](std::shared_ptr<TwoPInts<double>>& p) {
        return std::dynamic_pointer_cast<InCore4indexTPI<double>>(p)->pointer();
      }
    );

    size_t NB = basisSet.nBasis;
    size_t MB = basisSet2.nBasis;
    size_t NB2 = NB * MB;
    size_t NB3 = NB2 * MB;
    bool sameBasis = (&basisSet == &basisSet2);



    auto computeAndPlace = [&](size_t sh1, size_t sh2, size_t sh3, size_t sh4,
                               size_t b1s, size_t b2s, size_t b3s, size_t b4s,
                               size_t  n1, size_t  n2, size_t  n3, size_t  n4,
                               libint2::Engine& engine) {

      engine.compute2<
        libint2::Operator::coulomb, libint2::BraKet::xx_xx, 1>(
        basisSet.shells[sh1],
        basisSet.shells[sh2],
        basisSet2.shells[sh3],
        basisSet2.shells[sh4]
      );

      const auto& results = engine.results();

      // Get atomic centers for each basis function
      std::vector<size_t> ac;
      ac.push_back(basisSet.mapSh2Cen[sh1]);
      ac.push_back(basisSet.mapSh2Cen[sh2]);
      ac.push_back(basisSet2.mapSh2Cen[sh3]);
      ac.push_back(basisSet2.mapSh2Cen[sh4]);

      // Get constants
      size_t NB = basisSet.nBasis;
      size_t MB = basisSet2.nBasis;
      size_t NB2 = NB * NB;
      size_t NB3 = NB2 * MB;

      // Place shell quartet into persistent storage with
      // permutational symmetry
      for ( auto iC = 0, itot = 0; iC < 4; iC++ )
      for ( auto iXYZ = 0; iXYZ < 3; iXYZ++, itot++) {

        // Libint internal screening
        if ( results[itot] == nullptr ) continue;

        for(size_t i = 0ul, bf1 = b1s, ijkl = 0ul ; i < n1; ++i, bf1++) 
        for(size_t j = 0ul, bf2 = b2s             ; j < n2; ++j, bf2++) 
        for(size_t k = 0ul, bf3 = b3s             ; k < n3; ++k, bf3++) 
        for(size_t l = 0ul, bf4 = b4s             ; l < n4; ++l, bf4++, ++ijkl) {


          // Because this is just placement, not computation, the if statements
          //   in this hot loop shouldn't significantly increase time, but we
          //   should measure it eventually. Can refactor into if statements on
          //   the outside at the cost of ugly, duplicated code.
          //
          // This redundancy checking is here because of the required summation
          //   instead of simple assignment in the non-gradient ERI case

          // (12 | 34)
          eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf3*NB2 + bf4*NB3] += results[itot][ijkl];

          if ( sh3 != sh4 )
            // (12 | 43)
            eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf4*NB2 + bf3*NB3] += results[itot][ijkl];

          if ( sh1 != sh2 ) {
            // (21 | 34)
            eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf3*NB2 + bf4*NB3] += results[itot][ijkl];
            if ( sh3 != sh4 )
              // (21 | 43)
              eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf4*NB2 + bf3*NB3] += results[itot][ijkl];
          } // sh1/2
          
          if( sameBasis ) {
            if ( sh1 != sh3 || sh2 != sh4 ) {
              // (34 | 12)
              eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf1*NB2 + bf2*NB3] += results[itot][ijkl];

              if ( sh3 != sh4 )
                // (43 | 12)
                eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf1*NB2 + bf2*NB3] += results[itot][ijkl];

              if ( sh1 != sh2 ) {
                // (34 | 21)
                eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf2*NB2 + bf1*NB3] += results[itot][ijkl];

                if ( sh3 != sh4 )
                  // (43 | 21)
                  eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf2*NB2 + bf1*NB3] += results[itot][ijkl];
              } // sh1/2
            } // sh1/3 or sh2/4
          }

        } // ijkl loop

      } // ic/xyz loop

    };

    doERIInCore(1, eris, computeAndPlace, basisSet, basisSet2, op, options);

#ifdef __DEBUGERI__
    std::cout << "Two-Electron Integral Derivatives (ERIs)" << std::endl;
    for(auto iGrad = 0; iGrad < eris.size(); iGrad++ )
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << eris[iGrad][i + j*NB  + k*NB2 + l*NB3] << std::endl;
    };
#endif


  };

  template<>
  void GradInts<TwoPInts,double>::computeAOInts(BasisSet& basisSet,
    Molecule& mol, EMPerturbation& pert, OPERATOR op,
    const HamiltonianOptions &options)
  {
    computeAOInts(basisSet, basisSet, mol, pert, op, options);
  }

  template<>
  void GradInts<TwoPInts,dcomplex>::computeAOInts(BasisSet& basisSet,
    Molecule& mol, EMPerturbation& pert, OPERATOR op,
    const HamiltonianOptions &options)
  {
    CErr("Complex integral gradients not yet implemented!");
  }

  template<>
  void GradInts<TwoPInts,dcomplex>::computeAOInts(BasisSet& basisSet,
    BasisSet& basis2, Molecule& mol, EMPerturbation& pert, OPERATOR op,
    const HamiltonianOptions &options)
  {
    CErr("Complex integral gradients not yet implemented!");
  }


}; // namespace ChronusQ

