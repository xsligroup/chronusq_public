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

#include <cqlinalg.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/onepints.hpp>
#include <integrals.hpp>

#include <util/threads.hpp>
#include <chrono>
#include <utility>
#include <set>
#include <algorithm>

#include <libcint.hpp>

//#define __DEBUGERI__
//#define CD_PROGRESS


namespace ChronusQ {

  inline size_t compoundToSquare(size_t I, size_t NB) {
    size_t p,q;
    q = static_cast<size_t>(sqrt(2*I + 0.25) - 0.5);
    p = I - q * (q+1) / 2;
    return p + q * NB;
  }

  inline size_t squareToCompound(size_t I, size_t NB) {
    size_t p,q;
    p = I % NB;
    q = I / NB;
    return p + q * (q+1) / 2;
  }

  inline size_t toSquare(size_t p, size_t q, size_t NB) {
    return p + q * NB;
  }

  inline std::pair<size_t, size_t> anaSquare(size_t I, size_t NB) {
    size_t p,q;
    p = I % NB;
    q = I / NB;
    return std::make_pair(p,q);
  }

  inline size_t toCompound(size_t p, size_t q) {
    return p + q * (q+1) / 2;
  }

  inline std::pair<size_t, size_t> anaCompound(size_t I) {
    size_t p,q;
    q = static_cast<size_t>(sqrt(2*I + 0.25) - 0.5);
    p = I - q * (q+1) / 2;
    return std::make_pair(p,q);
  }

  template <typename IntsT>
  void InCoreRITPI<IntsT>::saveRawERI3J() {
    if (not saveRawERI_) return;

    size_t NB3 = NBRI * this->nBasis()*(this->nBasis() + 1) / 2;
    if (rawERI3J_) {
      if (CQMemManager::get().getSize(rawERI3J_) == NB3)
        return;
      CQMemManager::get().free(rawERI3J_);
    }
    try { rawERI3J_ = CQMemManager::get().calloc<IntsT>(NB3); }
    catch(...) {
      std::cout << std::fixed;
      std::cout << "Insufficient memory for the full RI-ERI tensor ("
      << (NB3/1e9) * sizeof(double) << " GB)" << std::endl;
      std::cout << std::endl << CQMemManager::get() << std::endl;
      CErr();
    }
    std::copy_n(ERI3J, NB3, rawERI3J_);
  }

  /**
   *  \brief Construct L=(P|Q)^{-1/2}
   */
  template <>
  void InCoreRITPI<dcomplex>::halfInverse2CenterERI(cqmatrix::Matrix<dcomplex> &S) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreRITPI<double>::halfInverse2CenterERI(cqmatrix::Matrix<double> &twocenterERI) {

    auto topERI3Trans = tick();

    double *S = twocenterERI.pointer();
    size_t NBRI = twocenterERI.dimension();

#ifdef __DEBUGERI__
    prettyPrintSmart(std::cout, "Cholesky S", S, NBRI, NBRI, NBRI);
#endif

    // S = L L^H -> L^H
    int INFO = lapack::potrf(lapack::Uplo::Upper,NBRI,S,NBRI);

    if (INFO)
      CErr("Error in Cholesky decomposition of auxiliary basis potential matrix. "
           "Error code : " + std::to_string(INFO));

    // Zero out LowerTriangular
    #pragma omp parallel for
    for (size_t i = 0; i < NBRI; i++)
      std::fill_n(S + i*NBRI + i + 1, NBRI - i - 1, 0.0);

    auto dur2CCholesky = tock(topERI3Trans);
    std::cout << "  RI-ERI3-Transformation-Cholesky duration = " << dur2CCholesky << " s " << std::endl;

    auto topTriInv = tick();

    // S = L^H^-1
    INFO = lapack::trtri(lapack::Uplo::Upper,lapack::Diag::NonUnit,NBRI,S,NBRI);

    if (INFO)
      CErr("Error in inverse of Cholesky decomposed auxiliary basis potential matrix.");

    auto durTriInv = tock(topTriInv);
    std::cout << "  RI-ERI3-Transformation-TriInv duration   = " << durTriInv << " s " << std::endl;
  }

  /**
   *  \brief Contract (P|Q)^{-1/2} with (Q|ij) to form L(Q|ij)
   *         Save L(Q|ij) in ERI3J
   */
  template <>
  void InCoreRITPI<dcomplex>::contract2CenterERI() {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreRITPI<double>::contract2CenterERI() {

    auto topGemm = tick();

    double *S = twocenterERI_->pointer();

    size_t NB2   = NB*(NB+1)/2;
    size_t NB3   = NB2*NBRI;
    // S^{-1/2}(Q|ij)
    auto ijK = CQMemManager::get().calloc<double>(NB3);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBRI,NB2,NBRI,double(1.),S,NBRI,pointer(),NBRI,double(0.),ijK,NBRI);

    auto durGemm = tock(topGemm);
    std::cout << "  RI-ERI3-Transformation-Gemm duration     = " << durGemm << " s " << std::endl;

#ifdef __DEBUGERI__
    // Debug output of the ERIs
    auto TempERI4 = CQMemManager::get().calloc<double>(NB2*NB2);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NB2,NB2,NBRI,double(1.),ijK,NB2,ijK,NB2,double(0.),TempERI4,NB2);
    std::cout << "Two-Electron Integrals (ERIs)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << TempERI4[i + j*NB  + k*NB2 + l*NB2*NB] << std::endl;
    };
    CQMemManager::get().free<double>(TempERI4);
#endif

    auto topCopy = tick();

    for (size_t pq = 0; pq < NB2; pq++) {
      auto pqAna = anaCompound(pq);
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.first, pqAna.second, NB));
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.second, pqAna.first, NB));
    }

    CQMemManager::get().free(ijK);

    auto durCopy = tock(topCopy);
    std::cout << "  RI-ERI3-Transformation-Copy duration     = " << durCopy << " s " << std::endl;

  } // InCoreRIERI<double>::contract2CenterERI


  /**
   *  \brief Compute three-center ERI (Q|ij) where Q is the auxiliary basis,
   *         and i,j are regular AO basis.
   *         Using Libint2 over the CGTO basis.
   */
  template <>
  void InCoreAuxBasisRIERI<dcomplex>::compute3CenterERI(BasisSet&, BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreAuxBasisRIERI<double>::compute3CenterERI(
      BasisSet &basisSet, BasisSet &auxBasisSet) {
    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      std::max(basisSet.maxPrim, auxBasisSet.maxPrim),
      std::max(basisSet.maxL, auxBasisSet.maxL),0);
    engines[0].set_precision(0.);
    engines[0].set(libint2::BraKet::xs_xx);
    const auto& unitshell = libint2::Shell::unit();

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

    // Allocate and zero out ERIs
    size_t NB    = basisSet.nBasis;
    this->setNRIBasis( auxBasisSet.nBasis );
    size_t NB2   = NB*NB;
    size_t NBRI2 = NBRI*NBRI;
    size_t NBNBRI= NB*NBRI;
    size_t NB3   = NB2*NBRI;



    auto topERI3 = tick();
    std::fill_n(ERI3J,NB3,0.);
    InCoreRITPI<double> &eri3j = *this;

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      size_t n1,n2,n3,i,j,k,ijk,bf1,bf2,bf3;
      size_t nShellDF = auxBasisSet.nShell;
      size_t nShell   = basisSet.nShell;
      // The outer loop runs over all DFBasis shells
      for(auto s1=0ul, bf1_s=0ul, s123=0ul; s1 < nShellDF; bf1_s+=n1, s1++) {

        n1 = auxBasisSet.shells[s1].size(); // Size of DFShell 1

      for(auto s2=0ul, bf2_s=0ul; s2 < nShell; bf2_s+=n2, s2++) {

        n2 = basisSet.shells[s2].size(); // Size of Shell 2

      for(auto s3=s2, bf3_s=bf2_s; s3 < nShell; bf3_s+=n3, s3++, s123++) {

        n3 = basisSet.shells[s3].size(); // Size of Shell 3

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s123 % nthreads != thread_id ) continue;
        #endif

        // Evaluate ERI3 for shell quartet
        engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(
          auxBasisSet.shells[s1],
          unitshell,
          basisSet.shells[s2],
          basisSet.shells[s3]
        );
        const auto *buff =  buf_vec[0] ;
        if(buff == nullptr) continue;

        // Place shell triplet into persistent storage
        for(i = 0ul, bf1 = bf1_s, ijk = 0ul  ; i < n1; ++i, bf1++)
          for(j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
            if (s2 == s3) {
              for(k = j, bf3 = bf2, ijk += j; k < n3; ++k, bf3++, ++ijk) {
                // (Q|12) -> RI-J
                eri3j.pointer()[bf1 + toCompound(bf2, bf3) * NBRI] = buff[ijk];
              }; // ijk loop
            } else {
              for(k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++, ++ijk) {
                // (Q|12) -> RI-J
                eri3j.pointer()[bf1 + toCompound(bf2, bf3) * NBRI] = buff[ijk];
              }; // ijk loop
            }

      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    auto durERI3 = tock(topERI3);
    std::cout << "  Libint-RI-ERI3 duration   = " << durERI3 << " s " << std::endl;

  }; // InCoreAuxBasisRIERI<double>::compute3CenterERI


  /**
   *  \brief Compute two-center ERI (P|Q) where P,Q are indices of the auxiliary basis.
   *         Using Libint2 over the CGTO basis.
   */
  template <>
  void InCoreAuxBasisRIERI<dcomplex>::compute2CenterERI(BasisSet&, dcomplex*) const {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreAuxBasisRIERI<double>::compute2CenterERI(
      BasisSet &auxBasisSet, double *S) const {

    std::fill_n(S,NBRI*NBRI,0.);

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      auxBasisSet.maxPrim, auxBasisSet.maxL,0);
    engines[0].set_precision(0.);
    engines[0].set(libint2::BraKet::xs_xs);
    const auto& unitshell = libint2::Shell::unit();

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    auto topERI2 = tick();
    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      size_t n1,n2,i,j,ij,bf1,bf2;
      size_t nShellDF = auxBasisSet.nShell;
      // The outer loop runs over all DFBasis shells
      for(auto s1=0ul, bf1_s=0ul, s12=0ul; s1 < nShellDF; bf1_s+=n1, s1++) {

        n1 = auxBasisSet.shells[s1].size(); // Size of DFShell 1

      for(auto s2=0ul, bf2_s=0ul; s2 < nShellDF; bf2_s+=n2, s2++, s12++) {

        n2 = auxBasisSet.shells[s2].size(); // Size of DFShell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        // Evaluate ERI2 for shell quartet
        engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xs_xs, 0>(
          auxBasisSet.shells[s1],
          unitshell,
          auxBasisSet.shells[s2],
          unitshell
        );
        const auto *buff =  buf_vec[0] ;
        if(buff == nullptr) continue;

        // Place shell doublet into persistent storage
        for(i = 0ul, bf1 = bf1_s, ij = 0ul ; i < n1; ++i, bf1++)
        for(j = 0ul, bf2 = bf2_s           ; j < n2; ++j, bf2++, ++ij) {
          // (1 | 2)
          S[bf1 + bf2*NBRI ] = buff[ij];
        }; // ij loop
      }; // s2
      }; // s1
    }; // omp region
    auto durERI2 = tock(topERI2);
    std::cout << "  Libint-RI-ERI2 duration   = " << durERI2 << " s " << std::endl;

  }; // InCoreAuxBasisRIERI<double>::compute2CenterERI


  /**
   *  \brief Compute and store the Auxiliary basis RI
   *  3-index ERI tensor using Libint2 over the CGTO basis.
   */ 
  template <>
  void InCoreAuxBasisRIERI<dcomplex>::computeAOInts(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreAuxBasisRIERI<double>::computeAOInts(BasisSet &basisSet, Molecule&,
      EMPerturbation&, OPERATOR op, const HamiltonianOptions &options) {

    if (op != ELECTRON_REPULSION)
      CErr("Only Electron repulsion integrals in InCoreAuxBasisRIERI<double>",std::cout);
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in InCoreAuxBasisRIERI<double>",std::cout);

    auto topLibintRI = tick();

    compute3CenterERI(basisSet, *auxBasisSet_);
    saveRawERI3J();

    twocenterERI_ = std::make_shared<cqmatrix::Matrix<double>>(NBRI);

    compute2CenterERI(*auxBasisSet_, twocenterERI_->pointer());
    if (saveRawERI_)
      rawERI2C_ = std::make_shared<cqmatrix::Matrix<double>>(*twocenterERI_);

    auto topERI3Trans = tick();
    halfInverse2CenterERI(*twocenterERI_);
    contract2CenterERI();

    auto durERI3Trans = tock(topERI3Trans);
    std::cout << "  RI-ERI3-Transformation duration = " << durERI3Trans << " s " << std::endl;

    auto durLibintRI = tock(topLibintRI);
    std::cout << "  Libint-RI duration   = " << durLibintRI << " s " << std::endl;

  }; // InCoreAuxBasisRIERI<double>::computeAOInts


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


  /**
   *  \brief Compute diagonal elements of ERI for Cholesky RI using libcint
   */
  void computeDiagonalPrebuilt4Index(double *diag,
      std::shared_ptr<InCore4indexTPI<double>> eri4I_){

    size_t NB = eri4I_->nBasis();

    auto topDiag = tick();

    size_t NB2 = NB*(NB+1)/2;
    #pragma omp parallel for
    for (size_t i = 0; i < NB2; i++) {
      size_t sI = compoundToSquare(i, NB);
      diag[i] = (*eri4I_)(sI,sI);
    }

    auto durDiag = tock(topDiag);
    std::cout << "  Cholesky-RI-Diagonal duration   = " << durDiag << " s " << std::endl;
  }


  /**
   *  \brief Compute diagonal elements of ERI for Cholesky RI using libcint
   */
  template <>
  std::map<std::pair<size_t, size_t>, dcomplex*>
  InCoreCholeskyRIERI<dcomplex>::computeDiagonalLibcint(
      BasisSet&, dcomplex*, bool) {

    std::map<std::pair<size_t, size_t>, dcomplex*> dummy;
    CErr("Only real GTOs are allowed",std::cout);
    return dummy;
  };
  template <>
  std::map<std::pair<size_t, size_t>, double*>
  InCoreCholeskyRIERI<double>::computeDiagonalLibcint(
      BasisSet &basisSet, double *diag, bool saveDiagBlocks) {
    
    std::map<std::pair<size_t, size_t>, double*> diagBlocks;

    if (eri4I_) {
      computeDiagonalPrebuilt4Index(diag, eri4I_);
      return diagBlocks;
    }

    auto topDiag = tick();

    if (saveDiagBlocks) {
      for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
        for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {

          size_t pqSize(basisSet.shells[P].size() * basisSet.shells[Q].size());
          diagBlocks[std::make_pair(P,Q)] = CQMemManager::get().calloc<double>(pqSize * pqSize);
        }
      }
    }

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel reduction(+:lc1ERI,lt1ERI)
    {
      size_t thread_id = GetThreadID();
      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
        for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {
          // Round Robbin work distribution
#ifdef _OPENMP
          if( PQ % nthreads != thread_id ) continue;
#endif

          if (saveDiagBlocks)
            buff = diagBlocks[std::make_pair(P,Q)];

          int shls[4]{int(P), int(Q), int(P), int(Q)};

          auto beginERI = tick();
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) {
            lt1ERI += tock(beginERI);
            lc1ERI++;
            continue;
          }
          lt1ERI += tock(beginERI);
          lc1ERI++;

          if (P == Q) {
            for (size_t p(0), pBegin(basisSet.mapSh2Bf[P]),
                 pSize(basisSet.shells[P].size());
                 p < pSize; p++) {
              for (size_t q(p); q < pSize; q++) {
                diag[toCompound(pBegin + p, pBegin + q)] =
                    buff[(p + q * pSize) * (1 + pSize * pSize)];
              }
            }

          } else {
            for (size_t p(0), pBegin(basisSet.mapSh2Bf[P]),
                 pSize(basisSet.shells[P].size());
                 p < pSize; p++) {
              for (size_t q(0), qBegin(basisSet.mapSh2Bf[Q]),
                   qSize(basisSet.shells[Q].size());
                   q < qSize; q++) {
                diag[toCompound(pBegin + p, qBegin + q)] =
                    buff[(p + q * pSize) * (1 + pSize * qSize)];
              }
            }
          }

#ifdef __DEBUGERI__
          size_t pqSize = basisSet.shells[P].size() * basisSet.shells[Q].size();
          prettyPrintSmart(std::cout, "ERI("+std::to_string(P)+","+std::to_string(Q)+")",
                           diagBlocks[std::make_pair(P,Q)], pqSize, pqSize, pqSize);
#endif
        }; // Q
      }; // P
    }; // omp region
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    auto durDiag = tock(topDiag);
    std::cout << "  Cholesky-RI-Diagonal duration   = " << durDiag << " s " << std::endl;

    return diagBlocks;

  }; // InCoreCholeskyRIERI<double>::computeDiagonalLibcint



  /**
   *  \brief Compute diagonal elements of ERI for Cholesky RI using libint2
   */
  template <>
  std::map<std::pair<size_t, size_t>, dcomplex*>
  InCoreCholeskyRIERI<dcomplex>::computeDiagonalLibint(
      BasisSet&, dcomplex*, bool) {
    std::map<std::pair<size_t, size_t>, dcomplex*> dummy;
    CErr("Only real GTOs are allowed",std::cout);
    return dummy;
  };
  template <>
  std::map<std::pair<size_t, size_t>, double*>
  InCoreCholeskyRIERI<double>::computeDiagonalLibint(
      BasisSet &basisSet, double *diag, bool saveDiagBlocks) {
    
    std::map<std::pair<size_t, size_t>, double*> diagBlocks;

    if (eri4I_) {
      computeDiagonalPrebuilt4Index(diag, eri4I_);
      return diagBlocks;
    }

    auto topDiag = tick();

    for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
      for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {

        if (saveDiagBlocks or basisSet.shells[P].ncontr() > 1 or basisSet.shells[Q].ncontr() > 1) {
          size_t pqSize(basisSet.shells[P].size() * basisSet.shells[Q].size());
          diagBlocks[std::make_pair(P,Q)] = CQMemManager::get().calloc<double>(pqSize * pqSize);
        }
      }
    }

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel reduction(+:lc1ERI,lt1ERI)
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
        for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {
          // Round Robbin work distribution
#ifdef _OPENMP
          if( PQ % nthreads != thread_id ) continue;
#endif

          if (basisSet.shells[P].ncontr() == 1 and basisSet.shells[Q].ncontr() == 1) {

            // Evaluate ERI for shell quartet
            auto beginERI = tick();
            engines[thread_id].compute2<
              libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                  basisSet.shells[P],
                  basisSet.shells[Q],
                  basisSet.shells[P],
                  basisSet.shells[Q]
            );
            lt1ERI += tock(beginERI);
            lc1ERI++;
            const auto *buff = buf_vec[0];
            if(buff == nullptr) continue;

            double *block = nullptr;
            if (saveDiagBlocks)
              block = diagBlocks[std::make_pair(P,Q)];

            if (P == Q) {
              for (size_t p(0), pBegin(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size());
                   p < pSize; p++) {
                for (size_t q(p); q < pSize; q++) {
                  diag[toCompound(pBegin + p, pBegin + q)] =
                      buff[(p * pSize + q) * (1 + pSize * pSize)];

                  if (saveDiagBlocks) {
                    size_t pqrs = (p * pSize + q) * pSize * pSize;
                    double *pqBlock = block + (p + q * pSize) * pSize * pSize;
                    for (size_t r(0); r < pSize; r++) {
                      for (size_t s(0); s < pSize; s++, pqrs++) {
                        pqBlock[r + s * pSize] = buff[pqrs];
                      }
                    }

                    if (p != q) {
                      pqrs = (q * pSize + p) * pSize * pSize;
                      pqBlock = block + (q + p * pSize) * pSize * pSize;
                      for (size_t r(0); r < pSize; r++) {
                        for (size_t s(0); s < pSize; s++, pqrs++) {
                          pqBlock[r + s * pSize] = buff[pqrs];
                        }
                      }
                    }

                  }

                }
              }

            } else {
              for (size_t p(0), pqrs(0), pBegin(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size());
                   p < pSize; p++) {
                for (size_t q(0), qBegin(basisSet.mapSh2Bf[Q]),
                     qSize(basisSet.shells[Q].size());
                     q < qSize; q++) {
                  diag[toCompound(pBegin + p, qBegin + q)] =
                      buff[(p * qSize + q) * (1 + pSize * qSize)];

                  if (saveDiagBlocks) {
                    double *pqBlock = block + (p + q * pSize) * pSize * qSize;
                    for (size_t r(0); r < pSize; r++) {
                      for (size_t s(0); s < qSize; s++, pqrs++) {
                        pqBlock[r + s * pSize] = buff[pqrs];
                      }
                    }

                  }

                }
              }
            }

          } else {

            double *block = diagBlocks[std::make_pair(P,Q)];

            size_t pContrSize = basisSet.shells[P].contr.size();
            size_t qContrSize = basisSet.shells[Q].contr.size();

            size_t pAMSize = shellPrims_[P][0].size();
            size_t qAMSize = shellPrims_[Q][0].size();

            size_t pqrsAMSize = pAMSize * qAMSize * pAMSize * qAMSize;

            const double *resPQRS;
            std::pair<size_t, double> counter_timer = libintGeneralContractionERI(
                P, Q, P, Q,
                basisSet, engines[thread_id],
                shellPrims_, coefBlocks_, workBlocks[thread_id],
                resPQRS);
            lt1ERI += counter_timer.second;
            lc1ERI += counter_timer.first;

            // Reorganize
            for (size_t SS = 0, pqrs = 0; SS < qContrSize; SS++)
              for (size_t s = 0; s < qAMSize; s++)
                for (size_t RR = 0; RR < pContrSize; RR++)
                  for (size_t r = 0; r < pAMSize; r++)
                    for (size_t QQ = 0; QQ < qContrSize; QQ++)
                      for (size_t q = 0; q < qAMSize; q++)
                        for (size_t PP = 0; PP < pContrSize; PP++)
                          for (size_t p = 0; p < pAMSize; p++, pqrs++) {
                            block[pqrs] = resPQRS[(s + qAMSize * (r + pAMSize * (q + qAMSize * p)))
                                + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + pContrSize * SS)))];
                          }

            if (P == Q) {
              for (size_t p(0), pBegin(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size());
                   p < pSize; p++) {
                for (size_t q(p); q < pSize; q++) {
                  diag[toCompound(pBegin + p, pBegin + q)] =
                      block[(p + q * pSize) * (1 + pSize * pSize)];
                }
              }

            } else {
              for (size_t p(0), pBegin(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size());
                   p < pSize; p++) {
                for (size_t q(0), qBegin(basisSet.mapSh2Bf[Q]),
                     qSize(basisSet.shells[Q].size());
                     q < qSize; q++) {
                  diag[toCompound(pBegin + p, qBegin + q)] =
                      block[(p + q * pSize) * (1 + pSize * qSize)];
                }
              }
            }
          }


        }; // Q
      }; // P
    }; // omp region
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    if (not saveDiagBlocks) {
      for (auto &kv : diagBlocks) {
        CQMemManager::get().free(kv.second);
      }
      diagBlocks.clear();
    }

    auto durDiag = tock(topDiag);
    std::cout << "  Cholesky-RI-Diagonal duration   = " << durDiag << " s " << std::endl;

    return diagBlocks;

  }; // InCoreCholeskyRIERI<double>::computeDiagonalLibint



  /**
   *  \brief Compute 3-index ERI vector Cholesky RI in a shell
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeTraditionalERIVectorByShellLibcint(
      BasisSet&, size_t, size_t, size_t, size_t, dcomplex*) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeTraditionalERIVectorByShellLibcint(
      BasisSet &basisSet, size_t r, size_t s,
      size_t R, size_t S, double *L) {

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel reduction(+:lc1ERI,lt1ERI)
    {
      size_t thread_id = GetThreadID();

      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll + cache_size*thread_id;
      int shls[4]{int(0), int(0), int(S), int(R)};

      for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
        for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {
          // Round Robbin work distribution
#ifdef _OPENMP
          if( PQ % nthreads != thread_id ) continue;
#endif

          shls[0] = int(Q);
          shls[1] = int(P);

          auto beginERI = tick();
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) {
            lt1ERI += tock(beginERI);
            lc1ERI++;
            continue;
          }
          lt1ERI += tock(beginERI);
          lc1ERI++;

          if (P == Q) {
            for (size_t pBegin(basisSet.mapSh2Bf[P]),
                 pSize(basisSet.shells[P].size()),
                 pEnd(pBegin + pSize),
                 p(pBegin),
                 rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                 p < pEnd; p++, rsp++) {
              for (size_t q(p),
                   rspq = rsp * pSize + p - pBegin;
                   q < pEnd; q++, rspq++) {
                L[toCompound(p,q)] = buff[rspq];
              }
            }
          } else {
            for (size_t p(basisSet.mapSh2Bf[P]),
                 pSize(basisSet.shells[P].size()),
                 pEnd(p + pSize),
                 rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                 p < pEnd; p++, rsp++) {
              for (size_t q(basisSet.mapSh2Bf[Q]),
                   qSize(basisSet.shells[Q].size()),
                   qEnd(q + qSize),
                   rspq = rsp * qSize;
                   q < qEnd; q++, rspq++) {
                L[toCompound(p,q)] = buff[rspq];
              }
            }
          }

        }; // Q
      }; // P
    }; // omp region
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

  }; // InCoreCholeskyRIERI<double>::computeTraditionalERIVectorByShellLibcint



  /**
   *  \brief Compute 3-index ERI vector Cholesky RI in a shell
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeTraditionalERIVectorByShellLibint(
      BasisSet&, size_t, size_t, size_t, size_t, dcomplex*) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeTraditionalERIVectorByShellLibint(
      BasisSet &basisSet, size_t r, size_t s,
      size_t R, size_t S, double *L) {

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t rSize = basisSet.shells[R].size();
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();
    size_t NB2 = basisSet.nBasis * (basisSet.nBasis + 1) / 2;

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel reduction(+:lc1ERI,lt1ERI)
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
        for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++) {
          // Round Robbin work distribution
#ifdef _OPENMP
          if( PQ % nthreads != thread_id ) continue;
#endif

          if (basisSet.shells[P].ncontr() == 1 and basisSet.shells[Q].ncontr() == 1
              and basisSet.shells[R].ncontr() == 1 and basisSet.shells[S].ncontr() == 1) {
            // Evaluate ERI for shell quartet
            auto beginERI = tick();
            engines[thread_id].compute2<
              libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                  basisSet.shells[R],
                  basisSet.shells[S],
                  basisSet.shells[P],
                  basisSet.shells[Q]
            );
            lt1ERI += tock(beginERI);
            lc1ERI++;
            const auto *buff =  buf_vec[0] ;
            if(buff == nullptr) continue;

            if (P == Q) {
              for (size_t pBegin(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size()),
                   pEnd(pBegin + pSize),
                   p(pBegin),
                   rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                   p < pEnd; p++, rsp++) {
                for (size_t q(p),
                     rspq = rsp * pSize + p - pBegin;
                     q < pEnd; q++, rspq++) {
                  L[toCompound(p,q)] = buff[rspq];
                }
              }
            } else {
              for (size_t p(basisSet.mapSh2Bf[P]),
                   pSize(basisSet.shells[P].size()),
                   pEnd(p + pSize),
                   rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                   p < pEnd; p++, rsp++) {
                for (size_t q(basisSet.mapSh2Bf[Q]),
                     qSize(basisSet.shells[Q].size()),
                     qEnd(q + qSize),
                     rspq = rsp * qSize;
                     q < qEnd; q++, rspq++) {
                  L[toCompound(p,q)] = buff[rspq];
                }
              }
            }

          } else {

            size_t pContrSize = basisSet.shells[P].contr.size();
            size_t qContrSize = basisSet.shells[Q].contr.size();
            size_t rContrSize = basisSet.shells[R].contr.size();

            size_t pAMSize = shellPrims_[P][0].size();
            size_t qAMSize = shellPrims_[Q][0].size();
            size_t rAMSize = shellPrims_[R][0].size();
            size_t sAMSize = shellPrims_[S][0].size();

            size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

            const double *resPQRS;
            std::pair<size_t, double> counter_timer = libintGeneralContractionERI(
                P, Q, R, S,
                basisSet, engines[thread_id],
                shellPrims_, coefBlocks_, workBlocks[thread_id],
                resPQRS);
            lt1ERI += counter_timer.second;
            lc1ERI += counter_timer.first;


            size_t rShift = r - rBegin;
            size_t RR = rShift / rAMSize, rr = rShift % rAMSize;

            size_t sShift = s - sBegin;
            size_t SS = sShift / sAMSize, ss = sShift % sAMSize;

            if (P == Q) {
              for (size_t q(0),
                   qBegin(basisSet.mapSh2Bf[Q]),
                   qSize(basisSet.shells[Q].size());
                   q < qSize; q++) {

                size_t QQ = q / qAMSize, qq = q % qAMSize;

                for (size_t p(0); p <= q; p++) {

                  size_t PP = p / pAMSize, pp = p % pAMSize;

                  L[toCompound(p + qBegin, q + qBegin)] =
                      resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                              + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                }
              }

            } else {
              for (size_t q(0),
                   qBegin(basisSet.mapSh2Bf[Q]),
                   qSize(basisSet.shells[Q].size());
                   q < qSize; q++) {

                size_t QQ = q / qAMSize, qq = q % qAMSize;

                for (size_t p(0),
                     pBegin(basisSet.mapSh2Bf[P]),
                     pSize(basisSet.shells[P].size());
                     p < pSize; p++) {

                  size_t PP = p / pAMSize, pp = p % pAMSize;

                  L[toCompound(p + pBegin, q + qBegin)] =
                      resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                              + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                }
              }
            }

          }

        }; // Q
      }; // P
    }; // omp region
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

  }; // InCoreCholeskyRIERI<double>::computeTraditionalERIVectorByShellLibint



  /**
   *  \brief Compute 3-index ERI vector Cholesky RI
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeTraditionalERIVector(BasisSet&, size_t, dcomplex*) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeTraditionalERIVector(
      BasisSet &basisSet, size_t pivot, double *L) {

    auto pivotPair = anaCompound(pivot);
    size_t r = pivotPair.first;
    size_t s = pivotPair.second;

    if (eri4I_) {
      size_t pivotSq = toSquare(r, s, NB);
      size_t NB2 = NB*(NB+1)/2;
      #pragma omp parallel for
      for (size_t i = 0; i < NB2; i++) {
        L[i] = (*eri4I_)(compoundToSquare(i, NB), pivotSq);
      }
      return;
    }

    size_t R = std::distance(basisSet.mapSh2Bf.begin(),
        std::upper_bound(basisSet.mapSh2Bf.begin(),
            basisSet.mapSh2Bf.end(), r)) - 1;
    size_t S = std::distance(basisSet.mapSh2Bf.begin(),
        std::upper_bound(basisSet.mapSh2Bf.begin(),
            basisSet.mapSh2Bf.end(), s)) - 1;

    if (libcint_)
      computeTraditionalERIVectorByShellLibcint(basisSet, r, s, R, S, L);
    else
      computeTraditionalERIVectorByShellLibint(basisSet, r, s, R, S, L);

  }; // InCoreCholeskyRIERI<double>::computeTraditionalERIVector


  /**
   *  \brief Allocate, compute and store the Cholesky RI
   *         3-index ERI tensor using Libint2 over the CGTO basis.
   *         Using standard index-by-index algorithm.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeCD_Traditional(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeCD_Traditional(BasisSet &basisSet) {

    std::cout << "Traditional Cholesky Decomposition:" << std::endl;
    std::cout << bannerMid << std::endl;

    auto topCholesky = tick();

    size_t NB2 = NB*(NB+1)/2;

    double *diag = CQMemManager::get().calloc<double>(NB2);

    if (libcint_)
      computeDiagonalLibcint(basisSet, diag);
    else
      computeDiagonalLibint(basisSet, diag);

    // Temporary cholesky factor
    std::vector<double*> L, allocs;

    pivots_.clear();
    size_t NBRI = 0;

    auto beginERIvec = topCholesky;
    auto beginCDalg = topCholesky;
    double cumeERIvec = 0, cumeCDalg = 0;

    while (NBRI < NB2) {
      // Select the pivot
      size_t pivot = 0;
      double Dmax = diag[0];
      for (size_t P = 0; P < NB2; P++) {
          if (Dmax < diag[P]) {
              Dmax = diag[P];
              pivot = P;
          }
      }

      // Check to see if convergence reached
      if (Dmax < tau_) break;

#ifdef CD_PROGRESS
      auto pivotPQ = anaCompound(pivot);
      std::cout << "    Selected: (" << pivotPQ.first << "," << pivotPQ.second
                << ")\tDiag: " << Dmax << std::endl;
#endif

      // If here, we're trying to add this row
      pivots_.push_back(pivot);
      double L_QQ = sqrt(Dmax);

      // If here, we're really going to add this row
      if (NBRI % NB == 0)
        allocs.push_back(CQMemManager::get().calloc<double>(NB * NB2));
      L.push_back(allocs.back() + (NBRI % NB) * NB2);

      beginERIvec = tick();
      computeTraditionalERIVector(basisSet, pivot, L[NBRI]);
      cumeERIvec += tock(beginERIvec);

      beginCDalg = tick();
      // [(m|Q) - L_m^P L_Q^P]
      for (size_t nBatch(NBRI / NB + 1), i(0); i < nBatch; i++)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, NB2, 1,
             NBRI - i * NB > NB ? NB : NBRI - i * NB,
             -1.0, allocs[i], NB2,
             allocs[i] + pivots_[NBRI], NB2,
             1.0, L[NBRI], NB2);

      // 1/L_QQ [(m|Q) - L_m^P L_Q^P]
      blas::scal(NB2, 1.0 / L_QQ, L[NBRI], 1);

      // Zero the upper triangle
      #pragma omp parallel for
      for (size_t P = 0; P < pivots_.size(); P++) {
          L[NBRI][pivots_[P]] = 0.0;
      }

      cumeCDalg += tock(beginCDalg);

      // Set the pivot factor
      L[NBRI][pivot] = L_QQ;

      // Update the Schur complement diagonal
      #pragma omp parallel for
      for (size_t P = 0; P < NB2; P++) {
          diag[P] -= L[NBRI][P] * L[NBRI][P];
      }

      // Force truly zero elements to zero
      #pragma omp parallel for
      for (size_t P = 0; P < pivots_.size(); P++) {
          diag[pivots_[P]] = 0.0;
      }

      NBRI++;
    }

    if (NBRI == 0)
      CErr("Cholesky decomposition threshold is greater than all TPI diagonal elements.");

    auto durCholesky = tock(topCholesky);
    std::cout << "  Cholesky-RI-Pivots-ERI count    = " << c1ERI << std::endl;
    std::cout << "  Cholesky-RI-Pivots-ERI duration = " << t1ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIvec duration     = " << cumeERIvec << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalg duration      = " << cumeCDalg << " s " << std::endl;
    std::cout << "  Cholesky-RI-misc duration       = "
              << durCholesky - cumeERIvec - cumeCDalg << " s " << std::endl;
    std::cout << "  Cholesky-RI-Traditional duration   = " << durCholesky << " s " << std::endl;

    this->setNRIBasis(NBRI);
    std::cout << "  Cholesky-RI auxiliary dimension = " << NBRI << std::endl;

    #pragma omp parallel for
    for (size_t Q = 0; Q < NBRI; Q++) {
      for (size_t ij = 0; ij < NB2; ij++) {
        auto ijP = anaCompound(ij);
        const double &v = L[Q][ij];
        (*this)(Q,ijP.first,ijP.second) = v;
        (*this)(Q,ijP.second,ijP.first) = v;
      }
    }

    CQMemManager::get().free(diag);
    for (double *p : allocs) {
      CQMemManager::get().free(p);
    }

  }; // InCoreCholeskyRIERI<double>::computeCD_Traditional

  /**
   *  \brief Compute the Mpq matrix in the span factor Cholesky pivots
   *         algorithm using Libcint over the CGTO basis.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeSpanFactorERIVectorLibcint(BasisSet&,
      const std::vector<size_t> &D, std::vector<SpanFactorShellPair>&,
      size_t lenQ, size_t qInd, dcomplex* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, dcomplex*> &diagBlocks,
      bool Mpq_only) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeSpanFactorERIVectorLibcint(BasisSet &basisSet,
      const std::vector<size_t> &D,
      std::vector<SpanFactorShellPair> &shellPair_Dindices,
      size_t lenQ, size_t qShellPairIndex, double* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      bool Mpq_only) {

    size_t lenD = D.size();
    size_t lenShellD = shellPair_Dindices.size();

    size_t R = shellPair_Dindices[qShellPairIndex].PQ.first;
    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t rSize = basisSet.shells[R].size();
    size_t S = shellPair_Dindices[qShellPairIndex].PQ.second;
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel for schedule(guided) reduction(+:lc1ERI,lt1ERI)
    for (size_t pShellPairIndex = 0; pShellPairIndex < lenShellD; pShellPairIndex++) {

      if (shellPair_Dindices[pShellPairIndex].evaluated
          or qShellPairIndex == pShellPairIndex)
        continue;

      size_t P = shellPair_Dindices[pShellPairIndex].PQ.first;
      size_t pBegin = basisSet.mapSh2Bf[P];
      size_t pSize = basisSet.shells[P].size();
      size_t Q = shellPair_Dindices[pShellPairIndex].PQ.second;
      size_t qBegin = basisSet.mapSh2Bf[Q];
      size_t qSize = basisSet.shells[Q].size();

      size_t thread_id = GetThreadID();

      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      int shls[4]{int(Q), int(P), int(S), int(R)};

      // Evaluate ERI for shell quartet
      auto beginERI = tick();
      if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) {
        lt1ERI += tock(beginERI);
        lc1ERI++;
        continue;
      }
      lt1ERI += tock(beginERI);
      lc1ERI++;

      if (Mpq_only) {
        for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

          size_t r = D[rs] % NB - rBegin;
          size_t s = D[rs] / NB - sBegin;

          for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

            size_t p = D[pq] % NB - pBegin;
            size_t q = D[pq] / NB - qBegin;

            size_t rspq = ((r * sSize + s) * pSize + p) * qSize + q;

            if (rs < lenQ) {
              ERIvecAlloc[pq + rs * lenD] = buff[rspq];
            }

            if (pq < lenQ) {
              ERIvecAlloc[rs + pq * lenD] = buff[rspq];
            }

          }; // pq
        }; // rs

      } else {
        bool needTranspose = shellPair_Dindices[pShellPairIndex].indices[0] < lenQ;

        for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

          size_t r = D[rs] % NB - rBegin;
          size_t s = D[rs] / NB - sBegin;

          for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

            size_t p = D[pq] % NB - pBegin;
            size_t q = D[pq] / NB - qBegin;

            size_t rspq = ((r * sSize + s) * pSize + p) * qSize + q;

            ERIvecAlloc[pq + rs * lenD] = buff[rspq];
            if (needTranspose)
              ERIvecAlloc[rs + pq * lenD] = buff[rspq];

          }; // pq
        }; // rs

      }; // Mpq_only

    }; // PQ
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    shellPair_Dindices[qShellPairIndex].evaluated = true;

  }; // InCoreCholeskyRIERI<double>::computeSpanFactorERIVectorLibcint


  /**
   *  \brief Compute the Mpq matrix in the span factor Cholesky pivots
   *         algorithm using Libint2 over the CGTO basis.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeSpanFactorERIVectorLibint(BasisSet&,
      const std::vector<size_t> &D, std::vector<SpanFactorShellPair>&,
      size_t lenQ, size_t qInd, dcomplex* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, dcomplex*> &diagBlocks,
      bool Mpq_only) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeSpanFactorERIVectorLibint(BasisSet &basisSet,
      const std::vector<size_t> &D,
      std::vector<SpanFactorShellPair> &shellPair_Dindices,
      size_t lenQ, size_t qShellPairIndex, double* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      bool Mpq_only) {

    size_t lenD = D.size();
    size_t lenShellD = shellPair_Dindices.size();

    size_t R = shellPair_Dindices[qShellPairIndex].PQ.first;
    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t rSize = basisSet.shells[R].size();
    size_t S = shellPair_Dindices[qShellPairIndex].PQ.second;
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel for schedule(guided) reduction(+:lc1ERI,lt1ERI)
    for (size_t pShellPairIndex = 0; pShellPairIndex < lenShellD; pShellPairIndex++) {

      if (shellPair_Dindices[pShellPairIndex].evaluated
          or qShellPairIndex == pShellPairIndex)
        continue;

      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      size_t P = shellPair_Dindices[pShellPairIndex].PQ.first;
      size_t pBegin = basisSet.mapSh2Bf[P];
      size_t pSize = basisSet.shells[P].size();
      size_t Q = shellPair_Dindices[pShellPairIndex].PQ.second;
      size_t qBegin = basisSet.mapSh2Bf[Q];
      size_t qSize = basisSet.shells[Q].size();

      if (basisSet.shells[P].ncontr() == 1 and basisSet.shells[Q].ncontr() == 1
          and basisSet.shells[R].ncontr() == 1 and basisSet.shells[S].ncontr() == 1) {

        // Evaluate ERI for shell quartet
        auto beginERI = tick();
        engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              basisSet.shells[R],
              basisSet.shells[S],
              basisSet.shells[P],
              basisSet.shells[Q]
        );
        lt1ERI += tock(beginERI);
        lc1ERI++;
        const auto *buff =  buf_vec[0] ;
        if(buff == nullptr) continue;

        if (Mpq_only) {
          for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

            size_t r = D[rs] % NB - rBegin;
            size_t s = D[rs] / NB - sBegin;

            for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

              size_t p = D[pq] % NB - pBegin;
              size_t q = D[pq] / NB - qBegin;

              size_t rspq = ((r * sSize + s) * pSize + p) * qSize + q;

              if (rs < lenQ) {
                ERIvecAlloc[pq + rs * lenD] = buff[rspq];
              }

              if (pq < lenQ) {
                ERIvecAlloc[rs + pq * lenD] = buff[rspq];
              }

            }; // pq
          }; // rs

        } else {
          bool needTranspose = shellPair_Dindices[pShellPairIndex].indices[0] < lenQ;

          for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

            size_t r = D[rs] % NB - rBegin;
            size_t s = D[rs] / NB - sBegin;

            for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

              size_t p = D[pq] % NB - pBegin;
              size_t q = D[pq] / NB - qBegin;

              size_t rspq = ((r * sSize + s) * pSize + p) * qSize + q;

              ERIvecAlloc[pq + rs * lenD] = buff[rspq];
              if (needTranspose)
                ERIvecAlloc[rs + pq * lenD] = buff[rspq];

            }; // pq
          }; // rs

        }; // Mpq_only

      } else {

        size_t pContrSize = basisSet.shells[P].contr.size();
        size_t qContrSize = basisSet.shells[Q].contr.size();
        size_t rContrSize = basisSet.shells[R].contr.size();

        size_t pAMSize = shellPrims_[P][0].size();
        size_t qAMSize = shellPrims_[Q][0].size();
        size_t rAMSize = shellPrims_[R][0].size();
        size_t sAMSize = shellPrims_[S][0].size();

        size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

        const double *resPQRS;
        std::pair<size_t, double> counter_timer = libintGeneralContractionERI(
            P, Q, R, S,
            basisSet, engines[thread_id],
            shellPrims_, coefBlocks_, workBlocks[thread_id],
            resPQRS);
        lt1ERI += counter_timer.second;
        lc1ERI += counter_timer.first;

        if (Mpq_only) {
          for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

            size_t r = D[rs] % NB - rBegin;
            size_t RR = r / rAMSize, rr = r % rAMSize;
            size_t s = D[rs] / NB - sBegin;
            size_t SS = s / sAMSize, ss = s % sAMSize;

            for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

              size_t p = D[pq] % NB - pBegin;
              size_t PP = p / pAMSize, pp = p % pAMSize;
              size_t q = D[pq] / NB - qBegin;
              size_t QQ = q / qAMSize, qq = q % qAMSize;

              size_t rspq = (ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)));

              if (rs < lenQ) {
                ERIvecAlloc[pq + rs * lenD] = resPQRS[rspq];
              }

              if (pq < lenQ) {
                ERIvecAlloc[rs + pq * lenD] = resPQRS[rspq];
              }

            }; // pq
          }; // rs

        } else {
          bool needTranspose = shellPair_Dindices[pShellPairIndex].indices[0] < lenQ;

          for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

            size_t r = D[rs] % NB - rBegin;
            size_t RR = r / rAMSize, rr = r % rAMSize;
            size_t s = D[rs] / NB - sBegin;
            size_t SS = s / sAMSize, ss = s % sAMSize;

            for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

              size_t p = D[pq] % NB - pBegin;
              size_t PP = p / pAMSize, pp = p % pAMSize;
              size_t q = D[pq] / NB - qBegin;
              size_t QQ = q / qAMSize, qq = q % qAMSize;

              size_t rspq = (ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)));

              ERIvecAlloc[pq + rs * lenD] = resPQRS[rspq];
              if (needTranspose)
                ERIvecAlloc[rs + pq * lenD] = resPQRS[rspq];

            }; // pq
          }; // rs

        }; // Mpq_only

      }; // if general contraction

    }; // PQ
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    shellPair_Dindices[qShellPairIndex].evaluated = true;

  }; // InCoreCholeskyRIERI<double>::computeSpanFactorERIVectorLibint


  /**
   *  \brief Compute the Mpq matrix in the span factor Cholesky pivots
   *         algorithm using Libint2 over the CGTO basis.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeSpanFactorERI(BasisSet&,
      const std::vector<size_t> &D, std::vector<SpanFactorShellPair>&,
      size_t lenQ, dcomplex* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, dcomplex*> &diagBlocks,
      bool Mpq_only) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeSpanFactorERI(BasisSet &basisSet,
      const std::vector<size_t> &D,
      std::vector<SpanFactorShellPair> &shellPair_Dindices,
      size_t lenQ, double* ERIvecAlloc,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      bool Mpq_only) {

    size_t lenShellD = shellPair_Dindices.size();
    size_t lenD = D.size();

    // Copy diagonal blocks
    if (not eri4I_) {
      // Determine the number of OpenMP threads
      size_t nthreads = GetNumThreads();

      #pragma omp parallel
      {
        size_t thread_id = GetThreadID();

        for (size_t qShellPairIndex = 0; qShellPairIndex < lenShellD; qShellPairIndex++) {

          // Round Robbin work distribution
          #ifdef _OPENMP
          if( qShellPairIndex % nthreads != thread_id ) continue;
          #endif

          if (shellPair_Dindices[qShellPairIndex].indices[0] >= lenQ)
            break;

          if (shellPair_Dindices[qShellPairIndex].evaluated)
            continue;

          size_t R = shellPair_Dindices[qShellPairIndex].PQ.first;
          size_t rBegin = basisSet.mapSh2Bf[R];
          size_t rSize = basisSet.shells[R].size();
          size_t S = shellPair_Dindices[qShellPairIndex].PQ.second;
          size_t sBegin = basisSet.mapSh2Bf[S];
          size_t sSize = basisSet.shells[S].size();

          double *block = diagBlocks[std::make_pair(R,S)];

          if (Mpq_only) {
            for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {
              if (rs >= lenQ)
                break;

              size_t r = D[rs] % NB - rBegin;
              size_t s = D[rs] / NB - sBegin;

              for (size_t pq : shellPair_Dindices[qShellPairIndex].indices) {

                size_t p = D[pq] % NB - rBegin;
                size_t q = D[pq] / NB - sBegin;

                size_t rspq = r + rSize * (s + sSize * (p + rSize * q));

                ERIvecAlloc[pq + rs * lenD] = block[rspq];

              }; // pq
            }; // rs

          } else {
            for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

              size_t r = D[rs] % NB - rBegin;
              size_t s = D[rs] / NB - sBegin;

              for (size_t pq : shellPair_Dindices[qShellPairIndex].indices) {

                size_t p = D[pq] % NB - rBegin;
                size_t q = D[pq] / NB - sBegin;

                size_t rspq = r + rSize * (s + sSize * (p + rSize * q));

                ERIvecAlloc[pq + rs * lenD] = block[rspq];

              }; // pq
            }; // rs

          }; // Mpq_only
        }
      }; // OMP region
    }

    for (size_t qShellPairIndex = 0; qShellPairIndex < lenShellD; qShellPairIndex++) {

      if (shellPair_Dindices[qShellPairIndex].indices[0] >= lenQ)
        break;

      if (not shellPair_Dindices[qShellPairIndex].evaluated) {

        if (eri4I_) {

          size_t lenD = D.size();

          if (Mpq_only) {
            #pragma omp parallel for schedule(guided)
            for (size_t pShellPairIndex = 0; pShellPairIndex < lenShellD; pShellPairIndex++) {

              if (shellPair_Dindices[pShellPairIndex].evaluated)
                continue;

              for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

                size_t Drs = D[rs];

                for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

                  size_t Dpq = D[pq];

                  if (rs < lenQ) {
                    ERIvecAlloc[pq + rs * lenD] = (*eri4I_)(Drs,Dpq);
                  }

                  if (pq < lenQ) {
                    ERIvecAlloc[rs + pq * lenD] = (*eri4I_)(Drs,Dpq);
                  }

                }; // pq
              }; // rs
            }; // PQ

          } else {
            #pragma omp parallel for schedule(guided)
            for (size_t pShellPairIndex = 0; pShellPairIndex < lenShellD; pShellPairIndex++) {

              if (shellPair_Dindices[pShellPairIndex].evaluated)
                continue;

              bool needTranspose = shellPair_Dindices[pShellPairIndex].indices[0] < lenQ;

              for (size_t rs : shellPair_Dindices[qShellPairIndex].indices) {

                size_t Drs = D[rs];

                for (size_t pq : shellPair_Dindices[pShellPairIndex].indices) {

                  size_t Dpq = D[pq];

                  ERIvecAlloc[pq + rs * lenD] = (*eri4I_)(Drs,Dpq);
                  if (needTranspose)
                    ERIvecAlloc[rs + pq * lenD] = (*eri4I_)(Drs,Dpq);

                }; // pq
              }; // rs
            }; // PQ

          }; // Mpq_only

          shellPair_Dindices[qShellPairIndex].evaluated = true;

        } else if (libcint_)
          computeSpanFactorERIVectorLibcint(basisSet, D, shellPair_Dindices, lenQ, qShellPairIndex, ERIvecAlloc, diagBlocks, Mpq_only);
        else
          computeSpanFactorERIVectorLibint(basisSet, D, shellPair_Dindices, lenQ, qShellPairIndex, ERIvecAlloc, diagBlocks, Mpq_only);

      }

    }; // RS

  }; // InCoreCholeskyRIERI<double>::computeSpanFactorERI


  /**
   *  \brief Compute the Cholesky RI pivots using the span factor algorithm
   *         and Libint2 over the CGTO basis.
   *         Algorithm refers to 10.1063/1.5083802
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeCDPivots_SpanFactor(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeCDPivots_SpanFactor(BasisSet &basisSet) {

    std::cout << "Span-Factor Pivots Determination:" << std::endl;
    std::cout << bannerMid << std::endl;

    auto topEffCDPivots = tick();

    size_t NB2 = NB*NB;


    // 2. Clear pivots_ and CD vectors
    pivots_.clear();
    double* L;


    // 1. Compute diagonal
    // 3. Select all diagonals greater than theoreshold
    double *diag = CQMemManager::get().calloc<double>(NB2);
    double *diagCompound = CQMemManager::get().calloc<double>(NB*(NB+1)/2);
    std::map<std::pair<size_t, size_t>, double*> diagBlocks;
    std::vector<size_t> D;

    if (libcint_)
      diagBlocks = computeDiagonalLibcint(basisSet, diagCompound, true);
    else
      diagBlocks = computeDiagonalLibint(basisSet, diagCompound, true);

    for (size_t q = 0; q < NB; q++ )
    for (size_t p = 0; p <= q; p++ ) {
      size_t pq = p + q*NB;
      size_t pqComp = toCompound(p, q);
      diag[pq] = diagCompound[pqComp];
      diag[q + p*NB] = diagCompound[pqComp];
      if (diag[pq] >= tau_)
        D.push_back(pq);
    }

    if (D.empty())
      CErr("Cholesky decomposition threshold is greater than all TPI diagonal elements.");

    CQMemManager::get().free(diagCompound);

    std::sort(D.begin(), D.end(),
        [&diag](size_t pq, size_t rs){ return diag[pq] > diag[rs]; });

    double Dmax = diag[D[0]];

    double Q_threshold = std::max(sigma_ * Dmax, tau_);

    size_t QsplitIndex = std::distance(D.begin(),
        std::upper_bound(D.begin(), D.end(), Q_threshold,
            [&diag](double v, size_t rs){ return v > diag[rs]; }));

    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        DindicesByShell = groupPivotsByShell(basisSet, D);

    std::vector<SpanFactorShellPair> shellPair_Dindices;

    shellPair_Dindices.reserve(DindicesByShell.size());

    for (auto &kv : DindicesByShell)
      shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

    std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

    std::vector<size_t> Dgroup(D.size());

    size_t qL = 0, qR = QsplitIndex;
    for (auto &sD : shellPair_Dindices) {
      for (size_t &pq : sD.indices) {
        if (diag[D[pq]] >= Q_threshold) {
          Dgroup[qL] = D[pq];
          pq = qL++;
        } else {
          Dgroup[qR] = D[pq];
          pq = qR++;
        }
      }
    }

    D.swap(Dgroup);
    Dgroup.clear();

#ifdef CD_PROGRESS
    std::cout << "    D size = " << D.size() << std::endl;
#endif

    L = CQMemManager::get().calloc<double>(D.size() * std::min(D.size(), maxQual_));

    auto beginERIvec = topEffCDPivots;
    auto beginCDalgMM = topEffCDPivots;
    auto beginCDalgMV = topEffCDPivots;
    auto beginShrink = topEffCDPivots;
    double cumeERIvec = 0, cumeCDalgMM = 0, cumeCDalgMV = 0, cumeShrink = 0;
    size_t shrinkCount = 0;

    double curERIvec = 0.0, curCDalgMM = 0.0, curShrink = 0.0;
    double curERIdur = 0.0;
    size_t curERIcount = 0;

    while (D.size() > 0) {

      Dmax = diag[D[0]];


      // 4. Select Q
      size_t lenQ = std::min(maxQual_, QsplitIndex);
      double inner_threshold = Q_threshold;
      auto it = std::lower_bound(shellPair_Dindices.begin(), shellPair_Dindices.end(),
          lenQ,
          [](SpanFactorShellPair &a, size_t max_qual) {
            return a.indices[0] < max_qual;
          });
      if (it != shellPair_Dindices.end()) {
        lenQ = std::min(it->indices[0], lenQ);
        inner_threshold = std::max(diag[D[it->indices[0]]], inner_threshold);
        inner_threshold = std::max(diag[D[lenQ]], inner_threshold);
      }


#ifdef __DEBUGERI__
    std::cout << "D: [ ";
    for (size_t d : D)
      std::cout << d << ", ";
    std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    lenQ: " << lenQ << std::endl;
#endif

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "L", L, D.size(), pivots_.size(), D.size());
#endif


      // 5. Build Mpq (New version build in inner loop)
      double *M = CQMemManager::get().calloc<double>(D.size() * lenQ);
      curERIdur = t1ERI;
      curERIcount = c1ERI;
      beginERIvec = tick();
      computeSpanFactorERI(basisSet, D, shellPair_Dindices, lenQ, M, diagBlocks, true);
      curERIvec = tock(beginERIvec);
      curERIdur = t1ERI - curERIdur;
      curERIcount = c1ERI - curERIcount;
      cumeERIvec += curERIvec;

#ifdef CD_PROGRESS
      std::cout << "    Current ERI count    = " << curERIcount << std::endl;
      std::cout << "    Current ERI duration = " << curERIdur << " s " << std::endl;
      std::cout << "    Current ERIvec = " << curERIvec << " s " << std::endl;
#endif

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq ERI", M, D.size(), lenQ, D.size());
#endif

      beginCDalgMM = tick();

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, D.size(), lenQ, pivots_.size(),
           -1.0, L, D.size(),
           L, D.size(),
           1.0, M, D.size());

      curCDalgMM = tock(beginCDalgMM);
      cumeCDalgMM += curCDalgMM;

#ifdef CD_PROGRESS
      std::cout << "    Current CDalgMM = " << curCDalgMM << " s " << std::endl;
#endif

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq", M, D.size(), lenQ, D.size());
#endif

      // 6. Inner loop, select new pivots in C
      auto topInner = tick();
      std::vector<size_t> C;
      double innerCDalgMV = 0.0;

      while (true) {

        size_t qShellPairIndex = 0;
        size_t Qmax = 0;

        while (Qmax < lenQ and
               std::find(C.begin(), C.end(), D[Qmax]) != C.end())
          Qmax++;

        if (Qmax >= lenQ)
          break;

#ifdef __DEBUGERI__
        std::cout << "Qmax: " << Qmax << std::endl;
#endif

        for (size_t i = Qmax + 1; i < lenQ; i++ )
          if (diag[D[i]] > diag[D[Qmax]] and
              std::find(C.begin(), C.end(), D[i]) == C.end())
            Qmax = i;

        if (diag[D[Qmax]] < inner_threshold)
          break;

#ifdef CD_PROGRESS
        std::cout << "    Selected: (" << D[Qmax] % NB << "," << D[Qmax] / NB
                  << ")\tDiag: " << diag[D[Qmax]] << std::endl;
#endif

        for (size_t sumPair = 0;
             qShellPairIndex < shellPair_Dindices.size();
             qShellPairIndex++ ) {
          sumPair += shellPair_Dindices[qShellPairIndex].indices.size();
          if (sumPair > Qmax)
            break;
        }

        double *Lq = L + (pivots_.size() + C.size()) * D.size();

        beginCDalgMV = tick();
        std::copy_n(M + Qmax * D.size(), D.size(), Lq);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, D.size(), 1, C.size(),
             -1.0, L + pivots_.size() * D.size(), D.size(),
             L + Qmax + pivots_.size() * D.size(), D.size(),
             1.0, Lq, D.size());

        blas::scal(D.size(), 1.0 / sqrt(diag[D[Qmax]]), Lq, 1);

        innerCDalgMV += tock(beginCDalgMV);

#ifdef __DEBUGERI__
        prettyPrintSmart(std::cout, "Lq", Lq, D.size(), 1, D.size());
#endif

        C.push_back(D[Qmax]);

#ifdef __DEBUGERI__
        std::cout << "C: [ ";
        for (size_t c : C)
          std::cout << c << ", ";
        std::cout << "]" << std::endl;
#endif

        #pragma omp parallel for
        for (size_t p = 0; p < D.size(); p++)
          diag[D[p]] -= Lq[p] * Lq[p];

#ifdef __DEBUGERI__
        std::cout << "Diag: [ ";
        for (size_t d : D)
          std::cout << diag[d] << ", ";
        std::cout << "]" << std::endl;
#endif
      }
      auto durInner = tock(topInner);
      cumeCDalgMV += innerCDalgMV;

#ifdef CD_PROGRESS
      std::cout << "    InnerLoop CDalgMV  = " << innerCDalgMV << " s " << std::endl;
      std::cout << "    InnerLoop duration = " << durInner << " s " << std::endl;
#endif

      CQMemManager::get().free(M);


      // 7. Merge C to pivots_
      std::copy_n(C.begin(), C.size(), std::back_inserter(pivots_));

#ifdef __DEBUGERI__
      std::cout << "B: [ ";
      for (size_t p : pivots_)
        std::cout << p << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    pivot Size = " << pivots_.size() << std::endl;
#endif


      // 3. Update D and L
      beginShrink = tick();

#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

      std::vector<size_t> sortArgs(D.size());
      std::iota(sortArgs.begin(), sortArgs.end(), 0);
      std::sort(sortArgs.begin(), sortArgs.end(),
                [&diag, &D](size_t i, size_t j){ return diag[D[i]] > diag[D[j]]; });

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      sortArgs.resize(std::distance(sortArgs.begin(),
          std::upper_bound(sortArgs.begin(), sortArgs.end(), tau_,
              [diag, D](double t, size_t j){ return t > diag[D[j]]; })));

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      if (sortArgs.size() == 0)
        break;

      size_t DpreSize = D.size();

      std::vector<size_t> Dpre;
      Dpre.swap(D);
      D.reserve(sortArgs.size());

      for (size_t i : sortArgs) {
        D.push_back(Dpre[i]);
      }

      Dmax = diag[D[0]];
      Q_threshold = std::max(sigma_ * Dmax, tau_);

      QsplitIndex = std::distance(D.begin(),
          std::upper_bound(D.begin(), D.end(), Q_threshold,
              [&diag](double v, size_t rs){ return v > diag[rs]; }));

      DindicesByShell = groupPivotsByShell(basisSet, D);

      shellPair_Dindices.clear();

      for (auto &kv : DindicesByShell)
        shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

      std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

      std::vector<size_t> sortArgsShell(sortArgs.size());
      Dgroup.resize(D.size());

      size_t qL = 0, qR = QsplitIndex;
      for (auto &sD : shellPair_Dindices) {
        for (size_t &pq : sD.indices) {
          if (diag[D[pq]] >= Q_threshold) {
            Dgroup[qL] = D[pq];
            sortArgsShell[qL] = sortArgs[pq];
            pq = qL++;
          } else {
            Dgroup[qR] = D[pq];
            sortArgsShell[qR] = sortArgs[pq];
            pq = qR++;
          }
        }
      }

      sortArgs.swap(sortArgsShell);
      D.swap(Dgroup);
      Dgroup.clear();

#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    D size = " << D.size() << std::endl;
#endif

      double *Lnew = CQMemManager::get().calloc<double>(D.size() * (pivots_.size() + std::min(D.size(), maxQual_)));

      #pragma omp parallel for
      for (size_t p = 0; p < pivots_.size(); p++) {
        double *Lp = L + p * DpreSize, *Lnp = Lnew + p * D.size();
        size_t j = 0;
        for (size_t i : sortArgs) {
          Lnp[j++] = Lp[i];
        }
      }

      CQMemManager::get().free(L);

      L = Lnew;

      shrinkCount++;
      curShrink = tock(beginShrink);
      cumeShrink += curShrink;

#ifdef CD_PROGRESS
      std::cout << "    Current Shrink = " << curShrink << " s " << std::endl;
#endif

    }



    CQMemManager::get().free(diag, L);
    for (auto &kv : diagBlocks) {
      CQMemManager::get().free(kv.second);
    }

    auto durEffCDPivots = tock(topEffCDPivots);
    std::cout << "  Cholesky-RI-Pivots-ERI count    = " << c1ERI << std::endl;
    std::cout << "  Cholesky-RI-Pivots-ERI duration = " << t1ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIvec duration     = " << cumeERIvec << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMM duration    = " << cumeCDalgMM << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMV duration    = " << cumeCDalgMV << " s " << std::endl;
    std::cout << "  Cholesky-RI-Shrink count        = " << shrinkCount << std::endl;
    std::cout << "  Cholesky-RI-Shrink duration     = " << cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-misc duration       = "
              << durEffCDPivots - cumeERIvec - cumeCDalgMM - cumeCDalgMV - cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-Span-Factor-Pivots duration = " << durEffCDPivots << " s " << std::endl;

  }; // InCoreCholeskyRIERI<double>::computeCDPivots_SpanFactor


  /**
   *  \brief Compute the Cholesky RI pivots using the span factor algorithm
   *         and Libint2 over the CGTO basis.
   *         Algorithm refers to 10.1063/1.5083802
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeCDPivots_DynamicERI(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeCDPivots_DynamicERI(BasisSet &basisSet) {

    std::cout << "Dynamic-ERI Pivots Determination:" << std::endl;
    std::cout << bannerMid << std::endl;

    auto topEffCDPivots = tick();

    size_t NB2 = NB*NB;


    // 2. Clear pivots_ and CD vectors
    pivots_.clear();
    double* L;


    // 1. Compute diagonal
    // 3. Select all diagonals greater than theoreshold
    double *diag = CQMemManager::get().calloc<double>(NB2);
    double *diagCompound = CQMemManager::get().calloc<double>(NB*(NB+1)/2);
    std::map<std::pair<size_t, size_t>, double*> diagBlocks;
    std::vector<size_t> D;

    if (libcint_)
      diagBlocks = computeDiagonalLibcint(basisSet, diagCompound, true);
    else
      diagBlocks = computeDiagonalLibint(basisSet, diagCompound, true);

    for (size_t q = 0; q < NB; q++ )
    for (size_t p = 0; p <= q; p++ ) {
      size_t pq = p + q*NB;
      size_t pqComp = toCompound(p, q);
      diag[pq] = diagCompound[pqComp];
      diag[q + p*NB] = diagCompound[pqComp];
      if (diag[pq] >= tau_)
        D.push_back(pq);
    }

    if (D.empty())
      CErr("Cholesky decomposition threshold is greater than all TPI diagonal elements.");

    CQMemManager::get().free(diagCompound);

    std::sort(D.begin(), D.end(),
        [&diag](size_t pq, size_t rs){ return diag[pq] > diag[rs]; });

    double Dmax = diag[D[0]];

    double Q_threshold = std::max(sigma_ * Dmax, tau_);

    size_t QsplitIndex = std::distance(D.begin(),
        std::upper_bound(D.begin(), D.end(), Q_threshold,
            [&diag](double v, size_t rs){ return v > diag[rs]; }));

    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        DindicesByShell = groupPivotsByShell(basisSet, D);

    std::vector<SpanFactorShellPair> shellPair_Dindices;

    shellPair_Dindices.reserve(DindicesByShell.size());

    for (auto &kv : DindicesByShell)
      shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

    std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

    size_t lenD = D.size();
#ifdef CD_PROGRESS
    std::cout << "    D size = " << lenD << std::endl;
#endif

    L = CQMemManager::get().calloc<double>(lenD * std::min(lenD, maxQual_));

    auto beginERIvec = topEffCDPivots;
    auto beginERIcopy = topEffCDPivots;
    auto beginERItranspose = topEffCDPivots;
    auto beginCDalgMM = topEffCDPivots;
    auto beginCDalgMV = topEffCDPivots;
    auto beginShrink = topEffCDPivots;
    double cumeERIvec = 0, cumeERIcopy = 0, cumeERItranspose = 0,
        cumeCDalgMM = 0, cumeCDalgMV = 0, cumeShrink = 0;
    size_t shrinkCount = 0;

    // 4. Select Q
    size_t lenQ = std::min(maxQual_, QsplitIndex);
    double inner_threshold = Q_threshold;
    if (lenQ < lenD)
      inner_threshold = std::max(diag[D[lenQ]], Q_threshold);

    std::vector<size_t> Dgroup(lenD);

    size_t qL = 0, qR = lenQ;
    for (auto &sD : shellPair_Dindices) {
      for (size_t &pq : sD.indices) {
        if (pq < lenQ) {
          Dgroup[qL] = D[pq];
          pq = qL++;
        } else {
          Dgroup[qR] = D[pq];
          pq = qR++;
        }
      }
    }

    D.swap(Dgroup);
    Dgroup.clear();



#ifdef __DEBUGERI__
    std::cout << "D: [ ";
    for (size_t d : D)
      std::cout << d << ", ";
    std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
    std::cout << "    lenQ: " << lenQ << std::endl;
#endif

    auto QendIt = std::lower_bound(shellPair_Dindices.begin(), shellPair_Dindices.end(),
        lenQ,
        [](SpanFactorShellPair &a, size_t max_qual) {
          return a.indices[0] < max_qual;
        });

    size_t nERIvec = 0;
    for (auto it = shellPair_Dindices.begin(); it < QendIt; it++) {
      nERIvec += it->indices.size();
    }

#ifdef CD_PROGRESS
    std::cout << "    ERI size = " << nERIvec << std::endl;
    std::cout << "    LevalEnd = " << lenQ << std::endl;
    std::cout << "    RevalEnd = " << nERIvec << std::endl;
    std::cout << "    SevalBegin = " << nERIvec << std::endl;
#endif

    double *ERIvecAlloc = CQMemManager::get().calloc<double>(lenD * nERIvec);
    double curERIvec = 0.0, curCDalgMM = 0.0, curShrink = 0.0;
    double curERIcopy = 0.0, curERItranspose = 0.0, curERIdur = 0.0;
    size_t curERIcount = 0;

    while (lenD > 0) {

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "L", L, lenD, pivots_.size(), lenD);
#endif

      // 5. Build Mpq (New version build in inner loop)
      double *M = CQMemManager::get().calloc<double>(lenD * lenQ);
      curERIdur = t1ERI;
      curERIcount = c1ERI;
      beginERIvec = tick();
      computeSpanFactorERI(basisSet, D, shellPair_Dindices, lenQ, ERIvecAlloc, diagBlocks, false);
      curERIvec = tock(beginERIvec);
      curERIdur = t1ERI - curERIdur;
      curERIcount = c1ERI - curERIcount;
      cumeERIvec += curERIvec;

#ifdef CD_PROGRESS
      std::cout << "    Current ERI count    = " << curERIcount << std::endl;
      std::cout << "    Current ERI duration = " << curERIdur << " s " << std::endl;
      std::cout << "    Current ERIvec = " << curERIvec << " s " << std::endl;
#endif

      beginERIcopy = tick();
      SetMat('N', lenD, lenQ, 1.0, ERIvecAlloc, lenD, M, lenD);
      curERIcopy = tock(beginERIcopy);
      cumeERIcopy += curERIcopy;

#ifdef CD_PROGRESS
      std::cout << "    Current ERIcopy = " << curERIcopy << " s " << std::endl;
#endif

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq ERI", M, lenD, lenQ, lenD);
#endif

      beginCDalgMM = tick();

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD, lenQ, pivots_.size(),
           -1.0, L, lenD,
           L, lenD,
           1.0, M, lenD);

      curCDalgMM = tock(beginCDalgMM);
      cumeCDalgMM += curCDalgMM;

#ifdef CD_PROGRESS
      std::cout << "    Current CDalgMM = " << curCDalgMM << " s " << std::endl;
#endif

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq", M, lenD, lenQ, lenD);
#endif

      // 6. Inner loop, select new pivots in C
      auto topInner = tick();
      std::vector<size_t> C;
      double innerCDalgMV = 0.0;

      while (true) {

        size_t qShellPairIndex = 0;
        size_t Qmax = 0;

        while (Qmax < lenQ and
               std::find(C.begin(), C.end(), D[Qmax]) != C.end())
          Qmax++;

        if (Qmax >= lenQ)
          break;

#ifdef __DEBUGERI__
        std::cout << "Qmax: " << Qmax << std::endl;
#endif

        for (size_t i = Qmax + 1; i < lenQ; i++ )
          if (diag[D[i]] > diag[D[Qmax]] and
              std::find(C.begin(), C.end(), D[i]) == C.end())
            Qmax = i;

        if (diag[D[Qmax]] < inner_threshold)
          break;

#ifdef CD_PROGRESS
        std::cout << "    Selected: (" << D[Qmax] % NB << "," << D[Qmax] / NB
                  << ")\tDiag: " << diag[D[Qmax]] << std::endl;
#endif

        for (size_t sumPair = 0;
             qShellPairIndex < shellPair_Dindices.size();
             qShellPairIndex++ ) {
          sumPair += shellPair_Dindices[qShellPairIndex].indices.size();
          if (sumPair > Qmax)
            break;
        }

        double *Lq = L + (pivots_.size() + C.size()) * lenD;

        beginCDalgMV = tick();
        std::copy_n(M + Qmax * lenD, lenD, Lq);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD, 1, C.size(),
             -1.0, L + pivots_.size() * lenD, lenD,
             L + Qmax + pivots_.size() * lenD, lenD,
             1.0, Lq, lenD);

        blas::scal(lenD, 1.0 / sqrt(diag[D[Qmax]]), Lq, 1);

        innerCDalgMV += tock(beginCDalgMV);

#ifdef __DEBUGERI__
        prettyPrintSmart(std::cout, "Lq", Lq, lenD, 1, lenD);
#endif

        C.push_back(D[Qmax]);

#ifdef __DEBUGERI__
        std::cout << "C: [ ";
        for (size_t c : C)
          std::cout << c << ", ";
        std::cout << "]" << std::endl;
#endif

        #pragma omp parallel for
        for (size_t p = 0; p < lenD; p++)
          diag[D[p]] -= Lq[p] * Lq[p];

#ifdef __DEBUGERI__
        std::cout << "Diag: [ ";
        for (size_t d : D)
          std::cout << diag[d] << ", ";
        std::cout << "]" << std::endl;
#endif
      }
      auto durInner = tock(topInner);
      cumeCDalgMV += innerCDalgMV;

#ifdef CD_PROGRESS
      std::cout << "    InnerLoop CDalgMV  = " << innerCDalgMV << " s " << std::endl;
      std::cout << "    InnerLoop duration = " << durInner << " s " << std::endl;
#endif

      CQMemManager::get().free(M);


      // 7. Merge C to pivots_
      std::copy_n(C.begin(), C.size(), std::back_inserter(pivots_));

#ifdef __DEBUGERI__
      std::cout << "B: [ ";
      for (size_t p : pivots_)
        std::cout << p << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    pivot Size = " << pivots_.size() << std::endl;
#endif


      // 3. Update D and L
      beginShrink = tick();

#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

      std::vector<size_t> sortArgs(lenD);
      std::iota(sortArgs.begin(), sortArgs.end(), 0);
      std::sort(sortArgs.begin(), sortArgs.end(),
                [&diag, &D](size_t i, size_t j){ return diag[D[i]] > diag[D[j]]; });

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      sortArgs.resize(std::distance(sortArgs.begin(),
          std::upper_bound(sortArgs.begin(), sortArgs.end(), tau_,
              [diag, D](double t, size_t j){ return t > diag[D[j]]; })));

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      size_t preLenD = lenD;
      lenD = sortArgs.size();

      if (lenD == 0)
        break;

      std::vector<size_t> Dpre;
      Dpre.swap(D);
      D.reserve(sortArgs.size());

      for (size_t i : sortArgs) {
        D.push_back(Dpre[i]);
      }

      Dmax = diag[D[0]];
      Q_threshold = std::max(sigma_ * Dmax, tau_);

      QsplitIndex = std::distance(D.begin(),
          std::upper_bound(D.begin(), D.end(), Q_threshold,
              [&diag](double v, size_t rs){ return v > diag[rs]; }));

      // 4. Select Q
      lenQ = std::min(maxQual_, QsplitIndex);
      inner_threshold = Q_threshold;
      if (lenQ < lenD)
        inner_threshold = std::max(diag[D[lenQ]], Q_threshold);

      DindicesByShell = groupPivotsByShell(basisSet, D);

      shellPair_Dindices.clear();

      for (auto &kv : DindicesByShell)
        shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

      for (auto &sD : shellPair_Dindices) {
        if (sortArgs[sD.indices[0]] < nERIvec) {
          sD.evaluated = true;
        }
      }

      std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

      QendIt = std::lower_bound(shellPair_Dindices.begin(), shellPair_Dindices.end(),
          lenQ,
          [](SpanFactorShellPair &a, size_t max_qual) {
            return a.indices[0] < max_qual;
          });

      std::sort(shellPair_Dindices.begin(), QendIt,
          [](SpanFactorShellPair &a, SpanFactorShellPair &b) {
            return a.evaluated == b.evaluated ? false : a.evaluated;
          });

      std::sort(QendIt, shellPair_Dindices.end(),
          [](SpanFactorShellPair &a, SpanFactorShellPair &b) {
            return a.evaluated == b.evaluated ? false : a.evaluated;
          });

      std::vector<size_t> sortArgsShell(sortArgs.size());
      Dgroup.resize(lenD);

      size_t preNERIvec = nERIvec;
      nERIvec = 0;
      for (auto it = shellPair_Dindices.begin(); it != shellPair_Dindices.end(); it++) {
        if (it < QendIt or it->evaluated)
          nERIvec += it->indices.size();
      }

      size_t qL = 0, qR = lenQ, LevalEnd = lenQ, RevalEnd = nERIvec, SevalBegin = nERIvec;
      bool preEvaluated = true, resumed = false;
      for (auto &sD : shellPair_Dindices) {
        if (not resumed and preEvaluated and not sD.evaluated) {
          LevalEnd = qL;
          RevalEnd = qR;
          preEvaluated = false;
        }
        if (not preEvaluated and sD.evaluated) {
          SevalBegin = qR;
          preEvaluated = true;
          resumed = true;
        }
        for (size_t &pq : sD.indices) {
          if (pq < lenQ) {
            Dgroup[qL] = D[pq];
            sortArgsShell[qL] = sortArgs[pq];
            pq = qL++;
          } else {
            Dgroup[qR] = D[pq];
            sortArgsShell[qR] = sortArgs[pq];
            pq = qR++;
          }
        }
      }

      sortArgs.swap(sortArgsShell);
      D.swap(Dgroup);
      Dgroup.clear();


#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    D size = " << lenD << std::endl;
#endif

      double *Lnew = CQMemManager::get().calloc<double>(lenD * (pivots_.size() + std::min(lenD, maxQual_)));

      #pragma omp parallel for
      for (size_t p = 0; p < pivots_.size(); p++) {
        double *Lp = L + p * preLenD, *Lnp = Lnew + p * lenD;
        size_t j = 0;
        for (size_t i : sortArgs) {
          Lnp[j++] = Lp[i];
        }
      }

      CQMemManager::get().free(L);
      L = Lnew;

      double *preERIvecAlloc = ERIvecAlloc;
      ERIvecAlloc = CQMemManager::get().calloc<double>(lenD * nERIvec);

      #pragma omp parallel for
      for (size_t p = 0; p < nERIvec; p++) {
        if (sortArgs[p] < preNERIvec) {
          double *ERIp = preERIvecAlloc + sortArgs[p] * preLenD;
          double *ERInp = ERIvecAlloc + p * lenD;
          size_t j = 0;
          for (size_t i : sortArgs) {
            ERInp[j++] = ERIp[i];
          }
        }
      }

      CQMemManager::get().free(preERIvecAlloc);

      shrinkCount++;
      curShrink = tock(beginShrink);
      cumeShrink += curShrink;



#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    lenQ: " << lenQ << std::endl;
      std::cout << "    ERI size = " << nERIvec << std::endl;
      std::cout << "    LevalEnd = " << LevalEnd << std::endl;
      std::cout << "    RevalEnd = " << RevalEnd << std::endl;
      std::cout << "    SevalBegin = " << SevalBegin << std::endl;
      std::cout << "    Current Shrink = " << curShrink << " s " << std::endl;
#endif

      // Transpose existing ERI elements
      beginERItranspose = tick();

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "ERI vecs before transpose", ERIvecAlloc, lenD, nERIvec, lenD);
#endif

      SetMat('C', lenQ - LevalEnd, LevalEnd,
             1.0, ERIvecAlloc + LevalEnd, lenD,
             ERIvecAlloc + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, LevalEnd,
             1.0, ERIvecAlloc + RevalEnd, lenD,
             ERIvecAlloc + RevalEnd * lenD, lenD);
      SetMat('C', lenQ - LevalEnd, RevalEnd - lenQ,
             1.0, ERIvecAlloc + LevalEnd + lenQ * lenD, lenD,
             ERIvecAlloc + lenQ + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, RevalEnd - lenQ,
             1.0, ERIvecAlloc + RevalEnd + lenQ * lenD, lenD,
             ERIvecAlloc + lenQ + RevalEnd * lenD, lenD);
      SetMat('C', lenQ - LevalEnd, nERIvec - SevalBegin,
             1.0, ERIvecAlloc + LevalEnd + SevalBegin * lenD, lenD,
             ERIvecAlloc + SevalBegin + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, nERIvec - SevalBegin,
             1.0, ERIvecAlloc + RevalEnd + SevalBegin * lenD, lenD,
             ERIvecAlloc + SevalBegin + RevalEnd * lenD, lenD);

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "ERI vecs after transpose", ERIvecAlloc, lenD, nERIvec, lenD);
#endif
      curERItranspose = tock(beginERItranspose);
      cumeERItranspose += curERItranspose;

#ifdef CD_PROGRESS
      std::cout << "    Current ERItrans = " << curERItranspose << " s " << std::endl;
#endif

    }



    CQMemManager::get().free(diag, L, ERIvecAlloc);
    for (auto &kv : diagBlocks) {
      CQMemManager::get().free(kv.second);
    }

    auto durEffCDPivots = tock(topEffCDPivots);
    std::cout << "  Cholesky-RI-Pivots-ERI count    = " << c1ERI << std::endl;
    std::cout << "  Cholesky-RI-Pivots-ERI duration = " << t1ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIvec duration     = " << cumeERIvec << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIcopy duration    = " << cumeERIcopy << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERItrans duration   = " << cumeERItranspose << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMM duration    = " << cumeCDalgMM << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMV duration    = " << cumeCDalgMV << " s " << std::endl;
    std::cout << "  Cholesky-RI-Shrink count        = " << shrinkCount << std::endl;
    std::cout << "  Cholesky-RI-Shrink duration     = " << cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-misc duration       = "
              << durEffCDPivots - cumeERIvec - cumeERIcopy - cumeERItranspose
                 - cumeCDalgMM - cumeCDalgMV - cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-Dynamic-ERI-Pivots duration = " << durEffCDPivots << " s " << std::endl;

  }; // InCoreCholeskyRIERI<double>::computeCDPivots_DynamicERI


  /**
   *  \brief Compute the Cholesky RI pivots using the span factor algorithm
   *         and Libint2 over the CGTO basis.
   *         Algorithm refers to 10.1063/1.5083802
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeCDPivots_SpanFactorReuse(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeCDPivots_SpanFactorReuse(BasisSet &basisSet) {

    std::cout << "Span-Factor-Reuse Pivots Determination:" << std::endl;
    std::cout << bannerMid << std::endl;

    auto topEffCDPivots = tick();

    size_t NB2 = NB*NB;


    // 2. Clear pivots_ and CD vectors
    pivots_.clear();
    double* L;


    // 1. Compute diagonal
    // 3. Select all diagonals greater than theoreshold
    double *diag = CQMemManager::get().calloc<double>(NB2);
    double *diagCompound = CQMemManager::get().calloc<double>(NB*(NB+1)/2);
    std::map<std::pair<size_t, size_t>, double*> diagBlocks;
    std::vector<size_t> D;

    if (libcint_)
      diagBlocks = computeDiagonalLibcint(basisSet, diagCompound, true);
    else
      diagBlocks = computeDiagonalLibint(basisSet, diagCompound, true);

    for (size_t q = 0; q < NB; q++ )
    for (size_t p = 0; p <= q; p++ ) {
      size_t pq = p + q*NB;
      size_t pqComp = toCompound(p, q);
      diag[pq] = diagCompound[pqComp];
      diag[q + p*NB] = diagCompound[pqComp];
      if (diag[pq] >= tau_)
        D.push_back(pq);
    }

    if (D.empty())
      CErr("Cholesky decomposition threshold is greater than all TPI diagonal elements.");

    CQMemManager::get().free(diagCompound);

    std::sort(D.begin(), D.end(),
        [&diag](size_t pq, size_t rs){ return diag[pq] > diag[rs]; });

    double Dmax = diag[D[0]];

    double Q_threshold = std::max(sigma_ * Dmax, tau_);

    size_t QsplitIndex = std::distance(D.begin(),
        std::upper_bound(D.begin(), D.end(), Q_threshold,
            [&diag](double v, size_t rs){ return v > diag[rs]; }));

    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        DindicesByShell = groupPivotsByShell(basisSet, D);

    std::vector<SpanFactorShellPair> shellPair_Dindices;

    shellPair_Dindices.reserve(DindicesByShell.size());

    for (auto &kv : DindicesByShell)
      shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

    std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

    size_t lenD = D.size();
#ifdef CD_PROGRESS
    std::cout << "    D size = " << lenD << std::endl;
#endif

    L = CQMemManager::get().calloc<double>(lenD * std::min(lenD, maxQual_));

    auto beginERIvec = topEffCDPivots;
    auto beginERIcopy = topEffCDPivots;
    auto beginERItranspose = topEffCDPivots;
    auto beginCDalgMM = topEffCDPivots;
    auto beginCDalgMV = topEffCDPivots;
    auto beginShrink = topEffCDPivots;
    double cumeERIvec = 0, cumeERIcopy = 0, cumeERItranspose = 0,
        cumeCDalgMM = 0, cumeCDalgMV = 0, cumeShrink = 0;
    size_t shrinkCount = 0;

    // 4. Select Q
    size_t lenQ = std::min(maxQual_, QsplitIndex);
    double inner_threshold = Q_threshold;
    if (lenQ < lenD)
      inner_threshold = std::max(diag[D[lenQ]], Q_threshold);

    std::vector<size_t> Dgroup(lenD);

    size_t qL = 0, qR = lenQ;
    for (auto &sD : shellPair_Dindices) {
      for (size_t &pq : sD.indices) {
        if (pq < lenQ) {
          Dgroup[qL] = D[pq];
          pq = qL++;
        } else {
          Dgroup[qR] = D[pq];
          pq = qR++;
        }
      }
    }

    D.swap(Dgroup);
    Dgroup.clear();



#ifdef __DEBUGERI__
    std::cout << "D: [ ";
    for (size_t d : D)
      std::cout << d << ", ";
    std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
    std::cout << "    lenQ: " << lenQ << std::endl;
#endif

    auto QendIt = std::lower_bound(shellPair_Dindices.begin(), shellPair_Dindices.end(),
        lenQ,
        [](SpanFactorShellPair &a, size_t max_qual) {
          return a.indices[0] < max_qual;
        });

    size_t nERIvec = 0;
    for (auto it = shellPair_Dindices.begin(); it < QendIt; it++) {
      nERIvec += it->indices.size();
    }

#ifdef CD_PROGRESS
    std::cout << "    ERI size = " << nERIvec << std::endl;
    std::cout << "    LevalEnd = " << lenQ << std::endl;
    std::cout << "    RevalEnd = " << nERIvec << std::endl;
    std::cout << "    SevalBegin = " << nERIvec << std::endl;
#endif

    double *ERIvecAlloc = CQMemManager::get().calloc<double>(lenD * nERIvec);
    double curERIvec = 0.0, curCDalgMM = 0.0, curShrink = 0.0;
    double curERIcopy = 0.0, curERItranspose = 0.0, curERIdur = 0.0;
    size_t curERIcount = 0;

    size_t preLenB = 0;
    size_t LevalEnd = 0, RevalEnd = lenQ, SevalBegin = nERIvec;

    while (lenD > 0) {

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "L", L, lenD, pivots_.size(), lenD);
#endif

      // Compute ERI
      curERIdur = t1ERI;
      curERIcount = c1ERI;
      beginERIvec = tick();
      computeSpanFactorERI(basisSet, D, shellPair_Dindices, lenQ, ERIvecAlloc, diagBlocks, false);
      curERIvec = tock(beginERIvec);
      curERIdur = t1ERI - curERIdur;
      curERIcount = c1ERI - curERIcount;
      cumeERIvec += curERIvec;

#ifdef CD_PROGRESS
      std::cout << "    Current ERI count    = " << curERIcount << std::endl;
      std::cout << "    Current ERI duration = " << curERIdur << " s " << std::endl;
      std::cout << "    Current ERIvec = " << curERIvec << " s " << std::endl;
#endif

      // Bring new ERI to preLenB level: 6 block GEMMs
      beginCDalgMM = tick();

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenQ - LevalEnd, lenQ - LevalEnd, preLenB,
           -1.0, L + LevalEnd, lenD,
           L + LevalEnd, lenD,
           1.0, ERIvecAlloc + LevalEnd + LevalEnd * lenD, lenD);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenQ - LevalEnd, SevalBegin - RevalEnd, preLenB,
           -1.0, L + LevalEnd, lenD,
           L + RevalEnd, lenD,
           1.0, ERIvecAlloc + LevalEnd + RevalEnd * lenD, lenD);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, SevalBegin - RevalEnd, lenQ - LevalEnd, preLenB,
           -1.0, L + RevalEnd, lenD,
           L + LevalEnd, lenD,
           1.0, ERIvecAlloc + RevalEnd + LevalEnd * lenD, lenD);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, SevalBegin - RevalEnd, SevalBegin - RevalEnd, preLenB,
           -1.0, L + RevalEnd, lenD,
           L + RevalEnd, lenD,
           1.0, ERIvecAlloc + RevalEnd + RevalEnd * lenD, lenD);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD - nERIvec, lenQ - LevalEnd, preLenB,
           -1.0, L + nERIvec, lenD,
           L + LevalEnd, lenD,
           1.0, ERIvecAlloc + nERIvec + LevalEnd * lenD, lenD);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD - nERIvec, SevalBegin - RevalEnd, preLenB,
           -1.0, L + nERIvec, lenD,
           L + RevalEnd, lenD,
           1.0, ERIvecAlloc + nERIvec + RevalEnd * lenD, lenD);

      // 5. Build Mpq (New version build in inner loop)

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq ERI", ERIvecAlloc, lenD, nERIvec, lenD);
#endif

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD, nERIvec, pivots_.size() - preLenB,
           -1.0, L + preLenB * lenD, lenD,
           L + preLenB * lenD, lenD,
           1.0, ERIvecAlloc, lenD);

      curCDalgMM = tock(beginCDalgMM);
      cumeCDalgMM += curCDalgMM;

#ifdef CD_PROGRESS
      std::cout << "    Current CDalgMM = " << curCDalgMM << " s " << std::endl;
#endif

      preLenB = pivots_.size();

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "Mpq", ERIvecAlloc, lenD, nERIvec, lenD);
#endif

      // 6. Inner loop, select new pivots in C
      auto topInner = tick();
      std::vector<size_t> C;
      double innerCDalgMV = 0.0;

      while (true) {

        size_t qShellPairIndex = 0;
        size_t Qmax = 0;

        while (Qmax < lenQ and
               std::find(C.begin(), C.end(), D[Qmax]) != C.end())
          Qmax++;

        if (Qmax >= lenQ)
          break;

#ifdef __DEBUGERI__
        std::cout << "Qmax: " << Qmax << std::endl;
#endif

        for (size_t i = Qmax + 1; i < lenQ; i++ )
          if (diag[D[i]] > diag[D[Qmax]] and
              std::find(C.begin(), C.end(), D[i]) == C.end())
            Qmax = i;

        if (diag[D[Qmax]] < inner_threshold)
          break;

#ifdef CD_PROGRESS
        std::cout << "    Selected: (" << D[Qmax] % NB << "," << D[Qmax] / NB
                  << ")\tDiag: " << diag[D[Qmax]] << std::endl;
#endif

        for (size_t sumPair = 0;
             qShellPairIndex < shellPair_Dindices.size();
             qShellPairIndex++ ) {
          sumPair += shellPair_Dindices[qShellPairIndex].indices.size();
          if (sumPair > Qmax)
            break;
        }

        double *Lq = L + (pivots_.size() + C.size()) * lenD;

        beginCDalgMV = tick();
        std::copy_n(ERIvecAlloc + Qmax * lenD, lenD, Lq);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, lenD, 1, C.size(),
             -1.0, L + pivots_.size() * lenD, lenD,
             L + Qmax + pivots_.size() * lenD, lenD,
             1.0, Lq, lenD);

        blas::scal(lenD, 1.0 / sqrt(diag[D[Qmax]]), Lq, 1);

        innerCDalgMV += tock(beginCDalgMV);

#ifdef __DEBUGERI__
        prettyPrintSmart(std::cout, "Lq", Lq, lenD, 1, lenD);
#endif

        C.push_back(D[Qmax]);

#ifdef __DEBUGERI__
        std::cout << "C: [ ";
        for (size_t c : C)
          std::cout << c << ", ";
        std::cout << "]" << std::endl;
#endif

        #pragma omp parallel for
        for (size_t p = 0; p < lenD; p++)
          diag[D[p]] -= Lq[p] * Lq[p];

#ifdef __DEBUGERI__
        std::cout << "Diag: [ ";
        for (size_t d : D)
          std::cout << diag[d] << ", ";
        std::cout << "]" << std::endl;
#endif
      }
      auto durInner = tock(topInner);
      cumeCDalgMV += innerCDalgMV;

#ifdef CD_PROGRESS
      std::cout << "    InnerLoop CDalgMV  = " << innerCDalgMV << " s " << std::endl;
      std::cout << "    InnerLoop duration = " << durInner << " s " << std::endl;
#endif


      // 7. Merge C to pivots_
      std::copy_n(C.begin(), C.size(), std::back_inserter(pivots_));

#ifdef __DEBUGERI__
      std::cout << "B: [ ";
      for (size_t p : pivots_)
        std::cout << p << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    pivot Size = " << pivots_.size() << std::endl;
#endif


      // 3. Update D and L
      beginShrink = tick();

#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

      std::vector<size_t> sortArgs(lenD);
      std::iota(sortArgs.begin(), sortArgs.end(), 0);
      std::sort(sortArgs.begin(), sortArgs.end(),
                [&diag, &D](size_t i, size_t j){ return diag[D[i]] > diag[D[j]]; });

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      sortArgs.resize(std::distance(sortArgs.begin(),
          std::upper_bound(sortArgs.begin(), sortArgs.end(), tau_,
              [diag, D](double t, size_t j){ return t > diag[D[j]]; })));

#ifdef __DEBUGERI__
      std::cout << "sortArgs: [ " << std::endl;
      for (size_t s : sortArgs)
        std::cout << s << ", " << diag[D[s]] << std::endl;
      std::cout << "]" << std::endl;
#endif

      size_t preLenD = lenD;
      lenD = sortArgs.size();

      if (lenD == 0)
        break;

      std::vector<size_t> Dpre;
      Dpre.swap(D);
      D.reserve(sortArgs.size());

      for (size_t i : sortArgs) {
        D.push_back(Dpre[i]);
      }

      Dmax = diag[D[0]];
      Q_threshold = std::max(sigma_ * Dmax, tau_);

      QsplitIndex = std::distance(D.begin(),
          std::upper_bound(D.begin(), D.end(), Q_threshold,
              [&diag](double v, size_t rs){ return v > diag[rs]; }));

      // 4. Select Q
      lenQ = std::min(maxQual_, QsplitIndex);
      inner_threshold = Q_threshold;
      if (lenQ < lenD)
        inner_threshold = std::max(diag[D[lenQ]], Q_threshold);

      DindicesByShell = groupPivotsByShell(basisSet, D);

      shellPair_Dindices.clear();

      for (auto &kv : DindicesByShell)
        shellPair_Dindices.push_back(SpanFactorShellPair{kv.first, kv.second, false});

      for (auto &sD : shellPair_Dindices) {
        if (sortArgs[sD.indices[0]] < nERIvec) {
          sD.evaluated = true;
        }
      }

      std::sort(shellPair_Dindices.begin(), shellPair_Dindices.end());

      QendIt = std::lower_bound(shellPair_Dindices.begin(), shellPair_Dindices.end(),
          lenQ,
          [](SpanFactorShellPair &a, size_t max_qual) {
            return a.indices[0] < max_qual;
          });

      std::sort(shellPair_Dindices.begin(), QendIt,
          [](SpanFactorShellPair &a, SpanFactorShellPair &b) {
            return a.evaluated == b.evaluated ? false : a.evaluated;
          });

      std::sort(QendIt, shellPair_Dindices.end(),
          [](SpanFactorShellPair &a, SpanFactorShellPair &b) {
            return a.evaluated == b.evaluated ? false : a.evaluated;
          });

      std::vector<size_t> sortArgsShell(sortArgs.size());
      Dgroup.resize(lenD);

      size_t preNERIvec = nERIvec;
      nERIvec = 0;
      for (auto it = shellPair_Dindices.begin(); it != shellPair_Dindices.end(); it++) {
        if (it < QendIt or it->evaluated)
          nERIvec += it->indices.size();
      }

      qL = 0, qR = lenQ, LevalEnd = lenQ, RevalEnd = nERIvec, SevalBegin = nERIvec;
      bool preEvaluated = true, resumed = false;
      for (auto &sD : shellPair_Dindices) {
        if (not resumed and preEvaluated and not sD.evaluated) {
          LevalEnd = qL;
          RevalEnd = qR;
          preEvaluated = false;
        }
        if (not preEvaluated and sD.evaluated) {
          SevalBegin = qR;
          preEvaluated = true;
          resumed = true;
        }
        for (size_t &pq : sD.indices) {
          if (pq < lenQ) {
            Dgroup[qL] = D[pq];
            sortArgsShell[qL] = sortArgs[pq];
            pq = qL++;
          } else {
            Dgroup[qR] = D[pq];
            sortArgsShell[qR] = sortArgs[pq];
            pq = qR++;
          }
        }
      }

      sortArgs.swap(sortArgsShell);
      D.swap(Dgroup);
      Dgroup.clear();


#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    D size = " << lenD << std::endl;
#endif

      double *Lnew = CQMemManager::get().calloc<double>(lenD * (pivots_.size() + std::min(lenD, maxQual_)));

      #pragma omp parallel for
      for (size_t p = 0; p < pivots_.size(); p++) {
        double *Lp = L + p * preLenD, *Lnp = Lnew + p * lenD;
        size_t j = 0;
        for (size_t i : sortArgs) {
          Lnp[j++] = Lp[i];
        }
      }

      CQMemManager::get().free(L);
      L = Lnew;

      double *preERIvecAlloc = ERIvecAlloc;
      ERIvecAlloc = CQMemManager::get().calloc<double>(lenD * nERIvec);

      #pragma omp parallel for
      for (size_t p = 0; p < nERIvec; p++) {
        if (sortArgs[p] < preNERIvec) {
          double *ERIp = preERIvecAlloc + sortArgs[p] * preLenD;
          double *ERInp = ERIvecAlloc + p * lenD;
          size_t j = 0;
          for (size_t i : sortArgs) {
            ERInp[j++] = ERIp[i];
          }
        }
      }

      CQMemManager::get().free(preERIvecAlloc);

      shrinkCount++;
      curShrink = tock(beginShrink);
      cumeShrink += curShrink;


#ifdef __DEBUGERI__
      std::cout << "D: [ ";
      for (size_t d : D)
        std::cout << d << ", ";
      std::cout << "]" << std::endl;
#endif

#ifdef CD_PROGRESS
      std::cout << "    lenQ: " << lenQ << std::endl;
      std::cout << "    ERI size = " << nERIvec << std::endl;
      std::cout << "    LevalEnd = " << LevalEnd << std::endl;
      std::cout << "    RevalEnd = " << RevalEnd << std::endl;
      std::cout << "    SevalBegin = " << SevalBegin << std::endl;
      std::cout << "    Current Shrink = " << curShrink << " s " << std::endl;
#endif

      // Transpose existing ERI elements
      beginERItranspose = tick();

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "ERI vecs before transpose", ERIvecAlloc, lenD, nERIvec, lenD);
#endif

      SetMat('C', lenQ - LevalEnd, LevalEnd,
             1.0, ERIvecAlloc + LevalEnd, lenD,
             ERIvecAlloc + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, LevalEnd,
             1.0, ERIvecAlloc + RevalEnd, lenD,
             ERIvecAlloc + RevalEnd * lenD, lenD);
      SetMat('C', lenQ - LevalEnd, RevalEnd - lenQ,
             1.0, ERIvecAlloc + LevalEnd + lenQ * lenD, lenD,
             ERIvecAlloc + lenQ + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, RevalEnd - lenQ,
             1.0, ERIvecAlloc + RevalEnd + lenQ * lenD, lenD,
             ERIvecAlloc + lenQ + RevalEnd * lenD, lenD);
      SetMat('C', lenQ - LevalEnd, nERIvec - SevalBegin,
             1.0, ERIvecAlloc + LevalEnd + SevalBegin * lenD, lenD,
             ERIvecAlloc + SevalBegin + LevalEnd * lenD, lenD);
      SetMat('C', SevalBegin - RevalEnd, nERIvec - SevalBegin,
             1.0, ERIvecAlloc + RevalEnd + SevalBegin * lenD, lenD,
             ERIvecAlloc + SevalBegin + RevalEnd * lenD, lenD);

#ifdef __DEBUGERI__
      prettyPrintSmart(std::cout, "ERI vecs after transpose", ERIvecAlloc, lenD, nERIvec, lenD);
#endif
      curERItranspose = tock(beginERItranspose);
      cumeERItranspose += curERItranspose;

#ifdef CD_PROGRESS
      std::cout << "    Current ERItrans = " << curERItranspose << " s " << std::endl;
#endif

    }



    CQMemManager::get().free(diag, L, ERIvecAlloc);
    for (auto &kv : diagBlocks) {
      CQMemManager::get().free(kv.second);
    }

    auto durEffCDPivots = tock(topEffCDPivots);
    std::cout << "  Cholesky-RI-Pivots-ERI count    = " << c1ERI << std::endl;
    std::cout << "  Cholesky-RI-Pivots-ERI duration = " << t1ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIvec duration     = " << cumeERIvec << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIcopy duration    = " << cumeERIcopy << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERItrans duration   = " << cumeERItranspose << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMM duration    = " << cumeCDalgMM << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalgMV duration    = " << cumeCDalgMV << " s " << std::endl;
    std::cout << "  Cholesky-RI-Shrink count        = " << shrinkCount << std::endl;
    std::cout << "  Cholesky-RI-Shrink duration     = " << cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-misc duration       = "
              << durEffCDPivots - cumeERIvec - cumeERIcopy - cumeERItranspose
                 - cumeCDalgMM - cumeCDalgMV - cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-Span-Factor-Reuse-Pivots duration = " << durEffCDPivots << " s " << std::endl;

  }; // InCoreCholeskyRIERI<double>::computeCDPivots_SpanFactorReuse



  /**
   * @brief The BisectMemManager class
   *        A memory manager class specially designed for computing
   *        Cholesky decomposition pivots.
   *        When the memory requirement of each vector is less than
   *        half of the original, number of available memory vectors
   *        double.
   * @warning Invalid for memVecLen > 2^32
   */
  template <typename T>
  class BisectMemManager {
  protected:
    std::vector<T*> rawPtrs_;
    std::vector<std::vector<size_t>> memPtrs_;
    std::set<T*> freePtrs_;
    size_t initialMemVecLen_;
    size_t initialNVec_;
    size_t memBlockSize_;
    size_t bisectCount_ = 0;
    size_t allocCount_ = 0;

    size_t blockIndex(T *ptr) const {
      size_t I = 0;
      for (T *r : rawPtrs_) {
        if (ptr >= r and ptr < r + memBlockSize_)
          return I;
        I++;
      }
      return I;
    }

    const std::array<std::pair<size_t,size_t>,5> masks{
      {{0xaaaaaaaa, 0x55555555},
      {0xcccccccc, 0x33333333},
      {0xf0f0f0f0, 0xf0f0f0f},
      {0xff00ff00, 0xff00ff},
      {0xffff0000, 0xffff}}
    };

    size_t reverseLastBisectCountBits(size_t n) const {
      size_t m = 0;
      while ((1<<m) < bisectCount_) {
        n = ((n & masks[m].first) >> (1<<m)) | ((n & masks[m].second) << (1<<m));
        m++;
      }
      n >>= (1<<m) - bisectCount_;
      return n;
    }

    std::pair<size_t,size_t> patchIndex(T *ptr) const {
      size_t b = blockIndex(ptr);
      if (not bisectCount_) return {b, 0};
      size_t shift = (ptr - rawPtrs_[b]) % initialMemVecLen_;
      shift /= currentVecLen();
      return {b, reverseLastBisectCountBits(shift)};
    }

  public:
    // Constructor
    BisectMemManager() = delete;
    BisectMemManager(size_t memVecLen, size_t nVec):
        initialMemVecLen_(memVecLen-1), initialNVec_(nVec) {
      // Find the smallest power of 2 not less than memVecLen
      // Invalid for memVecLen > 2^32
      initialMemVecLen_ |= initialMemVecLen_ >> 1;
      initialMemVecLen_ |= initialMemVecLen_ >> 2;
      initialMemVecLen_ |= initialMemVecLen_ >> 4;
      initialMemVecLen_ |= initialMemVecLen_ >> 8;
      initialMemVecLen_ |= initialMemVecLen_ >> 16;
      ++initialMemVecLen_;
      memBlockSize_ = initialMemVecLen_ * initialNVec_;
    }

    size_t currentVecLen() const {
      return initialMemVecLen_ >> bisectCount_;
    }

    size_t getAllocCount() const {
      return allocCount_;
    }

    T* malloc() {

      if (not freePtrs_.empty()) {
        T *ptr = *freePtrs_.begin();
        freePtrs_.erase(freePtrs_.begin());
        return ptr;
      }

      allocCount_++;

      for (size_t I = 0; I < memPtrs_.size(); I++) {

        for (size_t J = 0; J < memPtrs_[I].size(); J++) {

          size_t &nAlloc = memPtrs_[I][J];
          if (nAlloc < initialNVec_) {
            T* ptr = rawPtrs_[I] + initialMemVecLen_ * (nAlloc++);
            if (bisectCount_) {
              ptr += currentVecLen() * reverseLastBisectCountBits(J);
            }
            return ptr;
          }

        }

        if (memPtrs_[I].size() < (1 << bisectCount_)) {
          memPtrs_[I].push_back({1});
          return rawPtrs_[I] + currentVecLen() * reverseLastBisectCountBits(memPtrs_[I].size() - 1);
        }

      }

      T* ptr = CQMemManager::get().calloc<T>(memBlockSize_);
      rawPtrs_.push_back(ptr);
      memPtrs_.push_back({1});

#ifdef CD_PROGRESS
      std::cout << "    Malloc a block of size " << memBlockSize_ << std::endl;
#endif

      return ptr;

    }

    void free(T* &ptr) {

      freePtrs_.insert(ptr);
      ptr = nullptr;

    }

    void shrink(std::vector<size_t> keeps) {

      if (keeps.size() == 0)
        return;

      for (size_t I = 0; I < memPtrs_.size(); I++) {

        for (size_t J = 0; J < memPtrs_[I].size(); J++) {

          size_t nAlloc = memPtrs_[I][J];

          if (not nAlloc)
            continue;

          T *pStart = rawPtrs_[I] + currentVecLen() * reverseLastBisectCountBits(J);
          T *pEnd = pStart + initialMemVecLen_ * nAlloc;

          #pragma omp parallel for
          for (T *ptr = pStart; ptr < pEnd; ptr += initialMemVecLen_) {

            size_t K = 0;
            for (size_t L : keeps)
              ptr[K++] = ptr[L];

          }

        }

      }

      while (keeps.size() <= currentVecLen() / 2) {
        bisectCount_++;

#ifdef CD_PROGRESS
        std::cout << "    Shrink by half." << std::endl;
#endif

      }

    }

    ~BisectMemManager() {
      for (T* p : rawPtrs_) {
        CQMemManager::get().free(p);
      }
    }

  }; // class BisectMemManager


  struct DynamicERIVec {
    std::pair<size_t,size_t> pq;
    double diag;
    double *vec;

    bool operator<(const DynamicERIVec &rhs) const {
      return diag < rhs.diag;
    }
  };


  struct DynamicShellPair {
    std::pair<size_t,size_t> PQ;
    double maxDiag;
    bool evaluated;
    size_t beginIndex;
    std::vector<DynamicERIVec> vecs;

    bool operator<(const DynamicShellPair &rhs) const {
      return maxDiag < rhs.maxDiag;
    }
  };


  /**
   *  \brief Copy existing element cases in computing the ERI vectors
   *         for dynamic-all CD algorithm.
   */
  inline bool computeDynamicAllERIVectorCopyCases(
      DynamicShellPair &candidatesI,
      std::vector<DynamicERIVec> &raw_vec_RS,
      std::vector<std::pair<size_t, double*>> &rsToEval,
      size_t NB, size_t RSstart,
      size_t R, size_t rBegin, size_t rSize,
      size_t S, size_t sBegin, size_t sSize,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      std::shared_ptr<InCore4indexTPI<double>> eri4I) {

    size_t PQstart = candidatesI.beginIndex;
    std::vector<DynamicERIVec> &raw_vec_PQ = candidatesI.vecs;

    if (candidatesI.evaluated) {
      // 4.2. For shell pair already evaluated, copy ERI elements

      // (rs|pq) = (pq|rs)
      for (size_t J = 0; J < raw_vec_PQ.size(); J++) {
        if (raw_vec_PQ[J].vec) {
          for (auto &K : rsToEval) {
            K.second[PQstart + J] = raw_vec_PQ[J].vec[RSstart + K.first];
          }
        }
      }

      return true;

    } else if (eri4I) {

      for (auto &K : rsToEval) {
        size_t r = raw_vec_RS[K.first].pq.first;
        size_t s = raw_vec_RS[K.first].pq.second;
        size_t rs = r + s * NB;

        for (size_t J = 0; J < raw_vec_PQ.size(); J++) {
          size_t p = raw_vec_PQ[J].pq.first;
          size_t q = raw_vec_PQ[J].pq.second;

          K.second[PQstart + J] = (*eri4I)(rs,p + q * NB);
#ifdef __DEBUGERI__
          std::cout << "ERI(" << r << "," << s << "|" << p << "," << q << "):"
                    << (*eri4I)(rs,p + q * NB) << std::endl;
#endif
        }
      }

      return true;

    } else if (PQstart == RSstart) {
      // 4.3. For diagonal blocks, copy ERI elements

      double *block = diagBlocks[std::make_pair(R,S)];
      for (auto &K : rsToEval) {
        size_t r = raw_vec_RS[K.first].pq.first - rBegin;
        size_t s = raw_vec_RS[K.first].pq.second - sBegin;
        size_t rs = (r + s * rSize) * rSize * sSize;

        for (auto &J : rsToEval) {
          size_t p = raw_vec_RS[J.first].pq.first - rBegin;
          size_t q = raw_vec_RS[J.first].pq.second - sBegin;

          K.second[PQstart + J.first] = block[rs + p + q * rSize];
#ifdef __DEBUGERI__
          std::cout << "ERI(" << r+rBegin << "," << s+sBegin
                    << "|" << p+rBegin << "," << q+sBegin << "):"
                    << std::setprecision(4)
                    << block[rs + p + q * rSize] << std::endl;
#endif
        }
      }

      return true;

    }

    return false;
  }


  /**
   *  \brief Compute the ERI vectors for dynamic-all CD algorithm with
   *         Libcint over the CGTO basis.
   */
  void computeDynamicAllERIVectorLibcint(
      BasisSet &basisSet,
      std::vector<DynamicShellPair> &candidates,
      DynamicShellPair &shellPair,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      BisectMemManager<double> &dynamicMem,
      std::vector<std::pair<size_t, double*>> &rsToEval,
      int *atm, int nAtoms, int *bas, int nShells,
      double *env, double *buffAll, int buffN4, double *cacheAll, int cache_size,
      double tau, size_t NB, std::shared_ptr<InCore4indexTPI<double>> eri4I,
      size_t &c1ERI, double &t1ERI) {

    size_t R = shellPair.PQ.first;
    size_t S = shellPair.PQ.second;

    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t rSize = basisSet.shells[R].size();
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();

    size_t RSstart = shellPair.beginIndex;

    std::vector<DynamicERIVec> &raw_vec_RS = shellPair.vecs;

    rsToEval.clear();

    // 4.1. Initialize raw vectors in shell
    for (size_t K = 0; K < raw_vec_RS.size(); K++) {
      if (raw_vec_RS[K].diag >= tau) {
        raw_vec_RS[K].vec = dynamicMem.malloc();
        rsToEval.push_back(std::make_pair(K, raw_vec_RS[K].vec));
      }
    }

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel for schedule(guided) reduction(+:lc1ERI,lt1ERI)
    for (size_t I = 0; I < candidates.size(); I++) {

      auto &candidatesI = candidates[I];

      if (candidatesI.maxDiag < tau)
        continue;

      size_t PQstart = candidatesI.beginIndex;

      std::vector<DynamicERIVec> &raw_vec_PQ = candidatesI.vecs;

      if (computeDynamicAllERIVectorCopyCases(
            candidatesI, raw_vec_RS, rsToEval, NB, RSstart,
            R, rBegin, rSize, S, sBegin, sSize, diagBlocks, eri4I)) {

      } else {
        // 4.4. For shell pair never evaluated, compute ERI elements

        size_t P = candidatesI.PQ.first;
        size_t Q = candidatesI.PQ.second;

        size_t pBegin = basisSet.mapSh2Bf[P];
        size_t pSize = basisSet.shells[P].size();
        size_t qBegin = basisSet.mapSh2Bf[Q];
        size_t qSize = basisSet.shells[Q].size();


        size_t thread_id = GetThreadID();

        double *buff = buffAll + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;

        int shls[4]{int(Q), int(P), int(S), int(R)};

        auto beginERI = tick();
        if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) {
          lt1ERI += tock(beginERI);
          lc1ERI++;
          continue;
        }
        lt1ERI += tock(beginERI);
        lc1ERI++;

        for (auto &K : rsToEval) {
          size_t r = raw_vec_RS[K.first].pq.first - rBegin;
          size_t s = raw_vec_RS[K.first].pq.second - sBegin;

#ifdef __DEBUGERI__
          std::cout << "(" << r+rBegin << "," << s+sBegin << "|:" << std::endl;
#endif

          for (size_t J = 0; J < raw_vec_PQ.size(); J++) {
            size_t p = raw_vec_PQ[J].pq.first - pBegin;
            size_t q = raw_vec_PQ[J].pq.second - qBegin;

            K.second[PQstart + J] =
                buff[((r * sSize + s) * pSize + p) * qSize + q];

#ifdef __DEBUGERI__
            std::cout << "  |" << p+pBegin << "," << q+qBegin << "):" << K.second[PQstart + J] << std::endl;
#endif
          }
        }


#ifdef __DEBUGERI__
        size_t pqSize = pSize * qSize;
        size_t rsSize = rSize * sSize;
        prettyPrintSmart(std::cout, "ERI", buff, pqSize, rsSize, pqSize);
#endif

      }
    }
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    // 4.5. Mark the shell pair as ERI evaluated.
    shellPair.evaluated = true;
  }; // computeDynamicAllERIVectorLibcint


  /**
   *  \brief Compute the ERI vectors for dynamic-all CD algorithm with
   *         Libint2 over the CGTO basis.
   */
  void computeDynamicAllERIVectorLibint(
      BasisSet &basisSet,
      std::vector<DynamicShellPair> &candidates,
      DynamicShellPair &shellPair,
      std::map<std::pair<size_t, size_t>, double*> &diagBlocks,
      BisectMemManager<double> &dynamicMem,
      std::vector<std::pair<size_t, double*>> &rsToEval,
      std::vector<libint2::Engine> &engines,
      std::vector<std::vector<libint2::Shell>> &shellPrims_,
      std::vector<double*> &coefBlocks_,
      std::vector<double*> &workBlocks,
      double tau, size_t NB, std::shared_ptr<InCore4indexTPI<double>> eri4I,
      size_t &c1ERI, double &t1ERI) {

    size_t R = shellPair.PQ.first;
    size_t S = shellPair.PQ.second;

    size_t rBegin = basisSet.mapSh2Bf[R];
    size_t rSize = basisSet.shells[R].size();
    size_t sBegin = basisSet.mapSh2Bf[S];
    size_t sSize = basisSet.shells[S].size();

    size_t RSstart = shellPair.beginIndex;

    std::vector<DynamicERIVec> &raw_vec_RS = shellPair.vecs;

    rsToEval.clear();

    // 4.1. Initialize raw vectors in shell
    for (size_t K = 0; K < raw_vec_RS.size(); K++) {
      if (raw_vec_RS[K].diag >= tau) {
        raw_vec_RS[K].vec = dynamicMem.malloc();
        rsToEval.push_back(std::make_pair(K, raw_vec_RS[K].vec));
      }
    }

    size_t lc1ERI = 0;
    double lt1ERI = 0.0;
    #pragma omp parallel for schedule(guided) reduction(+:lc1ERI,lt1ERI)
    for (size_t I = 0; I < candidates.size(); I++) {

      auto &candidatesI = candidates[I];

      if (candidatesI.maxDiag < tau)
        continue;

      size_t PQstart = candidatesI.beginIndex;

      std::vector<DynamicERIVec> &raw_vec_PQ = candidatesI.vecs;

      if (computeDynamicAllERIVectorCopyCases(
            candidatesI, raw_vec_RS, rsToEval, NB, RSstart,
            R, rBegin, rSize, S, sBegin, sSize, diagBlocks, eri4I)) {

      } else {
        // 4.4. For shell pair never evaluated, compute ERI elements

        size_t P = candidatesI.PQ.first;
        size_t Q = candidatesI.PQ.second;

        size_t pBegin = basisSet.mapSh2Bf[P];
        size_t pSize = basisSet.shells[P].size();
        size_t qBegin = basisSet.mapSh2Bf[Q];
        size_t qSize = basisSet.shells[Q].size();


        size_t thread_id = GetThreadID();
        const auto& buf_vec = engines[thread_id].results();

        if (basisSet.shells[P].ncontr() == 1 and basisSet.shells[Q].ncontr() == 1
            and basisSet.shells[R].ncontr() == 1 and basisSet.shells[S].ncontr() == 1) {

          // Evaluate ERI for shell quartet
          auto beginERI = tick();
          engines[thread_id].compute2<
            libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                basisSet.shells[R],
                basisSet.shells[S],
                basisSet.shells[P],
                basisSet.shells[Q]
          );
          lt1ERI += tock(beginERI);
          lc1ERI++;
          const auto *buff = buf_vec[0];
          if(buff == nullptr) continue;

          for (auto &K : rsToEval) {
            size_t r = raw_vec_RS[K.first].pq.first - rBegin;
            size_t s = raw_vec_RS[K.first].pq.second - sBegin;

#ifdef __DEBUGERI__
            std::cout << "(" << r+rBegin << "," << s+sBegin << "|:" << std::endl;
#endif

            for (size_t J = 0; J < raw_vec_PQ.size(); J++) {
              size_t p = raw_vec_PQ[J].pq.first - pBegin;
              size_t q = raw_vec_PQ[J].pq.second - qBegin;

              K.second[PQstart + J] =
                  buff[((r * sSize + s) * pSize + p) * qSize + q];

#ifdef __DEBUGERI__
              std::cout << "  |" << p+pBegin << "," << q+qBegin << "):" << K.second[PQstart + J] << std::endl;
#endif
            }
          }



#ifdef __DEBUGERI__
          size_t pqSize = pSize * qSize;
          size_t rsSize = rSize * sSize;
          prettyPrintSmart(std::cout, "ERI", buf_vec[0], pqSize, rsSize, pqSize);
#endif
        } else {

          size_t pContrSize = basisSet.shells[P].contr.size();
          size_t qContrSize = basisSet.shells[Q].contr.size();
          size_t rContrSize = basisSet.shells[R].contr.size();

          size_t pAMSize = shellPrims_[P][0].size();
          size_t qAMSize = shellPrims_[Q][0].size();
          size_t rAMSize = shellPrims_[R][0].size();
          size_t sAMSize = shellPrims_[S][0].size();

          size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

          const double *resPQRS;
          std::pair<size_t, double> counter_timer = libintGeneralContractionERI(
              P, Q, R, S,
              basisSet, engines[thread_id],
              shellPrims_, coefBlocks_, workBlocks[thread_id],
              resPQRS);
          lt1ERI += counter_timer.second;
          lc1ERI += counter_timer.first;

          for (auto &K : rsToEval) {
            size_t r = raw_vec_RS[K.first].pq.first - rBegin;
            size_t RR = r / rAMSize, rr = r % rAMSize;
            size_t s = raw_vec_RS[K.first].pq.second - sBegin;
            size_t SS = s / sAMSize, ss = s % sAMSize;

#ifdef __DEBUGERI__
            std::cout << "(" << r+rBegin << "," << s+sBegin << "|:" << std::endl;
#endif

            for (size_t J = 0; J < raw_vec_PQ.size(); J++) {
              size_t p = raw_vec_PQ[J].pq.first - pBegin;
              size_t PP = p / pAMSize, pp = p % pAMSize;
              size_t q = raw_vec_PQ[J].pq.second - qBegin;
              size_t QQ = q / qAMSize, qq = q % qAMSize;

              K.second[PQstart + J] = resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];

#ifdef __DEBUGERI__
              std::cout << "  |" << p+pBegin << "," << q+qBegin << "):" << K.second[PQstart + J] << std::endl;
#endif
            }
          }

        }
      }
    }
    c1ERI += lc1ERI;
    t1ERI += lt1ERI;

    // 4.5. Mark the shell pair as ERI evaluated.
    shellPair.evaluated = true;
  }; // computeDynamicAllERIVectorLibint


  /**
   *  \brief Compute the Cholesky RI pivots that agree with standard algorithm
   *         and Libint2 over the CGTO basis.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeCDPivots_DynamicAll(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeCDPivots_DynamicAll(BasisSet &basisSet) {

    auto topDynCDPivots = tick();

    std::cout << "Dynamic-All Pivots Determination:" << std::endl;
    std::cout << bannerMid << std::endl;

    size_t NB2 = NB*(NB+1)/2;


    // 1. Clear pivots_ and Cholesky vectors
    pivots_.clear();
    double* L;


    // 2. Initialize candidates

    // 2.1. Compute diagonal ERI elements
    // 2.2. Select all diagonals greater than the threshold
    double *diag = CQMemManager::get().calloc<double>(NB2);

    std::map<std::pair<size_t, size_t>, double*> diagBlocks;

    std::vector<size_t> D;

    if (libcint_)
      diagBlocks = computeDiagonalLibcint(basisSet, diag, true);
    else
      diagBlocks = computeDiagonalLibint(basisSet, diag, true);

    for (size_t pq = 0; pq < NB2; pq++) {
      if (diag[pq] >= tau_)
        D.push_back(compoundToSquare(pq, NB));
    }

    if (D.empty())
      CErr("Cholesky decomposition threshold is greater than all TPI diagonal elements.");

#ifdef __DEBUGERI__
    std::cout << "Raw Diag: [ " << std::endl;
    for (size_t Q = 0; Q < NB; Q++)
      for (size_t P = 0; P <= Q; P++)
        std::cout << "(" << P << "," << Q
                  << "):\t" << diag[P + Q*NB] << ", " << std::endl;
    std::cout << "]" << std::endl;
#endif

    // 2.3. Group basis function pairs by shell
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        DindicesByShell = groupPivotsByShell(basisSet, D);


    // 2.4. Initialize candidate container
    std::vector<DynamicShellPair> candidates(DindicesByShell.size());

    size_t candidateSize = 0;
    size_t I = 0;
    for (auto &kv : DindicesByShell) {

      candidates[I].PQ.first = kv.first.first;
      candidates[I].PQ.second = kv.first.second;

      std::vector<DynamicERIVec> &raw_vec = candidates[I].vecs;

      double maxDiag = 0.0;
      for (size_t pInd : kv.second) {
        size_t pq = D[pInd];
        maxDiag = std::max(maxDiag, diag[squareToCompound(pq, NB)]);
        raw_vec.push_back({{pq % NB, pq / NB}, diag[squareToCompound(pq, NB)], nullptr});
      }

      candidates[I].maxDiag = maxDiag;
      candidates[I].evaluated = false;
      candidates[I].beginIndex = candidateSize;


      I++;
      candidateSize += raw_vec.size();
    }

#ifdef CD_PROGRESS
      std::cout << "  Initial candidate size = " << candidateSize << std::endl;
#endif


    // 2.5. Free memory
    CQMemManager::get().free(diag);
    DindicesByShell.clear();


    // Preparation for Cholesky loop

    std::vector<std::pair<size_t, double*>> rsToEval;
    std::vector<std::vector<double*>> ptrsToFreeThreads(GetNumThreads());

    // Memory control: avoid frequent calling of malloc and free.
    BisectMemManager<double> dynamicMem(candidateSize, NB);

    std::vector<size_t> preserve;
    preserve.reserve(candidateSize);

    size_t Lalloc = std::min(candidateSize, maxQual_);
    L = CQMemManager::get().calloc<double>(candidateSize * Lalloc);

    auto beginERIvec = topDynCDPivots;
    auto beginCDalg = topDynCDPivots;
    auto beginShrink = topDynCDPivots;
    double cumeERIvec = 0, cumeCDalg = 0, cumeShrink = 0;

    size_t shrinkCost = 0, extraCost = 0, preShrinkIter = 0;
    size_t preCandSize = candidateSize, curCandSize = 0;
    size_t shrinkCount = 0;

    while (true) {

      // 3. Locate a pivot shell
      auto shellPairIter = std::max_element(candidates.begin(), candidates.end());

      if (shellPairIter->maxDiag < tau_)
        break;

      std::vector<DynamicERIVec> &raw_vec_RS = shellPairIter->vecs;

      // 4. Evaluate ERI elements
      if (not shellPairIter->evaluated) {

        beginERIvec = tick();

        if (libcint_)
          computeDynamicAllERIVectorLibcint(basisSet, candidates, *shellPairIter,
              diagBlocks, dynamicMem, rsToEval,
              atm, nAtoms, bas, nShells, env, buffAll, buffN4, cacheAll, cache_size,
              tau_, NB, eri4I_, c1ERI, t1ERI);
        else
          computeDynamicAllERIVectorLibint(basisSet, candidates, *shellPairIter,
              diagBlocks, dynamicMem, rsToEval,
              engines, shellPrims_, coefBlocks_, workBlocks,
              tau_, NB, eri4I_, c1ERI, t1ERI);

        cumeERIvec += tock(beginERIvec);

      }

      // 5. Compute the CD vector
      auto rsIter = std::max_element(raw_vec_RS.begin(), raw_vec_RS.end());

#ifdef CD_PROGRESS
      std::cout << "    Selected: (" << rsIter->pq.first << "," << rsIter->pq.second
                << ")\tDiag: " << rsIter->diag << std::endl;
#endif

      size_t rsIndex = shellPairIter->beginIndex + std::distance(raw_vec_RS.begin(), rsIter);

      double *Lq = L + pivots_.size() * candidateSize;

      beginCDalg = tick();
      std::copy_n(rsIter->vec, candidateSize, Lq);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::Trans, candidateSize, 1, pivots_.size(),
           -1.0, L, candidateSize,
           L + rsIndex, candidateSize,
           1.0, Lq, candidateSize);
      blas::scal(candidateSize, 1.0 / sqrt(rsIter->diag), Lq, 1);
      cumeCDalg += tock(beginCDalg);

      rsIter->diag = 0.0;
      dynamicMem.free(rsIter->vec);
      pivots_.push_back(rsIter->pq.first + rsIter->pq.second * NB);

      // 6. Update diagonal elements
      curCandSize = 0;
      #pragma omp parallel for schedule(guided) reduction(+:curCandSize)
      for (size_t I = 0; I < candidates.size(); I++) {

        auto &candidatesI = candidates[I];

        if (candidatesI.maxDiag < tau_)
          continue;

        size_t RSstart = candidatesI.beginIndex;
        std::vector<DynamicERIVec> &raw_vec = candidatesI.vecs;

        size_t thread_id = GetThreadID();

        double maxDiag = 0.0;
        for (size_t J = 0; J < raw_vec.size(); J++) {
          double &diagJ = raw_vec[J].diag;
          if (diagJ >= tau_) {
            double dDiag = Lq[RSstart + J];
            diagJ -= dDiag * dDiag;
          }
          if (diagJ >= tau_) {
            curCandSize++;
          } else if (raw_vec[J].vec) {
            ptrsToFreeThreads[thread_id].push_back(raw_vec[J].vec);
            raw_vec[J].vec = nullptr;
          }
          maxDiag = std::max(maxDiag, diagJ);
        }
        candidatesI.maxDiag = maxDiag;
      }

      for (auto &ptrsToFree : ptrsToFreeThreads) {
        for (double *p : ptrsToFree)
          dynamicMem.free(p);
        ptrsToFree.clear();
      }

#ifdef __DEBUGERI__
      std::cout << "Diag: [ " << std::endl;
      for (auto &kv : candidates)
        for (auto &raw : std::get<5>(kv))
          if (std::get<2>(raw) >= tau_)
            std::cout << "(" << std::get<0>(raw) << "," << std::get<1>(raw)
                      << "):\t" << std::get<2>(raw) << ", " << std::endl;
      std::cout << "]" << std::endl;
#endif


      // 7. Shrink Cholesky and ERI vectors
      // 7.1. Evaluate shrink costs
      shrinkCost = (dynamicMem.getAllocCount() + pivots_.size()) * curCandSize;
      extraCost += (pivots_.size() - 1) * (preCandSize - curCandSize);

      if (pivots_.size() >= Lalloc
          or extraCost >= shrinkCost and (pivots_.size() - preShrinkIter) >= minShrinkCycle_) {

        beginShrink = tick();

        preserve.clear();
        size_t I = 0;
        for (auto &kv : candidates) {
          auto &raw_vec = kv.vecs;
          if (kv.maxDiag < tau_)
            I += raw_vec.size();
          else
            for (auto &rs : raw_vec) {
              if (rs.diag >= tau_) {
                preserve.push_back(I);
              }
              I++;
            }
        }


        if (preserve.size() == 0)
          break;


        candidates.erase(
            std::remove_if(candidates.begin(), candidates.end(),
                           [this](DynamicShellPair &a){
                             return a.maxDiag < tau_;
                           }),
            candidates.end());

        #pragma omp parallel for
        for (size_t I = 0; I < candidates.size(); I++) {
          auto &raw_vec = candidates[I].vecs;
          raw_vec.erase(
              std::remove_if(raw_vec.begin(), raw_vec.end(),
                             [this](DynamicERIVec &a){
                               return a.diag < tau_;
                             }),
              raw_vec.end());
        }

        candidateSize = 0;
        for (auto &kv : candidates) {
          kv.beginIndex = candidateSize;
          candidateSize += kv.vecs.size();
        }

        dynamicMem.shrink(preserve);

        Lalloc = pivots_.size() + std::min(candidateSize, maxQual_);
        double *Lnew = CQMemManager::get().calloc<double>(candidateSize * Lalloc);

        #pragma omp parallel for
        for (size_t p = 0; p < pivots_.size(); p++) {
          double *Lp = L + p * preCandSize, *Lnp = Lnew + p * candidateSize;
          size_t j = 0;
          for (size_t i : preserve) {
            Lnp[j++] = Lp[i];
          }
        }

        CQMemManager::get().free(L);

        L = Lnew;

        shrinkCount++;
        preShrinkIter = pivots_.size();
        preCandSize = curCandSize;
        extraCost = 0;

        cumeShrink += tock(beginShrink);

      }

    }

    CQMemManager::get().free(L);
    for (auto &kv : diagBlocks) {
      CQMemManager::get().free(kv.second);
    }


#ifdef __DEBUGERI__
    std::cout << "pivots: [ ";
    for (size_t p : pivots_)
      std::cout << p << ", ";
    std::cout << "]" << std::endl;
    std::vector<size_t> sortedPivots(pivots_);
    std::sort(sortedPivots.begin(), sortedPivots.end());
    std::cout << "sorted pivots: [ ";
    for (size_t p : sortedPivots)
      std::cout << p << ", ";
    std::cout << "]" << std::endl;
#endif

    auto durDynCDPivots = tock(topDynCDPivots);
    std::cout << "  Cholesky-RI-Pivots-ERI count    = " << c1ERI << std::endl;
    std::cout << "  Cholesky-RI-Pivots-ERI duration = " << t1ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-ERIvec duration     = " << cumeERIvec << " s " << std::endl;
    std::cout << "  Cholesky-RI-CDalg duration      = " << cumeCDalg << " s " << std::endl;
    std::cout << "  Cholesky-RI-Shrink count        = " << shrinkCount << std::endl;
    std::cout << "  Cholesky-RI-Shrink duration     = " << cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-misc duration       = "
              << durDynCDPivots - cumeERIvec - cumeCDalg - cumeShrink << " s " << std::endl;
    std::cout << "  Cholesky-RI-Dynamic-All-Pivots duration = " << durDynCDPivots << " s " << std::endl;

  }; // InCoreCholeskyRIERI<double>::computeCDPivots_DynamicAll


  /**
   *  \brief Group the pivot indices by Shell pairs
   */
  template <typename IntsT>
  std::map<std::pair<size_t,size_t>, std::vector<size_t>>
  InCoreCholeskyRIERI<IntsT>::groupPivotsByShell(
      BasisSet &basisSet, const std::vector<size_t> &pivots) {

    std::map<std::pair<size_t,size_t>, std::vector<size_t>> pivotsByShell;

    for (size_t i = 0; i < pivots.size(); i++) {

      size_t pivot = pivots[i];

      size_t r = pivot % basisSet.nBasis;
      size_t s = pivot / basisSet.nBasis;

      if ( r > s ) std::swap(r, s);

      size_t R = std::distance(basisSet.mapSh2Bf.begin(),
          std::upper_bound(basisSet.mapSh2Bf.begin(),
              basisSet.mapSh2Bf.end(), r)) - 1;
      size_t S = std::distance(basisSet.mapSh2Bf.begin(),
          std::upper_bound(basisSet.mapSh2Bf.begin(),
              basisSet.mapSh2Bf.end(), s)) - 1;

      pivotsByShell[std::make_pair(R, S)].push_back(i);

    }

    return pivotsByShell;

  }; // InCoreCholeskyRIERI<IntsT>::groupPivotsByShell


  /**
   *  \brief Compute the Cholesky RI 3-index ERI tensor
   *         using Libcint over the CGTO basis for pivot RI algorithm.
   *         Pivots are assumed to be already computed and stored
   *         in pivots_.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computePivotRI3indexERILibcint(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computePivotRI3indexERILibcint(BasisSet &basisSet) {

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        pivotIndicesByShell = groupPivotsByShell(basisSet, pivots_);

    size_t pivotShellSize = pivotIndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivotShells;
    pivotShells.reserve(pivotShellSize);
    for (auto &shell_pivot : pivotIndicesByShell) {
      pivotShells.push_back(shell_pivot.first);
    }

    size_t lc2ERI = 0;
    double lt2ERI = 0.0;
    #pragma omp parallel reduction(+:lc2ERI,lt2ERI)
    {
      size_t thread_id = GetThreadID();

      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll + cache_size*thread_id;

      for (size_t I = 0, IPQ = 0; I < pivotShellSize; I++) {

        const auto &RSpair = pivotShells[I];
        auto &shell_pivot = pivotIndicesByShell[RSpair];

        size_t R = RSpair.first;
        size_t S = RSpair.second;

        size_t rBegin = basisSet.mapSh2Bf[R];
        size_t sBegin = basisSet.mapSh2Bf[S];
        size_t rSize = basisSet.shells[R].size();
        size_t sSize = basisSet.shells[S].size();
        size_t rEnd = rBegin + rSize;
        size_t sEnd = sBegin + sSize;

        int shls[4]{int(0), int(0), int(S), int(R)};

        for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
          for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++, IPQ++) {

            // Round Robbin work distribution
            #ifdef _OPENMP
            if( IPQ % nthreads != thread_id ) continue;
            #endif

            auto PQpair = std::make_pair(P,Q);
            bool hasPivot = pivotIndicesByShell.find(PQpair) != pivotIndicesByShell.end();
            if (hasPivot and PQpair < RSpair)
              continue;

            shls[0] = int(Q);
            shls[1] = int(P);

            auto beginERI = tick();
            if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) {
              lt2ERI += tock(beginERI);
              lc2ERI++;
              continue;
            }
            lt2ERI += tock(beginERI);
            lc2ERI++;

            for (auto &pivot_index : shell_pivot) {

              auto rs = anaSquare(pivots_[pivot_index], basisSet.nBasis);
              size_t r = rs.first;
              size_t s = rs.second;

              if (P == Q) {
                for (size_t pBegin(basisSet.mapSh2Bf[P]),
                     pSize(basisSet.shells[P].size()),
                     pEnd(pBegin + pSize),
                     p(pBegin),
                     rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                     p < pEnd; p++, rsp++) {
                  for (size_t q(p),
                       rspq = rsp * pSize + p - pBegin;
                       q < pEnd; q++, rspq++) {
                    this->pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                  }
                }
              } else {
                for (size_t p(basisSet.mapSh2Bf[P]),
                     pSize(basisSet.shells[P].size()),
                     pEnd(p + pSize),
                     rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                     p < pEnd; p++, rsp++) {
                  for (size_t q(basisSet.mapSh2Bf[Q]),
                       qSize(basisSet.shells[Q].size()),
                       qEnd(q + qSize),
                       rspq = rsp * qSize;
                       q < qEnd; q++, rspq++) {
                    this->pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                  }
                }
              }
            }

            if (hasPivot and PQpair != RSpair) {

              size_t pBegin = basisSet.mapSh2Bf[P];
              size_t qBegin = basisSet.mapSh2Bf[Q];
              size_t pSize = basisSet.shells[P].size();
              size_t qSize = basisSet.shells[Q].size();

              const std::vector<size_t> &PQpivotIndicesInShell = pivotIndicesByShell[PQpair];

              for (auto &pivot_index : PQpivotIndicesInShell) {

                auto pq = anaSquare(pivots_[pivot_index], basisSet.nBasis);
                size_t p = pq.first;
                size_t q = pq.second;

                for (size_t r(rBegin); r < rEnd; r++) {
                  for (size_t s(R == S ? r : sBegin);
                       s < sEnd; s++) {
                    this->pointer()[pivot_index + toCompound(r,s) * NBRI] =
                        buff[(((r-rBegin) * sSize + (s-sBegin)) * pSize + (p - pBegin)) * qSize + (q - qBegin)];
                  }
                }
              }

            }

          }; // Q
        }; // P

      }
    }; // omp region
    c2ERI += lc2ERI;
    t2ERI += lt2ERI;

  }; // InCoreCholeskyRIERI<double>::computePivotRI3indexERILibcint


  /**
   *  \brief Compute the Cholesky RI 3-index ERI tensor
   *         using Libint over the CGTO basis for pivot RI algorithm.
   *         Pivots are assumed to be already computed and stored
   *         in pivots_.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computePivotRI3indexERILibint(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computePivotRI3indexERILibint(BasisSet &basisSet) {

    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();

    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        pivotIndicesByShell = groupPivotsByShell(basisSet, pivots_);

    size_t pivotShellSize = pivotIndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivotShells;
    pivotShells.reserve(pivotShellSize);
    for (auto &shell_pivot : pivotIndicesByShell) {
      pivotShells.push_back(shell_pivot.first);
    }

    size_t lc2ERI = 0;
    double lt2ERI = 0.0;
    #pragma omp parallel reduction(+:lc2ERI,lt2ERI)
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      for (size_t I = 0, IPQ = 0; I < pivotShellSize; I++) {

        const auto &RSpair = pivotShells[I];
        auto &shell_pivot = pivotIndicesByShell[RSpair];

        size_t R = RSpair.first;
        size_t S = RSpair.second;

        size_t rBegin = basisSet.mapSh2Bf[R];
        size_t sBegin = basisSet.mapSh2Bf[S];
        size_t rSize = basisSet.shells[R].size();
        size_t sSize = basisSet.shells[S].size();
        size_t rEnd = rBegin + rSize;
        size_t sEnd = sBegin + sSize;

        for (size_t P(0), PQ(0); P < basisSet.nShell; P++) {
          for (size_t Q = P; Q < basisSet.nShell; Q++, PQ++, IPQ++) {

            // Round Robbin work distribution
            #ifdef _OPENMP
            if( IPQ % nthreads != thread_id ) continue;
            #endif

            auto PQpair = std::make_pair(P,Q);
            bool hasPivot = pivotIndicesByShell.find(PQpair) != pivotIndicesByShell.end();
            if (hasPivot and PQpair < RSpair)
              continue;

            if (basisSet.shells[P].ncontr() == 1 and basisSet.shells[Q].ncontr() == 1
                and basisSet.shells[R].ncontr() == 1 and basisSet.shells[S].ncontr() == 1) {
              // Evaluate ERI for shell quartet
              auto beginERI = tick();
              engines[thread_id].compute2<
                libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                    basisSet.shells[R],
                    basisSet.shells[S],
                    basisSet.shells[P],
                    basisSet.shells[Q]
              );
              lt2ERI += tock(beginERI);
              lc2ERI++;
              const auto *buff =  buf_vec[0] ;
              if(buff == nullptr) continue;

              for (auto &pivot_index : shell_pivot) {

                auto rs = anaSquare(pivots_[pivot_index], basisSet.nBasis);
                size_t r = rs.first;
                size_t s = rs.second;

                if (P == Q) {
                  for (size_t pBegin(basisSet.mapSh2Bf[P]),
                       pSize(basisSet.shells[P].size()),
                       pEnd(pBegin + pSize),
                       p(pBegin),
                       rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                       p < pEnd; p++, rsp++) {
                    for (size_t q(p),
                         rspq = rsp * pSize + p - pBegin;
                         q < pEnd; q++, rspq++) {
                      this->pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                    }
                  }
                } else {
                  for (size_t p(basisSet.mapSh2Bf[P]),
                       pSize(basisSet.shells[P].size()),
                       pEnd(p + pSize),
                       rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                       p < pEnd; p++, rsp++) {
                    for (size_t q(basisSet.mapSh2Bf[Q]),
                         qSize(basisSet.shells[Q].size()),
                         qEnd(q + qSize),
                         rspq = rsp * qSize;
                         q < qEnd; q++, rspq++) {
                      this->pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                    }
                  }
                }
              }

              if (hasPivot and PQpair != RSpair) {

                size_t pBegin = basisSet.mapSh2Bf[P];
                size_t qBegin = basisSet.mapSh2Bf[Q];
                size_t pSize = basisSet.shells[P].size();
                size_t qSize = basisSet.shells[Q].size();

                const std::vector<size_t> &PQpivotIndicesInShell = pivotIndicesByShell[PQpair];

                for (auto &pivot_index : PQpivotIndicesInShell) {

                  auto pq = anaSquare(pivots_[pivot_index], basisSet.nBasis);
                  size_t p = pq.first;
                  size_t q = pq.second;

                  for (size_t r(rBegin); r < rEnd; r++) {
                    for (size_t s(R == S ? r : sBegin);
                         s < sEnd; s++) {
                      this->pointer()[pivot_index + toCompound(r,s) * NBRI] =
                          buff[(((r-rBegin) * sSize + (s-sBegin)) * pSize + (p - pBegin)) * qSize + (q - qBegin)];
                    }
                  }
                }

              }

            } else {

              size_t pContrSize = basisSet.shells[P].contr.size();
              size_t qContrSize = basisSet.shells[Q].contr.size();
              size_t rContrSize = basisSet.shells[R].contr.size();

              size_t pAMSize = shellPrims_[P][0].size();
              size_t qAMSize = shellPrims_[Q][0].size();
              size_t rAMSize = shellPrims_[R][0].size();
              size_t sAMSize = shellPrims_[S][0].size();

              size_t pqrsAMSize = pAMSize * qAMSize * rAMSize * sAMSize;

              const double *resPQRS;
              std::pair<size_t, double> counter_timer = libintGeneralContractionERI(
                  P, Q, R, S,
                  basisSet, engines[thread_id],
                  shellPrims_, coefBlocks_, workBlocks[thread_id],
                  resPQRS);
              lt2ERI += counter_timer.second;
              lc2ERI += counter_timer.first;

              for (auto &pivot_index : shell_pivot) {

                auto rs = anaSquare(pivots_[pivot_index], basisSet.nBasis);

                size_t r = rs.first - rBegin;
                size_t RR = r / rAMSize, rr = r % rAMSize;

                size_t s = rs.second - sBegin;
                size_t SS = s / sAMSize, ss = s % sAMSize;

                if (P == Q) {
                  for (size_t q(0),
                       qBegin(basisSet.mapSh2Bf[Q]),
                       qSize(basisSet.shells[Q].size());
                       q < qSize; q++) {

                    size_t QQ = q / qAMSize, qq = q % qAMSize;

                    for (size_t p(0); p <= q; p++) {

                      size_t PP = p / pAMSize, pp = p % pAMSize;

                      this->pointer()[pivot_index + toCompound(p + qBegin, q + qBegin) * NBRI] =
                          resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                    }
                  }

                } else {
                  for (size_t q(0),
                       qBegin(basisSet.mapSh2Bf[Q]),
                       qSize(basisSet.shells[Q].size());
                       q < qSize; q++) {

                    size_t QQ = q / qAMSize, qq = q % qAMSize;

                    for (size_t p(0),
                         pBegin(basisSet.mapSh2Bf[P]),
                         pSize(basisSet.shells[P].size());
                         p < pSize; p++) {

                      size_t PP = p / pAMSize, pp = p % pAMSize;

                      this->pointer()[pivot_index + toCompound(p + pBegin, q + qBegin) * NBRI] =
                          resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                    }
                  }
                }
              }

              if (hasPivot and PQpair != RSpair) {

                size_t pBegin = basisSet.mapSh2Bf[P];
                size_t qBegin = basisSet.mapSh2Bf[Q];

                const std::vector<size_t> &PQpivotIndicesInShell = pivotIndicesByShell[PQpair];

                for (auto &pivot_index : PQpivotIndicesInShell) {

                  auto pq = anaSquare(pivots_[pivot_index], basisSet.nBasis);

                  size_t p = pq.first - pBegin;
                  size_t PP = p / pAMSize, pp = p % pAMSize;

                  size_t q = pq.second - qBegin;
                  size_t QQ = q / qAMSize, qq = q % qAMSize;

                  for (size_t s(0); s < sSize; s++) {

                    size_t SS = s / sAMSize, ss = s % sAMSize;

                    for (size_t r(0), rEnd(R == S ? s : rSize - 1);
                         r <= rEnd; r++) {

                      size_t RR = r / rAMSize, rr = r % rAMSize;

                      this->pointer()[pivot_index + toCompound(r + rBegin, s + sBegin) * NBRI] =
                          resPQRS[(ss + sAMSize * (rr + rAMSize * (qq + qAMSize * pp)))
                                  + pqrsAMSize * (PP + pContrSize * (QQ + qContrSize * (RR + rContrSize * SS)))];
                    }
                  }
                }

              }

            }

          }; // Q
        }; // P

      }
    }; // omp region
    c2ERI += lc2ERI;
    t2ERI += lt2ERI;

  }; // InCoreCholeskyRIERI<double>::computePivotRI3indexERILibint


  template <>
  void InCoreCholeskyRIERI<dcomplex>::extractTwoCenterSubsetFrom3indexERI(
      const std::vector<size_t> &pivots, size_t NBRI, size_t NB,
      const double* eri3J, size_t LD3J, double* S, size_t LDS, bool upperTriOnly) {
    CErr("Only real GTOs are allowed",std::cout);
  }
  template <>
  void InCoreCholeskyRIERI<double>::extractTwoCenterSubsetFrom3indexERI(
      const std::vector<size_t> &pivots, size_t NBRI, size_t NB,
      const double* eri3J, size_t LD3J, double* S, size_t LDS, bool upperTriOnly) {

    auto topLibintPivot2Index = tick();

    // Only build Upper triangular part of S for Cholesky decomposition
    size_t qMax = pivots.size();
    #pragma omp parallel for
    for (size_t Q = 0; Q < qMax; Q++) {
      const double *ptr = eri3J + squareToCompound(pivots[Q],NB) * LD3J;
      for (size_t P = upperTriOnly ? Q : 0; P < NBRI; P++) {

        S[Q + P*LDS] = ptr[P];

      }
    }

    #ifdef __DEBUGERI__
    prettyPrintSmart(std::cout, "S", S, NBRI, NBRI, LDS);
    #endif

    auto durLibintPivot2Index = tock(topLibintPivot2Index);
    std::cout << "  Cholesky-RI-PivotRI-2index duration = " << durLibintPivot2Index << " s " << std::endl;
  }; // InCoreCholeskyRIERI<double>::extractTwoCenterSubsetFrom3indexERI


  /**
   *  \brief Allocate, compute and store the Cholesky RI
   *         3-index ERI tensor using Libint2 over the CGTO basis.
   *         Using shell-by-shell pivot RI algorithm.
   *         Pivots are assumed to be already computed and stored
   *         in pivots_.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computePivotRI(BasisSet&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computePivotRI(BasisSet &basisSet) {

    std::cout << std::endl << "Build 3-index RIERI tensor:" << std::endl;
    std::cout << bannerMid << std::endl;

    size_t NB2 = NB*NB, NBC = NB*(NB+1)/2;

    auto topLibintPivotRI = tick();

    setNRIBasis(pivots_.size());
    clear();

    if (eri4I_) {

      #pragma omp parallel for
      for (size_t P = 0; P < NBRI; P++) {
        size_t pivot_index = pivots_[P];
        double *Lrs = eri4I_->pointer() + pivot_index * NB2;

        for (size_t ij = 0; ij < NBC; ij++) {
          this->pointer()[P + ij * NBRI] = Lrs[compoundToSquare(ij, NB)];
        }
      }

    } else if (libcint_) {

      computePivotRI3indexERILibcint(basisSet);

    } else {

      computePivotRI3indexERILibint(basisSet);

    }; // if (eri4I_)

#ifdef __DEBUGERI__
    std::cout << "Ppq:" << std::endl;
    std::cout << "[";
    for (size_t Q = 0; Q < NBRI; Q++) {
      size_t r = pivots_[Q] % basisSet.nBasis;
      size_t s = pivots_[Q] / basisSet.nBasis;
      std::cout << "(" << r << "," << s << "):[";
      for (size_t p = 0; p < NB; p++) {
        std::cout << "[";
        for (size_t q = 0; q < NB; q++) {
          if (p < q)
            std::cout << this->pointer()[Q + toCompound(p,q) * NBRI] << "\t";
          else
            std::cout << this->pointer()[Q + toCompound(q,p) * NBRI] << "\t";
        }
        std::cout << "]" << std::endl;
      }
//      for (size_t ij = 0; ij < NBC; ij++)
//        std::cout << this->pointer()[Q + ij * NBRI] << "\t";
      std::cout << "]" << std::endl;
    }
    std::cout << "]" << std::endl;
#endif


    auto durLibintPivot3Index = tock(topLibintPivotRI);
    std::cout << "  Cholesky-RI-PivotRI-ERI count       = " << c2ERI << std::endl;
    std::cout << "  Cholesky-RI-PivotRI-ERI duration    = " << t2ERI << " s " << std::endl;
    std::cout << "  Cholesky-RI-PivotRI-3index duration = " << durLibintPivot3Index << " s " << std::endl;

    saveRawERI3J();

    twocenterERI_ = std::make_shared<cqmatrix::Matrix<double>>(NBRI);
    double *S = twocenterERI_->pointer();
    extractTwoCenterSubsetFrom3indexERI(pivots_, NBRI, NB, pointer(), NBRI, S, NBRI);
    if (saveRawERI_)
      rawERI2C_ = std::make_shared<cqmatrix::Matrix<double>>(*twocenterERI_);

#ifdef __DEBUGERI__
    std::cout << "Pivots:" << std::endl;
    for (size_t pivot : pivots_)
      std::cout << pivot << std::endl;
#endif

    auto topERI3Trans = tick();
    halfInverse2CenterERI(*twocenterERI_);
    contract2CenterERI();

    auto durERI3Trans = tock(topERI3Trans);
    std::cout << "  RI-ERI3-Transformation duration = " << durERI3Trans << " s " << std::endl;

    auto durLibintPivotRI = tock(topLibintPivotRI);
    std::cout << "  Cholesky-RI-PivotRI duration = " << durLibintPivotRI << " s " << std::endl;

  }; // InCoreCholeskyRIERI<double>::computePivotRI


  /**
   *  \brief Allocate, compute and store the Cholesky RI
   *  3-index ERI tensor using Libint2 over the CGTO basis.
   */
  template <>
  void InCoreCholeskyRIERI<dcomplex>::computeAOInts(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void InCoreCholeskyRIERI<double>::computeAOInts(BasisSet &basisSet, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {


    std::cout << std::left;
    std::cout << std::endl << "ERI Cholesky Decomposition";
    std::cout << ":" << std::endl << BannerTop << std::endl << std::endl;

    if (op != ELECTRON_REPULSION)
      CErr("Only Electron repulsion integrals in InCoreCholeskyRIERI<double>",std::cout);
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in InCoreCholeskyRIERI<double>",std::cout);

    if (options.Libcint) {
      libcint_ = true;
      if (not generalContraction_) {
        generalContraction_ = true;
        std::cout << "Warning: Forced general contraction algorithm for libcint." << std::endl << std::endl;
      }
    }

    std::cout << "Parameters and options:" << std::endl;
    std::cout << bannerMid << std::endl;
    const int fieldNameWidth(40);
    std::cout << "  " << std::setw(fieldNameWidth) << "Algorithm:";
    switch (alg_) {
    case CHOLESKY_ALG::TRADITIONAL:
      std::cout << "Traditional";
      break;
    case CHOLESKY_ALG::DYNAMIC_ALL:
      std::cout << "Dynamic-All";
      break;
    case CHOLESKY_ALG::SPAN_FACTOR:
      std::cout << "Span-Factor";
      break;
    case CHOLESKY_ALG::DYNAMIC_ERI:
      std::cout << "Dynamic-ERI";
      break;
    case CHOLESKY_ALG::SPAN_FACTOR_REUSE:
      std::cout << "Span-Factor-Reuse";
      break;
    }
    std::cout << std::endl;
    std::cout << "  " << std::setw(fieldNameWidth) << "Threshold:"
        << tau_ << std::endl;
    std::cout << "  " << std::setw(fieldNameWidth) << "ERI library:"
        << (libcint_ ? "Libcint" : "Libint2") << std::endl;
    std::cout << "  " << std::setw(fieldNameWidth) << "General contraction:"
        << (generalContraction_ ? "True" : "False") << std::endl;
    std::cout << "  " << std::setw(fieldNameWidth) << "Have already computed 4-index ERI:"
        << (eri4I_ ? "True" : "False") << std::endl;
    std::cout << "  " << std::setw(fieldNameWidth) << "Build 4-index ERI:"
        << (build4I_ ? "True" : "False") << std::endl;
    switch (alg_) {
    case CHOLESKY_ALG::DYNAMIC_ALL:
      std::cout << "  " << std::setw(fieldNameWidth) << "Min shrink cycle:"
          << minShrinkCycle_ << std::endl;
      break;
    case CHOLESKY_ALG::SPAN_FACTOR:
    case CHOLESKY_ALG::DYNAMIC_ERI:
    case CHOLESKY_ALG::SPAN_FACTOR_REUSE:
      std::cout << "  " << std::setw(fieldNameWidth) << "Sigma:"
          << sigma_ << std::endl;
      std::cout << "  " << std::setw(fieldNameWidth) << "Max qualification:"
          << maxQual_ << std::endl;
      break;
    default:
      break;
    }
    std::cout << std::endl;


    if (build4I_ and not eri4I_) {
      std::cout << "Building 4-index ERI:" << std::endl;
      std::cout << bannerMid << std::endl;
      auto top4I = tick();

      eri4I_ = std::make_shared<InCore4indexTPI<double>>(NB);
      eri4I_->computeAOInts(basisSet, mol, emPert, op, options);

      auto dur4I = tock(top4I);
      std::cout << "  Cholesky-4-Index duration   = " << dur4I << " s " << std::endl << std::endl;
    }

    auto topCholeskyRI = tick();

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();


    if (not libcint_) {
      engines.resize(nthreads);

      // Initialize the first engine for the integral evaluation
      engines[0] = libint2::Engine(libint2::Operator::coulomb,
        basisSet.maxPrim, basisSet.maxL,0);
      engines[0].set_precision(0.);

      // Copy over the engines to other threads if need be
      for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];
    }

    BasisSet groupedBasisSet(basisSet);
    if (generalContraction_) {

      groupedBasisSet = basisSet.groupGeneralContractionBasis();

      maxNcontrAMSize_ = std::max_element(groupedBasisSet.shells.begin(),
                                          groupedBasisSet.shells.end(),
                                          [](libint2::Shell &a, libint2::Shell &b) {
                                            return a.size() < b.size();
                                          })->size();

      if (options.Libcint) {
        if (basisSet.forceCart)
          CErr("Libcint + cartesian GTO NYI.");

        nAtoms = mol.nAtoms;
        nShells = groupedBasisSet.nShell;

        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
        bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
        env = CQMemManager::get().template malloc<double>(groupedBasisSet.getLibcintEnvLength(mol));

        groupedBasisSet.setLibcintEnv(mol, atm, bas, env);

        // Get threads result buffer
        buffN4 = maxNcontrAMSize_*maxNcontrAMSize_*maxNcontrAMSize_*maxNcontrAMSize_;
        buffAll = CQMemManager::get().malloc<double>(buffN4*nthreads);

        cache_size = 0;
        for (int i = 0; i < nShells; i++) {
          size_t n;
          int shls[4]{i,i,i,i};
          if (groupedBasisSet.forceCart) {
            n = int2e_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          } else {
            n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          }
          cache_size = std::max(cache_size, n);
        }
        cacheAll = CQMemManager::get().malloc<double>(cache_size*nthreads);

      } else {

        // Compute the mappings from primitives to CGTOs
        BasisSet primitives(basisSet.uncontractBasis());
        OnePInts<double> overlap(primitives.nBasis);

        overlap.computeAOInts(primitives, mol, emPert, OVERLAP, options);

        double *mapPrim2Cont = CQMemManager::get().malloc<double>(primitives.nBasis*basisSet.nBasis);
        basisSet.makeMapPrim2Cont(overlap.pointer(), mapPrim2Cont);

#ifdef __DEBUGERI__
        prettyPrintSmart(std::cout, "mapPrim2Cont", mapPrim2Cont,
                         basisSet.nBasis, basisSet.nPrimitive, basisSet.nBasis);
#endif

        // Clear objects
        for (double *p : coefBlocks_) {
          if (p) CQMemManager::get().free(p);
        }
        coefBlocks_.clear();
        coefBlocks_.resize(groupedBasisSet.nShell, nullptr);
        shellPrims_.clear();
        shellPrims_.reserve(groupedBasisSet.nShell);

        for (size_t P(0); P < groupedBasisSet.nShell; P++) {
          // Gather contraction coefficients
          libint2::Shell &shellP = groupedBasisSet.shells[P];
          size_t pContrSize = shellP.ncontr();
          size_t pAMSize = shellP.contr[0].size();
          if (pContrSize == 1) {
            maxNprimAMSize_ = std::max(maxNprimAMSize_, pAMSize);
            shellPrims_.emplace_back(1, shellP);
            coefBlocks_[P] = CQMemManager::get().malloc<double>(1);
            coefBlocks_[P][0] = 1.0;
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

          maxNprimAMSize_ = std::max(maxNprimAMSize_, pNprim * pAMSize);

          size_t pBegin = groupedBasisSet.mapSh2Bf[P];
          coefBlocks_[P] = CQMemManager::get().malloc<double>(pContrSize * pNprim);
          for (size_t c = 0; c < pContrSize; c++) {
            for (size_t i = 0; i < pNprim; i++) {
              coefBlocks_[P][i + c * pNprim]
                  = mapPrim2Cont[pBegin + c * pAMSize +
                                  groupedBasisSet.primitives[primsP[i]] * basisSet.nBasis];
            }
          }
          shellPrims_.push_back(std::move(primsP));
        }
        CQMemManager::get().free(mapPrim2Cont);

        workBlocks.resize(nthreads, nullptr);

        size_t maxL = basisSet.maxL;
        maxAMSize_ = basisSet.forceCart ?
                        (maxL + 1) * (maxL + 2) / 2 :
                        (2 * maxL + 1);
        size_t primAllocSize =  maxNcontrAMSize_ * maxNcontrAMSize_
                              * maxNcontrAMSize_ * maxNcontrAMSize_
                              + maxNprimAMSize_ * maxNcontrAMSize_
                              * maxNcontrAMSize_ * maxNcontrAMSize_
                              + maxNprimAMSize_ * maxNcontrAMSize_
                              * maxNcontrAMSize_ * maxAMSize_
                              + maxNprimAMSize_ * maxNcontrAMSize_
                              * maxAMSize_ * maxAMSize_
                              + maxNprimAMSize_ * maxAMSize_
                              * maxAMSize_ * maxAMSize_;
        for (size_t i = 0; i < nthreads; i++) {
          workBlocks[i] = CQMemManager::get().malloc<double>(primAllocSize);
        }

      }
    }

    switch (alg_) {
    case CHOLESKY_ALG::TRADITIONAL:
      computeCD_Traditional(groupedBasisSet);
      break;

    case CHOLESKY_ALG::DYNAMIC_ALL:
      computeCDPivots_DynamicAll(groupedBasisSet);
      std::cout << "  Cholesky-RI auxiliary dimension = " << pivots_.size() << std::endl;
      computePivotRI(groupedBasisSet);
      break;

    case CHOLESKY_ALG::SPAN_FACTOR:
      computeCDPivots_SpanFactor(groupedBasisSet);
      std::cout << "  Cholesky-RI auxiliary dimension = " << pivots_.size() << std::endl;
      computePivotRI(groupedBasisSet);
      break;

    case CHOLESKY_ALG::DYNAMIC_ERI:
      computeCDPivots_DynamicERI(groupedBasisSet);
      std::cout << "  Cholesky-RI auxiliary dimension = " << pivots_.size() << std::endl;
      computePivotRI(groupedBasisSet);
      break;

    case CHOLESKY_ALG::SPAN_FACTOR_REUSE:
      computeCDPivots_SpanFactorReuse(groupedBasisSet);
      std::cout << "  Cholesky-RI auxiliary dimension = " << pivots_.size() << std::endl;
      computePivotRI(groupedBasisSet);
      break;
    }

    if (generalContraction_) {

      for (double *p : coefBlocks_) {
        if (p) CQMemManager::get().free(p);
      }
      coefBlocks_.clear();

      if (options.Libcint) {
        CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);

      } else {

        for (size_t i = 0; i < nthreads; i++) {
          CQMemManager::get().free(workBlocks[i]);
        }
      }

    }

    auto durCholeskyRI = tock(topCholeskyRI);
    std::cout << std::endl << "Cholesky-RI duration = " << durCholeskyRI << " s " << std::endl;

    std::cout << std::endl << BannerEnd << std::endl;

  }; // InCoreCholeskyRIERI<double>::computeAOInts

}; // namespace ChronusQ

