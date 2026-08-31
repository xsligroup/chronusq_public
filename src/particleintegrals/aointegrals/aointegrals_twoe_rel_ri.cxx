/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <libcint.hpp>
#include <util/timer.hpp>

#define _REPORT_INCORE_INTEGRAL_TIMINGS

namespace ChronusQ {

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


  template <>
  void IncoreDCRITPIList<double>::computeCholeskyDCERI_CINT(
      BasisSet &originalBasisSet, Molecule &molecule_) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(LLLL_);
    size_t NBRI = aux.nRIBasis();

    if (NBRI == 0)
      CErr("Cholesky decomposition of (LL|LL) has not performed but CD-(SS|LL) is requested.");

    // Group pivots into shells
    const std::vector<size_t> &pivots = aux.pivots();
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
    pivotIndicesByShell = InCoreCholeskyRIERI<double>::groupPivotsByShell(basisSet_, pivots);

    size_t pivotShellSize = pivotIndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivotShells;
    pivotShells.reserve(pivotShellSize);
    for (auto &shell_pivot : pivotIndicesByShell) {
      pivotShells.push_back(shell_pivot.first);
    }

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
      return a.size() < b.size();
    })->size();

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    CQMemManager &memManager_ = CQMemManager::get();
    int *atm = memManager_.template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = memManager_.template malloc<int>(nShells * BAS_SLOTS);
    double *env = memManager_.template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      n = int2e_ipvip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
    }


    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    buffN4 *= 9;

    double *buffAll = memManager_.malloc<double>(buffN4*nthreads);
    double *cacheAll = memManager_.malloc<double>(cache_size*nthreads);

    /* Dirac-Coulomb Integrals */
    // Dirac-Coulomb ∇_i∇_j(ij|kl)

    for (size_t i = 0; i < 4; i++) {
      set_asymm_term(i, std::make_shared<InCoreAsymmRITPI<double>>(NB, LLLL_));
      asymm_term(i)->malloc();
      asymm_term(i)->clear();
    }

    int AxBx = 0;
    int AxBy = 1;
    int AxBz = 2;
    int AyBx = 3;
    int AyBy = 4;
    int AyBz = 5;
    int AzBx = 6;
    int AzBy = 7;
    int AzBz = 8;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
    auto topERIDC = tick();
#endif
    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,i,j,bf1,bf2;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) {

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

          n2 = basisSet_.shells[s2].size(); // Size of Shell 2

          for (size_t I = 0, IPQ = 0; I < pivotShellSize; I++, s1234++) {

            // Round Robbin work distribution
            #ifdef _OPENMP
            if( s1234 % nthreads != thread_id ) continue;
            #endif

            const auto &RSpair = pivotShells[I];
            auto &shell_pivot = pivotIndicesByShell[RSpair];

            size_t R = RSpair.first;
            size_t S = RSpair.second;

            size_t rBegin = basisSet_.mapSh2Bf[R];
            size_t sBegin = basisSet_.mapSh2Bf[S];
            size_t rSize = basisSet_.shells[R].size();
            size_t sSize = basisSet_.shells[S].size();
            size_t rEnd = rBegin + rSize;
            size_t sEnd = sBegin + sSize;

            size_t nQuad = n1*n2*rSize*sSize;

            shls[0] = int(s1);
            shls[1] = int(s2);
            shls[2] = int(R);
            shls[3] = int(S);

            if (basisSet_.forceCart) {
              if(int2e_ipvip1_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
            } else {
              if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
            }

            for (auto &pivot_index : shell_pivot) {

              auto rs = anaSquare(pivots[pivot_index], NB);
              size_t r = rs.first;
              size_t s = rs.second;
              size_t rsShift = ((r-rBegin) + (s-sBegin) * rSize) * n1 * n2;

              for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) {
                size_t iBegin = (s1 == s2 ? j : 0);
                for(i = iBegin, bf1 = bf1_s + iBegin ; i < n1; ++i, bf1++) {

                  size_t ijkl = i + j * n1 + rsShift;

                  // ∇A∙∇B(ij|kl)
                  auto dAdotdB = buff[AxBx*nQuad+ijkl] + buff[AyBy*nQuad+ijkl] + buff[AzBz*nQuad+ijkl];
                  // ∇Ax∇B(ijkl)
                  auto dAcrossdB_x =  buff[AyBz*nQuad+ijkl] - buff[AzBy*nQuad+ijkl];
                  auto dAcrossdB_y = -buff[AxBz*nQuad+ijkl] + buff[AzBx*nQuad+ijkl];
                  auto dAcrossdB_z =  buff[AxBy*nQuad+ijkl] - buff[AyBx*nQuad+ijkl];

                  auto QJI = pivot_index + toCompound(bf2,bf1) * NBRI;

                  // ∇A∙∇B(ij|kl) followed by ∇Ax∇B(ij|kl) X, Y, and Z
                  // (ij|kl)
                  components_[0]->pointer()[QJI] = dAdotdB;
                  components_[1]->pointer()[QJI] = dAcrossdB_x;
                  components_[2]->pointer()[QJI] = dAcrossdB_y;
                  components_[3]->pointer()[QJI] = dAcrossdB_z;
                  
                }
              }
            }


          }; // PQ
        }; // s2
      }; // s1

    }; // omp region

#ifdef  _REPORT_INCORE_INTEGRAL_TIMINGS
auto durERIDC = tock(topERIDC);
std::cout << "Libcint-ERI-Dirac-Coulomb duration   = " << durERIDC << std::endl;
#endif


    size_t NB2comp    = NB*(NB+1)/2;
    size_t NBcompNBRI = NB2comp*NBRI;
    auto ijK = CQMemManager::get().malloc<double>(NBcompNBRI);
    for (size_t i = 0; i < 4; i++) {
      // S^{-1/2}(Q|ij)
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                 NBRI,NB2comp,NBRI,double(1.),aux.twoIndexERI()->pointer(),NBRI,
                 components_[i]->pointer(),NBRI,double(0.),ijK,NBRI);

      for (size_t pq = 0; pq < NB2comp; pq++) {
        auto pqAna = anaCompound(pq);
        std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI],
                  components_[i]->pointer()+NBRI*toSquare(pqAna.second, pqAna.first, NB));
        if (i == 0)
          std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI],
                    components_[i]->pointer()+NBRI*toSquare(pqAna.first, pqAna.second, NB));
        else
          std::transform(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI],
                         components_[i]->pointer()+NBRI*toSquare(pqAna.first, pqAna.second, NB),
                         [](const double &a){ return -a; });
      }
    }

    CQMemManager::get().free(ijK);

    // Report error
    // if (is4CRI) {
    //   std::shared_ptr<Integrals<double>> aoints_double =
    //       std::dynamic_pointer_cast<Integrals<double>>(aoints);
    //   std::shared_ptr<InCoreRelERI<double>> relERI_double =
    //       std::dynamic_pointer_cast<InCoreRelERI<double>>(aoints_double->TPI);
    //
    //   double largeValueThreshold = 1e-8; ///< Threshold for elements computing percentage error
    //   OPTOPT( largeValueThreshold = input.getData<double>("INTS/SCHWARZ"); )
    //
    //   relERI_double->RI_direct_error(*basis, mol, emPert, ELECTRON_REPULSION, ssOptions.hamiltonianOptions, largeValueThreshold);
    // }

  } // IncoreDCRITPIList<double>::computeCholeskyDCERI_CINT

  template <>
  void IncoreDCRITPIList<dcomplex>::computeCholeskyDCERI_CINT(
      BasisSet &originalBasisSet, Molecule &molecule_) {
    CErr("GIAO integral evaluation is NOT implemented in class InCoreAsymmRITPI.");
  } // IncoreDCRITPIList<dcomplex>::computeCholeskyDCERI_CINT



  template <>
  void IncoreGauntRITPIList<double>::computeCholeskyGauntERI_CINT(
      BasisSet &originalBasisSet, Molecule &molecule_) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(LLLL_);
    size_t NBRI = aux.nRIBasis();

    if (NBRI == 0)
      CErr("Cholesky decomposition of (LL|LL) has not performed but CD-Gaunt is requested.");

    // Group pivots into shells
    const std::vector<size_t> &pivots = aux.pivots();
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
        pivotIndicesByShell = InCoreCholeskyRIERI<double>::groupPivotsByShell(basisSet_, pivots);

    size_t pivotShellSize = pivotIndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivotShells;
    pivotShells.reserve(pivotShellSize);
    for (auto &shell_pivot : pivotIndicesByShell) {
      pivotShells.push_back(shell_pivot.first);
    }

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    CQMemManager &memManager_ = CQMemManager::get();
    int *atm = memManager_.template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = memManager_.template malloc<int>(nShells * BAS_SLOTS);
    double *env = memManager_.template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      n = int2e_ip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
    }


    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    buffN4 *= 3;

    double *buffAll = memManager_.malloc<double>(buffN4*nthreads);
    double *cacheAll = memManager_.malloc<double>(cache_size*nthreads);

    /* (NABLA i j|R12 |k l) Integrals */
    // ∇_i(ij|kl)

    for (size_t i = 0; i < 3; i++) {
      set_asymm_term(i, std::make_shared<InCoreAsymmRITPI<double>>(NB, LLLL_));
      asymm_term(i)->malloc();
      asymm_term(i)->clear();
    }

    int Ax = 0;
    int Ay = 1;
    int Az = 2;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
    auto topERIDC = tick();
#endif
#pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,i,j,bf1,bf2;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll + buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) {

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

        for(size_t s2(0), bf2_s(0); s2 < nShells; bf2_s+=n2, s2++) {

          n2 = basisSet_.shells[s2].size(); // Size of Shell 2

          for (size_t I = 0, IPQ = 0; I < pivotShellSize; I++, s1234++) {

            // Round Robbin work distribution
#ifdef _OPENMP
            if( s1234 % nthreads != thread_id ) continue;
#endif

            const auto &RSpair = pivotShells[I];
            auto &shell_pivot = pivotIndicesByShell[RSpair];

            size_t R = RSpair.first;
            size_t S = RSpair.second;

            size_t rBegin = basisSet_.mapSh2Bf[R];
            size_t sBegin = basisSet_.mapSh2Bf[S];
            size_t rSize = basisSet_.shells[R].size();
            size_t sSize = basisSet_.shells[S].size();
            size_t rEnd = rBegin + rSize;
            size_t sEnd = sBegin + sSize;

            size_t nQuad = n1*n2*rSize*sSize;

            shls[0] = int(s1);
            shls[1] = int(s2);
            shls[2] = int(R);
            shls[3] = int(S);

            if (basisSet_.forceCart) {
              if(int2e_ip1_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
            } else {
              if(int2e_ip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
            }

            for (auto &pivot_index : shell_pivot) {

              auto rs = anaSquare(pivots[pivot_index], NB);
              size_t r = rs.first;
              size_t s = rs.second;
              size_t rsShift = ((r-rBegin) + (s-sBegin) * rSize) * n1 * n2;

              for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) {
                size_t iBegin = (s1 == s2 ? j : 0);
                for(i = iBegin, bf1 = bf1_s + iBegin ; i < n1; ++i, bf1++) {

                  size_t ijkl = i + j * n1 + rsShift;

                  auto QJI = pivot_index + (bf1 + bf2 * NB) * NBRI;

                  // ∇A(ij|kl) X, Y, and Z
                  components_[0]->pointer()[QJI] = buff[Ax*nQuad+ijkl];
                  components_[1]->pointer()[QJI] = buff[Ay*nQuad+ijkl];
                  components_[2]->pointer()[QJI] = buff[Az*nQuad+ijkl];

                }
              }
            }


          }; // PQ
        }; // s2
      }; // s1

    }; // omp region

#ifdef  _REPORT_INCORE_INTEGRAL_TIMINGS
    auto durERIDC = tock(topERIDC);
    std::cout << "Libcint-ERI-Dirac-Coulomb duration   = " << durERIDC << std::endl;
#endif


    size_t NB2     = NB*NB;
    size_t NB2NBRI = NB2*NBRI;
    auto ijK = CQMemManager::get().malloc<double>(NB2NBRI);
    for (size_t i = 0; i < 3; i++) {
      // S^{-1/2}(Q|∇ij)
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                 NBRI,NB2,NBRI,double(1.),aux.twoIndexERI()->pointer(),NBRI,
                 components_[i]->pointer(),NBRI,double(0.),ijK,NBRI);
      std::copy_n(ijK, NB2NBRI, components_[i]->pointer());
    }

    CQMemManager::get().free(ijK);

  } // IncoreGauntRITPIList<double>::computeCholeskyGauntERI_CINT

  template <>
  void IncoreGauntRITPIList<dcomplex>::computeCholeskyGauntERI_CINT(
      BasisSet &originalBasisSet, Molecule &molecule_) {
    CErr("GIAO integral evaluation is NOT implemented in class InCoreAsymmRITPI.");
  } // IncoreGauntRITPIList<dcomplex>::computeCholeskyGauntERI_CINT


  template <>
  void InCoreRelERI<double>::RI_direct_error(BasisSet &originalBasisSet, Molecule &molecule_,
      EMPerturbation&, OPERATOR, const HamiltonianOptions &hamiltonianOptions, double largeValueThreshold) const {


    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    CQMemManager &memManager_ = CQMemManager::get();
    int *atm = memManager_.template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = memManager_.template malloc<int>(nShells * BAS_SLOTS);
    double *env = memManager_.template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      if (basisSet_.forceCart) {
        n = int2e_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        if(hamiltonianOptions.DiracCoulomb) {
          n = int2e_ipvip1_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gaunt) {
          n = int2e_ip1ip2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.DiracCoulombSSSS) {
          n = int2e_ipvip1ipvip2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gauge) {
          n = int2e_gauge_r1_sps1sps2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
      } else {
        n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        if(hamiltonianOptions.DiracCoulomb) {
          n = int2e_ipvip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gaunt) {
          n = int2e_ip1ip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.DiracCoulombSSSS) {
          n = int2e_ipvip1ipvip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gauge) {
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
    }


    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    if(hamiltonianOptions.DiracCoulombSSSS)
      buffN4 *= 81;
    else if (hamiltonianOptions.Gauge)
      buffN4 *= 16;
    else if (hamiltonianOptions.DiracCoulomb or hamiltonianOptions.Gaunt)
      buffN4 *= 9;

    double *buffAll, *buffAll2;
    if(hamiltonianOptions.Gauge) buffAll = memManager_.malloc<double>(buffN4*nthreads*2);
    else buffAll = memManager_.malloc<double>(buffN4*nthreads);
    if(hamiltonianOptions.Gauge) buffAll2 = memManager_.malloc<double>(buffN4*nthreads*2);
    else buffAll2 = memManager_.malloc<double>(buffN4*nthreads);

    double *cacheAll = memManager_.malloc<double>(cache_size*nthreads);

    std::cout<<"Using Libcint "<<std::endl;


    std::shared_ptr<InCoreRITPI<double>> riLLLL
        = std::dynamic_pointer_cast<InCoreRITPI<double>>(LLLL_);

    if (riLLLL) {

    size_t NBRI = riLLLL->nRIBasis();

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
    auto topERI4 = tick();
#endif

    double maxError = 0.0, meanAbsError = 0.0, maxElem = 0.0, meanAbsElem = 0.0, maxPercent = 0.0, meanAbsPercent = 0.0;
    #pragma omp parallel reduction( max: maxError) reduction(+: meanAbsError) reduction( max: maxElem) reduction(+: meanAbsElem) reduction( max: maxPercent) reduction(+: meanAbsPercent)
    {
      int thread_id = GetThreadID();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll+buffN4*thread_id;
      double *buff2 = buffAll2+buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells;
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

        shls[0] = int(s1);
        shls[1] = int(s2);
        shls[2] = int(s3);
        shls[3] = int(s4);

        if (basisSet_.forceCart) {
          if(int2e_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        } else {
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        }

//        std::cout << "ERI4: " << s1 << " " << s2 << " " << s3 << " " << s4 << std::endl;
//
//        std::cout << "buff1" << std::endl;
//        ijkl = 0ul;
//        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
//          for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
//            for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
//              for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++)
//                std::cout << "(" << bf1 << "," << bf2 << "|" << bf3 << "," << bf4 << ") = " << buff[ijkl++] << std::endl;
//        prettyPrintSmart(std::cout,"Buff", buff,n1 * n2,n3 * n4,n1 * n2);

        // compute RI reformated integrals
        ijkl = 0ul;
        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
          for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
            for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
              for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++)
                blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                           1,1,NBRI,1.,riLLLL->pointer() + NBRI*(bf1 + bf2 * NB),NBRI,
                           riLLLL->pointer() + NBRI*(bf3 + bf4 * NB),NBRI,
                           0.,&buff2[ijkl++],1);
//        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
//                   n1*n2,n3*n4,NBRI,1.,riLLLL->pointer() + NBRI*(bf1_s + bf2_s * NB),NBRI,
//                   riLLLL->pointer() + NBRI*(bf3_s + bf4_s * NB),NBRI,0.,buff2,n1*n2);

//        std::cout << "buff2" << std::endl;
//        ijkl = 0ul;
//        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
//          for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
//            for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
//              for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++)
//                std::cout << "(" << bf1 << "," << bf2 << "|" << bf3 << "," << bf4 << ") = " << buff2[ijkl++] << std::endl;
//        prettyPrintSmart(std::cout,"Buff2", buff2,n1 * n2,n3 * n4,n1 * n2);

        blas::axpy(n1*n2*n3*n4,-1.,buff,1,buff2,1);

        // compute max element
        auto minmax = std::minmax_element(buff, buff+n1*n2*n3*n4);
        maxElem = std::max(maxElem, std::max(-*minmax.first, *minmax.second));
        // compute mean absolute element in buff
        double blockMeanAbsElem = std::accumulate(buff, buff + n1*n2*n3*n4, 0., [](auto a, auto b) { return a + std::abs(b); });

        // compute max error
        minmax = std::minmax_element(buff2, buff2+n1*n2*n3*n4);
        maxError = std::max(maxError, std::max(-*minmax.first, *minmax.second));
        // compute mean absolute error in buff2
        double blockMeanAbsError = std::accumulate(buff2, buff2 + n1*n2*n3*n4, 0., [](auto a, auto b) { return a + std::abs(b); });

        // compute percent error
        double blockMeanPercent = 0.0;
        size_t nQuad = n1*n2*n3*n4;
        for (size_t i = 0; i < nQuad; i++)
          if (std::abs(buff[i]) > largeValueThreshold) { // ignore small elements (relative error is meaningless
            double percent = std::abs(buff2[i]) / std::abs(buff[i]);
            blockMeanPercent += percent;
            maxPercent = std::max(maxPercent, percent);
          }

        int multiplier = 8;
        if (s1 == s2) multiplier /= 2;
        if (s3 == s4) multiplier /= 2;
        if (s1 == s3 and s2 == s4) multiplier /= 2;
        meanAbsElem += multiplier * blockMeanAbsElem;
        meanAbsError += multiplier * blockMeanAbsError;
        meanAbsPercent += multiplier * blockMeanPercent;

      }; // s4
      }; // s3
      }; // s2
      }; // s1

    }; // omp region

    std::cout << "LLLL element MAX = " << maxElem << std::endl;
    std::cout << "LLLL error   MAX = " << maxError << std::endl;
    std::cout << "LLLL percent MAX = " << maxPercent << std::endl;
    std::cout << "LLLL element MAE = " << meanAbsElem/NB4 << std::endl;
    std::cout << "LLLL error   MAE = " << meanAbsError/NB4 << std::endl;
    std::cout << "LLLL percent MAE = " << meanAbsPercent/NB4 << std::endl;

#ifdef _REPORT_INCORE_CONTRACTION_TIMINGS
    auto durERI4 = tock(topERI4);
    //std::cout << "L = "<< basisSet_.shells[s1].contr[0].l<<" "<<basisSet_.shells[s2].contr[0].l<<" "
    //	               << basisSet_.shells[s3].contr[0].l<<" "<<basisSet_.shells[s4].contr[0].l<<std::endl;
    std::cout << "Libcint-ERI4 duration   = " << durERI4 << std::endl;
#endif

#ifdef __DEBUGERI__
    // Debug output of the ERIs
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "Libcint ERI (ab|cd)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)(i, j, k, l) << std::endl;
    };
#endif // __DEBUGERI__

    } // if (riLLLL)


    std::shared_ptr<IncoreDCRITPIList<double>> riSSLL
        = std::dynamic_pointer_cast<IncoreDCRITPIList<double>>(DC_terms());

    /* Dirac-Coulomb Integrals */
    if(hamiltonianOptions.DiracCoulomb and riSSLL) { // Dirac-Coulomb ∇_i∇_j(ij|kl)

      size_t NBRI = riSSLL->nRIBasis();
      riLLLL = riSSLL->LLLL_term();

      int AxBx = 0;
      int AxBy = 1;
      int AxBz = 2;
      int AyBx = 3;
      int AyBy = 4;
      int AyBz = 5;
      int AzBx = 6;
      int AzBy = 7;
      int AzBz = 8;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
      auto topERIDC = tick();
#endif

      std::vector<double> maxErrors(nthreads*4, 0.0), sumAbsErrors(nthreads*4, 0.0);
      std::vector<double> maxElems(nthreads*4, 0.0), sumAbsElems(nthreads*4, 0.0);
      std::vector<double> maxPercents(nthreads*4, 0.0), sumAbsPercents(nthreads*4, 0.0);
#pragma omp parallel
      {
        int thread_id = GetThreadID();

        size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
        int shls[4];
        double *buff = buffAll + buffN4*thread_id;
        double *buff2 = buffAll2 + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
        double expo1,expo2;

        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) {

          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
          expo1 = basisSet_.shells[s1].alpha[0]; // Exponent of Shell 1

          for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

            n2 = basisSet_.shells[s2].size(); // Size of Shell 2
            expo2 = basisSet_.shells[s2].alpha[0]; // Exponent of Shell 2

            for(size_t s3(0), bf3_s(0); s3 < nShells ; bf3_s+=n3, s3++) {

              n3 = basisSet_.shells[s3].size(); // Size of Shell 3

              for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {

                n4 = basisSet_.shells[s4].size(); // Size of Shell 4

                // Round Robbin work distribution
#ifdef _OPENMP
                if( s1234 % nthreads != thread_id ) continue;
#endif

                shls[0] = int(s1);
                shls[1] = int(s2);
                shls[2] = int(s3);
                shls[3] = int(s4);

                if (basisSet_.forceCart) {
                  if(int2e_ipvip1_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
                } else {
                  if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
                }

                ijkl = 0ul;
                auto nQuad = n1*n2*n3*n4;
                for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
                  for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
                    for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
                      for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {

#ifdef __DEBUGERI__

                        std::cout << std::scientific << std::setprecision(16);
  	  std::cout <<"Libcint ∇A∙∇B(ij|kl)"<<std::endl;
  	  std::cout<<buff[AxBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AxBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AxBz*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBz*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBz*nQuad+ijkl]<<std::endl;

#endif

                        // ∇A∙∇B(ij|kl)
                        auto dAdotdB = buff[AxBx*nQuad+ijkl] + buff[AyBy*nQuad+ijkl] + buff[AzBz*nQuad+ijkl];
                        // ∇Ax∇B(ijkl)
                        auto dAcrossdB_x =  buff[AyBz*nQuad+ijkl] - buff[AzBy*nQuad+ijkl];
                        auto dAcrossdB_y = -buff[AxBz*nQuad+ijkl] + buff[AzBx*nQuad+ijkl];
                        auto dAcrossdB_z =  buff[AxBy*nQuad+ijkl] - buff[AyBx*nQuad+ijkl];

//                        auto IJKL = bf1 + bf2*NB + bf3*NB2 + bf4*NB3;
//                        auto IJLK = bf1 + bf2*NB + bf4*NB2 + bf3*NB3;
//                        auto JIKL = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
//                        auto JILK = bf2 + bf1*NB + bf4*NB2 + bf3*NB3;

                        // ∇A∙∇B(ij|kl) followed by ∇Ax∇B(ij|kl) X, Y, and Z
                        // (ij|kl)
//                        (*this)[0].pointer()[IJKL] = dAdotdB;
//                        (*this)[1].pointer()[IJKL] = dAcrossdB_x;
//                        (*this)[2].pointer()[IJKL] = dAcrossdB_y;
//                        (*this)[3].pointer()[IJKL] = dAcrossdB_z;
//                        // (ij|lk)
//                        (*this)[0].pointer()[IJLK] = dAdotdB;
//                        (*this)[1].pointer()[IJLK] = dAcrossdB_x;
//                        (*this)[2].pointer()[IJLK] = dAcrossdB_y;
//                        (*this)[3].pointer()[IJLK] = dAcrossdB_z;
//                        // (ji|kl)
//                        (*this)[0].pointer()[JIKL] = dAdotdB;
//                        (*this)[1].pointer()[JIKL] = -dAcrossdB_x;
//                        (*this)[2].pointer()[JIKL] = -dAcrossdB_y;
//                        (*this)[3].pointer()[JIKL] = -dAcrossdB_z;
//                        // (ji|lk)
//                        (*this)[0].pointer()[JILK] = dAdotdB;
//                        (*this)[1].pointer()[JILK] = -dAcrossdB_x;
//                        (*this)[2].pointer()[JILK] = -dAcrossdB_y;
//                        (*this)[3].pointer()[JILK] = -dAcrossdB_z;

                        // (ij|kl)
                        buff2[ijkl] = dAdotdB;
                        buff2[ijkl +     nQuad] = dAcrossdB_x;
                        buff2[ijkl + 2 * nQuad] = dAcrossdB_y;
                        buff2[ijkl + 3 * nQuad] = dAcrossdB_z;
                        ijkl++;

                      }; // ijkl loop

                // compute RI reformated integrals
                ijkl = 0ul;
                for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
                  for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
                    for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
                      for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(0)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riLLLL->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(1)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riLLLL->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +     nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(2)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riLLLL->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 2 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(3)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riLLLL->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 3 * nQuad], 1);
                        ijkl++;
                      }

                blas::axpy(4*nQuad,-1.,buff2,1,buff,1);

                double expoScalar = 4 * expo1 * expo2;
                for (size_t t = 0; t < 4; t++) {
                  // compute max element
                  auto minmax = std::minmax_element(buff2+t*nQuad, buff2+(t+1)*nQuad);
                  maxElems[thread_id + t*nthreads] = std::max(maxElems[thread_id + t*nthreads], std::max(-*minmax.first, *minmax.second)/expoScalar);
                  // compute mean absolute element in buff2
                  double sumAbsElem = std::accumulate(buff2+t*nQuad, buff2+(t+1)*nQuad, 0., [](double a, double b) { return a + std::abs(b); })/expoScalar;

                  // compute max error
                  minmax = std::minmax_element(buff+t*nQuad, buff+(t+1)*nQuad);
                  maxErrors[thread_id + t*nthreads] = std::max(maxErrors[thread_id + t*nthreads], std::max(-*minmax.first, *minmax.second)/expoScalar);
                  // compute mean absolute error in buff
                  double sumAbsError = std::accumulate(buff+t*nQuad, buff+(t+1)*nQuad, 0., [](double a, double b) { return a + std::abs(b); })/expoScalar;

                  // compute percent error
                  double sumAbsPercent = 0.0;
                  for (size_t i = 0; i < nQuad; i++)
                    if (std::abs(buff2[i+t*nQuad]) > largeValueThreshold * expoScalar) { // ignore small elements (relative error is meaningless
                      double percent = std::abs(buff[i+t*nQuad]) / std::abs(buff2[i+t*nQuad]);
                      sumAbsPercent += percent;
                      maxPercents[thread_id + t*nthreads] = std::max(maxPercents[thread_id + t*nthreads], percent);
                    }

                  int multiplier = 4;
                  if (s1 == s2) multiplier /= 2;
                  if (s3 == s4) multiplier /= 2;
                  sumAbsElems[thread_id + t*nthreads] += multiplier * sumAbsElem;
                  sumAbsErrors[thread_id + t*nthreads] += multiplier * sumAbsError;
                  sumAbsPercents[thread_id + t*nthreads] += multiplier * sumAbsPercent;
                }
              }; // s4
            }; // s3
          }; // s2
        }; // s1

      }; // omp region
      for (size_t t = 0; t < 4; t++) {
        double scale_denom = 4*SpeedOfLight()*SpeedOfLight();
        double maxElem = *std::max_element(&maxElems[t*nthreads],&maxElems[(t+1)*nthreads]);
        std::cout << "SSLL-" << t << " element MAX = " << maxElem << "  scale by 1/(2mc)^2 = " << maxElem/scale_denom << std::endl;
        double maxError = *std::max_element(&maxErrors[t*nthreads],&maxErrors[(t+1)*nthreads]);
        std::cout << "SSLL-" << t << " error   MAX = " << maxError << "  scale by 1/(2mc)^2 = " << maxError/scale_denom << std::endl;
        double maxPercent = *std::max_element(&maxPercents[t*nthreads],&maxPercents[(t+1)*nthreads]);
        std::cout << "SSLL-" << t << " percent MAX = " << maxPercent << std::endl;
        double meanAbsElem = std::accumulate(&sumAbsElems[t*nthreads],&sumAbsElems[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSLL-" << t << " element MAE = " << meanAbsElem << "  scale by 1/(2mc)^2 = " << meanAbsElem/scale_denom << std::endl;
        double meanAbsError = std::accumulate(&sumAbsErrors[t*nthreads],&sumAbsErrors[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSLL-" << t << " error   MAE = " << meanAbsError << "  scale by 1/(2mc)^2 = " << meanAbsError/scale_denom << std::endl;
        double meanAbsPercent = std::accumulate(&sumAbsPercents[t*nthreads],&sumAbsPercents[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSLL-" << t << " percent MAE = " << meanAbsPercent << std::endl;
      }

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
      auto durERIDC = tock(topERIDC);
      std::cout << "Libcint-ERI-Dirac-Coulomb duration   = " << durERIDC << std::endl;
#endif


#if 0
      std::cout << std::scientific << std::setprecision(16);
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
#endif

    } // Dirac-Coulomb ∇_i∇_j(ij|kl)









    std::shared_ptr<IncoreSSSSRITPIList<double>> riSSSS
        = std::dynamic_pointer_cast<IncoreSSSSRITPIList<double>>(SSSS_terms());

    /* Dirac-Coulomb (SSSS) Integrals */
    if(hamiltonianOptions.DiracCoulombSSSS and riSSSS) { // Dirac-Coulomb ∇_i∇_j∇_k∇_l(ij|kl)

      riSSLL = riSSSS->DC_term();
      size_t NBRI = riSSSS->nRIBasis();

      auto nERIRef = 0;
      if(hamiltonianOptions.DiracCoulomb) nERIRef +=4;
      if(hamiltonianOptions.Gaunt) nERIRef += 19;

      int AxBxCxDx = 0;
      int AxBxCxDy = 1;
      int AxBxCxDz = 2;
      int AxBxCyDx = 3;
      int AxBxCyDy = 4;
      int AxBxCyDz = 5;
      int AxBxCzDx = 6;
      int AxBxCzDy = 7;
      int AxBxCzDz = 8;

      int AxByCxDx = 9;
      int AxByCxDy = 10;
      int AxByCxDz = 11;
      int AxByCyDx = 12;
      int AxByCyDy = 13;
      int AxByCyDz = 14;
      int AxByCzDx = 15;
      int AxByCzDy = 16;
      int AxByCzDz = 17;

      int AxBzCxDx = 18;
      int AxBzCxDy = 19;
      int AxBzCxDz = 20;
      int AxBzCyDx = 21;
      int AxBzCyDy = 22;
      int AxBzCyDz = 23;
      int AxBzCzDx = 24;
      int AxBzCzDy = 25;
      int AxBzCzDz = 26;



      int AyBxCxDx = 27;
      int AyBxCxDy = 28;
      int AyBxCxDz = 29;
      int AyBxCyDx = 30;
      int AyBxCyDy = 31;
      int AyBxCyDz = 32;
      int AyBxCzDx = 33;
      int AyBxCzDy = 34;
      int AyBxCzDz = 35;

      int AyByCxDx = 36;
      int AyByCxDy = 37;
      int AyByCxDz = 38;
      int AyByCyDx = 39;
      int AyByCyDy = 40;
      int AyByCyDz = 41;
      int AyByCzDx = 42;
      int AyByCzDy = 43;
      int AyByCzDz = 44;

      int AyBzCxDx = 45;
      int AyBzCxDy = 46;
      int AyBzCxDz = 47;
      int AyBzCyDx = 48;
      int AyBzCyDy = 49;
      int AyBzCyDz = 50;
      int AyBzCzDx = 51;
      int AyBzCzDy = 52;
      int AyBzCzDz = 53;



      int AzBxCxDx = 54;
      int AzBxCxDy = 55;
      int AzBxCxDz = 56;
      int AzBxCyDx = 57;
      int AzBxCyDy = 58;
      int AzBxCyDz = 59;
      int AzBxCzDx = 60;
      int AzBxCzDy = 61;
      int AzBxCzDz = 62;

      int AzByCxDx = 63;
      int AzByCxDy = 64;
      int AzByCxDz = 65;
      int AzByCyDx = 66;
      int AzByCyDy = 67;
      int AzByCyDz = 68;
      int AzByCzDx = 69;
      int AzByCzDy = 70;
      int AzByCzDz = 71;

      int AzBzCxDx = 72;
      int AzBzCxDy = 73;
      int AzBzCxDz = 74;
      int AzBzCyDx = 75;
      int AzBzCyDy = 76;
      int AzBzCyDz = 77;
      int AzBzCzDx = 78;
      int AzBzCzDy = 79;
      int AzBzCzDz = 80;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
      auto topERIDCSSSS = tick();
#endif

      std::vector<double> maxErrors(nthreads*16, 0.0), sumAbsErrors(nthreads*16, 0.0);
      std::vector<double> maxElems(nthreads*16, 0.0), sumAbsElems(nthreads*16, 0.0);
      std::vector<double> maxPercents(nthreads*16, 0.0), sumAbsPercents(nthreads*16, 0.0);
#pragma omp parallel
      {
        int thread_id = GetThreadID();

        size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
        int shls[4];
        double *buff = buffAll + buffN4*thread_id;
        double *buff2 = buffAll2 + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
        double expo1,expo2,expo3,expo4;

        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) {

          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
          expo1 = basisSet_.shells[s1].alpha[0]; // Exponent of Shell 1

          for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
            //for(size_t s2(0), bf2_s(0); s2 < nShells; bf2_s+=n2, s2++) {

            n2 = basisSet_.shells[s2].size(); // Size of Shell 2
            expo2 = basisSet_.shells[s2].alpha[0]; // Exponent of Shell 2

            for(size_t s3(0), bf3_s(0); s3 < nShells; bf3_s+=n3, s3++) {
              //for(size_t s3(0), bf3_s(0); s3 < nShells; bf3_s+=n3, s3++) {

              n3 = basisSet_.shells[s3].size(); // Size of Shell 3
              expo3 = basisSet_.shells[s3].alpha[0]; // Exponent of Shell 3

              for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {
                //for(size_t s4(0), bf4_s(0); s4 < nShells; bf4_s+=n4, s4++, s1234++) {

                n4 = basisSet_.shells[s4].size(); // Size of Shell 4
                expo4 = basisSet_.shells[s4].alpha[0]; // Exponent of Shell 4

                // Round Robbin work distribution
#ifdef _OPENMP
                if( s1234 % nthreads != thread_id ) continue;
#endif

                shls[0] = int(s1);
                shls[1] = int(s2);
                shls[2] = int(s3);
                shls[3] = int(s4);

                if (basisSet_.forceCart) {
                  if(int2e_ipvip1ipvip2_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
                } else {
                  if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
                }

                ijkl = 0ul;
                auto nQuad = n1 * n2 * n3 * n4;
                for (l = 0ul, bf4 = bf4_s; l < n4; ++l, bf4++)
                  for (k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++)
                    for (j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
                      for (i = 0ul, bf1 = bf1_s; i < n1; ++i, bf1++) {

                        // (∇A∙∇B)(∇C∙∇D)(ij|kl)
                        auto dAdotdBdCdotdD =
                            buff[AxBxCxDx * nQuad + ijkl] + buff[AxBxCyDy * nQuad + ijkl] + buff[AxBxCzDz * nQuad + ijkl]
                            + buff[AyByCxDx * nQuad + ijkl] + buff[AyByCyDy * nQuad + ijkl] + buff[AyByCzDz * nQuad + ijkl]
                            + buff[AzBzCxDx * nQuad + ijkl] + buff[AzBzCyDy * nQuad + ijkl] + buff[AzBzCzDz * nQuad + ijkl];

                        // (∇Ax∇B)(∇C∙∇D)(ijkl)
                        auto dAcrossdB_xdCdotdD = buff[AyBzCxDx * nQuad + ijkl] - buff[AzByCxDx * nQuad + ijkl]
                                                  + buff[AyBzCyDy * nQuad + ijkl] - buff[AzByCyDy * nQuad + ijkl]
                                                  + buff[AyBzCzDz * nQuad + ijkl] - buff[AzByCzDz * nQuad + ijkl];

                        auto dAcrossdB_ydCdotdD = -buff[AxBzCxDx * nQuad + ijkl] + buff[AzBxCxDx * nQuad + ijkl]
                                                  - buff[AxBzCyDy * nQuad + ijkl] + buff[AzBxCyDy * nQuad + ijkl]
                                                  - buff[AxBzCzDz * nQuad + ijkl] + buff[AzBxCzDz * nQuad + ijkl];

                        auto dAcrossdB_zdCdotdD = buff[AxByCxDx * nQuad + ijkl] - buff[AyBxCxDx * nQuad + ijkl]
                                                  + buff[AxByCyDy * nQuad + ijkl] - buff[AyBxCyDy * nQuad + ijkl]
                                                  + buff[AxByCzDz * nQuad + ijkl] - buff[AyBxCzDz * nQuad + ijkl];

                        // (∇A∙∇B)(∇Cx∇D)(ijkl)
                        auto dAdotdBdCcrossdD_x = buff[AxBxCyDz * nQuad + ijkl] - buff[AxBxCzDy * nQuad + ijkl]
                                                  + buff[AyByCyDz * nQuad + ijkl] - buff[AyByCzDy * nQuad + ijkl]
                                                  + buff[AzBzCyDz * nQuad + ijkl] - buff[AzBzCzDy * nQuad + ijkl];

                        auto dAdotdBdCcrossdD_y = -buff[AxBxCxDz * nQuad + ijkl] + buff[AxBxCzDx * nQuad + ijkl]
                                                  - buff[AyByCxDz * nQuad + ijkl] + buff[AyByCzDx * nQuad + ijkl]
                                                  - buff[AzBzCxDz * nQuad + ijkl] + buff[AzBzCzDx * nQuad + ijkl];

                        auto dAdotdBdCcrossdD_z = buff[AxBxCxDy * nQuad + ijkl] - buff[AxBxCyDx * nQuad + ijkl]
                                                  + buff[AyByCxDy * nQuad + ijkl] - buff[AyByCyDx * nQuad + ijkl]
                                                  + buff[AzBzCxDy * nQuad + ijkl] - buff[AzBzCyDx * nQuad + ijkl];

                        // (∇Ax∇B)(∇Cx∇D)(ijkl)
                        auto dAcrossdB_xdCcrossdD_x = buff[AyBzCyDz * nQuad + ijkl] - buff[AzByCyDz * nQuad + ijkl]
                                                      - buff[AyBzCzDy * nQuad + ijkl] + buff[AzByCzDy * nQuad + ijkl];

                        auto dAcrossdB_xdCcrossdD_y = buff[AyBzCzDx * nQuad + ijkl] - buff[AzByCzDx * nQuad + ijkl]
                                                      - buff[AyBzCxDz * nQuad + ijkl] + buff[AzByCxDz * nQuad + ijkl];

                        auto dAcrossdB_xdCcrossdD_z = buff[AyBzCxDy * nQuad + ijkl] - buff[AzByCxDy * nQuad + ijkl]
                                                      - buff[AyBzCyDx * nQuad + ijkl] + buff[AzByCyDx * nQuad + ijkl];

                        auto dAcrossdB_ydCcrossdD_x = buff[AzBxCyDz * nQuad + ijkl] - buff[AxBzCyDz * nQuad + ijkl]
                                                      - buff[AzBxCzDy * nQuad + ijkl] + buff[AxBzCzDy * nQuad + ijkl];

                        auto dAcrossdB_ydCcrossdD_y = buff[AzBxCzDx * nQuad + ijkl] - buff[AxBzCzDx * nQuad + ijkl]
                                                      - buff[AzBxCxDz * nQuad + ijkl] + buff[AxBzCxDz * nQuad + ijkl];

                        auto dAcrossdB_ydCcrossdD_z = buff[AzBxCxDy * nQuad + ijkl] - buff[AxBzCxDy * nQuad + ijkl]
                                                      - buff[AzBxCyDx * nQuad + ijkl] + buff[AxBzCyDx * nQuad + ijkl];

                        auto dAcrossdB_zdCcrossdD_x = buff[AxByCyDz * nQuad + ijkl] - buff[AyBxCyDz * nQuad + ijkl]
                                                      - buff[AxByCzDy * nQuad + ijkl] + buff[AyBxCzDy * nQuad + ijkl];

                        auto dAcrossdB_zdCcrossdD_y = buff[AxByCzDx * nQuad + ijkl] - buff[AyBxCzDx * nQuad + ijkl]
                                                      - buff[AxByCxDz * nQuad + ijkl] + buff[AyBxCxDz * nQuad + ijkl];

                        auto dAcrossdB_zdCcrossdD_z = buff[AxByCxDy * nQuad + ijkl] - buff[AyBxCxDy * nQuad + ijkl]
                                                      - buff[AxByCyDx * nQuad + ijkl] + buff[AyBxCyDx * nQuad + ijkl];

                        //std::cout << std::scientific << std::setprecision(16);
                        //std::cout << "(" << bf1 << "," << bf2 << "|" << bf3 << "," << bf4 << ")  ";
                        //std::cout << dAcrossdB_xdCdotdD << std::endl;


                        auto IJKL = bf1 + bf2 * NB + bf3 * NB2 + bf4 * NB3;
                        auto IJLK = bf1 + bf2 * NB + bf4 * NB2 + bf3 * NB3;
                        auto JIKL = bf2 + bf1 * NB + bf3 * NB2 + bf4 * NB3;
                        auto JILK = bf2 + bf1 * NB + bf4 * NB2 + bf3 * NB3;
                        auto KLIJ = bf3 + bf4 * NB + bf1 * NB2 + bf2 * NB3;
                        auto LKIJ = bf4 + bf3 * NB + bf1 * NB2 + bf2 * NB3;
                        auto KLJI = bf3 + bf4 * NB + bf2 * NB2 + bf1 * NB3;
                        auto LKJI = bf4 + bf3 * NB + bf2 * NB2 + bf1 * NB3;

                        //auto KLIJ = bf1 + bf2*NB + bf3*NB2 + bf4*NB3;
                        //auto LKIJ = bf1 + bf2*NB + bf4*NB2 + bf3*NB3;
                        //auto KLJI = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
                        //auto LKJI = bf2 + bf1*NB + bf4*NB2 + bf3*NB3;
                        //auto IJKL = bf3 + bf4*NB + bf1*NB2 + bf2*NB3;
                        //auto IJLK = bf4 + bf3*NB + bf1*NB2 + bf2*NB3;
                        //auto JIKL = bf3 + bf4*NB + bf2*NB2 + bf1*NB3;
                        //auto JILK = bf4 + bf3*NB + bf2*NB2 + bf1*NB3;

//
//                        // (∇A∙∇B)(∇C∙∇D)(ij|kl)
//                        // (ij|kl)
//                        (*this)[nERIRef].pointer()[IJKL] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[IJLK] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[JIKL] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[JILK] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[KLIJ] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[LKIJ] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[KLJI] = dAdotdBdCdotdD;
//                        (*this)[nERIRef].pointer()[LKJI] = dAdotdBdCdotdD;
//
//                        // (∇Ax∇B)_x(∇C∙∇D)(ijkl)
//                        (*this)[nERIRef + 1].pointer()[IJKL] = dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 1].pointer()[IJLK] = dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 1].pointer()[JIKL] = -dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 1].pointer()[JILK] = -dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 1].pointer()[KLIJ] = dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 1].pointer()[LKIJ] = -dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 1].pointer()[KLJI] = dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 1].pointer()[LKJI] = -dAdotdBdCcrossdD_x;
//
//                        // (∇Ax∇B)_y(∇C∙∇D)(ijkl)
//                        (*this)[nERIRef + 2].pointer()[IJKL] = dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 2].pointer()[IJLK] = dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 2].pointer()[JIKL] = -dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 2].pointer()[JILK] = -dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 2].pointer()[KLIJ] = dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 2].pointer()[LKIJ] = -dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 2].pointer()[KLJI] = dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 2].pointer()[LKJI] = -dAdotdBdCcrossdD_y;
//
//                        // (∇Ax∇B)_z(∇C∙∇D)(ijkl)
//                        (*this)[nERIRef + 3].pointer()[IJKL] = dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 3].pointer()[IJLK] = dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 3].pointer()[JIKL] = -dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 3].pointer()[JILK] = -dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 3].pointer()[KLIJ] = dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 3].pointer()[LKIJ] = -dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 3].pointer()[KLJI] = dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 3].pointer()[LKJI] = -dAdotdBdCcrossdD_z;
//
//
//
//                        // (∇A∙∇B)(∇Cx∇D)_x(ijkl)
//                        (*this)[nERIRef + 4].pointer()[IJKL] = dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 4].pointer()[IJLK] = -dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 4].pointer()[JIKL] = dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 4].pointer()[JILK] = -dAdotdBdCcrossdD_x;
//                        (*this)[nERIRef + 4].pointer()[KLIJ] = dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 4].pointer()[LKIJ] = dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 4].pointer()[KLJI] = -dAcrossdB_xdCdotdD;
//                        (*this)[nERIRef + 4].pointer()[LKJI] = -dAcrossdB_xdCdotdD;
//
//                        // (∇A∙∇B)(∇Cx∇D)_y(ijkl)
//                        (*this)[nERIRef + 5].pointer()[IJKL] = dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 5].pointer()[IJLK] = -dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 5].pointer()[JIKL] = dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 5].pointer()[JILK] = -dAdotdBdCcrossdD_y;
//                        (*this)[nERIRef + 5].pointer()[KLIJ] = dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 5].pointer()[LKIJ] = dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 5].pointer()[KLJI] = -dAcrossdB_ydCdotdD;
//                        (*this)[nERIRef + 5].pointer()[LKJI] = -dAcrossdB_ydCdotdD;
//
//                        // (∇A∙∇B)(∇Cx∇D)_z(ijkl)
//                        (*this)[nERIRef + 6].pointer()[IJKL] = dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 6].pointer()[IJLK] = -dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 6].pointer()[JIKL] = dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 6].pointer()[JILK] = -dAdotdBdCcrossdD_z;
//                        (*this)[nERIRef + 6].pointer()[KLIJ] = dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 6].pointer()[LKIJ] = dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 6].pointer()[KLJI] = -dAcrossdB_zdCdotdD;
//                        (*this)[nERIRef + 6].pointer()[LKJI] = -dAcrossdB_zdCdotdD;
//
//
//
//                        // (∇Ax∇B)_x(∇Cx∇D)_x(ijkl)
//                        (*this)[nERIRef + 7].pointer()[IJKL] = dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[JILK] = dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_x;
//                        (*this)[nERIRef + 7].pointer()[LKJI] = dAcrossdB_xdCcrossdD_x;
//
//                        // (∇Ax∇B)_x(∇Cx∇D)_y(ijkl)
//                        (*this)[nERIRef + 8].pointer()[IJKL] = dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 8].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 8].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 8].pointer()[JILK] = dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 8].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 8].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 8].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 8].pointer()[LKJI] = dAcrossdB_ydCcrossdD_x;
//
//                        // (∇Ax∇B)_x(∇Cx∇D)_z(ijkl)
//                        (*this)[nERIRef + 9].pointer()[IJKL] = dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 9].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 9].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 9].pointer()[JILK] = dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 9].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 9].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 9].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 9].pointer()[LKJI] = dAcrossdB_zdCcrossdD_x;
//
//
//
//                        // (∇Ax∇B)_y(∇Cx∇D)_x(ijkl)
//                        (*this)[nERIRef + 10].pointer()[IJKL] = dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 10].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 10].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 10].pointer()[JILK] = dAcrossdB_ydCcrossdD_x;
//                        (*this)[nERIRef + 10].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 10].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 10].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_y;
//                        (*this)[nERIRef + 10].pointer()[LKJI] = dAcrossdB_xdCcrossdD_y;
//
//                        // (∇Ax∇B)_y(∇Cx∇D)_y(ijkl)
//                        (*this)[nERIRef + 11].pointer()[IJKL] = dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[JILK] = dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_y;
//                        (*this)[nERIRef + 11].pointer()[LKJI] = dAcrossdB_ydCcrossdD_y;
//
//                        // (∇Ax∇B)_y(∇Cx∇D)_z(ijkl)
//                        (*this)[nERIRef + 12].pointer()[IJKL] = dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 12].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 12].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 12].pointer()[JILK] = dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 12].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 12].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 12].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 12].pointer()[LKJI] = dAcrossdB_zdCcrossdD_y;
//
//
//
//                        // (∇Ax∇B)_z(∇Cx∇D)_x(ijkl)
//                        (*this)[nERIRef + 13].pointer()[IJKL] = dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 13].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 13].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 13].pointer()[JILK] = dAcrossdB_zdCcrossdD_x;
//                        (*this)[nERIRef + 13].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 13].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 13].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_z;
//                        (*this)[nERIRef + 13].pointer()[LKJI] = dAcrossdB_xdCcrossdD_z;
//
//                        // (∇Ax∇B)_z(∇Cx∇D)_y(ijkl)
//                        (*this)[nERIRef + 14].pointer()[IJKL] = dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 14].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 14].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 14].pointer()[JILK] = dAcrossdB_zdCcrossdD_y;
//                        (*this)[nERIRef + 14].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 14].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 14].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_z;
//                        (*this)[nERIRef + 14].pointer()[LKJI] = dAcrossdB_ydCcrossdD_z;
//
//                        // (∇Ax∇B)_z(∇Cx∇D)_z(ijkl)
//                        (*this)[nERIRef + 15].pointer()[IJKL] = dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[JILK] = dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_z;
//                        (*this)[nERIRef + 15].pointer()[LKJI] = dAcrossdB_zdCcrossdD_z;


                        // (ij|kl)
                        buff2[ijkl] = dAdotdBdCdotdD;
                        buff2[ijkl +     nQuad] = dAcrossdB_xdCdotdD;
                        buff2[ijkl + 2 * nQuad] = dAcrossdB_ydCdotdD;
                        buff2[ijkl + 3 * nQuad] = dAcrossdB_zdCdotdD;
                        buff2[ijkl + 4 * nQuad] = dAdotdBdCcrossdD_x;
                        buff2[ijkl + 5 * nQuad] = dAdotdBdCcrossdD_y;
                        buff2[ijkl + 6 * nQuad] = dAdotdBdCcrossdD_z;
                        buff2[ijkl + 7 * nQuad] = dAcrossdB_xdCcrossdD_x;
                        buff2[ijkl + 8 * nQuad] = dAcrossdB_xdCcrossdD_y;
                        buff2[ijkl + 9 * nQuad] = dAcrossdB_xdCcrossdD_z;
                        buff2[ijkl +10 * nQuad] = dAcrossdB_ydCcrossdD_x;
                        buff2[ijkl +11 * nQuad] = dAcrossdB_ydCcrossdD_y;
                        buff2[ijkl +12 * nQuad] = dAcrossdB_ydCcrossdD_z;
                        buff2[ijkl +13 * nQuad] = dAcrossdB_zdCcrossdD_x;
                        buff2[ijkl +14 * nQuad] = dAcrossdB_zdCcrossdD_y;
                        buff2[ijkl +15 * nQuad] = dAcrossdB_zdCcrossdD_z;
                        ijkl++;

                      }; // ijkl loop

                // compute RI reformated integrals
                ijkl = 0ul;
                for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
                  for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
                    for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
                      for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {
                        // (∇A∙∇B)(∇C∙∇D)(ij|kl)
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(0)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(0)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl], 1);
                        // (∇Ax∇B)(∇C∙∇D)(ijkl)
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(1)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(0)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +     nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(2)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(0)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 2 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(3)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(0)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 3 * nQuad], 1);
                        // (∇A∙∇B)(∇Cx∇D)(ijkl)
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(0)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(1)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 4 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(0)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(2)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 5 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(0)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(3)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 6 * nQuad], 1);
                        // (∇Ax∇B)(∇Cx∇D)(ijkl)
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(1)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(1)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 7 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(1)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(2)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 8 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(1)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(3)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl + 9 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(2)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(1)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +10 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(2)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(2)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +11 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(2)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(3)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +12 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(3)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(1)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +13 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(3)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(2)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +14 * nQuad], 1);
                        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                                   1, 1, NBRI, 1., riSSLL->asymm_term(3)->pointer() + NBRI * (bf1 + bf2 * NB), NBRI,
                                   riSSLL->asymm_term(3)->pointer() + NBRI * (bf3 + bf4 * NB), NBRI,
                                   0., &buff[ijkl +15 * nQuad], 1);
                        ijkl++;
                      }

                blas::axpy(16*nQuad,-1.,buff2,1,buff,1);

                double expoScalar = 16 * expo1 * expo2 * expo3 * expo4;
                for (size_t t = 0; t < 16; t++) {
                  // compute max element
                  auto minmax = std::minmax_element(buff2+t*nQuad, buff2+(t+1)*nQuad);
                  maxElems[thread_id + t*nthreads] = std::max(maxElems[thread_id + t*nthreads], std::max(-*minmax.first, *minmax.second)/expoScalar);
                  // compute mean absolute element in buff2
                  double sumAbsElem = std::accumulate(buff2+t*nQuad, buff2+(t+1)*nQuad, 0., [](double a, double b) { return a + std::abs(b); })/expoScalar;

                  // compute max error
                  minmax = std::minmax_element(buff+t*nQuad, buff+(t+1)*nQuad);
                  maxErrors[thread_id + t*nthreads] = std::max(maxErrors[thread_id + t*nthreads], std::max(-*minmax.first, *minmax.second)/expoScalar);
                  // compute mean absolute error in buff
                  double sumAbsError = std::accumulate(buff+t*nQuad, buff+(t+1)*nQuad, 0., [](double a, double b) { return a + std::abs(b); })/expoScalar;

                  // compute percent error
                  double sumAbsPercent = 0.0;
                  for (size_t i = 0; i < nQuad; i++)
                    if (std::abs(buff2[i+t*nQuad]) > largeValueThreshold * expoScalar) { // ignore small elements (relative error is meaningless
                      double percent = std::abs(buff[i+t*nQuad]) / std::abs(buff2[i+t*nQuad]);
                      sumAbsPercent += percent;
                      maxPercents[thread_id + t*nthreads] = std::max(maxPercents[thread_id + t*nthreads], percent);
                    }

                  int multiplier = 4;
                  if (s1 == s2) multiplier /= 2;
                  if (s3 == s4) multiplier /= 2;
                  sumAbsElems[thread_id + t*nthreads] += multiplier * sumAbsElem;
                  sumAbsErrors[thread_id + t*nthreads] += multiplier * sumAbsError;
                  sumAbsPercents[thread_id + t*nthreads] += multiplier * sumAbsPercent;
                }
              }; // s4
            }; // s3
          }; // s2
        }; // s1

      }; // omp region
      for (size_t t = 0; t < 16; t++) {
        double scale_denom = 4*SpeedOfLight()*SpeedOfLight();
        scale_denom *= scale_denom;
        double maxElem = *std::max_element(&maxElems[t*nthreads],&maxElems[(t+1)*nthreads]);
        std::cout << "SSSS-" << t << " element MAX = " << maxElem << "  scale by 1/(2mc)^2 = " << maxElem/scale_denom << std::endl;
        double maxError = *std::max_element(&maxErrors[t*nthreads],&maxErrors[(t+1)*nthreads]);
        std::cout << "SSSS-" << t << " error   MAX = " << maxError << "  scale by 1/(2mc)^2 = " << maxError/scale_denom << std::endl;
        double maxPercent = *std::max_element(&maxPercents[t*nthreads],&maxPercents[(t+1)*nthreads]);
        std::cout << "SSSS-" << t << " percent MAX = " << maxPercent << std::endl;
        double meanAbsElem = std::accumulate(&sumAbsElems[t*nthreads],&sumAbsElems[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSSS-" << t << " element MAE = " << meanAbsElem << "  scale by 1/(2mc)^2 = " << meanAbsElem/scale_denom << std::endl;
        double meanAbsError = std::accumulate(&sumAbsErrors[t*nthreads],&sumAbsErrors[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSSS-" << t << " error   MAE = " << meanAbsError << "  scale by 1/(2mc)^2 = " << meanAbsError/scale_denom << std::endl;
        double meanAbsPercent = std::accumulate(&sumAbsPercents[t*nthreads],&sumAbsPercents[(t+1)*nthreads], 0.0) / NB4;
        std::cout << "SSSS-" << t << " percent MAE = " << meanAbsPercent << std::endl;
      }


#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS
      auto durERIDCSSSS = tock(topERIDCSSSS);
      std::cout << "Libcint-ERI-Dirac-Coulomb-SSSS duration   = " << durERIDCSSSS << std::endl;
#endif

#if 0
      std::cout << std::scientific << std::setprecision(16);

      std::cout << "(∇A∙∇B)(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+1](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+2](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+3](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+4](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+5](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+6](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+7](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+8](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+9](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+10](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+11](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+12](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+13](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+14](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+15](i, j, k, l) << std::endl;
      };



#endif


    } // Dirac-Coulomb (SSSS) ∇_i∇_j∇_k∇_l(ij|kl)




    memManager_.free(cacheAll, buffAll, buffAll2, env, bas, atm);



  } // InCoreRelERI<double>::RI_direct_error

  template <>
  void InCoreRelERI<dcomplex>::RI_direct_error(BasisSet &originalBasisSet, Molecule &molecule_,
      EMPerturbation&, OPERATOR, const HamiltonianOptions &hamiltonianOptions, double largeValueThreshold) const {
    CErr("GIAO integral evaluation is NOT implemented in class InCoreRelERI.");
  } // InCoreRelERI<dcomplex>::RI_direct_error

}; // namespace ChronusQ
