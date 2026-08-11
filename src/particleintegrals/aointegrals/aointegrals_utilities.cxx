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
#include <lapack.hh>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>

#include <util/threads.hpp>
#include <chrono>
// Debug directives
//#define _DEBUGORTHO
//#define __DEBUGERI__

//#define _DEBUG_SCHWARZ_GRAD
//#define _VERIFY_SCHWARZ_GRAD_FD


namespace ChronusQ {



  /**
   *  \brief Allocate and evaluate the Schwarz bounds over the
   *  CGTO shell pairs.
   */ 
//  template <>
//  void AOIntegrals<dcomplex>::computeSchwarz() {
//    CErr("Only real GTOs are allowed",std::cout);
//  };
  template <typename IntsT>
  void DirectTPI<IntsT>::computeSchwarz() {

    if( schwarz() != nullptr ) CQMemManager::get().free(schwarz());
    if( schwarz2() != nullptr ) CQMemManager::get().free(schwarz2());

    // Allocate the schwarz tensor
    size_t nShell = basisSet().nShell;
    schwarz() = CQMemManager::get().malloc<double>(nShell*nShell);
    if (&basisSet() != &basisSet2())
      schwarz2() = CQMemManager::get().malloc<double>(basisSet2().nShell*basisSet2().nShell);

    // Define the libint2 integral engine
    libint2::Engine engine(libintOperator(),
      std::max(basisSet().maxPrim, basisSet2().maxPrim),
      std::max(basisSet().maxL,basisSet2().maxL), 0);

    if (this->kernel() == Kernel::ShortRangeErfc)
      engine.set_params(this->rangeSeparationParameter());

    engine.set_precision(0.); // Don't screen prims during evaluation

    const auto &buf_vec = engine.results();

    auto topSch = std::chrono::high_resolution_clock::now();
  
    size_t n1,n2;
    for(auto s1(0ul); s1 < basisSet().nShell; s1++) {
      n1 = basisSet().shells[s1].size(); // Size shell 1
    for(auto s2(0ul); s2 <= s1; s2++) {
      n2 = basisSet().shells[s2].size(); // Size shell 2



      // Evaluate the shell quartet (s1 s2 | s1 s2)
      engine.compute(
        basisSet().shells[s1],
        basisSet().shells[s2],
        basisSet().shells[s1],
        basisSet().shells[s2]
      );

      if(buf_vec[0] == nullptr) continue;

      // Allocate space to hold the diagonals
      double* diags = CQMemManager::get().malloc<double>(n1*n2);

      for(auto i(0), ij(0); i < n1; i++)
      for(auto j(0); j < n2; j++, ij++)
        diags[i + j*n1] = buf_vec[0][ij*n1*n2 + ij];


      schwarz()[s1 + s2*basisSet().nShell] =
        std::sqrt(lapack::lange(lapack::Norm::Inf,n1,n2,diags,n1));

      // Free up space
      CQMemManager::get().free(diags);

    } // loop s2
    } // loop s1

    if (&basisSet() != &basisSet2()) {
      // compute (rs|rs)
      for(auto s1(0ul); s1 < basisSet2().nShell; s1++) {
        n1 = basisSet2().shells[s1].size(); // Size shell 1
      for(auto s2(0ul); s2 <= s1; s2++) {
        n2 = basisSet2().shells[s2].size(); // Size shell 2



        // Evaluate the shell quartet (s1 s2 | s1 s2)
        engine.compute(
          basisSet2().shells[s1],
          basisSet2().shells[s2],
          basisSet2().shells[s1],
          basisSet2().shells[s2]
        );

        if(buf_vec[0] == nullptr) continue;

        // Allocate space to hold the diagonals
        double* diags = CQMemManager::get().malloc<double>(n1*n2);

        for(auto i(0), ij(0); i < n1; i++)
        for(auto j(0); j < n2; j++, ij++)
          diags[i + j*n1] = buf_vec[0][ij*n1*n2 + ij];


        schwarz2()[s1 + s2*basisSet2().nShell] =
          std::sqrt(lapack::lange(lapack::Norm::Inf,n1,n2,diags,n1));

        // Free up space
        CQMemManager::get().free(diags);

      } // loop s2
      } // loop s1
      // done computing (rs|rs)
    }

    auto botSch = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> durSch = botSch - topSch;

    HerMat('L',basisSet().nShell,schwarz(),basisSet().nShell);
    if (&basisSet() != &basisSet2())
      HerMat('L',basisSet2().nShell,schwarz2(),basisSet2().nShell);

#if 0
    prettyPrintSmart(std::cout,"Schwarz",schwarz,basisSet_.nShell,
      basisSet_.nShell,basisSet_.nShell);
#endif

  }; // DirectERI<double>::computeSchwarz
  template void DirectTPI<double>::computeSchwarz();
  template void DirectTPI<dcomplex>::computeSchwarz();

// =====================================================================
// computeSchwarz + computeSchwarzGrad -- FD-based rigorous R
//
// R_NM = max_d [ Q_{∂_d s_N · s_M} + Q_{s_N · ∂_d s_M} ]
//
// where each Q is computed by central differences on the diagonal ERI
// of the corresponding (displaced) shell pair, and the underlying
// Schwarz factor uses the element-wise max-element norm of the
// diagonal (not the L∞-induced row-sum norm).
//
// Compile-time flags:
//   _DEBUG_SCHWARZ_GRAD          -- R/Q sanity, per-(L1,L2) breakdown
//   _VERIFY_SCHWARZ_GRAD_FD      -- self-consistency check (production
//                                   R is FD R, so this should now show
//                                   R/R_FD == 1.0 to roundoff)
// =====================================================================

template <typename IntsT>
void DirectTPI<IntsT>::computeSchwarzGrad() {

  if (schwarzGrad()  != nullptr) CQMemManager::get().free(schwarzGrad());
  if (schwarzGrad2() != nullptr) CQMemManager::get().free(schwarzGrad2());

  // Ensure ordinary Schwarz exists (we use it for R/Q sanity stats only)
  if (schwarz() == nullptr) computeSchwarz();

  auto topT = std::chrono::high_resolution_clock::now();

  const auto eriOperator = libintOperator();
  const bool isShortRangeErfc = this->kernel() == Kernel::ShortRangeErfc;
  const double omega = this->rangeSeparationParameter();

  auto computeR_oneBasis = [eriOperator, isShortRangeErfc, omega]
    (BasisSet& bs, double*& R_out) {

    const size_t nShell = bs.nShell;
    R_out = CQMemManager::get().malloc<double>(nShell * nShell);
    std::memset(R_out, 0, nShell * nShell * sizeof(double));

    // FD displacement. h=1e-4 gives ~1e-8 truncation in R^2, ~1e-4
    // relative in R — well below screening thresholds.
    const double h = 1e-4;

    // Per-thread engines.
    const size_t nThreads = GetNumThreads();
    std::vector<libint2::Engine> engines(nThreads);
    engines[0] = libint2::Engine(eriOperator, bs.maxPrim, bs.maxL, 0);
    if (isShortRangeErfc)
      engines[0].set_params(omega);
    engines[0].set_precision(0.);
    for (size_t t = 1; t < nThreads; t++) engines[t] = engines[0];

    // Element-wise max of the diagonal (the *correct* shell-level Schwarz
    // quantity). NOT the L∞-induced row-sum norm.
    auto diag_max = [](const double* buf, size_t na, size_t nb) -> double {
      const size_t nbra = na * nb;
      double mx = 0.0;
      for (size_t i = 0; i < na; i++)
        for (size_t j = 0; j < nb; j++) {
          const size_t ij = i*nb + j;
          mx = std::max(mx, std::abs(buf[ij * (nbra + 1)]));
        }
      return mx;
    };

    auto displaced = [&](size_t s, int d, double sign) {
      libint2::Shell sh = bs.shells[s];
      sh.O[d] += sign * h;
      return sh;
    };

    // Q_{∂_d s1 · s2} via central differences on s1.
    auto Q_fd_left = [&](libint2::Engine& eng,
                         size_t s1, size_t s2, int d) -> double {
      auto s1p = displaced(s1, d, +1.0);
      auto s1m = displaced(s1, d, -1.0);
      const auto& m = bs.shells[s2];
      const size_t na = bs.shells[s1].size();
      const size_t nb = m.size();
      const size_t nelem = na*nb*na*nb;
      const auto& buf = eng.results();

      std::vector<double> A(nelem), B(nelem), C(nelem);
      eng.compute(s1p, m, s1p, m); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, A.begin());
      eng.compute(s1p, m, s1m, m); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, B.begin());
      eng.compute(s1m, m, s1m, m); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, C.begin());

      std::vector<double> comb(nelem);
      const double inv4h2 = 1.0 / (4.0 * h * h);
      for (size_t k = 0; k < nelem; k++)
        comb[k] = (A[k] - 2.0*B[k] + C[k]) * inv4h2;
      return std::sqrt(diag_max(comb.data(), na, nb));
    };

    // Q_{s1 · ∂_d s2} via central differences on s2.
    auto Q_fd_right = [&](libint2::Engine& eng,
                          size_t s1, size_t s2, int d) -> double {
      auto s2p = displaced(s2, d, +1.0);
      auto s2m = displaced(s2, d, -1.0);
      const auto& n = bs.shells[s1];
      const size_t na = n.size();
      const size_t nb = bs.shells[s2].size();
      const size_t nelem = na*nb*na*nb;
      const auto& buf = eng.results();

      std::vector<double> A(nelem), B(nelem), C(nelem);
      eng.compute(n, s2p, n, s2p); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, A.begin());
      eng.compute(n, s2p, n, s2m); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, B.begin());
      eng.compute(n, s2m, n, s2m); if (!buf[0]) return 0.0;
      std::copy_n(buf[0], nelem, C.begin());

      std::vector<double> comb(nelem);
      const double inv4h2 = 1.0 / (4.0 * h * h);
      for (size_t k = 0; k < nelem; k++)
        comb[k] = (A[k] - 2.0*B[k] + C[k]) * inv4h2;
      return std::sqrt(diag_max(comb.data(), na, nb));
    };

    // R is symmetric in (s1,s2). Compute upper triangle, mirror.
    #pragma omp parallel for schedule(dynamic)
    for (size_t s1 = 0; s1 < nShell; s1++) {
      const size_t tid = GetThreadID();
      auto& eng = engines[tid];
      for (size_t s2 = 0; s2 <= s1; s2++) {
        double R = 0.0;
        for (int d = 0; d < 3; d++) {
          const double Qa = Q_fd_left (eng, s1, s2, d);
          const double Qb = Q_fd_right(eng, s1, s2, d);
          R = std::max(R, Qa + Qb);
        }
        R_out[s1 + s2*nShell] = R;
        R_out[s2 + s1*nShell] = R;
      }
    }
  };

  computeR_oneBasis(basisSet(), schwarzGrad());
  if (&basisSet() != &basisSet2())
    computeR_oneBasis(basisSet2(), schwarzGrad2());

  auto botT = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> durT = botT - topT;

#ifdef _DEBUG_SCHWARZ_GRAD
  {
    const size_t nShell = basisSet().nShell;
    std::cout << "\n=== computeSchwarzGrad sanity (FD-based) ===\n";
    std::cout << "Shells:  " << nShell << "\n";
    std::cout << "Setup:   " << std::fixed << std::setprecision(4)
              << durT.count() << " s\n";

    double minRatio = 1e30, maxRatio = 0.0, sumRatio = 0.0;
    size_t nPairs = 0;
    constexpr int LMAX = 5;
    std::array<std::array<double, LMAX+1>, LMAX+1> ratioMin{}, ratioMax{}, ratioSum{};
    std::array<std::array<size_t, LMAX+1>, LMAX+1> ratioCnt{};
    for (int i = 0; i <= LMAX; i++)
      for (int j = 0; j <= LMAX; j++) { ratioMin[i][j] = 1e30; ratioMax[i][j] = 0.0; }

    for (size_t s1 = 0; s1 < nShell; s1++)
    for (size_t s2 = 0; s2 <= s1; s2++) {
      const double Q = schwarz()[s1 + s2*nShell];
      const double R = schwarzGrad()[s1 + s2*nShell];
      if (Q > 1e-15) {
        const double r = R/Q;
        minRatio = std::min(minRatio, r);
        maxRatio = std::max(maxRatio, r);
        sumRatio += r;
        nPairs++;
        int L1 = basisSet().shells[s1].contr[0].l;
        int L2 = basisSet().shells[s2].contr[0].l;
        if (L1 <= LMAX && L2 <= LMAX) {
          int a = std::min(L1, L2), b = std::max(L1, L2);
          ratioMin[a][b] = std::min(ratioMin[a][b], r);
          ratioMax[a][b] = std::max(ratioMax[a][b], r);
          ratioSum[a][b] += r;
          ratioCnt[a][b]++;
        }
      }
    }
    std::cout << "\nR/Q ratio stats (over " << nPairs << " unique pairs):\n";
    std::cout << "  min:  " << minRatio << "\n";
    std::cout << "  max:  " << maxRatio << "\n";
    std::cout << "  mean: " << sumRatio/nPairs << "\n";

    std::cout << "\nR/Q ratio by (L_min,L_max):\n";
    std::cout << "  (L1 L2)   count       min          mean         max\n";
    for (int a = 0; a <= LMAX; a++)
    for (int b = a; b <= LMAX; b++) {
      if (ratioCnt[a][b] == 0) continue;
      std::cout << "  (" << a << " " << b << ")   "
                << std::setw(7) << ratioCnt[a][b] << "   "
                << std::setw(11) << ratioMin[a][b] << "   "
                << std::setw(11) << ratioSum[a][b]/ratioCnt[a][b] << "   "
                << std::setw(11) << ratioMax[a][b] << "\n";
    }

    double maxAsym = 0.0;
    for (size_t s1 = 0; s1 < nShell; s1++)
    for (size_t s2 = 0; s2 < nShell; s2++) {
      const double r12 = schwarzGrad()[s1 + s2*nShell];
      const double r21 = schwarzGrad()[s2 + s1*nShell];
      maxAsym = std::max(maxAsym, std::abs(r12 - r21));
    }
    std::cout << "\nMax |R(NM) - R(MN)|: "
              << std::scientific << maxAsym << "\n";

    size_t nRzero = 0;
    for (size_t k = 0; k < nShell*nShell; k++)
      if (schwarzGrad()[k] == 0.0) nRzero++;
    std::cout << "R exact zeros: " << nRzero << " / " << nShell*nShell;
    if (nRzero > 0) std::cout << "   *** suspicious -- R should be nonzero for all pairs ***";
    std::cout << "\n================================\n" << std::endl;
  }
#endif

#ifdef _VERIFY_SCHWARZ_GRAD_FD
  // Self-consistency: production R is FD R, so R/R_FD should == 1.0
  // to roundoff for every sampled pair.
  {
    const size_t nShell = basisSet().nShell;
    const double h = 1e-4;
    const double tol = 1e-6;

    const size_t s1_stride = std::max<size_t>(1, nShell / 60);

    libint2::Engine eng_fd(libintOperator(),
                           basisSet().maxPrim, basisSet().maxL, 0);
    if (this->kernel() == Kernel::ShortRangeErfc)
      eng_fd.set_params(this->rangeSeparationParameter());
    eng_fd.set_precision(0.);
    const auto& buf_fd = eng_fd.results();

    auto displaced = [&](size_t s, int d, double sign) {
      libint2::Shell sh = basisSet().shells[s];
      sh.O[d] += sign * h;
      return sh;
    };
    auto diag_max = [&](const double* buf, size_t na, size_t nb) -> double {
      const size_t nbra = na * nb;
      double mx = 0.0;
      for (size_t i = 0; i < na; i++)
        for (size_t j = 0; j < nb; j++) {
          const size_t ij = i*nb + j;
          mx = std::max(mx, std::abs(buf[ij*(nbra+1)]));
        }
      return mx;
    };
    auto Q_fd_left = [&](size_t s1, size_t s2, int d) -> double {
      auto s1p = displaced(s1, d, +1.0);
      auto s1m = displaced(s1, d, -1.0);
      const auto& m = basisSet().shells[s2];
      const size_t na = basisSet().shells[s1].size();
      const size_t nb = m.size();
      std::vector<double> bufA(na*nb*na*nb), bufB(na*nb*na*nb), bufC(na*nb*na*nb);
      eng_fd.compute(s1p, m, s1p, m); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufA.begin());
      eng_fd.compute(s1p, m, s1m, m); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufB.begin());
      eng_fd.compute(s1m, m, s1m, m); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufC.begin());
      std::vector<double> comb(na*nb*na*nb);
      for (size_t k = 0; k < na*nb*na*nb; k++)
        comb[k] = (bufA[k] - 2.0*bufB[k] + bufC[k]) / (4.0*h*h);
      return std::sqrt(diag_max(comb.data(), na, nb));
    };
    auto Q_fd_right = [&](size_t s1, size_t s2, int d) -> double {
      auto s2p = displaced(s2, d, +1.0);
      auto s2m = displaced(s2, d, -1.0);
      const auto& n = basisSet().shells[s1];
      const size_t na = n.size();
      const size_t nb = basisSet().shells[s2].size();
      std::vector<double> bufA(na*nb*na*nb), bufB(na*nb*na*nb), bufC(na*nb*na*nb);
      eng_fd.compute(n, s2p, n, s2p); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufA.begin());
      eng_fd.compute(n, s2p, n, s2m); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufB.begin());
      eng_fd.compute(n, s2m, n, s2m); if (buf_fd[0] == nullptr) return 0.0;
      std::copy_n(buf_fd[0], na*nb*na*nb, bufC.begin());
      std::vector<double> comb(na*nb*na*nb);
      for (size_t k = 0; k < na*nb*na*nb; k++)
        comb[k] = (bufA[k] - 2.0*bufB[k] + bufC[k]) / (4.0*h*h);
      return std::sqrt(diag_max(comb.data(), na, nb));
    };

    double minRatioFD = 1e30;
    size_t worst_s1=0, worst_s2=0;
    int    worst_L1=0, worst_L2=0;
    double worst_R = 0.0, worst_R_FD = 0.0;
    size_t nViolations = 0;
    size_t nSampled = 0;

    std::cout << "\n=== _VERIFY_SCHWARZ_GRAD_FD ===\n";
    std::cout << "h = " << h << "  rel tol = " << tol << "\n";

    for (size_t s1 = 0; s1 < nShell; s1 += s1_stride) {
      for (size_t s2 = 0; s2 < nShell; s2++) {
        double R_FD = 0.0;
        for (int d = 0; d < 3; d++) {
          double Qa = Q_fd_left (s1, s2, d);
          double Qb = Q_fd_right(s1, s2, d);
          R_FD = std::max(R_FD, Qa + Qb);
        }
        const double R_an = schwarzGrad()[s1 + s2*nShell];
        nSampled++;
        if (R_FD > 0.0) {
          double ratio = R_an / R_FD;
          if (ratio < minRatioFD) {
            minRatioFD = ratio;
            worst_s1 = s1; worst_s2 = s2;
            worst_L1 = basisSet().shells[s1].contr[0].l;
            worst_L2 = basisSet().shells[s2].contr[0].l;
            worst_R = R_an; worst_R_FD = R_FD;
          }
          if (ratio < 1.0 - tol) nViolations++;
        }
      }
    }
    std::cout << "Sampled pairs:     " << nSampled << "\n";
    std::cout << "FD-rigor violations (R/R_FD < " << (1.0-tol) << "):  " << nViolations;
    if (nViolations > 0) std::cout << "   *** R IS NOT RIGOROUS ***";
    std::cout << "\n";
    std::cout << "min R/R_FD ratio:  " << std::scientific << minRatioFD << "\n";
    std::cout << "worst pair:        s1=" << worst_s1 << "  s2=" << worst_s2
              << "  (L1=" << worst_L1 << ", L2=" << worst_L2 << ")\n";
    std::cout << "  R_analytic = " << worst_R << "   R_FD = " << worst_R_FD << "\n";
    std::cout << "================================\n" << std::endl;
  }
#endif

};

template void DirectTPI<double>::computeSchwarzGrad();
template void DirectTPI<dcomplex>::computeSchwarzGrad();

}; // namespace ChronusQ

