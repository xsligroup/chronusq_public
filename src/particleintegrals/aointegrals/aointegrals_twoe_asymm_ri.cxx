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
#include <util/timer.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>

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
  void InCoreAsymmRITPI<double>::computeOneCholeskyPartialTPILibint(
      BasisSet &auxBasis, BasisSet &completeBasis) {
    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_ ? aux1_ : aux2_);
    size_t NBRI = aux.nRIBasis();
    
    auto beginSelectIndex = tick();
    // Group pivots into shells
    const std::vector<size_t> &pivots = aux.pivots();
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
    pivotIndicesByShell = InCoreCholeskyRIERI<double>::groupPivotsByShell(auxBasis, pivots);

    size_t pivotShellSize = pivotIndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivotShells;
    pivotShells.reserve(pivotShellSize);
    for (auto &shell_pivot : pivotIndicesByShell) {
      pivotShells.push_back(shell_pivot.first);
    }

    
    // Initialize libint engine
    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
                                 std::max(auxBasis.maxPrim, completeBasis.maxPrim),
                                 std::max(auxBasis.maxL, completeBasis.maxL),0);
    engines[0].set_precision(0.);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    // Compute TPI elements
    #pragma omp parallel
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      for (size_t I = 0, IPQ = 0; I < pivotShellSize; I++) {

        const auto &RSpair = pivotShells[I];
        auto &shell_pivot = pivotIndicesByShell[RSpair];

        size_t R = RSpair.first;
        size_t S = RSpair.second;

        size_t rBegin = auxBasis.mapSh2Bf[R];
        size_t sBegin = auxBasis.mapSh2Bf[S];
        size_t rSize = auxBasis.shells[R].size();
        size_t sSize = auxBasis.shells[S].size();
        size_t rEnd = rBegin + rSize;
        size_t sEnd = sBegin + sSize;

        for (size_t P(0), PQ(0); P < completeBasis.nShell; P++) {
          for (size_t Q = P; Q < completeBasis.nShell; Q++, PQ++, IPQ++) {

            // Round Robbin work distribution
            #ifdef _OPENMP
            if( IPQ % nthreads != thread_id ) continue;
            #endif

            // Evaluate ERI for shell quartet
            engines[thread_id].compute2<
            libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                auxBasis.shells[R],
                auxBasis.shells[S],
                completeBasis.shells[P],
                completeBasis.shells[Q]
                );
            const auto *buff =  buf_vec[0] ;
            if(buff == nullptr) continue;

            for (auto &pivot_index : shell_pivot) {

              auto rs = anaSquare(pivots[pivot_index], auxBasis.nBasis);
              size_t r = rs.first;
              size_t s = rs.second;

              if (P == Q) {
                for (size_t pBegin(completeBasis.mapSh2Bf[P]),
                    pSize(completeBasis.shells[P].size()),
                    pEnd(pBegin + pSize),
                    p(pBegin),
                    rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                    p < pEnd; p++, rsp++) {
                  for (size_t q(p),
                      rspq = rsp * pSize + p - pBegin;
                      q < pEnd; q++, rspq++) {
                    pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                  }
                }
              } else {
                for (size_t p(completeBasis.mapSh2Bf[P]),
                    pSize(completeBasis.shells[P].size()),
                    pEnd(p + pSize),
                    rsp = ((r-rBegin) * sSize + (s-sBegin)) * pSize;
                    p < pEnd; p++, rsp++) {
                  for (size_t q(completeBasis.mapSh2Bf[Q]),
                      qSize(completeBasis.shells[Q].size()),
                      qEnd(q + qSize),
                      rspq = rsp * qSize;
                      q < qEnd; q++, rspq++) {
                    pointer()[pivot_index + toCompound(p,q) * NBRI] = buff[rspq];
                  }
                }
              }
            }

          }; // Q
        }; // P

      }
    }; // omp region

    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;

    auto beginBuildPartialTPI = tick();
    size_t NB2   = completeBasis.nBasis*(completeBasis.nBasis+1)/2;
    size_t NB3   = NB2*NBRI;
    // S^{-1/2}(Q|ij)
    auto ijK = memManager().malloc<double>(NB3);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
               NBRI,NB2,NBRI,double(1.),aux.twoIndexERI()->pointer(),NBRI,
               pointer(),NBRI,double(0.),ijK,NBRI);

    for (size_t pq = 0; pq < NB2; pq++) {
      auto pqAna = anaCompound(pq);
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.first, pqAna.second, completeBasis.nBasis));
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.second, pqAna.first, completeBasis.nBasis));
    }

    memManager().free(ijK);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;

  } //InCoreAsymmRITPI<double>::computeOneCholeskyPartialTPILibint

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeOneCholeskyPartialTPILibint(BasisSet&, BasisSet&) {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeTwoCholeskyPartialTPILibint(BasisSet &basisSet1, BasisSet &basisSet2) {
    
    const InCoreCholeskyRIERI<double> &aux1 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_);
    size_t NBRI1 = aux1.nRIBasis();
    
    auto beginSelectIndex = tick();

    // Group pivots1 into shells
    const std::vector<size_t> &pivots1 = aux1.pivots();
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
    pivot1IndicesByShell = InCoreCholeskyRIERI<double>::groupPivotsByShell(basisSet1, pivots1);

    size_t pivot1ShellSize = pivot1IndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivot1Shells;
    pivot1Shells.reserve(pivot1ShellSize);
    for (auto &shell_pivot : pivot1IndicesByShell) {
      pivot1Shells.push_back(shell_pivot.first);
    }
    
    
    const InCoreCholeskyRIERI<double> &aux2 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
    size_t NBRI2 = aux2.nRIBasis();
    
    // Group pivots2 into shells
    const std::vector<size_t> &pivots2 = aux2.pivots();
    std::map<std::pair<size_t,size_t>, std::vector<size_t>>
    pivot2IndicesByShell = InCoreCholeskyRIERI<double>::groupPivotsByShell(basisSet2, pivots2);

    size_t pivot2ShellSize = pivot2IndicesByShell.size();
    std::vector<std::pair<size_t,size_t>> pivot2Shells;
    pivot2Shells.reserve(pivot2ShellSize);
    for (auto &shell_pivot : pivot2IndicesByShell) {
      pivot2Shells.push_back(shell_pivot.first);
    }
    

    // Initialize libint engine
    // Determine the number of OpenMP threads
    size_t nthreads = GetNumThreads();
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
                                 std::max(basisSet1.maxPrim, basisSet2.maxPrim),
                                 std::max(basisSet1.maxL, basisSet2.maxL),0);
    engines[0].set_precision(0.);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

    
    // Compute TPI elements
    #pragma omp parallel
    {
      size_t thread_id = GetThreadID();

      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();

      for (size_t I = 0, IJ = 0; I < pivot1ShellSize; I++) {

        const auto &RSpair = pivot1Shells[I];
        auto &shell_pivot1 = pivot1IndicesByShell[RSpair];

        size_t R = RSpair.first;
        size_t S = RSpair.second;

        size_t rBegin = basisSet1.mapSh2Bf[R];
        size_t sBegin = basisSet1.mapSh2Bf[S];
        size_t rSize = basisSet1.shells[R].size();
        size_t sSize = basisSet1.shells[S].size();
        size_t rEnd = rBegin + rSize;
        size_t sEnd = sBegin + sSize;

        for (size_t J = 0; J < pivot2ShellSize; J++, IJ++) {

          // Round Robbin work distribution
          #ifdef _OPENMP
          if( IJ % nthreads != thread_id ) continue;
          #endif

          const auto &PQpair = pivot2Shells[J];
          auto &shell_pivot2 = pivot2IndicesByShell[PQpair];

          size_t P = PQpair.first;
          size_t Q = PQpair.second;

          size_t pBegin = basisSet2.mapSh2Bf[P];
          size_t qBegin = basisSet2.mapSh2Bf[Q];
          size_t pSize = basisSet2.shells[P].size();
          size_t qSize = basisSet2.shells[Q].size();
          size_t pEnd = pBegin + pSize;
          size_t qEnd = qBegin + qSize;
          

          // Evaluate ERI for shell quartet
          engines[thread_id].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              basisSet1.shells[R],
              basisSet1.shells[S],
              basisSet2.shells[P],
              basisSet2.shells[Q]
              );
          const auto *buff =  buf_vec[0] ;
          if(buff == nullptr) continue;

          for (auto &pivot1_index : shell_pivot1) {

            auto rsPair = anaSquare(pivots1[pivot1_index], basisSet1.nBasis);
            size_t r = rsPair.first;
            size_t s = rsPair.second;
            size_t rs = ((r-rBegin) * sSize + (s-sBegin)) * pSize;

            for (auto &pivot2_index : shell_pivot2) {

              auto pqPair = anaSquare(pivots2[pivot2_index], basisSet2.nBasis);
              size_t p = pqPair.first;
              size_t q = pqPair.second;
              size_t rspq = (rs + (p-pBegin)) * qSize + (q-qBegin);
              
              pointer()[pivot1_index + pivot2_index * NBRI1] = buff[rspq];
              
            }
          }
          
        }

      }
    }; // omp region
    
    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;

    auto beginBuildPartialTPI = tick();
    double *SCR = memManager().template malloc<double>(NBRI1 * NBRI2);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               NBRI1, NBRI2, NBRI2, 1., pointer(), NBRI1, aux2.twoIndexERI()->pointer(), NBRI2, 0., SCR, NBRI1);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
               NBRI1, NBRI2, NBRI1, 1., aux1.twoIndexERI()->pointer(), NBRI1, SCR, NBRI1, 0., pointer(), NBRI1);
    memManager().free(SCR);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;
  } //InCoreAsymmRITPI<double>::computeTwoCholeskyPartialTPILibint

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeTwoCholeskyPartialTPILibint(BasisSet&, BasisSet&) {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  };

  template <>
  void InCoreAsymmRITPI<double>::computeOneCholeskyPartialTPIPrebuilt4Index(BasisSet &basisSet1, BasisSet &basisSet2){
    
    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_ ? aux1_ : aux2_);
    const std::vector<size_t>& pivots = aux.getPivots();

    size_t NBRI = aux.nRIBasis();
    size_t NB1 = basisSet1.nBasis;
    size_t NB2 = basisSet2.nBasis;
    size_t NB2_Squared = NB2 * NB2;

    // Created temporary object to select important elements
    double *ERI3I = aux.memManager().template malloc<double>(NB2_Squared * NBRI);

    auto beginSelectIndex = tick();

    #pragma omp parallel for
    for (size_t P = 0; P < NBRI; P++) {
      size_t pivot_index = pivots[P];
      for (size_t i = 0; i < NB2; i++) 
        for (size_t j = 0; j < NB2; j++) 
          ERI3I[P + (i+j*NB2) * NBRI] = aux1_ ?
          (*eri4I_)(pivot_index/NB1, pivot_index%NB1, i, j) : (*eri4I_)(i, j, pivot_index/NB1, pivot_index%NB1);
    }

    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;

    auto beginBuildPartialTPI = tick();

    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBRI,NB2_Squared,NBRI,1.,aux.twoIndexERI()->pointer(),NBRI,ERI3I,NBRI,0.,partialTPI_,NBRI);

    aux.memManager().free(ERI3I);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;
  }


  template <>
  void InCoreAsymmRITPI<dcomplex>::computeOneCholeskyPartialTPIPrebuilt4Index(BasisSet &basisSet1, BasisSet &basisSet2){
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeTwoCholeskyPartialTPIPrebuilt4Index(BasisSet &basisSet1, BasisSet &basisSet2){
    
    const InCoreCholeskyRIERI<double> &aux1 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_);
    const InCoreCholeskyRIERI<double> &aux2 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
    
    size_t NBRI1 = aux1.nRIBasis();
    size_t NBRI2 = aux2.nRIBasis();

    const std::vector<size_t>& pivots1 = aux1.getPivots();
    const std::vector<size_t>& pivots2 = aux2.getPivots();
    
    size_t NB1 = basisSet1.nBasis;
    size_t NB2 = basisSet2.nBasis;

    // Created temporary objects to select important elements
    double *NS = aux1.memManager().template malloc<double>(NBRI1 * NBRI2);
    double *SCR1 = aux1.memManager().template malloc<double>(NBRI1 * NBRI2);

    auto beginSelectIndex = tick();

    #pragma omp parallel for
    for (size_t P1 = 0; P1 < NBRI1; P1++) {
      size_t pivot_index1 = pivots1[P1];
      for (size_t P2 = 0; P2 < NBRI2; P2++) {
        size_t pivot_index2 = pivots2[P2];
        NS[P1 + P2* NBRI1] = (*eri4I_)(pivot_index1/NB1, pivot_index1%NB1, pivot_index2/NB2, pivot_index2%NB2); 
      }
    }

    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;

    auto beginBuildPartialTPI = tick();

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans, NBRI1, NBRI2, NBRI2, 1., NS, NBRI1, aux2.twoIndexERI()->pointer(), NBRI2, 0., SCR1, NBRI1);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans, NBRI1, NBRI2, NBRI1, 1., aux1.twoIndexERI()->pointer(), NBRI1, SCR1, NBRI1, 0., partialTPI_, NBRI1);

    aux1.memManager().free(NS, SCR1);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;
  }
  
  template <>
  void InCoreAsymmRITPI<dcomplex>::computeTwoCholeskyPartialTPIPrebuilt4Index(BasisSet &basisSet1, BasisSet &basisSet2){
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeAOInts(BasisSet &basisSet, BasisSet &basisSet2,
      Molecule& mol, EMPerturbation& emPert, OPERATOR, const HamiltonianOptions& options) {
    
    std::cout << "\nCalculating (ee|pp) Integrals Using Asymmetric Cholesky Decomposition: \n" << std::endl;
    
    if (build4I_ and not eri4I_) {
      std::cout << "     * Building full 4-index ERI for (ee|pp) per user's request" << std::endl;

      auto top4I = tick();
      eri4I_ = std::make_shared<InCore4indexTPI<double>>(this->memManager(),basisSet.nBasis,basisSet2.nBasis);
      eri4I_->computeAOInts(basisSet, basisSet2, mol, emPert, EP_ATTRACTION, options);
      auto dur4I = tock(top4I);
      std::cout << "       4-Index (ee|pp) evaluation duration   = " << dur4I << " s " << std::endl << std::endl;
    }

    if (aux1_) {
      std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux1 = 
          std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_); 
      if (aux2_) {
        std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux2 = 
            std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
        if (cd_aux1 and cd_aux2){
          if(cd_aux1->getPivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {-1., 1.};
            cd_aux1->computeAOInts(basisSet, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if(cd_aux2->getPivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {1., ProtMassPerE};
            cd_aux2->computeAOInts(basisSet2, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if (!partialTPI_) malloc();

          std::cout<< "     * Using elec and prot aux basis" << std::endl;
          auto topCDEP = tick();
          if (eri4I_) {
            std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
            computeTwoCholeskyPartialTPIPrebuilt4Index(basisSet, basisSet2);
          } else {
            std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
            computeTwoCholeskyPartialTPILibint(basisSet, basisSet2);
          }
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }
      } else {
        if (cd_aux1) {
          if(cd_aux1->getPivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {-1., 1.};
            cd_aux1->computeAOInts(basisSet, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if (!partialTPI_) malloc();
          
          std::cout<< "     * Using elec aux basis" << std::endl;
          auto topCDEP = tick();
          if (eri4I_){
            std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
            computeOneCholeskyPartialTPIPrebuilt4Index(basisSet, basisSet2);
          } else {
            std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
            computeOneCholeskyPartialTPILibint(basisSet, basisSet2);
          }
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }

      }

    } else {
      if (aux2_) {
        std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux2 = 
            std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
        if (cd_aux2) {
          if(cd_aux2->getPivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {1., ProtMassPerE};
            cd_aux2->computeAOInts(basisSet2, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if (!partialTPI_) malloc();

          std::cout<< "     * Using prot aux basis" << std::endl; 
          auto topCDEP = tick();
          if(eri4I_){
            std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
            computeOneCholeskyPartialTPIPrebuilt4Index(basisSet2, basisSet);
          }else{
            std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
            computeOneCholeskyPartialTPILibint(basisSet2, basisSet);
          }
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }

      } else {
        CErr("No aux available in IncoreAsymmRITPI::computeAOInts");
      }

    }
    
    if(printError_) printError(basisSet, basisSet2, mol, emPert);

    std::cout << std::endl << BannerEnd << std::endl;

  } // InCoreAsymmRITPI<double>::computeAOInts

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeAOInts(
      BasisSet&, BasisSet&, Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("GIAO integral evaluation is NOT implemented in class InCoreAsymmRITPI.");
  } // InCoreAsymmRITPI<dcomplex>::computeAOInts

}; // namespace ChronusQ
