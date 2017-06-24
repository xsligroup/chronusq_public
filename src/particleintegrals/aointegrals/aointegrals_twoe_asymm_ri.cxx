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
  void InCoreAsymmRITPI<double>::computeOneCholeskyRawSubTPILibint(
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
  } //InCoreAsymmRITPI<double>::computeOneCholeskyRawSubTPILibint

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeOneCholeskyRawSubTPILibint(BasisSet&, BasisSet&) {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeOneCholeskyPartialTPI() {
    auto beginBuildPartialTPI = tick();

    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_ ? aux1_ : aux2_);
    size_t NBRI = aux.nRIBasis();
    size_t NBcomplete = asymmCDalg_ == ASYMM_CD_ALG::INT1_AUX ? sNB : NB;
    size_t NB2   = NBcomplete*(NBcomplete+1)/2;
    size_t NB3   = NB2*NBRI;
    // S^{-1/2}(Q|ij)
    auto ijK = CQMemManager::get().malloc<double>(NB3);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
               NBRI,NB2,NBRI,double(1.),aux.twoIndexERI()->pointer(),NBRI,
               pointer(),NBRI,double(0.),ijK,NBRI);

    for (size_t pq = 0; pq < NB2; pq++) {
      auto pqAna = anaCompound(pq);
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.first, pqAna.second, NBcomplete));
      std::copy(&ijK[pq*NBRI], &ijK[pq*NBRI+NBRI], pointer()+NBRI*toSquare(pqAna.second, pqAna.first, NBcomplete));
    }

    CQMemManager::get().free(ijK);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;

  } //InCoreAsymmRITPI<double>::computeOneCholeskyPartialTPI

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeOneCholeskyPartialTPI() {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeTwoCholeskyRawSubTPILibint(BasisSet &basisSet1, BasisSet &basisSet2) {
    
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
  } //InCoreAsymmRITPI<double>::computeTwoCholeskyRawSubTPILibint

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeTwoCholeskyRawSubTPILibint(BasisSet&, BasisSet&) {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  };

  template <>
  void InCoreAsymmRITPI<double>::computeTwoCholeskyPartialTPI(){
    auto beginBuildPartialTPI = tick();

    const InCoreCholeskyRIERI<double> &aux1 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_);
    size_t NBRI1 = aux1.nRIBasis();
    const InCoreCholeskyRIERI<double> &aux2 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
    size_t NBRI2 = aux2.nRIBasis();

    double *SCR = CQMemManager::get().template malloc<double>(NBRI1 * NBRI2);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               NBRI1, NBRI2, NBRI2, 1., pointer(), NBRI1, aux2.twoIndexERI()->pointer(), NBRI2, 0., SCR, NBRI1);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
               NBRI1, NBRI2, NBRI1, 1., aux1.twoIndexERI()->pointer(), NBRI1, SCR, NBRI1, 0., pointer(), NBRI1);
    CQMemManager::get().free(SCR);

    double durBuildPartialTPI = tock(beginBuildPartialTPI);
    std::cout<< "  Cholesky-Asymm-Build-PartialTPI duration = " << durBuildPartialTPI << " s " << std::endl;
  } //InCoreAsymmRITPI<double>::computeTwoCholeskyPartialTPI

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeTwoCholeskyPartialTPI() {
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  };

  template <>
  void InCoreAsymmRITPI<double>::computeOneCholeskyRawSubTPIPrebuilt4Index(){
    
    const InCoreCholeskyRIERI<double> &aux =
        *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_ ? aux1_ : aux2_);
    const std::vector<size_t>& pivots = aux.pivots();

    size_t NBRI = aux.nRIBasis();
    size_t NB1 = NB;
    size_t NB2 = sNB;
    size_t NB2_Squared = NB2 * NB2;
    // Use compound index to select elements from pre-built 4-index
    size_t NB2_NTT = NB2*(NB2+1)/2;

    // Select Compound Index and save in pointer()
    auto beginSelectIndex = tick();
    if (aux1_){
      #pragma omp parallel for
      for (size_t P = 0; P < NBRI; P++) {
        size_t pivot_index = pivots[P];
        for (size_t pq = 0; pq < NB2_NTT; pq++) {
          auto pqAna = anaCompound(pq);
          pointer()[P + pq*NBRI] = (*eri4I_)(pivot_index/NB1, pivot_index%NB1, pqAna.first, pqAna.second);
        }
      }  
    }else{
      #pragma omp parallel for
      for (size_t P = 0; P < NBRI; P++) {
        size_t pivot_index = pivots[P];
        for (size_t pq = 0; pq < NB2_NTT; pq++) {
          auto pqAna = anaCompound(pq);
          pointer()[P + pq*NBRI] = (*eri4I_)(pqAna.first, pqAna.second, pivot_index/NB1, pivot_index%NB1);
        }
      }
    }
    
    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;
  }


  template <>
  void InCoreAsymmRITPI<dcomplex>::computeOneCholeskyRawSubTPIPrebuilt4Index(){
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeTwoCholeskyRawSubTPIPrebuilt4Index(){
    
    const InCoreCholeskyRIERI<double> &aux1 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_);
    const InCoreCholeskyRIERI<double> &aux2 = *std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
    
    size_t NBRI1 = aux1.nRIBasis();
    size_t NBRI2 = aux2.nRIBasis();

    const std::vector<size_t>& pivots1 = aux1.pivots();
    const std::vector<size_t>& pivots2 = aux2.pivots();
    
    size_t NB1 = NB;
    size_t NB2 = sNB;

    auto beginSelectIndex = tick();

    #pragma omp parallel for
    for (size_t P1 = 0; P1 < NBRI1; P1++) {
      size_t pivot_index1 = pivots1[P1];
      for (size_t P2 = 0; P2 < NBRI2; P2++) {
        size_t pivot_index2 = pivots2[P2];
        pointer()[P1 + P2* NBRI1] = (*eri4I_)(pivot_index1/NB1, pivot_index1%NB1, pivot_index2/NB2, pivot_index2%NB2);
      }
    }

    double durSelectIndex = tock(beginSelectIndex);
    std::cout<< "  Cholesky-Asymm-Select-PartialTPI-Index duration = " << durSelectIndex << " s " << std::endl;
  }
  
  template <>
  void InCoreAsymmRITPI<dcomplex>::computeTwoCholeskyRawSubTPIPrebuilt4Index(){
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::prebuilt4Index(BasisSet &basisSet, BasisSet &basisSet2,
      Molecule& mol, EMPerturbation& emPert, OPERATOR, const HamiltonianOptions& options) {

    if (build4I_ and not eri4I_) {
      std::cout << "     * Building full 4-index ERI for (ee|pp) per user's request" << std::endl;

      auto top4I = tick();
      eri4I_ = std::make_shared<InCore4indexTPI<double>>(basisSet.nBasis,basisSet2.nBasis);
      eri4I_->computeAOInts(basisSet, basisSet2, mol, emPert, EP_ATTRACTION, options);
      auto dur4I = tock(top4I);
      std::cout << "       4-Index (ee|pp) evaluation duration   = " << dur4I << " s " << std::endl << std::endl;
    }

  } // InCoreAsymmRITPI<double>::prebuilt4Index

  template <>
  void InCoreAsymmRITPI<dcomplex>::prebuilt4Index(BasisSet&, BasisSet&,
      Molecule& mol, EMPerturbation& emPert, OPERATOR, const HamiltonianOptions& options){
    CErr("Complex TPI for InCoreAsymmRITPI NYI",std::cout);
  }

  template <>
  void InCoreAsymmRITPI<double>::computeAOInts(BasisSet &basisSet, BasisSet &basisSet2,
      Molecule& mol, EMPerturbation& emPert, OPERATOR, const HamiltonianOptions& options) {
    
    std::cout << "\nCalculating (ee|pp) Integrals Using Asymmetric Cholesky Decomposition: \n" << std::endl;

    prebuilt4Index(basisSet, basisSet2, mol, emPert, EP_ATTRACTION, options);

    if (aux1_) {
      std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux1 = 
          std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1_); 
      if (aux2_) {
        std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux2 = 
            std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
        if (asymmCDalg_ != ASYMM_CD_ALG::CONNECTOR 
            and asymmCDalg_ != ASYMM_CD_ALG::COMBINEAUXBASIS
            and asymmCDalg_ != ASYMM_CD_ALG::COMBINEMATRIX)
          CErr("Auxiliary bases are set for both integrals, "
               "but requested algorithm is only one-side (INT1_AUX or INT2_AUX)");

        if (cd_aux1 and cd_aux2){

          // If cd_aux1 and/or cd_aux2 is created but has not been computed, 
          // (which cou be the case if the user forces to use two-aux asymm algorithms, but we dont have existing aux basis/bases)
          // Compute them first 
          if(cd_aux1->pivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {-1., 1.};
            cd_aux1->computeAOInts(basisSet, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }
          if(cd_aux2->pivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {1., ProtMassPerE};
            cd_aux2->computeAOInts(basisSet2, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          std::cout<< "     * Using elec and prot aux basis" << std::endl;
          auto topCDEP = tick();

          /*************************************************************
            Connector Method (Double One-Component RI in the paper)

            // partialTPI_ to be built:
              A_{\alpha \Gamma} = \sum_{\beta \in B_e \sum_{\Theta \in B_n}
                                  (K^{-1})_{\alpha \beta}      (\beta \vert \Theta)      K^{-1})__{\Theta \Gamma}

          *************************************************************/
          if (asymmCDalg_ == ASYMM_CD_ALG::CONNECTOR) {
            if (!partialTPI_) malloc();
            // Compute raw two-index TPI (\beta \vert \Theta) 
            if (eri4I_) {
              std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
              computeTwoCholeskyRawSubTPIPrebuilt4Index();
            } else {
              std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
              computeTwoCholeskyRawSubTPILibint(basisSet, basisSet2);
            }
            // Compute full A_{\alpha \Gamma} = (K^{-1})_{\alpha \beta}  (\beta \vert \Theta)  K^{-1})__{\Theta \Gamma}
            computeTwoCholeskyPartialTPI();



          } else if (asymmCDalg_ == ASYMM_CD_ALG::COMBINEAUXBASIS) {
          /*************************************************************
            CombineAuxBasis Method (Two-Component RI in the paper)

            // No partialTPI_ to be built, but changing aux1_ and aux2_
                aux1_: L_{p q, \kappa} &= \sum_{\lambda \in B_c} ( p q \vert \lambda ) ( K^{-T}_{\lambda, \kappa }  
                aux2_: L_{R S, \kappa} &= \sum_{\lambda \in B_c} ( R S \vert \lambda ) ( K^{-T}_{\lambda, \kappa }  

                And we get K^{-T} from a Cholesky decompostion of the two-component J matrix defined below:
                J = \begin{pmatrix}
                    (\alpha \vert \beta)  &  (\alpha \vert \Gamma)\\
                    (\Gamma \vert \alpha) &  (\Gamma \vert \Theta)
                    \end{pmatrix}

          *************************************************************/
            InCoreAsymmRITPI<double> asymmAux1Comp2(aux1_, basisSet2.nBasis, ASYMM_CD_ALG::INT1_AUX, build4I_);
            asymmAux1Comp2.malloc();
            InCoreAsymmRITPI<double> asymmAux2Comp1(basisSet.nBasis, aux2_, ASYMM_CD_ALG::INT2_AUX, build4I_);
            asymmAux2Comp1.malloc();

            // \kappa,\lambda \in B_c    is the union of     \alpha,beta \in B_e      and       \Gamma,Theta in \B_n
            // To build raw 3-index (p q \vert lambda) and (R S \vert lambda), we already have (p q \vert \beta) and (R S \vert \Gamma) from symm CD 
            // We just need to build (p q \vert \Gamma) in asymmAux2Comp1
            //                   and (R S \vert \beta)  in asymmAux1Comp2 (same procedure as single one-component RI)
            if (build4I_) {
              asymmAux1Comp2.prebuilt4Index(basisSet, basisSet2, mol, emPert, EP_ATTRACTION, options);
              asymmAux1Comp2.computeOneCholeskyRawSubTPIPrebuilt4Index();
              asymmAux2Comp1.prebuilt4Index(basisSet, basisSet2, mol, emPert, EP_ATTRACTION, options);
              asymmAux2Comp1.computeOneCholeskyRawSubTPIPrebuilt4Index();
            } else {
              asymmAux1Comp2.computeOneCholeskyRawSubTPILibint(basisSet, basisSet2);
              asymmAux2Comp1.computeOneCholeskyRawSubTPILibint(basisSet2, basisSet);
            }

            // Allocate J matrix (dimension is NBRI1+NBRI2)
            size_t combineNBRI = aux1_->nRIBasis() + aux2_->nRIBasis();
            cqmatrix::Matrix<double> twocenterERI(combineNBRI);

            auto asymmCopyBegin = tick();
            // Copy (\alpha \vert \beta) to upper left corner of square matrix J
            SetMat('N', aux1_->nRIBasis(), aux1_->nRIBasis(),
                   1.0, aux1_->rawERI2C()->pointer(), aux1_->nRIBasis(),
                   twocenterERI.pointer(), combineNBRI);
            // Copy (\Gamma \vert \Theta) to lower right corner of square matrix J
            SetMat('N', aux2_->nRIBasis(), aux2_->nRIBasis(),
                   1.0, aux2_->rawERI2C()->pointer(), aux2_->nRIBasis(),
                   twocenterERI.pointer() + aux1_->nRIBasis() * (combineNBRI + 1), combineNBRI);
            double durAsymmCopy = tock(asymmCopyBegin);
            std::cout<< "  Cholesky-Asymm-TwoIndex-Copy duration = " << durAsymmCopy << " s " << std::endl;

            // Fill (\alpha \vert \Gamma) in upper right corner of square matrix J
            // Achieved by only selecting \alpha in B_e in all (p q \vert \Gamma) elements
            InCoreCholeskyRIERI<double>::extractTwoCenterSubsetFrom3indexERI(
                cd_aux1->pivots(), aux2_->nRIBasis(), NB,
                asymmAux2Comp1.pointer(), aux2_->nRIBasis(),
                twocenterERI.pointer() + aux1_->nRIBasis() * combineNBRI, combineNBRI, false);

            // Lower triangle is not necessary
//            InCoreCholeskyRIERI<double>::extractTwoCenterSubsetFrom3indexERI(
//                cd_aux2->pivots(), aux1_->nRIBasis(), sNB,
//                asymmAux1Comp2.pointer(), aux1_->nRIBasis(),
//                twocenterERI.pointer() + aux1_->nRIBasis(), combineNBRI, false);

            // Calculate K^{-T} by calling cholesky decomposition on J
            InCoreRITPI<double>::halfInverse2CenterERI(twocenterERI);

            size_t maxNB = std::max(NB, sNB);
            auto ijK = CQMemManager::get().malloc<double>(combineNBRI * maxNB * (maxNB+1)/2);
            auto ijK1 = CQMemManager::get().malloc<double>(combineNBRI * maxNB * (maxNB+1)/2);
            auto ijK2 = CQMemManager::get().malloc<double>(combineNBRI * maxNB * (maxNB+1)/2);
            auto Scra = CQMemManager::get().malloc<double>(combineNBRI);

            size_t NB2 = NB * (NB+1)/2;
            
            auto asymmGemm1Begin = tick();
            
            // Compute L_{p q, \kappa} &= \sum_{\lambda \in B_c} ( p q \vert \lambda ) ( K^{-T}_{\lambda, \kappa } in two steps

            // 1. K^{-T} from 0 to NBRI1 column  * (p q \vert alpha) , save in ijK
            blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                       combineNBRI,NB2,aux1_->nRIBasis(),double(1.),twocenterERI.pointer(),combineNBRI,
                       aux1_->rawERI3J(),aux1_->nRIBasis(),double(0.),ijK,combineNBRI);
            // 2. K^{-T} from NBRI1 column to the end * (p q \vert Gamma)  , add to ijK
            blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                       combineNBRI,NB2,aux2_->nRIBasis(),double(1.),twocenterERI.pointer()+aux1_->nRIBasis(),combineNBRI,
                       asymmAux2Comp1.pointer(),aux2_->nRIBasis(),double(1.),ijK,combineNBRI);

            double durAsymmGemm1 = tock(asymmGemm1Begin);
            std::cout<< "  Cholesky-Asymm-Gemm1 duration = " << durAsymmGemm1 << " s " << std::endl;
            
            // Expand to full index 
            auto asymmConvertIndex1Begin = tick();
            std::shared_ptr<InCoreRITPI<double>> combineAux1 =
                std::make_shared<InCoreRITPI<double>>(NB, combineNBRI);
            for (size_t pq = 0; pq < NB2; pq++) {
              auto pqAna = anaCompound(pq);
              std::copy_n(&ijK[pq*combineNBRI], combineNBRI,
                          combineAux1->pointer()+combineNBRI*toSquare(pqAna.first, pqAna.second, NB));
              std::copy_n(&ijK[pq*combineNBRI], combineNBRI,
                          combineAux1->pointer()+combineNBRI*toSquare(pqAna.second, pqAna.first, NB));
            }
            double durConvertIndex1 = tock(asymmConvertIndex1Begin);
            std::cout<< "  Cholesky-Asymm-ConvertIndex1 duration = " << durConvertIndex1 << " s " << std::endl;
            
            auto asymmGemm2Begin = tick();
            NB2 = sNB * (sNB+1)/2;

            // Compute L_{R S, \kappa} &= \sum_{\lambda \in B_c} ( R S \vert \lambda ) ( K^{-T}_{\lambda, \kappa } in two steps

            // 1. K^{-T} from 0 to NBRI1 column  * (R S \vert alpha) , save in ijK
            blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                       combineNBRI,NB2,aux1_->nRIBasis(),double(1.),twocenterERI.pointer(),combineNBRI,
                       asymmAux1Comp2.pointer(),aux1_->nRIBasis(),double(0.),ijK,combineNBRI);
            // 2. K^{-T} from NBRI1 column to the end * (R S \vert Gamma)  , add to ijK           
            blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                       combineNBRI,NB2,aux2_->nRIBasis(),
                       double(1.),twocenterERI.pointer() + aux1_->nRIBasis(),combineNBRI,
                       aux2_->rawERI3J(),aux2_->nRIBasis(),double(1.),ijK,combineNBRI);
            double durAsymmGemm2 = tock(asymmGemm2Begin);
            std::cout<< "  Cholesky-Asymm-Gemm2 duration = " << durAsymmGemm2 << " s " << std::endl;

            // Expand to full index 
            auto asymmConvertIndex2Begin = tick();
            std::shared_ptr<InCoreRITPI<double>> combineAux2 =
                std::make_shared<InCoreRITPI<double>>(sNB, combineNBRI);
            for (size_t pq = 0; pq < NB2; pq++) {
              auto pqAna = anaCompound(pq);
              std::copy_n(&ijK[pq*combineNBRI], combineNBRI,
                          combineAux2->pointer()+combineNBRI*toSquare(pqAna.first, pqAna.second, sNB));
              std::copy_n(&ijK[pq*combineNBRI], combineNBRI,
                          combineAux2->pointer()+combineNBRI*toSquare(pqAna.second, pqAna.first, sNB));
            }
            double durConvertIndex2 = tock(asymmConvertIndex2Begin);
            std::cout<< "  Cholesky-Asymm-ConvertIndex2 duration = " << durConvertIndex2 << " s " << std::endl;

            CQMemManager::get().free(ijK, ijK1, ijK2, Scra);

            aux1_->clearRawERI();
            aux2_->clearRawERI();

            // Swap out aux1_ and aux2_ to be the combineNBRI 3-index tensors
            aux1_ = combineAux1;
            aux2_ = combineAux2;

          } else if (asymmCDalg_ == ASYMM_CD_ALG::COMBINEMATRIX) {
            CErr("COMBINEMATRIX algorithm for asymm-RI NYI");
          }
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }
      } else {

        /*************************************************************
          ElecAux Method (Single One-Component RI in the paper)

            // partialTPI_ to be built:
                L_{R S,\alpha} = \sum_{\beta \in B_e} ( R S \vert \beta ) K^{-T} _{\beta \alpha} 

        *************************************************************/
        if (asymmCDalg_ != ASYMM_CD_ALG::INT1_AUX)
          CErr("Only auxiliary basis for integral 1 available, "
               "but requested algorithm is not INT1_AUX");
        if (cd_aux1) {
          if(cd_aux1->pivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {-1., 1.};
            cd_aux1->computeAOInts(basisSet, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if (!partialTPI_) malloc();
          
          std::cout<< "     * Using elec aux basis" << std::endl;
          // Compute raw three-index TPI ( R S \vert \beta) 
          auto topCDEP = tick();
          if (eri4I_){
            std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
            computeOneCholeskyRawSubTPIPrebuilt4Index();
          } else {
            std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
            computeOneCholeskyRawSubTPILibint(basisSet, basisSet2);
          }
          // Compute  K^{-T} ( R S \vert \beta )
          computeOneCholeskyPartialTPI();
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }

      }

    } else {
      
      /*************************************************************
          ProtAux Method (Single One-Component RI in the paper)

          // partialTPI_ to be built:
              L_{p q, \Theta} = \sum_{\Gamma \in B_n} ( p q \vert \Gamma ) K^{-T} _{\Gamma \Theta} 

      *************************************************************/
      if (asymmCDalg_ != ASYMM_CD_ALG::INT2_AUX)
        CErr("Only auxiliary basis for integral 2 available, "
             "but requested algorithm is not INT2_AUX");
      if (aux2_) {
        std::shared_ptr<InCoreCholeskyRIERI<double>> cd_aux2 = 
            std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2_);
        if (cd_aux2) {
          if(cd_aux2->pivots().empty()) {
            HamiltonianOptions temp_opt;
            temp_opt.particle = {1., ProtMassPerE};
            cd_aux2->computeAOInts(basisSet2, mol, emPert, ELECTRON_REPULSION, temp_opt);
          }

          if (!partialTPI_) malloc();

          std::cout<< "     * Using prot aux basis" << std::endl; 
          auto topCDEP = tick();
          // Compute raw three-index TPI ( p q \vert \Gamma)
          if(eri4I_){
            std::cout<< "     * Computing PartialTPI for (ee|pp) with prebuilt 4-index (ee|pp)\n" << std::endl;
            computeOneCholeskyRawSubTPIPrebuilt4Index();
          }else{
            std::cout<< "     * Computing PartialTPI for (ee|pp) on the fly\n" << std::endl;
            computeOneCholeskyRawSubTPILibint(basisSet2, basisSet);
          }
          // Compute K^{-T} ( p q \vert \Gamma)
          computeOneCholeskyPartialTPI();
          auto durCDEP = tock(topCDEP);
          std::cout<< "  Cholesky-Asymm-Total duration = " << durCDEP << " s " << std::endl;

        } else {
          CErr("Aux-basis InCoreAsymmRITPI NYI");
        }

      } else {
        CErr("No aux available in IncoreAsymmRITPI::computeAOInts");
      }

    }
    
    if(reportError_) reportError(basisSet, basisSet2, mol, emPert);

    std::cout << std::endl << BannerEnd << std::endl;

  } // InCoreAsymmRITPI<double>::computeAOInts

  template <>
  void InCoreAsymmRITPI<dcomplex>::computeAOInts(
      BasisSet&, BasisSet&, Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("GIAO integral evaluation is NOT implemented in class InCoreAsymmRITPI.");
  } // InCoreAsymmRITPI<dcomplex>::computeAOInts

}; // namespace ChronusQ
