/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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
#ifdef CQ_ENABLE_SPARSE
#include <posthartreefock.hpp>
#include <newcibuilder/dascibuilder.hpp>
#include <newperturb.hpp>
#include <mointstransformer/impl.hpp>
#include <util/matout.hpp>

#define _TANH_REGULARIZATION

namespace ChronusQ {


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::buildHX(std::pair<MatsT,MatsT>& HX,
    DistributedSparseVectors<MatsT>& sigmaQ,
    size_t I) {

    #ifdef _DEBUG_GVVPT
      std::cout << " In Pure buildHX: " << I << std::endl;
    #endif

    HX = {MatsT(0.), MatsT(0.)};
    const MatsT E0_I = MatsT(E0_[I]);
    const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
    const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

    auto localViewSigmaIQ = bracategoricalSpace.createLocalCIVectorsView(sigmaQ, 0ul, 1ul);
    auto& nonZeroIQVec = sigmaQ.getVecs().getCol(I);
    const auto CBegin = localViewSigmaIQ.localCategoryBegin();
    const auto CEnd   = localViewSigmaIQ.localCategoryEnd();
    const MatsT* ptrFock_I = &(*fockDiag_)(0,I);

    #pragma omp parallel default(shared)
    {
      std::vector<size_t> occ(nCorrE);
      std::vector<uint64_t> dets;
      MatsT localHX1 = MatsT(0.);
      MatsT localHX2 = MatsT(0.);

      #pragma omp for schedule(dynamic)
      for (auto i = CBegin; i < CEnd; ++i) {
        const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
        const auto nEs = cat.SpaceOccupations();
        const size_t catOffset = localViewSigmaIQ.getCatOffset(i);
        const size_t nextCatOffset = (i + 1 < CEnd)
          ? localViewSigmaIQ.getCatOffset(i + 1) : localViewSigmaIQ.localLength();

        auto startIt = std::lower_bound(nonZeroIQVec.begin(), nonZeroIQVec.end(),
          std::make_pair(catOffset, MatsT(0.)),
          [](const auto& a, const auto& b){ return a.first < b.first; });
        if (startIt == nonZeroIQVec.end() || startIt->first >= nextCatOffset)
          continue;
        // Direct addressing: visitDeterminants(addr, addr+1, ...) allocates a
        // determinant-string buffer on every call, i.e. once per non-zero.
        const auto addresser = cat.template addresser<uint64_t>();
        dets.assign(nEs.size(), 0ul);

        for (auto it = startIt; it != nonZeroIQVec.end() && it->first < nextCatOffset; ++it) {
          addresser.addressToBitStrings(it->first - catOffset, dets);
          determinantsToOccs(dets, nEs, nOrbs, occ);

          MatsT H0qq = MatsT(0.);
          for (size_t orbind : occ)
            H0qq += ptrFock_I[orbind];
          // Get the denominator:
          auto sigma_iq = it->second;
#ifdef _TANH_REGULARIZATION
          MatsT Xqi = (std::abs(E0_I - H0qq) > 1e-12)
            ? std::tanh(H0qq - E0_I) * (sigma_iq / (E0_I - H0qq)) : -sigma_iq;
#else
          MatsT Xqi = (std::abs(E0_I - H0qq) > 1e-12)
            ? (sigma_iq / (E0_I - H0qq)) : MatsT(0.);
#endif
          localHX1 += SmartConj(sigma_iq) * Xqi;
          localHX2 += SmartConj(Xqi) * sigma_iq;
        }
      }
      #pragma omp critical
      {
        HX.first  += localHX1;
        HX.second += localHX2;
      }
    }

  } // DasPerturb::buildHIJ

   
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::buildHX(std::pair<MatsT,MatsT>& HX,
    DistributedSparseVectors<MatsT>& sigmaQ,
    size_t I, size_t J) {

#ifdef _DEBUG_GVVPT
      std::cout << "Mixed buildHX: " << I << " | " << J << std::endl;
#endif

    HX = {MatsT(0.), MatsT(0.)};
    const MatsT E0_I = MatsT(E0_[I]);
    const MatsT E0_J = MatsT(E0_[J]);
    const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
    const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

    auto localViewSigmaIQ = bracategoricalSpace.createLocalCIVectorsView(sigmaQ, 0ul, 1ul);
    auto& nonZeroIQVec = sigmaQ.getVecs().getCol(I);
    auto& nonZeroJQVec = sigmaQ.getVecs().getCol(J);
    const auto CBegin = localViewSigmaIQ.localCategoryBegin();
    const auto CEnd   = localViewSigmaIQ.localCategoryEnd();

    const MatsT* ptrFock_I = &(*fockDiag_)(0,I);
    const MatsT* ptrFock_J = &(*fockDiag_)(0,J);

    #pragma omp parallel default(shared)
    {
      std::vector<size_t> occ(nCorrE);
      std::vector<uint64_t> dets;
      MatsT localHX1 = MatsT(0.);
      MatsT localHX2 = MatsT(0.);

      #pragma omp for schedule(dynamic)
      for (auto i = CBegin; i < CEnd; ++i) {
        const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
        const auto nEs = cat.SpaceOccupations();
        const size_t catOffset = localViewSigmaIQ.getCatOffset(i);
        const size_t nextCatOffset = (i + 1 < CEnd)
          ? localViewSigmaIQ.getCatOffset(i + 1) : localViewSigmaIQ.localLength();

        auto startIt = std::lower_bound(nonZeroIQVec.begin(), nonZeroIQVec.end(),
          std::make_pair(catOffset, MatsT(0.)),
          [](const auto& a, const auto& b){ return a.first < b.first; });
        if (startIt == nonZeroIQVec.end() || startIt->first >= nextCatOffset)
          continue;
        auto jt = std::lower_bound(nonZeroJQVec.begin(), nonZeroJQVec.end(),
          std::make_pair(startIt->first, MatsT(0.)),
          [](const auto& a, const auto& b){ return a.first < b.first; });

        // Direct addressing: visitDeterminants(addr, addr+1, ...) allocates a
        // determinant-string buffer on every call, i.e. once per non-zero.
        const auto addresser = cat.template addresser<uint64_t>();
        dets.assign(nEs.size(), 0ul);

        for (auto it = startIt; it != nonZeroIQVec.end() && it->first < nextCatOffset; ++it) {
          while (jt != nonZeroJQVec.end() && jt->first < it->first) {
            ++jt;
          }
          if (jt != nonZeroJQVec.end() && jt->first == it->first) {
            addresser.addressToBitStrings(it->first - catOffset, dets);
            determinantsToOccs(dets, nEs, nOrbs, occ);

            MatsT H0qq_I = MatsT(0.);
            MatsT H0qq_J = MatsT(0.);
            for (size_t orbind : occ) {
              H0qq_I += ptrFock_I[orbind];
              H0qq_J += ptrFock_J[orbind];
            }
            // Get the denominator:
            auto sigma_iq = it->second;
            auto sigma_jq = jt->second;
#ifdef _TANH_REGULARIZATION
            MatsT Xqi = (std::abs(E0_I - H0qq_I) > 1e-12)
              ? std::tanh(H0qq_I - E0_I) * (sigma_iq / (E0_I - H0qq_I)) : -sigma_iq;
            MatsT Xqj = (std::abs(E0_J - H0qq_J) > 1e-12)
              ? std::tanh(H0qq_J - E0_J) * (sigma_jq / (E0_J - H0qq_J)) : -sigma_jq;
#else
            MatsT Xqi = (std::abs(E0_I - H0qq_I) > 1e-12)
              ? (sigma_iq / (E0_I - H0qq_I)) : MatsT(0.);
            MatsT Xqj = (std::abs(E0_J - H0qq_J) > 1e-12)
              ? (sigma_jq / (E0_J - H0qq_J)) : MatsT(0.);
#endif
            localHX1 += SmartConj(sigma_iq) * Xqj;
            localHX2 += SmartConj(Xqi) * sigma_jq;
          }
        }
      }
      #pragma omp critical
      {
        HX.first  += localHX1;
        HX.second += localHX2;
      }
    }

  } // DasPerturb::buildXPQ


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT, IntsT>::formEffectiveHSparse(std::shared_ptr<cqmatrix::Matrix<MatsT>>& effH) {

    const size_t PSize = Target_States_.size();
    const size_t SSize = PTopts.SECONDARYROOTS;
    auto PSdim = PSize + SSize;

    // Obtain H_PP for the CI States:
    auto H_PP = RefMCWfn_->StateEnergy;
    double eps = PTopts.EPS;
    std::pair<MatsT,MatsT> HX;
    auto* CIVector = dynamic_cast<DistributedSparseVectors<MatsT>*>(RefMCWfn_->CIVectors.get());

    // Build all <Q|H|I> in one multi-vector sigma build. Rebuilding H_QJ inside
    // the J loop cost O(PSdim^2) sigma builds where PSdim distinct vectors
    // exist; one nVec = PSdim pass also shares the excitation-list traversal
    // and the MPI broadcast of the CI vectors across all columns.
    auto H_Q = PTFactory_->braCategoricalSpace()->
      constructDistributedSparseCIVectors<MatsT>(this->comm, PSdim);
    std::vector<dcomplex> zeroCurEigvalues(PSdim, dcomplex(0., 0.));
    ptBuilder_->buildSigma(PSdim, *CIVector, 0ul, *H_Q, 0ul, zeroCurEigvalues.data(), eps);

    // Build Effective Hamiltonian. effH is Hermitian and only its lower
    // triangle is referenced by HermitianEigen(...,'L',...) below, so the
    // J > I half (the old Case IV, effH(I,J) = conj(effH(J,I))) is skipped.
    for (size_t I = 0; I < PSdim; ++I) {
      for (size_t J = 0; J <= I; ++J) {

        // Case I: H_pp (I = J)
        if (I == J and I < PSize) {
          buildHX(HX, *H_Q, I);
          (*effH)(I,J) += 0.5 * ( HX.first + HX.second );
        }

        // Case II: H_pp' (J < I < PSize)
        else if (I < PSize) {
          buildHX(HX, *H_Q, I, J);
          (*effH)(I,J) += 0.5 * ( HX.first + HX.second );
        }

        // Case III: H_SP (I >= PSize > J)
        else if (J < PSize) {
          buildHX(HX, *H_Q, I, J);
          (*effH)(I,J) += HX.first;
        }

        // I, J both in the secondary space: not built, as before.
      }
    }
    H_Q.reset();

    // Reduce the effH:
    MPIAllReduce(effH->pointer(), PSdim * PSdim, effH->pointer(), this->comm);
    for (size_t I = 0; I < PSdim; ++I) (*effH)(I,I) += H_PP[I];

  } // DasPerturb::formEffectiveH


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT, IntsT>::diagEffHSparse() {

    const bool stateavg = PTopts.STATEAVERAGE;
    const size_t PSize = Target_States_.size();
    const size_t SSize = PTopts.SECONDARYROOTS;
    auto PSdim = PSize + SSize;

    if (PSdim > RefMCWfn_->NStates) 
      { CErr("PSpace and SSpace are ill-defined w.r.t CI Space"); }

    // Form the Zeroth-order Hamiltonian and the Zeroth-Order Energies:
    if  (!stateavg) { 
      formFockDiag(); 
      computeZeroEnergy();
    } else { 
      formFockDiagStateAvg();
      computeSAZeroEnergy();
    }

#ifdef _DEBUG_GVVPT
    std::cout << " Zero-Order Energy : " << std::endl;
    for (double e0 : E0_) std::cout << e0 << "   ";
#endif

    // Build Effective GVVPT2 Hamiltonian:
    auto EffH = std::make_shared<cqmatrix::Matrix<MatsT>>(PSize+SSize, PSize+SSize);
    EffH->clear();
    formEffectiveHSparse(EffH);

#ifdef _DEBUG_GVVPT
    prettyPrintSmart(std::cout, "Eff H-Matrix: ", EffH->pointer(), PSdim, PSdim, PSdim, 1, 8);
#endif

    // See diagEffH: a non-finite element makes syev leave evals untouched.
    for (size_t IJ = 0; IJ < PSdim * PSdim; ++IJ)
      if (not std::isfinite(std::abs(EffH->pointer()[IJ])))
        CErr("Non-finite element in the GVVPT2 effective Hamiltonian");

    std::vector<double> evals(PSdim);
    int info = HermitianEigen('V', 'L', PSdim, EffH->pointer(), PSdim, evals.data());
    if (info != 0)
      CErr("Diagonalization of the GVVPT2 effective Hamiltonian failed");
    E2_ = std::move(evals);

#ifdef _DEBUG_GVVPT
    prettyPrintSmart(std::cout, "Eff H EigenValues: ", E2_.data(), PSdim, 1, PSdim);
#endif

  } // DasPerturb::diagEffH


} // ChronusQ
#endif
