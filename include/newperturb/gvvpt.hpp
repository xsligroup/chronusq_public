/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

#include <posthartreefock.hpp>
#include <newcibuilder/dascibuilder.hpp>
#include <newperturb.hpp>
#include <mointstransformer/impl.hpp>
#include <util/matout.hpp>

// #define _DEBUG_GVVPT
#define _TANH_REGULARIZATION

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::buildHX(std::pair<MatsT,MatsT>& HX, 
    const std::shared_ptr<DistributedVectors<MatsT>>& sigmaIQ, 
    size_t I) {

    HX = {MatsT(0.), MatsT(0.)};
    const MatsT E0_I = MatsT(E0_[I]);

    const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
    const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

    auto localViewSigmaIQ = bracategoricalSpace.createLocalCIVectorsView(*sigmaIQ, 0ul, 1ul);
    const auto CBegin = localViewSigmaIQ.localCategoryBegin();
    const auto CEnd   = localViewSigmaIQ.localCategoryEnd(); 
    const MatsT* ptrFock_I = &(*fockDiag_)(0,I);
    
    for (auto i = CBegin; i < CEnd; ++i) {
      const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
      const auto nEs = cat.SpaceOccupations();
      MatsT *catDiagH = localViewSigmaIQ.getCategoryPointer(i);
      const size_t nDets = cat.nDeterminants();
    	const size_t nDetsPerThread = std::ceil(double(nDets) / GetNumThreads());
      #pragma omp parallel default(shared)
      {
        std::vector<size_t> occ(nCorrE);
        MatsT localHX1 = MatsT(0.);
        MatsT localHX2 = MatsT(0.);
        auto detCatGen = cat.template generator<uint64_t>();
     	 	size_t iBegin = nDetsPerThread * GetThreadID();
     	 	size_t iEnd   = std::min(nDets, iBegin + nDetsPerThread);

        detCatGen.visitDeterminants(iBegin, iEnd,
     	    [&] (size_t addr, const auto& dets)  { 
       	  determinantsToOccs(dets, nEs, nOrbs, occ);
          MatsT sigma_iq = catDiagH[addr];
          MatsT H0qq = MatsT(0.);
          for (size_t orbind : occ) 
            H0qq += ptrFock_I[orbind];
//          H0qq = 0.5 * (H0qq + E0_I) 
//            + (0.5 * std::sqrt((H0qq - E0_I) + (4.0 * sigma_iq)));

         
          MatsT Xqi = (sigma_iq / (E0_I - H0qq));
          localHX1 += SmartConj(sigma_iq) * Xqi;
          localHX2 += SmartConj(Xqi) * sigma_iq;
        });
        #pragma omp critical
        {
          HX.first += localHX1;
          HX.second += localHX2;
        }
      }
    }


  } // DasPerturb::buildHIJ

  
  /*
  *
  * HX.first -> H_(IQ) X_(QJ)
  * HX.second -> X*(QI) * H_(QJ)
  *
  */
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::buildHX(std::pair<MatsT,MatsT>& HX, 
    const std::shared_ptr<DistributedVectors<MatsT>>& sigmaIQ, 
    const std::shared_ptr<DistributedVectors<MatsT>>& sigmaJQ, 
    size_t I, size_t J) {

    HX = {MatsT(0.), MatsT(0.)};
    const MatsT E0_I = MatsT(E0_[I]);
    const MatsT E0_J = MatsT(E0_[J]);
    const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
    const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

    auto localViewSigmaIQ = bracategoricalSpace.createLocalCIVectorsView(*sigmaIQ, 0ul, 1ul);
    auto localViewSigmaJQ = bracategoricalSpace.createLocalCIVectorsView(*sigmaJQ, 0ul, 1ul);
    const auto CBegin = localViewSigmaIQ.localCategoryBegin();
    const auto CEnd   = localViewSigmaIQ.localCategoryEnd(); 
    const MatsT* ptrFock_I = &(*fockDiag_)(0,I);
    const MatsT* ptrFock_J = &(*fockDiag_)(0,J);

    for (auto i = CBegin; i < CEnd; ++i) {
      const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
      const auto nEs = cat.SpaceOccupations();
      MatsT *catDiagHIQ = localViewSigmaIQ.getCategoryPointer(i);
      MatsT *catDiagHJQ = localViewSigmaJQ.getCategoryPointer(i);
      const size_t nDets = cat.nDeterminants();
    	const size_t nDetsPerThread = std::ceil(double(nDets) / GetNumThreads());

      #pragma omp parallel default(shared)
      {
        std::vector<size_t> occ(nCorrE);
        auto detCatGen = cat.template generator<uint64_t>();
     	 	size_t iBegin = nDetsPerThread * GetThreadID();
     	 	size_t iEnd   = std::min(nDets, iBegin + nDetsPerThread);
        MatsT localHX1 = MatsT(0.);
        MatsT localHX2 = MatsT(0.);

        detCatGen.visitDeterminants(iBegin, iEnd,
          [&] (size_t addr, const auto& dets) {
          determinantsToOccs(dets, nEs, nOrbs, occ);
          MatsT sigma_iq = catDiagHIQ[addr];
          MatsT sigma_jq = catDiagHJQ[addr];
          MatsT H0qq_I = MatsT(0.);
          MatsT H0qq_J = MatsT(0.);
          for (size_t orbind : occ) {
            H0qq_I += ptrFock_I[orbind];
            H0qq_J += ptrFock_J[orbind];
          }
//          H0qq_I = 0.5 * (H0qq_I + E0_I) 
//            + (0.5 * std::sqrt((H0qq_I - E0_I) + (4.0 * sigma_iq)));
//          H0qq_J = 0.5 * (H0qq_J + E0_J) 
//            + (0.5 * std::sqrt((H0qq_J - E0_J) + (4.0 * sigma_jq)));

          MatsT Xqi = (sigma_iq / (E0_I - H0qq_I));
          MatsT Xqj = (sigma_jq / (E0_J - H0qq_J));

          localHX1 += SmartConj(sigma_iq) * Xqj;
          localHX2 += SmartConj(Xqi) * sigma_jq;
        });
        #pragma omp critical
        {
          HX.first += localHX1;
          HX.second += localHX2;
        }
      }
    }

  } // DasPerturb::buildXPQ


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT, IntsT>::formEffectiveH(std::shared_ptr<cqmatrix::Matrix<MatsT>>& effH) {

    const size_t PSize = Target_States_.size();
    const size_t SSize = PTopts.SECONDARYROOTS;
    const size_t PSdim = PSize + SSize;

    // Obtain H_PP for the CI States:
    auto H_PP = RefMCWfn_->StateEnergy;
    auto H_QI = PTFactory_->braCategoricalSpace()->
      constructDistributedCIVectors<MatsT>(this->comm, 1ul);
    std::pair<MatsT,MatsT> HX;
    auto* CIVector = dynamic_cast<DistributedVectors<MatsT>*>(RefMCWfn_->CIVectors.get());

    for (size_t I = 0; I < PSdim; ++I) {
      ptBuilder_->buildSigma(1ul, *CIVector, I, *H_QI, 0ul);
      for (size_t J = 0; J <= I; ++J) {
        
        // Case I: H_pp (I = J)
        if (I == J and I < PSize) {
          buildHX(HX, H_QI, I);
          (*effH)(I,J) += 0.5 * ( HX.first + HX.second );
        }

        // Case II: H_pp'
        else if (I != J and I < PSize) {
          auto H_QJ = PTFactory_->braCategoricalSpace()->
            constructDistributedCIVectors<MatsT>(this->comm, 1ul);
          ptBuilder_->buildSigma(1ul, *CIVector, J, *H_QJ, 0ul);
          buildHX(HX, H_QI, H_QJ, I, J);
          (*effH)(I,J) += 0.5 * ( HX.first + HX.second );
          (*effH)(J,I) += (*effH)(I,J);
          H_QJ.reset();
        }

        // Case III: H_SP
        else if (I >= PSize and J < PSize) {
          auto H_QJ = PTFactory_->braCategoricalSpace()->
            constructDistributedCIVectors<MatsT>(this->comm, 1ul);
          ptBuilder_->buildSigma(1ul, *CIVector, J, *H_QJ, 0ul);
          buildHX(HX, H_QI, H_QJ, I, J);
          (*effH)(I,J) += HX.first;
          (*effH)(J,I) += (*effH)(I,J);
          H_QJ.reset();           
        } 
        //else if (I < PSize and J >= PSize) {
        //  auto H_QJ = PTFactory_->braCategoricalSpace()->
        //    constructDistributedCIVectors<MatsT>(this->comm, 1ul);
        //  ptBuilder_->buildSigma(1ul, *CIVector, J, *H_QJ, 0ul);
        //  buildHX(HX, H_QJ, H_QI, J, I);
        //  (*effH)(I,J) += SmartConj(HX.first);
        //  H_QJ.reset();           
        //}
      }
      H_QI->clear();
    }
    
    // Reduce the EffH on all nodes:
    MPIAllReduce(effH->pointer(), PSdim * PSdim, effH->pointer(), this->comm);
    for (size_t I = 0; I < PSdim; ++I) (*effH)(I,I) += H_PP[I];

  } // DasPerturb::formEffectiveH


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT, IntsT>::diagEffH() {

    const bool stateavg = PTopts.STATEAVERAGE;
    const size_t PSize = Target_States_.size();
    const size_t SSize = PTopts.SECONDARYROOTS;
    const size_t PSdim = PSize + SSize;

    if (PSdim > RefMCWfn_->NStates) 
    { 
      std::cout << "Ref NStates: " << RefMCWfn_->NStates << std::endl;
      std::cout << "PSDim : " << PSdim << std::endl;
      CErr("PSpace and SSpace are ill-defined w.r.t CI Space"); 
    }

    // Form the Zeroth-order Hamiltonian and the Zeroth-Order Energies:
    if  (!stateavg) { 
      formFockDiag(); 
      computeZeroEnergy();
    }
    else {
      //semiCanonicalize();
      formFockDiagStateAvg();
      computeSAZeroEnergy();
    }

#ifdef _DEBUG_GVVPT
    std::cout << " Zero-Order Energy : " << std::endl;
    for (double e0 : E0_) std::cout << e0 << "   ";
#endif

    auto EffH = std::make_shared<cqmatrix::Matrix<MatsT>>(PSdim, PSdim);
    EffH->clear();
    formEffectiveH(EffH);

#ifdef _DEBUG_GVVPT
    prettyPrintSmart(std::cout, "Eff H-Matrix: ", EffH->pointer(), PSdim, PSdim, PSdim, 1, 8);
#endif

    // A single NaN/Inf amplitude poisons the whole effective Hamiltonian and
    // makes syev bail without touching the eigenvalues, which would otherwise
    // be reported as a run with no correlation energy at all.
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
