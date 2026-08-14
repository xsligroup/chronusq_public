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

#include <posthartreefock.hpp>
#include <newcibuilder/dascibuilder.hpp>
#include <newperturb.hpp>
#include <mointstransformer/impl.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  /*
   * Compute Zero-Order Energy for ENPT2
   *
   */
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::CIEnergy() {
  
    this->E0_ = RefMCWfn_->StateEnergy;
    std::cout << " Zero-Order Energies (ENPT2): \n";
    for (size_t i = 0; i < E0_.size(); ++i) {
      std::cout << std::fixed << std::setw(12) <<  std::setprecision(10) 
      << E0_[i] << " | ";
      if ((i + 1) % 4 == 0) std::cout << std::endl;
    }
    for (auto &i: E0_) i -= (this->reference()->molecule().nucRepEnergy + this->coreEnergy);
    std::cout << std::endl << bannerTop << std::endl;


  }  // ENPT2 Zeroth-Order Energy

  /*
  * Cache the state-independent integrals entering the perturber diagonal
  *
  *   <Q|H|Q> = sum_t h_tt + 1/2 sum_tu A_tu,   A_tu = (tt|uu) - (tu|ut)
  *
  */
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::buildDiagIntCache() {

    if (not ttuu2eCache_.empty()) return;
    
    const auto& mrptMOSpace = this->corrSpace;
    const size_t nCorrO = mrptMOSpace.nCorrO;
    const auto& hCore_tt = *(this->moints->template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
    const auto& antiSymmetricERI_ttuu = *(this->moints->template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));
    hDiagCache_.resize(nCorrO);
    ttuu2eCache_.assign(nCorrO * nCorrO, 0.0);

    for (size_t t = 0ul; t < nCorrO; ++t)
      hDiagCache_[t] = std::real(hCore_tt(t, 0));

    // Symmetrize explicitly so that a row and a column gather give identical
    // results; the diagonal stays zero.
    for (size_t u = 0ul; u < nCorrO; ++u)
    for (size_t t = 0ul; t < u; ++t) {
      const double a = 0.5 * (std::real(antiSymmetricERI_ttuu(t,u)) +
                              std::real(antiSymmetricERI_ttuu(u,t)));
      ttuu2eCache_[t + u * nCorrO] = a;
      ttuu2eCache_[u + t * nCorrO] = a;
    }

  } // DasPerturb::buildDiagIntCache

	/*
	* Accumulates the sum over CAS Kets: 
	* 	L.H.S = 1 / (E_KK - E_0) 
	* 	where K: Perturber Space, I: CAS space
	*/
	template <typename MatsT, typename IntsT>
	void DasPerturb<MatsT,IntsT>::computeAmplitudes(size_t& state_index, 
    DistributedVectors<MatsT>& diagPTH, double& ENPT2) {


    const auto& mrptMOSpace = this->corrSpace;
    const size_t nCorrO = mrptMOSpace.nCorrO;
		const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    buildDiagIntCache();
    const double* hDiag = hDiagCache_.data();
    const double* twoE  = ttuu2eCache_.data();
    const double E0 = E0_[state_index];
    const double levelShift = PTopts.LEVELSHIFT;

    // main loop
  	const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
  	const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

  	auto diagHLocalView = bracategoricalSpace.createLocalCIVectorsView(diagPTH, 0ul, 1ul); 

  	for (auto i = diagHLocalView.localCategoryBegin(); i < diagHLocalView.localCategoryEnd(); ++i) {
		  const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
 	   	MatsT *catDiagH = diagHLocalView.getCategoryPointer(i); 
    	const auto nEs = cat.SpaceOccupations();
    	const size_t nDets = cat.nDeterminants();
    	const size_t nDetsPerThread = std::ceil(double(nDets) / GetNumThreads());
    	#pragma omp parallel default(shared)
    	{
     	 	// Each thread gets its own private objects
     	 	std::vector<size_t> occ(nCorrE);
     	 	auto detCatGen = cat.template generator<uint64_t>();
     	 	size_t iBegin = nDetsPerThread * GetThreadID();
     	 	size_t iEnd   = std::min(nDets, iBegin + nDetsPerThread);

     	 	double localE2 = 0.0;
     	 	detCatGen.visitDeterminants(iBegin, iEnd,
     	     	[&] (size_t addr, const auto& dets)  {
     	    determinantsToOccs(dets, nEs, nOrbs, occ);
          // <Q|H|Q>: occ is ascending, so the inner gather runs over a single
          // column of the (real, symmetric) antisymmetrized ERI matrix.
          double diag = 0.0;
          for (size_t a = 0ul; a < nCorrE; ++a) {
            const size_t t = occ[a];
            const double* aCol = twoE + t * nCorrO;
            diag += hDiag[t];
            for (size_t b = 0ul; b < a; ++b) diag += aCol[occ[b]];
          }
          const double D = E0 - diag;
          const double denom = D / (D * D + levelShift);
          MatsT sigma_k = catDiagH[addr];
          localE2 += std::real(SmartConj(sigma_k) * sigma_k) * denom;
          catDiagH[addr] = sigma_k * denom;
     	  });
     	 	#pragma omp atomic
     	 	ENPT2 += localE2;
    	} // parallel region
		}


	}	//DasPerturb::computeLHS


  /*
  *
  * Main compute for EN2 PT: Computes R.H.S -> sigma_k = \Sum_i <K|H|i>c_i
  * Gets denominator from computeAmplitudes
  * Does dot product to obtain [sum_k (1/denominator_k) sigma_k ^ 2]
  *
  */
	template <typename MatsT, typename IntsT>
	void DasPerturb<MatsT,IntsT>::computeEN2(size_t& state_index) {

    auto diagPT = PTFactory_->braCategoricalSpace()->
      constructDistributedCIVectors<MatsT>(this->comm,1ul);
    double ENPT2 = 0.0;
    diagPT->clear(0ul, 1ul);
    auto* CIVector = dynamic_cast<DistributedVectors<MatsT>*>(RefMCWfn_->CIVectors.get());

    // Build PT Sigma
    ProgramTimer::tick("Build PT Sigma");
    ptBuilder_->buildSigma(1ul, *CIVector, state_index, *diagPT, 0ul);
    ProgramTimer::tock("Build PT Sigma");

    // Compute amplitudes (fused with EN2 accumulation)
    ProgramTimer::tick("Compute Amplitudes");
    computeAmplitudes(state_index, *diagPT, ENPT2);
#ifdef CQ_ENABLE_MPI
    ENPT2 = MPIAllReduce(ENPT2, this->comm);
#endif
    ProgramTimer::tock("Compute Amplitudes");

    // Store PT2 corrections:
    std::cout << std::fixed << std::setprecision(10) << 
    "State Index: " << state_index << "  |  " <<
    "State-Specific ENPT2 correlation Energy: " << ENPT2 << 
    "\n" << std::endl;

    E2_.push_back( std::real(ENPT2) + E0_[state_index] + 
        this->coreEnergy + this->reference()->molecule().nucRepEnergy );


  }


  /*
  *
  * SPARSE IMPLEMENTSTION
	* Accumulates the sum over CAS Kets: 
	* 	L.H.S = 1 / (E_KK - E_0) 
	* 	where K: Perturber Space, I: CAS space
	*/
#ifdef CQ_ENABLE_SPARSE
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::computeSparseAmplitudes(
      size_t& state_index,
      DistributedSparseVectors<MatsT>& diagPTH,
      double& ENPT2) {
    
    const auto& mrptMOSpace = this->corrSpace;
    const size_t nCorrO = mrptMOSpace.nCorrO;
    const auto& bracategoricalSpace = *PTFactory_->braCategoricalSpace();
    buildDiagIntCache();
    const double* hDiag = hDiagCache_.data();
    const double* twoE  = ttuu2eCache_.data();
    double kappa = PTopts.LEVELSHIFT;
    const bool regularize = (kappa > 0.0);

    const auto& nOrbs = PTFactory_->nOrbitalsInEachSpace();
    const size_t nCorrE = PTFactory_->nTotalCorrElectrons();

    const double E0 = E0_[state_index];
    auto diagHLocalView = bracategoricalSpace.createLocalCIVectorsView(diagPTH, 0ul, 1ul);
    auto& stateVec = diagPTH.getVecs().getCol(0);
    const auto catBegin = diagHLocalView.localCategoryBegin();
    const auto catEnd   = diagHLocalView.localCategoryEnd();
    
    #pragma omp parallel default(shared) reduction(+:ENPT2)
    {
      std::vector<size_t> occ(nCorrE); // allocated once per thread
      std::vector<uint64_t> dets;      // reused determinant-string buffer
      #pragma omp for schedule(dynamic)
      for (auto i = catBegin; i < catEnd; ++i) {
        const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*bracategoricalSpace.getCategory(i));
        const auto nEs  = cat.SpaceOccupations();
        const size_t catOffset     = diagHLocalView.getCatOffset(i);
        const size_t nextCatOffset = (i + 1 < catEnd)
            ? diagHLocalView.getCatOffset(i + 1)
            : diagHLocalView.localLength();

        auto startIt = std::lower_bound(stateVec.begin(), stateVec.end(),
            std::make_pair(catOffset, MatsT(0.)),
            [](const auto& a, const auto& b){ return a.first < b.first; });
        if (startIt == stateVec.end() || startIt->first >= nextCatOffset)
          continue;
        // Address the determinants directly instead of going through
        // visitDeterminants(addr, addr+1, ...): the generator allocates a
        // fresh string buffer on every call, i.e. once per non-zero.
        const auto addresser = cat.template addresser<uint64_t>();
        dets.assign(nEs.size(), 0ul);

        for (auto it = startIt; it != stateVec.end() && it->first < nextCatOffset; ++it) {

          addresser.addressToBitStrings(it->first - catOffset, dets);
          determinantsToOccs(dets, nEs, nOrbs, occ);

          // <Q|H|Q>: occ is ascending, so the inner gather runs over a single
          // column of the (real, symmetric) antisymmetrized ERI matrix.
          double diag = 0.0;
          for (size_t a = 0ul; a < nCorrE; ++a) {
            const size_t t = occ[a];
            const double* aCol = twoE + t * nCorrO;
            diag += hDiag[t];
            for (size_t b = 0ul; b < a; ++b) diag += aCol[occ[b]];
          }

          double D = E0 - diag;
          double regShift = regularize ? (1.0 - std::exp(-kappa * D * D)) : 1.0;
          MatsT sigma_k = it->second;
          ENPT2 += std::real(SmartConj(sigma_k) * sigma_k) *
                           ((regShift * regShift) / D);
          // Update amplitude in-place
          it->second = (sigma_k * regShift) / D;

        }
      }
    }
  }

  /*
  * SPARSE IMPLEMENTATION:
  * Main compute for EN2 PT: Computes R.H.S -> sigma_k = \Sum_i <K|H|i>c_i
  * Gets denominator from computeAmplitudes
  * Does dot product to obtain [sum_k (1/denominator_k) sigma_k ^ 2]
  *
  */
  template <typename MatsT, typename IntsT>
	void DasPerturb<MatsT,IntsT>::computeEN2Sparse(size_t& state_index) {

    auto diagPT = PTFactory_->braCategoricalSpace()->
      constructDistributedSparseCIVectors<MatsT>(this->comm, 1ul);
    diagPT->clear(0ul, 1ul);
    auto* CIVector = dynamic_cast<DistributedSparseVectors<MatsT>*>(RefMCWfn_->CIVectors.get());
    std::unique_ptr<DistributedSparseVectors<MatsT>> convertedCIVector;
    if (!CIVector) {
      auto* tmp = dynamic_cast<DistributedVectors<MatsT>*>(RefMCWfn_->CIVectors.get());
      if (tmp) {
        convertedCIVector = std::make_unique<DistributedSparseVectors<MatsT>>(*tmp);
        CIVector = convertedCIVector.get();
      } else {
        CErr("CIVectors is neither DistributedSparseVectors nor DistributedVectors");
      }
    }

    // Build PT Sigma
    ProgramTimer::tick("Build PT Sigma");
    std::vector<dcomplex> zeroCurEigvalues(1, (0., 0.));
    double eps = PTopts.EPS;
    ptBuilder_->buildSigma(1ul, *CIVector, state_index, *diagPT, 0ul, zeroCurEigvalues.data(), eps);
    ProgramTimer::tock("Build PT Sigma");

    size_t localNDets = diagPT->getVecs().nonZeros(0ul, 1ul);
    size_t NDets = MPIAllReduce(localNDets, this->comm);  // sum across all MPI ranks
		std::cout << "Number of Non-Zero Perturber determinants:     " << NDets << std::endl;

    // Compute amplitudes
    ProgramTimer::tick("Compute Amplitudes");
    double ENPT2 = 0.;
    computeSparseAmplitudes(state_index, *diagPT, ENPT2);
    ENPT2 = MPIAllReduce(ENPT2, this->comm);
    ProgramTimer::tock("Compute Amplitudes");

    // Store PT2 corrections:
    std::cout << std::fixed << std::setprecision(10) << 
    "State Index: " << state_index << "  |  " <<
    "State-Specific ENPT2 correlation Energy: " << ENPT2 << 
    "\n" << std::endl;
   
    E2_.push_back( std::real(ENPT2) + E0_[state_index] + 
        this->coreEnergy + this->reference()->molecule().nucRepEnergy );
  
  }
#endif

} //DasPerturb :: ChronusQ	
