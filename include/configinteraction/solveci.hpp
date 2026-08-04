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

#include <configinteraction.hpp>
#include <cqlinalg/eig.hpp>
#include <itersolver.hpp>

// #define _DEBUG_CISOLVER_IMPL
namespace ChronusQ {

namespace {

// davidson guess
template <typename MatsT>
void davidsonGuess(size_t nGuess, const DistributedVectors<MatsT>& diagH, 
    SolverVectors<MatsT>& Guess,
   std::vector<std::pair<double, size_t>> energyRefs ) {
  
  std::cout << "  * use unit vector guess based on diagonal elements:" << std::endl;
   
  std::vector<size_t> nGuessMinDiagHIndices;
  std::vector<MatsT> nGuessMinDiagHValues;

  diagH.getKIndicesAndValues(nGuess, 0ul, nGuessMinDiagHIndices, nGuessMinDiagHValues, energyRefs,
      [](const MatsT& a, const MatsT& b) { return std::real(a) < std::real(b); }
  );
  
  Guess.clear();

  for(auto i = 0ul; i < nGuess; i++) { 
    std::cout << "    " << std::setw(9) << std::left << nGuessMinDiagHIndices[i]
              << std::setw(40) << std::left << nGuessMinDiagHValues[i] << std::endl;
    Guess.set(nGuessMinDiagHIndices[i], i, MatsT(1.));
  }
  

} // davidsonGuess

#ifdef CQ_ENABLE_SPARSE
// davidson guess w/o building diagH
template <typename MatsT>
void davidsonGuessNoDiagH(size_t nGuess,
    SolverVectors<MatsT>& Guess, const CategoricalSpace& categoricalSpace, std::shared_ptr<NewCIBuilder<MatsT>> ciBuilder) {

  std::cout << "  * use unit vector guess based on on-the-fly diagonal elements:" << std::endl;
  
  DistributedSparseVectors<MatsT> *GuessPtr = nullptr;
  size_t GuessShift = 0ul;
  tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(Guess,
    [&] (auto& GuessRef, size_t shift) {
      GuessPtr = &GuessRef;
      GuessShift = shift;
    }
  );


  std::vector<size_t> nGuessMinDiagHIndices;
  std::vector<MatsT> nGuessMinDiagHValues;
  auto const comp = +[](const MatsT& a, const MatsT& b) { return std::real(a) < std::real(b); };

  ciBuilder->largestDiagH(*GuessPtr, categoricalSpace, nGuess, nGuessMinDiagHIndices, nGuessMinDiagHValues, comp);

  Guess.clear();

  for(auto i = 0ul; i < nGuess; i++) {
    Guess.set(nGuessMinDiagHIndices[i], i, MatsT(1.));
    std::cout << "    " << std::setw(9) << std::left << nGuessMinDiagHIndices[i]
              << std::setw(40) << std::left << nGuessMinDiagHValues[i] << std::endl;
  }


} // davidsonGuess w/o building diagH
#endif


// davidson preconditioner
template <typename MatsT>
void davidsonPreconditioner(MatsT* localR, MatsT* localS, 
  const MatsT* localDiagonals, size_t localLD, const dcomplex * curEigenvalues, 
  size_t nVec) {
  
  const double small = 1e-12;

  #pragma omp parallel for schedule(static) collapse(2) default(shared)
  for (auto i = 0ul; i < localLD; ++i) {  
    for (auto j = 0ul; j < nVec; ++j) {
      size_t ij = i + j * localLD;
      double P = std::real(localDiagonals[i]) - std::real(curEigenvalues[j]);
      if(std::abs(P) > small) localS[ij] = localR[ij] / P;
    }
  }

} // davidsonPreconditioner


/*
 * Steps:
 * 1. convert R to sparse local view (which has category info
 * 2.  iterate inside the iterator loop, get the category i of iter.row()
 * 3. using the detstobits routine get the integral info for the address it.row() (see buildDiagH for info, here addr = it.row()
 */

// davidson sparse preconditioner w/o Hamiltonian diagonal build
#ifdef CQ_ENABLE_SPARSE
template <typename MatsT>
void davidsonSparsePreconditioner(DistributedSparseVectors<MatsT>& localR, size_t Rshift, DistributedSparseVectors<MatsT>& localS, size_t Sshift, 
  const CategoricalSpace& categoricalSpace, const DeterminantFactory& detFactory, const std::shared_ptr<IntegralsCollection>&  moints,  
  const dcomplex * curEigenvalues,
  size_t nVec, double C_eps) {

  const double small = 1e-12;

  if(&localR.getVecs() != &localS.getVecs()) { 
    CErr("S and R in davidsonPreconditioner must be the same");
  }
  
  auto RLocalView = categoricalSpace.createLocalCIVectorsView(localR, Rshift, nVec);

  const auto& nOrbs = detFactory.nOrbitalsInEachSpace();
  const size_t nCorrE = detFactory.nTotalCorrElectrons();

  const auto& hCore_tt = *(moints->template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
  const auto& antiSymmetricERI_ttuu = *(moints->template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));

  for (auto j = 0ul; j < nVec; j++) {

    #pragma omp parallel for schedule(static) default(shared)
    for(auto it = (RLocalView.getVecsPtr())->begin(Rshift + j); it != (RLocalView.getVecsPtr())->end(Rshift + j); it++) {
      
      size_t iCat = RLocalView.localCategoryBegin();
      size_t adjustedOffset = iCat == RLocalView.localCategoryEnd() - 1? RLocalView.localLength() : RLocalView.getCatOffest(iCat + 1);

      while(it->first >= adjustedOffset) {
	iCat++;
	adjustedOffset = iCat == RLocalView.localCategoryEnd() - 1? RLocalView.localLength() : RLocalView.getCatOffest(iCat + 1);
      }

      const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*categoricalSpace.getCategory(iCat));
      const auto nEs = cat.SpaceOccupations();
      std::vector<size_t> occ(nCorrE);
      auto detCatGen = cat.template generator<uint64_t>();
      auto addresser = detCatGen.addresser();
      std::vector<uint64_t> braDetStrings(detCatGen.nSpaces(), 0ul);

      addresser.addressToBitStrings(it->first - RLocalView.getCatOffest(iCat), braDetStrings);

      determinantsToOccs(braDetStrings, nEs, nOrbs, occ);
      MatsT tmp = MatsT(0.);
      for(const auto& t : occ) {
        tmp += hCore_tt(t, 0);
        for(const auto& u : occ) {
          tmp += MatsT(0.5) * antiSymmetricERI_ttuu(t, u);
        }
      }

      MatsT localDiagonal = tmp;
      double P = std::real(localDiagonal) - std::real(curEigenvalues[j]);
      if(std::abs(P) > small) it->second /= P;
    }
  }

  //Get rid of small coefficients (relative to the reference = epsilon * norm)
  std::vector<double> norms(nVec);
  for(size_t j = 0; j < nVec; j++) {
    norms[j] = localS.norm2F(j + Sshift, 1);
  }

  //prune
  #pragma omp parallel for schedule(static) default(shared)
  for(size_t j = 0; j < nVec; j++) {
    localS.getVecs().prune(C_eps * norms[j], Sshift + j);
  }

} // davidsonSparsePreconditioner
#endif
} // namespace

/*
 * \brief solve CI Hamiltonian for N States 
 * 
 * Options:
 * 1. build full matrix and then direct diagonalize it
 * 2. pass functions of build sigma to davidson iterative solver
 *
 */ 
template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::solveCI() {
  const size_t nDet = nDeterminants(); // total number of determinants
  const size_t nRoots   = this->NStates; // total number of CI roots to solve
  const auto& ketCategoricalSpace = *detFactory->ketCategoricalSpace();

  if (ciSettings.ciAlg == CI_FULL_MATRIX) {
    
    std::cout << "  Diagonalize CI Full Hamiltonian Matrix ... \n" << std::endl;
    dcomplex * ciEigenvalues = CQMemManager::get().malloc<dcomplex>(nDet); // storage to save the ciEigenvalues
    auto distributedFullH = ketCategoricalSpace.constructDistributedCIVectors<MatsT>(this->comm, nDet); 
    std::shared_ptr<RawVectors<MatsT>> fullH;
    
    MatsT *dummy = nullptr, *fullH_ptr = nullptr;
    
    ProgramTimer::tick("Full Matrix");
    ciBuilder->buildFullH(*distributedFullH);

    if (MPISize() > 1) {
      fullH = std::make_shared<RawVectors<MatsT>>(this->comm, nDet, nDet);
      distributedFullH->setRawVectors(0ul, *fullH, 0ul, nDet);
      distributedFullH = nullptr;
      if (MPIRank(this->comm) == 0) fullH_ptr = fullH->getPtr();
    } else {
      fullH_ptr = distributedFullH->getLocalPtr(); 
    }

    RawVectors<MatsT> ciEigenvectors(this->comm, nDet, nDet);
    
    if (MPIRank(this->comm) == 0) {
#ifdef _DEBUG_CISOLVER_IMPL
      std::cout << " <0|H|0> = "  << std::setprecision(16) 
                << fullH_ptr[0] + this->reference().molecule().nucRepEnergy + this->coreEnergy
                << std::endl;
      prettyPrintSmart(std::cout,"HH ciEigenvalues", fullH_ptr, nDet, nDet, nDet);
#endif
    
//    GeneralEigen('N', 'V', nDet, fullH_ptr, nDet, ciEigenvalues, dummy, 1, ciEigenvectors.getPtr(), nDet);
      HermitianEigen('V', 'L', nDet, fullH_ptr, nDet, ciEigenvalues);
      std::copy_n(fullH_ptr,nDet*nDet,ciEigenvectors.getPtr());

#ifdef _DEBUG_CISOLVER_IMPL
      ciEigenvectors.print(std::cout,"HH Eigenvectors");
      prettyPrintSmart(std::cout,"HH ciEigenvalues", ciEigenvalues, nDet, 1, nDet);
#endif
    }
    
    ProgramTimer::tock("Full Matrix");
    
    // copy over ciEigenvalues and eigenvectors
    std::shared_ptr<DistributedVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(CIVectors);

    CIVectorsCast->fromRawVectors(0ul, ciEigenvectors, 0ul, nRoots); 
    if (MPIRank(this->comm) == 0) {
      for (auto i = 0ul; i < nRoots; i++) {
        this->StateEnergy[i] = std::real(ciEigenvalues[i]);
      }
    } 
    CQMemManager::get().free(ciEigenvalues);

  } else if (ciSettings.ciAlg == CI_DAVIDSON and (not ciSettings.SparseDavidson)) {
   
    // release the memory of current CI Vecs to have more memory for intermediates 
    CIVectors = nullptr;
    
    // build diagonal H
    auto diagH = ketCategoricalSpace.constructDistributedCIVectors<MatsT>(this->comm, 1ul); 
    
    ciBuilder->buildDiagH(*diagH, ketCategoricalSpace);
    
#ifdef _DEBUG_CISOLVER_IMPL
    diagH->print(std::cout, "HH Diagonal H ", 0ul, 1ul);
#endif
      
    // set davidson parameters
    size_t kG = ciSettings.nDavidsonGuess;
    size_t m  = std::max(ciSettings.maxDavidsonSpace, kG);
    size_t nG = kG * nRoots;

    dcomplex * curEigenvalues = CQMemManager::get().malloc<dcomplex>(nG);
    
    using LinearTrans_t = typename IterDiagonalizer<MatsT>::LinearTrans_t;
      
    // define linear transformation
    double totalSigmaBuildTime = 0.;
    size_t nTotalSigmaBuilt = 0;

    LinearTrans_t func = [&] (size_t nVec, SolverVectors<MatsT>& V, SolverVectors<MatsT>& AV) {
      
      auto sigmaBuild = tick();
      
      DistributedVectors<MatsT> *VPtr = nullptr, *AVPtr = nullptr;
      size_t VShift = 0ul, AVShift = 0ul;
      tryDowncastReferenceTo<DistributedVectors<MatsT>>(V,
          [&] (auto& VRef, size_t shift) {
            VPtr = &VRef;
            VShift = shift;
          }
      );
      tryDowncastReferenceTo<DistributedVectors<MatsT>>(AV,
          [&] (auto& AVRef, size_t shift) {
            AVPtr = &AVRef;
            AVShift = shift;
          }
      );

      //V.print(std::cout, "V", 0ul, nVec);
	  ciBuilder->buildSigma(nVec, *VPtr, VShift, *AVPtr, AVShift);
      //AV.print(std::cout, "AV", 0ul, nVec);
      
      auto durationSigmaBuild = tock(sigmaBuild);
      std::cout << "      * Sigma Build Time: " << durationSigmaBuild << " s"
                << std::endl << std::endl;
      nTotalSigmaBuilt += nVec;
      totalSigmaBuildTime += durationSigmaBuild;
    }; 
      
    // define preconditioner
    LinearTrans_t PC = [&] (size_t nVec, SolverVectors<MatsT> &S, SolverVectors<MatsT> &R) {
      davidsonPreconditioner(tryGetDistributedVectorsLocalPointer(S),
          tryGetDistributedVectorsLocalPointer(R),
          tryGetDistributedVectorsLocalPointer(*diagH), diagH->localLength(),
          curEigenvalues, nVec);
    };
      
    std::function<std::shared_ptr<SolverVectors<MatsT>>(size_t)> distributedDASCIVecsGen = 
        [&] (size_t nVec) {
           return ketCategoricalSpace.constructDistributedCIVectors<MatsT>(this->comm, nVec); 
        };
    
    std::cout << "  Use Davidson Diagonalization ... \n" << std::endl;

    Davidson<MatsT> davidson(this->comm, nDet, 40, ciSettings.maxCIIter,
                             ciSettings.ciVectorConv, nRoots, func, PC,
                             distributedDASCIVecsGen);
     
    davidson.setM(m);
    davidson.setkG(kG);
    davidson.setEigForT(curEigenvalues);
    davidson.setHerm(true);
    std::vector<std::pair<double, size_t>> energyWindows;
    if (not ciSettings.energyRefs.empty()) {
      for (auto& [energy, energy_roots] : ciSettings.energyRefs)
      {
        energyWindows.emplace_back(energy, energy_roots * kG);
      } 
    } else {
        energyWindows.emplace_back(0.0, nG);
    }
    davidson.setGuess(nG, [&] (size_t nGuess, SolverVectors<MatsT> &Guess, size_t N) {
      davidsonGuess(std::min(nGuess, N), *diagH, Guess, energyWindows);
    });

    if (!ciSettings.energyRefs.empty()) {
      davidson.setEnergySpecific(ciSettings.energyRefs);
    }

    davidson.run();
    
    // copy over eigenvalues and eigenvectors
    auto VR = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(davidson.VR());
    if (VR) {
      CIVectors = VR;
    } else {
      CErr("Davidson is not using DistributedVectors");
    }
    
    std::cout << std::endl << "  * Built " << nTotalSigmaBuilt 
              << " Sigma vectors using " << totalSigmaBuildTime << " s, on average " 
              << totalSigmaBuildTime / nTotalSigmaBuilt << " s per vector" << std::endl;
    
    if (MPIRank(this->comm) == 0) {
      const auto davidsonEig = davidson.eigVal();
      for (auto i = 0ul; i < nRoots; i++) {
  	    this->StateEnergy[i] = std::real(davidsonEig[i]);
      }
    }

    CQMemManager::get().free(curEigenvalues);
  } 

#ifdef CQ_ENABLE_SPARSE 
  else if (ciSettings.ciAlg == CI_DAVIDSON and ciSettings.SparseDavidson) {

    // release the memory of current CI Vecs to have more memory for intermediates
    CIVectors = nullptr;
    std::cout << "Using SPARSE CI Davidson" << std::endl;
    // set davidson parameters
    size_t kG = ciSettings.nDavidsonGuess;
    size_t m  = std::max(ciSettings.maxDavidsonSpace, kG);
    size_t nG = kG * nRoots;

    dcomplex * curEigenvalues = CQMemManager::get().malloc<dcomplex>(nG);

    using LinearTrans_t = typename IterDiagonalizer<MatsT>::LinearTrans_t;
    using SubAFormer_t = typename Davidson<MatsT>::SubAFormer_t;






    SubAFormer_t formSubA = [&] (size_t nVec, size_t nTot, SolverVectors<MatsT>& V, std::vector<MatsT>& VAV) {

      DistributedSparseVectors<MatsT> *VPtr = nullptr;
      size_t VShift = 0ul;
      tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(V,
          [&] (auto& VRef, size_t shift) {
            VPtr = &VRef;
            VShift = shift;
          }
      );

      //Drop null CI coefficients
      //VPtr->getVecs().prune(ciSettings.SparseDavidsonEps);

      VPtr->getVecs().prune(0.0);
      std::cout << "      * Forming subspace matrix" << std::endl;

      ciBuilder->formSubA(nVec, nTot, *VPtr, VShift, VAV);
      //Reduce dot products from all threads into VAV[0]
      //
      for(size_t iVec = 0; iVec < nVec; iVec++)
        for(size_t jVec = 0; jVec < nTot; jVec++)
	  for (auto k = 1ul; k < GetNumThreads(); ++k)
	    VAV[jVec + iVec * nTot] += VAV[jVec + iVec * nTot + k * nTot * nVec];
      
      //Root rank should have all the dot products to form the subspace matrix
      if (MPISize(this->comm) > 1) {
	std::vector<MatsT> reducedVAV(nTot * nVec, MatsT(0.));
        MPIAllReduce(VAV.data(), nTot * nVec, reducedVAV.data(), this->comm);

      //copy reduced results back into the buffer
      for(size_t iVec = 0; iVec < nVec; iVec++)
        for(size_t jVec = 0; jVec < nTot; jVec++)
           VAV[jVec + iVec * nTot] = reducedVAV[jVec + iVec * nTot];

      }
    };




    // define linear transformation
    double totalSigmaBuildTime = 0.;
    size_t nTotalSigmaBuilt = 0;

    LinearTrans_t func = [&] (size_t nVec, SolverVectors<MatsT>& V, SolverVectors<MatsT>& AV) {

      auto sigmaBuild = tick();

      DistributedSparseVectors<MatsT> *VPtr = nullptr, *AVPtr = nullptr;
      size_t VShift = 0ul, AVShift = 0ul;
      tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(V,
          [&] (auto& VRef, size_t shift) {
            VPtr = &VRef;
            VShift = shift;
          }
      );
      tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(AV,
          [&] (auto& AVRef, size_t shift) {
            AVPtr = &AVRef;
            AVShift = shift;
          }
      );

      //Drop null CI coefficients
      VPtr->getVecs().prune(0.0);

      //Get the density of the CI vector
      std::cout << std::endl;
      for (size_t iVec = 0; iVec < nVec; iVec++) {
        size_t myNonZeros = VPtr->getVecs().nonZeros(VShift + iVec);
        size_t nonZeros = MPIAllReduce(myNonZeros, this->comm);
        std::cout << "        Root " << iVec << " has CI vector density: " << nonZeros / (double)(VPtr->length()) * 100 << " %" << std::endl;
      }

      ciBuilder->buildSigma(nVec, *VPtr, VShift, *AVPtr, AVShift, curEigenvalues, ciSettings.SparseDavidsonEps);

      auto durationSigmaBuild = tock(sigmaBuild);
      std::cout << "      * Sigma Build Time: " << durationSigmaBuild << " s"
                << std::endl << std::endl;
      nTotalSigmaBuilt += nVec;
      totalSigmaBuildTime += durationSigmaBuild;
    };






    // define preconditioner
    LinearTrans_t PC = [&] (size_t nVec, SolverVectors<MatsT> &S, SolverVectors<MatsT> &R) {


      DistributedSparseVectors<MatsT> *SPtr = nullptr, *RPtr = nullptr;
      size_t SShift = 0ul, RShift = 0ul;

      tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(S,
      [&] (auto& SRef, size_t shift) {
        SPtr = &SRef;
        SShift = shift;
      }
      );
      tryDowncastReferenceTo<DistributedSparseVectors<MatsT>>(R,
      [&] (auto& RRef, size_t shift) {
        RPtr = &RRef;
        RShift = shift;
      }
      );
      
      davidsonSparsePreconditioner(*SPtr, tryGetDistributedSparseVectorsLocalPointer(S), 
          *RPtr, tryGetDistributedSparseVectorsLocalPointer(R), 
          ketCategoricalSpace, *detFactory, this->moints,  
          curEigenvalues, nVec, ciSettings.SparseDavidsonEps);
    };

    std::function<std::shared_ptr<SolverVectors<MatsT>>(size_t)> distributedDASCIVecsGen =
        [&] (size_t nVec) {
           return ketCategoricalSpace.constructDistributedSparseCIVectors<MatsT>(this->comm, nVec);
        };

    std::cout << "  Use Davidson Diagonalization ... \n" << std::endl;

    Davidson<MatsT> davidson(this->comm, nDet, 40, ciSettings.maxCIIter,
                             ciSettings.ciVectorConv, nRoots, formSubA, func, PC,
                             distributedDASCIVecsGen, ciSettings.checkEigenValue, ciSettings.checkEigenVector, ciSettings.checkResidue);

    davidson.setM(m);
    davidson.setkG(kG);
    davidson.setEigForT(curEigenvalues);
    davidson.setHerm(true);
    davidson.setGuess(nG, [&] (size_t nGuess, SolverVectors<MatsT> &Guess, size_t N) {
      davidsonGuessNoDiagH(std::min(nGuess, N), Guess, ketCategoricalSpace, ciBuilder);
    });

    davidson.run();

    // copy over eigenvalues and eigenvectors
    auto VR = std::dynamic_pointer_cast<DistributedSparseVectors<MatsT>>(davidson.VR());
    if (VR) {
      CIVectors = VR;
    } else {
      CErr("Davidson is not using DistributedSparseVectors");
    }

    std::cout << std::endl << "  * Built " << nTotalSigmaBuilt
              << " Sigma vectors using " << totalSigmaBuildTime << " s, on average "
              << totalSigmaBuildTime / nTotalSigmaBuilt << " s per vector" << std::endl;

    if (MPIRank(this->comm) == 0) {
      const auto davidsonEig = davidson.eigVal();
      for (auto i = 0ul; i < nRoots; i++) {
  	    this->StateEnergy[i] = std::real(davidsonEig[i]);
      }
    }

    CQMemManager::get().free(curEigenvalues);
  } 
#endif

  else{
    CErr("Haven't Implement Other Diagonalization yet");
  }
    
  // add other parts of the energy 
  if (MPIRank(this->comm) == 0) {
    double EOther = this->reference()->molecule().nucRepEnergy + this->coreEnergy;
    for (auto i = 0ul; i < nRoots; i++) this->StateEnergy[i] += EOther;
  }
  MPIBCast(&(this->StateEnergy[0]), nRoots, 0, this->comm);
} // CISolver::solveCI

} // namespace ChronusQ
