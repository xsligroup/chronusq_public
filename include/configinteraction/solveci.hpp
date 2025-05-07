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
    SolverVectors<MatsT>& Guess) { 
  
  std::cout << "  * use unit vector guess based on diagonal elements:" << std::endl;
   
  std::vector<size_t> nGuessMinDiagHIndices;
  std::vector<MatsT> nGuessMinDiagHValues;

  diagH.getKIndicesAndValues(nGuess, 0ul, nGuessMinDiagHIndices, nGuessMinDiagHValues, 
      [](const MatsT& a, const MatsT& b) { return std::real(a) < std::real(b); }
  );
  
  Guess.clear();

  for(auto i = 0ul; i < nGuess; i++) { 
    Guess.set(nGuessMinDiagHIndices[i], i, MatsT(1.));
    std::cout << "    " << std::setw(9) << std::left << nGuessMinDiagHIndices[i]
              << std::setw(40) << std::left << nGuessMinDiagHValues[i] << std::endl;
  }
  

} // davidsonGuess


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
    CIVectors->fromRawVectors(0ul, ciEigenvectors, 0ul, nRoots); 
    if (MPIRank(this->comm) == 0) {
      for (auto i = 0ul; i < nRoots; i++) {
        this->StateEnergy[i] = std::real(ciEigenvalues[i]);
      }
    } 
    CQMemManager::get().free(ciEigenvalues);

  } else if (ciSettings.ciAlg == CI_DAVIDSON) {
   
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
    davidson.setGuess(nG, [&] (size_t nGuess, SolverVectors<MatsT> &Guess, size_t N) {
      davidsonGuess(std::min(nGuess, N), *diagH, Guess);
    });

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
  } else{
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
