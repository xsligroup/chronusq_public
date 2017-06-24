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

#include <mcscf.hpp>
#include <cibuilder/casci/impl.hpp>
#include <cibuilder/rasci/impl.hpp>
#include <cqlinalg/eig.hpp>
#include <itersolver.hpp>

// #define _DEBUG_CISOLVER_IMPL

namespace ChronusQ {

  // Different type
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  CISolver<MatsT,IntsT>::CISolver(const CISolver<MatsU,IntsT> & other) {};
  
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  CISolver<MatsT,IntsT>::CISolver(CISolver<MatsU,IntsT> && other) {};
  
  /*
   * \brief solve CI Hamiltonian for N States 
   * 
   * Options:
   * 1. build full matrix and then direct diagonalize it
   * 2. pass functions of build sigma to davidson iterative solver
   *
   */ 

  template <typename MatsT, typename IntsT>
  void CISolver<MatsT, IntsT>::solveCI(MCWaveFunction<MatsT, IntsT> & mcwfn,EMPerturbation& pert) {
    
    size_t NDet = mcwfn.NDet;
    size_t nR   = mcwfn.NStates;
    auto & StateEnergy = mcwfn.StateEnergy;
    auto & CIVecs = mcwfn.CIVecs;

    if (alg_ == CI_FULL_MATRIX) {
      
      std::cout << "  Diagonalize CI Full Hamitonian Matrix ... \n" << std::endl;
      dcomplex * Energy = CQMemManager::get().malloc<dcomplex>(NDet);
      MatsT * fullH     = CQMemManager::get().malloc<MatsT>(NDet*NDet); 
      MatsT * EigVec    = CQMemManager::get().malloc<MatsT>(NDet*NDet);
      MatsT * dummy     = nullptr;
      
      ProgramTimer::tick("Full Matrix");
      mcwfn.ciBuilder->buildFullH(mcwfn, fullH);
      ProgramTimer::tock("Full Matrix");
      
#ifdef _DEBUG_CISOLVER_IMPL
      prettyPrintSmart(std::cout,"HH full H", fullH, NDet, NDet, NDet);
      std::cout << " <0|H|0> = "  << std::setprecision(16) 
                << *fullH + mcwfn.reference().molecule().nucRepEnergy + mcwfn.InactEnergy
                << std::endl;
#endif


      if( pert_has_type(pert,Magnetic) ){
        std::cout << "Magnetic field detected. GeneralEigen will be used." << std::endl;
        GeneralEigen('N','V', NDet, fullH, NDet, Energy, dummy, 1, EigVec, NDet);
      } else {
        HermetianEigen('V', 'L', NDet, fullH, NDet, Energy);
        std::copy_n(fullH,NDet*NDet,EigVec);
      }
      
#ifdef _DEBUG_CISOLVER_IMPL
      prettyPrintSmart(std::cout,"HH Eigenvectors", EigVec, NDet, NDet, NDet);
      prettyPrintSmart(std::cout,"HH Eigenvalues", Energy, NDet, 1, NDet);
#endif
    
      // copy over eigenvalues and eigenvectors
      for (auto i = 0ul; i < nR; i++) {
        StateEnergy[i] = std::real(Energy[i]);
        std::copy_n(EigVec+i*NDet, NDet, CIVecs[i]);
      }
      CQMemManager::get().free(Energy, fullH, EigVec);

    } else if (alg_ == CI_DAVIDSON) {
      
      // build diagonal H
      MatsT  * diagH  = CQMemManager::get().malloc<MatsT>(NDet); 
      mcwfn.ciBuilder->buildDiagH(mcwfn, diagH);
      
#ifdef _DEBUG_CISOLVER_IMPL
      prettyPrintSmart(std::cout,"HH Diagonal H ", diagH, NDet, 1, NDet);
#endif
      
      // set davidson parameters
      size_t kG = nDavidsonGuess_;
      size_t m  = std::max(maxDavidsonSpace_, kG);
      size_t nG = kG*nR;

      dcomplex * curEig = CQMemManager::get().malloc<dcomplex>(nG);
      
      using LinearTrans_t = typename IterDiagonalizer<MatsT>::LinearTrans_t;
      
      // define linear transformation
      LinearTrans_t func = [&] (size_t nVec, SolverVectors<MatsT> &V, SolverVectors<MatsT> &AV) {
          
#ifdef CQ_ENABLE_MPI
	    // disable MPI for now
	    ROOT_ONLY(mcwfn.comm);
#endif
        
        ProgramTimer::tick("Sigma");
	    mcwfn.ciBuilder->buildSigma(mcwfn, nVec, tryGetRawVectorsPointer(V), tryGetRawVectorsPointer(AV));

        ProgramTimer::tock("Sigma");
      
      }; 
      
      // define preconditioner
      LinearTrans_t PC = [&] (size_t nVec, SolverVectors<MatsT> &S, SolverVectors<MatsT> &R) {

#ifdef CQ_ENABLE_MPI
	    // disable MPI for now
	    ROOT_ONLY(mcwfn.comm);
#endif
	    this->davidsonPC(nVec, NDet, tryGetRawVectorsPointer(S), tryGetRawVectorsPointer(R), diagH, curEig);
      }; 
      
      std::cout << "  Use Davidson Diagonalization ... \n" << std::endl;

      Davidson<MatsT> davidson(mcwfn.comm, NDet, 5, maxIter_,
        vectorConv_, nR, func, PC);
       
      davidson.setM(m);
      davidson.setkG(kG);
      if( pert_has_type(pert,Magnetic) ){
        std::cout << "Magnetic field detected. GeneralEigen will be used." << std::endl;
      }else
        davidson.setHerm(true);
      davidson.setEigForT(curEig);
      davidson.setGuess(nG, [&] (size_t nGuess, SolverVectors<MatsT> &Guess, size_t N) {
#ifdef CQ_ENABLE_MPI
	    // disable MPI for now
	    ROOT_ONLY(mcwfn.comm);
#endif
          this->davidsonGS(nGuess, N, diagH, tryGetRawVectorsPointer(Guess));
      });

      if (!energyRefs_.empty()) {
        davidson.setEnergySpecific(energyRefs_);
      }

      davidson.run();
	  
      if (MPIRank(mcwfn.comm) == 0) {
        // copy over eigenvalues and eigenvectors
        auto davidsonEig = davidson.eigVal();
        auto davidsonVec = tryGetRawVectorsPointer(*davidson.VR());
	    for (auto i = 0ul; i < nR; i++) {
	      StateEnergy[i] = std::real(davidsonEig[i]);
          std::copy_n(davidsonVec+i*NDet, NDet, CIVecs[i]);
	    }
      }

      CQMemManager::get().free(curEig, diagH);

    } else{
      CErr("Haven't Inplement Other Diagonalization yet");
    }
      
    // add other parts of the energy 
     double EOther = mcwfn.reference().molecule().nucRepEnergy + mcwfn.InactEnergy + mcwfn.EFieldNuc;
     for (auto i = 0ul; i < nR; i++) StateEnergy[i] += EOther;
  
  } // CISolver::solveCI
  
  
  // davidson guess
  template <typename MatsT, typename IntsT>
  void CISolver<MatsT,IntsT>::davidsonGS(size_t nGuess, size_t N, MatsT * HDiag, MatsT * Guess) { 
    
    std::cout << "  * use unit vector guess based on diagonal elements:" << std::endl;
     
    std::vector<size_t> guessIndices;        
    guessIndices.reserve(nGuess);

    std::vector<size_t> indx(N, 0);
    std::iota(indx.begin(), indx.end(), 0);    

    std::stable_sort(indx.begin(), indx.end(), 
      [&] (size_t i , size_t j) {
        return std::real(HDiag[i]) < std::real(HDiag[j]);
      }
    );

    // for energy specific Davidson
    if (energyRefs_.empty()) {

      std::copy_n(indx.begin(), nGuess, std::back_inserter(guessIndices));

    } else {

      size_t lowGuess = nGuess;
      for (auto & pair: energyRefs_) {
        lowGuess -= pair.second * nDavidsonGuess_;
      } 
      std::vector<size_t>::iterator curIterBegin = indx.begin() + lowGuess;
      std::copy(indx.begin(), curIterBegin, std::back_inserter(guessIndices)); 
      double Eoffset = std::real(HDiag[indx[0]]);

      for (auto & pair: energyRefs_) {
        double curERef =  pair.first + Eoffset;
        size_t curNGuess = nDavidsonGuess_ * pair.second;
        curIterBegin = std::lower_bound(curIterBegin, indx.end(), curERef,
                       [&HDiag](size_t i, double x){ return std::real(HDiag[i]) < x; });
        if (curIterBegin <= indx.end() - curNGuess) {
          std::copy_n(curIterBegin, curNGuess, std::back_inserter(guessIndices));
          curIterBegin += curNGuess;
        } else {
          CErr("No enough element above the reference energy to select.");
        }
      }
    }

    std::fill_n(Guess, nGuess*N, MatsT(0.));
    for(auto i = 0ul; i < nGuess; i++) { 
      Guess[i*N+guessIndices[i]] = 1.0;
      std::cout << "    " << std::setw(9) << std::left << guessIndices[i]
                << std::setw(40) << std::left << HDiag[guessIndices[i]] << std::endl;
    }

  } // CISolver::davidsonGS
  
  // davidson preconditioner
  template <typename MatsT, typename IntsT>
  void CISolver<MatsT,IntsT>::davidsonPC(size_t nVec, size_t N, MatsT * R, MatsT * S, 
    MatsT * HDiag, dcomplex * curEig) {
    
    auto D = HDiag;
    auto L = curEig;
    double small = 1e-12;
    double P; 
    
    for (auto i = 0ul; i < nVec; i++) {
      auto Ri = R + i*N; 
      auto Si = S + i*N; 
      for (auto j = 0ul; j < N; j++) { 
        P = std::real(D[j]) - std::real(L[i]);
	    if(std::abs(P) > small) Si[j] = Ri[j] / P;
      }
    }

  } // CISolver::davidsonPC

}; // namespace ChronusQ
