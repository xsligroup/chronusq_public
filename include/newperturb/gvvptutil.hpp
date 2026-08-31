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
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/eig.hpp>

#define DEBUG_SAVG

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::semiCanonicalize() {

    std::cout << " *** Begin semi-canonicalizing orbitals *** " << std::endl;
    computeSAOneRDM();
    const auto& refSpace = RefMCWfn_->corrSpace;
    size_t nInact = refSpace.nInact;
    size_t nCorrO = refSpace.nCorrO;
    size_t nVirtO = refSpace.nSVirt;
    size_t nMO = nInact + nCorrO + nVirtO;
    size_t ndim = (RefMCWfn_->StateEnergy).size();
    MatsT* dummy = nullptr;
    fockDiag_ = std::make_shared<cqmatrix::Matrix<MatsT>>(nMO,ndim);
    
    auto & hCore = *(this->moints->template getIntegral<OnePInts, MatsT>
                                                    ("hCore_Correlated_Space")); 
    auto & moERI = *(this->moints->template getIntegral<InCore4indexTPI, MatsT>
                                                    ("ERI_Correlated_Space"));

    if (nInact > 0) {
      // Begin with Core Fock Block:
      auto fockCore = std::make_shared<cqmatrix::Matrix<MatsT>>(nInact);
      auto fockCoreVec = std::make_shared<cqmatrix::Matrix<MatsT>>(nInact);

      for (size_t i = 0; i < nInact; ++i)
      for (size_t j = 0; j < nInact; ++j) {
        (*fockCore)(i,j) = hCore(i,j);
        for (size_t p = 0; p < nCorrO + nInact; ++p)
        for (size_t q = 0; q < nCorrO + nInact; ++q) {
          if (p >= nInact and q >= nInact)
            (*fockCore)(i,j) += (*oneRDMSA_)(p-nInact, q-nInact) * (moERI(i,j,p,q)
                - moERI(i,q,p,j));
          else if (p == q) 
            (*fockCore)(i,j) += (moERI(i,j,p,q) - moERI(i,q,p,j));
        }
      }
      // Diagonalize the Fock Core Matrix:
      std::vector<MatsT> diagFockCore(nInact);
//      GeneralEigen('N','V',nInact,fockCore->pointer(),nInact,diagFockCore.data(),dummy,1,fockCoreVec->pointer(),nInact);
      HermitianEigen('V','L',nInact,fockCore->pointer(),nInact,diagFockCore.data());
      for (const auto& val: diagFockCore)
        std::cout << val << " | ";
      std::cout << std::endl;
      
      for (size_t i = 0; i < nInact; ++i) {
        (*fockDiag_)(i,0) = diagFockCore[i]; }

    
    }

    if (nCorrO > 0) {
      // Build Active Fock Block:
      auto fockAct = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
      auto fockActVec = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);

      for (size_t t = nInact; t < nInact+nCorrO; ++t)
      for (size_t u = nInact; u < nInact+nCorrO; ++u) {
        (*fockAct)(t-nInact,u-nInact) = hCore(t,u);
        for (size_t p = 0; p < nInact+nCorrO; ++p)
        for (size_t q = 0; q < nInact+nCorrO; ++q) {
          if (p >= nInact and q >= nInact)
            (*fockAct)(t-nInact,u-nInact) += (*oneRDMSA_)(p-nInact, q-nInact) * (moERI(t,u,p,q)
                  - moERI(t,q,p,u));
          else if (p == q)
            (*fockAct)(t-nInact,u-nInact) += (moERI(t,u,p,q) - moERI(t,q,p,u));
        }
      }
      
      // Diagonalize Active Fock:
      std::vector<MatsT> diagFockAct(nCorrO);
//      GeneralEigen('N','V',nCorrO,fockAct->pointer(),nCorrO,diagFockAct.data(),dummy,1,fockActVec->pointer(),nCorrO);
      HermitianEigen('V','L',nCorrO,fockAct->pointer(),nCorrO,diagFockAct.data());
      
      for (const auto& val: diagFockAct)
        std::cout << val << " | ";
      std::cout << std::endl;

      for (size_t i = 0; i < nCorrO; ++i) {
        (*fockDiag_)(i+nInact,0) = diagFockAct[i]; }
    
    }

    if (nVirtO > 0) {
      // Build Virtual Fock Block:
      auto fockVirt = std::make_shared<cqmatrix::Matrix<MatsT>>(nVirtO);
      auto fockVirtVec = std::make_shared<cqmatrix::Matrix<MatsT>>(nVirtO);
      
      for (size_t a = nInact + nCorrO; a < nInact + nCorrO + nVirtO; ++a)
      for (size_t b = nInact + nCorrO; b < nInact + nCorrO + nVirtO; ++b) {
        (*fockVirt)(a - nInact - nCorrO,b - nInact - nCorrO) = hCore(a,b);
        for (size_t p = 0; p < nInact+nCorrO; ++p)
        for (size_t q = 0; q < nInact+nCorrO; ++q) {
          if (p >= nInact and q >= nInact)
            (*fockVirt)(a - nInact - nCorrO,b - nInact - nCorrO) 
              += (*oneRDMSA_)(p-nInact, q-nInact) * (moERI(a,b,p,q)
              - moERI(a,q,p,b));
          else if (p == q)
            (*fockVirt)(a - nInact - nCorrO,b - nInact - nCorrO) 
              += (moERI(a,b,p,q) - moERI(a,q,p,b));
        }
      }
      // Diagonalize virtual fock block:
      std::vector<MatsT> diagFockVirt(nVirtO);
//      GeneralEigen('N','V',nVirtO,fockVirt->pointer(),nVirtO,diagFockVirt.data(),dummy,1,fockVirtVec->pointer(),nVirtO);
      HermitianEigen('V','L',nVirtO,fockVirt->pointer(),nVirtO,diagFockVirt.data());
      for (const auto& val: diagFockVirt)
        std::cout << val << " | ";
      std::cout << std::endl;
    
      for (size_t i = 0; i < nVirtO; ++i) {
        (*fockDiag_)(i+nInact+nCorrO,0) = diagFockVirt[i]; }

    }
    
    // Repeated copy: To make this just one array!!!
    for (size_t state_index = 1ul; state_index < ndim; ++state_index)
    for (size_t i = 0ul; i < nMO; ++i)
      (*fockDiag_)(i, state_index) = (*fockDiag_)(i, 0);
    std::cout << " *** Formed State-Averaged Fock Diagonals with Semicanonicalized Orbitals *** " <<
      std::endl;

  } // DasPerturb::semiCanonicalize


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::computeSAOneRDM() {
   
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    const size_t ndim = Target_States_.size();
    size_t nCorrO = refMOSpace.nCorrO;
    oneRDMSA_ = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    oneRDMSA_->clear();

    for (size_t i = 0ul; i < ndim; ++i)
      *oneRDMSA_ += (1.0/ndim) * (*RefMCWfn_->oneRDM[Target_States_[i]]);
#ifdef DEBUG_SAVG
    std::cout << "Complete SARDM" << std::endl;
#endif

  }
  
  
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::formFockDiag() {

    std::cout << " *** Build Fock Matrix *** " << std::endl;
    const auto& mrptMOSpace = this->corrSpace;
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    const auto& h1e = *(this->moints->template getIntegral<OnePInts,MatsT>("hCore_Correlated_Space"));
    const auto& ERI = *(this->moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));
    const size_t ndim = (RefMCWfn_->StateEnergy).size();

    size_t nMO = mrptMOSpace.nCorrO;
    size_t nInact = refMOSpace.nInact;
    size_t nCorrO = refMOSpace.nCorrO;
    fockDiag_ = std::make_shared<cqmatrix::Matrix<MatsT>>(nMO, ndim);
#ifdef _DEBUG_GVVPT
    std::cout << nMO << " | " << nInact << " | " << nCorrO << std::endl;
#endif
     
    for (size_t state_index = 0; state_index < ndim; ++state_index)
    for (size_t p = 0ul; p < nMO; ++p) {
      (*fockDiag_)(p, state_index) = h1e(p,p); 
      for (auto k = 0ul; k < nInact; ++k) {
        (*fockDiag_)(p, state_index) += ERI(p,p,k,k) - ERI(p,k,k,p);
      }
      for (auto k = nInact; k < nInact + nCorrO; ++k)
      for (auto l = nInact; l < nInact + nCorrO; ++l) {
        (*fockDiag_)(p, state_index) += (*RefMCWfn_->oneRDM[state_index])(k-nInact,l-nInact) * 
          (ERI(p,p,k,l) - ERI(p,l,k,p));
      }
    }

#ifdef _DEBUG_GVVPT
    std::cout << "DEBUG FOCK DIAG: only ground state" << std::endl;
    for (size_t i = 0; i < nMO; ++i)
      std::cout << (*fockDiag_)(i, 0) << std::endl;
    std::cout << "END DEBUG: FOCK DIAG" << std::endl;
#endif


  } // genFock


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::formFockDiagStateAvg() {

    std::cout << " *** Build Fock Matrix State-Averaged *** " << std::endl;
    computeSAOneRDM();
    const auto& mrptMOSpace = this->corrSpace;
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    const auto& h1e = *(this->moints->template getIntegral<OnePInts,MatsT>("hCore_Correlated_Space"));
    const auto& ERI = *(this->moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));
    const size_t ndim = (RefMCWfn_->StateEnergy).size();

    size_t nMO = mrptMOSpace.nCorrO;
    size_t nInact = refMOSpace.nInact;
    size_t nCorrO = refMOSpace.nCorrO;

    fockDiag_ = std::make_shared<cqmatrix::Matrix<MatsT>>(nMO, ndim);
#ifdef _DEBUG_GVVPT
    std::cout << nMO << " | " << nInact << " | " << nCorrO << std::endl;
#endif

    // The state-averaged Fock is state independent: build one column, copy it.
    for (size_t p = 0ul; p < nMO; ++p) {
      (*fockDiag_)(p, 0) = h1e(p,p);
      for (auto k = 0ul; k < nInact; ++k) {
        (*fockDiag_)(p, 0) += ERI(p,p,k,k) - ERI(p,k,k,p);
      }
      for (auto k = nInact; k < nInact + nCorrO; ++k)
      for (auto l = nInact; l < nInact + nCorrO; ++l) {
        (*fockDiag_)(p, 0) += (*oneRDMSA_)(k-nInact,l-nInact) *
          (ERI(p,p,k,l) - ERI(p,l,k,p));
      }
    }
    for (size_t state_index = 1ul; state_index < ndim; ++state_index)
    for (size_t p = 0ul; p < nMO; ++p)
      (*fockDiag_)(p, state_index) = (*fockDiag_)(p, 0);

#ifdef _DEBUG_GVVPT
    std::cout << "DEBUG FOCK DIAG: only ground state" << std::endl;
    for (size_t i = 0; i < nMO; ++i)
      std::cout << (*fockDiag_)(i, 0) << std::endl;
    std::cout << "END DEBUG: FOCK DIAG" << std::endl;
#endif


  } // SAFock


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::computeZeroEnergy() {
  
    E0_ = RefMCWfn_->StateEnergy;
    std::fill(E0_.begin(), E0_.end(), 0.0);
    const size_t ndim = E0_.size(); 

    std::cout << " *** Compute Zero-Order Energy of Reference *** " << std::endl;
    const auto& mrptMOSpace = this->corrSpace;
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    size_t nMO = mrptMOSpace.nCorrO;
    size_t nInact = refMOSpace.nInact;
    size_t nCorrO = refMOSpace.nCorrO;

    for (size_t state_index = 0; state_index < ndim; ++state_index) 
    for (size_t i = 0ul; i < nInact+nCorrO; i++) {
      if (i < nInact) {
        E0_[state_index] += std::real((*fockDiag_)(i, state_index));
      }
      else {
        E0_[state_index] += std::real((*fockDiag_)(i, state_index) * 
          (*RefMCWfn_->oneRDM[state_index])(i-nInact,i-nInact)); 
      }
    }

#ifdef _DEBUG_GVVPT
    std::cout << "DEBUG ZERO-ENERGY" << std::endl;
    for (size_t i = 0; i < ndim; ++i)
      std::cout << E0_[i] << std::endl;
    std::cout << "END DEBUG: ZERO-ENERGY" << std::endl;
#endif

  } // computeE0


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::computeSAZeroEnergy() {
  
    E0_ = RefMCWfn_->StateEnergy;
    std::fill(E0_.begin(), E0_.end(), 0.0);
    const size_t ndim = E0_.size(); 

    std::cout << " *** Compute Zero-Order Energy of Reference State-Avg *** " << std::endl;
    const auto& mrptMOSpace = this->corrSpace;
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    size_t nMO = mrptMOSpace.nCorrO;
    size_t nInact = refMOSpace.nInact;
    size_t nCorrO = refMOSpace.nCorrO;

    for (size_t state_index = 0; state_index < ndim; ++state_index) 
    for (size_t i = 0ul; i < nInact+nCorrO; i++) {
      if (i < nInact) {
        E0_[state_index] += std::real((*fockDiag_)(i, state_index));
      }
      else {
        E0_[state_index] += std::real((*fockDiag_)(i, state_index) * 
          (*oneRDMSA_)(i-nInact,i-nInact)); 
      }
    }

#ifdef _DEBUG_GVVPT
    std::cout << "DEBUG ZERO-ENERGY" << std::endl;
    for (size_t i = 0; i < ndim; ++i)
      std::cout << E0_[i] << std::endl;
    std::cout << "END DEBUG: ZERO-ENERGY" << std::endl;
#endif

  } // computeE0


} // DasPerturb::ChronusQ

