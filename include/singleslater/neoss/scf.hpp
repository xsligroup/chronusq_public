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
#pragma once

#include <singleslater.hpp>
#include <singleslater/neoss.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::printProperties() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    applyToEach([](NEOSS<MatsT,IntsT>::SubSSPtr& ss) {
      ss->printMOInfo(std::cout);
      ss->printSpin(std::cout);
      ss->printMiscProperties(std::cout);
    });
    this->printMultipoles(std::cout);
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> NEOSS<MatsT,IntsT>::getFock() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> focks;
    applyToEach([&focks](SubSSPtr& ss) {
      for( auto& X: ss->getFock() )
        focks.push_back(X);
    });
    return focks;
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> NEOSS<MatsT,IntsT>::getOnePDM() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> dens;
    applyToEach([&dens](SubSSPtr& ss) {
      for( auto& X: ss->getOnePDM() )
        dens.push_back(X);
    });
    return dens;
  };

  template <typename MatsT, typename IntsT>
  std::vector<cqmatrix::Matrix<MatsT>> NEOSS<MatsT,IntsT>::getOnePDMOrtho() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    std::vector<cqmatrix::Matrix<MatsT>> dens;
    applyToEach([&dens](SubSSPtr& ss) {
      for( auto& X: ss->getOnePDMOrtho() )
        dens.push_back(X);
    });
    return dens;
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::setOnePDMOrtho(cqmatrix::Matrix<MatsT> *tempOnePDMOrtho) {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;

      auto neoMap = getSubsystemMap();
      auto neoSubsystemOrder = getOrder();
      // Loop over all subsystems
      size_t i = 0;
      for( auto& neoSubsystemLabel: neoSubsystemOrder ) {
        auto &neoSubsystem = neoMap[neoSubsystemLabel];
        neoSubsystem->setOnePDMOrtho(&tempOnePDMOrtho[i]);
        if(neoSubsystem->nC == 1) {
          if (neoSubsystem->iCS) i++;
          else i+=2;
        }
        else i++;
      }
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::setOnePDMAO(cqmatrix::Matrix<MatsT> *tempOnePDMAO) {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;

    auto neoMap = getSubsystemMap();
    auto neoSubsystemOrder = getOrder();
    // Loop over all subsystems
    size_t i = 0;
    for( auto& neoSubsystemLabel: neoSubsystemOrder ) {
      auto &neoSubsystem = neoMap[neoSubsystemLabel];
      neoSubsystem->setOnePDMAO(&tempOnePDMAO[i]);
      if(neoSubsystem->nC == 1) {
        if (neoSubsystem->iCS) i++;
        else i+=2;
      }
      else i++;
    }
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::ortho2aoDen() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    applyToEach([this](SubSSPtr& ss) {
      ss->ortho2aoDen();
    });
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::ortho2aoMOs() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    applyToEach([this](SubSSPtr& ss) {
      ss->ortho2aoMOs();
    });
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::ao2orthoDen() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    applyToEach([this](SubSSPtr& ss) {
      ss->ao2orthoDen();
    });
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> NEOSS<MatsT, IntsT>::getOrtho() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho;
    applyToEach([&ortho](SubSSPtr& ss) {
      for( auto& X: ss->getOrtho() )
        ortho.push_back(X);
    });
    return ortho;
  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT, IntsT>::setDenEqCoeff(bool val) {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;
    applyToEach([&val](NEOSS<MatsT,IntsT>::SubSSPtr& ss) {
      ss->setDenEqCoeff(val);
    });
  };

  template<typename MatsT, typename IntsT>
  void NEOSS<MatsT, IntsT>::initializeSCF() {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;

    this->moCoefficients.clear();
    this->moEigenvalues.clear();

    // Setup MO reference vector
    applyToEach([this](SubSSPtr& ss) {
      for( auto& m: ss->mo )
        this->moCoefficients.emplace_back(m);
    });

    // Setup Eigenvalue vector
    applyToEach([this](SubSSPtr& ss) {
      this->moEigenvalues.push_back(ss->eps1);
      if( ss->nC == 1 and !ss->iCS )
        this->moEigenvalues.push_back(ss->eps2);
    });

  }

  template<typename MatsT, typename IntsT>
  void NEOSS<MatsT, IntsT>::runSCF(EMPerturbation& pert) {
    using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;

    // Initialize properties
    //applyToEach([&pert](SubSSPtr& ss){
    //  ss->getNewOrbitals();
    //  ss->computeProperties(pert);
    //});

    // Setup MO reference vector
    std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> moRefs;
    applyToEach([&moRefs](SubSSPtr& ss) {
      for( auto& m: ss->mo )
        moRefs.emplace_back(m);
    });

    // Setup Eigenvalue vector
    std::vector<double*> epsVec;
    applyToEach([&epsVec](SubSSPtr& ss) {
      epsVec.push_back(ss->eps1);
      if( ss->nC == 1 and !ss->iCS )
        epsVec.push_back(ss->eps2);
    });

    this->orbitalModifier->runOrbitalModifier(pert, moRefs, epsVec);

    applyToEach([](SubSSPtr& ss) {
      ss->ao2orthoFock();
      ss->MOFOCK();
    });

    saveCurrentState();

#ifdef CQ_ENABLE_MPI
  // Broadcast the updated MOs to all MPI processes
    applyToEach([&](SubSSPtr& ss) {
      if( MPISize(ss->comm) > 1 ) {
        std::cerr  << "  *** Scattering the AO-MOs ***\n";
        size_t Nmo = ss->mo[0].dimension();
        MPIBCast(ss->mo[0].pointer(),Nmo*Nmo,0,ss->comm);
        if( ss->nC == 1 and not ss->iCS )
          MPIBCast(ss->mo[1].pointer(),Nmo*Nmo,0,ss->comm);

        std::cerr  << "  *** Scattering EPS ***\n";
        MPIBCast(ss->eps1,Nmo,0,ss->comm);
        if( ss->nC == 1 and not ss->iCS )
          MPIBCast(ss->eps2,Nmo,0,ss->comm);

        std::cerr  << "  *** Scattering FOCK ***\n";
        size_t fockDim = ss->fockMatrix->dimension();
        for(MatsT *mat : ss->fockMatrix->SZYXPointers())
          MPIBCast(mat,fockDim*fockDim,0,ss->comm);

        std::cerr  << "  *** Scattering the 1PDM ***\n";
        size_t denDim = ss->onePDM->dimension();
        for(auto p : ss->onePDM->SZYXPointers())
          MPIBCast(p,denDim*denDim,0,ss->comm);
      }
    });

#endif

  };

}; // namespace ChronusQ

