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
#include <singleslater/multiparticless.hpp>
#include <fockbuilder/rofock.hpp>
#include <util/matout.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::printProperties() {
    applyToEachLabeled([](const std::string& label, SubSSPtr& ss) {
      std::cout << "\n";
      std::cout << "Subsystem: " << label
                << " (charge = " << ss->particle.charge << ")\n";
      ss->printMOInfo(std::cout);
      ss->printSpin(std::cout);
      ss->printMiscProperties(std::cout);
    });

    this->printMultipoles(std::cout);
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> MultiParticleSS<MatsT,IntsT>::getFock() {
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> focks;
    applyToEach([&focks](SubSSPtr& ss) {
      for( auto& X: ss->getFock() )
        focks.push_back(X);
    },
    this->scfControls.NEOSubSystemOpt);
    return focks;
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> MultiParticleSS<MatsT,IntsT>::getOnePDM() {
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> dens;
    applyToEach([&dens](SubSSPtr& ss) {
      for( auto& X: ss->getOnePDM() )
        dens.push_back(X);
    },
    this->scfControls.NEOSubSystemOpt);
    return dens;
  };

  template <typename MatsT, typename IntsT>
  std::vector<cqmatrix::Matrix<MatsT>> MultiParticleSS<MatsT,IntsT>::getOnePDMOrtho() {
    std::vector<cqmatrix::Matrix<MatsT>> dens;
    applyToEach([&dens](SubSSPtr& ss) {
      for( auto& X: ss->getOnePDMOrtho() )
        dens.push_back(X);
    },
    this->scfControls.NEOSubSystemOpt);
    return dens;
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::setOnePDMOrtho(cqmatrix::Matrix<MatsT> *tempOnePDMOrtho) {
    size_t i = 0;
    applyToEach([&](SubSSPtr& ss) {
      ss->setOnePDMOrtho(&tempOnePDMOrtho[i]);
      if(ss->nC == 1) {
        if (ss->iCS) i++;
        else i += 2;
      }
      else i++;
    });
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::setOnePDMAO(cqmatrix::Matrix<MatsT> *tempOnePDMAO) {
    size_t i = 0;
    applyToEach([&](SubSSPtr& ss) {
      ss->setOnePDMAO(&tempOnePDMAO[i]);
      if(ss->nC == 1) {
        if (ss->iCS) i++;
        else i += 2;
      }
      else i++;
    });
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::ortho2aoDen() {
    applyToEach([](SubSSPtr& ss) {
      ss->ortho2aoDen();
    },
    this->scfControls.NEOSubSystemOpt);
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::ortho2aoMOs() {
    applyToEach([](SubSSPtr& ss) {
      ss->ortho2aoMOs();
    },
    this->scfControls.NEOSubSystemOpt);
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::ao2orthoDen() {
    applyToEach([](SubSSPtr& ss) {
      ss->ao2orthoDen();
    },
    this->scfControls.NEOSubSystemOpt);
  };

  template <typename MatsT, typename IntsT>
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> MultiParticleSS<MatsT, IntsT>::getOrtho() {
    std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho;
    applyToEach([&ortho](SubSSPtr& ss) {
      for( auto& X: ss->getOrtho() )
        ortho.push_back(X);
    },
    this->scfControls.NEOSubSystemOpt);
    return ortho;
  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT, IntsT>::setDenEqCoeff(bool val) {
    applyToEach([&val](SubSSPtr& ss) {
      ss->setDenEqCoeff(val);
    },
    this->scfControls.NEOSubSystemOpt);
  };

  template<typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT, IntsT>::initializeSCF() {
    this->moCoefficients.clear();
    this->moEigenvalues.clear();

    auto baseBuilderOf = [](SubSSPtr& ss) -> FockBuilder<MatsT,IntsT>* {
      auto* fb = ss->fockBuilder.get();
      if (auto* ipfb = dynamic_cast<InterParticleFockBase<MatsT,IntsT>*>(fb)) {
        fb = ipfb->getIntraParticleUpstream();
      }
      return fb;
    };

    // Setup MO reference vector
    applyToEach([this, &baseBuilderOf](SubSSPtr& ss) {
      bool iRO = (dynamic_cast<ROFock<MatsT, IntsT>*>(baseBuilderOf(ss)) != nullptr);
      if( iRO ) {
        this->moCoefficients.emplace_back(ss->mo[0]);
      } else {
        for( auto& m: ss->mo ) {
          this->moCoefficients.emplace_back(m);
        }
      }
    },
    this->scfControls.NEOSubSystemOpt);

    // Setup Eigenvalue vector
    applyToEach([this, &baseBuilderOf](SubSSPtr& ss) {
      this->moEigenvalues.push_back(ss->eps1);
      bool iRO = (dynamic_cast<ROFock<MatsT, IntsT>*>(baseBuilderOf(ss)) != nullptr);
      if( ss->nC == 1 and not(ss->iCS or iRO) )
        this->moEigenvalues.push_back(ss->eps2);
    },
    this->scfControls.NEOSubSystemOpt);

    // Symmetric per-subsystem ERI3J redistribution
    applyToEachLabeled([this](const std::string& label, SubSSPtr& ss) {
      if (auto tpi = std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(ss->aoints_->TPI)) {
        if (auto eri3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(tpi->eri3j())) {
          if (tpi->redistribute()) {
            std::cout << "Redistributing " << label
                      << " ERI3J to be split over auxiliary basis functions" << std::endl;
            eri3j->redistributeToSplitNBRI();
          }
        }
      }
    });

    // Asymmetric cross-pair ERI3J redistribution
    for (size_t i = 0; i < order_.size(); ++i) {
      for (size_t j = i + 1; j < order_.size(); ++j) {
        const auto& labelA = order_[i];
        const auto& labelB = order_[j];
        
        // Skip if no interaction between this pair
        if (!interIntegrals.at(labelA).count(labelB)) {
          continue;
        }

        auto ep_ints = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(
            interIntegrals.at(labelA).at(labelB).second);
        if (!ep_ints) continue;              // not an asymmetric RI object (e.g. 4-index) — skip
        if (!ep_ints->redistribute()) continue;

        std::cout << "Redistributing asymmetric ERI3J for [" << labelA << "-" << labelB
                  << "] to be split over auxiliary basis functions" << std::endl;

        auto asymmAlg = ep_ints->asymmCDalg();
        if (asymmAlg == ASYMM_CD_ALG::INT1_AUX or asymmAlg == ASYMM_CD_ALG::INT2_AUX) {
          if (auto asymm3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ep_ints->partialTPI())) {
            asymm3j->redistributeToSplitNBRI();
          } else {
            CErr("Failed to cast partial TPI to a DistributedERI3J object for [" + labelA + "-" + labelB + "]");
          }
        } else if (asymmAlg == ASYMM_CD_ALG::COMBINEAUXBASIS) {
          auto aux1 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ep_ints->getAux1()->eri3j());
          auto aux2 = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(ep_ints->getAux2()->eri3j());
          if (aux1 and aux2) {
            aux1->redistributeToSplitNBRI();
            aux2->redistributeToSplitNBRI();
          } else {
            CErr("Failed to cast aux1 or aux2 to a DistributedERI3J object for [" + labelA + "-" + labelB + "]");
          }
        } else {
          CErr("Unsupported asymmetric CD algorithm for ERI3J redistribution for [" + labelA + "-" + labelB + "]");
        }
      }
    }
  };

  template<typename MatsT, typename IntsT>
  bool MultiParticleSS<MatsT, IntsT>::secondSCF() {
    if(!this->scfControls.NEOStepwiseOpt) return false;

    // Stepwise is only defined for exactly 2 subsystems (one electron + one proton)
    if(order_.size() != 2)
      CErr("NEO stepwise optimization requires exactly two subsystems (electron + single proton subsystem)");

    // Check for stepwise convergence
    bool energyConverged = std::abs(this->totalEnergy - this->lastE) < this->scfControls.eneConvTol;
#ifdef CQ_ENABLE_MPI
    // Broadcast whether or not we're converged to ensure that all
    // MPI processes exit the NEO-SCF simultaneously
    if(MPISize(this->comm) > 1) MPIBCast(energyConverged, 0, this->comm);
#endif

    // We'll also do one final energyOnly step for populating all properties
    if(energyConverged && !this->scfControls.NEOSubSystemOpt.size())
      return false;

    if(energyConverged) {
      std::cout << "          Converged both subsystems, running one energy only iteration" << std::endl;
      this->scfControls.NEOSubSystemOpt.clear();
      this->scfControls.energyOnly = true;
      initializeSCF();
      return true;
    }

    this->lastE = this->totalEnergy;
    std::cout << "          Converged " << this->scfControls.NEOSubSystemOpt[0] << std::endl;

    // Swap active subsystem. Determine the two labels generically: electron is "E",
    // the other is whatever the single proton subsystem is labeled.
    std::string electronLabel = "E";
    std::string protonLabel;
    for(auto& l : order_) if(l != electronLabel) protonLabel = l;

    std::map<std::string,std::string> swaps = {
      {electronLabel, protonLabel}, {protonLabel, electronLabel}
    };
    this->scfControls.NEOSubSystemOpt[0] = swaps.at(this->scfControls.NEOSubSystemOpt[0]);
    initializeSCF();
    return true;
  };

}; // namespace ChronusQ

