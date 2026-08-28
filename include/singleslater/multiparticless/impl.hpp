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

#include <singleslater/multiparticless.hpp>
#include <cerr.hpp>
#include <util/timer.hpp>
#include <iostream>
#include <unordered_set>

namespace ChronusQ {

  // Add a subsystem to the NEO object
  //
  // label  Unique name for this subsystem
  // ss     SingleSlater object representing this subsystem
  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::addSubsystem(
    const std::string label,
    std::shared_ptr<SingleSlater<MatsT,IntsT>> ss)
  {

    if(label.empty())           CErr("Quantum subsystem label cannot be empty");
    if(!ss)                     CErr("Cannot add a null quantum subsystem");
    if(subsystems.count(label)) CErr("Duplicate quantum subsystem label: " + label);
    if(!ss->fockBuilder)        CErr("Subsystem " + label + " has no base FockBuilder");

    subsystems.emplace(label, ss);
    order_.push_back(label);

    // Add the base fockBuilder to make sure it has a persistent lifetime
    fockBuilders.insert({label, {ss->fockBuilder}});

    // Create new maps for the new subsystem
    interCoulomb.try_emplace(label);
    interIntegrals.try_emplace(label);
    gradInterInts.try_emplace(label);
    interFockBuilders.try_emplace(label);
    needsFullXCEnergyUpdate = true;

    // Edit the printout for NEO SCF jobs to be electronic reference + NEO.
    // For example: "C-GHF + NEO"
    if (label == "E") {
      this->refShortName_ = ss->refShortName_ + " + NEO";
      this->refLongName_  = ss->refLongName_  + " + NEO";
    }

  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::addInteraction(const std::string label1, const std::string label2,
    const std::shared_ptr<TwoPInts<IntsT>>& ints,
    bool contractSecond) {

    if(label1 == label2)          CErr("Interparticle interaction requires two different subsystems.");
    if(!ints)                     CErr("Cannot add null interparticle integrals for " + label1 + "-" + label2);
    if(!subsystems.count(label1)) CErr("Unknown subsystem label: " + label1);
    if(!subsystems.count(label2)) CErr("Unknown subsystem label: " + label2);
    if(interIntegrals.at(label1).count(label2)) 
      CErr("Interaction already exists between " + label1 + " and " + label2);

    
    auto& ss1 = subsystems.at(label1);
    auto& ss2 = subsystems.at(label2);
    std::shared_ptr<TPIContractions<MatsT,IntsT>> contraction1;
    std::shared_ptr<TPIContractions<MatsT,IntsT>> contraction2;
    std::shared_ptr<DirectTPI<IntsT>> directIntegrals;
    
    // Create contractions for each subsystem
    if(auto direct = std::dynamic_pointer_cast<DirectTPI<IntsT>>(ints)) {
      directIntegrals = direct;
      contraction1 = std::make_shared<GTODirectTPIContraction<MatsT,IntsT>>(direct);
      contraction2 = std::make_shared<GTODirectTPIContraction<MatsT,IntsT>>(direct);
    } else if(auto incore = std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(ints)) {
      contraction1 = std::make_shared<InCore4indexTPIContraction<MatsT,IntsT>>(incore);
      contraction2 = std::make_shared<InCore4indexTPIContraction<MatsT,IntsT>>(incore);
    } else if(auto ri = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(ints)) {
      if (ri->isDistributed()) {
        contraction1 = std::make_shared<DistributedAsymmRITPIContraction<MatsT,IntsT>>(ri);
        contraction2 = std::make_shared<DistributedAsymmRITPIContraction<MatsT,IntsT>>(ri);
      } else {
        contraction1 = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(ri);
        contraction2 = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(ri);
      }
    } else {
      CErr("Unsupported interparticle integral type for " + label1 + "-" + label2);
    }
    contraction1->contractSecond = contractSecond;
    contraction2->contractSecond = !contractSecond;
    contraction1->isCross = true;
    contraction2->isCross = true;

    // Generate InterParticleFockBuilder for each subsystem
    auto fock1 = std::make_shared<InterParticleFockBuilder<MatsT,IntsT>>(ss1->aoints_->options_);
    auto fock2 = std::make_shared<InterParticleFockBuilder<MatsT,IntsT>>(ss2->aoints_->options_);

    fock1->setAux(ss2.get());
    fock1->setContraction(contraction1);
    fock2->setAux(ss1.get());
    fock2->setContraction(contraction2);

    const bool useBatchedDirect = directIntegrals and
                                  ss1->basisSet().basisType == REAL_GTO and
                                  ss2->basisSet().basisType == REAL_GTO;

    if(useBatchedDirect) {
      // Direct contractions done in a batched way to avoid overhead when number of subsystems get large
      const bool subsystemOrderMatchesIntegralOrder = not contractSecond;
      BasisSet& firstIntegralBasis = subsystemOrderMatchesIntegralOrder ?
        directIntegrals->basisSet() : directIntegrals->basisSet2();
      BasisSet& secondIntegralBasis = subsystemOrderMatchesIntegralOrder ?
        directIntegrals->basisSet2() : directIntegrals->basisSet();

      if(&firstIntegralBasis != &ss1->basisSet() or
         &secondIntegralBasis != &ss2->basisSet())
        CErr("Direct interparticle integral basis order is inconsistent for " +
          label1 + "-" + label2);

      // Each subsystem accumulates its batched Coulomb contribution into one persistent matrix
      auto setupBatchedCoulombMatrix = [&](const std::string& label, size_t nBasis) {
        auto [matrix, created] = batchedDirectCoulombMatrices.try_emplace(label, cqmatrix::Matrix<MatsT>(nBasis));
        if(created) matrix->second.clear();
      };

      setupBatchedCoulombMatrix(label1, ss1->basisSet().nBasis);
      setupBatchedCoulombMatrix(label2, ss2->basisSet().nBasis);
      batchedDirectPairs.push_back({label1, label2, directIntegrals, subsystemOrderMatchesIntegralOrder});
    } else {
      // In-core and RI fockbuilders are setup in a recursive chain.
      auto matrix1 = interCoulomb.at(label1).emplace(label2,
        cqmatrix::Matrix<MatsT>(ss1->basisSet().nBasis));
      auto matrix2 = interCoulomb.at(label2).emplace(label1,
        cqmatrix::Matrix<MatsT>(ss2->basisSet().nBasis));

      fock1->setOutput(&matrix1.first->second);
      fock1->setUpstream(fockBuilders.at(label1).back().get());
      fock2->setOutput(&matrix2.first->second);
      fock2->setUpstream(fockBuilders.at(label2).back().get());

      fockBuilders.at(label1).push_back(fock1);
      fockBuilders.at(label2).push_back(fock2);
      ss1->fockBuilder = fockBuilders.at(label1).back();
      ss2->fockBuilder = fockBuilders.at(label2).back();

      // Store in maps for look-up
      // Note that the direct pairs are not in these maps. They are in stored in batchedDirectPairs 
      interFockBuilders.at(label1).emplace(label2, fock1);
      interFockBuilders.at(label2).emplace(label1, fock2);
    }

    // Store the interparticle integrals
    interIntegrals.at(label1).emplace(label2, std::make_pair(contractSecond, ints));
    interIntegrals.at(label2).emplace(label1, std::make_pair(!contractSecond, ints));

  }


  template <typename MatsT, typename IntsT>
  bool MultiParticleSS<MatsT,IntsT>::isBatchedPair(const std::string& first,
    const std::string& second) const {

    for(const auto& interaction : batchedDirectPairs)
      if((interaction.firstSubsystemLabel == first and interaction.secondSubsystemLabel == second) or
         (interaction.firstSubsystemLabel == second and interaction.secondSubsystemLabel == first))
        return true;
    return false;

  }


  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::formBatchedDirectInterparticleCoulomb(
    const std::vector<std::string>& targets, bool increment) {

    if(batchedDirectPairs.empty()) return;

    if(increment)
      CErr("Incremental Fock build is not supported by the batched direct interparticle contraction.");

    // An empty target list means forming Interparticle Coulomb Fock for every subsystem Fock.
    const auto& labels = targets.empty() ? order_ : targets;
    const std::unordered_set<std::string> requested(labels.begin(), labels.end());

    std::unordered_set<std::string> formed;
    std::vector<DirectInterparticleJContraction<MatsT,IntsT>> directInterparticleJContractions;
    directInterparticleJContractions.reserve(batchedDirectPairs.size());

    for(const auto& pair : batchedDirectPairs) {
      const bool formFirstCoulomb  = requested.count(pair.firstSubsystemLabel);
      const bool formSecondCoulomb = requested.count(pair.secondSubsystemLabel);
      // Skip if not requested
      if(not formFirstCoulomb and not formSecondCoulomb) continue;

      auto& firstSubsystem  = subsystems.at(pair.firstSubsystemLabel);
      auto& secondSubsystem = subsystems.at(pair.secondSubsystemLabel);

      if(formFirstCoulomb)    formed.insert(pair.firstSubsystemLabel);
      if(formSecondCoulomb)   formed.insert(pair.secondSubsystemLabel);
      MatsT* firstCoulomb = formFirstCoulomb ?
        batchedDirectCoulombMatrices.at(pair.firstSubsystemLabel).pointer() : nullptr;
      MatsT* secondCoulomb = formSecondCoulomb ?
        batchedDirectCoulombMatrices.at(pair.secondSubsystemLabel).pointer() : nullptr;
      const double scale = 2. * firstSubsystem->particle.charge * secondSubsystem->particle.charge;

      directInterparticleJContractions.push_back({
        pair.integrals,
        pair.subsystemOrderMatchesIntegralOrder,
        firstSubsystem->onePDM->S().pointer(),
        secondSubsystem->onePDM->S().pointer(),
        firstCoulomb,
        secondCoulomb,
        scale,
        formFirstCoulomb,
        formSecondCoulomb
      });
    }

    for(const auto& label : formed)
      batchedDirectCoulombMatrices.at(label).clear();

    if(directInterparticleJContractions.empty()) return;

    time_point directStart;
    if(this->scfControls.printContractionTiming) directStart = tick();
    batchedDirectInterparticleJContraction.JContract(this->comm, directInterparticleJContractions);
    const double directSeconds = this->scfControls.printContractionTiming ? tock(directStart) : 0.;

    if(MPIRank(this->comm) != 0) return;
    for(const auto& label : formed) {
      auto& subsystem = subsystems.at(label);
      auto& coulomb = batchedDirectCoulombMatrices.at(label);
      *subsystem->twoeH += coulomb;
      *subsystem->fockMatrix += coulomb;
    }

    if(this->scfControls.printContractionTiming) {
      std::cout << "        Batched-Direct-Interparticle duration = "
                << directSeconds << " s\n";
    }
  }


  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::addGradientIntegrals(
    std::string label1, std::string label2,
    std::shared_ptr<GradInts<TwoPInts,IntsT>> ints, bool contractSecond) {

    if(!ints)
      CErr("Cannot add null gradient integrals for " + label1 + "-" + label2);

    // Batched-direct pairs evaluate their derivative integrals inside BatchedDirectInterparticleJContraction::GradJContract, 
    // The setup below is InterParticleFockBuilder is not needed
    if(isBatchedPair(label1,label2)) return;

    if(!interFockBuilders.at(label1).count(label2))
      CErr("Cannot add gradient integrals before adding interaction " + label1 + "-" + label2);
    if(gradInterInts.at(label1).count(label2))
      CErr("Gradient interaction already exists between " + label1 + " and " + label2);
  
    gradInterInts.at(label1).emplace(label2, std::make_pair(contractSecond, ints));
    gradInterInts.at(label2).emplace(label1, std::make_pair(!contractSecond, ints));
  
    interFockBuilders.at(label1).at(label2)->setGradientIntegrals(ints.get());
    interFockBuilders.at(label2).at(label1)->setGradientIntegrals(ints.get());

  };


  template <typename MatsT, typename IntsT>
  std::vector<double> MultiParticleSS<MatsT,IntsT>::getGrad(EMPerturbation& pert,
    bool equil, bool saveInts, double xHFX) {

    const size_t nAtoms = this->molecule().nAtoms;
    const size_t nGrad = 3*nAtoms;
    std::vector<double> gradient(nGrad);

    // Include the classical nuclear repulsion exactly once.
    for(size_t iGrad = 0; iGrad < nGrad; iGrad++)
      gradient[iGrad] = this->molecule().nucRepForce[iGrad/3][iGrad%3];

    // Form derivative integrals for each unique inter-particle interaction
    // (skipped entirely on the batched-direct path).
    for(size_t i = 0; i < order_.size(); i++) {
      for(size_t j = i + 1; j < order_.size(); j++) {
        const auto& first = order_[i];
        const auto& second = order_[j];

        if(isBatchedPair(first, second)) continue;
        if(!gradInterInts.at(first).count(second))
          CErr("Missing gradient integrals for " + first + "-" + second);

        auto& firstBasis = subsystems.at(first)->basisSet();
        auto& secondBasis = subsystems.at(second)->basisSet();
        const auto& gradientIntegrals = gradInterInts.at(first).at(second);

        HamiltonianOptions options;
        options.OneEScalarRelativity = false;
        if(firstBasis.basisType == COMPLEX_GIAO or secondBasis.basisType == COMPLEX_GIAO) {
          if(firstBasis.basisType != secondBasis.basisType)
            CErr("Basis types for inter-particle gradient integrals must match");
          options.basisType = COMPLEX_GIAO;
        }

        // contractSecond records which side of the stored asymmetric
        //   integral is contracted; compute derivatives in that same order.
        if(not gradientIntegrals.first)
          gradientIntegrals.second->computeAOInts(firstBasis,secondBasis,
            this->molecule(),pert,EP_ATTRACTION,options);
        else
          gradientIntegrals.second->computeAOInts(secondBasis,firstBasis,
            this->molecule(),pert,EP_ATTRACTION,options);
      }
    }

    XCTerms xcTerms;
    std::unordered_set<std::string> unifiedIntraXC;
    if(hasInterXC()) {
      xcTerms = buildXCTerms();
      for(const auto& part : xcTerms.subsystems)
        if(part.formIntraXC) unifiedIntraXC.insert(part.label);
    }

    // Add each subsystem gradient, removing its duplicate nuclear term.
    applyToEachLabeled([&](const std::string& label, SubSSPtr& ss) {
      std::vector<double> localGradient;
      auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss);

      // Intra-XC owned by the unified driver is added below. Independent KS
      //   subsystems continue to use their established gradient path.
      const bool ownsXCGrad = ks and (not ks->functionals.empty() or (ks->intParam.useGauXC and this->gauxcUtils));
      if(ownsXCGrad and not unifiedIntraXC.count(label))
        localGradient = ks->getGrad(pert,equil,saveInts,xHFX);
      else {
        double subsystemXHFX = 1.;
        const ExchCXX::HybCoeffs *rshCoefficients = nullptr;
        if(label == "E" and ks and this->gauxcUtils and this->gauxcUtils->isRangeSeparatedHybrid()) {
          rshCoefficients = &this->gauxcUtils->hybridCoefficients;
          subsystemXHFX = rshCoefficients->alpha;
        } else if(ks and not ks->functionals.empty()) {
          subsystemXHFX = ks->functionals.back()->xHFX;
        }
        localGradient = ss->SingleSlater<MatsT,IntsT>::getGrad(pert,equil,saveInts,subsystemXHFX);
        if(rshCoefficients) {
          auto shortRangeGradient = ks->getShortRangeExchangeGrad(pert,rshCoefficients->beta);
          std::transform(localGradient.begin(),localGradient.end(),shortRangeGradient.begin(),localGradient.begin(),std::plus<double>());
        }
      }

      for(size_t iGrad = 0; iGrad < nGrad; iGrad++)
        gradient[iGrad] += localGradient[iGrad] -
          this->molecule().nucRepForce[iGrad/3][iGrad%3];
    });

    std::vector<DirectInterparticleGradJContraction<MatsT,IntsT>> directInterparticleGradJContractions;
    directInterparticleGradJContractions.reserve(batchedDirectPairs.size());

    for(const auto& pair : batchedDirectPairs) {
      auto& firstSubsystem  = subsystems.at(pair.firstSubsystemLabel);
      auto& secondSubsystem = subsystems.at(pair.secondSubsystemLabel);

      directInterparticleGradJContractions.push_back({
        pair.integrals,
        pair.subsystemOrderMatchesIntegralOrder,
        firstSubsystem->onePDM->S().pointer(),
        secondSubsystem->onePDM->S().pointer(),
        2. * firstSubsystem->particle.charge * secondSubsystem->particle.charge});
    }

    std::vector<double> batchedGradContribution(nGrad, 0.);
    if(not directInterparticleGradJContractions.empty())
      batchedGradContribution = batchedDirectInterparticleJContraction.GradJContract(
        this->comm, nGrad, directInterparticleGradJContractions);

    for(size_t iGrad = 0; iGrad < nGrad; iGrad++)
      gradient[iGrad] += batchedGradContribution[iGrad];

    // All inter-XC terms and unified intra-XC terms share one grid pass.
    if(not xcTerms.subsystems.empty()) {
      auto xcGradient = formXCGradient(pert,xcTerms);
      std::transform(gradient.begin(),gradient.end(),xcGradient.begin(),
        gradient.begin(),std::plus<double>());
    }

#ifdef CQ_HAS_D3
    if(this->d3Utils) {
      this->d3Utils->evaluate(this->molecule());
      for(size_t iGrad = 0; iGrad < nGrad; iGrad++)
        gradient[iGrad] += this->d3Utils->result().gradient[iGrad];
    }
#endif

    if(MPIRank(this->comm) == 0) {
      std::cout << "Total MultiParticleSS Gradient:" << std::endl;
      std::cout << std::setprecision(12);
      for(size_t iAt = 0; iAt < nAtoms; iAt++) {
        std::cout << " Gradient@I = " << iAt << ":";
        for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
          std::cout << "  " << gradient[3*iAt+iXYZ];
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    return gradient;

  };

  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::printSetup(std::ostream& out) {

    if(order_.size() != subsystems.size())
      CErr("MultiParticleSS ordering does not match subsystem storage.");

    auto findBatchedPair = [&](const std::string& first,
                               const std::string& second)
      -> const BatchedDirectPair* {
      for(const auto& interaction : batchedDirectPairs) {
        const bool sameOrder = interaction.firstSubsystemLabel == first and
          interaction.secondSubsystemLabel == second;
        const bool reverseOrder = interaction.firstSubsystemLabel == second and
          interaction.secondSubsystemLabel == first;
        if(sameOrder or reverseOrder) return &interaction;
      }
      return nullptr;
    };

    std::unordered_set<std::string> seen;
    size_t nDirectedChainBuilders = 0;
    size_t nDirectedBatchedPartners = 0;
    size_t nFockBuilderChainNodes = 0;

    out << "\n";
    out << "============================================================\n";
    out << " MultiParticleSS Setup\n";
    out << "============================================================\n";
    out << " Number of subsystems: " << order_.size() << "\n\n";

    for(const auto& label : order_) {

      if(!seen.emplace(label).second)
        CErr("Duplicate label in MultiParticleSS ordering: " + label);

      if(!subsystems.count(label))
        CErr("Ordered subsystem does not exist: " + label);

      if(!fockBuilders.count(label) ||
        !interIntegrals.count(label) ||
        !interCoulomb.count(label) ||
        !interFockBuilders.count(label))
        CErr("Incomplete storage for subsystem: " + label);

      const auto& ss = subsystems.at(label);
      const auto& chain = fockBuilders.at(label);
      const auto& interactions = interFockBuilders.at(label);

      if(chain.empty())
        CErr("Empty FockBuilder chain for subsystem: " + label);

      // Batched-direct partners of this subsystem, in registration order.
      std::vector<std::string> batchedPartners;
      for(const auto& interaction : batchedDirectPairs) {
        if(interaction.firstSubsystemLabel == label)
          batchedPartners.push_back(interaction.secondSubsystemLabel);
        else if(interaction.secondSubsystemLabel == label)
          batchedPartners.push_back(interaction.firstSubsystemLabel);
      }
      const size_t batchedCount = batchedPartners.size();

      // Check the size of chain builders (for incore/RI interactions)
      if(chain.size() != interactions.size() + 1)
        CErr("Incorrect FockBuilder chain length for subsystem: " + label);

      // Check that every registered partner without a chain builder MUST be a batched-direct pair, 
      // and the number of such partners must match batchedDirectPairs. 
      size_t batchedFromRegistry = 0;
      for(const auto& registered : interIntegrals.at(label)) {
        if(interactions.count(registered.first)) continue;
        if(!findBatchedPair(label, registered.first))
          CErr("Interaction " + label + "-" + registered.first +" has neither an InterParticleFockBuilder nor a batched direct pair");
        ++batchedFromRegistry;
      }

      if(batchedFromRegistry != batchedCount)
        CErr("Batched direct partner count does not match the registered interactions for subsystem: " + label);

      if(batchedCount and not batchedDirectCoulombMatrices.count(label))
        CErr("Missing batched direct Coulomb matrix for subsystem: " + label);

      if(ss->fockBuilder.get() != chain.back().get())
        CErr("Subsystem does not point to the end of its FockBuilder chain: " +
            label);

      for(size_t i = 1; i < chain.size(); ++i) {

        auto interBuilder =
          dynamic_cast<InterParticleFockBuilder<MatsT,IntsT>*>(
            chain[i].get());

        if(!interBuilder)
          CErr("Non-interparticle builder found inside interaction chain for " +
              label);

        if(interBuilder->getUpstream() != chain[i - 1].get())
          CErr("Broken upstream FockBuilder chain for subsystem: " + label);
      }

      out << " Subsystem " << label << "\n";
      out << "   Basis functions: "
          << ss->basisSet().nBasis << "\n";
      out << "   Particle charge: "
          << ss->particle.charge << "\n";
      out << "   Interactions: "
          << interIntegrals.at(label).size()
          << " (recursive chain: " << interactions.size()
          << ", batched direct: " << batchedCount << ")\n";
      out << "   Recursive FockBuilder chain: base";

      for(size_t i = 1; i < chain.size(); ++i) {

        std::string partner = "<unknown>";

        for(const auto& interaction : interactions) {
          if(interaction.second.get() == chain[i].get()) {
            partner = interaction.first;
            break;
          }
        }

        out << " -> " << partner;
      }

      out << "\n";

      out << "   Batched direct partners: ";
      if(batchedPartners.empty())
        out << "none";
      else
        for(size_t i = 0; i < batchedPartners.size(); ++i)
          out << (i ? ", " : "") << label << " <-> " << batchedPartners[i];

      out << "\n\n";

      nDirectedChainBuilders += interactions.size();
      nDirectedBatchedPartners += batchedCount;
      nFockBuilderChainNodes += chain.size();
    }

    size_t nPairInteractions = 0;
    size_t nBatchedPairs = 0;

    out << " Pair interactions\n";

    for(size_t i = 0; i < order_.size(); ++i) {
      for(size_t j = i + 1; j < order_.size(); ++j) {

        const auto& label1 = order_[i];
        const auto& label2 = order_[j];

        const bool has12 =
          interIntegrals.at(label1).count(label2);

        const bool has21 =
          interIntegrals.at(label2).count(label1);

        if(has12 != has21)
          CErr("Interaction exists in only one direction for " +
              label1 + "-" + label2);

        if(!has12)
          continue;

        const auto& interaction12 =
          interIntegrals.at(label1).at(label2);

        const auto& interaction21 =
          interIntegrals.at(label2).at(label1);

        if(interaction12.second != interaction21.second)
          CErr("Directional interaction entries do not share integrals for " +
              label1 + "-" + label2);

        if(interaction12.first == interaction21.first)
          CErr("Directional contractSecond flags are inconsistent for " +
              label1 + "-" + label2);

        const auto* batchedPair = findBatchedPair(label1, label2);
        if(batchedPair) {
          if(interCoulomb.at(label1).count(label2) or
             interCoulomb.at(label2).count(label1))
            CErr("Batched direct interaction has pair-local Coulomb matrices for " + label1 + "-" + label2);
          if(batchedPair->integrals != interaction12.second)
            CErr("Batched direct interaction does not share the registered integrals for " + label1 + "-" + label2);
          ++nBatchedPairs;
        } else if(!interCoulomb.at(label1).count(label2) or !interCoulomb.at(label2).count(label1)) {
          CErr("Missing Coulomb matrix for " + label1 + "-" + label2);
        }

        if(batchedPair) {
          if(interFockBuilders.at(label1).count(label2) || interFockBuilders.at(label2).count(label1))
            CErr("Batched direct interaction has a registered InterParticleFockBuilder for " +
                label1 + "-" + label2);
        } else if(!interFockBuilders.at(label1).count(label2) || !interFockBuilders.at(label2).count(label1)) {
          CErr("Missing InterParticleFockBuilder for " + label1 + "-" + label2);
        }

        out << "   " << label1 << " <-> " << label2
            << "   contractSecond(" << label1 << ") = "
            << std::boolalpha << interaction12.first
            << ", contractSecond(" << label2 << ") = "
            << interaction21.first << ", executor = "
            << (batchedPair ? "batched-direct" : "pair-chain") << "\n";

        ++nPairInteractions;
      }
    }

    const size_t nChainPairs = nPairInteractions - nBatchedPairs;

    if(nDirectedChainBuilders != 2 * nChainPairs)
      CErr("Directional interaction count does not match chain pair count.");
    if(nDirectedBatchedPartners != 2 * nBatchedPairs)
      CErr("Directional batched partner count does not match batched pair count.");
    if(nBatchedPairs != batchedDirectPairs.size())
      CErr("Batched direct interaction count does not match pair count.");

    const size_t nPossiblePairs =
      order_.size() * (order_.size() - 1) / 2;

    out << "\n";
    out << " Pair interactions: "
        << nPairInteractions << " / "
        << nPossiblePairs << " possible\n";

    out << " Recursive-chain interparticle pairs: "
        << nChainPairs << "\n";

    out << " Batched direct interparticle pairs: "
        << nBatchedPairs << "\n";

    out << " Directional interparticle builders: "
        << nDirectedChainBuilders << "\n";

    out << " Recursive FockBuilder chain nodes: "
        << nFockBuilderChainNodes << "\n";

    out << " Total FockBuilder objects: "
        << subsystems.size() + nDirectedChainBuilders << "\n";

    out << " Setup validation: PASSED\n";
    out << "============================================================\n\n";
  }

} // namespace ChronusQ
