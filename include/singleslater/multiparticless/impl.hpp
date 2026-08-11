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
    
    // Create contractions for each subsystem
    if(auto direct = std::dynamic_pointer_cast<DirectTPI<IntsT>>(ints)) {
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

    // Generate the output coulomb matrix for each subsystem
    auto matrix1 = interCoulomb.at(label1).emplace(label2, cqmatrix::Matrix<MatsT>(ss1->basisSet().nBasis));
    auto matrix2 = interCoulomb.at(label2).emplace(label1, cqmatrix::Matrix<MatsT>(ss2->basisSet().nBasis));

    // Generate InterParticleFockBuilder for each subsystem
    auto fock1 = std::make_shared<InterParticleFockBuilder<MatsT,IntsT>>(ss1->aoints_->options_);
    auto fock2 = std::make_shared<InterParticleFockBuilder<MatsT,IntsT>>(ss2->aoints_->options_);

    fock1->setAux(ss2.get());
    fock1->setOutput(&matrix1.first->second);
    fock1->setContraction(contraction1);
    fock1->setUpstream(fockBuilders[label1].back().get());

    fock2->setAux(ss1.get());
    fock2->setOutput(&matrix2.first->second);
    fock2->setContraction(contraction2);
    fock2->setUpstream(fockBuilders[label2].back().get());

    // Store the interparticle integrals and fockBuilders 
    interIntegrals.at(label1).emplace(label2, std::make_pair(contractSecond, ints));
    interIntegrals.at(label2).emplace(label1, std::make_pair(!contractSecond, ints));
    interFockBuilders.at(label1).emplace(label2, fock1);
    interFockBuilders.at(label2).emplace(label1, fock2);

    // Add the newly created interParticleFockBuilders to the list of fockBuilders for each subsystem
    fockBuilders.at(label1).push_back(fock1);
    fockBuilders.at(label2).push_back(fock2);

    // Set the fockBuilder for each subsystem to the newly created interParticleFockBuilders
    ss1->fockBuilder = fockBuilders[label1].back();
    ss2->fockBuilder = fockBuilders[label2].back();

  }


  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::addGradientIntegrals(
    std::string label1, std::string label2,
    std::shared_ptr<GradInts<TwoPInts,IntsT>> ints, bool contractSecond) {

    if(!ints)
      CErr("Cannot add null gradient integrals for " + label1 + "-" + label2);
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

    // Form derivative integrals for each unique inter-particle interaction.
    for(size_t i = 0; i < order_.size(); i++) {
      for(size_t j = i + 1; j < order_.size(); j++) {
        const auto& first = order_[i];
        const auto& second = order_[j];

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
  void MultiParticleSS<MatsT,IntsT>::buildOrbitalModifierOptions() {
    
    CErr("buildOrbitalModifierOptions not implemented for MultiParticleSS");

  }



  template <typename MatsT, typename IntsT>
  void MultiParticleSS<MatsT,IntsT>::printSetup(std::ostream& out) {

    if(order_.size() != subsystems.size())
      CErr("MultiParticleSS ordering does not match subsystem storage.");

    std::unordered_set<std::string> seen;
    size_t nDirectedInteractions = 0;
    size_t nTotalFockBuilders = 0;

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

      if(chain.size() != interactions.size() + 1)
        CErr("Incorrect FockBuilder chain length for subsystem: " + label);

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
          << interactions.size() << "\n";
      out << "   FockBuilder chain: base";

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

      out << "\n\n";

      nDirectedInteractions += interactions.size();
      nTotalFockBuilders += chain.size();
    }

    size_t nPairInteractions = 0;

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

        if(!interCoulomb.at(label1).count(label2) ||
          !interCoulomb.at(label2).count(label1))
          CErr("Missing Coulomb matrix for " + label1 + "-" + label2);

        if(!interFockBuilders.at(label1).count(label2) ||
          !interFockBuilders.at(label2).count(label1))
          CErr("Missing InterParticleFockBuilder for " +
              label1 + "-" + label2);

        out << "   " << label1 << " <-> " << label2
            << "   contractSecond(" << label1 << ") = "
            << std::boolalpha << interaction12.first
            << ", contractSecond(" << label2 << ") = "
            << interaction21.first << "\n";

        ++nPairInteractions;
      }
    }

    if(nDirectedInteractions != 2 * nPairInteractions)
      CErr("Directional interaction count does not match pair count.");

    const size_t nPossiblePairs =
      order_.size() * (order_.size() - 1) / 2;

    out << "\n";
    out << " Pair interactions: "
        << nPairInteractions << " / "
        << nPossiblePairs << " possible\n";

    out << " Directional interparticle builders: "
        << nDirectedInteractions << "\n";

    out << " Total FockBuilders: "
        << nTotalFockBuilders << "\n";

    out << " Setup validation: PASSED\n";
    out << "============================================================\n\n";
  }

} // namespace ChronusQ
