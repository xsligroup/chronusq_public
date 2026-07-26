/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#pragma once

#include <singleslater/kohnsham.hpp>

#include <unordered_set>

namespace ChronusQ {

  // Build the XC terms for all subsystems.
  template <typename MatsT, typename IntsT>
  typename MultiParticleSS<MatsT,IntsT>::XCTerms
  MultiParticleSS<MatsT,IntsT>::buildXCTerms() const {

    return buildXCTermsFor(order_);
  }

  /**
   *  Build the density inputs, functionals and VXC outputs for one shared-grid
   *  evaluation.
   */
  template <typename MatsT, typename IntsT>
  typename MultiParticleSS<MatsT,IntsT>::XCTerms
  MultiParticleSS<MatsT,IntsT>::buildXCTermsFor(
      const std::vector<std::string>& targets) const {

    XCTerms xcTerms;

    // Start from the requested output Fock matrices.
    std::unordered_set<std::string> targetSet;
    for(const auto& label : targets) {
      if( !subsystems.count(label) )
        CErr("Unknown MultiParticleSS XC target " + label);
      targetSet.insert(label);
    }
    // Before the first targeted build, include every functional for energy but
    //   still form VXC only for the requested targets. Subsequent builds can
    //   retain energies whose subsystem densities did not change.
    const bool updateAllEnergies = needsFullXCEnergyUpdate;

    std::unordered_map<std::string,size_t> compactIndex;
    auto addSubsystem = [&](const std::string& label) {
      auto iter = compactIndex.find(label);
      if( iter != compactIndex.end() ) return iter->second;

      size_t index = xcTerms.subsystems.size();
      compactIndex.emplace(label, index);
      xcTerms.subsystems.push_back({label, targetSet.count(label) != 0, false});
      return index;
    };

    // Targets activate their inter-particle functionals.
    for(size_t i = 0; i < interFunctionals.size(); ++i) {
      const auto& pair = interFunctionals[i];
      if( pair.functionals.empty() ) continue;

      bool containsTarget = targetSet.count(pair.first) or targetSet.count(pair.second);
      if( not containsTarget and not updateAllEnergies ) continue;

      size_t first = addSubsystem(pair.first);
      size_t second = addSubsystem(pair.second);

      // In-house EPC energy is integrated from the electronic side.
      if( subsystems.at(pair.first)->particle.charge >= 0 and
          subsystems.at(pair.second)->particle.charge >= 0 )
        CErr("Inter-particle XC energy requires one electronic subsystem");

      xcTerms.interTerms.push_back({i, first, second});
    }

    // Intra-XC belongs only to its subsystem target after the full energy update.
    for(size_t p = 0; p < order_.size(); ++p) {
      const auto& label = order_[p];
      auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(subsystems.at(label));
      if( not ks or ks->doVXC_ ) continue;

      bool hasIntraXC = not ks->functionals.empty();
      if( intParam.useGauXC ) {
        if( not this->gauxcUtils )
          CErr("GauXC MultiParticleSS requested without GauXCUtils");
        const auto& intra = this->gauxcUtils->multiparticle_functional_spec.intra_functionals;
        hasIntraXC = p < intra.size() and not intra[p].empty();
      }
      if( not hasIntraXC ) continue;

      bool isTarget = targetSet.count(label);
      if( not isTarget and not updateAllEnergies ) continue;

      auto& subsystem = xcTerms.subsystems[addSubsystem(label)];
      subsystem.formIntraXC = true;
    }

    return xcTerms;
  }

} // namespace ChronusQ
