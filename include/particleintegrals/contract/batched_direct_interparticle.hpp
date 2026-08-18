/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#pragma once

#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <util/mpi.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace ChronusQ {

  /**
   * Runtime data for one unordered, distinguishable-particle Coulomb contractiom.
   */
  template <typename MatsT, typename IntsT>
  struct DirectInterparticleJContraction {
    std::shared_ptr<DirectTPI<IntsT>> integrals;
    bool subsystemOrderMatchesIntegralOrder = true; // Whether the (first, second) subsystem order matches the (basisSet(), basisSet2()) order of the shared DirectTPI.

    const MatsT* firstDensity = nullptr;
    const MatsT* secondDensity = nullptr;

    MatsT* firstCoulomb = nullptr;
    MatsT* secondCoulomb = nullptr;

    double interactionScale = 1.;                   // Scaling prefactor to be determined by subsystem charge
    bool formFirstCoulomb = false;                  // Compute 1st subsystem's Couloumb matrix
    bool formSecondCoulomb = false;                 // Compute 2nd subsystem's Couloumb matrix
  };

  /**
   * Runtime data for one gradient Coulomb interaction.
   */
  template <typename MatsT, typename IntsT>
  struct DirectInterparticleGradJContraction {
    std::shared_ptr<DirectTPI<IntsT>> integrals;
    bool subsystemOrderMatchesIntegralOrder = true;

    const MatsT* firstDensity = nullptr;
    const MatsT* secondDensity = nullptr;

    double interactionScale = 1.;
  };

  /**
   * Batched, direct Coulomb contraction for distinguishable-particle pairs.
   */
  template <typename MatsT, typename IntsT>
  class BatchedDirectInterparticleJContraction {
    public:

      /**
       * Evaluate direct Coulomb interactions as one batched contraction.
       * Each interaction is traversed once and can digest both subsystem densities.
       * (compute integral once, both contract with both densities to compute both J contributions)
       *
       */
      void JContract(
        MPI_Comm comm,
        const std::vector<DirectInterparticleJContraction<MatsT,IntsT>>& directInterparticleJContractions);

      /**
       * Batched derivative Coulomb contraction, evaluated in trace mode.
       *   dE/dR_I += 0.5 * interactionScale * sum_{mu nu, lambda sigma} Re D_first(mu,nu) 
       *                  * d^I(mu nu | lambda sigma) Re D_second(lambda,sigma)
       *
       * Returns the summed molecular gradient contribution of every listed interaction
       */
      std::vector<double> GradJContract(
        MPI_Comm comm,
        size_t gradientComponentCount,
        const std::vector<DirectInterparticleGradJContraction<MatsT,IntsT>>& directInterparticleGradJContractions);
  };

} // namespace ChronusQ
