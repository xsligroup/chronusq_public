#pragma once

#include <memory>
#include <string>
#include <integrals.hpp>
#include <singleslater/base.hpp>

namespace ChronusQ {

  // Helper struct to store information about a quantum subsystem
  struct QuantumSubsystem {
    std::string inputLabel;        ///< Quantum subsystem input label (quantum subsystem symbol, specified by user in the GEOM section)
    std::string label;             ///< Quantum subsystem label (runtime label, different from input label in the case of distinguished particles within the same quantum subsystem)

    std::string qmSection;         ///< Quantum subsystem QM Reference Section
    std::string basisSection;      ///< Quantum subsystem Basis Section
    std::string dfbasisSection;    ///< Quantum subsystem DF Basis Section
    std::string guessBasisSection; ///< Quantum subsystem Guess Basis Section
    std::string intsSection;       ///< Quantum subsystem Integrals Section (optional)

    std::shared_ptr<BasisSet> basis = nullptr;
    std::shared_ptr<BasisSet> dfbasis = nullptr;
    std::shared_ptr<BasisSet> guessBasis = nullptr;

    IntegralOptions integralOptions;
    std::shared_ptr<IntegralsBase> integrals = nullptr;

    size_t particleIndex = 0;      ///< Particle index within the quantum subsystem (for distinguished particles within the same quantum subsystem)
    Particle particle;             ///< Particle properties for the quantum subsystem
    std::vector<size_t> atomIndices; ///< Molecule atom indices carrying this quantum subsystem's basis
    size_t nQuantumParticles = 0;   ///< Number of quantum particles in the quantum subsystem // TODOAL: maybe simplify all these variables later?
    SingleSlaterOptions ssOptions;  ///< SingleSlater options for the quantum subsystem
  };
  // Helper struct to store information about pair interaction between two quantum subsystems
  struct QuantumPairInteraction {
    std::string labelA;       ///< First quantum subsystem label
    std::string labelB;       ///< Second quantum subsystem label
    std::string inputLabelA;  ///< First quantum subsystem input label
    std::string inputLabelB;  ///< Second quantum subsystem input label

    std::string intsSection;  ///< Pair integral section (optional)

    IntegralOptions integralOptions;
    std::shared_ptr<IntegralsBase> integrals = nullptr;
  };

  std::vector<QuantumSubsystem> buildQuantumSubsystems(
    std::ostream& out, CQInputFile& input, Molecule& mol, const std::string& inputLabel);

  QuantumPairInteraction buildQuantumPairInteraction(std::ostream& out, CQInputFile& input, Molecule& mol, const QuantumSubsystem& sysA, const QuantumSubsystem& sysB);

  void printQuantumSubsystem(std::ostream& out, CQInputFile& input, const QuantumSubsystem& sys);
  void printQuantumSetup(std::ostream& out, CQInputFile& input, const std::vector<QuantumSubsystem>& quantumSubsystems, const std::vector<QuantumPairInteraction>& quantumPairInteractions);

  // Build 2-body integrals for all quantum subsystems and pair interactions
  void buildTwoBodyIntegrals(
      std::ostream& out, Molecule& mol,
      std::vector<QuantumSubsystem>& quantumSubsystems,
      std::vector<QuantumPairInteraction>& quantumPairInteractions);

  Particle getParticleForSubsystem(CQInputFile& input, const Molecule& mol, const QuantumSubsystem& sys);

}
