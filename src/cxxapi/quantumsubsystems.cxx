#include <quantumsubsystems.hpp>
#include <cxxapi/input.hpp>
#include <cxxapi/options.hpp>
#include <basisset.hpp>
#include <basisset/remove_linear_dep_shells.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <unordered_map>

namespace ChronusQ {

  std::vector<QuantumSubsystem> buildQuantumSubsystems(
    std::ostream& out, CQInputFile& input, Molecule& mol, const std::string& inputLabel) {
  
    // Helper function to ensure backward compatibility with legacy NEO input file format
    auto pick = [&](const std::string& canonical, const std::string& legacy) {
      return inputLabel == "QP" and not input.containsSection(canonical) and input.containsSection(legacy) ? legacy : canonical;
    };
  
    const std::string qmSection = pick(inputLabel + "QM","PROTQM");
  
    if(not input.containsSection(qmSection))
      CErr("Missing quantum subsystem section [" + qmSection + "]");
  
    auto invalidKeywords = CQQUANTUMSUBSYSTEMQM_VALID(input.getSection(qmSection));
    printInvalidKeys(invalidKeywords,qmSection);
  
    const bool distinguishable = input.containsData(qmSection + "/DISTINGUISHABLE") ?
      input.getData<bool>(qmSection + "/DISTINGUISHABLE") : false;

    std::vector<size_t> qAtomIndices;
    std::vector<size_t> basisAtomIndices;
    bool hasGhostCenters = false;
    for(size_t iAtm = 0; iAtm < mol.atoms.size(); ++iAtm) {
      const auto& atom = mol.atoms[iAtm];
      const bool belongsToSubsystem = atom.quantum and atom.quantumLabel == inputLabel;
      const bool isQuantumGhost = atom.atomicNumber == 0 and atom.atomicMass == atomicReference["GH-0"].atomicMass;
      if(belongsToSubsystem)
        qAtomIndices.push_back(iAtm);
      // GH carries quantum basis functions without adding a quantum particle.
      if(belongsToSubsystem or isQuantumGhost)
        basisAtomIndices.push_back(iAtm);
      if(isQuantumGhost)
        hasGhostCenters = true;
    }

    if(distinguishable and hasGhostCenters)
      CErr("Ghost centers are not supported for distinguishable quantum subsystem " + inputLabel);
  
    // Count number of particles specified in this quantum subsystem
    const size_t nPart = qAtomIndices.size();
    // Number of quantum subsystems to actually build (1 for indistinguishable particles, nPart for distinguishable particles)
    const size_t nSys = distinguishable and nPart > 1 ? nPart : 1;
  
    std::vector<QuantumSubsystem> systems;
    systems.reserve(nSys);
  
    for(size_t i = 0; i < nSys; ++i) {
      QuantumSubsystem sys;
  
      sys.inputLabel = inputLabel;
      sys.label = distinguishable ? inputLabel + std::to_string(i) : inputLabel;
      sys.particleIndex = distinguishable ? i : 0;
      sys.nQuantumParticles = distinguishable ? 1 : nPart;

      sys.qmSection = qmSection;
      sys.basisSection = pick(inputLabel + "BASIS","PBASIS");
      sys.dfbasisSection = pick(inputLabel + "DFBASIS","PDFBASIS");
      sys.guessBasisSection = pick(inputLabel + "GUESSBASIS","PGUESSBASIS");
      sys.intsSection = pick(inputLabel + "INTS","PINTS");
  
      if(not input.containsSection(sys.basisSection))
        CErr("Missing quantum subsystem basis section [" + sys.basisSection + "]");
  
      sys.particle = getParticleForSubsystem(input,mol,sys);
  
      // Distinguishable particles receive separate one-center basis objects.
      sys.atomIndices = distinguishable ? std::vector<size_t>{qAtomIndices[i]} : basisAtomIndices;
      sys.basis = CQBasisSetOptions(out,input,mol,sys.basisSection,sys.atomIndices);
      read_option_and_remove_linear_dependency(input,*sys.basis,mol,out);
  
      sys.dfbasis = input.containsSection(sys.dfbasisSection) ?
        CQBasisSetOptions(out,input,mol,sys.dfbasisSection,sys.atomIndices) : nullptr;
  
      sys.guessBasis = input.containsSection(sys.guessBasisSection) ?
        CQBasisSetOptions(out,input,mol,sys.guessBasisSection,sys.atomIndices) : nullptr;
  
      sys.integralOptions = getIntegralOptions(out,input,sys.basis,sys.dfbasis,sys.guessBasis,sys.intsSection);
  
      systems.push_back(std::move(sys));
    }
  
    return systems;
  }


  QuantumPairInteraction buildQuantumPairInteraction(
    std::ostream& out, CQInputFile& input, Molecule& mol, const QuantumSubsystem& sysA, const QuantumSubsystem& sysB) {

    QuantumPairInteraction pair;

    pair.labelA = sysA.label;
    pair.labelB = sysB.label;
    pair.inputLabelA = sysA.inputLabel;
    pair.inputLabelB = sysB.inputLabel;

    const std::string sectionAB = pair.inputLabelA + pair.inputLabelB + "INTS";
    const std::string sectionBA = pair.inputLabelB + pair.inputLabelA + "INTS";

    if(input.containsSection(sectionAB))
      pair.intsSection = sectionAB;
    else if(input.containsSection(sectionBA))
      pair.intsSection = sectionBA;
    else if(((pair.inputLabelA == "E" and pair.inputLabelB == "QP") or
             (pair.inputLabelB == "QP" and pair.inputLabelA == "E")) and
            input.containsSection("EPINTS"))
      pair.intsSection = "EPINTS";
    else
      pair.intsSection = sectionAB;

    pair.integralOptions = getIntegralOptions(out,input,nullptr,nullptr,nullptr,pair.intsSection, true);

    return pair;
  }


  void printQuantumSubsystem(std::ostream& out, CQInputFile& input, const QuantumSubsystem& sys) {

    auto yn    = [](bool b) { return b ? "yes" : "no"; };
    auto ptr   = [](const auto& p) { return p ? "built" : "null"; };
    auto sec   = [&](const std::string& s) { return std::string("[") + s + ", exists=" + yn(input.containsSection(s)) + "]"; };
    auto idx   = sys.inputLabel == sys.label ? std::string("-") : std::to_string(sys.particleIndex);
  
    out << "    " << sys.label << "  inputLabel=" << sys.inputLabel
        << ", particleIndex=" << idx
        << ", charge=" << sys.particle.charge
        << ", mass=" << sys.particle.mass << " me\n";
  
    out << "      sections: QM=" << sec(sys.qmSection)
        << ", BASIS=" << sec(sys.basisSection)
        << ", DFBASIS=" << sec(sys.dfbasisSection)
        << ", GUESSBASIS=" << sec(sys.guessBasisSection)
        << ", INTS=" << sec(sys.intsSection) << "\n";
  
    out << "      objects : basis=" << ptr(sys.basis)
        << ", dfbasis=" << ptr(sys.dfbasis)
        << ", guessBasis=" << ptr(sys.guessBasis)
        << ", integrals=" << ptr(sys.integrals) << "\n";
  }

  void printQuantumSetup(
    std::ostream& out, CQInputFile& input,
    const std::vector<QuantumSubsystem>& quantumSubsystems,
    const std::vector<QuantumPairInteraction>& quantumPairInteractions) {
  
    auto yn  = [](bool b) { return b ? "yes" : "no"; };
    auto ptr = [](const auto& p) { return p ? "built" : "null"; };
  
    auto findSys = [&](const std::string& label) -> const QuantumSubsystem* {
      for(const auto& sys : quantumSubsystems)
        if(sys.label == label) return &sys;
      return nullptr;
    };
  
    if(not quantumSubsystems.empty()) {
      out << "\n  Multicomponent quantum setup:\n";
      out << "\n  Quantum subsystems:\n";
      for(const auto& sys : quantumSubsystems)
        printQuantumSubsystem(out,input,sys);
    }
  
    if(not quantumPairInteractions.empty()) {
      out << "\n  Quantum pair interactions:\n";
  
      for(const auto& pair : quantumPairInteractions) {
        const auto* sysA = findSys(pair.labelA);
        const auto* sysB = findSys(pair.labelB);
  
        out << "    " << pair.labelA << "-" << pair.labelB
            << "  inputLabels=" << pair.inputLabelA << "-" << pair.inputLabelB
            << ", INTS=[" << pair.intsSection << ", exists=" << yn(input.containsSection(pair.intsSection)) << "]";
  
        if(sysA and sysB)
          out << ", chargeProduct=" << sysA->particle.charge * sysB->particle.charge;
  
        out << ", integrals=" << ptr(pair.integrals) << "\n";
      }
    }
  
    out << std::endl;
  }


  void buildTwoBodyIntegrals(
    std::ostream& out, Molecule& mol,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions) {
  
    out << BannerTop << std::endl;
    out << "\n  Building 2-body integral objects for quantum subsystems and pair interactions\n\n";
  
    std::unordered_map<std::string, QuantumSubsystem*> subsystemMap;
  
    for(auto& sys : quantumSubsystems)
      subsystemMap[sys.label] = &sys;
  
    // 1. Build all symmetric integral objects.
    for(auto& sys : quantumSubsystems)
      sys.integrals = sys.integralOptions.buildSymmIntegral(out, mol, sys.basis, sys.dfbasis, sys.label);
  
    // 2. Build all pair/asymmetric interaction integral objects.
    for(auto& pair : quantumPairInteractions) {
      auto* sysA = subsystemMap.at(pair.labelA);
      auto* sysB = subsystemMap.at(pair.labelB);
      pair.integrals = pair.integralOptions.buildAsymmIntegral(out,mol,
          sysA->basis, sysA->dfbasis, sysB->basis,
          sysA->integralOptions, sysB->integralOptions,
          sysA->integrals, sysB->integrals,
          pair.labelA, pair.labelB
      );
      // Set the particles for cross integrals, important for sign flip in GIAO.
      pair.integrals->options_.particle  = sysA->particle; // Bra particle
      pair.integrals->options_.particle2 = sysB->particle; // Ket particle
    }
  
    out << BannerEnd << std::endl;
  }

  Particle getParticleForSubsystem(
    CQInputFile &input, const Molecule &mol, const QuantumSubsystem &sys ) {
  
    // Electronic subsystem
    if(sys.label == "E")  return Particle{-1.0, 1.0};

    Particle p;
    // Non-electronic quantum subsystem: default from geometry label
    bool found = false;
    double massAMU = 0.0;
    double charge = 0.0;
  
    for(const auto &atom : mol.atoms) {
      if(!atom.quantum) continue;
      if(atom.quantumLabel != sys.inputLabel) continue;
      massAMU = atom.atomicMass;
      charge  = static_cast<double>(atom.atomicNumber);
      found = true;
      break;
    }
  
    if(!found)
      CErr("Cannot infer particle properties for quantum subsystem " + sys.label);
  
    double mass = 0.0;
  
    if(massAMU == atomicReference["H-1"].atomicMass)
      mass = ProtMassPerE;
    else if(massAMU == atomicReference["H-2"].atomicMass)
      mass = DeutMassPerE;
    else if(massAMU == atomicReference["H-3"].atomicMass)
      mass = TritMassPerE;
    else
      mass = massAMU * AUPerAMU - charge;
  
    p = Particle{charge, mass};
    
    // User override from QM section
    const std::string sec = sys.qmSection;
  
    if(input.containsData(sec + "/PARTICLECHARGE"))
      p.charge = input.getData<double>(sec + "/PARTICLECHARGE");

    if(input.containsData(sec + "/PARTICLEMASS"))
      p.mass = input.getData<double>(sec + "/PARTICLEMASS");
  
    return p;
  }

}
