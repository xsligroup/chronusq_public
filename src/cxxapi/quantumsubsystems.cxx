#include <quantumsubsystems.hpp>
#include <cxxapi/input.hpp>
#include <cxxapi/options.hpp>
#include <basisset.hpp>
#include <basisset/remove_linear_dep_shells.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <unordered_map>
#include <cctype>
#include <algorithm>

namespace ChronusQ {

  // Strip the trailing particle index to get the base label that distinguishable particles share for input-section lookup. 
  //   For example, if user has explicit "QP0","QP1","QP2" labels, they are supposed to provide "QP0QM", "QP0BASIS", "QP1QM", 
  //   "QP1BASIS" ... which could be a lot of work preparing the input.
  //   This helper allows extraction of base labels "QP0","QP1","QP2" -> "QP" such that they all read [QP*]
  static std::string baseLabel(const std::string& label) {
    size_t end = label.size();
    while(end > 0 and std::isdigit(static_cast<unsigned char>(label[end-1]))) --end;
    return label.substr(0, end);
  }

  static bool isGhostAtom(const Atom& atom) {
    return atom.atomicNumber == 0 and atom.atomicMass == atomicReference["GH-0"].atomicMass;
  }

  std::vector<QuantumSubsystem> buildQuantumSubsystems(
    std::ostream& out, CQInputFile& input, Molecule& mol, const std::string& inputLabel) {

    const std::string base = baseLabel(inputLabel);

    // Resolve an input section in following order: 
    //    exact per-particle label        (inputLabel + suffix)
    //    the shared base label           (inputLabel's base label + suffix)
    //    legacy label                    ("PROT"/"P" + suffix)
    // For example, for inputLabel QP0 and we are looking for section "BASIS"
    //    we will first look for "QP0BASIS", then "QPBASIS", then "PBASIS"
    auto resolveSection = [&](const std::string& suffix, const std::string& legacy) -> std::string {
      const std::string exact = inputLabel + suffix;   // e.g. "QP0BASIS"
      if(input.containsSection(exact)) return exact;
      const std::string shared = base + suffix;        // e.g. "QPBASIS"
      if(base != inputLabel and input.containsSection(shared)) return shared;
      if(base == "QP" and input.containsSection(legacy)) return legacy;  // e.g. "PBASIS"
      return exact;  // default; a required-but-missing section errors downstream
    };

    const std::string qmSection = resolveSection("QM","PROTQM");

    if(not input.containsSection(qmSection))
      CErr("Missing quantum subsystem section [" + qmSection + "]");
  
    auto invalidKeywords = CQQUANTUMSUBSYSTEMQM_VALID(input.getSection(qmSection));
    printInvalidKeys(invalidKeywords,qmSection);
  
    const bool distinguishable = input.containsData(qmSection + "/DISTINGUISHABLE") ?
      input.getData<bool>(qmSection + "/DISTINGUISHABLE") : false;

    const auto quantumLabels = mol.getQuantumSystemLabels();
    // Built a comma-separated list of valid labels for the ghost center error message.
    std::string validLabelsList;
    for(size_t i = 0; i < quantumLabels.size(); ++i) validLabelsList += (i ? ", " : "") + quantumLabels[i];

    // Collect this subsystem's real particles and its basis centers (real particles + labeled ghosts) in GEOMETRY ORDER
    std::vector<size_t> qAtomIndices;      // indices of real particles 
    std::vector<size_t> basisAtomIndices;  // indices of real particles + ghost centers, geometry order
    bool hasGhostCenters = false;
    for(size_t iAtm = 0; iAtm < mol.atoms.size(); ++iAtm) {
      const auto& atom = mol.atoms[iAtm];
      if(isGhostAtom(atom)) {
        // Reject an unlabeled ghost atom
        if(atom.quantumLabel.empty())
          CErr("Ghost center (GEOM atom " + std::to_string(iAtm + 1) + ") has no subsystem label. "
               "Label every ghost with one of [" + validLabelsList + "] to indicate which subsystem it belongs to.");
        // Reject a ghost atom whose label matches no valid subsystems
        if(std::find(quantumLabels.begin(),quantumLabels.end(),atom.quantumLabel) == quantumLabels.end())
          CErr("Ghost center (GEOM atom " + std::to_string(iAtm + 1) + ") is labeled '" +
               atom.quantumLabel + "', which does not belong to valid quantum subsystem [" + validLabelsList + "].");
        // Add to this subsystem's basis only when the label matches.
        if(atom.quantumLabel == inputLabel) {
          basisAtomIndices.push_back(iAtm);
          hasGhostCenters = true;
        }
      }
      else if(atom.quantum and atom.quantumLabel == inputLabel) {
        qAtomIndices.push_back(iAtm);
        basisAtomIndices.push_back(iAtm);
      }
    }

    // Number of real particles specified in this quantum subsystem
    const size_t nPart = qAtomIndices.size();
    // Subsystems to build: DISTINGUISHABLE auto-expands the nPart particles into one subsystem each 
    //   For example: QP on 2 atoms -> QP0, QP1
    const size_t nSys = distinguishable and nPart > 1 ? nPart : 1;
    const bool autoExpand = nSys > 1;

    // Auto-expanded distinguishable particles are anonymous (QP -> QP0, QP1),
    //   so a ghost cannot know which expanded particle owns it. Error out here.
    // Require user to use explicit per-particle labels + a label-tagged ghost instead.
    if(autoExpand and hasGhostCenters)
      CErr("Ghost basis is not supported with auto-expanded distinguishable particles "
           "(DISTINGUISHABLE=TRUE on a multi-particle label). Give each particle an "
           "explicit label (e.g. " + base + "0, " + base + "1) and tag each ghost "
           "with its owner (e.g. GH x y z " + base + "0).");
    std::vector<QuantumSubsystem> systems;
    systems.reserve(nSys);
  
    for(size_t i = 0; i < nSys; ++i) {
      QuantumSubsystem sys;
  
      sys.inputLabel = inputLabel;
      sys.label = autoExpand ? inputLabel + std::to_string(i) : inputLabel;
      sys.particleIndex = autoExpand ? i : 0;
      sys.nQuantumParticles = autoExpand ? 1 : nPart;

      sys.qmSection = qmSection;
      sys.basisSection = resolveSection("BASIS","PBASIS");
      sys.dfbasisSection = resolveSection("DFBASIS","PDFBASIS");
      sys.guessBasisSection = resolveSection("GUESSBASIS","PGUESSBASIS");
      sys.intsSection = resolveSection("INTS","PINTS");
  
      if(not input.containsSection(sys.basisSection))
        CErr("Missing quantum subsystem basis section [" + sys.basisSection + "]");
  
      sys.particle = getParticleForSubsystem(input,mol,sys);

      // Basis centers for this subsystem:
      //   autoExpand: one center per particle, no ghosts (guarded above).
      //   otherwise : real particle(s) + ghost centers, in geometry order.
      if(autoExpand)
        sys.atomIndices = std::vector<size_t>{qAtomIndices[i]};
      else
        sys.atomIndices = basisAtomIndices;
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

    const std::string baseA = baseLabel(pair.inputLabelA);
    const std::string baseB = baseLabel(pair.inputLabelB);

    // Resolve the pair-integral section, trying the exact labels then the shared
    //   base labels (either order), so one [EQPINTS]/[QPQPINTS] covers every
    //   expanded pair (E-QP0, E-QP1, QP0-QP1, ...). Legacy [EPINTS] still honored.
    const std::vector<std::string> candidates = {
      pair.inputLabelA + pair.inputLabelB + "INTS",   // e.g. "EQP0INTS"
      pair.inputLabelB + pair.inputLabelA + "INTS",
      baseA + baseB + "INTS",                          // e.g. "EQPINTS"
      baseB + baseA + "INTS",
    };
    pair.intsSection = candidates.front();
    bool resolved = false;
    for(const auto& sec : candidates)
      if(input.containsSection(sec)) { pair.intsSection = sec; resolved = true; break; }
    if(not resolved and
       ((baseA == "E" and baseB == "QP") or (baseA == "QP" and baseB == "E")) and
       input.containsSection("EPINTS"))
      pair.intsSection = "EPINTS";

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
      if(isGhostAtom(atom)) continue; // ghosts carry no particle identity (charge/mass)
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
      mass = ProtMassPerE();
    else if(massAMU == atomicReference["H-2"].atomicMass)
      mass = DeutMassPerE();
    else if(massAMU == atomicReference["H-3"].atomicMass)
      mass = TritMassPerE();
    else
      mass = massAMU * AUPerAMU() - charge;
  
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
