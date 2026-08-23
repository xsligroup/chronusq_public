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

#include <chronusq_sys.hpp>
#include <singleslater.hpp>
#include <matrix.hpp>

#include <dft.hpp>
#include <singleslater/kohnsham.hpp>
#include <fockbuilder/neofock.hpp>
#include <fockbuilder/interparticlefock.hpp>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/contract/batched_direct_interparticle.hpp>


namespace ChronusQ {

  // Pure virtual class for only interface functions
  struct MultiParticleSSBase {
    virtual std::shared_ptr<SingleSlaterBase> getSubSSBase(std::string label) = 0;
    virtual std::vector<std::string> getLabels() = 0;
    virtual void setSubSetup() = 0;
    virtual void saveSubsystemReferenceTypes() = 0;
  };

  template <typename MatsT, typename IntsT>
  class MultiParticleSS: virtual public SingleSlater<MatsT,IntsT>, public MultiParticleSSBase {

    template <typename MatsU, typename IntsU>
    friend class MultiParticleSS;

    template <typename F>
    void applyToEach(F func) {
      for(auto& label: order_) {
        func(subsystems.at(label));
      }
    }

    template <typename F>
    void applyToEachLabeled(F func) {
      for(auto& label: order_) {
        func(label, subsystems.at(label));
      }
    }

    template <typename F>
    void applyToEach(F func, std::vector<std::string> subset) {
      if(subset.size())
      {
        for(auto & label: subset)
        {
          func(subsystems.at(label));
        }
      }
      else
      {
        applyToEach(func);
      }
    }

    protected:

      // typedefs
      using SubSSPtr = std::shared_ptr<SingleSlater<MatsT,IntsT>>;

      template <typename T>
      using LabeledMap = std::unordered_map<std::string, T>;

      // Active XC terms, required density inputs and requested VXC outputs.
      struct XCTerms {
        struct Subsystem {
          std::string label;
          bool formVXC = false;
          bool formIntraXC = false;
        };
        struct InterTerm {
          size_t functionalIndex;
          size_t firstIndex, secondIndex;
        };
        std::vector<Subsystem> subsystems;
        std::vector<InterTerm> interTerms;
      };

      XCTerms buildXCTerms() const;
      XCTerms buildXCTermsFor(const std::vector<std::string>&) const;
      void formXCInHouse(EMPerturbation&, const XCTerms&);
      void formXCGauXC(EMPerturbation&, const XCTerms&);
      std::vector<double> formXCGradient(EMPerturbation&, const XCTerms&);
      std::vector<double> formXCGradientInHouse(EMPerturbation&, const XCTerms&);
      std::vector<double> formXCGradientGauXC(EMPerturbation&, const XCTerms&);

      //
      // MAIN STORAGE
      //
      std::unordered_map<std::string,SubSSPtr> subsystems;

      // Ordering of subsystems
      std::vector<std::string> order_;

      // Subsystem coulomb interactions
      // First map is the "external" particle and second map is the integrated
      //   particle.
      // EXAMPLE: interCoulomb["electron"]["proton"] is the coulomb matrix for
      //   the electronic subsystem (in the electronic basis) coming from the
      //   protonic coulomb potential.
      LabeledMap<LabeledMap<cqmatrix::Matrix<MatsT>>> interCoulomb;

      // Two particle integral objects (same storage scheme as above)
      // Boolean is contractSecond
      LabeledMap<LabeledMap<std::pair<bool, std::shared_ptr<TwoPInts<IntsT>>>>> interIntegrals;
      // Two particle gradient integrals
      LabeledMap<LabeledMap<std::pair<bool, std::shared_ptr<GradInts<TwoPInts,IntsT>>>>> gradInterInts;

      // Storage for full chain of FockBuilders (determines lifetime) for each subsystem
      LabeledMap<std::vector<std::shared_ptr<FockBuilder<MatsT,IntsT>>>> fockBuilders;

      // Direct access to specific inter-particle FockBuilders.
      // Note direct builders are registered here. They are handled by BatchedDirectInterparticleJContraction
      LabeledMap<LabeledMap<std::shared_ptr<InterParticleFockBuilder<MatsT,IntsT>>>> interFockBuilders;

      // Batched AO-direct contraction for inter-particle interactions
      struct BatchedDirectPair {
        std::string firstSubsystemLabel;
        std::string secondSubsystemLabel;
        std::shared_ptr<DirectTPI<IntsT>> integrals;
        // Whether the (first, second) subsystem order matches the (basisSet(), basisSet2()) order of the shared DirectTPI.
        bool subsystemOrderMatchesIntegralOrder = true;
      };
      std::vector<BatchedDirectPair> batchedDirectPairs;
      LabeledMap<cqmatrix::Matrix<MatsT>> batchedDirectCoulombMatrices;
      BatchedDirectInterparticleJContraction<MatsT,IntsT> batchedDirectInterparticleJContraction;

      void formBatchedDirectInterparticleCoulomb(const std::vector<std::string>& targets, bool increment);

      // True when this (unordered) subsystem pair is handled by the batched
      //   direct interparticle contraction instead of a recursive FockBuilder
      //   chain node. Such pairs have no entry in interFockBuilders.
      bool isBatchedPair(const std::string& first, const std::string& second) const;

      LabeledMap<SingleSlaterOptions> subsystemGuessOptions;

      // Storage for subsystem dipoles (bare particle dipoles moments, not including classical nuclear contributions)
      LabeledMap<std::array<double,3>> subsystemDipole;


      struct InterXCFunctional {
        std::string first, second;
        std::vector<std::shared_ptr<DFTFunctional>> functionals;
        double energy = 0.;
      };

      // Inter-particle correlation (e.g. EPC) functional 
      std::vector<InterXCFunctional> interFunctionals;

      // Update all XC energy contributions from all functionals
      //    Flag needed for when we want to evaluate VXC of some target subsystem and update the relevant XC energy,
      //    but not touch the energies of other unaffected subsystems
      bool needsFullXCEnergyUpdate = true;

      // Numerical integration controls for the inter-particle XC grid
      IntegrationParam intParam;

      double lastE = 0.0;       ///< Energy storage if doing a stepwise optimization

    public:

      // Main constructor
      // XXX: Does NOT construct electronic or protonic wavefunctions
      template <typename... Args>
      MultiParticleSS(MPI_Comm c, Molecule &mol, BasisSet &basis,
                std::shared_ptr<Integrals<IntsT>> aoi, Args... args) :
        SingleSlater<MatsT,IntsT>(c,mol,basis,aoi,args...),
        WaveFunctionBase(c,mol,basis,args...),
        QuantumBase(c,args...) {};

      // TODOAL: Implement these semantics
      MultiParticleSS(const MultiParticleSS&) = delete;
      MultiParticleSS(MultiParticleSS&&) = delete;
      MultiParticleSS& operator=(const MultiParticleSS&) = delete;
      MultiParticleSS& operator=(MultiParticleSS&&) = delete;

      // TODOAL: Add type conversion constructors for other SingleSlater types


      // Add a subsystem to the NEO object
      //
      // label  Unique name for this subsystem
      // ss     SingleSlater object representing this subsystem
      void addSubsystem(
        const std::string label,
        std::shared_ptr<SingleSlater<MatsT,IntsT>> ss);

      void setSubsystemGuessOptions(const std::string& label,
        const SingleSlaterOptions& options) {
        subsystemGuessOptions[label] = options;
      }

      // Add pair-wise interaction between two subsystems
      void addInteraction(const std::string label1, const std::string label2,
        const std::shared_ptr<TwoPInts<IntsT>>& ints,
        bool contractSecond);

      void addGradientIntegrals(std::string label1, std::string label2,
        std::shared_ptr<GradInts<TwoPInts,IntsT>> ints, bool contractSecond);

      void setSubSSTPIContraction(std::shared_ptr<TwoPInts<double>> tpi);

      // Register inter-particle correlation functionals for a subsystem pair.
      void setInterFunctionals(const std::string& label1, const std::string& label2,
        const std::vector<std::shared_ptr<DFTFunctional>>& funcs) {
        if(!subsystems.count(label1)) CErr("Unknown subsystem label: " + label1);
        if(!subsystems.count(label2)) CErr("Unknown subsystem label: " + label2);
        for(auto& pair : interFunctionals) {
          bool matches = (pair.first == label1 and pair.second == label2) or
                         (pair.first == label2 and pair.second == label1);
          if( !matches ) continue;
          pair.functionals = funcs;
          pair.energy = 0.;
          needsFullXCEnergyUpdate = true;
          return;
        }
        interFunctionals.push_back({label1, label2, funcs, 0.});
        needsFullXCEnergyUpdate = true;
      }

      // Whether any inter-particle correlation functionals are present
      bool hasInterXC() const {
        return std::any_of(interFunctionals.begin(), interFunctionals.end(),
          [](const auto& pair) { return not pair.functionals.empty(); });
      }

      void setIntParam(const IntegrationParam& ip) {
        intParam = ip;
        needsFullXCEnergyUpdate = true;
      }

      void setOrder(std::vector<std::string> labels) {
        order_ = labels;
        needsFullXCEnergyUpdate = true;
      }

      // Getters
			std::vector<std::string> getLabels() override {
				return order_;
			}
			
			std::pair<bool,std::shared_ptr<TwoPInts<IntsT>>> getCrossTPIs(std::string label1,std::string label2){
		  	return interIntegrals.at(label1).at(label2);
			}

      template <template <typename, typename> class T>
      std::shared_ptr<T<MatsT,IntsT>> getSubsystem(std::string label) {
        return std::dynamic_pointer_cast<T<MatsT,IntsT>>(subsystems.at(label));
      }
      std::shared_ptr<SingleSlaterBase> getSubSSBase(std::string label) override {
        return std::dynamic_pointer_cast<SingleSlaterBase>(subsystems.at(label));
      }
      std::shared_ptr<SingleSlater<MatsT,IntsT>> getSubSS(std::string label){
        return std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(subsystems.at(label));
      }
      template <template <typename, typename> class T>
      std::vector<std::shared_ptr<T<MatsT,IntsT>>> getAllSubsystems() {
        std::vector<std::shared_ptr<T<MatsT,IntsT>>> results;
        applyToEach([&](SubSSPtr& ss){
          results.push_back(std::dynamic_pointer_cast<T<MatsT,IntsT>>(ss));
        });
        return results;
      }
      std::vector<std::shared_ptr<SingleSlaterBase>> getAllSubSSBase() {
        std::vector<std::shared_ptr<SingleSlaterBase>> results;
        applyToEach([&](SubSSPtr& ss){
          results.push_back(std::dynamic_pointer_cast<SingleSlaterBase>(ss));
        });
        return results;
      }

      const std::unordered_map<std::string, SubSSPtr>& getSubsystemMap() {
        return subsystems;
      }

      const std::vector<std::string>& getOrder() {
        return order_;
      }

      const std::unordered_map<std::string, std::array<double,3>>& getSubsystemDipoles() const {
        return subsystemDipole;
      }

      std::shared_ptr<InterParticleFockBuilder<MatsT,IntsT>> getInterFockBuilder(std::string label1, std::string label2) const {
        return interFockBuilders.at(label1).at(label2);
      }

      void saveSubsystemReferenceTypes(std::string prefix) {
        ROOT_ONLY(this->comm);
        if( !this->savFile.exists() ) return;

        const std::string root = prefix + "MULTISS/";
        applyToEachLabeled([&](const std::string& label, SubSSPtr&){
          int refType = subsystemGuessOptions.at(label).refOptions.refType;
          this->savFile.safeWriteData(root + label + "/SCF/REFTYPE", &refType, {1});
        });
      }

      void saveSubsystemReferenceTypes() override {
        saveSubsystemReferenceTypes("");
      }

      void saveCurrentState(bool saveMO = true, std::string prefix = "") override {
        const std::string root = prefix + "MULTISS/";

        // Pass-through to each subsystem
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          ss->saveCurrentState(saveMO, root + label + "/");
        });
        saveSubsystemReferenceTypes(prefix);
        ROOT_ONLY(this->comm);

        if( !this->savFile.exists() )
          CErr("savFile does not exist!");

        // Per-subsystem energies, keyed by label
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          this->savFile.safeWriteData(root + label + "/ENERGY",
                                      &ss->totalEnergy, {1});
        });

        this->savFile.safeWriteData(root + "NUC_REP_ENERGY",
                                    &this->molecule().nucRepEnergy, {1});
        this->savFile.safeWriteData(root + "TOTAL_ENERGY",
                                    &this->totalEnergy, {1});

        // Aggregate molecular multipoles
        this->savFile.safeWriteData(root + "LEN_ELECTRIC_DIPOLE",
                                    &this->elecDipole[0], {3});
        this->savFile.safeWriteData(root + "LEN_ELECTRIC_QUADRUPOLE",
                                    &this->elecQuadrupole[0][0], {3, 3});
        this->savFile.safeWriteData(root + "LEN_ELECTRIC_OCTUPOLE",
                                    &this->elecOctupole[0][0][0], {3, 3, 3});

        // Per-subsystem bare dipoles
        if( not subsystemDipole.empty() ) {
          applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
            this->savFile.safeWriteData(root + label + "/DIPOLE",
                                        &subsystemDipole.at(label)[0], {3});
          });
        }
      }

      void initializeSCF() override;

      void formGuess(EMPerturbation &pert, const SingleSlaterOptions& ssopt) override {
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          const auto it = subsystemGuessOptions.find(label);
          ss->formGuess(pert, it == subsystemGuessOptions.end() ? ssopt : it->second);
        });
      }

      void formFockForTargets(EMPerturbation& emPert, const std::vector<std::string>& targets, 
          bool increment = false, double xHFX = 1.) {
        // Each subsystem builds its own Fock (J + inter-particle Coulomb via its
        //   FockBuilder chain, plus its own intra-particle VXC if it is a KohnSham
        //   object with doVXC_ enabled and it's an active XC target)
        applyToEach([&](SubSSPtr& ss){ ss->formFock(emPert, increment, xHFX); }, targets);

        formBatchedDirectInterparticleCoulomb(targets, increment);

        if( !hasInterXC() ) return;
        // Form XC on the shared grid in a single pass for target subsystems (default=all subsystems)
        auto xcTerms = targets.empty() ? buildXCTerms() : buildXCTermsFor(targets);
        formXC(emPert, xcTerms);
      }

      virtual void formFock(EMPerturbation& emPert, bool increment = false, double xHFX = 1.) override {
        formFockForTargets(emPert, this->scfControls.NEOSubSystemOpt, increment, xHFX);
      }

      // Unified intra- and inter-particle XC driver.
      void formXC(EMPerturbation&, const XCTerms&);

      void formCoreH(EMPerturbation& emPert, bool save) override {
        applyToEach([&](SubSSPtr& ss){ ss->formCoreH(emPert, save); });
      }

      virtual void formDensity() override {
        applyToEach([&](SubSSPtr& ss){ ss->formDensity(); },
          this->scfControls.NEOSubSystemOpt);
      }

      virtual void printOrbitalPopulation(std::ostream& out) {
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          out << bannerTop << std::endl;
          out << "Subsystem " << label << " MO Occupation:" << std::endl;
          ss->printOrbitalPopulation(out);
        });
        out << bannerTop << std::endl;
      }

      void setupRangeSeparatedHybridExchange() override {
        if (not this->gauxcUtils or not this->gauxcUtils->isRangeSeparatedHybrid())
          return;

        // Propagate shared settings (including gauxcUtils) to subsystems first
        setSubSetup();

        // Rebuild the electronic subsystem's short-range erfc exchange on the current geometry
        auto electronic = subsystems.find("E");
        if(electronic == subsystems.end())
          CErr("NEO range-separated hybrid exchange requires an electronic Kohn-Sham subsystem labeled E.");
        auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(electronic->second);
        if(not ks or ks->nC != 1)
          CErr("NEO range-separated hybrid exchange is currently implemented only for RKS/UKS Kohn-Sham electronic subsystems.");

        ks->setupRangeSeparatedHybridExchange();
      }

      // Propagate options that were set by value in the *Options functions
      void setSubSetup() override {
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          ss->scfControls = this->scfControls;
          if(const auto it = subsystemGuessOptions.find(label);
             it != subsystemGuessOptions.end()) {
            ss->scfControls.guess = it->second.scfControls.guess;
            ss->scfControls.guessBasis = it->second.scfControls.guessBasis;
          }
          ss->savFile = this->savFile;
          ss->fchkFileName = this->fchkFileName;
          ss->gauxcUtils = this->gauxcUtils;
          ss->scrBinFileName = this->scrBinFileName;
        });
      }

      void printSetup(std::ostream& out);

      std::vector<double> getGrad(EMPerturbation&, bool, bool, double xHFX = 1.) override;

      virtual void checkIdempotency(std::string system="") {
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss){
          ss->checkIdempotency(label);
        });
      }

      virtual void formEWDM(bool equil = false) {
        applyToEach([&](SubSSPtr& ss){ ss->formEWDM(equil); });
      }
      // Functions for OrbitalModifier
      virtual void runSCF(EMPerturbation&) override;
      virtual void buildOrbitalModifierOptions() override;
      void printProperties() override;
      virtual std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getOnePDM() override;
      virtual std::vector<cqmatrix::Matrix<MatsT>> getOnePDMOrtho() override;
      virtual std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getFock() override;
      virtual void setOnePDMOrtho(cqmatrix::Matrix<MatsT>*) override;
      virtual void setOnePDMAO(cqmatrix::Matrix<MatsT>*) override;
      virtual std::vector<std::shared_ptr<Orthogonalization<MatsT>>> getOrtho() override;
      virtual double getTotalEnergy() override { return this->totalEnergy; };
      virtual void setDenEqCoeff(bool val);
      virtual void ortho2aoMOs();
      virtual void ao2orthoDen();
      virtual void ortho2aoDen() override;

      virtual bool secondSCF() override;

      // Properties
      using QuantumBase::computeEnergy;
      void computeEnergy() override {

        this->totalEnergy = 0.;
        applyToEach([&](SubSSPtr& ss){
          ss->computeEnergy();
          // Each subsystem's energy already includes its own intra-particle XC
          //   energy (via KohnSham::computeEnergy) and its share of the
          //   inter-particle Coulomb energy. It also includes the classical
          //   nuclear repulsion energy, which must be removed here and added
          //   back once below.
          this->totalEnergy += ss->totalEnergy - this->molecule().nucRepEnergy;
        });

        // Add each inter-particle correlation energy exactly once.
        for(const auto& pair : interFunctionals)
          this->totalEnergy += pair.energy;

        this->totalEnergy += this->molecule().nucRepEnergy;

      }
    
      void computeMultipole(EMPerturbation& emPert,
        const std::vector<PROPERTY>& properties = {}) {
    
        // Zero the aggregate molecular multipoles
        for (auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          this->elecDipole[iXYZ] = 0.;
          for (auto jXYZ = 0; jXYZ < 3; jXYZ++) {
            this->elecQuadrupole[iXYZ][jXYZ] = 0.;
            for (auto kXYZ = 0; kXYZ < 3; kXYZ++)
              this->elecOctupole[iXYZ][jXYZ][kXYZ] = 0.;
          }
        }
    
        subsystemDipole.clear();
    
        // Accumulate each subsystem's multipole moments together
        applyToEachLabeled([&](const std::string& label, SubSSPtr& ss) {
          ss->computeMultipole(emPert, properties);
    
          for (auto iXYZ = 0; iXYZ < 3; iXYZ++) {
            this->elecDipole[iXYZ] += ss->elecDipole[iXYZ];
            for (auto jXYZ = 0; jXYZ < 3; jXYZ++) {
              this->elecQuadrupole[iXYZ][jXYZ] += ss->elecQuadrupole[iXYZ][jXYZ];
              for (auto kXYZ = 0; kXYZ < 3; kXYZ++)
                this->elecOctupole[iXYZ][jXYZ][kXYZ] += ss->elecOctupole[iXYZ][jXYZ][kXYZ];
            }
          }
    
          // Compute and store the bare density moment for each subsystem by removing the classical nuclear term
          std::array<double,3> bare = {ss->elecDipole[0], ss->elecDipole[1], ss->elecDipole[2]};
          for (auto& atom : this->molecule().atoms) {
            if (atom.quantum) continue;
            for (int iXYZ = 0; iXYZ < 3; iXYZ++)
              bare[iXYZ] -= atom.nucCharge * atom.coord[iXYZ];
          }
          subsystemDipole[label] = bare;
        });
    
        // Remove the over-counted classical nuclear contribution from the aggregate molecular multipoles
        const double overCount = static_cast<double>(subsystems.size() - 1);
        for (auto& atom : this->molecule().atoms) {
          if (atom.quantum) continue;
    
          for (int iXYZ = 0; iXYZ < 3; iXYZ++)
            this->elecDipole[iXYZ] -= overCount * atom.nucCharge * atom.coord[iXYZ];
    
          for (size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            for (size_t jXYZ = 0; jXYZ < 3; jXYZ++)
              this->elecQuadrupole[iXYZ][jXYZ] -=
                overCount * atom.nucCharge * atom.coord[iXYZ] * atom.coord[jXYZ];
    
          for (size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            for (size_t jXYZ = 0; jXYZ < 3; jXYZ++)
              for (size_t kXYZ = 0; kXYZ < 3; kXYZ++)
                this->elecOctupole[iXYZ][jXYZ][kXYZ] -=
                  overCount * atom.nucCharge * atom.coord[iXYZ] * atom.coord[jXYZ] * atom.coord[kXYZ];
        }
      };

      void computeSpin() override {
        applyToEach([](SubSSPtr& ss){ ss->computeSpin(); });      
      }

      void methodSpecificProperties() override {
        applyToEach([](SubSSPtr& ss){ ss->methodSpecificProperties(); });      
      }

      // Disable NR/stability for now
      MatsT* getNRCoeffs() override {
        CErr("NR NYI for NEO!");
        return nullptr;
      }

      std::pair<double,MatsT*> getStab() override {
        CErr("NR NYI for NEO!");
        return {0., nullptr};
      }
  };

}

#include <singleslater/multiparticless/impl.hpp>
#include <singleslater/multiparticless/scf.hpp>
#include <singleslater/multiparticless/xcterms.hpp>
#include <singleslater/multiparticless/vxc.hpp>
#include <singleslater/multiparticless/gradient.hpp>
