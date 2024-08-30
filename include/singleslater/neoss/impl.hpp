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

#include <singleslater/neoss.hpp>
#include <singleslater/neoss/scf.hpp>
#include <singleslater/neoss/gradient.hpp> 
#include <cerr.hpp>

namespace ChronusQ {

#define TRY_REF(_REF, input, output) \
  if( output == nullptr ) \
    if( auto casted = std::dynamic_pointer_cast<_REF<MatsU,IntsT>>(input) ) \
      output = std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>( \
        std::make_shared< _REF<MatsT,IntsT> >(*casted) \
      );

  // Helper function to handle determining type of subsystem
  template <typename MatsT, typename IntsT, typename MatsU>
  std::shared_ptr<SingleSlater<MatsT,IntsT>> makeNewSS(const std::shared_ptr<SingleSlater<MatsU,IntsT>>& old_ss) {
    std::shared_ptr<SingleSlater<MatsT,IntsT>> new_ss = nullptr;

    TRY_REF(KohnSham, old_ss, new_ss);
    TRY_REF(HartreeFock, old_ss, new_ss);

    if( new_ss == nullptr )
      CErr("Unrecognized reference in constructing NEOSS!");

    return new_ss;
  };

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  NEOSS<MatsT,IntsT>::NEOSS(const NEOSS<MatsU,IntsT>& other, int dummy) :
    SingleSlater<MatsT,IntsT>(dynamic_cast<const SingleSlater<MatsU,IntsT>&>(other), dummy),
    WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other)),
    QuantumBase(dynamic_cast<const QuantumBase&>(other)),
    order_(other.order_), functionals(other.functionals)
  {
    // Loop over all old subsystems
    for( auto& x: other.subsystems ) {

      // Easier to read names
      auto& label = x.first;
      auto& old_sys = x.second;

      // Find a non-NEO fockbuilder
      auto baseFock = old_sys->fockBuilder.get();
      if( auto p = dynamic_cast<NEOFockBuilder<MatsU,IntsT>*>(baseFock) ) {
        baseFock = p->getNonNEOUpstream();
      }

      // Assume just non-relativistic Fock here
      // TODO: Make this general to all fockbuilders
      auto newFock = std::make_shared<FockBuilder<MatsT,IntsT>>(*baseFock);

      // Get new SingleSlater object
      SubSSPtr new_sys = makeNewSS<MatsT>(old_sys);
      new_sys->fockBuilder = newFock;

      // If we have EPC functionals and this is a nuclear subsystem, add the
      // functionals (only works for electron/proton systems)
      if( old_sys->particle.charge > 0 ) {
        if( auto ks = std::dynamic_pointer_cast<KohnSham<MatsU,IntsT>>(old_sys) ) {
          auto new_ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(new_sys);
          new_ks->functionals.insert(new_ks->functionals.end(), functionals.begin(), functionals.end());
        }
      }

      // Add the new subsystem to the current object
      addSubsystem(label, new_sys, other.interIntegrals.at(label));
    }

    // Copy gradient integrals
    for( auto& gradSys1: other.gradInterInts ) {
      for( auto& gradSys2: gradSys1.second ) {

        auto label1 = gradSys1.first;
        auto label2 = gradSys2.first;
        auto intPtr = gradSys2.second.second;
        auto contractSecond = gradSys2.second.first;

        addGradientIntegrals(label1, label2, intPtr, contractSecond);
      }
    }

  }; // Other type copy constructor

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  NEOSS<MatsT,IntsT>::NEOSS(NEOSS<MatsU,IntsT>&& other, int dummy) :
    SingleSlater<MatsT,IntsT>(dynamic_cast<SingleSlater<MatsU,IntsT>&&>(other), dummy),
    WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(other)),
    QuantumBase(dynamic_cast<QuantumBase&&>(other)),
    order_(std::move(other.order_)), functionals(std::move(other.functionals))
  {
    CErr("NEOSS move constructor not fully implemented");
  }; // Other type move constructor


  // Add a subsystem to the NEO object
  //
  // label  Unique name for this subsystem
  // ss     SingleSlater object representing this subsystem
  // ints   Two particle integrals for interacting with all other subsystems
  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::addSubsystem(
    std::string label,
    std::shared_ptr<SingleSlater<MatsT,IntsT>> ss,
    const LabeledMap<std::pair<bool,std::shared_ptr<TwoPInts<IntsT>>>>& ints)
  {

    // Check that we have at least the same number of integrals as other systems
    if(ints.size() < subsystems.size())
      CErr("Subsystem and integral number mismatch in addSubsystem!");

    auto NB = ss->basisSet().nBasis;

    bool this_nuclear = ss->particle.charge > 0;
    auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss);

    // Add a FockBuilder and Coulomb matrix:
    //   - For the new system's interaction with each other system
    //   - For each other system's interaction with the new system

    // Add the base fockBuilder to make sure it has a persistent lifetime
    fockBuilders.insert({label, {ss->fockBuilder}});
    // New coulomb matrices to be added to the _new_ system
    std::unordered_map<std::string, cqmatrix::Matrix<MatsT>> newCoulombs;

    // Loop over other subsystems
    for( auto& x: subsystems ) {

      // first = label;
      // second = SingleSlater;
      auto other_NB = x.second->basisSet().nBasis;

      // Add a new coulomb matrix to the new system
      newCoulombs.insert({x.first, cqmatrix::Matrix<MatsT>(NB)});
      // Add a new coulomb matrix to the other system
      interCoulomb.at(x.first).insert({label, cqmatrix::Matrix<MatsT>(other_NB)});

      HamiltonianOptions this_options = ss->aoints_->options_;
      HamiltonianOptions other_options = x.second->aoints_->options_;

      // New fock builders
      auto this_newFock = std::make_shared<NEOFockBuilder<MatsT,IntsT>>(this_options);
      auto other_newFock = std::make_shared<NEOFockBuilder<MatsT,IntsT>>(other_options);

      this_newFock->setAux(x.second.get());
      this_newFock->setOutput(&newCoulombs.at(x.first));
      this_newFock->setUpstream(fockBuilders[label].back().get());

      other_newFock->setAux(ss.get());
      other_newFock->setOutput(&interCoulomb.at(x.first).at(label));
      other_newFock->setUpstream(fockBuilders[x.first].back().get());

      // Add integrals to other system (mostly for consistency)
      auto& contractSecond = ints.at(x.first).first;
      auto& tpi = ints.at(x.first).second;
      interIntegrals[x.first].insert({label, {not contractSecond, tpi}});

      // Contractions
      std::shared_ptr<TPIContractions<MatsT,IntsT>> this_cont;
      std::shared_ptr<TPIContractions<MatsT,IntsT>> other_cont;
      if( auto tpi_t = std::dynamic_pointer_cast<DirectTPI<IntsT>>(tpi) ) {
        this_cont = std::make_shared<GTODirectTPIContraction<MatsT,IntsT>>(tpi_t);
        other_cont = std::make_shared<GTODirectTPIContraction<MatsT,IntsT>>(tpi_t);
      }
      else if( auto tpi_t = std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(tpi) ) {
        this_cont = std::make_shared<InCore4indexTPIContraction<MatsT,IntsT>>(tpi_t);
        other_cont = std::make_shared<InCore4indexTPIContraction<MatsT,IntsT>>(tpi_t);
      }
      else if( auto tpi_t = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(tpi)) {
        this_cont = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(tpi_t);
        other_cont = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(tpi_t);
      }
      else {
        CErr("Invalid TwoPInts type for NEO!");
      }
      // Set contractSecond
      this_cont->contractSecond = contractSecond;
      other_cont->contractSecond = not contractSecond;

      // Set printContractionTiming
      this_cont->printContractionTiming = ss->TPI->printContractionTiming;
      other_cont->printContractionTiming = ss->TPI->printContractionTiming;

      this_newFock->setContraction(this_cont);
      other_newFock->setContraction(other_cont);

      fockBuilders[label].push_back(this_newFock);
      fockBuilders[x.first].push_back(other_newFock);

      // KS specific fock builders
      // XXX: This assumes that there is no correlation functional between
      //      nuclear subsystems
      auto other_ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(x.second);
      bool other_nuclear = x.second->particle.charge > 0;
      std::shared_ptr<NEOKohnShamBuilder<MatsT,IntsT>> this_newKSBuilder = nullptr;
      std::shared_ptr<NEOKohnShamBuilder<MatsT,IntsT>> other_newKSBuilder = nullptr;
      if( other_ks && other_nuclear && !this_nuclear ) {
        this_newKSBuilder = std::make_shared<NEOKohnShamBuilder<MatsT,IntsT>>(this_options);
        other_newKSBuilder = std::make_shared<NEOKohnShamBuilder<MatsT,IntsT>>(other_options);

        this_newKSBuilder->setIntParam(other_ks->intParam);
        other_newKSBuilder->setIntParam(other_ks->intParam);

        this_newKSBuilder->setFunctionals(other_ks->functionals);
        other_newKSBuilder->setFunctionals(other_ks->functionals);

        // Get EPC out of the functional list of the nuclear subsystem
        other_ks->functionals.clear();
      }
      else if( ks && this_nuclear && !other_nuclear ) {
        this_newKSBuilder = std::make_shared<NEOKohnShamBuilder<MatsT,IntsT>>(this_options);
        other_newKSBuilder = std::make_shared<NEOKohnShamBuilder<MatsT,IntsT>>(other_options);

        this_newKSBuilder->setIntParam(ks->intParam);
        other_newKSBuilder->setIntParam(ks->intParam);

        this_newKSBuilder->setFunctionals(ks->functionals);
        other_newKSBuilder->setFunctionals(ks->functionals);

        // Get EPC out of the functional list of the nuclear subsystem
        this->functionals.insert(functionals.end(),
          ks->functionals.begin(), ks->functionals.end());
        ks->functionals.clear();
      }

      if( this_newKSBuilder ) {
        this_newKSBuilder->setAux(x.second.get());
        this_newKSBuilder->setUpstream(fockBuilders[label].back().get());

        other_newKSBuilder->setAux(ss.get());
        other_newKSBuilder->setUpstream(fockBuilders[x.first].back().get());

        fockBuilders[label].push_back(this_newKSBuilder);
        fockBuilders[x.first].push_back(other_newKSBuilder);
      }

    }

    // Add the other system's coulomb matrices to the coulomb matrix storage
    interCoulomb.insert({label, std::move(newCoulombs)});

    // Add the single slater object to the list of systems
    subsystems[label] = ss;

    // Add the integrals to this system
    interIntegrals.insert({label, ints});

    // Update all FockBuilders used to the most recent version
    for( auto& x: subsystems ) {
      x.second->fockBuilder = fockBuilders.at(x.first).back();
    }

  };


  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::addGradientIntegrals(
    std::string label1, std::string label2,
    std::shared_ptr<GradInts<TwoPInts,IntsT>> ints, bool contractSecond) {

    gradInterInts.insert({label1, {}});
    gradInterInts.insert({label2, {}});

    gradInterInts[label1].insert({label2, {contractSecond, ints}});
    gradInterInts[label2].insert({label1, {not contractSecond, ints}});

    // TODO: Make this work with nested systems. Right now there is no
    //   guarantee that it will be placed on the right fockBuilder for more
    //   than two systems. May require labeling of nested fockBuilders.
    auto setGradInts = [&](std::shared_ptr<FockBuilder<MatsT,IntsT>>& fock) {
      if(auto neofock = std::dynamic_pointer_cast<NEOFockBuilder<MatsT,IntsT>>(fock)) {
        neofock->setGradientIntegrals(ints.get());
      } else if(auto neoks = std::dynamic_pointer_cast<NEOKohnShamBuilder<MatsT,IntsT>>(fock)){
        if(auto neofock = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(neoks->getUpstream()) )
          neofock->setGradientIntegrals(ints.get());
        else
          CErr("Upstream FockBuilder incorrectly set in setGradInts!");
      } else {
        CErr("Can't set gradient integrals on a non-NEOFockBuilder");
      }
    };

    setGradInts(fockBuilders[label1].back());
    setGradInts(fockBuilders[label2].back());

  };


  template <typename MatsT, typename IntsT>
  std::vector<double> NEOSS<MatsT,IntsT>::getGrad(EMPerturbation& pert,
    bool equil, bool saveInts, double xHFX) {

    // Constants and return value
    size_t nAtoms = this->molecule().nAtoms;
    size_t nGrad = 3*nAtoms;
    std::vector<double> gradient(nGrad);

    // Initialize with classical nuclear repulsion
    for( auto iGrad = 0; iGrad < nGrad; iGrad++ )
      gradient[iGrad] = this->molecule().nucRepForce[iGrad/3][iGrad%3];

    // Determine the names of the systems to loop over (so we can loop
    //   over system1 < system2)
    std::vector<std::string> systemList;
    if( order_.size() == subsystems.size() ) {
      systemList.insert(systemList.end(), order_.begin(), order_.end());
    }
    else {
      for( auto& system: subsystems )
        systemList.push_back(system.first);
    }

    // Compute required integrals (intrasystem integrals are handled with
    //   the subsystem's getGrad call)
    for( auto iSys = 0; iSys < systemList.size(); iSys++ ) {
      for( auto jSys = 0; jSys < iSys; jSys++ ) {

        auto& sys1Label = systemList[iSys];
        auto& sys2Label = systemList[jSys];

        auto& sys1Basis = subsystems[sys1Label]->basisSet();
        auto& sys2Basis = subsystems[sys2Label]->basisSet();

        bool sys1Left = not gradInterInts[sys1Label][sys2Label].first;
        auto& gradInt12 = gradInterInts[sys1Label][sys2Label].second;

        HamiltonianOptions options;
        options.OneEScalarRelativity = false;

        if( sys1Left )
          gradInt12->computeAOInts(sys1Basis, sys2Basis, this->molecule(),
            pert, EP_ATTRACTION, options);
        else
          gradInt12->computeAOInts(sys2Basis, sys1Basis, this->molecule(),
            pert, EP_ATTRACTION, options);
      }
    }

    // Obtain the non-xc part of each subsystem gradient
    applyToEach([&](SubSSPtr& ss) {
      if (auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss)) 
        xHFX = ks->functionals.size() != 0 ? ks->functionals.back()->xHFX : 1. ;
      auto localGrad = ss->SingleSlater<MatsT, IntsT>::getGrad(pert, equil, saveInts, xHFX);
      //auto localGrad = ss->getGrad(pert, equil, saveInts, xHFX);
      for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) 
        gradient[iGrad] += localGrad[iGrad] - this->molecule().nucRepForce[iGrad/3][iGrad%3];
    });

    if(auto pss = getSubsystem<KohnSham>("Protonic")){
      for(size_t ic = 0; ic < nAtoms; ic++){
        for(size_t XYZ = 0; XYZ < 3; XYZ++){
          this->EXCGradient[ic][XYZ] = 0.0;
          this->EPCGradientE[ic][XYZ] = 0.0;
          this->EPCGradientP[ic][XYZ] = 0.0;
        }
      } 
      // Obtain the xc part of the gradient (ee_xc + epc)
      formEXCGradient();
      for( auto iGrad = 0; iGrad < nGrad; iGrad++ ){
          gradient[iGrad] += this->EXCGradient[iGrad/3][iGrad%3];
          gradient[iGrad] += this->EPCGradientE[iGrad/3][iGrad%3];
          gradient[iGrad] += this->EPCGradientP[iGrad/3][iGrad%3];
      } 
    }

    // Add additional printout for debugging NEO-Ehrenfest
    auto printGrad = [&](std::string name, std::vector<double>& vecgrad) {
      std::cout << name << std::endl;
      std::cout << std::setprecision(12);
      for( auto iAt = 0; iAt < nAtoms; iAt++ ) {
        std::cout << " Gradient@I = " << iAt << ":";
        for( auto iCart = 0; iCart < 3; iCart++ ) {
          std::cout << "  " << vecgrad[iAt*3 + iCart];
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    };
    //printGrad("Total NEO Gradient:", gradient);

    return gradient;

  };

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::buildOrbitalModifierOptions() {
    // Modify SCFControls
    this->scfControls.printLevel    = this->printLevel;
    this->scfControls.refLongName_  = this->refLongName_;
    this->scfControls.refShortName_ = this->refShortName_;

    // Initialize ModifyOrbitalOptions
    OrbitalModifierDrivers<MatsT> modOrbOpt;

    // Register functions
    modOrbOpt.printProperties   = [this]() { this->printProperties(); };
    modOrbOpt.saveCurrentState  = [this]() { this->saveCurrentState(); };
    modOrbOpt.formFock          = [this](EMPerturbation& pert) { this->formFock(pert,false,1.); };
    modOrbOpt.computeProperties = [this](EMPerturbation& pert) { this->computeProperties(pert); };
    modOrbOpt.computeEnergy     = [this](EMPerturbation& pert) { this->computeEnergy(pert); };
    modOrbOpt.formDensity       = [this]() { this->formDensity(); };
    modOrbOpt.getFock           = [this]() { return this->getFock(); };
    modOrbOpt.getOnePDM         = [this]() { return this->getOnePDM(); };
    modOrbOpt.getOrtho          = [this]() { return this->getOrtho(); };
    modOrbOpt.setDenEqCoeff     = [this](bool val) { return this->setDenEqCoeff(val); };
    modOrbOpt.getTotalEnergy    = [this]() { return this->getTotalEnergy(); };

    // Make OrbitalModifier based on scfControls
    if( this->scfControls.scfAlg == _CONVENTIONAL_SCF ) {
      // Conventional SCF
      this->orbitalModifier = std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(
          std::make_shared<ConventionalSCF<MatsT>>(this->scfControls, this->comm, modOrbOpt));
    } else if( this->scfControls.scfAlg == _NEWTON_RAPHSON_SCF ) {
      // Newton-Raphson SCF
      CErr("Newton-Raphson SCF NYI for NEO methods!");
    } else {
      // SKIP SCF
      // this->scfControls.doExtrap = false;
      // this->orbitalModifier       = std::dynamic_pointer_cast<OrbitalModifier<MatsT>>(
      //    std::make_shared<SkipSCF<MatsT>>(this->scfControls, this->comm, modOrbOpt));
    }

  }

  template <typename MatsT, typename IntsT>
  void NEOSS<MatsT,IntsT>::setSubSSTPIContraction(std::shared_ptr<TwoPInts<double>> tpi){

      if( auto tpi_t = std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(tpi)){

        std::shared_ptr<InCoreAsymmRITPIContraction<MatsT,IntsT>> elec_cont = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(tpi_t);
        std::shared_ptr<InCoreAsymmRITPIContraction<MatsT,IntsT>> prot_cont = std::make_shared<InCoreAsymmRITPIContraction<MatsT,IntsT>>(tpi_t);

        elec_cont->contractSecond = false;
        prot_cont->contractSecond = true;

        if(auto elec_fock = std::dynamic_pointer_cast<NEOFockBuilder<MatsT,IntsT>>(this->fockBuilders["Electronic"].back()) ){
          elec_fock->setContraction(elec_cont);
          if(auto prot_fock = std::dynamic_pointer_cast<NEOFockBuilder<MatsT,IntsT>>(this->fockBuilders["Protonic"].back()) ){
            prot_fock->setContraction(prot_cont);
          }else{
            CErr("Back of prot fockbuilders is NOT a NEOFockBuildber object");
          }
        }else{
          CErr("Back of elec fockbuilders is NOT a NEOFockBuildber object");
        }
        
        

      }else{
        CErr("Expecting an asymmetric ERI in NEOSS<MatsT,IntsT>::setsetSubSSTPI");
      }


  }

} // namespace ChronusQ
