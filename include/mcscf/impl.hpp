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

#include <mcwavefunction.hpp>
#include <mcscf.hpp>
#include <mcscf/print.hpp>
#include <mcscf/rdm.hpp>
#include <mcscf/cisolver.hpp>
#include <util/matout.hpp>
#include <orbitalrotation.hpp>

// #define DEBUG_MCSCF_IMPL

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::run(EMPerturbation & externalPert) {
    
    ProgramTimer::tick("MCSCF Total");
    // Create combined perturbation
    EMPerturbation pert;
    // Add on the MCSCF Perturbation if present
    for( auto& field : this->mcscfPert.fields )
      pert.addField( field );
    // Finally add any additional Perturbations
    for( auto& field : externalPert.fields )
      pert.addField( field );

    // allocating memeory
    this->alloc(); 
    
    // Initial Printing
    this->printMCSCFHeader(pert);
    
    // MCSCF Intial CI solution
    std::cout << "Cycle 0:\n" << std::endl;
    FormattedLine(std::cout, "AO to MO Intergral Transformation ...");
    
    ProgramTimer::tick("Solve CI");
    
    ProgramTimer::tick("Integral Trans");
    if (!this->readCI or this->settings.doSCF)
      MCWaveFunction<MatsT,IntsT>::transformInts(pert);
    ProgramTimer::tock("Integral Trans");
    
    std::cout << std::left << std::setprecision(10); 
    FormattedLine(std::cout, "Inactive Energy:", this->InactEnergy);

    ProgramTimer::tick("Diagonalization");
    if (!this->readCI)
      this->ciSolver->solveCI(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), pert);
    ProgramTimer::tock("Diagonalization");
    
    ProgramTimer::tock("Solve CI");

    // Initial 1RDM construction
    MCWaveFunction<MatsT, IntsT>::computeOneRDM();
    // MCSCF Cycles 
    
    if(this->settings.doSCF) {

      this->printStateEnergy();
      
      std::vector<double> EPrev = std::vector<double>(this->NStates);
      std::fill_n(EPrev.begin(), this->NStates, 0.);
      double EDiff      = 0.;
      bool   converged  = false;
      
      for (auto iter = 0ul; iter < settings.maxSCFIter; iter++) {
        
        ProgramTimer::tick("Orbital Rotation");
        // Exam Energy
        if(this->StateAverage) {
          EDiff = 0.;
          std::vector<double> EDiff2 = std::vector<double>(this->NStates);
          for (auto i = 0ul; i < this->NStates; i++)
            EDiff2[i] = this->StateEnergy[i] - EPrev[i];
          
          EDiff = *std::max_element(EDiff2.begin(), EDiff2.end(), 
            [&] (double a, double b) { return std::abs(a) < std::abs(b); });

        } else {
          EDiff = this->StateEnergy.back() - EPrev.back(); 
        }
        
        std::cout << "  (Maximum) Energy Difference = " << std::setw(18) 
                  <<  std::right << EDiff << std::left << std::endl; 
        if(std::abs(EDiff) <= settings.scfEnergyConv) converged = true; 
        
        // compute RDMs
        if (this->StateAverage) {
          this->computeOneRDM(); 
          this->computeTwoRDM(); 
        } else {
          this->computeOneRDM(this->NStates - 1); 
          this->computeTwoRDM(this->NStates - 1);
        }
        
        // this->print1RDMs();
        
        // compute gradient 
        double orbitalGradientNorm = 
          moRotator->computeOrbGradient(pert, *oneRDMSOI, *twoRDMSOI);
        
        // Gradient Convergence exam
        std::cout << "  Orbital Gradient Residue   = " << std::setw(18) 
                  << std::right << orbitalGradientNorm << std::left << std::endl;
        
        if(converged and orbitalGradientNorm < settings.scfGradientConv) break; 
          
        converged = false;
        
        // start of the new cycle 
        std::cout << "\n\nCycle " << iter+1 << ":\n" << std::endl;
        
        // compute hession diagonal and rotate orbitals
        FormattedLine(std::cout, "Performing Orbital Rotation ...");
        
        moRotator->rotateMO(pert, *oneRDMSOI, *twoRDMSOI);
        
        ProgramTimer::tock("Orbital Rotation");
        
        this->mointsTF->clearAllCache();
        this->moints->clear();

        // print energy and update EPrev
        std::copy_n(this->StateEnergy.begin(), this->NStates, EPrev.begin());

        ProgramTimer::tick("Solve CI");
        
        // Re-transform intgrals and solve new CI
        FormattedLine(std::cout, "Redo AO to MO Intergral Transformation ...");
        ProgramTimer::tick("Integral Trans");
        MCWaveFunction<MatsT,IntsT>::transformInts(pert);
        ProgramTimer::tock("Integral Trans");

        std::cout << std::left << std::setprecision(10); 
        FormattedLine(std::cout, "Inactive Energy:", this->InactEnergy);
        ProgramTimer::tick("Diagonalization");
        this->ciSolver->solveCI(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), pert);
        ProgramTimer::tock("Diagonalization");
        
        ProgramTimer::tock("Solve CI");
      
        this->printStateEnergy();
      
        saveCurrentStates();
      
      } // SCF Iteration

      ROOT_ONLY(this->comm);

      if(not converged) 
        CErr("\n MCSCF failed to converged in " + std::to_string(settings.maxSCFIter) + " cycles !");
    
      // compute 1RDMs
      if (this->StateAverage) this->computeOneRDM(); 
      else this->computeOneRDM(this->NStates - 1); 
      
      // generate IVOs as needed
      ProgramTimer::tick("Gen IVOs");
      if (this->settings.doIVOs) moRotator->generateIVOs(pert, *oneRDMSOI);   
      ProgramTimer::tock("Gen IVOs");

    } // doSCF

    // Final printing and save states for restart
    std::cout << "\n\nMCSCF Complete!" << std::endl;
    std::cout << bannerEnd << std::endl;

    this->printMCSCFFooter();

    // property calculation
    ProgramTimer::tick("Property Eval");

    // dipole moment
    if( this->multipoleMoment )
      MCWaveFunction<MatsT,IntsT>::computeMultipole();

    // mulliken analysis
    if (this->PopulationAnalysis) {
      std::cout<<"\n\nPopulation analysis in mcscf."<<std::endl;
      MCWaveFunction<MatsT,IntsT>::populationAnalysis();
    }

    if (this->printRDMs==0 && this->SpinAnalysis) {
      std::cout<<"\n\nSpin analysis in mcscf."<<std::endl;
      MCWaveFunction<MatsT,IntsT>::spinAnalysis();
    }    

    // oscillator strength
    if (this->NosS1) {


      this->osc_str = CQMemManager::get().malloc<double>(this->NosS1*this->NStates);
      for (size_t s1 = 0ul; s1 < this->NosS1; s1++)
      for (size_t s2 = 0ul; s2 < this->NStates; s2++){
//        if (s2 < this->NosS1) this->osc_str[s2+s1*this->NStates] = 0.;
        if (s2 <= s1) this->osc_str[s2+s1*this->NStates] = 0.;
        else this->osc_str[s2+s1*this->NStates] = 
                MCWaveFunction<MatsT,IntsT>::oscillator_strength(s2,s1);
      }

    }

    ProgramTimer::tock("Property Eval");

    saveCurrentStates(true);

    ProgramTimer::tock("MCSCF Total");
 
  }; //MCSCF::run
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::saveCurrentStates( bool saveProp ) {
    
    MCWaveFunction<MatsT, IntsT>::saveCurrentStates(saveProp);
    
    // only save MO when doing orbital rotation
    if (settings.doSCF and this->savFile.exists()) {
      auto mo_dim = this->reference().mo[0].dimension();
      this->savFile.safeWriteData("SCF/MO1", this->reference().mo[0].pointer(), {mo_dim, mo_dim});
    }

  }; // MCSCF::saveCurrentStates


  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::alloc() {
    
    MCWaveFunction<MatsT,IntsT>::alloc();
    
    ciSolver = std::make_shared<CISolver<MatsT,IntsT>>(settings.ciAlg, 
      settings.maxCIIter, settings.ciVectorConv,
      settings.maxDavidsonSpace, settings.nDavidsonGuess,
      settings.energyRefs);
    
    if (this->settings.doSCF) {
      
      this->cacheHalfTransTPI_ = true;
      
      size_t nCorrO = this->MOPartition.nCorrO;
      oneRDMSOI = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
      twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);
      
      moRotator = std::make_shared<OrbitalRotation<MatsT, IntsT>>(
        dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), settings.ORSettings);
    }
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::dealloc() {
    oneRDMSOI = nullptr;
    twoRDMSOI = nullptr;
    ciSolver  = nullptr;
    moRotator = nullptr;
  }
}; // namespace ChronusQ

