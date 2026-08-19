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
#include <mcscf/base/impl.hpp>
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
    if (!mcwfn_->readCI or this->settings->doSCF)
      this->transformInts(pert);
    ProgramTimer::tock("Integral Trans");
    
    std::cout << std::left << std::setprecision(10); 
    FormattedLine(std::cout, "Inactive Energy:", mcwfn_->InactEnergy);

    ProgramTimer::tick("Diagonalization");
    if (!mcwfn_->readCI)
      this->ciSolver->solveCI(*mcwfn_, pert);
    ProgramTimer::tock("Diagonalization");
    
    ProgramTimer::tock("Solve CI");

    // Initial 1RDM construction
    this->computeOneRDM();
    // MCSCF Cycles 
    
    if(this->settings->doSCF) {

      this->printStateEnergy();
      
      std::vector<double> EPrev = std::vector<double>(this->NStates);
      std::fill_n(EPrev.begin(), this->NStates, 0.);
      double EDiff      = 0.;
      bool   converged  = false;
      
      for (auto iter = 0ul; iter < settings->maxSCFIter; iter++) {
        
        ProgramTimer::tick("Orbital Rotation");
        // Exam Energy
        if(this->StateAverage) {
          EDiff = 0.;
          std::vector<double> EDiff2 = std::vector<double>(this->NStates);
          for (auto i = 0ul; i < this->NStates; i++)
            EDiff2[i] = mcwfn_->StateEnergy->at(i) - EPrev[i];
          
          EDiff = *std::max_element(EDiff2.begin(), EDiff2.end(), 
            [&] (double a, double b) { return std::abs(a) < std::abs(b); });

        } else {
          EDiff = mcwfn_->StateEnergy->back() - EPrev.back(); 
        }
        
        std::cout << "  (Maximum) Energy Difference = " << std::setw(18) 
                  <<  std::right << EDiff << std::left << std::endl; 
        if(std::abs(EDiff) <= settings->scfEnergyConv) converged = true; 
        
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
        
        if(converged and orbitalGradientNorm < settings->scfGradientConv) break; 
          
        converged = false;
        
        // start of the new cycle 
        std::cout << "\n\nCycle " << iter+1 << ":\n" << std::endl;
        
        // compute hession diagonal and rotate orbitals
        FormattedLine(std::cout, "Performing Orbital Rotation ...");
        
        moRotator->rotateMO(pert, *oneRDMSOI, *twoRDMSOI);
        
        ProgramTimer::tock("Orbital Rotation");
        
        mcwfn_->mointsTF->clearAllCache();
        mcwfn_->moints->clear();

        // print energy and update EPrev
        std::copy_n(mcwfn_->StateEnergy->begin(), this->NStates, EPrev.begin());

        ProgramTimer::tick("Solve CI");
        
        // Re-transform intgrals and solve new CI
        FormattedLine(std::cout, "Redo AO to MO Intergral Transformation ...");
        ProgramTimer::tick("Integral Trans");
        this->transformInts(pert);
        ProgramTimer::tock("Integral Trans");

        std::cout << std::left << std::setprecision(10); 
        FormattedLine(std::cout, "Inactive Energy:", mcwfn_->InactEnergy);
        ProgramTimer::tick("Diagonalization");
        this->ciSolver->solveCI(*mcwfn_, pert);
        ProgramTimer::tock("Diagonalization");
        
        ProgramTimer::tock("Solve CI");
      
        this->printStateEnergy();
      
        saveCurrentStates();
      
      } // SCF Iteration

      ROOT_ONLY(this->comm);

      if(not converged) 
        CErr("\n MCSCF failed to converged in " + std::to_string(settings->maxSCFIter) + " cycles !");
    
      // compute 1RDMs
      if (this->StateAverage) this->computeOneRDM(); 
      else this->computeOneRDM(this->NStates - 1); 
      
      // generate IVOs as needed
      ProgramTimer::tick("Gen IVOs");
      if (this->settings->doIVOs) moRotator->generateIVOs(pert, *oneRDMSOI);   
      ProgramTimer::tock("Gen IVOs");

    } // doSCF

    // Final printing and save states for restart
    std::cout << "\n\nMCSCF Complete!" << std::endl;
    std::cout << bannerEnd << std::endl;

    if(this->settings->NatOrbs)
    {
      this->formNaturalOrbitals();
      if(this->settings->NatOrbRediag)
      {
        // Clear previous transformed integrals
        mcwfn_->moints->clear();      

        // Retransform integrals in new NO basis
        this->transformInts(pert);
        // If we read in the previous CI Vector, use Davison by default
        if(this->ciSolver->getAlg()==SKIP)
          this->ciSolver->switchAlgorithm(CI_DAVIDSON);
        this->ciSolver->solveCI(*mcwfn_, pert);
        this->computeOneRDM();
      }
    }

    this->printMCSCFFooter();

    // property calculation
    ProgramTimer::tick("Property Eval");

    // dipole moment
    if( this->settings->multipoleMoment )
      this->computeMultipole();

    // mulliken analysis
    if (this->settings->PopulationAnalysis) {
      std::cout<<"\n\nPopulation analysis in mcscf."<<std::endl;
      this->populationAnalysis();
    }

    if (mcwfn_->printRDMs==0 && this->settings->SpinAnalysis) {
      std::cout<<"\n\nSpin analysis in mcscf."<<std::endl;
      this->spinAnalysis();
    }    

    // oscillator strength
    if (this->settings->NosS1) {


      mcwfn_->osc_str.resize(this->settings->NosS1*this->NStates);
      for (size_t s1 = 0ul; s1 < this->settings->NosS1; s1++)
      for (size_t s2 = 0ul; s2 < this->NStates; s2++){
        if (s2 <= s1) mcwfn_->osc_str[s2+s1*this->NStates] = 0.;
        else {mcwfn_->osc_str[s2+s1*this->NStates] = 
          this->oscillator_strength(s2,s1);}
      }

    }

    ProgramTimer::tock("Property Eval");

    saveCurrentStates(true);

    ProgramTimer::tock("MCSCF Total");
 
  }; //MCSCF::run
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::saveCurrentStates( bool saveProp ) {
    
    ROOT_ONLY(this->comm);

    mcwfn_->saveCurrentStates(saveProp);
    
    // only save MO when doing orbital rotation
    if (settings->doSCF and this->savFile.exists()) {
      auto mo_dim = mcwfn_->reference().mo[0].nRows();
      this->savFile.safeWriteData("SCF/MO1", mcwfn_->reference().mo[0].pointer(), {mo_dim, mo_dim});
    }

  }; // MCSCF::saveCurrentStates


  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::alloc() {
    
    mcwfn_->savFile = this->savFile;
    mcwfn_->alloc();
    if(mcwfn_->readCI) mcwfn_->ReadGuessCIVector();
    
    ciSolver = std::make_shared<CISolver<MatsT,IntsT>>(settings->ciAlg, 
      settings->maxCIIter, settings->ciVectorConv,
      settings->maxDavidsonSpace, settings->nDavidsonGuess,
      settings->energyRefs);
    
    if (this->settings->doSCF) {
      
      mcwfn_->cacheHalfTransTPI_ = true;
      
      size_t nCorrO = mcwfn_->MOPartition.nCorrO;
      oneRDMSOI = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
      twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);

      // Pass these pointers to the MCWaveFunction Object
      mcwfn_->oneRDMSOI = oneRDMSOI;
      mcwfn_->twoRDMSOI = twoRDMSOI;

      // If doing state-averaging, we need to inform the MCWaveFunction of this primarily
      // for CubeGen of the average density
      if(this->StateAverage)
      {
        mcwfn_->StateAverage = true;
        mcwfn_->SAWeight = this->SAWeight;
      }
      
      moRotator = std::make_shared<OrbitalRotation<MatsT, IntsT>>(*mcwfn_,settings->ORSettings);
    }
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::dealloc() {
    oneRDMSOI = nullptr;
    twoRDMSOI = nullptr;
    ciSolver  = nullptr;
    moRotator = nullptr;
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::transformInts(EMPerturbation & pert) 
  {
    mcwfn_->transformInts(pert);
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::computeMultipole() 
  {
    mcwfn_->computeMultipole();
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::populationAnalysis() 
  {
    mcwfn_->populationAnalysis();
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::spinAnalysis() 
  {
    mcwfn_->spinAnalysis();
  }

  template <typename MatsT, typename IntsT>
  double MCSCF<MatsT,IntsT>::oscillator_strength(size_t s1, size_t s2) 
  {
    return mcwfn_->oscillator_strength(s1,s2);
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::formNaturalOrbitals() 
  {
    mcwfn_->formNaturalOrbs(this->settings->NatOrbs-1);
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::runCube(std::vector<std::shared_ptr<CubeGen>> cu)
  {
    mcwfn_->runCube(cu);
  }

}; // namespace ChronusQ

