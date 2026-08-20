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

#include <configinteraction.hpp>
#include <configinteraction/print.hpp>
#include <configinteraction/rdm.hpp>
#include <configinteraction/solveci.hpp>
#include <util/matout.hpp>
#include <orbitalrotation.hpp>

// #define DEBUG_ConfigInteraction_IMPL
// #define DEBUG_SCALAR_REL_INTS
// #define _DIFF_RDM

namespace ChronusQ {
  
template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::run(EMPerturbation & pert) {
  
  ProgramTimer::tick("Configuration Interaction Total");
   
  std::cout << "configuration interaction starts" << std::endl;

  // TODO: estimate and allocate memory
  this->alloc();

  // initialize DAS space and categories
  this->initialization();
  
  // Initial Printing
  this->printCIHeader();
  
  // detFactory->output(std::cout, "Configuration Interactions");

  // ConfigInteraction Initial CI solution
  std::cout << "Cycle 0:\n" << std::endl;
  FormattedLine(std::cout, "AO to MO Integral Transformation ...");
    
  ProgramTimer::tick("NEW Solve CI");
  
  ProgramTimer::tick("Integral Trans");
  this->prepareMOIntegrals(pert, *detFactory, ciSettings.doSCF, not ciSettings.doSCF);
  ProgramTimer::tock("Integral Trans");
  

  std::cout << std::left << std::setprecision(10); 
  FormattedLine(std::cout, "Core Energy:", this->coreEnergy);

  ProgramTimer::tick("Diagonalization");
  this->solveCI();
  ProgramTimer::tock("Diagonalization");
  
  ProgramTimer::tock("NEW Solve CI");

  // Initial 1RDM construction
  this->computeOneRDM();
  
  // SCF Cycles 
  if(this->ciSettings.doSCF) {

    this->printStateEnergy();
    
    std::vector<double> EPrev = std::vector<double>(this->NStates);
    std::fill_n(EPrev.begin(), this->NStates, 0.);
    double EDiff      = 0.;
    bool   converged  = false;
    
    for (auto iter = 0ul; iter < ciSettings.maxSCFIter; iter++) {
      
      ProgramTimer::tick("Orbital Rotation");
      // Exam Energy
      if (this->StateAverage) {
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
      if(std::abs(EDiff) <= ciSettings.scfEnergyConv) converged = true;
      
      // compute RDMs
      computeRDMsForOrbitalRotations();
      
      // this->print1RDMs();
      
      // compute gradient 
      auto orbGradStart = tick();
      double orbitalGradientNorm = 
        moRotator->computeOrbGradient(pert, *oneRDMSOI, *twoRDMSOI);
      auto durationOrbGrad = tock(orbGradStart);
      
      // Gradient Convergence exam
      // MPI Orbital Gradient only valid on Root so bcast (needed for convergence check)
      MPIBCast(&(orbitalGradientNorm), 1, 0, this->comm);
      std::cout << "  Orbital Gradient Residue   = " << std::setw(18) 
                  << std::right << orbitalGradientNorm << std::left << std::endl;
      std::cout << "  Orbital Gradient Duration  = " << std::setw(18)
                   << std::right << durationOrbGrad  << " s " << std::left << std::endl;
      
      if(converged and orbitalGradientNorm < ciSettings.scfGradientConv) break;
        
      converged = false;
      
      // start of the new cycle 
      std::cout << "\n\nCycle " << iter+1 << ":\n" << std::endl;
      
      // compute hession diagonal and rotate orbitals
      FormattedLine(std::cout, "Performing Orbital Rotation ...");
      
      moRotator->rotateMO(pert, *oneRDMSOI, *twoRDMSOI);
      // communicate new orbitals
      auto& mo = this->reference()->mo[0];
      auto mo_pointer = mo.pointer();
      size_t bcast_size = mo.nRows() * mo.nColumns();
      MPIBCast(mo_pointer,bcast_size,0,this->comm);
      
      ProgramTimer::tock("Orbital Rotation");
      
      this->mointsTF->clearAllCache();

      // print energy and update EPrev
      std::copy_n(this->StateEnergy.begin(), this->NStates, EPrev.begin());

      ProgramTimer::tick("NEW Solve CI");
      
      // Re-transform intgrals and solve new CI
      FormattedLine(std::cout, "Redo AO to MO Intergral Transformation ...");
      ProgramTimer::tick("Integral Trans");
      this->prepareMOIntegrals(pert, *detFactory, true, false);
      ProgramTimer::tock("Integral Trans");

      std::cout << std::left << std::setprecision(10); 
      FormattedLine(std::cout, "Core Energy:", this->coreEnergy);
      ProgramTimer::tick("Diagonalization");
      this->solveCI();
      ProgramTimer::tock("Diagonalization");
      
      ProgramTimer::tock("NEW Solve CI");
    
      this->printStateEnergy();
    
      saveCurrentStates();
    
    } // SCF Iteration
    
    if(not converged) 
      CErr("\n ConfigInteraction failed to converged in " + std::to_string(ciSettings.maxSCFIter) + " cycles !");
  
    // compute 1RDMs
    if (this->StateAverage) this->computeOneRDM(); 
    else this->computeOneRDM(this->NStates - 1); 
    
    // generate IVOs as needed
    ProgramTimer::tick("Gen IVOs");
    if (this->ciSettings.doIVOs) moRotator->generateIVOs(pert, *oneRDMSOI);
    ProgramTimer::tock("Gen IVOs");

  } // doSCF

  // Final printing and save states for restart
  std::cout << "\n\nConfiguration Interaction Complete!" << std::endl;
  std::cout << bannerEnd << std::endl;
  
  this->printCIFooter();

  // // property calculation

  // Mulliken analysis
  if (this->PopulationAnalysis) {
    std::cout<<"\n\nPopulation analysis in ConfigInt."<<std::endl;
    PostHartreeFock<MatsT,IntsT>::populationAnalysis();
  }

  // save OnePDMS
  if (this->saveOnePDMS) {
    std::cout<<"\n\nSaving One PDMs in ConfigInt."<<std::endl;
    PostHartreeFock<MatsT,IntsT>::saveOnePDMs();
  }

  // Spin analysis
  if (this->printRDMs==0 && this->SpinAndAngularAnalysis) {
      std::cout<<"\n\nSpin analysis in ConfigInt."<<std::endl;
        PostHartreeFock<MatsT,IntsT>::spinAndAngularAnalysis();
  }   

  // Compute excited state transtition dipole moments:
  ProgramTimer::tick("Property Eval");
  if (this->printTransDipole) {
    this->printAllTransitionDipoleMom(); }

  if (this->printDipole) {
    this->printAllStateSpecificDipoleMom(); }

  // oscillator strength
  if (this->osc_str) {
    
    auto &ref = *this->reference();
    this->osc_str_array.reserve(this->NosS1*this->NStates);
    for (size_t s1 = 0ul; s1 < this->NosS1; s1++)
    for (size_t s2 = this->NosS1; s2 < this->NStates; s2++) {
      if (this->osc_str_order == 0 && ref.nC < 4) {
        this->osc_str_array 
          .push_back(PostHartreeFock<MatsT,IntsT>::oscillator_strength(s2, s1));
      }
      else if (this->osc_str_order == 2 && ref.nC < 4) {
        this->osc_str_array 
          .push_back(PostHartreeFock<MatsT,IntsT>::secondorder_oscillator_strength(s2, s1));
      }
      else if (this->osc_str_order == 0 && ref.nC == 4) {
        this->osc_str_array 
          .push_back(PostHartreeFock<MatsT,IntsT>::oscillator_strength4C(s2, s1));
      }
      else { 
        CErr("OSCISTRENGTH ORDER NYI!");
      }
    }

#ifdef _DIFF_RDM   
    if (this->NStates > 1) PostHartreeFock<MatsT,IntsT>::OneRDMDiff();
#endif

  }
  ProgramTimer::tock("Property Eval");
  saveCurrentStates();

  ProgramTimer::tock("Configuration Interaction Total");

} //ConfigInteraction::run

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::saveCurrentStates() {
  
  ROOT_ONLY(this->comm);

  PostHartreeFock<MatsT,IntsT>::saveCurrentStates();
  
  // only save MO when doing orbital rotation
  if (ciSettings.doSCF and this->savFile.exists()) {
    auto mo_dim = this->reference()->mo[0].nRows();
    this->savFile.safeWriteData("SCF/MO1", this->reference()->mo[0].pointer(), {mo_dim, mo_dim});
  }

} // ConfigInteraction::saveCurrentStates

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::alloc() {
  // TODO: need to estimate memory usage
}

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::initialization() {
  
  PostHartreeFock<MatsT,IntsT>::alloc();
   
  // build detFactory
  ProgramTimer::tick("Determinant Factory");
  std::cout << "Initializing Determinant Factory" << std::endl;

  // TODO: add automate mechanism of break down large space to smaller spaces
  // Initialize a DetFactory object with the reference active spaces
  detFactory = std::make_shared<DeterminantFactory>(this->comm,
                                                    this->corrSpace.nCorrE, ciSettings.activeSpaces);
  
  // make sure the active space partitioning in CategoricalSpace is in reference to
  // detFactory's active spaces
  // newCategoricalSpace is an object of CategoricalSpace class.
  // std::cout << "Construct Categorical Space" << std::endl;
  auto buildCat = tick();
  auto newCategoricalSpace = detFactory->buildEmptyDetsSpace();

  // build and add user defined categories based on the reference occupation number
  for (auto const & refOccupation: ciSettings.refOcc) {
    newCategoricalSpace->addReferenceCategory(refOccupation);
  }
  // build new categories based on excitation operator
  newCategoricalSpace->expandCategory(ciSettings.maxInterSpaceEx);
  auto durationBuildCat = tock(buildCat);
  // std::cout << "Construct Categorical Space Done: " << durationBuildCat << " s"<<std::endl;
  // std::cout << "Total Number of Categories = "<<newCategoricalSpace->nCategories()<<std::endl;
  std::cout << "Total Number of Determinants = "<<newCategoricalSpace->nDeterminants()<<std::endl;

  // Check and modify NStates to Ndets if NStates > Ndets
  if (this->NStates > newCategoricalSpace->nDeterminants()) {
    this->NStates = newCategoricalSpace->nDeterminants();
    std::cout << "Requested Number of States > NDeterminants :: Modifying NStates to NDeterminants!" 
      << "\nNew NRoots = " << this->NStates << "\n" << std::endl;
  }

  // Check Number of Threads:
  size_t nThreads = GetNumThreads();
  if (nThreads > newCategoricalSpace->nDeterminants()) {
    SetNumThreads(newCategoricalSpace->nDeterminants());
    std::cout << "WARNING: Number of Threads set to -> " << GetNumThreads() <<
      std::endl;
  }

  // For MPI
  newCategoricalSpace->initializeDistributedCatMap(this->comm);
  detFactory->setKetCategoricalSpace(newCategoricalSpace);
  detFactory->setBraCategoricalSpace(newCategoricalSpace);

  // based on Bra and Ket categories, figures out non-zero excitation maps
  // in terms of active spaces between categories.
  std::cout << "Computing Graph" << std::endl;
  auto buildGraph = tick();
  detFactory->generateComputingGraph();
  auto durationBuildGraph = tock(buildGraph);
  std::cout << "Computing Graph Done: " << durationBuildGraph << " s"<<std::endl;

  auto computeExList = tick();
  detFactory->estimateMemoryRequirement();
  detFactory->computeExcitationList();
  auto durationComputExList = tock(computeExList);
  std::cout << "Computing Excitation List Done: " << durationComputExList << " s"<<std::endl;

  std::cout << "DAS initialization done in Determinant Factory!" << std::endl;
  ProgramTimer::tock("Determinant Factory");
  
  // allocate CI vector
  size_t NS = this->NStates;


  if(ciSettings.SparseDavidson) {
#ifdef CQ_ENABLE_SPARSE
      CIVectors = newCategoricalSpace->constructDistributedSparseCIVectors<MatsT>(this->comm, NS);
#else
      CErr("ENABLE_SPARSE was turned off during compilation! Recompile with CQ_ENABLE_SPARSE");   
#endif
  }
  else {
    CIVectors = newCategoricalSpace->constructDistributedCIVectors<MatsT>(this->comm, NS);
  }

  auto dasciBuilder = std::make_shared<DASCIBuilder<MatsT>>(this->comm, this->moints, *detFactory);
  dasciBuilder->setSigma2eContractionAlgorithm(ciSettings.ciSigma2eContAlg);

  ciBuilder = dasciBuilder;

  if (this->ciSettings.doSCF) {
    size_t nCorrO = this->corrSpace.nCorrO;
    oneRDMSOI = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);     
    moRotator = std::make_shared<NewOrbitalRotation<MatsT, IntsT>>(
      dynamic_cast<PostHartreeFock<MatsT,IntsT>&>(*this), ciSettings.ORSettings);
  }
}

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::dealloc() {
  oneRDMSOI = nullptr;
  twoRDMSOI = nullptr;
  ciBuilder = nullptr;
  detFactory = nullptr;
  CIVectors = nullptr;
  // moRotator = nullptr;
}


} // namespace ChronusQ

