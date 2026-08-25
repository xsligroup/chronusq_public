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

#include <posthartreefock.hpp> 
#include <newperturb.hpp>
#include <newcibuilder/dascibuilder.hpp>
#include <detstringmanager.hpp>
#include <mointstransformer/impl.hpp>

#include <util/matout.hpp>
// #define _DEBUG_PTCATBUILD

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::run(EMPerturbation & pert) {

    ProgramTimer::tick("MRPT Total");
    // Initialize MRPT Partitions:  
    std::cout << BannerTop << std::endl; 
    std::cout << "\n*** Initialize MRPT2 ***" << std::endl;
    PTInitialize();
    printMRPTHeader();

    // Perform AO -> MO Integral transforms:
    ProgramTimer::tick("Integral Trans");
    if (PTopts.SPINFREE) {
      auto& HOp = RefMCWfn_->reference()->fockBuilder->hamiltonianOptions_;
      HOp.SpinFree = PTopts.SPINFREE;
    }
    this->prepareMOIntegrals(pert, *PTFactory_, true, PTopts.ENPT);
    ProgramTimer::tock("Integral Trans");
    
    // Build CI Builder for PT2:
    ptBuilder_ = std::make_shared<DASCIBuilder<MatsT>>(this->comm, this->moints, *PTFactory_);
  	ptBuilder_->setSigma2eContractionAlgorithm(RefMCWfn_->ciSettings.ciSigma2eContAlg);
    std::cout << std::endl << bannerTop << std::endl;

    if (PTopts.ENPT) {
      CIEnergy();
      for (auto& state_index : Target_States_)
#ifdef CQ_ENABLE_SPARSE
        computeEN2Sparse(state_index);
#else
        computeEN2(state_index);
#endif
    }
    
    else if (PTopts.GVVPT) {
      bool sparseFlag = RefMCWfn_->ciSettings.SparseDavidson;
      if (sparseFlag) {
#ifdef CQ_ENABLE_SPARSE
        diagEffHSparse();
#endif
      }
      else {
        // test semicanonicalization:
        diagEffH(); 
      }   
    }
  
    else CErr("Invalid MRPT flavor requested!!");

    // Final MRPT results:
    printMRPTFooter();
    // Save MRPT2 Energies to bin:
    saveCurrentStates();
    // Deallocate Space
    dealloc();

    ProgramTimer::tock("MRPT Total");

  } // DasPerturb::run()


  /**
   *
   *  \brief Initialize mcwfn for Perturb calculation
   *         based on reference mcwfn.
   *         
   */
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::PTInitialize() {

    // Get Reference MO Space and set PT MO Space
    ProgramTimer::tick("PT Category Build");
    const auto& refMOSpace = RefMCWfn_->corrSpace;
    this->setupCorrelatedMOSpace(refMOSpace.nMO-refMOSpace.nNegMO-
    PTopts.FROZENCORE-PTopts.FROZENVIRTUAL, refMOSpace.nCorrE+
    refMOSpace.nInact-PTopts.FROZENCORE, PTopts.FROZENCORE, 
    PTopts.FROZENVIRTUAL);
    const auto& mrptMOSpace = this->corrSpace;
    const size_t nInact = refMOSpace.nInact - mrptMOSpace.nFCore;
    const size_t nCorrO = refMOSpace.nCorrO;
    const size_t nSVirt = refMOSpace.nSVirt - mrptMOSpace.nFVirt;
    const size_t nCorrE = mrptMOSpace.nCorrE;

    // Define PT-DAS Space creation specific variables:
    size_t offsetInact = mrptMOSpace.nNegMO + mrptMOSpace.nFCore;
    size_t offsetSVirt = offsetInact + nInact + nCorrO;

    // Declare PT Space build objects:
    std::vector<size_t> nCoreDASOrb, nVirtDASOrb;
    std::vector<ActiveSpaceParameters> ptSpaces;
    std::vector<std::vector<size_t>> Occ(1);

    // Efficient splitting of Core and Virtual Spaces for
    // DAS space building:
    size_t nDasCoreSpace, nDasVirtSpace;
    if (nInact > 0) { 
      nDasCoreSpace = (nInact / 10) + 1; }
    else { 
      nDasCoreSpace = 0; }
    if (nSVirt > 0) {
      nDasVirtSpace = (nSVirt / 10) + 1; }
    else { 
      nDasVirtSpace = 0; }

    if ( nDasCoreSpace == 0 && nDasVirtSpace == 0)
      CErr("PT2 Space is Null! Please recheck your orbital spaces!!");

    // Initialize the core and virtual DAS Spaces that will enter the 
    // PT Space
    if (nInact > 0) { 
      nCoreDASOrb.resize(nDasCoreSpace, nInact/nDasCoreSpace);
      for (size_t i=0; i < nInact % nDasCoreSpace; i++) 
        nCoreDASOrb[i] += 1; }
    if (nSVirt > 0) {
      nVirtDASOrb.resize(nDasVirtSpace, nSVirt/nDasVirtSpace);
      for (size_t i=0; i < nSVirt % nDasVirtSpace; i++) 
        nVirtDASOrb[i] += 1; }

    // Build PT occupations according to the Core, active and virtual occs:
    // use active space occs 
    ptSpaces.resize(RefMCWfn_->ciSettings.activeSpaces.size() + 
      nCoreDASOrb.size() + nVirtDASOrb.size());
    ptSpaces.clear();

    // Build the PT spaces and occupations:
    size_t DASGroup = 0;
    if (nCoreDASOrb.size() > 0) {
      for (auto &i: nCoreDASOrb) {
        size_t maxHole = std::min(i,size_t(2));
        maxHole = std::min(nCorrO-nCorrE, maxHole);
        ptSpaces.push_back(ActiveSpaceParameters({offsetInact, 
          i, i, i - maxHole, i, DASGroup, 0, int(maxHole), "DASCore", false}));
        offsetInact += i; }
      Occ[0] = nCoreDASOrb; }
    
    for (auto &i: RefMCWfn_->ciSettings.activeSpaces) {
      i.iDASGroup += DASGroup + 1;
      ptSpaces.push_back(i); }
    Occ[0].insert(Occ[0].end(), 
      RefMCWfn_->ciSettings.refOcc[0].begin(),
      RefMCWfn_->ciSettings.refOcc[0].end());
    
    if (nVirtDASOrb.size() > 0) {
      DASGroup = RefMCWfn_->ciSettings.activeSpaces.back().iDASGroup + 1;
      for (auto &i: nVirtDASOrb) {
        size_t maxPart = std::min(i,size_t(2));
        ptSpaces.push_back(ActiveSpaceParameters({offsetSVirt, 
          i, 0, 0, maxPart, DASGroup, int(maxPart), 0, "DASVirt", false})); 
        offsetSVirt += i; }
      Occ[0].insert(Occ[0].end(), nVirtDASOrb.size(), 0); }


    // Initialize a DetFactory object with the Ref+PT active spaces
    auto tempFactory = std::make_shared<DeterminantFactory>(this->comm,
      nCorrE, ptSpaces);
    auto tempCategoricalSpace = tempFactory->buildEmptyDetsSpace();
    for (auto const &Occupation: Occ) 
      tempCategoricalSpace->addReferenceCategory(Occupation);
    tempCategoricalSpace->expandCategory(-1);
    
    // Initialize PT and Ref factories
    PTFactory_ = std::make_shared<DeterminantFactory>(this->comm,
      nCorrE, ptSpaces);
    auto ptcategoricalSpace = PTFactory_->buildEmptyDetsSpace();
    auto refcategoricalSpace = PTFactory_->buildEmptyDetsSpace();
 
    // Adding reference categories to pt and ref categorical spaces:
    for (size_t i = 0ul; i < tempCategoricalSpace->nCategories(); ++i) {
      auto cat = tempCategoricalSpace->getCategory(i);
      auto spaceOcc = cat->SpaceOccupations();
      size_t tmpCoreDAS = std::accumulate(spaceOcc.begin(), 
        spaceOcc.begin() + nCoreDASOrb.size(), 0);
      size_t tmpVirtDAS = std::accumulate(spaceOcc.begin() + spaceOcc.size() -
        nVirtDASOrb.size(), spaceOcc.end(), 0);
      if (tmpCoreDAS == nInact && tmpVirtDAS == 0) {
        refcategoricalSpace->addReferenceCategory(spaceOcc);}
      else {
        ptcategoricalSpace->addReferenceCategory(spaceOcc);}

    }
    tempFactory.reset();
    
#ifdef _DEBUG_PTCATBUILD
    refcategoricalSpace->output(std::cout);
#endif
    // Restore Nthreads:
    size_t nThreads = RefMCWfn_->ciSettings.nThreads_;
    if (GetNumThreads() < nThreads) {
      SetNumThreads(nThreads);
      std::cout << "WARNING: Number of threads restored to: " << GetNumThreads() << std::endl;
    }

    // Build the Computation graphs and excitation lists for <Ref|H|PT>
    // calculations
    refcategoricalSpace->initializeDistributedCatMap(this->comm);
    this->PTFactory_->setKetCategoricalSpace(refcategoricalSpace);
    if (ptcategoricalSpace->nCategories() != 0) {
      ptcategoricalSpace->initializeDistributedCatMap(this->comm);
      this->PTFactory_->setBraCategoricalSpace(ptcategoricalSpace);}
    else { CErr("Perturbation Space is ill-defined!!"); }

    this->PTFactory_->generateComputingGraph(true);
    this->PTFactory_->estimateMemoryRequirement();
    this->PTFactory_->computeExcitationList();
    std::cout << "*** Initialized MRPT2 *** " << std::endl;
     
    ProgramTimer::tock("PT Category Build");


  } // DasPerturb::PTInitialize


  // Save current states
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::saveCurrentStates() {
    
    ROOT_ONLY(this->comm);
    if(this->savFile.exists()) {
    
      size_t NS = this->E2_.size();
      this->savFile.safeWriteData("MRPT2/SOI", &NS, {1});
      this->savFile.safeWriteData("MRPT2/ZERO_ENERGY", E0_.data(), {NS});
      this->savFile.safeWriteData("MRPT2/MRPT2_ENERGY", E2_.data(), {NS});
    }
  
  }

  // Allocate and deallocate memory
  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::alloc() {

    PostHartreeFock<MatsT,IntsT>::alloc();

  }; // DasPerturb::alloc()

  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::dealloc() {

    // Dealloc shared ptrs and vectors:
    PTFactory_.reset();
    ptBuilder_.reset();
    RefMCWfn_.reset();
    Target_States_.clear();
    fockDiag_.reset();
    E0_.clear();
    E2_.clear();

  }


} // namespace ChronusQ

// Include implementation files:
#include <newperturb/print.hpp>
#include <newperturb/enpt.hpp>
#include <newperturb/gvvptsparse.hpp>
#include <newperturb/gvvpt.hpp>
#include <newperturb/gvvptutil.hpp>
