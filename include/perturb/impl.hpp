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
#include <perturb.hpp>
#include <cibuilder/rasci.hpp>
#include <detstringmanager.hpp>
#include <mointstransformer/impl.hpp>

#include <util/matout.hpp>

//#define _DEBUG_PT2_impl

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::run(EMPerturbation & pert) {

    ProgramTimer::tick("PERTURB Total");
    std::cout << "\n\nLL PT2 Module" << std::endl;

    // Initialize MCWaveFunction object
    PTInitialize();

    printPTHeader();

    // transform integrals
    ProgramTimer::tick("Integral Trans");
    this->transformInts(pert);
    ProgramTimer::tock("Integral Trans");

    // eXtend Multi-State: form intermediate states
    if(PTopts.extendMS) rotateCIVecs();

    // get 1RDM from ref
    refMCwfn->computeOneRDM(); // TODO: now computing for all MCSCF states...
#ifdef _DEBUG_PT2_impl
    refMCwfn->print1RDMs();
#endif
    std::cout<<"computed 1RDM."<<std::endl;


    // use averaged 1RDM or not
    if (PTopts.saFock) {
      auto & oneRDM = refMCwfn->oneRDMSOI;
      std::cout<<"build state-averaged Fock: "<<std::endl;
      buildFock(*ptFock_, *oneRDM);

      if (PTopts.doIter) {
        std::cout<<"build diag LHS"<<std::endl;
        buildDiag(); // build diag for preconditioner using averaged Fock
      }
    }

    size_t nStates = this->NStates;
    // solve CASPT2 for each reference state
    for (auto i = 0ul; i < nStates; i++) {
      std::cout<<"\n\nSingle state: "<<SoI_[i]<<std::endl;
      runSingleState(i); 
    }

    if (nStates > 1) { // Multi-state calculations

      if (not PTopts.doGVV) {
        // build H_eff full matrix
        std::shared_ptr<cqmatrix::Matrix<MatsT>> H_eff =
                            std::make_shared<cqmatrix::Matrix<MatsT>>(nStates);
        buildHeff(*H_eff); 
        // diagonalize H_eff
        diagHeff(*H_eff);
      } else { // sigma build for GVVPT2 type
        // CASCI
        CErr("GVV NYI");
      }

    }

    else { // single state perturbed energy

      MatsT * HV = CQMemManager::get().malloc<MatsT>(this->NDet);
      computeHV(HV, this->CIVecs[0]);

      this->StateEnergy[0] = refMCwfn->StateEnergy[SoI_[0]]
                        + std::real(computeCHV(this->CIVecs[0], HV))
                        - std::real(computeShiftCorrection(0));
      CQMemManager::get().free(HV);
    }

    this->saveCurrentStates();
    this->printPTFooter();

    ProgramTimer::tock("PERTURB Total");

  } // PERTURB::run()

  /**
   * 
   * \brief solve CASPT2 for each reference state
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::runSingleState(size_t i) {

    auto & oneRDM = refMCwfn->oneRDM[SoI_[i]];

    // build ptFock for state-specific case
    if (not PTopts.saFock) {
      buildFock(*ptFock_, oneRDM);
      if (PTopts.doIter) buildDiag();
    }

    // compute zero-order energy for each state with state-specific or sa Fock
    if (not PTopts.extendMS) E0_[i] = computeZeroE(*ptFock_, oneRDM);
    
    // build RHS
    MatsT * RHS = CQMemManager::get().malloc<MatsT>(SDsize);
    buildRHS(RHS, i);

    // solve the linear equation
    solve(RHS, i);

    // Copy over RHS to ptV
    std::copy_n(RHS, SDsize, this->CIVecs[i]+CIsize);
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"LL PT2 RHS from solver", RHS, SDsize, 1, SDsize);
#endif


    CQMemManager::get().free(RHS);

  } // PERTURB::run()
    

  /**
   *
   *  \brief Initialize mcwfn for Perturb calculation
   *         based on reference mcwfn.
   *         
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::PTInitialize() {

    auto & mopart_ref = refMCwfn->MOPartition;
    // 1RDMs setting in ref
    if(PTopts.saFock) {
      refMCwfn->StateAverage = true;
      refMCwfn->SAWeight = std::vector<double>(refMCwfn->NStates);
      // setting equal weight for state-averaged 1RDM/Fock
      for (auto i = 0ul; i < refMCwfn->NStates; i++) {
        if(std::count(SoI_.begin(), SoI_.end(), i))
          refMCwfn->SAWeight[i] = 1./ this->NStates;
        else refMCwfn->SAWeight[i] = 0.;
        std::cout<<"state "<<i<<" weight: "<<refMCwfn->SAWeight[i]<<std::endl;
      }
      if ( not refMCwfn->oneRDMSOI)
        refMCwfn->oneRDMSOI = std::make_shared<cqmatrix::Matrix<MatsT>>(mopart_ref.nCorrO);
    }

    this->FourCompNoPair = refMCwfn->FourCompNoPair;
    auto & mopart = this->MOPartition;
    // set RAS scheme
    mopart.scheme = RAS;
    std::vector<size_t> nActO(3, 0);
    nActO[0] = mopart_ref.nInact - PTopts.frozenCore; // RAS 1
    nActO[1] = mopart_ref.nCorrO; // CAS space
    nActO[2] = mopart_ref.nElecMO - PTopts.frozenCore - PTopts.frozenVirtual
                - nActO[0] - nActO[1]; // RAS 3
    mopart.mxHole = std::min(2ul,nActO[0]);
    mopart.mxHole = std::min(std::accumulate(nActO.begin(), nActO.end(), 0)
                - (refMCwfn->referenceWaveFunction().nO - PTopts.frozenCore), mopart.mxHole);
    mopart.mxElec = std::min(2ul,nActO[2]);
    mopart.mxElec = std::min(refMCwfn->referenceWaveFunction().nO - PTopts.frozenCore,
                    mopart.mxElec);
    this->partitionMOSpace(nActO, refMCwfn->referenceWaveFunction().nO - PTopts.frozenCore);

    // select virtuals based on input
    if (not PTopts.selectVirtual.empty()) {
      std::vector<char> inputOrbIndices(mopart.orbIndices);
      size_t virOffset = mopart_ref.nInact + mopart_ref.nCorrO + mopart_ref.nNegMO;
      std::fill(inputOrbIndices.begin()+virOffset, inputOrbIndices.end(), 'N');
//      print_orbIndices(inputOrbIndices);

      set_orbital_index(inputOrbIndices, PTopts.selectVirtual, '3');
//      std::cout<<"nFVirt: "<<mopart.nFVirt<<std::endl;
//      std::cout<<"frozenvirtual: "<<PTopts.frozenVirtual<<std::endl;

      fill_default_index(inputOrbIndices, virOffset, 'S', mopart.nFVirt);
      mopart.orbIndices = inputOrbIndices;
      print_orbIndices(mopart.orbIndices);
      this->setActiveSpaceAndReOrder();

    }



#ifdef _DEBUG_PT2_impl
    std::cout<<"nActOs[0]: "<<mopart.nActOs[0]<<", [1]: "<<mopart.nActOs[1]
        <<", [2]: "<<mopart.nActOs[2]<<std::endl;
    std::cout<<"nCorrE: "<<mopart.nCorrE<<std::endl;
    std::cout<<"nCorrO: "<<mopart.nCorrO<<std::endl;
    std::cout<<"nMO: "<<mopart.nMO<<std::endl;
    std::cout<<"nElecMO: "<<mopart.nElecMO<<std::endl;
    std::cout<<"mxHole: "<<mopart.mxHole<<std::endl;
    std::cout<<"mxElec: "<<mopart.mxElec<<std::endl;
    std::cout<<"nInact: "<<mopart.nInact<<std::endl;
    std::cout<<"nFVirt: "<<mopart.nFVirt<<std::endl;
#endif

    CIsize = refMCwfn->NDet;
    SDsize = this->NDet - refMCwfn->NDet;

#ifdef _DEBUG_PT2_impl
    std::cout<<"NDet: "<<this->NDet<<std::endl;
    std::cout<<"CIsize: "<<CIsize<<std::endl;
    std::cout<<"SDsize: "<<SDsize<<std::endl;
#endif

    alloc();
       

  } // PERTURB::PTInitialize

  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::saveCurrentStates() {

    MCWaveFunction<MatsT, IntsT>::saveCurrentStates(false);

  } // PERTURB::saveCurrentStates


  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::alloc() {

    MCWaveFunction<MatsT,IntsT>::alloc();

    size_t nCorrO = this->MOPartition.nCorrO;

    E0_ = CQMemManager::get().malloc<MatsT>(this->NStates);
    ptFock_ = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);

    if (PTopts.doFull)
      // LHS = F - E_0  dimension: SDsize x SDsize
      LHS_ = std::make_shared<cqmatrix::Matrix<MatsT>>(SDsize);

    if (PTopts.doIter) diagLHS = CQMemManager::get().malloc<MatsT>(SDsize);
    
    std::cout<<"Perturb allocation finished."<<std::endl;

  }; // PERTURB::alloc()

  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::dealloc() {

    refMCwfn = nullptr;
    CQMemManager::get().free(E0_);
    E0_ = nullptr;
    LHS_ = nullptr;
    ptFock_ = nullptr;
    if (diagLHS) {
      CQMemManager::get().free(diagLHS);
      diagLHS = nullptr;
    }
  }

} // namespace ChronusQ

// Other implementation files
#include <perturb/print.hpp>
#include <perturb/sspt2.hpp>
#include <perturb/mspt2.hpp>
#include <perturb/solver.hpp>


