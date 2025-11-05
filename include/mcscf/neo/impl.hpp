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
#include <mcscf/neo/print.hpp>
#include <mcscf/neo/rdm.hpp>
#include <mcscf/neo/property.hpp>
#include <mcscf/neo/cube.hpp>
#include <orbitalrotation.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void NEOMCSCF<MatsT, IntsT>::transformMultipleInts(EMPerturbation & pert)
  {

    // Precompute field - nuclear moment contributions 
    this->precompute_NucEField(pert);

    // SMG 09/11/13
    // This mostly follows the structure defined in MCWavefunction defined in
    // /include/mcwavefunction/moints.hpp

    eeTF->setMORanges(*(this->ewfn_));

    // SMG 09/13/23
    // This is major hacky to get the protonic wavefunction to 
    // appropriately partition the MO space
    PPTF->setMORanges(*(this->pwfn_));

    // If the field has changed, we need to rebuild the AOHCore cache
    const double FIELD_DIFF_EPSILON = 1e-15;
    if( pert_has_type(pert,Electric) ) {
      std::array<double, 3> dip_field = pert.getDipoleAmp(Electric);
      double diff = 0.0;
      for(auto iXYZ = 0;    iXYZ < 3;     iXYZ++){
        diff += std::pow(dip_field[iXYZ] - this->old_dip_field[iXYZ], 2.0);
        this->old_dip_field[iXYZ] = dip_field[iXYZ];
      }
      diff = std::sqrt(diff);
      this->field_changed = (diff > FIELD_DIFF_EPSILON);
    }
    if(this->field_changed)
    {
      eeTF->clearAllCache(); 
      PPTF->clearAllCache();
    }

    size_t nTEOrb = this->ewfn_->MOPartition.nMO;
    size_t nCorrEO = this->ewfn_->MOPartition.nCorrO;
    size_t nInactE = this->ewfn_->MOPartition.nInact;
    size_t nInactE2 = nInactE * nInactE;

    size_t nTPOrb = this->pwfn_->MOPartition.nMO;
    size_t nCorrPO = this->pwfn_->MOPartition.nCorrO;
    size_t nInactP = this->pwfn_->MOPartition.nInact;
    size_t nInactP2 = nInactP * nInactP;

    size_t enCorrO = nCorrEO;
    size_t PnCorrO = nCorrPO;

    // TODO: Calculated the core energy, and appropriately transformed 
    // one electron integrals

    // Storage for the transformed matrices
    OnePInts<MatsT> ehCore_tu(enCorrO);
    OnePInts<MatsT> PhCore_tu(PnCorrO);
    OnePInts<MatsT> ehCoreP_tu(enCorrO);
    OnePInts<MatsT> PhCoreP_tu(PnCorrO);

    // Actual storage for the TEI's moving forward
    std::shared_ptr<InCore4indexTPI<MatsT>> eERI_tuwv = this->moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space");
    std::shared_ptr<InCore4indexTPI<MatsT>> PERI_tuwv = this->moints->template getIntegral<InCore4indexTPI,MatsT>("PRI_Correlated_Space");
    std::shared_ptr<InCore4indexTPI<MatsT>> ePERI     = this->moints->template getIntegral<InCore4indexTPI,MatsT>("eP_Correlated_Space");

    // Temporary storage for Two particle integrals
    // Lifetime is this function scope
    InCore4indexTPI<MatsT> PP_ijkl(nInactP);
    InCore4indexTPI<MatsT> PP_pqrs(nTPOrb);
    InCore4indexTPI<MatsT> eP_eCore_PActive(nInactE,PnCorrO);
    InCore4indexTPI<MatsT> eP_PCore_eActive(nCorrEO,nInactP);
    InCore4indexTPI<MatsT> eP_eCore_PCore(nInactE,nInactP);

    // Calculate the core energy
    MatsT * h1e_ii = CQMemManager::get().template malloc<MatsT>(nInactE);
    MatsT * h1p_ii = CQMemManager::get().template malloc<MatsT>(nInactP);
    MatsT * GDejj_ii = CQMemManager::get().template malloc<MatsT>(nInactE);
    MatsT * GDpjj_ii = CQMemManager::get().template malloc<MatsT>(nInactP);

    // Electronic Hamiltonian
    // Transform the oei's and get appropriate G[D] for the core orbitals
    eeTF->transformHCore(pert,h1e_ii,"ii",true);
    eeTF->transformGD(pert,'i',GDejj_ii,"jj",true,true,"WithInactive-i");
    //eeTF->transformGD(pert,'i',GDejj_ii,"jj",true);
    //PPTF->transformGD(pert,'i',GDpjj_ii,"jj",true);

    // Calcualte the electron only core energy
    MatsT ECore=0.;
    for(auto i = 0ul; i < nInactE; i++)
    {
      ECore+=h1e_ii[i]+0.5*GDejj_ii[i]; // + GDepjj_ii[i];
    }
    
    // Account for double occupation of core orbitals
    ECore *= 2.0;

    // Proton-Proton Hamiltonian
    // Note here were don't use transformGD since that assumes closed shell
    // but the protonic wavefunction is high spin
    //PPTF->transformGD(pert,'i',GDpjj_ii,"jj",true,true,"WithInactive-i");
    if(nInactP)
    {
      PPTF->transformHCore(pert,h1p_ii,"ii",true);
      PPTF->transformTPI(pert,PP_ijkl.pointer(),"ijkl");
      std::fill_n(GDpjj_ii,nInactP,MatsT(0.0));
      for(size_t i = 0; i < nInactP; i++)
      {
        for(size_t j = 0; j < nInactP; j++)
        {
          GDpjj_ii[i] += PP_ijkl(i,i,j,j)-PP_ijkl(i,j,i,j);
        }
      }
      for(auto i = 0ul; i < nInactP; i++)
      {
        ECore+=h1p_ii[i]+0.5*GDpjj_ii[i]; 
      }

      // Transform the (ij|tu) integrals for core Protonic Hamiltonian
      // This is an issue right now in that we just transform all the integrals to then
      // pick back out the pieces we care about
      if(nInactP)
      {
        PPTF->transformTPI(pert,PP_pqrs.pointer(),"pqrs");
      }
    }

    // Do the asymmetric integral transforms for the two different
    // core-active partitions
    if(nInactP && nInactE)
    {

      std::cout << "------------------------------------------------------"  << std::endl;
      std::cout << "Detected both inactive electrons and inactive protons."  << std::endl;
      std::cout << "This is okay, but will fully transform the two        "  << std::endl;
      std::cout << "particle (ee|PP) integral set! This may be painfully  "  << std::endl;
      std::cout << "slow pending your basis size!                         "  << std::endl;
      std::cout << "------------------------------------------------------"  << std::endl;

      ePTF->transformAsymmTPI(pert,eP_eCore_PCore.pointer(),eeTF,PPTF,"ijkl","ijkl",false);
      for(size_t i = 0; i < nInactE; i++)
      {
        for(size_t I = 0; I < nInactP; I++)
        {
          ECore -= 2.*eP_eCore_PCore(i,i,I,I);
        }
      }
    }
    this->InactEnergy = std::real(ECore);

    // Transform the one particle integrals in the correlated orbitals
    eeTF->transformHCore(pert,ehCore_tu.pointer(),"tu",false,'i');
    // Don't use the clever core tricks here since again that assumed double occupation
    // Not sure why this needs to be done, but the 1 particle integrals are NOT trasnformed
    // correctly in the presence of a field if the PPTF cache is not cleared
    PPTF->clearAllCache();
    PPTF->transformHCore(pert,PhCore_tu.pointer(),"tu");
    //ehCore_tu.output(std::cout,"Electron Core Hamiltoniai",true);
    //PhCore_tu.output(std::cout,"Proton Core Hamiltonian",true);

    // Transform the TPI's for the e-e and p-p subsystems
    if(!eERI_tuwv)
    {
      eERI_tuwv = std::make_shared<InCore4indexTPI<MatsT>>(nCorrEO);
      eeTF->transformTPI(pert,eERI_tuwv->pointer(),"tuvw");
      this->moints->addIntegral("ERI_Correlated_Space",eERI_tuwv);
    }
    if(!PERI_tuwv)
    {
      PERI_tuwv = std::make_shared<InCore4indexTPI<MatsT>>(nCorrPO);
      PPTF->transformTPI(pert,PERI_tuwv->pointer(),"tuvw");
      this->moints->addIntegral("PRI_Correlated_Space",PERI_tuwv);
    }
    // Transform the (ee|PP) integrals
    if(!ePERI)
    {
      ePERI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrEO,nCorrPO);
      ePTF->transformAsymmTPI(pert,ePERI->pointer(),eeTF,PPTF,"tuvw","tuvw");
      this->moints->addIntegral("eP_Correlated_Space",ePERI);
    }

    if(nInactE)
    {
      ePTF->transformAsymmTPI(pert,eP_eCore_PActive.pointer(),eeTF,PPTF,"ijkl","tuvw",false);
    }
    if(nInactP)
    {
      ePTF->transformAsymmTPI(pert,eP_PCore_eActive.pointer(),eeTF,PPTF,"tuvw","ijkl",false);
    }


    // Fold the TPI's into the OPI's
#pragma omp parallel for schedule(static) collapse(2) default(shared)
    for(auto u = 0ul; u < enCorrO; u++)
    for(auto t = 0ul; t < enCorrO; t++){
      MatsT tmp = 0.0;
      for(auto v = 0ul; v < enCorrO; v++)
      {
        tmp += 0.5 * eERI_tuwv->operator()(t,v,v,u);
      }
      if(nInactP)
      {
        for(auto I = 0; I < nInactP; I++)
        {
          tmp += eP_PCore_eActive(t,u,I,I);
        }
      }
      ehCoreP_tu(t,u) = ehCore_tu(t,u) - tmp;
    }


    // Fold the TPI's into the OPI's
#pragma omp parallel for schedule(static) collapse(2) default(shared)
    for(auto u = 0ul; u < PnCorrO; u++)
    for(auto t = 0ul; t < PnCorrO; t++){
      MatsT tmp = 0.0;
      for(auto v = 0ul; v < PnCorrO; v++)
      {
        tmp += 0.5 * PERI_tuwv->operator()(t,v,v,u);
      }
      if(nInactE)
      {
        for(auto i = 0; i < nInactE; i++)
        {
          tmp += 2.0*eP_eCore_PActive(i,i,t,u);
        }
      }
      if(nInactP)
      {
        for(size_t i = 0; i < nInactP; i++)
        {
          tmp -= PP_pqrs(i,i,t+nInactP,u+nInactP) - PP_pqrs(i,t+nInactP,i,u+nInactP);
        }
      }
      PhCoreP_tu(t,u) = PhCore_tu(t,u) - tmp;
    }

    // Save the two particle integrals
    this->moints->addIntegral("hCoreP_Correlated_Space",
          std::make_shared<OnePInts<MatsT>>(ehCoreP_tu));
    this->moints->addIntegral("PhCoreP_Correlated_Space",
          std::make_shared<OnePInts<MatsT>>(PhCoreP_tu));

    // Free temporary intermediates
    CQMemManager::get().free(h1e_ii,h1p_ii,GDejj_ii,GDpjj_ii);

    return;
  }

  template<typename MatsT, typename IntsT>
  void NEOMCSCF<MatsT,IntsT>::alloc()
  {
    // Instead of letting MCWavefunction do the allocation, we'll just do it all here
    //MCWaveFunction<MatsT,IntsT>::alloc();

    this->CIVecs = std::vector<MatsT*>(this->NStates);

    // oneRDM is the electronic subsystem
    // using the base MCWaveFunction oneRDM member
    this->oneRDM.reserve(this->NStates); 
    size_t enCorrO = this->ewfn_->MOPartition.nCorrO;

    // PoneRDM is protonoic subsystem
    this->PoneRDM.reserve(this->NStates);
    size_t PnCorrO = this->pwfn_->MOPartition.nCorrO;

    this->ciBuilder = std::make_shared<NEOCASCI<MatsT,IntsT>>();
    NEOCIBuilder = std::static_pointer_cast<NEOCASCI<MatsT,IntsT>>(this->ciBuilder);

    try {
      for (auto i = 0ul; i < this->NStates; i++) {
        this->CIVecs[i] = CQMemManager::get().template malloc<MatsT>(this->NDet);
        this->oneRDM.emplace_back(cqmatrix::Matrix<MatsT>(enCorrO)); 
        this->PoneRDM.emplace_back(cqmatrix::Matrix<MatsT>(PnCorrO)); 
      }
    } catch (...) {
      CErr("Not enough Memory to allocate CIVector for the specified number of determiants");
    }

    // computing list for detstring 
    std::cout << std::endl;
    FormattedLine(std::cout, "Compute Excitation List(s) ...");
    this->ewfn_->detStr->computeList();
    this->ewfn_->detStrBeta->computeList();
    this->pwfn_->detStr->computeList();
    this->pwfn_->detStrBeta->computeList();

    // Read in previous CI vectors, useful for hacking in initial guess CI vectors
    // to RTCASCI Symplectic split hamiltonian
    if(this->readCI)
    {
      //this->ReadGuessCIVector({this->savFile.fName()});
      this->ReadGuessCIVector();
    }

    // End what would normally be called by MCWaveFunction

    // ciSolver is a member of the MCSCF class
    this->ciSolver = std::make_shared<CISolver<MatsT,IntsT>>(this->settings.ciAlg,
      this->settings.maxCIIter, this->settings.ciVectorConv,
      this->settings.maxDavidsonSpace, this->settings.nDavidsonGuess,
      this->settings.energyRefs);
  }

  template<typename MatsT, typename IntsT>
  void NEOMCSCF<MatsT,IntsT>::formNaturalOrbitals()
  {
    this->pwfn_->formNaturalOrbs(this->PoneRDM[this->NatOrbs-1]);
    this->ewfn_->formNaturalOrbs(this->oneRDM[this->NatOrbs-1]);
  }


}; // namespace ChronusQ
