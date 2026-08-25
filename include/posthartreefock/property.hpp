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
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>
#include <physcon.hpp>

#include <configinteraction.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

 /*
  * \brief Perform mulliken analysis for each state
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. SingleSlater->populationAnalysis()
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::populationAnalysis(size_t i) {
    
    ROOT_ONLY(this->comm);

    std::cout << std::endl << "Population Analysis for State " << i+1 << ": ";

    // transform oneRDM to AO basis
    rdm2pdm(*this->oneRDM[i]);

    ref_->populationAnalysis();
    ref_->printMiscProperties(std::cout);

  }; // PostHartreeFock::populationAnalysis

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::populationAnalysis() {
    
    ROOT_ONLY(this->comm);

    for (auto i = 0ul; i < this->NStates; i++) {
      populationAnalysis(i);
    }

  }; // PostHartreeFock::populationAnalysis




 /*
  * \brief Store the AO representation of the PDM in bin
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. save SS->onePDM
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::saveOnePDMs(size_t i) {
    
    ROOT_ONLY(this->comm);

    std::cout << "Saving PDM for State " << i+1 << " as requested." << std::endl;

    // transform oneRDM to AO basis
    rdm2pdm(*this->oneRDM[i]);

    // Convert to AO basis and update PDM in ref
    std::string rdmStr = "POSTHF/RDM-"+std::to_string(i+1);
    savFile.safeWriteData(rdmStr, *ref_->onePDM);

  }; // PostHartreeFock::saveOnePDMs(i)

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::saveOnePDMs() {
   
    for (auto i : this->saveOnePDM_states) {

      PostHartreeFock::saveOnePDMs(i);

    }

  }; // PostHartreeFock::saveOnePDMs()

  /* Active space spin overlap matrices
 
     S_kl = C_k^\dagger S C_l

     stored as

     Saa, Sab, Sba, Sbb

     where 
     k,l are spin components;
     C_k is the k-spin component MO coefficients
     S_kl is the kl-spin overlap matrix
     S is the overlap matrix

  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::populateSpinOverlap() {

    size_t nC = ref_->nC;
    size_t NBasis = ref_->nAlphaOrbital();
    size_t nCorrO = this->corrSpace.nCorrO;
    size_t nInact = this->corrSpace.nInact;

    if (nC == 1 || nC == 2) {
      cqmatrix::Matrix<MatsT> CMOa(NBasis, nCorrO);
      cqmatrix::Matrix<MatsT> CMOb(NBasis, nCorrO);
      cqmatrix::Matrix<MatsT> scratch(nCorrO, NBasis);
      if(nC == 1) {
        for (auto i = 0; i < NBasis; i++) {
          for (auto j = 0; j < nCorrO; j++) {
            CMOa(i,j) = (ref_->moCoefficients[0].get())(i,j+nInact);
            CMOb(i,j) = (ref_->moCoefficients[1].get())(i,j+nInact);
          }
        }
      } else if(nC == 2) {
        for (auto i = 0; i < NBasis; i++) {
          for (auto j = 0; j < nCorrO; j++) {
            CMOa(i,j) = (ref_->moCoefficients[0].get())(i,j+nInact);
            CMOb(i,j) = (ref_->moCoefficients[0].get())(i+NBasis,j+nInact);
          }
        }
      }
      this->spin_overlap.clear();
      this->spin_overlap.reserve(4);
      for (auto i = 0ul; i < 4; i++) {
        this->spin_overlap.emplace_back(std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO)); 
      }
      // Computing aa/ab overlaps
      // Ca^\dagger * S
       blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOa.pointer(), NBasis, 
                 ref_->aoints_->overlap->matrix().pointer(), NBasis,
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);
      // Saa = Ca^\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 this->spin_overlap[0]->pointer(), nCorrO);
      // Sab = Ca^\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 this->spin_overlap[1]->pointer(), nCorrO);
      // Computing bb/ba overlaps
      // Cb^\dagger * S
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOb.pointer(), NBasis, 
                 ref_->aoints_->overlap->matrix().pointer(), NBasis,
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);
      // Sba = Cb\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 this->spin_overlap[2]->pointer(), nCorrO);
      // Sbb = Cb\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 this->spin_overlap[3]->pointer(), nCorrO);
    } else {
      CErr("4-Component spin expectation values NYI - Instance 2");
    }

  };


  /*
  * \brief Spin and angular analysis for each state
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. Compute spin/angular properties with 2RDM and print
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAndAngularAnalysis(size_t i) {
    size_t nCorrO = this->corrSpace.nCorrO;   
    size_t nCoreO = this->corrSpace.nFCore + this->corrSpace.nInact;
    auto tempOneRDM = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDM = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);

    twoRDM->clear(); 
    // make a different 2rdm object later
    
    std::cout << std::endl << "Spin and Angular Analysis for State " << i+1 << ": " << std::endl;

    computeTDM(i, i, tempOneRDM);
    // tempOneRDM->output(std::cout, "1RDM in MO basis for state " + std::to_string(i+1), true);
    rdm2pdm(*tempOneRDM);
    
    twoRDM = computeFull2RDM(i); // Compute 2RDM for the selected state in MO basis
    // twoRDM->output(std::cout, "2RDM in MO basis for state " + std::to_string(i+1), true);
    //print out the active space
    
    const size_t activeEndOff = this->corrSpace.nNegMO + this->corrSpace.nFCore
      + this->corrSpace.nInact + this->corrSpace.nCorrO;
    // std::cout << "Active space end offset: " << activeEndOff << std::endl;

    // Transform the selected state 1RDM back to AO space bc computeSpinAndAngularProperties expects AO basis
    ref_->computeSpinAndAngularProperties(twoRDM.get(), activeEndOff); 
    ref_->printSpin(std::cout, false);
    ref_->printAngularProperties(std::cout, false);

  }; // PostHartreeFock::spinAndAngularAnalysis

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAndAngularAnalysis() {

    for (auto i = 0ul; i < this->NStates; i++) {

      PostHartreeFock::spinAndAngularAnalysis(i);

    }

  }; 

  /*
  * brief: Compute Oscillator Strength in the MO basis
  * Arguments: Ground state s1, target state s2
  * Prints out Osc strength
  *  
  * Formula: 2/3 * (E2 - E1) (sum_pq <psi_0|(e . r)_pq|psi_f>)^2 
  *    
  */ 
  template <typename MatsT, typename IntsT>
  double PostHartreeFock<MatsT,IntsT>::oscillator_strength(size_t s2, size_t s1) {

    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;;
    size_t nCorrO = corrSpace.nCorrO;
    size_t nCoreO = corrSpace.nInact + corrSpace.nFCore;
    
    // compute transition density matrix for specific state
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpAOTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    auto tmpAOTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    
    computeTDM(s1, s2, tmpTDM1);
    rdm2pdm(*tmpTDM1, 1., true);
    if (ref_->nC == 1) *tmpAOTDM1 = (ref_->onePDM->S());
    else if (ref_->nC == 2) ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
    else {
      CErr("NYI");
    }
    
    computeTDM(s2, s1, tmpTDM2);
    rdm2pdm(*tmpTDM2, 1., true);
    if (ref_->nC == 1) *tmpAOTDM2 = (ref_->onePDM->S());
    else if (ref_->nC == 2) ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM2);
    else {
      CErr("NYI");
    }
   
    double f_lenGauge = 0.;
    double f_velGauge = 0.; 

    if (MPIRank(this->comm) == 0){

      const std::array<std::string,3> dipoleList =
        { "X","Y","Z" };

      // Obtain Dipole integrals (Length Gauge: AO basis)
      VectorInts<IntsT> AOdipole(nAO, 1, true);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (ref_->nC == 1)
          AOdipole[iXYZ] = (*ref_->aoints_->lenElectric)[dipoleList[iXYZ]];
        else if (ref_->nC == 2)
          AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->lenElectric)[dipoleList[iXYZ]]
                           ->template spatialToSpinBlock<IntsT>());
      }
      // fED2 ---> Elec Dipole: Calculate sum(<0|D_ab|n>^2) in AO basis:
      MatsT eD_lenGauge = MatsT(0.);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        eD_lenGauge += abs(blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1))
            * abs(blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[iXYZ]->pointer(),1));
      }

      // Replace AOdipole len Gauge integrals with vel gauge:
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (ref_->nC == 1)
          AOdipole[iXYZ] = (*ref_->aoints_->velElectric)[dipoleList[iXYZ]];
        else if (ref_->nC == 2)
          AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->velElectric)[dipoleList[iXYZ]]
                           ->template spatialToSpinBlock<IntsT>());
      }
      // fED2 ---> Elec Dipole: Calculate sum(<0|D_ab|n>^2) in AO basis:
      MatsT eD_velGauge = MatsT(0.);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        eD_velGauge += (-blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[iXYZ]->pointer(),1));
      }


      // oscillator strength f = 2/3 (E2 - E1) eD.
      f_lenGauge = (2./ 3.) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(eD_lenGauge);
      f_velGauge = (2./ 3.) * (1 / (StateEnergy[s2] - StateEnergy[s1])) * std::real(eD_velGauge);

      // output
      std::cout << "Excited State: " << std::setw(3) << std::right << s2+1
                << " to state: " << std::setw(3) << std::right << s1+1 << ":";
      std::cout << std::setw(15) << std::right << "E(Eh) = "
                << std::setprecision(8) << std::fixed << (StateEnergy[s2] - StateEnergy[s1]);
      std::cout << std::setw(15) << std::right << "    f(0)(Len Gauge) = "
                << std::setprecision(12) << std::fixed << f_lenGauge; 
      std::cout << std::setw(15) << std::right << "  ||    f(0)(Vel Gauge) = "
                << std::setprecision(12) << std::fixed << f_velGauge << std::endl;
    }

    MPIBCast(f_lenGauge, 0, this->comm);
    MPIBCast(f_velGauge, 0, this->comm);
    
    return f_lenGauge;

  } // PostHartreeFock::oscillator_strength


  /*
  *        Adding multipolar contributions: Magnetic dipole
  *        and electric quad, octupole moments to oscillator strengths.
  *        (This is in AO Basis: not tested)
  *
  * \brief Compute oscillator strength for MC wavefunction
  *         using AO electric dipole, AO magnetic dipole, 
  *         AO electric quad and MO coefficients and MO TDM.
  *         Only for 1C and 2C
  *         s1: initial state
  *         s2: final state
  */ 

  template <typename MatsT, typename IntsT>
  double PostHartreeFock<MatsT,IntsT>::secondorder_oscillator_strength(size_t s2, size_t s1) {

    if (ref_->nC == 4) CErr("Second-Order Moments NYI for 4CDC-DAS!!");

    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;;
    size_t nCorrO = corrSpace.nCorrO;
    size_t nCoreO = corrSpace.nInact + corrSpace.nFCore;
    double alpha = std::pow(1/SpeedOfLight(), 2);
    double excEn = StateEnergy[s2] - StateEnergy[s1];

    // compute transition density matrix for specific state
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpAOTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    auto tmpAOTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);

    computeTDM(s1, s2, tmpTDM1);
    rdm2pdm(*tmpTDM1, 1., true);
    if(ref_->nC == 1){
      *tmpAOTDM1 = (ref_->onePDM->S());
    }
    else {
      ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
    }

    computeTDM(s2, s1, tmpTDM2);
    rdm2pdm(*tmpTDM2, 1., true);
    if(ref_->nC == 1){
      *tmpAOTDM2 = (ref_->onePDM->S());
    }
    else {
      ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM2);
    }
    
    double f = 0.; 
    if (MPIRank(this->comm) == 0){

      const std::array<std::string,3> dipoleList =
        { "X","Y","Z" };
      const std::array<std::string,6> quadrupoleList =
        { "XX","XY","XZ","YY","YZ","ZZ" };
      const std::array<std::string,10> octupoleList =
        { "XXX","XXY","XXZ","XYY","XYZ","XZZ","YYY",
          "YYZ","YZZ","ZZZ" };
      const std::array<std::string,9> quadrupoleListAsymm =
        { "XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ" };

      // Obtain Dipole integrals (AO basis)
      VectorInts<IntsT> AOdipole(nAO, 1, true);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (ref_->nC == 1)
          AOdipole[iXYZ] = (*ref_->aoints_->velElectric)[dipoleList[iXYZ]];
        else if (ref_->nC == 2)
          AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->velElectric)[dipoleList[iXYZ]]
                           ->template spatialToSpinBlock<IntsT>());
      }
      
      // Obtain Magnetic Dipole intgrals (AO basis)
      VectorInts<IntsT> AO_Magnetic_Dipole(nAO, 1, true);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (ref_->nC == 1)
          AO_Magnetic_Dipole[iXYZ] = (*ref_->aoints_->magnetic)[dipoleList[iXYZ]];
        else if (ref_->nC == 2)
          AO_Magnetic_Dipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->magnetic)[dipoleList[iXYZ]]
                                  ->template spatialToSpinBlock<IntsT>());
        }
       
      // Obtain Magnetic Quadpole intgrals (AO basis)
      VectorInts<IntsT> AO_Magnetic_Quadpole(nAO, 2, false);
      for(auto iXYZ = 0; iXYZ < 9; iXYZ++) {
        if (ref_->nC == 1)
          AO_Magnetic_Quadpole[iXYZ] = (*ref_->aoints_->magnetic)[quadrupoleListAsymm[iXYZ]];
        else if (ref_->nC == 2){
          AO_Magnetic_Quadpole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->magnetic)[quadrupoleListAsymm[iXYZ]]
                                  ->template spatialToSpinBlock<IntsT>());}
        }

      // Obtain electric quad integrals (AO basis) 
      VectorInts<IntsT> AO_Elec_Quadpole(nAO, 2, true);
      for(auto iXYZ = 0; iXYZ < 6; iXYZ++) {
        if (ref_->nC == 1)
           AO_Elec_Quadpole[iXYZ] = (*ref_->aoints_->velElectric)[quadrupoleList[iXYZ]];
        else if (ref_->nC == 2)
           AO_Elec_Quadpole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->velElectric)[quadrupoleList[iXYZ]]
                                ->template spatialToSpinBlock<IntsT>());
      }

      // Obtain electric octupole integrals (AO basis) 
      VectorInts<IntsT> AO_Elec_Octupole(nAO, 3, true);
      for(auto iXYZ = 0; iXYZ < 10; iXYZ++) {
        if (ref_->nC == 1)
           AO_Elec_Octupole[iXYZ] = (*ref_->aoints_->velElectric)[octupoleList[iXYZ]];
        else if (ref_->nC == 2)
           AO_Elec_Octupole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->velElectric)[octupoleList[iXYZ]]
                                ->template spatialToSpinBlock<IntsT>());
      }

      // fED2 ---> Elec Dipole: Calculate sum(<0|D_ab|n>^2) in AO basis:
      MatsT eD = MatsT(0.);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        eD += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[iXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1));
      }

      // fMD2 ---> Magnetic Dipole: Calculate sum(<0|M_ab|n>^2) in AO basis:
      MatsT mD = MatsT(0.);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        mD += (-blas::dotu(nAO*nAO,tmpAOTDM2->pointer(),1,(AO_Magnetic_Dipole)[iXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO,tmpAOTDM1->pointer(),1,(AO_Magnetic_Dipole)[iXYZ]->pointer(),1));
      }
      // fEQ2 ---> Elec Quadpole: Calculate sum(<0|Q_ab|n>^2) - (1/3)(sum<0|Q_aa|n>)^2 in AO basis:
      MatsT eQ = MatsT(0.);
      MatsT tempQ_1 = MatsT(0.);
      MatsT tempQ_2 = MatsT(0.);

      for(auto iXYZ = 0; iXYZ < 6; iXYZ++){
        if(iXYZ == 0 || iXYZ == 3 || iXYZ == 5){
          eQ += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1))
              * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1));
          tempQ_1 += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1));
          tempQ_2 += (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1));
        
        } else{
          eQ += 2. * (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1))
              * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Quadpole[iXYZ]->pointer(),1));
        }
      }
      
      eQ = eQ - ((1./3.) * (tempQ_1 * tempQ_2));
            
      // fED x EO ---> Elec Octupole: Calculate sum(<0|D_b|n><0|eO_aab|n>):
      MatsT eO = MatsT(0.);
      std::map<int, std::vector<int>> indexMap = {
        {0, {0, 3, 5}},
        {1, {1, 6, 8}},
        {2, {2, 7, 9}}
      };
      for(auto eXYZ = 0; eXYZ < 3; eXYZ++){
          eO += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[eXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Octupole[indexMap[eXYZ][0]]->pointer(),1));
          eO += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[eXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Octupole[indexMap[eXYZ][1]]->pointer(),1));
          eO += (-blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[eXYZ]->pointer(),1))
            * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Elec_Octupole[indexMap[eXYZ][2]]->pointer(),1));
      } 

      // fED x MQ ---> Mag Quadpole: Calculate sum(levicivita_aby <0|D_b|n><0|mQ_ya|n>:
      MatsT mQ = MatsT(0.);
      indexMap = {
        {0, {5, 7}},
        {1, {6, 2}},
        {2, {1, 3}}
      };

      for(auto eXYZ = 0; eXYZ < 3; eXYZ++){
        mQ += (blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[eXYZ]->pointer(),1))
             * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1, AO_Magnetic_Quadpole[indexMap[eXYZ][0]]->pointer(),1));
        mQ -= (blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole[eXYZ]->pointer(),1))
             * (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AO_Magnetic_Quadpole[indexMap[eXYZ][1]]->pointer(),1));
      }

      const double f0     =  (2.0  / 3.0)  * (1.0 / excEn) * std::real(eD);
      const double fMD2   =  (1.0  / 6.0)  * alpha * excEn * std::real(mD);
      const double fEQ2   =  (1.0  / 20.0) * alpha * excEn * std::real(eQ);
      const double fEDxMQ =  (1.0  / 9.0)  * alpha * excEn * std::real(mQ);
      const double fEDxEO = -(2.0  / 45.0) * alpha * excEn * std::real(eO);
      const double f2 = f0 + fMD2 + fEQ2 + fEDxMQ + fEDxEO;

      std::cout << std::fixed << std::setprecision(8)
          << "Exc. State: " << (s2 + 1)
          << " | Ground State: " << (s1 + 1)
          << " | E(Eh) = " << excEn
          << std::endl;
      std::cout << BannerTop << std::endl;
      std::cout << std::fixed << std::setprecision(10)
          << "  f(0) = "   << f0
          << " |     f(2) = "   << f2
          << " |    f(EQ) = " << fEQ2
          << std::endl;

      std::cout << std::fixed << std::setprecision(10)
          << "f(MD2) = "   << fMD2
          << " | f(EDxMQ) = " << fEDxMQ
          << " | f(EDxEO) = " << fEDxEO
          << std::endl;
      std::cout << BannerTop << std::endl;

    }                                                                     
    MPIBCast(f, 0, this->comm);                                           
                                                                          
    return f;
  } //PostHartreeFock:Full 2nd-Order Osc Strength


  /*
   *  Precompute 4C AO Dipole Integral when required:
   *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::compute4CAODipole() {
  
    if (!AODipole4C_)
      AODipole4C_ = ref_->aoints_->lenElectric->gather4CDipole();

  }

  /*
  *
  * brief: Compute Oscillator Strength in the MO basis
  * Arguments: Ground state s1, target state s2
  * Prints out Osc strength
  *  
  * Formula: 2/3 * (E2 - E1) (sum_pq <psi_0|(e . r)_pq|psi_f>)^2 
  *    
  */ 
  template <typename MatsT, typename IntsT>
  double PostHartreeFock<MatsT,IntsT>::oscillator_strength4C(size_t s2, size_t s1) {

    if constexpr (std::is_same_v<MatsT, double>) {
      CErr("Incorrect function!! Should be 4-Component CI Calculation");
      return 0.;
    }

    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;;
    size_t nCorrO = corrSpace.nCorrO;
    
    // compute transition density matrix for specific state
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpAOTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    auto tmpAOTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);

    computeTDM(s1, s2, tmpTDM1);
    rdm2pdm(*tmpTDM1, 1., true);
    ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
    computeTDM(s2, s1, tmpTDM2);
    rdm2pdm(*tmpTDM2, 1., true);
    ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM2);

    double f = 0.;

    if (MPIRank(this->comm) == 0){

      const std::array<std::string,3> dipoleList =
        { "X","Y","Z" };
      this->compute4CAODipole();

      // fED2 ---> Elec Dipole: Calculate sum(<0|D_ab|n>^2) in AO basis:
      std::complex<double> eD(0., 0.);
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        
        cqmatrix::Matrix<MatsT> AOdipole = (*AODipole4C_)[iXYZ]
                                          .template spinGather<MatsT>();

        eD += (blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole.pointer(),1)) *
              (blas::dotu(nAO*nAO, tmpAOTDM2->pointer(),1,AOdipole.pointer(),1));
      }

      // oscillator strength f = 2/3 (E2 - E1) eD.
      f = (2./ 3.) * (StateEnergy[s2] - StateEnergy[s1]) * std::abs(eD);

      // output
      std::cout << "Excited State: " << std::setw(3) << std::right << s2+1
                << " to state: " << std::setw(3) << std::right << s1+1 << ":";
      std::cout << std::setw(15) << std::right << "E(Eh) = "
                << std::setprecision(8) << std::fixed << (StateEnergy[s2] - StateEnergy[s1]);
      std::cout << std::setw(15) << std::right << "  f(0) = "
                << std::setprecision(12) << std::fixed << f << std::endl;
    }

    MPIBCast(f, 0, this->comm);
    
    return f;
  } //PostHartreeFock:Full 2nd-Order Osc Strength

  template <typename MatsT, typename IntsT>
  std::vector<MatsT> PostHartreeFock<MatsT,IntsT>::computeStateSpecificDipoleMom(size_t i) {
  
    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;
    size_t nCorrO = corrSpace.nCorrO;
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpAOTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    std::vector<MatsT> SSDipole(3);  
    
    if (MPIRank(this->comm) == 0) {
      computeTDM(i, i, tmpTDM1);
      rdm2pdm(*tmpTDM1, 1., true);
      const std::array<std::string,3> dipoleList =
        { "X","Y","Z" };
      
      if (ref_->nC == 1) { 
        *tmpAOTDM1 = (ref_->onePDM->S());
        VectorInts<IntsT> AOdipole(nAO, 1, true);
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          AOdipole[iXYZ] = (*ref_->aoints_->lenElectric)[dipoleList[iXYZ]];
          SSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1);
        }
      }

      else if (ref_->nC == 2 ) { 
        ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
        VectorInts<IntsT> AOdipole(nAO, 1, true);
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->lenElectric)[dipoleList[iXYZ]]
                           ->template spatialToSpinBlock<IntsT>());
          SSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1);
        }
      }

      else if (ref_->nC == 4 ) {
        if constexpr (std::is_same_v<MatsT, dcomplex>) {
          ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
          for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
            auto AOdipole = (*AODipole4C_)[iXYZ].template spinGather<MatsT>();
            SSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(), 1, AOdipole.pointer(), 1);
          }
        } 
        else {
          CErr("4C requires complex MatsT");
          return {};
        }
      }
      
      else {
        CErr("Wrong Number of Components!!");
        return {};
      }

    }
    MPIBCast(SSDipole.data(), 3, 0, this->comm);
    return SSDipole;

  } // PostHartreeFock<MatsT,IntsT>::computeStateSpecificDipoleMom(size_t i)


  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::printAllStateSpecificDipoleMom() {

    size_t NStates = this->NStates;
    if (MPIRank(this->comm) == 0) {
      
      std::cout << "\n State Specific Dipole Moments (Length Gauge): \n" << std::endl;
      auto tdmTensor = std::make_shared<cqmatrix::Matrix<MatsT>>(3, NStates);
      for (size_t i = 0ul; i < NStates; ++i) {
      
      // fED2 ---> Elec Dipole: Calculate sum(<0|D_ab|n>^2) in AO basis:
      auto DipMoment = computeStateSpecificDipoleMom(i);
      for (size_t iXYZ = 0; iXYZ < 3; ++iXYZ) {
        (*tdmTensor)(iXYZ, i) = DipMoment[iXYZ]; }
      std::cout << "State: " << i << std::endl;
      std::cout << "ED(x) | ED(y) | ED(z)   \n" << DipMoment[0] << " | " << 
        DipMoment[1] << " | " << DipMoment[2] <<  "\n" << BannerTop << std::endl;
      } 

      if (savFile.exists()) {
        savFile.safeWriteData("POSTHF/STATESPECIFIC_DIPOLEMOMENTS", tdmTensor->pointer(), {NStates, 3});
      }

    }

  } //PostHartreeFock<MatsT,IntsT>::printAllStateSpecificDipoleMom()



  template <typename MatsT, typename IntsT>
  std::vector<MatsT> PostHartreeFock<MatsT,IntsT>::computeTransitionDipoleMom(size_t i, size_t j) {

    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;
    size_t nCorrO = corrSpace.nCorrO;
    
    // compute transition density matrix for specific state
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpAOTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
    std::vector<MatsT> TSDipole(3);  
    if (MPIRank(this->comm) == 0) {
    
      computeTDM(i, j, tmpTDM1);
      rdm2pdm(*tmpTDM1, 1., true);
      const std::array<std::string,3> dipoleList =
        { "X","Y","Z" };
      
      if (ref_->nC == 1) { 
        *tmpAOTDM1 = (ref_->onePDM->S());
        VectorInts<IntsT> AOdipole(nAO, 1, true);
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          AOdipole[iXYZ] = (*ref_->aoints_->lenElectric)[dipoleList[iXYZ]];
          TSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1);
        }
      }

      else if (ref_->nC == 2 ) { 
        ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
        VectorInts<IntsT> AOdipole(nAO, 1, true);
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->lenElectric)[dipoleList[iXYZ]]
                           ->template spatialToSpinBlock<IntsT>());
          TSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(),1,AOdipole[iXYZ]->pointer(),1);
        }
      }

      else if (ref_->nC == 4 ) {
        if constexpr (std::is_same_v<MatsT, dcomplex>) {
          ref_->onePDM-> template spinGather<MatsT>(*tmpAOTDM1);
          for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
            auto AOdipole = (*AODipole4C_)[iXYZ].template spinGather<MatsT>();
            TSDipole[iXYZ] = blas::dotu(nAO*nAO, tmpAOTDM1->pointer(), 1, AOdipole.pointer(), 1);
          }
        } 
        else {
          CErr("4C requires complex MatsT");
          return {};
        }
      }
      
      else {
        CErr("Wrong Number of Components!!");
        return {};
      }
    }
    MPIBCast(TSDipole.data(), 3, 0, this->comm);
    return TSDipole;

  } // PostHartreeFock<MatsT,IntsT>::computeTransitionDipoleMom


  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::printAllTransitionDipoleMom() {

    size_t NStates = this->NStates;
    size_t n = (NStates * (NStates - 1)) / 2;
    if (n == 0) { 
      CErr("Only one electronic state from CI! Please check calculation."); 
    }
   
    if (MPIRank(this->comm) == 0) {
      std::cout << "\n Transition Dipole Moments between states (Length Gauge) \n" << std::endl; 
      auto tdmTensor = std::make_shared<cqmatrix::Matrix<MatsT>>(3, n);
      for (size_t i = 0ul; i < NStates; ++i) {
        for (size_t j = 0ul; j < i; ++j) {
          
          size_t compIdx = (( i * (i - 1)) / 2) + j;
          auto tDipMoment = computeTransitionDipoleMom(i, j);
          for (size_t iXYZ = 0; iXYZ < 3; ++iXYZ) {
            (*tdmTensor)(iXYZ, compIdx) = tDipMoment[iXYZ]; }
          
          std::cout << "States: " << i << " | " << j << std::endl;
          std::cout << std::fixed << std::setprecision(14);
          std::cout << "Excitation Energy (in au): " << StateEnergy[i] - StateEnergy[j] << std::endl; 
          std::cout << "ED(x) | ED(y) | ED(z) \n " << tDipMoment[0] << " | " << 
            tDipMoment[1] << " | " << tDipMoment[2] << "\n" << BannerTop << std::endl;
        }
      }

      if (savFile.exists()) {
        savFile.safeWriteData("POSTHF/TRANSITION_DIPOLEMOMENTS", tdmTensor->pointer(), {n, 3});
      }

    }    

  } // PostHartreeFock<MatsT,IntsT>::printAllTransitionDipoleMom


/*
 * Adding Difference RDM analysis 
 *
 */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::OneRDMDiff() {

    if (this->NStates == 1) CErr("NStates should be greater than one");

    size_t nCorrO = corrSpace.nCorrO;
    auto RDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpRDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpRDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto diagRDMDiff = std::make_shared<cqmatrix::Matrix<MatsT>>(
       nCorrO,  this->NStates-this->NosS1);

    computeTDM(0,0,RDM1);
    blas::scal(nCorrO*nCorrO, MatsT(1./this->NosS1), RDM1->pointer(),1);
    for (auto i = 1; i < this->NosS1; i++) {
      computeTDM(i,i,tmpRDM1);
      blas::scal(nCorrO*nCorrO, MatsT(1./this->NosS1), tmpRDM1->pointer(), 1);
      blas::axpy(nCorrO*nCorrO, MatsT(1.), tmpRDM1->pointer(), 1, RDM1->pointer(), 1);
    }
    
    for (size_t i = 0; i < this->NStates-this->NosS1; i++) {
      computeTDM(this->NosS1+i, this->NosS1+i, tmpRDM2);
      std::cout << std::fixed << std::setprecision(2); 
      std::cout << "\n\nDiff RDM of State: " << (this->NosS1+i+1) << std::endl;
      
      for (size_t j = 0; j < nCorrO; j++) {
        (*diagRDMDiff)(j,i) = (*tmpRDM2)(j,j) - (*RDM1)(j,j);
        std::cout << " " << std::real((*diagRDMDiff)(j,i)) << " ";
      }
      
    }

    if (savFile.exists()) 
      savFile.safeWriteData("POSTHF/RDMDIFF", diagRDMDiff->pointer(), 
          {this->NStates-1, nCorrO});

    RDM1.reset();
    tmpRDM1.reset();
    tmpRDM2.reset();
    diagRDMDiff.reset();

  }
                       

}; // namespace ChronusQ
