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
  std::vector<cqmatrix::Matrix<MatsT>> PostHartreeFock<MatsT,IntsT>::spinOverlap() {

    size_t nC = ref_->nC;
    size_t NBasis = ref_->nAlphaOrbital();
    size_t nCorrO = this->corrSpace.nCorrO;
    size_t nInact = this->corrSpace.nInact;

    std::vector<cqmatrix::Matrix<MatsT>> spin_overlap;

    cqmatrix::Matrix<MatsT> S = ref_->aoints_->overlap->matrix();
    cqmatrix::Matrix<MatsT> CMOa(NBasis, nCorrO);
    cqmatrix::Matrix<MatsT> CMOb(NBasis, nCorrO);
    cqmatrix::Matrix<MatsT> scratch(nCorrO, NBasis);
   
    if(nC == 1) {
      SetMat('N', NBasis, NBasis, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nInact, NBasis, CMOa.pointer(), NBasis);
      SetMat('N', NBasis, NBasis, MatsT(1.0), ref_->moCoefficients[1].get().pointer()+NBasis*nInact, NBasis, CMOb.pointer(), NBasis);
    }

    else if(nC == 2) {
      SetMat('N', NBasis*nC, NBasis*nC, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nC*nInact, NBasis*nC, CMOa.pointer(), NBasis);
      SetMat('N', NBasis*nC, NBasis*nC, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nC*nInact+NBasis, NBasis*nC, CMOb.pointer(), NBasis);
    }

    else {
      CErr("4-Component spin expectation values NYI - Instance 1");
    }

    if (nC == 1 || nC == 2) {
      spin_overlap = std::vector<cqmatrix::Matrix<MatsT>>(4, cqmatrix::Matrix<MatsT>(nCorrO));

      // Computing aa/ab overlaps
      // Ca^\dagger * S
       blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOa.pointer(), NBasis, 
                 S.pointer(), NBasis, 
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);

      // Saa = Ca^\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[0].pointer(), nCorrO);

      // Sab = Ca^\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[1].pointer(), nCorrO);

      // Computing bb/ba overlaps
      // Cb^\dagger * S
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOb.pointer(), NBasis, 
                 S.pointer(), NBasis, 
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);

      // Sba = Cb\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[2].pointer(), nCorrO);

      // Sbb = Cb\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[3].pointer(), nCorrO);
    }

    else {
      CErr("4-Component spin expectation values NYI - Instance 2");
    }

    return spin_overlap;

  };



  /*
  * \brief Spin analysis for each state
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. SingleSlater->populationAnalysis()
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAnalysis(size_t i, std::vector<cqmatrix::Matrix<MatsT>>* spin_overlap) {

    MatsT ssq;
    MatsT sz = MatsT(0.0);
    MatsT ssq0 = MatsT(0.0);
    MatsT ssqx = MatsT(0.0);
    MatsT ssqy = MatsT(0.0);
    MatsT ssqz = MatsT(0.0);
    size_t nCorrO = this->corrSpace.nCorrO;   
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);

    twoRDMSOI->clear(); 
    
    std::cout << std::endl << "Spin Analysis for State " << i+1 << ": " << std::endl;

    ConfigurationInteraction<MatsT,IntsT>* ci = dynamic_cast<ConfigurationInteraction<MatsT,IntsT>*>(this);
    if (ci == nullptr) {
      CErr("Dyanmic cast from PostHartreeFock to ConfigurationInteraction failed.");
    }

    ci->compute2TDM(i, i, twoRDMSOI);

    #pragma omp declare reduction (+ : std::complex<double> : omp_out += omp_in) initializer(omp_priv = std::complex<double>(0, 0))
    #pragma omp parallel for collapse(4) schedule(static) default(shared) reduction(+:ssq0, sz, ssqx, ssqy, ssqz)
    for (auto p = 0ul; p < nCorrO; p++) {
      for (auto q = 0ul; q < nCorrO; q++) {
        for (auto r = 0ul; r < nCorrO; r++) {
          for (auto s = 0ul; s < nCorrO; s++) {
          
            MatsT actual2RDM;

            if (r == 0 and s == 0) {
              // sum_k <S_k^2> = (3/4) sum_pq ( gamma_pq (Saa+Sbb)_pq )
              ssq0 += (*this->oneRDM[i])(p,q) * (((*spin_overlap)[0])(p,q) + ((*spin_overlap)[3])(p,q));
   
              // <Sz> = (1/2) sum_pq ( gamma_pq (Saa-Sbb)_pq )
              sz += (*this->oneRDM[i])(p,q) * (((*spin_overlap)[0])(p,q) - ((*spin_overlap)[3])(p,q));
            }

            if (q == r) {
              actual2RDM = (*twoRDMSOI)(p,q,r,s) - (*this->oneRDM[i])(p,s);
            }

            else {
              actual2RDM = (*twoRDMSOI)(p,q,r,s);
            }

            // <S_k(1)S_k(2)> = sgn(k) (1/4) sum_pqrs ( Gamma_pqrs S_pq S_rs )
            ssqx += (((*spin_overlap)[1])(p,q) + ((*spin_overlap)[2])(p,q)) * actual2RDM * (((*spin_overlap)[1])(r,s) + ((*spin_overlap)[2])(r,s));
            ssqy += (((*spin_overlap)[1])(p,q) - ((*spin_overlap)[2])(p,q)) * actual2RDM * (((*spin_overlap)[1])(r,s) - ((*spin_overlap)[2])(r,s));
            ssqz += (((*spin_overlap)[0])(p,q) - ((*spin_overlap)[3])(p,q)) * actual2RDM * (((*spin_overlap)[0])(r,s) - ((*spin_overlap)[3])(r,s));

          }
        }
      }
    }

    ssq = MatsT(0.25)*ssqx + MatsT(0.25)*ssqz - MatsT(0.25)*ssqy + MatsT(0.75)*ssq0;

    std::cout << "<Sz> : " << MatsT(0.5)*sz << std::endl;
    std::cout << "<S^2> : " << ssq << std::endl << std::endl; 

    // for debugging
    #if 0 
    std::cout << "State : " << i+1 << "   Sx^2 : " << MatsT(0.25)*ssqx << std::endl;
    std::cout << "State : " << i+1 << "   Sy^2 : " << MatsT(-0.25)*ssqy << std::endl;
    std::cout << "State : " << i+1 << "   Sz^2 : " << MatsT(0.25)*ssqz << std::endl;
    std::cout << "State : " << i+1 << "   S0^2 : " << MatsT(0.75)*ssq0 << std::endl;
    #endif



  }; // PostHartreeFock::spinAnalysis

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAnalysis() {

    spin_overlap = spinOverlap();

    for (auto i = 0ul; i < this->NStates; i++) {

      PostHartreeFock::spinAnalysis(i, &spin_overlap);

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

      // oscillator strength f = 2/3 (E2 - E1) eD.
      f =  (2./ 3.) * (1 / (StateEnergy[s2] - StateEnergy[s1])) * std::real(eD);
      // Magnetic dipole contri. to oscillator strength f = 1/6 * alpha^2 * (E2 - E1) * eM.  
      f += (1./6.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(mD);
      // Electric quadpole contri. to oscillator strength f = 1/20 * alpha^2 * (E2 - E1) * eQ.
      f += (1./20.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(eQ);
      // Magnetic quadpole contri. to oscillator strength f = 1/9 * alpha^2 * (E2 - E1) * mQ.
      f += (1./9.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(mQ);
      // Electric Octupole contri. to oscillator strength f = (-2/45) * alpha^2 * (E2 - E1) * eO
      f -= (2./45.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(eO);

      // output
      std::cout << "Excited State: " << std::setw(3) << std::right << s2+1
                << " to state: " << std::setw(3) << std::right << s1+1 << ":";
      std::cout << std::setw(15) << std::right << "E(Eh) = "
                << std::setprecision(8) << std::fixed << (StateEnergy[s2] - StateEnergy[s1]);
      std::cout << std::setw(30) << std::right << "f(2) = "
                << std::setprecision(15) << std::fixed << f << std::endl;

      //Individual Contribution:
      std::cout << bannerTop << std::endl; 
      std::cout << std::setw(12) << "f(0)" << std::setw(8) << "||"
                << std::setw(10) << "f(MD2)"   << std::setw(8) << "||"
                << std::setw(10) << "f(EQ2)"   << std::setw(8) << "||"
                << std::setw(10) << "f(EDxMQ)" << std::setw(8) << "||"
                << std::setw(10) << "f(EDxEO)" << std::setw(8) << "||"
                << std::endl;
     
      std::cout << bannerTop << std::endl; 
      std::cout << std::fixed << std::setprecision(10)
                << std::setw(12) << std::fixed << std::setprecision(8) << (2./ 3.) * (1/ (StateEnergy[s2] - StateEnergy[s1])) * std::real(eD)
                << std::setw(20) << std::fixed << std::setprecision(10) << (1./6.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(mD)
                << std::setw(20) << std::fixed << std::setprecision(10) << (1./20.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(eQ)
                << std::setw(20) << std::fixed << std::setprecision(10) << (1./9.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(mQ)
                << std::setw(20) << std::fixed << std::setprecision(10) << (-2./45.) * std::pow(1/SpeedOfLight(), 2) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(eO)
                << "\n" << std::endl;                                    
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

  } // PostHartreeFock::4Coscillator_strength
                       
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
