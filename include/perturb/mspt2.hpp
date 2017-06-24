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

#include <perturb.hpp>
#include <cqlinalg.hpp>

#include <itersolver.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

//#define _DEBUG_PT2_impl

namespace ChronusQ {

  // To help change type dcomplex to MatsT
  template <typename MatsT>
  MatsT toMatsT(dcomplex a) {
    return a;
  }
  template <> double toMatsT(dcomplex a) {
    return std::real(a);
  }
  template dcomplex toMatsT(dcomplex);



  /**
   * 
   *  \brief rotate the CI vectors by diagonaling zero-order H_eff
   *         Only used in Extended MS
   *         ! ref is still using CAS settings
   *         Unitary transformation: U^ f U = d
   *         f = <M|f|N>
   *         the rotated CI vectors are saved in ref.CIVecs
   *
   */  
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::rotateCIVecs() {

    // build state averaged Fock
    refMCwfn->computeOneRDM(); // TODO: now computing for all MCSCF states...
#ifdef _DEBUG_PT2_impl
    refMCwfn->print1RDMs();
#endif
    std::cout<<"computed 1RDM."<<std::endl;

    auto & oneRDM = refMCwfn->oneRDMSOI;
    std::cout<<"!! XMS -- CAUTIOUS !!!"<<std::endl;
    std::cout<<" More tests needed. -- Lixin" <<std::endl;
    std::cout<<"XMS: build state-averaged Fock"<<std::endl;
    buildFock(*ptFock_, *oneRDM);

    // rotate the CI Vectors
    auto & ptFock = *ptFock_;
    auto & mopart = this->MOPartition;
    size_t nStates = this->NStates;
    auto & nActO = mopart.nActOs;

    // build zero-order H_eff
    cqmatrix::Matrix<MatsT> MfN(nStates);
    cqmatrix::Matrix<MatsT> TDMscr(mopart.nCorrO);
    MfN.clear();
    for (auto m = 0ul; m < nStates; m++)
    for (auto n = 0ul; n < nStates; n++) {
      if (m == n) // same state: <m|f|m> = zero-order E for m state
        MfN(m, n) = computeZeroE(ptFock, refMCwfn->oneRDM[SoI_[m]]);
      else { // <m|f|n> = \sum_uv ptFock_uv * TDM_uv
        TDMscr.clear();
        refMCwfn->ciBuilder->computeTDM(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*refMCwfn),
                        refMCwfn->CIVecs[SoI_[m]], refMCwfn->CIVecs[SoI_[n]], TDMscr);

        for (auto u = nActO[0]; u < nActO[0]+nActO[1]; u++)
        for (auto v = nActO[0]; v < nActO[0]+nActO[1]; v++)
          MfN(m, n) += ptFock(u, v) * TDMscr(u - nActO[0], v - nActO[0]);
      }
      std::cout<<"MfN m "<<m<<", n "<<n<<": "<< std::setprecision(10)<<MfN(m, n)<<std::endl;
    }
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"PT2 full H_eff", MfN.pointer(), nStates, nStates, nStates);
#endif

    // diagonalize zero-order H_eff
    //MatsT * U = CQMemManager::get().malloc<MatsT>(nStates*nStates);
    dcomplex * Heff_diag = CQMemManager::get().malloc<dcomplex>(nStates); 
    MatsT * dummy = nullptr;
    // TODO: GeneralEigen VS HermetianEigen?
    //GeneralEigen('N','V',nStates,MfN.pointer(),nStates,Heff_diag,dummy,1,U,nStates);
    HermetianEigen('V','L',nStates,MfN.pointer(),nStates,Heff_diag);

#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"PT2 H_eff_0 eigenvectors", U, nStates, nStates, nStates);
    prettyPrintSmart(std::cout,"PT2 H_eff_0 diagonal", Heff_diag, nStates, 1, nStates);
    //testing
    MatsT * UU = CQMemManager::get().malloc<MatsT>(nStates*nStates);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nStates,nStates,
               nStates,MatsT(1.),U,nStates,U,nStates,MatsT(0.),UU,nStates);
    prettyPrintSmart(std::cout,"PT2 U^U",UU, nStates, nStates, nStates);
    CQMemManager::get().free(UU);
#endif

    MatsT * CIVref = CQMemManager::get().malloc<MatsT>(nStates*CIsize);
    MatsT * CIVpt = CQMemManager::get().malloc<MatsT>(nStates*CIsize);
    std::fill_n(CIVref, nStates*CIsize, MatsT(0.));
    // save diagonalized H_eff into E0_
    for (auto i = 0ul; i < nStates; i++) {
      E0_[i] = toMatsT<MatsT>(Heff_diag[i]); // not including ptFock frozen orbital part
      std::copy_n(refMCwfn->CIVecs[SoI_[i]], CIsize, CIVref+i*CIsize);
//      std::cout<<"E0 state "<<i<<": "<< std::setprecision(12)<<E0_[i]<<std::endl;
    }
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"PT2 CAS CIVref", CIVref, CIsize, nStates, CIsize);
#endif

    // rotate CI vectors
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,CIsize,nStates,
               nStates,MatsT(1.),CIVref,CIsize,MfN.pointer(),nStates,MatsT(0.),CIVpt,CIsize);
//               nStates,MatsT(1.),CIVref,CIsize,U,nStates,MatsT(0.),CIVpt,CIsize);
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"PT2 rotated CIV", CIVpt, CIsize, nStates, CIsize);
#endif

    for (auto i = 0ul; i < nStates; i++) {
      std::copy_n(CIVpt+i*CIsize, CIsize,  refMCwfn->CIVecs[SoI_[i]]); // update rotated CI vectors
#ifdef _DEBUG_PT2_impl
      std::cout<<"state "<<i<<std::endl;
      prettyPrintSmart(std::cout,"PT2 CI Vecs after rotation", refMCwfn->CIVecs[i], CIsize,1,CIsize);
#endif
    }
    //CQMemManager::get().free(U, Heff_diag,CIVref,CIVpt);
    CQMemManager::get().free(Heff_diag,CIVref,CIVpt);

  } // PERTURB::rotateCIVecs

  /**
   *  
   *  \brief build H_eff for multi-state PT2
   *         when H_eff is in primary space (nStates X nStates)
   *         H_eff(M,N) = <M|H|N> + 1/2 (<M|HX|N> + <M|X^ H|N>)
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildHeff(cqmatrix::Matrix<MatsT> & H) {

    auto PT2HeffSt = tick();
//    std::cout<<"H_eff build starts:"<<std::endl;
    size_t nStates = this->NStates;

    H.clear();    
    std::vector<MatsT*> HV(nStates);
    for (auto i = 0ul; i < nStates; i++) {
      HV[i] = CQMemManager::get().malloc<MatsT>(this->NDet);
      computeHV(HV[i], this->CIVecs[i]);
    }
    for (auto m = 0ul; m < nStates; m++)
    for (auto n = 0ul; n < nStates; n++) {
      if (m == n) {
        H(m, m) += refMCwfn->StateEnergy[SoI_[m]] + 0.5 * (computeCHV(this->CIVecs[m], HV[m])
                + SmartConj(computeCHV(this->CIVecs[m], HV[m])))
                - std::real(computeShiftCorrection(m));
      } else {
        H(m, n) += 0.5 * (computeCHV(this->CIVecs[m], HV[n])
                + SmartConj(computeCHV(this->CIVecs[n], HV[m])));
      }  
    }
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"LL PT2 H_eff --", H.pointer(), nStates, nStates, nStates);
#endif

    double PT2Heffdur = tock(PT2HeffSt);
    std::cout << "\nPT2 Heff build - DURATION = " << std::setprecision(8)
                << PT2Heffdur << " s." << std::endl;

  } // PERTURB::buildHeff

  /**
   * 
   *  \brief full diagonalization of H_eff to get the perturbed state energies
   *  
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::diagHeff(cqmatrix::Matrix<MatsT> & H) {

    ProgramTimer::tick("Heff diag");
//    std::cout<<"H_eff diagonalization starts:"<<std::endl;
    size_t nStates = this->NStates;

    // GeneralEigen and HermetianEigen give the same results, but not the same eigvec.
//    MatsT * eigvec = CQMemManager::get().malloc<MatsT>(nStates*nStates);
    double * energy = CQMemManager::get().malloc<double>(nStates);
//    dcomplex * energy = CQMemManager::get().malloc<dcomplex>(nStates);
//    MatsT * dummy = nullptr;
    HermetianEigen('V','L',nStates,H.pointer(),nStates,energy);
//    GeneralEigen('N','V',nStates,H.pointer(),nStates,energy,dummy,1,eigvec,nStates);
#ifdef _DEBUG_PT2_impl
    MatsT * test = CQMemManager::get().malloc<MatsT>(nStates*nStates);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nStates,nStates,
               nStates,MatsT(1.),H.pointer(),nStates,H.pointer(),nStates,MatsT(0.),test,nStates);
    prettyPrintSmart(std::cout,"test H_eff diag vector", test,nStates,nStates,nStates);
    CQMemManager::get().free(test);
#endif

    for (auto i = 0ul; i < nStates; i++) {
      this->StateEnergy[i] = std::real(energy[i]);
//      std::cout<<"state: "<<i<<" energy raw: "<< std::setprecision(12)<<energy[i]<<std::endl;
    }
    CQMemManager::get().free(energy);
//    CQMemManager::get().free(eigvec, energy);    
    ProgramTimer::tock("Heff diag");
  } //PERTURB::diagHeff


} // namespace ChronusQ


