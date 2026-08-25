/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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

#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>
//#include <coupledcluster/EOMCC.hpp>
#include <util/math.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <functional>
#include <util/timer.hpp>
#include <coupledcluster/MBExpansion.hpp>
#include <itersolver/davidson.hpp>
#include <itersolver.hpp>

namespace ChronusQ{

  template <typename MatsT>
  EOMRCCSD<MatsT>::EOMRCCSD(const SafeFile &savFile,
                                CCIntermediates<MatsT> &intermediatesRref,
                                const EOMSettings &eomSettings,
                                const CoupledClusterSettings &ccSettings):
      EOMCCBase<MatsT>(savFile, intermediatesRref, eomSettings,ccSettings),
      fockMatrix_ta_RHF(intermediatesRref.fockMatrix),
      Moints(intermediatesRref.moInts),
      muMatrix_RHF(intermediatesRref.muMatrix),
      vLabel_RHF(intermediatesRref.vLabel), oLabel_RHF(intermediatesRref.oLabel),
      T1_RHF(intermediatesRref.T->get_tensor("OneBody")), T2_RHF(intermediatesRref.T->get_tensor("TwoBody")),
      tau_RHF(intermediatesRref.tau),
      F_ae_RHF(intermediatesRref.F_ae),
      F_mi_RHF(intermediatesRref.F_mi),
      F_me_RHF(intermediatesRref.F_me),
      W_mnij_RHF(intermediatesRref.W_mnij),
      W_abef_RHF(intermediatesRref.W_abef),
      W_mbej_RHF_baab(intermediatesRref.W_mbej_baab),
      W_mbej_RHF_baba(intermediatesRref.W_mbej_baba),
      W_mnie_RHF(intermediatesRref.W_mnie),
      W_amef_RHF(intermediatesRref.W_amef),
      W_mbij_RHF(intermediatesRref.W_mbij),
      W_abei_RHF(intermediatesRref.W_abei),
      G_ae_RHF(intermediatesRref.G_ae),
      G_mi_RHF(intermediatesRref.G_mi),
      D_ai_RHF(intermediatesRref.D_ai),
      D_abij_RHF(intermediatesRref.D_abij) {
    
    TAManager &TAmanager = TAManager::get();

    // without L, we don't need D
    if (not (eomSettings.oscillator_strength or ccSettings.crcc or ccSettings.computeDipole)) {
      if(intermediatesRref.D_ai)   TAmanager.free("vo", std::move(intermediatesRref.D_ai), true);
      if(intermediatesRref.D_abij) TAmanager.free("vvoo", std::move(intermediatesRref.D_abij), true);
    }

    nV_RHF = TAmanager.getRange(vLabel_RHF).extent();
    nO_RHF = TAmanager.getRange(oLabel_RHF).extent();
    nOVshift_RHF = nV_RHF * nO_RHF;

    this->Hbar_dim = nOVshift_RHF + nOVshift_RHF*(nOVshift_RHF+1)/2;

    this->tensor_builder_.push_back(std::string({intermediatesRref.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediatesRref.oLabel}));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->tensor_builder_.push_back(std::string({intermediatesRref.vLabel,intermediatesRref.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediatesRref.oLabel,intermediatesRref.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));
  }




  template <typename MatsT>
  void EOMRCCSD<MatsT>::initializeEOMCC() {
    TAManager &TAmanager = TAManager::get();

    if (not W_mnie_RHF.is_initialized()){
      W_mnie_RHF = TAmanager.malloc<MatsT>("ooov");
    }

    if (not W_amef_RHF.is_initialized()){
      W_amef_RHF = TAmanager.malloc<MatsT>("vovv");
    }

    if (not W_mbij_RHF.is_initialized()){
      W_mbij_RHF = TAmanager.malloc<MatsT>("ovoo");
    }

    if (not W_abei_RHF.is_initialized()){
      W_abei_RHF = TAmanager.malloc<MatsT>("vvvo");
    }

  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formF_ae() {
    F_ae_RHF("a,e") -= 0.5 * this->T1_RHF("a,m") * F_me_RHF("m,e");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formF_mi() {
    F_mi_RHF("m,i") += 0.5 * this->T1_RHF("e,i") * F_me_RHF("m,e");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_mnij() {
    W_mnij_RHF("m,n,i,j") += 0.5 * tau_RHF("e,f,i,j") * conj(Moints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_abef() {
    W_abef_RHF("a,b,e,f") += 0.5 * tau_RHF("a,b,m,n") * conj(Moints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_mbej() {
    W_mbej_RHF_baab("m,b,e,j") += 0.5 * this->T2_RHF("f,b,j,n") * conj(Moints["vvoo"]("f,e,m,n"));
    W_mbej_RHF_baba("m,b,e,j") -= 0.5 * this->T2_RHF("f,b,j,n") * conj(Moints["vvoo"]("e,f,m,n"));
    W_mbej_RHF_baba("m,b,e,j") += this->T2_RHF("b,f,j,n") * conj(Moints["vvoo"]("e,f,m,n"));
    W_mbej_RHF_baba("m,b,e,j") -= 0.5 * this->T2_RHF("b,f,j,n") * conj(Moints["vvoo"]("f,e,m,n"));
  } 

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_mnie() {
    W_mnie_RHF("m,n,i,e") = conj(Moints["vooo"]("e,i,n,m")) + this->T1_RHF("f,i") * conj(Moints["vvoo"]("f,e,m,n"));
  } 

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_amef() {
    W_amef_RHF("a,m,e,f") = conj(Moints["vvvo"]("e,f,a,m")) - this->T1_RHF("a,n") * conj(Moints["vvoo"]("f,e,m,n"));
  }
  
  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_mbij() {
    W_mbij_RHF("m,b,i,j") = Moints["vooo"]("b,m,j,i");
    W_mbij_RHF("m,b,i,j") += F_me_RHF("m,e") * this->T2_RHF("e,b,i,j");
    W_mbij_RHF("m,b,i,j") -= this->T1_RHF("b,n") * W_mnij_RHF("m,n,i,j");
    W_mbij_RHF("m,b,i,j") += conj(Moints["vvvo"]("f,e,b,m")) * tau_RHF("e,f,i,j");

    W_mbij_RHF("m,b,i,j") += 2.0 * conj(Moints["vooo"]("e,i,n,m")) * this->T2_RHF("b,e,j,n");
    W_mbij_RHF("m,b,i,j") -= conj(Moints["vooo"]("e,i,m,n")) * this->T2_RHF("b,e,j,n");
    W_mbij_RHF("m,b,i,j") -= conj(Moints["vooo"]("e,i,n,m")) * this->T2_RHF("e,b,j,n");
    W_mbij_RHF("m,b,i,j") -= conj(Moints["vooo"]("e,j,m,n")) * this->T2_RHF("e,b,i,n");

    TArray tmpRHF = TAManager::get().malloc<MatsT>("ovvo");
    tmpRHF("m,b,e,j") = Moints["voov"]("b,m,j,e");
    tmpRHF("m,b,e,j") += 2.0 * this->T2_RHF("f,b,n,j") * conj(Moints["vvoo"]("e,f,m,n"));
    tmpRHF("m,b,e,j") -= this->T2_RHF("b,f,n,j") * conj(Moints["vvoo"]("e,f,m,n"));
    tmpRHF("m,b,e,j") -= this->T2_RHF("f,b,n,j") * conj(Moints["vvoo"]("f,e,m,n"));
    W_mbij_RHF("m,b,i,j") += this->T1_RHF("e,i") * tmpRHF("m,b,e,j");

    tmpRHF("m,b,e,i") = Moints["vovo"]("b,m,e,i");
    tmpRHF("m,b,e,i") -= this->T2_RHF("b,f,n,i") * conj(Moints["vvoo"]("f,e,m,n"));
    W_mbij_RHF("m,b,i,j") += this->T1_RHF("e,j") * tmpRHF("m,b,e,i");
    TAManager::get().free("ovvo", std::move(tmpRHF));
  } 
  template <typename MatsT>
  void EOMRCCSD<MatsT>::formW_abei() {
    W_abei_RHF("a,b,e,i") = Moints["vvvo"]("a,b,e,i");
    W_abei_RHF("a,b,e,i") -= F_me_RHF("m,e") * this->T2_RHF("a,b,m,i");
    W_abei_RHF("a,b,e,i") += this->T1_RHF("f,i") * W_abef_RHF("a,b,e,f");
    W_abei_RHF("a,b,e,i") += conj(Moints["vooo"]("e,i,m,n")) * tau_RHF("a,b,m,n");

    W_abei_RHF("a,b,e,i") += 2.0 * conj(Moints["vvvo"]("e,f,a,m")) * this->T2_RHF("b,f,i,m");
    W_abei_RHF("a,b,e,i") -= conj(Moints["vvvo"]("f,e,a,m")) * this->T2_RHF("b,f,i,m");
    W_abei_RHF("a,b,e,i") -= conj(Moints["vvvo"]("e,f,a,m")) * this->T2_RHF("f,b,i,m");
    W_abei_RHF("a,b,e,i") -= conj(Moints["vvvo"]("f,e,b,m")) * this->T2_RHF("a,f,m,i");

    TArray tmpRHF = TAManager::get().malloc<MatsT>("ovvo");
    tmpRHF("m,b,e,i") = Moints["voov"]("b,m,i,e");
    tmpRHF("m,b,e,i") += 2.0 * this->T2_RHF("f,b,n,i") * conj(Moints["vvoo"]("e,f,m,n"));
    tmpRHF("m,b,e,i") -= this->T2_RHF("b,f,n,i") * conj(Moints["vvoo"]("e,f,m,n"));
    tmpRHF("m,b,e,i") -= this->T2_RHF("f,b,n,i") * conj(Moints["vvoo"]("f,e,m,n"));
    W_abei_RHF("a,b,e,i") -= this->T1_RHF("a,m") * tmpRHF("m,b,e,i");

    tmpRHF("m,a,e,i") = Moints["vovo"]("a,m,e,i");
    tmpRHF("m,a,e,i") -= this->T2_RHF("a,f,n,i") * conj(Moints["vvoo"]("f,e,m,n"));
    W_abei_RHF("a,b,e,i") -= this->T1_RHF("b,m") * tmpRHF("m,a,e,i");
    TAManager::get().free("ovvo", std::move(tmpRHF));
  }  

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formEOMIntermediates() {

    formF_ae();
    formF_mi();
    formW_mnij();
    formW_abef();
    formW_mbej();
    formW_mnie();
    formW_amef();
    formW_mbij();
    formW_abei();
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formR1_tilde_RHF(const TArray &R1, const TArray &R2, TArray &tildeR1) const {
    tildeR1("a,i") = F_ae_RHF("a,c") * R1("c,i");
    tildeR1("a,i") -= F_mi_RHF("k,i") * R1("a,k");
    tildeR1("a,i") += 2.0 * F_me_RHF("k,c") * R2("a,c,i,k");
    tildeR1("a,i") -= F_me_RHF("k,c") * R2("a,c,k,i");
    tildeR1("a,i") += 2.0 * W_mbej_RHF_baba("k,a,c,i") * R1("c,k");
    tildeR1("a,i") += W_mbej_RHF_baab("k,a,c,i") * R1("c,k");
    tildeR1("a,i") += 2.0 * W_amef_RHF("a,k,c,d") * R2("c,d,i,k");
    tildeR1("a,i") -= W_amef_RHF("a,k,c,d") * R2("d,c,i,k");
    tildeR1("a,i") -= 2.0 * W_mnie_RHF("k,l,i,d") * R2("a,d,k,l");
    tildeR1("a,i") += W_mnie_RHF("k,l,i,d") * R2("a,d,l,k");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formR2_tilde_RHF(const TArray &R1, const TArray &R2, TArray &tildeR2) const {
    tildeR2("a,b,i,j") = F_ae_RHF("b,e") * R2("a,e,i,j");
    tildeR2("a,b,i,j") += F_ae_RHF("a,e") * R2("e,b,i,j");
    tildeR2("a,b,i,j") -= F_mi_RHF("k,j") * R2("a,b,i,k");
    tildeR2("a,b,i,j") -= F_mi_RHF("k,i") * R2("a,b,k,j");
    tildeR2("a,b,i,j") += W_mnij_RHF("k,l,i,j") * R2("a,b,k,l");
    tildeR2("a,b,i,j") += W_abef_RHF("a,b,e,f") * R2("e,f,i,j");

    tildeR2("a,b,i,j") += 2.0 * W_mbej_RHF_baba("k,b,c,j") * R2("a,c,i,k");
    tildeR2("a,b,i,j") += W_mbej_RHF_baab("k,b,c,j") * R2("a,c,i,k");
    tildeR2("a,b,i,j") -= W_mbej_RHF_baba("k,b,c,j") * R2("c,a,i,k");

    tildeR2("a,b,i,j") += W_mbej_RHF_baab("k,a,c,j") * R2("c,b,i,k");

    tildeR2("a,b,i,j") += W_mbej_RHF_baab("k,b,c,i") * R2("c,a,j,k");

    tildeR2("a,b,i,j") += 2.0 * W_mbej_RHF_baba("k,a,c,i") * R2("b,c,j,k");
    tildeR2("a,b,i,j") += W_mbej_RHF_baab("k,a,c,i") * R2("b,c,j,k");
    tildeR2("a,b,i,j") -= W_mbej_RHF_baba("k,a,c,i") * R2("c,b,j,k");

    tildeR2("a,b,i,j") += W_abei_RHF("a,b,c,j") * R1("c,i"); // Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += W_abei_RHF("b,a,c,i") * R1("c,j");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") -= W_mbij_RHF("k,a,j,i") * R1("b,k");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") -= W_mbij_RHF("k,b,i,j") * R1("a,k");// Sign is different between Tianyuan's and literature

    TAManager &TAmanager = TAManager::get();
    //[Asthana:2019:4102] Break equations A8, A9 and A10 to save computational cost
    TArray tmp1 = TAmanager.malloc<MatsT>("vv");
    tmp1("b,d") = 2.0 * W_amef_RHF("b,k,d,c") * R1("c,k");
    tmp1("b,d") -= W_amef_RHF("b,k,c,d") * R1("c,k");
    tildeR2("a,b,i,j") += tmp1("b,d") * this->T2_RHF("a,d,i,j");
    tildeR2("a,b,i,j") += tmp1("a,d") * this->T2_RHF("d,b,i,j");
    TArray tmp2 = TAmanager.malloc<MatsT>("oo");
    tmp2("l,j") = 2.0 * W_mnie_RHF("l,k,j,c") * R1("c,k");
    tmp2("l,j") -= W_mnie_RHF("k,l,j,c") * R1("c,k");
    tildeR2("a,b,i,j") -= tmp2("l,j") * this->T2_RHF("a,b,i,l");
    tildeR2("a,b,i,j") -= tmp2("l,i") * this->T2_RHF("a,b,l,j");

    //reuse tmp2 container
    tmp2("l,j") = 2.0 * conj(this->Moints.at("vvoo")("d,c,k,l")) * R2("c,d,j,k");
    tmp2("l,j") -= conj(this->Moints.at("vvoo")("c,d,k,l")) * R2("c,d,j,k");
    tildeR2("a,b,i,j") -= tmp2("l,j") * this->T2_RHF("a,b,i,l");
    tildeR2("a,b,i,j") -= tmp2("l,i") * this->T2_RHF("a,b,l,j");

    //reuse tmp1 container
    tmp1("b,e") = 2.0 * conj(this->Moints.at("vvoo")("e,d,k,l")) * R2("b,d,k,l");
    tmp1("b,e") -= conj(this->Moints.at("vvoo")("e,d,l,k")) * R2("b,d,k,l");
    tildeR2("a,b,i,j") -= tmp1("b,e") * this->T2_RHF("a,e,i,j");
    tildeR2("a,b,i,j") -= tmp1("a,e") * this->T2_RHF("e,b,i,j");

    TAmanager.free("vv", std::move(tmp1));
    TAmanager.free("oo", std::move(tmp2));
  }


  template <typename MatsT>
  void EOMRCCSD<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {
    const TArray &V1 =  V.get_tensor("OneBody");
    const TArray &V2 =  V.get_tensor("TwoBody");
    TArray &HV1 = HV.get_tensor("OneBody");
    TArray &HV2 = HV.get_tensor("TwoBody");
    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR1_tilde_RHF(V1, V2, HV1);
        formR2_tilde_RHF(V1, V2, HV2);
        break;
      case EOMCCEigenVecType::LEFT:
        TAManager &TAmanager = TAManager::get();
        TArray G_ae_RHF = TAmanager.malloc<MatsT>("vv");
        updateG_ae_RHF(V2, G_ae_RHF);
        TArray G_mi_RHF = TAmanager.malloc<MatsT>("oo");
        updateG_mi_RHF(V2, G_mi_RHF);
        formL1_tilde_RHF(V1, V2, G_ae_RHF, G_mi_RHF, HV1);
        formL2_tilde_RHF(V1, V2, G_ae_RHF, G_mi_RHF, HV2);
        TAmanager.free("vv", std::move(G_ae_RHF));
        TAmanager.free("oo", std::move(G_mi_RHF));
        break;
    }
  }


  template <typename MatsT>
  void EOMRCCSD<MatsT>::buildRightZeroBody(size_t nVec) {
    // first of all, check if the imaginary component of excitation energies are under the threshold
    for (size_t i=0; i<nVec; i++) {
      if (abs(std::imag(this->theta[i])) > this->eomSettings.energy_imag_tol) {
        CErr("EOM-RCCSD: Imaginary component of at least one VEE is greater than tolerance (default: 1.0e-12).");
      }
    }

    // Assumes MBExpansionSet R_ type

    std::shared_ptr<MBExpansionSet<MatsT>> VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);

    std::cout << "buildRightZeroBody:" << std::endl;
    for (size_t i = 0; i < nVec; i++) {
      MatsT r0_1_RHF = 2.0 * F_me_RHF("i,a").dot(VR->get(i).get_tensor("OneBody")("a,i")).get();

      TA::get_default_world().gop.fence();
      MatsT r0_2_RHF = 2.0 * conj(this->Moints["vvoo"]("a,b,i,j")).dot(VR->get(i).get_tensor("TwoBody")("a,b,i,j")).get();
      r0_2_RHF -= conj(this->Moints["vvoo"]("b,a,i,j")).dot(VR->get(i).get_tensor("TwoBody")("a,b,i,j")).get();
      TA::get_default_world().gop.fence();
      if constexpr (std::is_same_v<MatsT, double>) {
        VR->get(i).zeroBody() = (r0_1_RHF + r0_2_RHF) / std::real(this->theta[i]);
      } else {
        VR->get(i).zeroBody() = (r0_1_RHF + r0_2_RHF) / this->theta[i];
      }
    }
  }

  template <typename MatsT>
  inline size_t EOMRCCSD<MatsT>::toCompoundS_RHF(size_t a, size_t i) const {
    return a + i * nV_RHF;
  }

  template <typename MatsT>
  inline size_t EOMRCCSD<MatsT>::toCompoundD_RHF(size_t a, size_t b, size_t i, size_t j) const {
    size_t ai = a * nO_RHF + i, bj = b * nO_RHF + j;
    if (ai > bj)
      CErr("ai > bj in toCompoundD_RHF, which should not happen for RHF case. a: "
        + std::to_string(a) + ", b: " + std::to_string(b)
        + ", i: " + std::to_string(i) + ", j: " + std::to_string(j));
    return bj*(bj+1)/2 + ai;
  }

  template <typename MatsT>
  EOMRCCSD<MatsT>::~EOMRCCSD() {

    //if (theta) CQMemManager::get().free(theta);

    TAManager &TAmanager = TAManager::get();

    if (Rho_ij_RHF) TAmanager.free("oo", std::move(Rho_ij_RHF), true);
    if (Rho_ab_RHF) TAmanager.free("vv", std::move(Rho_ab_RHF), true);
    if (Rho_ia_RHF) TAmanager.free("ov", std::move(Rho_ia_RHF), true);
    if (Rho_ai_RHF) TAmanager.free("vo", std::move(Rho_ai_RHF), true);

  }

  // for CCS_guess
  template <typename MatsT>
  void EOMRCCSD<MatsT>::fillGuess(MatsT *guess_vec, size_t n_vec) const{
    CErr("EOMRCCSD does not support CCS guess.");
  }


   /// for full_diagonaization routine
  template <typename MatsT>
  void EOMRCCSD<MatsT>::buildDiag(MatsT * diag, const std::vector<double> &) const {

    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_RHF).extent();
    size_t nO = TAmanager.getRange(oLabel_RHF).extent();

    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));
    MatsT * diag2 = diag + nOVshift_RHF;

    TA::foreach_inplace( F_ae_RHF, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t a = x[0], b = x[0];
        for (size_t i = 0; i < nO; ++i)
          diag[toCompoundS_RHF(a,i)] += v;

        for (size_t b = a; b < nV; ++b)
          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }

        for (size_t j = 0; j < nO; j++)
          for (size_t a = 0; a <= b; a++)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }
      }

    });

    TA::foreach_inplace( F_mi_RHF, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t i = x[0], j = x[0];
        for (size_t a = 0; a < nV; ++a)
          diag[toCompoundS_RHF(a,i)] -= v;

        for (size_t b = 0; b < nV; ++b)
          for (size_t a = 0; a <= b; a++)
            for (size_t j = (a < b? 0 : i); j < nO; ++j) {
              diag2[toCompoundD_RHF(a, b, i, j)] -= v;
            }

        for (size_t b = 0; b < nV; b++)
          for (size_t a = 0; a <= b; a++)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] -= v;
            }
      }

    });

    TA::foreach_inplace( W_mbej_RHF_baab, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[3] or lobound[1] != lobound[2])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          x[3] = x[0];
          x[2] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], a = x[1], b = x[1], j = x[0];
          diag[toCompoundS_RHF(a,i)] += v;

          for (size_t b = a; b < nV; ++b)
            for (size_t j = (a < b? 0 : i); j < nO; ++j) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }

          for (size_t b = a; b < nV; ++b)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }

          for (size_t a = 0; a <= b; ++a)
            for (size_t j = (a < b? 0 : i); j < nO; ++j) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }

          for (size_t a = 0; a <= b; ++a)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }
        }
      }

    });
    TA::foreach_inplace( W_mbej_RHF_baba, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[3] or lobound[1] != lobound[2])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          x[3] = x[0];
          x[2] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], a = x[1], b = x[1], j = x[0];
          diag[toCompoundS_RHF(a,i)] += v;

          for (size_t b = a; b < nV; ++b)
            for (size_t j = (a < b? 0 : i); j < nO; ++j) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }

          for (size_t a = 0; a <= b; ++a)
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_mnij_RHF, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
      return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], j = x[1];

          for (size_t b = 0; b < nV; ++b)
            for (size_t a = 0; a <= b; ++a) {
              if (a == b and i > j) continue;
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_abef_RHF, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < std::min(x[1] + 1, static_cast<std::size_t>(upbound[0])); ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t a = x[0], b = x[1];

          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < nO; ++i) {
              if (a == b and i > j) continue;
              diag2[toCompoundD_RHF(a, b, i, j)] += v;
            }
        }
      }

    });


    TA::foreach_inplace(T2_RHF, this->Moints.at("vvoo"),
                        [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& T2tile, const TA::Tensor<MatsT>& Vtile) {
      const auto& lobound = T2tile.range().lobound();
      const auto& upbound = T2tile.range().upbound();

      std::vector<std::size_t> x{0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              size_t a = x[0], b = x[1], i = x[2], j = x[3];
              MatsT v = SmartConj(T2tile[x]) * Vtile[x];

              for (size_t a = 0; a <= b; ++a) {
                if (a == b and i > j) continue;
                diag2[toCompoundD_RHF(a, b, i, j)] -= v;
              }

              for (size_t b = a; b < nV; ++b) {
                if (a == b and i > j) continue;
                diag2[toCompoundD_RHF(a, b, i, j)] -= v;
              }

              if (a <= b) {
                for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++)
                  diag2[toCompoundD_RHF(a, b, i, j)] -= v;

                for (size_t j = (a < b? 0 : i); j < nO; ++j)
                  diag2[toCompoundD_RHF(a, b, i, j)] -= v;
              }

            }
    });

    TA::get_default_world().gop.fence();

    MatsT *diag_copy = CQMemManager::get().malloc<MatsT>(this->Hbar_dim);
    std::copy_n(diag, this->Hbar_dim, diag_copy);
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));
    MPIAllReduce(diag_copy, this->Hbar_dim, diag, MPI_COMM_WORLD);
    CQMemManager::get().free(diag_copy);

    TA::get_default_world().gop.fence();
    MBExpansion<MatsT> diag_TA(this->tensor_builder_, false, MBTensorSymmetry::RCCSD);
    diag_TA.fromRaw(diag, false);

  }


  template <typename MatsT>
  typename Davidson<MatsT>::VecsGen_t EOMRCCSD<MatsT>::EmptyDavidsonVectorBuilder(){
      // Algorithm with implicit Hbar matrix
      typename Davidson<MatsT>::VecsGen_t vecsGenEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSet<MatsT>>(this->tensor_builder_, nVec, this->savFile_, MBTensorSymmetry::RCCSD);
        }; // implicit vecsGenerator

        return vecsGenEOM;
      } else {
        CErr("EOMRCCSD: EmptyDavidsonVectorBuilder is only implemented for implicit Hbar type.");
      }
      return vecsGenEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMRCCSD<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
//      typename Davidson<MatsT>::LinearTrans_t funcEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        this->funcEOM = [this, &eigenVecType]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSet<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;
          try {
            V_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t i = 0; i < nVec; i++) {
            const MBExpansion<MatsT> &ViRHF =  V_ptr->get(i + Vshift);
            MBExpansion<MatsT> &AViRHF = AV_ptr->get(i + AVshift);
            buildSigma(ViRHF, AViRHF, eigenVecType);
            AViRHF.enforceSymmetry();
          }
        }; // implicit sigmaBuilder
      } else {
        CErr("EOMRCCSD: DavidsonResidualBuilder with implicit Hbar type is the only one implemented.");
      }
      return this->funcEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMRCCSD<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){

      double PCsmall = this->eomSettings.davidson_preCond_small;

//      typename Davidson<MatsT>::LinearTrans_t PCEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        this->PCEOM = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          AV.set_data(0, nVec, V, 0);

          MBExpansionSet<MatsT> *AV_ptr = nullptr;
          size_t AVshift = 0;

          try {
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t iVec = 0; iVec < nVec; iVec++) {

            MBExpansion<MatsT> &curB = AV_ptr->get(iVec + AVshift);
            MatsT curEigI = 0.0;
            if constexpr (std::is_same_v<MatsT, double>) {
              curEigI = curEig[iVec].real();
            } else {
              curEigI = curEig[iVec];
            }

            TA::foreach_inplace(curB.get_tensor("OneBody"), [iVec, curEigI, eomDiag, this, PCsmall](TA::Tensor<MatsT> &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              MatsT denom = 0.0;
              std::vector<std::size_t> x{0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                  denom = curEigI - eomDiag[toCompoundS_RHF(x[0], x[1])];
                  if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                }
            });

            MatsT *diagD = eomDiag + this->intermediates_.nVir * this->intermediates_.nOcc;
            TA::foreach_inplace(curB.get_tensor("TwoBody"), [iVec, curEigI, diagD, this, PCsmall](TA::Tensor<MatsT> &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              MatsT denom = 0.0;
              std::vector<std::size_t> x{0, 0, 0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                  size_t a = x[0], b = x[1];
                  if (a > b)
                    std::swap(a, b);
                  for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                    for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                      size_t i = x[2], j = x[3];
                      if (i > j)
                        std::swap(i, j);
                      denom = curEigI - diagD[toCompoundD_RHF(a, b, i, j)];
                      if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                    }
                }
            });
            TA::get_default_world().gop.fence();

            curB.enforceSymmetry();
          }
        }; // implicit preConditioner

      } else {
        CErr("EOMRCCSD: DavidsonPreconditionerBuilder with implicit Hbar type is the only one implemented.");
      }
      return this->PCEOM;
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::initializeDensity() {
    TAManager &TAmanager = TAManager::get();

    if (not Rho_ij_RHF.is_initialized()){
      Rho_ij_RHF = TAmanager.malloc<MatsT>("oo");
    }

    if (not Rho_ab_RHF.is_initialized()){
      Rho_ab_RHF = TAmanager.malloc<MatsT>("vv");
    }

    if (not Rho_ia_RHF.is_initialized()){
      Rho_ia_RHF = TAmanager.malloc<MatsT>("ov");
    }

    if (not Rho_ai_RHF.is_initialized()){
      Rho_ai_RHF = TAmanager.malloc<MatsT>("vo");
    }
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::buildDensity_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){
    formRho_ij_RHF(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
    formRho_ab_RHF(t1, t2, r0, r1, r2, l0, l1, l2);
    formRho_ia_RHF(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
    formRho_ai_RHF(t1, t2, r0, r1, r2, l0, l1, l2);
  }


  template <typename MatsT>
  void EOMRCCSD<MatsT>::formRho_ij_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame) {


    if (isSame){
      TA::foreach_inplace(Rho_ij_RHF, [](TA::Tensor<MatsT> &tile){

        const auto& lobound = tile.range().lobound();
        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
            if(x[0] == x[1])
              tile[x] = 1.0;
            else
              tile[x] = 0.0;
      });

      this->Rho_ij_RHF("i,j") -= r0 * l1("e,i") * t1("e,j");

    }
    else{
      this->Rho_ij_RHF("i,j") = - r0 * l1("e,i") * t1("e,j");
    }


    this->Rho_ij_RHF("i,j") -= 2.0 * r0 * l2("f,e,i,m") * t2("f,e,j,m");
    this->Rho_ij_RHF("i,j") += r0 * l2("e,f,i,m") * t2("f,e,j,m");
    this->Rho_ij_RHF("i,j") -= l1("e,i") * r1("e,j");
    this->Rho_ij_RHF("i,j") -=  2.0 * l2("f,e,i,m") * r2("f,e,j,m");
    this->Rho_ij_RHF("i,j") +=  l2("e,f,i,m") * r2("f,e,j,m");
    TArray tmp = TAManager::get().malloc<MatsT>("ov");
    tmp("i,f") = 2.0 * l2("f,e,i,m") * r1("e,m");
    tmp("i,f") -= l2("e,f,i,m") * r1("e,m");
    this->Rho_ij_RHF("i,j") -= tmp("i,f") * t1("f,j");
    this->Rho_ij_RHF("i,j") = this->Rho_ij_RHF("j,i");

    TAManager::get().free("ov", std::move(tmp));
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formRho_ab_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2) {

    this->Rho_ab_RHF("a,b") = r0 * l1("b,m") * t1("a,m");
    this->Rho_ab_RHF("a,b") += 2.0 * r0 * l2("e,b,m,n") * t2("e,a,m,n");
    this->Rho_ab_RHF("a,b") -= r0 * l2("e,b,n,m") * t2("e,a,m,n");
    this->Rho_ab_RHF("a,b") += l1("b,m") * r1("a,m");
    this->Rho_ab_RHF("a,b") += 2.0 * l2("e,b,m,n") * r2("e,a,m,n");
    this->Rho_ab_RHF("a,b") -= l2("e,b,n,m") * r2("e,a,m,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ov");
    tmp("n,b") = 2.0 * l2("e,b,m,n") * r1("e,m");
    tmp("n,b") -= l2("e,b,n,m") * r1("e,m");
    this->Rho_ab_RHF("a,b") += tmp("n,b") * t1("a,n");
    this->Rho_ab_RHF("a,b") = this->Rho_ab_RHF("b,a");

    TAManager::get().free("ov", std::move(tmp));
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formRho_ai_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2){
    this->Rho_ai_RHF("a,i") = r0 * l1("a,i");
    this->Rho_ai_RHF("a,i") += 2.0 * l2("a,e,i,m") * r1("e,m");
    this->Rho_ai_RHF("a,i") -= l2("e,a,i,m") * r1("e,m");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formRho_ia_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){

    if(isSame){
      this->Rho_ia_RHF("i,a") = t1("a,i");
      this->Rho_ia_RHF("i,a") += 2.0 * r0 * l1("e,m") * t2("a,e,i,m");
    }
    else{
      this->Rho_ia_RHF("i,a") = 2.0 * r0 * l1("e,m") * t2("a,e,i,m");
    }
    this->Rho_ia_RHF("i,a") -= r0 * l1("e,m") * t2("a,e,m,i");

    TAManager &TAmanager = TAManager::get();

    this->Rho_ia_RHF("i,a") += l0 * r1("a,i");
    TArray tmp1 = TAmanager.malloc<MatsT>("oo");
    tmp1("m,i") = l1("e,m") * t1("e,i");
    this->Rho_ia_RHF("i,a") += - r0 * tmp1("m,i") * t1("a,m");

    this->Rho_ia_RHF("i,a") += - tmp1("m,i") * r1("a,m");



    tmp1("m,i") = 2.0 * l2("e,f,m,n") * t2("e,f,i,n");
    tmp1("m,i") -= l2("e,f,m,n") * t2("f,e,i,n");
    this->Rho_ia_RHF("i,a") -= r0 * tmp1("m,i") * t1("a,m");

    this->Rho_ia_RHF("i,a") -= tmp1("m,i") * r1("a,m");

    TArray tmp2 = TAmanager.malloc<MatsT>("vv");
    tmp2("a,e") = 2.0 * l2("e,f,m,n") * t2("a,f,m,n");
    tmp2("a,e") -= l2("e,f,n,m") * t2("a,f,m,n");
    this->Rho_ia_RHF("i,a") -= r0 * tmp2("a,e") * t1("e,i");

    this->Rho_ia_RHF("i,a") += 2.0 * l1("e,m") * r2("a,e,i,m");
    this->Rho_ia_RHF("i,a") -= l1("e,m") * r2("e,a,i,m");


    tmp1("m,i") = l1("e,m") * r1("e,i");
    this->Rho_ia_RHF("i,a") -= tmp1("m,i") * t1("a,m");


    this->Rho_ia_RHF("i,a") -= tmp2("a,e") * r1("e,i");


    tmp2("a,e") = 2.0 * l2("e,f,m,n") * r2("a,f,m,n");
    tmp2("a,e") -= l2("e,f,n,m") * r2("a,f,m,n");
    this->Rho_ia_RHF("i,a") -= tmp2("a,e") * t1("e,i");


    tmp1("m,i") = 2.0 * l2("e,f,m,n") * r2("e,f,i,n");
    tmp1("m,i") -= l2("f,e,m,n") * r2("e,f,i,n");
    this->Rho_ia_RHF("i,a") -= tmp1("m,i") * t1("a,m");

    TArray tmp3 = TAmanager.malloc<MatsT>("ov");
    tmp3("m,e") = 2.0 * l2("e,f,m,n") * r1("f,n");
    tmp3("m,e") -= l2("f,e,m,n") * r1("f,n");


    tmp1("m,i") = tmp3("m,e") * t1("e,i");
    this->Rho_ia_RHF("i,a") -= tmp1("m,i") * t1("a,m");

    this->Rho_ia_RHF("i,a") += 2.0 * tmp3("m,e") * t2("a,e,i,m");
    this->Rho_ia_RHF("i,a") -= tmp3("m,e") * t2("a,e,m,i");

    TAmanager.free("oo", std::move(tmp1));
    TAmanager.free("vv", std::move(tmp2));
    TAmanager.free("ov", std::move(tmp3));
  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMRCCSD<MatsT>::calcTransitionDipole(
    const MatsT r0, const TArray& r1, const TArray& r2,
    const MatsT l0,const TArray& l1, const TArray& l2, bool isSame) {

    std::array<MatsT, 3> mu;

    buildDensity_RHF(T1_RHF, T2_RHF, r0, r1, r2, l0, l1, l2, isSame);
    for (size_t j = 0; j < 3; j++) {
      mu[j]  = dot(muMatrix_RHF[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij_RHF("i,j")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(muMatrix_RHF[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia_RHF("i,a")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(muMatrix_RHF[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai_RHF("a,i")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(muMatrix_RHF[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab_RHF("a,b")).get();
      TA::get_default_world().gop.fence();
      mu[j] *= 2.0; // factor of 2 for spin summation in RHF
    }
    TA::get_default_world().gop.fence();

    if (isSame) // Add frozen core contribution to the transition dipole moment if the left and right states are the same
      for (size_t j = 0; j < 3; j++) {
        mu[j] += this->intermediates_.Mu_fzc[j];
      }

    return mu;

  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMRCCSD<MatsT>::calcGround2ExcitedTransitionDipole(size_t i) {

    // Prepare left ground state amplitudes
    MatsT L0_ = this->Lg_->zeroBody();
    TArray &L1_ = this->Lg_->get_tensor("OneBody");
    TArray &L2_ = this->Lg_->get_tensor("TwoBody");

    MBExpansionSet<MatsT> &Reom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);
    return calcTransitionDipole(Reom.get(i).zeroBody(),
                                Reom.get(i).get_tensor("OneBody"),
                                Reom.get(i).get_tensor("TwoBody"),
                                L0_, L1_, L2_);
  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMRCCSD<MatsT>::calcExcited2GroundTransitionDipole(size_t i) {

    // Prepare right ground state amplitudes
    TAManager &TAmanager = TAManager::get();
    TArray Rg1 = TAmanager.malloc<MatsT>("vo");
    Rg1("a,i") = 0.0 * Rg1("a,i");
    TArray Rg2 = TAmanager.malloc<MatsT>("vvoo");
    Rg2("a,b,i,j") = 0.0 * Rg2("a,b,i,j");

    MBExpansionSet<MatsT> &Leom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->L_);

    std::array<MatsT, 3> mu;
    mu = calcTransitionDipole(1.0, Rg1, Rg2,
                              Leom.get(i).zeroBody(),
                              Leom.get(i).get_tensor("OneBody"),
                              Leom.get(i).get_tensor("TwoBody"));
    TAmanager.free("vo", std::move(Rg1));
    TAmanager.free("vvoo", std::move(Rg2));

    return mu;
  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMRCCSD<MatsT>::calcExcited2ExcitedTransitionDipole(size_t i, size_t j) {

    MBExpansionSet<MatsT> &Leom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->L_);
    MBExpansionSet<MatsT> &Reom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);

    return calcTransitionDipole(Reom.get(j).zeroBody(),
                                Reom.get(j).get_tensor("OneBody"),
                                Reom.get(j).get_tensor("TwoBody"),
                                Leom.get(i).zeroBody(),
                                Leom.get(i).get_tensor("OneBody"),
                                Leom.get(i).get_tensor("TwoBody"),
                                i == j);
  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMRCCSD<MatsT>::calcGroundDipole() {

    // Prepare left ground state amplitudes
    MatsT L0_ = this->Lg_->zeroBody();
    TArray &L1_ = this->Lg_->get_tensor("OneBody");
    TArray &L2_ = this->Lg_->get_tensor("TwoBody");

    // Prepare right ground state amplitudes
    TAManager &TAmanager = TAManager::get();
    TArray Rg1 = TAmanager.malloc<MatsT>("vo");
    Rg1("a,i") = 0.0 * Rg1("a,i");
    TArray Rg2 = TAmanager.malloc<MatsT>("vvoo");
    Rg2("a,b,i,j") = 0.0 * Rg2("a,b,i,j");

    std::array<MatsT, 3> mu;
    mu = calcTransitionDipole(1.0, Rg1, Rg2, L0_, L1_, L2_, true);
    TAmanager.free("vo", std::move(Rg1));
    TAmanager.free("vvoo", std::move(Rg2));

    return mu;

  }

  template <typename MatsT>
  dcomplex EOMRCCSD<MatsT>::calcOscillatorStrength(size_t i){
    std::array<MatsT, 3> mu_g2x, mu_x2g;

    mu_g2x = calcGround2ExcitedTransitionDipole(i);
    mu_x2g = calcExcited2GroundTransitionDipole(i);

    MatsT DS = mu_g2x[0] * mu_x2g[0] + mu_g2x[1] * mu_x2g[1] + mu_g2x[2] * mu_x2g[2];

    dcomplex f = 2./3 * this->theta[i] * DS;
    return f;

  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::initializeLambda() {

    if (this->Lg_ == nullptr) this->initializeGroundStateLambda();
    this->Lg_->zeroBody() = 1.0;
    this->Lg_->get_tensor("OneBody")("a,i") = conj(T1_RHF("a,i"));
    this->Lg_->get_tensor("TwoBody")("a,b,i,j") = conj(T2_RHF("a,b,i,j"));

    TAManager &TAmanager = TAManager::get();
    if (not G_ae_RHF.is_initialized()){
      G_ae_RHF = TAmanager.malloc<MatsT>("vv");
    }
    if (not G_mi_RHF.is_initialized()){
      G_mi_RHF = TAmanager.malloc<MatsT>("oo");
    }

  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::updateG_ae_RHF(const TArray &L2, TArray &G_ae_RHF) const {
    G_ae_RHF("a,e") = - 2.0 * T2_RHF("e,f,m,n") * L2("a,f,m,n");
    G_ae_RHF("a,e") += T2_RHF("e,f,n,m") * L2("a,f,m,n");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::updateG_mi_RHF(const TArray &L2, TArray &G_mi_RHF) const {
    G_mi_RHF("m,i") = 2.0 * T2_RHF("e,f,m,n") * L2("e,f,i,n");
    G_mi_RHF("m,i") -= T2_RHF("f,e,m,n") * L2("e,f,i,n");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formL1_tilde_RHF(const TArray &L1, const TArray &L2, const TArray &G_ae_RHF, const TArray &G_mi_RHF,
                                          TArray &tildeL1) const {

    tildeL1("a,i") = F_ae_RHF("e,a") * L1("e,i");
    tildeL1("a,i") -= F_mi_RHF("i,m") * L1("a,m");
    tildeL1("a,i") += 2.0 * L1("e,m") * W_mbej_RHF_baba("i,e,a,m");
    tildeL1("a,i") += L1("e,m") * W_mbej_RHF_baab("i,e,a,m");
    tildeL1("a,i") += 2.0 * L2("e,f,i,m") * W_abei_RHF("e,f,a,m");
    tildeL1("a,i") -= L2("f,e,i,m") * W_abei_RHF("e,f,a,m");
    tildeL1("a,i") -= 2.0 * L2("a,e,m,n") * W_mbij_RHF("i,e,m,n");
    tildeL1("a,i") += L2("a,e,n,m") * W_mbij_RHF("i,e,m,n");
    tildeL1("a,i") -= 2.0 * G_ae_RHF("e,f") * W_amef_RHF("e,i,f,a");
    tildeL1("a,i") += G_ae_RHF("e,f") * W_amef_RHF("e,i,a,f");
    tildeL1("a,i") -= 2.0 * G_mi_RHF("m,n") * W_mnie_RHF("m,i,n,a");
    tildeL1("a,i") += G_mi_RHF("m,n") * W_mnie_RHF("i,m,n,a");

  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::formL2_tilde_RHF(const TArray &L1, const TArray &L2, const TArray &G_ae_RHF, const TArray &G_mi_RHF,
                                          TArray &tildeL2) const {

    tildeL2("a,b,i,j") = L2("a,e,i,j") * F_ae_RHF("e,b");
    tildeL2("a,b,i,j") += L2("e,b,i,j") * F_ae_RHF("e,a");

    tildeL2("a,b,i,j") -= L2("a,b,i,m") * F_mi_RHF("j,m");
    tildeL2("a,b,i,j") -= L2("a,b,m,j") * F_mi_RHF("i,m");

    tildeL2("a,b,i,j") += L2("a,b,m,n") * W_mnij_RHF("i,j,m,n");

    tildeL2("a,b,i,j") += L2("e,f,i,j") * W_abef_RHF("e,f,a,b");

    tildeL2("a,b,i,j") += L1("e,i") * W_amef_RHF("e,j,a,b");
    tildeL2("a,b,i,j") += L1("e,j") * W_amef_RHF("e,i,b,a");

    tildeL2("a,b,i,j") -= L1("a,m") * W_mnie_RHF("i,j,m,b");
    tildeL2("a,b,i,j") -= L1("b,m") * W_mnie_RHF("j,i,m,a");

    tildeL2("a,b,i,j") += 2.0 * L2("a,e,i,m") * W_mbej_RHF_baba("j,e,b,m");
    tildeL2("a,b,i,j") += L2("a,e,i,m") * W_mbej_RHF_baab("j,e,b,m");
    tildeL2("a,b,i,j") -= L2("e,a,i,m") * W_mbej_RHF_baba("j,e,b,m");
    tildeL2("a,b,i,j") += L2("e,a,j,m") * W_mbej_RHF_baab("i,e,b,m");
    tildeL2("a,b,i,j") += L2("e,b,i,m") * W_mbej_RHF_baab("j,e,a,m");
    tildeL2("a,b,i,j") += 2.0 * L2("b,e,j,m") * W_mbej_RHF_baba("i,e,a,m");
    tildeL2("a,b,i,j") += L2("b,e,j,m") * W_mbej_RHF_baab("i,e,a,m");
    tildeL2("a,b,i,j") -= L2("e,b,j,m") * W_mbej_RHF_baba("i,e,a,m");

    tildeL2("a,b,i,j") += L1("a,i") * F_me_RHF("j,b");
    tildeL2("a,b,i,j") += L1("b,j") * F_me_RHF("i,a");

    tildeL2("a,b,i,j") += conj(this->Moints["vvoo"]("a,e,i,j")) * G_ae_RHF("b,e");
    tildeL2("a,b,i,j") += conj(this->Moints["vvoo"]("e,b,i,j")) * G_ae_RHF("a,e");

    tildeL2("a,b,i,j") -= conj(this->Moints["vvoo"]("a,b,i,m")) * G_mi_RHF("m,j");
    tildeL2("a,b,i,j") -= conj(this->Moints["vvoo"]("a,b,m,j")) * G_mi_RHF("m,i");
  }

  template <typename MatsT>
  void EOMRCCSD<MatsT>::runLambda() {

    if (this->eomSettings.restart_l or this->eomSettings.restart_r) { // Restart from file, read in amplitudes

      TA::get_default_world().gop.fence();
      bool lg_amp_exist = false;
      if (MPIRank() == 0) lg_amp_exist = this->savFile_.exists("/CC/LAMBDA_AMPLITUDE");
      MPIBCast(lg_amp_exist, 0, MPI_COMM_WORLD);

      if (lg_amp_exist) {
        if (this->Lg_ == nullptr) this->initializeGroundStateLambda();
        std::cout << "Reading Lambda amplitudes..." << std::endl;
        size_t size = this->Lg_->length();
        MatsT * lg_amp = CQMemManager::get().malloc<MatsT>(size);
        TA::get_default_world().gop.fence();
        if (MPIRank() == 0) this->savFile_.readData("/CC/LAMBDA_AMPLITUDE", lg_amp);
        MPIBCast(lg_amp, size, 0, MPI_COMM_WORLD);
        TA::get_default_world().gop.fence();
        this->Lg_->fromRaw(lg_amp, false);
        this->Lg_->zeroBody() = 1.0;
        TA::get_default_world().gop.fence();
        CQMemManager::get().free(lg_amp);
        std::cout << "Reading Lambda amplitudes finished." << std::endl;

        return;
      }

    }

    auto lambda_start = tick();

    initializeLambda();
    MBExpansion<MatsT>& Lg = *(this->Lg_);
    TArray &L1_ = Lg.get_tensor("OneBody");
    TArray &L2_ = Lg.get_tensor("TwoBody");
    MBExpansion<MatsT> L_old(Lg);

    std::shared_ptr<DIISTA<MatsT>> ldiis = nullptr;
    if(this->ccSettings_.useDIIS){
      ldiis = std::make_shared<DIISTA<MatsT>>(this->ccSettings_.nDIIS);
    }

    MatsT pseudoEnergy = 0.0;

    std::cout << "Start Lambda iterations!" << std::endl;

    std::cout << std::endl;

    std::cout << std::setw(18) << std::left <<  "Lambda Iterations";
    std::cout << std::setw(34) << std::left << "Pseudo Energy (Eh)";
    std::cout << std::setw(19) << std::right << "\u0394PE (Eh)";
    std::cout << std::setw(19) << std::right << "|\u0394L|";
    std::cout << '\n';
    std::cout << std::setw(18) << std::left <<  "  -------------";
    std::cout << std::setw(34) << std::left << "-----------------";
    std::cout << std::setw(18) << std::right << "--------";
    std::cout << std::setw(18) << std::right << "----";
    std::cout << '\n' << std::endl;

    for (auto iter = 0; iter < this->ccSettings_.maxiter; iter++){
      L_old = Lg;

      updateG_ae_RHF(L_old.get_tensor("TwoBody"), G_ae_RHF);
      updateG_mi_RHF(L_old.get_tensor("TwoBody"), G_mi_RHF);

      formL1_tilde_RHF(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody"), G_ae_RHF, G_mi_RHF, L1_);
      L1_("a,i") += F_me_RHF("i,a");
      L1_("a,i") = L_old.get_tensor("OneBody")("a,i") + L1_("a,i") * D_ai_RHF("a,i");

      formL2_tilde_RHF(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody"), G_ae_RHF, G_mi_RHF, L2_);
      L2_("a,b,i,j") += conj(this->Moints["vvoo"]("a,b,i,j"));
      L2_("a,b,i,j") = L_old.get_tensor("TwoBody")("a,b,i,j") + L2_("a,b,i,j") * D_abij_RHF("a,b,i,j");

      if(this->ccSettings_.useDIIS){

        // give solution vector to diis
        ldiis->WriteVector(Lg);

        // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
        L_old.scale(-1.0);
        L_old.axpy(1.0, Lg);

        //set error vector in DIIS
        ldiis->WriteErrorVector(L_old);

        // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
        ldiis->Extrapolate(Lg);
      } else {
        L_old.axpy(-1.0, Lg);
      }
      
      MatsT PE_old = pseudoEnergy;
      MatsT PEOneBody = 2.0 * this->fockMatrix_ta_RHF["vo"]("a,i").dot(L1_("a,i")).get();
      MatsT PETwoBodyL2 = 2.0 * this->Moints["vvoo"]("a,b,i,j").dot(L2_("a,b,i,j")).get();
      PETwoBodyL2 -= this->Moints["vvoo"]("b,a,i,j").dot(L2_("a,b,i,j")).get();
      MatsT PETwoBodyL1 = 2.0 * this->Moints["vvoo"]("c,d,k,l").dot(L1_("c,k") * L1_("d,l")).get();
      PETwoBodyL1 -= this->Moints["vvoo"]("d,c,k,l").dot(L1_("c,k") * L1_("d,l")).get();
      TA::get_default_world().gop.fence();    
      pseudoEnergy = PEOneBody + PETwoBodyL2 + PETwoBodyL1;

      double dPE = std::abs(pseudoEnergy - PE_old);
      double dL = L_old.norm();

      std::cout << std::setprecision(12) << std::fixed;
      std::cout << "  Iteration "  << std::setw(6) << std::left << iter;
      std::cout << std::setw(34) << std::left << std::fixed << pseudoEnergy;
      std::cout << std::setw(18) << std::right << std::fixed << dPE;
      std::cout << std::setw(18) << std::right << std::fixed << dL;
      std::cout << std::endl;

      if ( dPE < this->ccSettings_.eConv and dL < this->ccSettings_.tConv) {
        std::cout << "\n  Lambda iteration converged in " << iter << " steps." << std::endl;
        std::cout << "\n  Lambda Completed: Iteration total time "<< std::setw(10) << std::right
                  << std::setprecision(6) << tock(lambda_start) << " s" << std::endl;

        if(this->eomSettings.save_l){
          size_t size = Lg.length();
          MatsT * lg_amp = CQMemManager::get().malloc<MatsT>(size);
          TA::get_default_world().gop.fence();
          Lg.toRaw(lg_amp, false);
          TA::get_default_world().gop.fence();
          if(this->savFile_.exists()){
            this->savFile_.safeWriteData("/CC/LAMBDA_AMPLITUDE", lg_amp, {size});
          }
          CQMemManager::get().free(lg_amp);
        }

        std::cout << bannerEnd << std::endl;
        break;
      }

      if (iter == this->ccSettings_.maxiter - 1){
        CErr(std::string("Lambda iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }
    }
  }
  
}
