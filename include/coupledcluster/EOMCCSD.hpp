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
  EOMCCSD<MatsT>::EOMCCSD(const SafeFile &savFile,
                                CCIntermediates<MatsT> &intermediates,
                                const EOMSettings &eomSettings,
                                const CoupledClusterSettings &ccSettings):
      EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
      vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
      T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")),
      tau(intermediates.tau),
      F_ae(intermediates.F_ae),
      F_mi(intermediates.F_mi),
      F_me(intermediates.F_me),
      W_mnij(intermediates.W_mnij),
      W_abef(intermediates.W_abef),
      W_mbej(intermediates.W_mbej),
      W_mnie(intermediates.W_mnie),
      W_amef(intermediates.W_amef),
      W_mbij(intermediates.W_mbij),
      W_abei(intermediates.W_abei),
      G_ae(intermediates.G_ae),
      G_mi(intermediates.G_mi),
      D_ai(intermediates.D_ai),
      D_abij(intermediates.D_abij),
      Rho_ij(intermediates.Rho_ij),
      Rho_ab(intermediates.Rho_ab),
      Rho_ia(intermediates.Rho_ia),
      Rho_ai(intermediates.Rho_ai) {
    
    TAManager &TAmanager = TAManager::get();

    // without L, we don't need D
    if (not (eomSettings.oscillator_strength or ccSettings.crcc or ccSettings.computeDipole)) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    nV_ = TAmanager.getRange(vLabel_).extent();
    nO_ = TAmanager.getRange(oLabel_).extent();
    nOVshift_ = nV_ * nO_;
    nO2shift_ = nO_ * (nO_ - 1) / 2;
    nV2shift_ = nV_ * (nV_ - 1) / 2;

    this->Hbar_dimension_offsets.emplace("OneBody", nOVshift_);
    this->Hbar_dimension_offsets.emplace("TwoBody", nO2shift_*nV2shift_);

    this->Hbar_dim = nOVshift_ + nO2shift_ * nV2shift_;

    this->outOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);
    abIndices_.clear();
    abIndices_.resize(nV_, std::vector<size_t>(nV_, this->outOfBound_));
    ijIndices_.clear();
    ijIndices_.resize(nO_, std::vector<size_t>(nO_, outOfBound_));
    nOVshift_ = nO_ * nV_;
    nO2shift_ = nO_ * (nO_ - 1) / 2 ;
    nV2shift_ = nV_ * (nV_ - 1) / 2 ;
    size_t idx = 0;
    for (size_t b = 0; b < nV_; b++) {
      for (size_t a = 0; a < std::min(b, nV_); a++) {
        abIndices_[a][b] = idx;
        abIndices_[b][a] = idx++;
      }
    }
    idx = 0;
    for (size_t j = 0; j < nO_; j++) {
      for (size_t i = 0; i < std::min(j, nO_); i++) {
        ijIndices_[i][j] = idx;
        ijIndices_[j][i] = idx++;
      }
    }

    this->tensor_builder_.push_back(std::string({intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->tensor_builder_.push_back(std::string({intermediates.vLabel,intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));
  }




  template <typename MatsT>
  void EOMCCSD<MatsT>::initializeEOMCC() {
    TAManager &TAmanager = TAManager::get();

    if (not W_mnie.is_initialized()){
      W_mnie = TAmanager.malloc<MatsT>("ooov");
    }

    if (not W_amef.is_initialized()){
      W_amef = TAmanager.malloc<MatsT>("vovv");
    }

    if (not W_mbij.is_initialized()){
      W_mbij = TAmanager.malloc<MatsT>("ovoo");
    }

    if (not W_abei.is_initialized()){
      W_abei = TAmanager.malloc<MatsT>("vvvo");
    }

  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formF_ae() {
    F_ae("a,e") -= 0.5 * this->T1_("a,m") * F_me("m,e");
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formF_mi() {
    F_mi("m,i") += 0.5 * this->T1_("e,i") * F_me("m,e");
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_mnij() {
    W_mnij("m,n,i,j") += 0.25 * tau("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_abef() {
    W_abef("a,b,e,f") += 0.25 * tau("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_mbej() {
    W_mbej("m,b,e,j") -= 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_mnie() {
    W_mnie("m,n,i,e") = - conj(this->antiSymMoints["vooo"]("e,i,m,n")) + this->T1_("f,i") * conj(this->antiSymMoints["vvoo"]("f,e,m,n"));
  } 

  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_amef() {
    W_amef("a,m,e,f") = conj(this->antiSymMoints["vvvo"]("e,f,a,m")) - this->T1_("a,n") * conj(this->antiSymMoints["vvoo"]("e,f,n,m"));
  } 
  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_mbij() {
    W_mbij("m,b,i,j") = - this->antiSymMoints["vooo"]("b,m,i,j") - F_me("m,e") * this->T2_("b,e,i,j");
    W_mbij("m,b,i,j") += - this->T1_("b,n") * W_mnij("m,n,i,j");
    W_mbij("m,b,i,j") += - 0.5 * conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * tau("e,f,i,j");
    W_mbij("m,b,i,j") += - conj(this->antiSymMoints["vooo"]("e,i,m,n")) * this->T2_("b,e,j,n");
    W_mbij("m,b,i,j") += conj(this->antiSymMoints["vooo"]("e,j,m,n")) * this->T2_("b,e,i,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,j") = - this->antiSymMoints["vovo"]("b,m,e,j") - this->T2_("b,f,n,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    W_mbij("m,b,i,j") += this->T1_("e,i") * tmp("m,b,e,j");
    W_mbij("m,b,i,j") += - this->T1_("e,j") * tmp("m,b,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
  } 
  template <typename MatsT>
  void EOMCCSD<MatsT>::formW_abei() {
    W_abei("a,b,e,i") = this->antiSymMoints["vvvo"]("a,b,e,i") - F_me("m,e") * this->T2_("a,b,m,i");
    W_abei("a,b,e,i") += this->T1_("f,i") * W_abef("a,b,e,f");
    W_abei("a,b,e,i") +=  0.5 * conj(this->antiSymMoints["vooo"]("e,i,m,n")) * tau("a,b,m,n");
    W_abei("a,b,e,i") += conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * this->T2_("a,f,m,i");
    W_abei("a,b,e,i") += -conj(this->antiSymMoints["vvvo"]("e,f,a,m")) * this->T2_("b,f,m,i");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,i") = - this->antiSymMoints["vovo"]("b,m,e,i") - this->T2_("b,f,n,i") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    W_abei("a,b,e,i") += - this->T1_("a,m") * tmp("m,b,e,i");
    W_abei("a,b,e,i") +=  this->T1_("b,m") * tmp("m,a,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
  }  

  template <typename MatsT>
  void EOMCCSD<MatsT>::formEOMIntermediates() {
    TAManager &TAmanager = TAManager::get();

    // If DFCCSD was used for ground state, initiate necessary slices
    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      // Build slices of ERI needed to build intermediates
      this->antiSymMoints["vooo"] = TAmanager.malloc<MatsT>("vooo");
      this->antiSymMoints["vooo"]("a,n,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["boo"]("Q,n,j");
      this->antiSymMoints["vooo"]("a,n,i,j") -= this->antiSymMoints["vooo"]("a,n,j,i");

      this->antiSymMoints["vvoo"] = TAmanager.malloc<MatsT>("vvoo");
      this->antiSymMoints["vvoo"]("a,b,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvoo"]("a,b,i,j") -= this->antiSymMoints["vvoo"]("a,b,j,i");

      this->antiSymMoints["vovo"] = TAmanager.malloc<MatsT>("vovo");
      this->antiSymMoints["vovo"]("a,m,e,i")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["boo"]("Q,m,i");
      this->antiSymMoints["vovo"]("a,m,e,i") -= this->riMoints["bvo"]("Q,a,i") * this->riMoints["bov"]("Q,m,e");

      this->antiSymMoints["vvvo"] = TAmanager.malloc<MatsT>("vvvo");
      this->antiSymMoints["vvvo"]("a,b,e,j")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvvo"]("a,b,e,j") -= this->antiSymMoints["vvvo"]("b,a,e,j");
    }

    formF_ae();
    formF_mi();
    formW_mnij();
    formW_abef();
    formW_mbej();
    formW_mnie();
    formW_amef();
    formW_mbij();
    formW_abei();

    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      // "vvoo" slice still needed in formR2_tilde and buildRightZeroBody, don't clear here
      TAmanager.free("vooo",std::move(this->antiSymMoints["vooo"]));
      TAmanager.free("vovo",std::move(this->antiSymMoints["vovo"]));
      TAmanager.free("vvvo",std::move(this->antiSymMoints["vvvo"]));
      TA::get_default_world().gop.fence();
    }
  } 

  template <typename MatsT>
  void EOMCCSD<MatsT>::formR1_tilde(const TArray &R1, const TArray &R2, TArray &tildeR1) const {
    tildeR1("a,i") = F_ae("a,c") * R1("c,i");
    tildeR1("a,i") += - F_mi("k,i") * R1("a,k");
    tildeR1("a,i") += F_me("k,c") * R2("a,c,i,k");
    tildeR1("a,i") += W_mbej("k,a,c,i") * R1("c,k");
    tildeR1("a,i") += 0.5 * W_amef("a,k,c,d") * R2("c,d,i,k");
    tildeR1("a,i") += - 0.5 * W_mnie("k,l,i,d") * R2("a,d,k,l");
  }  

  template <typename MatsT>
  void EOMCCSD<MatsT>::formR2_tilde(const TArray &R1, const TArray &R2, TArray &tildeR2) const {
    tildeR2("a,b,i,j") = F_ae("b,e") * R2("a,e,i,j");
    tildeR2("a,b,i,j") += - F_ae("a,e") * R2("b,e,i,j");
    tildeR2("a,b,i,j") += - F_mi("k,j") * R2("a,b,i,k");
    tildeR2("a,b,i,j") += F_mi("k,i") * R2("a,b,j,k");
    tildeR2("a,b,i,j") += 0.5 * W_mnij("k,l,i,j") * R2("a,b,k,l");
    tildeR2("a,b,i,j") += 0.5 * W_abef("a,b,e,f") * R2("e,f,i,j");

    tildeR2("a,b,i,j") += W_mbej("k,b,c,j") * R2("a,c,i,k");
    tildeR2("a,b,i,j") += - W_mbej("k,a,c,j") * R2("b,c,i,k");
    tildeR2("a,b,i,j") += - W_mbej("k,b,c,i") * R2("a,c,j,k");
    tildeR2("a,b,i,j") += W_mbej("k,a,c,i") * R2("b,c,j,k");

    tildeR2("a,b,i,j") += W_abei("a,b,c,j") * R1("c,i"); // Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += - W_abei("a,b,c,i") * R1("c,j");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += -W_mbij("k,a,j,i") * R1("b,k");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += W_mbij("k,b,j,i") * R1("a,k");// Sign is different between Tianyuan's and literature

    TAManager &TAmanager = TAManager::get();
    //[Asthana:2019:4102] Break equations A8, A9 and A10 to save computational cost
    TArray tmp1 = TAmanager.malloc<MatsT>("vv");
    tmp1("b,d") = W_amef("b,k,d,c") * R1("c,k");
    tildeR2("a,b,i,j") += tmp1("b,d") * this->T2_("a,d,i,j");
    tildeR2("a,b,i,j") += - tmp1("a,d") * this->T2_("b,d,i,j");
    TArray tmp2 = TAmanager.malloc<MatsT>("oo");
    tmp2("l,j") = W_mnie("l,k,j,c") * R1("c,k");
    tildeR2("a,b,i,j") += - tmp2("l,j") * this->T2_("a,b,i,l");
    tildeR2("a,b,i,j") += tmp2("l,i") * this->T2_("a,b,j,l");

    //reuse tmp2 container
    tmp2("l,j") = - 0.5 * conj(this->antiSymMoints.at("vvoo")("d,c,k,l")) * R2("c,d,j,k");

    tildeR2("a,b,i,j") += tmp2("l,j") * this->T2_("a,b,i,l");
    tildeR2("a,b,i,j") += - tmp2("l,i") * this->T2_("a,b,j,l");
    //reuse tmp1 container
    tmp1("b,e") = 0.5 * conj(this->antiSymMoints.at("vvoo")("e,d,k,l")) * R2("b,d,k,l");

    tildeR2("a,b,i,j") += - tmp1("b,e") * this->T2_("a,e,i,j");

    tildeR2("a,b,i,j") += tmp1("a,e") * this->T2_("b,e,i,j");

    TAmanager.free("vv", std::move(tmp1));
    TAmanager.free("oo", std::move(tmp2));
  }


  template <typename MatsT>
  void EOMCCSD<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {
    const TArray &V1 =  V.get_tensor("OneBody");
    const TArray &V2 =  V.get_tensor("TwoBody");
    TArray &HV1 = HV.get_tensor("OneBody");
    TArray &HV2 = HV.get_tensor("TwoBody");
    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR1_tilde(V1, V2, HV1);
        formR2_tilde(V1, V2, HV2);
        break;
      case EOMCCEigenVecType::LEFT:
        TAManager &TAmanager = TAManager::get();
        TArray G_ae = TAmanager.malloc<MatsT>("vv");
        updateG_ae(V2, G_ae);
        TArray G_mi = TAmanager.malloc<MatsT>("oo");
        updateG_mi(V2, G_mi);
        formL1_tilde(V1, V2, G_ae, G_mi, HV1);
        formL2_tilde(V1, V2, G_ae, G_mi, HV2);
        TAmanager.free("vv", std::move(G_ae));
        TAmanager.free("oo", std::move(G_mi));
        break;
    }
  }


  template <typename MatsT>
  void EOMCCSD<MatsT>::buildRightZeroBody(size_t nVec) {
    // Assumes MBExpansionSet R_ type

    std::shared_ptr<MBExpansionSet<MatsT>> VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);

    for (size_t i = 0; i < nVec; i++) {
      MatsT r0_1 = F_me("i,a").dot(VR->get(i).get_tensor("OneBody")("a,i"));
      TA::get_default_world().gop.fence();    
      MatsT r0_2 = conj(this->antiSymMoints["vvoo"]("a,b,i,j")).dot(VR->get(i).get_tensor("TwoBody")("a,b,i,j"));
      TA::get_default_world().gop.fence();
      if constexpr (std::is_same_v<MatsT, double>) {
        VR->get(i).zeroBody() = (r0_1 + 0.25 * r0_2) / std::real(this->theta[i]);
      } else {
        VR->get(i).zeroBody() = (r0_1 + 0.25 * r0_2) / this->theta[i];
      }
    }
  }

  template <typename MatsT>
  inline double EOMCCSD<MatsT>::signD(size_t a, size_t b, size_t i, size_t j) const {
    if (a == b or i == j) {
      a = 0;
      b = 0;
      i = 0;
      j = 0; 
      return 0.0;
    }
    double sign = 1.0;
    if (a > b) {
      std::swap(a,b);
      sign *= -1.0;
    }
    if (i > j) {
      std::swap(i,j);
      sign *= -1.0;
    }
    return sign;
  }

  template <typename MatsT>
  inline size_t EOMCCSD<MatsT>::toCompoundS(size_t a, size_t i) const {
    if (a >= nV_ or i >= nO_)
      return outOfBound_;
    return a + i * nV_;
  }

  template <typename MatsT>
  inline size_t EOMCCSD<MatsT>::toCompoundD(size_t a, size_t b, size_t i, size_t j) const {
    size_t ab = abIndices_[a][b], ij = ijIndices_[i][j];
    if (ab == outOfBound_ or ij == outOfBound_)
      return outOfBound_;
    return ab + ij * nV2shift_;
  }

  template <typename MatsT>
  inline size_t EOMCCSD<MatsT>::toCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const {
    size_t ai = toCompoundS(a,i), bj = toCompoundS(b,j);
    if (ai == outOfBound_ or bj == outOfBound_)
      return outOfBound_;
    return ai + bj * ldH;
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT>::toCompoundSD(size_t e, size_t m,
                                                                         size_t a, size_t b, size_t i, size_t j, size_t ldH) const {
    size_t em = toCompoundS(e,m), abij = toCompoundD(a,b,i,j);
    if (em == outOfBound_ or abij == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    return std::make_pair(em + abij * ldH, sign);
  }


  template <typename MatsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT>::toCompoundDS(size_t a, size_t b, size_t i, size_t j,
                                                                         size_t e, size_t m, size_t ldH) const {
    size_t em = toCompoundS(e,m), abij = toCompoundD(a,b,i,j);
    if (em == outOfBound_ or abij == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    return std::make_pair(abij + em * ldH, sign);
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT>::toCompoundDD(size_t a, size_t b, size_t i, size_t j,
                                                                         size_t c, size_t d, size_t k, size_t l, size_t ldH) const {
    size_t abij = toCompoundD(a,b,i,j), cdkl = toCompoundD(c,d,k,l);
    if (abij == outOfBound_ or cdkl == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    sign *= signD(c,d,k,l);
    return std::make_pair(abij + cdkl * ldH, sign);
  }

  template <typename MatsT>
  EOMCCSD<MatsT>::~EOMCCSD() {

    //if (theta) CQMemManager::get().free(theta);
    
    TAManager &TAmanager = TAManager::get();

    if (Rho_ij) TAmanager.free("oo", std::move(Rho_ij), true);
    if (Rho_ab) TAmanager.free("vv", std::move(Rho_ab), true);
    if (Rho_ia) TAmanager.free("ov", std::move(Rho_ia), true);
    if (Rho_ai) TAmanager.free("vo", std::move(Rho_ai), true);

  }

//for full_diagonalization
   template <typename MatsT>
  cqmatrix::Matrix<MatsT> EOMCCSD<MatsT>::buildHbar(bool includeGroundState) const {
  
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    cqmatrix::Matrix<MatsT> fullMat(includeGroundState ? this->Hbar_dim + 1 : this->Hbar_dim);
    fullMat.clear();

    MatsT * Hbar = fullMat.pointer();
    size_t ldH = fullMat.nRows();
    size_t nCol = fullMat.nColumns();

    MatsT * Hbar0S = nullptr;
    MatsT * Hbar0D = nullptr;
    MatsT * HbarSS = Hbar;
    MatsT * HbarSD = HbarSS + nOVshift_ * ldH;
    MatsT * HbarDS = HbarSS + nOVshift_;
    MatsT * HbarDD = HbarSD + nOVshift_;
    if (includeGroundState) {
      Hbar0S = Hbar + ldH;
      Hbar0D = Hbar0S + nOVshift_ * ldH;
      HbarSS += 1 + ldH;
      HbarSD += 1 + ldH;
      HbarDS += 1 + ldH;
      HbarDD += 1 + ldH;
    }


    TA::foreach_inplace( F_ae, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") = F_ae("a,c") * R1("c,i");
          size_t a = x[0], c = x[1];

          for (size_t i = 0; i < nO; ++i) {
            size_t idx = toCompoundSS(a,i,c,i,ldH);
            if (isInBound(idx))
              HbarSS[idx] += v;
          }

          // tildeR2("a,b,i,j") = F_ae("b,e") * R2("a,e,i,j");
          // tildeR2("a,b,i,j") += - F_ae("a,e") * R2("b,e,i,j");
          size_t b = x[0], e = x[1];

          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i)
              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = toCompoundDD(a,b,i,j,a,e,i,j,ldH);
                if (isInBound(idx_sgn.first))
                  HbarDD[idx_sgn.first] += idx_sgn.second * v;
              }

        }
    });

    TA::foreach_inplace( F_mi, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") += - F_mi("k,i") * R1("a,k");
          size_t k = x[0], i = x[1];

          for (size_t a = 0; a < nV; ++a) {
            size_t idx = toCompoundSS(a,i,a,k,ldH);
            if (isInBound(idx))
              HbarSS[idx] -= v;
          }

          // tildeR2("a,b,i,j") += - F_mi("k,j") * R2("a,b,i,k");
          // tildeR2("a,b,i,j") += F_mi("k,i") * R2("a,b,j,k");
          size_t j = x[1];

          for (size_t b = 0; b < nV; ++b)
            for (size_t a = 0; a < b; ++a)
              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = toCompoundDD(a,b,i,j,a,b,i,k,ldH);
                if (isInBound(idx_sgn.first))
                  HbarDD[idx_sgn.first] -= idx_sgn.second * v;
              }
        }
    });

    TA::foreach_inplace( F_me, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") += F_me("k,c") * R2("a,c,i,k");
          size_t k = x[0], c = x[1];

          if (includeGroundState) {
            size_t idx = toCompoundS(c,k);
            if (isInBound(idx))
              Hbar0S[idx * ldH] = v;
          }

          for (size_t a = 0; a < nV; ++a)
            for (size_t i = 0; i < nO; ++i) {
              auto idx_sgn = toCompoundSD(a,i,a,c,i,k,ldH);
              if (isInBound(idx_sgn.first))
                HbarSD[idx_sgn.first] += idx_sgn.second * v;
            }
        }
    });

    if (includeGroundState)
      TA::foreach_inplace( this->antiSymMoints["vvoo"], [&](TA::Tensor<MatsT>& tile) {

        const auto& lobound = tile.range().lobound();
        if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
          return;

        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0,0,0,0};
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
                MatsT v = tile[x];

                size_t a = x[0], b = x[1], i = x[2], j = x[3];

                size_t idx = toCompoundD(a,b,i,j);
                if (isInBound(idx))
                  Hbar0D[idx * ldH] = SmartConj(v);

              }

      });

    TA::foreach_inplace( W_mbej, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += W_mbej("k,a,c,i") * R1("c,k");
              size_t k = x[0], a = x[1], c = x[2], i = x[3];
              size_t idx = toCompoundSS(a,i,c,k,ldH);
              if (isInBound(idx))
                HbarSS[idx] += v;

              // tildeR2("a,b,i,j") += W_mbej("k,b,c,j") * R2("a,c,i,k");
              // tildeR2("a,b,i,j") += - W_mbej("k,a,c,j") * R2("b,c,i,k");
              // tildeR2("a,b,i,j") += - W_mbej("k,b,c,i") * R2("a,c,j,k");
              // tildeR2("a,b,i,j") += W_mbej("k,a,c,i") * R2("b,c,j,k");
              for (size_t b = 0; b < nV; ++b)
                for (size_t j = 0; j < nO; ++j) {
                  auto idx_sgn = toCompoundDD(a,b,i,j,b,c,j,k,ldH);
                  if (isInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_mnij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += 0.5 * W_mnij("k,l,i,j") * R2("a,b,k,l");
              size_t k = x[0], l = x[1], i = x[2], j = x[3];

              for (size_t b = 0; b < nV; ++b)
                for (size_t a = 0; a < b; ++a) {
                  auto idx_sgn = toCompoundDD(a,b,i,j,a,b,k,l,ldH);
                  if (isInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_abef, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += 0.5 * W_abef("a,b,e,f") * R2("e,f,i,j");
              size_t a = x[0], b = x[1], e = x[2], f = x[3];

              for (size_t j = 0; j < nO; ++j)
                for (size_t i = 0; i < j; ++i) {
                  auto idx_sgn = toCompoundDD(a,b,i,j,e,f,i,j,ldH);
                  if (isInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_abei, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += W_abei("a,b,c,j") * R1("c,i");
              // tildeR2("a,b,i,j") += - W_abei("a,b,c,i") * R1("c,j");
              size_t a = x[0], b = x[1], c = x[2], j = x[3];

              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = toCompoundDS(a,b,i,j,c,i,ldH);
                if (isInBound(idx_sgn.first))
                  HbarDS[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_mbij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[2] < lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < std::min(x[2], static_cast<std::size_t>(upbound[3])); ++x[3]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += -W_mbij("k,a,j,i") * R1("b,k");
              // tildeR2("a,b,i,j") += W_mbij("k,b,j,i") * R1("a,k");
              size_t k = x[0], b = x[1], j = x[2], i = x[3];

              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = toCompoundDS(a,b,i,j,a,k,ldH);
                if (isInBound(idx_sgn.first))
                  HbarDS[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_amef, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += 0.5 * W_amef("a,k,c,d") * R2("c,d,i,k");
              size_t a = x[0], k = x[1], c = x[2], d = x[3];

              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = toCompoundSD(a,i,c,d,i,k,ldH);
                if (isInBound(idx_sgn.first))
                  HbarSD[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_mnie, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += - 0.5 * W_mnie("k,l,i,d") * R2("a,d,k,l");
              size_t k = x[0], l = x[1], i = x[2], d = x[3];

              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = toCompoundSD(a,i,a,d,k,l,ldH);
                if (isInBound(idx_sgn.first))
                  HbarSD[idx_sgn.first] -= idx_sgn.second * v;
              }

            }
    });

    TArray WT_ckabij = TAmanager.malloc<MatsT>("vovvoo");
    WT_ckabij("c,k,a,b,i,j") = 0.5 * W_amef("b,k,d,c") * T2_("a,d,i,j");
    WT_ckabij("c,k,a,b,i,j") -= 0.5 * W_mnie("l,k,j,c") * T2_("a,b,i,l");
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( WT_ckabij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") += W_amef("b,k,d,c") * T2_("a,d,i,j") * R1("c,k");
                  // tildeR2("a,b,i,j") -= W_amef("a,k,d,c") * T2_("b,d,i,j") * R1("c,k");
                  // tildeR2("a,b,i,j") -= W_mnie("l,k,j,c") * T2_("a,b,i,l") * R1("c,k");
                  // tildeR2("a,b,i,j") += W_mnie("l,k,i,c") * T2_("a,b,j,l") * R1("c,k");
                  size_t c = x[0], k = x[1], a = x[2], b = x[3], i = x[4], j = x[5];

                  auto idx_sgn = toCompoundDS(a,b,i,j,c,k,ldH);
                  if (isInBound(idx_sgn.first))
                    HbarDS[idx_sgn.first] += idx_sgn.second * v;

                }
    });
    TAmanager.free("vovvoo", std::move(WT_ckabij), true);

    TArray TV_abidck = TAmanager.malloc<MatsT>("vvovvo");
    TV_abidck("a,b,i,d,c,k") = T2_("a,b,i,l") * conj(this->antiSymMoints.at("vvoo")("d,c,k,l"));
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( TV_abidck, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[4] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < std::min(x[3], static_cast<std::size_t>(upbound[4])); ++x[4])
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") -= 0.5 * T2_("a,b,i,l") * conj(V("d,c,k,l")) * R2("c,d,j,k");
                  // tildeR2("a,b,i,j") += 0.5 * T2_("a,b,j,l") * conj(V("d,c,k,l")) * R2("c,d,i,k");
                  size_t a = x[0], b = x[1], i = x[2], d = x[3], c = x[4], k = x[5];

                  for (size_t j = 0; j < nO; ++j) {
                    auto idx_sgn = toCompoundDD(a,b,i,j,c,d,j,k,ldH);
                    if (isInBound(idx_sgn.first))
                      HbarDD[idx_sgn.first] -= idx_sgn.second * v;
                  }

                }
    });
    TAmanager.free("vvovvo", std::move(TV_abidck), true);

    TArray TV_aijdkl = TAmanager.malloc<MatsT>("voovoo");
    TV_aijdkl("a,i,j,d,k,l") = T2_("e,a,i,j") * conj(this->antiSymMoints.at("vvoo")("e,d,k,l"));
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( TV_aijdkl, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[1] > lobound[2] or lobound[4] > lobound[5])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
          for(x[1] = lobound[1]; x[1] < std::min(x[2], static_cast<std::size_t>(upbound[1])); ++x[1])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5])
                for(x[4] = lobound[4]; x[4] < std::min(x[5], static_cast<std::size_t>(upbound[4])); ++x[4]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") += 0.5 * T2_("e,a,i,j") * conj(V("e,d,k,l")) * R2("b,d,k,l");
                  // tildeR2("a,b,i,j") -= 0.5 * T2_("e,b,i,j") * conj(V("e,d,k,l")) * R2("a,d,k,l");
                  size_t a = x[0], i = x[1], j = x[2], d = x[3], k = x[4], l = x[5];

                  for (size_t b = 0; b < nV; ++b) {
                    auto idx_sgn = toCompoundDD(a,b,i,j,b,d,k,l,ldH);
                    if (isInBound(idx_sgn.first))
                      HbarDD[idx_sgn.first] += idx_sgn.second * v;
                  }

                }
    });
    TAmanager.free("voovoo", std::move(TV_aijdkl), true);

    TA::get_default_world().gop.fence();

    MatsT *Hbar_copy = CQMemManager::get().malloc<MatsT>(ldH*nCol);
    std::copy_n(Hbar, ldH*nCol, Hbar_copy);
    std::fill_n(Hbar, ldH*nCol, MatsT(0.0));
    MPIAllReduce(Hbar_copy, ldH*nCol, Hbar, MPI_COMM_WORLD);
    CQMemManager::get().free(Hbar_copy);

    TA::get_default_world().gop.fence();

    return fullMat;

  }
  
  // for CCS_guess 
  template <typename MatsT>
  void EOMCCSD<MatsT>::fillGuess(MatsT *guess_vec, size_t n_vec) const{
    
    if (n_vec > nOVshift_) {
      CErr("EOMCCSD: asking for more roots than the single excitatin space dimension.");
    } 

  
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    cqmatrix::Matrix<MatsT> fullMat(nOVshift_);
    fullMat.clear();

    MatsT * Hbar = fullMat.pointer();
    size_t ldH = fullMat.nRows();
    size_t nCol = fullMat.nColumns();

    MatsT * HbarSS = Hbar;


    TA::foreach_inplace( F_ae, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") = F_ae("a,c") * R1("c,i");
          size_t a = x[0], c = x[1];

          for (size_t i = 0; i < nO; ++i) {
            size_t idx = toCompoundSS(a,i,c,i,ldH);
            if (isInBound(idx))
              HbarSS[idx] += v;
          }

        }
    });

    TA::foreach_inplace( F_mi, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") += - F_mi("k,i") * R1("a,k");
          size_t k = x[0], i = x[1];

          for (size_t a = 0; a < nV; ++a) {
            size_t idx = toCompoundSS(a,i,a,k,ldH);
            if (isInBound(idx))
              HbarSS[idx] -= v;
          }

        }
    });


    TA::foreach_inplace( W_mbej, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += W_mbej("k,a,c,i") * R1("c,k");
              size_t k = x[0], a = x[1], c = x[2], i = x[3];
              size_t idx = toCompoundSS(a,i,c,k,ldH);
              if (isInBound(idx))
                HbarSS[idx] += v;


            }
    });


    TA::get_default_world().gop.fence();

    MatsT *Hbar_copy = CQMemManager::get().malloc<MatsT>(ldH*nCol);
    std::copy_n(Hbar, ldH*nCol, Hbar_copy);
    std::fill_n(Hbar, ldH*nCol, MatsT(0.0));
    MPIAllReduce(Hbar_copy, ldH*nCol, Hbar, MPI_COMM_WORLD);
    CQMemManager::get().free(Hbar_copy);

    TA::get_default_world().gop.fence();


    dcomplex * theta = CQMemManager::get().malloc<dcomplex>(nOVshift_);
    MatsT * VR    = CQMemManager::get().malloc<MatsT>(nOVshift_ * nOVshift_);
    MatsT * dummy = nullptr;
    if (MPIRank() == 0) GeneralEigen('N', 'V', nOVshift_, fullMat.pointer(), nOVshift_, theta, dummy, 1, VR, nOVshift_);
    std::copy_n(VR, nOVshift_ * n_vec, guess_vec);
    CQMemManager::get().free(theta);
    CQMemManager::get().free(VR);

  }


   /// for full_diagonaization routine
  template <typename MatsT>
  void EOMCCSD<MatsT>::buildDiag(MatsT * diag, const std::vector<double> &eps) const {

//    for (size_t a = 0; a < NV; ++a)
//      for (size_t i = 0; i < NO; ++i) {
//        // F_ae[a][a] - F_mi[i][i] + W_mbej[i][a][a][i]
//      }
//
//    for (size_t a = 0; a < NV; ++a)
//      for (size_t b = 0; b < a; ++b)
//        for (size_t i = 0; i < NO; ++i)
//          for (size_t j = 0; j < i; ++j) {
//            //   F_ae[b][b] + F_ae[a][a] - F_mi[j][j] - F_mi[i][i]
//            // + W_mnij[i][j][i][j] + W_abef[a][b][a][b]
//            // + W_mbej[j][b][b][j] + W_mbej[j][a][a][j]
//            // + W_mbej[i][b][b][i] + W_mbej[i][a][a][i]
//            // - sum{l} T2[a][b][i][l] * conj(V[a][b][i][l])
//            // - sum{l} T2[a][b][l][j] * conj(V[a][b][l][j])
//            // - sum{e} T2[a][e][i][j] * conj(V[a][e][i][j])
//            // - sum{e} T2[e][b][i][j] * conj(V[e][b][i][j])
//          }
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));
    MatsT * diag2 = diag + nOVshift_;

    TA::foreach_inplace( F_ae, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t a = x[0], b = x[0];
        for (size_t i = 0; i < nO; ++i)
          diag[toCompoundS(a,i)] += v;

        for (size_t b = a + 1; b < nV; ++b)
          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }

        for (size_t a = 0; a < b; ++a)
          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }
      }

    });

    TA::foreach_inplace( F_mi, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t i = x[0], j = x[0];
        for (size_t a = 0; a < nV; ++a)
          diag[toCompoundS(a,i)] -= v;

        for (size_t b = 0; b < nV; ++b)
          for (size_t a = 0; a < b; ++a)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[toCompoundD(a, b, i, j)] -= v;
            }

        for (size_t b = 0; b < nV; ++b)
          for (size_t a = 0; a < b; ++a)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] -= v;
            }
      }

    });

    TA::foreach_inplace( W_mbej, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
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
          diag[toCompoundS(a,i)] += v;

          for (size_t b = a + 1; b < nV; ++b)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }

          for (size_t b = a + 1; b < nV; ++b)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }

          for (size_t a = 0; a < b; ++a)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }

          for (size_t a = 0; a < b; ++a)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_mnij, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
      return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], j = x[1];

          for (size_t b = 0; b < nV; ++b)
            for (size_t a = 0; a < b; ++a) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_abef, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t a = x[0], b = x[1];

          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[toCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    
    TA::foreach_inplace(T2_, this->antiSymMoints.at("vvoo"),
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

              if (i < j) {

                for (size_t a = 0; a < b; ++a)
                  diag2[toCompoundD(a, b, i, j)] -= v;

                for (size_t b = a + 1; b < nV; ++b)
                  diag2[toCompoundD(a, b, i, j)] -= v;

              }

              if (a < b) {

                for (size_t i = 0; i < j; ++i)
                  diag2[toCompoundD(a, b, i, j)] -= v;

                for (size_t j = i + 1; j < nO; ++j)
                  diag2[toCompoundD(a, b, i, j)] -= v;

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

  }


  template <typename MatsT>
  typename Davidson<MatsT>::VecsGen_t EOMCCSD<MatsT>::EmptyDavidsonVectorBuilder(){
      // Algorithm with implicit Hbar matrix
      typename Davidson<MatsT>::VecsGen_t vecsGenEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSet<MatsT>>(this->tensor_builder_, nVec, this->savFile_);
        }; // implicit vecsGenerator

        return vecsGenEOM;
      }
      fullMat = nullptr;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or 
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        std::cout << "  *** Start building the full matrix for explicit diagonalization ***" << std::endl;

        auto beginBuildHbar = tick();
        fullMat = std::make_shared<cqmatrix::Matrix<MatsT>>(buildHbar(false));
        std::cout << "    * Build Hbar spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginBuildHbar) << " s." << std::endl;
      }

      // Algorithm for debug, comparing implicit and explicit
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        typename Davidson<MatsT>::VecsGen_t vecsGenRaw = [/*&Hbar_dim,*/ this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<RawVectors<MatsT>>(
              MPI_COMM_WORLD, this->Hbar_dim, nVec
              );
        };

        typename Davidson<MatsT>::VecsGen_t vecsGenDebug =
            [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSetDebug<MatsT>>(
              this->tensor_builder_, nVec, this->savFile_, MPI_COMM_WORLD
              );
        };
        return vecsGenDebug;
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return vecsGenEOM;
      return vecsGenEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMCCSD<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
//      typename Davidson<MatsT>::LinearTrans_t funcEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
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
            const MBExpansion<MatsT> &Vi = V_ptr->get(i + Vshift);
            MBExpansion<MatsT> &AVi = AV_ptr->get(i + AVshift);
            //buildSigma(Vi.get_tensor("OneBody"), Vi.get_tensor("TwoBody"), AVi.get_tensor("OneBody"), AVi.get_tensor("TwoBody"), eigenVecType);
            buildSigma(Vi, AVi, eigenVecType);
            AVi.enforceSymmetry();
          }
//V.print(std::cout, "R-imp", 0, nVec);
//AV.print(std::cout, "H*R-imp", 0, nVec);
        }; // implicit sigmaBuilder
      }
//      typename Davidson<MatsT>::LinearTrans_t funcRaw;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        this->funcRaw = [this, &eigenVecType]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {
          ROOT_ONLY(MPI_COMM_WORLD);

          size_t N = this->Hbar_dim;
          auto V_ptr = tryGetRawVectorsPointer(V);
          auto AV_ptr = tryGetRawVectorsPointer(AV);

          switch(eigenVecType) {
            case EOMCCEigenVecType::RIGHT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
                         N,nVec,N,MatsT(1.),fullMat->pointer(),N,V_ptr,N,MatsT(0.),AV_ptr,N);
              break;
            case EOMCCEigenVecType::LEFT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                         nVec,N,N,MatsT(1.),V_ptr,N,fullMat->pointer(),N,MatsT(0.),AV_ptr,nVec);
              IMatCopy('T', nVec, N, 1.0, AV_ptr,nVec, N);
              break;
          }
//V.print(std::cout, "R-exp", 0, nVec);
//AV.print(std::cout, "H*R-exp", 0, nVec);
        }; // explicit sigmaBuilder
      }

      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
 //       typename Davidson<MatsT>::LinearTrans_t
        this->funcDebug = [this]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSetDebug<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<MatsT> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<MatsT> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<MatsT> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<MatsT> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::funcDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          this->funcEOM(nVec, V_EOM, AV_EOM);
          this->funcRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::funcDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };
        return this->funcDebug;

      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) return this->funcEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return this->funcRaw;

      return this->funcEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMCCSD<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){
      
      double PCsmall = this->eomSettings.davidson_preCond_small;

//      typename Davidson<MatsT>::LinearTrans_t PCEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
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
                  denom = curEigI - eomDiag[toCompoundS(x[0], x[1])];
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
                  if (x[0] == x[1])
                    continue;
                  size_t a = x[0], b = x[1];
                  for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                    for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                      if (x[2] == x[3])
                        continue;
                      size_t i = x[2], j = x[3];
                      signD(a,b,i,j);
                      denom = curEigI - diagD[toCompoundD(a, b, i, j)];
                      if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                    }
                }
            });
            TA::get_default_world().gop.fence();

            curB.enforceSymmetry();
          }
        }; // implicit preConditioner

      }
//      typename Davidson<MatsT>::LinearTrans_t PCRaw;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or 
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {

        this->PCRaw = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          ROOT_ONLY(MPI_COMM_WORLD);

          //            prettyPrintSmart(std::cout, "eomDiag", eomDiag, Hbar_dim, 1, Hbar_dim);

          MatsT denom = 0.0, nom = 0.0;
          const MatsT *Vptr = tryGetRawVectorsPointer(V);
          MatsT *AVptr = tryGetRawVectorsPointer(AV);

          // Scale by inverse diagonals
          for (size_t i = 0; i < nVec; i++) {
            MatsT curEigI = 0.0;
            if constexpr (std::is_same_v<MatsT, double>) {
              curEigI = curEig[i].real();
            } else {
              curEigI = curEig[i];
            }
            for(auto k = 0ul; k < this->Hbar_dim; k++ ) {
              nom = Vptr[k];
              denom = curEigI - eomDiag[k];
              if (std::abs(denom) >= PCsmall)
                AVptr[k] = nom / denom;
              else
                AVptr[k] = nom;
            }
            Vptr += this->Hbar_dim;
            AVptr += this->Hbar_dim;
          }

        }; // explicit preConditioner
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
//        typename Davidson<MatsT>::LinearTrans_t
        this->PCDebug = [this]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSetDebug<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<MatsT> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<MatsT> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<MatsT> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<MatsT> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::PCDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          this->PCEOM(nVec, V_EOM, AV_EOM);
          this->PCRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::PCDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };
        return this->PCDebug;
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) return this->PCEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return this->PCRaw;

      return this->PCEOM;
  }
}
