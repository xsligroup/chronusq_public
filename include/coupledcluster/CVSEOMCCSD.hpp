/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <util/math.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <functional>
#include <util/timer.hpp>

namespace ChronusQ{
  template <typename MatsT>
  CVSEOMCCSD<MatsT>::CVSEOMCCSD(const SafeFile &savFile,
                                      CCIntermediates<MatsT> &intermediates,
                                      const EOMSettings &eomSettings,
                                      const CoupledClusterSettings &ccSettings):
    EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
    vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
    hLabel_(intermediates.hLabel), lLabel_(intermediates.lLabel),
    rLabel_(intermediates.rLabel), cLabel_(intermediates.cLabel),
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
    reuse_tmps_(intermediates.sigmaOps),
    tmps_(intermediates.tempOps),
    T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")){

    TAManager &TAmanager = TAManager::get();

    // without L, we don't need D
    if (eomSettings.oscillator_strength == false) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    assignCVSIndices();

    size_t offset = 0;
    this->tensor_builder_.push_back(std::string({intermediates.rLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.cLabel}));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->Hbar_dimension_offsets.emplace("ai", offset);
    offset += nOVshift_;
    this->tensor_builder_.push_back(std::string({intermediates.rLabel,intermediates.rLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.cLabel,intermediates.cLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody_core"));
    this->Hbar_dimension_offsets.emplace("abij", offset);
    offset += nV2shift_ * nCVSOCore_ * (nCVSOCore_ - 1) / 2 ;
    this->tensor_builder_.push_back(std::string({intermediates.rLabel,intermediates.rLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.cLabel,intermediates.hLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody_val"));
    this->Hbar_dimension_offsets.emplace("abiJ", offset);

    b.cc = TAmanager.toBlockRange("cc");
    b.ch = TAmanager.toBlockRange("ch");
    b.cr = TAmanager.toBlockRange("cr");
    b.hc = TAmanager.toBlockRange("hc");
    b.hh = TAmanager.toBlockRange("hh");
    b.hr = TAmanager.toBlockRange("hr");
    b.oc = TAmanager.toBlockRange("oc");
    b.oh = TAmanager.toBlockRange("oh");
    b.rc = TAmanager.toBlockRange("rc");
    b.rr = TAmanager.toBlockRange("rr");
    b.co = TAmanager.toBlockRange("co");
    b.ho = TAmanager.toBlockRange("ho");
    b.rh = TAmanager.toBlockRange("rh");
    b.cccc = TAmanager.toBlockRange("cccc");
    b.ccch = TAmanager.toBlockRange("ccch");
    b.cchh = TAmanager.toBlockRange("cchh");
    b.cccr = TAmanager.toBlockRange("cccr");
    b.chcc = TAmanager.toBlockRange("chcc");
    b.chch = TAmanager.toBlockRange("chch");
    b.chcr = TAmanager.toBlockRange("chcr");
    b.crcc = TAmanager.toBlockRange("crcc");
    b.crch = TAmanager.toBlockRange("crch");
    b.crrc = TAmanager.toBlockRange("crrc");
    b.crrh = TAmanager.toBlockRange("crrh");
    b.hrrc = TAmanager.toBlockRange("hrrc");
    b.hrrh = TAmanager.toBlockRange("hrrh");
    b.occr = TAmanager.toBlockRange("occr");
    b.ochr = TAmanager.toBlockRange("ochr");
    b.rcrr = TAmanager.toBlockRange("rcrr");
    b.rhrr = TAmanager.toBlockRange("rhrr");
    b.rrcc = TAmanager.toBlockRange("rrcc");
    b.rrch = TAmanager.toBlockRange("rrch");
    b.rrco = TAmanager.toBlockRange("rrco");
    b.rrho = TAmanager.toBlockRange("rrho");
    b.rrhh = TAmanager.toBlockRange("rrhh");
    b.rroc = TAmanager.toBlockRange("rroc");
    b.rroh = TAmanager.toBlockRange("rroh");
    b.rrrc = TAmanager.toBlockRange("rrrc");
    b.rrrh = TAmanager.toBlockRange("rrrh");
    b.rrrr = TAmanager.toBlockRange("rrrr");




    }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::assignCVSIndices() {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    nCVSOActive_ = nO;
    nCVSOCore_ = this->eomSettings.cvs_core.size();
    nCVSOValance_ = nCVSOActive_ - nCVSOCore_;

    nCVSVActive_ = nV;
    nCVSVContinuum_ = this->eomSettings.external_virtual.size();
    if (nCVSVContinuum_ != nCVSVActive_) {
      CErr("CVSEOMCCSD: CVS does not support virtual valence space.");
    }

    nOVshift_ = nCVSOCore_ * nCVSVContinuum_;
    nO2shift_ = nCVSOCore_ * (nCVSOCore_ - 1) / 2 + nCVSOCore_ * nCVSOValance_;
    nV2shift_ = nCVSVContinuum_ * (nCVSVContinuum_ - 1) / 2 ;

    this->Hbar_dim = nOVshift_ + nO2shift_ * nV2shift_;
    CVSoutOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);

    CVSabIndices_.clear();
    CVSabIndices_.resize(nV, std::vector<size_t>(nV, CVSoutOfBound_));

    size_t idx = 0;
    for (size_t b = 0; b < nCVSVContinuum_; b++) {
      for (size_t a = 0; a < b; a++) {
        CVSabIndices_[a][b] = idx;
        CVSabIndices_[b][a] = idx++;
      }
    }


    CVSijIndices_.clear();
    CVSijIndices_.resize(nO, std::vector<size_t>(nO, CVSoutOfBound_));

    idx = 0;
    for (size_t j = 0; j < nCVSOCore_; j++) {
      for (size_t i = 0; i < j; i++) {
        CVSijIndices_[i][j] = idx;
        CVSijIndices_[j][i] = idx++;
      }
    }
    idx = 0;
    for (size_t j = nCVSOCore_; j < nCVSOActive_; j++) {
      for (size_t i = 0; i < nCVSOCore_; i++) {
        CVSijIndices_[i][j] = idx;
        CVSijIndices_[j][i] = idx++;
      }
    }

  }

  //template <typename MatsT>
  //inline double CVSEOMCCSD<MatsT>::signD(size_t a, size_t b, size_t i, size_t j) const {
  //  if (a == b or i == j) {
  //    a = 0;
  //    b = 0;
  //    i = 0;
  //    j = 0; 
  //    return 0.0;
  //  }
  //  double sign = 1.0;
  //  if (a > b) {
  //    std::swap(a,b);
  //    sign *= -1.0;
  //  }
  //  if (i > j) {
  //    std::swap(i,j);
  //    sign *= -1.0;
  //  }
  //  return sign;
  //}

  // in this particular implimentation, a and i must include frozen space, as they are the index from HF calculations un-ordered by space
  template <typename MatsT>
  inline size_t CVSEOMCCSD<MatsT>::toCompoundS(size_t a, size_t i) const {
    if (a >= nCVSVContinuum_ or i >= nCVSOCore_)
      return CVSoutOfBound_;
    return a + i * nCVSVContinuum_;
  }

  // in this particular implimentation, a and i must include frozen space, as they are the index from HF calculations un-ordered by space
  template <typename MatsT>
  inline size_t CVSEOMCCSD<MatsT>::toCompoundD(size_t a, size_t b, size_t i, size_t j) const {
    size_t ab = CVSabIndices_[a][b], ij = CVSijIndices_[i][j];
    if (ab == CVSoutOfBound_ or ij == CVSoutOfBound_)
      return CVSoutOfBound_;
    return ab + ij * nV2shift_;
  }

  template <typename MatsT>
  inline size_t CVSEOMCCSD<MatsT>::toCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const {
    size_t ai = toCompoundS(a,i), bj = toCompoundS(b,j);
    if (ai == CVSoutOfBound_ or bj == CVSoutOfBound_)
      return CVSoutOfBound_;
    return ai + bj * ldH;
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::initializeEOMCC() {
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
  void CVSEOMCCSD<MatsT>::formF_ae() {
    F_ae("a,e") -= 0.5 * this->T1_("a,m") * F_me("m,e");
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formF_mi() {
    F_mi("m,i") += 0.5 * this->T1_("e,i") * F_me("m,e");
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_mnij() {
    W_mnij("m,n,i,j") += 0.25 * tau("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_abef() {
    W_abef("a,b,e,f") += 0.25 * tau("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_mbej() {
    W_mbej("m,b,e,j") -= 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_mnie() {
    W_mnie("m,n,i,e") = - conj(this->antiSymMoints["vooo"]("e,i,m,n")) + this->T1_("f,i") * conj(this->antiSymMoints["vvoo"]("f,e,m,n"));
  } 

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_amef() {
    W_amef("a,m,e,f") = conj(this->antiSymMoints["vvvo"]("e,f,a,m")) - this->T1_("a,n") * conj(this->antiSymMoints["vvoo"]("e,f,n,m"));
  } 
  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formW_mbij() {
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
  void CVSEOMCCSD<MatsT>::formW_abei() {
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
  void CVSEOMCCSD<MatsT>::formEOMIntermediates() {
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


//  //for full_diagonalization
//   template <typename MatsT>
//  cqmatrix::Matrix<MatsT> CVSEOMCCSD<MatsT>::buildHbar(bool includeGroundState) const{} 

  /// for full_diagonaization and Davidson guess
  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::buildDiag(MatsT * diag, const std::vector<double> &eps) const {

    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));

    TAManager &TAmanager = TAManager::get();
    size_t nV = nCVSVContinuum_;
    size_t nC = nCVSOCore_;
    size_t nVa = nCVSOValance_;
    size_t nO = TAmanager.getRange(oLabel_).extent();

    MatsT * diag2 = diag + nOVshift_;
    MatsT * diag3 = diag + nOVshift_ + nC*(nC-1)*nV*(nV-1)/4;


    for (auto a = 0; a < nV; a++){
      for (auto i = 0; i < nC; i++){
        diag[toCompoundS(a,i)] = eps[a+nO] - eps[i];
      }
    }

    for (auto b = 0; b < nV; b++){
      for (auto a = 0; a < b; a++){
        for (auto j = 0; j < nC; j++){
          for (auto i = 0; i < j; i++){
            diag2[toCompoundD(a,b,i,j)] = eps[a+nO] + eps[b+nO] - eps[i] - eps[j];
          }
        }
      }
    }

    for (auto b = 0; b < nV; b++){
      for (auto a = 0; a < b; a++){
        for (auto j = nC; j < nO; j++){
          for (auto i = 0; i < nC; i++){
            diag3[toCompoundD(a,b,i,j)] = eps[a+nO] + eps[b+nO] - eps[i] - eps[j];
          }
        }
      }
    }

    TA::get_default_world().gop.fence();
    MatsT *diag_copy = CQMemManager::get().malloc<MatsT>(this->Hbar_dim);
    std::copy_n(diag, this->Hbar_dim, diag_copy);
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));
    MPIAllReduce(diag_copy, this->Hbar_dim, diag, MPI_COMM_WORLD);
    CQMemManager::get().free(diag_copy);

  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {

    const TArray &V1 = V.get_tensor("OneBody");
    TArray &HV1 = HV.get_tensor("OneBody");
    const TArray &V2 = V.get_tensor("TwoBody_core");
    TArray &HV2 = HV.get_tensor("TwoBody_core");
    const TArray &V2Val = V.get_tensor("TwoBody_val");
    TArray &HV2Val = HV.get_tensor("TwoBody_val");

    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR1_tilde(V1, V2, V2Val, HV1);
        formR2_tilde(V1, V2, V2Val, HV2);
        formR2Val_tilde(V1, V2, V2Val, HV2Val);
        break;
      case EOMCCEigenVecType::LEFT:
        updateG_ae(V2, V2Val, G_ae);
        updateG_mi(V2, V2Val, G_mi);
        formL1_tilde(V1,V2,V2Val,G_ae,G_mi,HV1);
        formL2_tilde(V1,V2,V2Val,G_ae,G_mi,HV2);
        formL2Val_tilde(V1,V2,V2Val,G_ae,G_mi,HV2Val);
	      break;
    }
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::buildRightZeroBody(size_t nVec) {

    std::shared_ptr<MBExpansionSet<MatsT>> VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);

    for (size_t i = 0; i < nVec; i++) {
      MatsT r0_1 = F_me("i,a").block(b.cr).dot(VR->get(i).get_tensor("OneBody")("a,i"));
      MatsT r0_2 = conj(this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrcc)).dot(VR->get(i).get_tensor("TwoBody_core")("a,b,i,j"));
      MatsT r0_3 = conj(this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrch)).dot(VR->get(i).get_tensor("TwoBody_val")("a,b,i,j"));
      if constexpr (std::is_same_v<MatsT, double>) {
        VR->get(i).zeroBody() = (r0_1 + 0.25 * r0_2 + 0.5 * r0_3) / std::real(this->theta[i]);
      } else {
        VR->get(i).zeroBody() = (r0_1 + 0.25 * r0_2 + 0.5 * r0_3) / this->theta[i];
      }
    }
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formR1_tilde(const TArray &R1, const TArray &R2, const TArray &R2Val, TArray &tildeR1) const {

    tildeR1("a,i") = F_ae("a,c").block(b.rr) * R1("c,i");
    tildeR1("a,i") -= R1("a,j") * F_mi("j,i").block(b.cc);
    tildeR1("a,i") += R1("b,j") * W_mbej("j,a,b,i").block(b.crrc);

    tildeR1("a,i") += R2("a,b,i,j") * F_me("j,b").block(b.cr);
    tildeR1("a,i") += R2Val("a,b,i,j") * F_me("j,b").block(b.hr);

    tildeR1("a,i") -= 0.5 * R2("a,b,j,k") * W_mnie("j,k,i,b").block(b.cccr);
    tildeR1("a,i") -= R2Val("a,b,j,k") * W_mnie("j,k,i,b").block(b.chcr); // factor of two applied for second uirtual block

    tildeR1("a,i") += 0.5 * R2("b,c,i,j") * W_amef("a,j,b,c").block(b.rcrr);
    tildeR1("a,i") += 0.5 * R2Val("b,c,i,j") * W_amef("a,j,b,c").block(b.rhrr);

  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formR2_tilde(const TArray &R1, const TArray &R2, const TArray &R2Val, TArray &tildeR2) const {
    TAManager &TAmanager = TAManager::get();

    TArray tmp = TAmanager.malloc<dcomplex>("rrcc");

    // W_mbij * R2 
    tmp("b,a,i,j") = R1("b,m") * W_mbij("m,a,i,j").block(b.crcc);
    tildeR2("a,b,i,j") = tmp("b,a,i,j");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");

    //H_abk;ijc ---2
    TArray tmp_co = TAmanager.malloc<dcomplex>("co");
    tmp_co("i,m") = R1("e,n") * W_mnie("m,n,i,e").block(b.occr);
    //tmp_co("i,m") = R1("e,n") * W_mnie("m,n,i,e").block(b.occr).block(TAmanager.toRange("occr"));
    tmp("a,b,j,i") = tmp_co("i,m") * this->T2_("a,b,j,m").block(b.rrco);
    tildeR2("a,b,i,j") += tmp("a,b,j,i");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");

    // W_abei * R2 
    tmp("a,b,j,i") = R1("e,j") * W_abei("a,b,e,i").block(b.rrrc);
    tildeR2("a,b,i,j") -= tmp("a,b,j,i");
    tildeR2("a,b,i,j") += tmp("a,b,i,j");

    //H_abk;ijc ---1
    TArray tmp_uv = TAmanager.malloc<dcomplex>("rr");
    tmp_uv("a,e") = R1("f,m") * W_amef("a,m,e,f").block(b.rcrr);
    tmp("a,b,i,j") = tmp_uv("a,e") * this->T2_("b,e,i,j").block(b.rrcc);
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");
    tildeR2("a,b,i,j") += tmp("b,a,i,j");

    // F_ae * R2 
    tmp("a,b,j,i") = R2("a,b,j,m") * F_mi("m,i").block(b.cc);
    tildeR2("a,b,i,j") += tmp("a,b,j,i");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");
    tmp("a,b,j,i") = R2Val("a,b,j,m") * F_mi("m,i").block(b.hc);
    tildeR2("a,b,i,j") += tmp("a,b,j,i");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");

    // H_abk,icd * R2
    tmp_co("i,m") = R2("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rroc));
    tmp_co("i,m") += R2Val("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rroh));
    tmp("a,b,j,i") = 0.5 * tmp_co("i,m") * this->T2_("a,b,j,m").block(b.rrco);
    tildeR2("a,b,i,j") += tmp("a,b,j,i");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");

    // F_mi * R2
    tmp("a,b,i,j") = R2("a,e,i,j") * F_ae("b,e").block(b.rr);
    tildeR2("a,b,i,j") += tmp("a,b,i,j");
    tildeR2("a,b,i,j") -= tmp("b,a,i,j");

    // H_akl,ijd * R2
    tmp_uv("a,e") = R2("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rrcc));
    tmp_uv("a,e") += 2.0 * R2Val("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rrch)); // factor of two applied for remaining virtual contribution
    tmp("a,b,i,j") = 0.5 * tmp_uv("a,e") * this->T2_("b,e,i,j").block(b.rrcc);
    tildeR2("a,b,i,j") += tmp("a,b,i,j");
    tildeR2("a,b,i,j") -= tmp("b,a,i,j");

    // W_mbej * R2
    tmp("a,b,j,i") = R2("a,e,m,j") * W_mbej("m,b,e,i").block(b.crrc);
    tildeR2("a,b,i,j") += tmp("a,b,j,i");
    tildeR2("a,b,i,j") -= tmp("b,a,j,i");
    tildeR2("a,b,i,j") -= tmp("a,b,i,j");
    tildeR2("a,b,i,j") += tmp("b,a,i,j");
    tmp("a,b,j,i") = R2Val("a,e,j,m") * W_mbej("m,b,e,i").block(b.hrrc);
    tildeR2("a,b,i,j") -= tmp("a,b,j,i");
    tildeR2("a,b,i,j") += tmp("b,a,j,i");
    tildeR2("a,b,i,j") += tmp("a,b,i,j");
    tildeR2("a,b,i,j") -= tmp("b,a,i,j");

    // W_mnij * R2
    tildeR2("a,b,i,j") += 0.5 * R2("a,b,m,n") * W_mnij("m,n,i,j").block(b.cccc);
    tildeR2("a,b,i,j") += R2Val("a,b,m,n") * W_mnij("m,n,i,j").block(b.chcc); // factor of two applied for remaining uirtual contribution

    // W_abef * R2
    tildeR2("a,b,i,j") += 0.5 * R2("e,f,i,j") * W_abef("a,b,e,f").block(b.rrrr);
    
    TAmanager.free("rrcc", std::move(tmp));
    TAmanager.free("rr", std::move(tmp_uv));
    TAmanager.free("co", std::move(tmp_co));
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formR2Val_tilde(const TArray &R1, const TArray &R2, const TArray &R2Val, TArray &tildeR2Val) const {
    TAManager &TAmanager = TAManager::get();

    TArray tmp = TAmanager.malloc<dcomplex>("rrch");

    // W_mbij * R2
    tmp("a,b,i,j") = R1("b,m") * W_mbij("m,a,i,j").block(b.crch);
    tildeR2Val("a,b,i,j") = tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") -= tmp("b,a,i,j");

    //H_abk;ijc ---2
    TArray tmp_co = TAmanager.malloc<dcomplex>("co");
    tmp_co("i,m") = R1("e,n") * W_mnie("m,n,i,e").block(b.occr);
    tildeR2Val("a,b,i,j") += tmp_co("i,m") * this->T2_("a,b,j,m").block(b.rrho);
    TArray tmp_Vo = TAmanager.malloc<dcomplex>("ho");
    tmp_Vo("i,m") = R1("e,n") * W_mnie("m,n,i,e").block(b.ochr);
    tildeR2Val("a,b,i,j") -= tmp_Vo("j,m") * this->T2_("a,b,i,m").block(b.rrco);

    // W_abei * R2
    tildeR2Val("a,b,i,j") += R1("e,i") * W_abei("a,b,e,j").block(b.rrrh);

    //H_abk;ijc ---1
    TArray tmp_uv = TAmanager.malloc<dcomplex>("rr");
    tmp_uv("a,e") = R1("f,m") * W_amef("a,m,e,f").block(b.rcrr);
    tmp("a,b,i,j") = tmp_uv("a,e") * this->T2_("b,e,i,j").block(b.rrch);
    tildeR2Val("a,b,i,j") -= tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") += tmp("b,a,i,j");

    // F_ae * R2
    tildeR2Val("a,b,i,j") -= R2Val("a,b,m,j") * F_mi("m,i").block(b.cc);
    tildeR2Val("a,b,i,j") -= R2("a,b,i,m") * F_mi("m,j").block(b.ch);
    tildeR2Val("a,b,i,j") -= R2Val("a,b,i,m") * F_mi("m,j").block(b.hh);

    // H_abk,icd * R2
    tmp_co("i,m") = R2("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rroc));
    tmp_co("i,m") += R2Val("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rroh));
    tildeR2Val("a,b,i,j") += 0.5 * tmp_co("i,m") * this->T2_("a,b,j,m").block(b.rrho);
    tmp_Vo("i,m") = -1.0 * R2Val("e,f,n,i") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rroc));
    tildeR2Val("a,b,i,j") -= 0.5 * tmp_Vo("j,m") * this->T2_("a,b,i,m").block(b.rrco);

    // F_mi * R2
    tmp("a,b,i,j") = R2Val("a,e,i,j") * F_ae("b,e").block(b.rr);
    tildeR2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") -= tmp("b,a,i,j");

    // H_akl,ijd * R2
    tmp_uv("a,e") = R2("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rrcc));
    tmp_uv("a,e") += 2.0 * R2Val("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n").block(b.rrch)); // factor of two applied for remaining uirtual contribution
    tmp("a,b,i,j") = 0.5 * tmp_uv("a,e") * this->T2_("b,e,i,j").block(b.rrch);
    tildeR2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") -= tmp("b,a,i,j");

    // W_mbej * R2
    tmp("a,b,i,j") = R2Val("a,e,m,j") * W_mbej("m,b,e,i").block(b.crrc);
    tildeR2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") -= tmp("b,a,i,j");
    tmp("a,b,i,j") = R2("a,e,m,i") * W_mbej("m,b,e,j").block(b.crrh);
    tildeR2Val("a,b,i,j") -= tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") += tmp("b,a,i,j");
    tmp("a,b,i,j") = R2Val("a,e,i,m") * W_mbej("m,b,e,j").block(b.hrrh);
    tildeR2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeR2Val("a,b,i,j") -= tmp("b,a,i,j");

    // W_mnij * R2
    tildeR2Val("a,b,i,j") += 0.5 * R2("a,b,m,n") * W_mnij("m,n,i,j").block(b.ccch);
    tildeR2Val("a,b,i,j") += R2Val("a,b,m,n") * W_mnij("m,n,i,j").block(b.chch); // factor of 2 applied for remaining uirtual contribution

    // W_abef * R2
    tildeR2Val("a,b,i,j") += 0.5 * R2Val("e,f,i,j") * W_abef("a,b,e,f").block(b.rrrr);

    TAmanager.free("rrch", std::move(tmp));
    TAmanager.free("rr", std::move(tmp_uv));
    TAmanager.free("ho", std::move(tmp_Vo));
    TAmanager.free("co", std::move(tmp_co));
  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::updateG_ae(const TArray &L2, const TArray &L2Val, TArray &G_ae) const {
    G_ae("a,e")  = - 0.5 * this->T2_("e,f,m,n").block(b.rrcc) * L2("a,f,m,n");
    G_ae("a,e") -=         this->T2_("e,f,m,n").block(b.rrch) * L2Val("a,f,m,n"); // factor of two applied for missing valence contribution
  }  

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::updateG_mi(const TArray &L2, const TArray &L2Val, TArray &G_mi) const{
    G_mi("m,i").block(b.oc)  = -0.5 * this->T2_("e,f,n,m").block(b.rrco) * L2("e,f,i,n"); // all equations have at least one factor of -1
    G_mi("m,i").block(b.oc) -=  0.5 * this->T2_("e,f,n,m").block(b.rrho) * L2Val("e,f,i,n");
    G_mi("m,i").block(b.oh)  =  0.5 * this->T2_("e,f,n,m").block(b.rrco) * L2Val("e,f,n,i"); // factor of -1 due to swapping of i and n in L2Val
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formL1_tilde(const TArray &L1, const TArray &L2, const TArray &L2Val,
                                             const TArray &G_ae, const TArray &G_mi, TArray &tildeL1) const {
	  
    tildeL1("a,i")  = F_ae("e,a").block(b.rr) * L1("e,i");
    tildeL1("a,i") -= F_mi("i,m").block(b.cc) * L1("a,m");
    tildeL1("a,i") += L1("e,m") * W_mbej("i,e,a,m").block(b.crrc);

    tildeL1("a,i") += 0.5 * L2("e,f,i,m") * W_abei("e,f,a,m").block(b.rrrc);
    tildeL1("a,i") += 0.5 * L2Val("e,f,i,m") * W_abei("e,f,a,m").block(b.rrrh);

    tildeL1("a,i") -= 0.5 * L2("a,e,m,n") * W_mbij("i,e,m,n").block(b.crcc);
    tildeL1("a,i") -= L2Val("a,e,m,n") * W_mbij("i,e,m,n").block(b.crch); // factor of 2 applied for extra ualence term

    tildeL1("a,i") -= G_ae("e,f") * W_amef("e,i,f,a").block(b.rcrr);

    tildeL1("a,i") -= G_mi("m,n").block(b.oc) * W_mnie("m,i,n,a").block(b.occr);
    tildeL1("a,i") -= G_mi("m,n").block(b.oh) * W_mnie("m,i,n,a").block(b.ochr);

  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formL2_tilde(const TArray &L1, const TArray &L2,  const TArray &L2Val,
                                             const TArray &G_ae, const TArray &G_mi, TArray &tildeL2) const {

    TAManager &TAmanager = TAManager::get();

    TArray tmp = TAmanager.malloc<dcomplex>("rrcc");
    tmp("a,b,i,j")  = L2("a,e,i,j") * F_ae("e,b").block(b.rr);
    tildeL2("a,b,i,j")  = tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("b,a,i,j");

    tmp("a,b,i,j") = L2("a,b,i,m") * F_mi("j,m").block(b.cc);
    tildeL2("a,b,i,j") -= tmp("a,b,i,j");
    tildeL2("a,b,i,j") += tmp("a,b,j,i");

    tmp("a,b,i,j") = L2Val("a,b,i,m") * F_mi("j,m").block(b.ch);
    tildeL2("a,b,i,j") -= tmp("a,b,i,j");
    tildeL2("a,b,i,j") += tmp("a,b,j,i");

    tildeL2("a,b,i,j") += 0.5 * L2("a,b,m,n") * W_mnij("i,j,m,n").block(b.cccc);
    tildeL2("a,b,i,j") += L2Val("a,b,m,n") * W_mnij("i,j,m,n").block(b.ccch); // factor of two applied for missing valence contribution

    tildeL2("a,b,i,j") += 0.5 * L2("e,f,i,j") * W_abef("e,f,a,b").block(b.rrrr);

    tmp("a,b,i,j") = L1("e,i") * W_amef("e,j,a,b").block(b.rcrr);
    tildeL2("a,b,i,j") += tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("a,b,j,i");

    tmp("a,b,i,j") = L1("a,m") * W_mnie("i,j,m,b").block(b.cccr);
    tildeL2("a,b,i,j") -= tmp("a,b,i,j");
    tildeL2("a,b,i,j") += tmp("b,a,i,j");

    tmp("a,b,i,j") = L2("a,e,i,m") * W_mbej("j,e,b,m").block(b.crrc);
    tildeL2("a,b,i,j") += tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("a,b,j,i");
    tildeL2("a,b,i,j") -= tmp("b,a,i,j");
    tildeL2("a,b,i,j") += tmp("b,a,j,i");

    tmp("a,b,i,j") = L2Val("a,e,i,m") * W_mbej("j,e,b,m").block(b.crrh);
    tildeL2("a,b,i,j") += tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("a,b,j,i");
    tildeL2("a,b,i,j") -= tmp("b,a,i,j");
    tildeL2("a,b,i,j") += tmp("b,a,j,i");

    tmp("a,b,i,j") = L1("a,i") * F_me("j,b").block(b.cr);
    tildeL2("a,b,i,j") += tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("a,b,j,i");
    tildeL2("a,b,i,j") -= tmp("b,a,i,j");
    tildeL2("a,b,i,j") += tmp("b,a,j,i");

    tmp("a,b,i,j") = conj(this->antiSymMoints["vvoo"]("a,e,i,j").block(b.rrcc)) * G_ae("b,e");
    tildeL2("a,b,i,j") += tmp("a,b,i,j");
    tildeL2("a,b,i,j") -= tmp("b,a,i,j");

    tmp("a,b,i,j") = conj(this->antiSymMoints["vvoo"]("a,b,m,i").block(b.rroc)) * G_mi("m,j").block(b.oc);
    tildeL2("a,b,i,j") += tmp("a,b,i,j"); // signs flipped on these contributions to account for
    tildeL2("a,b,i,j") -= tmp("a,b,j,i"); // the flip in i and m in the MO term

    TAmanager.free("rrcc", std::move(tmp));
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formL2Val_tilde(const TArray &L1, const TArray &L2, const TArray &L2Val, 
                                                const TArray &G_ae, const TArray &G_mi, TArray &tildeL2Val) const {
    TAManager &TAmanager = TAManager::get();

    TArray tmp = TAmanager.malloc<dcomplex>("rrch");
    tmp("a,b,i,j") = L2Val("a,e,i,j") * F_ae("e,b").block(b.rr);
    tildeL2Val("a,b,i,j")  =  tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") -=  tmp("b,a,i,j");

    tildeL2Val("a,b,i,j") -= L2("a,b,i,m") * F_mi("j,m").block(b.hc);
    tildeL2Val("a,b,i,j") -= L2Val("a,b,i,m") * F_mi("j,m").block(b.hh);

    tildeL2Val("a,b,i,j") -= L2Val("a,b,m,j") * F_mi("i,m").block(b.cc); // factor of -1 for swapping m and j

    tildeL2Val("a,b,i,j") += 0.5 * L2("a,b,m,n") * W_mnij("i,j,m,n").block(b.chcc);
    tildeL2Val("a,b,i,j") += L2Val("a,b,m,n") * W_mnij("i,j,m,n").block(b.chch); // factor of two applied for missing valence contribution

    tildeL2Val("a,b,i,j") += 0.5 * L2Val("e,f,i,j") * W_abef("e,f,a,b").block(b.rrrr);

    tildeL2Val("a,b,i,j") += L1("e,i") * W_amef("e,j,a,b").block(b.rhrr);

    tmp("a,b,i,j") = L1("a,m") * W_mnie("i,j,m,b").block(b.chcr);
    tildeL2Val("a,b,i,j") -= tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") += tmp("b,a,i,j");

    tmp("a,b,i,j") = L2("a,e,i,m") * W_mbej("j,e,b,m").block(b.hrrc);
    tildeL2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") -= tmp("b,a,i,j");

    tmp("a,b,i,j") = L2Val("a,e,i,m") * W_mbej("j,e,b,m").block(b.hrrh);
    tildeL2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") -= tmp("b,a,i,j");

    tmp("a,b,i,j") = -L2Val("a,e,m,j") * W_mbej("i,e,b,m").block(b.crrc);// factor of -1 for swapping m and j
    tildeL2Val("a,b,i,j") -= tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") += tmp("b,a,i,j");

    tmp("a,b,i,j") = L1("a,i") * F_me("j,b").block(b.hr);
    tildeL2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") -= tmp("b,a,i,j");

    tmp("a,b,i,j") = conj(this->antiSymMoints["vvoo"]("a,e,i,j").block(b.rrch)) * G_ae("b,e");
    tildeL2Val("a,b,i,j") += tmp("a,b,i,j");
    tildeL2Val("a,b,i,j") -= tmp("b,a,i,j");

    tildeL2Val("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,b,m,i").block(b.rroc)) * G_mi("m,j").block(b.oh); // sign flipped to account for swap of i and m in MO 
    tildeL2Val("a,b,i,j") -= conj(this->antiSymMoints["vvoo"]("a,b,m,j").block(b.rroh)) * G_mi("m,i").block(b.oc);

    TAmanager.free("rrch", std::move(tmp));
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::runLambda() {
    TAManager &TAmanager = TAManager::get();

    initializeLambda();
    MBExpansion<MatsT>& Lg = *(this->Lg_);
    TArray &L1_ = Lg.get_tensor("OneBody");
    TArray &L2_ = Lg.get_tensor("TwoBody_core");
    TArray &L2Val_ = Lg.get_tensor("TwoBody_val");
    MBExpansion<MatsT> L_old(Lg);

    std::shared_ptr<DIISTA<MatsT>> ldiis = nullptr;
    if(this->ccSettings_.useDIIS){
      std::cout << " CVS Lambda DIIS NYI, performing optimization without " << std::endl;
      //ldiis = std::make_shared<DIISTA<MatsT>>(ccSettings_.nDIIS);
    }

    MatsT pseudoEnergy = 0.0;

    std::cout << std::endl;
    std::cout << "CVS Lambda Iterations" << std::endl;
    std::cout << std::endl;
    std::cout << "    ";
    std::cout << " iter";
    std::cout << "                  PE";
    std::cout << "                 dPE";
    std::cout << "                  dL";
    std::cout << std::endl;

    for (auto iter = 0; iter < this->ccSettings_.maxiter; iter++){
      L_old = Lg;

      updateG_ae(L_old.get_tensor("TwoBody_core"), L_old.get_tensor("TwoBody_val"),  G_ae);
      updateG_mi(L_old.get_tensor("TwoBody_core"), L_old.get_tensor("TwoBody_val"),  G_mi);
      std::cout << std::setw(20) << std::setprecision(12);

      formL1_tilde(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody_core"), L_old.get_tensor("TwoBody_val"), G_ae, G_mi, L1_);
      L1_("a,i") += F_me("i,a").block(b.cr);
      L1_("a,i") = L_old.get_tensor("OneBody")("a,i") + L1_("a,i") * D_ai("a,i").block(b.rc);

      formL2_tilde(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody_core"), L_old.get_tensor("TwoBody_val"), G_ae, G_mi, L2_);
      L2_("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrcc));
      L2_("a,b,i,j") = L_old.get_tensor("TwoBody_core")("a,b,i,j") + L2_("a,b,i,j") * D_abij("a,b,i,j").block(b.rrcc);

      formL2Val_tilde(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody_core"), L_old.get_tensor("TwoBody_val"), G_ae, G_mi, L2Val_);
      L2Val_("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrch));
      L2Val_("a,b,i,j") = L_old.get_tensor("TwoBody_val")("a,b,i,j") + L2Val_("a,b,i,j") * D_abij("a,b,i,j").block(b.rrch);

      L_old.axpy(-1.0, Lg);

      MatsT PE_old = pseudoEnergy;
      MatsT PEOneBody = this->fockMatrix_ta["vo"]("a,i").block(b.rc).dot(L1_("a,i"));
      MatsT PETwoBodyL2 = this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrcc).dot(L2_("a,b,i,j"));
      MatsT PETwoBodyL2Val = this->antiSymMoints["vvoo"]("a,b,i,j").block(b.rrch).dot(L2Val_("a,b,i,j"));
      MatsT PETwoBodyL1 = this->antiSymMoints["vvoo"]("c,d,k,l").block(b.rrcc).dot(L1_("c,k") * L1_("d,l"));

      pseudoEnergy = PEOneBody + 0.25 * PETwoBodyL2 + 0.5 * (PETwoBodyL2Val + PETwoBodyL1);

      double dL = L_old.norm();

      std::cout << "    ";
      std::cout << std::setw(5) << iter;
      std::cout << std::setw(20) << std::setprecision(12) << std::fixed << pseudoEnergy;
      std::cout << std::setw(20) << std::setprecision(12) << std::fixed << pseudoEnergy - PE_old;
      std::cout << std::setw(20) << std::setprecision(12) << std::fixed << dL;
      std::cout << std::endl;
      if ( std::abs(pseudoEnergy - PE_old) < this->ccSettings_.eConv and dL < this->ccSettings_.tConv) {
        std::cout << "CVS Lambda iteration converged in " << iter << " steps." << std::endl;
	break;
      }

      if (iter == this->ccSettings_.maxiter-1){
        CErr(std::string("CVS Lambda iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }

    }

    finalizeLambda();
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::buildCVSLeftIntermediates(){}

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::cleanNonCVSIntermediates(){
    TAManager &TAmanager = TAManager::get();

    // I don't know why, but freeing muMatrix causes segfault
    //for (auto it = muMatrix.begin(); it != muMatrix.end(); it++)
    //  if (it->first.substr(1) != "rr"){
    //    TAmanager.free(it->first.substr(1), std::move(it->second), true);
    //    muMatrix.erase(it);
    //  }

    TAmanager.free("rroo", std::move(tau), true);
   
    // I don't know why, but freeing T1 T2 causes segfault
    //TAmanager.free("or", std::move(this->T1_), true);
    //TAmanager.free("oorr", std::move(this->T2_), true);

  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::finalizeLambda() {
    // L needs G_ae and G_mi
    //TAManager &TAmanager = TAManager::get();
    //if (G_ae) TAmanager.free("rr", std::move(G_ae), true);
    //if (G_mi) TAmanager.free("oo", std::move(G_mi), true);
  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::initializeLambda() {
    TAManager &TAmanager = TAManager::get();
    
    TArray &LG1_rc = this->Lg_->get_tensor("OneBody");
    TArray &LG2_rrcc = this->Lg_->get_tensor("TwoBody_core");
    TArray &LG2_rrch = this->Lg_->get_tensor("TwoBody_val");
    this->Lg_->zeroBody() = 1.0;
    LG1_rc("a,i") = conj(this->T1_("a,i").block(b.rc));
    LG2_rrcc("a,b,i,j") = conj(this->T2_("a,b,i,j").block(b.rrcc));
    LG2_rrch("a,b,i,j") = conj(this->T2_("a,b,i,j").block(b.rrch));
    TA::get_default_world().gop.fence();

    if (not G_ae.is_initialized()){
      G_ae = TAmanager.malloc<MatsT>("rr");
    }

    if (not G_mi.is_initialized()){
      G_mi = TAmanager.malloc<MatsT>("oo");
    }
  }

  // for CCS_guess 
  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::fillGuess(MatsT *guess_vec, size_t n_vec) const{
 
    if (n_vec > nOVshift_) {
      CErr("CVSEOMCCSD: asking for more roots than the single excitation space dimension.");
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


  template <typename MatsT>
  typename Davidson<MatsT>::VecsGen_t CVSEOMCCSD<MatsT>::EmptyDavidsonVectorBuilder(){
    // Algorithm with implicit Hbar matrix
    typename Davidson<MatsT>::VecsGen_t vecsGenEOM;
    if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
      vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
        return std::make_shared<MBExpansionSet<MatsT>>(this->tensor_builder_, nVec, this->savFile_);
      }; // implicit vecsGenerator

      return vecsGenEOM;
    } else {
      CErr("CVSEOMCCSD: only IMPLICIT algorithm is implemented.");
    }
    return vecsGenEOM;

  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t CVSEOMCCSD<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
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
          const MBExpansion<MatsT> &Vi = V_ptr->get(i + Vshift);
          MBExpansion<MatsT> &AVi = AV_ptr->get(i + AVshift);
          buildSigma(Vi, AVi, eigenVecType);
          TA::get_default_world().gop.fence();
          AVi.enforceSymmetry();
        }

      }; // implicit sigmaBuilder
    } else {
      CErr("CVSEOMCCSD: only IMPLICIT algorithm is implemented.");
    }
    return this->funcEOM;
  }


  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t CVSEOMCCSD<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){

      double PCsmall = this->eomSettings.davidson_preCond_small;
      size_t n_c = nCVSOCore_;

      this->PCEOM = [this, eomDiag, curEig, PCsmall, n_c]( size_t nVec, SolverVectors<MatsT> &V,
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
          if constexpr (std::is_same_v<MatsT, dcomplex>) {
            curEigI = curEig[iVec];
          } else {
            curEigI = std::real(curEig[iVec]);
          }

          TA::foreach_inplace(curB.get_tensor("OneBody"), [iVec, curEigI, eomDiag, this, PCsmall,n_c](TA::Tensor<MatsT> &tile){
          const auto& lobound = tile.range().lobound();
          const auto& upbound = tile.range().upbound();

          MatsT denom = 0.0;
          std::vector<std::size_t> x{0, 0};
          for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
            for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
              if (x[0] == x[1])
                continue;
              denom = curEigI - eomDiag[this->Hbar_dimension_offsets.at("ai")+toCompoundS(x[0], x[1])];
              if (std::abs(denom) >= PCsmall) tile[x] /= denom;
            }
          });
          TA::get_default_world().gop.fence();
          TA::foreach_inplace(curB.get_tensor("TwoBody_core"), [iVec, curEigI, eomDiag, this, PCsmall, n_c](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            MatsT denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                    if (x[0] == x[1])
                      continue;
                    if (x[2] == x[3])
                      continue;
                    size_t a = x[0], b = x[1], i = x[2], j = x[3];
                    //signD(a,b,i,j);
                    denom = curEigI - eomDiag[this->Hbar_dimension_offsets.at("abij")+toCompoundD(a, b, i, j)];
                    if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                  }
              }
          });
          TA::get_default_world().gop.fence();
          TA::foreach_inplace(curB.get_tensor("TwoBody_val"), [iVec, curEigI, eomDiag, this, PCsmall, n_c](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            MatsT denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                    if (x[0] == x[1])
                      continue;
                    size_t a = x[0], b = x[1], i = x[2], j = x[3];
                    //signD(a,b,i,j);
                    denom = curEigI - eomDiag[this->Hbar_dimension_offsets.at("abiJ")+toCompoundD(a, b, i, j+n_c)];
                    if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                  }
              }
          });

          curB.enforceSymmetry();
        }
      }; // implicit preConditioner
      return this->PCEOM;
  }

  template <typename MatsT>
  CVSEOMCCSD<MatsT>::~CVSEOMCCSD(){
  } 

}
