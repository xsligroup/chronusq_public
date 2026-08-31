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

#include <chronusq_sys.hpp>
#include <util/threads.hpp>
#include <coupledcluster.hpp>

namespace ChronusQ{


  template <typename MatsT>
  CCSDT<MatsT>::CCSDT(const SafeFile &savFile,
                          CCIntermediates<MatsT> &intermediates,
                          const CoupledClusterSettings &ccSettings):
      CCBase<MatsT>(savFile, intermediates, ccSettings),
      T3_(intermediates.T->get_tensor("ThreeBody")),
      tau_(intermediates.tau),
      tilde_tau_(intermediates.tilde_tau),
      Fae_(intermediates.F_ae),
      Fmi_(intermediates.F_mi),
      Fme_(intermediates.F_me),
      Wmnie_(intermediates.W_mnie), // not in CCSD.hpp
      Wamef_(intermediates.W_amef), // not in CCSD.hpp
      Wmnij_(intermediates.W_mnij),
      Wabef_(intermediates.W_abef),
      Wmbej_(intermediates.W_mbej),
      Wmbij_(intermediates.W_mbij), // not in CCSD.hpp
      Wabei_(intermediates.W_abei), // not in CCSD.hpp
      eps(intermediates.eps){}

  template <typename MatsT>
  void CCSDT<MatsT>::run() {
    runConventional();
  }

  template <typename MatsT>
  void CCSDT<MatsT>::initIntermediates() {

    TAManager &TAmanager = TAManager::get();
    if (not tau_.is_initialized()){
      tau_ = TAmanager.malloc<MatsT>("vvoo");
    }
    if (not tilde_tau_.is_initialized()){
      tilde_tau_ = TAmanager.malloc<MatsT>("vvoo");
    }
    if (not Fae_.is_initialized()){
      Fae_ = TAmanager.malloc<MatsT>("vv");
    }
    if (not Fmi_.is_initialized()){
      Fmi_ = TAmanager.malloc<MatsT>("oo");
    }
    if (not Fme_.is_initialized()){
      Fme_ = TAmanager.malloc<MatsT>("ov");
    }
    if (not Wmnie_.is_initialized()){
      Wmnie_ = TAmanager.malloc<MatsT>("ooov");
    }
    if (not Wamef_.is_initialized()){
      Wamef_ = TAmanager.malloc<MatsT>("vovv");
    }
    if (not Wmnij_.is_initialized()){
      Wmnij_ = TAmanager.malloc<MatsT>("oooo");
    }
    if (not Wabef_.is_initialized()){
      Wabef_ = TAmanager.malloc<MatsT>("vvvv");
    }
    if (not Wmbej_.is_initialized()){
      Wmbej_ = TAmanager.malloc<MatsT>("ovvo");
    }
    if (not Wmbij_.is_initialized()){
      Wmbij_ = TAmanager.malloc<MatsT>("ovoo");
    }
    if (not Wabei_.is_initialized()){
      Wabei_ = TAmanager.malloc<MatsT>("vvvo");
    }

  }
  
  template <typename MatsT>
  void CCSDT<MatsT>::initAmplitudes() {
    if(this->ccSettings_.restart){
      size_t size = this->T_.length();
      MatsT * t_amp = CQMemManager::get().malloc<MatsT>(size);
      TA::get_default_world().gop.fence();
      if (MPIRank() == 0) this->savFile_.readData("/CC/T_AMPLITUDE", t_amp);
      if (MPIRank() == 0) this->savFile_.readData("/CC/CORRELATION_ENERGY", &this->CorrE);
      MPIBCast(t_amp, size, 0, MPI_COMM_WORLD);
      MPIBCast(&this->intermediates_.E_ref, 1, 0, MPI_COMM_WORLD);
      MPIBCast(&this->CorrE               , 1, 0, MPI_COMM_WORLD);
      MPIBCast(t_amp, size, 0, MPI_COMM_WORLD);
      TA::get_default_world().gop.fence();
      this->T_.fromRaw(t_amp, false);
      TA::get_default_world().gop.fence();
      CQMemManager::get().free(t_amp);
    } else {
      this->T_.scale(0.0);
    }
  }

  template <typename MatsT>
  void CCSDT<MatsT>::printAnalysis() {
    std::cout << bannerTop << std::endl;
    std::cout << "Coupled Cluster Wave Function Analysis:" << std::endl;
    std::cout << std::setw(22) << "   max(|t1|)  " << std::setprecision(4) << abs_max(this->T1_).get() << std::endl;
    std::cout << std::setw(22) << "   max(|t2|)  " << std::setprecision(4) << abs_max(this->T2_).get() << std::endl;
    std::cout << std::setw(22) << "   max(|t3|)  " << std::setprecision(4) << abs_max(this->T3_).get() << std::endl;
    double normT1Sq = squared_norm(this->T1_).get();
    double normT2Sq = squared_norm(this->T2_).get();
    double normT3Sq = squared_norm(this->T3_).get();
    std::cout << std::setw(22) << "   sum of t1 weights  " << std::setprecision(4) << normT1Sq << std::endl;
    std::cout << std::setw(22) << "   sum of t2 weights  " << std::setprecision(4) << normT2Sq << std::endl;
    std::cout << std::setw(22) << "   sum of t3 weights  " << std::setprecision(4) << normT3Sq << std::endl;
    std::cout << std::setw(22) << "   T1 diagnostic  " << std::setprecision(4) << sqrt(normT1Sq / this->intermediates_.nOcc) << std::endl;
    std::cout << std::setw(22) << "   T2 diagnostic  " << std::setprecision(4) << sqrt(normT2Sq / this->intermediates_.nOcc) << std::endl;
    std::cout << std::setw(22) << "   T3 diagnostic  " << std::setprecision(4) << sqrt(normT3Sq / this->intermediates_.nOcc) << std::endl;
    std::cout << bannerEnd << std::endl;
  }

  /**
   * Build Eq. III(d.1) and Eq. III(d.2)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tau_and_tilde_tau() {
    tau_("a,b,i,j") = 0.5 * this->T1_("a,i") * this->T1_("b,j");
    tau_("a,b,i,j") -= tau_("b,a,i,j");
    tau_("a,b,i,j") -= tau_("a,b,j,i");

    tilde_tau_("a,b,i,j") = 0.5 * tau_("a,b,i,j") + this->T2_("a,b,i,j");
    tau_("a,b,i,j") += this->T2_("a,b,i,j");
  }

  /**
   * Build Eq. III(a.1)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Fae() {
    Fae_("a,e") = this->fockMatrix_ta["vv"]("a,e");
    Fae_("a,e") -= 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("a,m");
    Fae_("a,e") += this->T1_("f,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
    Fae_("a,e") -= 0.5 * tilde_tau_("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.2)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Fmi() {
    Fmi_("m,i") = this->fockMatrix_ta["oo"]("m,i");
    Fmi_("m,i") += 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("e,i");
    Fmi_("m,i") -= this->T1_("e,n") * conj(this->antiSymMoints["vooo"]("e,i,m,n"));
    Fmi_("m,i") += 0.5 * tilde_tau_("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n")); // TODO: better performance if we have this->antiSymMoints["oovv"]
  }

  /**
   * Build Eq. III(a.3)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Fme() {
    Fme_("m,e") = this->fockMatrix_ta["ov"]("m,e");
    Fme_("m,e") += this->T1_("f,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.4)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Wmnij() {
    Wmnij_("m,n,i,j") = -this->T1_("e,j") * conj(this->antiSymMoints["vooo"]("e,i,m,n")); // TODO: better performance if we have this->antiSymMoints["ooov"]
    Wmnij_("m,n,i,j") -= Wmnij_("m,n,j,i");
    Wmnij_("m,n,i,j") += this->antiSymMoints["oooo"]("m,n,i,j");
    Wmnij_("m,n,i,j") += 0.25 * tau_("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.5)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Wabef() {
    Wabef_("a,b,e,f") = -this->T1_("b,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
    Wabef_("a,b,e,f") -= Wabef_("b,a,e,f");
    Wabef_("a,b,e,f") += this->antiSymMoints["vvvv"]("a,b,e,f");
    Wabef_("a,b,e,f") += 0.25 * tau_("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.6)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::build_tilde_Wmbej() {
    Wmbej_("m,b,e,j") = -this->antiSymMoints["vovo"]("b,m,e,j"); // TODO: better performance if we have this->antiSymMoints["voov"]

    Wmbej_("m,b,e,j") -= this->T1_("f,j") * conj(this->antiSymMoints["vvvo"]("e,f,b,m")); // TODO: better performance if we have this->antiSymMoints["ovvv"]
    Wmbej_("m,b,e,j") -= this->T1_("b,n") * conj(this->antiSymMoints["vooo"]("e,j,m,n"));

    TArray tmp = TAManager::get().malloc<MatsT>("vvoo");
    tmp("f,b,j,n") = 0.5 * this->T2_("f,b,j,n");
    tmp("f,b,j,n") += this->T1_("f,j") * this->T1_("b,n");
    Wmbej_("m,b,e,j") -= tmp("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    TAManager::get().free("vvoo", std::move(tmp));
  }

  template <typename MatsT>
  void CCSDT<MatsT>::formW_mnie() {
    Wmnie_("m,n,i,e") = - conj(this->antiSymMoints["vooo"]("e,i,m,n")) + this->T1_("f,i") * conj(this->antiSymMoints["vvoo"]("f,e,m,n"));
  }

  template <typename MatsT>
  void CCSDT<MatsT>::formW_amef() {
    Wamef_("a,m,e,f") = conj(this->antiSymMoints["vvvo"]("e,f,a,m")) - this->T1_("a,n") * conj(this->antiSymMoints["vvoo"]("e,f,n,m"));
  }

  template <typename MatsT>
  void CCSDT<MatsT>::formW_mbij() {
    Wmbij_("m,b,i,j") = - this->antiSymMoints["vooo"]("b,m,i,j") - Fme_("m,e") * this->T2_("b,e,i,j");
    // modify Wmnij_ on the fly
    TArray tmp_mnij = TAManager::get().malloc<MatsT>("oooo");
    tmp_mnij("m,n,i,j") = Wmnij_("m,n,i,j") + 0.25 * tau_("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    Wmbij_("m,b,i,j") += - this->T1_("b,n") * tmp_mnij("m,n,i,j");
    TAManager::get().free("oooo", std::move(tmp_mnij));
    // end modify Wmnij_
    Wmbij_("m,b,i,j") += - 0.5 * conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * tau_("e,f,i,j");
    Wmbij_("m,b,i,j") += - conj(this->antiSymMoints["vooo"]("e,i,m,n")) * this->T2_("b,e,j,n");
    Wmbij_("m,b,i,j") += conj(this->antiSymMoints["vooo"]("e,j,m,n")) * this->T2_("b,e,i,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,j") = - this->antiSymMoints["vovo"]("b,m,e,j") - this->T2_("b,f,n,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    Wmbij_("m,b,i,j") += this->T1_("e,i") * tmp("m,b,e,j");
    Wmbij_("m,b,i,j") += - this->T1_("e,j") * tmp("m,b,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
    if (this->include_t3_) {
      // T3 part
      Wmbij_("m,b,i,j") += 0.5 * this->T3_("e,b,f,i,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    }
  }

  template <typename MatsT>
  void CCSDT<MatsT>::formW_abei() {
    Wabei_("a,b,e,i") = this->antiSymMoints["vvvo"]("a,b,e,i") - Fme_("m,e") * this->T2_("a,b,m,i");
    // modify Wabef on the fly
    TArray tmp_abef = TAManager::get().malloc<MatsT>("vvvv");
    tmp_abef("a,b,e,f") = Wabef_("a,b,e,f") + 0.25 * tau_("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    Wabei_("a,b,e,i") += this->T1_("f,i") * tmp_abef("a,b,e,f");
    TAManager::get().free("vvvv", std::move(tmp_abef));
    // end modify Wabef
    Wabei_("a,b,e,i") +=  0.5 * conj(this->antiSymMoints["vooo"]("e,i,m,n")) * tau_("a,b,m,n");
    Wabei_("a,b,e,i") += conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * this->T2_("a,f,m,i");
    Wabei_("a,b,e,i") += -conj(this->antiSymMoints["vvvo"]("e,f,a,m")) * this->T2_("b,f,m,i");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,i") = - this->antiSymMoints["vovo"]("b,m,e,i") - this->T2_("b,f,n,i") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    Wabei_("a,b,e,i") += - this->T1_("a,m") * tmp("m,b,e,i");
    Wabei_("a,b,e,i") +=  this->T1_("b,m") * tmp("m,a,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
    if (this->include_t3_) {
      // T3 part
      Wabei_("a,b,e,i") -= 0.5 * this->T3_("a,b,f,m,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    }
  }

  /**
   * Build Eq. III(a.1-6)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::buildIntermediates() {
    build_tau_and_tilde_tau();
    build_tilde_Fae();
    build_tilde_Fmi();
    build_tilde_Fme();
    build_tilde_Wmnij();
    build_tilde_Wabef();
    build_tilde_Wmbej();
    formW_mnie();
    formW_amef();
    formW_mbij();
    formW_abei();
  }

  /**
   * Build Eq. I(a)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::updateT1(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old) {
    this->T1_("a,i") = this->fockMatrix_ta["vo"]("a,i");
    this->T1_("a,i") += Fae_("a,e") * T1_old("e,i");
    this->T1_("a,i") -= T1_old("a,m") * Fmi_("m,i");
    this->T1_("a,i") += Fme_("m,e") * T2_old("a,e,i,m");
    this->T1_("a,i") -= T1_old("e,m") * this->antiSymMoints["vovo"]("a,m,e,i");
    this->T1_("a,i") += 0.5 * T2_old("e,f,i,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m")); // TODO: better performance if we have this->antiSymMoints["vvov"]
    this->T1_("a,i") -= 0.5 * T2_old("a,e,m,n") * conj(this->antiSymMoints["vooo"]("e,i,n,m")); // TODO: better performance if we have this->antiSymMoints["ovoo"]
    if (this->include_t3_) {
      // T3 part
      this->T1_("a,i") += 0.25 * T3_old("a,e,f,i,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    }

    this->T1_("a,i") = T1_old("a,i") + this->T1_("a,i") * this->Dai_("a,i");
  }

  /**
   * Build Eq. I(b)
   */
  template <typename MatsT>
  void CCSDT<MatsT>::updateT2(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old) {

    TAManager &TAmanager = TAManager::get();

    this->T2_("a,b,i,j") = this->antiSymMoints["vvoo"]("a,b,i,j");

    TArray TMPbe = TAmanager.malloc<MatsT>("vv");
    TMPbe("b,e") = Fae_("b,e");
    TMPbe("b,e") -= 0.5 * T1_old("b,m") * Fme_("m,e");
    TArray Pabij = TAmanager.malloc<MatsT>("vvoo");
    Pabij("a,b,i,j") = T2_old("a,e,i,j") * TMPbe("b,e");
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    this->T2_("a,b,i,j") -= Pabij("b,a,i,j");
    TAmanager.free("vv", std::move(TMPbe));

    TArray TMPmj = TAmanager.malloc<MatsT>("oo");
    TMPmj("m,j") = Fmi_("m,j");
    TMPmj("m,j") += 0.5 * T1_old("e,j") * Fme_("m,e");
    Pabij("a,b,i,j") = T2_old("a,b,i,m") * TMPmj("m,j");
    this->T2_("a,b,i,j") -= Pabij("a,b,i,j");
    this->T2_("a,b,i,j") += Pabij("a,b,j,i");
    TAmanager.free("oo", std::move(TMPmj));

    this->T2_("a,b,i,j") += 0.5 * tau_("a,b,m,n") * Wmnij_("m,n,i,j");
    this->T2_("a,b,i,j") += 0.5 * tau_("e,f,i,j") * Wabef_("a,b,e,f");

    TArray TMPmbij = TAmanager.malloc<MatsT>("ovoo");
    TMPmbij("m,b,i,j") = T1_old("e,i") * this->antiSymMoints["vovo"]("b,m,e,j");
    Pabij("a,b,i,j") = T1_old("a,m") * TMPmbij("m,b,i,j");
    Pabij("a,b,i,j") += T2_old("a,e,i,m") * Wmbej_("m,b,e,j");
    Pabij("a,b,i,j") -= Pabij("a,b,j,i");
    Pabij("a,b,i,j") -= Pabij("b,a,i,j");
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    TAmanager.free("ovoo", std::move(TMPmbij));

    Pabij("a,b,i,j") = T1_old("e,i") * this->antiSymMoints["vvvo"]("a,b,e,j");
    if (this->include_t3_) {
      // T3 part
      Pabij("a,b,i,j") -= 0.5 * Wmnie_("m,n,i,f") * T3_old("a,b,f,m,j,n");
    }
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    this->T2_("a,b,i,j") -= Pabij("a,b,j,i");

    Pabij("a,b,i,j") = -T1_old("b,m") * this->antiSymMoints["vooo"]("a,m,i,j");
    if (this->include_t3_) {
      // T3 part
      Pabij("a,b,i,j") += 0.5 * Wamef_("a,n,e,f") * T3_old("e,b,f,i,j,n");
    }
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    this->T2_("a,b,i,j") -= Pabij("b,a,i,j");

    TAmanager.free("vvoo", std::move(Pabij));
    if (this->include_t3_) {
      // T3 part
      this->T2_("a,b,i,j") += Fme_("m,e") * T3_old("a,b,e,i,j,m");
    }

    this->T2_("a,b,i,j") = T2_old("a,b,i,j") + this->T2_("a,b,i,j") * this->Dabij_("a,b,i,j");
  }

  template <typename MatsT>
  void CCSDT<MatsT>::updateT3(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old) {

    TAManager &TAmanager = TAManager::get();

    // <ijkacb| (Hbar T_2)_C |0>
    this->T3_("a,b,c,i,j,k")  = 0.25 * Wmbij_("m,b,i,j") * T2_old("a,c,k,m");
    this->T3_("a,b,c,i,j,k") += 0.25 * T2_old("e,b,i,j") * Wabei_("a,c,e,k");
    // subtract double counting of <ijkabc| Hbar1 1/2 T2^2 |0> from W_mbij and W_abei
    this->T3_("a,b,c,i,j,k") -= 0.25 * T2_old("e,b,i,j") * Fme_("m,e") * T2_old("a,c,k,m");

    TArray tmp_mi = TAmanager.malloc<MatsT>("oo");
    tmp_mi("m,i") = Fmi_("m,i") + 0.5 * Fme_("m,e") * T1_old("e,i");
    this->T3_("a,b,c,i,j,k") -= 1.0/12.0 * T3_old("a,b,c,i,j,m") * tmp_mi("m,k");
    TAmanager.free("oo",std::move(tmp_mi));

    TArray tmp_ae = TAmanager.malloc<MatsT>("vv");
    tmp_ae("a,e") = Fae_("a,e") - 0.5 * T1_old("a,m") * Fme_("m,e");
    this->T3_("a,b,c,i,j,k") += 1.0/12.0 * tmp_ae("a,e") * T3_old("e,b,c,i,j,k");
    TAmanager.free("vv",std::move(tmp_ae));

    TArray tmp_abef = TAmanager.malloc<MatsT>("vvvv");
    tmp_abef("a,b,e,f") = Wabef_("a,b,e,f") + 0.25 * tau_("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    this->T3_("a,b,c,i,j,k") += 1.0/24.0 * tmp_abef("a,b,e,f") * T3_old("e,f,c,i,j,k");
    TAManager::get().free("vvvv", std::move(tmp_abef));

    TArray tmp_mnij = TAmanager.malloc<MatsT>("oooo");
    tmp_mnij("m,n,i,j") = Wmnij_("m,n,i,j") + 0.25 * tau_("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    this->T3_("a,b,c,i,j,k") += 1.0/24.0 * T3_old("a,b,c,i,m,n") * tmp_mnij("m,n,j,k");
    TAmanager.free("oooo",std::move(tmp_mnij));

    TArray tmp_mbej = TAmanager.malloc<MatsT>("ovvo");
    tmp_mbej("m,b,e,j") = Wmbej_("m,b,e,j") - 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    this->T3_("a,b,c,i,j,k") += 0.25 * tmp_mbej("m,a,e,i") * T3_old("e,b,c,m,j,k");
    TAmanager.free("ovvo",std::move(tmp_mbej));

    TArray Pabcijk = TAmanager.malloc<MatsT>("vvvooo");
    // Antisymmetrize abc
    Pabcijk("a,b,c,i,j,k")  = this->T3_("a,b,c,i,j,k");
    Pabcijk("a,b,c,i,j,k") -= this->T3_("b,a,c,i,j,k");
    Pabcijk("a,b,c,i,j,k") -= this->T3_("c,b,a,i,j,k");
    Pabcijk("a,b,c,i,j,k") -= this->T3_("a,c,b,i,j,k");
    Pabcijk("a,b,c,i,j,k") += this->T3_("c,a,b,i,j,k");
    Pabcijk("a,b,c,i,j,k") += this->T3_("b,c,a,i,j,k");
    // Antisymmetrize ijk
    this->T3_("a,b,c,i,j,k")  = Pabcijk("a,b,c,i,j,k");
    this->T3_("a,b,c,i,j,k") -= Pabcijk("a,b,c,j,i,k");
    this->T3_("a,b,c,i,j,k") -= Pabcijk("a,b,c,k,j,i");
    this->T3_("a,b,c,i,j,k") -= Pabcijk("a,b,c,i,k,j");
    this->T3_("a,b,c,i,j,k") += Pabcijk("a,b,c,k,i,j");
    this->T3_("a,b,c,i,j,k") += Pabcijk("a,b,c,j,k,i");
    TAmanager.free("vvvooo",std::move(Pabcijk));
    
    double denomshift_ = this->ccSettings_.denomshift;

    TA::foreach_inplace( this->T3_, [ this, denomshift_ ](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
            for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
              for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                  for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){
                    size_t a = x[0]+this->nO_, b = x[1]+this->nO_, c = x[2]+this->nO_;
                    size_t i = x[3], j = x[4], k = x[5];
                    tile[x] /= ((eps[i] + eps[j] + eps[k] - eps[a] - eps[b] - eps[c]) - denomshift_);
                  }
    });
    TA::get_default_world().gop.fence();

    this->T3_("a,b,c,i,j,k") += T3_old("a,b,c,i,j,k");

  }

  template <typename MatsT>
  void CCSDT<MatsT>::runConventional(){

    auto cc_start = tick();

    std::shared_ptr<DIISTA<MatsT> > diis = nullptr;
    if(this->ccSettings_.useDIIS){
      diis = std::make_shared<DIISTA<MatsT>>(this->ccSettings_.nDIIS);
    }

    this->CorrE = 0.0; // make public and not initialized in coupledcluster.hpp

    initIntermediates();

    initAmplitudes();

    if (this->ccSettings_.skipCC && this->ccSettings_.restart) return;

    std::vector<std::string> tmp;
    tmp.push_back(std::string({this->vLabel_}));
    tmp.push_back(std::string({this->oLabel_}));
    tmp.push_back(std::string("OneBody"));
    tmp.push_back(std::string({this->vLabel_,this->vLabel_}));
    tmp.push_back(std::string({this->oLabel_,this->oLabel_}));
    tmp.push_back(std::string("TwoBody"));
    tmp.push_back(std::string({this->vLabel_,this->vLabel_,this->vLabel_}));
    tmp.push_back(std::string({this->oLabel_,this->oLabel_,this->oLabel_}));
    tmp.push_back(std::string("ThreeBody"));
    MBExpansion<MatsT> T_old(tmp);

    std::cout << std::setw(18) << std::left <<  "  CC Iterations";
    std::cout << std::setw(22) << std::right << "Corr. Energy (Eh)";
    std::cout << std::setw(19) << std::right << "\u0394Ec (Eh)";
    std::cout << std::setw(19) << std::right << "|\u0394T|";
    std::cout << '\n';
    std::cout << std::setw(18) << std::left <<  "  -------------";
    std::cout << std::setw(22) << std::right << "-----------------";
    std::cout << std::setw(18) << std::right << "--------";
    std::cout << std::setw(18) << std::right << "----";
    std::cout << '\n' << std::endl;

    if (this->ccSettings_.CCSDinit == true){
      this->include_t3_ = false;
    }
    else
      this->include_t3_ = true;

    for (auto iter = 0; iter < this->ccSettings_.maxiter; iter++){

      T_old = this->T_;

      buildIntermediates();

      updateT1(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"), T_old.get_tensor("ThreeBody"));
      updateT2(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"), T_old.get_tensor("ThreeBody"));
      if (this->include_t3_) {
        updateT3(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"), T_old.get_tensor("ThreeBody"));
      }

      if (this->ccSettings_.useDIIS){
        // diis
        this->doDIIS(T_old, diis);
      } else {
        T_old.axpy(-1, this->T_);
      }

      MatsT Eold = this->CorrE;
      this->getCorrEnergy();
      this->intermediates_.E_cc = this->intermediates_.E_ref + std::real(this->CorrE);
      double dE = std::abs(this->CorrE - Eold);

      double dT = T_old.norm();

      std::cout << std::setprecision(12) << std::fixed;
      std::cout << "  Iteration "  << std::setw(6) << std::left << iter;
      std::cout << std::setw(22) << std::right << std::fixed << this->CorrE;
      std::cout << std::setw(18) << std::right << std::fixed << dE;
      std::cout << std::setw(18) << std::right << std::fixed << std::abs(dT);
      std::cout << std::endl;

      if (dE < this->ccSettings_.eConv and dT < this->ccSettings_.tConv) {
        if (this->include_t3_) {
          std::cout << "\n  CC Completed: Corr. E is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->CorrE << " Eh" << std::endl;
          std::cout << "\n  CC Completed: Total E is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE << " Eh" << std::endl;
          std::cout << "\n  CC Completed: Iteration total time "<< std::setw(10) << std::right
                    << std::setprecision(6) << tock(cc_start) << " s" << std::endl;

          if(this->ccSettings_.save){
            size_t size = this->T_.length();
            MatsT * t_amp = CQMemManager::get().malloc<MatsT>(size);
            TA::get_default_world().gop.fence();
            this->T_.toRaw(t_amp, false);
            TA::get_default_world().gop.fence();
            if(this->savFile_.exists()){
              this->savFile_.safeWriteData("/CC/T_AMPLITUDE", t_amp, {size});
            }
            CQMemManager::get().free(t_amp);
          }

          if (this->savFile_.exists()) {
            this->savFile_.safeWriteData("/CC/REFERENCE_ENERGY",&this->intermediates_.E_ref, {1});
            this->savFile_.safeWriteData("/CC/CORRELATION_ENERGY",&this->CorrE, {1});
          }
          
          std::cout << BannerEnd << std::endl;

          printAnalysis();

          std::cout << BannerEnd << std::endl;
          
          break;
        }
        else {
          std::cout << "  CCSD initial guess converged, engaging T3 updates" << std::endl;
          this->include_t3_ = true;
//          this->ccSettings_.useDIIS = false;
        }
      }

      if(iter == this->ccSettings_.maxiter - 1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }
    }

  }




  template <typename MatsT>
  CCSDT<MatsT>::~CCSDT() {

    TAManager &TAmanager = TAManager::get();

  }

  template <typename MatsT>
  void CCSDT<MatsT>::cleanMemory(){
    TAManager &TAmanager = TAManager::get();
//    if(Fae_) TAmanager.free("vv", std::move(Fae_), true);
//    if(Fmi_) TAmanager.free("oo", std::move(Fmi_), true);
    if(tilde_tau_) TAmanager.free("vvoo", std::move(tilde_tau_), true);
  }

  template <typename MatsT>
  size_t CCSDT<MatsT>::estimate_mem_peak() const {
    // will need further checking as the TA objects are dynamically allocated and freed at runtime
    TAManager &TAmanager = TAManager::get();

    size_t nDIIS = this->ccSettings_.useDIIS ? this->ccSettings_.nDIIS : 0;
    size_t count = 0;
    //                                               1         2       3       4         5           6     7      8         9
    count += 9 * TAmanager.elem_per_TA("oo");     // muX_oo,   muY_oo, muZ_oo, coreH_oo, fock_oo,    Fmi_, TMPmj, moDen_oo, tmp_mi
    count += 6 * TAmanager.elem_per_TA("ov");     // muX_ov,   muY_ov, muZ_ov, coreH_ov, fock_ov,    Fme_
    count += 7 * TAmanager.elem_per_TA("vo");     // muX_vo,   muY_vo, muZ_vo, coreH_vo, fock_vo,    Dai,  T1
    count += 8 * TAmanager.elem_per_TA("vv");     // muX_vv,   muY_vv, muZ_vv, coreH_vv, fock_vv,    Fae_, TMPbe, tmp_ae
    count += 3 * TAmanager.elem_per_TA("oooo");   // ERI_oooo, Wmnij_, tmp_mnij
    count += 1 * TAmanager.elem_per_TA("ooov");   // Wmnie_
    count += 2 * TAmanager.elem_per_TA("ovoo");   // TMPmbij,  Wmbij_
    count += 1 * TAmanager.elem_per_TA("vooo");   // ERI_vooo
    count += 3 * TAmanager.elem_per_TA("ovvo");   // Wmbej_,   tmp,    tmp_mbej
    count += 1 * TAmanager.elem_per_TA("vovo");   // ERI_vovo
    count += 7 * TAmanager.elem_per_TA("vvoo");   // ERI_vvoo, Dabij,  T2,     tau_,     tilde_tau_, tmp,  Pabij,
    count += 2 * TAmanager.elem_per_TA("vvvo");   // ERI_vvvo, Wabei_,
    count += 1 * TAmanager.elem_per_TA("vovv");   // Wamef_
    count += 3 * TAmanager.elem_per_TA("vvvv");   // ERI_vvvv, W_abef, tmp_abef
    count += 3 * TAmanager.elem_per_TA("vvvooo"); // T3,       Pabcijk, ????

    if (nDIIS) {
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vo");   // T1 DIIS copy?
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vvoo"); // T2 DIIS copy?
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vvvooo"); // T3 DIIS copy?
    }

    return count * sizeof(MatsT);
  }

};


