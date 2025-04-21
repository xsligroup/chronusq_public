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
  DFCCSD<MatsT>::DFCCSD(const SafeFile &savFile,
                          CCIntermediates<MatsT> &intermediates,
                          const CoupledClusterSettings &ccSettings):
      CCBase<MatsT>(savFile, intermediates, ccSettings),
      tau_(intermediates.tau),
      tilde_tau_(intermediates.tilde_tau),
      Fae_(intermediates.F_ae),
      Fmi_(intermediates.F_mi),
      Fme_(intermediates.F_me),
      Wmnij_(intermediates.W_mnij),
      Wabef_(intermediates.W_abef),
      Wmbej_(intermediates.W_mbej){}

  template <typename MatsT>
  void DFCCSD<MatsT>::run() {
    runConventional();
  }

  template <typename MatsT>
  void DFCCSD<MatsT>::initIntermediates() {
    TAManager &TAmanager = TAManager::get();

    // DF-specific intermediates -- T1_transformed F_N
    if (not F_ov.is_initialized()){
      F_ov = TAmanager.malloc<MatsT>("ov");
    }
    if (not F_oo.is_initialized()){
      F_oo = TAmanager.malloc<MatsT>("oo");
    }
    if (not F_vv.is_initialized()){
      F_vv = TAmanager.malloc<MatsT>("vv");
    }
    if (not F_vo.is_initialized()){
      F_vo = TAmanager.malloc<MatsT>("vo");
    }

    // DF-specific intermediates -- T1_transformed B(Q)_pq
    if (not B_oo.is_initialized()){
      B_oo = TAmanager.malloc<MatsT>("boo");
    }
    if (not B_vv.is_initialized()){
      B_vv = TAmanager.malloc<MatsT>("bvv");
    }
    if (not B_vo.is_initialized()){
      B_vo = TAmanager.malloc<MatsT>("bvo");
    }

  }
  
  template <typename MatsT>
  void DFCCSD<MatsT>::initAmplitudes() {
    if(this->ccSettings_.restart){
      size_t size = this->T_.length();
      dcomplex * t_amp = CQMemManager::get().malloc<dcomplex>(size);
      TA::get_default_world().gop.fence();
      if (MPIRank() == 0) this->savFile_.readData("/CC/T_AMPLITUDE", t_amp);
      if (MPIRank() == 0) this->savFile_.readData("/CC/REFERENCE_ENERGY",   &this->intermediates_.E_ref);
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
  void DFCCSD<MatsT>::t1TransformIntegrals(){
    // Dressed 3-index integrals
    // Dressed Ws can be constructed from dressed 3-index integrals on the fly (see update_T1 and update_T2)
//#define DEBUG_DFCCSD
#ifdef DEBUG_DFCCSD
    std::cout << "I am in B_oo" << std::endl;
#endif
    // Qo2 + Qo2v
    B_oo("Q,m,i")  = this->riMoInts["boo"]("Q,m,i");
    B_oo("Q,m,i") += this->riMoInts["bov"]("Q,m,e") * this->T1_("e,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in B_vv" << std::endl;
#endif
    // Qv2
    B_vv("Q,a,e")  = this->riMoInts["bvv"]("Q,a,e");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in B_vo" << std::endl;
#endif
    // Qvo + Qo2v + Qov2
    B_vo("Q,a,i")  = this->riMoInts["bvo"]("Q,a,i");
    B_vo("Q,a,i") -= this->T1_("a,m") * B_oo("Q,m,i");
    B_vo("Q,a,i") += B_vv("Q,a,e") * this->T1_("e,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in the last B_vv" << std::endl;
#endif
    // Qo2v2
    B_vv("Q,a,e") -= this->T1_("a,m") * this->riMoInts["bov"]("Q,m,e");

    // Dressed Fock
    // Qov
    TArray tmp_Q = TAManager::get().malloc<MatsT>("b");
    tmp_Q("Q") = this->riMoInts["bov"]("Q,m,e") * this->T1_("e,m");
    // Qo2v
    TArray X_Qoo = TAManager::get().malloc<MatsT>("boo");
    X_Qoo("Q,m,i") = this->riMoInts["bov"]("Q,m,e") * this->T1_("e,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in F_ov" << std::endl;
#endif
    // ov + Q + Qo
    F_ov("m,e")  = this->fockMatrix_ta["ov"]("m,e");
    F_ov("m,e") += tmp_Q("Q") * this->riMoInts["bov"]("Q,m,e");
    F_ov("m,e") -= X_Qoo("Q,m,n") * this->riMoInts["bov"]("Q,n,e");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in F_oo" << std::endl;
#endif
    // o2 + o2v + Qo2 + Qo3
    F_oo("m,i")  = this->fockMatrix_ta["oo"]("m,i");
    F_oo("m,i") += this->fockMatrix_ta["ov"]("m,e") * this->T1_("e,i");
    F_oo("m,i") += tmp_Q("Q") * B_oo("Q,m,i");
    F_oo("m,i") -= X_Qoo("Q,m,n") * B_oo("Q,n,i");
    TAManager::get().free("boo", std::move(X_Qoo));

    // Qov2
    TArray X_Qvo = TAManager::get().malloc<MatsT>("bvo");
    X_Qvo("Q,a,i") = B_vv("Q,a,f") * this->T1_("f,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in F_vv" << std::endl;
#endif
    // v2 + ov2 + Qv2 + Qov2
    F_vv("a,e")  = this->fockMatrix_ta["vv"]("a,e");
    F_vv("a,e") -= this->T1_("a,m") * this->fockMatrix_ta["ov"]("m,e");
    F_vv("a,e") += tmp_Q("Q") * B_vv("Q,a,e");
    F_vv("a,e") -= X_Qvo("Q,a,n") * this->riMoInts["bov"]("Q,n,e");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in F_vo" << std::endl;
#endif
    // ov + ov2 + o2v + o2v2 + Qov + Qo2v
    F_vo("a,i")  = this->fockMatrix_ta["vo"]("a,i");
    F_vo("a,i") += this->fockMatrix_ta["vv"]("a,e") * this->T1_("e,i");
    F_vo("a,i") -= this->T1_("a,m") * this->fockMatrix_ta["oo"]("m,i");
    F_vo("a,i") -= this->T1_("a,m") * this->fockMatrix_ta["ov"]("m,e") * this->T1_("e,i");
    F_vo("a,i") += tmp_Q("Q") * B_vo("Q,a,i");
    TAManager::get().free("b", std::move(tmp_Q));
    F_vo("a,i") -= X_Qvo("Q,a,m") * B_oo("Q,m,i");
    TAManager::get().free("bvo", std::move(X_Qvo));
  }

  template <typename MatsT>
  void DFCCSD<MatsT>::runConventional(){

    auto cc_start = tick();

    std::shared_ptr<DIISTA<MatsT> > diis = nullptr;
    if(this->ccSettings_.useDIIS){
      diis = std::make_shared<DIISTA<MatsT>>(this->ccSettings_.nDIIS);
    }

    this->CorrE = 0.0; // make public and not initialized in coupledcluster.hpp

    initIntermediates();

    initAmplitudes();
    if (this->ccSettings_.skipCC && this->ccSettings_.restart) return;

    MBExpansion<MatsT> T_old(this->T_);

    std::cout << std::setw(18) << std::left <<  "  CC Iterations";
    std::cout << std::setw(34) << std::left << "Corr. Energy (Eh)";
    std::cout << std::setw(19) << std::right << "\u0394Ec (Eh)";
    std::cout << std::setw(19) << std::right << "|\u0394T|";
    std::cout << std::endl;
    std::cout << std::setw(18) << std::left <<  "  -------------";
    std::cout << std::setw(34) << std::left << "-----------------";
    std::cout << std::setw(18) << std::right << "--------";
    std::cout << std::setw(18) << std::right << "----";
    std::cout << std::endl << std::endl;

    for (auto iter = 0; iter < this->ccSettings_.maxiter; iter++){

      T_old = this->T_;

      // build T1-transformed integrals
#ifdef DEBUG_DFCCSD
      std::cout << "I am now entering t1TransformIntegrals" << std::endl;
#endif
      t1TransformIntegrals();

#ifdef DEBUG_DFCCSD
      std::cout << "I am now entering updateT1" << std::endl;
#endif
      updateT1(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"));
#ifdef DEBUG_DFCCSD
      std::cout << "I am now entering updateT2" << std::endl;
#endif
      updateT2(T_old.get_tensor("TwoBody"));
#ifdef DEBUG_DFCCSD
      std::cout << "I am now entering DIIS" << std::endl;
#endif
      if (this->ccSettings_.useDIIS){
      // diis
        this->doDIIS(T_old, diis);
      } else {
        T_old.axpy(-1, this->T_);
      }

      MatsT Eold = this->CorrE;

      // Compute ERI_vvoo on the fly to reduce persistent o2v2 storage
      this->antiSymMoints["vvoo"] = TAManager::get().malloc<MatsT>("vvoo");
      this->antiSymMoints["vvoo"]("p,r,q,s")  = this->riMoInts["bvo"]("L,p,q") * this->riMoInts["bvo"]("L,r,s");
      this->antiSymMoints["vvoo"]("p,q,r,s") -= this->antiSymMoints["vvoo"]("p,q,s,r");
      this->getCorrEnergy();
      TAManager::get().free("vvoo", std::move(this->antiSymMoints["vvoo"]));

      this->intermediates_.E_cc = this->intermediates_.E_ref + std::real(this->CorrE);
      double dE = std::abs(this->CorrE - Eold);

      double dT = T_old.norm();

      std::cout << std::setprecision(12) << std::fixed;
      std::cout << "  Iteration "  << std::setw(6) << std::left << iter;
      std::cout << std::setw(34) << std::left << std::fixed << this->CorrE;
      std::cout << std::setw(18) << std::right << std::fixed << dE;
      std::cout << std::setw(18) << std::right << std::fixed << std::abs(dT);
      std::cout << std::endl;

      if (dE < this->ccSettings_.eConv and dT < this->ccSettings_.tConv) {

        std::cout << std::endl << "  CC Completed: Corr. E is "<< std::setw(18) << std::right
                  << std::setprecision(12) << this->CorrE << " Eh" << std::endl;
        std::cout << std::endl << "  CC Completed: Total E is "<< std::setw(18) << std::right
                  << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE << " Eh" << std::endl;
        std::cout << std::endl << "  CC Completed: Iteration total time "<< std::setw(10) << std::right
                  << std::setprecision(6) << tock(cc_start) << " s" << std::endl;

        if(this->ccSettings_.save){
          size_t size = this->T_.length();
          dcomplex * t_amp = CQMemManager::get().malloc<dcomplex>(size);
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

        // run CCSD(T) *after* saving amplitudes and corrE
        if (this->ccSettings_.pertT3) {
          std::cout << "  CCSD(T) with DF-CCSD is in progress, check back again soon!" << std::endl;
          CErr(std::string("DF-CCSD(T) NYI..."));
        }

        std::cout << BannerEnd << std::endl;

        break;
      }

      if(iter == this->ccSettings_.maxiter - 1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }
    }
  }

  /**
   * Build Eq. III(d.1) and Eq. III(d.2)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tau_and_tilde_tau() {
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
  void DFCCSD<MatsT>::build_tilde_Fae() {
    Fae_("a,e") = this->fockMatrix_ta["vv"]("a,e");
    Fae_("a,e") -= 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("a,m");
    Fae_("a,e") += this->T1_("f,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
    Fae_("a,e") -= 0.5 * tilde_tau_("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
//    Fae_("a,e") = Fae_("a,e") - this->fockMatrix_ta["vv_diag"]("a,e");
  }

  /**
   * Build Eq. III(a.2)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tilde_Fmi() {
    Fmi_("m,i") = this->fockMatrix_ta["oo"]("m,i");
    Fmi_("m,i") += 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("e,i");
    Fmi_("m,i") -= this->T1_("e,n") * conj(this->antiSymMoints["vooo"]("e,i,m,n"));
    Fmi_("m,i") += 0.5 * tilde_tau_("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n")); // TODO: better performance if we have this->antiSymMoints["oovv"]
//    Fmi_("m,i") = Fmi_("m,i") - this->fockMatrix_ta["oo_diag"]("m,i");
  }

  /**
   * Build Eq. III(a.3)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tilde_Fme() {
    Fme_("m,e") = this->fockMatrix_ta["ov"]("m,e");
    Fme_("m,e") += this->T1_("f,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.4)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tilde_Wmnij() {
    Wmnij_("m,n,i,j") = -this->T1_("e,j") * conj(this->antiSymMoints["vooo"]("e,i,m,n")); // TODO: better performance if we have this->antiSymMoints["ooov"]
    Wmnij_("m,n,i,j") -= Wmnij_("m,n,j,i");
    Wmnij_("m,n,i,j") += this->antiSymMoints["oooo"]("m,n,i,j");
    Wmnij_("m,n,i,j") += 0.25 * tau_("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.5)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tilde_Wabef() {
    Wabef_("a,b,e,f") = -this->T1_("b,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
    Wabef_("a,b,e,f") -= Wabef_("b,a,e,f");
    Wabef_("a,b,e,f") += this->antiSymMoints["vvvv"]("a,b,e,f");
    Wabef_("a,b,e,f") += 0.25 * tau_("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.6)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::build_tilde_Wmbej() {
    Wmbej_("m,b,e,j") = -this->antiSymMoints["vovo"]("b,m,e,j"); // TODO: better performance if we have this->antiSymMoints["voov"]
    Wmbej_("m,b,e,j") -= this->T1_("f,j") * conj(this->antiSymMoints["vvvo"]("e,f,b,m")); // TODO: better performance if we have this->antiSymMoints["ovvv"]
    Wmbej_("m,b,e,j") -= this->T1_("b,n") * conj(this->antiSymMoints["vooo"]("e,j,m,n"));
    TArray tmp = TAManager::get().malloc<MatsT>("vvoo");
    tmp("f,b,j,n") = 0.5 * this->T2_("f,b,j,n");
    tmp("f,b,j,n") += this->T1_("f,j") * this->T1_("b,n");
    Wmbej_("m,b,e,j") -= tmp("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    TAManager::get().free("vvoo", std::move(tmp));
  }

  /**
   * Build Eq. III(a.1-6)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::buildIntermediates() {
    TAManager &TAmanager = TAManager::get();

    // Literature intermediates, may need for EOM
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
    if (not Wmnij_.is_initialized()){
      Wmnij_ = TAmanager.malloc<MatsT>("oooo");
    }
    if (not Wabef_.is_initialized()){
      Wabef_ = TAmanager.malloc<MatsT>("vvvv");
    }
    if (not Wmbej_.is_initialized()){
      Wmbej_ = TAmanager.malloc<MatsT>("ovvo");
    }

    // Build slices of ERI needed to build intermediates
    this->antiSymMoints["oooo"] = TAmanager.malloc<MatsT>("oooo");
    this->antiSymMoints["oooo"]("m,n,i,j")  = this->riMoInts["boo"]("Q,m,i") * this->riMoInts["boo"]("Q,n,j");
    this->antiSymMoints["oooo"]("m,n,i,j") -= this->antiSymMoints["oooo"]("m,n,j,i");

    this->antiSymMoints["vooo"] = TAmanager.malloc<MatsT>("vooo");
    this->antiSymMoints["vooo"]("a,n,i,j")  = this->riMoInts["bvo"]("Q,a,i") * this->riMoInts["boo"]("Q,n,j");
    this->antiSymMoints["vooo"]("a,n,i,j") -= this->antiSymMoints["vooo"]("a,n,j,i");

    this->antiSymMoints["vvoo"] = TAmanager.malloc<MatsT>("vvoo");
    this->antiSymMoints["vvoo"]("a,b,i,j")  = this->riMoInts["bvo"]("Q,a,i") * this->riMoInts["bvo"]("Q,b,j");
    this->antiSymMoints["vvoo"]("a,b,i,j") -= this->antiSymMoints["vvoo"]("a,b,j,i");

    this->antiSymMoints["vovo"] = TAmanager.malloc<MatsT>("vovo");
    this->antiSymMoints["vovo"]("a,m,e,i")  = this->riMoInts["bvv"]("Q,a,e") * this->riMoInts["boo"]("Q,m,i");
    this->antiSymMoints["vovo"]("a,m,e,i") -= this->riMoInts["bvo"]("Q,a,i") * this->riMoInts["bov"]("Q,m,e");

    this->antiSymMoints["vvvo"] = TAmanager.malloc<MatsT>("vvvo");
    this->antiSymMoints["vvvo"]("a,b,e,j")  = this->riMoInts["bvv"]("Q,a,e") * this->riMoInts["bvo"]("Q,b,j");
    this->antiSymMoints["vvvo"]("a,b,e,j") -= this->antiSymMoints["vvvo"]("b,a,e,j");

    this->antiSymMoints["vvvv"] = TAmanager.malloc<MatsT>("vvvv");
    this->antiSymMoints["vvvv"]("a,b,e,f")  = this->riMoInts["bvv"]("Q,a,e") * this->riMoInts["bvv"]("Q,b,f");
    this->antiSymMoints["vvvv"]("a,b,e,f") -= this->antiSymMoints["vvvv"]("a,b,f,e");

    build_tau_and_tilde_tau();
    build_tilde_Fae();
    build_tilde_Fmi();
    build_tilde_Fme();
    build_tilde_Wmnij();
    build_tilde_Wabef();
    build_tilde_Wmbej();

    TAmanager.free("oooo",std::move(this->antiSymMoints["oooo"]));
    TAmanager.free("vooo",std::move(this->antiSymMoints["vooo"]));
    TAmanager.free("vvoo",std::move(this->antiSymMoints["vvoo"]));
    TAmanager.free("vovo",std::move(this->antiSymMoints["vovo"]));
    TAmanager.free("vvvo",std::move(this->antiSymMoints["vvvo"]));
    TAmanager.free("vvvv",std::move(this->antiSymMoints["vvvv"]));
    TA::get_default_world().gop.fence();
  }

  /**
   * Build Eq. I(a)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::updateT1(const TArray &T1_old, const TArray &T2_old) {
#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T1: Fvo" << std::endl;
#endif
    // ov
    this->T1_("a,i")  = F_vo("a,i");
#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T1: Fov*T2" << std::endl;
#endif
    // o2v2
    this->T1_("a,i") += F_ov("m,e") * T2_old("e,a,m,i");
#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T1: Wooov*T2 + Wvovv*T2" << std::endl;
#endif
    // Qo2v2 + Qo2v + Qov2
    TArray tmp_Qvo = TAManager::get().malloc<MatsT>("bvo");
    tmp_Qvo("Q,a,i") = this->riMoInts["bov"]("Q,m,e") * T2_old("a,e,i,m");
    this->T1_("a,i") -= B_oo("Q,m,i") * tmp_Qvo("Q,a,m");
    this->T1_("a,i") += B_vv("Q,a,e") * tmp_Qvo("Q,e,i");
    TAManager::get().free("bvo", std::move(tmp_Qvo));

    this->T1_("a,i") = T1_old("a,i") + this->T1_("a,i") * this->Dai_("a,i");
  }

  /**
   * Build Eq. I(b)
   */
  template <typename MatsT>
  void DFCCSD<MatsT>::updateT2(const TArray &T2_old) {
#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: Wvvoo" << std::endl;
#endif
    // 2*o2v2
    this->T2_("a,b,i,j")  = B_vo("Q,a,i") * B_vo("Q,b,j");
    this->T2_("a,b,i,j") -= B_vo("Q,a,j") * B_vo("Q,b,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: Wvvvv * T2" << std::endl;
#endif

/*
 *  Stephen comments:
 *  Currently, the loop-based PPL term is correctly implemented
 *  but only works for serial run. MPI-parallelization requires
 *  this->T2_ assignment to be run separately and is not trivial
 *  in TiledArray (at least to me). So PPL is run as usual, i.e.
 *  by building the W_vvvv intermediate on the fly.
 */
//    // Do W_vvvv contractions for a<=b as a loop
//    // Define TileRanges
//    size_t o_tr = (this->intermediates_.nOcc + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
//    size_t v_tr = (this->intermediates_.nVir + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
//    size_t Q_tr = (this->intermediates_.nRI + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
//
//    auto compute_o2v4_ab = [&,this](size_t &a, size_t &b){
//      // Define blocks
//      typedef std::vector<size_t> block;
//      block   Ba_lower = {   0,    a,    0};
//      block   Ba_upper = {Q_tr,  a+1, v_tr};
//      block   Bb_lower = {   0,    b,    0};
//      block   Bb_upper = {Q_tr,  b+1, v_tr};
//      block res2_lower = {   a,    b,    0,    0};
//      block res2_upper = { a+1,  b+1, o_tr, o_tr};
//
//      // Get appropriate chunks of RI integrals
//      auto Ba = B_vv("Q,a,e").block(Ba_lower,Ba_upper);
//      auto Bb = B_vv("Q,b,f").block(Bb_lower,Bb_upper);
//
//      // Construct small vvvv object
//      TArray W_vvvv;
//      W_vvvv("a,b,e,f") = Ba * Bb;
//
//      // Add contribution
//      this->T2_("a,b,i,j").block(res2_lower,res2_upper) += W_vvvv("a,b,e,f") * T2_old("e,f,i,j");
//
//      // If a!=b, also add (b,a,i,j) contribution
//      if (a != b) {
//        res2_lower = {   b,   a,    0,    0};
//        res2_upper = { b+1, a+1, o_tr, o_tr};
//
//        TArray W_vvvv2;
//        W_vvvv2("b,a,f,e") = Bb * Ba;
//
//        this->T2_("b,a,i,j").block(res2_lower, res2_upper) += W_vvvv2("b,a,f,e") * T2_old("f,e,i,j");
//      }
//
//    };
//
//    // Loop over a>=b
//    for (size_t a = 0; a < v_tr; ++a)
//      for (size_t b = 0; b <= a; ++b){
//#ifdef DEBUG_DFCCSD
//        std::cout << "  Running Wvvvv at a = " << a << " and b = " << b;
//#endif
//        compute_o2v4_ab(a,b);
//    };

    // Construct W_vvvv
    // Qv4 + o2v4
    TArray W_vvvv = TAManager::get().malloc<MatsT>("vvvv");
    W_vvvv("a,b,e,f")  = B_vv("Q,a,e") * B_vv("Q,b,f");
    this->T2_("a,b,i,j") += W_vvvv("a,b,e,f") * T2_old("e,f,i,j");
    TAManager::get().free("vvvv", std::move(W_vvvv));

    // Qo2v2
    TArray W_oovv = TAManager::get().malloc<MatsT>("oovv");
    W_oovv("m,n,e,f")  = this->riMoInts["bov"]("Q,m,e") * this->riMoInts["bov"]("Q,n,f");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: (Woooo + Woovv * T2) * T2" << std::endl;
#endif
    // Qo4 + Qo4v2 + Qo4v2
    TArray W_oooo = TAManager::get().malloc<MatsT>("oooo");
    W_oooo("m,n,i,j")  = B_oo("Q,m,i") * B_oo("Q,n,j");
    W_oooo("m,n,i,j") += 0.5 * W_oovv("m,n,e,f") * T2_old("e,f,i,j");
    this->T2_("a,b,i,j") += W_oooo("m,n,i,j") * T2_old("a,b,m,n");
    TAManager::get().free("oooo", std::move(W_oooo));

#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: (Wvoov + Woovv * T2) * T2" << std::endl;
#endif
    // 2*Qo2v2 + 2*o3v3 + 4*o2v2
    TArray W_voov = TAManager::get().malloc<MatsT>("voov");
    W_voov("a,m,i,e")  = B_vo("Q,a,i") * this->riMoInts["bov"]("Q,m,e");
    W_voov("a,m,i,e") -= B_oo("Q,m,i") * B_vv("Q,a,e");
    W_voov("a,m,i,e") += 0.5 * (W_oovv("m,n,e,f") - W_oovv("m,n,f,e")) * T2_old("f,a,n,i");
    TArray tmp_t2 = TAManager::get().malloc<MatsT>("vvoo");
    tmp_t2("a,b,i,j") = W_voov("a,m,i,e") * T2_old("e,b,m,j");
    TAManager::get().free("voov", std::move(W_voov));
    this->T2_("a,b,i,j") += tmp_t2("a,b,i,j");
    this->T2_("a,b,i,j") -= tmp_t2("a,b,j,i");
    this->T2_("a,b,i,j") -= tmp_t2("b,a,i,j");
    this->T2_("a,b,i,j") += tmp_t2("b,a,j,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: (Foo + Woovv * T2) * T2" << std::endl;
#endif
    // 2*o3v2 + 2*o2v2
    F_oo("m,i") += W_oovv("m,n,e,f") * T2_old("e,f,i,n");
    tmp_t2("a,b,i,j") = F_oo("m,i") * T2_old("a,b,m,j");
    this->T2_("a,b,i,j") -= tmp_t2("a,b,i,j");
    this->T2_("a,b,i,j") += tmp_t2("a,b,j,i");

#ifdef DEBUG_DFCCSD
    std::cout << "I am in update_T2: (Fvv + Woovv * T2) * T2" << std::endl;
#endif
    // 2*o2v3 + 2*o2v2
    F_vv("a,e") -= W_oovv("m,n,e,f") * T2_old("a,f,m,n");
    TAManager::get().free("oovv", std::move(W_oovv));
    tmp_t2("a,b,i,j") = F_vv("a,e") * T2_old("e,b,i,j");
    this->T2_("a,b,i,j") += tmp_t2("a,b,i,j");
    this->T2_("a,b,i,j") -= tmp_t2("b,a,i,j");
    TAManager::get().free("vvoo", std::move(tmp_t2));

    this->T2_("a,b,i,j") = T2_old("a,b,i,j") + this->T2_("a,b,i,j") * this->Dabij_("a,b,i,j");
  }

  template <typename MatsT>
  DFCCSD<MatsT>::~DFCCSD() {

    cleanMemory();


  }
  
  template <typename MatsT>
  void DFCCSD<MatsT>::cleanMemory(){
    TAManager &TAmanager = TAManager::get();
    //Free all the DF-specific intermediates here!
    if (F_ov) TAmanager.free("ov", std::move(F_ov), true);
    if (F_oo) TAmanager.free("oo", std::move(F_oo), true);
    if (F_vv) TAmanager.free("vv", std::move(F_vv), true);
    if (F_vo) TAmanager.free("vo", std::move(F_vo), true);

    if (B_oo) TAmanager.free("boo", std::move(B_oo), true);
    if (B_vv) TAmanager.free("bvv", std::move(B_vv), true);
    if (B_vo) TAmanager.free("bvo", std::move(B_vo), true);

//    if (Fae_) TAmanager.free("vv", std::move(Fae_), true);
//    if (Fmi_) TAmanager.free("oo", std::move(Fmi_), true);
    if (tilde_tau_) TAmanager.free("vvoo", std::move(tilde_tau_), true);

    TAmanager.discard_cache();

  }
};
