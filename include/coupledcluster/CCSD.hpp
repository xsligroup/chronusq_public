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

  template <typename MatsT, typename IntsT>
  CCSD<MatsT,IntsT>::CCSD(const SafeFile &savFile,
                          CCIntermediates<MatsT> &intermediates,
                          const CoupledClusterSettings &ccSettings):
      savFile_(savFile),
      vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
      intermediates_(intermediates), ccSettings_(ccSettings),
      fockMatrix_ta(intermediates.fockMatrix),
      antiSymMoints(intermediates.antiSymMoInts),
      T_(*intermediates.T),
      T1_(intermediates.T->oneBody()),
      T2_(intermediates.T->twoBody()),
      tau_(intermediates.tau),
      tilde_tau_(intermediates.tilde_tau),
      Fae_wDiag_(intermediates.F_ae),
      Fmi_wDiag_(intermediates.F_mi),
      Fme_(intermediates.F_me),
      Wmnij_(intermediates.W_mnij),
      Wabef_(intermediates.W_abef),
      Wmbej_(intermediates.W_mbej),
      Dai_(intermediates.D_ai),
      Dabij_(intermediates.D_abij){}

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::run() {
    runConventional();
  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initIntermediates() {
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
    if (not Fae_wDiag_.is_initialized()){
      Fae_wDiag_ = TAmanager.malloc<MatsT>("vv");
    }
    if (not Fmi_.is_initialized()){
      Fmi_ = TAmanager.malloc<MatsT>("oo");
    }
    if (not Fmi_wDiag_.is_initialized()){
      Fmi_wDiag_ = TAmanager.malloc<MatsT>("oo");
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

  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initIntermediates_xll() {
    TAManager &TAmanager = TAManager::get();
    if (not this->A_1.is_initialized()){
      this->A_1 = TAmanager.malloc<MatsT>("oo");
    }
    if (not this->A_2.is_initialized()){
      this->A_2 = TAmanager.malloc<MatsT>("vv");
    }
    if (not this->A_3.is_initialized()){
      this->A_3 = TAmanager.malloc<MatsT>("ov");
    }
    if (not this->A_4.is_initialized()){
      this->A_4 = TAmanager.malloc<MatsT>("oovo");
    }
  
    if (not this->B_1.is_initialized()){
      this->B_1 = TAmanager.malloc<MatsT>("oo");
    }
    if (not this->B_2.is_initialized()){
      this->B_2 = TAmanager.malloc<MatsT>("vv");
    }
    if (not this->B_3.is_initialized()){
      this->B_3 = TAmanager.malloc<MatsT>("ovoo");
    }
    if (not this->B_5.is_initialized()){
      this->B_5 = TAmanager.malloc<MatsT>("oooo");
    }
    if (not this->B_6.is_initialized()){
      this->B_6 = TAmanager.malloc<MatsT>("ovvo");
    }
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::initAmplitudes() {
    T_.scale(0.0);
  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::doDIIS(EOMCCSDVector<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis ){

      // give solution vector to diis
      diis->WriteVector(T_);

      // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
      T_old.scale(-1.0);
      T_old.axpy(1.0, T_);

      //set error vector in DIIS
      diis->WriteErrorVector(T_old);

      // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
      diis->Extrapolate(T_);

  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::printBanner(double Eref) const {

    std::cout << BannerTop << std::endl;
    std::cout << "Coupled Cluster (CC) Settings:" << std::endl << std::endl;

    std::cout << std::setw(45) << std::left << "  Type:"
    << "Singles and Doubles" << std::endl;

    std::cout << std::setw(45) << std::left << "  Occupied Orbitals:"
    << TAManager::get().getRange(oLabel_).extent() << std::endl;
    std::cout << std::setw(45) << std::left << "  Virtual Orbitals:"
    << TAManager::get().getRange(vLabel_).extent() << std::endl;
    std::cout << std::setw(45) << std::left << "  Reference energy:"
              << std::setprecision(10) << std::fixed << Eref << std::endl;

    std::cout << std::setprecision(6) << std::scientific;

    std::cout << std::setw(45) << std::left << "  Energy Convergence Tolerance:"
    << this->ccSettings_.eConv << std::endl;
    std::cout << std::setw(45) << std::left << "  Amplitude Convergence Tolerance:"
    << this->ccSettings_.tConv << std::endl;

    std::cout << std::setw(45) << std::left << "  Direct Inversion of Iterative Subspace:";

    if ( this->ccSettings_.useDIIS ) {
      std::cout << "On" << std::endl;
      std::cout << std::left << "    * DIIS will track up to " << this->ccSettings_.nDIIS
      << " previous iterations" << std::endl;
    }
    else
      std::cout << "Off" << std::endl;

    std::cout << std::endl;

    std::pair<double, char> mem_postfix = memSize(intermediates_.estimate_mem_peak(ccSettings_.useDIIS ? ccSettings_.nDIIS : 0));
    std::cout << std::setw(45) << std::left << "  Estimated TiledArray memory requirement: " << std::fixed << std::setprecision(1)
    << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    std::cout << BannerMid << std::endl << std::endl;

  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::runConventional(){

    auto cc_start = tick();
    TAManager &TAmanager = TAManager::get();

    std::shared_ptr<DIISTA<MatsT> > diis = nullptr;
    if(this->ccSettings_.useDIIS){
      diis = std::make_shared<DIISTA<MatsT>>(this->ccSettings_.nDIIS);
    }
//#define XLL_ALG
#ifdef XLL_ALG
    initIntermediates_xll();
#else
    initIntermediates();
#endif
    initAmplitudes();

    EOMCCSDVector<MatsT> T_old(vLabel_, oLabel_);

#ifdef XLL_ALG
    TArray T2_old_copy = TAmanager.malloc<MatsT>("vvoo");
#endif

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

      T_old = T_;

#ifdef XLL_ALG
      T2_old_copy("p,q,r,s") = this->T2_("p,q,r,s");

      // Xiaolin's algorithm
      formA_1();
      formA_2();
      formA_3();
      formA_4();

      formB_1();
      formB_2();
      formB_3();
      formB_5();
      formB_6();

      updateT1_xxl(T_old.oneBody(), T2_old_copy);
      //T2_old_copy will be modified in updateT2_xxl()!
      updateT2_xxl(T_old.oneBody(), T2_old_copy);
#else

      buildIntermediates();
#ifdef DEBUG_CCSD
      std::cout << "tau_:" << tau_ << std::endl;
      std::cout << "tilde_tau_:" << tilde_tau_ << std::endl;
      std::cout << "Fae_:" << Fae_ << std::endl;
      std::cout << "Fmi_:" << Fmi_ << std::endl;
      std::cout << "Fme_:" << Fme_ << std::endl;
      std::cout << "Wmnij_:" << Wmnij_ << std::endl;
      std::cout << "Wabef_:" << Wabef_ << std::endl;
      std::cout << "Wmbej_:" << Wmbej_ << std::endl;
#endif

      updateT1(T_old.oneBody(), T_old.twoBody());
      updateT2(T_old.oneBody(), T_old.twoBody());
#endif
#ifdef DEBUG_CCSD
      std::cout << "T1_:" << T1_ << std::endl;
      std::cout << "T2_:" << T2_ << std::endl;
#endif
      
      if (this->ccSettings_.useDIIS){
      // diis
        doDIIS(T_old, diis);
      } else {
        T_old.axpy(-1, T_);
      }

      MatsT Eold = CorrE;
      getCorrEnergy();
      double dE = std::abs(CorrE - Eold);

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
                  << std::setprecision(12) << intermediates_.E_ref + CorrE << " Eh" << std::endl;
        std::cout << std::endl << "  CC Completed: Iteration total time "<< std::setw(10) << std::right
                  << std::setprecision(6) << tock(cc_start) << " s" << std::endl;
        
        if (savFile_.exists()) {
          double realCoreE = std::real(this->CorrE);
          savFile_.safeWriteData("/CC/CORRELATION_ENERGY",&realCoreE, {1});
        }
        std::cout << bannerEnd << std::endl;

        printAnalysis();
        
        std::cout << BannerEnd << std::endl;
        
        break;
      }

      if(iter == ccSettings_.maxiter - 1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(ccSettings_.maxiter) + " steps." );
      }
    }

#ifdef XLL_ALG
    TAmanager.free("vvoo", std::move(T2_old_copy));
#endif

  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::getCorrEnergy() {
    MatsT CorrEOneBody = fockMatrix_ta["ov"]("i,a").dot(T1_("a,i")).get();
    MatsT CorrETwoBodyT2 = conj(antiSymMoints["vvoo"]("c,d,k,l")).dot(T2_("c,d,k,l")).get();
    MatsT CorrETwoBodyT1 = conj(antiSymMoints["vvoo"]("c,d,k,l")).dot(T1_("c,k") * T1_("d,l")).get();
    TA::get_default_world().gop.fence();
    this->CorrE = CorrEOneBody + 0.25 * (CorrETwoBodyT2 + 2.0 * CorrETwoBodyT1);
  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::printAnalysis() {
    std::cout << bannerTop << std::endl;
    std::cout << "Coupled Cluster Wave Function Analysis:" << std::endl;
    std::cout << std::setw(22) << "   max(|t1|)  " << std::setprecision(4) << abs_max(T1_).get() << std::endl;
    std::cout << std::setw(22) << "   max(|t2|)  " << std::setprecision(4) << abs_max(T2_).get() << std::endl;
    double normT1Sq = squared_norm(T1_).get();
    double normT2Sq = squared_norm(T2_).get();
    std::cout << std::setw(22) << "   sum of t1 weights  " << std::setprecision(4) << normT1Sq << std::endl;
    std::cout << std::setw(22) << "   sum of t2 weights  " << std::setprecision(4) << normT2Sq << std::endl;
    std::cout << std::setw(22) << "   T1 diagnostic  " << std::setprecision(4) << sqrt(normT1Sq / intermediates_.nOcc) << std::endl;
    std::cout << std::setw(22) << "   T2 diagnostic  " << std::setprecision(4) << sqrt(normT2Sq / intermediates_.nOcc) << std::endl;
    std::cout << bannerEnd << std::endl;
  }

  /**
   * Build Eq. III(d.1) and Eq. III(d.2)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tau_and_tilde_tau() {
    tau_("a,b,i,j") = 0.5 * T1_("a,i") * T1_("b,j");
    tau_("a,b,i,j") -= tau_("b,a,i,j");
    tau_("a,b,i,j") -= tau_("a,b,j,i");

    tilde_tau_("a,b,i,j") = 0.5 * tau_("a,b,i,j") + T2_("a,b,i,j");
    tau_("a,b,i,j") += T2_("a,b,i,j");
  }

  /**
   * Build Eq. III(a.1)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Fae() {
    Fae_wDiag_("a,e") = fockMatrix_ta["vv"]("a,e");
#ifdef DEBUG_CCSD
    std::cout << "Fae_1:" << Fae_ << std::endl;
#endif

    Fae_wDiag_("a,e") -= 0.5 * fockMatrix_ta["ov"]("m,e") * T1_("a,m");
#ifdef DEBUG_CCSD
    std::cout << "Fae_2:" << Fae_ << std::endl;
#endif
    Fae_wDiag_("a,e") += T1_("f,m") * conj(antiSymMoints["vvvo"]("e,f,a,m"));
#ifdef DEBUG_CCSD
    std::cout << "Fae_3:" << Fae_ << std::endl;
#endif

    Fae_wDiag_("a,e") -= 0.5 * tilde_tau_("a,f,m,n") * conj(antiSymMoints["vvoo"]("e,f,m,n"));
#ifdef DEBUG_CCSD
    std::cout << "Fae_4:" << Fae_ << std::endl;
#endif

    Fae_("a,e") = Fae_wDiag_("a,e") - fockMatrix_ta["vv_diag"]("a,e");
  }

  /**
   * Build Eq. III(a.2)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Fmi() {
    Fmi_wDiag_("m,i") = fockMatrix_ta["oo"]("m,i");

    Fmi_wDiag_("m,i") += 0.5 * fockMatrix_ta["ov"]("m,e") * T1_("e,i");
    Fmi_wDiag_("m,i") -= T1_("e,n") * conj(antiSymMoints["vooo"]("e,i,m,n"));

    Fmi_wDiag_("m,i") += 0.5 * tilde_tau_("e,f,i,n") * conj(antiSymMoints["vvoo"]("e,f,m,n")); // TODO: better performance if we have antiSymMoints["oovv"]

    Fmi_("m,i") = Fmi_wDiag_("m,i") - fockMatrix_ta["oo_diag"]("m,i");
  }

  /**
   * Build Eq. III(a.3)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Fme() {
    Fme_("m,e") = fockMatrix_ta["ov"]("m,e");
    Fme_("m,e") += T1_("f,n") * conj(antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.4)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Wmnij() {
    Wmnij_("m,n,i,j") = -T1_("e,j") * conj(antiSymMoints["vooo"]("e,i,m,n")); // TODO: better performance if we have antiSymMoints["ooov"]
    Wmnij_("m,n,i,j") -= Wmnij_("m,n,j,i");

    Wmnij_("m,n,i,j") += antiSymMoints["oooo"]("m,n,i,j");

    Wmnij_("m,n,i,j") += 0.25 * tau_("e,f,i,j") * conj(antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.5)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Wabef() {
    Wabef_("a,b,e,f") = -T1_("b,m") * conj(antiSymMoints["vvvo"]("e,f,a,m"));
    Wabef_("a,b,e,f") -= Wabef_("b,a,e,f");

    Wabef_("a,b,e,f") += antiSymMoints["vvvv"]("a,b,e,f");

    Wabef_("a,b,e,f") += 0.25 * tau_("a,b,m,n") * conj(antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.6)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::build_tilde_Wmbej() {
    Wmbej_("m,b,e,j") = -antiSymMoints["vovo"]("b,m,e,j"); // TODO: better performance if we have antiSymMoints["voov"]

    Wmbej_("m,b,e,j") -= T1_("f,j") * conj(antiSymMoints["vvvo"]("e,f,b,m")); // TODO: better performance if we have antiSymMoints["ovvv"]
    Wmbej_("m,b,e,j") -= T1_("b,n") * conj(antiSymMoints["vooo"]("e,j,m,n"));

    TArray tmp = TAManager::get().malloc<MatsT>("vvoo");
    tmp("f,b,j,n") = 0.5 * T2_("f,b,j,n");
    tmp("f,b,j,n") += T1_("f,j") * T1_("b,n");
    Wmbej_("m,b,e,j") -= tmp("f,b,j,n") * conj(antiSymMoints["vvoo"]("e,f,m,n"));
    TAManager::get().free("vvoo", std::move(tmp));
  }

  /**
   * Build Eq. III(a.1-6)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::buildIntermediates() {
    build_tau_and_tilde_tau();
    build_tilde_Fae();
    build_tilde_Fmi();
    build_tilde_Fme();
    build_tilde_Wmnij();
    build_tilde_Wabef();
    build_tilde_Wmbej();
  }

  /**
   * Build Eq. I(a)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT1(const TArray T1_old, const TArray T2_old) {
    T1_("a,i") = fockMatrix_ta["vo"]("a,i");
    T1_("a,i") += Fae_("a,e") * T1_old("e,i");
    T1_("a,i") -= T1_old("a,m") * Fmi_("m,i");
    T1_("a,i") += Fme_("m,e") * T2_old("a,e,i,m");
    T1_("a,i") -= T1_old("e,m") * antiSymMoints["vovo"]("a,m,e,i");
    T1_("a,i") += 0.5 * T2_old("e,f,i,m") * conj(antiSymMoints["vvvo"]("e,f,a,m")); // TODO: better performance if we have antiSymMoints["vvov"]
    T1_("a,i") -= 0.5 * T2_old("a,e,m,n") * conj(antiSymMoints["vooo"]("e,i,n,m")); // TODO: better performance if we have antiSymMoints["ovoo"]

    T1_("a,i") = T1_("a,i") * Dai_("a,i");
  }

  /**
   * Build Eq. I(b)
   */
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT2(const TArray T1_old, const TArray T2_old) {
    TAManager &TAmanager = TAManager::get();

    T2_("a,b,i,j") = antiSymMoints["vvoo"]("a,b,i,j");

    TArray TMPbe = TAmanager.malloc<MatsT>("vv");
    TMPbe("b,e") = Fae_("b,e");
    TMPbe("b,e") -= 0.5 * T1_old("b,m") * Fme_("m,e");
    TArray Pabij = TAmanager.malloc<MatsT>("vvoo");
    Pabij("a,b,i,j") = T2_old("a,e,i,j") * TMPbe("b,e");
    T2_("a,b,i,j") += Pabij("a,b,i,j");
    T2_("a,b,i,j") -= Pabij("b,a,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_1:" << T2_ << std::endl;
#endif
    TAmanager.free("vv", std::move(TMPbe));

    TArray TMPmj = TAmanager.malloc<MatsT>("oo");
    TMPmj("m,j") = Fmi_("m,j");
    TMPmj("m,j") += 0.5 * T1_old("e,j") * Fme_("m,e");
    Pabij("a,b,i,j") = T2_old("a,b,i,m") * TMPmj("m,j");
    T2_("a,b,i,j") -= Pabij("a,b,i,j");
    T2_("a,b,i,j") += Pabij("a,b,j,i");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_2:" << T2_ << std::endl;
#endif
    TAmanager.free("oo", std::move(TMPmj));

    T2_("a,b,i,j") += 0.5 * tau_("a,b,m,n") * Wmnij_("m,n,i,j");
    T2_("a,b,i,j") += 0.5 * tau_("e,f,i,j") * Wabef_("a,b,e,f");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_3:" << T2_ << std::endl;
#endif

    TArray TMPmbij = TAmanager.malloc<MatsT>("ovoo");
    TMPmbij("m,b,i,j") = T1_old("e,i") * antiSymMoints["vovo"]("b,m,e,j");
    Pabij("a,b,i,j") = T1_old("a,m") * TMPmbij("m,b,i,j");
    Pabij("a,b,i,j") += T2_old("a,e,i,m") * Wmbej_("m,b,e,j");
    Pabij("a,b,i,j") -= Pabij("a,b,j,i");
    Pabij("a,b,i,j") -= Pabij("b,a,i,j");
    T2_("a,b,i,j") += Pabij("a,b,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_4:" << T2_ << std::endl;
#endif
    TAmanager.free("ovoo", std::move(TMPmbij));

    Pabij("a,b,i,j") = T1_old("e,i") * antiSymMoints["vvvo"]("a,b,e,j");
    T2_("a,b,i,j") += Pabij("a,b,i,j");
    T2_("a,b,i,j") -= Pabij("a,b,j,i");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_5:" << T2_ << std::endl;
#endif

    Pabij("a,b,i,j") = -T1_old("a,m") * antiSymMoints["vooo"]("b,m,i,j"); // TODO: change sign
    T2_("a,b,i,j") -= Pabij("a,b,i,j");
    T2_("a,b,i,j") += Pabij("b,a,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_6:" << T2_ << std::endl;
#endif

    T2_("a,b,i,j") = T2_("a,b,i,j") * Dabij_("a,b,i,j");
    TAmanager.free("vvoo", std::move(Pabij));
  }

  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT1_xxl(const TArray T1_old, const TArray T2_old) {
  
    this->T1_("a,i") =   fockMatrix_ta["vo"]("a,i");
    this->T1_("a,i") +=  - this->antiSymMoints["vovo"]("a,k,c,i") * T1_old("c,k");
    this->T1_("a,i") +=  0.5  *   conj(this->antiSymMoints["vvvo"]("d,c,a,k")) * T2_old("c,d,k,i");  
    this->T1_("a,i") += A_1("k,i") * T1_old("a,k");   
    this->T1_("a,i") += A_2("a,c") * T1_old("c,i");
    this->T1_("a,i") += A_3("k,c") * T2_old("c,a,k,i");
    this->T1_("a,i") += A_4("k,l,c,i") * T2_old("c,a,k,l");
  
    this->T1_("a,i") = this->T1_("a,i") * Dai_("a,i");
  
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_1() {
  
    A_1("k,i") =   conj(this->antiSymMoints["vooo"]("c,i,k,l")) * this->T1_("c,l"); 
    A_1("k,i") +=  - 0.5  *   conj(this->antiSymMoints["vvoo"]("c,d,l,k")) * this->T2_("c,d,l,i"); 
    A_1("k,i") += -fockMatrix_ta["oo"]("k,i") + fockMatrix_ta["oo_diag"]("k,i");
    TA::TArray<MatsT> A_11 = TAManager::get().malloc<MatsT>("ov");
    A_11("k,c") = - conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T1_("d,l");
    A_11("k,c") += -fockMatrix_ta["ov"]("k,c");
    A_1("k,i") += A_11("k,c") * this->T1_("c,i");
    TAManager::get().free("ov", std::move(A_11));
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_2() {
    A_2("a,c") = - conj(this->antiSymMoints["vvvo"]("d,c,a,k")) * this->T1_("d,k");
    A_2("a,c") +=   fockMatrix_ta["vv"]("a,c") - fockMatrix_ta["vv_diag"]("a,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_3() {
    A_3("k,c") =    conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T1_("d,l");
    A_3("k,c") +=   fockMatrix_ta["ov"]("k,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formA_4() {
    A_4("k,l,c,i") = - 0.5  *   conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T1_("d,i");
    A_4("k,l,c,i") += - 0.5  *   conj(this->antiSymMoints["vooo"]("c,i,k,l"));
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::updateT2_xxl(const TArray T1_old, TArray T2_old) {
    this->T2_("a,b,i,j") =    this->antiSymMoints["vvoo"]("a,b,i,j");
    this->T2_("a,b,i,j") +=   B_1("k,i") * T2_old("a,b,k,j");
    this->T2_("a,b,i,j") += -B_1("k,j") * T2_old("a,b,k,i");
    this->T2_("a,b,i,j") +=   B_2("b,c") * T2_old("c,a,i,j");
    this->T2_("a,b,i,j") += -B_2("a,c") * T2_old("c,b,i,j");
    this->T2_("a,b,i,j") +=   B_3("k,a,i,j") * T1_old("b,k");
    this->T2_("a,b,i,j") += -B_3("k,b,i,j") * T1_old("a,k");
    this->T2_("a,b,i,j") += B_5("k,l,i,j") * T2_old("a,b,k,l");
    this->T2_("a,b,i,j") +=   B_6("k,a,c,j") * T2_old("c,b,k,i");
    this->T2_("a,b,i,j") += -B_6("k,b,c,j") * T2_old("c,a,k,i");
    this->T2_("a,b,i,j") +=   B_6("k,b,c,i") * T2_old("c,a,k,j");
    this->T2_("a,b,i,j") += -B_6("k,a,c,i") * T2_old("c,b,k,j");
    this->T2_("a,b,i,j") +=  this->antiSymMoints["vvvo"]("a,b,c,j") * T1_old("c,i");
    this->T2_("a,b,i,j") += -this->antiSymMoints["vvvo"]("a,b,c,i") * T1_old("c,j");
    //reuse the T2_old container
    T2_old("c,d,i,j") += T1_old("c,i") * T1_old("d,j");
    T2_old("c,d,i,j") += - T1_old("c,j") * T1_old("d,i");
    this->T2_("a,b,i,j") += 0.5 * this->antiSymMoints["vvvv"]("a,b,c,d") * T2_old("c,d,i,j");
  
    this->T2_("a,b,i,j") = this->T2_("a,b,i,j") * Dabij_("a,b,i,j");
  
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_1() {
    B_1("k,i") =    conj(this->antiSymMoints["vooo"]("c,i,k,l")) * this->T1_("c,l");
    B_1("k,i") +=  0.5  *   conj(this->antiSymMoints["vvoo"]("d,c,k,l")) * this->T2_("d,c,l,i");
    B_1("k,i") += -fockMatrix_ta["oo"]("k,i") + fockMatrix_ta["oo_diag"]("k,i");
    TA::TArray<MatsT> B_11 = TAManager::get().malloc<MatsT>("ov");
    B_11("k,c") = conj(this->antiSymMoints["vvoo"]("d,c,k,l")) * this->T1_("d,l");
    B_11("k,c") += -fockMatrix_ta["ov"]("k,c");
    B_1("k,i") += B_11("k,c") * this->T1_("c,i");
    TAManager::get().free("ov", std::move(B_11));
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_2() {
    B_2("b,c") =    conj(this->antiSymMoints["vvvo"]("d,c,b,k")) * this->T1_("d,k");
    B_2("b,c") +=  - 0.5  *   conj(this->antiSymMoints["vvoo"]("c,d,l,k")) * this->T2_("d,b,l,k"); 
    B_2("b,c") += -fockMatrix_ta["vv"]("b,c") + fockMatrix_ta["vv_diag"]("b,c");
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_3() {
    TAManager &TAmanager = TAManager::get();
    B_3("k,a,i,j") = 0.5  *   conj(this->antiSymMoints["vvvo"]("c,d,a,k")) * this->T2_("d,c,i,j"); 
    B_3("k,a,i,j") +=    this->antiSymMoints["vooo"]("a,k,j,i");
    TA::TArray<MatsT> B_31 = TAmanager.malloc<MatsT>("ov");
    B_31("k,c") = 0.0 * B_31("k,c");
    TA::TArray<MatsT> B_32 = TAmanager.malloc<MatsT>("oooo");
    B_32("k,l,i,j") = 0.0 * B_32("k,l,i,j");
    TA::TArray<MatsT> B_321 = TAmanager.malloc<MatsT>("oovo");
    B_321("k,l,c,j") = 0.0 * B_321("k,l,c,j");
    TA::TArray<MatsT> B_33 = TAmanager.malloc<MatsT>("ovvo");
    B_33("k,a,c,j") = 0.0 * B_33("k,a,c,j");
    B_321("k,l,c,j") += - 0.25   *  conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T1_("d,j");
    B_321("k,l,c,j") += - 0.5  *   conj(this->antiSymMoints["vooo"]("c,j,k,l"));
    B_31("k,c") +=    conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T1_("d,l");
    B_31("k,c") +=   fockMatrix_ta["ov"]("k,c");
    B_32("k,l,i,j") += -0.25   *  conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T2_("c,d,i,j");
    B_32("k,l,i,j") += - 0.5  *   this->antiSymMoints["oooo"]("k,l,i,j");
    B_32("k,l,i,j") +=   B_321("k,l,c,j") * this->T1_("c,i");
    B_32("k,l,i,j") += -B_321("k,l,c,i") * this->T1_("c,j");
    TAmanager.free("oovo", std::move(B_321));
    B_33("k,a,c,j") +=  0.5  *   conj(this->antiSymMoints["vvvo"]("d,c,a,k")) * this->T1_("d,j");
    B_33("k,a,c,j") +=   - this->antiSymMoints["vovo"]("a,k,c,j");
    TA::TArray<MatsT> B_34 = TAmanager.malloc<MatsT>("oovo");
    B_34("l,k,c,j") = 0.0 * B_34("l,k,c,j");
    B_34("l,k,c,j") +=   - conj(this->antiSymMoints["vvoo"]("c,d,l,k")) * this->T1_("d,j"); 
    B_34("l,k,c,j") +=   - conj(this->antiSymMoints["vooo"]("c,j,l,k"));  
    B_3("k,a,i,j") += B_31("k,c") * this->T2_("c,a,i,j");
    B_3("k,a,i,j") += B_32("k,l,i,j") * this->T1_("a,l");
    B_3("k,a,i,j") +=   B_33("k,a,c,j") * this->T1_("c,i");
    B_3("k,a,i,j") += -B_33("k,a,c,i") * this->T1_("c,j");
    B_3("k,a,i,j") +=   B_34("l,k,c,j") * this->T2_("c,a,l,i");
    B_3("k,a,i,j") += -B_34("l,k,c,i") * this->T2_("c,a,l,j");
    TAmanager.free("ov", std::move(B_31), true);
    TAmanager.free("oooo", std::move(B_32), true);
    TAmanager.free("ovvo", std::move(B_33), true);
    TAmanager.free("oovo", std::move(B_34), true);
  }
  
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_5() {
    B_5("k,l,i,j") = 0.25   *  conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T2_("c,d,i,j");
    B_5("k,l,i,j") +=  0.5  *   this->antiSymMoints["oooo"]("k,l,i,j");
    TA::TArray<MatsT> B_51 = TAManager::get().malloc<MatsT>("oovo");
    B_51("k,l,c,i") = 0.0 * B_51("k,l,c,i");
    B_51("k,l,c,i") += 0.25   *  conj(this->antiSymMoints["vvoo"]("d,c,k,l")) * this->T1_("d,i"); 
    B_51("k,l,c,i") += - 0.5  *   conj(this->antiSymMoints["vooo"]("c,i,k,l"));
    B_5("k,l,i,j") +=   B_51("k,l,c,i") * this->T1_("c,j");
    B_5("k,l,i,j") += -B_51("k,l,c,j") * this->T1_("c,i");
    TAManager::get().free("oovo", std::move(B_51), true);
  }
  template <typename MatsT, typename IntsT>
  void CCSD<MatsT,IntsT>::formB_6() {
    B_6("k,a,c,j") = - conj(this->antiSymMoints["vvvo"]("d,c,a,k")) * this->T1_("d,j");
    B_6("k,a,c,j") += - 0.5  *   conj(this->antiSymMoints["vvoo"]("c,d,k,l")) * this->T2_("d,a,l,j");
    B_6("k,a,c,j") +=  this->antiSymMoints["vovo"]("a,k,c,j");
  }

  template <typename MatsT, typename IntsT>
  CCSD<MatsT,IntsT>::~CCSD() {

    TAManager &TAmanager = TAManager::get();

    if (Fae_) TAmanager.free("vv", std::move(Fae_), true);
    if (Fmi_) TAmanager.free("oo", std::move(Fmi_), true);
    if (A_1) TAmanager.free("oo", std::move(A_1), true);
    if (A_2) TAmanager.free("vv", std::move(A_2), true);
    if (A_3) TAmanager.free("ov", std::move(A_3), true);
    if (A_4) TAmanager.free("oovo", std::move(A_4), true);
    if (B_1) TAmanager.free("oo", std::move(B_1), true);
    if (B_2) TAmanager.free("vv", std::move(B_2), true);
    if (B_3) TAmanager.free("ovoo", std::move(B_3), true);
    if (B_5) TAmanager.free("oooo", std::move(B_5), true);
    if (B_6) TAmanager.free("ovvo", std::move(B_6), true);

  }
  
};


