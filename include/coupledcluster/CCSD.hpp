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
  CCSD<MatsT>::CCSD(const SafeFile &savFile,
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
      Wmbej_(intermediates.W_mbej),
      eps(intermediates.eps){}

  template <typename MatsT>
  void CCSD<MatsT>::run() {
    runConventional();
  }

  template <typename MatsT>
  void CCSD<MatsT>::initIntermediates() {
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
  
  template <typename MatsT>
  void CCSD<MatsT>::initAmplitudes() {
    if(this->ccSettings_.restart){
      size_t size = this->T_.length();
      MatsT * t_amp = CQMemManager::get().malloc<MatsT>(size);
      TA::get_default_world().gop.fence();
      if (MPIRank() == 0) this->savFile_.readData("/CC/T_AMPLITUDE", t_amp);
      if (MPIRank() == 0) this->savFile_.readData("/CC/REFERENCE_ENERGY",   &this->intermediates_.E_ref);
      if (MPIRank() == 0) this->savFile_.readData("/CC/CORRELATION_ENERGY", &this->CorrE);
      MPIBCast(t_amp, size, 0, MPI_COMM_WORLD);
      MPIBCast(&this->intermediates_.E_ref, 1, 0, MPI_COMM_WORLD);
      MPIBCast(&this->CorrE               , 1, 0, MPI_COMM_WORLD);
      TA::get_default_world().gop.fence();
      this->T_.fromRaw(t_amp, false);
      TA::get_default_world().gop.fence();
      CQMemManager::get().free(t_amp);
    } else {
      this->T_.scale(0.0);
    }

  }

  template <typename MatsT>
  void CCBase<MatsT>::doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis ){

      // give solution vector to diis
      diis->WriteVector(this->T_);

      // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
      T_old.scale(-1.0);
      T_old.axpy(1.0, this->T_);

      //set error vector in DIIS
      diis->WriteErrorVector(T_old);

      // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
      diis->Extrapolate(this->T_);

  }


  template <typename MatsT>
  void CCBase<MatsT>::printBanner(double Eref) const {

    std::cout << BannerTop << std::endl;
    std::cout << "Coupled Cluster (CC) Settings:" << std::endl << std::endl;;

    if ( this->ccSettings_.cctype == CC_TYPE::CCSD ) {
      std::cout << std::setw(45) << std::left << "  Type:"
                << "Singles and Doubles" << std::endl;
    } else if ( this->ccSettings_.cctype == CC_TYPE::CCSDT ) {
      std::cout << std::setw(45) << std::left << "  Type:"
                << "Singles, Doubles, and Triples" << std::endl;
    }
    std::cout << std::setw(45) << std::left << "  Triples correction:";
    if ( this->ccSettings_.pertT3 ) {
      std::cout << "(T) with loop over "
                << (this->ccSettings_.loop_abc? "abc" : "ijk") << std::endl;
    }
    else if ( this->ccSettings_.crcc ) {
      std::cout << "CR-CC(2,3) with loop over"
                << (this->ccSettings_.loop_abc? "abc" : "ijk") << std::endl;
    }
    else {
      std::cout << "None" << std::endl;
    }

    std::cout << std::setw(45) << std::left << "  Occupied Orbitals:"
    << TAManager::get().getRange(this->oLabel_).extent() << std::endl;
    std::cout << std::setw(45) << std::left << "  Virtual Orbitals:"
    << TAManager::get().getRange(this->vLabel_).extent() << std::endl;
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

    std::cout << std::setw(45) << std::left << "  Denominator shift:"
    << this->ccSettings_.denomshift << std::endl;

    std::cout << std::endl;


    std::pair<double, char> mem_postfix = memSize(estimate_mem_peak());
    std::cout << std::setw(45) << std::left << "  Estimated TiledArray memory requirement: " << std::fixed << std::setprecision(1)
    << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    std::cout << BannerMid << std::endl << std::endl;

  }


  template <typename MatsT>
  size_t CCSD<MatsT>::estimate_mem_peak() const {
    // will need further checking as the TA objects are dynamically allocated and freed at runtime
    TAManager &TAmanager = TAManager::get();

    size_t nDIIS = this->ccSettings_.useDIIS ? this->ccSettings_.nDIIS : 0;
    size_t count = 0;
    //                                             1         2       3       4         5           6     7      8
    count += 8 * TAmanager.elem_per_TA("oo");   // muX_oo,   muY_oo, muZ_oo, coreH_oo, fock_oo,    Fmi_, TMPmj, moDen_oo
    count += 6 * TAmanager.elem_per_TA("ov");   // muX_ov,   muY_ov, muZ_ov, coreH_ov, fock_ov,    Fme_
    count += 7 * TAmanager.elem_per_TA("vo");   // muX_vo,   muY_vo, muZ_vo, coreH_vo, fock_vo,    Dai,  T1
    count += 7 * TAmanager.elem_per_TA("vv");   // muX_vv,   muY_vv, muZ_vv, coreH_vv, fock_vv,    Fae_, TMPbe
    count += 2 * TAmanager.elem_per_TA("oooo"); // ERI_oooo, Wmnij_
    count += 1 * TAmanager.elem_per_TA("ovoo"); // TMPmbij
    count += 1 * TAmanager.elem_per_TA("vooo"); // ERI_vooo
    count += 1 * TAmanager.elem_per_TA("ovvo"); // Wmbej_,   tmp
    count += 1 * TAmanager.elem_per_TA("vovo"); // ERI_vovo
    count += 7 * TAmanager.elem_per_TA("vvoo"); // ERI_vvoo, Dabij,  T2,     tau_,     tilde_tau_, tmp,  Pabij,
    count += 1 * TAmanager.elem_per_TA("vvvo"); // ERI_vvvo,
    count += 2 * TAmanager.elem_per_TA("vvvv"); // ERI_vvvv, W_abef

    if (nDIIS) {
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vo");   // T1 DIIS copy?
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vvoo"); // T2 DIIS copy?
    }

    // CCSD(T) additional objects
    if (this->ccSettings_.pertT3) {
      size_t count_pertt3 = 0;

      if (this->ccSettings_.loop_abc) {
        // containers are t3o3 if looping over abc
        count_pertt3 += 3 * TAmanager.elem_per_TA("tttooo"); // vt1, vt2, tmp
      } else {
        // containers are v3t3 if looping over ijk
        count_pertt3 += 3 * TAmanager.elem_per_TA("vvvttt"); // vt1, vt2, tmp
      }

      // account for MPI parallelization
      auto &global_world = TA::get_default_world();
      const auto size = global_world.size();
      if (this->ccSettings_.triplesMPI && size > 1) count_pertt3 = size * count_pertt3;
      count += count_pertt3;
    }

    // CR-CC(2,3) additional objects
    if (this->ccSettings_.crcc) {
      size_t count_crcc = 0;

      if (this->ccSettings_.loop_abc) {
        // containers are t3o3 if looping over abc
        count_crcc += 6 * TAmanager.elem_per_TA("tttooo"); // m3, l3a, l3b, l3c, l3d, tmp
      } else {
        // containers are v3t3 if looping over ijk
        count_crcc += 6 * TAmanager.elem_per_TA("vvvttt"); // m3, l3a, l3b, l3c, l3d, tmp
      }

      // account for MPI parallelization
      auto &global_world = TA::get_default_world();
      const auto size = global_world.size();
      if (this->ccSettings_.triplesMPI && size > 1) count_crcc = size * count_crcc;
      count += count_crcc;
    }

    return count * sizeof(MatsT);
  }



  template <typename MatsT>
  void CCSD<MatsT>::runConventional(){

    auto cc_start = tick();
    TAManager &TAmanager = TAManager::get();

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

      updateT1(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"));
      updateT2(T_old.get_tensor("OneBody"), T_old.get_tensor("TwoBody"));
#ifdef DEBUG_CCSD
      std::cout << "T1_:" << this->T1_ << std::endl;
      std::cout << "T2_:" << this->T2_ << std::endl;
#endif

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

        // run CCSD(T) *after* saving amplitudes and corrE
        if (this->ccSettings_.pertT3) {
          auto ccsdpt_start = tick();
          std::cout << std::setw(18) << std::left <<  "  -------------";
          std::cout << std::setw(34) << std::left << "-----------------";
          std::cout << std::setw(18) << std::right << "--------";
          std::cout << std::setw(18) << std::right << "----";
          std::cout << std::endl << std::endl;

//#define DEBUG_PERTT3
#ifdef DEBUG_PERTT3
          std::cout << "  Running (T) routine - building T3-sized objects" <<std::endl;

          runPertT3_straight(this->T1_,this->T2_);

          // Print results
          std::cout << std::endl << "  (T) Completed: CCSD    energy is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: (T) correction is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->PertT3Energy << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: CCSD(T) energy is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE + this->PertT3Energy << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: CCSD(T) total time "<< std::setw(10) << std::right
                    << std::setprecision(6) << tock(ccsdpt_start) << " s" << std::endl;

          this->PertT3Energy = 0.0;
          ccsdpt_start = tick();
#endif

          std::cout << "  Running (T) routine - batched algorithm" <<std::endl;
          if (this->ccSettings_.loop_abc) {
            runPertT3abc(this->T1_,this->T2_);
          } else {
            runPertT3ijk(this->T1_,this->T2_);
          }

          // Print results
          std::cout << std::endl << "  (T) Completed: CCSD    energy is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: (T) correction is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->PertT3Energy << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: CCSD(T) energy is "<< std::setw(18) << std::right
                    << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE + this->PertT3Energy << " Eh" << std::endl;
          std::cout << std::endl << "  (T) Completed: CCSD(T) total time "<< std::setw(10) << std::right
                    << std::setprecision(6) << tock(ccsdpt_start) << " s" << std::endl;

          if (this->savFile_.exists()) {
            this->savFile_.safeWriteData("/CC/CCSD(T)_CORRECTION",&this->PertT3Energy, {1});
          }

        }

        std::cout << bannerEnd << std::endl;

        printAnalysis();
        
        std::cout << BannerEnd << std::endl;

        break;
      }

      if(iter == this->ccSettings_.maxiter - 1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }
    }

  }

  template <typename MatsT>
  void CCBase<MatsT>::getCorrEnergy() {
    MatsT CorrEOneBody = this->fockMatrix_ta["ov"]("i,a").dot(this->T1_("a,i"));
    MatsT CorrETwoBodyT2 = conj(this->antiSymMoints["vvoo"]("c,d,k,l")).dot(this->T2_("c,d,k,l"));
    MatsT CorrETwoBodyT1 = conj(this->antiSymMoints["vvoo"]("c,d,k,l")).dot(this->T1_("c,k") * this->T1_("d,l"));
    TA::get_default_world().gop.fence();
    this->CorrE = CorrEOneBody + 0.25 * (CorrETwoBodyT2 + 2.0 * CorrETwoBodyT1);
  }

  template <typename MatsT>
  void CCSD<MatsT>::printAnalysis() {
    std::cout << bannerTop << std::endl;
    std::cout << "Coupled Cluster Wave Function Analysis:" << std::endl;
    std::cout << std::setw(22) << "   max(|t1|)  " << std::setprecision(4) << abs_max(this->T1_).get() << std::endl;
    std::cout << std::setw(22) << "   max(|t2|)  " << std::setprecision(4) << abs_max(this->T2_).get() << std::endl;
    double normT1Sq = squared_norm(this->T1_).get();
    double normT2Sq = squared_norm(this->T2_).get();
    std::cout << std::setw(22) << "   sum of t1 weights  " << std::setprecision(4) << normT1Sq << std::endl;
    std::cout << std::setw(22) << "   sum of t2 weights  " << std::setprecision(4) << normT2Sq << std::endl;
    std::cout << std::setw(22) << "   T1 diagnostic  " << std::setprecision(4) << sqrt(normT1Sq / this->intermediates_.nOcc) << std::endl;
    std::cout << std::setw(22) << "   T2 diagnostic  " << std::setprecision(4) << sqrt(normT2Sq / this->intermediates_.nOcc) << std::endl;
    std::cout << bannerEnd << std::endl;
  }

  /**
   * Build Eq. III(d.1) and Eq. III(d.2)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tau_and_tilde_tau() {
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
  void CCSD<MatsT>::build_tilde_Fae() {
    Fae_("a,e") = this->fockMatrix_ta["vv"]("a,e");
#ifdef DEBUG_CCSD
    std::cout << "Fae_1:" << Fae_ << std::endl;
#endif

    Fae_("a,e") -= 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("a,m");
#ifdef DEBUG_CCSD
    std::cout << "Fae_2:" << Fae_ << std::endl;
#endif
    Fae_("a,e") += this->T1_("f,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
#ifdef DEBUG_CCSD
    std::cout << "Fae_3:" << Fae_ << std::endl;
#endif

    Fae_("a,e") -= 0.5 * tilde_tau_("a,f,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
#ifdef DEBUG_CCSD
    std::cout << "Fae_4:" << Fae_ << std::endl;
#endif
  }

  /**
   * Build Eq. III(a.2)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tilde_Fmi() {
    Fmi_("m,i") = this->fockMatrix_ta["oo"]("m,i");

    Fmi_("m,i") += 0.5 * this->fockMatrix_ta["ov"]("m,e") * this->T1_("e,i");
    Fmi_("m,i") -= this->T1_("e,n") * conj(this->antiSymMoints["vooo"]("e,i,m,n"));

    Fmi_("m,i") += 0.5 * tilde_tau_("e,f,i,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n")); // TODO: better performance if we have this->antiSymMoints["oovv"]
  }

  /**
   * Build Eq. III(a.3)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tilde_Fme() {
    Fme_("m,e") = this->fockMatrix_ta["ov"]("m,e");
    Fme_("m,e") += this->T1_("f,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.4)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tilde_Wmnij() {
    Wmnij_("m,n,i,j") = -this->T1_("e,j") * conj(this->antiSymMoints["vooo"]("e,i,m,n")); // TODO: better performance if we have this->antiSymMoints["ooov"]
    Wmnij_("m,n,i,j") -= Wmnij_("m,n,j,i");

    Wmnij_("m,n,i,j") += this->antiSymMoints["oooo"]("m,n,i,j");

    Wmnij_("m,n,i,j") += 0.25 * tau_("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.5)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tilde_Wabef() {
    Wabef_("a,b,e,f") = -this->T1_("b,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m"));
    Wabef_("a,b,e,f") -= Wabef_("b,a,e,f");

    Wabef_("a,b,e,f") += this->antiSymMoints["vvvv"]("a,b,e,f");

    Wabef_("a,b,e,f") += 0.25 * tau_("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.6)
   */
  template <typename MatsT>
  void CCSD<MatsT>::build_tilde_Wmbej() {
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
  void CCSD<MatsT>::buildIntermediates() {
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
  template <typename MatsT>
  void CCSD<MatsT>::updateT1(const TArray T1_old, const TArray T2_old) {
    this->T1_("a,i") = this->fockMatrix_ta["vo"]("a,i");
    this->T1_("a,i") += Fae_("a,e") * T1_old("e,i");
    this->T1_("a,i") -= T1_old("a,m") * Fmi_("m,i");
    this->T1_("a,i") += Fme_("m,e") * T2_old("a,e,i,m");
    this->T1_("a,i") -= T1_old("e,m") * this->antiSymMoints["vovo"]("a,m,e,i");
    this->T1_("a,i") += 0.5 * T2_old("e,f,i,m") * conj(this->antiSymMoints["vvvo"]("e,f,a,m")); // TODO: better performance if we have this->antiSymMoints["vvov"]
    this->T1_("a,i") -= 0.5 * T2_old("a,e,m,n") * conj(this->antiSymMoints["vooo"]("e,i,n,m")); // TODO: better performance if we have this->antiSymMoints["ovoo"]

    this->T1_("a,i") = T1_old("a,i") + this->T1_("a,i") * this->Dai_("a,i");
  }

  /**
   * Build Eq. I(b)
   */
  template <typename MatsT>
  void CCSD<MatsT>::updateT2(const TArray T1_old, const TArray T2_old) {

    TAManager &TAmanager = TAManager::get();

    this->T2_("a,b,i,j") = this->antiSymMoints["vvoo"]("a,b,i,j");

    TArray TMPbe = TAmanager.malloc<MatsT>("vv");
    TMPbe("b,e") = Fae_("b,e");
    TMPbe("b,e") -= 0.5 * T1_old("b,m") * Fme_("m,e");
    TArray Pabij = TAmanager.malloc<MatsT>("vvoo");
    Pabij("a,b,i,j") = T2_old("a,e,i,j") * TMPbe("b,e");
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    this->T2_("a,b,i,j") -= Pabij("b,a,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_1:" << this->T2_ << std::endl;
#endif
    TAmanager.free("vv", std::move(TMPbe));

    TArray TMPmj = TAmanager.malloc<MatsT>("oo");
    TMPmj("m,j") = Fmi_("m,j");
    TMPmj("m,j") += 0.5 * T1_old("e,j") * Fme_("m,e");
    Pabij("a,b,i,j") = T2_old("a,b,i,m") * TMPmj("m,j");
    this->T2_("a,b,i,j") -= Pabij("a,b,i,j");
    this->T2_("a,b,i,j") += Pabij("a,b,j,i");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_2:" << this->T2_ << std::endl;
#endif
    TAmanager.free("oo", std::move(TMPmj));

    this->T2_("a,b,i,j") += 0.5 * tau_("a,b,m,n") * Wmnij_("m,n,i,j");
    this->T2_("a,b,i,j") += 0.5 * tau_("e,f,i,j") * Wabef_("a,b,e,f");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_3:" << this->T2_ << std::endl;
#endif

    TArray TMPmbij = TAmanager.malloc<MatsT>("ovoo");
    TMPmbij("m,b,i,j") = T1_old("e,i") * this->antiSymMoints["vovo"]("b,m,e,j");
    Pabij("a,b,i,j") = T1_old("a,m") * TMPmbij("m,b,i,j");
    Pabij("a,b,i,j") += T2_old("a,e,i,m") * Wmbej_("m,b,e,j");
    Pabij("a,b,i,j") -= Pabij("a,b,j,i");
    Pabij("a,b,i,j") -= Pabij("b,a,i,j");
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_4:" << this->T2_ << std::endl;
#endif
    TAmanager.free("ovoo", std::move(TMPmbij));

    Pabij("a,b,i,j") = T1_old("e,i") * this->antiSymMoints["vvvo"]("a,b,e,j");
    this->T2_("a,b,i,j") += Pabij("a,b,i,j");
    this->T2_("a,b,i,j") -= Pabij("a,b,j,i");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_5:" << this->T2_ << std::endl;
#endif

    Pabij("a,b,i,j") = -T1_old("a,m") * this->antiSymMoints["vooo"]("b,m,i,j"); // TODO: change sign
    this->T2_("a,b,i,j") -= Pabij("a,b,i,j");
    this->T2_("a,b,i,j") += Pabij("b,a,i,j");
#ifdef DEBUG_CCSD
    std::cout << "rhs_T2_6:" << this->T2_ << std::endl;
#endif

    TAmanager.free("vvoo", std::move(Pabij));
    this->T2_("a,b,i,j") = T2_old("a,b,i,j") + this->T2_("a,b,i,j") * this->Dabij_("a,b,i,j");
  }

  template <typename MatsT>
  void CCSD<MatsT>::runPertT3ijk(const TArray &T1_, const TArray &T2_) {

    TAManager &TAmanager = TAManager::get();

    std::cout << "  Initializing..." << std::endl;

    size_t o_tr = (this->intermediates_.nOcc + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
    size_t v_tr = (this->intermediates_.nVir + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);

    std::cout << "  Occ. TileRange: " << o_tr << std::endl;
    std::cout << "  Vir. TileRange: " << v_tr << std::endl;

    // lambda to compute VT2(C)
    auto compute_vt2 = [&](size_t &i, size_t &j, size_t &k, TArray &t3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      // grab relevant blocks
      TArray tmp = TAmanager.malloc<MatsT>("vvvttt");

      // t_{abc}^{ijk} contribution
      // t2_{eb}^{ij}
      block t21_lower = {0,0,i,j};
      block t21_upper = {v_tr,v_tr,i+1,j+1};
      // v_{ac}^{ek}
      block v1_lower = {0,0,0,k};
      block v1_upper = {v_tr,v_tr,v_tr,k+1};
      //t2_{ac}^{mk}
      block t22_lower = {0,0,0,k};
      block t22_upper = {v_tr,v_tr,o_tr,k+1};
      // v_{bm}^{ij}
      block v2_lower = {0,0,i,j};
      block v2_upper = {v_tr,o_tr,i+1,j+1};

      auto t2_ebij = T2_("e,b,i,j").block(t21_lower,t21_upper);
      auto v_acek  = this->antiSymMoints["vvvo"]("a,c,e,k").block(v1_lower,v1_upper);
      auto t2_acmk = T2_("a,c,m,k").block(t22_lower,t22_upper);
      auto v_bmij  = this->antiSymMoints["vooo"]("b,m,i,j").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k") = v_acek * t2_ebij + t2_acmk * v_bmij;

      // t_{abc}^{kji} contribution
      // t2_{eb}^{kj}
      t21_lower = {0,0,k,j};
      t21_upper = {v_tr,v_tr,k+1,j+1};
      // v_{ac}^{ei}
      v1_lower = {0,0,0,i};
      v1_upper = {v_tr,v_tr,v_tr,i+1};
      //t2_{ac}^{mi}
      t22_lower = {0,0,0,i};
      t22_upper = {v_tr,v_tr,o_tr,i+1};
      // v_{bm}^{kj}
      v2_lower = {0,0,k,j};
      v2_upper = {v_tr,o_tr,k+1,j+1};

      auto t2_ebkj = T2_("e,b,k,j").block(t21_lower,t21_upper);
      auto v_acei  = this->antiSymMoints["vvvo"]("a,c,e,i").block(v1_lower,v1_upper);
      auto t2_acmi = T2_("a,c,m,i").block(t22_lower,t22_upper);
      auto v_bmkj  = this->antiSymMoints["vooo"]("b,m,k,j").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k") -= v_acei * t2_ebkj + t2_acmi * v_bmkj;

      // t_{abc}^{ikj} contribution
      // t2_{eb}^{ik}
      t21_lower = {0,0,i,k};
      t21_upper = {v_tr,v_tr,i+1,k+1};
      // v_{ac}^{ej}
      v1_lower = {0,0,0,j};
      v1_upper = {v_tr,v_tr,v_tr,j+1};
      //t2_{ac}^{mj}
      t22_lower = {0,0,0,j};
      t22_upper = {v_tr,v_tr,o_tr,j+1};
      // v_{bm}^{ik}
      v2_lower = {0,0,i,k};
      v2_upper = {v_tr,o_tr,i+1,k+1};

      auto t2_ebik = T2_("e,b,i,k").block(t21_lower,t21_upper);
      auto v_acej  = this->antiSymMoints["vvvo"]("a,c,e,j").block(v1_lower,v1_upper);
      auto t2_acmj = T2_("a,c,m,j").block(t22_lower,t22_upper);
      auto v_bmik  = this->antiSymMoints["vooo"]("b,m,i,k").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k") -= v_acej * t2_ebik + t2_acmj * v_bmik;

      // apply A(ac/b)
      t3 = tmp.clone(); // deep copy for safety
      t3("a,b,c,i,j,k") -= tmp("b,a,c,i,j,k");
      t3("a,b,c,i,j,k") -= tmp("a,c,b,i,j,k");

      TAmanager.free("vvvttt", std::move(tmp));
    };

    // lambda to compute VT1(DC)
    auto compute_vt1 = [&](size_t &i, size_t &j, size_t &k, TArray &t3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      TArray tmp    = TAmanager.malloc<MatsT>("vvvttt");

      // t_{abc}^{ijk} contribution
      // t1_{c}_{k}
      block t1_lower = {0,k};
      block t1_upper = {v_tr,k+1};
      // v_{ab}^{ij}
      block v_lower = {0,0,i,j};
      block v_upper = {v_tr,v_tr,i+1,j+1};

      auto v_abij = this->antiSymMoints["vvoo"]("a,b,i,j").block(v_lower,v_upper);
      auto t1_ck  = T1_("c,k").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") = v_abij * t1_ck;

      // t_{abc}^{kji} contribution
      // t1_{c}_{i}
      t1_lower = {0,i};
      t1_upper = {v_tr,i+1};
      // v_{ab}^{kj}
      v_lower = {0,0,k,j};
      v_upper = {v_tr,v_tr,k+1,j+1};

      auto v_abkj = this->antiSymMoints["vvoo"]("a,b,k,j").block(v_lower,v_upper);
      auto t1_ci  = T1_("c,i").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") -= v_abkj * t1_ci;

      // t_{abc}^{ikj} contribution
      // t1_{c}_{j}
      t1_lower = {0,j};
      t1_upper = {v_tr,j+1};
      // v_{ab}^{ik}
      v_lower = {0,0,i,k};
      v_upper = {v_tr,v_tr,i+1,k+1};

      auto v_abik = this->antiSymMoints["vvoo"]("a,b,i,k").block(v_lower,v_upper);
      auto t1_cj  = T1_("c,j").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") -= v_abik * t1_cj;

      // apply A(ab/c)
      t3  = tmp.clone(); // deep copy for safety
      t3("a,b,c,i,j,k") -= tmp("c,b,a,i,j,k");
      t3("a,b,c,i,j,k") -= tmp("a,c,b,i,j,k");

      // Return (VT1)*
      t3("a,b,c,i,j,k") = conj(t3("a,b,c,i,j,k"));

      TAmanager.free("vvvttt", std::move(tmp));
    };

    // MPI parallelization from MPQC
    // split global_world
    auto &global_world = TA::get_default_world();

    const auto rank = global_world.rank();
    const auto size = global_world.size();

    madness::World *tmp_ptr;
    std::shared_ptr<madness::World> world_ptr;

    if (this->ccSettings_.triplesMPI && size > 1) {
      SafeMPI::Group group = global_world.mpi.comm().Get_group().Incl(1, &rank);
      SafeMPI::Intracomm comm = global_world.mpi.comm().Create(group);
      world_ptr = std::make_shared<madness::World>(comm);
      tmp_ptr = world_ptr.get();
    } else {
      tmp_ptr = &global_world;
    }
    auto &this_world = *tmp_ptr;
    global_world.gop.fence();

    // global iteration
    size_t global_iter = 0;

    // this will set TA::get_default_world() from this point on to &this_world
    TA::set_default_world(this_world);

    // loop over i,j,k
    std::cout << "  Starting loop over i>=j>=k tiles" << std::endl;
    std::cout << "  Running " << o_tr*(o_tr+1)*(o_tr+2)/6 << " iterations on "
              << size << " MPI ranks " << std::endl;

    // parse loop cut
    size_t end_loop = o_tr;
    if (this->ccSettings_.triples_end) end_loop = this->ccSettings_.triples_end;

    // loop from i=triples_begin and j=k=0 to i=j=k=end_loop-1
    for (size_t i = this->ccSettings_.triples_begin; i < end_loop; ++i)
      for (size_t j = 0; j <= i ; ++j)
        for (size_t k = 0; k <= j; ++k) {
          auto ijk_start = tick();

          global_iter++;

          // Distribute ijk term for MPI-parallel version
          // round-robin distribute the loop
          if (this->ccSettings_.triplesMPI && global_iter % size != rank) continue;

          TArray vt1_ijk = TAmanager.malloc<MatsT>("vvvttt");
          TArray vt2_ijk = TAmanager.malloc<MatsT>("vvvttt");

          // <ijkabc| (VT2)_C |0>
          compute_vt2(i,j,k,vt2_ijk);

          // <0| (VT2)*_C + (VT1)*_DC | ijkabc)
          compute_vt1(i,j,k,vt1_ijk);
          vt1_ijk("a,b,c,i,j,k") += conj(vt2_ijk("a,b,c,i,j,k")).set_world(this_world);

          // divide by MP denominator
          TA::foreach_inplace(vt1_ijk, [this, &i, &j, &k](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){
                        // offset i,j,k of blocked TiledArray object, otherwise wrong eps[x] index
                        size_t i_off = i*this->ccSettings_.blksize;
                        size_t j_off = j*this->ccSettings_.blksize;
                        size_t k_off = k*this->ccSettings_.blksize;

                        size_t a = x[0]+this->intermediates_.nOcc;
                        size_t b = x[1]+this->intermediates_.nOcc;
                        size_t c = x[2]+this->intermediates_.nOcc;
                        size_t ii = x[3] + i_off, jj = x[4] + j_off, kk = x[5] + k_off;

                        tile[x] /= this->intermediates_.eps[ii] + this->intermediates_.eps[jj] + this->intermediates_.eps[kk]
                                 - this->intermediates_.eps[a]  - this->intermediates_.eps[b]  - this->intermediates_.eps[c];
                      }
          });

          MatsT tmp_en = vt1_ijk("a,b,c,i,j,k").dot(vt2_ijk("a,b,c,i,j,k"));
          if ((i==j) && (j==k)) {
            // ijk -> 3!, abc -> 3!
            tmp_en *= 1.0/36.0;
          }
          else if ((i==j) != (j==k)) {
            // ijk -> 2!, abc -> 3!
            tmp_en *= 1.0/12.0;
          }
          else {
            // ijk -> 1!, abc -> 3!
            tmp_en *= 1.0/6.0;
          }
          this->PertT3Energy += tmp_en;

          TAmanager.free("vvvttt", std::move(vt1_ijk));
          TAmanager.free("vvvttt", std::move(vt2_ijk));

          std::cout << "  i: " << i << " j: " << j << " k: " << k << " done in " << tock(ijk_start) << "s"
                    << " from global_iter " << global_iter << " in rank " << rank << std::endl;
        }

    this_world.gop.fence();
    global_world.gop.fence();

    TA::set_default_world(global_world);

    // sum over contribution if MPI parallel is run
    if (this->ccSettings_.triplesMPI && size >1) global_world.gop.sum(this->PertT3Energy);

  }

  template <typename MatsT>
  void CCSD<MatsT>::runPertT3abc(const TArray &T1_, const TArray &T2_) {

    TAManager &TAmanager = TAManager::get();

    std::cout << "  Initializing..." << std::endl;

    size_t o_tr = (this->intermediates_.nOcc + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
    size_t v_tr = (this->intermediates_.nVir + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);

    std::cout << "  Occ. TileRange: " << o_tr << std::endl;
    std::cout << "  Vir. TileRange: " << v_tr << std::endl;

    // lambda to compute VT2(C)
    auto compute_vt2 = [&](size_t &a, size_t &b, size_t &c, TArray &t3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      // grab relevant blocks
      TArray tmp     = TAmanager.malloc<MatsT>("tttooo");

      // t_{abc}^{ijk} contribution
      // t2_{eb}^{ij}
      block t21_lower = {0,b,0,0};
      block t21_upper = {v_tr,b+1,o_tr,o_tr};
      // v_{ac}^{ek}
      block v1_lower = {a,c,0,0};
      block v1_upper = {a+1,c+1,v_tr,o_tr};
      //t2_{ac}^{mk}
      block t22_lower = {a,c,0,0};
      block t22_upper = {a+1,c+1,o_tr,o_tr};
      // v_{bm}^{ij}
      block v2_lower = {b,0,0,0};
      block v2_upper = {b+1,o_tr,o_tr,o_tr};

      auto t2_ebij = T2_("e,b,i,j").block(t21_lower,t21_upper);
      auto v_acek  = this->antiSymMoints["vvvo"]("a,c,e,k").block(v1_lower,v1_upper);
      auto t2_acmk = T2_("a,c,m,k").block(t22_lower,t22_upper);
      auto v_bmij  = this->antiSymMoints["vooo"]("b,m,i,j").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k")  = v_acek * t2_ebij + t2_acmk * v_bmij;

      // t_{bac}^{ijk} contribution
      // t2_{ea}^{ij}
      t21_lower = {0,a,0,0};
      t21_upper = {v_tr,a+1,o_tr,o_tr};
      // v_{bc}^{ek}
      v1_lower = {b,c,0,0};
      v1_upper = {b+1,c+1,v_tr,o_tr};
      //t2_{bc}^{mk}
      t22_lower = {b,c,0,0};
      t22_upper = {b+1,c+1,o_tr,o_tr};
      // v_{am}^{ij}
      v2_lower = {a,0,0,0};
      v2_upper = {a+1,o_tr,o_tr,o_tr};

      // grab relevant blocks
      auto t2_eaij = T2_("e,a,i,j").block(t21_lower,t21_upper);
      auto v_bcek  = this->antiSymMoints["vvvo"]("b,c,e,k").block(v1_lower,v1_upper);
      auto t2_bcmk = T2_("b,c,m,k").block(t22_lower,t22_upper);
      auto v_amij  = this->antiSymMoints["vooo"]("a,m,i,j").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k") -= v_bcek * t2_eaij + t2_bcmk * v_amij;

      // t_{acb}^{ijk} contribution
      // t2_{ec}^{ij}
      t21_lower = {0,c,0,0};
      t21_upper = {v_tr,c+1,o_tr,o_tr};
      // v_{ab}^{ek}
      v1_lower = {a,b,0,0};
      v1_upper = {a+1,b+1,v_tr,o_tr};
      //t2_{ab}^{mk}
      t22_lower = {a,b,0,0};
      t22_upper = {a+1,b+1,o_tr,o_tr};
      // v_{mc}^{ij}
      v2_lower = {c,0,0,0};
      v2_upper = {c+1,o_tr,o_tr,o_tr};

      // grab relevant blocks
      auto t2_ecij = T2_("e,c,i,j").block(t21_lower,t21_upper);
      auto v_abek  = this->antiSymMoints["vvvo"]("a,b,e,k").block(v1_lower,v1_upper);
      auto t2_abmk = T2_("a,b,m,k").block(t22_lower,t22_upper);
      auto v_cmij  = this->antiSymMoints["vooo"]("c,m,i,j").block(v2_lower,v2_upper);

      tmp("a,b,c,i,j,k") -= v_abek * t2_ecij + t2_abmk * v_cmij;

      // apply A(ij/k)
      t3 = tmp.clone(); // deep copy for safety
      t3("a,b,c,i,j,k") -= tmp("a,b,c,k,j,i");
      t3("a,b,c,i,j,k") -= tmp("a,b,c,i,k,j");

      TAmanager.free("tttooo", std::move(tmp));
    };

    // lambda to compute VT1(DC)
    auto compute_vt1 = [&](size_t &a, size_t &b, size_t &c, TArray &t3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      TArray tmp    = TAmanager.malloc<MatsT>("tttooo");

      // t_{abc}^{ijk} contribution
      // t1_{c}_{k}
      block t1_lower = {c,0};
      block t1_upper = {c+1,o_tr};
      // v_{ab}^{ij}
      block v_lower = {a,b,0,0};
      block v_upper = {a+1,b+1,o_tr,o_tr};

      auto v_abij = this->antiSymMoints["vvoo"]("a,b,i,j").block(v_lower,v_upper);
      auto t1_ck  = T1_("c,k").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") = v_abij * t1_ck;

      // t_{cba}^{ijk} contribution
      // t1_{a}_{k}
      t1_lower = {a,0};
      t1_upper = {a+1,o_tr};
      // v_{cb}^{ij}
      v_lower = {c,b,0,0};
      v_upper = {c+1,b+1,o_tr,o_tr};

      auto v_cbij = this->antiSymMoints["vvoo"]("c,b,i,j").block(v_lower,v_upper);
      auto t1_ak  = T1_("a,k").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") -= v_cbij * t1_ak;

      // t_{acb}^{ijk} contribution
      // t1_{b}_{k}
      t1_lower = {b,0};
      t1_upper = {b+1,o_tr};
      // v_{ac}^{ij}
      v_lower = {a,c,0,0};
      v_upper = {a+1,c+1,o_tr,o_tr};

      auto v_acij = this->antiSymMoints["vvoo"]("a,c,i,j").block(v_lower,v_upper);
      auto t1_bk  = T1_("b,k").block(t1_lower,t1_upper);

      tmp("a,b,c,i,j,k") -= v_acij * t1_bk;

      // apply A(ij/k)
      t3  = tmp.clone(); // deep copy for safety
      t3("a,b,c,i,j,k") -= tmp("a,b,c,k,j,i");
      t3("a,b,c,i,j,k") -= tmp("a,b,c,i,k,j");

      // Return (VT1)*
      t3("a,b,c,i,j,k") = conj(t3("a,b,c,i,j,k"));

      TAmanager.free("tttooo", std::move(tmp));
    };

    // MPI parallelization from MPQC
    // split global_world
    auto &global_world = TA::get_default_world();

    const auto rank = global_world.rank();
    const auto size = global_world.size();

    madness::World *tmp_ptr;
    std::shared_ptr<madness::World> world_ptr;

    if (this->ccSettings_.triplesMPI && size > 1) {
      SafeMPI::Group group = global_world.mpi.comm().Get_group().Incl(1, &rank);
      SafeMPI::Intracomm comm = global_world.mpi.comm().Create(group);
      world_ptr = std::make_shared<madness::World>(comm);
      tmp_ptr = world_ptr.get();
    } else {
      tmp_ptr = &global_world;
    }
    auto &this_world = *tmp_ptr;
    global_world.gop.fence();

    // global iteration
    size_t global_iter = 0;

    // this will set TA::get_default_world() from this point on to &this_world
    TA::set_default_world(this_world);

    // loop over a,b,c
    std::cout << "  Starting loop over a>=b>=c tiles" << std::endl;
    std::cout << "  Running " << v_tr*(v_tr+1)*(v_tr+2)/6 << " iterations on "
              << size << " MPI ranks " << std::endl;

    // parse loop cut
    size_t end_loop = v_tr;
    if (this->ccSettings_.triples_end) end_loop = this->ccSettings_.triples_end;

    // loop from a=triples_begin and b=c=0 to a=b=c=end_loop-1
    for (size_t a = this->ccSettings_.triples_begin; a < end_loop; ++a)
      for (size_t b = 0; b <= a ; ++b)
        for (size_t c = 0; c <= b; ++c) {
          auto abc_start = tick();

          global_iter++;

          // Distribute abc term for MPI-parallel version
          // round-robin distribute the loop
          if (this->ccSettings_.triplesMPI && global_iter % size != rank) continue;

          TArray vt1_abc = TAmanager.malloc<MatsT>("tttooo");
          TArray vt2_abc = TAmanager.malloc<MatsT>("tttooo");

          // <ijkabc| (VT2)_C |0>
          compute_vt2(a,b,c,vt2_abc);

          // <0| (VT2)*_C + (VT1)*_DC | ijkabc)
          compute_vt1(a,b,c,vt1_abc);
          vt1_abc("a,b,c,i,j,k") += conj(vt2_abc("a,b,c,i,j,k")).set_world(this_world);

          // divide by MP denominator
          TA::foreach_inplace(vt1_abc, [this, &a, &b, &c](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){
                        // offset a,b,c of blocked TiledArray object, otherwise wrong eps[x] index
                        size_t a_off = a*this->ccSettings_.blksize;
                        size_t b_off = b*this->ccSettings_.blksize;
                        size_t c_off = c*this->ccSettings_.blksize;

                        size_t aa = x[0]+this->intermediates_.nOcc+a_off;
                        size_t bb = x[1]+this->intermediates_.nOcc+b_off;
                        size_t cc = x[2]+this->intermediates_.nOcc+c_off;
                        size_t i = x[3], j = x[4], k = x[5];

                        tile[x] /= this->intermediates_.eps[i]  + this->intermediates_.eps[j]  + this->intermediates_.eps[k]
                                 - this->intermediates_.eps[aa] - this->intermediates_.eps[bb] - this->intermediates_.eps[cc];
                      }
          });

          MatsT tmp_en = vt1_abc("a,b,c,i,j,k").dot(vt2_abc("a,b,c,i,j,k"));
          if ((a==b) && (b==c)) {
            // abc -> 3!, ijk -> 3!
            tmp_en *= 1.0/36.0;
          }
          else if ((a==b) != (b==c)) {
            // abc -> 2!, ijk -> 3!
            tmp_en *= 1.0/12.0;
          }
          else {
            // abc -> 1!, ijk -> 3!
            tmp_en *= 1.0/6.0;
          }
          this->PertT3Energy += tmp_en;

          TAmanager.free("tttooo", std::move(vt1_abc));
          TAmanager.free("tttooo", std::move(vt2_abc));

          std::cout << "  a: " << a << " b: " << b << " c: " << c << " done in " << tock(abc_start) << "s"
                    << " from global_iter " << global_iter << " in rank " << rank << std::endl;
        }

    this_world.gop.fence();
    global_world.gop.fence();

    TA::set_default_world(global_world);

    // sum over contribution if MPI parallel is run
    if (this->ccSettings_.triplesMPI && size >1) global_world.gop.sum(this->PertT3Energy);

  }


  template <typename MatsT>
  CCSD<MatsT>::~CCSD() {

    TAManager &TAmanager = TAManager::get();

    cleanMemory();

  }

  template <typename MatsT>
  void CCSD<MatsT>::cleanMemory(){
    TAManager &TAmanager = TAManager::get();
//    if(Fae_) TAmanager.free("vv", std::move(Fae_), true);
//    if(Fmi_) TAmanager.free("oo", std::move(Fmi_), true);
    if(tilde_tau_) TAmanager.free("vvoo", std::move(tilde_tau_), true);
  }
  
};


