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
  void EOMCCSD<MatsT>::initializeLambda() {

    if (this->Lg_ == nullptr) this->initializeGroundStateLambda();
    this->Lg_->zeroBody() = 1.0;
    this->Lg_->get_tensor("OneBody")("a,i") = conj(T1_("a,i"));
    this->Lg_->get_tensor("TwoBody")("a,b,i,j") = conj(T2_("a,b,i,j"));

    TAManager &TAmanager = TAManager::get();
    if (not G_ae.is_initialized()){
      G_ae = TAmanager.malloc<MatsT>("vv");
    }
    if (not G_mi.is_initialized()){
      G_mi = TAmanager.malloc<MatsT>("oo");
    }

  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::updateG_ae(const TArray &L2, TArray &G_ae) const {
    G_ae("a,e") = - 0.5 * T2_("e,f,m,n") * L2("a,f,m,n");
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::updateG_mi(const TArray &L2, TArray &G_mi) const {
    G_mi("m,i") = 0.5 * T2_("e,f,m,n") * L2("e,f,i,n");
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formL1_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                                          TArray &tildeL1) const {

    tildeL1("a,i") = F_ae("e,a") * L1("e,i");
    tildeL1("a,i") += - F_mi("i,m") * L1("a,m");
    tildeL1("a,i") += L1("e,m") * W_mbej("i,e,a,m");
    tildeL1("a,i") += 0.5 * L2("e,f,i,m") * W_abei("e,f,a,m");
    tildeL1("a,i") += - 0.5 * L2("a,e,m,n") * W_mbij("i,e,m,n");
    tildeL1("a,i") += - G_ae("e,f") * W_amef("e,i,f,a");
    tildeL1("a,i") += - G_mi("m,n") * W_mnie("m,i,n,a");

  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formL2_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                                          TArray &tildeL2) const {

    tildeL2("a,b,i,j") = L2("a,e,i,j") * F_ae("e,b");
    tildeL2("a,b,i,j") += - L2("b,e,i,j") * F_ae("e,a");

    tildeL2("a,b,i,j") += - L2("a,b,i,m") * F_mi("j,m");
    tildeL2("a,b,i,j") += L2("a,b,j,m") * F_mi("i,m");

    tildeL2("a,b,i,j") += 0.5 * L2("a,b,m,n") * W_mnij("i,j,m,n");

    tildeL2("a,b,i,j") += 0.5 * L2("e,f,i,j") * W_abef("e,f,a,b");

    tildeL2("a,b,i,j") += L1("e,i") * W_amef("e,j,a,b");
    tildeL2("a,b,i,j") += - L1("e,j") * W_amef("e,i,a,b");

    tildeL2("a,b,i,j") += - L1("a,m") * W_mnie("i,j,m,b");
    tildeL2("a,b,i,j") += L1("b,m") * W_mnie("i,j,m,a");

    tildeL2("a,b,i,j") += L2("a,e,i,m") * W_mbej("j,e,b,m");
    tildeL2("a,b,i,j") += - L2("a,e,j,m") * W_mbej("i,e,b,m");
    tildeL2("a,b,i,j") += - L2("b,e,i,m") * W_mbej("j,e,a,m");
    tildeL2("a,b,i,j") += L2("b,e,j,m") * W_mbej("i,e,a,m");

    tildeL2("a,b,i,j") += L1("a,i") * F_me("j,b");
    tildeL2("a,b,i,j") += - L1("a,j") * F_me("i,b");
    tildeL2("a,b,i,j") += - L1("b,i") * F_me("j,a");
    tildeL2("a,b,i,j") += L1("b,j") * F_me("i,a");

    tildeL2("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,e,i,j")) * G_ae("b,e");
    tildeL2("a,b,i,j") += - conj(this->antiSymMoints["vvoo"]("b,e,i,j")) * G_ae("a,e");

    tildeL2("a,b,i,j") += - conj(this->antiSymMoints["vvoo"]("a,b,i,m")) * G_mi("m,j");
    tildeL2("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,b,j,m")) * G_mi("m,i");

  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::runLambda() {

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

      updateG_ae(L_old.get_tensor("TwoBody"), G_ae);
      updateG_mi(L_old.get_tensor("TwoBody"), G_mi);

      formL1_tilde(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody"), G_ae, G_mi, L1_);
      L1_("a,i") += F_me("i,a");
      L1_("a,i") = L_old.get_tensor("OneBody")("a,i") + L1_("a,i") * D_ai("a,i");

      formL2_tilde(L_old.get_tensor("OneBody"), L_old.get_tensor("TwoBody"), G_ae, G_mi, L2_);
      L2_("a,b,i,j") += conj(this->antiSymMoints["vvoo"]("a,b,i,j"));
      L2_("a,b,i,j") = L_old.get_tensor("TwoBody")("a,b,i,j") + L2_("a,b,i,j") * D_abij("a,b,i,j");
      
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
      MatsT PEOneBody = this->fockMatrix_ta["vo"]("a,i").dot(L1_("a,i")).get();
      MatsT PETwoBodyL2 = this->antiSymMoints["vvoo"]("a,b,i,j").dot(L2_("a,b,i,j")).get();
      MatsT PETwoBodyL1 = this->antiSymMoints["vvoo"]("c,d,k,l").dot(L1_("c,k") * L1_("d,l")).get();
      TA::get_default_world().gop.fence();    
      pseudoEnergy = PEOneBody + 0.25 * (PETwoBodyL2 + 2.0 * PETwoBodyL1);

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

  template <typename MatsT>
  void EOMCCSD<MatsT>::runCR(const MatsT &CorrE){
    auto crcc_start = tick();
    std::cout << "  Entering CR-CC(2,3) routine..." << std::endl;

    // Initialize L1, L2
    //MBExpansionSet<MatsT> &Lg_eom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->Lg_);
    TArray &L1_ = this->Lg_->get_tensor("OneBody");
    TArray &L2_ = this->Lg_->get_tensor("TwoBody");

    TAManager &TAmanager = TAManager::get();

    builddenom(T2_);

    if (this->ccSettings_.loop_abc) {
      runCRbatchabc(T1_, T2_, L1_, L2_);
    } else {
      runCRbatchijk(T1_, T2_, L1_, L2_);
    }

    std::cout << "\n  (2,3) Completed: CCSD        energy is "<< std::setw(18) << std::right
              << std::setprecision(12) << this->intermediates_.E_ref + CorrE << " Eh" << std::endl;
    std::cout << "\n  (2,3) Completed: MM(2,3)A correction is "<< std::setw(18) << std::right
              << std::setprecision(12) << CRCC23Energy_A << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: MM(2,3)B correction is "<< std::setw(18) << std::right
              << std::setprecision(12) << CRCC23Energy_B << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: MM(2,3)C correction is "<< std::setw(18) << std::right
              << std::setprecision(12) << CRCC23Energy_C << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: MM(2,3)D correction is "<< std::setw(18) << std::right
              << std::setprecision(12) << CRCC23Energy_D << " Eh" << std::endl;
    std::cout << "\n  (2,3) Completed: CR-CC(2,3)A  energy is "<< std::setw(18) << std::right
              << std::setprecision(12) << this->intermediates_.E_ref + CorrE + CRCC23Energy_A << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: CR-CC(2,3)B  energy is "<< std::setw(18) << std::right
              << std::setprecision(12) << this->intermediates_.E_ref + CorrE + CRCC23Energy_B << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: CR-CC(2,3)C  energy is "<< std::setw(18) << std::right
              << std::setprecision(12) << this->intermediates_.E_ref + CorrE + CRCC23Energy_C << " Eh" << std::endl;
    std::cout << "  (2,3) Completed: CR-CC(2,3)D  energy is "<< std::setw(18) << std::right
              << std::setprecision(12) << this->intermediates_.E_ref + CorrE + CRCC23Energy_D << " Eh" << std::endl;
    std::cout << "\n  (2,3) Completed: CR-CC(2,3) total time "<< std::setw(10) << std::right
              << std::setprecision(6) << tock(crcc_start) << " s" << std::endl;

    if (this->savFile_.exists()) {
      std::vector<MatsT> CREnergies = {CRCC23Energy_A,
                                       CRCC23Energy_B,
                                       CRCC23Energy_C,
                                       CRCC23Energy_D};
      this->savFile_.safeWriteData("/CC/CR-CC(2,3)_CORRECTION",CREnergies.data(), {4});
    }

  }

  // Consider making denominator as vectors to minimize T3-sized object usage
  template <typename MatsT>
  void EOMCCSD<MatsT>::builddenom(const TArray &T2){
    auto timer = tick();
    std::cout << "  Building denominator ..." << std::endl;

    TAManager &TAmanager = TAManager::get();
    size_t nO = TAmanager.getRange(oLabel_).extent();
    size_t nV = TAmanager.getRange(vLabel_).extent();

    // Initialize arrays
    size_t nOrb_ = nO + nV;
    size_t nO2 = nO * nO;
    size_t nOV = nO * nV;
    size_t nV2 = nV * nV;

    // Just making sure
    hbar1.clear();
    hbar2.clear();
    hbar3.clear();

    // Resize and fill with zeroes
    hbar1.resize(nOrb_, 0.0);
    hbar2.resize(nO2+nOV+nV2, 0.0);
    hbar3.resize(nO*nV*nOrb_, 0.0);

    // Build hbar1 part of denom
    // Populate hbar1[o]
    TA::foreach_inplace( F_mi, [&](TA::Tensor<MatsT> &tile){
      const auto &lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto &upbound = tile.range().upbound();
        std::size_t x[2] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          size_t i = x[0];
          hbar1[i] = tile[x];
        }
      }
    });

    // Populate hbar1[v]
    TA::foreach_inplace( F_ae, [&](TA::Tensor<MatsT> &tile){
      const auto &lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto &upbound = tile.range().upbound();
        std::size_t x[2] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          size_t a = nO + x[0];
          hbar1[a] = tile[x];
        }
      }
    });
    TA::get_default_world().gop.fence();
    MatsT *hbar1_copy = CQMemManager::get().malloc<MatsT>(nOrb_);
    std::copy_n(hbar1.data(), nOrb_, hbar1_copy);
    std::fill_n(hbar1.data(), nOrb_, MatsT(0.0));
    MPIAllReduce(hbar1_copy, nOrb_, hbar1.data(), MPI_COMM_WORLD);
    CQMemManager::get().free(hbar1_copy);
    TA::get_default_world().gop.fence();


    // Build hbar2 part of denom
    // Populate hbar2[oo]
    TA::foreach_inplace( W_mnij, [&](TA::Tensor<MatsT> &tile){
      const auto &lobound = tile.range().lobound();
      if (lobound[0] == lobound[2] && lobound[1] == lobound[3]) {
        const auto &upbound = tile.range().upbound();
        std::size_t x[4] = {0, 0, 0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
            x[2] = x[0]; //oooo, m->i
            x[3] = x[1]; //oooo, n->j
            size_t ij = x[0] * nO + x[1];
            hbar2[ij] = tile[x];
        }
      }
    });

    // Populate hbar2[ov]
    TA::foreach_inplace( W_mbej, [&](TA::Tensor<MatsT> &tile){
      const auto &lobound = tile.range().lobound();
      if (lobound[0] == lobound[3] && lobound[1] == lobound[2]) {
        const auto &upbound = tile.range().upbound();
        std::size_t x[4] = {0, 0, 0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
            x[3] = x[0]; //ovvo, m->j
            x[2] = x[1]; //ovvo, b->e
            size_t ia = nO2 + x[0] * nV + x[1];
            hbar2[ia] = tile[x];
        }
      }
    });

    // Populate hbar2[vv]
    TA::foreach_inplace( W_abef, [&](TA::Tensor<MatsT> &tile){
      const auto &lobound = tile.range().lobound();
      if (lobound[0] == lobound[2] && lobound[1] == lobound[3]) {
        const auto &upbound = tile.range().upbound();
        std::size_t x[4] = {0, 0, 0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
            x[2] = x[0]; //vvvv, a->e
            x[3] = x[1]; //vvvv, b->f
            size_t ab = (nO2 + nOV) + x[0] * nV + x[1];
            hbar2[ab] = tile[x];
          }
      }
    });
    TA::get_default_world().gop.fence();

    MatsT *hbar2_copy = CQMemManager::get().malloc<MatsT>(nO2+nOV+nV2);
    std::copy_n(hbar2.data(), nO2+nOV+nV2, hbar2_copy);
    std::fill_n(hbar2.data(), nO2+nOV+nV2, MatsT(0.0));
    MPIAllReduce(hbar2_copy, nO2+nOV+nV2, hbar2.data(), MPI_COMM_WORLD);
    CQMemManager::get().free(hbar2_copy);

    TA::get_default_world().gop.fence();


    // Build hbar3 part of denom
    {
      // Populate hbar3[ovo]
      TArray tmp_VT2 = TAmanager.malloc<MatsT>("ovo");
      TArray tmp_ones = TAmanager.malloc_fresh<MatsT>("v");
      tmp_ones.fill(1.0);
      // Consider using TA::expressions::einsum, but have to be careful with indices
      tmp_VT2("i,b,j") = conj(this->antiSymMoints["vvoo"]("e,b,i,j")) * T2("e,b,i,j") * tmp_ones("e");
      TA::get_default_world().gop.fence();

      TA::foreach_inplace(tmp_VT2, [&](TA::Tensor<MatsT> &tile) {
        const auto &lobound = tile.range().lobound();
        const auto &upbound = tile.range().upbound();
        std::size_t x[3] = {0, 0, 0};
        for (x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for (x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
            for (x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]) {
              size_t iaj = (x[0] * nV + x[1]) * nO + x[2];
              hbar3[iaj] = tile[x];
            }
      });
      TAManager::get().free("ovo", std::move(tmp_VT2));
      TAManager::get().free("v", std::move(tmp_ones));
    }

    {
      // Populate hbar3[ovv]
      TArray tmp_VT2 = TAmanager.malloc<MatsT>("ovv");
      TArray tmp_ones = TAmanager.malloc_fresh<MatsT>("o");
      tmp_ones.fill(1.0);
      tmp_VT2("i,a,b") = -conj(this->antiSymMoints["vvoo"]("a,b,i,m")) * T2("a,b,i,m") * tmp_ones("m");
      TA::get_default_world().gop.fence();

      TA::foreach_inplace( tmp_VT2, [&](TA::Tensor<MatsT> &tile){
        const auto &lobound = tile.range().lobound();
        const auto &upbound = tile.range().upbound();
        std::size_t x[3] = {0, 0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
          for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
            for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]){
              size_t iab = nO2*nV + (x[0] * nV + x[1]) * nV + x[2];
              hbar3[iab] = tile[x];
            }
      });
      TA::get_default_world().gop.fence();
      TAManager::get().free("ovv", std::move(tmp_VT2));
      TAManager::get().free("o", std::move(tmp_ones));
    }
    TA::get_default_world().gop.fence();

    MatsT *hbar3_copy = CQMemManager::get().malloc<MatsT>(nOV*nOrb_);
    std::copy_n(hbar3.data(), nOV*nOrb_, hbar3_copy);
    std::fill_n(hbar3.data(), nOV*nOrb_, MatsT(0.0));
    MPIAllReduce(hbar3_copy, nOV*nOrb_, hbar3.data(), MPI_COMM_WORLD);
    CQMemManager::get().free(hbar3_copy);
    TA::get_default_world().gop.fence();


    std::cout << "  Finished building denominator in " << std::setprecision(6) << tock(timer) << " s" << std::endl;
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::runCRbatchijk(const TArray &T1_, const TArray &T2_, const TArray &L1_, const TArray &L2_) {

    TAManager &TAmanager = TAManager::get();

    std::cout << "  Initializing..." << std::endl;

    size_t o_tr = (this->intermediates_.nOcc + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
    size_t v_tr = (this->intermediates_.nVir + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);

    std::cout << "  Occ. TileRange: " << o_tr << std::endl;
    std::cout << "  Vir. TileRange: " << v_tr << std::endl;

    // lambda to compute <ijkabc|Hbar|0>
    auto compute_m3 = [&](size_t &i, size_t &j, size_t &k, TArray &m3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      // grab relevant blocks
      TArray tmp = TAmanager.malloc<MatsT>("vvvttt");

      // m_{abc}^{ijk} contribution
      // w_{mb}^{ij}
      block wm_lower = {0,0,i,j};
      block wm_upper = {o_tr,v_tr,i+1,j+1};
      // v_{ac}^{ek}
      block we_lower = {0,0,0,k};
      block we_upper = {v_tr,v_tr,v_tr,k+1};
      // t2_{eb}^{ij}
      block t2e_lower = {0,0,i,j};
      block t2e_upper = {v_tr,v_tr,i+1,j+1};
      //t2_{ac}^{km}
      block t2m_lower = {0,0,k,0};
      block t2m_upper = {v_tr,v_tr,k+1,o_tr};

      auto w_mbij  = W_mbij("m,b,i,j").block(wm_lower,wm_upper);
      auto w_acek  = W_abei("a,c,e,k").block(we_lower,we_upper);
      auto t2_ebij = T2_("e,b,i,j").block(t2e_lower,t2e_upper);
      auto t2_ackm = T2_("a,c,k,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k")  = w_mbij * t2_ackm;
      tmp("a,b,c,i,j,k") += t2_ebij * w_acek;
      tmp("a,b,c,i,j,k") -= t2_ebij * F_me("m,e") * t2_ackm;

      // t_{abc}^{kji} contribution
      // w_{mb}^{kj}
      wm_lower = {0,0,k,j};
      wm_upper = {o_tr,v_tr,k+1,j+1};
      // v_{ac}^{ei}
      we_lower = {0,0,0,i};
      we_upper = {v_tr,v_tr,v_tr,i+1};
      // t2_{eb}^{kj}
      t2e_lower = {0,0,k,j};
      t2e_upper = {v_tr,v_tr,k+1,j+1};
      //t2_{ac}^{im}
      t2m_lower = {0,0,i,0};
      t2m_upper = {v_tr,v_tr,i+1,o_tr};

      auto w_mbkj  = W_mbij("m,b,k,j").block(wm_lower,wm_upper);
      auto w_acei  = W_abei("a,c,e,i").block(we_lower,we_upper);
      auto t2_ebkj = T2_("e,b,k,j").block(t2e_lower,t2e_upper);
      auto t2_acim = T2_("a,c,i,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k") -= w_mbkj * t2_acim;
      tmp("a,b,c,i,j,k") -= t2_ebkj * w_acei;
      tmp("a,b,c,i,j,k") += t2_ebkj * F_me("m,e") * t2_acim;

      // t_{abc}^{ikj} contribution
      // w_{mb}^{ik}
      wm_lower = {0,0,i,k};
      wm_upper = {o_tr,v_tr,i+1,k+1};
      // v_{ac}^{ej}
      we_lower = {0,0,0,j};
      we_upper = {v_tr,v_tr,v_tr,j+1};
      // t2_{eb}^{ik}
      t2e_lower = {0,0,i,k};
      t2e_upper = {v_tr,v_tr,i+1,k+1};
      //t2_{ac}^{jm}
      t2m_lower = {0,0,j,0};
      t2m_upper = {v_tr,v_tr,j+1,o_tr};

      auto w_mbik = W_mbij("m,b,i,k").block(wm_lower,wm_upper);
      auto w_acej = W_abei("a,c,e,j").block(we_lower,we_upper);
      auto t2_ebik = T2_("e,b,i,k").block(t2e_lower,t2e_upper);
      auto t2_acjm = T2_("a,c,j,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k") -= w_mbik * t2_acjm;
      tmp("a,b,c,i,j,k") -= t2_ebik * w_acej;
      tmp("a,b,c,i,j,k") += t2_ebik * F_me("m,e") * t2_acjm;

      // apply A(ac/b)
      m3 = tmp.clone(); // deep copy for safety
      m3("a,b,c,i,j,k") -= tmp("b,a,c,i,j,k");
      m3("a,b,c,i,j,k") -= tmp("a,c,b,i,j,k");

      TAmanager.free("vvvttt", std::move(tmp));
    };

    // lambda to compute <0|(1+L1+L2) Hbar|ijkabc>
    auto compute_l3 = [&](size_t &i, size_t &j, size_t &k, TArray &l3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      TArray tmp = TAmanager.malloc<MatsT>("vvvttt");

      // l_{abc}^{ijk} contribution
      // w_{ij}^{mb} * l_{km}^{ac}
      block wm_lower = {i,j,0,0}, wm_upper = {i+1,j+1,o_tr,v_tr};
      block lm_lower = {0,0,k,0}, lm_upper = {v_tr,v_tr,k+1,o_tr};
      // w_{ek}^{ac} * l_{ij}^{eb}
      block we_lower = {0,k,0,0}, we_upper = {v_tr,k+1,v_tr,v_tr};
      block le_lower = {0,0,i,j}, le_upper = {v_tr,v_tr,i+1,j+1};
      // w_{ij}^{ac} * l_{k}^{b}
      block wdc_lower = {0,0,i,j}, wdc_upper = {v_tr,v_tr,i+1,j+1};
      block l1_lower  = {0,k},     l1_upper  = {v_tr,k+1};
      // f_{k}^{b} * l_{ij}^{ac}
      block fdc_lower = {k,0},     fdc_upper = {k+1,v_tr};
      block ldc_lower = {0,0,i,j}, ldc_upper = {v_tr,v_tr,i+1,j+1};

      auto w_mbij  = W_mnie("i,j,m,b").block(wm_lower,wm_upper);
      auto l2_ackm = L2_("a,c,k,m").block(lm_lower,lm_upper);
      auto w_acek  = W_amef("e,k,a,c").block(we_lower,we_upper);
      auto l2_ebij = L2_("e,b,i,j").block(le_lower,le_upper);
      auto w_acij  = conj(this->antiSymMoints["vvoo"]("a,c,i,j").block(wdc_lower,wdc_upper));
      auto l1_bk   = L1_("b,k").block(l1_lower,l1_upper);
      auto f_kb    = F_me("k,b").block(fdc_lower,fdc_upper);
      auto l2_acij = L2_("a,c,i,j").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k")  = w_mbij * l2_ackm;
      tmp("a,b,c,i,j,k") += l2_ebij * w_acek;
      tmp("a,b,c,i,j,k") -= w_acij * l1_bk;
      tmp("a,b,c,i,j,k") -= l2_acij * f_kb;

      // l_{abc}^{kji} contribution
      // w_{kj}^{mb} * l_{im}^{ac}
      wm_lower = {k,j,0,0}, wm_upper = {k+1,j+1,o_tr,v_tr};
      lm_lower = {0,0,i,0}, lm_upper = {v_tr,v_tr,i+1,o_tr};
      // w_{ei}^{ac} * l_{kj}^{eb}
      we_lower = {0,i,0,0}, we_upper = {v_tr,i+1,v_tr,v_tr};
      le_lower = {0,0,k,j}, le_upper = {v_tr,v_tr,k+1,j+1};
      // w_{kj}^{ac} * l_{i}^{b}
      wdc_lower = {0,0,k,j}, wdc_upper = {v_tr,v_tr,k+1,j+1};
      l1_lower  = {0,i},     l1_upper  = {v_tr,i+1};
      // f_{i}^{b} * l_{kj}^{ac}
      fdc_lower = {i,0},     fdc_upper = {i+1,v_tr};
      ldc_lower = {0,0,k,j}, ldc_upper = {v_tr,v_tr,k+1,j+1};

      auto w_mbkj  = W_mnie("k,j,m,b").block(wm_lower,wm_upper);
      auto l2_acim = L2_("a,c,i,m").block(lm_lower,lm_upper);
      auto w_acei  = W_amef("e,i,a,c").block(we_lower,we_upper);
      auto l2_ebkj = L2_("e,b,k,j").block(le_lower,le_upper);
      auto w_ackj  = conj(this->antiSymMoints["vvoo"]("a,c,k,j").block(wdc_lower,wdc_upper));
      auto l1_bi   = L1_("b,i").block(l1_lower,l1_upper);
      auto f_ib    = F_me("i,b").block(fdc_lower,fdc_upper);
      auto l2_ackj = L2_("a,c,k,j").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k") -= w_mbkj * l2_acim;
      tmp("a,b,c,i,j,k") -= l2_ebkj * w_acei;
      tmp("a,b,c,i,j,k") += w_ackj * l1_bi;
      tmp("a,b,c,i,j,k") += l2_ackj * f_ib;

      // l_{abc}^{ikj} contribution
      // w_{ik}^{mb} * l_{jm}^{ac}
      wm_lower = {i,k,0,0}, wm_upper = {i+1,k+1,o_tr,v_tr};
      lm_lower = {0,0,j,0}, lm_upper = {v_tr,v_tr,j+1,o_tr};
      // w_{ej}^{ac} * l_{ik}^{eb}
      we_lower = {0,j,0,0}, we_upper = {v_tr,j+1,v_tr,v_tr};
      le_lower = {0,0,i,k}, le_upper = {v_tr,v_tr,i+1,k+1};
      // w_{ik}^{ac} * l_{j}^{b}
      wdc_lower = {0,0,i,k}, wdc_upper = {v_tr,v_tr,i+1,k+1};
      l1_lower  = {0,j},     l1_upper  = {v_tr,j+1};
      // f_{j}^{b} * l_{ik}^{ac}
      fdc_lower = {j,0},     fdc_upper = {j+1,v_tr};
      ldc_lower = {0,0,i,k}, ldc_upper = {v_tr,v_tr,i+1,k+1};

      auto w_mbik  = W_mnie("i,k,m,b").block(wm_lower,wm_upper);
      auto l2_acjm = L2_("a,c,j,m").block(lm_lower,lm_upper);
      auto w_acej  = W_amef("e,j,a,c").block(we_lower,we_upper);
      auto l2_ebik = L2_("e,b,i,k").block(le_lower,le_upper);
      auto w_acik  = conj(this->antiSymMoints["vvoo"]("a,c,i,k").block(wdc_lower,wdc_upper));
      auto l1_bj   = L1_("b,j").block(l1_lower,l1_upper);
      auto f_jb    = F_me("j,b").block(fdc_lower,fdc_upper);
      auto l2_acik = L2_("a,c,i,k").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k") -= w_mbik * l2_acjm;
      tmp("a,b,c,i,j,k") -= l2_ebik * w_acej;
      tmp("a,b,c,i,j,k") += w_acik * l1_bj;
      tmp("a,b,c,i,j,k") += l2_acik * f_jb;

      // apply A(ac/b)
      l3  = tmp.clone(); // deep copy for safety
      l3("a,b,c,i,j,k") -= tmp("b,a,c,i,j,k");
      l3("a,b,c,i,j,k") -= tmp("a,c,b,i,j,k");
      
      TAmanager.free("vvvttt", std::move(tmp));
    };

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

          TArray m3_ijk   = TAmanager.malloc<MatsT>("vvvttt");
          TArray l3_ijk_d = TAmanager.malloc<MatsT>("vvvttt");

          compute_m3(i,j,k,m3_ijk);
          compute_l3(i,j,k,l3_ijk_d);

          size_t nO = this->intermediates_.nOcc;
          size_t nV = this->intermediates_.nVir;
          size_t nO2 = nO * nO;
          size_t nOV = nO * nV;

          TArray l3_ijk_a = TAmanager.malloc<MatsT>("vvvttt");
          TArray l3_ijk_b = TAmanager.malloc<MatsT>("vvvttt");
          TArray l3_ijk_c = TAmanager.malloc<MatsT>("vvvttt");

          l3_ijk_a = l3_ijk_d.clone(); // another v3t3 object
          l3_ijk_b = l3_ijk_d.clone(); // another v3t3 object
          l3_ijk_c = l3_ijk_d.clone(); // another v3t3 object

          // offset i,j,k of blocked TiledArray object, otherwise wrong eps[x] index
          size_t i_off = i*this->ccSettings_.blksize;
          size_t j_off = j*this->ccSettings_.blksize;
          size_t k_off = k*this->ccSettings_.blksize;

          // A denominator: <ijkabc|F|ijkabc>
          TA::foreach_inplace(l3_ijk_a, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t a = x[0]+nO, b = x[1]+nO, c = x[2]+nO;
                        size_t ii = x[3]+i_off, jj = x[4]+j_off, kk = x[5]+k_off;

                        // D3 construction
                        // <ijkabc| Fock |ijkabc> = f(v) - f(o), invert sign
                        MatsT denom = this->intermediates_.eps[ii] + this->intermediates_.eps[jj] + this->intermediates_.eps[kk]
                                    - this->intermediates_.eps[a] - this->intermediates_.eps[b] - this->intermediates_.eps[c];
                        tile[x] /= denom;
                      }
          });

          // B denominator: <ijkabc|H1|ijkabc>
          TA::foreach_inplace(l3_ijk_b, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t a = x[0]+nO, b = x[1]+nO, c = x[2]+nO;
                        size_t ii = x[3]+i_off, jj = x[4]+j_off, kk = x[5]+k_off;

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[ii] + hbar1[jj] + hbar1[kk]
                                    - hbar1[a] - hbar1[b] - hbar1[c];
                        tile[x] /= denom;
                      }
          });

          // C denominator: <ijkabc|H1+H2|ijkabc>
          TA::foreach_inplace(l3_ijk_c, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t a = x[0]+nO, b = x[1]+nO, c = x[2]+nO;
                        size_t ii = x[3]+i_off, jj = x[4]+j_off, kk = x[5]+k_off;

                        size_t ij = ii * nO + jj;
                        size_t ik = ii * nO + kk;
                        size_t jk = jj * nO + kk;

                        size_t ia = nO2 + ii * nV + x[0];
                        size_t ib = nO2 + ii * nV + x[1];
                        size_t ic = nO2 + ii * nV + x[2];
                        size_t ja = nO2 + jj * nV + x[0];
                        size_t jb = nO2 + jj * nV + x[1];
                        size_t jc = nO2 + jj * nV + x[2];
                        size_t ka = nO2 + kk * nV + x[0];
                        size_t kb = nO2 + kk * nV + x[1];
                        size_t kc = nO2 + kk * nV + x[2];

                        size_t ab = (nO2 + nOV) + x[0] * nV + (x[1]);
                        size_t ac = (nO2 + nOV) + x[0] * nV + (x[2]);
                        size_t bc = (nO2 + nOV) + x[1] * nV + (x[2]);

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[ii] + hbar1[jj] + hbar1[kk]
                                    - hbar1[a] - hbar1[b] - hbar1[c];
                        // <ijkabc| Hbar2 |ijkabc> = h(oo) + h(vv) + h(ov), invert sign
                        denom -= hbar2[ij] + hbar2[ik] + hbar2[jk];
                        denom -= hbar2[ab] + hbar2[ac] + hbar2[bc];
                        denom -= hbar2[ia] + hbar2[ib] + hbar2[ic];
                        denom -= hbar2[ja] + hbar2[jb] + hbar2[jc];
                        denom -= hbar2[ka] + hbar2[kb] + hbar2[kc];
                        tile[x] /= denom;
                      }
          });

          // D denominator: <ijkabc|H1+H2+H3|ijkabc
          TA::foreach_inplace(l3_ijk_d, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t a = x[0]+nO, b = x[1]+nO, c = x[2]+nO;
                        size_t ii = x[3]+i_off, jj = x[4]+j_off, kk = x[5]+k_off;

                        size_t ij = ii * nO + jj;
                        size_t ik = ii * nO + kk;
                        size_t jk = jj * nO + kk;

                        size_t ia = nO2 + ii * nV + x[0];
                        size_t ib = nO2 + ii * nV + x[1];
                        size_t ic = nO2 + ii * nV + x[2];
                        size_t ja = nO2 + jj * nV + x[0];
                        size_t jb = nO2 + jj * nV + x[1];
                        size_t jc = nO2 + jj * nV + x[2];
                        size_t ka = nO2 + kk * nV + x[0];
                        size_t kb = nO2 + kk * nV + x[1];
                        size_t kc = nO2 + kk * nV + x[2];

                        size_t ab = (nO2 + nOV) + x[0] * nV + (x[1]);
                        size_t ac = (nO2 + nOV) + x[0] * nV + (x[2]);
                        size_t bc = (nO2 + nOV) + x[1] * nV + (x[2]);

                        size_t iaj = (ii * nV + x[0]) * nO + jj;
                        size_t ibj = (ii * nV + x[1]) * nO + jj;
                        size_t icj = (ii * nV + x[2]) * nO + jj;
                        size_t iak = (ii * nV + x[0]) * nO + kk;
                        size_t ibk = (ii * nV + x[1]) * nO + kk;
                        size_t ick = (ii * nV + x[2]) * nO + kk;
                        size_t jak = (jj * nV + x[0]) * nO + kk;
                        size_t jbk = (jj * nV + x[1]) * nO + kk;
                        size_t jck = (jj * nV + x[2]) * nO + kk;

                        size_t iab = nO2*nV + (ii * nV + x[0]) * nV + x[1];
                        size_t iac = nO2*nV + (ii * nV + x[0]) * nV + x[2];
                        size_t ibc = nO2*nV + (ii * nV + x[1]) * nV + x[2];
                        size_t jab = nO2*nV + (jj * nV + x[0]) * nV + x[1];
                        size_t jac = nO2*nV + (jj * nV + x[0]) * nV + x[2];
                        size_t jbc = nO2*nV + (jj * nV + x[1]) * nV + x[2];
                        size_t kab = nO2*nV + (kk * nV + x[0]) * nV + x[1];
                        size_t kac = nO2*nV + (kk * nV + x[0]) * nV + x[2];
                        size_t kbc = nO2*nV + (kk * nV + x[1]) * nV + x[2];

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[ii] + hbar1[jj] + hbar1[kk]
                                    - hbar1[a] - hbar1[b] - hbar1[c];
                        // <ijkabc| Hbar2 |ijkabc> = h(oo) + h(vv) + h(ov), invert sign
                        denom -= hbar2[ij] + hbar2[ik] + hbar2[jk];
                        denom -= hbar2[ab] + hbar2[ac] + hbar2[bc];
                        denom -= hbar2[ia] + hbar2[ib] + hbar2[ic];
                        denom -= hbar2[ja] + hbar2[jb] + hbar2[jc];
                        denom -= hbar2[ka] + hbar2[kb] + hbar2[kc];
                        // <ijkabc| Hbar3 |ijkabc> = h(ovv) - h(oov), invert sign
                        denom += hbar3[iaj] + hbar3[ibj] + hbar3[icj];
                        denom += hbar3[iak] + hbar3[ibk] + hbar3[ick];
                        denom += hbar3[jak] + hbar3[jbk] + hbar3[jck];
                        denom -= hbar3[iab] + hbar3[iac] + hbar3[ibc];
                        denom -= hbar3[jab] + hbar3[jac] + hbar3[jbc];
                        denom -= hbar3[kab] + hbar3[kac] + hbar3[kbc];
                        tile[x] /= denom;
                      }
          });

          MatsT tmp_en_a = l3_ijk_a("a,b,c,i,j,k").dot(m3_ijk("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();
          MatsT tmp_en_b = l3_ijk_b("a,b,c,i,j,k").dot(m3_ijk("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();
          MatsT tmp_en_c = l3_ijk_c("a,b,c,i,j,k").dot(m3_ijk("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();
          MatsT tmp_en_d = l3_ijk_d("a,b,c,i,j,k").dot(m3_ijk("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();
          if ((i==j) && (j==k)) {
            // ijk -> 3!, abc -> 3!
            tmp_en_a *= 1.0/36.0;
            tmp_en_b *= 1.0/36.0;
            tmp_en_c *= 1.0/36.0;
            tmp_en_d *= 1.0/36.0;
          }
          else if ((i==j) || (j==k) || (i==k)) {
            // ijk -> 2!, abc -> 3!
            tmp_en_a *= 1.0/12.0;
            tmp_en_b *= 1.0/12.0;
            tmp_en_c *= 1.0/12.0;
            tmp_en_d *= 1.0/12.0;
          }
          else {
            // ijk -> 1!, abc -> 3!
            tmp_en_a *= 1.0/6.0;
            tmp_en_b *= 1.0/6.0;
            tmp_en_c *= 1.0/6.0;
            tmp_en_d *= 1.0/6.0;
          }
          this->CRCC23Energy_A += tmp_en_a;
          this->CRCC23Energy_B += tmp_en_b;
          this->CRCC23Energy_C += tmp_en_c;
          this->CRCC23Energy_D += tmp_en_d;

          TAmanager.free("vvvttt", std::move(m3_ijk));
          TAmanager.free("vvvttt", std::move(l3_ijk_a));
          TAmanager.free("vvvttt", std::move(l3_ijk_b));
          TAmanager.free("vvvttt", std::move(l3_ijk_c));
          TAmanager.free("vvvttt", std::move(l3_ijk_d));

          std::cout << "  i: " << i << " j: " << j << " k: " << k << " done in " << tock(ijk_start) << "s"
                    << " from global_iter " << global_iter << " in rank " << rank << std::endl;
        }

    this_world.gop.fence();
    global_world.gop.fence();

    TA::set_default_world(global_world);

    // sum over contribution if MPI parallel is run
    if (this->ccSettings_.triplesMPI && size >1) {
    global_world.gop.sum(this->CRCC23Energy_A);
    global_world.gop.sum(this->CRCC23Energy_B);
    global_world.gop.sum(this->CRCC23Energy_C);
    global_world.gop.sum(this->CRCC23Energy_D);
    }

  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::runCRbatchabc(const TArray &T1_, const TArray &T2_, const TArray &L1_, const TArray &L2_) {

    TAManager &TAmanager = TAManager::get();

    std::cout << "  Initializing..." << std::endl;

    size_t o_tr = (this->intermediates_.nOcc + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);
    size_t v_tr = (this->intermediates_.nVir + this->ccSettings_.blksize - 1)/(this->ccSettings_.blksize);

    std::cout << "  Occ. TileRange: " << o_tr << std::endl;
    std::cout << "  Vir. TileRange: " << v_tr << std::endl;

    // lambda to compute <ijkabc|Hbar|0>
    auto compute_m3 = [&](size_t &a, size_t &b, size_t &c, TArray &m3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      // grab relevant blocks
      TArray tmp = TAmanager.malloc<MatsT>("tttooo");

      // m_{abc}^{ijk} contribution
      // w_{mb}^{ij}
      block wm_lower = {0,b,0,0};
      block wm_upper = {o_tr,b+1,o_tr,o_tr};
      // v_{ac}^{ek}
      block we_lower = {a,c,0,0};
      block we_upper = {a+1,c+1,v_tr,o_tr};
      // t2_{eb}^{ij}
      block t2e_lower = {0,b,0,0};
      block t2e_upper = {v_tr,b+1,o_tr,o_tr};
      //t2_{ac}^{km}
      block t2m_lower = {a,c,0,0};
      block t2m_upper = {a+1,c+1,o_tr,o_tr};

      auto w_mbij  = W_mbij("m,b,i,j").block(wm_lower,wm_upper);
      auto w_acek  = W_abei("a,c,e,k").block(we_lower,we_upper);
      auto t2_ebij = T2_("e,b,i,j").block(t2e_lower,t2e_upper);
      auto t2_ackm = T2_("a,c,k,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k")  = w_mbij * t2_ackm;
      tmp("a,b,c,i,j,k") += t2_ebij * w_acek;
      tmp("a,b,c,i,j,k") -= t2_ebij * F_me("m,e") * t2_ackm;

      // t_{bac}^{ijk} contribution
      // w_{ma}^{ij}
      wm_lower = {0,a,0,0};
      wm_upper = {o_tr,a+1,o_tr,o_tr};
      // w_{bc}^{ek}
      we_lower = {b,c,0,0};
      we_upper = {b+1,c+1,v_tr,o_tr};
      // t2_{eb}^{ij}
      t2e_lower = {0,a,0,0};
      t2e_upper = {v_tr,a+1,o_tr,o_tr};
      //t2_{ac}^{km}
      t2m_lower = {b,c,0,0};
      t2m_upper = {b+1,c+1,o_tr,o_tr};

      auto w_maij  = W_mbij("m,a,i,j").block(wm_lower,wm_upper);
      auto w_bcek  = W_abei("b,c,e,k").block(we_lower,we_upper);
      auto t2_eaij = T2_("e,a,i,j").block(t2e_lower,t2e_upper);
      auto t2_bckm = T2_("b,c,k,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k") -= w_maij * t2_bckm;
      tmp("a,b,c,i,j,k") -= t2_eaij * w_bcek;
      tmp("a,b,c,i,j,k") += t2_eaij * F_me("m,e") * t2_bckm;

      // t_{acb}^{ijk} contribution
      // w_{mc}^{ij}
      wm_lower = {0,c,0,0};
      wm_upper = {o_tr,c+1,o_tr,o_tr};
      // v_{ab}^{ek}
      we_lower = {a,b,0,0};
      we_upper = {a+1,b+1,v_tr,o_tr};
      // t2_{ec}^{ij}
      t2e_lower = {0,c,0,0};
      t2e_upper = {v_tr,c+1,o_tr,o_tr};
      //t2_{ab}^{km}
      t2m_lower = {a,b,0,0};
      t2m_upper = {a+1,b+1,o_tr,o_tr};

      auto w_mcij  = W_mbij("m,c,i,j").block(wm_lower,wm_upper);
      auto w_abek  = W_abei("a,b,e,k").block(we_lower,we_upper);
      auto t2_ecij = T2_("e,c,i,j").block(t2e_lower,t2e_upper);
      auto t2_abkm = T2_("a,b,k,m").block(t2m_lower,t2m_upper);

      tmp("a,b,c,i,j,k") -= w_mcij * t2_abkm;
      tmp("a,b,c,i,j,k") -= t2_ecij * w_abek;
      tmp("a,b,c,i,j,k") += t2_ecij * F_me("m,e") * t2_abkm;

      // apply A(ij/k)
      m3 = tmp.clone(); // deep copy for safety
      m3("a,b,c,i,j,k") -= tmp("a,b,c,k,j,i");
      m3("a,b,c,i,j,k") -= tmp("a,b,c,i,k,j");

      TAmanager.free("tttooo", std::move(tmp));
    };

    // lambda to compute <0|(1+L1+L2) Hbar|ijkabc>
    auto compute_l3 = [&](size_t &a, size_t &b, size_t &c, TArray &l3) {

      // for blocking purposes
      typedef std::vector<size_t> block;

      TArray tmp = TAmanager.malloc<MatsT>("tttooo");

      // l_{abc}^{ijk} contribution
      // w_{ij}^{mb} * l_{km}^{ac}
      block wm_lower = {0,0,0,b}, wm_upper = {o_tr,o_tr,o_tr,b+1};
      block lm_lower = {a,c,0,0}, lm_upper = {a+1,c+1,o_tr,o_tr};
      // w_{ek}^{ac} * l_{ij}^{eb}
      block we_lower = {0,0,a,c}, we_upper = {v_tr,o_tr,a+1,c+1};
      block le_lower = {0,b,0,0}, le_upper = {v_tr,b+1,o_tr,o_tr};
      // w_{ij}^{ac} * l_{k}^{b}
      block wdc_lower = {a,c,0,0}, wdc_upper = {a+1,c+1,o_tr,o_tr};
      block l1_lower  = {b,0},     l1_upper  = {b+1,o_tr};
      // f_{k}^{b} * l_{ij}^{ac}
      block fdc_lower = {0,b},     fdc_upper = {o_tr,b+1};
      block ldc_lower = {a,c,0,0}, ldc_upper = {a+1,c+1,o_tr,o_tr};

      auto w_mbij  = W_mnie("i,j,m,b").block(wm_lower,wm_upper);
      auto l2_ackm = L2_("a,c,k,m").block(lm_lower,lm_upper);
      auto w_acek  = W_amef("e,k,a,c").block(we_lower,we_upper);
      auto l2_ebij = L2_("e,b,i,j").block(le_lower,le_upper);
      auto w_acij  = conj(this->antiSymMoints["vvoo"]("a,c,i,j").block(wdc_lower,wdc_upper));
      auto l1_bk   = L1_("b,k").block(l1_lower,l1_upper);
      auto f_kb    = F_me("k,b").block(fdc_lower,fdc_upper);
      auto l2_acij = L2_("a,c,i,j").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k")  = w_mbij * l2_ackm;
      tmp("a,b,c,i,j,k") += l2_ebij * w_acek;
      tmp("a,b,c,i,j,k") -= w_acij * l1_bk;
      tmp("a,b,c,i,j,k") -= l2_acij * f_kb;

      // l_{bac}^{ijk} contribution
      // w_{ij}^{ma} * l_{km}^{bc}
      wm_lower = {0,0,0,a}, wm_upper = {o_tr,o_tr,o_tr,a+1};
      lm_lower = {b,c,0,0}, lm_upper = {b+1,c+1,o_tr,o_tr};
      // w_{ek}^{bc} * l_{ij}^{ea}
      we_lower = {0,0,b,c}, we_upper = {v_tr,o_tr,b+1,c+1};
      le_lower = {0,a,0,0}, le_upper = {v_tr,a+1,o_tr,o_tr};
      // w_{ij}^{bc} * l_{k}^{a}
      wdc_lower = {b,c,0,0}, wdc_upper = {b+1,c+1,o_tr,o_tr};
      l1_lower  = {a,0},     l1_upper  = {a+1,o_tr};
      // f_{k}^{a} * l_{ij}^{bc}
      fdc_lower = {0,a},     fdc_upper = {o_tr,a+1};
      ldc_lower = {b,c,0,0}, ldc_upper = {b+1,c+1,o_tr,o_tr};

      auto w_maij  = W_mnie("i,j,m,a").block(wm_lower,wm_upper);
      auto l2_bckm = L2_("b,c,k,m").block(lm_lower,lm_upper);
      auto w_bcek  = W_amef("e,k,b,c").block(we_lower,we_upper);
      auto l2_eaij = L2_("e,a,i,j").block(le_lower,le_upper);
      auto w_bcij  = conj(this->antiSymMoints["vvoo"]("b,c,i,j").block(wdc_lower,wdc_upper));
      auto l1_ak   = L1_("a,k").block(l1_lower,l1_upper);
      auto f_ka    = F_me("k,a").block(fdc_lower,fdc_upper);
      auto l2_bcij = L2_("b,c,i,j").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k") -= w_maij * l2_bckm;
      tmp("a,b,c,i,j,k") -= l2_eaij * w_bcek;
      tmp("a,b,c,i,j,k") += w_bcij * l1_ak;
      tmp("a,b,c,i,j,k") += l2_bcij * f_ka;

      // l_{acb}^{ijk} contribution
      // w_{ij}^{mb} * l_{km}^{ac}
      wm_lower = {0,0,0,c}, wm_upper = {o_tr,o_tr,o_tr,c+1};
      lm_lower = {a,b,0,0}, lm_upper = {a+1,b+1,o_tr,o_tr};
      // w_{ek}^{ac} * l_{ij}^{eb}
      we_lower = {0,0,a,b}, we_upper = {v_tr,o_tr,a+1,b+1};
      le_lower = {0,c,0,0}, le_upper = {v_tr,c+1,o_tr,o_tr};
      // w_{ij}^{ac} * l_{k}^{b}
      wdc_lower = {a,b,0,0}, wdc_upper = {a+1,b+1,o_tr,o_tr};
      l1_lower  = {c,0},     l1_upper  = {c+1,o_tr};
      // f_{k}^{b} * l_{ij}^{ac}
      fdc_lower = {0,c},     fdc_upper = {o_tr,c+1};
      ldc_lower = {a,b,0,0}, ldc_upper = {a+1,b+1,o_tr,o_tr};

      auto w_mcij  = W_mnie("i,j,m,c").block(wm_lower,wm_upper);
      auto l2_abkm = L2_("a,b,k,m").block(lm_lower,lm_upper);
      auto w_abek  = W_amef("e,k,a,b").block(we_lower,we_upper);
      auto l2_ecij = L2_("e,c,i,j").block(le_lower,le_upper);
      auto w_abij  = conj(this->antiSymMoints["vvoo"]("a,b,i,j").block(wdc_lower,wdc_upper));
      auto l1_ck   = L1_("c,k").block(l1_lower,l1_upper);
      auto f_kc    = F_me("k,c").block(fdc_lower,fdc_upper);
      auto l2_abij = L2_("a,b,i,j").block(ldc_lower,ldc_upper);

      tmp("a,b,c,i,j,k") -= w_mcij * l2_abkm;
      tmp("a,b,c,i,j,k") -= l2_ecij * w_abek;
      tmp("a,b,c,i,j,k") += w_abij * l1_ck;
      tmp("a,b,c,i,j,k") += l2_abij * f_kc;

      // apply A(ij/k)
      l3  = tmp.clone(); // deep copy for safety
      l3("a,b,c,i,j,k") -= tmp("a,b,c,k,j,i");
      l3("a,b,c,i,j,k") -= tmp("a,b,c,i,k,j");

      TAmanager.free("tttooo", std::move(tmp));
    };

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
    std::cout << "  Running " << v_tr*(v_tr+1)*(v_tr+2)/6 << " iterations among "
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

          TArray m3_abc   = TAmanager.malloc<MatsT>("tttooo");
          TArray l3_abc_d = TAmanager.malloc<MatsT>("tttooo");

          compute_m3(a,b,c,m3_abc);
          compute_l3(a,b,c,l3_abc_d);

          size_t nO = this->intermediates_.nOcc;
          size_t nV = this->intermediates_.nVir;
          size_t nO2 = nO * nO;
          size_t nOV = nO * nV;

          TArray l3_abc_a = TAmanager.malloc<MatsT>("tttooo");
          TArray l3_abc_b = TAmanager.malloc<MatsT>("tttooo");
          TArray l3_abc_c = TAmanager.malloc<MatsT>("tttooo");

          l3_abc_a = l3_abc_d.clone(); // another t3o3 object
          l3_abc_b = l3_abc_d.clone(); // another t3o3 object
          l3_abc_c = l3_abc_d.clone(); // another t3o3 object

          // offset a,b,c of blocked TiledArray object, otherwise wrong eps[x] index
          size_t a_off = a*this->ccSettings_.blksize;
          size_t b_off = b*this->ccSettings_.blksize;
          size_t c_off = c*this->ccSettings_.blksize;

          // A denominator: <ijkabc|F|ijkabc>
          TA::foreach_inplace(l3_abc_a, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t aa = x[0]+nO+a_off, bb = x[1]+nO+b_off, cc = x[2]+nO+c_off;
                        size_t i = x[3], j = x[4], k = x[5];

                        // D3 construction
                        // <ijkabc| Fock |ijkabc> = f(v) - f(o), invert sign
                        MatsT denom = this->intermediates_.eps[i] + this->intermediates_.eps[j] + this->intermediates_.eps[k]
                                    - this->intermediates_.eps[aa] - this->intermediates_.eps[bb] - this->intermediates_.eps[cc];
                        tile[x] /= denom;
                      }
          });

          // B denominator: <ijkabc|H1|ijkabc>
          TA::foreach_inplace(l3_abc_b, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t aa = x[0]+nO+a_off, bb = x[1]+nO+b_off, cc = x[2]+nO+c_off;
                        size_t i = x[3], j = x[4], k = x[5];

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[i] + hbar1[j] + hbar1[k]
                                    - hbar1[aa] - hbar1[bb] - hbar1[cc];
                        tile[x] /= denom;
                      }
          });

          // C denominator: <ijkabc|H1+H2|ijkabc>
          TA::foreach_inplace(l3_abc_c, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t aa = x[0]+nO+a_off, bb = x[1]+nO+b_off, cc = x[2]+nO+c_off;
                        size_t i = x[3], j = x[4], k = x[5];

                        size_t ij = x[3] * nO + x[4];
                        size_t ik = x[3] * nO + x[5];
                        size_t jk = x[4] * nO + x[5];

                        size_t ia = nO2 + x[3] * nV + (x[0] + a_off);
                        size_t ib = nO2 + x[3] * nV + (x[1] + b_off);
                        size_t ic = nO2 + x[3] * nV + (x[2] + c_off);
                        size_t ja = nO2 + x[4] * nV + (x[0] + a_off);
                        size_t jb = nO2 + x[4] * nV + (x[1] + b_off);
                        size_t jc = nO2 + x[4] * nV + (x[2] + c_off);
                        size_t ka = nO2 + x[5] * nV + (x[0] + a_off);
                        size_t kb = nO2 + x[5] * nV + (x[1] + b_off);
                        size_t kc = nO2 + x[5] * nV + (x[2] + c_off);

                        size_t ab = (nO2 + nOV) + (x[0] + a_off) * nV + (x[1] + b_off);
                        size_t ac = (nO2 + nOV) + (x[0] + a_off) * nV + (x[2] + c_off);
                        size_t bc = (nO2 + nOV) + (x[1] + b_off) * nV + (x[2] + c_off);

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[i] + hbar1[j] + hbar1[k]
                                    - hbar1[aa] - hbar1[bb] - hbar1[cc];
                        // <ijkabc| Hbar2 |ijkabc> = h(oo) + h(vv) + h(ov), invert sign
                        denom -= hbar2[ij] + hbar2[ik] + hbar2[jk];
                        denom -= hbar2[ab] + hbar2[ac] + hbar2[bc];
                        denom -= hbar2[ia] + hbar2[ib] + hbar2[ic];
                        denom -= hbar2[ja] + hbar2[jb] + hbar2[jc];
                        denom -= hbar2[ka] + hbar2[kb] + hbar2[kc];
                        tile[x] /= denom;
                      }
          });

          // D denominator: <ijkabc|H1+H2+H3|ijkabc
          TA::foreach_inplace(l3_abc_d, [&, this](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0,};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){

                        size_t aa = x[0]+nO+a_off, bb = x[1]+nO+b_off, cc = x[2]+nO+c_off;
                        size_t i = x[3], j = x[4], k = x[5];

                        size_t ij = x[3] * nO + x[4];
                        size_t ik = x[3] * nO + x[5];
                        size_t jk = x[4] * nO + x[5];

                        size_t ia = nO2 + x[3] * nV + (x[0] + a_off);
                        size_t ib = nO2 + x[3] * nV + (x[1] + b_off);
                        size_t ic = nO2 + x[3] * nV + (x[2] + c_off);
                        size_t ja = nO2 + x[4] * nV + (x[0] + a_off);
                        size_t jb = nO2 + x[4] * nV + (x[1] + b_off);
                        size_t jc = nO2 + x[4] * nV + (x[2] + c_off);
                        size_t ka = nO2 + x[5] * nV + (x[0] + a_off);
                        size_t kb = nO2 + x[5] * nV + (x[1] + b_off);
                        size_t kc = nO2 + x[5] * nV + (x[2] + c_off);

                        size_t ab = (nO2 + nOV) + (x[0] + a_off) * nV + (x[1] + b_off);
                        size_t ac = (nO2 + nOV) + (x[0] + a_off) * nV + (x[2] + c_off);
                        size_t bc = (nO2 + nOV) + (x[1] + b_off) * nV + (x[2] + c_off);

                        size_t iaj = (x[3] * nV + (x[0] + a_off)) * nO + x[4];
                        size_t ibj = (x[3] * nV + (x[1] + b_off)) * nO + x[4];
                        size_t icj = (x[3] * nV + (x[2] + c_off)) * nO + x[4];
                        size_t iak = (x[3] * nV + (x[0] + a_off)) * nO + x[5];
                        size_t ibk = (x[3] * nV + (x[1] + b_off)) * nO + x[5];
                        size_t ick = (x[3] * nV + (x[2] + c_off)) * nO + x[5];
                        size_t jak = (x[4] * nV + (x[0] + a_off)) * nO + x[5];
                        size_t jbk = (x[4] * nV + (x[1] + b_off)) * nO + x[5];
                        size_t jck = (x[4] * nV + (x[2] + c_off)) * nO + x[5];

                        size_t iab = nO2*nV + (x[3] * nV + (x[0] + a_off)) * nV + (x[1] + b_off);
                        size_t iac = nO2*nV + (x[3] * nV + (x[0] + a_off)) * nV + (x[2] + c_off);
                        size_t ibc = nO2*nV + (x[3] * nV + (x[1] + b_off)) * nV + (x[2] + c_off);
                        size_t jab = nO2*nV + (x[4] * nV + (x[0] + a_off)) * nV + (x[1] + b_off);
                        size_t jac = nO2*nV + (x[4] * nV + (x[0] + a_off)) * nV + (x[2] + c_off);
                        size_t jbc = nO2*nV + (x[4] * nV + (x[1] + b_off)) * nV + (x[2] + c_off);
                        size_t kab = nO2*nV + (x[5] * nV + (x[0] + a_off)) * nV + (x[1] + b_off);
                        size_t kac = nO2*nV + (x[5] * nV + (x[0] + a_off)) * nV + (x[2] + c_off);
                        size_t kbc = nO2*nV + (x[5] * nV + (x[1] + b_off)) * nV + (x[2] + c_off);

                        // D3 construction
                        // <ijkabc| Hbar1 |ijkabc> = h(v) - h(o), invert sign
                        MatsT denom = hbar1[i] + hbar1[j] + hbar1[k]
                                    - hbar1[aa] - hbar1[bb] - hbar1[cc];
                        // <ijkabc| Hbar2 |ijkabc> = h(oo) + h(vv) + h(ov), invert sign
                        denom -= hbar2[ij] + hbar2[ik] + hbar2[jk];
                        denom -= hbar2[ab] + hbar2[ac] + hbar2[bc];
                        denom -= hbar2[ia] + hbar2[ib] + hbar2[ic];
                        denom -= hbar2[ja] + hbar2[jb] + hbar2[jc];
                        denom -= hbar2[ka] + hbar2[kb] + hbar2[kc];
                        // <ijkabc| Hbar3 |ijkabc> = h(ovv) - h(oov), invert sign
                        denom += hbar3[iaj] + hbar3[ibj] + hbar3[icj];
                        denom += hbar3[iak] + hbar3[ibk] + hbar3[ick];
                        denom += hbar3[jak] + hbar3[jbk] + hbar3[jck];
                        denom -= hbar3[iab] + hbar3[iac] + hbar3[ibc];
                        denom -= hbar3[jab] + hbar3[jac] + hbar3[jbc];
                        denom -= hbar3[kab] + hbar3[kac] + hbar3[kbc];
                        tile[x] /= denom;
                      }
          });

          MatsT tmp_en_a = l3_abc_a("a,b,c,i,j,k").dot(m3_abc("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();    
          MatsT tmp_en_b = l3_abc_b("a,b,c,i,j,k").dot(m3_abc("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();    
          MatsT tmp_en_c = l3_abc_c("a,b,c,i,j,k").dot(m3_abc("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();    
          MatsT tmp_en_d = l3_abc_d("a,b,c,i,j,k").dot(m3_abc("a,b,c,i,j,k"));
          TA::get_default_world().gop.fence();    
          if ((a==b) && (b==c)) {
            // abc -> 3!, ijk -> 3!
            tmp_en_a *= 1.0/36.0;
            tmp_en_b *= 1.0/36.0;
            tmp_en_c *= 1.0/36.0;
            tmp_en_d *= 1.0/36.0;
          }
          else if ((a==b) || (b==c) || (a==c)) {
            // abc -> 2!, ijk -> 3!
            tmp_en_a *= 1.0/12.0;
            tmp_en_b *= 1.0/12.0;
            tmp_en_c *= 1.0/12.0;
            tmp_en_d *= 1.0/12.0;
          }
          else {
            // abc -> 1!, ijk -> 3!
            tmp_en_a *= 1.0/6.0;
            tmp_en_b *= 1.0/6.0;
            tmp_en_c *= 1.0/6.0;
            tmp_en_d *= 1.0/6.0;
          }
          this->CRCC23Energy_A += tmp_en_a;
          this->CRCC23Energy_B += tmp_en_b;
          this->CRCC23Energy_C += tmp_en_c;
          this->CRCC23Energy_D += tmp_en_d;

          TAmanager.free("tttooo", std::move(m3_abc));
          TAmanager.free("tttooo", std::move(l3_abc_a));
          TAmanager.free("tttooo", std::move(l3_abc_b));
          TAmanager.free("tttooo", std::move(l3_abc_c));
          TAmanager.free("tttooo", std::move(l3_abc_d));

          std::cout << "  a: " << a << " b: " << b << " c: " << c << " done in " << tock(abc_start) << "s"
                    << " from global_iter " << global_iter << " in rank " << rank << std::endl;
        }

    this_world.gop.fence();
    global_world.gop.fence();

    TA::set_default_world(global_world);

    // sum over contribution if MPI parallel is run
    if (this->ccSettings_.triplesMPI && size >1) {
      global_world.gop.sum(this->CRCC23Energy_A);
      global_world.gop.sum(this->CRCC23Energy_B);
      global_world.gop.sum(this->CRCC23Energy_C);
      global_world.gop.sum(this->CRCC23Energy_D);
    }
  }

}; // namespace
