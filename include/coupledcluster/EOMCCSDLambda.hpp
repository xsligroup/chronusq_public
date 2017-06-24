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
#include <util/math.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <functional>
#include <util/timer.hpp>

namespace ChronusQ{

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::initilizeLambda() {

    Lg_.zeroBody() = 1.0;
    Lg_.oneBody()("a,i") = conj(T1_("a,i"));
    Lg_.twoBody()("a,b,i,j") = conj(T2_("a,b,i,j"));

    TAManager &TAmanager = TAManager::get();
    if (not G_ae.is_initialized()){
      G_ae = TAmanager.malloc<MatsT>("vv");
    }
    if (not G_mi.is_initialized()){
      G_mi = TAmanager.malloc<MatsT>("oo");
    }

  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::updateG_ae(const TArray &L2, TArray &G_ae) const {
    G_ae("a,e") = - 0.5 * T2_("e,f,m,n") * L2("a,f,m,n");
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::updateG_mi(const TArray &L2, TArray &G_mi) const {
    G_mi("m,i") = 0.5 * T2_("e,f,m,n") * L2("e,f,i,n");
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formL1_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                                          TArray &tildeL1, bool subtractDiagonal) const {

    tildeL1("a,i") = F_ae("e,a") * L1("e,i");
    tildeL1("a,i") += - F_mi("i,m") * L1("a,m");
    tildeL1("a,i") += L1("e,m") * W_mbej("i,e,a,m");
    tildeL1("a,i") += 0.5 * L2("e,f,i,m") * W_abei("e,f,a,m");
    tildeL1("a,i") += - 0.5 * L2("a,e,m,n") * W_mbij("i,e,m,n");
    tildeL1("a,i") += - G_ae("e,f") * W_amef("e,i,f,a");
    tildeL1("a,i") += - G_mi("m,n") * W_mnie("m,i,n,a");

    if (subtractDiagonal) {
      tildeL1("a,i") -= fockMatrix_ta["vv_diag"]("e,a") * L1("e,i");
      tildeL1("a,i") += fockMatrix_ta["oo_diag"]("i,m") * L1("a,m");
    }

  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formL2_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                                          TArray &tildeL2, bool subtractDiagonal) const {

    if (subtractDiagonal) {
      tildeL2("a,b,i,j") = L2("a,e,i,j") * (F_ae("e,b") - fockMatrix_ta["vv_diag"]("e,b"));
      tildeL2("a,b,i,j") += - L2("b,e,i,j") * (F_ae("e,a") - fockMatrix_ta["vv_diag"]("e,a"));

      tildeL2("a,b,i,j") += - L2("a,b,i,m") * (F_mi("j,m") - fockMatrix_ta["oo_diag"]("j,m"));
      tildeL2("a,b,i,j") += L2("a,b,j,m") * (F_mi("i,m") - fockMatrix_ta["oo_diag"]("i,m"));

    } else {
      tildeL2("a,b,i,j") = L2("a,e,i,j") * F_ae("e,b");
      tildeL2("a,b,i,j") += - L2("b,e,i,j") * F_ae("e,a");

      tildeL2("a,b,i,j") += - L2("a,b,i,m") * F_mi("j,m");
      tildeL2("a,b,i,j") += L2("a,b,j,m") * F_mi("i,m");

    }

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

    tildeL2("a,b,i,j") += conj(antiSymMoints["vvoo"]("a,e,i,j")) * G_ae("b,e");
    tildeL2("a,b,i,j") += - conj(antiSymMoints["vvoo"]("b,e,i,j")) * G_ae("a,e");

    tildeL2("a,b,i,j") += - conj(antiSymMoints["vvoo"]("a,b,i,m")) * G_mi("m,j");
    tildeL2("a,b,i,j") += conj(antiSymMoints["vvoo"]("a,b,j,m")) * G_mi("m,i");

  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::runLambda() {

    auto lambda_start = tick();

    initilizeLambda();
    TArray &L1_ = Lg_.oneBody();
    TArray &L2_ = Lg_.twoBody();
    EOMCCSDVector<MatsT> L_old(vLabel_, oLabel_);

    std::shared_ptr<DIISTA<MatsT>> ldiis = nullptr;
    if(ccSettings_.useDIIS){
      ldiis = std::make_shared<DIISTA<MatsT>>(ccSettings_.nDIIS);
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

    for (auto iter = 0; iter < ccSettings_.maxiter; iter++){
      L_old = Lg_;

      updateG_ae(L_old.twoBody(), G_ae);
      updateG_mi(L_old.twoBody(), G_mi);

      formL1_tilde(L_old.oneBody(), L_old.twoBody(), G_ae, G_mi, L1_, true);
      L1_("a,i") += F_me("i,a");
      L1_("a,i") = L1_("a,i") * D_ai("a,i");

      formL2_tilde(L_old.oneBody(), L_old.twoBody(), G_ae, G_mi, L2_, true);
      L2_("a,b,i,j") += conj(antiSymMoints["vvoo"]("a,b,i,j"));
      L2_("a,b,i,j") = L2_("a,b,i,j") * D_abij("a,b,i,j");
      
      if(ccSettings_.useDIIS){

        // give solution vector to diis
        ldiis->WriteVector(Lg_);

        // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
        L_old.scale(-1.0);
        L_old.axpy(1.0, Lg_);

        //set error vector in DIIS
        ldiis->WriteErrorVector(L_old);

        // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
        ldiis->Extrapolate(Lg_);
      } else {
        L_old.axpy(-1.0, Lg_);
      }
      
      MatsT PE_old = pseudoEnergy;
      MatsT PEOneBody = fockMatrix_ta["vo"]("a,i").dot(L1_("a,i")).get();
      MatsT PETwoBodyL2 = antiSymMoints["vvoo"]("a,b,i,j").dot(L2_("a,b,i,j")).get();
      MatsT PETwoBodyL1 = antiSymMoints["vvoo"]("c,d,k,l").dot(L1_("c,k") * L1_("d,l")).get();
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

      if ( dPE < ccSettings_.eConv and dL < ccSettings_.tConv) {
        std::cout << "\n  Lambda iteration converged in " << iter << " steps." << std::endl;
        std::cout << "\n  Lambda Completed: Iteration total time "<< std::setw(10) << std::right
                  << std::setprecision(6) << tock(lambda_start) << " s" << std::endl;
        std::cout << bannerEnd << std::endl;
        break;
      }

      if (iter == ccSettings_.maxiter - 1){
        CErr(std::string("Lambda iterations didn't converge in ") + std::to_string(ccSettings_.maxiter) + " steps." );
      }
    }
  }

}; // namespace