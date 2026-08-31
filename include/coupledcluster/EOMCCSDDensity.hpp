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

#include <physcon.hpp>
#include <coupledcluster/EOMCCSD.hpp>

namespace ChronusQ{

  template <typename MatsT>
  void EOMCCSD<MatsT>::initializeDensity() {
    TAManager &TAmanager = TAManager::get();

    if (not Rho_ij.is_initialized()){
      Rho_ij = TAmanager.malloc<MatsT>("oo");
    }

    if (not Rho_ab.is_initialized()){
      Rho_ab = TAmanager.malloc<MatsT>("vv");
    }

    if (not Rho_ia.is_initialized()){
      Rho_ia = TAmanager.malloc<MatsT>("ov");
    }

    if (not Rho_ai.is_initialized()){
      Rho_ai = TAmanager.malloc<MatsT>("vo");
    }
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::buildDensity(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){
      formRho_ij(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
      formRho_ab(t1, t2, r0, r1, r2, l0, l1, l2);
      formRho_ia(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
      formRho_ai(t1, t2, r0, r1, r2, l0, l1, l2);
  }


  template <typename MatsT>
  void EOMCCSD<MatsT>::formRho_ij(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame) {


    if (isSame){
      TA::foreach_inplace(Rho_ij, [](TA::Tensor<MatsT> &tile){

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

      this->Rho_ij("i,j") -= r0 * l1("e,i") * t1("e,j");

    }
    else{
      this->Rho_ij("i,j") = - r0 * l1("e,i") * t1("e,j");
    }


    this->Rho_ij("i,j") += - 0.5 * r0 * l2("f,e,i,m") * t2("f,e,j,m");
    this->Rho_ij("i,j") += - l1("e,i") * r1("e,j");
    this->Rho_ij("i,j") +=  - 0.5 * l2("f,e,i,m") * r2("f,e,j,m");
    TArray tmp = TAManager::get().malloc<MatsT>("ov");
    tmp("i,f") = l2("f,e,i,m") * r1("e,m");
    this->Rho_ij("i,j") += - tmp("i,f") * t1("f,j");
    this->Rho_ij("i,j") = this->Rho_ij("j,i");

    TAManager::get().free("ov", std::move(tmp));
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formRho_ab(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2) {

    this->Rho_ab("a,b") = r0 * l1("b,m") * t1("a,m");
    this->Rho_ab("a,b") += 0.5 * r0 * l2("e,b,m,n") * t2("e,a,m,n");
    this->Rho_ab("a,b") += l1("b,m") * r1("a,m");
    this->Rho_ab("a,b") += 0.5 * l2("e,b,m,n") * r2("e,a,m,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ov");
    tmp("n,b") = l2("e,b,m,n") * r1("e,m");
    this->Rho_ab("a,b") += tmp("n,b") * t1("a,n");
    this->Rho_ab("a,b") = this->Rho_ab("b,a");

    TAManager::get().free("ov", std::move(tmp));
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formRho_ai(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2){
    this->Rho_ai("a,i") = r0 * l1("a,i");
    this->Rho_ai("a,i") += l2("a,e,i,m") * r1("e,m");
  }

  template <typename MatsT>
  void EOMCCSD<MatsT>::formRho_ia(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){

    if(isSame){
      this->Rho_ia("i,a") = t1("a,i");
      this->Rho_ia("i,a") += r0 * l1("e,m") * t2("a,e,i,m");
    }
    else{
      this->Rho_ia("i,a") = r0 * l1("e,m") * t2("a,e,i,m");
    }

    TAManager &TAmanager = TAManager::get();

    this->Rho_ia("i,a") += l0 * r1("a,i");
    TArray tmp1 = TAmanager.malloc<MatsT>("oo");
    tmp1("m,i") = l1("e,m") * t1("e,i");
    this->Rho_ia("i,a") += - r0 * tmp1("m,i") * t1("a,m");

    this->Rho_ia("i,a") += - tmp1("m,i") * r1("a,m");



    tmp1("m,i") = l2("e,f,m,n") * t2("e,f,i,n");
    this->Rho_ia("i,a") += - 0.5 * r0 * tmp1("m,i") * t1("a,m");

    this->Rho_ia("i,a") += - 0.5 * tmp1("m,i") * r1("a,m");

    TArray tmp2 = TAmanager.malloc<MatsT>("vv");
    tmp2("a,e") = l2("e,f,m,n") * t2("a,f,m,n");
    this->Rho_ia("i,a") += - 0.5 * r0 * tmp2("a,e") * t1("e,i");

    this->Rho_ia("i,a") += l1("e,m") * r2("a,e,i,m");


    tmp1("m,i") = l1("e,m") * r1("e,i");
    this->Rho_ia("i,a") += - tmp1("m,i") * t1("a,m");


    this->Rho_ia("i,a") += - 0.5 * tmp2("a,e") * r1("e,i");


    tmp2("a,e") = l2("e,f,m,n") * r2("a,f,m,n");
    this->Rho_ia("i,a") += - 0.5 * tmp2("a,e") * t1("e,i");


    tmp1("m,i") = l2("e,f,m,n") * r2("e,f,i,n");
    this->Rho_ia("i,a") += - 0.5 * tmp1("m,i") * t1("a,m");

    TArray tmp3 = TAmanager.malloc<MatsT>("ov");
    tmp3("m,e") = l2("e,f,m,n") * r1("f,n");


    tmp1("m,i") = tmp3("m,e") * t1("e,i");
    this->Rho_ia("i,a") += - tmp1("m,i") * t1("a,m");

    this->Rho_ia("i,a") += tmp3("m,e") * t2("a,e,i,m");

    TAmanager.free("oo", std::move(tmp1));
    TAmanager.free("vv", std::move(tmp2));
    TAmanager.free("ov", std::move(tmp3));

  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMCCSD<MatsT>::calcTransitionDipole(
    const MatsT r0, const TArray& r1, const TArray& r2,
    const MatsT l0,const TArray& l1, const TArray& l2, bool isSame) {

    std::array<MatsT, 3> mu;

    buildDensity(T1_, T2_, r0, r1, r2, l0, l1, l2, isSame);
    for (size_t j = 0; j < 3; j++) {
      mu[j]  = dot(this->muMatrix[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij("i,j")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia("i,a")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai("a,i")).get();
      TA::get_default_world().gop.fence();
      mu[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab("a,b")).get();
      TA::get_default_world().gop.fence();
    }
    TA::get_default_world().gop.fence();

    if (isSame) // Add frozen core contribution to the transition dipole moment if the left and right states are the same
      for (size_t j = 0; j < 3; j++) {
        mu[j] += this->intermediates_.Mu_fzc[j];
      }

    return mu;

  }

  template <typename MatsT>
  std::array<MatsT, 3> EOMCCSD<MatsT>::calcGround2ExcitedTransitionDipole(size_t i) {

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
  std::array<MatsT, 3> EOMCCSD<MatsT>::calcExcited2GroundTransitionDipole(size_t i) {

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
  std::array<MatsT, 3> EOMCCSD<MatsT>::calcExcited2ExcitedTransitionDipole(size_t i, size_t j) {

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
  std::array<MatsT, 3> EOMCCSD<MatsT>::calcGroundDipole() {

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
  dcomplex EOMCCSD<MatsT>::calcOscillatorStrength(size_t i){
    std::array<MatsT, 3> mu_g2x, mu_x2g;

    mu_g2x = calcGround2ExcitedTransitionDipole(i);
    mu_x2g = calcExcited2GroundTransitionDipole(i);

    MatsT DS = mu_g2x[0] * mu_x2g[0] + mu_g2x[1] * mu_x2g[1] + mu_g2x[2] * mu_x2g[2];

    dcomplex f = 2./3 * this->theta[i] * DS;
    return f;

  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::initializeDensity() {
    TAManager &TAmanager = TAManager::get();

    if (not Rho_ij.is_initialized()){
      Rho_ij = TAmanager.malloc<MatsT>("oo");
    }

    if (not Rho_ab.is_initialized()){
      Rho_ab = TAmanager.malloc<MatsT>("vv");
    }

    if (not Rho_ia.is_initialized()){
      Rho_ia = TAmanager.malloc<MatsT>("ov");
    }

    if (not Rho_ai.is_initialized()){
      Rho_ai = TAmanager.malloc<MatsT>("vo");
    }
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::buildDensity(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val, bool isSame){

      formRho_ij(t1, t2, r0, r1, r2, r2val, l0, l1, l2, l2val, isSame);
      formRho_ab(t1, t2, r0, r1, r2, r2val, l0, l1, l2, l2val);
      formRho_ia(t1, t2, r0, r1, r2, r2val, l0, l1, l2, l2val, isSame);
      formRho_ai(t1, t2, r0, r1, r2, r2val, l0, l1, l2, l2val);
  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formRho_ai(const TArray& t1_ai, const TArray& t2_abij, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2Val) {
    this->Rho_ai("a,i").block(b.rc)  = r0 * l1("a,i");
    this->Rho_ai("a,i").block(b.rc) += l2("a,e,i,m") * r1("e,m");
    this->Rho_ai("a,i").block(b.rh)  = -1.0 * l2Val("a,e,m,i") * r1("e,m"); // factor of -1 for swapping m and i
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formRho_ab(const TArray& t1_ai, const TArray& t2_abij, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2Val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2Val) {
    TAManager &TAmanager = TAManager::get();

    this->Rho_ab("a,b").block(b.rr) = r0 * l1("b,m") * t1_ai("a,m").block(b.rc);
    this->Rho_ab("a,b").block(b.rr) += 0.5 * r0 * l2("e,b,m,n") * t2_abij("e,a,m,n").block(b.rrcc);
    this->Rho_ab("a,b").block(b.rr) += r0 * l2Val("e,b,m,n") * t2_abij("e,a,m,n").block(b.rrch); // factor of two applied for missing valence contribution
    this->Rho_ab("a,b").block(b.rr) += l1("b,m") * r1("a,m");
    this->Rho_ab("a,b").block(b.rr) += 0.5 * l2("e,b,m,n") * r2("e,a,m,n");
    this->Rho_ab("a,b").block(b.rr) += l2Val("e,b,m,n") * r2Val("e,a,m,n"); // factor of two applied for missing valence contribution
    TArray tmp = TAmanager.malloc<dcomplex>("cr");
    tmp("n,b") = l2("e,b,m,n") * r1("e,m");
    this->Rho_ab("a,b").block(b.rr) += tmp("n,b") * t1_ai("a,n").block(b.rc);
    TArray tmpVal = TAmanager.malloc<dcomplex>("hr");
    tmpVal("n,b") = l2Val("e,b,m,n") * r1("e,m");
    this->Rho_ab("a,b").block(b.rr) += tmpVal("n,b") * t1_ai("a,n").block(b.rh);
    this->Rho_ab("a,b").block(b.rr) = this->Rho_ab("b,a").block(b.rr);

    TAManager::get().free("cr", std::move(tmp));
    TAManager::get().free("hr", std::move(tmpVal));
  }


  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formRho_ij(const TArray& t1_ai, const TArray& t2_abij, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2Val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2Val, bool isSame) {
    TAManager &TAmanager = TAManager::get();

    if (isSame){
      TA::foreach_inplace(Rho_ij, [](TA::Tensor<MatsT> &tile){

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

      this->Rho_ij("i,j").block(b.cc) -= r0 * l1("e,i") * t1_ai("e,j").block(b.rc);
      this->Rho_ij("i,j").block(b.ch) -= r0 * l1("e,i") * t1_ai("e,j").block(b.rh);

    }
    else{
      this->Rho_ij("i,j").block(b.cc) = - r0 * l1("e,i") * t1_ai("e,j").block(b.rc);
      this->Rho_ij("i,j").block(b.ch) = - r0 * l1("e,i") * t1_ai("e,j").block(b.rh);
    }


    this->Rho_ij("i,j").block(b.cc) += - 0.5 * r0 * l2("f,e,i,m") * t2_abij("f,e,j,m").block(b.rrcc);
    this->Rho_ij("i,j").block(b.cc) += - 0.5 * r0 * l2Val("f,e,i,m") * t2_abij("f,e,j,m").block(b.rrch);
    this->Rho_ij("i,j").block(b.ch) +=   0.5 * r0 * l2("f,e,i,m") * t2_abij("f,e,m,j").block(b.rrch); // factor of -1 due to swapping m and j 
    this->Rho_ij("i,j").block(b.ch) += - 0.5 * r0 * l2Val("f,e,i,m") * t2_abij("f,e,j,m").block(b.rrhh);
    this->Rho_ij("i,j").block(b.hc)  =   0.5 * r0 * l2Val("f,e,m,i") * t2_abij("f,e,j,m").block(b.rrcc); // factor of -1 due to swapping m and j
    this->Rho_ij("i,j").block(b.hh)  = - 0.5 * r0 * l2Val("f,e,m,i") * t2_abij("f,e,m,j").block(b.rrch); // factor of -1 due to swapping m and i and swapping j and m

    this->Rho_ij("i,j").block(b.cc) += - l1("e,i") * r1("e,j");

    this->Rho_ij("i,j").block(b.cc) += - 0.5 * l2("f,e,i,m") * r2("f,e,j,m");
    this->Rho_ij("i,j").block(b.cc) += - 0.5 * l2Val("f,e,i,m") * r2Val("f,e,j,m");
    this->Rho_ij("i,j").block(b.ch) +=   0.5 * l2("f,e,i,m") * r2Val("f,e,m,j"); // factor of -1 due to swapping m and j
    this->Rho_ij("i,j").block(b.hc) +=   0.5 * l2Val("f,e,m,i") * r2("f,e,j,m"); // factor of -1 due to swapping m and i 
    this->Rho_ij("i,j").block(b.hh) += - 0.5 * l2Val("f,e,m,i") * r2Val("f,e,m,j"); // factor of -1 due to swapping m and j and swapping m and i

    TArray tmp = TAmanager.malloc<dcomplex>("cr");
    tmp("i,f") = l2("f,e,i,m") * r1("e,m");
    TArray tmpVal = TAmanager.malloc<dcomplex>("hr");
    tmpVal("i,f") = -1.0 * l2Val("f,e,m,i") * r1("e,m"); // factor of -1 due to swapping m and i
    this->Rho_ij("i,j").block(b.cc) += - tmp("i,f") * t1_ai("f,j").block(b.rc);
    this->Rho_ij("i,j").block(b.ch) += - tmp("i,f") * t1_ai("f,j").block(b.rh);
    this->Rho_ij("i,j").block(b.hc) += - tmpVal("i,f") * t1_ai("f,j").block(b.rc);
    this->Rho_ij("i,j").block(b.hh) += - tmpVal("i,f") * t1_ai("f,j").block(b.rh);

    //this->Rho_ij("i,j").block(b.cc) = this->Rho_ij("j,i").block(b.cc);
    //this->Rho_ij("i,j").block(b.hh) = this->Rho_ij("j,i").block(b.hh);

    TAManager::get().free("cr", std::move(tmp));
    TAManager::get().free("hr", std::move(tmpVal));
  }

  template <typename MatsT>
  void CVSEOMCCSD<MatsT>::formRho_ia(const TArray& t1_ai, const TArray& t2_abij, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2Val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2Val, bool isSame) {
    TAManager &TAmanager = TAManager::get();

    if(isSame){
      this->Rho_ia("i,a").block(b.cr) = t1_ai("a,i").block(b.rc);
      this->Rho_ia("i,a").block(b.hr) = t1_ai("a,i").block(b.rh);
      this->Rho_ia("i,a").block(b.cr) += r0 * l1("e,m") * t2_abij("a,e,i,m").block(b.rrcc);
      this->Rho_ia("i,a").block(b.hr) -= r0 * l1("e,m") * t2_abij("a,e,m,i").block(b.rrch); // factor of -1 due to swapping m and i 
    }
    else{
      this->Rho_ia("i,a").block(b.cr) =        r0 * l1("e,m") * t2_abij("a,e,i,m").block(b.rrcc);
      this->Rho_ia("i,a").block(b.hr) = -1.0 * r0 * l1("e,m") * t2_abij("a,e,m,i").block(b.rrch); // factor of -1 due to swapping m and i
    }

    this->Rho_ia("i,a").block(b.cr) += l0 * r1("a,i");

    TArray tmp1cc = TAmanager.malloc<dcomplex>("cc");
    tmp1cc("m,i") = l1("e,m") * t1_ai("e,i").block(b.rc);
    TArray tmp1cV = TAmanager.malloc<dcomplex>("ch");
    tmp1cV("m,i") = l1("e,m") * t1_ai("e,i").block(b.rh);

    this->Rho_ia("i,a").block(b.cr) += - r0 * tmp1cc("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.hr) += - r0 * tmp1cV("m,i") * t1_ai("a,m").block(b.rc);

    this->Rho_ia("i,a").block(b.cr) += - tmp1cc("m,i") * r1("a,m");
    this->Rho_ia("i,a").block(b.hr) += - tmp1cV("m,i") * r1("a,m");

    tmp1cc("m,i")  = l2("e,f,m,n") * t2_abij("e,f,i,n").block(b.rrcc);
    tmp1cc("m,i") += l2Val("e,f,m,n") * t2_abij("e,f,i,n").block(b.rrch);
    tmp1cV("m,i")  = -1.0 * l2("e,f,m,n") * t2_abij("e,f,n,i").block(b.rrch); // factor of -1 applied for swapping i and n
    tmp1cV("m,i") += l2Val("e,f,m,n") * t2_abij("e,f,i,n").block(b.rrhh);
    TArray tmp1Vc = TAmanager.malloc<dcomplex>("hc");
    tmp1Vc("m,i") = -1.0 * l2Val("e,f,n,m") * t2_abij("e,f,i,n").block(b.rrcc); // factor of -1 applied for swapping n and m 
    TArray tmp1VV = TAmanager.malloc<dcomplex>("hh");
    tmp1VV("m,i") = l2Val("e,f,n,m") * t2_abij("e,f,n,i").block(b.rrch); // factor of -1 applied for swapping n and m and swapping n and i

    this->Rho_ia("i,a").block(b.cr) += - 0.5 * r0 * tmp1cc("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.cr) += - 0.5 * r0 * tmp1Vc("m,i") * t1_ai("a,m").block(b.rh); // factor of -1 applied for swapping m and i
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * r0 * tmp1cV("m,i") * t1_ai("a,m").block(b.rc); 
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * r0 * tmp1VV("m,i") * t1_ai("a,m").block(b.rh); 

    this->Rho_ia("i,a").block(b.cr) += - 0.5 * tmp1cc("m,i") * r1("a,m");
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * tmp1cV("m,i") * r1("a,m");

    TArray tmp2 = TAmanager.malloc<dcomplex>("rr");
    tmp2("a,e")  = l2("e,f,m,n") * t2_abij("a,f,m,n").block(b.rrcc);
    tmp2("a,e") += 2.0 * l2Val("e,f,m,n") * t2_abij("a,f,m,n").block(b.rrch); // factor of two applied for missing ualence contribution
    this->Rho_ia("i,a").block(b.cr) += - 0.5 * r0 * tmp2("a,e") * t1_ai("e,i").block(b.rc);
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * r0 * tmp2("a,e") * t1_ai("e,i").block(b.rh);

    this->Rho_ia("i,a").block(b.cr) += l1("e,m") * r2("a,e,i,m");
    this->Rho_ia("i,a").block(b.hr) -= l1("e,m") * r2Val("a,e,m,i"); // factor of -1 applied for swapping i and m 

    tmp1cc("m,i") = l1("e,m") * r1("e,i");
    this->Rho_ia("i,a").block(b.cr) += - tmp1cc("m,i") * t1_ai("a,m").block(b.rc);

    this->Rho_ia("i,a").block(b.cr) += - 0.5 * tmp2("a,e") * r1("e,i");

    tmp2("a,e")  = l2("e,f,m,n") * r2("a,f,m,n");
    tmp2("a,e") += 2.0 * l2Val("e,f,m,n") * r2Val("a,f,m,n"); // factor of two applied for missing ualence contribution
    this->Rho_ia("i,a").block(b.cr) += - 0.5 * tmp2("a,e") * t1_ai("e,i").block(b.rc);
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * tmp2("a,e") * t1_ai("e,i").block(b.rh);

    tmp1cc("m,i")  = l2("e,f,m,n") * r2("e,f,i,n");
    tmp1cc("m,i") += l2Val("e,f,m,n") * r2Val("e,f,i,n");
    tmp1cV("m,i")  = -1.0 * l2("e,f,m,n") * r2Val("e,f,n,i"); // factor of -1 for swapping i and n
    tmp1Vc("m,i")  = -1.0 * l2Val("e,f,n,m") * r2("e,f,i,n"); // factor of -1 for swapping m and n
    tmp1VV("m,i")  = l2Val("e,f,n,m") * r2Val("e,f,n,i"); // factor of -1 for swapping m and n and n and i 

    this->Rho_ia("i,a").block(b.cr) += - 0.5 * tmp1cc("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.cr) += - 0.5 * tmp1Vc("m,i") * t1_ai("a,m").block(b.rh); // factor of -1 for swapping m and i 
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * tmp1cV("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.hr) += - 0.5 * tmp1VV("m,i") * t1_ai("a,m").block(b.rh);

    TArray tmp3cu = TAmanager.malloc<dcomplex>("cr");
    TArray tmp3Vu = TAmanager.malloc<dcomplex>("hr");
    tmp3cu("m,e") = l2("e,f,m,n") * r1("f,n");
    tmp3Vu("m,e") = l2Val("e,f,n,m") * -1.0 * r1("f,n"); // factor of -1 for swapping m and n

    tmp1cc("m,i") = tmp3cu("m,e") * t1_ai("e,i").block(b.rc);
    tmp1cV("m,i") = tmp3cu("m,e") * t1_ai("e,i").block(b.rh);
    tmp1Vc("m,i") = tmp3Vu("m,e") * t1_ai("e,i").block(b.rc);
    tmp1VV("m,i") = tmp3Vu("m,e") * t1_ai("e,i").block(b.rh);

    this->Rho_ia("i,a").block(b.cr) += - tmp1cc("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.cr) += - tmp1Vc("m,i") * t1_ai("a,m").block(b.rh);
    this->Rho_ia("i,a").block(b.hr) += - tmp1cV("m,i") * t1_ai("a,m").block(b.rc);
    this->Rho_ia("i,a").block(b.hr) += - tmp1VV("m,i") * t1_ai("a,m").block(b.rh);

    this->Rho_ia("i,a").block(b.cr) += tmp3cu("m,e") * t2_abij("a,e,i,m").block(b.rrcc);
    this->Rho_ia("i,a").block(b.cr) += tmp3Vu("m,e") * t2_abij("a,e,i,m").block(b.rrch);
    this->Rho_ia("i,a").block(b.hr) -= tmp3cu("m,e") * t2_abij("a,e,m,i").block(b.rrch); // factor of -1 for swapping m and i 
    this->Rho_ia("i,a").block(b.hr) += tmp3Vu("m,e") * t2_abij("a,e,i,m").block(b.rrhh);

    TAManager::get().free("cc", std::move(tmp1cc));
    TAManager::get().free("ch", std::move(tmp1cV));
    TAManager::get().free("hc", std::move(tmp1Vc));
    TAManager::get().free("hh", std::move(tmp1VV));

    TAManager::get().free("rr", std::move(tmp2));

    TAManager::get().free("cr", std::move(tmp3cu));
    TAManager::get().free("hr", std::move(tmp3Vu));
  }

  template <typename MatsT>
  dcomplex CVSEOMCCSD<MatsT>::calcOscillatorStrength(size_t i){
    std::array<MatsT, 3> mu_g2x, mu_x2g;


    MBExpansionSet<MatsT> &Reom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);
    MBExpansionSet<MatsT> &Leom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->L_);

    MatsT L0_ = this->Lg_->zeroBody();
    TArray &L1_ = this->Lg_->get_tensor("OneBody");
    TArray &L2_ = this->Lg_->get_tensor("TwoBody_core");
    TArray &L2Val_ = this->Lg_->get_tensor("TwoBody_val");

    buildDensity(T1_, T2_,
                 Reom.get(i).zeroBody(),
                 Reom.get(i).get_tensor("OneBody"),
                 Reom.get(i).get_tensor("TwoBody_core"),
                 Reom.get(i).get_tensor("TwoBody_val"),
                 L0_, L1_, L2_, L2Val_);
    for (size_t j = 0; j < 3; j++) {
      mu_g2x[j]  = dot(this->muMatrix[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij("i,j")).get();
      TA::get_default_world().gop.fence();    
      mu_g2x[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia("i,a")).get();
      TA::get_default_world().gop.fence();    
      mu_g2x[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai("a,i")).get();
      TA::get_default_world().gop.fence();    
      mu_g2x[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab("a,b")).get();
      TA::get_default_world().gop.fence();   

    }

    TAManager &TAmanager = TAManager::get();
    TArray Rg1 = TAmanager.malloc<MatsT>("vc");
    Rg1("a,i") = 0.0 * Rg1("a,i");
    TArray Rg2 = TAmanager.malloc<MatsT>("vvcc");
    Rg2("a,b,i,j") = 0.0 * Rg2("a,b,i,j");
    TArray Rg2Val = TAmanager.malloc<MatsT>("vvch");
    Rg2Val("a,b,i,j") = 0.0 * Rg2Val("a,b,i,j");
    buildDensity(T1_, T2_, 1.0, Rg1, Rg2, Rg2Val,
                 Leom.get(i).zeroBody(),
                 Leom.get(i).get_tensor("OneBody"),
                 Leom.get(i).get_tensor("TwoBody_core"),
                 Leom.get(i).get_tensor("TwoBody_val"));
    for (size_t j = 0; j < 3; j++) {
      mu_x2g[j]  = dot(this->muMatrix[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij("i,j")).get();
      TA::get_default_world().gop.fence();    
      mu_x2g[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia("i,a")).get();
      TA::get_default_world().gop.fence();    
      mu_x2g[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai("a,i")).get();
      TA::get_default_world().gop.fence();    
      mu_x2g[j] += dot(this->muMatrix[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab("a,b")).get();
      TA::get_default_world().gop.fence();    
    }
    TAmanager.free("vo", std::move(Rg1));
    TAmanager.free("vvcc", std::move(Rg2));
    TAmanager.free("vvch", std::move(Rg2Val));


    MatsT DS = mu_g2x[0] * mu_x2g[0] + mu_g2x[1] * mu_x2g[1] + mu_g2x[2] * mu_x2g[2];

    dcomplex f = 2./3 * this->theta[i] * DS;
    return f;

  }

}; // namespace
