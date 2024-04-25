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

#include <physcon.hpp>
#include <coupledcluster.hpp>

namespace ChronusQ{

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::initilizeDensity() {
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

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::buildDensity(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){

      formRho_ij(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
      formRho_ab(t1, t2, r0, r1, r2, l0, l1, l2);
      formRho_ia(t1, t2, r0, r1, r2, l0, l1, l2, isSame);
      formRho_ai(t1, t2, r0, r1, r2, l0, l1, l2);
  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formRho_ij(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame) {


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

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formRho_ab(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2) {

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

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formRho_ai(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2){
    this->Rho_ai("a,i") = r0 * l1("a,i");
    this->Rho_ai("a,i") += l2("a,e,i,m") * r1("e,m");
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formRho_ia(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame){

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

  template <typename MatsT, typename IntsT>
  MatsT EOMCCSD<MatsT,IntsT>::calcOscillatorStrength(size_t i){
    std::array<MatsT, 3> mu_g2x, mu_x2g;
    
    EOMCCSDVectorSet<MatsT> &Reom = *std::dynamic_pointer_cast<EOMCCSDVectorSet<MatsT>>(R_);
    EOMCCSDVectorSet<MatsT> &Leom = *std::dynamic_pointer_cast<EOMCCSDVectorSet<MatsT>>(L_);

    MatsT L0_ = Lg_.zeroBody();
    TArray &L1_ = Lg_.oneBody();
    TArray &L2_ = Lg_.twoBody();

    buildDensity(T1_, T2_,
                 Reom.get(i).zeroBody(),
                 Reom.get(i).oneBody(),
                 Reom.get(i).twoBody(),
                 L0_, L1_, L2_);
    for (size_t j = 0; j < 3; j++) {
      mu_g2x[j]  = dot(muMatrix[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij("i,j")).get();
      mu_g2x[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia("i,a")).get();
      mu_g2x[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai("a,i")).get();
      mu_g2x[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab("a,b")).get();
    }
    TA::get_default_world().gop.fence();

    TAManager &TAmanager = TAManager::get();
    TArray Rg1 = TAmanager.malloc<MatsT>("vo");
    Rg1("a,i") = 0.0 * Rg1("a,i");
    TArray Rg2 = TAmanager.malloc<MatsT>("vvoo");
    Rg2("a,b,i,j") = 0.0 * Rg2("a,b,i,j");
    buildDensity(T1_, T2_, 1.0, Rg1, Rg2,
                 Leom.get(i).zeroBody(),
                 Leom.get(i).oneBody(),
                 Leom.get(i).twoBody());
    for (size_t j = 0; j < 3; j++) {
      mu_x2g[j]  = dot(muMatrix[static_cast<char>('X' + j) + std::string("oo")]("i,j"), Rho_ij("i,j")).get();
      mu_x2g[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("ov")]("i,a"), Rho_ia("i,a")).get();
      mu_x2g[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("vo")]("a,i"), Rho_ai("a,i")).get();
      mu_x2g[j] += dot(muMatrix[static_cast<char>('X' + j) + std::string("vv")]("a,b"), Rho_ab("a,b")).get();
    }
    TA::get_default_world().gop.fence();
    TAmanager.free("vo", std::move(Rg1));
    TAmanager.free("vvoo", std::move(Rg2));


    MatsT DS = mu_g2x[0] * mu_x2g[0] + mu_g2x[1] * mu_x2g[1] + mu_g2x[2] * mu_x2g[2];

    MatsT f = 2./3 * theta[i] * DS;
    return f;

  }

}; // namespace