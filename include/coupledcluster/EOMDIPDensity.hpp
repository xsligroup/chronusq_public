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
#include <coupledcluster.hpp>

namespace ChronusQ{

  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::initializeDensity() {
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
  void EOMDIP_3h1p<MatsT>::buildDensity(const TArray& t1, const TArray& t2, const TArray& r2, const TArray& r3, const TArray& l2, const TArray& l3) {


    TAManager &TAmanager = TAManager::get();
    MatsT scalar_0 = dot(l2("l,k"), r2("l,k"));
    TA::get_default_world().gop.fence();    
    MatsT scalar_1 = dot(l3("a,m,l,k"), r3("a,m,l,k"));
    TA::get_default_world().gop.fence();    

    // tempOps["oo_0"] += 0.500000 r3("I,a,l,k,i") l3("I,a,l,k,j") 
    // flops: o4v1L1: 1, o2v0: 1 | mem: o2v0: 2, 
    tempOps["oo_0"] = TAmanager.malloc<MatsT>("oo");
    tempOps["oo_0"]("i,j") = 0.500000 * r3("a,l,k,i") * l3("a,l,k,j");

    // Rho_ij += -0.500000 l3(a,l,k,j) r3(a,l,k,i) 
    // flops: o2v0: 1 | mem: o2v0: 1, 
    Rho_ij("i,j") -= tempOps["oo_0"]("i,j");

    // Rho_ai += -0.500000 l3(b,l,k,j) r3(b,l,k,i) t1(a,j) 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    Rho_ai("a,i") -= tempOps["oo_0"]("i,j") * t1("a,j");
    TAmanager.free("oo", std::move(tempOps["oo_0"]));

    // tempOps["vo_1"] += 0.500000 r2("I,l,k") l3("I,a,l,k,j") 
    // flops: o3v1L1: 1, o1v1: 1 | mem: o1v1: 2, 
    tempOps["vo_1"] = TAmanager.malloc<MatsT>("vo");
    tempOps["vo_1"]("a,j") = 0.500000 * r2("l,k") * l3("a,l,k,j");

    // Rho_ij += -0.500000 l3(a,l,k,j) r2(l,k) t1(a,i) 
    // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2, 
    Rho_ij("i,j") -= tempOps["vo_1"]("a,j") * t1("a,i");

    // Rho_ia += 0.500000 l3(a,k,j,i) r2(k,j) 
    // flops: o1v1: 1 | mem: o1v1: 1, 
    Rho_ia("i,a") += tempOps["vo_1"]("a,i");

    // Rho_ai += -0.500000 l3(b,l,k,j) r2(l,k) t2(a,b,j,i) 
    // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2, 
    Rho_ai("a,i") -= tempOps["vo_1"]("b,j") * t2("a,b,j,i");

    // Rho_ai += -0.500000 l3(b,l,k,j) r2(l,k) t1(a,j) t1(b,i) 
    // flops: o2v1: 2, o1v1: 1 | mem: o1v1: 2, o2v0: 1, 
    Rho_ai("a,i") -= t1("b,i") * tempOps["vo_1"]("b,j") * t1("a,j");

    // Rho_ab += 0.500000 l3(a,k,j,i) r2(k,j) t1(b,i) 
    // flops: o1v2: 1, o0v2: 1 | mem: o0v2: 2, 
    Rho_ab("a,b") += tempOps["vo_1"]("a,i") * t1("b,i");
    TAmanager.free("vo", std::move(tempOps["vo_1"]));

    // tempOps["oo_2"] += 1.000000 r2("I,k,i") l2("I,k,j") 
    // flops: o3v0L1: 1, o2v0: 1 | mem: o2v0: 2, 
    tempOps["oo_2"] = TAmanager.malloc<MatsT>("oo");
    tempOps["oo_2"]("i,j") = r2("k,i") * l2("k,j");

    // Rho_ij += -1.000000 l2(k,j) r2(k,i) 
    // flops: o2v0: 1 | mem: o2v0: 1, 
    Rho_ij("i,j") -= tempOps["oo_2"]("i,j");

    // Rho_ai += -1.000000 l2(k,j) r2(k,i) t1(a,j) 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    Rho_ai("a,i") -= tempOps["oo_2"]("i,j") * t1("a,j");
    TAmanager.free("oo", std::move(tempOps["oo_2"]));

    // Rho_ij += 0.500000 d(i,j) l2(l,k) r2(l,k) 
    // flops: o2v0: 1, o0v0: 1 | mem: o2v0: 1, o0v0: 1, 
    Rho_ij("i,j") += 0.500000 * scalar_0 * Id_oo("i,j");

    // Rho_ij += 0.166667 d(i,j) l3(a,m,l,k) r3(a,m,l,k) 
    // flops: o2v0: 1, o0v0: 1 | mem: o2v0: 1, o0v0: 1, 
    Rho_ij("i,j") += 0.166667 * scalar_1 * Id_oo("i,j");

    // Rho_ai += 0.500000 l2(k,j) r3(a,k,j,i) 
    // flops: o3v1L1: 1, o1v1: 1 | mem: o1v1: 2, 
    Rho_ai("a,i") += 0.500000 * l2("k,j") * r3("a,k,j,i");

    // Rho_ai += 0.166667 l3(b,l,k,j) r3(b,l,k,j) t1(a,i) 
    // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1, 
    Rho_ai("a,i") += 0.166667 * scalar_1 * t1("a,i");

    // Rho_ai += 0.500000 l2(k,j) r2(k,j) t1(a,i) 
    // flops: o1v1: 1, o0v0: 1 | mem: o1v1: 1, o0v0: 1, 
    Rho_ai("a,i") += 0.500000 * scalar_0 * t1("a,i");

    // Rho_ai += -0.166667 l3(b,l,k,j) r3(a,l,k,j) t1(b,i) 
    // flops: o4v1L1: 2, o1v1: 1 | mem: o4v0L1: 1, o1v1: 2, 
    Rho_ai("a,i") -= 0.166667 * l3("b,l,k,j") * t1("b,i") * r3("a,l,k,j");

    // Rho_ai += -0.500000 l3(b,l,j,k) r2(l,i) t2(a,b,j,k) 
    // flops: o4v1L1: 1, o3v2: 1, o1v1: 1 | mem: o3v1: 1, o1v1: 2, 
    Rho_ai("a,i") -= 0.500000 * l3("b,l,j,k") * r2("l,i") * t2("a,b,j,k");

    // Rho_ab += 0.166667 l3(a,k,j,i) r3(b,k,j,i) 
    // flops: o3v2L1: 1, o0v2: 1 | mem: o0v2: 2, 
    Rho_ab("a,b") += 0.166667 * l3("a,k,j,i") * r3("b,k,j,i");

  }

  template <typename MatsT>
  dcomplex EOMDIP_3h1p<MatsT>::calcOscillatorStrength(size_t i){
    std::array<MatsT, 3> mu_g2x, mu_x2g;
    
    MBExpansionSet<MatsT> &Reom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->R_);
    MBExpansionSet<MatsT> &Leom = *std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(this->L_);

    TArray &L2_ = this->Lg_->get_tensor("OneBody");
    TArray &L3_ = this->Lg_->get_tensor("TwoBody");

    buildDensity(T1_, T2_,
                 Reom.get(i).get_tensor("OneBody"),
                 Reom.get(i).get_tensor("TwoBody"),
                 L2_, L3_);
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
    TArray Rg1 = TAmanager.malloc<MatsT>("oo");
    Rg1("i,j") = 0.0 * Rg1("i,j");
    TArray Rg2 = TAmanager.malloc<MatsT>("vooo");
    Rg2("a,i,j,k") = 0.0 * Rg2("a,i,j,k");
    buildDensity(T1_, T2_, Rg1, Rg2,
                 Leom.get(i).get_tensor("OneBody"),
                 Leom.get(i).get_tensor("TwoBody"));
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
    TAmanager.free("vvoo", std::move(Rg2));


    MatsT DS = mu_g2x[0] * mu_x2g[0] + mu_g2x[1] * mu_x2g[1] + mu_g2x[2] * mu_x2g[2];

    dcomplex f = 2./3 * this->theta[i] * DS;
    return f;

  }


}; // namespace
