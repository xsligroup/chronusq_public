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
#include <coupledcluster/MBExpansion.hpp>
#include <itersolver/davidson.hpp>
#include <itersolver.hpp>
#include <coupledcluster/EOMDIP_3h1p_full_Hbar.hpp>
#include <coupledcluster/EOMDIP_3h1p_guess.hpp>
#include <coupledcluster/EOMDIPDensity.hpp>

namespace ChronusQ{

  template <typename MatsT>
  EOMDIP_3h1p<MatsT>::EOMDIP_3h1p(const SafeFile &savFile,
                                CCIntermediates<MatsT> &intermediates,
                                const EOMSettings &eomSettings,
                                const CoupledClusterSettings &ccSettings):
    EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
    vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
    tempPerm_oo(intermediates.tempPerm_oo),
    F_me(intermediates.F_me),
    F_ae(intermediates.F_ae),
    F_mi(intermediates.F_mi),
    tau(intermediates.tau),
    W_mnij(intermediates.W_mnij),
    W_mbej(intermediates.W_mbej),
    //tempPerm_vooo(intermediates.tempPerm_vooo),
    sigmaOps(intermediates.sigmaOps),
    tempOps(intermediates.tempOps),
    T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")) {

    TAManager &TAmanager = TAManager::get();

    // without L, we don't need D
    if (eomSettings.oscillator_strength == false) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    nV_ = TAmanager.getRange(vLabel_).extent();
    nO_ = TAmanager.getRange(oLabel_).extent();
    nO2shift_ = nO_ * (nO_ - 1) / 2;
    nO3shift_ = nO_ * (nO_ - 1) * (nO_ - 2) / 6;
   
    this->Hbar_dimension_offsets.emplace("TwoBody", nO2shift_);
    this->Hbar_dimension_offsets.emplace("ThreeBody", nO3shift_*nV_);
//#define DEBUG_DIP 
    this->Hbar_dim = nO2shift_ + nV_ * nO3shift_;
    this->outOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);
    
    ijIndices_.clear();
    ijIndices_.resize(nO_, std::vector<size_t>(nO_, outOfBound_));
    ijkIndices_.clear();
    ijkIndices_.resize(nO_, std::vector<std::vector<size_t>>(nO_, std::vector<size_t>(nO_, outOfBound_)));
    
    size_t idx = 0;
    for (size_t j = 0; j < nO_; j++) {
      for (size_t i = 0; i < std::min(j, nO_); i++) {
        ijIndices_[i][j] = idx;
        ijIndices_[j][i] = idx++;
      }
    }

    idx = 0;
    for (size_t k = 0; k < nO_; k++) {
      for (size_t j = 0; j < std::min(k, nO_); j++) {
        for (size_t i = 0; i < std::min(j, nO_); i++) {
          ijkIndices_[i][j][k] = idx;
          ijkIndices_[j][k][i] = idx;
          ijkIndices_[k][i][j] = idx;
          ijkIndices_[j][i][k] = idx;
          ijkIndices_[i][k][j] = idx;
          ijkIndices_[k][j][i] = idx++;
        }
      }
    }

    this->tensor_builder_.push_back(std::string(""));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->tensor_builder_.push_back(std::string({intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));
    }

//  template <typename MatsT>
//  void EOMDIP_3h1p<MatsT>::buildCCSDEnergy(){ 
//    // get CCSD ground state energy
//
//    // H += -0.500000 <j,i||a,b> t1(a,i) t1(b,j)
//    // flops: o0v0: 1 | mem: o0v0: 1,
//    ccsd_energy = -0.500000 * scalar_4;
//
//    // ccsd_energy += 0.250000 <j,i||a,b> t2(a,b,j,i)
//    // flops: o0v0: 1 | mem: o0v0: 1,
//    ccsd_energy += 0.250000 * scalar_3;
//
//    // ccsd_energy += -0.500000 <j,i||j,i>
//    // flops: o0v0: 1 | mem: o0v0: 1,
//    ccsd_energy -= 0.500000 * scalar_2;
//
//    // ccsd_energy += 1.000000 f(i,a) t1(a,i)
//    // flops: o0v0: 1 | mem: o0v0: 1,
//    ccsd_energy += scalar_1;
//
//    // ccsd_energy += 1.000000 f(i,i)
//    // flops: o0v0: 1 | mem: o0v0: 1,
//    ccsd_energy += scalar_0;
//
//    }





  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::initializeEOMCC(){ 
    TAManager &TAmanager = TAManager::get();
//    if (not F_ae.is_initialized()){
//      F_ae = TAmanager.malloc<MatsT>("vv");
//    }
//    if (not F_mi.is_initialized()){
//      F_mi = TAmanager.malloc<MatsT>("oo");
//    }
//    if (not W_mnij.is_initialized()){
//      W_mnij = TAmanager.malloc<MatsT>("vvoo");
//    }
//    if (not W_mbej.is_initialized()){
//      W_mbej = TAmanager.malloc<MatsT>("ovvo");
//    }
    if (not tempPerm_oo.is_initialized()){
      tempPerm_oo = TAmanager.malloc<MatsT>("oo");
    }
    //if (not tempPerm_vooo.is_initialized()){
    //  tempPerm_vooo = TAmanager.malloc<MatsT>("vooo");
    //}
    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      this->antiSymMoints["oooo"] = TAmanager.malloc<MatsT>("oooo");
      this->antiSymMoints["oooo"]("m,n,i,j")  = this->riMoints["boo"]("Q,m,i") * this->riMoints["boo"]("Q,n,j");
      this->antiSymMoints["oooo"]("m,n,i,j") -= this->antiSymMoints["oooo"]("m,n,j,i");
      this->antiSymMoints["vooo"] = TAmanager.malloc<MatsT>("vooo");
      this->antiSymMoints["vooo"]("a,n,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["boo"]("Q,n,j");
      this->antiSymMoints["vooo"]("a,n,i,j") -= this->antiSymMoints["vooo"]("a,n,j,i");
      this->antiSymMoints["vovo"] = TAmanager.malloc<MatsT>("vovo");
      this->antiSymMoints["vovo"]("a,n,e,j")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["boo"]("Q,n,j");
      this->antiSymMoints["vovo"]("a,n,e,j") -= this->riMoints["bvo"]("Q,a,j") * this->riMoints["bov"]("Q,n,e");
      this->antiSymMoints["vvoo"] = TAmanager.malloc<MatsT>("vvoo");
      this->antiSymMoints["vvoo"]("a,b,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvoo"]("a,b,i,j") -= this->antiSymMoints["vvoo"]("a,b,j,i");
    }
    if (not Id_oo.is_initialized()){
      Id_oo = TAmanager.malloc<MatsT>("oo");
    }
    //if (not Id_oooo.is_initialized()){
    //  Id_oooo = TAmanager.malloc<MatsT>("oooo");
    //}
    TA::foreach_inplace(Id_oo, [&](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0,};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          if(x[0]==x[1]) 
            tile[x] = 1.0;
          else 
            tile[x] = 0.0;
    });
    TA::get_default_world().gop.fence();    

    //Id_oooo("i,j,k,l") = Id_oo("i,k") * Id_oo("j,l");

    sigmaOps.emplace(std::make_pair("vooo_0", TAmanager.malloc<MatsT>("vooo")));
    sigmaOps.emplace(std::make_pair("vo_1", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("oooo_4", TAmanager.malloc<MatsT>("oooo")));
    sigmaOps.emplace(std::make_pair("vooo_5", TAmanager.malloc<MatsT>("vooo")));
    sigmaOps.emplace(std::make_pair("vvoo_7", TAmanager.malloc<MatsT>("vvoo")));
    sigmaOps.emplace(std::make_pair("vo_8", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vv_9", TAmanager.malloc<MatsT>("vv")));
    sigmaOps.emplace(std::make_pair("oo_10", TAmanager.malloc<MatsT>("oo")));
    sigmaOps.emplace(std::make_pair("vo_12", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_13", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vooo_14", TAmanager.malloc<MatsT>("vooo")));
    sigmaOps.emplace(std::make_pair("oooo_16", TAmanager.malloc<MatsT>("oooo")));
    sigmaOps.emplace(std::make_pair("vv_19", TAmanager.malloc<MatsT>("vv")));
    sigmaOps.emplace(std::make_pair("vo_20", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_21", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("oo_22", TAmanager.malloc<MatsT>("oo")));
    sigmaOps.emplace(std::make_pair("vo_23", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("oo_24", TAmanager.malloc<MatsT>("oo")));
    sigmaOps.emplace(std::make_pair("vo_25", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_15", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("oooo_17", TAmanager.malloc<MatsT>("oooo")));
    sigmaOps.emplace(std::make_pair("vo_31", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("oo_33", TAmanager.malloc<MatsT>("oo")));
    sigmaOps.emplace(std::make_pair("vooo_30", TAmanager.malloc<MatsT>("vooo")));
    sigmaOps.emplace(std::make_pair("vo_32", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_34", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_35", TAmanager.malloc<MatsT>("vo")));
    sigmaOps.emplace(std::make_pair("vo_36", TAmanager.malloc<MatsT>("vo")));

    TA::get_default_world().gop.fence();    
  }

  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::formEOMIntermediates() {

    //scalar_0 = dot(this->fockMatrix_ta["oo"]("o0,o1"), Id_oo("o0,o1")).get();
    //TA::get_default_world().gop.fence();    
    //scalar_1 = dot(this->fockMatrix_ta["ov"]("k,a"), this->T1_("a,k")).get();
    //TA::get_default_world().gop.fence();    
    //scalar_2 = dot(this->antiSymMoints["oooo"]("o0,o1,o2,o3"), Id_oooo("o0,o1,o2,o3")).get();
    //TA::get_default_world().gop.fence();    
    //scalar_3 = dot(conj(this->antiSymMoints["vvoo"]("a,b,l,k")), this->T2_("a,b,l,k")).get();
    //TA::get_default_world().gop.fence();    
    //scalar_4 = dot(conj(this->antiSymMoints["vvoo"]("a,b,l,k")) * this->T1_("a,k"), this->T1_("b,l")).get();
    //TA::get_default_world().gop.fence();    

    //buildCCSDEnergy(); 

    // sigmaOps["vooo_0"] += 1.000000 this->T1_("a,j") conj(this->antiSymMoints["vvoo"]("a,b,l,k")) 
    // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_0"]("b,j,l,k") = this->T1_("a,j") * conj(this->antiSymMoints["vvoo"]("a,b,l,k"));

    // sigmaOps["vo_1"] += 1.000000 conj(this->antiSymMoints["vvoo"]("a,b,l,k")) this->T1_("a,k") 
    // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_1"]("b,l") = conj(this->antiSymMoints["vvoo"]("a,b,l,k")) * this->T1_("a,k");

    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      this->antiSymMoints["vvvo"] = TAManager::get().malloc<MatsT>("vvvo");
      this->antiSymMoints["vvvo"]("a,b,e,j")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvvo"]("a,b,e,j") -= this->antiSymMoints["vvvo"]("b,a,e,j");
    };
    // sigmaOps["oooo_4"] += 1.000000 this->T2_("b,c,j,k") conj(this->antiSymMoints["vvoo"]("b,c,m,l")) 
    // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2, 
    sigmaOps["oooo_4"]("j,k,m,l") = this->T2_("b,c,j,k") * conj(this->antiSymMoints["vvoo"]("b,c,m,l"));

    // sigmaOps["vvoo_7"] += 1.000000 conj(this->antiSymMoints["vvvo"]("b,c,a,l")) this->T1_("b,k") 
    // flops: o2v3: 1, o2v2: 1 | mem: o2v2: 2, 
    sigmaOps["vvoo_7"]("a,c,l,k") = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,k");

    // sigmaOps["vo_8"] += 1.000000 this->T2_("b,c,k,l") conj(this->antiSymMoints["vvvo"]("b,c,a,l")) 
    // flops: o2v3: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_8"]("a,k") = this->T2_("b,c,k,l") * conj(this->antiSymMoints["vvvo"]("b,c,a,l"));

    // sigmaOps["vv_9"] += 1.000000 conj(this->antiSymMoints["vvoo"]("b,c,m,l")) this->T2_("b,a,m,l") 
    // flops: o2v3: 1, o0v2: 1 | mem: o0v2: 2, 
    sigmaOps["vv_9"]("c,a") = conj(this->antiSymMoints["vvoo"]("b,c,m,l")) * this->T2_("b,a,m,l");

    // sigmaOps["oo_10"] += 1.000000 this->T2_("a,b,j,k") conj(this->antiSymMoints["vvoo"]("a,b,l,k")) 
    // flops: o3v2: 1, o2v0: 1 | mem: o2v0: 2, 
    sigmaOps["oo_10"]("j,l") = this->T2_("a,b,j,k") * conj(this->antiSymMoints["vvoo"]("a,b,l,k"));

    // sigmaOps["vo_12"] += 1.000000 conj(this->antiSymMoints["vvoo"]("b,c,m,l")) this->T2_("b,c,k,m") this->T1_("a,l") 
    // flops: o3v2: 1, o2v1: 1, o1v1: 1 | mem: o1v1: 2, o2v0: 1, 
    sigmaOps["vo_12"]("a,k") = conj(this->antiSymMoints["vvoo"]("b,c,m,l")) * this->T2_("b,c,k,m") * this->T1_("a,l");

    // sigmaOps["vo_13"] += 1.000000 this->T2_("b,a,m,l") conj(this->antiSymMoints["vooo"]("b,k,m,l")) 
    // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_13"]("a,k") = this->T2_("b,a,m,l") * conj(this->antiSymMoints["vooo"]("b,k,m,l"));

    // sigmaOps["vooo_14"] += 1.000000 this->T2_("b,a,j,k") this->fockMatrix_ta["ov"]("l,b") 
    // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_14"]("a,j,k,l") = this->T2_("b,a,j,k") * this->fockMatrix_ta["ov"]("l,b");

    // sigmaOps["oooo_16"] += 1.000000 conj(this->antiSymMoints["vooo"]("a,j,l,k")) this->T1_("a,i") 
    // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2, 
    sigmaOps["oooo_16"]("l,k,j,i") = conj(this->antiSymMoints["vooo"]("a,j,l,k")) * this->T1_("a,i");

    // sigmaOps["vv_19"] += 1.000000 conj(this->antiSymMoints["vvvo"]("b,c,a,l")) this->T1_("b,l") 
    // flops: o1v3: 1, o0v2: 1 | mem: o0v2: 2, 
    sigmaOps["vv_19"]("a,c") = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,l");

    // sigmaOps["vo_20"] += 1.000000 this->antiSymMoints["vovo"]("a,l,b,k") this->T1_("b,l") 
    // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_20"]("a,k") = this->antiSymMoints["vovo"]("a,l,b,k") * this->T1_("b,l");

    // sigmaOps["vo_21"] += 1.000000 this->fockMatrix_ta["ov"]("l,b") this->T2_("b,a,k,l") 
    // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_21"]("a,k") = this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,k,l");

    // sigmaOps["oo_22"] += 1.000000 conj(this->antiSymMoints["vooo"]("a,j,l,k")) this->T1_("a,k") 
    // flops: o3v1: 1, o2v0: 1 | mem: o2v0: 2, 
    sigmaOps["oo_22"]("l,j") = conj(this->antiSymMoints["vooo"]("a,j,l,k")) * this->T1_("a,k");

    // sigmaOps["vo_23"] += 1.000000 this->fockMatrix_ta["vv"]("a,b") this->T1_("b,k") 
    // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_23"]("a,k") = this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,k");

    // sigmaOps["oo_24"] += 1.000000 this->fockMatrix_ta["ov"]("k,a") this->T1_("a,j") 
    // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2, 
    sigmaOps["oo_24"]("k,j") = this->fockMatrix_ta["ov"]("k,a") * this->T1_("a,j");

    // sigmaOps["vo_25"] += 1.000000 this->fockMatrix_ta["oo"]("l,k") this->T1_("a,l") 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_25"]("a,k") = this->fockMatrix_ta["oo"]("l,k") * this->T1_("a,l");

    // sigmaOps["vo_15"] += 1.000000 this->T2_("c,a,m,l") tempOp1_vooo("c,i,m,l") 
    // flops: o3v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_15"]("a,i") = this->T2_("c,a,m,l") * sigmaOps["vooo_0"]("c,i,m,l");

    // sigmaOps["oooo_17"] += 1.000000 this->T1_("b,i") tempOp1_vooo("b,j,l,k") 
    // flops: o4v1: 1, o4v0: 1 | mem: o4v0: 2, 
    sigmaOps["oooo_17"]("i,j,l,k") = this->T1_("b,i") * sigmaOps["vooo_0"]("b,j,l,k");

    // sigmaOps["vo_31"] += 1.000000 tempOp2_vo("c,m") this->T2_("c,a,k,m") 
    // flops: o2v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_31"]("a,k") = sigmaOps["vo_1"]("c,m") * this->T2_("c,a,k,m");

    // sigmaOps["oo_33"] += 1.000000 this->T1_("b,j") tempOp2_vo("b,l") 
    // flops: o2v1: 1, o2v0: 1 | mem: o2v0: 2, 
    sigmaOps["oo_33"]("j,l") = this->T1_("b,j") * sigmaOps["vo_1"]("b,l");

    // sigmaOps["vooo_3"] += 1.000000 this->T2_("b,c,j,k") conj(this->antiSymMoints["vvvo"]("b,c,a,l"))
    // flops: o3v3: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,j,k,l") = 0.5 * this->T2_("b,c,j,k") * conj(this->antiSymMoints["vvvo"]("b,c,a,l"));

    // sigmaOps["vooo_30"] += 1.000000 this->T1_("a,l") tempOp5_oooo("i,j,m,l") 
    // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,i,j,m") += 0.5 * this->T1_("a,l") * sigmaOps["oooo_4"]("i,j,m,l");

    // sigmaOps["vooo_29"] += 1.000000 tempOp18_oooo("j,k,m,l") this->T1_("a,l") 
    // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,j,k,m") -= sigmaOps["oooo_17"]("j,k,m,l") * this->T1_("a,l");
    
    // sigmaOps["vooo_26"] += 1.000000 this->T2_("c,a,i,j") tempOp2_vo("c,m") 
    // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,i,j,m") += this->T2_("c,a,i,j") * sigmaOps["vo_1"]("c,m");

    // sigmaOps["vooo_27"] += 1.000000 this->T1_("c,i") tempOp8_vvoo("a,c,l,j") 
    // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,j,k,m") -= this->T1_("c,j") * sigmaOps["vvoo_7"]("a,c,m,k");

    // sigmaOps["vooo_18"] += 1.000000 this->antiSymMoints["oooo"]("m,l,j,k") this->T1_("a,l") 
    // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_30"]("a,j,k,m") += this->antiSymMoints["oooo"]("m,l,j,k") * this->T1_("a,l");

    if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT)
       TAManager::get().free("vvvo", std::move(this->antiSymMoints["vvvo"]), true);

    // sigmaOps["vooo_5"] += 1.000000 conj(this->antiSymMoints["vooo"]("b,k,m,l")) this->T2_("b,a,j,l") 
    // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_5"]("a,m,j,k") = conj(this->antiSymMoints["vooo"]("b,j,m,l")) * this->T2_("b,a,k,l");

    // sigmaOps["vooo_28"] += 1.000000 this->T1_("a,l") tempOp17_oooo("m,l,j,k") 
    // flops: o4v1: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_5"]("a,m,j,k") += this->T1_("a,l") * sigmaOps["oooo_16"]("m,l,j,k");

    // sigmaOps["vooo_11"] += 1.000000 this->T1_("b,j") this->antiSymMoints["vovo"]("a,l,b,k") 
    // flops: o3v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_5"]("a,m,j,k") -= this->T1_("b,j") * this->antiSymMoints["vovo"]("a,m,b,k");

    // sigmaOps["vooo_6"] += 1.000000 this->T2_("c,a,k,l") tempOp1_vooo("c,j,m,l") 
    // flops: o4v2: 1, o3v1: 1 | mem: o3v1: 2, 
    sigmaOps["vooo_5"]("a,m,j,k") += this->T2_("c,a,j,l") * sigmaOps["vooo_0"]("c,k,m,l");

    // sigmaOps["vvoo_2"] += 1.000000 this->T2_("b,a,k,l") conj(this->antiSymMoints["vvoo"]("b,c,m,l")) 
    // flops: o3v3: 1, o2v2: 1 | mem: o2v2: 2, 
    sigmaOps["vvoo_7"]("a,c,m,k") += this->T2_("b,a,k,l") * conj(this->antiSymMoints["vvoo"]("b,c,m,l"));

    // sigmaOps["vo_32"] += 1.000000 this->T1_("c,i") tempOp20_vv("a,c") 
    // flops: o1v2: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_32"]("a,i") = this->T1_("c,i") * sigmaOps["vv_19"]("a,c");

    // sigmaOps["vo_34"] += 1.000000 tempOp23_oo("m,k") this->T1_("a,m") 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_34"]("a,k") = sigmaOps["oo_22"]("m,k") * this->T1_("a,m");

    // sigmaOps["vo_35"] += 1.000000 tempOp25_oo("l,k") this->T1_("a,l") 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_35"]("a,k") = sigmaOps["oo_24"]("l,k") * this->T1_("a,l");

    // sigmaOps["vo_36"] += 1.000000 this->T1_("a,m") tempOp34_oo("i,m") 
    // flops: o2v1: 1, o1v1: 1 | mem: o1v1: 2, 
    sigmaOps["vo_36"]("a,i") = this->T1_("a,m") * sigmaOps["oo_33"]("i,m");
}




  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::formR_tilde(const TArray &R2, const TArray &R3, TArray &sigmar_R2, TArray &sigmar_R3) const {
    TAManager &TAmanager = TAManager::get();
#ifndef DEBUG_DIP

    // tempOps["vooo_4"] += 1.000000 r3_Lvooo("I,c,j,k,l") sigmaOps["vvoo_7"]("a,c,l,i")
    // flops: o4v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_4"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_4"]("a,j,k,i") = R3("c,j,k,l") * sigmaOps["vvoo_7"]("a,c,l,i");

    // sigmar_R3 += -1.000000 <l,a||b,c> r3(c,j,k,l) t1(b,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") = tempOps["vooo_4"]("a,j,k,i");

    // sigmar_R3 += -1.000000 P(j,k) <l,a||b,c> r3(c,i,j,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_4"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_4"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_4"]));

    // tempOps["vooo_5"] += 1.000000 r3_Lvooo("I,b,i,j,l") this->antiSymMoints["vovo"]("a,l,b,k")
    // flops: o4v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_5"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_5"]("a,i,j,k") = R3("b,i,j,l") * this->antiSymMoints["vovo"]("a,l,b,k");

    // sigmar_R3 += 1.000000 <l,a||b,i> r3(b,j,k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_5"]("a,j,k,i");

    // sigmar_R3 += 1.000000 P(j,k) <l,a||b,k> r3(b,i,j,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_5"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_5"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_5"]));

    // tempOps["vooo_6"] += 0.500000 sigmaOps["oooo_16"]("m,l,i,k") r3_Lvooo("I,a,j,m,l")
    // flops: o5v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_6"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_6"]("a,i,k,j") = 0.500000 * sigmaOps["oooo_16"]("m,l,i,k") * R3("a,j,m,l");

    // sigmar_R3 += 0.500000 P(j,k) <m,l||b,i> r3(a,j,m,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_6"]("a,i,k,j");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_6"]("a,i,j,k");

    // sigmar_R3 += -0.500000 P(i,k) <m,l||b,j> r3(a,i,m,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_6"]("a,j,k,i");
    sigmar_R3("a,i,j,k") += tempOps["vooo_6"]("a,j,i,k");

    // sigmar_R3 += 0.500000 P(i,j) <m,l||b,k> r3(a,i,m,l) t1(b,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_6"]("a,k,j,i");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_6"]("a,k,i,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_6"]));


    // tempOps["vooo_9"] += 1.000000 conj(this->antiSymMoints["vooo"]("b,k,m,l")) r3_Lvooo("I,b,i,j,m") this->T1_("a,l")
    // flops: o5v1L1: 1, o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, o4v0L1: 1,
    tempOps["vooo_9"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_9"]("a,k,i,j") = conj(this->antiSymMoints["vooo"]("b,k,m,l")) * R3("b,i,j,m") * this->T1_("a,l");

    // sigmar_R3 += -1.000000 P(j,k) <m,l||b,k> r3(b,i,j,m) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_9"]("a,k,i,j");
    sigmar_R3("a,i,j,k") += tempOps["vooo_9"]("a,j,i,k");

    // sigmar_R3 += -1.000000 <m,l||b,i> r3(b,j,k,m) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_9"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_9"]));

    // tempOps["vooo_10"] += 1.000000 sigmaOps["vooo_0"]("c,k,m,l") r3_Lvooo("I,c,i,j,m") this->T1_("a,l")
    // flops: o5v1L1: 1, o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2, o4v0L1: 1,
    tempOps["vooo_10"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_10"]("a,k,i,j") = sigmaOps["vooo_0"]("c,k,m,l") * R3("c,i,j,m") * this->T1_("a,l");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,c> r3(c,i,j,m) t1(b,k) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_10"]("a,k,i,j");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_10"]("a,j,i,k");

    // sigmar_R3 += 1.000000 <m,l||b,c> r3(c,j,k,m) t1(b,i) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_10"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_10"]));


    // tempOps["vooo_15"] += 0.500000 r3_Lvooo("I,a,k,m,l") this->antiSymMoints["oooo"]("m,l,i,j")
    // flops: o5v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_15"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_15"]("a,k,i,j") = 0.500000 * R3("a,k,m,l") * this->antiSymMoints["oooo"]("m,l,i,j");

    // sigmar_R3 += 0.500000 <m,l||i,j> r3(a,k,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_15"]("a,k,i,j");

    // sigmar_R3 += 0.500000 P(i,j) <m,l||j,k> r3(a,i,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_15"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_15"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_15"]));

    // tempOps["vooo_16"] += 0.250000 r3_Lvooo("I,a,i,m,l") sigmaOps["oooo_4"]("j,k,m,l")
    // flops: o5v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_16"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_16"]("a,i,j,k") = 0.250000 * R3("a,i,m,l") * sigmaOps["oooo_4"]("j,k,m,l");

    // sigmar_R3 += 0.250000 P(i,j) <m,l||b,c> r3(a,i,m,l) t2(b,c,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_16"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_16"]("a,j,i,k");

    // sigmar_R3 += 0.250000 <m,l||b,c> r3(a,k,m,l) t2(b,c,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_16"]("a,k,i,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_16"]));

    // tempOps["vooo_17"] += 0.500000 r3_Lvooo("I,a,k,m,l") sigmaOps["oooo_17"]("i,j,m,l")
    // flops: o5v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_17"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_17"]("a,k,i,j") = 0.500000 * R3("a,k,m,l") * sigmaOps["oooo_17"]("i,j,m,l");

    // sigmar_R3 += -0.500000 <m,l||b,c> r3(a,k,m,l) t1(b,j) t1(c,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_17"]("a,k,i,j");

    // sigmar_R3 += -0.500000 P(i,j) <m,l||b,c> r3(a,i,m,l) t1(b,k) t1(c,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_17"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_17"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_17"]));

    // tempOps["vooo_19"] += 0.500000 conj(this->antiSymMoints["vvoo"]("b,c,m,l")) r3_Lvooo("I,c,k,m,l") this->T2_("b,a,i,j")
    // flops: o3v2L1: 2, o3v1L1: 1 | mem: o3v1L1: 2, o1v1L1: 1,
    tempOps["vooo_19"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_19"]("a,k,i,j") = 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,m,l")) * R3("c,k,m,l") * this->T2_("b,a,i,j");

    // sigmar_R3 += -0.500000 P(i,j) <m,l||b,c> r3(c,i,m,l) t2(b,a,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_19"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_19"]("a,j,i,k");

    // sigmar_R3 += -0.500000 <m,l||b,c> r3(c,k,m,l) t2(b,a,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_19"]("a,k,i,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_19"]));

    // tempOps["vooo_20"] += 0.500000 conj(this->antiSymMoints["vooo"]("b,i,m,l")) r2_Loo("I,m,l") this->T2_("b,a,j,k")
    // flops: o3v2L1: 1, o3v1L1: 2 | mem: o3v1L1: 2, o1v1L1: 1,
    tempOps["vooo_20"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_20"]("a,i,j,k") = 0.500000 * conj(this->antiSymMoints["vooo"]("b,i,m,l")) * R2("m,l") * this->T2_("b,a,j,k");

    // sigmar_R3 += -0.500000 <m,l||b,i> r2(m,l) t2(b,a,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_20"]("a,i,j,k");

    // sigmar_R3 += -0.500000 P(j,k) <m,l||b,k> r2(m,l) t2(b,a,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_20"]("a,k,i,j");
    sigmar_R3("a,i,j,k") += tempOps["vooo_20"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_20"]));

    // tempOps["vooo_21"] += 0.500000 r2_Loo("I,m,l") sigmaOps["vooo_0"]("c,i,m,l") this->T2_("c,a,j,k")
    // flops: o3v2L1: 1, o3v1L1: 2 | mem: o3v1L1: 2, o1v1L1: 1,
    tempOps["vooo_21"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_21"]("a,i,j,k") = 0.500000 * R2("m,l") * sigmaOps["vooo_0"]("c,i,m,l") * this->T2_("c,a,j,k");

    // sigmar_R3 += 0.500000 P(j,k) <m,l||b,c> r2(m,l) t1(b,k) t2(c,a,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_21"]("a,k,i,j");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_21"]("a,j,i,k");

    // sigmar_R3 += 0.500000 <m,l||b,c> r2(m,l) t1(b,i) t2(c,a,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_21"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_21"]));

    // tempOps["vooo_23"] += 1.000000 r2_Loo("I,i,l") sigmaOps["vooo_11"]("a,k,l,j")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_23"] = TAmanager.malloc<MatsT>("vooo");
    // tempOps["vooo_25"] += 1.000000 r2_Loo("I,i,m") sigmaOps["vooo_5"]("a,m,j,k")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    // tempOps["vooo_24"] += 1.000000 r2_Loo("I,i,m") sigmaOps["vooo_28"]("a,m,j,k")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_23"]("a,i,j,k") =- R2("i,m") * sigmaOps["vooo_5"]("a,m,j,k");
    // tempOps["vooo_26"] += 1.000000 r2_Loo("I,i,m") sigmaOps["vooo_6"]("a,j,k,m")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,

    // sigmar_R3 += 1.000000 P(i,k) <l,a||b,j> r2(i,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_23"]("a,i,k,j");
    sigmar_R3("a,i,j,k") += tempOps["vooo_23"]("a,k,i,j");

    // sigmar_R3 += -1.000000 P(j,k) <l,a||b,i> r2(j,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_23"]("a,j,k,i");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_23"]("a,k,j,i");

    // sigmar_R3 += -1.000000 P(i,j) <l,a||b,k> r2(i,l) t1(b,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_23"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_23"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_23"]));

    // tempOps["vooo_34"] += 1.000000 sigmaOps["oo_33"]("i,m") r3_Lvooo("I,a,j,k,m")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_34"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_34"]("a,i,j,k") = sigmaOps["oo_33"]("i,m") * R3("a,j,k,m");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,c> r3(a,i,j,m) t1(b,l) t1(c,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_34"]("a,k,i,j");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_34"]("a,j,i,k");

    // sigmar_R3 += 1.000000 <m,l||b,c> r3(a,j,k,m) t1(b,l) t1(c,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_34"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_34"]));

    // tempOps["vooo_35"] += 0.500000 r2_Loo("I,i,m") sigmaOps["vooo_30"]("a,j,k,m")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_35"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_35"]("a,i,j,k") = R2("i,m") * sigmaOps["vooo_30"]("a,j,k,m");

    // sigmar_R3 += 0.500000 <m,l||b,c> r2(k,m) t1(a,l) t2(b,c,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_35"]("a,k,i,j");

    // sigmar_R3 += 0.500000 P(i,j) <m,l||b,c> r2(i,m) t1(a,l) t2(b,c,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_35"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_35"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_35"]));

    // tempOps["vooo_44"] += 1.000000 sigmaOps["vooo_14"]("a,j,k,l") r2_Loo("I,i,l")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_44"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_44"]("a,j,k,i") = sigmaOps["vooo_14"]("a,j,k,l") * R2("i,l");

    // sigmar_R3 += -1.000000 P(i,j) f(l,b) r2(i,l) t2(b,a,j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_44"]("a,j,k,i");
    sigmar_R3("a,i,j,k") += tempOps["vooo_44"]("a,i,k,j");

    // sigmar_R3 += -1.000000 f(l,b) r2(k,l) t2(b,a,i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_44"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_44"]));

    // tempOps["vooo_37"] += 1.000000 this->antiSymMoints["vooo"]("a,l,i,j") r2_Loo("I,k,l")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_37"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_37"]("a,i,j,k") = this->antiSymMoints["vooo"]("a,l,i,j") * R2("k,l");

    // sigmar_R3 += -1.000000 <l,a||i,j> r2(k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_37"]("a,i,j,k");

    // sigmar_R3 += -1.000000 P(i,j) <l,a||j,k> r2(i,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_37"]("a,j,k,i");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_37"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_37"]));

    // tempOps["vooo_38"] += 1.000000 r3_Lvooo("I,a,i,j,m") sigmaOps["oo_22"]("m,k")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_38"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_38"]("a,i,j,k") = R3("a,i,j,m") * sigmaOps["oo_22"]("m,k");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,k> r3(a,i,j,m) t1(b,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_38"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_38"]("a,i,k,j");

    // sigmar_R3 += 1.000000 <m,l||b,i> r3(a,j,k,m) t1(b,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_38"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_38"]));

    // tempOps["vooo_41"] += 0.500000 sigmaOps["oo_10"]("k,m") r3_Lvooo("I,a,i,j,m")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_41"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_41"]("a,k,i,j") = 0.500000 * sigmaOps["oo_10"]("k,m") * R3("a,i,j,m");

    // sigmar_R3 += -0.500000 <m,l||b,c> r3(a,j,k,m) t2(b,c,i,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_41"]("a,i,j,k");

    // sigmar_R3 += -0.500000 P(j,k) <m,l||b,c> r3(a,i,j,m) t2(b,c,k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_41"]("a,k,i,j");
    sigmar_R3("a,i,j,k") += tempOps["vooo_41"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tempOps["vooo_41"]));

    // tempOps["vooo_42"] += 1.000000 r3_Lvooo("I,a,i,j,l") sigmaOps["oo_24"]("l,k")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_42"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_42"]("a,i,j,k") = R3("a,i,j,l") * sigmaOps["oo_24"]("l,k");

    // sigmar_R3 += -1.000000 P(j,k) f(l,b) r3(a,i,j,l) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_42"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_42"]("a,i,k,j");

    // sigmar_R3 += -1.000000 f(l,b) r3(a,j,k,l) t1(b,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_42"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_42"]));

    // tempOps["vooo_45"] += 1.000000 r3_Lvooo("I,a,i,j,l") this->fockMatrix_ta["oo"]("l,k")
    // flops: o4v1L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    tempOps["vooo_45"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_45"]("a,i,j,k") = R3("a,i,j,l") * this->fockMatrix_ta["oo"]("l,k");

    // sigmar_R3 += -1.000000 P(j,k) f(l,k) r3(a,i,j,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_45"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_45"]("a,i,k,j");

    // sigmar_R3 += -1.000000 f(l,i) r3(a,j,k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_45"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_45"]));


    // tempOps["vooo_50"] += 1.000000 r2_Loo("I,j,k") sigmaOps["vo_36"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_50"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_50"]("a,j,k,i") = R2("j,k") * sigmaOps["vo_36"]("a,i");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,c> r2(i,j) t1(b,l) t1(c,k) t1(a,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_50"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_50"]("a,i,k,j");

    // sigmar_R3 += 1.000000 <m,l||b,c> r2(j,k) t1(b,l) t1(c,i) t1(a,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_50"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_50"]));

    // tempOps["vooo_51"] += 1.000000 r2_Loo("I,i,j") sigmaOps["vo_35"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_51"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_51"]("a,i,j,k") = R2("i,j") * sigmaOps["vo_35"]("a,k");

    // sigmar_R3 += -1.000000 P(j,k) f(l,b) r2(i,j) t1(b,k) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_51"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_51"]("a,i,k,j");

    // sigmar_R3 += -1.000000 f(l,b) r2(j,k) t1(b,i) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_51"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_51"]));

    // tempOps["vooo_52"] += 1.000000 r2_Loo("I,j,k") sigmaOps["vo_34"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_52"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_52"]("a,j,k,i") = R2("j,k") * sigmaOps["vo_34"]("a,i");

    // sigmar_R3 += 1.000000 <m,l||b,i> r2(j,k) t1(b,l) t1(a,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_52"]("a,j,k,i");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,k> r2(i,j) t1(b,l) t1(a,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_52"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_52"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_52"]));

    // tempOps["vooo_53"] += 1.000000 r2_Loo("I,i,j") sigmaOps["vo_32"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_53"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_53"]("a,i,j,k") = R2("i,j") * sigmaOps["vo_32"]("a,k");

    // sigmar_R3 += 1.000000 <l,a||b,c> r2(j,k) t1(b,l) t1(c,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_53"]("a,j,k,i");

    // sigmar_R3 += 1.000000 P(j,k) <l,a||b,c> r2(i,j) t1(b,l) t1(c,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_53"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_53"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_53"]));

    // tempOps["vooo_54"] += 1.000000 r2_Loo("I,j,k") sigmaOps["vo_25"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_54"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_54"]("a,j,k,i") = R2("j,k") * sigmaOps["vo_25"]("a,i");

    // sigmar_R3 += -1.000000 f(l,i) r2(j,k) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_54"]("a,j,k,i");

    // sigmar_R3 += -1.000000 P(j,k) f(l,k) r2(i,j) t1(a,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_54"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_54"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_54"]));

    // tempOps["vooo_55"] += 1.000000 r2_Loo("I,j,k") sigmaOps["vo_23"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_55"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_55"]("a,j,k,i") = R2("j,k") * sigmaOps["vo_23"]("a,i");

    // sigmar_R3 += 1.000000 P(j,k) f(a,b) r2(i,j) t1(b,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_55"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_55"]("a,i,k,j");

    // sigmar_R3 += 1.000000 f(a,b) r2(j,k) t1(b,i)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_55"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_55"]));

    // tempOps["vooo_56"] += 1.000000 r2_Loo("I,j,k") sigmaOps["vo_21"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_56"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_56"]("a,j,k,i") = R2("j,k") * sigmaOps["vo_21"]("a,i");

    // sigmar_R3 += -1.000000 P(j,k) f(l,b) r2(i,j) t2(b,a,k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_56"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_56"]("a,i,k,j");

    // sigmar_R3 += -1.000000 f(l,b) r2(j,k) t2(b,a,i,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_56"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_56"]));

    // tempOps["vooo_57"] += 1.000000 r2_Loo("I,i,j") sigmaOps["vo_20"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_57"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_57"]("a,i,j,k") = R2("i,j") * sigmaOps["vo_20"]("a,k");

    // sigmar_R3 += 1.000000 <l,a||b,i> r2(j,k) t1(b,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_57"]("a,j,k,i");

    // sigmar_R3 += 1.000000 P(j,k) <l,a||b,k> r2(i,j) t1(b,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_57"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_57"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_57"]));

    // tempOps["vooo_58"] += 0.500000 r2_Loo("I,i,j") sigmaOps["vo_12"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_58"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_58"]("a,i,j,k") = 0.500000 * R2("i,j") * sigmaOps["vo_12"]("a,k");

    // sigmar_R3 += 0.500000 P(j,k) <m,l||b,c> r2(i,j) t1(a,l) t2(b,c,k,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_58"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_58"]("a,i,k,j");

    // sigmar_R3 += 0.500000 <m,l||b,c> r2(j,k) t1(a,l) t2(b,c,i,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_58"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_58"]));

    // tempOps["vooo_59"] += 0.500000 r2_Loo("I,i,j") sigmaOps["vo_15"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_59"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_59"]("a,i,j,k") = 0.500000 * R2("i,j") * sigmaOps["vo_15"]("a,k");

    // sigmar_R3 += 0.500000 P(j,k) <m,l||b,c> r2(i,j) t1(b,k) t2(c,a,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_59"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_59"]("a,i,k,j");

    // sigmar_R3 += 0.500000 <m,l||b,c> r2(j,k) t1(b,i) t2(c,a,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_59"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_59"]));

    // tempOps["vooo_60"] += 0.500000 r2_Loo("I,i,j") sigmaOps["vo_8"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_60"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_60"]("a,i,j,k") = 0.500000 * R2("i,j") * sigmaOps["vo_8"]("a,k");

    // sigmar_R3 += -0.500000 P(j,k) <l,a||b,c> r2(i,j) t2(b,c,k,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_60"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_60"]("a,i,k,j");

    // sigmar_R3 += -0.500000 <l,a||b,c> r2(j,k) t2(b,c,i,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_60"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_60"]));

    // tempOps["vooo_61"] += 1.000000 r2_Loo("I,i,j") sigmaOps["vo_31"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_61"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_61"]("a,i,j,k") = R2("i,j") * sigmaOps["vo_31"]("a,k");

    // sigmar_R3 += 1.000000 <m,l||b,c> r2(j,k) t1(b,l) t2(c,a,i,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_61"]("a,j,k,i");

    // sigmar_R3 += 1.000000 P(j,k) <m,l||b,c> r2(i,j) t1(b,l) t2(c,a,k,m)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_61"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_61"]("a,i,k,j");
    TAmanager.free("vooo", std::move(tempOps["vooo_61"]));

    // tempOps["vooo_62"] += 1.000000 r2_Loo("I,i,j") this->fockMatrix_ta["vo"]("a,k")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_62"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_62"]("a,i,j,k") = R2("i,j") * this->fockMatrix_ta["vo"]("a,k");

    // sigmar_R3 += 1.000000 P(j,k) f(a,k) r2(i,j)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_62"]("a,i,j,k");
    sigmar_R3("a,i,j,k") -= tempOps["vooo_62"]("a,i,k,j");

    // sigmar_R3 += 1.000000 f(a,i) r2(j,k)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") += tempOps["vooo_62"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_62"]));

    // tempOps["vooo_63"] += 0.500000 r2_Loo("I,j,k") sigmaOps["vo_13"]("a,i")
    // flops: o3v1L1: 2 | mem: o3v1L1: 2,
    tempOps["vooo_63"] = TAmanager.malloc<MatsT>("vooo");
    tempOps["vooo_63"]("a,j,k,i") = 0.500000 * R2("j,k") * sigmaOps["vo_13"]("a,i");

    // sigmar_R3 += -0.500000 P(j,k) <m,l||b,k> r2(i,j) t2(b,a,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_63"]("a,i,j,k");
    sigmar_R3("a,i,j,k") += tempOps["vooo_63"]("a,i,k,j");

    // sigmar_R3 += -0.500000 <m,l||b,i> r2(j,k) t2(b,a,m,l)
    // flops: o3v1: 1 | mem: o3v1: 1,
    sigmar_R3("a,i,j,k") -= tempOps["vooo_63"]("a,j,k,i");
    TAmanager.free("vooo", std::move(tempOps["vooo_63"]));
#endif

    // sigmar_R2 += 0.500000 <l,k||i,j> r2(l,k)
    // flops: o4v0L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
    sigmar_R2("i,j") = 0.500000 * this->antiSymMoints["oooo"]("l,k,i,j") * R2("l,k");

    // sigmar_R2 += 0.250000 <l,k||a,b> r2(l,k) t2(a,b,i,j)
    // flops: o2v2L1: 2, o2v0L1: 1 | mem: o0v2L1: 1, o2v0L1: 2,
    sigmar_R2("i,j") += 0.250000 * conj(this->antiSymMoints["vvoo"]("a,b,l,k")) * R2("l,k") * this->T2_("a,b,i,j");

    // sigmar_R2 += 1.000000 f(k,a) r3(a,i,j,k)
    // flops: o3v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
    sigmar_R2("i,j") += this->fockMatrix_ta["ov"]("k,a") * R3("a,i,j,k");

    // sigmar_R2 += -0.500000 <l,k||a,b> r2(l,k) t1(a,j) t1(b,i)
    // flops: o4v0L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
    sigmar_R2("i,j") -= 0.500000 * sigmaOps["oooo_17"]("i,j,l,k") * R2("l,k");

    // sigmar_R2 += -0.500000 P(i,j) <l,k||a,b> r3(b,i,l,k) t1(a,j)
    // flops: o4v1L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = 0.500000 * sigmaOps["vooo_0"]("b,j,l,k") * R3("b,i,l,k");
    sigmar_R2("i,j") -= tempPerm_oo("i,j");
    sigmar_R2("i,j") += tempPerm_oo("j,i");

    // sigmar_R2 += -1.000000 P(i,j) f(k,j) r2(i,k)
    // flops: o3v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = this->fockMatrix_ta["oo"]("k,j") * R2("i,k");
    sigmar_R2("i,j") -= tempPerm_oo("i,j");
    sigmar_R2("i,j") += tempPerm_oo("j,i");

    // sigmar_R2 += -1.000000 P(i,j) f(k,a) r2(i,k) t1(a,j)
    // flops: o3v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = sigmaOps["oo_24"]("k,j") * R2("i,k");
    sigmar_R2("i,j") -= tempPerm_oo("i,j");
    sigmar_R2("i,j") += tempPerm_oo("j,i");

//    // sigmar_R2 += 0.250000 <l,k||a,b> r2(i,j) t2(a,b,l,k)
//    // flops: o2v0L1: 1, o0v0: 1 | mem: o2v0L1: 1, o0v0: 1,
//    sigmar_R2("i,j") += 0.250000 * scalar_3 * R2("i,j");

//    // sigmar_R2 += -0.500000 <l,k||a,b> r2(i,j) t1(a,k) t1(b,l)
//    // flops: o2v0L1: 1, o0v0: 1 | mem: o2v0L1: 1, o0v0: 1,
//    sigmar_R2("i,j") -= 0.500000 * scalar_4 * R2("i,j");

    // sigmar_R2 += 1.000000 P(i,j) <l,k||a,b> r2(i,l) t1(a,k) t1(b,j)
    // flops: o3v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = sigmaOps["oo_33"]("j,l") * R2("i,l");
    sigmar_R2("i,j") += tempPerm_oo("i,j");
    sigmar_R2("i,j") -= tempPerm_oo("j,i");

    // sigmar_R2 += 1.000000 P(i,j) <l,k||a,j> r2(i,l) t1(a,k)
    // flops: o3v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = sigmaOps["oo_22"]("l,j") * R2("i,l");
    sigmar_R2("i,j") += tempPerm_oo("i,j");
    sigmar_R2("i,j") -= tempPerm_oo("j,i");

    // sigmar_R2 += 0.500000 P(i,j) <l,k||a,j> r3(a,i,l,k)
    // flops: o4v1L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = 0.500000 * conj(this->antiSymMoints["vooo"]("a,j,l,k")) * R3("a,i,l,k");
    sigmar_R2("i,j") += tempPerm_oo("i,j");
    sigmar_R2("i,j") -= tempPerm_oo("j,i");

//    // sigmar_R2 += 1.000000 f(k,a) r2(i,j) t1(a,k)
//    // flops: o2v0L1: 1, o0v0: 1 | mem: o2v0L1: 1, o0v0: 1,
//    sigmar_R2("i,j") += scalar_1 * R2("i,j");

    // sigmar_R2 += -0.500000 P(i,j) <l,k||a,b> r2(i,l) t2(a,b,j,k)
    // flops: o3v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = 0.500000 * sigmaOps["oo_10"]("j,l") * R2("i,l");
    sigmar_R2("i,j") -= tempPerm_oo("i,j");
    sigmar_R2("i,j") += tempPerm_oo("j,i");

//    // sigmar_R2 += 1.000000 f(k,k) r2(i,j)
//    // flops: o2v0L1: 1, o0v0: 1 | mem: o2v0L1: 1, o0v0: 1,
//    sigmar_R2("i,j") += scalar_0 * R2("i,j");

    // sigmar_R2 += -1.000000 <l,k||a,b> r3(b,i,j,l) t1(a,k)
    // flops: o3v1L1: 1, o2v0L1: 1 | mem: o2v0L1: 2,
    sigmar_R2("i,j") -= sigmaOps["vo_1"]("b,l") * R3("b,i,j,l");

    // sigmar_R2 += 0.500000 P(i,j) <l,k||a,j> r2(l,k) t1(a,i)
    // flops: o4v0L1: 1, o2v0L1: 2 | mem: o2v0L1: 2,
    tempPerm_oo("i,j") = 0.500000 * sigmaOps["oooo_16"]("l,k,j,i") * R2("l,k");
    sigmar_R2("i,j") += tempPerm_oo("i,j");
    sigmar_R2("i,j") -= tempPerm_oo("j,i");

//    // sigmar_R2 += -0.500000 <l,k||l,k> r2(i,j)
//    // flops: o2v0L1: 1, o0v0: 1 | mem: o2v0L1: 1, o0v0: 1,
//    sigmar_R2("i,j") -= 0.500000 * scalar_2 * R2("i,j");

#ifndef DEBUG_DIP
    // sigmar_R3 += 1.000000 f(a,b) r3(b,i,j,k)
    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    sigmar_R3("a,i,j,k") += this->fockMatrix_ta["vv"]("a,b") * R3("b,i,j,k");

//    // sigmar_R3 += 1.000000 f(l,l) r3(a,i,j,k)
//    // flops: o3v1L1: 1, o0v0: 1 | mem: o3v1L1: 1, o0v0: 1,
//    sigmar_R3("a,i,j,k") += scalar_0 * R3("a,i,j,k");
//
//    // sigmar_R3 += 1.000000 f(l,b) r3(a,i,j,k) t1(b,l)
//    // flops: o3v1L1: 1, o0v0: 1 | mem: o3v1L1: 1, o0v0: 1,
//    sigmar_R3("a,i,j,k") += scalar_1 * R3("a,i,j,k");
//
//    // sigmar_R3 += -0.500000 <m,l||m,l> r3(a,i,j,k)
//    // flops: o3v1L1: 1, o0v0: 1 | mem: o3v1L1: 1, o0v0: 1,
//    sigmar_R3("a,i,j,k") -= 0.500000 * scalar_2 * R3("a,i,j,k");

    // sigmar_R3 += 1.000000 <l,a||b,c> r3(c,i,j,k) t1(b,l)
    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    sigmar_R3("a,i,j,k") -= sigmaOps["vv_19"]("a,c") * R3("c,i,j,k");

    // sigmar_R3 += 1.000000 <m,l||b,c> r3(c,i,j,k) t1(b,l) t1(a,m)
    // flops: o4v1L1: 2, o3v1L1: 1 | mem: o3v1L1: 2, o4v0L1: 1,
    sigmar_R3("a,i,j,k") += sigmaOps["vo_1"]("c,m") * R3("c,i,j,k") * this->T1_("a,m");

//    // sigmar_R3 += 0.250000 <m,l||b,c> r3(a,i,j,k) t2(b,c,m,l)
//    // flops: o3v1L1: 1, o0v0: 1 | mem: o3v1L1: 1, o0v0: 1,
//    sigmar_R3("a,i,j,k") += 0.250000 * scalar_3 * R3("a,i,j,k");

    // sigmar_R3 += -0.500000 <m,l||b,c> r3(c,i,j,k) t2(b,a,m,l)
    // flops: o3v2L1: 1, o3v1L1: 1 | mem: o3v1L1: 2,
    sigmar_R3("a,i,j,k") -= 0.500000 * sigmaOps["vv_9"]("c,a") * R3("c,i,j,k");

//    // sigmar_R3 += -0.500000 <m,l||b,c> r3(a,i,j,k) t1(b,l) t1(c,m)
//    // flops: o3v1L1: 1, o0v0: 1 | mem: o3v1L1: 1, o0v0: 1,
//    sigmar_R3("a,i,j,k") -= 0.500000 * scalar_4 * R3("a,i,j,k");

    // sigmar_R3 += -1.000000 f(l,b) r3(b,i,j,k) t1(a,l)
    // flops: o4v1L1: 2, o3v1L1: 1 | mem: o3v1L1: 2, o4v0L1: 1,
    sigmar_R3("a,i,j,k") -= this->fockMatrix_ta["ov"]("l,b") * R3("b,i,j,k") * this->T1_("a,l");
#endif
  }



  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::formL_tilde(const TArray &L2, const TArray &L3, TArray &sigmal_L2, TArray &sigmal_L3) const{} 


  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {
    const TArray &V2 =  V.get_tensor("OneBody");
    const TArray &V3 =  V.get_tensor("TwoBody");
    TArray &HV2 = HV.get_tensor("OneBody");
    TArray &HV3 = HV.get_tensor("TwoBody");
    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR_tilde(V2, V3, HV2, HV3);
        break;
      case EOMCCEigenVecType::LEFT:
        formL_tilde(V2, V3, HV2, HV3);
        break;
    }
  }


  // DIP does not have ground state
  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::buildRightZeroBody(size_t nVec) {}

  template <typename MatsT>
  inline size_t EOMDIP_3h1p<MatsT>::toCompoundS(size_t i, size_t j) const {
    if (i >= nO_ or j >= nO_)
      return outOfBound_;
    return ijIndices_[i][j];
  }

  //template <typename MatsT>
  //inline double EOMDIP_3h1p<MatsT>::signS(size_t i, size_t j) const {
  //  return signD(i,j);
  //}
  
  //template <typename MatsT>
  //inline double EOMDIP_3h1p<MatsT>::signD(size_t a, size_t i, size_t j, size_t k) const {
  //  return signT(i,j,k);
  //}


  template <typename MatsT>
  inline double EOMDIP_3h1p<MatsT>::signD(size_t i, size_t j) const {
    if (i == j) {
      i = 0;
      j = 0;
      return 0.0;
    }
    double sign = 1.0;
    if (i > j) {
      std::swap(i,j);
      sign *= -1.0;
    }
    return sign;
  }
  template <typename MatsT>
  inline double EOMDIP_3h1p<MatsT>::signT(size_t i, size_t j, size_t k) const {
    if (i == j or i == k or j == k) {
      i = 0;
      j = 0;
      k = 0;
      return 0.0;
    }
    double sign = 1.0;
    if (i > j) {
      std::swap(i,j);
      sign *= -1.0;
    }
    if (j > k) {
      std::swap(j,k);
      sign *= -1.0;
    }
    if (i > j) {
      std::swap(i,j);
      sign *= -1.0;
    }
    return sign;
  }
  template <typename MatsT>
  inline size_t EOMDIP_3h1p<MatsT>::toCompoundD(size_t a, size_t i, size_t j, size_t k) const {
    size_t ijk = ijkIndices_[i][j][k];
    if (ijk == outOfBound_)
      return outOfBound_;
    return a + ijk * nV_;
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMDIP_3h1p<MatsT>::toCompoundSS(size_t i, size_t j, size_t k, size_t l, size_t ldH) const {
    size_t ij = toCompoundS(i,j), kl = toCompoundS(k,l);
    if (ij == outOfBound_ or kl == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(i,j);
    sign *= signD(k,l);
    return std::make_pair(ij + kl * ldH, sign);
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMDIP_3h1p<MatsT>::toCompoundSD(size_t i, size_t j,
                                                                         size_t a, size_t k, size_t l, size_t m, size_t ldH) const {
    size_t ij = toCompoundS(i,j), aklm = toCompoundD(a,k,l,m);
    if (ij == outOfBound_ or aklm == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(i,j);
    sign *= signT(k,l,m);
    return std::make_pair(ij + aklm * ldH, sign);
  }


  template <typename MatsT>
  inline std::pair<size_t, double> EOMDIP_3h1p<MatsT>::toCompoundDS(size_t a, size_t i, size_t j, size_t k,
                                                                         size_t l, size_t m, size_t ldH) const {
    size_t lm = toCompoundS(l,m), aijk = toCompoundD(a,i,j,k);
    if (lm == outOfBound_ or aijk == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(l,m);
    sign *= signT(i,j,k);
    return std::make_pair(aijk + lm * ldH, sign);
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMDIP_3h1p<MatsT>::toCompoundDD(size_t a, size_t i, size_t j, size_t k,
                                                                         size_t b, size_t l, size_t m, size_t n, size_t ldH) const {
    size_t aijk = toCompoundD(a,i,j,k), blmn = toCompoundD(b,l,m,n);
    if (aijk == outOfBound_ or blmn == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signT(l,m,n);
    sign *= signT(i,j,k);
    return std::make_pair(aijk + blmn * ldH, sign);
  }

  template <typename MatsT>
  EOMDIP_3h1p<MatsT>::~EOMDIP_3h1p(){ 
    TAManager &TAmanager = TAManager::get();
    TAmanager.free("oo", std::move(tempPerm_oo), true);
    //TAmanager.free("vooo", std::move(tempPerm_vooo), true);
    TAmanager.free("oo", std::move(Id_oo), true);
    //TAmanager.free("oooo", std::move(Id_oooo), true);
    TAmanager.free("vooo",std::move(sigmaOps["vooo_0"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_1"]), true);
    TAmanager.free("oooo",std::move(sigmaOps["oooo_4"]), true);
    TAmanager.free("vooo",std::move(sigmaOps["vooo_5"]), true);
    TAmanager.free("vvoo",std::move(sigmaOps["vvoo_7"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_8"]), true);
    TAmanager.free("vv",std::move(sigmaOps["vv_9"]), true);
    TAmanager.free("oo",std::move(sigmaOps["oo_10"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_12"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_13"]), true);
    TAmanager.free("vooo",std::move(sigmaOps["vooo_14"]), true);
    TAmanager.free("oooo",std::move(sigmaOps["oooo_16"]), true);
    TAmanager.free("vv",std::move(sigmaOps["vv_19"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_20"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_21"]), true);
    TAmanager.free("oo",std::move(sigmaOps["oo_22"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_23"]), true);
    TAmanager.free("oo",std::move(sigmaOps["oo_24"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_25"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_15"]), true);
    TAmanager.free("oooo",std::move(sigmaOps["oooo_17"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_31"]), true);
    TAmanager.free("oo",std::move(sigmaOps["oo_33"]), true);
    TAmanager.free("vooo",std::move(sigmaOps["vooo_30"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_32"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_34"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_35"]), true);
    TAmanager.free("vo",std::move(sigmaOps["vo_36"]), true);


  }
  //for full_diagonalization
   template <typename MatsT>
  cqmatrix::Matrix<MatsT> EOMDIP_3h1p<MatsT>::buildHbar(bool includeGroundState) const{
    if (includeGroundState) CErr("DIP should not include Ground State in the Hamiltonian.");

    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    cqmatrix::Matrix<MatsT> fullMat(this->Hbar_dim);
    fullMat.clear();

    MatsT * Hbar = fullMat.pointer();
    size_t ldH = fullMat.nRows();
    size_t nCol = fullMat.nColumns();

    MatsT * HbarSS = Hbar;
    MatsT * HbarSD = HbarSS + nO2shift_ * ldH;
    MatsT * HbarDS = HbarSS + nO2shift_;
    MatsT * HbarDD = HbarSD + nO2shift_;

    TArray H_oooo = TAmanager.malloc<MatsT>("oooo");  
    TArray H_oovooo = TAmanager.malloc<MatsT>("oovooo");  
    TArray H_vooooo = TAmanager.malloc<MatsT>("vooooo");  
    TArray H_vooovooo = TAmanager.malloc<MatsT>("vooovooo");  
 
#ifdef DEBUG_DIP
    std::cout<<"fock"<<std::endl<<this->fockMatrix_ta["oo"]<<this->fockMatrix_ta["ov"]<<this->fockMatrix_ta["vo"]<<this->fockMatrix_ta["vv"]<<std::endl;
    std::cout<<"eri"<<std::endl<<this->antiSymMoints["oooo"]<<this->antiSymMoints["vooo"]<<this->antiSymMoints["vovo"]<<this->antiSymMoints["vvoo"]<<this->antiSymMoints["vvvo"]<<this->antiSymMoints["vvvv"]<<std::endl;
    std::cout<<"t1"<<std::endl<<this->T1_<<std::endl;
    std::cout<<"t2"<<std::endl<<this->T2_<<std::endl;
    buildHbarSS(H_oooo);
#else
    buildHbarTA(H_oooo, H_oovooo, H_vooooo, H_vooovooo);
#endif
    //std::cout<<"H_oooo"<<std::endl<<H_oooo<<std::endl;
    //std::cout<<"H_oovooo"<<std::endl<<H_oovooo<<std::endl;
    //std::cout<<"H_vooooo"<<std::endl<<H_vooooo<<std::endl;
    //std::cout<<"H_vooovooo"<<std::endl<<H_vooovooo<<std::endl;

    TA::foreach_inplace( H_oooo, [this, &ldH, HbarSS ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          if (x[0] == x[1])
            continue;
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
              if (x[2] == x[3])
                continue;
              MatsT v = tile[x];

              size_t i = x[0], j = x[1], k = x[2], l = x[3];
              if (i<j && k<l){
                auto idx_sgn = toCompoundSS(i,j,k,l,ldH);
//std::cout<<"SS "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                if (isInBound(idx_sgn.first)) {
                  HbarSS[idx_sgn.first] = idx_sgn.second * v;
                }
              }

            }
        }
    });
    TA::foreach_inplace( H_oovooo, [this, &ldH, HbarSD ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          if (x[0] == x[1])
            continue;
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) 
              for(x[4] = lobound[4]; x[4] != upbound[4]; ++x[4]) 
                for(x[5] = lobound[5]; x[5] != upbound[5]; ++x[5]) {
                  if (x[3] == x[4] or x[4] == x[5] or x[3] == x[5])
                    continue;
                  MatsT v = tile[x];
    
                  //size_t i = x[0], j = x[1], k = x[2], l = x[3], m = x[4], b = x[5];
                  size_t i = x[0], j = x[1], b = x[2], k = x[3], l = x[4], m = x[5];
                  if (i<j && k<l && l<m){
                    auto idx_sgn = toCompoundSD(i,j,b,k,l,m,ldH);
//std::cout<<"SD "<<i<<" "<<j<<" "<<b<<" "<<k<<" "<<l<<" "<<m<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                    if (isInBound(idx_sgn.first)) {
                      HbarSD[idx_sgn.first] = idx_sgn.second * v;
                    }
                  }
            }
        }
    });
    TA::foreach_inplace( H_vooooo, [this, &ldH, HbarDS ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) 
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2]) 
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
              if (x[1] == x[2] or x[2] == x[3] or x[1] == x[3])
                continue;
              for(x[4] = lobound[4]; x[4] != upbound[4]; ++x[4])
                for(x[5] = lobound[5]; x[5] != upbound[5]; ++x[5]) {
                  if (x[4] == x[5])
                    continue;
                  MatsT v = tile[x];

                  //size_t i = x[0], j = x[1], k = x[2], a = x[3], l = x[4], m = x[5];
                  size_t a = x[0], i = x[1], j = x[2], k = x[3], l = x[4], m = x[5];
                  if (i<j && j<k && l<m){
                    auto idx_sgn = toCompoundDS(a,i,j,k,l,m,ldH);
//std::cout<<"DS "<<a<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                    if (isInBound(idx_sgn.first)) {
                      HbarDS[idx_sgn.first] = idx_sgn.second * v;
                    }
                  }
                }
        }
    });
    TA::foreach_inplace( H_vooovooo, [this, &ldH, HbarDD ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0, 0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) 
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2]) 
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
              if (x[1] == x[2] or x[2] == x[3] or x[1] == x[3])
                continue;
              for(x[4] = lobound[4]; x[4] != upbound[4]; ++x[4])
                for(x[5] = lobound[5]; x[5] != upbound[5]; ++x[5]) 
                  for(x[6] = lobound[6]; x[6] != upbound[6]; ++x[6])
                    for(x[7] = lobound[7]; x[7] != upbound[7]; ++x[7]) {
                      if (x[5] == x[6] or x[6] == x[7] or x[5] == x[7]) 
                        continue;
                      MatsT v = tile[x];

                      //size_t i = x[0], j = x[1], k = x[2], a = x[3], l = x[4], m = x[5], n = x[6], b = x[7];
                      size_t a = x[0], i = x[1], j = x[2], k = x[3], b = x[4], l = x[5], m = x[6], n = x[7];
                      if (i<j && j<k && l<m && m<n){
                        auto idx_sgn = toCompoundDD(a,i,j,k,b,l,m,n,ldH);
//std::cout<<"DD "<<a<<" "<<i<<" "<<j<<" "<<k<<" "<<b<<" "<<l<<" "<<m<<" "<<n<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                        if (isInBound(idx_sgn.first)) {
                          HbarDD[idx_sgn.first] = idx_sgn.second * v;
                        }
                      }
                    }
            }
    });

    TA::get_default_world().gop.fence();
    TAmanager.free("oooo", std::move(H_oooo));
    TAmanager.free("vooooo", std::move(H_vooooo));
    TAmanager.free("oovooo", std::move(H_oovooo));
    TAmanager.free("vooovooo", std::move(H_vooovooo));
    return fullMat;
  }
  //for full_diagonalization
   template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::buildHbarSS(TArray & H_oooo) const{
    TAManager &TAmanager = TAManager::get();
    // %%%% Prepare Scalars %%%%
//    dcomplex scalar_0, scalar_1, scalar_2, scalar_3, scalar_4;
    
    TArray tempPerm_oooo = TAmanager.malloc<MatsT>("oooo");  
    
    // %%%% Assign Scalars %%%%
//    scalar_0 = dot(this->fockMatrix_ta["oo"]("o0,o1"), Id_oo("o0,o1")).get();
//    TA::get_default_world().gop.fence();    
//    scalar_1 = dot(this->fockMatrix_ta["ov"]("m,a"), this->T1_("a,m")).get();
//    TA::get_default_world().gop.fence();    
//    scalar_2 = dot(this->antiSymMoints["oooo"]("o0,o1,o2,o3"), Id_oo("o0,o2"), Id_oo("o1,o3")).get();
//    TA::get_default_world().gop.fence();    
//    scalar_3 = dot(conj(this->antiSymMoints["vvoo"]("a,b,n,m")), this->T2_("a,b,n,m")).get();
//    TA::get_default_world().gop.fence();    
//    scalar_4 = dot(conj(this->antiSymMoints["vvoo"]("a,b,n,m")) * this->T1_("a,m"), this->T1_("b,n")).get();
//    TA::get_default_world().gop.fence();    


    // tempOps["oooo_0"] += 1.000000 conj(this->antiSymMoints["vvoo"]("a,b,l,m")) this->T2_("a,b,j,m") Id_oo("i,k") 
    // flops: o3v2: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    tempOps["oooo_0"] = TAmanager.malloc<MatsT>("oooo");
    tempOps["oooo_0"]("l,j,i,k") = 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T2_("a,b,j,m") * Id_oo("i,k");

    // H_oooo += -0.500000 d(i,k) <l,m||a,b> t2(a,b,j,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_0"]("l,j,i,k");

    // H_oooo += -0.500000 d(j,l) <k,m||a,b> t2(a,b,i,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_0"]("k,i,j,l");

    // H_oooo += 0.500000 d(i,l) <k,m||a,b> t2(a,b,j,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_0"]("k,j,i,l");

    // H_oooo += 0.500000 d(j,k) <l,m||a,b> t2(a,b,i,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_0"]("l,i,j,k");
    TAmanager.free("oooo", std::move(tempOps["oooo_0"]));

    // tempOps["oooo_1"] += 1.000000 conj(this->antiSymMoints["vvoo"]("a,b,l,m")) this->T1_("a,m") this->T1_("b,j") Id_oo("i,k") 
    // flops: o2v2: 1, o4v0: 2, o2v1: 1 | mem: o4v0: 2, o1v1: 1, o2v0: 1, 
    tempOps["oooo_1"] = TAmanager.malloc<MatsT>("oooo");
    tempOps["oooo_1"]("l,j,i,k") = conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T1_("a,m") * this->T1_("b,j") * Id_oo("i,k");

    // H_oooo += 1.000000 d(i,k) <l,m||a,b> t1(a,m) t1(b,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_1"]("l,j,i,k");

    // H_oooo += 1.000000 d(j,l) <k,m||a,b> t1(a,m) t1(b,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_1"]("k,i,j,l");

    // H_oooo += -1.000000 d(j,k) <l,m||a,b> t1(a,m) t1(b,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_1"]("l,i,j,k");

    // H_oooo += -1.000000 d(i,l) <k,m||a,b> t1(a,m) t1(b,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_1"]("k,j,i,l");
    TAmanager.free("oooo", std::move(tempOps["oooo_1"]));

    // tempOps["oooo_2"] += 1.000000 conj(this->antiSymMoints["vooo"]("a,j,l,m")) this->T1_("a,m") Id_oo("i,k") 
    // flops: o3v1: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    tempOps["oooo_2"] = TAmanager.malloc<MatsT>("oooo");
    tempOps["oooo_2"]("l,j,i,k") = conj(this->antiSymMoints["vooo"]("a,j,l,m")) * this->T1_("a,m") * Id_oo("i,k");

    // H_oooo += -1.000000 d(j,k) <l,m||a,i> t1(a,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_2"]("l,i,j,k");

    // H_oooo += 1.000000 d(i,k) <l,m||a,j> t1(a,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_2"]("l,j,i,k");

    // H_oooo += 1.000000 d(j,l) <k,m||a,i> t1(a,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_2"]("k,i,j,l");

    // H_oooo += -1.000000 d(i,l) <k,m||a,j> t1(a,m) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_2"]("k,j,i,l");
    TAmanager.free("oooo", std::move(tempOps["oooo_2"]));

    // tempOps["oooo_3"] += 1.000000 Id_oo("j,l") Id_oo("i,k") 
    // flops: o4v0: 2 | mem: o4v0: 2, 
//    tempOps["oooo_3"] = TAmanager.malloc<MatsT>("oooo");
//    tempOps["oooo_3"]("j,l,i,k") = Id_oo("j,l") * Id_oo("i,k");
//
//    // H_oooo += 1.000000 d(j,l) d(i,k) f(m,a) t1(a,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") += scalar_1 * tempOps["oooo_3"]("j,l,i,k");
//
//    // H_oooo += -1.000000 d(i,l) d(j,k) f(m,a) t1(a,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") -= scalar_1 * tempOps["oooo_3"]("i,l,j,k");
//
//    // H_oooo += -1.000000 d(i,l) d(j,k) f(m,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") -= scalar_0 * tempOps["oooo_3"]("i,l,j,k");
//
//    // H_oooo += 1.000000 d(j,l) d(i,k) f(m,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") += scalar_0 * tempOps["oooo_3"]("j,l,i,k");
//
//    // H_oooo += 0.500000 d(i,l) d(j,k) <n,m||a,b> t1(a,m) t1(b,n) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") += 0.500000 * scalar_4 * tempOps["oooo_3"]("i,l,j,k");
//
//    // H_oooo += -0.500000 d(j,l) d(i,k) <n,m||a,b> t1(a,m) t1(b,n) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") -= 0.500000 * scalar_4 * tempOps["oooo_3"]("j,l,i,k");
//
//    // H_oooo += 0.250000 d(j,l) d(i,k) <n,m||a,b> t2(a,b,n,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") += 0.250000 * scalar_3 * tempOps["oooo_3"]("j,l,i,k");
//
//    // H_oooo += -0.250000 d(i,l) d(j,k) <n,m||a,b> t2(a,b,n,m) 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") -= 0.250000 * scalar_3 * tempOps["oooo_3"]("i,l,j,k");
//
//    // H_oooo += 0.500000 d(i,l) d(j,k) <n,m||n,m> 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") += 0.500000 * scalar_2 * tempOps["oooo_3"]("i,l,j,k");
//
//    // H_oooo += -0.500000 d(j,l) d(i,k) <n,m||n,m> 
//    // flops: o4v0: 1, o0v0: 1 | mem: o4v0: 1, o0v0: 1, 
//    H_oooo("i,j,k,l") -= 0.500000 * scalar_2 * tempOps["oooo_3"]("j,l,i,k");
//    TAmanager.free("oooo", std::move(tempOps["oooo_3"]));

    // tempOps["oooo_4"] += 1.000000 this->fockMatrix_ta["ov"]("l,a") this->T1_("a,j") Id_oo("i,k") 
    // flops: o4v0: 2, o2v1: 1 | mem: o4v0: 2, o2v0: 1, 
    tempOps["oooo_4"] = TAmanager.malloc<MatsT>("oooo");
    tempOps["oooo_4"]("l,j,i,k") = this->fockMatrix_ta["ov"]("l,a") * this->T1_("a,j") * Id_oo("i,k");

    // H_oooo += -1.000000 d(i,k) f(l,a) t1(a,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_4"]("l,j,i,k");

    // H_oooo += 1.000000 d(j,k) f(l,a) t1(a,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_4"]("l,i,j,k");

    // H_oooo += -1.000000 d(j,l) f(k,a) t1(a,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_4"]("k,i,j,l");

    // H_oooo += 1.000000 d(i,l) f(k,a) t1(a,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_4"]("k,j,i,l");
    TAmanager.free("oooo", std::move(tempOps["oooo_4"]));

    // tempOps["oooo_5"] += 1.000000 this->fockMatrix_ta["oo"]("l,j") Id_oo("i,k") 
    // flops: o4v0: 2 | mem: o4v0: 2, 
    tempOps["oooo_5"] = TAmanager.malloc<MatsT>("oooo");
    tempOps["oooo_5"]("l,j,i,k") = this->fockMatrix_ta["oo"]("l,j") * Id_oo("i,k");

    // H_oooo += 1.000000 d(i,l) f(k,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_5"]("k,j,i,l");

    // H_oooo += -1.000000 d(j,l) f(k,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_5"]("k,i,j,l");

    // H_oooo += -1.000000 d(i,k) f(l,j) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") -= tempOps["oooo_5"]("l,j,i,k");

    // H_oooo += 1.000000 d(j,k) f(l,i) 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += tempOps["oooo_5"]("l,i,j,k");
    TAmanager.free("oooo", std::move(tempOps["oooo_5"]));

    // H_oooo += -1.000000 <k,l||a,b> t1(a,j) t1(b,i) 
    // flops: o3v2: 1, o4v1: 1, o4v0: 1 | mem: o3v1: 1, o4v0: 2, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T1_("a,j") * this->T1_("b,i");

    // H_oooo += 0.500000 <k,l||a,b> t2(a,b,i,j) 
    // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") += 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T2_("a,b,i,j");

    // H_oooo += 1.000000 P(i,j) <k,l||a,j> t1(a,i) 
    // flops: o4v1: 1, o4v0: 2 | mem: o4v0: 2, 
    tempPerm_oooo("i,j,k,l") = conj(this->antiSymMoints["vooo"]("a,j,k,l")) * this->T1_("a,i");
    H_oooo("i,j,k,l") += tempPerm_oooo("i,j,k,l");
    H_oooo("i,j,k,l") -= tempPerm_oooo("j,i,k,l");

    // H_oooo += 1.000000 <k,l||i,j> 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += this->antiSymMoints["oooo"]("k,l,i,j");

    TA::get_default_world().gop.fence();
    TAmanager.free("oooo", std::move(tempPerm_oooo));
  
  } 

   /// for full_diagonaization and Davidson routine
  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::buildDiag(MatsT * diag, const std::vector<double> &eps) const {
    TAManager &TAmanager = TAManager::get();
    size_t n_v = TAmanager.getRange(vLabel_).extent();
    size_t n_o = TAmanager.getRange(oLabel_).extent();
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));

    // Option 1, use orbital energies to provide the guess    
    for (auto i = 0; i < n_o; i++){
      for (auto j = 0; j < i; j++){
        diag[toCompoundS(i,j)] = - eps[i] - eps[j];
      }
    }
    MatsT * diag2 = diag + nO2shift_;
    for (auto a = 0; a < n_v; a++){
      for (auto i = 0; i < n_o; i++){
        for (auto j = 0; j < i; j++){
          for (auto k = 0; k < j; k++){
            diag2[toCompoundD(a,i,j,k)] = eps[a+n_o] - eps[i] - eps[j] - eps[k];
          }
        }
      }
    }
 
    // Option 2: use actual diagonals    
    //for (size_t i = 0; i < NO; ++i)
    //  for (size_t j = 0; j < i; ++j) {
    //    // - F_mi[i][i] - F_mi[j][j] + W_mnij[i][j][j][i]
    //  }
    //
    //for (size_t a = 0; a < NV; ++a)
    //  for (size_t i = 0; i < NO; ++i)
    //    for (size_t j = 0; j < i; ++j) 
    //      for (size_t k = 0; k < j; ++k) {
    //        // - F_mi[i][i] - F_mi[j][j] - F_mi[k][k] + F_ae[a][a]
    //        // + W_mnij[i][j][j][i] + W_mnij[j][k][k][j] + W_mnij[i][k][k][i]
    //        // + W_mbej[i][a][a][i] + W_mbej[j][a][a][j] + W_mbej[k][a][a][k]
    //        // - sum{c} T2[c][a][j][k] * conj(V[c][a][j][k])
    //        // - sum{c} T2[c][a][i][k] * conj(V[c][a][i][k])
    //        // - sum{c} T2[c][a][i][j] * conj(V[c][a][i][j])
    //      }

//    size_t nV = TAmanager.getRange(vLabel_).extent();
//    size_t nO = TAmanager.getRange(oLabel_).extent();
//
//    std::fill_n(diag, this->Hbar_dim, ccsd_energy);
//    MatsT * diag2 = diag + nO2shift_;
//
//    TA::foreach_inplace( F_mi, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
//      const auto& lobound = tile.range().lobound();
//      if (lobound[0] != lobound[1]) return;
//
//      const auto& upbound = tile.range().upbound();
//
//      std::size_t x[] = {0, 0};
//      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
//        x[1] = x[0];
//        MatsT v = tile[x];
//
//        size_t i = x[0], j = x[0], k = x[0];
//
//        for (size_t j = 0; j < i; ++j) {
//          diag[toCompoundS(i,j)] -= v;
//        }
//
//        for (size_t i = j + 1; i < nO; ++i) {
//          diag[toCompoundS(i,j)] -= v;
//        }
//#ifndef DEBUG_DIP
//        for (size_t j = 0; j < i; ++j)
//          for (size_t k = 0; k < j; ++k)
//            for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] -= v;
//            }
//
//        for (size_t i = j + 1; i < nO; ++i)
//          for (size_t k = 0; k < j; ++k)
//            for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] -= v;
//            }
//
//        for (size_t j = k + 1; j < nO; ++j)
//          for (size_t i = j + 1; i < nO; ++i)
//            for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] -= v;
//            }
//#endif            
//      }
//
//    });
//#ifndef DEBUG_DIP
//    TA::foreach_inplace( F_ae, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
//      const auto& lobound = tile.range().lobound();
//      if (lobound[0] != lobound[1]) return;
//
//      const auto& upbound = tile.range().upbound();
//
//      std::size_t x[] = {0, 0};
//      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
//        x[1] = x[0];
//        MatsT v = tile[x];
//
//        size_t a = x[0];
//
//        for (size_t i = 0; i < nO; ++i)
//          for (size_t j = 0; j < i; ++j)
//            for (size_t k = 0; k < j; ++k)
//              diag2[toCompoundD(a, i, j, k)] -= v;
//
//      }
//
//    });
//#endif            
//
//    TA::foreach_inplace( W_mnij, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
//
//      const auto& lobound = tile.range().lobound();
//      if (lobound[0] != lobound[3] or lobound[1] != lobound[2] or lobound[0] > lobound[1])
//      return;
//
//      const auto& upbound = tile.range().upbound();
//
//      std::size_t x[] = {0,0,0,0};
//      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
//        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0]) {
//          x[3] = x[0];
//          x[2] = x[1];
//          MatsT v = tile[x];
//
//          size_t i = x[0], j = x[1];
//
//          if (i > j) {
//            diag[toCompoundS(i, j)] += v;
//          }
//#ifndef DEBUG_DIP
//          if (i > j) {
//            for (size_t k = 0; k < j; ++k)
//              for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] += v;
//            }
//          }
//          size_t k = x[1]; j = x[0];  i = x[0];
//          if (j > k) {
//            for (size_t i = j + 1; i < nO; ++i)
//              for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] += v;
//            }
//          }
//          if (i > k) {
//            for (size_t j = k + 1; j < i; ++j)
//              for (size_t a = 0; a < nV; ++a) {
//              diag2[toCompoundD(a, i, j, k)] += v;
//            }
//          }
//#endif
//        }
//      }
//
//    });
//
//
//#ifndef DEBUG_DIP
//    //        // + W_mbej[i][a][a][i] + W_mbej[j][a][a][j] + W_mbej[k][a][a][k]
//	TA::foreach_inplace( W_mbej, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
//
//      const auto& lobound = tile.range().lobound();
//      if (lobound[0] != lobound[3] or lobound[1] != lobound[2])
//      return;
//
//      const auto& upbound = tile.range().upbound();
//
//      std::size_t x[] = {0,0,0,0};
//      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
//        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
//          x[3] = x[0];
//          x[2] = x[1];
//          MatsT v = tile[x];
//
//          size_t i = x[0], j = x[0], k = x[0], a = x[1];
//
//          for (size_t j = 0; j < i; ++j) 
//            for (size_t k = 0; k < j; ++k)
//            diag2[toCompoundD(a, i, j, k)] += v;
//          for (size_t i = j + 1; i < nO; ++i)
//          	for (size_t k = 0; k < j; ++k)
//            diag2[toCompoundD(a, i, j, k)] += v;
//          for (size_t i = j + 1; i < nO; ++i)
//            for (size_t j = k + 1; j < i; ++j)
//            diag2[toCompoundD(a, i, j, k)] += v;
//        }
//      }
//
//    });
//
//    TA::foreach_inplace(T2_, this->antiSymMoints.at("vvoo"),
//                        [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& T2tile, const TA::Tensor<MatsT>& Vtile) {
//      const auto& lobound = T2tile.range().lobound();
//      const auto& upbound = T2tile.range().upbound();
//
//      std::vector<std::size_t> x{0, 0, 0, 0};
//      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
//        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
//          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
//            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
//              MatsT v = SmartConj(T2tile[x]) * Vtile[x];
//
//              size_t a = x[0];
//
//			  size_t j = x[2], k = x[3];
//              if (j > k) {
//              for (size_t i = j + 1; i < nO; ++i)
//                for (size_t b = a + 1; b < nV; ++b)
//                  diag2[toCompoundD(a, i, j, k)] -= v;
//              }
//
//              size_t i = x[2]; k = x[3];
//              if (i > k) {
//              for (size_t j = k + 1; j < i; ++j)
//                for (size_t b = a + 1; b < nV; ++b)
//                  diag2[toCompoundD(a, i, j, k)] -= v;
//              }
//
//              i = x[2]; j = x[3];
//              if (i > j) {
//              for (size_t k = 0; k < j; ++k)
//                for (size_t b = a + 1; b < nV; ++b)
//                  diag2[toCompoundD(a, i, j, k)] -= v;
//              }
//
//
//            }
//    });
//#endif
//    TA::get_default_world().gop.fence();
//    TA::get_default_world().gop.template reduce(diag, this->Hbar_dim, std::plus<MatsT>());
//

  }

//  template <typename MatsT>
//  void EOMDIP_3h1p<MatsT>::formF_ae() {
//    F_ae("a,e") -= 0.5 * this->T1_("a,m") * F_me("m,e");
//  }
//
//  template <typename MatsT>
//  void EOMDIP_3h1p<MatsT>::formF_mi() {
//    F_mi("m,i") += 0.5 * this->T1_("e,i") * F_me("m,e");
//  }
//
//  template <typename MatsT>
//  void EOMDIP_3h1p<MatsT>::formW_mnij() {
//    W_mnij("m,n,i,j") += 0.25 * tau("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
//  }
//
//  template <typename MatsT>
//  void EOMDIP_3h1p<MatsT>::formW_mbej() {
//    W_mbej("m,b,e,j") -= 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
//  }

  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::runLambda(){} 


  template <typename MatsT>
  typename Davidson<MatsT>::VecsGen_t EOMDIP_3h1p<MatsT>::EmptyDavidsonVectorBuilder(){
      // Algorithm with implicit Hbar matrix
      typename Davidson<MatsT>::VecsGen_t vecsGenEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSet<MatsT>>(this->tensor_builder_, nVec, this->savFile_);
        }; // implicit vecsGenerator

        return vecsGenEOM;
      }
      fullMat = nullptr;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or 
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        std::cout << "  *** Start building the full matrix for explicit diagonalization ***" << std::endl;

        auto beginBuildHbar = tick();
        fullMat = std::make_shared<cqmatrix::Matrix<MatsT>>(buildHbar(false));
        std::cout << "    * Build Hbar spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginBuildHbar) << " s." << std::endl;
      }

      // Algorithm for debug, comparing implicit and explicit
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        typename Davidson<MatsT>::VecsGen_t vecsGenRaw = [/*&Hbar_dim,*/ this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<RawVectors<MatsT>>(
              MPI_COMM_WORLD, this->Hbar_dim, nVec
              );
        };

        typename Davidson<MatsT>::VecsGen_t vecsGenDebug =
            [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSetDebug<MatsT>>(
              this->tensor_builder_, nVec, this->savFile_, MPI_COMM_WORLD
              );
        };
        return vecsGenDebug;
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return vecsGenEOM;
      return vecsGenEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMDIP_3h1p<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
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
            //buildSigma(Vi.get_tensor("OneBody"), Vi.get_tensor("TwoBody"), AVi.get_tensor("OneBody"), AVi.get_tensor("TwoBody"), eigenVecType);
            buildSigma(Vi, AVi, eigenVecType);
            TA::get_default_world().gop.fence();
            AVi.enforceSymmetry();
          }

//V.print(std::cout, "R-imp", 0, nVec);
//AV.print(std::cout, "H*R-imp", 0, nVec);
        }; // implicit sigmaBuilder
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        this->funcRaw = [this, &eigenVecType]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {
          ROOT_ONLY(MPI_COMM_WORLD);

          size_t N = this->Hbar_dim;
          auto V_ptr = tryGetRawVectorsPointer(V);
          auto AV_ptr = tryGetRawVectorsPointer(AV);

          switch(eigenVecType) {
            case EOMCCEigenVecType::RIGHT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
                         N,nVec,N,MatsT(1.),fullMat->pointer(),N,V_ptr,N,MatsT(0.),AV_ptr,N);
              break;
            case EOMCCEigenVecType::LEFT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                         nVec,N,N,MatsT(1.),V_ptr,N,fullMat->pointer(),N,MatsT(0.),AV_ptr,nVec);
              IMatCopy('T', nVec, N, 1.0, AV_ptr,nVec, N);
              break;
          }

//V.print(std::cout, "R-exp", 0, nVec);
//AV.print(std::cout, "H*R-exp", 0, nVec);
        }; // explicit sigmaBuilder
      }

      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        this->funcDebug = [this]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSetDebug<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<MatsT> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<MatsT> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<MatsT> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<MatsT> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::funcDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          this->funcEOM(nVec, V_EOM, AV_EOM);
          this->funcRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::funcDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };
        return this->funcDebug;

      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) return this->funcEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return this->funcRaw;

      return this->funcEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMDIP_3h1p<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){
      
      double PCsmall = this->eomSettings.davidson_preCond_small;

      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT or
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        this->PCEOM = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<MatsT> &V,
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
            if constexpr (std::is_same_v<MatsT, double>) {
              curEigI = curEig[iVec].real();
            } else {
              curEigI = curEig[iVec];
            }

            TA::foreach_inplace(curB.get_tensor("OneBody"), [iVec, curEigI, eomDiag, this, PCsmall](TA::Tensor<MatsT> &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              MatsT denom = 0.0;
              std::vector<std::size_t> x{0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                  if (x[0] == x[1])
                    continue;
                  denom = curEigI - eomDiag[toCompoundS(x[0], x[1])];
                  if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                }
            });
            TA::get_default_world().gop.fence();

            size_t nO = this->intermediates_.nOcc;
            MatsT *diagD = eomDiag + nO * (nO - 1) / 2;

            TA::foreach_inplace(curB.get_tensor("TwoBody"), [iVec, curEigI, diagD, this, PCsmall](TA::Tensor<MatsT> &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              MatsT denom = 0.0;
              std::vector<std::size_t> x{0, 0, 0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
                size_t a = x[0];
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
                  for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                    for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                      if (x[2] == x[3] || x[1] == x[2] || x[1] == x[3])
                        continue;
                      size_t i = x[1], j = x[2], k = x[3];
                      //signD(a,i,j,k);
                      denom = curEigI - diagD[toCompoundD(a, i, j, k)];
                      if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                    }
                }
            });
            TA::get_default_world().gop.fence();

            curB.enforceSymmetry();
          }
        }; // implicit preConditioner

      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT or 
          this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {

        this->PCRaw = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          ROOT_ONLY(MPI_COMM_WORLD);

          //            prettyPrintSmart(std::cout, "eomDiag", eomDiag, Hbar_dim, 1, Hbar_dim);

          MatsT denom = 0.0, nom = 0.0;
          const MatsT *Vptr = tryGetRawVectorsPointer(V);
          MatsT *AVptr = tryGetRawVectorsPointer(AV);

          // Scale by inverse diagonals
          for (size_t i = 0; i < nVec; i++) {
            MatsT curEigI = 0.0;
            if constexpr (std::is_same_v<MatsT, double>) {
              curEigI = curEig[i].real();
            } else {
              curEigI = curEig[i];
            }
            for(auto k = 0ul; k < this->Hbar_dim; k++ ) {
              nom = Vptr[k];
              denom = curEigI - eomDiag[k];
              if (std::abs(denom) >= PCsmall)
                AVptr[k] = nom / denom;
              else
                AVptr[k] = nom;
            }
            Vptr += this->Hbar_dim;
            AVptr += this->Hbar_dim;
          }

        }; // explicit preConditioner
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        this->PCDebug = [this]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSetDebug<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSetDebug<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<MatsT> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<MatsT> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<MatsT> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<MatsT> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::PCDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          this->PCEOM(nVec, V_EOM, AV_EOM);
          this->PCRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::PCDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };
        return this->PCDebug;
      }
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) return this->PCEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT) return this->PCRaw;

      return this->PCEOM;
  }
}
