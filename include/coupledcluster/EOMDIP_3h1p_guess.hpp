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
//#include <coupledcluster/MBExpansion.hpp>
//#include <itersolver/davidson.hpp>
//#include <itersolver.hpp>

namespace ChronusQ{
  //for full_diagonalization
  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::fillGuess(MatsT *guess_vec, size_t n_vec) const{

    TAManager &TAmanager = TAManager::get();
    // %%%% Prepare Scalars %%%%

    TArray tempPerm_oooo = TAmanager.malloc<MatsT>("oooo");

    /// ****** p†q ****** ///


    TArray Id_oo = TAmanager.malloc<MatsT>("oo");
    TArray H_oooo = TAmanager.malloc<MatsT>("oooo");
    
    TA::foreach_inplace(Id_oo, [&](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          if(x[0]==x[1]) 
            tile[x] = 1.0;
          else 
            tile[x] = 0.0;
    });
    TA::get_default_world().gop.fence();    

    // H_oooo += -1.000000 d(i,k) f(l,j) 
    // flops: o4v0: 2 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") -= Id_oo("i,k") * this->fockMatrix_ta["oo"]("l,j");

    // H_oooo += 1.000000 d(i,l) f(k,j) 
    // flops: o4v0: 2 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") += Id_oo("i,l") * this->fockMatrix_ta["oo"]("k,j");

    // H_oooo += 1.000000 d(j,k) f(l,i) 
    // flops: o4v0: 2 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") += Id_oo("j,k") * this->fockMatrix_ta["oo"]("l,i");

    // H_oooo += -1.000000 d(j,l) f(k,i) 
    // flops: o4v0: 2 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") -= Id_oo("j,l") * this->fockMatrix_ta["oo"]("k,i");

    // H_oooo += -1.000000 d(i,k) f(l,a) t1(a,j) 
    // flops: o4v0: 2, o2v1: 1 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= this->fockMatrix_ta["ov"]("l,a") * this->T1_("a,j") * Id_oo("i,k");

    // H_oooo += 1.000000 d(i,l) f(k,a) t1(a,j) 
    // flops: o4v0: 2, o2v1: 1 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += this->fockMatrix_ta["ov"]("k,a") * this->T1_("a,j") * Id_oo("i,l");

    // H_oooo += 1.000000 d(j,k) f(l,a) t1(a,i) 
    // flops: o4v0: 2, o2v1: 1 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += this->fockMatrix_ta["ov"]("l,a") * this->T1_("a,i") * Id_oo("j,k");

    // H_oooo += -1.000000 d(j,l) f(k,a) t1(a,i) 
    // flops: o4v0: 2, o2v1: 1 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= this->fockMatrix_ta["ov"]("k,a") * this->T1_("a,i") * Id_oo("j,l");

    // H_oooo += 1.000000 <k,l||i,j> 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += this->antiSymMoints["oooo"]("k,l,i,j");

    // H_oooo += 1.000000 d(i,k) <l,m||a,j> t1(a,m) 
    // flops: o3v1: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += conj(this->antiSymMoints["vooo"]("a,j,l,m")) * this->T1_("a,m") * Id_oo("i,k");

    // H_oooo += -1.000000 d(i,l) <k,m||a,j> t1(a,m) 
    // flops: o3v1: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vooo"]("a,j,k,m")) * this->T1_("a,m") * Id_oo("i,l");

    // H_oooo += 1.000000 P(i,j) <k,l||a,j> t1(a,i) 
    // flops: o4v1: 1, o4v0: 2 | mem: o4v0: 2, 
    tempPerm_oooo("i,j,k,l") = conj(this->antiSymMoints["vooo"]("a,j,k,l")) * this->T1_("a,i");
    H_oooo("i,j,k,l") += tempPerm_oooo("i,j,k,l");
    H_oooo("i,j,k,l") -= tempPerm_oooo("j,i,k,l");

    // H_oooo += -1.000000 d(j,k) <l,m||a,i> t1(a,m) 
    // flops: o3v1: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vooo"]("a,i,l,m")) * this->T1_("a,m") * Id_oo("j,k");

    // H_oooo += 1.000000 d(j,l) <k,m||a,i> t1(a,m) 
    // flops: o3v1: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += conj(this->antiSymMoints["vooo"]("a,i,k,m")) * this->T1_("a,m") * Id_oo("j,l");

    // H_oooo += -0.500000 d(i,k) <l,m||a,b> t2(a,b,j,m) 
    // flops: o3v2: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T2_("a,b,j,m") * Id_oo("i,k");

    // H_oooo += 0.500000 d(i,l) <k,m||a,b> t2(a,b,j,m) 
    // flops: o3v2: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,k,m")) * this->T2_("a,b,j,m") * Id_oo("i,l");

    // H_oooo += 0.500000 d(j,k) <l,m||a,b> t2(a,b,i,m) 
    // flops: o3v2: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") += 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T2_("a,b,i,m") * Id_oo("j,k");

    // H_oooo += -0.500000 d(j,l) <k,m||a,b> t2(a,b,i,m) 
    // flops: o3v2: 1, o4v0: 2 | mem: o4v0: 2, o2v0: 1, 
    H_oooo("i,j,k,l") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,k,m")) * this->T2_("a,b,i,m") * Id_oo("j,l");

    // H_oooo += 0.500000 <k,l||a,b> t2(a,b,i,j) 
    // flops: o4v2: 1, o4v0: 1 | mem: o4v0: 2, 
    H_oooo("i,j,k,l") += 0.500000 * conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T2_("a,b,i,j");

    // H_oooo += 1.000000 d(i,k) <l,m||a,b> t1(a,m) t1(b,j) 
    // flops: o2v2: 1, o4v0: 2, o2v1: 1 | mem: o4v0: 2, o1v1: 1, o2v0: 1, 
    H_oooo("i,j,k,l") += conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T1_("a,m") * this->T1_("b,j") * Id_oo("i,k");

    // H_oooo += -1.000000 d(i,l) <k,m||a,b> t1(a,m) t1(b,j) 
    // flops: o2v2: 1, o4v0: 2, o2v1: 1 | mem: o4v0: 2, o1v1: 1, o2v0: 1, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vvoo"]("a,b,k,m")) * this->T1_("a,m") * this->T1_("b,j") * Id_oo("i,l");

    // H_oooo += -1.000000 d(j,k) <l,m||a,b> t1(a,m) t1(b,i) 
    // flops: o2v2: 1, o4v0: 2, o2v1: 1 | mem: o4v0: 2, o1v1: 1, o2v0: 1, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T1_("a,m") * this->T1_("b,i") * Id_oo("j,k");

    // H_oooo += 1.000000 d(j,l) <k,m||a,b> t1(a,m) t1(b,i) 
    // flops: o2v2: 1, o4v0: 2, o2v1: 1 | mem: o4v0: 2, o1v1: 1, o2v0: 1, 
    H_oooo("i,j,k,l") += conj(this->antiSymMoints["vvoo"]("a,b,k,m")) * this->T1_("a,m") * this->T1_("b,i") * Id_oo("j,l");

    // H_oooo += -1.000000 <k,l||a,b> t1(a,j) t1(b,i) 
    // flops: o3v2: 1, o4v1: 1, o4v0: 1 | mem: o3v1: 1, o4v0: 2, 
    H_oooo("i,j,k,l") -= conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T1_("a,j") * this->T1_("b,i");


    TAmanager.free("oooo", std::move(tempPerm_oooo));
    cqmatrix::Matrix<MatsT> fullMat(nO2shift_);
    MatsT * raw = fullMat.pointer();
    TA::foreach_inplace(H_oooo, [raw, this](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              if(x[0] >= x[1] or x[2] >= x[3]) continue;
              size_t i = x[0];
              size_t j = x[1];
              size_t k = x[2];
              size_t l = x[3];
              size_t ij = i + j*(j-1)/2;
              size_t kl = k + l*(l-1)/2;
              size_t idx = ij + this->nO2shift_ * kl;
              raw[idx] = tile[x];
            }
    });
    TA::get_default_world().gop.fence();    

    MatsT *raw_copy = CQMemManager::get().malloc<MatsT>(nO2shift_*nO2shift_);
    std::copy_n(raw, nO2shift_*nO2shift_, raw_copy);
    std::fill_n(raw, nO2shift_*nO2shift_, MatsT(0.0));
    MPIAllReduce(raw_copy, nO2shift_*nO2shift_, raw, MPI_COMM_WORLD);
    CQMemManager::get().free(raw_copy);


    dcomplex * theta = CQMemManager::get().malloc<dcomplex>(nO2shift_);
    MatsT * VR    = CQMemManager::get().malloc<MatsT>(nO2shift_ * nO2shift_);
    MatsT * dummy = nullptr;
    if (MPIRank() == 0) GeneralEigen('N', 'V', nO2shift_, fullMat.pointer(), nO2shift_, theta, dummy, 1, VR, nO2shift_);
    std::copy_n(VR, nO2shift_ * n_vec, guess_vec);
    CQMemManager::get().free(theta);
    CQMemManager::get().free(VR);

  }

}

