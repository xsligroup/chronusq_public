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
#include <coupledcluster/MBExpansion.hpp>
#include <itersolver/davidson.hpp>
#include <itersolver.hpp>

namespace ChronusQ{
  //for full_diagonalization
  template <typename MatsT>
  void EOMDIP_3h1p<MatsT>::buildHbarTA(TArray & H_oooo, TArray & H_oovooo, TArray & H_vooooo, TArray & H_vooovooo) const{

    TAManager &TAmanager = TAManager::get();
    // %%%% Prepare Scalars %%%%
    MatsT scalar_0, scalar_1, scalar_2, scalar_3, scalar_4;
    
    // %%%% Assign Scalars %%%%
    TArray  Id_oooo = TAmanager.malloc<MatsT>("oooo");
    Id_oooo("i,j,k,l") = Id_oo("i,k") * Id_oo("j,l");
  
    scalar_0 = dot(this->fockMatrix_ta["oo"]("o0,o1"), Id_oo("o0,o1")).get();
    TA::get_default_world().gop.fence();    
    scalar_1 = dot(this->fockMatrix_ta["ov"]("m,a"), this->T1_("a,m")).get();
    TA::get_default_world().gop.fence();    
    scalar_2 = dot(this->antiSymMoints["oooo"]("o0,o1,o2,o3"), Id_oooo("o0,o1,o2,o3")).get();
    TA::get_default_world().gop.fence();    
    scalar_3 = dot(conj(this->antiSymMoints["vvoo"]("a,b,n,m")), this->T2_("a,b,n,m")).get();
    TA::get_default_world().gop.fence();    
    scalar_4 = dot(conj(this->antiSymMoints["vvoo"]("a,b,n,m")) * this->T1_("a,m"), this->T1_("b,n")).get();
    TA::get_default_world().gop.fence();    

    TArray Id_vv = TAmanager.malloc<MatsT>("vv");
    
    TA::foreach_inplace(Id_vv, [&](TA::Tensor<MatsT> &tile){
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
    
    TArray tempPerm_vooooo = TAmanager.malloc<MatsT>("vooooo");
    TArray tempPerm_oooo = TAmanager.malloc<MatsT>("oooo");
    TArray tempPerm_vvoooooo = TAmanager.malloc<MatsT>("vvoooooo");

    /// ****** p†q ****** ///



    // H_oooo += 1.000000 d(j,l) d(i,k) f(m,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") += scalar_0 * Id_oo("j,l") * Id_oo("i,k");

    // H_oooo += -1.000000 d(i,l) d(j,k) f(m,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") -= scalar_0 * Id_oo("i,l") * Id_oo("j,k");

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

    // H_oooo += 1.000000 d(j,l) d(i,k) f(m,a) t1(a,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") += scalar_1 * Id_oo("j,l") * Id_oo("i,k");

    // H_oooo += -1.000000 d(i,l) d(j,k) f(m,a) t1(a,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") -= scalar_1 * Id_oo("i,l") * Id_oo("j,k");

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

    // H_oooo += -0.500000 d(j,l) d(i,k) <n,m||n,m> 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") -= 0.500000 * scalar_2 * Id_oo("j,l") * Id_oo("i,k");

    // H_oooo += 0.500000 d(i,l) d(j,k) <n,m||n,m> 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") += 0.500000 * scalar_2 * Id_oo("i,l") * Id_oo("j,k");

//dcomplex temp,temp2;    
//TA::get_default_world().gop.fence();    
//TA::foreach_inplace(H_oooo, [&](TA::Tensor<MatsT> &tile){
//  std::vector<std::size_t> x{2, 3, 0, 1};
//  temp = tile[x];
//  std::vector<std::size_t> x2{3, 2, 0, 1};
//  temp2 = tile[x2];
//});
//TA::get_default_world().gop.fence();    
//std::cout<<"temp "<<temp<<" and "<< temp2 <<std::flush;
//TA::get_default_world().gop.fence();    
    
    // H_oooo += 1.000000 <k,l||i,j> 
    // flops: o4v0: 1 | mem: o4v0: 1, 
    H_oooo("i,j,k,l") += this->antiSymMoints["oooo"]("k,l,i,j");
//TA::get_default_world().gop.fence();    
//TA::foreach_inplace(H_oooo, [&](TA::Tensor<MatsT> &tile){
//  std::vector<std::size_t> x{2, 3, 0, 1};
//  temp = tile[x];
//  std::vector<std::size_t> x2{3, 2, 0, 1};
//  temp2 = tile[x2];
//});
//TA::get_default_world().gop.fence();    
//std::cout<<" vs "<<temp<<" and "<< temp2<< std::endl<<std::flush;
//TA::get_default_world().gop.fence();    
    

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

    // H_oooo += 0.250000 d(j,l) d(i,k) <n,m||a,b> t2(a,b,n,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") += 0.250000 * scalar_3 * Id_oo("j,l") * Id_oo("i,k");

    // H_oooo += -0.250000 d(i,l) d(j,k) <n,m||a,b> t2(a,b,n,m) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") -= 0.250000 * scalar_3 * Id_oo("i,l") * Id_oo("j,k");

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

    // H_oooo += -0.500000 d(j,l) d(i,k) <n,m||a,b> t1(a,m) t1(b,n) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") -= 0.500000 * scalar_4 * Id_oo("j,l") * Id_oo("i,k");

    // H_oooo += 0.500000 d(i,l) d(j,k) <n,m||a,b> t1(a,m) t1(b,n) 
    // flops: o4v0: 2, o0v0: 1 | mem: o4v0: 2, o0v0: 1, 
    H_oooo("i,j,k,l") += 0.500000 * scalar_4 * Id_oo("i,l") * Id_oo("j,k");

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




/// ****** H(a,i,j,k,b,l,m,n) ****** ///




    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,m) d(i,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_0 * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,m) d(j,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_0 * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,m) d(i,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_0 * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,m) d(k,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_0 * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,m) d(j,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_0 * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,m) d(k,l) f(o,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_0 * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,m) d(i,l) f(n,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,m") * Id_oo("i,l") * this->fockMatrix_ta["oo"]("n,k") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,m) d(j,l) f(n,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,m") * Id_oo("j,l") * this->fockMatrix_ta["oo"]("n,k") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,l) f(m,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,n") * Id_oo("i,l") * this->fockMatrix_ta["oo"]("m,k") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(i,m) f(l,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,n") * Id_oo("i,m") * this->fockMatrix_ta["oo"]("l,k") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,l) f(m,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,n") * Id_oo("j,l") * this->fockMatrix_ta["oo"]("m,k") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(j,m) f(l,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,n") * Id_oo("j,m") * this->fockMatrix_ta["oo"]("l,k") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,m) d(i,l) f(n,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,m") * Id_oo("i,l") * this->fockMatrix_ta["oo"]("n,j") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,m) d(k,l) f(n,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,m") * Id_oo("k,l") * this->fockMatrix_ta["oo"]("n,j") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,l) f(m,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,n") * Id_oo("i,l") * this->fockMatrix_ta["oo"]("m,j") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(i,m) f(l,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,n") * Id_oo("i,m") * this->fockMatrix_ta["oo"]("l,j") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,l) f(m,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,n") * Id_oo("k,l") * this->fockMatrix_ta["oo"]("m,j") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(k,m) f(l,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,n") * Id_oo("k,m") * this->fockMatrix_ta["oo"]("l,j") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,m) d(j,l) f(n,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,m") * Id_oo("j,l") * this->fockMatrix_ta["oo"]("n,i") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,m) d(k,l) f(n,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,m") * Id_oo("k,l") * this->fockMatrix_ta["oo"]("n,i") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,l) f(m,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,n") * Id_oo("j,l") * this->fockMatrix_ta["oo"]("m,i") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(j,m) f(l,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,n") * Id_oo("j,m") * this->fockMatrix_ta["oo"]("l,i") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,l) f(m,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,n") * Id_oo("k,l") * this->fockMatrix_ta["oo"]("m,i") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(k,m) f(l,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,n") * Id_oo("k,m") * this->fockMatrix_ta["oo"]("l,i") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(k,n) d(j,m) d(i,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += -1.000000 d(k,n) d(i,m) d(j,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += -1.000000 d(j,n) d(k,m) d(i,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += 1.000000 d(j,n) d(i,m) d(k,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += 1.000000 d(i,n) d(k,m) d(j,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += -1.000000 d(i,n) d(j,m) d(k,l) f(a,b) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * this->fockMatrix_ta["vv"]("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,m) d(i,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_1 * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,m) d(j,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_1 * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,m) d(i,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_1 * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,m) d(k,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_1 * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,m) d(j,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += scalar_1 * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,m) d(k,l) f(o,c) t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= scalar_1 * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,m) d(i,l) f(n,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,k") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,m) d(j,l) f(n,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,k") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,l) f(m,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,k") * Id_oo("j,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(i,m) f(l,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,k") * Id_oo("j,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,l) f(m,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,k") * Id_oo("i,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(j,m) f(l,c) t1(c,k) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,k") * Id_oo("i,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,m) d(i,l) f(n,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,j") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,m) d(k,l) f(n,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,j") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,l) f(m,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,j") * Id_oo("k,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(i,m) f(l,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,j") * Id_oo("k,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,l) f(m,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,j") * Id_oo("i,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(k,m) f(l,c) t1(c,j) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,j") * Id_oo("i,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,m) d(j,l) f(n,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,i") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,m) d(k,l) f(n,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("n,c") * this->T1_("c,i") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,l) f(m,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,i") * Id_oo("k,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(j,m) f(l,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,i") * Id_oo("k,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,l) f(m,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,i") * Id_oo("j,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(k,m) f(l,c) t1(c,i) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("l,c") * this->T1_("c,i") * Id_oo("j,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(k,n) d(j,m) d(i,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(k,n) d(i,m) d(j,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(j,n) d(k,m) d(i,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(j,n) d(i,m) d(k,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(i,n) d(k,m) d(j,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(i,n) d(j,m) d(k,l) f(o,b) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 1, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += this->fockMatrix_ta["ov"]("o,b") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += -0.500000 d(a,b) d(k,n) d(j,m) d(i,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_2 * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,n) d(i,m) d(j,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_2 * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,n) d(k,m) d(i,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_2 * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,n) d(i,m) d(k,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_2 * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,n) d(k,m) d(j,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_2 * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,n) d(j,m) d(k,l) <t,o||t,o> 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_2 * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,l) <m,n||j,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_vv("a,b") * Id_oo("i,l") * this->antiSymMoints["oooo"]("m,n,j,k");

    // H_vooovooo += -1.000000 d(a,b) d(i,m) <l,n||j,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_vv("a,b") * Id_oo("i,m") * this->antiSymMoints["oooo"]("l,n,j,k");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) <l,m||j,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_vv("a,b") * Id_oo("i,n") * this->antiSymMoints["oooo"]("l,m,j,k");

    // H_vooovooo += -1.000000 d(a,b) d(j,l) <m,n||i,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_vv("a,b") * Id_oo("j,l") * this->antiSymMoints["oooo"]("m,n,i,k");

    // H_vooovooo += 1.000000 d(a,b) d(j,m) <l,n||i,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_vv("a,b") * Id_oo("j,m") * this->antiSymMoints["oooo"]("l,n,i,k");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) <l,m||i,k> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_vv("a,b") * Id_oo("j,n") * this->antiSymMoints["oooo"]("l,m,i,k");

    // H_vooovooo += 1.000000 d(a,b) d(k,l) <m,n||i,j> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_vv("a,b") * Id_oo("k,l") * this->antiSymMoints["oooo"]("m,n,i,j");

    // H_vooovooo += -1.000000 d(a,b) d(k,m) <l,n||i,j> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_vv("a,b") * Id_oo("k,m") * this->antiSymMoints["oooo"]("l,n,i,j");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) <l,m||i,j> 
    // flops: o6v2: 2, o2v2: 1 | mem: o6v2: 2, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_vv("a,b") * Id_oo("k,n") * this->antiSymMoints["oooo"]("l,m,i,j");

    // H_vooovooo += 1.000000 d(j,m) d(i,l) <n,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,m") * Id_oo("i,l") * this->antiSymMoints["vovo"]("a,n,b,k");

    // H_vooovooo += -1.000000 d(i,m) d(j,l) <n,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,m") * Id_oo("j,l") * this->antiSymMoints["vovo"]("a,n,b,k");

    // H_vooovooo += -1.000000 d(j,n) d(i,l) <m,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,n") * Id_oo("i,l") * this->antiSymMoints["vovo"]("a,m,b,k");

    // H_vooovooo += 1.000000 d(j,n) d(i,m) <l,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,n") * Id_oo("i,m") * this->antiSymMoints["vovo"]("a,l,b,k");

    // H_vooovooo += 1.000000 d(i,n) d(j,l) <m,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,n") * Id_oo("j,l") * this->antiSymMoints["vovo"]("a,m,b,k");

    // H_vooovooo += -1.000000 d(i,n) d(j,m) <l,a||b,k> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,n") * Id_oo("j,m") * this->antiSymMoints["vovo"]("a,l,b,k");

    // H_vooovooo += -1.000000 d(k,m) d(i,l) <n,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,m") * Id_oo("i,l") * this->antiSymMoints["vovo"]("a,n,b,j");

    // H_vooovooo += 1.000000 d(i,m) d(k,l) <n,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,m") * Id_oo("k,l") * this->antiSymMoints["vovo"]("a,n,b,j");

    // H_vooovooo += 1.000000 d(k,n) d(i,l) <m,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,n") * Id_oo("i,l") * this->antiSymMoints["vovo"]("a,m,b,j");

    // H_vooovooo += -1.000000 d(k,n) d(i,m) <l,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,n") * Id_oo("i,m") * this->antiSymMoints["vovo"]("a,l,b,j");

    // H_vooovooo += -1.000000 d(i,n) d(k,l) <m,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("i,n") * Id_oo("k,l") * this->antiSymMoints["vovo"]("a,m,b,j");

    // H_vooovooo += 1.000000 d(i,n) d(k,m) <l,a||b,j> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("i,n") * Id_oo("k,m") * this->antiSymMoints["vovo"]("a,l,b,j");

    // H_vooovooo += 1.000000 d(k,m) d(j,l) <n,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,m") * Id_oo("j,l") * this->antiSymMoints["vovo"]("a,n,b,i");

    // H_vooovooo += -1.000000 d(j,m) d(k,l) <n,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,m") * Id_oo("k,l") * this->antiSymMoints["vovo"]("a,n,b,i");

    // H_vooovooo += -1.000000 d(k,n) d(j,l) <m,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("k,n") * Id_oo("j,l") * this->antiSymMoints["vovo"]("a,m,b,i");

    // H_vooovooo += 1.000000 d(k,n) d(j,m) <l,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("k,n") * Id_oo("j,m") * this->antiSymMoints["vovo"]("a,l,b,i");

    // H_vooovooo += 1.000000 d(j,n) d(k,l) <m,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= Id_oo("j,n") * Id_oo("k,l") * this->antiSymMoints["vovo"]("a,m,b,i");

    // H_vooovooo += -1.000000 d(j,n) d(k,m) <l,a||b,i> 
    // flops: o6v2: 2, o4v0: 1 | mem: o6v2: 2, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += Id_oo("j,n") * Id_oo("k,m") * this->antiSymMoints["vovo"]("a,l,b,i");

    // H_vooovooo += 1.000000 d(a,b) d(j,m) d(i,l) <n,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,k,n,o")) * this->T1_("c,o") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,m) d(j,l) <n,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,k,n,o")) * this->T1_("c,o") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(i,l) <m,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,k,m,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,m) <l,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,k,l,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(j,l) <m,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,k,m,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,m) <l,o||c,k> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,k,l,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 P(j,k) d(a,b) d(i,l) <m,n||c,k> t1(c,j) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,m,n")) * this->T1_("c,j") * Id_oo("i,l") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,k,j,l,m,n");

    // H_vooovooo += -1.000000 P(j,k) d(a,b) d(i,m) <l,n||c,k> t1(c,j) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,l,n")) * this->T1_("c,j") * Id_oo("i,m") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,k,j,l,m,n");

    // H_vooovooo += 1.000000 P(j,k) d(a,b) d(i,n) <l,m||c,k> t1(c,j) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,l,m")) * this->T1_("c,j") * Id_oo("i,n") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,k,j,l,m,n");

    // H_vooovooo += -1.000000 P(i,k) d(a,b) d(j,l) <m,n||c,k> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,m,n")) * this->T1_("c,i") * Id_oo("j,l") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,k,j,i,l,m,n");

    // H_vooovooo += 1.000000 P(i,k) d(a,b) d(j,m) <l,n||c,k> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,l,n")) * this->T1_("c,i") * Id_oo("j,m") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,k,j,i,l,m,n");

    // H_vooovooo += -1.000000 P(i,k) d(a,b) d(j,n) <l,m||c,k> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,k,l,m")) * this->T1_("c,i") * Id_oo("j,n") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,k,j,i,l,m,n");

    // H_vooovooo += -1.000000 d(j,m) d(i,l) <n,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,k,n,o")) * this->T1_("a,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(i,m) d(j,l) <n,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,k,n,o")) * this->T1_("a,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(j,n) d(i,l) <m,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,k,m,o")) * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(j,n) d(i,m) <l,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,k,l,o")) * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,m");

    // H_vooovooo += -1.000000 d(i,n) d(j,l) <m,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,k,m,o")) * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(i,n) d(j,m) <l,o||b,k> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,k,l,o")) * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,m");

    // H_vooovooo += -1.000000 d(a,b) d(k,m) d(i,l) <n,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,j,n,o")) * this->T1_("c,o") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,m) d(k,l) <n,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,j,n,o")) * this->T1_("c,o") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(i,l) <m,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,j,m,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,m) <l,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,j,l,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(k,l) <m,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,j,m,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,m) <l,o||c,j> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,j,l,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 P(i,j) d(a,b) d(k,l) <m,n||c,j> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,j,m,n")) * this->T1_("c,i") * Id_oo("k,l") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,j,i,k,l,m,n");

    // H_vooovooo += -1.000000 P(i,j) d(a,b) d(k,m) <l,n||c,j> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,j,l,n")) * this->T1_("c,i") * Id_oo("k,m") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,j,i,k,l,m,n");

    // H_vooovooo += 1.000000 P(i,j) d(a,b) d(k,n) <l,m||c,j> t1(c,i) 
    // flops: o6v2: 3, o6v0: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    tempPerm_vvoooooo("a,b,i,j,k,l,m,n") = conj(this->antiSymMoints["vooo"]("c,j,l,m")) * this->T1_("c,i") * Id_oo("k,n") * Id_vv("a,b");
    H_vooovooo("a,i,j,k,b,l,m,n") += tempPerm_vvoooooo("a,b,i,j,k,l,m,n");
    H_vooovooo("a,i,j,k,b,l,m,n") -= tempPerm_vvoooooo("a,b,j,i,k,l,m,n");

    // H_vooovooo += 1.000000 d(k,m) d(i,l) <n,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,j,n,o")) * this->T1_("a,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(i,m) d(k,l) <n,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,j,n,o")) * this->T1_("a,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(k,n) d(i,l) <m,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,j,m,o")) * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(k,n) d(i,m) <l,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,j,l,o")) * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,m");

    // H_vooovooo += 1.000000 d(i,n) d(k,l) <m,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,j,m,o")) * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(i,n) d(k,m) <l,o||b,j> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,j,l,o")) * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,m");

    // H_vooovooo += 1.000000 d(a,b) d(k,m) d(j,l) <n,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,i,n,o")) * this->T1_("c,o") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,m) d(k,l) <n,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,i,n,o")) * this->T1_("c,o") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(j,l) <m,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,i,m,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,m) <l,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,i,l,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(k,l) <m,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("c,i,m,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,m) <l,o||c,i> t1(c,o) 
    // flops: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("c,i,l,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(k,m) d(j,l) <n,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,i,n,o")) * this->T1_("a,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(j,m) d(k,l) <n,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,i,n,o")) * this->T1_("a,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(k,n) d(j,l) <m,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,i,m,o")) * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(k,n) d(j,m) <l,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,i,l,o")) * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,m");

    // H_vooovooo += -1.000000 d(j,n) d(k,l) <m,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vooo"]("b,i,m,o")) * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(j,n) d(k,m) <l,o||b,i> t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vooo"]("b,i,l,o")) * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,m");

    // H_vooovooo += 1.000000 d(k,n) d(j,m) d(i,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(k,n) d(i,m) d(j,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,n) d(k,m) d(i,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(j,n) d(i,m) d(k,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(i,n) d(k,m) d(j,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(i,n) d(j,m) d(k,l) <o,a||c,b> t1(c,o) 
    // flops: o6v2: 2, o4v2: 1, o1v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,o")) * this->T1_("c,o") * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(j,m) d(i,l) <n,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,k") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(i,m) d(j,l) <n,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,k") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(j,n) d(i,l) <m,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,k") * Id_oo("j,n") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(j,n) d(i,m) <l,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,k") * Id_oo("j,n") * Id_oo("i,m");

    // H_vooovooo += -1.000000 d(i,n) d(j,l) <m,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,k") * Id_oo("i,n") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(i,n) d(j,m) <l,a||c,b> t1(c,k) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,k") * Id_oo("i,n") * Id_oo("j,m");

    // H_vooovooo += 1.000000 d(k,m) d(i,l) <n,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,j") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(i,m) d(k,l) <n,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,j") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(k,n) d(i,l) <m,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,j") * Id_oo("k,n") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(k,n) d(i,m) <l,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,j") * Id_oo("k,n") * Id_oo("i,m");

    // H_vooovooo += 1.000000 d(i,n) d(k,l) <m,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,j") * Id_oo("i,n") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(i,n) d(k,m) <l,a||c,b> t1(c,j) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,j") * Id_oo("i,n") * Id_oo("k,m");

    // H_vooovooo += -1.000000 d(k,m) d(j,l) <n,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,i") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(j,m) d(k,l) <n,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,n")) * this->T1_("c,i") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(k,n) d(j,l) <m,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,i") * Id_oo("k,n") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(k,n) d(j,m) <l,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,i") * Id_oo("k,n") * Id_oo("j,m");

    // H_vooovooo += -1.000000 d(j,n) d(k,l) <m,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvvo"]("c,b,a,m")) * this->T1_("c,i") * Id_oo("j,n") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(j,n) d(k,m) <l,a||c,b> t1(c,i) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvvo"]("c,b,a,l")) * this->T1_("c,i") * Id_oo("j,n") * Id_oo("k,m");

    // H_vooovooo += 0.250000 d(a,b) d(k,n) d(j,m) d(i,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.250000 * scalar_3 * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.250000 d(a,b) d(k,n) d(i,m) d(j,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.250000 * scalar_3 * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -0.250000 d(a,b) d(j,n) d(k,m) d(i,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.250000 * scalar_3 * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 0.250000 d(a,b) d(j,n) d(i,m) d(k,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.250000 * scalar_3 * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 0.250000 d(a,b) d(i,n) d(k,m) d(j,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.250000 * scalar_3 * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -0.250000 d(a,b) d(i,n) d(j,m) d(k,l) <t,o||c,d> t2(c,d,t,o) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.250000 * scalar_3 * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,m) d(i,l) <n,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,k,o") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,m) d(j,l) <n,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,k,o") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,n) d(i,l) <m,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,k,o") * Id_oo("j,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,n) d(i,m) <l,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,k,o") * Id_oo("j,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,n) d(j,l) <m,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,k,o") * Id_oo("i,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,n) d(j,m) <l,o||c,d> t2(c,d,k,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,k,o") * Id_oo("i,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,m) d(i,l) <n,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,j,o") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,m) d(k,l) <n,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,j,o") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(k,n) d(i,l) <m,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,j,o") * Id_oo("k,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,n) d(i,m) <l,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,j,o") * Id_oo("k,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,n) d(k,l) <m,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,j,o") * Id_oo("i,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,n) d(k,m) <l,o||c,d> t2(c,d,j,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,j,o") * Id_oo("i,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(k,m) d(j,l) <n,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,i,o") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,m) d(k,l) <n,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T2_("c,d,i,o") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,n) d(j,l) <m,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,i,o") * Id_oo("k,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(k,n) d(j,m) <l,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,i,o") * Id_oo("k,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,n) d(k,l) <m,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T2_("c,d,i,o") * Id_oo("j,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,n) d(k,m) <l,o||c,d> t2(c,d,i,o) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T2_("c,d,i,o") * Id_oo("j,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,l) <m,n||c,d> t2(c,d,j,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,d,j,k") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,m) <l,n||c,d> t2(c,d,j,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T2_("c,d,j,k") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,n) <l,m||c,d> t2(c,d,j,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T2_("c,d,j,k") * Id_oo("i,n") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,l) <m,n||c,d> t2(c,d,i,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,d,i,k") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,m) <l,n||c,d> t2(c,d,i,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T2_("c,d,i,k") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,n) <l,m||c,d> t2(c,d,i,k) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T2_("c,d,i,k") * Id_oo("j,n") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,l) <m,n||c,d> t2(c,d,i,j) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,d,i,j") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(k,m) <l,n||c,d> t2(c,d,i,j) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T2_("c,d,i,j") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,n) <l,m||c,d> t2(c,d,i,j) 
    // flops: o6v2: 2, o4v2: 1, o6v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T2_("c,d,i,j") * Id_oo("k,n") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(k,n) d(j,m) d(i,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += 0.500000 d(k,n) d(i,m) d(j,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += 0.500000 d(j,n) d(k,m) d(i,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += -0.500000 d(j,n) d(i,m) d(k,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += -0.500000 d(i,n) d(k,m) d(j,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += 0.500000 d(i,n) d(j,m) d(k,l) <t,o||c,b> t2(c,a,t,o) 
    // flops: o6v2: 2, o4v2: 1, o2v3: 1, o2v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T2_("c,a,t,o") * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(j,m) d(i,l) <n,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,k,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(i,m) d(j,l) <n,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,k,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,n) d(i,l) <m,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,k,o") * Id_oo("j,n") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(j,n) d(i,m) <l,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,k,o") * Id_oo("j,n") * Id_oo("i,m");

    // H_vooovooo += 1.000000 d(i,n) d(j,l) <m,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,k,o") * Id_oo("i,n") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(i,n) d(j,m) <l,o||c,b> t2(c,a,k,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,k,o") * Id_oo("i,n") * Id_oo("j,m");

    // H_vooovooo += -1.000000 d(k,m) d(i,l) <n,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,j,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(i,m) d(k,l) <n,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,j,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(k,n) d(i,l) <m,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,j,o") * Id_oo("k,n") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(k,n) d(i,m) <l,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,j,o") * Id_oo("k,n") * Id_oo("i,m");

    // H_vooovooo += -1.000000 d(i,n) d(k,l) <m,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,j,o") * Id_oo("i,n") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(i,n) d(k,m) <l,o||c,b> t2(c,a,j,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,j,o") * Id_oo("i,n") * Id_oo("k,m");

    // H_vooovooo += 1.000000 d(k,m) d(j,l) <n,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,i,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,m) d(k,l) <n,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T2_("c,a,i,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(k,n) d(j,l) <m,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,i,o") * Id_oo("k,n") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(k,n) d(j,m) <l,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,i,o") * Id_oo("k,n") * Id_oo("j,m");

    // H_vooovooo += 1.000000 d(j,n) d(k,l) <m,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T2_("c,a,i,o") * Id_oo("j,n") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(j,n) d(k,m) <l,o||c,b> t2(c,a,i,o) 
    // flops: o6v2: 2, o3v3: 1, o4v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T2_("c,a,i,o") * Id_oo("j,n") * Id_oo("k,m");

    // H_vooovooo += -1.000000 d(i,l) <m,n||c,b> t2(c,a,j,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,n")) * this->T2_("c,a,j,k") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(i,m) <l,n||c,b> t2(c,a,j,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,n")) * this->T2_("c,a,j,k") * Id_oo("i,m");

    // H_vooovooo += -1.000000 d(i,n) <l,m||c,b> t2(c,a,j,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,m")) * this->T2_("c,a,j,k") * Id_oo("i,n");

    // H_vooovooo += 1.000000 d(j,l) <m,n||c,b> t2(c,a,i,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,n")) * this->T2_("c,a,i,k") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,m) <l,n||c,b> t2(c,a,i,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,n")) * this->T2_("c,a,i,k") * Id_oo("j,m");

    // H_vooovooo += 1.000000 d(j,n) <l,m||c,b> t2(c,a,i,k) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,m")) * this->T2_("c,a,i,k") * Id_oo("j,n");

    // H_vooovooo += -1.000000 d(k,l) <m,n||c,b> t2(c,a,i,j) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,n")) * this->T2_("c,a,i,j") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(k,m) <l,n||c,b> t2(c,a,i,j) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,n")) * this->T2_("c,a,i,j") * Id_oo("k,m");

    // H_vooovooo += -1.000000 d(k,n) <l,m||c,b> t2(c,a,i,j) 
    // flops: o6v2: 2, o4v3: 1 | mem: o6v2: 2, o4v2: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,m")) * this->T2_("c,a,i,j") * Id_oo("k,n");

    // H_vooovooo += -0.500000 d(a,b) d(k,n) d(j,m) d(i,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_4 * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(k,n) d(i,m) d(j,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_4 * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(j,n) d(k,m) d(i,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_4 * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(j,n) d(i,m) d(k,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_4 * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -0.500000 d(a,b) d(i,n) d(k,m) d(j,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * scalar_4 * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(a,b) d(i,n) d(j,m) d(k,l) <t,o||c,d> t1(c,o) t1(d,t) 
    // flops: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o0v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * scalar_4 * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,m) d(i,l) <n,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("j,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,m) d(j,l) <n,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("i,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(i,l) <m,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("j,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(i,m) <l,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("j,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(j,l) <m,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("i,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(j,m) <l,o||c,d> t1(c,o) t1(d,k) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,k") * Id_oo("i,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,m) d(i,l) <n,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("k,m") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,m) d(k,l) <n,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("i,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(i,l) <m,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("k,n") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(i,m) <l,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("k,n") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) d(k,l) <m,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("i,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,n) d(k,m) <l,o||c,d> t1(c,o) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,j") * Id_oo("i,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,m) d(j,l) <n,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("k,m") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,m) d(k,l) <n,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,n,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("j,m") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) d(j,l) <m,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("k,n") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,n) d(j,m) <l,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("k,n") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) d(k,l) <m,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,m,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("j,n") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,n) d(k,m) <l,o||c,d> t1(c,o) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o2v2: 1, o4v0: 1, o2v1: 1 | mem: o6v2: 2, o6v0: 1, o4v0: 1, o1v1: 1, o2v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,o")) * this->T1_("c,o") * this->T1_("d,i") * Id_oo("j,n") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += 0.500000 d(k,n) d(j,m) d(i,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += -0.500000 d(k,n) d(i,m) d(j,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += -0.500000 d(j,n) d(k,m) d(i,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += 0.500000 d(j,n) d(i,m) d(k,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += 0.500000 d(i,n) d(k,m) d(j,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += -0.500000 d(i,n) d(j,m) d(k,l) <t,o||c,b> t1(c,o) t1(a,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,o") * this->T1_("a,t") * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(a,b) d(i,l) <m,n||c,d> t1(c,k) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,k") * this->T1_("d,j") * Id_oo("i,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(i,m) <l,n||c,d> t1(c,k) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T1_("c,k") * this->T1_("d,j") * Id_oo("i,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(i,n) <l,m||c,d> t1(c,k) t1(d,j) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T1_("c,k") * this->T1_("d,j") * Id_oo("i,n") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,l) <m,n||c,d> t1(c,k) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,k") * this->T1_("d,i") * Id_oo("j,l") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(j,m) <l,n||c,d> t1(c,k) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T1_("c,k") * this->T1_("d,i") * Id_oo("j,m") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(j,n) <l,m||c,d> t1(c,k) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T1_("c,k") * this->T1_("d,i") * Id_oo("j,n") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(j,m) d(i,l) <n,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(i,m) d(j,l) <n,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,n) d(i,l) <m,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(j,n) d(i,m) <l,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,m");

    // H_vooovooo += 1.000000 d(i,n) d(j,l) <m,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(i,n) d(j,m) <l,o||c,b> t1(c,k) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,k") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,m");

    // H_vooovooo += -1.000000 d(a,b) d(k,l) <m,n||c,d> t1(c,j) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,j") * this->T1_("d,i") * Id_oo("k,l") * Id_vv("a,b");

    // H_vooovooo += 1.000000 d(a,b) d(k,m) <l,n||c,d> t1(c,j) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,d,l,n")) * this->T1_("c,j") * this->T1_("d,i") * Id_oo("k,m") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(a,b) d(k,n) <l,m||c,d> t1(c,j) t1(d,i) 
    // flops: o6v2: 2, o6v0: 1, o3v2: 1, o4v1: 1 | mem: o6v2: 2, o6v0: 1, o3v1: 1, o4v0: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,d,l,m")) * this->T1_("c,j") * this->T1_("d,i") * Id_oo("k,n") * Id_vv("a,b");

    // H_vooovooo += -1.000000 d(k,m) d(i,l) <n,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += 1.000000 d(i,m) d(k,l) <n,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(k,n) d(i,l) <m,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,l");

    // H_vooovooo += -1.000000 d(k,n) d(i,m) <l,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,m");

    // H_vooovooo += -1.000000 d(i,n) d(k,l) <m,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,l");

    // H_vooovooo += 1.000000 d(i,n) d(k,m) <l,o||c,b> t1(c,j) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,j") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,m");

    // H_vooovooo += 1.000000 d(k,m) d(j,l) <n,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += -1.000000 d(j,m) d(k,l) <n,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,n,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(k,n) d(j,l) <m,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,l");

    // H_vooovooo += 1.000000 d(k,n) d(j,m) <l,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,m");

    // H_vooovooo += 1.000000 d(j,n) d(k,l) <m,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("c,b,m,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,l");

    // H_vooovooo += -1.000000 d(j,n) d(k,m) <l,o||c,b> t1(c,i) t1(a,o) 
    // flops: o6v2: 2, o4v2: 1, o3v2: 2 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o3v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("c,b,l,o")) * this->T1_("c,i") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,m");

    // H_vooovooo += -0.500000 d(k,n) d(j,m) d(i,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooovooo += 0.500000 d(k,n) d(i,m) d(j,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooovooo += 0.500000 d(j,n) d(k,m) d(i,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooovooo += -0.500000 d(j,n) d(i,m) d(k,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooovooo += -0.500000 d(i,n) d(k,m) d(j,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooovooo += 0.500000 d(i,n) d(j,m) d(k,l) <t,o||c,b> t1(a,o) t1(c,t) 
    // flops: o6v2: 2, o4v2: 1, o2v2: 2, o1v2: 1 | mem: o6v2: 2, o4v2: 1, o2v2: 1, o0v2: 1, o1v1: 1, 
    H_vooovooo("a,i,j,k,b,l,m,n") += 0.500000 * conj(this->antiSymMoints["vvoo"]("c,b,t,o")) * this->T1_("c,t") * this->T1_("a,o") * Id_oo("i,n") * Id_oo("j,m") * Id_oo("k,l");




/// ****** H(i,j,b,l,m,n) ****** ///




    // H_oovooo += 1.000000 d(j,m) d(i,l) f(n,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("j,m") * Id_oo("i,l") * this->fockMatrix_ta["ov"]("n,b");

    // H_oovooo += -1.000000 d(i,m) d(j,l) f(n,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("i,m") * Id_oo("j,l") * this->fockMatrix_ta["ov"]("n,b");

    // H_oovooo += -1.000000 d(j,n) d(i,l) f(m,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("j,n") * Id_oo("i,l") * this->fockMatrix_ta["ov"]("m,b");

    // H_oovooo += 1.000000 d(j,n) d(i,m) f(l,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("j,n") * Id_oo("i,m") * this->fockMatrix_ta["ov"]("l,b");

    // H_oovooo += 1.000000 d(i,n) d(j,l) f(m,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("i,n") * Id_oo("j,l") * this->fockMatrix_ta["ov"]("m,b");

    // H_oovooo += -1.000000 d(i,n) d(j,m) f(l,b) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("i,n") * Id_oo("j,m") * this->fockMatrix_ta["ov"]("l,b");

    // H_oovooo += 1.000000 d(i,l) <m,n||b,j> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("i,l") * conj(this->antiSymMoints["vooo"]("b,j,m,n"));

    // H_oovooo += -1.000000 d(i,m) <l,n||b,j> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("i,m") * conj(this->antiSymMoints["vooo"]("b,j,l,n"));

    // H_oovooo += 1.000000 d(i,n) <l,m||b,j> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("i,n") * conj(this->antiSymMoints["vooo"]("b,j,l,m"));

    // H_oovooo += -1.000000 d(j,l) <m,n||b,i> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("j,l") * conj(this->antiSymMoints["vooo"]("b,i,m,n"));

    // H_oovooo += 1.000000 d(j,m) <l,n||b,i> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") += Id_oo("j,m") * conj(this->antiSymMoints["vooo"]("b,i,l,n"));

    // H_oovooo += -1.000000 d(j,n) <l,m||b,i> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_oovooo("i,j,b,l,m,n") -= Id_oo("j,n") * conj(this->antiSymMoints["vooo"]("b,i,l,m"));

    // H_oovooo += -1.000000 d(j,m) d(i,l) <n,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,n,k")) * this->T1_("a,k") * Id_oo("j,m") * Id_oo("i,l");

    // H_oovooo += 1.000000 d(i,m) d(j,l) <n,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,n,k")) * this->T1_("a,k") * Id_oo("i,m") * Id_oo("j,l");

    // H_oovooo += 1.000000 d(j,n) d(i,l) <m,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,m,k")) * this->T1_("a,k") * Id_oo("j,n") * Id_oo("i,l");

    // H_oovooo += -1.000000 d(j,n) d(i,m) <l,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,l,k")) * this->T1_("a,k") * Id_oo("j,n") * Id_oo("i,m");

    // H_oovooo += -1.000000 d(i,n) d(j,l) <m,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,m,k")) * this->T1_("a,k") * Id_oo("i,n") * Id_oo("j,l");

    // H_oovooo += 1.000000 d(i,n) d(j,m) <l,k||a,b> t1(a,k) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,l,k")) * this->T1_("a,k") * Id_oo("i,n") * Id_oo("j,m");

    // H_oovooo += -1.000000 d(i,l) <m,n||a,b> t1(a,j) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,m,n")) * this->T1_("a,j") * Id_oo("i,l");

    // H_oovooo += 1.000000 d(i,m) <l,n||a,b> t1(a,j) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,l,n")) * this->T1_("a,j") * Id_oo("i,m");

    // H_oovooo += -1.000000 d(i,n) <l,m||a,b> t1(a,j) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T1_("a,j") * Id_oo("i,n");

    // H_oovooo += 1.000000 d(j,l) <m,n||a,b> t1(a,i) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,m,n")) * this->T1_("a,i") * Id_oo("j,l");

    // H_oovooo += -1.000000 d(j,m) <l,n||a,b> t1(a,i) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") -= conj(this->antiSymMoints["vvoo"]("a,b,l,n")) * this->T1_("a,i") * Id_oo("j,m");

    // H_oovooo += 1.000000 d(j,n) <l,m||a,b> t1(a,i) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_oovooo("i,j,b,l,m,n") += conj(this->antiSymMoints["vvoo"]("a,b,l,m")) * this->T1_("a,i") * Id_oo("j,n");




/// ****** H(a,i,j,k,l,m) ****** ///




    // H_vooooo += 1.000000 d(j,m) d(i,l) f(a,k) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("j,m") * Id_oo("i,l") * this->fockMatrix_ta["vo"]("a,k");

    // H_vooooo += -1.000000 d(i,m) d(j,l) f(a,k) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("i,m") * Id_oo("j,l") * this->fockMatrix_ta["vo"]("a,k");

    // H_vooooo += -1.000000 d(k,m) d(i,l) f(a,j) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("k,m") * Id_oo("i,l") * this->fockMatrix_ta["vo"]("a,j");

    // H_vooooo += 1.000000 d(i,m) d(k,l) f(a,j) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("i,m") * Id_oo("k,l") * this->fockMatrix_ta["vo"]("a,j");

    // H_vooooo += 1.000000 d(k,m) d(j,l) f(a,i) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("k,m") * Id_oo("j,l") * this->fockMatrix_ta["vo"]("a,i");

    // H_vooooo += -1.000000 d(j,m) d(k,l) f(a,i) 
    // flops: o5v1: 2, o4v0: 1 | mem: o5v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("j,m") * Id_oo("k,l") * this->fockMatrix_ta["vo"]("a,i");

    // H_vooooo += -1.000000 d(j,m) d(i,l) f(n,k) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["oo"]("n,k") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(j,l) f(n,k) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["oo"]("n,k") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(k,m) d(i,l) f(n,j) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["oo"]("n,j") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(k,l) f(n,j) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["oo"]("n,j") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) d(j,l) f(n,i) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["oo"]("n,i") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) d(k,l) f(n,i) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["oo"]("n,i") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(j,m) d(i,l) f(a,b) t1(b,k) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,k") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) f(a,b) t1(b,k) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,k") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(k,m) d(i,l) f(a,b) t1(b,j) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,j") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) f(a,b) t1(b,j) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,j") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) d(j,l) f(a,b) t1(b,i) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,i") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) f(a,b) t1(b,i) 
    // flops: o5v1: 2, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,i") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(j,m) d(i,l) f(n,b) t2(b,a,k,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(j,l) f(n,b) t2(b,a,k,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(k,m) d(i,l) f(n,b) t2(b,a,j,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(k,l) f(n,b) t2(b,a,j,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) d(j,l) f(n,b) t2(b,a,i,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) d(k,l) f(n,b) t2(b,a,i,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T2_("b,a,i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(i,l) f(m,b) t2(b,a,j,k) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("m,b") * this->T2_("b,a,j,k") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) f(l,b) t2(b,a,j,k) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,j,k") * Id_oo("i,m");

    // H_vooooo += 1.000000 d(j,l) f(m,b) t2(b,a,i,k) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("m,b") * this->T2_("b,a,i,k") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) f(l,b) t2(b,a,i,k) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,i,k") * Id_oo("j,m");

    // H_vooooo += -1.000000 d(k,l) f(m,b) t2(b,a,i,j) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("m,b") * this->T2_("b,a,i,j") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) f(l,b) t2(b,a,i,j) 
    // flops: o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,i,j") * Id_oo("k,m");

    // H_vooooo += -1.000000 d(j,m) d(i,l) f(n,b) t1(b,k) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,k") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(j,l) f(n,b) t1(b,k) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,k") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(k,m) d(i,l) f(n,b) t1(b,j) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,j") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(k,l) f(n,b) t1(b,j) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,j") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) d(j,l) f(n,b) t1(b,i) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,i") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) d(k,l) f(n,b) t1(b,i) t1(a,n) 
    // flops: o5v1: 2, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += this->fockMatrix_ta["ov"]("n,b") * this->T1_("b,i") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(i,l) <m,a||j,k> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("i,l") * this->antiSymMoints["vooo"]("a,m,j,k");

    // H_vooooo += 1.000000 d(i,m) <l,a||j,k> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("i,m") * this->antiSymMoints["vooo"]("a,l,j,k");

    // H_vooooo += 1.000000 d(j,l) <m,a||i,k> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("j,l") * this->antiSymMoints["vooo"]("a,m,i,k");

    // H_vooooo += -1.000000 d(j,m) <l,a||i,k> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("j,m") * this->antiSymMoints["vooo"]("a,l,i,k");

    // H_vooooo += -1.000000 d(k,l) <m,a||i,j> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") += Id_oo("k,l") * this->antiSymMoints["vooo"]("a,m,i,j");

    // H_vooooo += 1.000000 d(k,m) <l,a||i,j> 
    // flops: o5v1: 2 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= Id_oo("k,m") * this->antiSymMoints["vooo"]("a,l,i,j");

    // H_vooooo += 1.000000 d(i,l) <m,n||j,k> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["oooo"]("m,n,j,k") * this->T1_("a,n") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) <l,n||j,k> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["oooo"]("l,n,j,k") * this->T1_("a,n") * Id_oo("i,m");

    // H_vooooo += -1.000000 d(j,l) <m,n||i,k> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["oooo"]("m,n,i,k") * this->T1_("a,n") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) <l,n||i,k> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["oooo"]("l,n,i,k") * this->T1_("a,n") * Id_oo("j,m");

    // H_vooooo += 1.000000 d(k,l) <m,n||i,j> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["oooo"]("m,n,i,j") * this->T1_("a,n") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) <l,n||i,j> t1(a,n) 
    // flops: o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["oooo"]("l,n,i,j") * this->T1_("a,n") * Id_oo("k,m");

    // H_vooooo += 1.000000 d(j,m) d(i,l) <n,a||b,k> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["vovo"]("a,n,b,k") * this->T1_("b,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) <n,a||b,k> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["vovo"]("a,n,b,k") * this->T1_("b,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 P(j,k) d(i,l) <m,a||b,k> t1(b,j) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,m,b,k") * this->T1_("b,j") * Id_oo("i,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += 1.000000 P(j,k) d(i,m) <l,a||b,k> t1(b,j) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,l,b,k") * this->T1_("b,j") * Id_oo("i,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += 1.000000 P(i,k) d(j,l) <m,a||b,k> t1(b,i) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,m,b,k") * this->T1_("b,i") * Id_oo("j,l");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += -1.000000 P(i,k) d(j,m) <l,a||b,k> t1(b,i) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,l,b,k") * this->T1_("b,i") * Id_oo("j,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += -1.000000 d(k,m) d(i,l) <n,a||b,j> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["vovo"]("a,n,b,j") * this->T1_("b,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) <n,a||b,j> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["vovo"]("a,n,b,j") * this->T1_("b,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 P(i,j) d(k,l) <m,a||b,j> t1(b,i) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,m,b,j") * this->T1_("b,i") * Id_oo("k,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += 1.000000 P(i,j) d(k,m) <l,a||b,j> t1(b,i) 
    // flops: o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = this->antiSymMoints["vovo"]("a,l,b,j") * this->T1_("b,i") * Id_oo("k,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += 1.000000 d(k,m) d(j,l) <n,a||b,i> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= this->antiSymMoints["vovo"]("a,n,b,i") * this->T1_("b,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) <n,a||b,i> t1(b,n) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += this->antiSymMoints["vovo"]("a,n,b,i") * this->T1_("b,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -0.500000 d(j,m) d(i,l) <o,n||b,k> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vooo"]("b,k,o,n")) * this->T2_("b,a,o,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += 0.500000 d(i,m) d(j,l) <o,n||b,k> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vooo"]("b,k,o,n")) * this->T2_("b,a,o,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 P(j,k) d(i,l) <m,n||b,k> t2(b,a,j,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,m,n")) * this->T2_("b,a,j,n") * Id_oo("i,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += -1.000000 P(j,k) d(i,m) <l,n||b,k> t2(b,a,j,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,l,n")) * this->T2_("b,a,j,n") * Id_oo("i,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += -1.000000 P(i,k) d(j,l) <m,n||b,k> t2(b,a,i,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,m,n")) * this->T2_("b,a,i,n") * Id_oo("j,l");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += 1.000000 P(i,k) d(j,m) <l,n||b,k> t2(b,a,i,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,l,n")) * this->T2_("b,a,i,n") * Id_oo("j,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += -1.000000 P(j,k) <l,m||b,k> t2(b,a,i,j) 
    // flops: o5v2: 1, o5v1: 2 | mem: o5v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,l,m")) * this->T2_("b,a,i,j");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += 0.500000 d(k,m) d(i,l) <o,n||b,j> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vooo"]("b,j,o,n")) * this->T2_("b,a,o,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += -0.500000 d(i,m) d(k,l) <o,n||b,j> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vooo"]("b,j,o,n")) * this->T2_("b,a,o,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 P(i,j) d(k,l) <m,n||b,j> t2(b,a,i,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,j,m,n")) * this->T2_("b,a,i,n") * Id_oo("k,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += -1.000000 P(i,j) d(k,m) <l,n||b,j> t2(b,a,i,n) 
    // flops: o4v2: 1, o5v1: 3 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,j,l,n")) * this->T2_("b,a,i,n") * Id_oo("k,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += -0.500000 d(k,m) d(j,l) <o,n||b,i> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vooo"]("b,i,o,n")) * this->T2_("b,a,o,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += 0.500000 d(j,m) d(k,l) <o,n||b,i> t2(b,a,o,n) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vooo"]("b,i,o,n")) * this->T2_("b,a,o,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 <l,m||b,i> t2(b,a,j,k) 
    // flops: o5v2: 1, o5v1: 1 | mem: o5v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vooo"]("b,i,l,m")) * this->T2_("b,a,j,k");

    // H_vooooo += -0.500000 d(j,m) d(i,l) <n,a||b,c> t2(b,c,k,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,k,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += 0.500000 d(i,m) d(j,l) <n,a||b,c> t2(b,c,k,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,k,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 0.500000 d(k,m) d(i,l) <n,a||b,c> t2(b,c,j,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,j,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += -0.500000 d(i,m) d(k,l) <n,a||b,c> t2(b,c,j,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,j,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -0.500000 d(k,m) d(j,l) <n,a||b,c> t2(b,c,i,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,i,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += 0.500000 d(j,m) d(k,l) <n,a||b,c> t2(b,c,i,n) 
    // flops: o5v1: 2, o2v3: 1, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T2_("b,c,i,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -0.500000 d(i,l) <m,a||b,c> t2(b,c,j,k) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T2_("b,c,j,k") * Id_oo("i,l");

    // H_vooooo += 0.500000 d(i,m) <l,a||b,c> t2(b,c,j,k) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T2_("b,c,j,k") * Id_oo("i,m");

    // H_vooooo += 0.500000 d(j,l) <m,a||b,c> t2(b,c,i,k) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T2_("b,c,i,k") * Id_oo("j,l");

    // H_vooooo += -0.500000 d(j,m) <l,a||b,c> t2(b,c,i,k) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T2_("b,c,i,k") * Id_oo("j,m");

    // H_vooooo += -0.500000 d(k,l) <m,a||b,c> t2(b,c,i,j) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T2_("b,c,i,j") * Id_oo("k,l");

    // H_vooooo += 0.500000 d(k,m) <l,a||b,c> t2(b,c,i,j) 
    // flops: o3v3: 1, o5v1: 2 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T2_("b,c,i,j") * Id_oo("k,m");

    // H_vooooo += 1.000000 d(j,m) d(i,l) <o,n||b,c> t1(b,n) t2(c,a,k,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,k,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) <o,n||b,c> t1(b,n) t2(c,a,k,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,k,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(k,m) d(i,l) <o,n||b,c> t1(b,n) t2(c,a,j,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,j,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) <o,n||b,c> t1(b,n) t2(c,a,j,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,j,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) d(j,l) <o,n||b,c> t1(b,n) t2(c,a,i,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,i,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) <o,n||b,c> t1(b,n) t2(c,a,i,o) 
    // flops: o5v1: 2, o2v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 2, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T2_("c,a,i,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(i,l) <m,n||b,c> t1(b,n) t2(c,a,j,k) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,n") * this->T2_("c,a,j,k") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) <l,n||b,c> t1(b,n) t2(c,a,j,k) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,n") * this->T2_("c,a,j,k") * Id_oo("i,m");

    // H_vooooo += -1.000000 d(j,l) <m,n||b,c> t1(b,n) t2(c,a,i,k) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,n") * this->T2_("c,a,i,k") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) <l,n||b,c> t1(b,n) t2(c,a,i,k) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,n") * this->T2_("c,a,i,k") * Id_oo("j,m");

    // H_vooooo += 1.000000 d(k,l) <m,n||b,c> t1(b,n) t2(c,a,i,j) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,n") * this->T2_("c,a,i,j") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) <l,n||b,c> t1(b,n) t2(c,a,i,j) 
    // flops: o5v1: 2, o3v2: 1, o2v2: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,n") * this->T2_("c,a,i,j") * Id_oo("k,m");

    // H_vooooo += 0.500000 d(j,m) d(i,l) <o,n||b,c> t1(b,k) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,k") * this->T2_("c,a,o,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -0.500000 d(i,m) d(j,l) <o,n||b,c> t1(b,k) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,k") * this->T2_("c,a,o,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 P(j,k) d(i,l) <m,n||b,c> t1(b,k) t2(c,a,j,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,k") * this->T2_("c,a,j,n") * Id_oo("i,l");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += 1.000000 P(j,k) d(i,m) <l,n||b,c> t1(b,k) t2(c,a,j,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,k") * this->T2_("c,a,j,n") * Id_oo("i,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += 1.000000 P(i,k) d(j,l) <m,n||b,c> t1(b,k) t2(c,a,i,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,k") * this->T2_("c,a,i,n") * Id_oo("j,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += -1.000000 P(i,k) d(j,m) <l,n||b,c> t1(b,k) t2(c,a,i,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,k") * this->T2_("c,a,i,n") * Id_oo("j,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += 1.000000 P(j,k) <l,m||b,c> t1(b,k) t2(c,a,i,j) 
    // flops: o5v2: 1, o5v1: 2, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T1_("b,k") * this->T2_("c,a,i,j");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += -0.500000 d(k,m) d(i,l) <o,n||b,c> t1(b,j) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,j") * this->T2_("c,a,o,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 0.500000 d(i,m) d(k,l) <o,n||b,c> t1(b,j) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,j") * this->T2_("c,a,o,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 P(i,j) d(k,l) <m,n||b,c> t1(b,j) t2(c,a,i,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,j") * this->T2_("c,a,i,n") * Id_oo("k,l");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += 1.000000 P(i,j) d(k,m) <l,n||b,c> t1(b,j) t2(c,a,i,n) 
    // flops: o4v2: 1, o5v1: 3, o3v2: 1 | mem: o5v1: 2, o3v1: 2, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,j") * this->T2_("c,a,i,n") * Id_oo("k,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += 0.500000 d(k,m) d(j,l) <o,n||b,c> t1(b,i) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,i") * this->T2_("c,a,o,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -0.500000 d(j,m) d(k,l) <o,n||b,c> t1(b,i) t2(c,a,o,n) 
    // flops: o5v1: 2, o3v2: 2, o3v1: 1 | mem: o5v1: 2, o3v1: 2, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,i") * this->T2_("c,a,o,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 <l,m||b,c> t1(b,i) t2(c,a,j,k) 
    // flops: o5v2: 1, o5v1: 1, o3v2: 1 | mem: o5v1: 2, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T1_("b,i") * this->T2_("c,a,j,k");

    // H_vooooo += 0.500000 d(j,m) d(i,l) <o,n||b,c> t1(a,n) t2(b,c,k,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,k,o") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -0.500000 d(i,m) d(j,l) <o,n||b,c> t1(a,n) t2(b,c,k,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,k,o") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -0.500000 d(k,m) d(i,l) <o,n||b,c> t1(a,n) t2(b,c,j,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,j,o") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 0.500000 d(i,m) d(k,l) <o,n||b,c> t1(a,n) t2(b,c,j,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,j,o") * this->T1_("a,n") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 0.500000 d(k,m) d(j,l) <o,n||b,c> t1(a,n) t2(b,c,i,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,i,o") * this->T1_("a,n") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -0.500000 d(j,m) d(k,l) <o,n||b,c> t1(a,n) t2(b,c,i,o) 
    // flops: o5v1: 2, o3v2: 1, o3v1: 1, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T2_("b,c,i,o") * this->T1_("a,n") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 0.500000 d(i,l) <m,n||b,c> t1(a,n) t2(b,c,j,k) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T2_("b,c,j,k") * this->T1_("a,n") * Id_oo("i,l");

    // H_vooooo += -0.500000 d(i,m) <l,n||b,c> t1(a,n) t2(b,c,j,k) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T2_("b,c,j,k") * this->T1_("a,n") * Id_oo("i,m");

    // H_vooooo += -0.500000 d(j,l) <m,n||b,c> t1(a,n) t2(b,c,i,k) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T2_("b,c,i,k") * this->T1_("a,n") * Id_oo("j,l");

    // H_vooooo += 0.500000 d(j,m) <l,n||b,c> t1(a,n) t2(b,c,i,k) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T2_("b,c,i,k") * this->T1_("a,n") * Id_oo("j,m");

    // H_vooooo += 0.500000 d(k,l) <m,n||b,c> t1(a,n) t2(b,c,i,j) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T2_("b,c,i,j") * this->T1_("a,n") * Id_oo("k,l");

    // H_vooooo += -0.500000 d(k,m) <l,n||b,c> t1(a,n) t2(b,c,i,j) 
    // flops: o4v2: 1, o5v1: 2, o4v1: 1 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= 0.500000 * conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T2_("b,c,i,j") * this->T1_("a,n") * Id_oo("k,m");

    // H_vooooo += 1.000000 d(j,m) d(i,l) <o,n||b,k> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vooo"]("b,k,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) <o,n||b,k> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vooo"]("b,k,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += 1.000000 P(j,k) d(i,l) <m,n||b,k> t1(b,j) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,m,n")) * this->T1_("b,j") * this->T1_("a,n") * Id_oo("i,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += -1.000000 P(j,k) d(i,m) <l,n||b,k> t1(b,j) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,l,n")) * this->T1_("b,j") * this->T1_("a,n") * Id_oo("i,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,k,j,l,m");

    // H_vooooo += -1.000000 P(i,k) d(j,l) <m,n||b,k> t1(b,i) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,m,n")) * this->T1_("b,i") * this->T1_("a,n") * Id_oo("j,l");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += 1.000000 P(i,k) d(j,m) <l,n||b,k> t1(b,i) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,k,l,n")) * this->T1_("b,i") * this->T1_("a,n") * Id_oo("j,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,k,j,i,l,m");

    // H_vooooo += -1.000000 d(k,m) d(i,l) <o,n||b,j> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vooo"]("b,j,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) <o,n||b,j> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vooo"]("b,j,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 P(i,j) d(k,l) <m,n||b,j> t1(b,i) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,j,m,n")) * this->T1_("b,i") * this->T1_("a,n") * Id_oo("k,l");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += -1.000000 P(i,j) d(k,m) <l,n||b,j> t1(b,i) t1(a,n) 
    // flops: o5v1: 3, o4v1: 2 | mem: o5v1: 2, o3v1: 1, o4v0: 1, 
    tempPerm_vooooo("a,i,j,k,l,m") = conj(this->antiSymMoints["vooo"]("b,j,l,n")) * this->T1_("b,i") * this->T1_("a,n") * Id_oo("k,m");
    H_vooooo("a,i,j,k,l,m") -= tempPerm_vooooo("a,i,j,k,l,m");
    H_vooooo("a,i,j,k,l,m") += tempPerm_vooooo("a,j,i,k,l,m");

    // H_vooooo += 1.000000 d(k,m) d(j,l) <o,n||b,i> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vooo"]("b,i,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) <o,n||b,i> t1(b,n) t1(a,o) 
    // flops: o5v1: 2, o3v1: 2, o2v1: 1 | mem: o5v1: 2, o3v1: 1, o1v1: 1, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vooo"]("b,i,o,n")) * this->T1_("b,n") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(j,m) d(i,l) <n,a||b,c> t1(b,n) t1(c,k) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,k") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) <n,a||b,c> t1(b,n) t1(c,k) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,k") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(k,m) d(i,l) <n,a||b,c> t1(b,n) t1(c,j) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,j") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) <n,a||b,c> t1(b,n) t1(c,j) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,j") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) d(j,l) <n,a||b,c> t1(b,n) t1(c,i) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,i") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) <n,a||b,c> t1(b,n) t1(c,i) 
    // flops: o5v1: 2, o1v3: 1, o3v1: 1, o1v2: 1 | mem: o5v1: 2, o3v1: 1, o0v2: 1, o1v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,n")) * this->T1_("b,n") * this->T1_("c,i") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(i,l) <m,a||b,c> t1(b,k) t1(c,j) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T1_("b,k") * this->T1_("c,j") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) <l,a||b,c> t1(b,k) t1(c,j) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,k") * this->T1_("c,j") * Id_oo("i,m");

    // H_vooooo += -1.000000 d(j,l) <m,a||b,c> t1(b,k) t1(c,i) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T1_("b,k") * this->T1_("c,i") * Id_oo("j,l");

    // H_vooooo += 1.000000 d(j,m) <l,a||b,c> t1(b,k) t1(c,i) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,k") * this->T1_("c,i") * Id_oo("j,m");

    // H_vooooo += 1.000000 d(k,l) <m,a||b,c> t1(b,j) t1(c,i) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvvo"]("b,c,a,m")) * this->T1_("b,j") * this->T1_("c,i") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(k,m) <l,a||b,c> t1(b,j) t1(c,i) 
    // flops: o5v1: 2, o2v3: 1, o3v2: 1 | mem: o5v1: 2, o2v2: 1, o3v1: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,j") * this->T1_("c,i") * Id_oo("k,m");

    // H_vooooo += 1.000000 d(j,m) d(i,l) <o,n||b,c> t1(b,n) t1(c,k) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,k") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("i,l");

    // H_vooooo += -1.000000 d(i,m) d(j,l) <o,n||b,c> t1(b,n) t1(c,k) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,k") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(k,m) d(i,l) <o,n||b,c> t1(b,n) t1(c,j) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,j") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) d(k,l) <o,n||b,c> t1(b,n) t1(c,j) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,j") * this->T1_("a,o") * Id_oo("i,m") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) d(j,l) <o,n||b,c> t1(b,n) t1(c,i) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,i") * this->T1_("a,o") * Id_oo("k,m") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) d(k,l) <o,n||b,c> t1(b,n) t1(c,i) t1(a,o) 
    // flops: o5v1: 2, o2v2: 1, o3v1: 1, o2v1: 2 | mem: o5v1: 2, o3v1: 1, o1v1: 2, o2v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,o,n")) * this->T1_("b,n") * this->T1_("c,i") * this->T1_("a,o") * Id_oo("j,m") * Id_oo("k,l");

    // H_vooooo += -1.000000 d(i,l) <m,n||b,c> t1(b,k) t1(c,j) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,k") * this->T1_("c,j") * this->T1_("a,n") * Id_oo("i,l");

    // H_vooooo += 1.000000 d(i,m) <l,n||b,c> t1(b,k) t1(c,j) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,k") * this->T1_("c,j") * this->T1_("a,n") * Id_oo("i,m");

    // H_vooooo += 1.000000 d(j,l) <m,n||b,c> t1(b,k) t1(c,i) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,k") * this->T1_("c,i") * this->T1_("a,n") * Id_oo("j,l");

    // H_vooooo += -1.000000 d(j,m) <l,n||b,c> t1(b,k) t1(c,i) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,k") * this->T1_("c,i") * this->T1_("a,n") * Id_oo("j,m");

    // H_vooooo += -1.000000 d(k,l) <m,n||b,c> t1(b,j) t1(c,i) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") -= conj(this->antiSymMoints["vvoo"]("b,c,m,n")) * this->T1_("b,j") * this->T1_("c,i") * this->T1_("a,n") * Id_oo("k,l");

    // H_vooooo += 1.000000 d(k,m) <l,n||b,c> t1(b,j) t1(c,i) t1(a,n) 
    // flops: o5v1: 2, o3v2: 1, o4v1: 2 | mem: o5v1: 2, o3v1: 2, o4v0: 1, 
    H_vooooo("a,i,j,k,l,m") += conj(this->antiSymMoints["vvoo"]("b,c,l,n")) * this->T1_("b,j") * this->T1_("c,i") * this->T1_("a,n") * Id_oo("k,m");



    TAmanager.free("oooo", std::move(tempPerm_oooo));
    TAmanager.free("vooooo", std::move(tempPerm_vooooo));
    TAmanager.free("vvoooooo", std::move(tempPerm_vvoooooo));
    
	TAmanager.free("vv", std::move(Id_vv));

  }
}

