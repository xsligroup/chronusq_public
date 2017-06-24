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

#include <singleslater/kohnsham.hpp>
#include <util/preprocessor.hpp>
#include <quantum/preprocessor.hpp>

#define KOHNSHAM_COLLECTIVE_OP(OP_MEMBER,OP_VEC_OP) \
  OP_MEMBER(this,other,isGGA_)\
  OP_MEMBER(this,other,functionals)\
  OP_MEMBER(this,other,intParam)\
  OP_MEMBER(this,other,XCEnergy);

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  KohnSham<MatsT,IntsT>::KohnSham(const KohnSham<MatsU,IntsT> &other, int dummy) :
    SingleSlater<MatsT,IntsT>(dynamic_cast<const SingleSlater<MatsU,IntsT>&>(other),dummy),
    QuantumBase(dynamic_cast<const QuantumBase&>(other)),
    WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other)),
    VXC(std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(*other.VXC)),
    doVXC_(other.doVXC_)
    { 
      KOHNSHAM_COLLECTIVE_OP(COPY_OTHER_MEMBER,COPY_OTHER_MEMBER_VEC_OP);
    };

  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  KohnSham<MatsT,IntsT>::KohnSham(KohnSham<MatsU,IntsT> &&other, int dummy) :
    SingleSlater<MatsT,IntsT>(dynamic_cast<SingleSlater<MatsU,IntsT>&&>(std::move(other)),dummy),
    QuantumBase(dynamic_cast<QuantumBase&&>(std::move(other))),
    WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(std::move(other))),
    VXC(std::move(other.VXC)), doVXC_(std::move(other.doVXC_))
    { 
      KOHNSHAM_COLLECTIVE_OP(MOVE_OTHER_MEMBER,MOVE_OTHER_MEMBER_VEC_OP);
    };

  template <typename MatsT, typename IntsT>
  KohnSham<MatsT,IntsT>::KohnSham(const KohnSham<MatsT,IntsT> &other) :
    KohnSham(other,0) { }

  template <typename MatsT, typename IntsT>
  KohnSham<MatsT,IntsT>::KohnSham(KohnSham<MatsT,IntsT> &&other) :
    KohnSham(std::move(other),0) { }



}; // namespace ChronusQ

#include <singleslater/kohnsham/vxc.hpp> // VXC build
#include <singleslater/kohnsham/gradient.hpp> // VXC build

