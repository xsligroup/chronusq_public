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

#include <corehbuilder/x2c.hpp>
#include <util/preprocessor.hpp>

// FIXME: For copy and move, this only populates the lists, not the
// explicit pointers
#define X2C_COLLECTIVE_OP(OP_OP,OP_VEC_OP) \
  /* Handle Operators */\
  OP_OP(IntsT,this,other,mapPrim2Cont);\
  OP_OP(IntsT,this,other,UK);\
  OP_OP(double,this,other,p);\
  OP_OP(MatsT,this,other,UL);\
  OP_OP(MatsT,this,other,US);


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  X2C<MatsT,IntsT>::X2C(const X2C<MatsT,IntsT> &other) :
    X2C(other,0) {}

  template <typename MatsT, typename IntsT>
  X2C<MatsT,IntsT>::X2C(X2C<MatsT,IntsT> &&other) :
    X2C(std::move(other),0) {}

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  X2C<MatsT,IntsT>::X2C(const X2C<MatsU,IntsT> &other, int dummy) :
    aoints_(other.aoints_), ssOptions_(other.ssOptions_),
    molecule_(other.molecule_), basisSet_(other.basisSet_),
    uncontractedBasis_(other.uncontractedBasis_),
    uncontractedInts_(other.uncontractedInts_),
    nPrimUse_(other.nPrimUse_),
    W(other.W ? std::make_shared<cqmatrix::Matrix<MatsT>>(*other.W) : nullptr) {

    X2C_COLLECTIVE_OP(COPY_OTHER_MEMBER_OP, COPY_OTHER_MEMBER_VEC_OP)

  }

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  X2C<MatsT,IntsT>::X2C(X2C<MatsU,IntsT> &&other, int dummy) :
    aoints_(other.aoints_), ssOptions_(other.ssOptions_),
    molecule_(other.molecule_), basisSet_(other.basisSet_),
    uncontractedBasis_(other.uncontractedBasis_),
    uncontractedInts_(other.uncontractedInts_),
    nPrimUse_(other.nPrimUse_),
    W(other.W ? std::make_shared<cqmatrix::Matrix<MatsT>>(*other.W) : nullptr) {

    X2C_COLLECTIVE_OP(MOVE_OTHER_MEMBER_OP, MOVE_OTHER_MEMBER_VEC_OP)

  }

  template <typename MatsT, typename IntsT>
  void X2C<MatsT,IntsT>::dealloc() {

    X2C_COLLECTIVE_OP(DEALLOC_OP_5, DEALLOC_VEC_OP_5)

  }

}; // namespace ChronusQ
