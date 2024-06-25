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
#include <realtime.hpp>

namespace ChronusQ {

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::genInitialState() {

  if (curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT *>>(
            vecManager);
    auto derived_ref =
        std::dynamic_pointer_cast<MCWaveFunction<MatsT, IntsT>>(reference_);
    if (derived_ref)
      vecManagerDerived->buildInitCIVec(derived_ref);
  } else if (curState.curStep == RealTimeAlgorithm::RTRungeKuttaOrderFour) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT *>>(
            vecManager);
    auto derived_ref =
        std::dynamic_pointer_cast<MCWaveFunction<MatsT, IntsT>>(reference_);
    if (derived_ref)
      vecManagerDerived->buildInitCIVec(derived_ref);
  }
};

template <typename oper_t>
template <typename MatsT, typename IntsT>
void RealTimeMultiSlaterVectorManagerSSO<oper_t>::buildInitCIVec(
    std::shared_ptr<MCWaveFunction<MatsT, IntsT>> ref) {
  auto NDet = this->get_vecSize_();

  auto derived_ref = dynamic_cast<MCWaveFunction<MatsT, IntsT> *>(ref.get());

  RTMS::fill(C_real_t, 0.0, NDet);
  if (this->initmethod == MSInitialState::LinearCombination) {
    for (auto state : this->init_detail) {
      RTMS::add(derived_ref->CIVecs[state.second - 1], C_real_t, NDet,
                state.first);
    }
  } else if (this->initmethod == MSInitialState::CustomCI) {
    for (auto det : this->init_detail) {
      C_real_t[det.second - 1] = det.first;
    }
  } else {
    CErr("Unclear how you'd like to create your RTCI initial State?");
  }
  double norm;
  RTMS::normalize(C_real_t, NDet, norm);
  // Norm check set to 1e-6 which I /think/ matches the default CIConv Tol
  if (abs(norm - 1.0) > 1e-6) {
    CErr("Norm of initial CI Vector is neq 1.0, check the input!");
  }

}; // RealTimeMultiSlaterVectorManagerSSO::buildInitCIVec

template <typename oper_t>
template <typename MatsT, typename IntsT>
void RealTimeMultiSlaterVectorManagerRK4<oper_t>::buildInitCIVec(
    std::shared_ptr<MCWaveFunction<MatsT, IntsT>> ref) {
  const auto NDet = this->get_vecSize_();
  auto derived_ref = dynamic_cast<MCWaveFunction<MatsT, IntsT> *>(ref.get());
  RTMS::fill(C_t, 0.0, NDet);
  if (this->initmethod == MSInitialState::LinearCombination) {
    for (auto state : this->init_detail) {
      RTMS::add(derived_ref->CIVecs[state.second - 1], C_t, NDet, state.first);
    }
  } else if (this->initmethod == MSInitialState::CustomCI) {
    for (auto det : this->init_detail) {
      C_t[det.second - 1] = det.first;
    }
  } else {
    CErr("Unclear how you'd like to create your RTCI initial State?");
  }
  dcomplex norm;
  RTMS::normalize(C_t, NDet, norm);
  // Norm check set to 1e-6 which I /think/ matches the default CIConv Tol
  if (abs(norm - dcomplex(1.0)) > 1e-6) {
    CErr("Norm of initial CI Vector is neq 1.0, check the input!");
  }
}; // RealTimeMultiSlaterVectorManagerRK4::buildInitCIVec

}; // namespace ChronusQ
