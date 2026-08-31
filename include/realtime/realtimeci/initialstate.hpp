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
#include <realtime.hpp>

namespace ChronusQ {

template <typename MatsT, typename IntsT>
void RealTimeCI<MatsT, IntsT>::genInitialState() {
  if (this->curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT>>(
            this->vecManager);

    RTMS::fill(vecManagerDerived->C_real_t, 0.0 );
    auto derived_C_real_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(
        vecManagerDerived->C_real_t);
    if (vecManagerDerived->initmethod == MSInitialState::LinearCombination) {
      for (auto state : vecManagerDerived->init_detail) {
        // RTMS::add(reference_->CIVecs[state.second - 1],
        // vecManagerDerived->C_real_t, , state.first);
        auto NDet = reference_->NDet;
        blas::axpy(NDet, state.first, reference_->CIVecs[state.second - 1], 1,
                   derived_C_real_t->getPtr(), 1);
      }
    } else if (vecManagerDerived->initmethod == MSInitialState::CustomCI) {
      for (auto det : vecManagerDerived->init_detail) {
        vecManagerDerived->C_real_t->set(det.second - 1, 0, det.first);
      }
    } else {
      CErr("Unclear how you'd like to create your RTCI initial State?");
    }
    double norm;
    RTMS::normalize(vecManagerDerived->C_real_t, norm);
    // Norm check set to 1e-6 which I /think/ matches the default CIConv Tol
    if (abs(norm - 1.0) > 1e-6) {
      CErr("Norm of initial CI Vector is neq 1.0, check the input!");
    }
  } else if (this->curState.curStep ==
             RealTimeAlgorithm::RTRungeKuttaOrderFour) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT>>(
            this->vecManager);
    RTMS::fill(vecManagerDerived->C_t, 0.0);
    auto derived_C_t =
        std::dynamic_pointer_cast<RawVectors<MatsT>>(vecManagerDerived->C_t);
    if (vecManagerDerived->initmethod == MSInitialState::LinearCombination) {
      for (auto state : vecManagerDerived->init_detail) {
        // RTMS::add(reference_->CIVecs[state.second - 1],
        // vecManagerDerived->C_t, state.first);
        auto NDet = reference_->NDet;
        blas::axpy(NDet, state.first, reference_->CIVecs[state.second - 1], 1,
                   derived_C_t->getPtr(), 1);
      }
    } else if (vecManagerDerived->initmethod == MSInitialState::CustomCI) {
      for (auto det : vecManagerDerived->init_detail) {
        // horrible lol
        vecManagerDerived->C_t->set(det.second - 1, 0, det.first);
      }
    } else {
      CErr("Unclear how you'd like to create your RTCI initial State?");
    }
    dcomplex norm;
    RTMS::normalize(vecManagerDerived->C_t, norm);
    // Norm check set to 1e-6 which I /think/ matches the default CIConv Tol
    if (abs(norm - dcomplex(1.0)) > 1e-6) {
      CErr("Norm of initial CI Vector is neq 1.0, check the input!");
    }
  }
};

}; // namespace ChronusQ
