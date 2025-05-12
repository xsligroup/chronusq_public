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
void RealTimeCI<MatsT, IntsT>::calculateDipole() {

  SingleSlater<MatsT, IntsT> *ss_ptr = &reference_->reference();
  size_t nCorrO = reference_->MOPartition.nCorrO;

  cqmatrix::Matrix<MatsT> oneRDM(nCorrO);
  if (this->curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT>>(
            this->vecManager);
    auto derived_C_real_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(
        vecManagerDerived->C_real_t);
    auto derived_C_imag_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(
        vecManagerDerived->C_imag_t);
    cqmatrix::Matrix<MatsT> oneRDM_r(nCorrO);
    cqmatrix::Matrix<MatsT> oneRDM_i(nCorrO);
    reference_->ciBuilder->computeOneRDM(*reference_,
                                         derived_C_real_t->getPtr(), oneRDM_r);
    reference_->ciBuilder->computeOneRDM(*reference_,
                                         derived_C_imag_t->getPtr(), oneRDM_i);
    oneRDM = oneRDM_r + oneRDM_i;
  } else if (this->curState.curStep ==
             RealTimeAlgorithm::RTRungeKuttaOrderFour) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT>>(
            this->vecManager);
    auto derived_C_t =
        std::dynamic_pointer_cast<RawVectors<MatsT>>(vecManagerDerived->C_t);
    reference_->ciBuilder->computeOneRDM(*reference_, derived_C_t->getPtr(),
                                         oneRDM);
  }
  // Convert to AO basis and update PDM in ref
  reference_->rdm2pdm(oneRDM);
  EMPerturbation emPert;
  ss_ptr->computeMultipole(emPert);
  std::copy(ss_ptr->elecDipole.begin(), ss_ptr->elecDipole.end(),
            this->Dipole.begin());
}
}; // namespace ChronusQ
