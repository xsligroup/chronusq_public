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
void RealTimeMultiSlater<MatsT, IntsT>::calculateDipole() {
  auto derived_ref =
      std::dynamic_pointer_cast<MCWaveFunction<MatsT, IntsT>>(reference_);

  SingleSlater<MatsT, IntsT> *ss_ptr = &derived_ref->reference();
  size_t nCorrO = derived_ref->MOPartition.nCorrO;

  cqmatrix::Matrix<MatsT> oneRDM(nCorrO);
  if (curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT *>>(
            vecManager);
    cqmatrix::Matrix<MatsT> oneRDM_r(nCorrO);
    cqmatrix::Matrix<MatsT> oneRDM_i(nCorrO);
    derived_ref->ciBuilder->computeOneRDM(
        *derived_ref, vecManagerDerived->C_real_t, oneRDM_r);
    derived_ref->ciBuilder->computeOneRDM(
        *derived_ref, vecManagerDerived->C_imag_t, oneRDM_i);
    oneRDM = oneRDM_r + oneRDM_i;
  } else if (curState.curStep == RealTimeAlgorithm::RTRungeKuttaOrderFour) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT *>>(
            vecManager);
    derived_ref->ciBuilder->computeOneRDM(*derived_ref, vecManagerDerived->C_t,
                                          oneRDM);
  }
  // Convert to AO basis and update PDM in ref
  derived_ref->rdm2pdm(oneRDM);
  EMPerturbation emPert;
  ss_ptr->computeMultipole(emPert);
  std::copy(ss_ptr->elecDipole.begin(), ss_ptr->elecDipole.end(),
            Dipole.begin());
}
}; // namespace ChronusQ
