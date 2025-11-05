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
#include <cibuilder/neo/impl.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void RealTimeNEOCI<MatsT, IntsT>::calculateDipole() {

    SingleSlater<MatsT,IntsT> * ess_ptr = &neomcref_->ewfn_->reference();
    SingleSlater<MatsT,IntsT> * pss_ptr = &neomcref_->pwfn_->reference();

    size_t enCorrO = neomcref_->ewfn_->MOPartition.nCorrO;
    size_t pnCorrO = neomcref_->pwfn_->MOPartition.nCorrO;

    cqmatrix::Matrix<MatsT> eoneRDM(enCorrO);
    cqmatrix::Matrix<MatsT> poneRDM(pnCorrO);

    if(this->curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator)
    {
      auto vecManagerDerived = std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT>>(this->vecManager);
      auto derived_C_real_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(vecManagerDerived->C_real_t);
      auto derived_C_imag_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(vecManagerDerived->C_imag_t);

      // electronic contribution
      cqmatrix::Matrix<MatsT> eoneRDM_r(enCorrO);
      cqmatrix::Matrix<MatsT> eoneRDM_i(enCorrO);
      neocibuilder_->computeOneRDM(*neomcref_,derived_C_real_t->getPtr(), eoneRDM_r);
      neocibuilder_->computeOneRDM(*neomcref_,derived_C_imag_t->getPtr(), eoneRDM_i);
      eoneRDM = eoneRDM_r + eoneRDM_i;

      // Nuclear contribution
      cqmatrix::Matrix<MatsT> poneRDM_r(pnCorrO);
      cqmatrix::Matrix<MatsT> poneRDM_i(pnCorrO);
      neocibuilder_->computePOneRDM(*neomcref_,derived_C_real_t->getPtr(), poneRDM_r);
      neocibuilder_->computePOneRDM(*neomcref_,derived_C_imag_t->getPtr(), poneRDM_i);
      poneRDM = poneRDM_r + poneRDM_i;
    }
    else if(this->curState.curStep == RealTimeAlgorithm::RTRungeKuttaOrderFour)
    {
      auto vecManagerDerived = std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT>>(this->vecManager);
      auto derived_C_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(vecManagerDerived->C_t);
      neocibuilder_->computeOneRDM(*neomcref_,derived_C_t->getPtr(),eoneRDM);
      neocibuilder_->computePOneRDM(*neomcref_,derived_C_t->getPtr(),poneRDM);
    }
    neomcref_->ewfn_->rdm2pdm(eoneRDM);
    neomcref_->pwfn_->rdm2pdm(poneRDM);
    EMPerturbation emPert;
    ess_ptr->computeMultipole(emPert,{ELECTRIC_DIPOLE});
    pss_ptr->computeMultipole(emPert,{ELECTRIC_DIPOLE});
    // Each of ess_ptr->elecDipole and pss_ptr->elecDipole contain one copy of 
    // the classical nuclear dipole moment, so subtract that out here and store the
    // total dipole moment in this->Dipole
    for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
    {
      this->Dipole[iXYZ] = ess_ptr->elecDipole[iXYZ] + pss_ptr->elecDipole[iXYZ] - classicalNucDipole[iXYZ];
      protDipole[iXYZ] = pss_ptr->elecDipole[iXYZ] - classicalNucDipole[iXYZ];
    }

}
}; // namespace ChronusQ
