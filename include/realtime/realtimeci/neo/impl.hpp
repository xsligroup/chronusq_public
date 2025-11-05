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
  void RealTimeNEOCI<MatsT, IntsT>::saveState(EMPerturbation & pert) {

    this->data.ProtDipole.push_back(protDipole);

    // Save all the usual RT Data
    RealTimeCI<MatsT,IntsT>::saveState(pert);

    // Save additional NEO specific quantities
    if (this->savFile.exists())
    {
      size_t nSteps = 0;
      size_t maxStep = (size_t)((this->intScheme.tMax + this->intScheme.deltaT / 4) / this->intScheme.deltaT);

      if ((this->curState.iStep + 1) % this->intScheme.iSave == 0 and this->curState.iStep != this->intScheme.restoreStep)
        nSteps = this->intScheme.iSave;
      else if (this->curState.iStep == maxStep) 
        nSteps = (this->curState.iStep - this->intScheme.restoreStep) % this->intScheme.iSave + 1;

      size_t lastPos = this->curState.iStep - nSteps + 1;
      size_t memLastPos = this->data.Time.size() - nSteps;

      if(nSteps != 0)
      {
        this->savFile.partialWriteData("RTNEW/PROT_LEN_ELEC_DIPOLE",&this->data.ProtDipole[memLastPos][0],
                                       {lastPos,0},{nSteps,3},{memLastPos,0},{this->data.Time.size(),3});
      }
    }
  }

  template <typename MatsT, typename IntsT>
  void RealTimeNEOCI<MatsT, IntsT>::createRTDataSets(size_t maxPoints) {

    // Create all the usual RT DataSets
    RealTimeCI<MatsT,IntsT>::createRTDataSets(maxPoints);

    if(this->restart)
      return;

    if (maxPoints == 0)
      maxPoints = (size_t)(((this->intScheme.tMax + this->intScheme.deltaT / 4) / this->intScheme.deltaT) + 1);

    // Containers for saving
    this->savFile.template createDataSet<double>("RTNEW/CLASSICAL_NUC_LEN_ELEC_DIPOLE",{3});
    this->savFile.template createDataSet<double>("RTNEW/PROT_LEN_ELEC_DIPOLE",{maxPoints,3});

    this->savFile.safeWriteData("RTNEW/CLASSICAL_NUC_LEN_ELEC_DIPOLE",&classicalNucDipole[0],{3});

  }

}; // namespace ChronusQ

#include <realtime/realtimeci/neo/cube.hpp>
#include <realtime/realtimeci/neo/dipole.hpp>
