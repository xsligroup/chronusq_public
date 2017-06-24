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

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <matrix.hpp>
#include <realtime.hpp>

#include <util/matout.hpp>
#include <util/timer.hpp>

namespace ChronusQ {

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::createRTDataSets(size_t maxPoints) {
  if (restart)
    return;

  if (maxPoints == 0)
    maxPoints =
        (size_t)(((intScheme.tMax + intScheme.deltaT / 4) / intScheme.deltaT) +
                 1);
  // maxPoints = intScheme.tMax / intScheme.deltaT + 1;

  savFile.createGroup("RT");

  savFile.createDataSet<double>("RT/TIME", {maxPoints});
  savFile.createDataSet<double>("RT/ENERGY", {maxPoints});
  savFile.createDataSet<double>("RT/LEN_ELEC_DIPOLE", {maxPoints, 3});
  savFile.createDataSet<double>("RT/LEN_ELEC_DIPOLE_FIELD", {maxPoints, 3});
  if (this->CIPopFreq != 0) {
    hsize_t nPop = (maxPoints / this->CIPopFreq);
    if (this->CIPopFreq != 1 && (maxPoints - 1) % this->CIPopFreq == 0)
      nPop += 1;
    if (maxPoints == 1)
      nPop = 1;
    hsize_t NStates = vecManager->get_vecSize_();
    savFile.createDataSet<double>("RT/CIPOPULATION", {nPop, NStates});
  }
  if (this->RealTimeCorrelationFunctionFreq != 0) {
    hsize_t maxRTCorrPts = (intScheme.tMax + intScheme.deltaT / 4 -
                            this->RealTimeCorrelationFunctionStart) /
                               intScheme.deltaT +
                           1;
    hsize_t nRTCorr = (maxRTCorrPts / this->RealTimeCorrelationFunctionFreq);
    if (this->RealTimeCorrelationFunctionFreq != 1 &&
        (maxRTCorrPts - 1) % this->RealTimeCorrelationFunctionFreq == 0)
      nRTCorr += 1;
    if (maxRTCorrPts == 1)
      nRTCorr = 1;
    savFile.createDataSet<dcomplex>("RT/REALTIMECORRELATIONFUNCTION",
                                    {nRTCorr});
  }

}; // RealTimeMultiSlater::createRTDataSets

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::restoreState() {
  CErr("NYI!");
  hsize_t maxPoints = intScheme.tMax / intScheme.deltaT + 1;

  if (savFile.getDims("RT/TIME")[0] != maxPoints)
    CErr("Mismatched requested and saved propagation length!");

  /*
  // Restore time dependent density
  try {
    savFile.readData("RT/TD_1PDM", *propagator_.onePDM);
    savFile.readData("RT/TD_1PDM_ORTHO", *propagator_.onePDMOrtho);
  } catch(...) { }

  // Find last time step that was checkpointed
  double* timeData = memManager_.template malloc<double>(maxPoints);
  savFile.readData("RT/TIME", timeData);
  int offset = *timeData < 1e-10 ? -1 : 0;
  size_t restoreStep = offset + std::distance( timeData,
    std::find_if( timeData+1, timeData+maxPoints,
      [](double x){ return x < 1e-10; }
    )
  );
  memManager_.free(timeData);

  if( printLevel > 0 ) {
    std::cout << "  *** Restoring from step " << restoreStep << " (";
    std::cout << std::setprecision(4) << restoreStep * intScheme.deltaT;
    std::cout << " AU) ***" << std::endl;
  }

  intScheme.restoreStep = restoreStep;
  */

}; // RealTimeMultiSlater::restoreState

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::saveState(EMPerturbation &pert_t) {
  data.Time.push_back(curState.xTime);
  data.Energy.push_back(this->total_energy);
  data.ElecDipole.push_back(Dipole);
  if (pert_t.fields.size() > 0)
    data.ElecDipoleField.push_back(pert_t.getDipoleAmp(Electric));

  // Write to file
  if (savFile.exists()) {

    hsize_t nSteps = 0;
    size_t maxStep =
        (size_t)((intScheme.tMax + intScheme.deltaT / 4) / intScheme.deltaT);

    if ((curState.iStep + 1) % intScheme.iSave == 0 and
        curState.iStep != intScheme.restoreStep)
      nSteps = intScheme.iSave;
    else if (curState.iStep == maxStep) {
      nSteps = (curState.iStep - intScheme.restoreStep) % intScheme.iSave + 1;
      for (auto a : data.Energy) {
        std::cout << a << std::endl;
      }
    }

    hsize_t lastPos = curState.iStep - nSteps + 1;
    hsize_t memLastPos = data.Time.size() - nSteps;

    if (nSteps != 0) {
      if (printLevel > 0)
        std::cout << "  *** Saving data to binary file ***" << std::endl;
      savFile.partialWriteData("RT/TIME", data.Time.data(), {lastPos}, {nSteps},
                               {memLastPos}, {data.Time.size()});
      savFile.partialWriteData("RT/ENERGY", data.Energy.data(), {lastPos},
                               {nSteps}, {memLastPos}, {data.Energy.size()});
      savFile.partialWriteData("RT/LEN_ELEC_DIPOLE", &data.ElecDipole[0][0],
                               {lastPos, 0}, {nSteps, 3}, {memLastPos, 0},
                               {data.Time.size(), 3});

      if (data.ElecDipoleField.size() > 0)
        savFile.partialWriteData(
            "RT/LEN_ELEC_DIPOLE_FIELD", &data.ElecDipoleField[0][0],
            {lastPos, 0}, {nSteps, 3}, {memLastPos, 0}, {data.Time.size(), 3});

      // savFile.safeWriteData("RT/TD_1PDM", *propagator_.onePDM);
      // if ( curState.curStep == PropagationStep::ModifiedMidpoint )
      //   savFile.safeWriteData("RT/TD_1PDM_ORTHO",*DOSav[0]);
      // else
      //   savFile.safeWriteData("RT/TD_1PDM_ORTHO",*propagator_.onePDMOrtho);
    }
  }
}; // RealTimeMultiSlater::saveState

}; // namespace ChronusQ
