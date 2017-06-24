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
void RealTimeMultiSlater<MatsT, IntsT>::doPropagation() {
  ProgramTimer::tick("Real Time Total");
  printRTHeader();
  // Generate the initial wavefunction to propogate
  curState.curStep = intScheme.integrationAlgorithm;
  genInitialState();
  if (savFile.exists())
    if (restart)
      restoreState();
  size_t maxStep =
      (size_t)((intScheme.tMax + intScheme.deltaT / 4) / intScheme.deltaT);
  bool Start(false);  // Start the RT iterations
  bool Finish(false); // Wrap up the RT iterations
  size_t RealTimeCorrelationFunctionFirstStep =
      (size_t)((this->RealTimeCorrelationFunctionStart + intScheme.deltaT / 2) /
               intScheme.deltaT);
  curState.xTime = intScheme.restoreStep * intScheme.deltaT;
  for (curState.iStep = intScheme.restoreStep; curState.iStep <= maxStep;
       curState.xTime += intScheme.deltaT, curState.iStep++) {
    ProgramTimer::tick("Real Time Iter");
    // Perturbation for the current time
    EMPerturbation pert_t = pert.getPert(curState.xTime);

    // "Start" if this is the first step
    Start = (curState.iStep == intScheme.restoreStep);

    // "Finish" if this is the last step
    Finish = (curState.iStep == maxStep);

    // Determine the step type for the current integration step
    // All steps are same
    curState.curStep = intScheme.integrationAlgorithm;

    // set stepsize
    curState.stepSize = intScheme.deltaT;

    // propogate the CI vector
    propagateWFN(Start, Finish);

    // calculate Dipole Moment
    calculateDipole();
    // Print progress line in the output file
    printRTStep();
    if (this->CIPopFreq != 0 && curState.iStep % this->CIPopFreq == 0) {
      CIPop();
    }
    if (this->RealTimeCorrelationFunctionFreq != 0 &&
        curState.iStep == RealTimeCorrelationFunctionFirstStep) {
      vecManager->allocateCorrelationFunctionMemory();
      // we dont alloc in memory.h cause we dont know if we need it there yet
      // (RealTimeCorrelationFunctionFirstStep is not set)
      RealTimeCorrelationFunction();
    }
    if (this->RealTimeCorrelationFunctionFreq != 0 &&
        curState.iStep > RealTimeCorrelationFunctionFirstStep &&
        (curState.iStep - RealTimeCorrelationFunctionFirstStep) %
                this->RealTimeCorrelationFunctionFreq ==
            0) {
      RealTimeCorrelationFunction();
    }

    // printing after the propagation forward in time
    // because the SSO method creates the imaginary component
    // at the current time during the forward propagation

    // Printing out real-time CI vector
    if (printCIVec)
      CErr("Implement printCIVec");

    // Save data
    saveState(pert_t);

    ProgramTimer::tock("Real Time Iter");

  } // Time loop

  ProgramTimer::tock("Real Time Total");

}; // RealTimeMultiSlater::doPropagation

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::propagateWFN(bool Start, bool Finish) {
  if (Start && Finish) {
    CErr("RealTimeMultiSlater::propagateWFN error. Can't be both starting and "
         "finishing propagation.");
  }
  ProgramTimer::tick("Propagate WFN");
  if (curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    propagateWFN_SSO(Start, Finish);
  } else if (curState.curStep == RealTimeAlgorithm::RTRungeKuttaOrderFour) {
    propagateWFN_RK4(Start, Finish);
  } else {
    CErr("Other propagation algorithms not yet implemented");
  }
  ProgramTimer::tock("Propagate WFN");
}; // RealTimeMultiSlater::propagateWFN
}; // namespace ChronusQ
