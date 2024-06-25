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
void RealTimeMultiSlater<MatsT, IntsT>::RealTimeCorrelationFunction() {
  dcomplex curr_corr(0.0, 0.0);
  if (curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT *>>(
            vecManager);
    size_t NDet = vecManagerDerived->get_vecSize_();
    // (Crealep, Cimagep) dot (Crealt, Cimagt)
    // (a + bi)\dagger dot (c + di)
    // (a dot c + b dot d) + (a dot d - b dot c) i
    double dot_val;
    RTMS::dot(vecManagerDerived->C_real_epsilon, vecManagerDerived->C_real_t,
              NDet, dot_val);
    curr_corr += dot_val;
    RTMS::dot(vecManagerDerived->C_imag_epsilon, vecManagerDerived->C_imag_t,
              NDet, dot_val);
    curr_corr += dot_val;
    RTMS::dot(vecManagerDerived->C_real_epsilon, vecManagerDerived->C_imag_t,
              NDet, dot_val);
    curr_corr += dcomplex(0.0, dot_val);
    RTMS::dot(vecManagerDerived->C_imag_epsilon, vecManagerDerived->C_real_t,
              NDet, dot_val);
    curr_corr -= dcomplex(0.0, dot_val);
  } else {
    CErr("NYI");
  }
  size_t RealTimeCorrelationFunctionFirstStep =
      (size_t)((this->RealTimeCorrelationFunctionStart + intScheme.deltaT / 2) /
               intScheme.deltaT);
  if (savFile.exists()) {
    hsize_t location = (curState.iStep - RealTimeCorrelationFunctionFirstStep) /
                       this->RealTimeCorrelationFunctionFreq;
    hsize_t lastPos = curState.iStep - RealTimeCorrelationFunctionFirstStep;
    hsize_t memLastPos = curState.iStep - RealTimeCorrelationFunctionFirstStep;
    savFile.partialWriteData("RT/REALTIMECORRELATIONFUNCTION", &curr_corr,
                             {location}, {1}, {0}, {1});
  }
}; // RealTimeMultiSlater:: RTCorrFunc

}; // namespace ChronusQ
