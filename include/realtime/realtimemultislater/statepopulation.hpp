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
void RealTimeMultiSlater<MatsT, IntsT>::CIPop() {
  std::vector<double> populations;
  std::vector<MatsT> overlaps;
  if (curState.curStep == RealTimeAlgorithm::RTSymplecticSplitOperator) {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerSSO<MatsT *>>(
            vecManager);
    auto derived_ref =
        std::dynamic_pointer_cast<MCWaveFunction<MatsT, IntsT>>(reference_);
    derived_ref->computeOverlaps(vecManagerDerived->C_real_t, overlaps);
    populations.resize(overlaps.size());
    for (auto ovlp_i = 0; ovlp_i < overlaps.size(); ovlp_i++) {
      populations[ovlp_i] += std::pow(std::abs(overlaps[ovlp_i]), 2);
    }
    overlaps.clear();
    derived_ref->computeOverlaps(vecManagerDerived->C_imag_t, overlaps);
    for (auto ovlp_i = 0; ovlp_i < overlaps.size(); ovlp_i++) {
      populations[ovlp_i] += std::pow(std::abs(overlaps[ovlp_i]), 2);
    }
  } else {
    auto vecManagerDerived =
        std::dynamic_pointer_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT *>>(
            vecManager);
    auto derived_ref =
        std::dynamic_pointer_cast<MCWaveFunction<MatsT, IntsT>>(reference_);
    derived_ref->computeOverlaps(vecManagerDerived->C_t, overlaps);
    populations.resize(overlaps.size());
    for (auto ovlp_i = 0; ovlp_i < overlaps.size(); ovlp_i++) {
      populations[ovlp_i] += std::pow(std::abs(overlaps[ovlp_i]), 2);
    }
  }

  if (savFile.exists()) {
    size_t fullDim = populations.size();
    hsize_t location = curState.iStep / this->CIPopFreq;
    savFile.partialWriteData("RT/CIPOPULATION", populations.data(),
                             {location, 0}, {1, fullDim}, {0, 0}, {1, fullDim});
  }

}; // RealTimeMultiSlater:: CIPop

}; // namespace ChronusQ
