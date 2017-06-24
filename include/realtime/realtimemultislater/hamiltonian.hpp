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
void RealTimeMultiSlater<MatsT, IntsT>::formHamiltonian(double t) {
  ProgramTimer::timeOp("Form Hamiltonian", [&]() {
    // Get perturbation for the current time and build a Fock matrix
    EMPerturbation pert_t = pert.getPert(t);

    // Add on the SCF Perturbation
    if (intScheme.includeSCFField)
      for (auto &field : scfPert.fields)
        pert_t.addField(field);

    auto dipAmp = pert_t.getDipoleAmp(Electric);
    double sum = 0.0;
    for (auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      sum += std::abs(dipAmp[iXYZ] - old_amp[iXYZ]);
    }
    if (sum < transform_threshold) {
       return;
    }
    for (auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      old_amp[iXYZ] =  dipAmp[iXYZ];
    }
    auto derived_ref =
        dynamic_cast<MCWaveFunction<MatsT, IntsT> *>(reference_.get());
    if (time_independent_ham) {
      EMPerturbation dummy_field;
      derived_ref->transformInts(dummy_field);
    } else {
      derived_ref->transformInts(pert_t);
    }
  });
};

}; // namespace ChronusQ
