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
void RealTimeMultiSlater<MatsT, IntsT>::propagateWFN_RK4(bool Start,
                                                         bool Finish) {
  // Ref: doi:10.1063/1.5126945
  // Could implement low storage RK4 if desired. 10.1016/j.jcp.2009.11.006
  // At minimum, I think this would involve 1 more sigma build
  // Variable time step low storage methods could enable larger time steps which
  // might be worth the expense example CK90 integrator used here
  // doi:10.1021/acs.jctc.2c00490 For now implement standard RK4 (Butcher's
  // tableau): 0   | 1/2 | 1/2 1/2 |  0   1/2 1   |  0    0    1
  // ------------------------
  //     | 1/6  1/3  1/3  1/6

  // |k1>        = -i H(t) |psi(t)>
  // |ktemp>      = |psi(t)> + 1/2dt |k1>
  // |k2>        = -i H(t+dt/2) |ktemp>
  // |ktemp>      = |psi(t)> + 1/2dt |k2>
  // |k3>        = -i H(t+dt/2) |ktemp>
  // |ktemp>      = |psi(t)> + dt |k3>
  // |k4>        = -i H(t+dt) |ktemp>
  // |psi(t+dt)> = |psi(t)> + dt/6|k1> + dt/3|k1> + dt/3|k2> + dt/6|k3>

  // If starting, copy psi from mcwfn obj into |psi> for psi(0+dt)
  // Otherwise copy psi(t+dt) to psi(t) do propagation to form psi(t+dt).
  // properties etc are evaluated on psi(t) copy psi(t+dt) and evaluate
  // properties If finishing, copy psi(t+dt) to psi(t), no propagation just
  // evaluate properties
  auto vecManagerDerived =
      dynamic_cast<RealTimeMultiSlaterVectorManagerRK4<MatsT *> *>(
          vecManager.get());
  auto derived_ref =
      dynamic_cast<MCWaveFunction<MatsT, IntsT> *>(reference_.get());
  size_t NDet = vecManagerDerived->get_vecSize_();
  if (Start) {
    // transform integrals for t
    this->formHamiltonian(curState.xTime);
  }
  if (!Start) {
    RTMS::copy(vecManagerDerived->C_tplusdt, vecManagerDerived->C_t, NDet);
  }

  this->total_energy = derived_ref->reference().molecule().nucRepEnergy +
                       derived_ref->InactEnergy;
  const dcomplex negative_i(0.0, -1.0);
  if (!Finish) {
    // transform integrals for t
    // this->formHamiltonian(curState.xTime); integrals already transformed for
    // current time

    // |k1>        = -i H(t) |psi(t)>
    this->buildSigma(vecManagerDerived->C_t, vecManagerDerived->k1,
                     curState.xTime);
    dcomplex dot_result;
    RTMS::dot(vecManagerDerived->C_t, vecManagerDerived->k1, NDet, dot_result);
    this->total_energy += std::real(dot_result);
    RTMS::scal(vecManagerDerived->k1, NDet, negative_i);

    // |ktemp>      = |psi(t)> + 1/2dt |k1>
    RTMS::copy(vecManagerDerived->C_t, vecManagerDerived->ktemp, NDet);
    RTMS::add(vecManagerDerived->k1, vecManagerDerived->ktemp, NDet,
              0.5 * curState.stepSize);

    // transform integrals for t + 0.5 dt
    this->formHamiltonian(curState.xTime + (0.5 * curState.stepSize));

    // |k2>        = -i H(t+dt/2) |ktemp>
    this->buildSigma(vecManagerDerived->ktemp, vecManagerDerived->k2,
                     curState.xTime + (0.5 * curState.stepSize));
    RTMS::scal(vecManagerDerived->k2, NDet, negative_i);

    // |ktemp>      = |psi(t)> + 1/2dt |k2>
    RTMS::copy(vecManagerDerived->C_t, vecManagerDerived->ktemp, NDet);
    RTMS::add(vecManagerDerived->k2, vecManagerDerived->ktemp, NDet,
              0.5 * curState.stepSize);

    // |k3>        = -i H(t+dt/2) |ktemp>
    this->buildSigma(vecManagerDerived->ktemp, vecManagerDerived->k3,
                     curState.xTime + (0.5 * curState.stepSize));
    RTMS::scal(vecManagerDerived->k3, NDet, negative_i);

    // |ktemp>      = |psi(t)> + dt |k3>
    RTMS::copy(vecManagerDerived->C_t, vecManagerDerived->ktemp, NDet);
    RTMS::add(vecManagerDerived->k3, vecManagerDerived->ktemp, NDet,
              curState.stepSize);

    // transform integrals for t + dt
    this->formHamiltonian(curState.xTime + (curState.stepSize));

    // |k4>        = -i H(t+dt) |ktemp>
    this->buildSigma(vecManagerDerived->ktemp, vecManagerDerived->k4,
                     curState.xTime + (curState.stepSize));
    RTMS::scal(vecManagerDerived->k4, NDet, negative_i);

    // |psi(t+dt)> = |psi(t)> + dt/6|k1> + dt/3|k1> + dt/3|k2> + dt/6|k3>
    RTMS::copy(vecManagerDerived->C_t, vecManagerDerived->C_tplusdt, NDet);
    RTMS::add(vecManagerDerived->k1, vecManagerDerived->C_tplusdt, NDet,
              curState.stepSize / 6.0);
    RTMS::add(vecManagerDerived->k2, vecManagerDerived->C_tplusdt, NDet,
              curState.stepSize / 3.0);
    RTMS::add(vecManagerDerived->k3, vecManagerDerived->C_tplusdt, NDet,
              curState.stepSize / 3.0);
    RTMS::add(vecManagerDerived->k4, vecManagerDerived->C_tplusdt, NDet,
              curState.stepSize / 6.0);
    double norm = 0.0;
    RTMS::normalize(vecManagerDerived->C_tplusdt, NDet, norm);
  } else {
    // if we are finished we need our energy!
    this->buildSigma(vecManagerDerived->C_t, vecManagerDerived->k1,
                     curState.xTime);
    dcomplex dot_result;
    RTMS::dot(vecManagerDerived->C_t, vecManagerDerived->k1, NDet, dot_result);
    this->total_energy += std::real(dot_result);
  }
}

}; // namespace ChronusQ
