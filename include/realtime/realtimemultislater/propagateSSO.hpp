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
void RealTimeMultiSlater<MatsT, IntsT>::propagateWFN_SSO(bool Start,
                                                         bool Finish) {
  // TODO figure out upcast for cc/dmrg wfn
  auto vecManagerDerived = dynamic_cast<RealTimeMultiSlaterVectorManagerSSO<double *> *>(vecManager.get());
  // auto derived_ref = dynamic_cast<MCWaveFunction<double,
  // double>*>(reference_.get());
  auto derived_ref = dynamic_cast<MCWaveFunction<MatsT, IntsT> *>(reference_.get());
  size_t NDet = vecManagerDerived->get_vecSize_();
  // C(t) = p(t) + q(t) i
  // Ref: doi:10.1021/acs.jctc.8b00381
  if (Start) {
    // transform integrals for t0
    this->formHamiltonian(curState.xTime);

    // dq(t0) = H(t0) p(t0)
    this->buildSigma(vecManagerDerived->C_real_t, vecManagerDerived->dC, curState.xTime);
    this->total_energy = derived_ref->reference().molecule().nucRepEnergy + derived_ref->InactEnergy;
    double dot_result;
    RTMS::dot(vecManagerDerived->C_real_t, vecManagerDerived->dC, NDet, dot_result);
    this->total_energy += dot_result;

    // q(t0 + 0.5 dt) = q(t0) - 1/2 dt dq(t0)
    RTMS::copy(vecManagerDerived->C_imag_t, vecManagerDerived->C_imag_tplushalfdt, NDet);
    RTMS::add(vecManagerDerived->dC, vecManagerDerived->C_imag_tplushalfdt, NDet, -0.5 * curState.stepSize);

    // transform integrals for t0 + 0.5 dt
    this->formHamiltonian(curState.xTime + (0.5 * curState.stepSize));

    // dp(t0 + 0.5 dt) = H(t0 + 0.5 dt) q(t0 + 0.5 dt)
    this->buildSigma(vecManagerDerived->C_imag_tplushalfdt, vecManagerDerived->dC, curState.xTime + (0.5 * curState.stepSize));

    // p(t0 + dt) =  p(t0) + dt dp(t0 + 0.5 dt)
    RTMS::copy(vecManagerDerived->C_real_t, vecManagerDerived->C_real_tplusdt, NDet);
    RTMS::add(vecManagerDerived->dC, vecManagerDerived->C_real_tplusdt, NDet, curState.stepSize);
    for (auto i = 0; i < NDet; i++) {
     std::cout << vecManagerDerived->C_imag_tplushalfdt[i] << std::endl;
    }
  } else if (Finish) {
    RTMS::copy(vecManagerDerived->C_imag_tplushalfdt, vecManagerDerived->C_imag_t, NDet);
    RTMS::copy(vecManagerDerived->C_real_tplusdt, vecManagerDerived->C_real_t, NDet);

    // transform integrals for t
    this->formHamiltonian(curState.xTime);

    // dq(t) = + H(t) p(t)
    this->buildSigma(vecManagerDerived->C_real_t, vecManagerDerived->dC, curState.xTime);
    this->total_energy = derived_ref->reference().molecule().nucRepEnergy + derived_ref->InactEnergy;
    double dot_result;
    RTMS::dot(vecManagerDerived->C_real_t, vecManagerDerived->dC, NDet, dot_result);
    this->total_energy += dot_result;

    // q(t) = q(t - 0.5 dt) - 0.5 dt dq(t)
    RTMS::add(vecManagerDerived->dC, vecManagerDerived->C_imag_t, NDet, -0.5 * curState.stepSize);

    this->buildSigma(vecManagerDerived->C_imag_t, vecManagerDerived->dC, curState.xTime);
    dot_result = 0.0;
    RTMS::dot(vecManagerDerived->C_imag_t, vecManagerDerived->dC, NDet, dot_result);
    this->total_energy += dot_result;
  } else {
    RTMS::copy(vecManagerDerived->C_imag_tplushalfdt, vecManagerDerived->C_imag_tminushalfdt, NDet);
    RTMS::copy(vecManagerDerived->C_real_tplusdt, vecManagerDerived->C_real_t, NDet);
    RTMS::fill(vecManagerDerived->C_imag_t, 0.0, NDet);

    // transform integrals for t
    this->formHamiltonian(curState.xTime);

    // dq(t) = + H(t) p(t)
    this->buildSigma(vecManagerDerived->C_real_t, vecManagerDerived->dC, curState.xTime);
    this->total_energy = derived_ref->reference().molecule().nucRepEnergy + derived_ref->InactEnergy;
    double dot_result;
    RTMS::dot(vecManagerDerived->C_real_t, vecManagerDerived->dC, NDet, dot_result);
    this->total_energy += dot_result;

    // q(t + 0.5 dt) = q(t - 0.5 dt) - dt dq(t)
    RTMS::add(vecManagerDerived->dC, vecManagerDerived->C_imag_tplushalfdt, NDet, -curState.stepSize);

    // transform integrals for t + 0.5 dt
    this->formHamiltonian(curState.xTime + (0.5 * curState.stepSize));

    // dp(t + 0.5 dt) = H(t + 0.5 dt) q(t + 0.5 dt)
    this->buildSigma(vecManagerDerived->C_imag_tplushalfdt, vecManagerDerived->dC, curState.xTime + (0.5 * curState.stepSize));

    // p(t + dt) =  p(t) + dt dp(t + 0.5 dt)
    RTMS::add(vecManagerDerived->dC, vecManagerDerived->C_real_tplusdt, NDet, curState.stepSize);

    // q(t) = 0.5 q(t - 0.5 dt) + 0.5 q(t + 0.5 dt)
    RTMS::add(vecManagerDerived->C_imag_tminushalfdt, vecManagerDerived->C_imag_t, NDet, 0.5);
    RTMS::add(vecManagerDerived->C_imag_tplushalfdt, vecManagerDerived->C_imag_t, NDet, 0.5);

    this->buildSigma(vecManagerDerived->C_imag_t, vecManagerDerived->dC, curState.xTime + (0.5 * curState.stepSize));
    dot_result = 0.0;
    RTMS::dot(vecManagerDerived->C_imag_t, vecManagerDerived->dC, NDet, dot_result);
    this->total_energy += dot_result;
  }
}

}; // namespace ChronusQ
