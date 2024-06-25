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
#include <realtime.hpp>

namespace ChronusQ {

namespace RTMS {
template <typename oper_t>
void copy(oper_t source, oper_t dest, size_t vecSize_) {
  std::copy_n(source, vecSize_, dest);
};

template <typename oper_t, typename T>
void add(oper_t source, oper_t dest, size_t vecSize_, T factor) {
  blas::axpy(vecSize_, factor, source, 1, dest, 1);
};

template <typename oper_t, typename T>
void fill(oper_t dest, T val, size_t vecSize_) {
  std::fill_n(dest, vecSize_, val);
};

template <typename oper_t, typename T>
void dot(oper_t source_1, oper_t source_2, size_t vecSize_, T &result) {
  result = blas::dot(vecSize_, source_1, 1, source_2, 1);
};

template <typename oper_t, typename T>
void scal(oper_t source, size_t vecSize_, T factor) {
  blas::scal(vecSize_, factor, source, 1);
};

template <typename oper_t, typename T>
void normalize(oper_t source, size_t vecSize_, T &result) {
  result = Normalize(vecSize_, source, 1);
};
} // namespace RTMS

class RealTimeMultiSlaterVectorManagerBase {
protected:
  size_t vecSize_;
  double inactiveEnergy_ = 0.0;

public:
  // For handling multislater initial wavefunctions
  MSInitialState initmethod;
  std::vector<std::pair<double, size_t>> init_detail;
  // Take a shared_ptr to the underlying MCWaveFunction object to get the CI
  // Vectors
  virtual void allocateMemory() = 0;
  virtual void allocateCorrelationFunctionMemory() = 0; // this is a separate function since it is optional
  virtual void cleanupMemory() = 0;

  void set_inactiveEnergy_(double in_E) { this->inactiveEnergy_ = in_E; };
  double get_inactiveEnergy_() { return this->inactiveEnergy_; }
  void set_vecSize_(int in_size) { this->vecSize_ = in_size; };
  size_t get_vecSize_() { return this->vecSize_; }
};

template <typename oper_t>
class RealTimeMultiSlaterVectorManagerSSO
    : public RealTimeMultiSlaterVectorManagerBase {
public:
  oper_t C_real_t;
  oper_t C_real_tplusdt;
  oper_t C_imag_tminushalfdt;
  oper_t C_imag_tplushalfdt;
  oper_t C_imag_t;
  oper_t dC;

  oper_t C_real_epsilon; // wave functions at t=\epsilon for the accumulation of
                         // the RT correlation function
  oper_t C_imag_epsilon;

  template <typename MatsT, typename IntsT>
  void buildInitCIVec(std::shared_ptr<MCWaveFunction<MatsT, IntsT>>);

  void allocateMemory() override {
    // allocating
    C_real_t = CQMemManager::get().malloc<double>(this->get_vecSize_());
    C_real_tplusdt = CQMemManager::get().malloc<double>(this->get_vecSize_());

    C_imag_tminushalfdt =
        CQMemManager::get().malloc<double>(this->get_vecSize_());
    C_imag_tplushalfdt =
        CQMemManager::get().malloc<double>(this->get_vecSize_());
    C_imag_t = CQMemManager::get().malloc<double>(this->get_vecSize_());

    dC = CQMemManager::get().malloc<double>(this->get_vecSize_());

    RTMS::fill(C_real_t, 0.0, this->get_vecSize_());
    RTMS::fill(C_real_tplusdt, 0.0, this->get_vecSize_());
    RTMS::fill(C_imag_tminushalfdt, 0.0, this->get_vecSize_());
    RTMS::fill(C_imag_tplushalfdt, 0.0, this->get_vecSize_());
    RTMS::fill(C_imag_t, 0.0, this->get_vecSize_()),
    RTMS::fill(dC, 0.0, this->get_vecSize_());
  }

  void allocateCorrelationFunctionMemory() override {
    C_real_epsilon = CQMemManager::get().malloc<double>(this->get_vecSize_());
    C_imag_epsilon = CQMemManager::get().malloc<double>(this->get_vecSize_());
    // std::fill_n(C_real_epsilon, this->get_vecSize_(), 0.0);
    const auto Nelem = this->get_vecSize_();
    RTMS::fill(C_imag_epsilon, 0.0, Nelem);
    RTMS::copy(C_real_t, C_real_epsilon, Nelem);
    RTMS::copy(C_imag_t, C_imag_epsilon, Nelem);
    // Real Time Correlation Function requires CONJ(C(epsilon) C(epsilon+t) so
    // we scale the imag part by negative 1 here
    // blas::axpy(this->get_vecSize_(), 1, C_imag_t, 1, C_imag_epsilon, 1);
  }

  void cleanupMemory() override {
    // free mem
    CQMemManager::get().free(C_real_t);
    CQMemManager::get().free(C_real_tplusdt);

    CQMemManager::get().free(C_imag_t);
    CQMemManager::get().free(C_imag_tminushalfdt);
    CQMemManager::get().free(C_imag_tplushalfdt);

    CQMemManager::get().free(dC);
    if (C_real_epsilon)
      CQMemManager::get().free(C_real_epsilon);
    if (C_imag_epsilon)
      CQMemManager::get().free(C_imag_epsilon);
  }
};

template <typename oper_t>
class RealTimeMultiSlaterVectorManagerRK4
    : public RealTimeMultiSlaterVectorManagerBase {
public:
  oper_t C_t;
  oper_t C_tplusdt;

  oper_t k1;
  oper_t k2;
  oper_t k3;
  oper_t k4;

  oper_t ktemp;

  oper_t C_epsilon; // wave functions at t=\epsilon for the accumulation of the
                    // RT correlation function

  template <typename MatsT, typename IntsT>
  void buildInitCIVec(std::shared_ptr<MCWaveFunction<MatsT, IntsT>>);

  void allocateMemory() override {
    // allocating
    C_t = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    C_tplusdt = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    k1 = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    k2 = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    k3 = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    k4 = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    ktemp = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());

    const auto Nelem = this->get_vecSize_();
    RTMS::fill(C_t, 0.0, Nelem);
    RTMS::fill(C_tplusdt, 0.0, Nelem);
    RTMS::fill(k1, 0.0, Nelem);
    RTMS::fill(k2, 0.0, Nelem);
    RTMS::fill(k3, 0.0, Nelem);
    RTMS::fill(k4, 0.0, Nelem);
    RTMS::fill(ktemp, 0.0, Nelem);
  }

  void allocateCorrelationFunctionMemory() override {
    C_epsilon = CQMemManager::get().malloc<dcomplex>(this->get_vecSize_());
    const auto Nelem = this->get_vecSize_();
    RTMS::copy(C_t, C_epsilon, Nelem);
  }

  void cleanupMemory() override {
    // free mem
    CQMemManager::get().free(C_t);
    CQMemManager::get().free(C_tplusdt);

    CQMemManager::get().free(k1);
    CQMemManager::get().free(k2);
    CQMemManager::get().free(k3);
    CQMemManager::get().free(k4);

    CQMemManager::get().free(ktemp);
    if (C_epsilon)
      CQMemManager::get().free(C_epsilon);
  }
};

}; // namespace ChronusQ
