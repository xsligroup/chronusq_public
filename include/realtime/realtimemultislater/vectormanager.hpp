/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <itersolver.hpp>
#include <realtime.hpp>

namespace ChronusQ {

namespace RTMS {
template <typename MatsT>
void copy(std::shared_ptr<SolverVectors<MatsT>> source,
          std::shared_ptr<SolverVectors<MatsT>> dest) {
  // should check that source + dest sizes are the same..
  if (source->size() * source->length() != dest->size() * dest->length()) {
    CErr("Mismatched sizes");
  }
  dest->set_data(0, source->size(), *source, 0);
};

template <typename MatsT, typename T>
void add(std::shared_ptr<SolverVectors<MatsT>> source,
         std::shared_ptr<SolverVectors<MatsT>> dest,
         T factor) {
  size_t shiftY = 0, nVec = 1, shiftX = 0;
  dest->axpy(shiftY, nVec, factor, *source, shiftX);
};

// template <typename MatsT, typename T>
// void add(MatsT* source, std::shared_ptr<SolverVectors<MatsT>> dest, size_t
// , T factor) {
//   size_t shiftY = 0, nVec = 1, shiftX = 0;
//   auto derived_dest= std::dynamic_pointer_cast<RawVectors<MatsT>>(dest);
//   blas::axpy(dest->size() * dest->length(), factor, source, 1,
//   derived_dest->getPtr(shiftY), 1);
// };

template <typename MatsT, typename T>
void fill(std::shared_ptr<SolverVectors<MatsT>> dest, T val) {
  // MatsT MatsT_val = (MatsT) val;
  for (int i = 0; i < dest->length(); i++) {
    for (int j = 0; j < dest->size(); j++) {
      dest->set(i, j, val);
    }
  }
};

template <typename MatsT>
void dot(std::shared_ptr<SolverVectors<MatsT>> source_1,
         std::shared_ptr<SolverVectors<MatsT>> source_2,
         MatsT &result) {
  size_t shiftA = 0, shiftB = 0;
  int64_t m = 1, n = 1, ldc = 1;
  bool conjA = true;
  source_1->dot_product(shiftA, *source_2, shiftB, m, n, &result, ldc, conjA);
};

template <typename MatsT, typename T>
void scal(std::shared_ptr<SolverVectors<MatsT>> source,
          T factor) {
  size_t shift = 0, nVec = 1;
  source->scale(factor, shift, nVec);
};

template <typename MatsT, typename T>
void normalize(std::shared_ptr<SolverVectors<MatsT>> source,
               T &result) {
  size_t shift = 0, nVec = 1;
  result = source->norm2F(shift, nVec);
  source->scale(1.0 / result, shift, nVec);
};
} // namespace RTMS

class RealTimeMultiSlaterVectorManagerBase {
protected:
  double inactiveEnergy_ = 0.0;
  MPI_Comm comm_;

public:
  RealTimeMultiSlaterVectorManagerBase() = delete;
  RealTimeMultiSlaterVectorManagerBase(
      const RealTimeMultiSlaterVectorManagerBase &) = delete;
  RealTimeMultiSlaterVectorManagerBase(
      RealTimeMultiSlaterVectorManagerBase &&) = delete;
  RealTimeMultiSlaterVectorManagerBase(MPI_Comm comm) : comm_(comm){};

  // For handling multislater initial wavefunctions
  MSInitialState initmethod;
  std::vector<std::pair<double, size_t>> init_detail;
  // Take a shared_ptr to the underlying MCWaveFunction object to get the CI
  // Vectors
  virtual void allocateMemory(bool nonhermitian = false) = 0;
  virtual void
  allocateCorrelationFunctionMemory() = 0; // this is a separate function since
                                           // it is optional
  virtual void cleanupMemory() = 0;

  void set_inactiveEnergy_(double in_E) { this->inactiveEnergy_ = in_E; };
  double get_inactiveEnergy_() { return this->inactiveEnergy_; }
};

template <typename MatsT>
class RealTimeMultiSlaterVectorManagerSSO
    : public RealTimeMultiSlaterVectorManagerBase {
public:
  RealTimeMultiSlaterVectorManagerSSO() = delete;
  RealTimeMultiSlaterVectorManagerSSO(
      const RealTimeMultiSlaterVectorManagerSSO &) = delete;
  RealTimeMultiSlaterVectorManagerSSO(RealTimeMultiSlaterVectorManagerSSO &&) =
      delete;
  RealTimeMultiSlaterVectorManagerSSO(MPI_Comm comm, std::function<std::shared_ptr<SolverVectors<MatsT>>(const size_t)>& in_vec_creator)
      : RealTimeMultiSlaterVectorManagerBase(comm), create_vecs(in_vec_creator){};

  std::shared_ptr<SolverVectors<MatsT>> C_real_t;
  std::shared_ptr<SolverVectors<MatsT>> C_real_tplusdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_tminushalfdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_tplushalfdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_t;
  std::shared_ptr<SolverVectors<MatsT>> dC;

  std::function< std::shared_ptr<SolverVectors<MatsT>>(const size_t )> create_vecs;

  std::shared_ptr<SolverVectors<MatsT>>
      C_real_epsilon; // wave functions at t=\epsilon for the accumulation of
                      // the RT correlation function
  std::shared_ptr<SolverVectors<MatsT>> C_imag_epsilon;

  void allocateMemory(bool nonhermitian = false) override {
    // allocating
    C_real_t = create_vecs(1);
    C_real_tplusdt = create_vecs(1);

    if (nonhermitian) {
      CErr("Nonhermitian SSO not implemented!");
    }

    C_imag_tminushalfdt = create_vecs(1);
    C_imag_tplushalfdt = create_vecs(1);
    C_imag_t = create_vecs(1);

    dC = create_vecs(1);
  }

  void allocateCorrelationFunctionMemory() override {
    C_real_epsilon = create_vecs(1);
    C_imag_epsilon = create_vecs(1);
    RTMS::copy(C_real_t, C_real_epsilon );
    RTMS::copy(C_imag_t, C_imag_epsilon );
    // Real Time Correlation Function requires CONJ(C(epsilon) C(epsilon+t) so
    // we scale the imag part by negative 1 here
    // blas::axpy(1, C_imag_t, 1, C_imag_epsilon, 1);
  }

  void cleanupMemory() override {
    // not needed
  }
};

template <typename MatsT>
class RealTimeMultiSlaterVectorManagerRK4
    : public RealTimeMultiSlaterVectorManagerBase {
public:
  RealTimeMultiSlaterVectorManagerRK4() = delete;
  RealTimeMultiSlaterVectorManagerRK4(
      const RealTimeMultiSlaterVectorManagerRK4 &) = delete;
  RealTimeMultiSlaterVectorManagerRK4(RealTimeMultiSlaterVectorManagerRK4 &&) =
      delete;
  RealTimeMultiSlaterVectorManagerRK4(MPI_Comm comm, std::function<std::shared_ptr<SolverVectors<MatsT>>(const size_t)>& in_vec_creator)
      : RealTimeMultiSlaterVectorManagerBase(comm), create_vecs(in_vec_creator) {};

  std::shared_ptr<SolverVectors<MatsT>> C_t;
  std::shared_ptr<SolverVectors<MatsT>> C_tplusdt;

  std::shared_ptr<SolverVectors<MatsT>> left_C_t;
  std::shared_ptr<SolverVectors<MatsT>> left_C_tplusdt;

  std::shared_ptr<SolverVectors<MatsT>> k1;
  std::shared_ptr<SolverVectors<MatsT>> k2;
  std::shared_ptr<SolverVectors<MatsT>> k3;
  std::shared_ptr<SolverVectors<MatsT>> k4;

  std::shared_ptr<SolverVectors<MatsT>> ktemp;

  std::shared_ptr<SolverVectors<MatsT>>
      C_epsilon; // wave functions at t=\epsilon for the accumulation of the
                 // RT correlation function

  std::function< std::shared_ptr<SolverVectors<MatsT>>(const size_t )> create_vecs;

  void allocateMemory(bool nonhermitian = false) override {
    C_t = create_vecs(1);
    C_tplusdt = create_vecs(1);
    if (nonhermitian) {
      left_C_t = create_vecs(1);
      left_C_tplusdt = create_vecs(1);
    }
    k1 = create_vecs(1); 
    k2 = create_vecs(1); 
    k3 = create_vecs(1); 
    k4 = create_vecs(1); 
    ktemp = create_vecs(1);
  }

  void allocateCorrelationFunctionMemory() override {
    C_epsilon = create_vecs(1);
    RTMS::copy(C_t, C_epsilon);
  }

  void cleanupMemory() override {
    // not needed
  }
};

}; // namespace ChronusQ
