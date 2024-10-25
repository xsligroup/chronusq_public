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
#include <itersolver.hpp>

namespace ChronusQ {

namespace RTMS {
template <typename MatsT>
void copy(std::shared_ptr<SolverVectors<MatsT>> source, std::shared_ptr<SolverVectors<MatsT>> dest, size_t vecSize_) {
  auto derived_source = std::dynamic_pointer_cast<RawVectors<MatsT>>(source);
  auto derived_dest= std::dynamic_pointer_cast<RawVectors<MatsT>>(dest);
  //should check that source + dest sizes are the same..
  if (source->size() * source->length() != dest->size() * dest->length()){
    CErr("Mismatched sizes");
  }
  std::copy_n(derived_source->getPtr(), source->size() * source->length(), derived_dest->getPtr());
};

template <typename MatsT, typename T>
void add(std::shared_ptr<SolverVectors<MatsT>> source, std::shared_ptr<SolverVectors<MatsT>> dest, size_t vecSize_, T factor) {
  size_t shiftY = 0, nVec = 1, shiftX = 0;
  //dest->axpy(shiftY, nVec, factor, *source, shiftX);
  auto derived_source= std::dynamic_pointer_cast<RawVectors<MatsT>>(source);
  auto derived_dest= std::dynamic_pointer_cast<RawVectors<MatsT>>(dest);
  blas::axpy(dest->size() * dest->length(), factor, derived_source->getPtr(shiftX), 1, derived_dest->getPtr(shiftY), 1);
};

template <typename MatsT, typename T>
void add(MatsT* source, std::shared_ptr<SolverVectors<MatsT>> dest, size_t vecSize_, T factor) {
  size_t shiftY = 0, nVec = 1, shiftX = 0;
  auto derived_dest= std::dynamic_pointer_cast<RawVectors<MatsT>>(dest);
  blas::axpy(dest->size() * dest->length(), factor, source, 1, derived_dest->getPtr(shiftY), 1);
};


template <typename MatsT, typename T>
void fill(std::shared_ptr<SolverVectors<MatsT>> dest, T val, size_t vecSize_) {
  auto derived_dest= std::dynamic_pointer_cast<RawVectors<MatsT>>(dest);
  MatsT MatsT_val = (MatsT) val;
  std::fill_n(derived_dest->getPtr(), dest->length() * dest->size(), MatsT_val);
};

template <typename MatsT>
void dot(std::shared_ptr<SolverVectors<MatsT>> source_1, std::shared_ptr<SolverVectors<MatsT>> source_2, size_t vecSize_, MatsT &result) {
    size_t shiftA = 0, shiftB =0;
    int64_t m =1, n = 1, ldc=1;
    bool conjA = true;
    auto derived_source_1= std::dynamic_pointer_cast<RawVectors<MatsT>>(source_1);
    auto derived_source_2= std::dynamic_pointer_cast<RawVectors<MatsT>>(source_2);
    MatsT *val;
    //derived_source_1->dot_product(shiftA, *derived_source_2, shiftB, m, n, val, ldc, conjA);
    //result = *val;
    result = blas::dot(vecSize_, derived_source_1->getPtr(), 1, derived_source_2->getPtr(), 1);
};

template <typename MatsT, typename T>
void scal(std::shared_ptr<SolverVectors<MatsT>> source, size_t vecSize_, T factor) {
  size_t shift=0, nVec=1;
  source->scale(factor, shift, nVec);
};

template <typename MatsT, typename T>
void normalize(std::shared_ptr<SolverVectors<MatsT>> source, size_t vecSize_, T &result) {
  size_t shift=0, nVec=1;
  result = source->norm2F(shift, nVec);
  auto derived_source= std::dynamic_pointer_cast<RawVectors<MatsT>>(source);
  source->scale(1.0/result, shift, nVec);
};
} // namespace RTMS

class RealTimeMultiSlaterVectorManagerBase {
protected:
  size_t vecSize_;
  double inactiveEnergy_ = 0.0;
  MPI_Comm comm_;

public:
    RealTimeMultiSlaterVectorManagerBase()                 = delete;
    RealTimeMultiSlaterVectorManagerBase(const RealTimeMultiSlaterVectorManagerBase &) = delete;
    RealTimeMultiSlaterVectorManagerBase(RealTimeMultiSlaterVectorManagerBase &&)      = delete;
    RealTimeMultiSlaterVectorManagerBase(MPI_Comm comm ) : comm_(comm)  {};

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

template <typename MatsT>
class RealTimeMultiSlaterVectorManagerSSO
    : public RealTimeMultiSlaterVectorManagerBase {
public:
    RealTimeMultiSlaterVectorManagerSSO()                 = delete;
    RealTimeMultiSlaterVectorManagerSSO(const RealTimeMultiSlaterVectorManagerSSO &) = delete;
    RealTimeMultiSlaterVectorManagerSSO(RealTimeMultiSlaterVectorManagerSSO &&)      = delete;
    RealTimeMultiSlaterVectorManagerSSO(MPI_Comm comm ) : RealTimeMultiSlaterVectorManagerBase(comm)  {};

  std::shared_ptr<SolverVectors<MatsT>> C_real_t;
  std::shared_ptr<SolverVectors<MatsT>> C_real_tplusdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_tminushalfdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_tplushalfdt;
  std::shared_ptr<SolverVectors<MatsT>> C_imag_t;
  std::shared_ptr<SolverVectors<MatsT>> dC;

  std::shared_ptr<SolverVectors<MatsT>> C_real_epsilon; // wave functions at t=\epsilon for the accumulation of
                         // the RT correlation function
  std::shared_ptr<SolverVectors<MatsT>> C_imag_epsilon;

  template <typename IntsT>
  void buildInitCIVec(std::shared_ptr<MCWaveFunction<MatsT, IntsT>>);

  void allocateMemory() override {
    // allocating
    C_real_t       = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);
    C_real_tplusdt = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);

    C_imag_tminushalfdt = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);
    C_imag_tplushalfdt = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);
    C_imag_t = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);

    dC       = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);

    // RTMS::fill(C_real_t, 0.0, this->get_vecSize_());
    // RTMS::fill(C_real_tplusdt, 0.0, this->get_vecSize_());
    // RTMS::fill(C_imag_tminushalfdt, 0.0, this->get_vecSize_());
    // RTMS::fill(C_imag_tplushalfdt, 0.0, this->get_vecSize_());
    // RTMS::fill(C_imag_t, 0.0, this->get_vecSize_()),
    // RTMS::fill(dC, 0.0, this->get_vecSize_());
  }

  void allocateCorrelationFunctionMemory() override {
    auto derived_C_real_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(C_real_t);
    auto derived_C_imag_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(C_imag_t);
    C_real_epsilon = std::make_shared<RawVectors<MatsT>>(*derived_C_real_t);
    C_imag_epsilon = std::make_shared<RawVectors<MatsT>>(*derived_C_imag_t);
    // std::fill_n(C_real_epsilon, this->get_vecSize_(), 0.0);
    const auto Nelem = this->get_vecSize_();
    // RTMS::fill(C_imag_epsilon, 0.0, Nelem);
    RTMS::copy(C_real_t, C_real_epsilon, Nelem);
    RTMS::copy(C_imag_t, C_imag_epsilon, Nelem);
    // Real Time Correlation Function requires CONJ(C(epsilon) C(epsilon+t) so
    // we scale the imag part by negative 1 here
    // blas::axpy(this->get_vecSize_(), 1, C_imag_t, 1, C_imag_epsilon, 1);
  }

  void cleanupMemory() override {
    // free mem
    //CQMemManager::get().free(C_real_t);
    //CQMemManager::get().free(C_real_tplusdt);

    //CQMemManager::get().free(C_imag_t);
    //CQMemManager::get().free(C_imag_tminushalfdt);
    //CQMemManager::get().free(C_imag_tplushalfdt);

    //CQMemManager::get().free(dC);
    //if (C_real_epsilon)
    //  CQMemManager::get().free(C_real_epsilon);
    //if (C_imag_epsilon)
    //  CQMemManager::get().free(C_imag_epsilon);
  }
};

template <typename MatsT>
class RealTimeMultiSlaterVectorManagerRK4
    : public RealTimeMultiSlaterVectorManagerBase {
public:
    RealTimeMultiSlaterVectorManagerRK4()                 = delete;
    RealTimeMultiSlaterVectorManagerRK4(const RealTimeMultiSlaterVectorManagerRK4 &) = delete;
    RealTimeMultiSlaterVectorManagerRK4(RealTimeMultiSlaterVectorManagerRK4 &&)      = delete;
    RealTimeMultiSlaterVectorManagerRK4(MPI_Comm comm ) : RealTimeMultiSlaterVectorManagerBase(comm)  {};

  std::shared_ptr<SolverVectors<MatsT>> C_t;
  std::shared_ptr<SolverVectors<MatsT>> C_tplusdt;

  std::shared_ptr<SolverVectors<MatsT>> k1;
  std::shared_ptr<SolverVectors<MatsT>> k2;
  std::shared_ptr<SolverVectors<MatsT>> k3;
  std::shared_ptr<SolverVectors<MatsT>> k4;

  std::shared_ptr<SolverVectors<MatsT>> ktemp;

  std::shared_ptr<SolverVectors<MatsT>> C_epsilon; // wave functions at t=\epsilon for the accumulation of the
                    // RT correlation function

  template <typename IntsT>
  void buildInitCIVec(std::shared_ptr<MCWaveFunction<MatsT, IntsT>>);

  void allocateMemory() override {
    // allocating
    C_t       = std::make_shared<RawVectors<MatsT>>(      this->comm_, this->get_vecSize_(), 1);
    C_tplusdt       = std::make_shared<RawVectors<MatsT>>(this->comm_, this->get_vecSize_(), 1);
    k1       = std::make_shared<RawVectors<MatsT>>(       this->comm_, this->get_vecSize_(), 1);
    k2       = std::make_shared<RawVectors<MatsT>>(       this->comm_, this->get_vecSize_(), 1);
    k3      = std::make_shared<RawVectors<MatsT>>(        this->comm_, this->get_vecSize_(), 1);
    k4       = std::make_shared<RawVectors<MatsT>>(       this->comm_, this->get_vecSize_(), 1);
    ktemp      = std::make_shared<RawVectors<MatsT>>(     this->comm_, this->get_vecSize_(), 1);

    const auto Nelem = this->get_vecSize_();
    // RTMS::fill(C_t, 0.0, Nelem);
    // RTMS::fill(C_tplusdt, 0.0, Nelem);
    // RTMS::fill(k1, 0.0, Nelem);
    // RTMS::fill(k2, 0.0, Nelem);
    // RTMS::fill(k3, 0.0, Nelem);
    // RTMS::fill(k4, 0.0, Nelem);
    // RTMS::fill(ktemp, 0.0, Nelem);
  }

  void allocateCorrelationFunctionMemory() override {
    auto derived_C_t = std::dynamic_pointer_cast<RawVectors<MatsT>>(C_t);
    C_epsilon = std::make_shared<RawVectors<MatsT>>(*derived_C_t);
    const auto Nelem = this->get_vecSize_();
  }

  void cleanupMemory() override {
    // free mem
    // CQMemManager::get().free(C_t);
    // CQMemManager::get().free(C_tplusdt);

    // CQMemManager::get().free(k1);
    // CQMemManager::get().free(k2);
    // CQMemManager::get().free(k3);
    // CQMemManager::get().free(k4);

    // CQMemManager::get().free(ktemp);
    // if (C_epsilon)
    //   CQMemManager::get().free(C_epsilon);
  }
};

}; // namespace ChronusQ
