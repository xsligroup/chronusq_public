/*
*  This file is part of the Chronus Quantum (ChronusQ) software package
*
*  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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

#include <coupledcluster/EOMCCSDVector.hpp>

namespace ChronusQ{

  template <typename MatsT>
  EOMCCSDVector<MatsT>::EOMCCSDVector(char vLabel, char oLabel): vLabel_(vLabel), oLabel_(oLabel) {
    initialize();
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>::EOMCCSDVector(char vLabel, char oLabel, const TArray &V1, const TArray &V2)
      : vLabel_(vLabel), oLabel_(oLabel) {
    initialize();
    assign(V1, V2);
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>::EOMCCSDVector(char vLabel, char oLabel, MatsT V0, const TArray &V1, const TArray &V2)
      : vLabel_(vLabel), oLabel_(oLabel) {
    initialize();
    assign(V0, V1, V2);
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>::EOMCCSDVector(const EOMCCSDVector<MatsT> &other) {
    initialize();
    operator=(other);
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>::EOMCCSDVector(EOMCCSDVector<MatsT> &&other) {
    swap(other);
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>::~EOMCCSDVector() {
    if (V1_) TAManager::get().free(std::string({vLabel_, oLabel_}), std::move(V1_));
    if (V2_) TAManager::get().free(std::string({vLabel_, vLabel_, oLabel_, oLabel_}), std::move(V2_));
  }

  template <typename MatsT>
  EOMCCSDVector<MatsT>& EOMCCSDVector<MatsT>::operator=(const EOMCCSDVector<MatsT>& other) {
    if (this != &other) {
      V0_ = other.V0_;
      V1_("p,q") = other.V1_("p,q");
      V2_("p,q,r,s") = other.V2_("p,q,r,s");
    }
    return *this;
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::swap(EOMCCSDVector<MatsT> &other) {
    std::swap(V0_, other.V0_);
    std::swap(V1_, other.V1_);
    std::swap(V2_, other.V2_);
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::assign(MatsT V0, const TArray &V1, const TArray &V2) {

    V0_ = V0;
    assign(V1, V2);
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::assign(const TArray &V1, const TArray &V2) {

    oneBody()("p,q") = V1("p,q");
    twoBody()("p,q,r,s") = V2("p,q,r,s");
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::initialize() {

    V0_ = 0.0;

    if (not V1_.is_initialized()){
      V1_ = TAManager::get().malloc<MatsT>(std::string({vLabel_, oLabel_}));
    }

    if (not V2_.is_initialized()){
      V2_ = TAManager::get().malloc<MatsT>(std::string({vLabel_, vLabel_, oLabel_, oLabel_}));
    }
  }

  template <typename MatsT>
  MatsT EOMCCSDVector<MatsT>::dot(const EOMCCSDVector<MatsT> &other, bool conjA) const {
    MatsT twoBodyDot = 0.0, oneBodyDot = 0.0, zeroBodyDot = 0.0;
    if (conjA) {
      twoBodyDot = V2_("a,b,i,j").inner_product(other.V2_("a,b,i,j")).get();
      oneBodyDot = V1_("a,i").inner_product(other.V1_("a,i")).get();
      zeroBodyDot = SmartConj(V0_) * other.V0_;

    } else {
      twoBodyDot = V2_("a,b,i,j").dot(other.V2_("a,b,i,j")).get();
      oneBodyDot = V1_("a,i").dot(other.V1_("a,i")).get();
      zeroBodyDot = V0_ * other.V0_;
    }
    TA::get_default_world().gop.fence();

    return zeroBodyDot + oneBodyDot + 0.25 * twoBodyDot;
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::axpy(MatsT alpha, const EOMCCSDVector<MatsT> &X) {
    V0_ += alpha * X.V0_;
    V1_("a,i") += alpha * X.V1_("a,i");
    V2_("a,b,i,j") += alpha * X.V2_("a,b,i,j");
  }

  template <typename MatsT>
  double EOMCCSDVector<MatsT>::norm() const {
    return std::sqrt(std::real(dot(*this)));
  }

  template <typename MatsT>
  double EOMCCSDVector<MatsT>::absmax() const {
    return std::max({std::abs(V0_), TA::abs_max(V1_).get(), TA::abs_max(V2_).get()});
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::scale(MatsT factor) {
    V0_ *= factor;
    V1_("a,i") = factor * V1_("a,i");
    V2_("a,b,i,j") = factor * V2_("a,b,i,j");
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::conjugate() {
    if (std::is_same<MatsT, double>::value) return;
    V0_ = std::conj(V0_);
    V1_("a,i") = conj(V1_("a,i"));
    V2_("a,b,i,j") = conj(V2_("a,b,i,j"));
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::normalize() {
    scale(1.0/norm());
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::enforceTwoBodySymmetry() {
    V2_("a,b,i,j") -= V2_("b,a,i,j");
    V2_("a,b,i,j") -= V2_("a,b,j,i");
    V2_("a,b,i,j") = 0.25 * V2_("a,b,i,j");
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::projectOut(const EOMCCSDVector<MatsT> &other, bool normalized) {
    MatsT coef = dot(other);
    if (not normalized)
      coef /= other.dot(other);
    axpy(-coef, other);
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::setElem(size_t idx, MatsT elem) {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    if (idx < nO * nV) {
      setOneBodyElem(idx % nV, idx / nV, elem);
    } else {
      idx -= nO * nV;
      const size_t nV2shift = nV * (nV - 1) / 2;
      size_t ab = idx % nV2shift, ij = idx / nV2shift;
      size_t b = static_cast<size_t>(sqrt(2*ab + 0.25) + 0.5);
      size_t a = ab - b * (b-1) / 2;
      size_t j = static_cast<size_t>(sqrt(2*ij + 0.25) + 0.5);
      size_t i = ij - j * (j-1) / 2;
      setTwoBodyElem(a, b, i, j, elem);
    }
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::setZeroBodyElem(MatsT elem) {
    V0_ = elem;
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::setOneBodyElem(size_t a, size_t i, MatsT elem) {
    for(auto it = oneBody().begin(); it != oneBody().end(); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      if (lobound[0] <= a and a < upbound[0]
          and lobound[1] <= i and i < upbound[1]) {
        std::vector<std::size_t> x{a,i};
        tile[x] = elem;
      }

    }
    TA::get_default_world().gop.fence();
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::setTwoBodyElem(size_t a, size_t b, size_t i, size_t j, MatsT elem) {
    for(auto it = twoBody().begin(); it != twoBody().end(); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      if (lobound[0] <= a and a < upbound[0]
          and lobound[1] <= b and b < upbound[1]
          and lobound[2] <= i and i < upbound[2]
          and lobound[3] <= j and j < upbound[3]) {
        std::vector<std::size_t> x{a,b,i,j};
        tile[x] = elem;
      }

      if (lobound[0] <= b and b < upbound[0]
          and lobound[1] <= a and a < upbound[1]
          and lobound[2] <= i and i < upbound[2]
          and lobound[3] <= j and j < upbound[3]) {
        std::vector<std::size_t> x{b,a,i,j};
        tile[x] = -elem;
      }

      if (lobound[0] <= a and a < upbound[0]
          and lobound[1] <= b and b < upbound[1]
          and lobound[2] <= j and j < upbound[2]
          and lobound[3] <= i and i < upbound[3]) {
        std::vector<std::size_t> x{a,b,j,i};
        tile[x] = -elem;
      }

      if (lobound[0] <= b and b < upbound[0]
          and lobound[1] <= a and a < upbound[1]
          and lobound[2] <= j and j < upbound[2]
          and lobound[3] <= i and i < upbound[3]) {
        std::vector<std::size_t> x{b,a,j,i};
        tile[x] = elem;
      }

    }
    TA::get_default_world().gop.fence();
  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::toRaw(MatsT *raw, bool includeZeroBody) const {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    std::fill_n(raw, length(includeZeroBody), 0.0);
    MatsT *raw1 = raw + (includeZeroBody ? 1 : 0);

    for(auto it = std::begin(V1_); it != std::end(V1_); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          raw1[x[0] + x[1] * nV] = tile[x];
    }

    const size_t oneBodyShift = nO * nV;
    const size_t nV2shift = nV * (nV - 1) / 2;

    for(auto it = std::begin(V2_); it != std::end(V2_); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = std::max(x[0]+1, static_cast<std::size_t>(lobound[1])); x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = std::max(x[2]+1, static_cast<std::size_t>(lobound[3])); x[3] < upbound[3]; ++x[3]) {
              raw1[oneBodyShift + x[0] + x[1] * (x[1] - 1) / 2 + (x[2] + x[3] * (x[3] - 1) / 2) * nV2shift] = tile[x];
            }
    }

    TA::get_default_world().gop.fence();

    TA::get_default_world().gop.template reduce(raw1, length(false), std::plus<MatsT>());

    if (includeZeroBody)
      raw[0] = V0_;

  }

  template <typename MatsT>
  template <typename IntsT>
  void EOMCCSDVector<MatsT>::toRaw(MatsT *raw, const EOMCCSD<MatsT,IntsT> &eom, bool includeZeroBody) const {

    std::fill_n(raw, length(includeZeroBody), 0.0);
    MatsT *raw1 = raw + (includeZeroBody ? 1 : 0);

    for(auto it = std::begin(V1_); it != std::end(V1_); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          size_t idx = eom.CVStoCompoundS(x[0], x[1]);
          if (eom.CVSisInBound(idx))
            raw1[idx] = tile[x];
        }
    }

    const size_t oneBodyShift = eom.CVSoneBodySize();
    MatsT *raw2 = raw1 + oneBodyShift;

    for(auto it = std::begin(V2_); it != std::end(V2_); ++it) {

      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = std::max(x[0]+1, static_cast<std::size_t>(lobound[1])); x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = std::max(x[2]+1, static_cast<std::size_t>(lobound[3])); x[3] < upbound[3]; ++x[3]) {
              size_t idx = eom.CVStoCompoundD(x[0], x[1], x[2], x[3]);
              if (eom.CVSisInBound(idx))
                raw2[idx] = tile[x];
            }
    }

    TA::get_default_world().gop.fence();

    TA::get_default_world().gop.template reduce(raw1, length(false), std::plus<MatsT>());

    if (includeZeroBody)
      raw[0] = V0_;

  }

  template <typename MatsT>
  template <typename IntsT>
  void EOMCCSDVector<MatsT>::fromRaw(const MatsT *raw, const EOMCCSD<MatsT,IntsT> &eom, bool hasZeroBody) {

    const MatsT *raw1 = raw;
    if (hasZeroBody) {
      raw1++;
      V0_ = raw[0];
    }

    TA::foreach_inplace(V1_, [&eom, raw1](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          size_t idx = eom.CVStoCompoundS(x[0], x[1]);
          tile[x] = eom.CVSisInBound(idx) ? raw1[idx] : 0.0;
        }
    });

    const size_t oneBodyShift = eom.CVSoneBodySize();
    const MatsT *raw2 = raw1 + oneBodyShift;

    TA::foreach_inplace(V2_, [&eom, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0,0};
      size_t a,b,i,j;
      double signABIJ;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              a = x[0];
              b = x[1];
              i = x[2];
              j = x[3];
              signABIJ = eom.signD(a,b,i,j);
              if (signABIJ == 0.0)
                tile[x] = 0.0;
              else {
                size_t idx = eom.CVStoCompoundD(a, b, i, j);
                tile[x] = eom.CVSisInBound(idx) ? signABIJ * raw2[idx] : 0.0;
              }
            }
    });

    TA::get_default_world().gop.fence();

  }

  template <typename MatsT>
  void EOMCCSDVector<MatsT>::fromRaw(const MatsT *raw, bool hasZeroBody) {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    const MatsT *raw1 = raw;
    if (hasZeroBody) {
      raw1++;
      V0_ = raw[0];
    }

    TA::foreach_inplace(V1_, [this, raw1, nV](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          tile[x] = raw1[x[0] + x[1] * nV];
    });

    const size_t oneBodyShift = nO * nV;
    const size_t nV2shift = nV * (nV - 1) / 2;
    const MatsT *raw2 = raw1 + oneBodyShift;

    TA::foreach_inplace(V2_, [raw2, nV2shift](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0,0};
      size_t a,b,i,j;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          double signAB = 1.0;
          a = x[0];
          b = x[1];
          if (a > b) {
            std::swap(a,b);
            signAB = -1.0;
          } else if (a == b) {
            for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
              for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                tile[x] = 0.0;
              }
            continue;
          } else {
            signAB = 1.0;
          }
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              double signABIJ = signAB;
              i = x[2];
              j = x[3];
              if (i > j) {
                std::swap(i,j);
                signABIJ *= -1.0;
              } else if (i == j) {
                tile[x] = 0.0;
                continue;
              }
              tile[x] = signABIJ * raw2[a + b * (b - 1) / 2 + (i + j * (j - 1) / 2) * nV2shift];
            }
        }
    });

    TA::get_default_world().gop.fence();

  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::initialize(size_t nVec) {
    vecs_.reserve(nVec);
    for (size_t i = 0; i < nVec; i++)
      vecs_.emplace_back(vLabel_, oLabel_);
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                               MatsT alpha, MatsT const *B, int64_t ldb,
                               MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const {

    if (transB != blas::Op::NoTrans)
      CErr("Transpose of B matrix NYI in EOMCCSDVectorSet::multiply_matrix");

    C.scale(beta, shiftC, n);
    
    tryDowncastReferenceTo<EOMCCSDVectorSet<MatsT>>(C,
        [&] (auto& CRef, size_t extraShiftC) {
          shiftC += extraShiftC;
          for (size_t j = 0; j < n; j++) {
            EOMCCSDVector<MatsT> &C_vec = CRef.get(j + shiftC);
            for (size_t i = 0; i < k; i++) {
              C_vec.axpy(alpha * B[i + j * ldb], get(i + shiftA));
            }
          }
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                                            int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA) const {

    tryDowncastReferenceTo<EOMCCSDVectorSet<MatsT>>(B,
        [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < m; i++) {
            const EOMCCSDVector<MatsT> &A_vec = get(i + shiftA);
            for (size_t j = 0; j < n; j++) {
              const EOMCCSDVector<MatsT> &B_vec = BRef.get(j + shiftB);
              C[i + j * ldc] = A_vec.dot(B_vec, conjA);
            }
          }
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable) {
    if (moveable) {
      swap_data(shiftA, nVec, const_cast<SolverVectors<MatsT>&>(B), shiftB);
      return;
    }

    tryDowncastReferenceTo<EOMCCSDVectorSet<MatsT>>(B,
        [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftA) = BRef.get(i + shiftB);
          }
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) {

    tryDowncastReferenceTo<EOMCCSDVectorSet<MatsT>>(B,
        [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftA).swap(BRef.get(i + shiftB));
          }
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::scale(MatsT scalar, size_t shiftA, size_t nVec) {
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "EOMCCSDVectorSet<MatsT>::scale");

    for (size_t i = 0; i < nVec; i++) {
      get(i + shiftA).scale(scalar);
    }
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::conjugate(size_t shiftA, size_t nVec) {
    if (std::is_same<MatsT, double>::value) return;
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "EOMCCSDVectorSet<MatsT>::conjugate");

    for (size_t i = 0; i < nVec; i++) {
      get(i + shiftA).conjugate();
    }
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) {

    tryDowncastReferenceTo<EOMCCSDVectorSet<MatsT>>(X,
        [&] (auto& XRef, size_t extraShiftX) {
          shiftX += extraShiftX;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftY).axpy(alpha, XRef.get(i + shiftX));
          }
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) {
    CErr("EOMCCSDVectorSet::trsm NYI.");
    abort();
  }

  template <typename MatsT>
  int EOMCCSDVectorSet<MatsT>::QR(size_t shift, size_t nVec, MatsT *R, int LDR) {
    CErr("EOMCCSDVectorSet::QR NYI.");
    abort();
  }

  template <typename MatsT>
  double EOMCCSDVectorSet<MatsT>::norm2F(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;

    this->sizeCheck(shift + nVec, "EOMCCSDVectorSet<MatsT>::norm2F");

    double norm = 0.0;

    for (size_t i = 0; i < nVec; i++) {
      double norm_i = get(i + shift).norm();
      norm += norm_i * norm_i;
    }

    return std::sqrt(norm);
  }

  template <typename MatsT>
  double EOMCCSDVectorSet<MatsT>::maxNormElement(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;

    this->sizeCheck(shift + nVec, "EOMCCSDVectorSet<MatsT>::maxNormElement");

    double absmax = 0.0;

    for (size_t i = 0; i < nVec; i++) {
      double absmax_i = get(i + shift).absmax();
      absmax = std::max(absmax, absmax_i);
    }

    return absmax;
  }

  template <typename MatsT>
  RawVectors<MatsT> EOMCCSDVectorSet<MatsT>::toRaw(
      MPI_Comm c, bool includeZeroBody, size_t shift, size_t nVec) const {
    if (nVec == 0) return RawVectors<MatsT>(c, length(includeZeroBody), nVec);

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shift, "EOMCCSDVectorSet<MatsT>::toRaw");
      nVec = size() - shift;
    } else
      this->sizeCheck(shift + nVec, "EOMCCSDVectorSet<MatsT>::toRaw");

//    print(std::cout, "EOMCCSDVectorSet::toRaw::EOM", shift, nVec);

    RawVectors<MatsT> raw(c, length(includeZeroBody), nVec);
    MatsT* rawPtr = raw.getPtr();
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(includeZeroBody));

    for (size_t i = 0; i < nVec; i++)
      get(shift + i).toRaw(rawPtr + i * length(includeZeroBody), includeZeroBody);

//    raw.print(std::cout, "EOMCCSDVectorSet::toRaw::Raw", 0, nVec);

    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

    return raw;

  }

  template <typename MatsT>
  template <typename IntsT>
  RawVectors<MatsT> EOMCCSDVectorSet<MatsT>::toRaw(
      MPI_Comm c, const EOMCCSD<MatsT,IntsT> &eom,
      bool includeZeroBody, size_t shift, size_t nVec) const {
    if (nVec == 0) return RawVectors<MatsT>(c, length(eom, includeZeroBody), nVec);

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shift, "EOMCCSDVectorSet<MatsT>::toRaw");
      nVec = size() - shift;
    } else
      this->sizeCheck(shift + nVec, "EOMCCSDVectorSet<MatsT>::toRaw");

//    print(std::cout, "EOMCCSDVectorSet::toRaw::EOM", shift, nVec);

    RawVectors<MatsT> raw(c, length(eom, includeZeroBody), nVec);
    MatsT* rawPtr = raw.getPtr();
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(eom, includeZeroBody));

    for (size_t i = 0; i < nVec; i++)
      get(shift + i).toRaw(rawPtr + i * length(eom, includeZeroBody), eom, includeZeroBody);

//    raw.print(std::cout, "EOMCCSDVectorSet::toRaw::Raw", 0, nVec);

    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

    return raw;

  }

  template <typename MatsT>
  template <typename IntsT>
  void EOMCCSDVectorSet<MatsT>::fromRaw(MPI_Comm c,
                                        const RawVectors<MatsT> &raw, const EOMCCSD<MatsT,IntsT> &eom,
                                        bool hasZeroBody, size_t shiftThis, size_t shiftRaw, size_t nVec) {
    if (nVec == 0) return;

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shiftThis, "EOMCCSDVectorSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw, "EOMCCSDVectorSet<MatsT>::fromRaw");
      nVec = std::min(size() - shiftThis, raw.size() - shiftRaw);
    } else {
      this->sizeCheck(shiftThis + nVec, "EOMCCSDVectorSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw + nVec, "EOMCCSDVectorSet<MatsT>::fromRaw");
    }

//    raw.print(std::cout, "EOMCCSDVectorSet::fromRaw::Raw", shiftRaw, nVec);

    MatsT* rawPtr = const_cast<MatsT*>(raw.getPtr(shiftRaw));
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(eom, hasZeroBody));

    TA::get_default_world().gop.template broadcast(rawPtr, nVec * length(eom, hasZeroBody), 0);

    for (size_t i = 0; i < nVec; i++)
      get(shiftThis + i).fromRaw(rawPtr + i * length(eom, hasZeroBody), eom, hasZeroBody);

//    print(std::cout, "EOMCCSDVectorSet::fromRaw::EOM", shiftThis, nVec);

    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

  }

  template <typename MatsT>
  void EOMCCSDVectorSet<MatsT>::fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw,
                                        bool hasZeroBody, size_t shiftThis, size_t shiftRaw, size_t nVec) {
    if (nVec == 0) return;

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shiftThis, "EOMCCSDVectorSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw, "EOMCCSDVectorSet<MatsT>::fromRaw");
      nVec = std::min(size() - shiftThis, raw.size() - shiftRaw);
    } else {
      this->sizeCheck(shiftThis + nVec, "EOMCCSDVectorSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw + nVec, "EOMCCSDVectorSet<MatsT>::fromRaw");
    }

    //    raw.print(std::cout, "EOMCCSDVectorSet::fromRaw::Raw", shiftRaw, nVec);

    MatsT* rawPtr = const_cast<MatsT*>(raw.getPtr(shiftRaw));
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(hasZeroBody));

    TA::get_default_world().gop.template broadcast(rawPtr, nVec * length(hasZeroBody), 0);

    for (size_t i = 0; i < nVec; i++)
      get(shiftThis + i).fromRaw(rawPtr + i * length(hasZeroBody), hasZeroBody);

    //    print(std::cout, "EOMCCSDVectorSet::fromRaw::EOM", shiftThis, nVec);

    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

  }

  template <typename MatsT>
  double EOMCCSDVectorSetDebug<MatsT>::compareDebug(size_t shift, size_t nVec) {
    if (nVec == 0) return 0.0;

    if (nVec == std::numeric_limits<size_t>::max()) {
      eomccSet_.sizeCheck(shift, "EOMCCSDVectorSetDebug<MatsT>::compareDebug");
      rawSet_.sizeCheck(shift, "EOMCCSDVectorSetDebug<MatsT>::compareDebug");
      nVec = size() - shift;
    } else {
      eomccSet_.sizeCheck(shift + nVec, "EOMCCSDVectorSetDebug<MatsT>::compareDebug");
      rawSet_.sizeCheck(shift + nVec, "EOMCCSDVectorSetDebug<MatsT>::compareDebug");
    }

    RawVectors<MatsT> raw = eomccSet_.toRaw(rawSet_.getMPIcomm(), false, shift, nVec);

    raw.axpy(0, nVec, -1.0, rawSet_, shift);

    return raw.norm2F(0, nVec);
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                                MatsT alpha, MatsT const *B, int64_t ldb,
                                                MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const {

    tryDowncastReferenceTo<EOMCCSDVectorSetDebug<MatsT>>(C,
        [&] (auto& C_debug, size_t extraShiftC) {
          shiftC += extraShiftC;
          eomccSet_.multiply_matrix(shiftA, transB, n, k, alpha, B, ldb, beta, C_debug.getEOMCCSet(), shiftC);

          TA::get_default_world().gop.fence();
          rawSet_.multiply_matrix(shiftA, transB, n, k, alpha, B, ldb, beta, C_debug.getRawSet(), shiftC);

          std::cout << "EOMCCSDVectorSetDebug::multiply_matrix error = "
                    << C_debug.compareDebug(shiftC, n) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                                                 int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA) const {
    if (m * n == 0) return;

    tryDowncastReferenceTo<EOMCCSDVectorSetDebug<MatsT>>(B,
        [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          eomccSet_.dot_product(shiftA, B_debug.eomccSet_, shiftB, m, n, C, ldc, conjA);

          MatsT *C_ref = CQMemManager::get().malloc<MatsT>(m * n);

          TA::get_default_world().gop.fence();
          rawSet_.dot_product(shiftA, B_debug.rawSet_, shiftB, m, n, C_ref, m, conjA);

          if (ldc == m)
            blas::axpy(m * n, -1.0, C, 1, C_ref, 1);
          else
            for (size_t i = 0; i < n; i++)
              blas::axpy(m, -1.0, C + i * ldc, 1, C_ref + i * m, 1);

          std::cout << "EOMCCSDVectorSetDebug::dot_product error = "
                    << blas::nrm2(m * n, C_ref, 1) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable) {
    if (moveable) {
      swap_data(shiftA, nVec, const_cast<SolverVectors<MatsT>&>(B), shiftB);
      return;
    }

    tryDowncastReferenceTo<EOMCCSDVectorSetDebug<MatsT>>(B,
        [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          eomccSet_.set_data(shiftA, nVec, B_debug.eomccSet_, shiftB);

          TA::get_default_world().gop.fence();
          rawSet_.set_data(shiftA, nVec, B_debug.rawSet_, shiftB);

          std::cout << "EOMCCSDVectorSetDebug::set_data error = "
                    << compareDebug(shiftA, nVec) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) {

    tryDowncastReferenceTo<EOMCCSDVectorSetDebug<MatsT>>(B,
        [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          
          eomccSet_.swap_data(shiftA, nVec, B_debug.eomccSet_, shiftB);

          TA::get_default_world().gop.fence();
          rawSet_.swap_data(shiftA, nVec, B_debug.rawSet_, shiftB);

          std::cout << "EOMCCSDVectorSetDebug::swap error = "
          << compareDebug(shiftA, nVec) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::scale(MatsT scalar, size_t shiftA, size_t nVec) {

    eomccSet_.scale(scalar, shiftA, nVec);

    TA::get_default_world().gop.fence();
    rawSet_.scale(scalar, shiftA, nVec);

    std::cout << "EOMCCSDVectorSetDebug::scale error = "
              << compareDebug(shiftA, nVec) << std::endl;
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::conjugate(size_t shiftA, size_t nVec) {

    eomccSet_.conjugate(shiftA, nVec);

    TA::get_default_world().gop.fence();
    rawSet_.conjugate(shiftA, nVec);

    std::cout << "EOMCCSDVectorSetDebug::conjugate error = "
              << compareDebug(shiftA, nVec) << std::endl;
  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) {

    tryDowncastReferenceTo<EOMCCSDVectorSetDebug<MatsT>>(X,
        [&] (auto& X_debug, size_t extraShiftX) {
          shiftX += extraShiftX;
          eomccSet_.axpy(shiftY, nVec, alpha, X_debug.eomccSet_, shiftX);

          TA::get_default_world().gop.fence();
          rawSet_.axpy(shiftY, nVec, alpha, X_debug.rawSet_, shiftX);

          std::cout << "EOMCCSDVectorSetDebug::axpy error = "
                    << compareDebug(shiftY, nVec) << std::endl;
        }
    );
  }

  template <typename MatsT>
  size_t EOMCCSDVectorSetDebug<MatsT>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                                                   size_t NRe, double eps) {
    double vecDiffBefore = compareDebug(shift, Mold + Mnew);

    std::cout << "EOMCCSDVectorSetDebug::GramSchmidt before error = "
              << vecDiffBefore << std::endl;

    size_t iOrtho = SolverVectors<MatsT>::GramSchmidt(shift, Mold, Mnew, NRe, eps);

    double vecDiffAfter = compareDebug(shift, Mold + Mnew);

    std::cout << "EOMCCSDVectorSetDebug::GramSchmidt after error = "
              << vecDiffAfter << std::endl;

    return iOrtho;

  }

  template <typename MatsT>
  void EOMCCSDVectorSetDebug<MatsT>::trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) {
    eomccSet_.trsm(shift, n, alpha, A, lda);

    TA::get_default_world().gop.fence();
    rawSet_.trsm(shift, n, alpha, A, lda);

    std::cout << "EOMCCSDVectorSetDebug::trsm error = "
              << compareDebug(shift, n) << std::endl;

  }

  template <typename MatsT>
  int EOMCCSDVectorSetDebug<MatsT>::QR(size_t shift, size_t nVec, MatsT *R, int LDR) {
    int iOrtho = eomccSet_.QR(shift, nVec, R, LDR);

    TA::get_default_world().gop.fence();
    int iOrthoRaw = rawSet_.QR(shift, nVec, R, LDR);

    std::cout << "EOMCCSDVectorSetDebug::QR error = "
              << compareDebug(shift, nVec) << std::endl;

    if (iOrtho != iOrthoRaw)
      std::cout << "EOMCCSDVectorSetDebug::QR iOrtho differs: eomccSet = "
                << iOrtho << ", rawSet = " << iOrthoRaw << std::endl;

    return iOrtho;
  }

  template <typename MatsT>
  double EOMCCSDVectorSetDebug<MatsT>::norm2F(size_t shift, size_t nVec) const {

    double norm = eomccSet_.norm2F(shift, nVec);

    TA::get_default_world().gop.fence();
    double normRaw = rawSet_.norm2F(shift, nVec);

    std::cout << "EOMCCSDVectorSetDebug::norm2F error = "
              << std::abs(norm - normRaw) << std::endl;

    return norm;
  }

  template <typename MatsT>
  double EOMCCSDVectorSetDebug<MatsT>::maxNormElement(size_t shift, size_t nVec) const {

    double absmax = eomccSet_.maxNormElement(shift, nVec);

    TA::get_default_world().gop.fence();
    double absmaxRaw = rawSet_.maxNormElement(shift, nVec);

    std::cout << "EOMCCSDVectorSetDebug::maxNormElement error = "
              << std::abs(absmax - absmaxRaw) << std::endl;

    return absmax;
  }

}
