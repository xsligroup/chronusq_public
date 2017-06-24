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

#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>
#include <util/math.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <functional>
#include <util/timer.hpp>

namespace ChronusQ{

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::assignCVSIndices() {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    CVSOindices_.clear();
    CVSOindices_.resize(nO, nO);

    size_t oIndex = nO;
    for (auto it = eomSettings.frozen_occupied.rbegin(); it != eomSettings.frozen_occupied.rend(); it++) {
      if (*it >= nO + nV)
        CErr("CVS-EOMCC: Orbital index in input EOMCC.FROZENOCCUPIED out of range");
      if (*it >= nO)
        CErr("CVS-EOMCC: Virtual orbital appears in input EOMCC.FROZENOCCUPIED");
      CVSOindices_[*it] = --oIndex;
    }

    oIndex = 0;
    for (size_t i : eomSettings.cvs_core) {
      if (i >= nO + nV)
        CErr("CVS-EOMCC: Orbital index in input EOMCC.CVSCORE out of range");
      if (i >= nO)
        CErr("CVS-EOMCC: Virtual orbital appears in input EOMCC.CVSCORE");
      if (CVSOindices_[i] != nO)
        CErr("CVS-EOMCC: Orbital appear in both EOMCC.CVSCORE and EOMCC.FROZENOCCUPIED inputs");
      CVSOindices_[i] = oIndex++;
    }

    for (size_t &i : CVSOindices_)
      if (i == nO)
        i = oIndex++;


    CVSVindices_.clear();
    CVSVindices_.resize(nV, nV);

    size_t vIndex = nV;
    for (auto it = eomSettings.frozen_virtual.rbegin(); it != eomSettings.frozen_virtual.rend(); it++) {
      if (*it >= nO + nV)
        CErr("CVS-EOMCC: Orbital index in input EOMCC.FROZENVRITUAL out of range");
      if (*it < nO)
        CErr("CVS-EOMCC: Occupied orbital appears in input EOMCC.FROZENVRITUAL");
      CVSVindices_[*it - nO] = --vIndex;
    }

    vIndex = 0;

    for (size_t i : eomSettings.cvs_virtual) {
      if (i >= nO + nV)
        CErr("CVS-EOMCC: Orbital index in input EOMCC.CVSVIRTUAL out of range");
      if (i < nO)
        CErr("CVS-EOMCC: Occupied orbital appears in input EOMCC.CVSVIRTUAL");
      if (CVSVindices_[i - nO] != nV)
        CErr("CVS-EOMCC: Orbital appear in both EOMCC.CVSVIRTUAL and EOMCC.FROZENVRITUAL inputs");
      CVSVindices_[i - nO] = vIndex++;
    }

    for (size_t &i : CVSVindices_)
      if (i == nV)
        i = vIndex++;

    nCVSOActive_ = nO - eomSettings.frozen_occupied.size();
    nCVSOCore_ = eomSettings.cvs_core.size() == 0? nCVSOActive_ : eomSettings.cvs_core.size();
    nCVSOValance_ = nCVSOActive_ - nCVSOCore_;

    nCVSVActive_ = nV - eomSettings.frozen_virtual.size();
    nCVSVContinuum_ = eomSettings.cvs_virtual.size() == 0? nCVSVActive_ : eomSettings.cvs_virtual.size();
    nCVSVValance_ = nCVSVActive_ - nCVSVContinuum_;


    nOVshift_ = nCVSOCore_ * nCVSVContinuum_;
    nO2shift_ = nCVSOCore_ * (nCVSOCore_ - 1) / 2 + nCVSOCore_ * nCVSOValance_;
    nV2shift_ = nCVSVContinuum_ * (nCVSVContinuum_ - 1) / 2 + nCVSVContinuum_ * nCVSVValance_;

    Hbar_dim = nOVshift_ + nO2shift_ * nV2shift_;
    CVSoutOfBound_ = (Hbar_dim + 1) * (Hbar_dim + 1);

    CVSabIndices_.clear();
    CVSabIndices_.resize(nV, std::vector<size_t>(nV, CVSoutOfBound_));

    std::vector<size_t> CVSVorbitals_(nV, 0);
    for (size_t a = 0; a < nV; a++) {
      CVSVorbitals_[CVSVindices_[a]] = a;
    }

    size_t idx = 0;
    for (size_t b = 0; b < nCVSVActive_; b++) {
      for (size_t a = 0; a < std::min(b, nCVSVContinuum_); a++) {
        size_t aa = CVSVorbitals_[a], bb = CVSVorbitals_[b];
        CVSabIndices_[aa][bb] = idx;
        CVSabIndices_[bb][aa] = idx++;
      }
    }


    CVSijIndices_.clear();
    CVSijIndices_.resize(nO, std::vector<size_t>(nO, CVSoutOfBound_));

    std::vector<size_t> CVSOorbitals_(nO, 0);
    for (size_t i = 0; i < nO; i++) {
      CVSOorbitals_[CVSOindices_[i]] = i;
    }

    idx = 0;
    for (size_t j = 0; j < nCVSOActive_; j++) {
      for (size_t i = 0; i < std::min(j, nCVSOCore_); i++) {
        size_t ii = CVSOorbitals_[i], jj = CVSOorbitals_[j];
        CVSijIndices_[ii][jj] = idx;
        CVSijIndices_[jj][ii] = idx++;
      }
    }

  }

  template <typename MatsT, typename IntsT>
  inline size_t EOMCCSD<MatsT,IntsT>::CVStoCompoundS(size_t a, size_t i) const {
    size_t aa = CVSVindices_[a], ii = CVSOindices_[i];
    if (aa >= nCVSVContinuum_ or ii >= nCVSOCore_)
      return CVSoutOfBound_;
    return aa + ii * nCVSVContinuum_;
  }

  template <typename MatsT, typename IntsT>
  inline size_t EOMCCSD<MatsT,IntsT>::CVStoCompoundD(size_t a, size_t b, size_t i, size_t j) const {
    size_t ab = CVSabIndices_[a][b], ij = CVSijIndices_[i][j];
    if (ab == CVSoutOfBound_ or ij == CVSoutOfBound_)
      return CVSoutOfBound_;
    return ab + ij * nV2shift_;
  }

  template <typename MatsT, typename IntsT>
  inline size_t EOMCCSD<MatsT,IntsT>::CVStoCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const {
    size_t ai = CVStoCompoundS(a,i), bj = CVStoCompoundS(b,j);
    if (ai == CVSoutOfBound_ or bj == CVSoutOfBound_)
      return CVSoutOfBound_;
    return ai + bj * ldH;
  }

  template <typename MatsT, typename IntsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT,IntsT>::CVStoCompoundSD(size_t e, size_t m,
                                                                         size_t a, size_t b, size_t i, size_t j, size_t ldH) const {
    size_t em = CVStoCompoundS(e,m), abij = CVStoCompoundD(a,b,i,j);
    if (em == CVSoutOfBound_ or abij == CVSoutOfBound_)
      return std::make_pair(CVSoutOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    return std::make_pair(em + abij * ldH, sign);
  }

  template <typename MatsT, typename IntsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT,IntsT>::CVStoCompoundDS(size_t a, size_t b, size_t i, size_t j,
                                                                         size_t e, size_t m, size_t ldH) const {
    size_t em = CVStoCompoundS(e,m), abij = CVStoCompoundD(a,b,i,j);
    if (em == CVSoutOfBound_ or abij == CVSoutOfBound_)
      return std::make_pair(CVSoutOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    return std::make_pair(abij + em * ldH, sign);
  }

  template <typename MatsT, typename IntsT>
  inline std::pair<size_t, double> EOMCCSD<MatsT,IntsT>::CVStoCompoundDD(size_t a, size_t b, size_t i, size_t j,
                                                                         size_t c, size_t d, size_t k, size_t l, size_t ldH) const {
    size_t abij = CVStoCompoundD(a,b,i,j), cdkl = CVStoCompoundD(c,d,k,l);
    if (abij == CVSoutOfBound_ or cdkl == CVSoutOfBound_)
      return std::make_pair(CVSoutOfBound_, 0.0);
    double sign = signD(a,b,i,j);
    sign *= signD(c,d,k,l);
    return std::make_pair(abij + cdkl * ldH, sign);
  }


  template <typename MatsT, typename IntsT>
  cqmatrix::Matrix<MatsT> EOMCCSD<MatsT,IntsT>::buildHbarCVS(bool includeGroundState) const {
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    cqmatrix::Matrix<MatsT> fullMat(includeGroundState ? Hbar_dim + 1 : Hbar_dim);
    fullMat.clear();

    MatsT * Hbar = fullMat.pointer();
    size_t ldH = fullMat.dimension();
    size_t nCol = fullMat.dimension();

    MatsT * Hbar0S = nullptr;
    MatsT * Hbar0D = nullptr;
    MatsT * HbarSS = Hbar;
    MatsT * HbarSD = HbarSS + nOVshift_ * ldH;
    MatsT * HbarDS = HbarSS + nOVshift_;
    MatsT * HbarDD = HbarSD + nOVshift_;
    if (includeGroundState) {
      Hbar0S = Hbar + ldH;
      Hbar0D = Hbar0S + nOVshift_ * ldH;
      HbarSS += 1 + ldH;
      HbarSD += 1 + ldH;
      HbarDS += 1 + ldH;
      HbarDD += 1 + ldH;
    }

    TA::foreach_inplace( F_ae, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") = F_ae("a,c") * R1("c,i");
          size_t a = x[0], c = x[1];

          for (size_t i = 0; i < nO; ++i) {
            size_t idx = CVStoCompoundSS(a,i,c,i,ldH);
            if (CVSisInBound(idx))
              HbarSS[idx] += v;
          }

          // tildeR2("a,b,i,j") = F_ae("b,e") * R2("a,e,i,j");
          // tildeR2("a,b,i,j") += - F_ae("a,e") * R2("b,e,i,j");
          size_t b = x[0], e = x[1];

          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i)
              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = CVStoCompoundDD(a,b,i,j,a,e,i,j,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarDD[idx_sgn.first] += idx_sgn.second * v;
              }

        }
    });

    TA::foreach_inplace( F_mi, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") += - F_mi("k,i") * R1("a,k");
          size_t k = x[0], i = x[1];

          for (size_t a = 0; a < nV; ++a) {
            size_t idx = CVStoCompoundSS(a,i,a,k,ldH);
            if (CVSisInBound(idx))
              HbarSS[idx] -= v;
          }

          // tildeR2("a,b,i,j") += - F_mi("k,j") * R2("a,b,i,k");
          // tildeR2("a,b,i,j") += F_mi("k,i") * R2("a,b,j,k");
          size_t j = x[1];

          for (size_t b = 0; b < nV; ++b)
            for (size_t a = 0; a < b; ++a)
              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = CVStoCompoundDD(a,b,i,j,a,b,i,k,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarDD[idx_sgn.first] -= idx_sgn.second * v;
              }
        }
    });

    TA::foreach_inplace( F_me, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          MatsT v = tile[x];

          // tildeR1("a,i") += F_me("k,c") * R2("a,c,i,k");
          size_t k = x[0], c = x[1];

          if (includeGroundState) {
            size_t idx = CVStoCompoundS(c,k);
            if (CVSisInBound(idx))
              Hbar0S[idx * ldH] = v;
          }

          for (size_t a = 0; a < nV; ++a)
            for (size_t i = 0; i < nO; ++i) {
              auto idx_sgn = CVStoCompoundSD(a,i,a,c,i,k,ldH);
              if (CVSisInBound(idx_sgn.first))
                HbarSD[idx_sgn.first] += idx_sgn.second * v;
            }
        }
    });

    if (includeGroundState)
      TA::foreach_inplace( antiSymMoints["vvoo"], [&](TA::Tensor<MatsT>& tile) {

        const auto& lobound = tile.range().lobound();
        if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
          return;

        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0,0,0,0};
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
                MatsT v = tile[x];

                size_t a = x[0], b = x[1], i = x[2], j = x[3];

                size_t idx = CVStoCompoundD(a,b,i,j);
                if (CVSisInBound(idx))
                  Hbar0D[idx * ldH] = SmartConj(v);

              }

      });

    TA::foreach_inplace( W_mbej, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += W_mbej("k,a,c,i") * R1("c,k");
              size_t k = x[0], a = x[1], c = x[2], i = x[3];
              size_t idx = CVStoCompoundSS(a,i,c,k,ldH);
              if (CVSisInBound(idx))
                HbarSS[idx] += v;

              // tildeR2("a,b,i,j") += W_mbej("k,b,c,j") * R2("a,c,i,k");
              // tildeR2("a,b,i,j") += - W_mbej("k,a,c,j") * R2("b,c,i,k");
              // tildeR2("a,b,i,j") += - W_mbej("k,b,c,i") * R2("a,c,j,k");
              // tildeR2("a,b,i,j") += W_mbej("k,a,c,i") * R2("b,c,j,k");
              for (size_t b = 0; b < nV; ++b)
                for (size_t j = 0; j < nO; ++j) {
                  auto idx_sgn = CVStoCompoundDD(a,b,i,j,b,c,j,k,ldH);
                  if (CVSisInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_mnij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += 0.5 * W_mnij("k,l,i,j") * R2("a,b,k,l");
              size_t k = x[0], l = x[1], i = x[2], j = x[3];

              for (size_t b = 0; b < nV; ++b)
                for (size_t a = 0; a < b; ++a) {
                  auto idx_sgn = CVStoCompoundDD(a,b,i,j,a,b,k,l,ldH);
                  if (CVSisInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_abef, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += 0.5 * W_abef("a,b,e,f") * R2("e,f,i,j");
              size_t a = x[0], b = x[1], e = x[2], f = x[3];

              for (size_t j = 0; j < nO; ++j)
                for (size_t i = 0; i < j; ++i) {
                  auto idx_sgn = CVStoCompoundDD(a,b,i,j,e,f,i,j,ldH);
                  if (CVSisInBound(idx_sgn.first))
                    HbarDD[idx_sgn.first] += idx_sgn.second * v;
                }

            }
    });

    TA::foreach_inplace( W_abei, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += W_abei("a,b,c,j") * R1("c,i");
              // tildeR2("a,b,i,j") += - W_abei("a,b,c,i") * R1("c,j");
              size_t a = x[0], b = x[1], c = x[2], j = x[3];

              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = CVStoCompoundDS(a,b,i,j,c,i,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarDS[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_mbij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[2] < lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < std::min(x[2], static_cast<std::size_t>(upbound[3])); ++x[3]) {
              MatsT v = tile[x];

              // tildeR2("a,b,i,j") += -W_mbij("k,a,j,i") * R1("b,k");
              // tildeR2("a,b,i,j") += W_mbij("k,b,j,i") * R1("a,k");
              size_t k = x[0], b = x[1], j = x[2], i = x[3];

              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = CVStoCompoundDS(a,b,i,j,a,k,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarDS[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_amef, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[2] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
            for(x[2] = lobound[2]; x[2] < std::min(x[3], static_cast<std::size_t>(upbound[2])); ++x[2]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += 0.5 * W_amef("a,k,c,d") * R2("c,d,i,k");
              size_t a = x[0], k = x[1], c = x[2], d = x[3];

              for (size_t i = 0; i < nO; ++i) {
                auto idx_sgn = CVStoCompoundSD(a,i,c,d,i,k,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarSD[idx_sgn.first] += idx_sgn.second * v;
              }

            }
    });

    TA::foreach_inplace( W_mnie, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              MatsT v = tile[x];

              // tildeR1("a,i") += - 0.5 * W_mnie("k,l,i,d") * R2("a,d,k,l");
              size_t k = x[0], l = x[1], i = x[2], d = x[3];

              for (size_t a = 0; a < nV; ++a) {
                auto idx_sgn = CVStoCompoundSD(a,i,a,d,k,l,ldH);
                if (CVSisInBound(idx_sgn.first))
                  HbarSD[idx_sgn.first] -= idx_sgn.second * v;
              }

            }
    });

    TArray WT_ckabij = TAmanager.malloc<MatsT>("vovvoo");
    WT_ckabij("c,k,a,b,i,j") = 0.5 * W_amef("b,k,d,c") * T2_("a,d,i,j");
    WT_ckabij("c,k,a,b,i,j") -= 0.5 * W_mnie("l,k,j,c") * T2_("a,b,i,l");
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( WT_ckabij, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4])
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") += W_amef("b,k,d,c") * T2_("a,d,i,j") * R1("c,k");
                  // tildeR2("a,b,i,j") -= W_amef("a,k,d,c") * T2_("b,d,i,j") * R1("c,k");
                  // tildeR2("a,b,i,j") -= W_mnie("l,k,j,c") * T2_("a,b,i,l") * R1("c,k");
                  // tildeR2("a,b,i,j") += W_mnie("l,k,i,c") * T2_("a,b,j,l") * R1("c,k");
                  size_t c = x[0], k = x[1], a = x[2], b = x[3], i = x[4], j = x[5];

                  auto idx_sgn = CVStoCompoundDS(a,b,i,j,c,k,ldH);
                  if (CVSisInBound(idx_sgn.first))
                    HbarDS[idx_sgn.first] += idx_sgn.second * v;

                }
    });
    TAmanager.free("vovvoo", std::move(WT_ckabij), true);

    TArray TV_abidck = TAmanager.malloc<MatsT>("vvovvo");
    TV_abidck("a,b,i,d,c,k") = T2_("a,b,i,l") * conj(antiSymMoints.at("vvoo")("d,c,k,l"));
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( TV_abidck, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] > lobound[1] or lobound[4] > lobound[3])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < std::min(x[3], static_cast<std::size_t>(upbound[4])); ++x[4])
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") -= 0.5 * T2_("a,b,i,l") * conj(V("d,c,k,l")) * R2("c,d,j,k");
                  // tildeR2("a,b,i,j") += 0.5 * T2_("a,b,j,l") * conj(V("d,c,k,l")) * R2("c,d,i,k");
                  size_t a = x[0], b = x[1], i = x[2], d = x[3], c = x[4], k = x[5];

                  for (size_t j = 0; j < nO; ++j) {
                    auto idx_sgn = CVStoCompoundDD(a,b,i,j,c,d,j,k,ldH);
                    if (CVSisInBound(idx_sgn.first))
                      HbarDD[idx_sgn.first] -= idx_sgn.second * v;
                  }

                }
    });
    TAmanager.free("vvovvo", std::move(TV_abidck), true);

    TArray TV_aijdkl = TAmanager.malloc<MatsT>("voovoo");
    TV_aijdkl("a,i,j,d,k,l") = T2_("e,a,i,j") * conj(antiSymMoints.at("vvoo")("e,d,k,l"));
    TA::get_default_world().gop.fence();

    TA::foreach_inplace( TV_aijdkl, [&](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[1] > lobound[2] or lobound[4] > lobound[5])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
          for(x[1] = lobound[1]; x[1] < std::min(x[2], static_cast<std::size_t>(upbound[1])); ++x[1])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5])
                for(x[4] = lobound[4]; x[4] < std::min(x[5], static_cast<std::size_t>(upbound[4])); ++x[4]) {
                  MatsT v = tile[x];

                  // tildeR2("a,b,i,j") += 0.5 * T2_("e,a,i,j") * conj(V("e,d,k,l")) * R2("b,d,k,l");
                  // tildeR2("a,b,i,j") -= 0.5 * T2_("e,b,i,j") * conj(V("e,d,k,l")) * R2("a,d,k,l");
                  size_t a = x[0], i = x[1], j = x[2], d = x[3], k = x[4], l = x[5];

                  for (size_t b = 0; b < nV; ++b) {
                    auto idx_sgn = CVStoCompoundDD(a,b,i,j,b,d,k,l,ldH);
                    if (CVSisInBound(idx_sgn.first))
                      HbarDD[idx_sgn.first] += idx_sgn.second * v;
                  }

                }
    });
    TAmanager.free("voovoo", std::move(TV_aijdkl), true);

    TA::get_default_world().gop.fence();

    TA::get_default_world().gop.template reduce(Hbar, ldH*nCol, std::plus<MatsT>());
    TA::get_default_world().gop.fence();

    return fullMat;

  }

};