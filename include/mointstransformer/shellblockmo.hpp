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

#include <util/timer.hpp>
#include <cqlinalg.hpp>
#include <matrix.hpp>

namespace ChronusQ {

/*
 * ShellBlockMO interface
 */ 
template <typename MatsT>
class ShellBlockMO {
 public:
  ShellBlockMO() = default;
  ShellBlockMO(const ShellBlockMO &) = delete;
  ShellBlockMO(ShellBlockMO &&) = delete;
  
  // compute density(nu, mu) = C(nu, q) x C(mu, p)^*
  void computeAODensityFromMO(
      size_t p, size_t q, 
      const cqmatrix::Matrix<MatsT>& shell_mu_block_mo,
      const cqmatrix::Matrix<MatsT>& shell_nu_block_mo,
      cqmatrix::Matrix<MatsT>& density, bool increment = false) const {
    
    size_t nMu = shell_mu_block_mo.nRows();
    size_t nNu = shell_nu_block_mo.nRows();
    density.resize(nNu, nMu);
    MatsT scale = increment ? MatsT(1.0) : MatsT(0.0);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
        nNu, nMu, 1, MatsT(1.), shell_nu_block_mo.pointer() + q * nNu, nNu,
        shell_mu_block_mo.pointer() + p * nMu, nMu, scale, density.pointer(), nNu);
  }
  
  // transform out(p, q) = C(mu, p)^* x in(mu, nu) x C(nu, q) 
  void transform(const cqmatrix::Matrix<MatsT>& in, 
      const cqmatrix::Matrix<MatsT>& shell_mu_block_mo,
      const cqmatrix::Matrix<MatsT>& shell_nu_block_mo,
      cqmatrix::Matrix<MatsT>& intermeidate,
      cqmatrix::Matrix<MatsT>& out, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const {
    
    size_t pOff = pOff_size.first;
    size_t qOff = qOff_size.first;
    size_t np = pOff_size.second;
    size_t nq = qOff_size.second;
    size_t nMu = shell_mu_block_mo.nRows();
    size_t nNu = shell_nu_block_mo.nRows();

    // intermeidate(nu, p) = ConjTrans(in(mu, nu)) x C(mu, p)
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
        nNu, np, nMu, MatsT(1.), in.pointer(), nMu, shell_mu_block_mo.pointer() + pOff * nMu, nMu, 
        MatsT(0.), intermeidate.pointer(), nNu);
    
    MatsT scale = increment ? MatsT(1.0) : MatsT(0.0);
    // out(p, q) = ConjTrans(intermeidate(nu, p)) x C(nu, q)
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
        np, nq, nNu, MatsT(1.), intermeidate.pointer(), nNu, shell_nu_block_mo.pointer() + qOff * nNu, nNu, 
        scale, out.pointer(), np);
  }

  // major interfaces
  virtual void genSymmDenLLMS(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& symmDenLLMS) const = 0;  

  virtual void transformLL(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size, 
      bool increment = false) const = 0;
}; // class ShellBlockMO

template <typename MatsT>
class OneCShellBlockMO: public ShellBlockMO<MatsT> {
 protected:
  std::vector<cqmatrix::Matrix<MatsT>> shell_block_mo_;
  mutable std::vector<cqmatrix::Matrix<MatsT>> interSCR_;
  mutable std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSCR_; 

 public:
  OneCShellBlockMO() = delete;
  OneCShellBlockMO(const OneCShellBlockMO &) = delete;
  OneCShellBlockMO(OneCShellBlockMO &&) = delete;
  OneCShellBlockMO(const cqmatrix::Matrix<MatsT>& mo,
      const std::vector<size_t>& shellSizes,
      size_t maxTransfromSize) {
    
    shell_block_mo_.reserve(shellSizes.size());
    const MatsT* mo_ptr = mo.pointer();
    size_t nAO = mo.nRows();
    size_t nMO = mo.nColumns();
    size_t maxShellSize = 0ul;
    for (const auto& shSize : shellSizes) {
      shell_block_mo_.emplace_back(shSize, nMO);
      SetMat('N', shSize, nMO, MatsT(1.), mo_ptr, nAO, shell_block_mo_.back().pointer(), shSize); 
      mo_ptr += shSize;
      maxShellSize = std::max(maxShellSize, shSize);
    }
    for (auto i = 0ul; i < GetNumThreads(); ++i) {
      interSCR_.emplace_back(maxShellSize, maxTransfromSize);
      pauliSCR_.emplace_back(maxShellSize, false, false);
    }
  }
  
  // major interfaces
  void genSymmDenLLMS(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& symmDenLLMS) const override {
    this->computeAODensityFromMO(p, q, shell_block_mo_[shell_mu], shell_block_mo_[shell_nu], symmDenLLMS.S());
    auto& pauliSCR = pauliSCR_[GetThreadID()];
    this->computeAODensityFromMO(p, q, shell_block_mo_[shell_nu], shell_block_mo_[shell_mu], pauliSCR.S());
    MatrixAXPY('T', MatsT(1.), pauliSCR.S(), symmDenLLMS.S()); 
  }
  
  void transformLL(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const override {
    this->transform(pauli.S(), shell_block_mo_[shell_mu], shell_block_mo_[shell_nu],
        interSCR_[GetThreadID()], mat, pOff_size, qOff_size, increment);
  }
  
}; // class OneCShellBlockMO

template <typename MatsT>
class TwoCShellBlockMO: public ShellBlockMO<MatsT> {
 protected:
  std::vector<cqmatrix::Matrix<MatsT>> shell_block_mo_;
  mutable std::vector<cqmatrix::Matrix<MatsT>> interSCR_;
  mutable std::vector<cqmatrix::Matrix<MatsT>> spinorSCR_;
  mutable std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSCR_; 

 public:
  TwoCShellBlockMO() = delete;
  TwoCShellBlockMO(const TwoCShellBlockMO &) = delete;
  TwoCShellBlockMO(TwoCShellBlockMO &&) = delete;
  TwoCShellBlockMO(const cqmatrix::Matrix<MatsT>& mo,
      const std::vector<size_t>& shellSizes, 
      size_t maxTransfromSize) {
    shell_block_mo_.reserve(shellSizes.size());
    const MatsT* mo_ptr = mo.pointer();
    size_t nAO = mo.nRows();
    size_t nAO1C = nAO / 2;
    size_t nMO = mo.nColumns();
    size_t maxShellSize = 0ul;
    for (const auto& shSize1C : shellSizes) {
      size_t shSize = shSize1C * 2;
      shell_block_mo_.emplace_back(shSize, nMO);
      // alpha part
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr, nAO, shell_block_mo_.back().pointer(), shSize); 
      // beta part
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr + nAO1C, nAO, shell_block_mo_.back().pointer() + shSize1C, shSize); 
      mo_ptr += shSize1C;
      maxShellSize = std::max(maxShellSize, shSize1C);
    } 
    for (auto i = 0ul; i < GetNumThreads(); ++i) {
      interSCR_.emplace_back(maxShellSize * 2, maxTransfromSize);
      spinorSCR_.emplace_back(maxShellSize * 2);  
      pauliSCR_.emplace_back(maxShellSize, false, false);
    }
  }
  
  // major interfaces
  void genSymmDenLLMS(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& symmDenLLMS) const override {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& pauliSCR = pauliSCR_[thread_id];
    
    size_t nNu = shell_block_mo_[shell_nu].nRows() / 2; 
    size_t nMu = shell_block_mo_[shell_mu].nRows() / 2;
    symmDenLLMS.resize(nNu, nMu);
    pauliSCR.resize(nMu, nNu);
    
    this->computeAODensityFromMO(p ,q, shell_block_mo_[shell_mu], shell_block_mo_[shell_nu], spinorSCR); 
    spinorSCR.spinScatter(symmDenLLMS, false, false);
    
    this->computeAODensityFromMO(p ,q, shell_block_mo_[shell_nu], shell_block_mo_[shell_mu], spinorSCR); 
    spinorSCR.spinScatter(pauliSCR, false, false);
    
    MatrixAXPY('T', MatsT(1.), pauliSCR.S(), symmDenLLMS.S());
  }
  
  void transformLL(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const override {
    
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& intermeidate = interSCR_[thread_id];
    
    spinorSCR.resize(2 * pauli.nRows(), 2 * pauli.nColumns());
    pauli.spinGather(spinorSCR);
    
    this->transform(spinorSCR, shell_block_mo_[shell_mu], shell_block_mo_[shell_nu],
        intermeidate, mat, pOff_size, qOff_size, increment);
  }

}; // class TwoCShellBlockMO

template <typename MatsT>
class FourCShellBlockMO: public ShellBlockMO<MatsT> {
 protected:
  std::vector<cqmatrix::Matrix<MatsT>> shell_block_large_mo_;
  std::vector<cqmatrix::Matrix<MatsT>> shell_block_small_mo_;
  mutable std::vector<cqmatrix::Matrix<MatsT>> interSCR_;
  mutable std::vector<cqmatrix::Matrix<MatsT>> spinorSCR_;
  mutable std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSCR_; 

 public:
  FourCShellBlockMO() = delete;
  FourCShellBlockMO(const FourCShellBlockMO &) = delete;
  FourCShellBlockMO(FourCShellBlockMO &&) = delete;
  FourCShellBlockMO(const cqmatrix::Matrix<MatsT>& mo,
      const std::vector<size_t>& shellSizes, 
      size_t maxTransfromSize) {
    const MatsT* mo_ptr = mo.pointer();
    size_t nAO = mo.nRows();
    size_t nAO1C = nAO / 4;
    size_t nAO2C = nAO / 2;
    size_t nMO = mo.nColumns(); 
    size_t maxShellSize = 0ul;
    for (const auto& shSize1C : shellSizes) {
      size_t shSize2C = shSize1C * 2;
      shell_block_large_mo_.emplace_back(shSize2C, nMO);
      shell_block_small_mo_.emplace_back(shSize2C, nMO);
      // alpha part
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr, nAO, shell_block_large_mo_.back().pointer(), shSize2C); 
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr + nAO1C, nAO, shell_block_small_mo_.back().pointer(), shSize2C); 
      // beta part
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr + nAO2C, nAO, shell_block_large_mo_.back().pointer() + shSize1C, shSize2C); 
      SetMat('N', shSize1C, nMO, MatsT(1.), mo_ptr + nAO1C + nAO2C, nAO, shell_block_small_mo_.back().pointer() + shSize1C, shSize2C); 
      mo_ptr += shSize1C;
      maxShellSize = std::max(maxShellSize, shSize1C);
    } 
    
    for (auto i = 0ul; i < GetNumThreads(); ++i) {
      interSCR_.emplace_back(maxShellSize * 2, maxTransfromSize);
      spinorSCR_.emplace_back(maxShellSize * 2);  
      pauliSCR_.emplace_back(maxShellSize, true, true);
    }
  }
  
  // major interfaces
  void genSymmDenLLMS(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& symmDenLLMS) const override {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& pauliSCR = pauliSCR_[thread_id];
    
    this->computeAODensityFromMO(p ,q, shell_block_large_mo_[shell_mu], shell_block_large_mo_[shell_nu], spinorSCR); 
    size_t nNu = spinorSCR.nRows() / 2;
    size_t nMu = spinorSCR.nColumns() / 2;
    symmDenLLMS.resize(nNu, nMu);
    spinorSCR.spinScatter(symmDenLLMS, false, false);
    
    this->computeAODensityFromMO(p ,q, shell_block_large_mo_[shell_nu], shell_block_large_mo_[shell_mu], spinorSCR); 
    pauliSCR.resize(nMu, nNu);
    spinorSCR.spinScatter(pauliSCR, false, false);
    
    MatrixAXPY('T', MatsT(1.), pauliSCR.S(), symmDenLLMS.S()); 
  } // genSymmDenLLMS

  void genSymmDenSS(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& symmDenSS) const {
    
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& pauliSCR = pauliSCR_[thread_id];
    
    this->computeAODensityFromMO(p ,q, shell_block_small_mo_[shell_mu], shell_block_small_mo_[shell_nu], spinorSCR); 
    size_t nNu = spinorSCR.nRows() / 2;
    size_t nMu = spinorSCR.nColumns() / 2;
    symmDenSS.resize(nNu, nMu);
    spinorSCR.spinScatter(symmDenSS, true, true);
    
    this->computeAODensityFromMO(p ,q, shell_block_small_mo_[shell_nu], shell_block_small_mo_[shell_mu], spinorSCR); 
    pauliSCR.resize(nMu, nNu);
    spinorSCR.spinScatter(pauliSCR, true, true);
    
    MatrixAXPY('T', MatsT(1.), pauliSCR.S(), symmDenSS.S()); 
    MatrixAXPY('T', MatsT(-1.), pauliSCR.Z(), symmDenSS.Z()); 
    MatrixAXPY('T', MatsT(-1.), pauliSCR.Y(), symmDenSS.Y()); 
    MatrixAXPY('T', MatsT(-1.), pauliSCR.X(), symmDenSS.X()); 
    
  } // genSymmDenSS
  
  void genDenLSpmDenSL(size_t p, size_t q, size_t shell_mu, size_t shell_nu,
      cqmatrix::PauliSpinorMatrices<MatsT>& denLSpmDenSL) const {

    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& pauliSCR = pauliSCR_[thread_id];
    
    this->computeAODensityFromMO(p ,q, shell_block_large_mo_[shell_mu], shell_block_small_mo_[shell_nu], spinorSCR); 
    size_t nNu = spinorSCR.nRows() / 2;
    size_t nMu = spinorSCR.nColumns() / 2;
    denLSpmDenSL.resize(nNu, nMu);
    spinorSCR.spinScatter(denLSpmDenSL, true, true);
    
    this->computeAODensityFromMO(p ,q, shell_block_small_mo_[shell_nu], shell_block_large_mo_[shell_mu], spinorSCR); 
    pauliSCR.resize(nMu, nNu);
    spinorSCR.spinScatter(pauliSCR, true, true);
    
    MatrixAXPY('T', MatsT(-1.), pauliSCR.S(), denLSpmDenSL.S()); 
    MatrixAXPY('T', MatsT(1.), pauliSCR.Z(), denLSpmDenSL.Z()); 
    MatrixAXPY('T', MatsT(1.), pauliSCR.Y(), denLSpmDenSL.Y()); 
    MatrixAXPY('T', MatsT(1.), pauliSCR.X(), denLSpmDenSL.X()); 
  }
  
  void transformLL(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const override {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& intermeidate = interSCR_[thread_id];
    spinorSCR.resize(2 * pauli.nRows(), 2 * pauli.nColumns());
    pauli.spinGather(spinorSCR);
    this->transform(spinorSCR, shell_block_large_mo_[shell_mu], shell_block_large_mo_[shell_nu],
        intermeidate, mat, pOff_size, qOff_size, increment);
  }

  void transformLS(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& intermeidate = interSCR_[thread_id];
    spinorSCR.resize(2 * pauli.nRows(), 2 * pauli.nColumns());
    pauli.spinGather(spinorSCR);
    this->transform(spinorSCR, shell_block_large_mo_[shell_mu], shell_block_small_mo_[shell_nu],
        intermeidate, mat, pOff_size, qOff_size, increment);
  }

  void transformSL(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& intermeidate = interSCR_[thread_id];
    spinorSCR.resize(2 * pauli.nRows(), 2 * pauli.nColumns());
    pauli.spinGather(spinorSCR);
    this->transform(spinorSCR, shell_block_small_mo_[shell_mu], shell_block_large_mo_[shell_nu],
        intermeidate, mat, pOff_size, qOff_size, increment);
  }

  void transformSS(const cqmatrix::PauliSpinorMatrices<MatsT>& pauli, 
      size_t shell_mu, size_t shell_nu, cqmatrix::Matrix<MatsT>& mat, 
      const std::pair<size_t, size_t>& pOff_size,
      const std::pair<size_t, size_t>& qOff_size,
      bool increment = false) const {
    size_t thread_id = GetThreadID();
    auto& spinorSCR = spinorSCR_[thread_id];
    auto& intermeidate = interSCR_[thread_id];
    spinorSCR.resize(2 * pauli.nRows(), 2 * pauli.nColumns());
    pauli.spinGather(spinorSCR);
    this->transform(spinorSCR, shell_block_small_mo_[shell_mu], shell_block_small_mo_[shell_nu],
        intermeidate, mat, pOff_size, qOff_size, increment);
  }
}; // class FourCShellBlockMO

} // namespace ChronusQ
