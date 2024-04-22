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

#include <chronusq_sys.hpp>
#include <hamiltonianoptions.hpp>
#include <libcint.hpp>

namespace ChronusQ {

class LibcintEngine {
 private:
  
  int nAtoms_ = 0;
  int nShells_ = 0;
  int *atom_ = nullptr;
  int *basis_ = nullptr;
  double *env_ = nullptr;
  
  std::vector<double*> cache_; 
  std::vector<size_t> shellBegins_;
  std::vector<size_t> shellSizes_; 
  size_t maxShellSize_ = 0;
   

 public:
  LibcintEngine() = delete;
  LibcintEngine(const LibcintEngine&) = delete;
  LibcintEngine(LibcintEngine &&) = delete;
  
  LibcintEngine(const BasisSet& basisSet, const Molecule& molecule): 
      nAtoms_(molecule.nAtoms), nShells_(basisSet.nShell) {
    
    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    atom_ = CQMemManager::get().malloc<int>(nAtoms_ * ATM_SLOTS);
    basis_ = CQMemManager::get().malloc<int>(nShells_ * BAS_SLOTS);
    env_ = CQMemManager::get().malloc<double>(basisSet.getLibcintEnvLength(molecule));
    basisSet.setLibcintEnv(molecule, atom_, basis_, env_);
    
    shellSizes_.reserve(basisSet.nShell);
    size_t bf = 0ul;
    for (const auto& sh : basisSet.shells) {
      shellBegins_.push_back(bf);
      bf += sh.size();
      shellSizes_.push_back(sh.size());
    }
    maxShellSize_ = *std::max_element(shellSizes_.begin(), shellSizes_.end());
    // std::cout << "libcint nShells_ =" << nShells_ << ", maxShellSize_ =" << maxShellSize_ << std::endl;
  }

  ~LibcintEngine() {
    if (atom_) CQMemManager::get().free(atom_);
    if (basis_) CQMemManager::get().free(basis_);
    if (env_) CQMemManager::get().free(env_);
    for (auto & c : cache_) {
      if (c) CQMemManager::get().free(c);
    }
  }
  
  // getters
  size_t nShells() const { return nShells_; }
  size_t maxShellSize() const { return maxShellSize_; }
  size_t shellSize(size_t i) const { return shellSizes_[i]; } 
  size_t shellBegin(size_t i) const { return shellBegins_[i]; }
  const std::vector<size_t>& shellBegins() const { return shellBegins_; }
  const std::vector<size_t>& shellSizes() const { return shellSizes_; }
  std::pair<size_t, size_t> shellOff(size_t i) const { return {shellBegins_[i], shellSizes_[i]}; }

  // Main interfaces
  void allocate_int2e_cache(const HamiltonianOptions& HOp) {
    size_t cache_size = 0;
    for (int i = 0; i < nShells_; i++) {
      int shells[4]{i, i, i, i};
      if (HOp.BareCoulomb or HOp.DiracCoulomb) {
        cache_size = std::max(cache_size, compute_int2e_sph(nullptr, shells));
      }  
      if (HOp.DiracCoulomb or HOp.DiracCoulombSSSS) {
        cache_size = std::max(cache_size, compute_int2e_ipvip1ipvip2_sph(nullptr, shells));
      } 
      if (HOp.DiracCoulomb) {
        cache_size = std::max(cache_size, compute_int2e_ipvip1_sph(nullptr, shells));
      } 
      if (HOp.Gaunt) {
        cache_size = std::max(cache_size, compute_int2e_ip1ip2_sph(nullptr, shells));
      } 
      if (HOp.Gauge) {
        cache_size = std::max(cache_size, compute_int2e_gauge_r1_ssp1sps2_sph(nullptr, shells));
        cache_size = std::max(cache_size, compute_int2e_gauge_r2_ssp1sps2_sph(nullptr, shells));
        cache_size = std::max(cache_size, compute_int2e_gauge_r1_ssp1ssp2_sph(nullptr, shells));
        cache_size = std::max(cache_size, compute_int2e_gauge_r2_ssp1ssp2_sph(nullptr, shells));
      } 
    }
    
    //std::cout << " cache_size = " << cache_size << std::endl;
    for (auto & c : cache_) {
      if (c) CQMemManager::get().free(c);
    }
    cache_.clear();
    for (auto i = 0ul; i < GetNumThreads(); ++i) {
      cache_.push_back(CQMemManager::get().malloc<double>(cache_size));
    }
  } // allocate_int2e_cache

  // define computing interfaces
#define DEFINE_CQ_LibcintEngine_INT2E_FUNC(function_name)                                    \
  size_t compute_##function_name(double *out, int *shells) const {                           \
    if (out == nullptr) {                                                                    \
      return function_name(nullptr, nullptr, shells, atom_, nAtoms_, basis_, nShells_, env_, \
        nullptr, nullptr);                                                                   \
    }                                                                                        \
    return function_name(out, nullptr, shells, atom_, nAtoms_, basis_, nShells_, env_,       \
        nullptr, cache_[GetThreadID()]);                                                     \
  }

  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_ipvip1ipvip2_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_ipvip1_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_ip1ip2_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_gauge_r1_ssp1sps2_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_gauge_r2_ssp1sps2_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_gauge_r1_ssp1ssp2_sph)
  DEFINE_CQ_LibcintEngine_INT2E_FUNC(int2e_gauge_r2_ssp1ssp2_sph)
  
}; // class LibcintEngine

} // namespace ChronusQ
