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

#include <mointstransformer.hpp>
#include <mointstransformer/shellblockmo.hpp>

#include <libcint/engine.hpp>
#include <util/timer.hpp>
#include <fockbuilder/fourcompfock/batchgd.hpp>
#include <cqlinalg.hpp>
#include <matrix.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasutil.hpp>


// #define _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING

namespace ChronusQ {

namespace {

inline double getMaxShBlkNorm(const cqmatrix::Matrix<double>& shBlkNorm,
    size_t s1, size_t s2, size_t s3, size_t s4) {
  double candidates[6]{shBlkNorm(s1, s2), shBlkNorm(s1, s3), shBlkNorm(s1, s4),
      shBlkNorm(s2, s3), shBlkNorm(s2, s4), shBlkNorm(s3, s4)};
  return *std::max_element(candidates, candidates + 6); 
}

} // namespace

class SchwarzIntegrals {
 private:
  
  // Schwarz Intgrals
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzERI = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzSSSS = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzGaunt = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzGauge = nullptr;
  
 public:
  SchwarzIntegrals() = delete;
  SchwarzIntegrals(const SchwarzIntegrals&) = delete;
  SchwarzIntegrals(SchwarzIntegrals&&) = delete;
  SchwarzIntegrals(const HamiltonianOptions& HOp, const LibcintEngine& cint) {
    computeSchwarzIntegrals(HOp, cint);
  } 
  ~SchwarzIntegrals() { dealloc(); }
  
  void dealloc() {
    SchwarzERI = nullptr;
    SchwarzSSSS = nullptr;
    SchwarzGaunt = nullptr;
    SchwarzGauge = nullptr;
  }
  
  double maxBareCoulomb(size_t s1, size_t s2, size_t s3, size_t s4) const {
    return SchwarzERI->operator()(s1, s2) * SchwarzERI->operator()(s3, s4);  
  }
  double maxDiracCoulomb(size_t s1, size_t s2, size_t s3, size_t s4) const {
    return SchwarzSSSS->operator()(s1, s2) * SchwarzERI->operator()(s3, s4);  
  }
  double maxDiracCoulombSSSS(size_t s1, size_t s2, size_t s3, size_t s4) const {
    return SchwarzSSSS->operator()(s1, s2) * SchwarzSSSS->operator()(s3, s4);  
  }
  double maxGaunt(size_t s1, size_t s2, size_t s3, size_t s4) const {
    return SchwarzGaunt->operator()(s1, s2) * SchwarzGaunt->operator()(s3, s4);  
  }
  double maxGauge(size_t s1, size_t s2, size_t s3, size_t s4) const {
    return SchwarzGauge->operator()(s1, s2) * SchwarzGauge->operator()(s3, s4);  
  }

  void computeSchwarzIntegrals(const HamiltonianOptions& HOp, const LibcintEngine& cint) {
    dealloc();
    size_t maxShellSize = cint.maxShellSize();
    size_t nShells = cint.nShells();
    size_t buffN4 = maxShellSize * maxShellSize * maxShellSize * maxShellSize;
    size_t nThreads = GetNumThreads();

    if (HOp.BareCoulomb or HOp.DiracCoulomb) {
      size_t nERI = 1;
      double *buffAll = CQMemManager::get().malloc<double>(nERI * buffN4 * nThreads);
      SchwarzERI = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzERI->clear(); 
      auto& SchwarzERIMat = *SchwarzERI;

      #pragma omp parallel
      {
        size_t thread_id = GetThreadID();
        size_t n1, n2;
        int shls[4];
        double *buff = buffAll + nERI * buffN4 * thread_id;
        for (size_t s1(0ul), s12(0ul); s1 < nShells; s1++) {
          n1 = cint.shellSize(s1);
          for (size_t s2(0ul); s2 <= s1; s2++, s12++) {
            // Round Robbin work distribution
            #ifdef _OPENMP
            if( s12 % nThreads != thread_id ) continue;
            #endif
            n2 = cint.shellSize(s2);
            shls[0] = s1; 
            shls[1] = s2;
            shls[2] = s1;
            shls[3] = s2;
            if (cint.compute_int2e_sph(buff, shls) == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            double result = std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff, nQuad));
            SchwarzERIMat(s2, s1) = result; 
            SchwarzERIMat(s1, s2) = result; 
          }
        }
      } // parallel region
      CQMemManager::get().free(buffAll);
    } // SchwarzERI
    
    if (HOp.DiracCoulomb or HOp.DiracCoulombSSSS) {
      size_t nERI = 81;
      double *buffAll = CQMemManager::get().malloc<double>(nERI * buffN4 * nThreads);
      SchwarzSSSS = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzSSSS->clear(); 
      auto& SchwarzSSSSMat = *SchwarzSSSS;
      double C2 = 1. / (4 * SpeedOfLight * SpeedOfLight);
      
      #pragma omp parallel
      {
        size_t thread_id = GetThreadID();
        size_t n1, n2;
        int shls[4];
        double *buff = buffAll + nERI * buffN4 * thread_id;
        for (size_t s1(0ul), s12(0ul); s1 < nShells; s1++) {
          n1 = cint.shellSize(s1);
          for (size_t s2(0ul); s2 <= s1; s2++, s12++) {
            // Round Robbin work distribution
            #ifdef _OPENMP
            if( s12 % nThreads != thread_id ) continue;
            #endif
            n2 = cint.shellSize(s2);
            shls[0] = s1; 
            shls[1] = s2;
            shls[2] = s1;
            shls[3] = s2;
            if (cint.compute_int2e_ipvip1ipvip2_sph(buff, shls) == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            double result = C2 * std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff, nQuad));
            SchwarzSSSSMat(s2, s1) = result; 
            SchwarzSSSSMat(s1, s2) = result; 
          }
        }
      } // parallel region
      CQMemManager::get().free(buffAll);
    } // SchwarzSSSS

    if (HOp.Gaunt) {
      size_t nERI = 9;
      double *buffAll = CQMemManager::get().malloc<double>(nERI * buffN4 * nThreads);
      SchwarzGaunt = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzGaunt->clear(); 
      auto& SchwarzGauntMat = *SchwarzGaunt;
      double C1 = 1. / (2 * SpeedOfLight);
      
      #pragma omp parallel
      {
        size_t thread_id = GetThreadID();
        size_t n1, n2;
        int shls[4];
        double *buff = buffAll + nERI * buffN4 * thread_id;
        for (size_t s1(0ul), s12(0ul); s1 < nShells; s1++) {
          n1 = cint.shellSize(s1);
          for (size_t s2(0ul); s2 < nShells; s2++, s12++) {
            // Round Robbin work distribution
            #ifdef _OPENMP
            if( s12 % nThreads != thread_id ) continue;
            #endif
            n2 = cint.shellSize(s2);
            shls[0] = s1; 
            shls[1] = s2;
            shls[2] = s1;
            shls[3] = s2;
            if (cint.compute_int2e_ip1ip2_sph(buff, shls) == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            SchwarzGauntMat(s1, s2) = C1 * std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff, nQuad)); 
          }
        }
      } // parallel region
      CQMemManager::get().free(buffAll);
    } // SchwarzGaunt

    if (HOp.Gauge) {
      size_t nERI = 16;
      double *buffAll = CQMemManager::get().malloc<double>(2 * nERI * buffN4 * nThreads);
      SchwarzGauge = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzGauge->clear(); 
      auto& SchwarzGaugeMat = *SchwarzGauge;
      double C1 = 1. / (2 * SpeedOfLight);
      
      #pragma omp parallel
      {
        size_t thread_id = GetThreadID();
        size_t n1, n2;
        int shls[4];
        double *buff1 = buffAll + nERI * buffN4 * thread_id;
        double *buff2 = buff1 + nERI * buffN4 * nThreads;
        for (size_t s1(0ul), s12(0ul); s1 < nShells; s1++) {
          n1 = cint.shellSize(s1);
          for (size_t s2(0ul); s2 < nShells; s2++, s12++) {
            // Round Robbin work distribution
            #ifdef _OPENMP
            if( s12 % nThreads != thread_id ) continue;
            #endif
            n2 = cint.shellSize(s2);
            shls[0] = s1; 
            shls[1] = s2;
            shls[2] = s1;
            shls[3] = s2;
            auto skip1 = cint.compute_int2e_gauge_r1_ssp1ssp2_sph(buff1, shls);
            auto skip2 = cint.compute_int2e_gauge_r2_ssp1ssp2_sph(buff2, shls);
            if (skip1 == 0 and skip2 == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            auto nBuff =  nQuad * nERI;
            for (auto i = 0ul; i < nBuff; i++) buff1[i] -= buff2[i];
            SchwarzGaugeMat(s1, s2) = C1 * std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff1, nQuad)); 
          }
        }
      } // parallel region
      CQMemManager::get().free(buffAll);
    } // SchwarzGauge
  } // computeSchwarzIntegrals

}; // class SchwarzIntegrals

template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::directTransformTPI(EMPerturbation & pert,
    MatsT* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes) {
   
  size_t poff = off_sizes[0].first;
  size_t qoff = off_sizes[1].first;
  size_t roff = off_sizes[2].first;
  size_t soff = off_sizes[3].first;
  size_t np = off_sizes[0].second;
  size_t nq = off_sizes[1].second;
  size_t nr = off_sizes[2].second;
  size_t ns = off_sizes[3].second;
  size_t npq = np * nq;  

  bool pqSymm = (poff == qoff) and (np == nq); 
  bool rsSymm = (roff == soff) and (nr == ns); 
  
  std::vector<std::pair<size_t,size_t>> rsPairs;
  for (auto s = 0ul; s < ns; s++) 
  for (auto r = 0ul; r < nr; r++) {
    rsPairs.push_back({r, s});
    if (pqSymm and rsSymm and r == s) break;
  } 
  
  std::fill_n(MOTPI, npq * rsPairs.size(), MatsT(0.));
  
  /************************************/
  /* Get env objects from ss          */
  /************************************/
  const auto& mo = ss_.mo[0];
  const auto& nC = ss_.nC;
  auto& HOp = ss_.fockBuilder->hamiltonianOptions_; 
  
  // hack bareCoulomb for 1C and 2C
  if (nC != 4) {
    HOp.BareCoulomb = true; 
  } 
  
  BasisSet basisSet = ss_.basisSet_.groupGeneralContractionBasis(); 
  size_t nShell = basisSet.nShell;
  size_t nShell2 = nShell * nShell;
  size_t nShell4 = nShell2 * nShell2;

  // Set up LibcintEngine  
  //std::cout << "Set up LibcintEngine" << std::endl;
  LibcintEngine cint(basisSet, ss_.molecule_);
  cint.allocate_int2e_cache(HOp);
  size_t maxShellSize = cint.maxShellSize();
  size_t maxShellSize2 = maxShellSize * maxShellSize;
  size_t maxShellSize4 = maxShellSize2 * maxShellSize2;

  // Get the sizes all density SCR
  size_t maxNDenSCR = 0ul;
  size_t n_s34 = 0ul;
  for (auto s3 = 0ul; s3 < nShell; s3++) 
  for (auto s4 = 0ul; s4 <= s3; s4++, n_s34++) { 
    maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
  } // (s3, s4) batches
  maxNDenSCR *= rsPairs.size();
  
#ifdef CQ_ENABLE_MPI
  size_t nNodes =  MPISize(comm_);
  std::vector<size_t> s12Assignment;
  size_t n_s12Assignment = 0ul; 
  if (MPIRank(comm_) == 0) {
    std::vector<size_t> nodeIdHeap(nNodes);
    std::vector<size_t> nodeLoads(nNodes, 0ul);
    std::iota(nodeIdHeap.begin(), nodeIdHeap.end(), 0ul);
    auto comp = [&] (size_t i, size_t j) { 
        return nodeLoads[i] > nodeLoads[j]; 
    };
    
    for (auto s1 = 0ul; s1 < nShell; s1++) 
    for (auto s2 = 0ul; s2 <= s1; s2++) {
      s12Assignment.push_back(nodeIdHeap[0]);
      nodeLoads[nodeIdHeap[0]] += cint.shellSize(s1) * cint.shellSize(s2);
      std::make_heap(nodeIdHeap.begin(), nodeIdHeap.end(), comp);
    }
    // for (auto iNode = 0ul; iNode < nNodes; iNode++) {
    //   std::cout << "iNode = " << iNode << ", loads = " << nodeLoads[iNode] << std::endl; 
    // } 
    n_s12Assignment = s12Assignment.size(); 
  }
  
  MPIBCast(n_s12Assignment, 0, comm_);
  if (MPIRank(comm_) != 0) s12Assignment.resize(n_s12Assignment);
  MPIBCast(s12Assignment.data(), n_s12Assignment, 0, comm_);
  // std::cout << "s12Assignment: " << std::endl; 
  // for (auto s1 = 0ul, s12 = 0ul; s1 < nShell; s1++) 
  // for (auto s2 = 0ul; s2 <= s1; s2++, s12++) {
  //   std::cout << "s1 = " << s1 << ", s2 = " << s2 << ", handdled by Node " << s12Assignment[s12] << std::endl;
  // }
#endif

  // Set up density generator
  std::shared_ptr<ShellBlockMO<MatsT>> shBlockMO = nullptr;
  if (nC == 1) {
    shBlockMO = std::make_shared<OneCShellBlockMO<MatsT>>(mo, cint.shellSizes(), np);
  } else if (nC == 2) {
    shBlockMO = std::make_shared<TwoCShellBlockMO<MatsT>>(mo, cint.shellSizes(), np);
  } else if (nC == 4) {  
    shBlockMO = std::make_shared<FourCShellBlockMO<MatsT>>(mo, cint.shellSizes(), np);
  }
  
  size_t nThreads = GetNumThreads();
  size_t LAThreads = GetLAThreads();
  SetLAThreads(1); // Turn off parallelism in LA functions
  
  /********************************/
  /* Allocate Caches              */
  /********************************/
  
  std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSpinorMSSCRs;
  std::vector<cqmatrix::Matrix<MatsT>> rsERISCRs; 
  for (auto i = 0ul; i < nThreads; ++i) { 
    for (auto j = 0ul; j < rsPairs.size(); ++j) {
      pauliSpinorMSSCRs.emplace_back(maxShellSize, false, false);
    }
    rsERISCRs.emplace_back(np, nq);
  }

  /***********************************/
  /* Prepare Schwarz Screening       */
  /***********************************/
  
  // compute Schwarz ERIs      
  // std::cout << "compute Schwarz Integrals" << std::endl;
  SchwarzIntegrals schwarzInts(HOp, cint);
  auto schwarzThreshold = std::dynamic_pointer_cast<DirectTPI<IntsT>>(ss_.aoints_->TPI)->threshSchwarz();  

  // compute Matrix Norms
  std::vector<cqmatrix::Matrix<double>> shBlkNormsSymmDenLLMS_rs;
  for (auto i = 0ul; i < rsPairs.size(); ++i) { 
    shBlkNormsSymmDenLLMS_rs.emplace_back(nShell); 
  }
  
  
#ifdef CQ_ENABLE_MPI
  size_t shBlkNorm_n = (rsPairs.size() + nNodes - 1) / nNodes; // get the ceilings
  size_t shBlkNorm_begin = shBlkNorm_n * MPIRank(comm_);
  size_t shBlkNorm_end = std::min(shBlkNorm_begin + shBlkNorm_n, rsPairs.size());
#else
  size_t shBlkNorm_n = rsPairs.size();
  size_t shBlkNorm_begin = 0ul;
  size_t shBlkNorm_end = rsPairs.size();
#endif        
  // std::cout << "compute shBlkNorm: nTask = " << shBlkNorm_n << std::endl;
  // std::cout << "shBlkNorm_begin = " << shBlkNorm_begin << ", shBlkNorm_end =" << shBlkNorm_end << std::endl;  
  
  #pragma omp parallel for
  for(auto i = shBlkNorm_begin; i < shBlkNorm_end; ++i) {
    const auto& [r, s] = rsPairs[i];
    auto& shBlkNorm = shBlkNormsSymmDenLLMS_rs[i];
    auto& pauli = pauliSpinorMSSCRs[GetThreadID()]; 

    // FIXME genSymmDenLLMS is ALWAYS symmetric
    for (auto s1 = 0ul; s1 < nShell; s1++) 
    for (auto s2 = 0ul; s2 <= s1; s2++) {
      shBlockMO->genSymmDenLLMS(r + roff, s + soff, s1, s2, pauli);
      double result = pauli.norm(lapack::Norm::Inf);
      shBlkNorm(s2, s1) = result;
      shBlkNorm(s1, s2) = result; 
    }
    // shBlkNorm.output(std::cout, "shBlkNorm[" + std::to_string(i) + "]", true); 
  } // [r, s]
  
#ifdef CQ_ENABLE_MPI
  ProgramTimer::tick("MOINTSTRANSFORM TPI TRANS MPI COMM");
  for(auto i = 0ul; i < rsPairs.size(); ++i) {
    int root = i / shBlkNorm_n;
    // std::cout << "Broadcast shBlkNorm " << i << ", handdled by Node " << root  << std::endl;
    auto& shBlkNorm = shBlkNormsSymmDenLLMS_rs[i];
    MPIBCast(shBlkNorm.pointer(), nShell2, root, comm_);
  }
  ProgramTimer::tock("MOINTSTRANSFORM TPI TRANS MPI COMM");
#endif

  auto maxShBlkNormsSymmDenLLMS_rs = shBlkNormsSymmDenLLMS_rs[0];
  
  #pragma omp parallel for
  for (auto s1 = 0ul; s1 < nShell; s1++) 
  for (auto s2 = 0ul; s2 <= s1; s2++) {
    for (auto i = 1ul; i < rsPairs.size(); ++i) {
      maxShBlkNormsSymmDenLLMS_rs(s1, s2) = std::max(
          maxShBlkNormsSymmDenLLMS_rs(s1, s2), shBlkNormsSymmDenLLMS_rs[i](s1, s2));
    } 
    maxShBlkNormsSymmDenLLMS_rs(s2, s1) = maxShBlkNormsSymmDenLLMS_rs(s1, s2);
  }

  // auto maxShBlkNormsSymmDenLLMS_pq = maxShBlkNormsSymmDenLLMS_rs;
  // if (poff != roff or np != nr or qoff != soff or nq != ns) {
  //   for (auto s1 = 0ul; s1 < nShell; s1++) 
  //   for (auto s2 = 0ul; s2 <= s1; s2++) {
  //     double result = 0.;
  //     auto& pauli = pauliSpinorMSSCRs[GetThreadID()];
  //     for (auto q = 0ul; q < nq; q++) 
  //     for (auto p = 0ul; p < np; p++) {
  //       shBlockMO->genSymmDenLLMS(p + poff, q + qoff, s1, s2, pauli);
  //       result = std::max(result,  pauli.norm(lapack::Norm::Inf));
  //       if (pqSymm and rsSymm and p == q) break; 
  //     }
  //     maxShBlkNormsSymmDenLLMS_pq(s1, s2) = result;
  //     maxShBlkNormsSymmDenLLMS_pq(s2, s1) = result;
  //   }
  // }
  // maxShBlkNormsSymmDenLLMS_rs.output(std::cout, "maxShBlkNormsSymmDenLLMS_rs", true);

  // Bare Coulomb
  if (HOp.BareCoulomb) {

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::vector<double> tInts_all(nThreads, 0.);
    std::vector<double> t1_2_all(nThreads, 0.);
    std::vector<double> t2_2_all(nThreads, 0.);
    std::vector<double> tDensity_all(nThreads, 0.);
    std::vector<double> tUpdate_all(nThreads, 0.);
#endif
    
    double *buffERIAll = CQMemManager::get().malloc<double>(maxShellSize4 * nThreads);
    size_t availbleMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availbleMem);
    MatsT *buffDensity =  CQMemManager::get().malloc<MatsT>(nDenSCR);
    
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;
    std::vector<std::vector<MatsT*>> s43DenPtrsAll; // 
    s34PairsAll.resize(nThreads);
    s43DenPtrsAll.resize(nThreads);
    
    for (auto i_s34 = 0ul; i_s34 < n_s34; ) { 
      
      /**************************************/
      /*   Form AO Densities                */
      /**************************************/ 
      
      // FIXME:try to do better parallellism here
      // try best to evenly distribute the workloads across different threads 
      std::vector<size_t> threadIdHeap(nThreads);
      std::vector<size_t> threadLoads(nThreads, 0ul);
      std::iota(threadIdHeap.begin(), threadIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) {
        return threadLoads[i] > threadLoads[j];
      };
      
      for (auto iThread = 0ul; iThread < nThreads; iThread++) {
        s34PairsAll[iThread].clear();
        s43DenPtrsAll[iThread].clear();
      }
      
      size_t nMem = 0ul, nMemOff = 0ul;
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell or nMem < nDenSCR; s3++) 
      for (auto s4 = 0ul; s4 <= s3; s4++, s34++) {
        if (s34 < i_s34) continue;
        size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
        nMem += nDen34SCR;
        if (nMem > nDenSCR) break;
        
        // assigning to a thread
        size_t iThread = threadIdHeap[0];
        threadLoads[iThread] += nDen34SCR;
        std::make_heap(threadIdHeap.begin(), threadIdHeap.end(), comp);
        
        s34PairsAll[iThread].push_back({s3, s4});
        s43DenPtrsAll[iThread].push_back(buffDensity + nMemOff);
        nMemOff += nDen34SCR;

        i_s34++;
      } // s34 assignment 
      
      // std::cout << "i_s34 = " << i_s34 << std::endl;
      // for (auto iThread = 0ul; iThread < nThreads; iThread++) {
      //   std::cout << "iThread = " << iThread << ", loads = " << threadLoads[iThread] 
      //             << ", nTasks = " << s34PairsAll[iThread].size() << std::endl; 
      // } 

      #pragma omp parallel 
      {
        int thread_id = GetThreadID();
        const auto& s34Pairs = s34PairsAll[thread_id]; 
        const auto& s43DenPtrs = s43DenPtrsAll[thread_id]; 
        auto& denSCR = pauliSpinorMSSCRs[thread_id]; 
        
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        auto& tDensity = tDensity_all[thread_id];
        auto topDensity = tick();
#endif

        for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
          MatsT* denPtr = s43DenPtrs[s34];
          const auto& [s3, s4] = s34Pairs[s34];
          size_t nsh34 = cint.shellSize(s3) * cint.shellSize(s4); 
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++, denPtr+=nsh34) {  
            const auto& [r, s] = rsPairs[iMat];
            shBlockMO->genSymmDenLLMS(r + roff, s + soff, s3, s4, denSCR);
            std::copy_n(denSCR.S().pointer(), nsh34, denPtr); 
          }
        }

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        tDensity += tock(topDensity);
#endif
 
      } // end of the parallel region
      
      for (auto s1 = 0ul, s12 = 0ul; s1 < nShell; s1++) 
      for (auto s2 = 0ul; s2 <= s1; s2++, s12++) { 
        
#ifdef CQ_ENABLE_MPI
      if (s12Assignment[s12] != MPIRank(comm_)) continue;
#endif        
        size_t n1 = cint.shellSize(s1);
        size_t n2 = cint.shellSize(s2);      
        
        // std::cout << " s1 =" << s1 << ", s2 = " << s2 << std::endl;

        /**************************************/
        /*  First Half Transformation         */
        /**************************************/ 
        #pragma omp parallel
        {
          auto topFirstHalf = tick();
          
          int thread_id = GetThreadID();
          const auto& s34Pairs = s34PairsAll[thread_id]; 
          const auto& s43DenPtrs = s43DenPtrsAll[thread_id]; 
          double *buff = buffERIAll + maxShellSize4 * thread_id;
          
          int shls[4]; 
          shls[0] = int(s1); 
          shls[1] = int(s2); 
          
          // initialize cache
          size_t ADLL12_loc_off = thread_id * rsPairs.size();
          for (auto iMat = ADLL12_loc_off; iMat < ADLL12_loc_off + rsPairs.size(); iMat++) {
            pauliSpinorMSSCRs[iMat].resize(n1, n2);  
            pauliSpinorMSSCRs[iMat].clear();  
          } 
          
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
          auto& tInts = tInts_all[thread_id];
          auto& t1_2 = t1_2_all[thread_id];
#endif
          
          for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
            const auto& [s3, s4] = s34Pairs[s34];
            
            auto maxBareCoulomb = schwarzInts.maxBareCoulomb(s1, s2, s3, s4);
            if (getMaxShBlkNorm(maxShBlkNormsSymmDenLLMS_rs, s1, s2, s3, s4) *
                maxBareCoulomb < schwarzThreshold) continue;
            
            shls[2] = int(s3); 
            shls[3] = int(s4);
            size_t n3 = cint.shellSize(s3);
            size_t n4 = cint.shellSize(s4);
            size_t nsh34 = n3 * n4;
            double s12_deg = (s1 == s2) ? 1.0 : 2.0;
            double s34_deg = (s3 == s4) ? 1.0 : 2.0;
            double s1234_deg = s12_deg * s34_deg * 0.5;
             
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto topInts = tick();
#endif

            if (cint.compute_int2e_sph(buff, shls) == 0) continue;
            auto nQuad = n1 * n2 * n3 * n4;
            for(auto i = 0ul; i < nQuad; i++) buff[i] *= s1234_deg;

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            tInts += tock(topInts);
            auto top1_2 = tick();
#endif
            
            MatsT* symmDLLMS43_ptr = s43DenPtrs[s34];
            for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++, symmDLLMS43_ptr+=nsh34) {
              if (getMaxShBlkNorm(shBlkNormsSymmDenLLMS_rs[iMat], s1, s2, s3, s4) 
                  * maxBareCoulomb < schwarzThreshold) continue;
              
              const auto& [r, s] = rsPairs[iMat];
              auto& ADLL12 = pauliSpinorMSSCRs[ADLL12_loc_off + iMat]; 
              auto& ADLLMS12 = ADLL12.S();
              
              // FIXME: make this as GEMM call
              for(auto l = 0ul, mnkl=0ul ; l < n4; ++l) 
              for(auto k = 0ul; k < n3; ++k) 
              for(auto n = 0ul; n < n2; ++n) 
              for(auto m = 0ul; m < n1; ++m, ++mnkl) {
                ADLLMS12(m, n) += buff[mnkl] * symmDLLMS43_ptr[l + k * n4]; 
              } // mnkl
            } // iMat
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            t1_2 += tock(top1_2);
#endif
          } // (s3, s4)
          // #pragma omp critical
          // std::cout << "thread_id = " << thread_id << ", time =" << tock(topFirstHalf) << std::endl;
        } // end of parallel region    
      
        /**************************************/
        /*  Second Half Transformation        */
        /************************************ */ 
        #pragma omp parallel
        { 
          int thread_id = GetThreadID();
          auto& rsERI = rsERISCRs[thread_id];
 
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
          auto& t2_2 = t2_2_all[thread_id];
          auto& tUpdate = tUpdate_all[thread_id];
#endif

          #pragma omp for
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++) {   

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto top2_2 = tick();
#endif
            auto& ADLL12 = pauliSpinorMSSCRs[iMat]; 
            auto& ADLLMS12 = ADLL12.S();
            for (auto iThread = 1ul; iThread < nThreads; iThread++) {
              ADLL12 += pauliSpinorMSSCRs[iMat + iThread * rsPairs.size()];    
            }

            
            shBlockMO->transformLL(ADLL12, s1, s2, rsERI, off_sizes[0], off_sizes[1]);  
            ADLLMS12.inplace_T();
            shBlockMO->transformLL(ADLL12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);  

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            t2_2 += tock(top2_2);
            auto topUpdate = tick();
#endif
            
            MatsT scale = MatsT(1.0);
            blas::axpy(npq, scale, rsERI.pointer(), 1, MOTPI + iMat * npq, 1); 
            
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            tUpdate += tock(topUpdate);
#endif
          } // iMat
        } // end of parallel region
      } // (s1, s2)    
    } // (s3, s4) Density Batching 
    
    CQMemManager::get().free(buffERIAll, buffDensity);
    
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    auto printTimings = [] (const std::string& section, 
        const std::vector<double> ts) {
        std::cout << std::setw(20) << section  << ": " 
                  << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                  << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                  << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;   
    };  

    std::cout << "\nTiming for fully direct transformation: " << std::endl;
    printTimings("t(Ints)", tInts_all);
    printTimings("t(Density)", tDensity_all);
    printTimings("t(1/2)", t1_2_all);
    printTimings("t(2/2)", t2_2_all);
    printTimings("t(Update)", tUpdate_all);
    std::cout << std::endl;   
#endif
  
  } // Bare Coulomb

  if (nC != 4) {
    HOp.BareCoulomb = false;
  }

  SetLAThreads(LAThreads);// Turn threads for LA back on

#ifdef CQ_ENABLE_MPI
  ProgramTimer::tick("MOINTSTRANSFORM TPI TRANS MPI COMM");
  for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++) { 
    MPIAllReduce(MOTPI + iMat * npq, npq, MOTPI + iMat * npq, comm_);
  }
  ProgramTimer::tock("MOINTSTRANSFORM TPI TRANS MPI COMM");
#endif
  
  // std::cout << "ERI Norm (with rs symmetry) = " << std::setprecision(16) << 
  //    lapack::lange(lapack::Norm::Fro, npq, rsPairs.size(), MOTPI, npq) << std::endl;
  
  // restore to no symmetry
  if (pqSymm and rsSymm) {
    // if pqSymm and rsSymm:  r <= s upper triagnle 
    for (int iMat = rsPairs.size() - 1; iMat >= 0; iMat--) {
      const auto& [r, s] = rsPairs[iMat];
      SetMat('N', np, nq, MatsT(1.), MOTPI + npq * iMat, np, MOTPI + npq * (r + s * nr), np);
    }
    for (const auto& [r, s] : rsPairs) {
      if (s != r) { 
        SetMat('C', np, nq, MatsT(1.), MOTPI + npq * (r + s * nr), np, MOTPI + npq * (s + r * nr), np);
      }
    }
  }
  
  // check Matrix Norm of ERI
  //std::cout << "ERI Norm = " << std::setprecision(16) << 
  //   lapack::lange(lapack::Norm::Fro, npq, nr * ns, MOTPI, npq) << std::endl;
  
  return;
} // directTransformTPI

} // namespace ChronusQ
