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


#define _MOINTSTRANSFORMER_TPI_FULL_DIRECT_SCALAR_TIMING

namespace ChronusQ {

class SchwarzIntegralsScalar {
 private:
  
  // Schwarz Intgrals
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzERI = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzSSSS = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzGaunt = nullptr;
  std::shared_ptr<cqmatrix::Matrix<double>> SchwarzGauge = nullptr;
  
 public:
  SchwarzIntegralsScalar() = delete;
  SchwarzIntegralsScalar(const SchwarzIntegralsScalar&) = delete;
  SchwarzIntegralsScalar(SchwarzIntegralsScalar&&) = delete;
  SchwarzIntegralsScalar(const HamiltonianOptions& HOp, const LibcintEngine& cint) {
    computeSchwarzIntegrals(HOp, cint);
  } 
  ~SchwarzIntegralsScalar() { dealloc(); }
  
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
      size_t nERI = 1;
      double *buffAll = CQMemManager::get().malloc<double>(nERI * buffN4 * nThreads);
      SchwarzSSSS = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzSSSS->clear(); 
      auto& SchwarzSSSSMat = *SchwarzSSSS;
      double C2 = 1. / (4 * SpeedOfLight() * SpeedOfLight());
      
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
            if (cint.compute_int2e_pp1pp2_sph(buff, shls) == 0) continue;
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
      size_t nERI = 1;
      double *buffAll = CQMemManager::get().malloc<double>(nERI * buffN4 * nThreads);
      SchwarzGaunt = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzGaunt->clear(); 
      auto& SchwarzGauntMat = *SchwarzGaunt;
      double C1 = 1. / (2 * SpeedOfLight());
      
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
            if (cint.compute_int2e_gaunt_ps1ps2_sph(buff, shls) == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            SchwarzGauntMat(s1, s2) = C1 * std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff, nQuad)); 
          }
        }
      } // parallel region

      // Sym 
      for (auto s1 = 0ul; s1 < nShells; s1++)
      for (auto s2 = 0ul; s2 <= s1; s2++) {
        double schwarzMax = std::max(SchwarzGauntMat(s1, s2), SchwarzGauntMat(s2, s1));
        SchwarzGauntMat(s1, s2) = schwarzMax;
        SchwarzGauntMat(s2, s1) = schwarzMax;
      }
      // SchwarzGauntMat.output(std::cout, "gaunt", true);
      CQMemManager::get().free(buffAll);
    } // SchwarzGaunt

    if (HOp.Gauge) {
      size_t nERI = 4;
      double *buffAll = CQMemManager::get().malloc<double>(2 * nERI * buffN4 * nThreads);
      SchwarzGauge = std::make_shared<cqmatrix::Matrix<double>>(nShells);
      SchwarzGauge->clear(); 
      auto& SchwarzGaugeMat = *SchwarzGauge;
      double C1 = 1. / (2 * SpeedOfLight());
      
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
            auto skip1 = cint.compute_int2e_gauge_r1_sp1sp2_sph(buff1, shls);
            auto skip2 = cint.compute_int2e_gauge_r2_sp1sp2_sph(buff2, shls);
            if (skip1 == 0 and skip2 == 0) continue;
            auto nQuad = n1 * n2 * n1 * n2;
            auto nBuff =  nQuad * nERI;
            for (auto i = 0ul; i < nBuff; i++) buff1[i] -= buff2[i];
            SchwarzGaugeMat(s1, s2) = C1 * std::sqrt(lapack::lange(lapack::Norm::Max, nQuad, nERI, buff1, nQuad)); 
          }
        }
      } // parallel region

      // Sym 
      for (auto s1 = 0ul; s1 < nShells; s1++)
      for (auto s2 = 0ul; s2 <= s1; s2++) {
        double schwarzMax = std::max(SchwarzGaugeMat(s1, s2), SchwarzGaugeMat(s2, s1));
        SchwarzGaugeMat(s1, s2) = schwarzMax;
        SchwarzGaugeMat(s2, s1) = schwarzMax;
      }
    //  SchwarzGaugeMat.output(std::cout, "gauge", true);
      CQMemManager::get().free(buffAll);
    } // SchwarzGauge
  } // computeSchwarzIntegrals

}; // class SchwarzIntegrals


template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::directTransformScalarTPIBatch(EMPerturbation & pert,
    MatsT* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes) const {
   
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
    // XSLI: why pqSymm is needed here? Hang only considered {p}={q}={r}={s}
    rsPairs.push_back({r, s});
    if (pqSymm and rsSymm and r == s) break;
  }
  
  std::fill_n(MOTPI, npq * rsPairs.size(), MatsT(0.));
  
  /************************************/
  /* Get env objects from ss          */
  /************************************/
  // XSLI: this will not work for UHF where mo[0] and mo[1] are needed
  const auto& mo = ss_.mo[0];
  const auto& nC = ss_.nC;
  auto& HOp = ss_.fockBuilder->hamiltonianOptions_; 
  
  // hack bareCoulomb for 1C and 2C
  if (nC != 4) {
    HOp.BareCoulomb = true; 
  } 
  // else {
  //   HOp.Gaunt = false;
  //   HOp.Gauge = false;
  // }

  // XSLI: why GeneralContractionBasis is needed here?
  BasisSet basisSet = ss_.basisSet_.groupGeneralContractionBasis(); 
  size_t nShell = basisSet.nShell;
  size_t nShell2 = nShell * nShell;
  size_t nShell4 = nShell2 * nShell2;

  // Set up LibcintEngine  
  LibcintEngine cint(basisSet, ss_.molecule_);
  cint.allocate_int2eScalar_cache(HOp);
  size_t maxShellSize = cint.maxShellSize();
  size_t maxShellSize2 = maxShellSize * maxShellSize;
  size_t maxShellSize4 = maxShellSize2 * maxShellSize2;

  #ifdef CQ_ENABLE_MPI
    size_t nNodes =  MPISize(comm_);
    std::vector<size_t> s12Assignment;
    size_t n_s12Assignment = 0ul; 
    if (MPIRank(comm_) == 0) {
      std::vector<size_t> nodeIdHeap(nNodes);
      std::vector<size_t> nodeLoads(nNodes, 0ul);
      // fill node ID from 0 to nNodes-1
      std::iota(nodeIdHeap.begin(), nodeIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) { 
          return nodeLoads[i] > nodeLoads[j]; 
      };

      // estimate loads based on the number of basis pairs in each shell pair.
      // cost of high angular momentum shell is not considered
      if (HOp.Gaunt or HOp.Gauge) {
        for (auto s1 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 < nShell; s2++) {
          s12Assignment.push_back(nodeIdHeap[0]);
          nodeLoads[nodeIdHeap[0]] += cint.shellSize(s1) * cint.shellSize(s2);
          std::make_heap(nodeIdHeap.begin(), nodeIdHeap.end(), comp);
        }
      } else if(HOp.BareCoulomb or HOp.DiracCoulomb or HOp.DiracCoulombSSSS) {
        for (auto s1 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 <= s1; s2++) {
          s12Assignment.push_back(nodeIdHeap[0]);
          nodeLoads[nodeIdHeap[0]] += cint.shellSize(s1) * cint.shellSize(s2);
          std::make_heap(nodeIdHeap.begin(), nodeIdHeap.end(), comp);
        }  
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

  // MO: for debug multi vectors could be one
  std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSpinorLLMSSCRs;
  std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSpinorSSSCRs;
  std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> pauliSpinorSLSCRs;
  std::vector<cqmatrix::Matrix<MatsT>> rsERISCRs; 
  for (auto i = 0ul; i < nThreads; ++i) { 
    pauliSpinorLLMSSCRs.emplace_back(maxShellSize, false, false);
    pauliSpinorSSSCRs.emplace_back(maxShellSize, false, false);
    pauliSpinorSLSCRs.emplace_back(maxShellSize, true, true);
    rsERISCRs.emplace_back(np, nq);
  }

  /***********************************/
  /* Prepare Schwarz Screening       */
  /***********************************/
  // compute Schwarz ERIs      
  SchwarzIntegralsScalar schwarzInts(HOp, cint);
  auto schwarzThreshold = std::dynamic_pointer_cast<DirectTPI<IntsT>>(ss_.aoints_->TPI)->threshSchwarz();  

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

  // compute Matrix Norms
  std::vector<cqmatrix::Matrix<double>> shBlkNormsSymmDenLLMS_rs;
  std::vector<cqmatrix::Matrix<double>> shBlkNormsSymmDenSS_rs;
  std::vector<cqmatrix::Matrix<double>> shBlkNormsSymmDenSL_rs;
  for (auto i = 0ul; i < rsPairs.size(); ++i) {
    shBlkNormsSymmDenLLMS_rs.emplace_back(nShell, nShell);
    if (HOp.DiracCoulomb or HOp.DiracCoulombSSSS) {
      shBlkNormsSymmDenSS_rs.emplace_back(nShell, nShell); 
    } if (HOp.Gaunt or HOp.Gauge) {
      shBlkNormsSymmDenSL_rs.emplace_back(nShell, nShell); 
    }
  } // rsPair
  // :)
  cqmatrix::Matrix<double> maxShBlkNormsSymmDenLLMS_rs(nShell, nShell);
  cqmatrix::Matrix<double> maxShBlkNormsSymmDenSS_rs(nShell, nShell);
  cqmatrix::Matrix<double> maxShBlkNormsSymmDenSL_rs(nShell, nShell);

  // BareCoulomb and DiracCoulomb
  maxShBlkNormsSymmDenLLMS_rs = shBlkNormsSymmDenLLMS_rs[0];
  #pragma omp parallel for
  for(auto i = shBlkNorm_begin; i < shBlkNorm_end; ++i) {
    const auto& [r, s] = rsPairs[i];

      auto& shBlkNorm = shBlkNormsSymmDenLLMS_rs[i];
      auto& pauli = pauliSpinorLLMSSCRs[GetThreadID()];

      for (auto s1 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 < nShell; s2++) {

          shBlockMO->genSymmDenLLMS(r + roff, s + soff, s1, s2, pauli);
          double result = pauli.norm(lapack::Norm::Inf);
          shBlkNorm(s2, s1) = result;
          shBlkNorm(s1, s2) = result;
        } // (s1, s2)
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
  
  maxShBlkNormsSymmDenLLMS_rs = shBlkNormsSymmDenLLMS_rs[0];
  #pragma omp parallel for
  for (auto s1 = 0ul; s1 < nShell; s1++) 
  for (auto s2 = 0ul; s2 <= s1; s2++) {
    for (auto i = 1ul; i < rsPairs.size(); ++i) {
      maxShBlkNormsSymmDenLLMS_rs(s1, s2) = std::max(
          maxShBlkNormsSymmDenLLMS_rs(s1, s2), shBlkNormsSymmDenLLMS_rs[i](s1, s2));
    } // rsPairs
    maxShBlkNormsSymmDenLLMS_rs(s2, s1) = maxShBlkNormsSymmDenLLMS_rs(s1, s2);
  } // (s1, s2)

  if (HOp.DiracCoulomb or HOp.DiracCoulombSSSS) {
    maxShBlkNormsSymmDenSS_rs = shBlkNormsSymmDenSS_rs[0];

    #pragma omp parallel for
    for(auto i = shBlkNorm_begin; i < shBlkNorm_end; ++i) {
      const auto& [r, s] = rsPairs[i];

        auto& shBlkNorm = shBlkNormsSymmDenSS_rs[i];
        auto& pauli = pauliSpinorSSSCRs[GetThreadID()];
        for (auto s1 = 0ul; s1 < nShell; s1++)
          for (auto s2 = 0ul; s2 <= s1; s2++) {
            shBlockMO->genSymmDenSSMS(r + roff, s + soff, s1, s2, pauli);
            double result = pauli.norm(lapack::Norm::Inf);
            shBlkNorm(s2, s1) = result;
            shBlkNorm(s1, s2) = result;
        } // (s1, s2)
    } // [r, s]
    #ifdef CQ_ENABLE_MPI
    for(auto i = 0ul; i < rsPairs.size(); ++i) {
      int root = i / shBlkNorm_n;
      // std::cout << "Broadcast shBlkNorm " << i << ", handdled by Node " << root  << std::endl;
      auto& shBlkNorm = shBlkNormsSymmDenSS_rs[i];
      MPIBCast(shBlkNorm.pointer(), nShell2, root, comm_);
    }
    #endif
    
    maxShBlkNormsSymmDenSS_rs = shBlkNormsSymmDenSS_rs[0];
    #pragma omp parallel for
    for (auto s1 = 0ul; s1 < nShell; s1++) 
    for (auto s2 = 0ul; s2 <= s1; s2++) {
      for (auto i = 1ul; i < rsPairs.size(); ++i) {
        maxShBlkNormsSymmDenSS_rs(s1, s2) = std::max(
            maxShBlkNormsSymmDenSS_rs(s1, s2), shBlkNormsSymmDenSS_rs[i](s1, s2));
      } // rsPairs
      maxShBlkNormsSymmDenSS_rs(s2, s1) = maxShBlkNormsSymmDenSS_rs(s1, s2);
    } // (s1, s2)
  } // DiracCoulomb or DiracCoulombSSSS

  if (HOp.Gaunt or HOp.Gauge) {
    maxShBlkNormsSymmDenSL_rs = shBlkNormsSymmDenSL_rs[0];

      #pragma omp parallel for
      for(auto i = shBlkNorm_begin; i < shBlkNorm_end; ++i) {
        const auto& [r, s] = rsPairs[i];

          auto& shBlkNorm = shBlkNormsSymmDenSL_rs[i];
          auto& pauli = pauliSpinorSLSCRs[GetThreadID()];

          // Get the max norm of s1 and s2 to keep symetry in screening
          for (auto s1 = 0ul; s1 < nShell; s1++)
            for (auto s2 = 0ul; s2 <= s1; s2++) {
              shBlockMO->genDenLSpmDenSL(r + roff, s + soff, s1, s2, pauli);
              double result12 = pauli.norm(lapack::Norm::Inf);
              shBlockMO->genDenLSpmDenSL(r + roff, s + soff, s2, s1, pauli);
              double result21 = pauli.norm(lapack::Norm::Inf);
              double result = std::max(result12, result21);
              shBlkNorm(s2, s1) = result;
              shBlkNorm(s1, s2) = result;
            }

        // shBlkNorm.output(std::cout, "shBlkNormSL[" + std::to_string(i) + "]", true); 
      } // [r, s]

    #ifdef CQ_ENABLE_MPI
    for(auto i = 0ul; i < rsPairs.size(); ++i) {
      int root = i / shBlkNorm_n;
      // std::cout << "Broadcast shBlkNorm " << i << ", handdled by Node " << root  << std::endl;
      auto& shBlkNorm = shBlkNormsSymmDenSL_rs[i];
      MPIBCast(shBlkNorm.pointer(), nShell2, root, comm_);
    }
    #endif

    maxShBlkNormsSymmDenSL_rs = shBlkNormsSymmDenSL_rs[0];
    
    #pragma omp parallel for
    for (auto s1 = 0ul; s1 < nShell; s1++) 
    for (auto s2 = 0ul; s2 <= s1; s2++) {
      for (auto i = 1ul; i < rsPairs.size(); ++i) {
        maxShBlkNormsSymmDenSL_rs(s1, s2) = std::max(
            maxShBlkNormsSymmDenSL_rs(s1, s2), shBlkNormsSymmDenSL_rs[i](s1, s2));
      } // rsPairs
      maxShBlkNormsSymmDenSL_rs(s2, s1) = maxShBlkNormsSymmDenSL_rs(s1, s2);
    } // (s1, s2)
   } // Gaunt or Gauge

 
  /***********************************/
  /*                                 */
  /* Start of Bare-Coulomb           */
  /*                                 */
  /***********************************/
  if (HOp.BareCoulomb) {
    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::cout << "  Transforming 2e-INTS: BareCoulomb ..." << std::endl;
    #endif
    auto startTime = tick();
    // Get the sizes all density SCR
    size_t maxNDenSCR = 0ul;
    size_t n_s34 = 0ul;
    for (auto s3 = 0ul; s3 < nShell; s3++)
      for (auto s4 = 0ul; s4 <= s3; s4++, n_s34++) {
        maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
      } // (s3, s4) batches
    maxNDenSCR *= rsPairs.size();

    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::vector<double> tInts_all(nThreads, 0.);
    std::vector<double> t1_2_all(nThreads, 0.);
    std::vector<double> t2_2_all(nThreads, 0.);
    std::vector<double> tDensity_all(nThreads, 0.);
    std::vector<double> tUpdate_all(nThreads, 0.);

    std::vector<size_t> nConSkipped(nThreads, 0);
    std::vector<size_t> nIntSkipped(nThreads, 0);
    #endif   
    
    double *buffERIAll = CQMemManager::get().malloc<double>(maxShellSize4 * nThreads);
    size_t availableMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availableMem);
    MatsT *buffDensity =  CQMemManager::get().malloc<MatsT>(nDenSCR);
    
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;
    std::vector<std::vector<MatsT*>> s43DenPtrsAll; // 
    s34PairsAll.resize(nThreads);
    s43DenPtrsAll.resize(nThreads);

    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorLLMSSCRs;
    for (auto i = 0ul; i < nThreads; ++i) {
      for (auto j = 0ul; j < rsPairs.size(); ++j) {
        s12SpinorLLMSSCRs.emplace_back(maxShellSize, false, false);
      }
    }

    // This is the top level loop over s34 pairs
    for (auto i_s34 = 0ul; i_s34 < n_s34; ) {
      
      /**************************************/
      /*   Form AO Densities                */
      /**************************************/ 
      
      // FIXME:try to do better parallelism here
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
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell and nMem < nDenSCR; s3++)
      for (auto s4 = 0ul; s4 <= s3; s4++, s34++) {
        if (s34 < i_s34) continue;
        size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
        nMem += nDen34SCR;
        // MO: This seems redudnat see s3 loop 5 lines up
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
        auto& denSCR = s12SpinorLLMSSCRs[thread_id]; 
        
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
            // get all rs density for a given s3s4 shell pair
            // C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, large component, assuming real integrals
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
            s12SpinorLLMSSCRs[iMat].resize(n1, n2);  
            s12SpinorLLMSSCRs[iMat].clear();  
          } 

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto& tInts = tInts_all[thread_id];
            auto& t1_2 = t1_2_all[thread_id];
  #endif
          
          for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
            const auto& [s3, s4] = s34Pairs[s34];
            
            auto maxBareCoulomb = schwarzInts.maxBareCoulomb(s1, s2, s3, s4);
            // MO: Screening
            if (getMaxShBlkNorm(maxShBlkNormsSymmDenLLMS_rs, s1, s2, s3, s4) *
                maxBareCoulomb < schwarzThreshold) {
                  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                  nIntSkipped[thread_id] ++;
                  #endif
                  continue;
                } 

            shls[2] = int(s3); 
            shls[3] = int(s4);
            size_t n3 = cint.shellSize(s3);
            size_t n4 = cint.shellSize(s4);
            size_t nsh34 = n3 * n4;
            double s12_deg = (s1 == s2) ? 1.0 : 2.0;
            double s34_deg = (s3 == s4) ? 1.0 : 2.0;
            // double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
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

              const auto& [r, s] = rsPairs[iMat];

              if (getMaxShBlkNorm(shBlkNormsSymmDenLLMS_rs[iMat], s1, s2, s3, s4) 
                  * maxBareCoulomb < schwarzThreshold) {
                    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                    nConSkipped[thread_id] ++;
                    #endif
                    continue;
                  }
              
              auto& ADLL12 = s12SpinorLLMSSCRs[ADLL12_loc_off + iMat]; 
              auto& ADLLMS12 = ADLL12.S();
              
              // FIXME: make this as GEMM call
              // First half transform
              for(auto l = 0ul, mnkl=0ul ; l < n4; ++l) 
              for(auto k = 0ul; k < n3; ++k) 
              for(auto n = 0ul; n < n2; ++n) 
              for(auto m = 0ul; m < n1; ++m, ++mnkl) {
                // This is because Hang only used 4-fold symmetry here. The symmetry between 12 and 34 is not considered yet.
                ADLLMS12(m, n) += buff[mnkl] * symmDLLMS43_ptr[l + k * n4];
              } // mnkl

              // MO: to make this a blas GEMV the resutling ADLLMS12 should be reshaped.
              // The raw eri has the shape of n1*n2 x n3*n4 and thus the reuslting matrix is flattend n1*n2 wich should be reshaped into n x m
              // blas::gemv(blas::Layout::ColMajor, blas::Op::NoTrans, n1*n2, n3*n4, MatsT(1.0), buff, 
              //   n3*n4, symmDLLMS43_ptr, 1, 0.0, ADLLMS12.pointer(), 1);

            } // iMat
  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t1_2 += tock(top1_2);
  #endif
          } // (s3, s4)
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
            auto& ADLL12 = s12SpinorLLMSSCRs[iMat]; 
            auto& ADLLMS12 = ADLL12.S();
            for (auto iThread = 1ul; iThread < nThreads; iThread++) {
              ADLL12 += s12SpinorLLMSSCRs[iMat + iThread * rsPairs.size()];    
            }

            // MO: the two transforms are to get the symetric MO
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

    double totalTime = tock(startTime);
    FormattedLine(std::cout, "Time to Bare-Coulomb Transform(s): ", totalTime); 
    
  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
      auto printTimings = [] (const std::string& section, 
          const std::vector<double> ts) {
          std::cout << std::setw(20) << section  << ": " 
                    << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                    << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                    << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;   
      };  

      size_t nIntSkippedAcc = std::accumulate(nIntSkipped.begin(),nIntSkipped.end(),0);
      size_t nConSkippedAcc = std::accumulate(nConSkipped.begin(),nConSkipped.end(),0);

      #ifdef CQ_ENABLE_MPI
      MPIAllReduce(&nIntSkippedAcc, 1, &nIntSkippedAcc, comm_);
      MPIAllReduce(&nConSkippedAcc, 1, &nConSkippedAcc, comm_);
      #endif

      std::cout << "\nTiming for bare Coulomb fully direct transformation: " << std::endl;
      std::cout << "Total time: " << totalTime << "s" << std::endl;
      std::cout << "Integrals skipped: " << nIntSkippedAcc << std::endl;
      std::cout << "Contractions skipped: " << nConSkippedAcc << std::endl;
      printTimings("t(Ints)", tInts_all);
      printTimings("t(Density)", tDensity_all);
      printTimings("t(1/2)", t1_2_all);
      printTimings("t(2/2)", t2_2_all);
      printTimings("t(Update)", tUpdate_all);
      std::cout << std::endl;   
  #endif
  
  } // Bare Coulomb
  /*********************************/
  /*                               */
  /* End of Bare-Coulomb           */
  /*                               */
  /*********************************/


  /******************************************/
  /*                                        */
  /* Start of Dirac-Coulomb C(2)            */
  /* includes DC-LLLL and DC-LLSS/SSLL      */
  /*                                        */
  /******************************************/

  // Check if MatsT is dcomplex at compile time aka only compile this blok of code for the dcomplex
  // Needed for comp with nr code
  
  if (HOp.DiracCoulomb) {
    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::cout << "  Transforming 2e-INTS: DiracCoulomb (LL|LL, LS|LS) ..." << std::endl;
    #endif
    auto startTime = tick();
    enum ERI_2ND_DERIV {
      AxBx,
      AxBy,
      AxBz,
      AyBx,
      AyBy,
      AyBz,
      AzBx,
      AzBy,
      AzBz
    };

    int nERI = 1;
    int nSave = 4;
    size_t NB  = maxShellSize*4;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;
    size_t NB4_2 = 2*NB4;
    size_t NB4_3 = 3*NB4;

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
      std::vector<double> tInts_all(nThreads, 0.);
      std::vector<double> t1_2_all(nThreads, 0.);
      std::vector<double> t2_2_all(nThreads, 0.);
      std::vector<double> tDensity_all(nThreads, 0.);
      std::vector<double> tUpdate_all(nThreads, 0.);

      std::vector<size_t> nConSkipped(nThreads, 0);
      std::vector<size_t> nIntSkipped(nThreads, 0);
  #endif

    // Allocate memory for raw ERI (2 buff)
    double *buffERIAll = CQMemManager::get().template malloc<double>(2*nERI*maxShellSize4*nThreads);
    // Allocate memory for assembled ERI, i.e., dot product and cross product (4 buff)
    int nBuff = 2;
    double *ERIBuffer = CQMemManager::get().template malloc<double>(nSave*nBuff*NB4*nThreads);

    // Get the sizes all density SCR
    size_t maxNDenSCR = 0ul;
    size_t n_s34 = 0ul;
    for (auto s3 = 0ul; s3 < nShell; s3++)
      for (auto s4 = 0ul; s4 <= s3; s4++, n_s34++) {
        maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
      } // (s3, s4) batches
    maxNDenSCR *= 5*rsPairs.size();

    size_t availableMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availableMem);

    // Allocate memory for density to be contracted with ERI
    MatsT *buffDensityLLMS =  CQMemManager::get().template malloc<MatsT>(nDenSCR/5);
    MatsT *buffDensitySSMS =  CQMemManager::get().template malloc<MatsT>(nDenSCR/5);

    // storage to save shell IDs of s3 and s4 for each thread
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;

    // storage to save s34 densities of all rs pairs for each thread
    std::vector<std::vector<MatsT*>> s43DenLLMSPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSSMSPtrsAll;

    s34PairsAll.resize(nThreads);
    s43DenLLMSPtrsAll.resize(nThreads);
    s43DenSSMSPtrsAll.resize(nThreads);

    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorLLSCRs;
    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorSSSCRs;
    for (auto i = 0ul; i < nThreads; ++i) {
      for (auto j = 0ul; j < rsPairs.size(); ++j) {
        s12SpinorLLSCRs.emplace_back(maxShellSize, false, false);
        s12SpinorSSSCRs.emplace_back(maxShellSize, false, false);
      }
    }

    for (auto i_s34 = 0ul; i_s34 < n_s34; ) {

      /**************************************/
      /*   Form AO Densities                */
      /**************************************/

      // FIXME:try to do better parallelism here
      // try best to evenly distribute the workloads across different threads
      std::vector<size_t> threadIdHeap(nThreads);
      std::vector<size_t> threadLoads(nThreads, 0ul);
      std::iota(threadIdHeap.begin(), threadIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) {
          return threadLoads[i] > threadLoads[j];
      };

      for (auto iThread = 0ul; iThread < nThreads; iThread++) {
        s34PairsAll[iThread].clear();
        s43DenLLMSPtrsAll[iThread].clear();
        s43DenSSMSPtrsAll[iThread].clear();
      }

      size_t nMem = 0ul, nMemOff = 0ul;
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell and nMem < nDenSCR; s3++)
      for (auto s4 = 0ul; s4 <= s3; s4++, s34++) {
        if (s34 < i_s34) continue;
        size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
        nMem += nDen34SCR*5;
        if (nMem > nDenSCR) break;

        // assigning to a thread
        size_t iThread = threadIdHeap[0];
        threadLoads[iThread] += nDen34SCR;
        std::make_heap(threadIdHeap.begin(), threadIdHeap.end(), comp);

        s34PairsAll[iThread].push_back({s3, s4});
        s43DenLLMSPtrsAll[iThread].push_back(buffDensityLLMS + nMemOff);
        s43DenSSMSPtrsAll[iThread].push_back(buffDensitySSMS + nMemOff);

        nMemOff += nDen34SCR;

        i_s34++;
      } // s34 assignment

  #pragma omp parallel
      {
        int thread_id = GetThreadID();
        const auto& s34Pairs = s34PairsAll[thread_id];
        const auto& s43DenLLMSPtrs = s43DenLLMSPtrsAll[thread_id];
        const auto& s43DenSSMSPtrs = s43DenSSMSPtrsAll[thread_id];
        auto& denLLSCR = s12SpinorLLSCRs[thread_id];
        auto& denSSSCR = s12SpinorSSSCRs[thread_id];

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        auto& tDensity = tDensity_all[thread_id];
        auto topDensity = tick();
  #endif

        for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
          MatsT* denPtrLLMS = s43DenLLMSPtrs[s34];
          MatsT* denPtrSSMS = s43DenSSMSPtrs[s34];
          const auto& [s3, s4] = s34Pairs[s34];
          size_t nsh34 = cint.shellSize(s3) * cint.shellSize(s4);
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++, denPtrLLMS+=nsh34,
                  denPtrSSMS+=nsh34) {
            const auto& [r, s] = rsPairs[iMat];

            // S: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, large component, assuming real integrals
            shBlockMO->genSymmDenLLMS(r + roff, s + soff, s3, s4, denLLSCR);
            std::copy_n(denLLSCR.S().pointer(), nsh34, denPtrLLMS);

            shBlockMO->genSymmDenSSMS(r + roff, s + soff, s3, s4, denSSSCR);
            std::copy_n(denSSSCR.S().pointer(), nsh34, denPtrSSMS);

          } // iMat
        } // s34Pairs

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
          tDensity += tock(topDensity);
  #endif

      } // end of the parallel region

      for (auto s1 = 0ul, s12 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 <= s1; s2++, s12++) {

  #ifdef CQ_ENABLE_MPI
            if (s12Assignment[s12] != MPIRank(comm_)) continue;
  #endif


          // std::cout << " s1 =" << s1 << ", s2 = " << s2 << std::endl;

          /**************************************/
          /*  First Half Transformation         */
          /**************************************/
  #pragma omp parallel
          {
            size_t n1 = cint.shellSize(s1);
            size_t n2 = cint.shellSize(s2);

            auto topFirstHalf = tick();

            int thread_id = GetThreadID();
            const auto& s34Pairs = s34PairsAll[thread_id];
            const auto& s43DenLLMSPtrs = s43DenLLMSPtrsAll[thread_id];
            const auto& s43DenSSMSPtrs = s43DenSSMSPtrsAll[thread_id];
            //const auto& s43DenSSMXPtrs = s43DenSSMXPtrsAll[thread_id];
            //const auto& s43DenSSMYPtrs = s43DenSSMYPtrsAll[thread_id];
            //const auto& s43DenSSMZPtrs = s43DenSSMZPtrsAll[thread_id];

            double *buff1 = buffERIAll + nERI * maxShellSize4 * thread_id;
            double *buff2 = buffERIAll + nERI * maxShellSize4 * nThreads + nERI * maxShellSize4 * thread_id;
            double *ERIBuffABmn   = &ERIBuffer[thread_id*nSave*NB4];
            double *ERIBuffCDkl   = &ERIBuffer[nThreads*nSave*NB4 + thread_id*nSave*NB4];

            int shls[4];

            // initialize storage for (𝜇𝜈|rs) where 𝜇∈s1 and 𝜈∈s2
            size_t AD12_loc_off = thread_id * rsPairs.size();
            for (auto iMat = AD12_loc_off; iMat < AD12_loc_off + rsPairs.size(); iMat++) {
              s12SpinorLLSCRs[iMat].resize(n1, n2);
              s12SpinorLLSCRs[iMat].clear();
              s12SpinorSSSCRs[iMat].resize(n1, n2);
              s12SpinorSSSCRs[iMat].clear();
            }

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              auto& tInts = tInts_all[thread_id];
            auto& t1_2 = t1_2_all[thread_id];
  #endif

            for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
              const auto& [s3, s4] = s34Pairs[s34];

              // get the max schwarz for both LLSS and SSLL
              double maxDiracCoulomb = std::max(schwarzInts.maxDiracCoulomb(s1, s2, s3, s4), schwarzInts.maxDiracCoulomb(s3, s4, s1, s2));

              // MO: Screening
              if (std::max(getMaxShBlkNorm(maxShBlkNormsSymmDenLLMS_rs, s1, s2, s3, s4), getMaxShBlkNorm(maxShBlkNormsSymmDenSS_rs, s1, s2, s3, s4))
                  * maxDiracCoulomb < schwarzThreshold) {
                  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                  nIntSkipped[thread_id] ++;
                  #endif
                  continue;
                 }
                 
              shls[0] = int(s1);
              shls[1] = int(s2);
              shls[2] = int(s3);
              shls[3] = int(s4);
              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              size_t n3 = cint.shellSize(s3);
              size_t n4 = cint.shellSize(s4);

              size_t nsh34 = n3 * n4;

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                auto topInts = tick();
  #endif

              // ∇A∙∇B(mn|kl); int2e_pp1 is the spin-free analogue of int2e_ipvip1
              // (derivatives on the bra pair only). int2e_pp1pp2 would be the SSSS
              // operator and is wrong here.
              bool skip1 = (cint.compute_int2e_pp1_sph(buff1, shls) == 0);
              auto nQuad = n1 * n2 * n3 * n4;
              //for(auto i = 0ul; i < nQuad; i++) buff[i] *= s1234_deg;

              // swap r1 and r2 so that we can compute ∇C∇D(mn|kl) integrals
              shls[0] = int(s3);
              shls[1] = int(s4);
              shls[2] = int(s1);
              shls[3] = int(s2);

              // MO: done using the same engine by (kl|mn)
              // ∇C∇D(mn|kl)
              // bool skip2 = (cint.compute_int2e_ipvip1_sph(buff2, shls) == 0);
              bool skip2 = (cint.compute_int2e_pp1_sph(buff2, shls) == 0);
              
              // FIXME: separate both buff builds
              // Only skip if both cint calls are empty
              if (skip1 && skip2) continue;

              double s34_deg = (s3 == s4) ? 0.5 : 1.0;
              double s12_deg = (s1 == s2) ? 0.5 : 1.0;

              for(auto l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
              for(auto k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k)
              for(auto n =   maxShellSize            ; n <   maxShellSize + n2; ++n)
              for(auto m = 0                                ; m <                  n1; ++m, ++mnkl) {

                /* Dirac-Coulomb */
                // ∇A∙∇B(mn|kl): int2e_pp1 returns a single component with the
                // dot product already contracted inside CINTgout2e_int2e_pp1
                auto dAdotdB = buff1[mnkl];
                // ∇Ax∇B(mn|kl)
                //auto dAcrossdB_x =  buff1[AyBz*nQuad+mnkl] - buff1[AzBy*nQuad+mnkl];
                //auto dAcrossdB_y = -buff1[AxBz*nQuad+mnkl] + buff1[AzBx*nQuad+mnkl];
                //auto dAcrossdB_z =  buff1[AxBy*nQuad+mnkl] - buff1[AyBx*nQuad+mnkl];

                auto MNKL = m + n*NB + k*NB2 + l*NB3;
                auto KLMN = k + l*NB + m*NB2 + n*NB3;

                // ∇A∙∇B(mn|kl) followed by ∇Ax∇B(mn|kl) X, Y, and Z
                // (mn|kl)
                ERIBuffABmn[       MNKL] =  s12_deg*s34_deg*dAdotdB;
                //ERIBuffABmn[   NB4+MNKL] =  s12_deg*s34_deg*dAcrossdB_x;
                //ERIBuffABmn[ NB4_2+MNKL] =  s12_deg*s34_deg*dAcrossdB_y;
                //ERIBuffABmn[ NB4_3+MNKL] =  s12_deg*s34_deg*dAcrossdB_z;

              } // ∇A∇B integral preparation loop

              for(auto n =   maxShellSize, mnkl = 0ul ; n <   maxShellSize + n2; ++n)
              for(auto m = 0                                 ; m <                  n1; ++m)
              for(auto l = 3*maxShellSize             ; l < 3*maxShellSize + n4; ++l)
              for(auto k = 2*maxShellSize             ; k < 2*maxShellSize + n3; ++k, ++mnkl) {
                /* Dirac-Coulomb */
                // ∇C∙∇D(mn|kl): single component, dot product already contracted
                auto dCdotdD = buff2[mnkl];
                // ∇Cx∇D(mn|kl)
                //auto dCcrossdD_x =  buff2[AyBz*nQuad+mnkl] - buff2[AzBy*nQuad+mnkl];
                //auto dCcrossdD_y = -buff2[AxBz*nQuad+mnkl] + buff2[AzBx*nQuad+mnkl];
                //auto dCcrossdD_z =  buff2[AxBy*nQuad+mnkl] - buff2[AyBx*nQuad+mnkl];

                auto MNKL = m + n*NB + k*NB2 + l*NB3;
                auto KLMN = k + l*NB + m*NB2 + n*NB3;

                // ∇C∙∇D(mn|kl) followed by ∇Cx∇D(mn|kl) X, Y, and Z
                // (mn|kl)
                ERIBuffCDkl[       MNKL] =  s12_deg*s34_deg*dCdotdD;
                //ERIBuffCDkl[   NB4+MNKL] =  s12_deg*s34_deg*dCcrossdD_x;
                //ERIBuffCDkl[ NB4_2+MNKL] =  s12_deg*s34_deg*dCcrossdD_y;
                //ERIBuffCDkl[ NB4_3+MNKL] =  s12_deg*s34_deg*dCcrossdD_z;

              } // ∇C∇D integral preparation loop

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                tInts += tock(topInts);
              auto top1_2 = tick();
  #endif

              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              n3 = cint.shellSize(s3);
              n4 = cint.shellSize(s4);
              
              MatsT* symmDLLMS43_ptr = s43DenLLMSPtrs[s34];
              MatsT* symmDSSMS43_ptr = s43DenSSMSPtrs[s34];
              //MatsT* symmDSSMX43_ptr = s43DenSSMXPtrs[s34];
              //MatsT* symmDSSMY43_ptr = s43DenSSMYPtrs[s34];
              //MatsT* symmDSSMZ43_ptr = s43DenSSMZPtrs[s34];


              for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++, symmDLLMS43_ptr+=nsh34,
                      symmDSSMS43_ptr+=nsh34) {

                const auto& [r, s] = rsPairs[iMat];

                // MO: Screening
                  if (std::max(getMaxShBlkNorm(shBlkNormsSymmDenLLMS_rs[iMat], s1, s2, s3, s4), getMaxShBlkNorm(shBlkNormsSymmDenSS_rs[iMat], s1, s2, s3, s4))
                      * maxDiracCoulomb < schwarzThreshold) {
                    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                    nConSkipped[thread_id] ++;
                    #endif
                    continue;
                   }

                auto& ADSS12 = s12SpinorSSSCRs[AD12_loc_off + iMat];
                auto& ADSSMS12 = ADSS12.S();
                //auto& ADSSMX12 = ADSS12.X();
                //auto& ADSSMY12 = ADSS12.Y();
                //auto& ADSSMZ12 = ADSS12.Z();

                auto& ADLL12 = s12SpinorLLSCRs[AD12_loc_off + iMat];
                auto& ADLLMS12 = ADLL12.S();


                for(auto m = 0ul, ms = 0ul; m < n1 and ms < n1; ++m, ++ms)
                for(auto n =   maxShellSize, ns = 0ul; n <   maxShellSize + n2 and ns < n2; ++n, ++ns)
                for(auto k = 2*maxShellSize, ks = 0ul; k < 2*maxShellSize + n3 and ks < n3; ++k, ++ks)
                for(auto l = 3*maxShellSize, ls = 0ul; l < 3*maxShellSize + n4 and ls < n4; ++l,++ls) {

                  auto MNKL = m + n*NB + k*NB2 + l*NB3;

                  auto DotPrdMNKL = MNKL;
                  //auto CrossXMNKL = MNKL+NB4;
                  //auto CrossYMNKL = MNKL+NB4_2;
                  //auto CrossZMNKL = MNKL+NB4_3;

                  ADSSMS12(ms,ns) += ERIBuffABmn[DotPrdMNKL]*symmDLLMS43_ptr[ls+ks*n4];
                  //ADSSMX12(ms,ns) += ERIBuffABmn[CrossXMNKL]*symmDLLMS43_ptr[ls+ks*n4];
                  //ADSSMY12(ms,ns) += ERIBuffABmn[CrossYMNKL]*symmDLLMS43_ptr[ls+ks*n4];
                  //ADSSMZ12(ms,ns) += ERIBuffABmn[CrossZMNKL]*symmDLLMS43_ptr[ls+ks*n4];
        
                  ADLLMS12(ms,ns) += ERIBuffCDkl[DotPrdMNKL] * symmDSSMS43_ptr[ls+ks*n4];
                                  //+ (ERIBuffCDkl[CrossXMNKL] * symmDSSMX43_ptr[ls+ks*n4]
                                  //+  ERIBuffCDkl[CrossYMNKL] * symmDSSMY43_ptr[ls+ks*n4] 
                                  //+  ERIBuffCDkl[CrossZMNKL] * symmDSSMZ43_ptr[ls+ks*n4] ) * dcomplex(0.,1.);

                } // mnkl
              } // iMat
  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                t1_2 += tock(top1_2);
  #endif
            } // (s3, s4)
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
              const auto& [r, s] = rsPairs[iMat];

              auto& ADLL12 = s12SpinorLLSCRs[iMat];
              auto& ADSS12 = s12SpinorSSSCRs[iMat];
              for (auto iThread = 1ul; iThread < nThreads; iThread++) {
                s12SpinorLLSCRs[iMat] += s12SpinorLLSCRs[iMat + iThread * rsPairs.size()];
                s12SpinorSSSCRs[iMat] += s12SpinorSSSCRs[iMat + iThread * rsPairs.size()];
              }

              // Second half tranformation 
              // TransfromLL. Contraction with the LL density
              shBlockMO->transformLL(ADLL12, s1, s2, rsERI, off_sizes[0], off_sizes[1]);
              ADLL12.S().inplace_T();
              shBlockMO->transformLL(ADLL12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);

              shBlockMO->transformSS(ADSS12, s1, s2, rsERI, off_sizes[0], off_sizes[1], true);
              ADSS12.S().inplace_T();
              shBlockMO->transformSS(ADSS12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                t2_2 += tock(top2_2);
              auto topUpdate = tick();
  #endif

              // Scaling of the eri. 2 * 1/4c^2
              MatsT C2 = 2./(4*SpeedOfLight()*SpeedOfLight());
              MatsT scale = C2;
              blas::axpy(npq, scale, rsERI.pointer(), 1, MOTPI + iMat * npq, 1);

  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                tUpdate += tock(topUpdate);
  #endif
            } // iMat
          } // end of parallel region
        } // (s1, s2)
    } // (s3, s4) Density Batching

    CQMemManager::get().free(buffERIAll, ERIBuffer, buffDensityLLMS,
             buffDensitySSMS);


    double totalTime = tock(startTime);
    FormattedLine(std::cout, "Time to DC(LL/LS) Transform(s): ", totalTime);
  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING

      auto printTimings = [] (const std::string& section,
          const std::vector<double> ts) {
          std::cout << std::setw(20) << section  << ": "
                    << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                    << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                    << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;
      };

      size_t nIntSkippedAcc = std::accumulate(nIntSkipped.begin(),nIntSkipped.end(),0);
      size_t nConSkippedAcc = std::accumulate(nConSkipped.begin(),nConSkipped.end(),0);

      #ifdef CQ_ENABLE_MPI
      MPIAllReduce(&nIntSkippedAcc, 1, &nIntSkippedAcc, comm_);
      MPIAllReduce(&nConSkippedAcc, 1, &nConSkippedAcc, comm_);
      #endif

      std::cout << "\nTiming for Dirac Coulomb fully direct transformation: " << std::endl;
      std::cout << "Total time: " << totalTime << "s" << std::endl;
      std::cout << "Integrals skipped: " << nIntSkippedAcc << std::endl;
      std::cout << "Contractions skipped: " << nConSkippedAcc << std::endl;
      printTimings("t(Ints)", tInts_all);
      printTimings("t(Density)", tDensity_all);
      printTimings("t(1/2)", t1_2_all);
      printTimings("t(2/2)", t2_2_all);
      printTimings("t(Update)", tUpdate_all);
      std::cout << std::endl;
  #endif

  } // DiracCoulomb

  /******************************************/
  /*                                        */
  /*   End of Dirac-Coulomb LL and C(2)-SS  */
  /*                                        */
  /******************************************/


  /****************************************/
  /*                                      */
  /* Start of Dirac-Coulomb C(4)-SSSS     */
  /*                                      */
  /****************************************/
  if (HOp.DiracCoulombSSSS) {
    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::cout << "  Transforming 2e-INTS: DiracCoulomb (SS|SS) ..." << std::endl;
    #endif
    auto startTime = tick();
    enum ERI_4TH_DERIV {
        AxBxCxDx,
        AxBxCxDy,
        AxBxCxDz,
        AxBxCyDx,
        AxBxCyDy,
        AxBxCyDz,
        AxBxCzDx,
        AxBxCzDy,
        AxBxCzDz,
        AxByCxDx,
        AxByCxDy,
        AxByCxDz,
        AxByCyDx,
        AxByCyDy,
        AxByCyDz,
        AxByCzDx,
        AxByCzDy,
        AxByCzDz,
        AxBzCxDx,
        AxBzCxDy,
        AxBzCxDz,
        AxBzCyDx,
        AxBzCyDy,
        AxBzCyDz,
        AxBzCzDx,
        AxBzCzDy,
        AxBzCzDz,
        AyBxCxDx,
        AyBxCxDy,
        AyBxCxDz,
        AyBxCyDx,
        AyBxCyDy,
        AyBxCyDz,
        AyBxCzDx,
        AyBxCzDy,
        AyBxCzDz,
        AyByCxDx,
        AyByCxDy,
        AyByCxDz,
        AyByCyDx,
        AyByCyDy,
        AyByCyDz,
        AyByCzDx,
        AyByCzDy,
        AyByCzDz,
        AyBzCxDx,
        AyBzCxDy,
        AyBzCxDz,
        AyBzCyDx,
        AyBzCyDy,
        AyBzCyDz,
        AyBzCzDx,
        AyBzCzDy,
        AyBzCzDz,
        AzBxCxDx,
        AzBxCxDy,
        AzBxCxDz,
        AzBxCyDx,
        AzBxCyDy,
        AzBxCyDz,
        AzBxCzDx,
        AzBxCzDy,
        AzBxCzDz,
        AzByCxDx,
        AzByCxDy,
        AzByCxDz,
        AzByCyDx,
        AzByCyDy,
        AzByCyDz,
        AzByCzDx,
        AzByCzDy,
        AzByCzDz,
        AzBzCxDx,
        AzBzCxDy,
        AzBzCxDz,
        AzBzCyDx,
        AzBzCyDy,
        AzBzCyDz,
        AzBzCzDx,
        AzBzCzDy,
        AzBzCzDz,
    };



    size_t NB  = maxShellSize*4;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;
    size_t NB4_2 = 2*NB4;
    size_t NB4_3 = 3*NB4;
    size_t NB4_4  = 4*NB4;
    size_t NB4_5  = 5*NB4;
    size_t NB4_6  = 6*NB4;
    size_t NB4_7  = 7*NB4;
    size_t NB4_8  = 8*NB4;
    size_t NB4_9  = 9*NB4;
    size_t NB4_10 =10*NB4;
    size_t NB4_11 =11*NB4;
    size_t NB4_12 =12*NB4;
    size_t NB4_13 =13*NB4;
    size_t NB4_14 =14*NB4;
    size_t NB4_15 =15*NB4;

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::vector<double> tInts_all(nThreads, 0.);
    std::vector<double> t1_2_all(nThreads, 0.);
    std::vector<double> t2_2_all(nThreads, 0.);
    std::vector<double> tDensity_all(nThreads, 0.);
    std::vector<double> tUpdate_all(nThreads, 0.);

    std::vector<size_t> nConSkipped(nThreads, 0);
    std::vector<size_t> nIntSkipped(nThreads, 0);
#endif
    int nERI = 1;
    int nSave = 16;
    // Allocate memory for raw ERI
    double *buffERIAll = CQMemManager::get().template malloc<double>(nERI*maxShellSize4*nThreads);
    // Allocate memory for assembled ERI, i.e., dot product and cross product
    double *ERIBuffer = CQMemManager::get().template malloc<double>(nSave*NB4*nThreads);

    size_t maxNDenSCR = 0ul;
    size_t n_s34 = 0ul;
    for (auto s3 = 0ul; s3 < nShell; s3++)
      for (auto s4 = 0ul; s4 <= s3; s4++, n_s34++) {
        maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
      } // (s3, s4) batches
    maxNDenSCR *= 4*rsPairs.size();

    size_t availableMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availableMem);

    // Allocate memory for density to be contracted with ERI
    MatsT *buffDensitySSMS =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    //MatsT *buffDensitySSMX =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    //MatsT *buffDensitySSMY =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    //MatsT *buffDensitySSMZ =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);

    // storage to save shell IDs of s3 and s4 for each thread
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;

    // storage to save s34 densities of all rs pairs for each thread
    std::vector<std::vector<MatsT*>> s43DenSSMSPtrsAll;
    //std::vector<std::vector<MatsT*>> s43DenSSMXPtrsAll;
    //std::vector<std::vector<MatsT*>> s43DenSSMYPtrsAll;
    //std::vector<std::vector<MatsT*>> s43DenSSMZPtrsAll;

    s34PairsAll.resize(nThreads);
    s43DenSSMSPtrsAll.resize(nThreads);
    //s43DenSSMXPtrsAll.resize(nThreads);
    //s43DenSSMYPtrsAll.resize(nThreads);
    //s43DenSSMZPtrsAll.resize(nThreads);

    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorSSSCRs;
    for (auto i = 0ul; i < nThreads; ++i) {
      for (auto j = 0ul; j < rsPairs.size(); ++j) {
        s12SpinorSSSCRs.emplace_back(maxShellSize, false, false);
      }
    }

    for (auto i_s34 = 0ul; i_s34 < n_s34; ) {

      /**************************************/
      /*   Form AO Densities                */
      /**************************************/

      // FIXME:try to do better parallelism here
      // try best to evenly distribute the workloads across different threads
      std::vector<size_t> threadIdHeap(nThreads);
      std::vector<size_t> threadLoads(nThreads, 0ul);
      std::iota(threadIdHeap.begin(), threadIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) {
          return threadLoads[i] > threadLoads[j];
      };

      for (auto iThread = 0ul; iThread < nThreads; iThread++) {
        s34PairsAll[iThread].clear();
        s43DenSSMSPtrsAll[iThread].clear();
        //s43DenSSMXPtrsAll[iThread].clear();
        //s43DenSSMYPtrsAll[iThread].clear();
        //s43DenSSMZPtrsAll[iThread].clear();
      }

      size_t nMem = 0ul, nMemOff = 0ul;
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell and nMem < nDenSCR; s3++)
        for (auto s4 = 0ul; s4 <= s3; s4++, s34++) {
          if (s34 < i_s34) continue;
          size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
          nMem += nDen34SCR*4;
          if (nMem > nDenSCR) break;

          // assigning to a thread
          size_t iThread = threadIdHeap[0];
          threadLoads[iThread] += nDen34SCR;
          std::make_heap(threadIdHeap.begin(), threadIdHeap.end(), comp);

          s34PairsAll[iThread].push_back({s3, s4});
          s43DenSSMSPtrsAll[iThread].push_back(buffDensitySSMS + nMemOff);
          //s43DenSSMXPtrsAll[iThread].push_back(buffDensitySSMX + nMemOff);
          //s43DenSSMYPtrsAll[iThread].push_back(buffDensitySSMY + nMemOff);
          //s43DenSSMZPtrsAll[iThread].push_back(buffDensitySSMZ + nMemOff);

          nMemOff += nDen34SCR;

          i_s34++;
        } // s34 assignment

#pragma omp parallel
      {
        int thread_id = GetThreadID();
        const auto& s34Pairs = s34PairsAll[thread_id];
        const auto& s43DenSSMSPtrs = s43DenSSMSPtrsAll[thread_id];
        //const auto& s43DenSSMXPtrs = s43DenSSMXPtrsAll[thread_id];
        //const auto& s43DenSSMYPtrs = s43DenSSMYPtrsAll[thread_id];
        //const auto& s43DenSSMZPtrs = s43DenSSMZPtrsAll[thread_id];
        auto& denSSSCR = s12SpinorSSSCRs[thread_id];

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        auto& tDensity = tDensity_all[thread_id];
        auto topDensity = tick();
#endif

        for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
          MatsT* denPtrSSMS = s43DenSSMSPtrs[s34];
          //MatsT* denPtrSSMX = s43DenSSMXPtrs[s34];
          //MatsT* denPtrSSMY = s43DenSSMYPtrs[s34];
          //MatsT* denPtrSSMZ = s43DenSSMZPtrs[s34];
          const auto& [s3, s4] = s34Pairs[s34];
          size_t nsh34 = cint.shellSize(s3) * cint.shellSize(s4);
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
               denPtrSSMS+=nsh34) {
            const auto& [r, s] = rsPairs[iMat];

            shBlockMO->genSymmDenSSMS(r + roff, s + soff, s3, s4, denSSSCR);
            // S: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSSSCR.S().pointer(), nsh34, denPtrSSMS);
            // X: C_s,s3xC_r,s4^* - C_r,s3^*xC_s,s4, small component, assuming real integrals
            //std::copy_n(denSSSCR.X().pointer(), nsh34, denPtrSSMX);
            // Y: C_s,s3xC_r,s4^* - C_r,s3^*xC_s,s4, small component, assuming real integrals
            //std::copy_n(denSSSCR.Y().pointer(), nsh34, denPtrSSMY);
            // Z: C_s,s3xC_r,s4^* - C_r,s3^*xC_s,s4, small component, assuming real integrals
            //std::copy_n(denSSSCR.Z().pointer(), nsh34, denPtrSSMZ);
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

          /**************************************/
          /*  First Half Transformation         */
          /**************************************/
#pragma omp parallel
          {
            size_t n1 = cint.shellSize(s1);
            size_t n2 = cint.shellSize(s2);

            auto topFirstHalf = tick();

            int thread_id = GetThreadID();
            const auto& s34Pairs = s34PairsAll[thread_id];
            const auto& s43DenSSMSPtrs = s43DenSSMSPtrsAll[thread_id];
            //const auto& s43DenSSMXPtrs = s43DenSSMXPtrsAll[thread_id];
            //const auto& s43DenSSMYPtrs = s43DenSSMYPtrsAll[thread_id];
            //const auto& s43DenSSMZPtrs = s43DenSSMZPtrsAll[thread_id];

            double *buff = buffERIAll + nERI*maxShellSize4 * thread_id;
            double *ERIBuffABCD   = &ERIBuffer[thread_id*nSave*NB4];

            int shls[4];
            shls[0] = int(s1);
            shls[1] = int(s2);

            // initialize storage for (𝜇𝜈|rs) where 𝜇∈s1 and 𝜈∈s2
            size_t AD12_loc_off = thread_id * rsPairs.size();
            for (auto iMat = AD12_loc_off; iMat < AD12_loc_off + rsPairs.size(); iMat++) {
              s12SpinorSSSCRs[iMat].resize(n1, n2);
              s12SpinorSSSCRs[iMat].clear();
            }

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto& tInts = tInts_all[thread_id];
          auto& t1_2 = t1_2_all[thread_id];
#endif

            for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
              const auto& [s3, s4] = s34Pairs[s34];

              // MO: Screening
              auto maxDiracCoulombSSSS = schwarzInts.maxDiracCoulombSSSS(s1, s2, s3, s4);
              if (getMaxShBlkNorm(maxShBlkNormsSymmDenSS_rs, s1, s2, s3, s4) *
                 maxDiracCoulombSSSS < schwarzThreshold) {
                  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                  nIntSkipped[thread_id] ++;
                  #endif
                  continue;
                 }

              shls[0] = int(s1);
              shls[1] = int(s2);
              shls[2] = int(s3);
              shls[3] = int(s4);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              auto topInts = tick();
#endif

              // if (cint.compute_int2e_ipvip1ipvip2_sph(buff, shls) == 0) continue;
              if (cint.compute_int2e_pp1pp2_sph(buff, shls) == 0) continue;
              
              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              size_t n3 = cint.shellSize(s3);
              size_t n4 = cint.shellSize(s4);
              size_t nsh34 = n3 * n4;
              auto nQuad = n1 * n2 * n3 * n4;

              for(auto l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
              for(auto k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k)
              for(auto n =   maxShellSize            ; n <   maxShellSize + n2; ++n)
              for(auto m = 0                         ; m <                  n1; ++m, ++mnkl) {

                // (∇A∙∇B)(∇C∙∇D)(mnkl): int2e_pp1pp2 returns a single component with
                // both dot products already contracted inside CINTgout2e_int2e_pp1pp2
                auto dAdotdBdCdotdD = buff[mnkl];

                // (∇Ax∇B)(∇C∙∇D)(mnkl)
                //auto dAcrossdB_xdCdotdD =  buff[AyBzCxDx*nQuad+mnkl] - buff[AzByCxDx*nQuad+mnkl]
                //                           + buff[AyBzCyDy*nQuad+mnkl] - buff[AzByCyDy*nQuad+mnkl]
                //                           + buff[AyBzCzDz*nQuad+mnkl] - buff[AzByCzDz*nQuad+mnkl];

                //auto dAcrossdB_ydCdotdD = -buff[AxBzCxDx*nQuad+mnkl] + buff[AzBxCxDx*nQuad+mnkl]
                //                          -buff[AxBzCyDy*nQuad+mnkl] + buff[AzBxCyDy*nQuad+mnkl]
                //                          -buff[AxBzCzDz*nQuad+mnkl] + buff[AzBxCzDz*nQuad+mnkl];

                //auto dAcrossdB_zdCdotdD =  buff[AxByCxDx*nQuad+mnkl] - buff[AyBxCxDx*nQuad+mnkl]
                //                           + buff[AxByCyDy*nQuad+mnkl] - buff[AyBxCyDy*nQuad+mnkl]
                //                           + buff[AxByCzDz*nQuad+mnkl] - buff[AyBxCzDz*nQuad+mnkl];

                //// (∇A∙∇B)(∇Cx∇D)(mnkl)
                //auto dAdotdBdCcrossdD_x =  buff[AxBxCyDz*nQuad+mnkl] - buff[AxBxCzDy*nQuad+mnkl]
                //                           + buff[AyByCyDz*nQuad+mnkl] - buff[AyByCzDy*nQuad+mnkl]
                //                           + buff[AzBzCyDz*nQuad+mnkl] - buff[AzBzCzDy*nQuad+mnkl];

                //auto dAdotdBdCcrossdD_y = -buff[AxBxCxDz*nQuad+mnkl] + buff[AxBxCzDx*nQuad+mnkl]
                //                          -buff[AyByCxDz*nQuad+mnkl] + buff[AyByCzDx*nQuad+mnkl]
                //                          -buff[AzBzCxDz*nQuad+mnkl] + buff[AzBzCzDx*nQuad+mnkl];

                //auto dAdotdBdCcrossdD_z =  buff[AxBxCxDy*nQuad+mnkl] - buff[AxBxCyDx*nQuad+mnkl]
                //                           + buff[AyByCxDy*nQuad+mnkl] - buff[AyByCyDx*nQuad+mnkl]
                //                           + buff[AzBzCxDy*nQuad+mnkl] - buff[AzBzCyDx*nQuad+mnkl];

                //// (∇Ax∇B)(∇Cx∇D)(mnkl)
                //auto dAcrossdB_xdCcrossdD_x =  buff[AyBzCyDz*nQuad+mnkl] - buff[AzByCyDz*nQuad+mnkl]
                //                               - buff[AyBzCzDy*nQuad+mnkl] + buff[AzByCzDy*nQuad+mnkl];

                //auto dAcrossdB_xdCcrossdD_y =  buff[AyBzCzDx*nQuad+mnkl] - buff[AzByCzDx*nQuad+mnkl]
                //                               - buff[AyBzCxDz*nQuad+mnkl] + buff[AzByCxDz*nQuad+mnkl];

                //auto dAcrossdB_xdCcrossdD_z =  buff[AyBzCxDy*nQuad+mnkl] - buff[AzByCxDy*nQuad+mnkl]
                //                               - buff[AyBzCyDx*nQuad+mnkl] + buff[AzByCyDx*nQuad+mnkl];

                //auto dAcrossdB_ydCcrossdD_x =  buff[AzBxCyDz*nQuad+mnkl] - buff[AxBzCyDz*nQuad+mnkl]
                //                               - buff[AzBxCzDy*nQuad+mnkl] + buff[AxBzCzDy*nQuad+mnkl];

                //auto dAcrossdB_ydCcrossdD_y =  buff[AzBxCzDx*nQuad+mnkl] - buff[AxBzCzDx*nQuad+mnkl]
                //                               - buff[AzBxCxDz*nQuad+mnkl] + buff[AxBzCxDz*nQuad+mnkl];

                //auto dAcrossdB_ydCcrossdD_z =  buff[AzBxCxDy*nQuad+mnkl] - buff[AxBzCxDy*nQuad+mnkl]
                //                               - buff[AzBxCyDx*nQuad+mnkl] + buff[AxBzCyDx*nQuad+mnkl];

                //auto dAcrossdB_zdCcrossdD_x =  buff[AxByCyDz*nQuad+mnkl] - buff[AyBxCyDz*nQuad+mnkl]
                //                               - buff[AxByCzDy*nQuad+mnkl] + buff[AyBxCzDy*nQuad+mnkl];

                //auto dAcrossdB_zdCcrossdD_y =  buff[AxByCzDx*nQuad+mnkl] - buff[AyBxCzDx*nQuad+mnkl]
                //                               - buff[AxByCxDz*nQuad+mnkl] + buff[AyBxCxDz*nQuad+mnkl];

                //auto dAcrossdB_zdCcrossdD_z =  buff[AxByCxDy*nQuad+mnkl] - buff[AyBxCxDy*nQuad+mnkl]
                //                               - buff[AxByCyDx*nQuad+mnkl] + buff[AyBxCyDx*nQuad+mnkl];

                auto MNKL = m + n*NB + k*NB2 + l*NB3;

                // (mn|kl)
                ERIBuffABCD[         MNKL] =  dAdotdBdCdotdD;
                //ERIBuffABCD[   NB4 + MNKL] =  dAcrossdB_xdCdotdD;
                //ERIBuffABCD[ NB4_2 + MNKL] =  dAcrossdB_ydCdotdD;
                //ERIBuffABCD[ NB4_3 + MNKL] =  dAcrossdB_zdCdotdD;
                //ERIBuffABCD[ NB4_4 + MNKL] =  dAdotdBdCcrossdD_x;
                //ERIBuffABCD[ NB4_5 + MNKL] =  dAdotdBdCcrossdD_y;
                //ERIBuffABCD[ NB4_6 + MNKL] =  dAdotdBdCcrossdD_z;
                //ERIBuffABCD[ NB4_7 + MNKL] =  dAcrossdB_xdCcrossdD_x;
                //ERIBuffABCD[ NB4_8 + MNKL] =  dAcrossdB_xdCcrossdD_y;
                //ERIBuffABCD[ NB4_9 + MNKL] =  dAcrossdB_xdCcrossdD_z;
                //ERIBuffABCD[NB4_10 + MNKL] =  dAcrossdB_ydCcrossdD_x;
                //ERIBuffABCD[NB4_11 + MNKL] =  dAcrossdB_ydCcrossdD_y;
                //ERIBuffABCD[NB4_12 + MNKL] =  dAcrossdB_ydCcrossdD_z;
                //ERIBuffABCD[NB4_13 + MNKL] =  dAcrossdB_zdCcrossdD_x;
                //ERIBuffABCD[NB4_14 + MNKL] =  dAcrossdB_zdCcrossdD_y;
                //ERIBuffABCD[NB4_15 + MNKL] =  dAcrossdB_zdCcrossdD_z;

              } // ∇A∇B∇C∇D (SSSS) integrals

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tInts += tock(topInts);
            auto top1_2 = tick();
#endif
              //MatsT* symmDLLMS43_ptr = s43DenLLMSPtrs[s34];
              MatsT* symmDSSMS43_ptr = s43DenSSMSPtrs[s34];
              //MatsT* symmDSSMX43_ptr = s43DenSSMXPtrs[s34];
              //MatsT* symmDSSMY43_ptr = s43DenSSMYPtrs[s34];
              //MatsT* symmDSSMZ43_ptr = s43DenSSMZPtrs[s34];

              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              n3 = cint.shellSize(s3);
              n4 = cint.shellSize(s4);

              for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
                      symmDSSMS43_ptr+=nsh34) {

                const auto& [r, s] = rsPairs[iMat];

                // MO: Screening
                if (getMaxShBlkNorm(shBlkNormsSymmDenSS_rs[iMat], s1, s2, s3, s4)
                   * maxDiracCoulombSSSS < schwarzThreshold) {
                    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                    nConSkipped[thread_id] ++;
                    #endif
                    continue;
                   }

                auto& ADSS12 = s12SpinorSSSCRs[AD12_loc_off + iMat];
                auto& ADSSMS12 = ADSS12.S();
                //auto& ADSSMX12 = ADSS12.X();
                //auto& ADSSMY12 = ADSS12.Y();
                //auto& ADSSMZ12 = ADSS12.Z();

                for(auto m = 0ul, ms = 0ul; m < n1 and ms < n1; ++m, ++ms)
                for(auto n =   maxShellSize, ns = 0ul; n <maxShellSize + n2 and ns < n2; ++n, ++ns)
                for(auto k = 2*maxShellSize, ks = 0ul; k < 2*maxShellSize + n3 and ks < n3; ++k, ++ks)
                for(auto l = 3*maxShellSize, ls = 0ul; l < 3*maxShellSize + n4 and ls < n4; ++l, ++ls) {

                  auto MNKL = m + n*NB + k*NB2 + l*NB3;
                  auto KLMN = k + l*NB + m*NB2 + n*NB3;

                  auto MNKLdAdotdBdCdotdD          = ERIBuffABCD[         MNKL];
                  //auto MNKLdAcrossdB_xdCdotdD      = ERIBuffABCD[   NB4 + MNKL];
                  //auto MNKLdAcrossdB_ydCdotdD      = ERIBuffABCD[ NB4_2 + MNKL];
                  //auto MNKLdAcrossdB_zdCdotdD      = ERIBuffABCD[ NB4_3 + MNKL];
                  //auto MNKLdAdotdBdCcrossdD_x      = ERIBuffABCD[ NB4_4 + MNKL];
                  //auto MNKLdAdotdBdCcrossdD_y      = ERIBuffABCD[ NB4_5 + MNKL];
                  //auto MNKLdAdotdBdCcrossdD_z      = ERIBuffABCD[ NB4_6 + MNKL];
                  //auto MNKLdAcrossdB_xdCcrossdD_x  = ERIBuffABCD[ NB4_7 + MNKL];
                  //auto MNKLdAcrossdB_xdCcrossdD_y  = ERIBuffABCD[ NB4_8 + MNKL];
                  //auto MNKLdAcrossdB_xdCcrossdD_z  = ERIBuffABCD[ NB4_9 + MNKL];
                  //auto MNKLdAcrossdB_ydCcrossdD_x  = ERIBuffABCD[NB4_10 + MNKL];
                  //auto MNKLdAcrossdB_ydCcrossdD_y  = ERIBuffABCD[NB4_11 + MNKL];
                  //auto MNKLdAcrossdB_ydCcrossdD_z  = ERIBuffABCD[NB4_12 + MNKL];
                  //auto MNKLdAcrossdB_zdCcrossdD_x  = ERIBuffABCD[NB4_13 + MNKL];
                  //auto MNKLdAcrossdB_zdCcrossdD_y  = ERIBuffABCD[NB4_14 + MNKL];
                  //auto MNKLdAcrossdB_zdCcrossdD_z  = ERIBuffABCD[NB4_15 + MNKL];

                  auto ScaleS3S4 = s3==s4 ? 0.5 : 1.0;
                  auto ScaleF = s1==s2? ScaleS3S4*0.5: ScaleS3S4;

                  // First half transform.
                  /* Equation 70 in the paper */
                  ADSSMS12(ms, ns) += ScaleF * (symmDSSMS43_ptr[ls + ks * n4] * MNKLdAdotdBdCdotdD);
                                   // + (symmDSSMZ43_ptr[ls + ks * n4] * MNKLdAdotdBdCcrossdD_z
                                   // + symmDSSMX43_ptr[ls + ks * n4] * MNKLdAdotdBdCcrossdD_x
                                   // + symmDSSMY43_ptr[ls + ks * n4] * MNKLdAdotdBdCcrossdD_y) *dcomplex(0., 1.));


                //  /* Equation 71 in the paper */
                //  ADSSMZ12(ms, ns) +=
                //          ScaleF * (symmDSSMS43_ptr[ls + ks * n4] * MNKLdAcrossdB_zdCdotdD * dcomplex(0., 1.)
                //                    - symmDSSMZ43_ptr[ls + ks * n4] * MNKLdAcrossdB_zdCcrossdD_z
                //                    - symmDSSMX43_ptr[ls + ks * n4] * MNKLdAcrossdB_zdCcrossdD_x
                //                    - symmDSSMY43_ptr[ls + ks * n4] * MNKLdAcrossdB_zdCcrossdD_y);

                //  /* Equation 72 in the paper */
                //  ADSSMX12(ms, ns) +=
                //          ScaleF * (symmDSSMS43_ptr[ls + ks * n4] * MNKLdAcrossdB_xdCdotdD * dcomplex(0., 1.)
                //                    - symmDSSMZ43_ptr[ls + ks * n4] * MNKLdAcrossdB_xdCcrossdD_z
                //                    - symmDSSMX43_ptr[ls + ks * n4] * MNKLdAcrossdB_xdCcrossdD_x
                //                    - symmDSSMY43_ptr[ls + ks * n4] * MNKLdAcrossdB_xdCcrossdD_y);

                //  /* Equation 73 in the paper */
                //  ADSSMY12(ms, ns) +=
                //          ScaleF * (symmDSSMS43_ptr[ls + ks * n4] * MNKLdAcrossdB_ydCdotdD * dcomplex(0., 1.)
                //                    - symmDSSMZ43_ptr[ls + ks * n4] * MNKLdAcrossdB_ydCcrossdD_z
                //                    - symmDSSMX43_ptr[ls + ks * n4] * MNKLdAcrossdB_ydCcrossdD_x
                //                    - symmDSSMY43_ptr[ls + ks * n4] * MNKLdAcrossdB_ydCcrossdD_y);

                  } // mnkl
              } // iMat
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t1_2 += tock(top1_2);
#endif
            } // (s3, s4)
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
              const auto& [r, s] = rsPairs[iMat];

              auto& ADSS12 = s12SpinorSSSCRs[iMat];
              for (auto iThread = 1ul; iThread < nThreads; iThread++) {
                s12SpinorSSSCRs[iMat] += s12SpinorSSSCRs[iMat + iThread * rsPairs.size()];
              }

              // Second half trnasformation. 
              shBlockMO->transformSS(ADSS12, s1, s2, rsERI, off_sizes[0], off_sizes[1]);
              ADSS12.S().inplace_T();
              shBlockMO->transformSS(ADSS12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);


#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t2_2 += tock(top2_2);
            auto topUpdate = tick();
#endif
              // Scaling with 2 * 1/16c^4.
              MatsT C4 = 2./(16*SpeedOfLight()*SpeedOfLight()*SpeedOfLight()*SpeedOfLight());
              MatsT scale = C4;
              blas::axpy(npq, scale, rsERI.pointer(), 1, MOTPI + iMat * npq, 1);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tUpdate += tock(topUpdate);
#endif
            } // iMat
          } // end of parallel region
        } // (s1, s2)
    } // (s3, s4) Density Batching

    CQMemManager::get().free(buffERIAll, ERIBuffer, //buffDensityLLMS,
             buffDensitySSMS);

  // std::cout << "After DC + SSSS ERI Norm = " << std::setprecision(16) << 
  // lapack::lange(lapack::Norm::Fro, npq, nr * ns, MOTPI, npq) << std::endl;

    double totalTime = tock(startTime);
    FormattedLine(std::cout, "Time to DC(SS) Transform(s): ", totalTime);
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    auto printTimings = [] (const std::string& section,
        const std::vector<double> ts) {
        std::cout << std::setw(20) << section  << ": "
                  << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                  << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                  << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;
    };

    size_t nIntSkippedAcc = std::accumulate(nIntSkipped.begin(),nIntSkipped.end(),0);
    size_t nConSkippedAcc = std::accumulate(nConSkipped.begin(),nConSkipped.end(),0);

    #ifdef CQ_ENABLE_MPI
    MPIAllReduce(&nIntSkippedAcc, 1, &nIntSkippedAcc, comm_);
    MPIAllReduce(&nConSkippedAcc, 1, &nConSkippedAcc, comm_);
    #endif

    std::cout << "\nTiming for DC SSSS fully direct transformation: " << std::endl;
    std::cout << "Total time: " << totalTime << "s" << std::endl;
    std::cout << "Integrals skipped: " << nIntSkippedAcc << std::endl;
    std::cout << "Contractions skipped: " << nConSkippedAcc << std::endl;
    printTimings("t(Ints)", tInts_all);
    printTimings("t(Density)", tDensity_all);
    printTimings("t(1/2)", t1_2_all);
    printTimings("t(2/2)", t2_2_all);
    printTimings("t(Update)", tUpdate_all);
    std::cout << std::endl;
#endif

  } // SSSS
  /*************************************/
  /*                                   */
  /* End of Dirac-Coulomb C(4)-SSSS    */
  /*                                   */
  /*************************************/


  /*******************/
  /*                 */
  /* Start of Gaunt  */
  /*                 */
  /*******************/
  if (HOp.Gaunt) {
    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::cout << "  Transforming 2e-INTS: Gaunt ..." << std::endl;
    #endif
    auto startTime = tick();
    enum Gaunt_2ND_DERIV_AC {
        AxCx,
        AxCy,
        AxCz,
        AyCx,
        AyCy,
        AyCz,
        AzCx,
        AzCy,
        AzCz
    };

    enum Gaunt_2ND_DERIV_BC {
        BxCx,
        BxCy,
        BxCz,
        ByCx,
        ByCy,
        ByCz,
        BzCx,
        BzCy,
        BzCz
    };


    size_t NB  = maxShellSize*4;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;
    size_t NB4_2 = 2*NB4;
    size_t NB4_3 = 3*NB4;
    size_t NB4_4  = 4*NB4;
    size_t NB4_5  = 5*NB4;
    size_t NB4_6  = 6*NB4;
    size_t NB4_7  = 7*NB4;
    size_t NB4_8  = 8*NB4;
    size_t NB4_9  = 9*NB4;
    size_t NB4_10 =10*NB4;
    size_t NB4_11 =11*NB4;
    size_t NB4_12 =12*NB4;

    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::vector<double> tInts_all(nThreads, 0.);
    std::vector<double> t1_2_all(nThreads, 0.);
    std::vector<double> t2_2_all(nThreads, 0.);
    std::vector<double> tDensity_all(nThreads, 0.);
    std::vector<double> tUpdate_all(nThreads, 0.);

    std::vector<size_t> nConSkipped(nThreads, 0);
    std::vector<size_t> nIntSkipped(nThreads, 0);
    #endif
    // Allocate memory for raw ERI
    int nERI = 1;
    int nSave = 13;
    double *buffERIAll = CQMemManager::get().template malloc<double>(nERI*maxShellSize4*nThreads);
    // Allocate memory for assembled ERI, i.e., dot product and cross product
    double *ERIBuffer = CQMemManager::get().template malloc<double>(2*nSave*NB4*nThreads);

    // Get the sizes all density SCR
    size_t maxNDenSCR = 0ul;
    size_t n_s34 = 0ul;
    for (auto s3 = 0ul; s3 < nShell; s3++)
      for (auto s4 = 0ul; s4 < nShell; s4++, n_s34++) {
        maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
      } // (s3, s4) batches
    maxNDenSCR *= 4*rsPairs.size();

    size_t availableMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availableMem);

    // Allocate memory for density to be contracted with ERI
    MatsT *buffDensitySLMS =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMX =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMY =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMZ =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);

    // storage to save shell IDs of s3 and s4 for each thread
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;

    // storage to save s34 densities of all rs pairs for each thread
    std::vector<std::vector<MatsT*>> s43DenSLMSPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMXPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMYPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMZPtrsAll;

    s34PairsAll.resize(nThreads);
    s43DenSLMSPtrsAll.resize(nThreads);
    s43DenSLMXPtrsAll.resize(nThreads);
    s43DenSLMYPtrsAll.resize(nThreads);
    s43DenSLMZPtrsAll.resize(nThreads);

    // the spin-free Gaunt kernel is spin-diagonal: every Pauli component of the
    // LS-/+SL density contributes with the same integral, so all four are needed
    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorLSSCRs;
    for (auto i = 0ul; i < nThreads; ++i) {
      for (auto j = 0ul; j < rsPairs.size(); ++j) {
        s12SpinorLSSCRs.emplace_back(maxShellSize, true, true);
      }
    }
    
    for (auto i_s34 = 0ul; i_s34 < n_s34; ) {

      /**************************************/
      /*   Form AO Densities                */
      /**************************************/

      // FIXME:try to do better parallelism here
      // try best to evenly distribute the workloads across different threads
      std::vector<size_t> threadIdHeap(nThreads);
      std::vector<size_t> threadLoads(nThreads, 0ul);
      std::iota(threadIdHeap.begin(), threadIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) {
          return threadLoads[i] > threadLoads[j];
      };

      for (auto iThread = 0ul; iThread < nThreads; iThread++) {
        s34PairsAll[iThread].clear();
        s43DenSLMSPtrsAll[iThread].clear();
        s43DenSLMXPtrsAll[iThread].clear();
        s43DenSLMYPtrsAll[iThread].clear();
        s43DenSLMZPtrsAll[iThread].clear();
      }

      size_t nMem = 0ul, nMemOff = 0ul;
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell and nMem < nDenSCR; s3++)
        for (auto s4 = 0ul; s4 < nShell; s4++, s34++) {
          if (s34 < i_s34) continue;
          size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
          nMem += nDen34SCR*4;
          if (nMem > nDenSCR) break;

          // assigning to a thread
          size_t iThread = threadIdHeap[0];
          threadLoads[iThread] += nDen34SCR;
          std::make_heap(threadIdHeap.begin(), threadIdHeap.end(), comp);

          s34PairsAll[iThread].push_back({s3, s4});
          s43DenSLMSPtrsAll[iThread].push_back(buffDensitySLMS + nMemOff);
          s43DenSLMXPtrsAll[iThread].push_back(buffDensitySLMX + nMemOff);
          s43DenSLMYPtrsAll[iThread].push_back(buffDensitySLMY + nMemOff);
          s43DenSLMZPtrsAll[iThread].push_back(buffDensitySLMZ + nMemOff);

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
        const auto& s43DenSLMSPtrs = s43DenSLMSPtrsAll[thread_id];
        const auto& s43DenSLMXPtrs = s43DenSLMXPtrsAll[thread_id];
        const auto& s43DenSLMYPtrs = s43DenSLMYPtrsAll[thread_id];
        const auto& s43DenSLMZPtrs = s43DenSLMZPtrsAll[thread_id];
        auto& denSLSCR = s12SpinorLSSCRs[thread_id];

        #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        auto& tDensity = tDensity_all[thread_id];
        auto topDensity = tick();

        #endif

        for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
          MatsT* denPtrSLMS = s43DenSLMSPtrs[s34];
          MatsT* denPtrSLMX = s43DenSLMXPtrs[s34];
          MatsT* denPtrSLMY = s43DenSLMYPtrs[s34];
          MatsT* denPtrSLMZ = s43DenSLMZPtrs[s34];
          const auto& [s3, s4] = s34Pairs[s34];
          size_t nsh34 = cint.shellSize(s3) * cint.shellSize(s4);
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
               denPtrSLMS+=nsh34, denPtrSLMX+=nsh34, denPtrSLMY+=nsh34,
               denPtrSLMZ+=nsh34) {
            const auto& [r, s] = rsPairs[iMat];

            shBlockMO->genDenLSpmDenSL(r + roff, s + soff, s3, s4, denSLSCR);

            // S: C_s,s3xC_r,s4^* - C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.S().pointer(), nsh34, denPtrSLMS);
            // X: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.X().pointer(), nsh34, denPtrSLMX);
            // Y: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.Y().pointer(), nsh34, denPtrSLMY);
            // Z: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.Z().pointer(), nsh34, denPtrSLMZ);

          } // iMat
        } // s34Pairs

        #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        tDensity += tock(topDensity);
        #endif

      } // end of the parallel region

      for (auto s1 = 0ul, s12 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 < nShell; s2++, s12++) {

          #ifdef CQ_ENABLE_MPI
          if (s12Assignment[s12] != MPIRank(comm_)) continue;
          #endif

          /**************************************/
          /*  First Half Transformation         */
          /**************************************/
          #pragma omp parallel
          {
            size_t n1 = cint.shellSize(s1);
            size_t n2 = cint.shellSize(s2);

            auto topFirstHalf = tick();

            int thread_id = GetThreadID();
            const auto& s34Pairs = s34PairsAll[thread_id];
            const auto& s43DenSLMSPtrs = s43DenSLMSPtrsAll[thread_id];
            const auto& s43DenSLMXPtrs = s43DenSLMXPtrsAll[thread_id];
            const auto& s43DenSLMYPtrs = s43DenSLMYPtrsAll[thread_id];
            const auto& s43DenSLMZPtrs = s43DenSLMZPtrsAll[thread_id];

            double *buff = buffERIAll + nERI*maxShellSize4 * thread_id;
            double *ERIBuffBC   = &ERIBuffer[thread_id*nSave*NB4];
            double *ERIBuffBD   = &ERIBuffer[nSave*NB4*nThreads + thread_id*nSave*NB4];

            int shls[4];
            shls[0] = int(s1);
            shls[1] = int(s2);

            // initialize storage for (𝜇𝜈|rs) where 𝜇∈s1 and 𝜈∈s2
            size_t AD12_loc_off = thread_id * rsPairs.size();
            for (auto iMat = AD12_loc_off; iMat < AD12_loc_off + rsPairs.size(); iMat++) {
              s12SpinorLSSCRs[iMat].resize(n1, n2);
              s12SpinorLSSCRs[iMat].clear();
            }

            #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto& tInts = tInts_all[thread_id];
            auto& t1_2 = t1_2_all[thread_id];
            #endif

            for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
              const auto& [s3, s4] = s34Pairs[s34];

              auto maxGaunt = schwarzInts.maxGaunt(s1, s2, s3, s4);

              if (getMaxShBlkNorm(maxShBlkNormsSymmDenSL_rs, s1, s2, s3, s4) *
                 maxGaunt < schwarzThreshold) {
                  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                  nIntSkipped[thread_id] ++;
                  #endif
                  continue;
                 }

              // MO: index this way to get the ∇B∇C
              shls[0] = int(s2);
              shls[1] = int(s1);
              shls[2] = int(s3);
              shls[3] = int(s4);

              #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              auto topInts = tick();
              #endif

              // MO: This is done by switching m and n so that the cint gives the right intgrals
              // ∇B∇C(mn|kl)
              if (cint.compute_int2e_gaunt_ps1ps2_sph(buff, shls) == 0) continue;
              
              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              size_t n3 = cint.shellSize(s3);
              size_t n4 = cint.shellSize(s4);

              size_t nsh34 = n3 * n4;
              auto nQuad = n1 * n2 * n3 * n4;

              for(auto l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
              for(auto k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k)
              for(auto m = 0                                ; m <                  n1; ++m)
              for(auto n =   maxShellSize            ; n <   maxShellSize + n2; ++n, ++mnkl) {

                auto dAdotdC = -buff[mnkl];
                // ∇Ax∇C(mn|kl)
                //auto dAcrossdC_x =  buff[AyCz*nQuad+mnkl] - buff[AzCy*nQuad+mnkl];
                //auto dAcrossdC_y = -buff[AxCz*nQuad+mnkl] + buff[AzCx*nQuad+mnkl];
                //auto dAcrossdC_z =  buff[AxCy*nQuad+mnkl] - buff[AyCx*nQuad+mnkl];
                // Change the index so that we do ∇B∙∇C(ij|kl) using the ∇B∇C engine
                auto MNKL = m + n*NB + k*NB2 + l*NB3;

                // (kl|mn)

                // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
                ERIBuffBC[      MNKL] = dAdotdC;
                //ERIBuffBC[  NB4+MNKL] = dAcrossdC_x;
                //ERIBuffBC[NB4_2+MNKL] = dAcrossdC_y;
                //ERIBuffBC[NB4_3+MNKL] = dAcrossdC_z;

                // ∇B_x∇C_x(mn|kl) - ∇B∙∇C(mn|kl)
                //ERIBuffBC[NB4_4+MNKL] = buff[BxCx*nQuad+mnkl] - dAdotdC;

                // ∇B_y∇C_x(mn|kl)
                //ERIBuffBC[NB4_5+MNKL] = buff[ByCx*nQuad+mnkl];

                // ∇B_z∇C_x(mn|kl)
                //ERIBuffBC[NB4_6+MNKL] = buff[BzCx*nQuad+mnkl];

                // ∇B_x∇C_y(mn|kl)
                //ERIBuffBC[NB4_7+MNKL] = buff[BxCy*nQuad+mnkl];

                // ∇B_y∇C_y(mn|kl) - ∇B∙∇C(mn|kl)
                //ERIBuffBC[NB4_8+MNKL] = buff[ByCy*nQuad+mnkl] - dAdotdC;

                // ∇B_z∇C_y(mn|kl)
                //ERIBuffBC[NB4_9+MNKL] = buff[BzCy*nQuad+mnkl];

                // ∇B_x∇C_z(mn|kl)
                //ERIBuffBC[NB4_10+MNKL] = buff[BxCz*nQuad+mnkl];

                // ∇B_y∇C_z(mn|kl)
                //ERIBuffBC[NB4_11+MNKL] = buff[ByCz*nQuad+mnkl];

                // ∇B_z∇C_z(mn|kl) - ∇B∙∇C(mn|kl)
                //ERIBuffBC[NB4_12+MNKL] = buff[BzCz*nQuad+mnkl] - dAdotdC;

              } // ∇B∇C integral preparation loop

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tInts += tock(topInts);
            auto top1_2 = tick();
#endif
              MatsT* symmDSLMS43_ptr = s43DenSLMSPtrs[s34];
              MatsT* symmDSLMX43_ptr = s43DenSLMXPtrs[s34];
              MatsT* symmDSLMY43_ptr = s43DenSLMYPtrs[s34];
              MatsT* symmDSLMZ43_ptr = s43DenSLMZPtrs[s34];

              for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
                      symmDSLMS43_ptr+=nsh34, symmDSLMX43_ptr+=nsh34,
                      symmDSLMY43_ptr+=nsh34, symmDSLMZ43_ptr+=nsh34) {                

                const auto& [r, s] = rsPairs[iMat];

                // MO: Screening
                if (getMaxShBlkNorm(shBlkNormsSymmDenSL_rs[iMat], s1, s2, s3, s4)
                    * maxGaunt < schwarzThreshold) {
                    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                    nConSkipped[thread_id] ++;
                    #endif
                    continue;
                    }

                auto& ADLS12 = s12SpinorLSSCRs[AD12_loc_off + iMat];
                auto& ADLSMS12 = ADLS12.S();
                auto& ADLSMX12 = ADLS12.X();
                auto& ADLSMY12 = ADLS12.Y();
                auto& ADLSMZ12 = ADLS12.Z();

                for(auto m = 0ul, ms = 0ul; m < n1 and ms < n1; ++m, ++ms)
                for(auto n =   maxShellSize, ns = 0ul; n <   maxShellSize + n2 and ns < n2; ++n, ++ns)
                for(auto k = 2*maxShellSize, ks = 0ul; k < 2*maxShellSize + n3 and ks < n3; ++k, ++ks)
                for(auto l = 3*maxShellSize, ls = 0ul; l < 3*maxShellSize + n4 and ls < n4; ++l, ++ls) {

                  auto MNKL = m + n * NB + k * NB2 + l * NB3;

                  //ERI0 :   ∇B∙∇C(mn|kl)
                  //ERI1 :   ∇Bx∇C(mn|kl)-X
                  //ERI2 :   ∇Bx∇C(mn|kl)-Y
                  //ERI3 :   ∇Bx∇C(mn|kl)-Z
                  //ERI4 :   ∇B_x∇C_x(mn|kl) - ∇B∙∇C(mn|kl)
                  //ERI5 :   ∇B_y∇C_x(mn|kl)
                  //ERI6 :   ∇B_z∇C_x(mn|kl)
                  //ERI7 :   ∇B_x∇C_y(mn|kl)
                  //ERI8 :   ∇B_y∇C_y(mn|kl) - ∇B∙∇C(mn|kl)
                  //ERI9 :   ∇B_z∇C_y(mn|kl)
                  //ERI10:   ∇B_x∇C_z(mn|kl)
                  //ERI11:   ∇B_y∇C_z(mn|kl)
                  //ERI12:   ∇B_z∇C_z(mn|kl) - ∇B∙∇C(mn|kl)

                  /*++++++++++++++++++++++++*/
                  /* Start of Gaunt (LS|LS) */
                  /*++++++++++++++++++++++++*/

                  size_t bf43 = ls + ks * n4;

                  // Spin-free reduction of the spinor expressions above. The
                  // spin-dependent engine contracts
                  //   MS <- -D_S (∇B∙∇C) + i D_c (∇Bx∇C)_c
                  //   Mi <- i D_S (∇Bx∇C)_i + D_j (∇B_j∇C_i - δ_ij ∇B∙∇C)
                  // Spin averaging kills (∇Bx∇C) and replaces the individual
                  // ∇B_j∇C_i by their isotropic part (1/3) δ_ij (∇B∙∇C), which is
                  // all the spin-free engine (int2e_gaunt_ps1ps2) provides. Hence
                  // the scalar component keeps the full kernel while each Pauli
                  // component keeps (1/3 - 1) = -2/3 of it.
                  //
                  // NOTE: the AO Fock engines (direct4C_libcint_spinfree.hpp and
                  // direct4C_libcint_coulombonly_spinfree.hpp) use -1 for all four
                  // components; that is what makes the old SSFock CASCI path
                  // disagree with its own SCF energy for spin-free Gaunt.
                  const double fcSFVec = 2. / 3.;

                  // 91 + 136
                  ADLSMS12(ms, ns) += -symmDSLMS43_ptr[bf43] * ERIBuffBC[MNKL];

                  // 92 x + 137 x
                  ADLSMX12(ms, ns) += -fcSFVec * symmDSLMX43_ptr[bf43] * ERIBuffBC[MNKL];

                  // 92 y + 137 y
                  ADLSMY12(ms, ns) += -fcSFVec * symmDSLMY43_ptr[bf43] * ERIBuffBC[MNKL];

                  // 92 z +137 z
                  ADLSMZ12(ms, ns) += -fcSFVec * symmDSLMZ43_ptr[bf43] * ERIBuffBC[MNKL];

                } // mnkl
              } // iMat
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t1_2 += tock(top1_2);
#endif

            } // (s3, s4)
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

              for (auto iThread = 1ul; iThread < nThreads; iThread++) {
                s12SpinorLSSCRs[iMat] += s12SpinorLSSCRs[iMat + iThread * rsPairs.size()];
              }

              auto& ADLS12 = s12SpinorLSSCRs[iMat];
              shBlockMO->transformLS(ADLS12, s1, s2, rsERI, off_sizes[0], off_sizes[1]);
              // (SL) half from the (LS) half by symmetry:
              // CSLMS = -[CLSMS]^T, CSLM{X,Y,Z} = [CLSM{X,Y,Z}]^T
              ADLS12.S().inplace_scaleT(MatsT(-1.), 'T');
              ADLS12.X().inplace_scaleT(MatsT(1.), 'T');
              ADLS12.Y().inplace_scaleT(MatsT(1.), 'T');
              ADLS12.Z().inplace_scaleT(MatsT(1.), 'T');
              shBlockMO->transformSL(ADLS12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t2_2 += tock(top2_2);
            auto topUpdate = tick();
#endif

              MatsT C2 = 2./(4*SpeedOfLight()*SpeedOfLight());
              // if gauge is requested scale the gaunt term by 1/2
              MatsT scale = (HOp.Gauge) ? C2 / 2. : C2;
              
              blas::axpy(npq, scale, rsERI.pointer(), 1, MOTPI + iMat * npq, 1);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tUpdate += tock(topUpdate);
#endif
            } // iMat
          } // end of parallel region
        } // (s1, s2)
      } // (s3, s4) Density Batching

    CQMemManager::get().free(buffERIAll, ERIBuffer, buffDensitySLMS,
             buffDensitySLMX, buffDensitySLMY, buffDensitySLMZ);

  //  std::cout << "After Gaunt ERI Norm = " << std::setprecision(16) << 
  //   lapack::lange(lapack::Norm::Fro, npq, nr * ns, MOTPI, npq) << std::endl;

    double totalTime = tock(startTime);
    FormattedLine(std::cout, "Time to Gaunt Transform(s): ", totalTime);
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    auto printTimings = [] (const std::string& section,
        const std::vector<double> ts) {
        std::cout << std::setw(20) << section  << ": "
                  << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                  << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                  << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;
    };

    size_t nIntSkippedAcc = std::accumulate(nIntSkipped.begin(),nIntSkipped.end(),0);
    size_t nConSkippedAcc = std::accumulate(nConSkipped.begin(),nConSkipped.end(),0);

    #ifdef CQ_ENABLE_MPI
    MPIAllReduce(&nIntSkippedAcc, 1, &nIntSkippedAcc, comm_);
    MPIAllReduce(&nConSkippedAcc, 1, &nConSkippedAcc, comm_);
    #endif

    std::cout << "\nTiming for Gaunt fully direct transformation: " << std::endl;
    std::cout << "Total time: " << totalTime << "s" << std::endl;
    std::cout << "Integrals skipped: " << nIntSkippedAcc << std::endl;
    std::cout << "Contractions skipped: " << nConSkippedAcc << std::endl;
    printTimings("t(Ints)", tInts_all);
    printTimings("t(Density)", tDensity_all);
    printTimings("t(1/2)", t1_2_all);
    printTimings("t(2/2)", t2_2_all);
    printTimings("t(Update)", tUpdate_all);
    std::cout << std::endl;
#endif

  } // Gaunt
  /*******************/
  /*                 */
  /*   End of Gaunt  */
  /*                 */
  /*******************/


  /*******************/
  /*                 */
  /* Start of Gauge  */
  /*                 */
  /*******************/
  if (HOp.Gauge) {
    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::cout << "  Transforming 2e-INTS: Gauge ..." << std::endl;
    #endif
    auto startTime = tick();
    enum Gauge_ERI {
        SxSx, // σ_x * σ_x     0
        SySx, // σ_y * σ_x     1
        SzSx, // σ_z * σ_x     2
        ISx,  // I   * σ_x     3
        SxSy, // σ_x * σ_y     4
        SySy, // σ_y * σ_y     5
        SzSy, // σ_z * σ_y     6
        ISy,  // I   * σ_y     7
        SxSz, // σ_x * σ_z     8
        SySz, // σ_y * σ_z     9
        SzSz, // σ_z * σ_z     10
        ISz,  // I   * σ_z     11
        SxI,  // σ_x * I       12
        SyI,  // σ_y * I       13
        SzI,  // σ_z * I       14
        II    // I   * I       15
      };


    size_t NB  = maxShellSize*4;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;
    size_t NB4_2 = 2*NB4;
    size_t NB4_3 = 3*NB4;
    size_t NB4_4  = 4*NB4;  
    size_t NB4_5  = 5*NB4;  
    size_t NB4_6  = 6*NB4;  
    size_t NB4_7  = 7*NB4;  
    size_t NB4_8  = 8*NB4;  
    size_t NB4_9  = 9*NB4;  
    size_t NB4_10 =10*NB4;  
    size_t NB4_11 =11*NB4;  
    size_t NB4_12 =12*NB4;  
    size_t NB4_13 =13*NB4;  
    size_t NB4_14 =14*NB4;  
    size_t NB4_15 =15*NB4;  

    size_t ss = 0     ;    // ERI00 (ss)(ij|kl)
    size_t sx = NB4   ;    // ERI01 (sσ)_x(ijkl)
    size_t sy = NB4_2 ;    // ERI02 (sσ)_y(ijkl)
    size_t sz = NB4_3 ;    // ERI03 (sσ)_z(ijkl)
    size_t xs = NB4_4 ;    // ERI04 (σs)_x(ijkl)
    size_t ys = NB4_5 ;    // ERI05 (σs)_y(ijkl)
    size_t zs = NB4_6 ;    // ERI06 (σs)_z(ijkl)
    size_t xx = NB4_7;     // ERI07 (σ_x σ_x)(ijkl)
    size_t xy = NB4_8;     // ERI08 (σ_x σ_y)(ijkl)
    size_t xz = NB4_9;     // ERI09 (σ_x σ_z)(ijkl)
    size_t yx = NB4_10;    // ERI10 (σ_y σ_x)(ijkl)
    size_t yy = NB4_11;    // ERI11 (σ_y σ_y)(ijkl)
    size_t yz = NB4_12;    // ERI12 (σ_y σ_z)(ijkl)
    size_t zx = NB4_13;    // ERI13 (σ_z σ_x)(ijkl)
    size_t zy = NB4_14;    // ERI14 (σ_z σ_y)(ijkl)
    size_t zz = NB4_15;    // ERI15 (σ_z σ_z)(ijkl)

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    std::vector<double> tInts_all(nThreads, 0.);
    std::vector<double> t1_2_all(nThreads, 0.);
    std::vector<double> t2_2_all(nThreads, 0.);
    std::vector<double> tDensity_all(nThreads, 0.);
    std::vector<double> tUpdate_all(nThreads, 0.);

    std::vector<size_t> nConSkipped(nThreads, 0);
    std::vector<size_t> nIntSkipped(nThreads, 0);
#endif

    // Allocate memory for raw ERI
    // int2e_gauge_r{1,2}_sp1ps2 have ncomp = 4 (3 spin-dependent + 1 spin-free);
    // nERI must match or libcint writes past each thread's slab of buffERIAll
    int nERI = 4;
    int nSave = 16;
    double *buffERIAll = CQMemManager::get().template malloc<double>(2*nERI*maxShellSize4*nThreads);
    // Allocate memory for assembled ERI, i.e., dot product and cross product
    double *ERIBuffer = CQMemManager::get().template malloc<double>(2*nSave*NB4*nThreads);

    // Get the sizes all density SCR
    size_t maxNDenSCR = 0ul;
    size_t n_s34 = 0ul;
    for (auto s3 = 0ul; s3 < nShell; s3++)
      for (auto s4 = 0ul; s4 < nShell; s4++, n_s34++) {
        maxNDenSCR += cint.shellSize(s3) * cint.shellSize(s4);
      } // (s3, s4) batches
    maxNDenSCR *= 4*rsPairs.size();

    size_t availableMem = CQMemManager::get().max_avail_allocatable<MatsT>(1, maxNDenSCR);
    size_t nDenSCR = std::min(maxNDenSCR, availableMem);

    // Allocate memory for density to be contracted with ERI
    MatsT *buffDensitySLMS =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMX =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMY =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);
    MatsT *buffDensitySLMZ =  CQMemManager::get().template malloc<MatsT>(nDenSCR/4);

    // storage to save shell IDs of s3 and s4 for each thread
    std::vector<std::vector<std::pair<size_t, size_t>>> s34PairsAll;

    // storage to save s34 densities of all rs pairs for each thread
    std::vector<std::vector<MatsT*>> s43DenSLMSPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMXPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMYPtrsAll;
    std::vector<std::vector<MatsT*>> s43DenSLMZPtrsAll;

    s34PairsAll.resize(nThreads);
    s43DenSLMSPtrsAll.resize(nThreads);
    s43DenSLMXPtrsAll.resize(nThreads);
    s43DenSLMYPtrsAll.resize(nThreads);
    s43DenSLMZPtrsAll.resize(nThreads);

    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> s12SpinorLSSCRs;
    for (auto i = 0ul; i < nThreads; ++i) {
      for (auto j = 0ul; j < rsPairs.size(); ++j) {
        s12SpinorLSSCRs.emplace_back(maxShellSize, true, true);
      }
    }
    
    for (auto i_s34 = 0ul; i_s34 < n_s34; ) {

      /**************************************/
      /*   Form AO Densities                */
      /**************************************/

      // FIXME:try to do better parallelism here
      // try best to evenly distribute the workloads across different threads
      std::vector<size_t> threadIdHeap(nThreads);
      std::vector<size_t> threadLoads(nThreads, 0ul);
      std::iota(threadIdHeap.begin(), threadIdHeap.end(), 0ul);
      auto comp = [&] (size_t i, size_t j) {
          return threadLoads[i] > threadLoads[j];
      };

      for (auto iThread = 0ul; iThread < nThreads; iThread++) {
        s34PairsAll[iThread].clear();
        s43DenSLMXPtrsAll[iThread].clear();
        s43DenSLMYPtrsAll[iThread].clear();
        s43DenSLMZPtrsAll[iThread].clear();
        s43DenSLMSPtrsAll[iThread].clear();
      }

      size_t nMem = 0ul, nMemOff = 0ul;
      for (auto s3 = 0ul, s34 = 0ul; s3 < nShell and nMem < nDenSCR; s3++)
        for (auto s4 = 0ul; s4 < nShell; s4++, s34++) {
          if (s34 < i_s34) continue;
          size_t nDen34SCR = cint.shellSize(s3) * cint.shellSize(s4) * rsPairs.size();
          nMem += nDen34SCR*4;
          if (nMem > nDenSCR) break;

          // assigning to a thread
          size_t iThread = threadIdHeap[0];
          threadLoads[iThread] += nDen34SCR;
          std::make_heap(threadIdHeap.begin(), threadIdHeap.end(), comp);

          s34PairsAll[iThread].push_back({s3, s4});
          s43DenSLMSPtrsAll[iThread].push_back(buffDensitySLMS + nMemOff);
          s43DenSLMXPtrsAll[iThread].push_back(buffDensitySLMX + nMemOff);
          s43DenSLMYPtrsAll[iThread].push_back(buffDensitySLMY + nMemOff);
          s43DenSLMZPtrsAll[iThread].push_back(buffDensitySLMZ + nMemOff);

          nMemOff += nDen34SCR;

          i_s34++;
        } // s34 assignment

#pragma omp parallel
      {
        int thread_id = GetThreadID();
        const auto& s34Pairs = s34PairsAll[thread_id];
        const auto& s43DenSLMSPtrs = s43DenSLMSPtrsAll[thread_id];
        const auto& s43DenSLMXPtrs = s43DenSLMXPtrsAll[thread_id];
        const auto& s43DenSLMYPtrs = s43DenSLMYPtrsAll[thread_id];
        const auto& s43DenSLMZPtrs = s43DenSLMZPtrsAll[thread_id];
        auto& denSLSCR = s12SpinorLSSCRs[thread_id];

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        auto& tDensity = tDensity_all[thread_id];
        auto topDensity = tick();
#endif

        for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
          MatsT* denPtrSLMS = s43DenSLMSPtrs[s34];
          MatsT* denPtrSLMX = s43DenSLMXPtrs[s34];
          MatsT* denPtrSLMY = s43DenSLMYPtrs[s34];
          MatsT* denPtrSLMZ = s43DenSLMZPtrs[s34];
          const auto& [s3, s4] = s34Pairs[s34];
          size_t nsh34 = cint.shellSize(s3) * cint.shellSize(s4);
          for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
               denPtrSLMS+=nsh34, denPtrSLMX+=nsh34, denPtrSLMY+=nsh34,
               denPtrSLMZ+=nsh34) {
            const auto& [r, s] = rsPairs[iMat];

            shBlockMO->genDenLSpmDenSL(r + roff, s + soff, s3, s4, denSLSCR);

            // S: C_s,s3xC_r,s4^* - C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.S().pointer(), nsh34, denPtrSLMS);
            // X: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.X().pointer(), nsh34, denPtrSLMX);
            // Y: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.Y().pointer(), nsh34, denPtrSLMY);
            // Z: C_s,s3xC_r,s4^* + C_r,s3^*xC_s,s4, small component, assuming real integrals
            std::copy_n(denSLSCR.Z().pointer(), nsh34, denPtrSLMZ);

          } // iMat
        } // s34Pairs

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
        tDensity += tock(topDensity);
#endif

      } // end of the parallel region

      for (auto s1 = 0ul, s12 = 0ul; s1 < nShell; s1++)
        for (auto s2 = 0ul; s2 < nShell; s2++, s12++) {

#ifdef CQ_ENABLE_MPI
          if (s12Assignment[s12] != MPIRank(comm_)) continue;
#endif

          /**************************************/
          /*  First Half Transformation         */
          /**************************************/
#pragma omp parallel
          {
            size_t n1 = cint.shellSize(s1);
            size_t n2 = cint.shellSize(s2);

            auto topFirstHalf = tick();

            int thread_id = GetThreadID();
            const auto& s34Pairs = s34PairsAll[thread_id];
            const auto& s43DenSLMSPtrs = s43DenSLMSPtrsAll[thread_id];
            const auto& s43DenSLMXPtrs = s43DenSLMXPtrsAll[thread_id];
            const auto& s43DenSLMYPtrs = s43DenSLMYPtrsAll[thread_id];
            const auto& s43DenSLMZPtrs = s43DenSLMZPtrsAll[thread_id];
            double *buff1 = buffERIAll + nERI * maxShellSize4 * thread_id;
            double *buff2 = buffERIAll + nERI * maxShellSize4 * nThreads + nERI * maxShellSize4 * thread_id;

            
            double *ERIBuffBC   = &ERIBuffer[thread_id*nSave*NB4];
            double *ERIBuffBD   = &ERIBuffer[nSave*NB4*nThreads + thread_id*nSave*NB4];

            int shls[4];
            shls[0] = int(s1);
            shls[1] = int(s2);

            // initialize storage for (𝜇𝜈|rs) where 𝜇∈s1 and 𝜈∈s2
            size_t AD12_loc_off = thread_id * rsPairs.size();
            for (auto iMat = AD12_loc_off; iMat < AD12_loc_off + rsPairs.size(); iMat++) {
              s12SpinorLSSCRs[iMat].resize(n1, n2);
              s12SpinorLSSCRs[iMat].clear();
            }

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
            auto& tInts = tInts_all[thread_id];
            auto& t1_2 = t1_2_all[thread_id];
#endif

            for (auto s34 = 0ul; s34 < s34Pairs.size(); s34++) {
              const auto& [s3, s4] = s34Pairs[s34];

              auto maxGauge = schwarzInts.maxGauge(s1, s2, s3, s4);

              if (getMaxShBlkNorm(maxShBlkNormsSymmDenSL_rs, s1, s2, s3, s4) *
                 maxGauge < schwarzThreshold) {
                  #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                  nIntSkipped[thread_id] ++;
                  #endif
                  continue;
                 }

              shls[0] = int(s1);
              shls[1] = int(s2);
              shls[2] = int(s3);
              shls[3] = int(s4);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              auto topInts = tick();
#endif

              //∇B∇C
              auto skip1 = cint.compute_int2e_gauge_r1_sp1ps2_sph(buff1, shls);
              auto skip2 = cint.compute_int2e_gauge_r2_sp1ps2_sph(buff2, shls);
              if (skip1 == 0 and skip2 == 0) continue;

              n1 = cint.shellSize(s1);
              n2 = cint.shellSize(s2);
              size_t n3 = cint.shellSize(s3);
              size_t n4 = cint.shellSize(s4);

              size_t nsh34 = n3 * n4;
              auto nQuad = n1 * n2 * n3 * n4;

              for(auto l = 3*maxShellSize, mnkl = 0ul; l < 3*maxShellSize + n4; ++l)
              for(auto k = 2*maxShellSize            ; k < 2*maxShellSize + n3; ++k)
              for(auto n =   maxShellSize            ; n <   maxShellSize + n2; ++n)
              for(auto m = 0                                ; m <                  n1; ++m, ++mnkl) {

                auto MNKL = m + n*NB + k*NB2 + l*NB3;
                //enum Gauge_ERI {
                //  SxSx,  // σ_x * σ_x     0
                //  SySx,  // σ_y * σ_x     1
                //  SzSx,  // σ_z * σ_x     2
                //  ISx,   // I   * σ_x     3
                //  SxSy,  // σ_x * σ_y     4
                //  SySy,  // σ_y * σ_y     5
                //  SzSy,  // σ_z * σ_y     6
                //  ISy,   // I   * σ_y     7
                //  SxSz,  // σ_x * σ_z     8
                //  SySz,  // σ_y * σ_z     9
                //  SzSz,  // σ_z * σ_z     10
                //  ISz,   // I   * σ_z     11
                //  SxI,   // σ_x * I       12
                //  SyI,   // σ_y * I       13
                //  SzI,   // σ_z * I       14
                //  II     // I   * I       15
                //};

                // auto ss = 0     ;    // ERI00 (ss)(ij|kl)
                // auto sx = NB4   ;    // ERI01 (sσ)_x(ijkl)
                // auto sy = NB4_2 ;    // ERI02 (sσ)_y(ijkl)
                // auto sz = NB4_3 ;    // ERI03 (sσ)_z(ijkl)
                // auto xs = NB4_4 ;    // ERI04 (σs)_x(ijkl)
                // auto ys = NB4_5 ;    // ERI05 (σs)_y(ijkl)
                // auto zs = NB4_6 ;    // ERI06 (σs)_z(ijkl)
                // auto xx = NB4_7;     // ERI07 (σ_x σ_x)(ijkl)
                // auto xy = NB4_8;     // ERI08 (σ_x σ_y)(ijkl)
                // auto xz = NB4_9;     // ERI09 (σ_x σ_z)(ijkl)
                // auto yx = NB4_10;    // ERI10 (σ_y σ_x)(ijkl)
                // auto yy = NB4_11;    // ERI11 (σ_y σ_y)(ijkl)
                // auto yz = NB4_12;    // ERI12 (σ_y σ_z)(ijkl)
                // auto zx = NB4_13;    // ERI13 (σ_z σ_x)(ijkl)
                // auto zy = NB4_14;    // ERI14 (σ_z σ_y)(ijkl)
                // auto zz = NB4_15;    // ERI15 (σ_z σ_z)(ijkl)

                // MNKL

                // (ss): the spin-free (I x I) piece is component 3 of the 4-component
                // sp1ps2 buffer, not component 15 of the 16-component spinor layout
                ERIBuffBC[MNKL] = buff1[3*nQuad+mnkl] - buff2[3*nQuad+mnkl];

                // (sigma . sigma): components 0, 1, 2 are the spin-free
                // (sigma_x sigma_x), (sigma_y sigma_y) and (sigma_z sigma_z) pieces
                ERIBuffBC[NB4+MNKL] = -(buff1[mnkl]           - buff2[mnkl])
                                      -(buff1[nQuad+mnkl]     - buff2[nQuad+mnkl])
                                      -(buff1[2*nQuad+mnkl]   - buff2[2*nQuad+mnkl]);
                
                // (sσ)_x
                // ERIBuffBC[NB4+MNKL] = buff1[3*nQuad+mnkl] - buff2[3*nQuad+mnkl];
                
                // (sσ)_y
                // ERIBuffBC[NB4_2+MNKL] = buff1[7*nQuad+mnkl] - buff2[7*nQuad+mnkl];
                
                // (sσ)_z
                // ERIBuffBC[NB4_3+MNKL] = buff1[11*nQuad+mnkl] - buff2[11*nQuad+mnkl];
                
                // (σs)_x
                // ERIBuffBC[NB4_4+MNKL] = buff1[12*nQuad+mnkl] - buff2[12*nQuad+mnkl];
                
                // (σs)_y
                // ERIBuffBC[NB4_5+MNKL] = buff1[13*nQuad+mnkl] - buff2[13*nQuad+mnkl];
                
                // (σs)_z
                // ERIBuffBC[NB4_6+MNKL] = buff1[14*nQuad+mnkl] - buff2[14*nQuad+mnkl];
                
                // σ_x*σ_x
                // ERIBuffBC[NB4_7+MNKL] = -(buff1[mnkl] - buff2[mnkl]);
                
                // σ_x*σ_y
                // ERIBuffBC[NB4_8+MNKL] = -(buff1[4*nQuad+mnkl] - buff2[4*nQuad+mnkl]);
                
                // σ_x*σ_z
                // ERIBuffBC[NB4_9+MNKL] = -(buff1[8*nQuad+mnkl] - buff2[8*nQuad+mnkl]);
                
                // σ_y*σ_x
                // ERIBuffBC[NB4_10+MNKL] = -(buff1[1*nQuad+mnkl] - buff2[1*nQuad+mnkl]);
                
                // σ_y*σ_y
                // ERIBuffBC[NB4_11+MNKL] = -(buff1[5*nQuad+mnkl] - buff2[5*nQuad+mnkl]);
                
                // σ_y*σ_z
                // ERIBuffBC[NB4_12+MNKL] = -(buff1[9*nQuad+mnkl] - buff2[9*nQuad+mnkl]);
                
                // σ_z*σ_x
                // ERIBuffBC[NB4_13+MNKL] = -(buff1[2*nQuad+mnkl] - buff2[2*nQuad+mnkl]);
                
                // σ_z*σ_y
                // ERIBuffBC[NB4_14+MNKL] = -(buff1[6*nQuad+mnkl] - buff2[6*nQuad+mnkl]);
                
                // σ_z*σ_z
                //ERIBuffBC[NB4_15+MNKL] = -(buff1[10*nQuad+mnkl] - buff2[10*nQuad+mnkl]);

              } // ∇B∇C integral preparation loop

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tInts += tock(topInts);
            auto top1_2 = tick();
#endif
              MatsT* symmDSLMS43_ptr = s43DenSLMSPtrs[s34];
              MatsT* symmDSLMX43_ptr = s43DenSLMXPtrs[s34];
              MatsT* symmDSLMY43_ptr = s43DenSLMYPtrs[s34];
              MatsT* symmDSLMZ43_ptr = s43DenSLMZPtrs[s34];

              for (auto iMat = 0ul; iMat < rsPairs.size(); iMat++,
                      symmDSLMS43_ptr+=nsh34, symmDSLMX43_ptr+=nsh34,
                      symmDSLMY43_ptr+=nsh34, symmDSLMZ43_ptr+=nsh34) {         
                               
                const auto& [r, s] = rsPairs[iMat];

                // MO: Screening
                if (getMaxShBlkNorm(shBlkNormsSymmDenSL_rs[iMat], s1, s2, s3, s4)
                    * maxGauge < schwarzThreshold) {
                    #ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
                    nConSkipped[thread_id] ++;
                    #endif
                    continue;
                    }

                auto& ADLS12 = s12SpinorLSSCRs[AD12_loc_off + iMat];
                auto& ADLSMS12 = ADLS12.S();
                auto& ADLSMX12 = ADLS12.X();
                auto& ADLSMY12 = ADLS12.Y();
                auto& ADLSMZ12 = ADLS12.Z();
                
                for(auto m = 0ul, ms = 0ul; m < n1 and ms < n1; ++m, ++ms)
                for(auto n =   maxShellSize, ns = 0ul; n <   maxShellSize + n2 and ns < n2; ++n, ++ns)
                for(auto k = 2*maxShellSize, ks = 0ul; k < 2*maxShellSize + n3 and ks < n3; ++k, ++ks)
                for(auto l = 3*maxShellSize, ls = 0ul; l < 3*maxShellSize + n4 and ls < n4; ++l, ++ls) {

                  auto MNKL = m + n * NB + k * NB2 + l * NB3;
                  auto KLMN = k + l * NB + m * NB2 + n * NB3;
                  auto MNLK = m + n * NB + l * NB2 + k * NB3;
                  auto LKNM = l + k * NB + n * NB2 + m * NB3;

                  /*++++++++++++++++++++++++*/
                  /* Start of Gauge (LL|SS) */
                  /*++++++++++++++++++++++++*/

                  size_t bf43 = ls + ks * n4;

                  // Spin-free reduction of the spinor expressions: the single-sigma
                  // (sigma s) / (s sigma) kernels are spin-dependent and vanish, and
                  // sigma_i sigma_j keeps only its isotropic part,
                  // (1/3) delta_ij (sigma . sigma). The scalar component therefore
                  // contracts with the (ss) kernel and each Pauli component with
                  // one third of the (sigma . sigma) kernel.
                  const double fcSFVec = 1. / 3.;

                  // 232 S + 233 S
                  ADLSMS12(ms, ns) += -( symmDSLMS43_ptr[bf43] * ERIBuffBC[MNKL] );
                  // 232 x + 233 x
                  ADLSMX12(ms, ns) += -fcSFVec * symmDSLMX43_ptr[bf43] * ERIBuffBC[MNKL + NB4];
                  // 232 y + 233 y
                  ADLSMY12(ms, ns) += -fcSFVec * symmDSLMY43_ptr[bf43] * ERIBuffBC[MNKL + NB4];
                  // 232 z + 233 z
                  ADLSMZ12(ms, ns) += -fcSFVec * symmDSLMZ43_ptr[bf43] * ERIBuffBC[MNKL + NB4];

                } // mnkl
              } // iMat
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t1_2 += tock(top1_2);
#endif

            } // (s3, s4)
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

              for (auto iThread = 1ul; iThread < nThreads; iThread++) {
                s12SpinorLSSCRs[iMat] += s12SpinorLSSCRs[iMat + iThread * rsPairs.size()];
              }

              auto& ADLS12 = s12SpinorLSSCRs[iMat];
              shBlockMO->transformLS(ADLS12, s1, s2, rsERI, off_sizes[0], off_sizes[1]);
              // same as Gaunt: CSLMS = -[CLSMS]^T, CSLM{X,Y,Z} = [CLSM{X,Y,Z}]^T
              ADLS12.S().inplace_scaleT(MatsT(-1.), 'T');
              ADLS12.X().inplace_scaleT(MatsT(1.), 'T');
              ADLS12.Y().inplace_scaleT(MatsT(1.), 'T');
              ADLS12.Z().inplace_scaleT(MatsT(1.), 'T');
              shBlockMO->transformSL(ADLS12, s2, s1, rsERI, off_sizes[0], off_sizes[1], true);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              t2_2 += tock(top2_2);
            auto topUpdate = tick();
#endif

              MatsT C2 = 2./(4*SpeedOfLight()*SpeedOfLight());
              MatsT scale = C2 / 2.;
              
              blas::axpy(npq, scale, rsERI.pointer(), 1, MOTPI + iMat * npq, 1);

#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
              tUpdate += tock(topUpdate);
#endif
            } // iMat
          } // end of parallel region
        } // (s1, s2)
      } // (s3, s4) Density Batching

    CQMemManager::get().free(buffERIAll, ERIBuffer, buffDensitySLMS,
             buffDensitySLMX, buffDensitySLMY, buffDensitySLMZ);

  //  std::cout << "After Gauge ERI Norm = " << std::setprecision(16) << 
  //   lapack::lange(lapack::Norm::Fro, npq, nr * ns, MOTPI, npq) << std::endl;

    double totalTime = tock(startTime);
    FormattedLine(std::cout, "Time to Gauge Transform(s): ", totalTime);
#ifdef _MOINTSTRANSFORMER_TPI_FULL_DIRECT_TIMING
    auto printTimings = [] (const std::string& section,
        const std::vector<double> ts) {
        std::cout << std::setw(20) << section  << ": "
                  << "average time = " << std::accumulate(ts.begin(), ts.end(), double(0.)) / ts.size() << " s"
                  << ", max time = " << *std::max_element(ts.begin(), ts.end()) << " s"
                  << ", mim time = " << *std::min_element(ts.begin(), ts.end()) << " s" << std::endl;
    };

    size_t nIntSkippedAcc = std::accumulate(nIntSkipped.begin(),nIntSkipped.end(),0);
    size_t nConSkippedAcc = std::accumulate(nConSkipped.begin(),nConSkipped.end(),0);

    #ifdef CQ_ENABLE_MPI
    MPIAllReduce(&nIntSkippedAcc, 1, &nIntSkippedAcc, comm_);
    MPIAllReduce(&nConSkippedAcc, 1, &nConSkippedAcc, comm_);
    #endif



    std::cout << "\nTiming for Gauge fully direct transformation: " << std::endl;
    std::cout << "Total time: " << totalTime << "s" << std::endl;
    std::cout << "Integrals skipped: " << nIntSkippedAcc << std::endl;
    std::cout << "Contractions skipped: " << nConSkippedAcc << std::endl;
    printTimings("t(Ints)", tInts_all);
    printTimings("t(Density)", tDensity_all);
    printTimings("t(1/2)", t1_2_all);
    printTimings("t(2/2)", t2_2_all);
    printTimings("t(Update)", tUpdate_all);
    std::cout << std::endl;
#endif
  } // Gauge
  /*******************/
  /*                 */
  /*   End of Gauge  */
  /*                 */
  /*******************/
  // } // Check if MatsT is dcomplex at compile time 

  if (nC != 4) {
    HOp.BareCoulomb = false;
  }
  // else {
  //   HOp.Gaunt = false;
  //   HOp.Gauge = false;
  // }

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

#if 1
  // restore to no symmetry

  if (pqSymm and rsSymm) {
    // if pqSymm and rsSymm:  r <= s upper triangle
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
#endif

  // check Matrix Norm of ERI
  // std::cout << "Total ERI Norm = " << std::setprecision(16) << 
  //   lapack::lange(lapack::Norm::Fro, npq, nr * ns, MOTPI, npq) << std::endl;
  
  return;
} // directTransformScalarTPIBatch

} // namespace ChronusQ
