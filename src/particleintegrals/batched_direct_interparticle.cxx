/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */

#include <particleintegrals/contract/batched_direct_interparticle.hpp>

#include <cerr.hpp>
#include <util/math.hpp>
#include <util/threads.hpp>

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
#include <iomanip>
#include <iostream>
#include <util/timer.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <iterator>
#include <limits>
#include <utility>

namespace ChronusQ {

  namespace {

    // Named tuning constants (kept out of the hot loops on purpose).
    constexpr size_t derivativeBufferCount = 12;       // 3 xyz x 4 shells
    constexpr size_t profileCounterCount = 3;          // candidate/evaluated/nonzero
    constexpr size_t parallelFoldThreshold = 1 << 20;  // serial vs OpenMP fold
    constexpr double machineEpsilon = std::numeric_limits<double>::epsilon();

    // Significant shell pair, flattened for the hot loops.
    struct ShellPairInfo {
      size_t firstShell;
      size_t secondShell;
      size_t firstBasisFunction;
      size_t secondBasisFunction;
      size_t firstShellSize;
      size_t secondShellSize;
    };

    // BasisSet pointer, with its triangular shell-pair workingData list.
    struct BasisInfo {
      BasisSet* basis = nullptr;
      std::vector<ShellPairInfo> significantShellPairs;
    };

    // Source density (to be contracted out) and its conservative shell-block norms.
    template <typename MatsT>
    struct DensityInfo {
      size_t basisInfoIndex = 0;
      const MatsT* density = nullptr;
      std::vector<double> shellBlockInfinityNorms;
    };

    // Coulomb output and its offsets in the packed buffers.
    template <typename MatsT>
    struct OutputInfo {
      size_t basisInfoIndex = 0;
      MatsT* output = nullptr;
      size_t rankOffset = 0;
      size_t threadOffset = 0;
      bool receivesInnerSideWrites = false;
    };

    // Interaction pair, oriented so the larger shell-pair list is the OUTER basis and the other is the INNER basis.
    template <typename MatsT>
    struct OrientedPair {
      size_t outerBasisInfoIndex = 0;
      size_t innerBasisInfoIndex = 0;
      const double* outerSchwarzBounds = nullptr;
      const double* innerSchwarzBounds = nullptr;

      const MatsT* outerDensity = nullptr;
      const MatsT* innerDensity = nullptr;
      size_t outerDensityInfoIndex = 0;
      size_t innerDensityInfoIndex = 0;

      bool formOuterOutput = false;
      bool formInnerOutput = false;
      size_t outerOutputIndex = 0;
      size_t innerOutputIndex = 0;

      double interactionScale = 1.;
      double schwarzThreshold = 0.;
    };

    // Interactions sharing an outer basis are visited by the same outer
    // shell-pair tasks (one task per outer shell pair).
    struct OuterBasisGroup {
      size_t outerBasisInfoIndex = 0;
      std::vector<size_t> interactionIndices;
    };

    struct ShellPairTask {
      size_t outerBasisGroupIndex;
      size_t outerShellPairIndex;
    };

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
    struct DirectProfile {
      size_t activeInteractions = 0;
      size_t workUnits = 0;
      size_t candidateQuartets = 0;
      size_t evaluatedQuartets = 0;
      size_t nonzeroQuartets = 0;
      double preparationSeconds = 0.;
      double engineSetupSeconds = 0.;
      double contractionSeconds = 0.;
      double threadFoldSeconds = 0.;
      double mpiReductionSeconds = 0.;
      double outputSeconds = 0.;
      double totalSeconds = 0.;
    };
#endif

    // --------------------------------------------------------------
    // Shared helpers (used by both the Coulomb and the gradient kernel)
    // --------------------------------------------------------------

    std::vector<ShellPairInfo> collectSignificantShellPairs(
      const BasisSet& basis) {

      std::vector<ShellPairInfo> shellPairs;

      for(size_t firstShell = 0; firstShell < basis.nShell; ++firstShell) {
        const auto significantPartners = basis.shellData.sigShellPair.find(firstShell);
        if(significantPartners == basis.shellData.sigShellPair.end())
          CErr("Missing significant shell-pair data in batched direct interparticle contraction.");

        for(const size_t secondShell : significantPartners->second) {
          shellPairs.push_back({
            firstShell,
            secondShell,
            basis.mapSh2Bf[firstShell],
            basis.mapSh2Bf[secondShell],
            basis.shells[firstShell].size(),
            basis.shells[secondShell].size()
          });
        }
      }

      return shellPairs;
    }

    template <typename MatsT>
    std::vector<double> formHermitianDensityShellBlockInfinityNorms(
      const BasisSet& basis, const MatsT* density) {

      if(density == nullptr)
        CErr("Null density in batched direct interparticle screening.");

      std::vector<double> shellBlockNorms(basis.nShell * basis.nShell, 0.);

      for(size_t firstShell = 0; firstShell < basis.nShell; ++firstShell) {
        const size_t firstBasisFunction = basis.mapSh2Bf[firstShell];
        const size_t firstShellSize = basis.shells[firstShell].size();

        for(size_t secondShell = 0; secondShell < basis.nShell; ++secondShell) {
          const size_t secondBasisFunction = basis.mapSh2Bf[secondShell];
          const size_t secondShellSize = basis.shells[secondShell].size();
          double infinityNorm = 0.;

          for(size_t firstFunction = 0; firstFunction < firstShellSize; ++firstFunction) {
            double rowSum = 0.;
            for(size_t secondFunction = 0; secondFunction < secondShellSize; ++secondFunction) {
              const size_t matrixIndex = firstBasisFunction + firstFunction + (secondBasisFunction + secondFunction) * basis.nBasis;
              const double densityMagnitude = std::abs(density[matrixIndex]);
              rowSum += densityMagnitude;
            }
            infinityNorm = std::max(infinityNorm, rowSum);
          }

          shellBlockNorms[firstShell + secondShell * basis.nShell] = infinityNorm;
        }
      }

      // The digest reads the transpose block. Use the larger row-sum norm in
      // either orientation so screening stays conservative.
      for(size_t firstShell = 0; firstShell < basis.nShell; ++firstShell)
      for(size_t secondShell = 0; secondShell < firstShell; ++secondShell) {
        const size_t firstIndex = firstShell + secondShell * basis.nShell;
        const size_t secondIndex = secondShell + firstShell * basis.nShell;
        const double blockNorm = std::max(shellBlockNorms[firstIndex], shellBlockNorms[secondIndex]);
        shellBlockNorms[firstIndex] = blockNorm;
        shellBlockNorms[secondIndex] = blockNorm;
      }

      return shellBlockNorms;
    }

    template <typename MatsT>
    void accumulateHermitianOutput(
      const MatsT* packedRawOutput, MatsT* output, size_t basisFunctionCount) {

      // The shell-pair digest stores one triangular contribution. Restore the
      // Hermitian Coulomb matrix after the rank reduction.
      for(size_t column = 0; column < basisFunctionCount; ++column)
      for(size_t row = 0; row < basisFunctionCount; ++row) {
        const size_t element = row + column * basisFunctionCount;
        const size_t transposeElement = column + row * basisFunctionCount;
        output[element] += MatsT(0.5) * (packedRawOutput[element] + SmartConj(packedRawOutput[transposeElement]));
      }
    }

    size_t registerBasis(std::vector<BasisInfo>& basisInfo, BasisSet& basis) {
      for(size_t index = 0; index < basisInfo.size(); ++index)
        if(basisInfo[index].basis == &basis) return index;

      const size_t index = basisInfo.size();
      basisInfo.push_back({&basis, collectSignificantShellPairs(basis)});
      return index;
    }

    template <typename MatsT>
    size_t registerDensity(std::vector<DensityInfo<MatsT>>& densityInfo,
      const std::vector<BasisInfo>& basisInfo, const size_t basisIndex,
      const MatsT* density, const bool computeShellBlockNorms = true) {

      if(density == nullptr)
        CErr("Null source density in batched direct interparticle contraction.");

      for(size_t index = 0; index < densityInfo.size(); ++index)
        if(densityInfo[index].basisInfoIndex == basisIndex and
           densityInfo[index].density == density)
          return index;

      const size_t index = densityInfo.size();
      densityInfo.push_back({basisIndex, density,
        computeShellBlockNorms ? formHermitianDensityShellBlockInfinityNorms(*basisInfo[basisIndex].basis, density) : std::vector<double>()});
      return index;
    }

    // Register a requested output (Coulomb or gradient record) and deduplicate
    // by pointer. receivesInnerSideWrites marks outputs that inner-side digests
    // can write, which need thread-local accumulation.
    template <typename OutputRecordT>
    size_t registerOutput(std::vector<OutputRecordT>& outputInfo,
      const size_t basisIndex, auto output, const bool receivesInnerSideWrites) {

      if(output == nullptr)
        CErr("Null requested output in batched direct interparticle contraction.");

      for(size_t index = 0; index < outputInfo.size(); ++index) {
        auto& registered = outputInfo[index];
        if(registered.output != output) continue;
        if(registered.basisInfoIndex != basisIndex)
          CErr("One batched direct interparticle output has inconsistent bases.");
        registered.receivesInnerSideWrites |= receivesInnerSideWrites;
        return index;
      }

      const size_t index = outputInfo.size();
      outputInfo.push_back({basisIndex, output, 0, 0, receivesInnerSideWrites});
      return index;
    }

    // Orient one pair: the subsystem with the larger shell-pair list becomes
    // the OUTER basis, the other the INNER basis. Also assigns the densities,
    // the requested-output flags, and the Schwarz bounds (nullptr when the
    // caller does not screen, e.g. the gradient path).
    template <typename MatsT>
    void orientPair(OrientedPair<MatsT>& pair,
      const bool subsystem1IsOuter,
      const size_t subsystem1BasisIndex, const size_t subsystem2BasisIndex,
      const MatsT* firstDensity, const MatsT* secondDensity,
      const bool formFirstOutput, const bool formSecondOutput,
      const double* subsystem1Schwarz, const double* subsystem2Schwarz) {

      pair.outerBasisInfoIndex = subsystem1IsOuter ? subsystem1BasisIndex : subsystem2BasisIndex;
      pair.innerBasisInfoIndex = subsystem1IsOuter ? subsystem2BasisIndex : subsystem1BasisIndex;
      pair.outerSchwarzBounds = subsystem1IsOuter ? subsystem1Schwarz : subsystem2Schwarz;
      pair.innerSchwarzBounds = subsystem1IsOuter ? subsystem2Schwarz : subsystem1Schwarz;
      pair.outerDensity = subsystem1IsOuter ? firstDensity : secondDensity;
      pair.innerDensity = subsystem1IsOuter ? secondDensity : firstDensity;
      pair.formOuterOutput = subsystem1IsOuter ? formFirstOutput : formSecondOutput;
      pair.formInnerOutput = subsystem1IsOuter ? formSecondOutput : formFirstOutput;
    }

    template <typename OrientedPairT>
    std::vector<OuterBasisGroup> buildOuterBasisGroups(
      const std::vector<OrientedPairT>& orientedPairs) {

      std::vector<OuterBasisGroup> groups;
      for(size_t interactionIndex = 0; interactionIndex < orientedPairs.size(); ++interactionIndex) {
        const size_t outerBasisInfoIndex = orientedPairs[interactionIndex].outerBasisInfoIndex;
        auto matchingGroup = std::find_if(
          groups.begin(), groups.end(),
          [&](const auto& group) {
            return group.outerBasisInfoIndex == outerBasisInfoIndex;
          });
        if(matchingGroup == groups.end()) {
          groups.push_back({outerBasisInfoIndex, {}});
          matchingGroup = std::prev(groups.end());
        }
        matchingGroup->interactionIndices.push_back(interactionIndex);
      }
      return groups;
    }

    std::vector<ShellPairTask> buildShellPairTasks(
      const std::vector<BasisInfo>& basisInfo,
      const std::vector<OuterBasisGroup>& outerBasisGroups) {

      std::vector<ShellPairTask> shellPairTasks;
      for(size_t groupIndex = 0; groupIndex < outerBasisGroups.size(); ++groupIndex) {
        const auto& outerShellPairs = basisInfo[outerBasisGroups[groupIndex].outerBasisInfoIndex].significantShellPairs;
        for(size_t outerPairIndex = 0; outerPairIndex < outerShellPairs.size(); ++outerPairIndex)
          shellPairTasks.push_back({groupIndex, outerPairIndex});
      }
      return shellPairTasks;
    }

    template <typename MatsT>
    void foldThreadLocalOutputs(
      const std::vector<MatsT>& threadOutputs, const size_t threadCount,
      const size_t threadLocalInnerElementCount,
      const std::vector<size_t>& threadLocalToRankMap,
      std::vector<MatsT>& rankOutput) {

      auto foldElement = [&](size_t threadLocalIndex) {
        MatsT value = MatsT(0.);
        for(size_t thread = 0; thread < threadCount; ++thread)
          value += threadOutputs[thread * threadLocalInnerElementCount + threadLocalIndex];
        rankOutput[threadLocalToRankMap[threadLocalIndex]] += value;
      };

      // Heuristic: only launch a second OpenMP region when the thread-local
      // working data is large enough to matter; the result is identical
      // either way.
      if(threadLocalInnerElementCount * threadCount >= parallelFoldThreshold) {
        #pragma omp parallel for schedule(static)
        for(long long threadLocalElement = 0; threadLocalElement < static_cast<long long>(threadLocalInnerElementCount); ++threadLocalElement)
          foldElement(static_cast<size_t>(threadLocalElement));
      } else {
        for(size_t threadLocalElement = 0; threadLocalElement < threadLocalInnerElementCount; ++threadLocalElement)
          foldElement(threadLocalElement);
      }
    }

    // --------------------------------------------------------------
    // Coulomb kernel: the prepared batch, its workingData, and the contract
    // --------------------------------------------------------------

    template <typename MatsT>
    struct CoulombBatch {
      std::vector<BasisInfo> basisInfo;
      std::vector<DensityInfo<MatsT>> densityInfo;
      std::vector<OutputInfo<MatsT>> outputInfo;
      std::vector<OrientedPair<MatsT>> orientedPairs;
      size_t maximumPrimitiveCount = 0;
      size_t maximumAngularMomentum = 0;
    };

    template <typename MatsT, typename IntsT>
    CoulombBatch<MatsT> prepareCoulombBatch(
      const std::vector<DirectInterparticleJContraction<MatsT,IntsT>>& directInterparticleJContractions) {

      CoulombBatch<MatsT> batch;
      auto& basisInfo = batch.basisInfo;
      auto& densityInfo = batch.densityInfo;
      auto& outputInfo = batch.outputInfo;
      auto& orientedPairs = batch.orientedPairs;
      orientedPairs.reserve(directInterparticleJContractions.size());

      for(const auto& jContraction : directInterparticleJContractions) {
        if(not jContraction.formFirstCoulomb and not jContraction.formSecondCoulomb)
          CErr("Batched direct interparticle interaction requests no Coulomb output.");
        if(not jContraction.integrals)
          CErr("Null DirectTPI in batched direct interparticle contraction.");
        if(jContraction.integrals->kernel() != TPI_KERNEL::Coulomb)
          CErr("Batched direct interparticle contraction supports the Coulomb TPI kernel only.");
        BasisSet& integralBasis1 = jContraction.integrals->basisSet();
        BasisSet& integralBasis2 = jContraction.integrals->basisSet2();
        if(integralBasis1.basisType != REAL_GTO or integralBasis2.basisType != REAL_GTO)
          CErr("Batched direct interparticle contraction supports real GTO bases only.");

        const bool integralBasesAreIdentical = &integralBasis1 == &integralBasis2;
        if(jContraction.integrals->schwarz() == nullptr or (not integralBasesAreIdentical and jContraction.integrals->schwarz2() == nullptr))
          jContraction.integrals->computeSchwarz();

        const double* integralBasis1Schwarz = jContraction.integrals->schwarz();
        const double* integralBasis2Schwarz = integralBasesAreIdentical ?
          integralBasis1Schwarz : jContraction.integrals->schwarz2();

        // Map the two subsystems onto the two integral bases. If the flag is
        // true, subsystem 1 is on integral basis 1 (and subsystem 2 on basis
        // 2); if false, the mapping is swapped.
        BasisSet& subsystem1Basis = jContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis1 : integralBasis2;
        BasisSet& subsystem2Basis = jContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis2 : integralBasis1;
        const double* subsystem1Schwarz = jContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis1Schwarz : integralBasis2Schwarz;
        const double* subsystem2Schwarz = jContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis2Schwarz : integralBasis1Schwarz;

        const size_t subsystem1BasisIndex = registerBasis(basisInfo, subsystem1Basis);
        const size_t subsystem2BasisIndex = registerBasis(basisInfo, subsystem2Basis);
        // Putting the larger shell-pair list outside gives E--P batches a
        // shared electronic workingData axis and leaves the small proton pairs inside.
        const bool subsystem1IsOuter = basisInfo[subsystem1BasisIndex].significantShellPairs.size() >= basisInfo[subsystem2BasisIndex].significantShellPairs.size();

        OrientedPair<MatsT> pair;
        pair.interactionScale = jContraction.interactionScale;
        pair.schwarzThreshold = jContraction.integrals->threshSchwarz();

        orientPair(pair, subsystem1IsOuter, subsystem1BasisIndex, subsystem2BasisIndex,
          jContraction.firstDensity, jContraction.secondDensity, jContraction.formFirstCoulomb, jContraction.formSecondCoulomb,
          subsystem1Schwarz, subsystem2Schwarz);

        if(pair.formOuterOutput) {
          pair.innerDensityInfoIndex = registerDensity(densityInfo, basisInfo, pair.innerBasisInfoIndex, pair.innerDensity);
          pair.outerOutputIndex = registerOutput(outputInfo, pair.outerBasisInfoIndex,
            subsystem1IsOuter ? jContraction.firstCoulomb : jContraction.secondCoulomb, false);
        }
        if(pair.formInnerOutput) {
          pair.outerDensityInfoIndex = registerDensity(densityInfo, basisInfo, pair.outerBasisInfoIndex, pair.outerDensity);
          pair.innerOutputIndex = registerOutput(outputInfo, pair.innerBasisInfoIndex,
            subsystem1IsOuter ? jContraction.secondCoulomb : jContraction.firstCoulomb, true);
        }

        batch.maximumPrimitiveCount = std::max(batch.maximumPrimitiveCount, std::max(subsystem1Basis.maxPrim, subsystem2Basis.maxPrim));
        batch.maximumAngularMomentum = std::max(batch.maximumAngularMomentum, std::max(subsystem1Basis.maxL, subsystem2Basis.maxL));
        orientedPairs.push_back(std::move(pair));
      }

      return batch;
    }

    template <typename MatsT>
    struct CoulombWorkingData {
      size_t rankBufferElementCount = 0;
      size_t threadLocalInnerElementCount = 0;
      std::vector<size_t> threadLocalToRankMap;
      std::vector<OuterBasisGroup> outerBasisGroups;
      std::vector<ShellPairTask> shellPairTasks;
    };

    template <typename MatsT>
    CoulombWorkingData<MatsT> packCoulombWorkingData(
      CoulombBatch<MatsT>& batch) {

      CoulombWorkingData<MatsT> workingData;
      auto& rankBufferElementCount = workingData.rankBufferElementCount;
      auto& threadLocalInnerElementCount = workingData.threadLocalInnerElementCount;

      // Assign every output its slice of the packed rank buffer, and the
      // inner-written outputs their slice of the thread-local buffer.
      for(auto& output : batch.outputInfo) {
        const size_t basisFunctionCount = batch.basisInfo[output.basisInfoIndex].basis->nBasis;
        if(basisFunctionCount != 0 and
           basisFunctionCount > std::numeric_limits<size_t>::max() / basisFunctionCount)
          CErr("Batched direct interparticle output size overflow.");

        const size_t matrixElementCount = basisFunctionCount * basisFunctionCount;
        if(rankBufferElementCount > std::numeric_limits<size_t>::max() - matrixElementCount)
          CErr("Batched direct interparticle packed output size overflow.");
        output.rankOffset = workingData.rankBufferElementCount;
        workingData.rankBufferElementCount += matrixElementCount;

        if(not output.receivesInnerSideWrites) continue;
        if(threadLocalInnerElementCount >
           std::numeric_limits<size_t>::max() - matrixElementCount)
          CErr("Batched direct interparticle thread output size overflow.");
        output.threadOffset = workingData.threadLocalInnerElementCount;
        workingData.threadLocalInnerElementCount += matrixElementCount;
      }

      // Map every thread-local element to its rank-buffer location.
      workingData.threadLocalToRankMap.assign(threadLocalInnerElementCount, 0);
      for(const auto& output : batch.outputInfo) {
        if(not output.receivesInnerSideWrites) continue;
        const size_t basisFunctionCount = batch.basisInfo[output.basisInfoIndex].basis->nBasis;
        const size_t matrixElementCount = basisFunctionCount * basisFunctionCount;
        for(size_t element = 0; element < matrixElementCount; ++element)
          workingData.threadLocalToRankMap[output.threadOffset + element] = output.rankOffset + element;
      }

      workingData.outerBasisGroups = buildOuterBasisGroups(batch.orientedPairs);
      workingData.shellPairTasks = buildShellPairTasks(batch.basisInfo, workingData.outerBasisGroups);
      return workingData;
    }

    template <typename MatsT>
    void contractCoulomb(MPI_Comm comm,
      const CoulombBatch<MatsT>& batch,
      const CoulombWorkingData<MatsT>& workingData
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      , DirectProfile& profile
#endif
      ) {

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      profile.activeInteractions = batch.orientedPairs.size();
      profile.workUnits = workingData.shellPairTasks.size();
#endif

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      // The profile appends a few counters to the packed buffer so they ride
      // the same collective.
      const size_t profileOffset = workingData.rankBufferElementCount;
      const size_t reductionElementCount = workingData.rankBufferElementCount + profileCounterCount;
#else
      const size_t reductionElementCount = workingData.rankBufferElementCount;
#endif
      if(reductionElementCount > static_cast<size_t>(std::numeric_limits<int>::max()))
        CErr("Batched direct interparticle reduction exceeds the MPI count limit.");

      const size_t threadCount = GetNumThreads();
      if(threadCount == 0)
        CErr("Batched direct interparticle contraction received zero OpenMP threads.");
      if(workingData.threadLocalInnerElementCount != 0 and threadCount > std::numeric_limits<size_t>::max() / workingData.threadLocalInnerElementCount)
        CErr("Batched direct interparticle thread-local output size overflow.");

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      const auto engineSetupStart = tick();
#endif
      const size_t laThreadCount = GetLAThreads();
      SetLAThreads(1);
      std::vector<libint2::Engine> engines(threadCount);
      engines[0] = libint2::Engine(libint2::Operator::coulomb, batch.maximumPrimitiveCount, batch.maximumAngularMomentum, 0);
      engines[0].set_precision(machineEpsilon);
      for(size_t thread = 1; thread < threadCount; ++thread)
        engines[thread] = engines[0];

      std::vector<MatsT> rankOutput(reductionElementCount, MatsT(0.));
      std::vector<MatsT> threadOutputs(threadCount * workingData.threadLocalInnerElementCount, MatsT(0.));
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      std::vector<size_t> threadCandidateQuartets(threadCount, 0);
      std::vector<size_t> threadEvaluatedQuartets(threadCount, 0);
      std::vector<size_t> threadNonzeroQuartets(threadCount, 0);
      profile.engineSetupSeconds = tock(engineSetupStart);
#endif

      const size_t mpiRank = MPIRank(comm);
      const size_t mpiSize = MPISize(comm);
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      const auto contractionStart = tick();
#endif

      // Phase 3: each task owns one outer shell pair. It evaluates compatible
      // inner shell pairs and digests both Coulomb directions from one ERI.
      #pragma omp parallel for schedule(dynamic, 1)
      for(long long shellPairTaskIndex = static_cast<long long>(mpiRank); shellPairTaskIndex < static_cast<long long>(workingData.shellPairTasks.size()); shellPairTaskIndex += static_cast<long long>(mpiSize)) {

        const size_t threadIndex = GetThreadID();
        auto& engine = engines[threadIndex];
        MatsT* threadLocalOutput = workingData.threadLocalInnerElementCount == 0 ? nullptr : threadOutputs.data() + threadIndex * workingData.threadLocalInnerElementCount;

        const auto& task = workingData.shellPairTasks[static_cast<size_t>(shellPairTaskIndex)];
        const auto& basisGroup = workingData.outerBasisGroups[task.outerBasisGroupIndex];
        const auto& outerBasisData = batch.basisInfo[basisGroup.outerBasisInfoIndex];
        const BasisSet& outerBasis = *outerBasisData.basis;
        const auto& outerShellPair = outerBasisData.significantShellPairs[task.outerShellPairIndex];

        const size_t outerFirstShell = outerShellPair.firstShell;
        const size_t outerSecondShell = outerShellPair.secondShell;
        const double outerPairDegeneracy = outerFirstShell == outerSecondShell ? 1. : 2.;

        for(const size_t pairIndex : basisGroup.interactionIndices) {
          const auto& interaction = batch.orientedPairs[pairIndex];
          MatsT* outerOutput = nullptr;
          if(interaction.formOuterOutput) {
            const auto& output = batch.outputInfo[interaction.outerOutputIndex];
            if(output.receivesInnerSideWrites)
              outerOutput = workingData.threadLocalInnerElementCount == 0 ? nullptr :
                threadLocalOutput + output.threadOffset;
            else {
              // Each task owns one outer shell-pair block, so an outer-only
              // output has exactly one writer per element.
              outerOutput = rankOutput.data() + output.rankOffset;
            }
          }
          MatsT* innerOutput = nullptr;
          if(interaction.formInnerOutput) {
            const auto& output = batch.outputInfo[interaction.innerOutputIndex];
            innerOutput = workingData.threadLocalInnerElementCount == 0 ? nullptr :
              threadLocalOutput + output.threadOffset;
          }
          const auto& innerBasisData = batch.basisInfo[interaction.innerBasisInfoIndex];
          const BasisSet& innerBasis = *innerBasisData.basis;
          const double outerSchwarzBound = interaction.outerSchwarzBounds[outerFirstShell + outerSecondShell * outerBasis.nShell];

          for(const auto& innerShellPair : innerBasisData.significantShellPairs) {
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
            ++threadCandidateQuartets[threadIndex];
#endif

            const size_t innerFirstShell = innerShellPair.firstShell;
            const size_t innerSecondShell = innerShellPair.secondShell;
            const double combinedSchwarzBound = outerSchwarzBound * interaction.innerSchwarzBounds[innerFirstShell + innerSecondShell * innerBasis.nShell];

            // Screen each requested direction with the density that feeds it.
            bool digestOuterOutput = false;
            if(interaction.formOuterOutput) {
              const auto& sourceDensity = batch.densityInfo[interaction.innerDensityInfoIndex];
              digestOuterOutput = combinedSchwarzBound * sourceDensity.shellBlockInfinityNorms[innerFirstShell + innerSecondShell * innerBasis.nShell] >= interaction.schwarzThreshold;
            }

            bool digestInnerOutput = false;
            if(interaction.formInnerOutput) {
              const auto& sourceDensity = batch.densityInfo[interaction.outerDensityInfoIndex];
              digestInnerOutput = combinedSchwarzBound * sourceDensity.shellBlockInfinityNorms[outerFirstShell + outerSecondShell * outerBasis.nShell] >= interaction.schwarzThreshold;
            }

            if(not digestOuterOutput and not digestInnerOutput) continue;
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
            ++threadEvaluatedQuartets[threadIndex];
#endif

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(outerBasis.shells[outerFirstShell], outerBasis.shells[outerSecondShell], innerBasis.shells[innerFirstShell], innerBasis.shells[innerSecondShell]);

            const double* integralBuffer = engine.results()[0];
            if(integralBuffer == nullptr) continue;
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
            ++threadNonzeroQuartets[threadIndex];
#endif

            // The same ERI buffer feeds both contractions:
            // J_outer += (outer|inner) * D_inner and
            // J_inner += (outer|inner) * D_outer.
            const double innerPairDegeneracy = innerFirstShell == innerSecondShell ? 1. : 2.;
            const double quartetScale = interaction.interactionScale * outerPairDegeneracy * innerPairDegeneracy;
            size_t integralIndex = 0;

            for(size_t outerFirstFunction = 0; outerFirstFunction < outerShellPair.firstShellSize; ++outerFirstFunction)
            for(size_t outerSecondFunction = 0; outerSecondFunction < outerShellPair.secondShellSize; ++outerSecondFunction)
            for(size_t innerFirstFunction = 0; innerFirstFunction < innerShellPair.firstShellSize; ++innerFirstFunction)
            for(size_t innerSecondFunction = 0; innerSecondFunction < innerShellPair.secondShellSize; ++innerSecondFunction, ++integralIndex) {

              const size_t outerFirstBasisFunction = outerShellPair.firstBasisFunction + outerFirstFunction;
              const size_t outerSecondBasisFunction = outerShellPair.secondBasisFunction + outerSecondFunction;
              const size_t innerFirstBasisFunction = innerShellPair.firstBasisFunction + innerFirstFunction;
              const size_t innerSecondBasisFunction = innerShellPair.secondBasisFunction + innerSecondFunction;
              const double scaledIntegral = quartetScale * integralBuffer[integralIndex];

              if(digestOuterOutput) {
                const size_t sourceDensityElement = innerSecondBasisFunction + innerFirstBasisFunction * innerBasis.nBasis;
                const size_t outputElement = outerFirstBasisFunction + outerSecondBasisFunction * outerBasis.nBasis;
                outerOutput[outputElement] += MatsT(scaledIntegral * std::real(interaction.innerDensity[sourceDensityElement]));
              }

              if(digestInnerOutput) {
                const size_t sourceDensityElement = outerSecondBasisFunction + outerFirstBasisFunction * outerBasis.nBasis;
                const size_t outputElement = innerFirstBasisFunction + innerSecondBasisFunction * innerBasis.nBasis;
                innerOutput[outputElement] += MatsT(scaledIntegral * std::real(interaction.outerDensity[sourceDensityElement]));
              }
            }
          }
        }
      }
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      profile.contractionSeconds = tock(contractionStart);
#endif

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      const auto threadFoldStart = tick();
      for(size_t thread = 0; thread < threadCount; ++thread) {
        rankOutput[profileOffset + 0] += MatsT(threadCandidateQuartets[thread]);
        rankOutput[profileOffset + 1] += MatsT(threadEvaluatedQuartets[thread]);
        rankOutput[profileOffset + 2] += MatsT(threadNonzeroQuartets[thread]);
      }
#endif

      // Phase 4: fold only the outputs that had multiple inner-side writers.
      foldThreadLocalOutputs(threadOutputs, threadCount, workingData.threadLocalInnerElementCount, workingData.threadLocalToRankMap, rankOutput);
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      profile.threadFoldSeconds = tock(threadFoldStart);
#endif

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      const auto mpiReductionStart = tick();
#endif
      std::vector<MatsT> reducedOutput(reductionElementCount, MatsT(0.));
      // Phase 5: reduce all output matrices in one packed MPI operation.
      MPIReduce(rankOutput.data(), static_cast<int>(reductionElementCount), reducedOutput.data(), 0, comm);
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      profile.mpiReductionSeconds = tock(mpiReductionStart);
#endif

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      const auto outputAccumulationStart = tick();
#endif
      if(mpiRank == 0) {
        for(const auto& registeredOutput : batch.outputInfo) {
          const BasisSet& outputBasis = *batch.basisInfo[registeredOutput.basisInfoIndex].basis;
          accumulateHermitianOutput(reducedOutput.data() + registeredOutput.rankOffset, registeredOutput.output, outputBasis.nBasis);
        }

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
        profile.candidateQuartets = static_cast<size_t>(std::llround(
          std::real(reducedOutput[profileOffset + 0])));
        profile.evaluatedQuartets = static_cast<size_t>(std::llround(
          std::real(reducedOutput[profileOffset + 1])));
        profile.nonzeroQuartets = static_cast<size_t>(std::llround(
          std::real(reducedOutput[profileOffset + 2])));
#endif
      }
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      profile.outputSeconds = tock(outputAccumulationStart);
      if(mpiRank == 0) {
        const auto previousFlags = std::cout.flags();
        const auto previousPrecision = std::cout.precision();
        std::cout << "Batched direct interparticle profile\n"
                  << "  interactions / shell-pair tasks = "
                  << profile.activeInteractions << " / " << profile.workUnits << "\n"
                  << "  candidate / evaluated / nonzero quartets = "
                  << profile.candidateQuartets << " / "
                  << profile.evaluatedQuartets << " / "
                  << profile.nonzeroQuartets << "\n"
                  << "  preparation / engine / contract / fold / reduce / output = "
                  << std::fixed << std::setprecision(6)
                  << profile.preparationSeconds << " / "
                  << profile.engineSetupSeconds << " / "
                  << profile.contractionSeconds << " / "
                  << profile.threadFoldSeconds << " / "
                  << profile.mpiReductionSeconds << " / "
                  << profile.outputSeconds << " s\n"
                  << "  total = " << profile.totalSeconds << " s\n";
        std::cout.flags(previousFlags);
        std::cout.precision(previousPrecision);
      }
#endif

    // Turn threads for LA back on.
    SetLAThreads(laThreadCount);
    }

    // --------------------------------------------------------------
    // Gradient kernel: the prepared batch, its workingData, and the contract
    // --------------------------------------------------------------

    template <typename MatsT>
    struct GradientBatch {
      std::vector<BasisInfo> basisInfo;
      std::vector<DensityInfo<MatsT>> densityInfo;
      std::vector<OrientedPair<MatsT>> orientedPairs;
      size_t maximumPrimitiveCount = 0;
      size_t maximumAngularMomentum = 0;
    };

    template <typename MatsT, typename IntsT>
    GradientBatch<MatsT> prepareGradientBatch(
      const std::vector<DirectInterparticleGradJContraction<MatsT,IntsT>>& directInterparticleGradJContractions) {

      GradientBatch<MatsT> batch;
      auto& basisInfo = batch.basisInfo;
      auto& densityInfo = batch.densityInfo;
      auto& orientedPairs = batch.orientedPairs;
      orientedPairs.reserve(directInterparticleGradJContractions.size());

      for(const auto& gradJContraction : directInterparticleGradJContractions) {
        if(not gradJContraction.integrals)
          CErr("Null DirectTPI in batched direct interparticle gradient contraction.");
        if(gradJContraction.integrals->kernel() != TPI_KERNEL::Coulomb)
          CErr("Batched direct interparticle contraction supports the Coulomb TPI kernel only.");
        BasisSet& integralBasis1 = gradJContraction.integrals->basisSet();
        BasisSet& integralBasis2 = gradJContraction.integrals->basisSet2();
        if(integralBasis1.basisType != REAL_GTO or integralBasis2.basisType != REAL_GTO)
          CErr("Batched direct interparticle gradient contraction supports real GTO bases only.");

        BasisSet& subsystem1Basis = gradJContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis1 : integralBasis2;
        BasisSet& subsystem2Basis = gradJContraction.subsystemOrderMatchesIntegralOrder ?
          integralBasis2 : integralBasis1;

        const size_t subsystem1BasisIndex = registerBasis(basisInfo, subsystem1Basis);
        const size_t subsystem2BasisIndex = registerBasis(basisInfo, subsystem2Basis);
        const bool subsystem1IsOuter = basisInfo[subsystem1BasisIndex].significantShellPairs.size() >= basisInfo[subsystem2BasisIndex].significantShellPairs.size();

        OrientedPair<MatsT> pair;
        pair.interactionScale = gradJContraction.interactionScale;
        // Trace mode always digests both sides of the interaction, so both densities are always required.
        orientPair(pair, subsystem1IsOuter, subsystem1BasisIndex, subsystem2BasisIndex,
          gradJContraction.firstDensity, gradJContraction.secondDensity, true, true,
          nullptr, nullptr); // No Horns 1991 screening as that calls finite difference computeSchwarzGrad() for electron-proton pair

        pair.innerDensityInfoIndex = registerDensity(densityInfo, basisInfo, pair.innerBasisInfoIndex, pair.innerDensity, false);
        pair.outerDensityInfoIndex = registerDensity(densityInfo, basisInfo, pair.outerBasisInfoIndex, pair.outerDensity, false);

        batch.maximumPrimitiveCount = std::max(batch.maximumPrimitiveCount, std::max(subsystem1Basis.maxPrim, subsystem2Basis.maxPrim));
        batch.maximumAngularMomentum = std::max(batch.maximumAngularMomentum, std::max(subsystem1Basis.maxL, subsystem2Basis.maxL));
        orientedPairs.push_back(std::move(pair));
      }

      return batch;
    }

    struct GradientWorkingData {
      std::vector<OuterBasisGroup> outerBasisGroups;
      std::vector<ShellPairTask> shellPairTasks;
    };

    template <typename MatsT>
    GradientWorkingData packGradientWorkingData(GradientBatch<MatsT>& batch) {

      GradientWorkingData workingData;
      workingData.outerBasisGroups = buildOuterBasisGroups(batch.orientedPairs);
      workingData.shellPairTasks = buildShellPairTasks(
        batch.basisInfo, workingData.outerBasisGroups);
      return workingData;
    }

    // TRACE mode (mirrors DirectGradContraction::directScaffoldGradImpl).
    // Accumulate the scalar gradient without saving the derivative F_I matrices
    template <typename MatsT>
    std::vector<double> contractGradient(MPI_Comm comm,
      const size_t gradientComponentCount,
      const GradientBatch<MatsT>& batch,
      const GradientWorkingData& workingData) {

      if(gradientComponentCount > static_cast<size_t>(std::numeric_limits<int>::max()))
        CErr("Batched direct interparticle gradient reduction exceeds the MPI count limit.");

      const size_t threadCount = GetNumThreads();
      if(threadCount == 0)
        CErr("Batched direct interparticle gradient contraction received zero OpenMP threads.");

      const size_t laThreadCount = GetLAThreads();
      SetLAThreads(1);
      std::vector<libint2::Engine> engines(threadCount);
      engines[0] = libint2::Engine(libint2::Operator::coulomb, batch.maximumPrimitiveCount, batch.maximumAngularMomentum, 1);
      engines[0].set_precision(machineEpsilon);
      for(size_t thread = 1; thread < threadCount; ++thread)
        engines[thread] = engines[0];

      // One scalar accumulator per Cartesian component per thread.
      std::vector<std::vector<double>> threadGradients(threadCount,std::vector<double>(gradientComponentCount, 0.));

      const size_t mpiRank = MPIRank(comm);
      const size_t mpiSize = MPISize(comm);

      // Phase 3: each task owns one outer shell pair and evaluates the twelve
      // first-derivative ERI buffers once, tracing every requested interaction.
      #pragma omp parallel for schedule(dynamic, 1)
      for(long long shellPairTaskIndex = static_cast<long long>(mpiRank); shellPairTaskIndex < static_cast<long long>(workingData.shellPairTasks.size()); shellPairTaskIndex += static_cast<long long>(mpiSize)) {

        const size_t threadIndex = GetThreadID();
        auto& engine = engines[threadIndex];
        double* threadGradient = threadGradients[threadIndex].data();

        const auto& task = workingData.shellPairTasks[static_cast<size_t>(shellPairTaskIndex)];
        const auto& basisGroup = workingData.outerBasisGroups[task.outerBasisGroupIndex];
        const auto& outerBasisData = batch.basisInfo[basisGroup.outerBasisInfoIndex];
        const BasisSet& outerBasis = *outerBasisData.basis;
        const auto& outerShellPair = outerBasisData.significantShellPairs[task.outerShellPairIndex];

        const size_t outerFirstShell = outerShellPair.firstShell;
        const size_t outerSecondShell = outerShellPair.secondShell;
        const size_t outerCenter1 = outerBasis.mapSh2Cen[outerFirstShell];
        const size_t outerCenter2 = outerBasis.mapSh2Cen[outerSecondShell];
        const double outerPairDegeneracy = outerFirstShell == outerSecondShell ? 1. : 2.;

        for(const size_t pairIndex : basisGroup.interactionIndices) {
          const auto& interaction = batch.orientedPairs[pairIndex];

          const auto& innerBasisData = batch.basisInfo[interaction.innerBasisInfoIndex];
          const BasisSet& innerBasis = *innerBasisData.basis;

          for(const auto& innerShellPair : innerBasisData.significantShellPairs) {

            const size_t innerFirstShell = innerShellPair.firstShell;
            const size_t innerSecondShell = innerShellPair.secondShell;
            const size_t innerCenter1 = innerBasis.mapSh2Cen[innerFirstShell];
            const size_t innerCenter2 = innerBasis.mapSh2Cen[innerSecondShell];
            const double innerPairDegeneracy = innerFirstShell == innerSecondShell ? 1. : 2.;
            const double quartetScale = 0.5 * interaction.interactionScale * outerPairDegeneracy * innerPairDegeneracy;

            engine.compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 1>(outerBasis.shells[outerFirstShell], outerBasis.shells[outerSecondShell], innerBasis.shells[innerFirstShell], innerBasis.shells[innerSecondShell]);

            const auto& resultBuffers = engine.results();
            if(resultBuffers.size() < derivativeBufferCount or
               resultBuffers[0] == nullptr)
              continue;

            for(size_t derivativeBuffer = 0; derivativeBuffer < derivativeBufferCount; ++derivativeBuffer) {
              const double* integralBuffer = resultBuffers[derivativeBuffer];
              if(integralBuffer == nullptr) continue;

              const size_t xyz = derivativeBuffer % 3;
              const size_t shellIndex = derivativeBuffer / 3;
              const size_t gradientCenter = shellIndex == 0 ? outerCenter1 : shellIndex == 1 ? outerCenter2 : shellIndex == 2 ? innerCenter1 : innerCenter2;
              const size_t gradientComponent = gradientCenter * 3 + xyz;

              double quartetTrace = 0.;
              size_t integralIndex = 0;
              for(size_t outerFirstFunction = 0; outerFirstFunction < outerShellPair.firstShellSize; ++outerFirstFunction)
              for(size_t outerSecondFunction = 0; outerSecondFunction < outerShellPair.secondShellSize; ++outerSecondFunction) {

                const size_t outerFirstBasisFunction = outerShellPair.firstBasisFunction + outerFirstFunction;
                const size_t outerSecondBasisFunction = outerShellPair.secondBasisFunction + outerSecondFunction;
                // Loop invariant in the inner shell pair.
                const double outerDensityElement = std::real(interaction.outerDensity[outerSecondBasisFunction + outerFirstBasisFunction * outerBasis.nBasis]);

                for(size_t innerFirstFunction = 0; innerFirstFunction < innerShellPair.firstShellSize; ++innerFirstFunction)
                for(size_t innerSecondFunction = 0; innerSecondFunction < innerShellPair.secondShellSize; ++innerSecondFunction, ++integralIndex) {

                  const size_t innerFirstBasisFunction = innerShellPair.firstBasisFunction + innerFirstFunction;
                  const size_t innerSecondBasisFunction = innerShellPair.secondBasisFunction + innerSecondFunction;
                  const double innerDensityElement = std::real(interaction.innerDensity[innerSecondBasisFunction + innerFirstBasisFunction * innerBasis.nBasis]);

                  quartetTrace += outerDensityElement * innerDensityElement * integralBuffer[integralIndex];
                }
              }

              threadGradient[gradientComponent] += quartetScale * quartetTrace;
            }
          }
        }
      }

      // Phase 4: sum the thread accumulators, then reduce over the
      // communicator so every rank holds the full gradient contribution.
      std::vector<double> rankGradient(gradientComponentCount, 0.);
      for(size_t thread = 0; thread < threadCount; ++thread)
      for(size_t gradient = 0; gradient < gradientComponentCount; ++gradient)
        rankGradient[gradient] += threadGradients[thread][gradient];

      std::vector<double> reducedGradient(gradientComponentCount, 0.);
      MPIAllReduce(rankGradient.data(), static_cast<int>(gradientComponentCount), reducedGradient.data(), comm);

      // Turn threads for LA back on.
      SetLAThreads(laThreadCount);

      return reducedGradient;
    }

  } // namespace

  template <typename MatsT, typename IntsT>
  void BatchedDirectInterparticleJContraction<MatsT,IntsT>::JContract(
    MPI_Comm comm,
    const std::vector<DirectInterparticleJContraction<MatsT,IntsT>>& directInterparticleJContractions) {

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
    DirectProfile profile;
    const auto totalStart = tick();
    const auto preparationStart = tick();
#endif

    // Phase 1: register + orient every pair.
    auto batch = prepareCoulombBatch(directInterparticleJContractions);
    if(batch.orientedPairs.empty()) return;

    // Phase 2: pack the output buffers and build the shell-pair task list.
    auto workingData = packCoulombWorkingData(batch);

#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
    profile.preparationSeconds = tock(preparationStart);
    profile.totalSeconds = tock(totalStart);
#endif

    // Phases 3-5: contraction, thread fold, packed reduce, Hermitize.
    contractCoulomb(comm, batch, workingData
#ifdef PROFILE_BATCHED_DIRECT_INTERPARTICLE
      , profile
#endif
    );
  }

  template <typename MatsT, typename IntsT>
  std::vector<double> BatchedDirectInterparticleJContraction<MatsT,IntsT>::GradJContract(
    MPI_Comm comm,
    const size_t gradientComponentCount,
    const std::vector<DirectInterparticleGradJContraction<MatsT,IntsT>>& directInterparticleGradJContractions) {

    // Phase 1: register + orient every pair.
    auto batch = prepareGradientBatch(directInterparticleGradJContractions);
    if(batch.orientedPairs.empty())
      return std::vector<double>(gradientComponentCount, 0.);

    // Phase 2: build the interaction grouping and the task list.
    auto workingData = packGradientWorkingData(batch);

    // Phases 3-4: derivative contraction in trace mode, thread sum, all-reduce.
    return contractGradient(comm, gradientComponentCount, batch, workingData);
  }

  template void
  BatchedDirectInterparticleJContraction<double,double>::JContract(MPI_Comm,
    const std::vector<DirectInterparticleJContraction<double,double>>&);

  template void
  BatchedDirectInterparticleJContraction<dcomplex,double>::JContract(MPI_Comm,
    const std::vector<DirectInterparticleJContraction<dcomplex,double>>&);

  template void
  BatchedDirectInterparticleJContraction<dcomplex,dcomplex>::JContract(MPI_Comm,
    const std::vector<DirectInterparticleJContraction<dcomplex,dcomplex>>&);

  template std::vector<double>
  BatchedDirectInterparticleJContraction<double,double>::GradJContract(MPI_Comm,
    size_t, const std::vector<DirectInterparticleGradJContraction<double,double>>&);

  template std::vector<double>
  BatchedDirectInterparticleJContraction<dcomplex,double>::GradJContract(MPI_Comm,
    size_t, const std::vector<DirectInterparticleGradJContraction<dcomplex,double>>&);

  template std::vector<double>
  BatchedDirectInterparticleJContraction<dcomplex,dcomplex>::GradJContract(MPI_Comm,
    size_t, const std::vector<DirectInterparticleGradJContraction<dcomplex,dcomplex>>&);

} // namespace ChronusQ
