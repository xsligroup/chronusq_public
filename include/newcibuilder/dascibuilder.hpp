/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#ifndef CIBUILDER_HPP
#error This file can only be included in newcibuilder.hpp
#endif

#include <particleintegrals/dasints.hpp>
#include <util/scratch.hpp>

namespace ChronusQ {

enum DASCISigma2eContAlg {
  NAIVE = 0,  // naive loop
  KNOWLESHANDY = 1, 
  OLSENROOS = 2,  
  FRISCHLI = 3,
  SMALLBLOCKHAMILTONIAN = 4,
  LARGEBLOCKHAMILTONIAN = 5,
}; // enum DASCISigma2eContAlg

/**
 *  \brief The DASCIBuilder Class.
 *  
 *  The Distributed Active Space Configuration Interaction Builder
 *
 *  Assume the every determinantal category are full
 *
 */
template <typename MatsT>
class DASCIBuilder: public NewCIBuilder<MatsT> {

protected:
  DASCISigma2eContAlg twoEcontAlg_ = DASCISigma2eContAlg::NAIVE;
  
  void estimateMemoryForIntermediates();

public:
  // Constructors

  // default constructors
  DASCIBuilder() = delete;
  DASCIBuilder(MPI_Comm comm,
    const std::shared_ptr<const IntegralsCollection> moints,
    const DeterminantFactory& detF):
    NewCIBuilder<MatsT>(comm, moints, detF) {};

  DASCIBuilder(const DASCIBuilder<MatsT>& other) = default;
  DASCIBuilder(DASCIBuilder<MatsT>&& other)      = default;
  ~DASCIBuilder() = default;
  
  // DASCIBuilder Interfacing Functions
  void setSigma2eContractionAlgorithm(std::string alg);
  
  void buildSigma1e(const LocalCIVectorsView<const MatsT>& C, 
      const LocalCIVectorsView<MatsT>& Sigma) const override;
  void buildSigma2e(const LocalCIVectorsView<const MatsT>& C, 
      const LocalCIVectorsView<MatsT>& Sigma) const override;
  
  void build1TDM(
      const LocalCIVectorsView<const MatsT>& CBra, 
      const LocalCIVectorsView<const MatsT>& CKet, 
      cqmatrix::Matrix<MatsT>& oneTDM,
      const double scale = 1.0) const override;
  
  void build2TDM(
      const LocalCIVectorsView<const MatsT>& CBra, 
      const LocalCIVectorsView<const MatsT>& CKet, 
      InCore4indexTPI<MatsT>& twoTDM,
      const double scale = 1.0) const override;
  
}; // class DASCIBuilder

// helper class
namespace DASCISigma2eBuilder {
  template <typename MatsT>
  void buildNaive(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      const std::vector<size_t>& KExOffs, const std::vector<size_t>& LExOffs, 
      std::shared_ptr<TensorLooper>& NonExLooper, const double symmFact);

  template <typename MatsT>
  void buildKnowlesHandy(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      const std::vector<size_t>& KExOffs, std::shared_ptr<TensorLooper>& KNonExLooper, 
      const std::vector<size_t>& LExOffs, std::shared_ptr<TensorLooper>& LNonExLooper,
      MatsT* omega, MatsT* lambda, MatsT* X, const double symmFactor);

  template <typename MatsT>
  void buildOlsenRoos(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      const std::vector<size_t>& KExOffs, std::shared_ptr<TensorLooper>& KNonExLooper, 
      const std::vector<size_t>& LExOffs, std::shared_ptr<TensorLooper>& LNonExLooper,
      MatsT* omega, MatsT* sSCR, MatsT* X, const double symmFactor);

  template <typename MatsT>
  void buildFrischLi(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      const std::vector<size_t>& KExOffs, std::shared_ptr<TensorLooper>& KNonExLooper, 
      const std::vector<size_t>& LExOffs, std::shared_ptr<TensorLooper>& LNonExLooper,
      MatsT* cSCR, MatsT* lambda, MatsT* X, const double symmFactor);

  template <typename MatsT>
  void buildSmallBlockHamiltonian(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      const std::vector<size_t>& KExOffs, std::shared_ptr<TensorLooper>& KNonExLooper, 
      const std::vector<size_t>& LExOffs, std::shared_ptr<TensorLooper>& LNonExLooper,
      MatsT* cSCR, MatsT* sSCR, MatsT* X, const double symmFactor);

  template <typename MatsT>
  void buildLargeBlockHamiltonian(
      size_t nVec, const MatsT* C, size_t LDC, MatsT* Sigma, size_t LDS,
      const DASTwoPInts<MatsT>& s2e, DoubleFullCD1eExListGenerator& double1eExListsGen,
      std::shared_ptr<TensorLooper>& KExLooper, std::shared_ptr<TensorLooper>& KNonExLooper, 
      std::shared_ptr<TensorLooper>& LExLooper, std::shared_ptr<TensorLooper>& LNonExLooper,
      MatsT* cSCR, MatsT* sSCR, MatsT* X, const double symmFactor);
}; // namespace DASCISigma2eBuilder

} // namespace ChronusQ
