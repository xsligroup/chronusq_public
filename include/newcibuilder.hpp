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
#include <matrix.hpp>
#include <posthartreefock.hpp>
#include <detfactory.hpp>
#include <detfactory/localcivectorsview.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <util/scratch.hpp>

#define CIBUILDER_HPP

//#define DEBUG_CIBuilder

namespace ChronusQ {
  
/* 
 * Brief Definition of CIBuilder Class
 *
 * For different contraction algorithm
 *
 */
template <typename MatsT>
class NewCIBuilder { 

 protected:
  MPI_Comm comm_;
  const std::shared_ptr<const IntegralsCollection> moints_; // integrals
  const DeterminantFactory& detFactory_;      // for CI vectors

  // for scratches
  std::unordered_map<std::string, size_t> nSCR_;
  
  void estimateLocalMemoryForMPIComm();

 public:
  
  // default Constructor
  NewCIBuilder() = delete;
  
  NewCIBuilder(MPI_Comm comm, 
      const std::shared_ptr<const IntegralsCollection> moints,
      const DeterminantFactory& detF):
      comm_(comm), moints_(moints), detFactory_(detF) { }


  NewCIBuilder(const NewCIBuilder<MatsT>&) = default;
  NewCIBuilder(NewCIBuilder<MatsT>&&)      = default;
  ~NewCIBuilder() = default;

  // Virtual CI Build Functions
  virtual void buildDiagH(DistributedVectors<MatsT>& diagH, 
      const CategoricalSpace& categoricalSpace) const;

  virtual void buildFullH(DistributedVectors<MatsT>& fullH) const;

  virtual void buildSigma(size_t nVec, const DistributedVectors<MatsT>& C,
      size_t CShift, DistributedVectors<MatsT>& Sigma, size_t SigmaShift) const;

#ifdef CQ_ENABLE_SPARSE
  virtual void largestDiagH(DistributedSparseVectors<MatsT>& C,
    const CategoricalSpace& categoricalSpace, size_t K, std::vector<size_t>& kIndices, std::vector<MatsT>& kValues, bool (*comp)(const MatsT&,const MatsT&)) const;

  virtual void buildSigma(size_t nVec, const DistributedSparseVectors<MatsT>& C,
      size_t CShift, DistributedSparseVectors<MatsT>& Sigma, size_t SigmaShift, dcomplex* curEigenvalues, double eps) const;
#endif

  virtual void buildSigma1e(const LocalCIVectorsView<const MatsT>& C, 
      const LocalCIVectorsView<MatsT>& Sigma) const = 0;
  virtual void buildSigma2e(const LocalCIVectorsView<const MatsT>& C,
      const LocalCIVectorsView<MatsT>& Sigma) const = 0;

#ifdef CQ_ENABLE_SPARSE
  virtual void buildSigma1e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC, const LocalCISparseVectorsView<MatsT>& Sigma, HashSparseMatrix<MatsT>& SCRSigma_hash, dcomplex* curEigenvalues, double eps) const = 0;
  virtual void buildSigma2e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC,
      const LocalCISparseVectorsView<MatsT>& Sigma, HashSparseMatrix<MatsT>& SCRSigma_hash, dcomplex* curEigenvalues, double eps) const = 0;




  virtual void formSubA(size_t nVec, size_t nTot, const DistributedSparseVectors<MatsT>& C,
      size_t CShift, std::vector<MatsT>& VAV) const;
    virtual void formSubA2e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC,
      std::vector<MatsT>& VAV, size_t iVec, size_t jVec) const = 0;

  virtual void formSubA1e(const LocalCISparseVectorsView<MatsT>& C, const LocalCISparseVectorsView<MatsT>& myC,
      std::vector<MatsT>& VAV, size_t iVec, size_t jVec) const = 0;
#endif



  // build 1TDM and 2TDM
  virtual void buildTDM(const DistributedVectors<MatsT>& CBra, 
      const DistributedVectors<MatsT>& CKet, 
      const std::pair<size_t, size_t>& stateIndices,  
      std::shared_ptr<cqmatrix::Matrix<MatsT>> oneTDM,
      bool oneTDMIncrement = false, double oneTDMScale = 1.0,
      std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM = nullptr,
      bool twoTDMIncrement = false, double twoTDMScale = 1.0
#ifdef CQ_ENABLE_MPI 
      , bool reduceOneTDM = true, bool reduceTwoTDM = true 
#endif
      ) const;
 
#ifdef CQ_ENABLE_SPARSE 
   virtual void buildTDM(const DistributedSparseVectors<MatsT>& CBra,
      const DistributedSparseVectors<MatsT>& CKet,
      const std::pair<size_t, size_t>& stateIndices,
      std::shared_ptr<cqmatrix::Matrix<MatsT>> oneTDM,
      bool oneTDMIncrement = false, double oneTDMScale = 1.0,
      std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM = nullptr,
      bool twoTDMIncrement = false, double twoTDMScale = 1.0
#ifdef CQ_ENABLE_MPI
      , bool reduceOneTDM = true, bool reduceTwoTDM = true
#endif
      ) const;
#endif

  virtual void build1TDM(
      const LocalCIVectorsView<const MatsT>& CBra, 
      const LocalCIVectorsView<const MatsT>& CKet, 
      cqmatrix::Matrix<MatsT>& oneTDM,
      const double scale = 1.0) const = 0;
  
  virtual void build2TDM(
      const LocalCIVectorsView<const MatsT>& CBra, 
      const LocalCIVectorsView<const MatsT>& CKet, 
      InCore4indexTPI<MatsT>& twoTDM,
      const double scale = 1.0) const = 0;

#ifdef CQ_ENABLE_SPARSE
  virtual void build1TDM(
      const LocalCISparseVectorsView<MatsT>& CBra,
      const LocalCISparseVectorsView<MatsT>& CKet,
      cqmatrix::Matrix<MatsT>& oneTDM,
      const double scale = 1.0) const = 0;

  virtual void build2TDM(
      const LocalCISparseVectorsView<MatsT>& CBra,
      const LocalCISparseVectorsView<MatsT>& CKet,
      InCore4indexTPI<MatsT>& twoTDM,
      const double scale = 1.0) const = 0;
#endif
}; // class CIBuilder

} // namespace ChronusQ

// Include declaration for specialization of CIBuilder
#include <newcibuilder/dascibuilder.hpp>
