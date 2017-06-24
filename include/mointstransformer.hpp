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

#include <particleintegrals/twopints/incoreritpi.hpp>
#include <fields.hpp>
#include <singleslater.hpp>
#include <matrix.hpp>
#include <mcwavefunction/base.hpp>

#include <type_traits>

// #define _DEBUG_MOINTSTRANSFORMER_CACHE

namespace ChronusQ {

/**
 *  \brief enum types that enable fast integral transformation  
 *         for cases where two or more index are same 
 */

enum TPI_TRANS_DELTA_TYPE {
  NO_KRONECKER_DELTA,     // outputs (pq|rs)
  // Coulomb Type
  KRONECKER_DELTA_PQ,     // outputs (pp|rs) as T(p,r,s)
  KRONECKER_DELTA_RS,     // outputs (pq|rr) as T(p,q,r)
  KRONECKER_DELTA_PQ_RS,  // outputs (pp|rr) as T(p,r)
  // eXchange Type
  KRONECKER_DELTA_PS,     // outputs (pq|rp) as T(p,q,r)
  KRONECKER_DELTA_RQ,     // outputs (pr|rs) as T(p,r,s)
  KRONECKER_DELTA_PS_RQ,  // outputs (pr|rp) as T(p,r)
};

/**
 *  \brief Templated class to handle AO to MO integral transformation
 *  for both one-body and two-body terms 
 *
 *  WARNING: for 1C only alpha part of the MO (mo[0]) is used 
 *
 */

template<typename MatsT, typename IntsT>
class MOIntsTransformer {

protected:

  TPI_TRANSFORMATION_ALG TPITransAlg_;
  MPI_Comm comm_;
  SingleSlater<MatsT,IntsT> & ss_;
  
  // storage of cache integral intermidates
  IntegralsCollection ints_cache_ = IntegralsCollection(); 
  
  // variables for moints type
  std::vector<std::set<char>> symbol_sets_;
  std::vector<std::pair<size_t, size_t>> mo_ranges_;
  
public:

  // Constructor
  MOIntsTransformer() = delete;
  MOIntsTransformer( const MOIntsTransformer & ) = default;
  MOIntsTransformer( MOIntsTransformer && ) = default;
  
  /**
   *  MOIntsTransformer Constructor. Constructs a MOIntsTransformer object
   *
   *  \param [in] ss  ... SingleSlater reference, which provides AO integrals
   *                      for transformation into MO basis via Direc or InCore
   *  \param [in] alg ... Algorithm for two particle integral transformation 
   *                      options see include/integrals.hpp
   */                      
  MOIntsTransformer(SingleSlater<MatsT,IntsT> & ss,
    TPI_TRANSFORMATION_ALG alg = DIRECT_N6): comm_(ss.comm),
    ss_(ss), TPITransAlg_(alg) {
      
      if (ss.nC == 4 and alg == INCORE_N5) {
        auto & hamiltonianOptions = ss.fockBuilder->getHamiltonianOptions();
        if (hamiltonianOptions.Gaunt or
            hamiltonianOptions.DiracCoulombSSSS or 
            hamiltonianOptions.Gauge) 
        CErr("MOIntsTransformer for above Hamiltonian Options is NYI !");
      }

      if (alg == DIRECT_N5) {
        CErr("DIRECT N5 MOIntsTransformer NYI !");   
      }
      
      // set default MO ranges as single slater
      setMORanges();

  }
 
  ~MOIntsTransformer() { clearAllCache(); };
  void clearAllCache() { ints_cache_.clear(); };
  
  /* 
   * helper functions to handles MO ranges through
   *   a map of index char sets with {offset, N}
   */

  // setters 
  void resetMORanges() {
      symbol_sets_.clear();
      mo_ranges_.clear();
  }
  void addMORanges(const std::set<char> &, const std::pair<size_t, size_t> &);
  void setMORanges(size_t nFrozenCore = 0, size_t nFrozenVirt = 0);  
  void setMORanges(const MCWaveFunctionBase & );  // for MCWavefunction 
  
  // parsers
  std::pair<size_t,size_t> parseMOType(const char);
  std::vector<std::pair<size_t,size_t>> parseMOType(const std::string &);
  char getUniqueSymbol(char type);
  
  // printers
  void printOffSizes(const std::vector<std::pair<size_t,size_t>> &);
  void printMORangesSummary();
  
  /** 
   * Major Intefaces to obtain different intgrals in MO basis
   * See details in include/mointstransformer/impl.hpp
   */
  
  // Methods to transform HCore 
  void transformHCore(EMPerturbation &, MatsT * MOHCore, 
    const std::string & moType = "pq", 
    bool deltaPQ = false, const char coreIndex = '\0');
  
  // Methods to transform G(D)
  void transformGD(EMPerturbation &, const cqmatrix::Matrix<MatsT> &, bool,
    MatsT *, const std::string & moType = "pq", bool deltaPQ = false, 
    bool cacheAOGD = false, const std::string & cacheId = "");
  void transformGD(EMPerturbation &, const char, 
    MatsT *, const std::string & moType = "pq", bool deltaPQ = false, 
    bool cacheAOGD = false, const std::string & cacheId = "");
  
  // Methods to transform TPI 
  void transformTPI(EMPerturbation & pert, MatsT* MOTPI, 
    const std::string & moType = "pqrs", 
    bool cacheIntermediates = true, bool withExchange = false,
    TPI_TRANS_DELTA_TYPE delta = NO_KRONECKER_DELTA);
  
  // MPI-enabled direct transform TPI
  void directTransformTPI(EMPerturbation & pert, MatsT* MOTPI,
    const std::string & moType = "pqrs");
  
  void directTransformTPI(EMPerturbation & pert,
    MatsT* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes);
 
 // Helper function used during transformation 
  std::shared_ptr<OnePInts<MatsT>> formAOHCore(EMPerturbation &, 
    bool cacheAOHCore = true, const char coreIndex = '\0');
  void subsetTransformHCore(EMPerturbation &, 
    const std::vector<std::pair<size_t,size_t>> &, 
    MatsT*, bool, const char coreIndex = '\0');
  void subsetTransformOPI(const std::vector<std::pair<size_t,size_t>> &,
    const OnePInts<MatsT> & AOOPI, MatsT*, bool); 
  
  std::shared_ptr<OnePInts<MatsT>> formAOGD(EMPerturbation &, 
    const cqmatrix::Matrix<MatsT> &, bool, bool cacheAOGD = false, 
    const std::string & cacheId = "");
  std::shared_ptr<cqmatrix::Matrix<MatsT>> formInactDen(const char coreIndex);
  void subsetTransformGD(EMPerturbation &, 
    const cqmatrix::Matrix<MatsT> &, bool, MatsT*, 
    const std::vector<std::pair<size_t,size_t>> &, 
    bool, bool cacheAOGD = false, const std::string & cacheId = ""); 
  
  void subsetTransformTPISSFockN6(EMPerturbation &, 
    const std::vector<std::pair<size_t,size_t>> &, 
    MatsT*, const std::string &, bool,
    TPI_TRANS_DELTA_TYPE delta = NO_KRONECKER_DELTA);
  
  void performSubsetTransformTPISSFockN6(EMPerturbation &, 
    const std::vector<std::pair<size_t,size_t>> &, 
    MatsT*, const std::string &, bool,
    TPI_TRANS_DELTA_TYPE delta = NO_KRONECKER_DELTA);
  
  std::shared_ptr<InCore4indexTPI<MatsT>> formAOTPIInCore(bool);
  void subsetTransformTPIInCoreN5(const std::vector<std::pair<size_t,size_t>> &, 
    MatsT*, bool, TPI_TRANS_DELTA_TYPE delta = NO_KRONECKER_DELTA);


}; // class MOIntsTransformer

} // namespace ChronusQ
