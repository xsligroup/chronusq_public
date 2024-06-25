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
#include <cerr.hpp>
#include <singleslater.hpp>
#include <mcwavefunction.hpp>
#include <orbitalrotation.hpp>
//#include <particleintegrals/twopints/incore4indextpi.hpp>

//#define _DEBUG_MCSCF

namespace ChronusQ {
  
  enum CIDiagonalizationAlgorithm {
    CI_FULL_MATRIX,
    CI_DAVIDSON,
    CI_GPLHR,
  }; // struct CIDiagonalizationAlgorithm
  
  // Settings 
  struct MCSCFSettings {
     
     // CISettings
     CIDiagonalizationAlgorithm ciAlg = CI_FULL_MATRIX; 
     
     // for davidson and gplhr
	 size_t maxCIIter        = 128;        
     double ciVectorConv     = 1.0e-6;    
     size_t maxDavidsonSpace = 50;
     size_t nDavidsonGuess   = 3;
     std::vector<std::pair<double, size_t>> energyRefs;

     // SCF Settings 
     bool doSCF           = false;
     bool doIVOs          = false;
     
     size_t maxSCFIter         = 0;
     double scfEnergyConv      = 1.0e-8; 
     double scfGradientConv    = 1.0e-4;
     OrbitalRotationSettings ORSettings;
  
     MCSCFSettings() {
       ORSettings.rotate_within_correlated = false;
     }
     
     MCSCFSettings(const MCSCFSettings &) = default;
     MCSCFSettings(MCSCFSettings &&) = default;

     void print(bool, size_t);
  }; // struct MCSCFSettings 
  
  template <typename MatsT, typename IntsT>
  class CISolver { 
     
  protected:
    
    CIDiagonalizationAlgorithm alg_ = CI_FULL_MATRIX;

    // for davidson and gplhr
    size_t maxIter_    = 128;       /// < Max Number of CI iteration 
    double vectorConv_ = 1.0e-6;    /// < Convergence criteria in terms of vector residue norm
    size_t maxDavidsonSpace_ = 50;  /// < Max davidson space in terms of n times of NRoots  
    size_t nDavidsonGuess_ = 3;     /// < number of guess in the intial davidson first a few iterations 
    std::vector<std::pair<double, size_t>> energyRefs_;

    void davidsonGS(size_t, size_t, MatsT *, MatsT *);
    void davidsonPC(size_t, size_t, MatsT *, MatsT *, MatsT *, dcomplex *);
  
  public:
    
    // default Constructor
    CISolver() = default;
    CISolver(CIDiagonalizationAlgorithm alg, size_t maxIter = 128,
      double vectorConv =  1.0e-6, size_t maxDSpace = 50, size_t nDGuess = 3,
      std::vector<std::pair<double, size_t>> eRefs = {}) {
      switchAlgorithm(alg, maxIter, vectorConv, maxDSpace, nDGuess, eRefs);
    };

    // typeconversion
	template <typename MatsU>
    CISolver(const CISolver<MatsU,IntsT> &);
	
	template <typename MatsU>
    CISolver(CISolver<MatsU,IntsT> &&);
    
	~CISolver() = default;
  
    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<CISolver<MatsU,IntsT>>
    convert(const std::shared_ptr<CISolver<MatsT,IntsT>>&);
    
    // switch Algrithm;
    void switchAlgorithm(CIDiagonalizationAlgorithm alg, 
      size_t maxIter = 128, double vectorConv = 1.0e-6, 
      size_t maxDSpace = 50, size_t nDGuess = 3,
      std::vector<std::pair<double, size_t>> eRefs = {}) {
        alg_ = alg;
        maxIter_ = maxIter;
        vectorConv_ = vectorConv;
        maxDavidsonSpace_ = maxDSpace;
        nDavidsonGuess_ = nDGuess;
        if (!eRefs.empty()) energyRefs_ = eRefs;
    }

    // solve CI
	void solveCI(MCWaveFunction<MatsT,IntsT> &,EMPerturbation &);
  
  }; // class CISolver

  template <typename MatsT, typename IntsT>
  class MCSCF : public MCWaveFunction<MatsT,IntsT> {
  
  protected:
    // Useful Typedefs
    typedef MatsT *                   oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;
  
  public:
    
    MCSCFSettings settings;
    std::shared_ptr<CISolver<MatsT,IntsT>>        ciSolver  = nullptr;
    std::shared_ptr<OrbitalRotation<MatsT,IntsT>> moRotator = nullptr;  

    // Reduced density Matrices (RDMs) only span over correlated space
    // SOI: state of interest, for orbital rotation
    // it's either state specific or state averaged RDM
    std::shared_ptr<cqmatrix::Matrix<MatsT>>    oneRDMSOI = nullptr;
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = nullptr;
    // Disable default, copy and move constructors
    MCSCF()              = delete;
    MCSCF(const MCSCF &) = delete;
    MCSCF(MCSCF &&)      = delete;

    // Constructors
    
    /**
     *  \brief MCSCF Constructor.
     *
     *  Stores references to a "reference" SingleSlater object
     *  and makes a copy of the reference into a complex
     *  SingleSlater object for the propagation.
     */ 
    template <typename MatsU>
    MCSCF(SingleSlater<MatsU,IntsT> & ref, size_t NS) : 
        MCWaveFunction<MatsT,IntsT>(ref, NS) { };  // MCSCF constructor
  
    ~MCSCF(){ dealloc(); }

    // MCSCF procedural functions
    void run(EMPerturbation &);       // From MCWaveFunctionBase
    
    // compute all RDM (and state average) if no inputs 
    void computeOneRDM();
    void computeOneRDM(size_t);
    void computeTwoRDM();
    void computeTwoRDM(size_t);

    void saveCurrentStates(bool prop = false);

    void printStateEnergy();
    void printMCSCFHeader(EMPerturbation &);
    void printMCSCFFooter();

    // Memory functions
    void alloc();
    void dealloc();

  }; // class MCSCF

}; // namespace ChronusQ



