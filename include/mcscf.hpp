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
#include <singleslater/neoss.hpp>
#include <mcscf/base.hpp>
#include <mcwavefunction.hpp>
#include <orbitalrotation.hpp>
//#include <particleintegrals/twopints/incore4indextpi.hpp>

//#define _DEBUG_MCSCF

namespace ChronusQ {
  
 
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

    // get Algorithm
    CIDiagonalizationAlgorithm getAlg(){ return alg_; }

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
	virtual void solveCI(MCWaveFunction<MatsT,IntsT> &,EMPerturbation &);
  
  }; // class CISolver

  template <typename MatsT, typename IntsT>
  class MCSCF : public MCSCFBase
  {
    // Potentially for compatibility reasons
    //template<typename MatsT, typename IntsT>
    friend class MCWaveFunction<MatsT,IntsT>;

    private:

      std::shared_ptr<MCWaveFunction<MatsT,IntsT>> mcwfn_;

    protected:

      std::shared_ptr<CISolver<MatsT,IntsT>>        ciSolver  = nullptr;
      std::shared_ptr<OrbitalRotation<MatsT,IntsT>> moRotator = nullptr;  
      std::shared_ptr<cqmatrix::Matrix<MatsT>>    oneRDMSOI = nullptr;
      std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = nullptr;

    public:

      MCSCF() = delete;
      MCSCF(const MCSCF&) = delete;
      MCSCF(MCSCF &&) = delete;

      MCSCF(std::shared_ptr<MCWaveFunction<MatsT,IntsT>>mcwfn):
        MCSCFBase(mcwfn->comm,mcwfn->NStates),
        mcwfn_(mcwfn)
        {
          this->NDet = mcwfn_->NDet;
          this->NStates = mcwfn_->NStates;
        };

      // MCSCF procedural functions
      void run(EMPerturbation & pert)override;
      void runCube(std::vector<std::shared_ptr<CubeGen>>) override;

      // Functions that generally are wrappers around MCWavefunction, 
      // but might be overwritten by classes which derive from MCSCF
      virtual void transformInts(EMPerturbation & pert);
      virtual void computeMultipole();
      virtual void populationAnalysis();
      virtual void spinAnalysis();
      virtual double oscillator_strength(size_t, size_t);
      virtual void formNaturalOrbitals();

      // compute all RDM (and state average) if no inputs 
      virtual void computeOneRDM();
      virtual void computeOneRDM(size_t);
      void computeTwoRDM();
      void computeTwoRDM(size_t);

      void saveCurrentStates(bool prop = false);

      void printStateEnergy();
      virtual void printMCSCFHeader(EMPerturbation &);
      virtual void printMCSCFFooter();

      // Memory functions
      virtual void alloc();
      void dealloc();

  }; // class MCSCF

}; // namespace ChronusQ
