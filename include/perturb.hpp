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
#include <util/files.hpp>
#include <itersolver.hpp>
#include <mcwavefunction.hpp>
#include <mcscf.hpp>

#define _DEBUG_PT2

namespace ChronusQ {

  struct PerturbOptions {

    bool doFull = false;   /// < full matrix build for left-hand-side
    bool doIter = true;  /// < use GMRES iterative solver
    bool saFock = true;  /// < use state-averaged density matrix to build Fock.
    bool extendMS = false;  /// < use extended multi-state scheme
    // TODO: GVVPT2
    bool doGVV = false;  /// < build H_eff within the whole primary space

    // level shift parameter
    double levelShift = 0.0;
    // imaginary shift
    double imaginaryShift = 0.0;

    // frozen orbital settings
    size_t frozenCore = 0;
    size_t frozenVirtual = 0;
    std::string selectVirtual;

    // for iterative solver
    size_t maxIter = 32;  /// < Max Number of iterations
    double convCrit = 1.0e-6;     /// < Convergence criteria


  };


  template <typename MatsT, typename IntsT>
  class PERTURB : public MCWaveFunction<MatsT,IntsT> {

  protected:
    // Useful Typedefs
    typedef MatsT *                   oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

    size_t CIsize;  /// Model space size: # CAS determinants
    size_t SDsize;  /// External space size (singles and doubles)

    oper_t E0_; /// < zero-order energies

    std::shared_ptr<MCSCF<MatsT,IntsT>> refMCwfn = nullptr; /// < MCSCF reference
    std::shared_ptr<cqmatrix::Matrix<MatsT>> LHS_ = nullptr; /// < LHS matrix
    std::shared_ptr<cqmatrix::Matrix<MatsT>> ptFock_ = nullptr; /// < Fock matrix
    oper_t diagLHS = nullptr;

    std::vector<MatsT> shifts = { MatsT(0.) }; /// < for GMRES setup


  public:

    PerturbOptions PTopts;
    std::vector<size_t> SoI_; /// < State of Interest

    // Disable default, copy and move constructors
    PERTURB()                = delete;
    PERTURB(const PERTURB &) = delete;
    PERTURB(PERTURB &&)      = delete;

    // Constructors

    /**
     *  \brief PERTURB Constructor
     *  
     *  Obtain references from a "reference" MCSCF object
     *
     */
    PERTURB(std::shared_ptr<MCSCF<MatsT,IntsT>> refMCwfn, std::vector<size_t> SoI) : 
        MCWaveFunction<MatsT,IntsT>(refMCwfn->reference(), SoI.size()), 
        SoI_(SoI), refMCwfn(refMCwfn) { 

      // MRPT with 1C reference not implemented
      if (refMCwfn->reference().nC==1)
        CErr("1c + MRPT not yet implemented");
      if (refMCwfn->MOPartition.scheme == RAS)
        CErr("RAS + MRPT not yet implemented");

    }

    ~PERTURB() { dealloc(); }

    std::shared_ptr<MCSCF<MatsT,IntsT>> mcwfnRef() { return refMCwfn; }
//    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> ptwfn() { return PTwfn; }

    // PERTURB procedural functions
    void run(EMPerturbation &);
    void runSingleState(size_t);

    // Initialize 1-order perturbed wavefunction
    void PTInitialize();

    // Single state PT2 functions
    void buildFock(cqmatrix::Matrix<MatsT> &, cqmatrix::Matrix<MatsT> &);
    MatsT computeZeroE(cqmatrix::Matrix<MatsT> &, cqmatrix::Matrix<MatsT> &);
    void buildFull(cqmatrix::Matrix<MatsT> &, MatsT);
    void buildRHS(oper_t, size_t);
    void computeHV(oper_t, oper_t);
    MatsT computeCHV(oper_t, oper_t);
    MatsT computeShift(size_t, MatsT);
    MatsT computeShiftCorrection(size_t);

    void solve(oper_t, size_t);
    // itersolve relevant
    void buildDiag();
    void PT2preCond(size_t, MatsT, SolverVectors<MatsT> &, SolverVectors<MatsT> &, oper_t);
    void buildAV(oper_t, oper_t, MatsT);
    void iterSolve(oper_t, MatsT);

    // Multi-state PT2 functions
    void rotateCIVecs();
    void buildHeff(cqmatrix::Matrix<MatsT> &);
    void diagHeff(cqmatrix::Matrix<MatsT> &);

    void saveCurrentStates();

    //printing functions
    void printPTHeader();
    void printPTFooter();

    // Memory functions
    void alloc();
    void dealloc();

  }; // class PERTURB

}; // namespace ChronusQ


