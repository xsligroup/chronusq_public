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
#include <wavefunction/base.hpp>
#include <integrals.hpp>
#include <fields.hpp>
#include <util/files.hpp>

// #define TEST_MOINTSTRANSFORMER

namespace ChronusQ {

  // Type of SingleSlater object
  enum RefType {
    isRawRef,  // non-specified
    isRRef,    // RHF/DFT
    isURef,    // UHF/DFT
    isRORef,   // ROHF/DFT
    isGRef,    // GHF/DFT
    isTwoCRef, // Two-component
    isX2CRef,  // X2C
    isFourCRef // Four-component
  };

  // A struct that stores reference information
  struct RefOptions {

    std::string RCflag = "REAL"; // Real or Complex

    RefType refType = isRRef;    // R/U/G/2c/X2C/4c

    bool isKSRef = false;        // HF or DFT
    bool isEPCRef = false;       // NEO-KS or not

    size_t nC = 1;               // number of component
    bool iCS = true;             // closed shell or not

    std::string funcName;        // DFT functional name
  };


  /**
   *  \brief A datastructure to hold the information
   *  pertaining to the control of the Kohn--Sham
   *  numerical integration.
   */
  struct IntegrationParam {
    double epsilon      = 1e-12; ///< Screening parameter
    size_t nAng         = 302;   ///< # Angular points
    size_t nRad         = 100;   ///< # Radial points
    size_t nRadPerBatch = 4;     ///< # Radial points / macro batch
  };

  /**
   *  The Single Slater guess types
   */ 
  enum SS_GUESS {
    CORE,
    SAD,
    TIGHT,
    RANDOM,
    READMO,
    READDEN,
    FCHKMO,
    // Specific Guess Options For NEO
    NEOTightProton,
    NEOConvergeClassical
  };
  /**
   *  The types of steps for the SCF
   */
  enum SCF_STEP { _CONVENTIONAL_SCF_STEP, _NEWTON_RAPHSON_STEP };

  /**
   *  SCF Algorithms
   */
  enum SCF_ALG { _CONVENTIONAL_SCF, _NEWTON_RAPHSON_SCF, _SKIP_SCF };

  enum DIIS_ALG {
    CDIIS,    ///< Commutator DIIS
    EDIIS,    ///< Energy DIIS
    CEDIIS,   ///< Commutator & Energy DIIS
    NONE = -1
  };

  enum NR_APPROX {
    FULL_NR,
    QUASI_BFGS,
    QUASI_SR1,
    GRAD_DESCENT
  };

  /**
   *  \brief A struct to hold the information pertaining to
   *  the control of an SCF procedure.
   *
   *  Holds information like convergence critera, DIIS settings, 
   *  max iterations, etc.
   */ 
  struct SCFControls {

    // Convergence criteria
    double denConvTol = 1e-8;  ///< Density convergence criteria
    double eneConvTol = 1e-10; ///< Energy convergence criteria
    double FDCConvTol = 1e-8; ///< Gradient convergence criteria

    // TODO: need to add logic to set this
    // Extrapolation flag for DIIS and damping
    bool doExtrap = true;     ///< Whether to extrapolate Fock matrix

    // Algorithm and step
    SCF_STEP  scfStep = _CONVENTIONAL_SCF_STEP;
    SCF_ALG   scfAlg  = _CONVENTIONAL_SCF;
    NR_APPROX nrAlg   = QUASI_BFGS;         ///< NR approximation(i.e. quasi-Newton)
    double nrTrust = 0.1;                   ///< Initial trust region for NR SCF 
    double nrLevelShift = 0.; 				///< Level shift for diagonal hessian

    // Guess Settings
    SS_GUESS guess = SAD;
    SS_GUESS prot_guess = NEOTightProton;

    // DIIS settings 
    DIIS_ALG diisAlg = CEDIIS; ///< Type of DIIS extrapolation 
    size_t nKeep     = 10;     ///< Number of matrices to use for DIIS
    double cediisSwitch = 0.05; ///< When to switch from EDIIS to CDIIS

    // Static Damping settings
    bool   doDamp         = false;           ///< Flag for turning on damping
    double dampStartParam = 0.7;            ///< Starting damping parameter
    double dampParam      = dampStartParam; ///< Current Damp parameter 
    double dampError      = 1e-3; ///< Energy oscillation to turn off damp

    // Incremental Fock build settings
    bool   doIncFock = false; ///< Whether to perform an incremental fock build
    size_t nIncFock  = 20;   ///< Restart incremental fock build after n steps

    // Misc control
    size_t maxSCFIter = 128; ///< Maximum SCF iterations.

    // Printing
    size_t printMOCoeffs = 0;
    size_t printLevel = 1;
    std::string refLongName_;
    std::string refShortName_;

  }; // SCFControls struct

  class SingleSlaterBase;

  struct SingleSlaterOptions {

    HamiltonianOptions hamiltonianOptions;

    RefOptions refOptions;

    IntegrationParam intParam;

    SCFControls scfControls;

    std::shared_ptr<SingleSlaterBase> buildSingleSlater(
        std::ostream &out, CQMemManager &mem,
        Molecule &mol, BasisSet &basis,
        std::shared_ptr<IntegralsBase> aoints) const;

  };

  /**
   *  \brief The SingleSlaterBase class. The abstraction of information
   *  relating to the SingleSlater class which are independent of storage
   *  type.
   *
   *  Specializes WaveFunctionBase interface.
   *
   *  See SingleSlater for further docs.
   */ 
  class SingleSlaterBase : virtual public WaveFunctionBase {

  protected:

  private:
  public:

    std::string refLongName_;  ///< Long form of the reference name
    std::string refShortName_; ///< Short form of the reference name

    // Save / Restart File
    SafeFile savFile;

    // Fchk File
    std::string fchkFileName;

    // Scratch Bin File
    std::string scrBinFileName;
       
    // Print Controls
    size_t printLevel; ///< Print Level

    // Current Timings
    double GDDur;

    // SCF Variables
    SCFControls    scfControls; ///< Controls for the SCF procedure

    // Pair function for SingleSlater MO swap
    std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;

    // Constructors (all defaulted)
    SingleSlaterBase(const SingleSlaterBase &) = default;
    SingleSlaterBase(SingleSlaterBase &&)      = default;

    SingleSlaterBase() = delete;

    SingleSlaterBase(MPI_Comm c, CQMemManager &mem, Molecule &mol, BasisSet &basis,
      size_t _nC, bool iCS, Particle p) : 
      WaveFunctionBase(c, mem,mol,basis,_nC,iCS,p), QuantumBase(c, mem,_nC,iCS,p),
      printLevel((MPIRank(c) == 0) ? 2 : 0) { };
      


    // Procedural Functions to be defined in all derived classes
      
    // In essence, all derived classes should be able to:
    //   Form a Fock matrix with the ability to increment
    virtual void formFock(EMPerturbation &, bool increment = false, double xHFX = 1.) = 0;

    // Function to build the modifyOrbitals object which determines which
    // algorithm is used
    virtual void buildModifyOrbitals() = 0;
    virtual void runModifyOrbitals(EMPerturbation&) = 0;

    //   Form an initial Guess (which populates the Fock, Density 
    //   and energy)
    virtual void formGuess(const SingleSlaterOptions&) = 0;

    //   Form the core Hamiltonian
    virtual void formCoreH(EMPerturbation&, bool) = 0;

    //   Save the current state of the wave function
    virtual void saveCurrentState() = 0;

    //   Print various matricies
    virtual void printFock(std::ostream& )     = 0;
    virtual void print1PDMOrtho(std::ostream&) = 0;
    virtual void printGD(std::ostream&)        = 0;
    virtual void printJ(std::ostream&)         = 0;
    virtual void printK(std::ostream&)         = 0;

    virtual void printFockTimings(std::ostream&) = 0;

#ifdef TEST_MOINTSTRANSFORMER
    virtual void MOIntsTransformationTest(EMPerturbation &pert) = 0;
#endif
  }; // class SingleSlaterBase

}; // namespace ChronusQ


