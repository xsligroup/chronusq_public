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

#include <cmath>
#include <physcon.hpp>
#include <cxxapi/input.hpp>


namespace ChronusQ {

  // Type of SingleSlater object
  enum RefType {
    isRawRef,  // non-specified
    isRRef,    // RHF/DFT
    isURef,    // UHF/DFT
    isRORef,   // ROHF/DFT
    isTwoCRef, // Two-component
    isFourCRef // Four-component
  };

  // A struct that stores reference information
  struct RefOptions {

    std::string RCflag = "REAL"; // Real or Complex

    RefType refType = isRRef;    // R/U/RO/2c/4c

    bool isKSRef = false;        // HF or DFT
    bool isEPCRef = false;       // NEO-KS or not
    bool isX2CRef = false;       // If user used legacy way to reference X2C

    size_t nC = 1;               // number of component
    bool iCS = true;             // closed shell or not

    std::string funcName;        // DFT functional name

    std::string refstring(RefType reftype) {
      switch (reftype) {
        case isRawRef:
          return "isRawRef";
        case isRRef:
          return "isRRef";
        case isURef:
          return "isURef";
        case isRORef:
          return "isRORef";
        case isTwoCRef:
          return "isTwoCRef";
        case isFourCRef:
          return "isFourCRef";
        default:
          return "";
      }
    } 
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
    bool   useGauXC     = false; ///< Use GauXC as the DFT engine
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
   *  Holds information like convergence criteria, DIIS settings,
   *  max iterations, etc.
   */
  struct SCFControls {

    // Convergence criteria
    double rmsdPConvTol = 1e-7; ///< RSMDP Density convergence criteria
    double maxdPConvTol = 1e-5; ///< MaxDP Density convergence criteria
    double eneConvTol = 1e-5; ///< Energy convergence criteria
    double smallEnergy= 1e-9; ///< Small energy threshold

    // TODO: need to add logic to set this
    // Extrapolation flag for DIIS and damping
    bool doExtrap = true;     ///< Whether to extrapolate Fock matrix

    bool energyOnly = false;  ///< Skip SCF

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
    DIIS_ALG diisAlg = CDIIS; ///< Type of DIIS extrapolation
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
    bool   printContractionTiming = false; ///< Whether to print contraction timing during SCF
    std::string refLongName_;
    std::string refShortName_;

    void parseSection(const InputMap &dict);

  }; // SCFControls struct


  enum class RestartAlgorithm {
    ForwardEuler,
    ModifiedMidpoint,
    ExplicitMagnus2
  };

  enum class PropagatorAlgorithm {
    Diagonalization,
    TaylorExpansion,
    ChebyshevExpansion
  };

  enum class FieldEnvelopeType {
    LinRamp,
    Gaussian,
    Step,
    PlaneWave
  };

  enum class RealTimeAlgorithm {
      Uninitialized,
      RTForwardEuler,
      RTModifiedMidpoint,
      RTExplicitMagnus2,
      RTSymplecticSplitOperator,
      RTRungeKuttaOrderFour
  };

  enum class MSInitialState {
    LinearCombination,
    CustomCI,
  };

  struct TDSCFOptions {

    RealTimeAlgorithm     integrationAlgorithm     = RealTimeAlgorithm::RTModifiedMidpoint; ///< Integration Algorithm
    RealTimeAlgorithm     protIntegrationAlgorithm = RealTimeAlgorithm::Uninitialized;      ///< Protonic Integration Algorithm
    RestartAlgorithm      restartAlgorithm         = RestartAlgorithm::ExplicitMagnus2;     ///< Restart Step
    PropagatorAlgorithm   propagatorAlgorithm      = PropagatorAlgorithm::Diagonalization;  ///< exp(-iF) Algorithm

    double tMax    = 0;  ///< Max simulation time in AU. Upon input, user can specify tMax or maxSteps
    size_t maxSteps= 0;    ///< Max number of steps. Upon input, user can specify tMax or maxSteps
    double deltaT  = 0.01; ///< Time-step in AU

    bool   doMD                = false;
    bool   includeTau          = false;
    size_t totalMDSteps        = 0; 
    size_t rtMaxStepsPerMDStep = 0;

    size_t iRestart  = 50;         ///< Restart MMUT every N steps
    size_t iSave     = 50;         ///< Save progress every N steps
    size_t iPrint    = 50;          ///< Print progress every N steps
    long int restoreFromStep = 0;    ///< Restore propagation from this step

    bool  includeSCFField = true;  ///< Whether to include the SCF field

    size_t rtGaunt = 1; /// < Calculate Gaunt every N steps
    size_t rtGauge = 1; /// < Calculate Gauge every N steps
    size_t rtBreit = 1; /// < Calculate Breit(gaunt and gauge) every N steps
    size_t Rtprintden = 0;
    size_t orbitalPopFreq = 0; ///< Print orbital population every 'orbitalPopFreq' steps during RT propagation

    bool saveOnePDM = false;   ///< Whether to save 1PDM in AO basis to bin file during RT propagation

    void parseSection(const InputMap &dict);
  };

  //struct SCFOptions {
  //};

  struct OrbitalModifierOptions {

    SCFControls scfControls;
    TDSCFOptions tdSCFOptions;
    //PrintControls printControls;

    OrbitalModifierOptions(SCFControls sc, TDSCFOptions tdsc) : scfControls(sc), tdSCFOptions(tdsc) {};

  }; // struct OrbitalModifierOptions



  /**
 *  The Single Slater guess types
 */
  enum SingleSlaterGuessType {
    CoreGuess,
    SADGuess,
    TightGuess,
    RandomGuess,
    ReadBin,
    ReadGaussFCHK,
    ClassicalGuess // For NEO
  };

  struct SingleSlaterGuessOptions {

    SingleSlaterGuessType electronicGuess = SADGuess;
    SingleSlaterGuessType nuclearGuess = TightGuess;
    std::vector<std::vector<size_t>> alphaElectronicMOSwap;
    std::vector<std::vector<size_t>> betaElectronicMOSwap;
    std::vector<std::vector<size_t>> alphaNuclearMOSwap;
    std::vector<std::vector<size_t>> betaNuclearMOSwap;

    void parseSection(const InputMap &dict);

  }; // struct GuessOptions

}; // namespace ChronusQ

