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

#include <util/math.hpp>
#include <detstringmanager.hpp>
#include <wavefunction/base.hpp>
#include <manybodywavefunction.hpp>
#include <manybodywavefunction/base.hpp>
#include <neworbitalrotation.hpp>

namespace ChronusQ {

  enum CIDiagonalizationAlgorithm {
    CI_FULL_MATRIX,
    CI_DAVIDSON,
    CI_GPLHR,
    SKIP,
  }; // struct CIDiagonalizationAlgorithm
 
  // Settings 
  struct MCSCFSettings {
     
     // CISettings
     CIDiagonalizationAlgorithm ciAlg = CI_FULL_MATRIX; 
     bool ciAlgUserSet = false;

     // Number of roots solved for in MCSCF
     size_t NStates = 1;
     
     // for davidson and gplhr
     size_t maxCIIter        = 128;
     double ciVectorConv     = 1.0e-6;    
     size_t maxDavidsonSpace = 50;
     size_t nDavidsonGuess   = 3;
     std::vector<std::pair<double, size_t>> energyRefs;

     // For Natural Orbitals 
     size_t NatOrbs = 0; // index for which root to for natural orbitals for
                         // Note this is saved as 1 indexed (i.e., NatOrbs=1 will
                         // form natural orbitals for the lowest energy root)
                         // Default value of 0 causes no Natural orbital formation
     // By default we want to re-express the CI vectors
     // in the new natural orbital basis
     bool NatOrbRediag = true;

     // Post MCSCF Analysis options (control flow is handled by MCSCF, 
     // the actual calculation is handled by MCWaveFunction)
     bool PopulationAnalysis = false; // default is do not do Mulliken analysis
     bool SpinAnalysis = false; // default is do not do Spin analysis
     size_t NosS1 = 0; // number of initial states s1 for oscillator strength
     bool multipoleMoment = false; // default is do not compute multipole moments


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

  /**
   *  \brief The MCSCFBase class. The abstraction of information
   *  relating to the MCSCF method 
   */
  class MCSCFBase {

  public:

    std::shared_ptr<MCSCFSettings> settings;

    SafeFile savFile;    ///< Data File, for restart
    MPI_Comm comm;

    // These are all inherent to MCWaveFunctionBase but we'll grab copies here
    // for convenience through the calculation
    size_t NDet;         ///  < Number of Determinants
    size_t NStates = 1;  ///  < Number of States
    std::shared_ptr<std::vector<double>> StateEnergy;

    // Specifics of the SCF in MCSCF calculation
    bool StateAverage    = false;
    std::vector<double> SAWeight;

    // Perturbation
    EMPerturbation mcscfPert;

    MCSCFBase() = delete;
    MCSCFBase(const MCSCFBase&) = delete;
    MCSCFBase(MCSCFBase&&)      = delete;

    /**
     *  MCWaveFunctionBase Constructor. Constructs a WaveFunctionBase object
     *
     *  \param [in] NS Number of States constructed by MCWaveFunction
     */
    MCSCFBase(MPI_Comm c, size_t NS):
      comm(c), NStates(NS) 
      {}; // MCWaveFunctionBase Constructor.

    ~MCSCFBase() { dealloc(); };

    void turnOnStateAverage(const std::vector<double> &);

    // Virtural Run function
    virtual void run(EMPerturbation &) = 0;
    virtual void runCube(std::vector<std::shared_ptr<CubeGen>>) = 0;

    // Post-processing functions
    //virtual void runCube(std::vector<std::shared_ptr<CubeGen>>) = 0;

    void dealloc () { }

  }; // class MCWaveFunctionBase

}; // namespace ChronusQ
