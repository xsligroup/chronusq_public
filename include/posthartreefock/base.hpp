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
#include <posthartreefock/base/space.hpp>

namespace ChronusQ {

/**
 *  \brief The PostHartreeFockBase class. The abstraction of information
 *  relating to the MCWaveFunction class which are independent of storage
 *  type.
 *
 *
 *  See WaveFunction for further docs.
 */
class PostHartreeFockBase {

public:

  SafeFile savFile;    ///< Data File, for restart
  MPI_Comm comm;

  std::shared_ptr<WaveFunctionBase> wfnRef_;

  CorrelatedMOSpace corrSpace;
  bool FourCompNoPair = true;
  
  size_t NStates = 1;  ///  < Number of States
  
  double coreEnergy;
  std::vector<double> StateEnergy;

  std::vector<cart_t> SExpectState;
  std::vector<double> SSqState;
  std::vector<cart_t> LExpectState;
  std::vector<double> LSqState;
  std::vector<cart_t> JExpectState;
  std::vector<double> JSqState;
  std::vector<double> SLState;

  bool StateAverage    = false;
  std::vector<double> SAWeight;

  bool SpinAndAngularAnalysis = false; // default is do not do spin/angular analysis
  bool PopulationAnalysis = false; // default is do not do Mulliken analysis
  bool osc_str = false; //default is do not do any oscillator strength
  size_t osc_str_order = 0; //default is to do oscillator strengths within dipole approximation
  size_t NosS1 = 1; // number of initial states s1 for oscillator strength
  std::vector<double> osc_str_array; // matrix to save oscillator strength

  bool printTransDipole = false;
  bool printDipole = false;
  bool saveOnePDMS = false;
  std::vector<size_t> saveOnePDM_states;

  // Options for CubeGen
  CubeGenOptions cubeOptsPostHF;

  // Print Settings
  bool printDetailedCICoeffs = false;
  size_t printMOCoeffs = 0;
  size_t printRDMs     = 0;
  double rdmCut        = 0.10;

  PostHartreeFockBase()                           = delete;
  PostHartreeFockBase(const PostHartreeFockBase &) = default;
  PostHartreeFockBase(PostHartreeFockBase &&)      = default;

  /**
   *  PostHartreeFockBase Constructor. Constructs a WaveFunctionBase object
   *
   *  \param [in] NS Number of States constructed by MCWaveFunction
   */
  PostHartreeFockBase(MPI_Comm c,
    std::shared_ptr<WaveFunctionBase> wfnRef, size_t NS):
    comm(c), wfnRef_(wfnRef), NStates(NS) {

    alloc();

  }; // PostHartreeFockBase Constructor.

  ~PostHartreeFockBase() { dealloc(); };

  void setupCorrelatedMOSpace(size_t nCorrO, size_t nCorrE, size_t nFCore, size_t nFVirt);
  void turnOnStateAverage(const std::vector<double> &);
  void setCorrelatedSpaceAndReOrder();

  // Virtural Run function
  virtual void run(EMPerturbation &) = 0;
  virtual void setMORanges() = 0;
  virtual void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>&, SpinType) = 0;

  // Post-processing functions
  virtual void runCube(std::vector<std::shared_ptr<CubeGen>>) = 0;

  void alloc() {
    this->StateEnergy.clear();
    this->StateEnergy.resize(this->NStates, 0.);

    this->SExpectState.clear();
    this->SExpectState.resize(this->NStates, {0., 0., 0.});
    this->SSqState.clear();
    this->SSqState.resize(this->NStates, 0.);
    this->LExpectState.clear();
    this->LExpectState.resize(this->NStates, {0., 0., 0.});
    this->LSqState.clear();
    this->LSqState.resize(this->NStates, 0.);
    this->JExpectState.clear();
    this->JExpectState.resize(this->NStates, {0., 0., 0.});
    this->JSqState.clear();
    this->JSqState.resize(this->NStates, 0.);
    this->SLState.clear();
    this->SLState.resize(this->NStates, 0.);
  }

  void dealloc () { }

}; // class PostHartreeFockBase

} // namespace ChronusQ
