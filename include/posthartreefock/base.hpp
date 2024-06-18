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

  bool StateAverage    = false;
  std::vector<double> SAWeight;

  bool SpinAnalysis = false; // default is do not do Spin analysis
  bool PopulationAnalysis = false; // default is do not do Mulliken analysis
  size_t NosS1 = 0; // number of initial states s1 for oscillator strength
  double * osc_str = nullptr; // matrix to save oscillator strength

  // Options for CubeGen
  CubeGenOptions cubeOptsPostHF;

  // Print Settings
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
  virtual void runCube(std::shared_ptr<CubeGen>) = 0;

  void alloc() {
    this->StateEnergy.clear();
    this->StateEnergy.resize(this->NStates, 0.);
  }

  void dealloc () { }

}; // class PostHartreeFockBase

} // namespace ChronusQ
