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

namespace ChronusQ {

  enum DetScheme {
    CAS,
    RAS,
    GAS,
    GENERIC_DET
  };

  struct MOSpacePartition {
    // TODO: make variables for cases beyound CAS
    //Parameters for space partition

    DetScheme scheme = CAS;

    size_t nMO=0;     /// < Total Number of Molecular Orbitals
    size_t nElecMO=0; /// < Total Number of Electronic Molecular Orbitals
    size_t nNegMO=0;  /// < Total Number of Negative Energy Molecular Orbitals
    size_t nInact=0;  /// < Number of Uncorrelated Core Orbitals
    size_t nFVirt=0;  /// < Number of Uncorrelated Virtual Orbtals
    size_t nFCore=0;  /// < Number of Frozen Core Orbitals (won't do rotations)
    size_t nDVirt=0;  /// < Number of Virtual being discarded

    size_t nCorrO=0;  /// < Total Number of Correlated Orbitals
    size_t nCorrE=0;  /// < Total Number of Correlated Electrons

    std::vector<char> orbIndices; /// defining types of orbitals

    std::vector<size_t> nActOs;  /// < Number of Correlated Orbitals in each Active Space
    std::vector<size_t> nActEs;  /// < Number of Correlated Electrons in each Active space

    // only for 1C
    size_t nCorrEA=0;                 /// < Number of Correlated Alpha Electrons
    size_t nCorrEB=0;                 /// < Number of Correlated Beta  Electrons
    std::vector<size_t> nActEAs; /// < Number of Correlated Alpha Electrons in each Active space
    std::vector<size_t> nActEBs; /// < Number of Correlated Beta  Electrons in each Active space

    // for RAS
    size_t mxHole=0;  /// < Maximum number of holes in RAS 1 space
    size_t mxElec=0;  /// < Maximum number of electrons in RAS 3 space
    std::vector<int> fCat;  /// < Category offset for RAS string

    MOSpacePartition() = default;
    MOSpacePartition(const MOSpacePartition &) = default;
    MOSpacePartition(MOSpacePartition &&)      = default;
  
  };

  // helper function to parse orbital index
  void set_orbital_index(std::vector<char> &orbindex, std::string inputstring, char C);
  void fill_default_index(std::vector<char> &orbindex, size_t iter, char C, size_t N);
  void print_orbIndices(std::vector<char> orbindex);


  /**
   *  \brief The MCWaveFunctionBase class. The abstraction of information
   *  relating to the MCWaveFunction class which are independent of storage
   *  type.
   *
   *
   *  See WaveFunction for further docs.
   */
  class MCWaveFunctionBase : public ManyBodyWavefunctionBase {

  public:

    typedef std::vector<std::vector<int>> int_matrix;

    SafeFile savFile;    ///< Data File, for restart
    MPI_Comm comm;

    MOSpacePartition MOPartition;

    size_t NDet;         ///  < Number of Determinants
    size_t NStates = 1;  ///  < Number of States
    bool FourCompNoPair = true; /// default as true if 4C

    // TODO: DetString now only works for CASCI String
    std::shared_ptr<DetStringManager> detStr     = nullptr;
    std::shared_ptr<DetStringManager> detStrBeta = nullptr; // only for 1C

    double InactEnergy;
    std::vector<double> StateEnergy;
    // Storage for Field-Nuclear dipole interactions
    // This is additive to the diagonal in CI theory, so it can be 
    // simply added to the total state energies on convergence
    double EFieldNuc = 0.0;

    // Flags for avoiding redundant integral transformations
    // If the applied field changed, we'll need to clear the AO Cache
    // for transforming the one particle integrals
    bool field_changed   = true;
    std::array<double, 3> old_dip_field;
    // Unused for now since the only place orbitals are currently changing
    // in MCSCF is during orbital optimization and the cache is just cleared
    // at that time
    bool orbital_changed = false;

    bool StateAverage    = false;
    std::vector<double> SAWeight;
    
    bool PopulationAnalysis = false; // default is do not do Mulliken analysis
    bool SpinAnalysis = false; // default is do not do Spin analysis
    size_t NosS1 = 0; // number of initial states s1 for oscillator strength
    double * osc_str = nullptr; // matrix to save oscillator strength
    bool multipoleMoment = false; // default is do not compute multipole moments

    // Length gauge electric multipoles
    std::vector<cart_t> elecDipoles;        ///< Electric Dipole in the length gauge
    std::vector<cartmat_t> elecQuadrupoles; ///< Electric Quadrupole in the length gauge
    std::vector<cartrk3_t> elecOctupoles;   ///< Electric Octupole in the length gauge

    bool readCI = false; ///< Read CI vectors and state energies from rstfiles

    // Perturbation
    EMPerturbation mcscfPert;

    // Options for CubeGen
    CubeGenOptions cubeOptsMC;

    // Print Settings
    size_t printMOCoeffs = 0;
    size_t printRDMs = 0;
    double rdmCut = 0.10;

    MCWaveFunctionBase() = delete;
    MCWaveFunctionBase(const MCWaveFunctionBase &) = default;
    MCWaveFunctionBase(MCWaveFunctionBase &&)      = default;

    /**
     *  MCWaveFunctionBase Constructor. Constructs a WaveFunctionBase object
     *
     *  \param [in] NS Number of States constructed by MCWaveFunction
     */
    MCWaveFunctionBase(MPI_Comm c, size_t NS):
      comm(c), NStates(NS) {

      alloc();

    }; // MCWaveFunctionBase Constructor.

    ~MCWaveFunctionBase() { dealloc(); };

    void partitionMOSpace(std::vector<size_t>, size_t);
    void turnOnStateAverage(const std::vector<double> &);
    std::vector<int> genfCat(std::vector<std::vector<size_t>>,
      std::vector<std::vector<size_t>>, size_t, size_t);
    std::vector<int> genfCat(std::vector<std::vector<size_t>>, size_t);

    // Virtural Run function
    virtual void run(EMPerturbation &) = 0;
    virtual void setMORanges() = 0;
    virtual WaveFunctionBase & referenceWaveFunction() = 0;
    virtual void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>&, SpinType) = 0;
    void setActiveSpaceAndReOrder();

    // Post-processing functions
    virtual void runCube(std::vector<std::shared_ptr<CubeGen>>, EMPerturbation &) = 0;

    void alloc() {
      this->StateEnergy.clear();
      this->StateEnergy.resize(this->NStates, 0.);
    }

    void dealloc () { }

  }; // class MCWaveFunctionBase

}; // namespace ChronusQ
