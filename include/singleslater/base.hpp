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
#include <cubegen.hpp>
#include <util/files.hpp>
#include <orbitalmodifieroptions.hpp>
#include <gauxcutils.hpp> 

// #define TEST_MOINTSTRANSFORMER

namespace ChronusQ {



  class SingleSlaterBase;

  struct SingleSlaterOptions {

    HamiltonianOptions hamiltonianOptions;

    RefOptions refOptions;

    IntegrationParam intParam;

    SCFControls scfControls;

    std::shared_ptr<SingleSlaterBase> buildSingleSlater(
        std::ostream &out,
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

    // GauXC utility class                                                      
    std::shared_ptr<GauXCUtils>  gauxcUtils;

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

    // Options for CubeGen
    CubeGenOptions cubeOptsSS;

    // Pair function for SingleSlater MO swap
    std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;

    // Constructors (all defaulted)
    SingleSlaterBase(const SingleSlaterBase &) = default;
    SingleSlaterBase(SingleSlaterBase &&)      = default;

    SingleSlaterBase() = delete;

    SingleSlaterBase(MPI_Comm c, Molecule &mol, BasisSet &basis,
      size_t _nC, bool iCS, Particle p) : 
      WaveFunctionBase(c,mol,basis,_nC,iCS,p), QuantumBase(c,_nC,iCS,p),
      printLevel((MPIRank(c) == 0) ? 2 : 0) { };
      


    // Procedural Functions to be defined in all derived classes
    virtual void initializeSCF() = 0;

    // In essence, all derived classes should be able to:
    //   Form a Fock matrix with the ability to increment
    virtual void formFock(EMPerturbation &, bool increment = false, double xHFX = 1.) = 0;
    // Function to build the orbitalModifier object which determines which
    // algorithm is used
    virtual void buildOrbitalModifierOptions() = 0;
    virtual void runSCF(EMPerturbation&) = 0;

    //   Form an initial Guess (which populates the Fock, Density 
    //   and energy)
    virtual void formGuess(const SingleSlaterOptions&) = 0;

    //   Form the core Hamiltonian
    virtual void formCoreH(EMPerturbation&, bool) = 0;

    //   Save the current state of the wave function
    virtual void saveCurrentState(bool saveMO = true) = 0;

    //   Print various matricies
    virtual void printFock(std::ostream& )     = 0;
    virtual void print1PDMOrtho(std::ostream&) = 0;
    virtual void printGD(std::ostream&)        = 0;
    virtual void printJ(std::ostream&)         = 0;
    virtual void printK(std::ostream&)         = 0;

    virtual void printFockTimings(std::ostream&) = 0;

    // Post-processing functions
    virtual void runCube(std::vector<std::shared_ptr<CubeGen>>, EMPerturbation &, std::string prefix="", std::shared_ptr<Molecule> = nullptr) = 0;

#ifdef TEST_MOINTSTRANSFORMER
    virtual void MOIntsTransformationTest(EMPerturbation &pert) = 0;
#endif
  }; // class SingleSlaterBase

}; // namespace ChronusQ


