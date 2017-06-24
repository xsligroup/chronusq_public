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
#include <quantum/base.hpp>

namespace ChronusQ {


  // Type of spin
  enum SpinType {
    isAlpha= 0,  // default orbital manifold
    isBeta = 1,  // beta manifold for unrestricted
  };

  /**
   *  \brief The WaveFunctionBase class. The abstraction of information
   *  relating to the WaveFunction class which are independent of storage
   *  type.
   *
   *  Specializes QuantumBase interface.
   *
   *  See WaveFunction for further docs.
   */ 
  class WaveFunctionBase : virtual public QuantumBase {

  public:
    // Member data

    Molecule &molecule_; ///< A reference of the Molecule
    BasisSet &basisSet_; ///< BasisSet for the GTO basis defintion

    size_t nO;  ///< Total number of occupied orbitals
    size_t nV;  ///< Total number of virtual orbitals
    size_t nOA; ///< Number of occupied alpha orbitals (nC == 1)
    size_t nOB; ///< Number of occupied beta orbitals  (nC == 1)
    size_t nVA; ///< Number of virtual alpha orbitals  (nC == 1)
    size_t nVB; ///< Number of virtual beta orbitals   (nC == 1)

    
    // Disable default constructor
    WaveFunctionBase() = delete;

    // Default copy and move constructors
    WaveFunctionBase(const WaveFunctionBase &) = default;
    WaveFunctionBase(WaveFunctionBase &&)      = default;



    /**
     *  WaveFunctionBase Constructor. Constructs a WaveFunctionBase object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] _nC  Number of spin components (1 and 2 are supported)
     *  \param [in] iCS  Whether or not to treat as closed shell
     */ 
    WaveFunctionBase(MPI_Comm c, Molecule &mol, BasisSet &basis,
      size_t _nC, bool iCS, Particle p) : 
      QuantumBase(c,_nC,iCS,p), 
      molecule_(mol), basisSet_(basis) 
      { }; // WaveFunctionBase ctor 

    // Member Functions
    Molecule& molecule() { return molecule_; }
    BasisSet& basisSet() { return basisSet_; }

    size_t nAlphaOrbital() const { return nOA + nVA; }
    size_t nBetaOrbital() const { return nOB + nVB; }
    size_t nOrbital() const { return nO + nV;}
    // Print Functions
    virtual void printMO(std::ostream&)  = 0;
    virtual void printEPS(std::ostream&) = 0;
    virtual void printMOInfo(std::ostream&, size_t a = 0) = 0;

    // MO swap functions
    virtual void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>&, SpinType sp) = 0;

    // Get geometric gradients
    virtual std::vector<double> getGrad(EMPerturbation&, bool equil,
      bool saveInts, double xHFX = 1.) = 0;

  }; // class WaveFunctionBase

}; // namespace ChronusQ

