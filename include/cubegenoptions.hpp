/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

namespace ChronusQ {

  // The standard options for obitals a user may want to plot
  // Currently unused 
  enum class MO_CLASSES {
    ALL,
    CUSTOM
  };

  // Standard options for which however many roots we may want to calculate
  // in an MCSCF calculation
  enum class CI_CUBE_ROOT_CLASSES {
    GS,
    ALL,
    AVERAGE,
    CUSTOM
  };

  struct CubeGenOptions {
 
    // Commonly-used options
    bool denCube = false; // Turn on density evaluation
    bool orbCube = false; // Turn on orbital evaluation

    // Options related to orbitals
    bool MagnitudeAndPhase = false; // MOs printed as magnitude/phase instead of real/imaginary
    MO_CLASSES whichMO; // Standard orbital choice options
    std::vector<size_t> custom_orb_request; // User-defined choice of orbitals

    // Options for MCSCF Calculations
    CI_CUBE_ROOT_CLASSES whichCIRoots = CI_CUBE_ROOT_CLASSES::GS; // By default only the lowest energy root generates cubes
    std::vector<size_t> custom_root_request; // User-defined set of roots to calculate

    // Miscellaneous options
    std::string cubeFileName = ""; //Preceding string for cube file names

    /**
     * @brief add an orbital in a custom list of MO's to generate
     * 
     * NOTE: input is 1 indexed, but internal storage is 0 indexed
    */
    void addMOtoList(size_t mo)
    {
      custom_orb_request.push_back(mo-1);
    }

    /**
     * @brief add a root in a custom list of roots to generate
     * 
     * NOTE: input is 1 indexed, but internal storage is 0 indexed
    */
    void addRoottoList(size_t root)
    {
      custom_root_request.push_back(root-1);
    }

  };

}

