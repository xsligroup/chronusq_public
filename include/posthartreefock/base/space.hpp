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

namespace ChronusQ {

/*
 *  Indix Mapping
 *
 *  'N' -> Negative Energy MO, only for 4C 
 *  'C' -> Frozen Core orbitals 
 *  'I' -> Inactive Core 
 *  'A' -> Active Space for Correlation Calculations 
 *  'S' -> Secondary Virtual 
 *  'V' -> Frozen Virtual Orbitals 
 */

struct CorrelatedMOSpace {
  
  size_t nMO = 0;     /// < Total Number of Molecular Orbitals
  size_t nElecMO = 0; /// < Total Number of Electronic Molecular Orbitals
  size_t nNegMO = 0;  /// < Total Number of Negative Energy Molecular Orbitals
  
  size_t nFCore = 0;  /// < Number of Frozen Core Orbitals    (No rotations)
  size_t nInact = 0;  /// < Number of Uncorrelated Core Orbitals
  size_t nSVirt = 0;  /// < Number of Uncorrelated (Secondary) Virtual Orbtals
  size_t nFVirt = 0;  /// < Number of Frozen Virtual Orbitals (No rotations)

  size_t nCorrO = 0;  /// < Total Number of Correlated Orbitals
  size_t nCorrE = 0;  /// < Total Number of Correlated Electrons

  std::vector<char> orbIndices; /// defining types of orbitals

  // only for 1C
  size_t nCorrEA = 0; /// < Number of Correlated Alpha Electrons
  size_t nCorrEB = 0; /// < Number of Correlated Beta  Electrons

  CorrelatedMOSpace() = default;
  CorrelatedMOSpace(const CorrelatedMOSpace &) = default;
  CorrelatedMOSpace(CorrelatedMOSpace &&)      = default;

}; // struct CorrelatedMOSpace

// ActiveSpaceParameters only include space restrictions.
// Electron occupations are defined as we build the determinant space
struct ActiveSpaceParameters {
  size_t MOOffset  = 0;   /// < orbital index offset
  size_t nOrbitals = 0;   /// < number of active orbitals
  size_t nElectrons= 0;   /// < this holds the reference electron occupation
  size_t minOcc    = 0;   /// < minimum electron occupation
  size_t maxOcc    = 0;   /// < maximum electron occupation
  size_t iDASGroup = 0;   /// < DAS Group ID, used to identify group excitation restrictions
  int    eLimit    = 0;   /// < max number of electrons excited into this DAS, restricted by group excitation.
  int    hLimit    = 0;   /// < max number of electrons excited out of this DAS, restricted by group excitation.
  std::string DASLabel = " ";    /// < a user defined label of the space
  bool intContraction = false;    /// < is this part of an internal contraction?
  SpinType spin = isAlpha;       /// < reserved for one-component case
  std::pair<size_t, size_t> range() const { return {MOOffset, nOrbitals}; };
}; // struct ActiveSpaceParameters
  
inline std::ostream& 
operator <<(std::ostream& os, const std::vector<ActiveSpaceParameters>& actSpaces) {
  size_t counter = 1ul;
  for (auto const & s: actSpaces) {
    os << "    - Space " << std::setw(3) << counter
       << ": NOrb = " << std::setw(3) << s.nOrbitals
       << ", AllowedOcc = " << std::setw(3) << s.minOcc
       << " ~ " << std::setw(3) << s.maxOcc
       << ", DAS Group = " << s.iDASGroup
       << ", MORange = " << std::setw(3) << s.MOOffset + 1
       << " ~ " << std::setw(3) << s.MOOffset + s.nOrbitals
       << std::endl;  
    counter++;
  }
  os << std::left << std::endl;
  return os;
} // printActiveSpaces

} // namespace ChronusQ
