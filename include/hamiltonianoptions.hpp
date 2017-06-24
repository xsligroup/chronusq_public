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
#include <basisset.hpp>
#include <particleintegrals.hpp>


namespace ChronusQ {

  enum KineticBalance {RKBPauli, RKBSpinor, UKBScalar};

  enum class TYPE_4C {All, SpinFree, SpinDependent};

  enum class APPROXIMATION_TYPE_4C {None, ThreeCenter, TwoCenter, OneCenter, AtomicMeanField};

  enum class X2C_TYPE {OFF, ONEE, TWOE, FOCK};

  /**
   * Type of screened nuclear spin–orbit approximation
   * BOETTGER:    Scaling factors proposed by Boettger, Phys. Rev. B 62, 7809 (2000)
   * DC,DCB:      Scaling factors proposed by Ehrman et al., J. Chem. Theory Comput. 19, 5785 (2023)
   * ROW_DEP_DCB: Row-dependent version of DCB
   */
  enum class SNSO_TYPE {BOETTGER, DC, DCB, ROW_DEP_DCB};

  struct ATOMIC_X2C_TYPE {
    bool isolateAtom;  ///< If atomic OEI feel only the basis origin nuclei potential
    bool diagonalOnly; ///< If only diagonal blocks of Hamiltonian are X2C corrected

    std::string toString() const {
      if (isolateAtom) {
        if (diagonalOnly)
          return "ALH";
        else
          return "ALU";
      } else {
        if (diagonalOnly)
          return "DLH";
        else
          return "DLU";
      }
    }
  };

  /**
   *  The particle types to evaluate integrals
   */
  struct Particle {
    double charge = -1.0; // particle charge
    double mass = 1.0;    // particle mass
  };


  struct HamiltonianOptions {

    // Number of components
//    size_t nComponents = 1; // 1, 2, or 4

    // Time-reversal symmetry
    bool KramersSymmetry = false; // Whether or not to use the Kramers' symmetry

    // Integral Options
    BASIS_FUNCTION_TYPE basisType = REAL_GTO; //GTO or GIAO
    bool finiteWidthNuc = false; // Use finite nuclei in integral evaluations
    bool Libcint = false; // Use Libcint library instead of Libint
    Particle particle; // Particle type

    // One-Component Options
    bool PerturbativeScalarRelativity = false; // Add perturbative scalar relativity
    bool PerturbativeSpinOrbit = false; // Add perturbative spin-orbit

    // Two-Component Options
    X2C_TYPE x2cType = X2C_TYPE::OFF; //Type of X2C
    bool OneEScalarRelativity = true; //scalar relativity
    bool OneESpinOrbit = true; //spin-orbit relativity
    bool SNSO = true; // Use Boetteger factor to scale one-electron spin-orbit
    SNSO_TYPE snsoType = SNSO_TYPE::BOETTGER; // Type of screened nuclear spin–orbit
    bool AtomicMeanField = false; // Use atomic mean field two-electron spin-orbit
    bool AtomicX2C = false; // Use atomic X2C
    ATOMIC_X2C_TYPE AtomicX2CType; // The type of atomic X2C

    // Four-Component Options
    KineticBalance kineticBalance = RKBPauli; // Choose the kinetic-balance condition; currently, only RKBPauli is implemented
    bool BareCoulomb = true; // Do bare Coulomb only in Restricted-Kinetic balance (RKB)
    bool DiracCoulomb = true; // Dirac-Coulomb without SSSS
    bool DiracCoulombSSSS = true; // SSSS to Dirac-Coulomb
    bool Gaunt = false; // Gaunt
    bool Gauge = false; // Gauge
    TYPE_4C DiracCoulombType = TYPE_4C::All; // Type of Dirac-Coulomb - All, Spin Free, Spin Dependent
    TYPE_4C SSSSType = TYPE_4C::All;         // Type of SSSS - All, Spin Free, Spin Dependent
    TYPE_4C GauntType = TYPE_4C::All;        // Type of Gaunt - All, Spin Free, Spin Dependent
    TYPE_4C GaugeType = TYPE_4C::All;        // Type of Gauge - All, Spin Free, Spin Dependent
                               
    APPROXIMATION_TYPE_4C DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::None; // Type of Dirac-Coulomb approximation - 3-Center, 2-Center, 1-Center, Atomic Mean Field
    APPROXIMATION_TYPE_4C SSSSApproximationType = APPROXIMATION_TYPE_4C::None;         // Type of SSSS approximation - 3-Center, 2-Center, 1-Center, Atomic Mean Field
    APPROXIMATION_TYPE_4C GauntApproximationType = APPROXIMATION_TYPE_4C::None;        // Type of Gaunt approximation - 3-Center, 2-Center, 1-Center, Atomic Mean Field
    APPROXIMATION_TYPE_4C GaugeApproximationType = APPROXIMATION_TYPE_4C::None;        // Type of Gauge approximation - 3-Center, 2-Center, 1-Center, Atomic Mean Field

    bool updateGaunt = false;      //Default False. True if open Gaunt and need to update. False if do not update gaunt term at current step 
    bool updateGauge = false;      //Default False. True if open Gauge and need to update. False if do not update gauge term at current step  

    bool includeTau = false;

  }; // struct HamiltonianOptions

  std::ostream& operator<<(std::ostream&, const HamiltonianOptions&);

}; // namespace ChronusQ

