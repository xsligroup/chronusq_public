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

#include <cerr.hpp>
#include <cmath>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace ChronusQ {

  /**
   *  Physical constants can only be set once 
   *  During testing we need either legacy or new constants so we need a way to reset these
   *  this should only be exposed to testing so its declared but no implementation exists except in testing 
   */
  namespace CQExclusivelyDuringTesting {

    void clearAllPhysicalConstants();

  }

  class PersistentConstant {
  public:
    explicit PersistentConstant(std::string name) : name_(std::move(name)) {
      registry().push_back(this);
    }

    PersistentConstant(const PersistentConstant &) = delete;
    PersistentConstant &operator=(const PersistentConstant &) = delete;

    double operator()(double input) {
      if (value_)
        CErr(std::string("Physical constant ") + name_ +
          " can only be set once.");
      value_ = input;
      return *value_;
    }

    double operator()() const {
      if (not value_)
        CErr(std::string("Physical constant ") + name_ +
          " was used before it was set.");
      return *value_;
    }

  private:
    // Every constant enrols itself here, so that a reset cannot quietly miss
    // one that was added after the reset was written.
    static std::vector<PersistentConstant *> &registry() {
      static std::vector<PersistentConstant *> constants;
      return constants;
    }

    void clear() { value_.reset(); }

    friend void CQExclusivelyDuringTesting::clearAllPhysicalConstants();

    std::string name_;
    std::optional<double> value_;
  };


  // Base physical constants (CODATA, via SciPy 1.18.0).
  // Regenerate with bin/gen_physcon_defaults.py
  inline constexpr double SpeedOfLight_SIDefault   = 299792458.;
  inline constexpr double SpeedOfLightDefault      = 137.03599917700001;
  inline constexpr double KgPerAMUDefault          = 1.6605390689199999e-27;
  inline constexpr double KgPerEDefault            = 9.1093837138999998e-31;
  inline constexpr double CouPerElDefault          = 1.6021766339999999e-19;
  inline constexpr double PlanckConstDefault       = 6.6260701499999998e-34;
  inline constexpr double AvogConstDefault         = 6.0221407599999999e+23;
  inline constexpr double BoltzmannConst_SIDefault = 1.3806490000000001e-23;
  inline constexpr double ProtMassPerEDefault      = 1836.1526734260001;
  inline constexpr double DeutMassPerEDefault      = 3670.4829676549998;
  inline constexpr double TritMassPerEDefault      = 5496.9215355099996;

  inline constexpr double VacElPermityDefault  = CouPerElDefault * CouPerElDefault * SpeedOfLightDefault / (2. * PlanckConstDefault * SpeedOfLight_SIDefault);

  // Persistent physical constants. Each starts empty and accepts one value.
  inline PersistentConstant SpeedOfLight("SpeedOfLight");
  inline PersistentConstant KgPerAMU("KgPerAMU");
  inline PersistentConstant KgPerE("KgPerE");
  inline PersistentConstant CouPerEl("CouPerEl");
  inline PersistentConstant PlanckConst("PlanckConst");
  inline PersistentConstant AvogConst("AvogConst");
  inline PersistentConstant ProtMassPerE("ProtMassPerE");
  inline PersistentConstant DeutMassPerE("DeutMassPerE");
  inline PersistentConstant TritMassPerE("TritMassPerE");
  inline PersistentConstant BoltzmannConst("BoltzmannConst");

  // Derived constants. These are initialized after the base constants and
  // cannot be set directly from the input file.
  inline PersistentConstant SpeedOfLight_CM("SpeedOfLight_CM");
  inline PersistentConstant AUPerAMU("AUPerAMU");
  inline PersistentConstant CouPerEl_ESU("CouPerEl_ESU");
  inline PersistentConstant MassEl_KG("MassEl_KG");
  inline PersistentConstant HBar("HBar");
  inline PersistentConstant Rotatory_CGS_Length("Rotatory_CGS_Length");
  inline PersistentConstant Rotatory_CGS_Vel("Rotatory_CGS_Vel");

  inline PersistentConstant AngPerBohr("AngPerBohr");
  inline PersistentConstant EBohrPerDebye("EBohrPerDebye");
  inline PersistentConstant EVPerHartree("EVPerHartree");
  inline PersistentConstant NMPerHartree("NMPerHartree");
  inline PersistentConstant FSPerAUTime("FSPerAUTime");
  inline PersistentConstant JPerHartree("JPerHartree");
  inline PersistentConstant VacElPermity("VacElPermity");

  inline PersistentConstant HBar_SI("HBar_SI");
  inline PersistentConstant BohrRadius_SI("BohrRadius_SI");
  inline PersistentConstant BoltzmannConst_SI("BoltzmannConst_SI");
  inline PersistentConstant SpeedOfLight_SI("SpeedOfLight_SI");


}; // namespace ChronusQ
