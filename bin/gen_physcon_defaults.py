#!/usr/bin/env python3
"""Print the default physical constant block for include/physcon.hpp.

base defaults values from CODATA (via scipy)
derived default values are constexpr expressions of the base defaults
"""
import scipy
import scipy.constants as sc

print(f"""  // Base physical constants (CODATA, via SciPy {scipy.__version__}).
  // Regenerate with bin/gen_physcon_defaults.py
  inline constexpr double SpeedOfLight_SIDefault   = {sc.c:.17g}.;
  inline constexpr double SpeedOfLightDefault      = {sc.value('inverse fine-structure constant'):.17g};
  inline constexpr double KgPerAMUDefault          = {sc.value('atomic mass constant'):.17g};
  inline constexpr double KgPerEDefault            = {sc.m_e:.17g};
  inline constexpr double CouPerElDefault          = {sc.e:.17g};
  inline constexpr double PlanckConstDefault       = {sc.h:.17g};
  inline constexpr double AvogConstDefault         = {sc.N_A:.17g};
  inline constexpr double BoltzmannConst_SIDefault = {sc.k:.17g};
  inline constexpr double ProtMassPerEDefault      = {sc.value('proton-electron mass ratio'):.17g};
  inline constexpr double DeutMassPerEDefault      = {sc.value('deuteron-electron mass ratio'):.17g};
  inline constexpr double TritMassPerEDefault      = {sc.value('triton-electron mass ratio'):.17g};
  
  inline constexpr double VacElPermityDefault  = CouPerElDefault * CouPerElDefault * SpeedOfLightDefault / (2. * PlanckConstDefault * SpeedOfLight_SIDefault);
  """)
