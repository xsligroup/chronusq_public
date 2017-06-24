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

#include <cmath>

namespace ChronusQ {

  // Physical Constants 
  // XXX: Could use some citations / more accurate values here
  constexpr double AngPerBohr    = 0.52917721092;
  constexpr double KgPerAMU      = 1.6650538921e-27;
  constexpr double CouPerEl      = 1.602176565e-19;
  constexpr double PlanckConst   = 6.62606957e-34;
  constexpr double AvogConst     = 6.02214129e+23;
  constexpr double EBohrPerDebye = 0.393430307;
  constexpr double EVPerHartree  = 27.211396132;
  constexpr double NMPerHartree  = 45.56335;
  constexpr double SpeedOfLight  = 137.035999074;
  //constexpr double SpeedOfLight  = 137.03599967994; // PySCF
  constexpr double FSPerAUTime   = 2.4188843265857e-2;
  constexpr double JPerHartree   = 4.35974434e-18;
  constexpr double ProtMassPerE  = 1836.15267343;
  constexpr double AUPerAMU      = 1822.888486217313;

  // Things in odd unit systems / derived
  constexpr double SpeedOfLight_CM = 2.99792458e+10; 
  constexpr double CouPerEl_ESU = 
    CouPerEl * SpeedOfLight_CM / 10.;
  constexpr double MassEl_KG = 
    1e4 * JPerHartree / SpeedOfLight_CM / SpeedOfLight_CM *
    SpeedOfLight * SpeedOfLight;
  constexpr double HBar = PlanckConst / 2. / M_PI;


  // Unit conversions
  constexpr double Rotatory_CGS_Length = 
    1e40 * CouPerEl_ESU * CouPerEl_ESU * HBar * AngPerBohr *
    1e7 * 1e-8 / (1e3 * MassEl_KG * SpeedOfLight_CM);   

  constexpr double Rotatory_CGS_Vel = 
    1e40 * CouPerEl_ESU * CouPerEl_ESU * HBar*HBar*HBar * 
    1e21 / ( MassEl_KG * MassEl_KG * SpeedOfLight_CM * 
        AngPerBohr * JPerHartree * 1e5 );
};

