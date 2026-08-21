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
#include <cxxapi/options.hpp>
#include <util/mpi.hpp>
#include <cerr.hpp>
#include <physcon.hpp>
#include <gauxc/physcon.hpp>

namespace ChronusQ {
  /**
   *  Check valid keywords in the section.
   */
  std::set<std::string> CQPHYSCON_VALID(const std::map<std::string, std::string>& inputSection) {
    const std::set<std::string> allowedKeywords = {
      "SPEEDOFLIGHT",
      "KGPERAMU",
      "KGPERE",
      "COUPEREL",
      "PLANCKCONST",
      "AVOGCONST",
      "PROTMASSPERE",
      "DEUTMASSPERE",
      "TRITMASSPERE",
      "BOLTZMANNCONST",
      "LEGACY_CONSTANTS"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  } // CQPHYSCON_VALID

  /**
   *  Initialize physical constants from input, legacy values, or defaults.
   *
   *  \param [in] out   Output device for data output
   *  \param [in] input Input file datastructure
   */
  void CQPhysConSetOptions(std::ostream &out, CQInputFile &input) {
    const bool legacyConstants =
      input.containsData("PHYSCON/LEGACY_CONSTANTS") and
      input.getData<bool>("PHYSCON/LEGACY_CONSTANTS");

    if (legacyConstants) {
      for (const auto &keyword : input.getDataInSection("PHYSCON"))
        if (keyword != "LEGACY_CONSTANTS")
          CErr("Cannot set PHYSCON/" + keyword +
            " when PHYSCON/LEGACY_CONSTANTS is true", out);

      // Frozen pre-CODATA-update values.
      SpeedOfLight(137.035999074);
      KgPerAMU(1.66053906892e-27);
      KgPerE(9.1093837139e-31);
      CouPerEl(1.602176565e-19);
      PlanckConst(6.62606957e-34);
      AvogConst(6.02214129e+23);
      ProtMassPerE(1836.152673426);
      DeutMassPerE(3670.482967655);
      TritMassPerE(5496.92153551);

      SpeedOfLight_CM(2.99792458e10);
      SpeedOfLight_SI(SpeedOfLight_CM() / 100);
      HBar_SI(PlanckConst() / (2. * M_PI));
      VacElPermity(8.8541878188e-12);
      JPerHartree(4.35974434e-18);
      EVPerHartree(27.211396132);
      NMPerHartree(45.56335);
      FSPerAUTime(2.4188843265857e-2);
      BohrRadius_SI(HBar_SI() * SpeedOfLight() / (KgPerE() * SpeedOfLight_SI()));
      AngPerBohr(0.52917721092);
      EBohrPerDebye(0.393430307);
      BoltzmannConst(3.166811563e-6);
    } else {
#define CQ_SET_PHYSICAL_CONSTANT(NAME, KEY) \
      if (input.containsData("PHYSCON/" KEY)) \
        NAME(input.getData<double>("PHYSCON/" KEY)); \
      else \
        NAME(NAME##Default)

      CQ_SET_PHYSICAL_CONSTANT(SpeedOfLight,   "SPEEDOFLIGHT");
      CQ_SET_PHYSICAL_CONSTANT(KgPerAMU,       "KGPERAMU");
      CQ_SET_PHYSICAL_CONSTANT(KgPerE,         "KGPERE");
      CQ_SET_PHYSICAL_CONSTANT(CouPerEl,       "COUPEREL");
      CQ_SET_PHYSICAL_CONSTANT(PlanckConst,    "PLANCKCONST");
      CQ_SET_PHYSICAL_CONSTANT(AvogConst,      "AVOGCONST");
      CQ_SET_PHYSICAL_CONSTANT(ProtMassPerE,   "PROTMASSPERE");
      CQ_SET_PHYSICAL_CONSTANT(DeutMassPerE,   "DEUTMASSPERE");
      CQ_SET_PHYSICAL_CONSTANT(TritMassPerE,   "TRITMASSPERE");

      VacElPermity(VacElPermityDefault);
      SpeedOfLight_CM(SpeedOfLight() * 100. * CouPerEl() * CouPerEl() /
        (2. * VacElPermity() * PlanckConst()));
      SpeedOfLight_SI(SpeedOfLight_CM() / 100);
      HBar_SI(PlanckConst() / (2. * M_PI));
      JPerHartree(KgPerE() * SpeedOfLight_SI() * SpeedOfLight_SI() / (SpeedOfLight() * SpeedOfLight()));
      EVPerHartree(JPerHartree() / CouPerEl());
      NMPerHartree(1e9 * PlanckConst() * SpeedOfLight_SI() / JPerHartree());
      FSPerAUTime(1e15 * HBar_SI() / JPerHartree());
      BohrRadius_SI(HBar_SI() * SpeedOfLight() / (KgPerE() * SpeedOfLight_SI()));
      AngPerBohr(1e10 * BohrRadius_SI());
      EBohrPerDebye((1e-21 / SpeedOfLight_SI()) / (CouPerEl() * BohrRadius_SI()));


      double BoltzmannConstDefault = BoltzmannConst_SIDefault / JPerHartree();
      CQ_SET_PHYSICAL_CONSTANT(BoltzmannConst, "BOLTZMANNCONST");

#undef CQ_SET_PHYSICAL_CONSTANT
    }
    
    AUPerAMU(KgPerAMU() / KgPerE());
    CouPerEl_ESU(CouPerEl() * SpeedOfLight_CM() / 10.);
    MassEl_KG(1e4 * JPerHartree() / SpeedOfLight_CM() /
      SpeedOfLight_CM() * SpeedOfLight() * SpeedOfLight());
    HBar(PlanckConst() / 2. / M_PI);
    Rotatory_CGS_Length(1e40 * CouPerEl_ESU() * CouPerEl_ESU() * HBar() *
      AngPerBohr() * 1e7 * 1e-8 /
      (1e3 * MassEl_KG() * SpeedOfLight_CM()));
    Rotatory_CGS_Vel(1e40 * CouPerEl_ESU() * CouPerEl_ESU() * HBar() *
      HBar() * HBar() * 1e21 /
      (MassEl_KG() * MassEl_KG() * SpeedOfLight_CM() * AngPerBohr() *
       JPerHartree() * 1e5));

    // Set external physical constants (dependencies, i.e. GauXC)
    GauXC::SpeedOfLight = SpeedOfLight();
    GauXC::RKB_factor = 1./(4.*SpeedOfLight()*SpeedOfLight());

    // Physical constant printing
    size_t width = 40;
    out << std::endl << "Physical Constants";
    out << ":" << std::endl << BannerTop << std::endl << std::endl;
    out << std::left << std::scientific << std::setprecision(10);
    out << "  " << std::setw(width) << "Speed of Light (a.u.): "
        << std::setw(width) << SpeedOfLight() << std::endl;
    out << "  " << std::setw(width) << "Speed of Light (m/s): "
        << std::setw(width) << SpeedOfLight_SI() << std::endl;
    out << "  " << std::setw(width) << "Angstrom per Bohr: "
        << std::setw(width) << AngPerBohr() << std::endl;
    out << "  " << std::setw(width) << "Kilograms per Atomic Mass Unit: "
        << std::setw(width) << KgPerAMU() << std::endl;
    out << "  " << std::setw(width) << "Electron mass in Kilograms: "
        << std::setw(width) << KgPerE() << std::endl;
    out << "  " << std::setw(width) << "Electron charge (Coulombs): "
        << std::setw(width) << CouPerEl() << std::endl;
    out << "  " << std::setw(width) << "Planck Constant: "
        << std::setw(width) << PlanckConst() << std::endl;
    out << "  " << std::setw(width) << "Avogadro Constant: "
        << std::setw(width) << AvogConst() << std::endl;
    out << "  " << std::setw(width) << "Boltzmann Constant: "
        << std::setw(width) << BoltzmannConst() << std::endl;
    out << "  " << std::setw(width) << "EBohrPerDebye: "
        << std::setw(width) << EBohrPerDebye() << std::endl;
    out << "  " << std::setw(width) << "Electronvolts per Hartree: "
        << std::setw(width) << EVPerHartree() << std::endl;
    out << "  " << std::setw(width) << "Nanometers per Hartree: "
        << std::setw(width) << NMPerHartree() << std::endl;
    out << "  " << std::setw(width) << "Joules per Hartree: "
        << std::setw(width) << JPerHartree() << std::endl;
    out << "  " << std::setw(width) << "Femtoseconds per Atomic Unit: "
        << std::setw(width) << FSPerAUTime() << std::endl;
    out << "  " << std::setw(width) << "Proton/Electron mass ratio: "
        << std::setw(width) << ProtMassPerE() << std::endl;
    out << "  " << std::setw(width) << "Deuteron/Electron mass ratio: "
        << std::setw(width) << DeutMassPerE() << std::endl;
    out << "  " << std::setw(width) << "Triton/Electron mass ratio: "
        << std::setw(width) << TritMassPerE() << std::endl;
    out << "  " << std::setw(width) << "Vacuum Electric Permittivity: "
        << std::setw(width) << VacElPermity() << std::endl;
  }
}; // namespace ChronusQ
