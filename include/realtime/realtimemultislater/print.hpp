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

#include <cxxapi/output.hpp>
#include <physcon.hpp>
#include <realtime.hpp>

namespace ChronusQ {

void RTMSFormattedLine(std::ostream &out, std::string s) {
  out << std::setw(38) << "  " + s << std::endl;
}

template <typename T>
void RTMSFormattedLine(std::ostream &out, std::string s, T v) {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <typename T, typename U>
void RTMSFormattedLine(std::ostream &out, std::string s, T v, U u) {
  out << std::setw(38) << "  " + s << v << u << std::endl;
}

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::printRTHeader() {

  // No printing if silent
  if (this->printLevel == 0)
    return;

  std::cout << BannerTop << std::endl;
  std::cout << "Real-Time Propagation Settings:" << std::endl << std::endl;

  std::cout << std::left << std::setprecision(7);
  std::string AUTime = " \u0127 / Eh";

  RTMSFormattedLine(std::cout, "* Simulation Parameters:");

  int nSteps = intScheme.tMax / intScheme.deltaT;
  RTMSFormattedLine(std::cout, "Simulation Time:", intScheme.tMax, AUTime);
  RTMSFormattedLine(std::cout, "Number of Steps:", nSteps);
  RTMSFormattedLine(std::cout, "Step Size:", intScheme.deltaT, AUTime);
  RTMSFormattedLine(std::cout, " ", intScheme.deltaT * FSPerAUTime, " fs");

  std::cout << std::endl;
  RTMSFormattedLine(std::cout, "* Integration Parameters:");

  std::string methString = "Runge Kutta Fourth Order";
  if (intScheme.integrationAlgorithm ==
      RealTimeAlgorithm::RTSymplecticSplitOperator)
    methString = "Symplectic Split Operator";

  RTMSFormattedLine(std::cout, "Electronic Integration:", methString);

  if (pert.fields.size() > 0) {
    std::cout << std::endl;
    RTMSFormattedLine(std::cout, "* Perturbation:\n");

    for (auto &field : pert.fields) {
      std::cout << std::setw(4) << " ";
      std::cout << "Field " << std::distance(&field, &pert.fields[0]) + 1
                << ":  ";
      auto amp = field->getAmp(0);
      if (dynamic_cast<TDDipoleField &>(*field).emFieldTyp == Electric)
        std::cout << "Electric";
      else
        std::cout << "Magnetic";

      std::cout << " ";

      if (amp.size() == 3)
        std::cout << "Dipole";

      std::cout << " Field\n";

      std::cout << std::setw(4) << " ";
      std::cout << std::setw(20) << " * Amplitude (AU)"
                << "{ ";
      for (auto i = 0; i < amp.size(); i++) {
        std::cout << amp[i];
        if (i != amp.size() - 1)
          std::cout << ", ";
      }
      std::cout << " }\n";

      std::cout << std::setw(4) << " ";
      try {
        StepField &env = dynamic_cast<StepField &>(*field->envelope);
        std::cout << std::setw(20) << " * Step Field";
        std::cout << std::setw(9) << "TON = " << std::setw(10) << env.tOn;
        std::cout << "   ";
        std::cout << std::setw(9) << "TOFF = " << std::setw(10) << env.tOff;
        std::cout << std::endl;
      } catch (...) {
      }
    }
  }

  std::cout << std::endl << BannerTop << std::endl;

  std::cout << std::endl << std::fixed << std::right;

  if (this->printLevel == 1) {
    std::cout << std::setprecision(4);
    std::cout << std::setw(11) << "Time (AU)"
              << " ";

    std::cout << std::setprecision(10);
    std::cout << std::setw(16) << "Energy (Eh)"
              << " ";

    std::cout << std::setprecision(8);
    std::cout << std::setw(16) << "Dipole (X)"
              << " ";
    std::cout << std::setw(16) << "Dipole (Y)"
              << " ";
    std::cout << std::setw(16) << "Dipole (Z)"
              << " ";

    std::cout << std::endl << bannerTop << std::endl << std::endl;
  }

  if (this->printLevel == -1)
    this->printLevel = 0;
};

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::printRTStep() {
  if (this->printLevel == 1) {
    printStepSummary();
  } else if (this->printLevel > 1) {
    printStepDetail();
  }
};

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::printStepSummary() {
  std::cout << std::fixed << std::right;
  std::cout << std::setprecision(4);
  std::cout << std::setw(11) << curState.xTime << " ";
  std::cout << std::setprecision(10);
  std::cout << std::setw(16) << this->totalEnergy() << " ";
  std::cout << std::setprecision(8);
  std::cout << std::setw(16) << Dipole[0] << " ";
  std::cout << std::setw(16) << Dipole[1] << " ";
  std::cout << std::setw(16) << Dipole[2] << " ";
  std::cout << std::endl;
};

template <typename MatsT, typename IntsT>
void RealTimeMultiSlater<MatsT, IntsT>::printStepDetail() {
  std::cout << bannerTop << "\n\n";
  std::cout << std::fixed << std::right;
  std::cout << "Step: " << std::setw(7) << curState.iStep << '\n';
  std::cout << std::setprecision(5) << "Time: ";
  std::cout << std::setw(11) << curState.xTime << " (au) | ";
  std::cout << std::setw(11) << curState.xTime * FSPerAUTime << " (fs)\n";
  std::cout << std::setprecision(12) << "Energy: ";
  std::cout << std::setw(24) << this->totalEnergy() << " (Hartree)\n";
  std::cout << std::setprecision(8) << "Dipole: ";
  std::cout << std::setw(16) << Dipole[0] / EBohrPerDebye << " ";
  std::cout << std::setw(16) << Dipole[1] / EBohrPerDebye << " ";
  std::cout << std::setw(16) << Dipole[2] / EBohrPerDebye << " (Debye)";
  std::cout << std::endl;
};

}; // namespace ChronusQ
