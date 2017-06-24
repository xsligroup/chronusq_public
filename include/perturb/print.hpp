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

#include <perturb.hpp>
#include <util/matout.hpp>
#include <cxxapi/output.hpp>
// #include <physcon.hpp>

namespace ChronusQ {

  /**
   * 
   *  /brief setup formatted line
   *
   */
  void PERTURBFormattedLine(std::ostream &out, std::string s) {
    out << std::setw(45) << "  " + s << std::endl;
  }

  template <typename T>
  void PERTURBFormattedLine(std::ostream &out, std::string s, T v) {
    out << std::setw(45) << "  " + s << v << std::endl;
  }

  template <typename T, typename U>
  void PERTURBFormattedLine(std::ostream &out, std::string s, T v, U u) {
    out << std::setw(45) << "  " + s << v << u << std::endl;
  }

  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::printPTHeader() {

//    auto & ref = dynamic_cast<MCSCF<MatsT,IntsT>&>(*refMCwfn);
    auto & mopart = this->MOPartition;

    std::cout << bannerTop << std::endl;
    std::cout << "PERTURB Settings:"  << std::endl << std::endl;

    std::cout << std::left << std::setprecision(3) << std::scientific;

    std::string job_title = "* Job:  ";
    if(this->NStates == 1) job_title += "Single state PT2";
    else if(PTopts.saFock) {
      job_title += "State-averaged ";
      if (PTopts.extendMS) job_title += "Extended ";
      job_title += "Multi-state PT2";
    } else job_title += "State-specific Multi-state PT2";
    PERTURBFormattedLine(std::cout,job_title);

    std::cout << std::endl;
    PERTURBFormattedLine(std::cout,"* MO Space Partition:");
    PERTURBFormattedLine(std::cout,"  Number of Frozen Core Orbitals:",     mopart.nInact);
    PERTURBFormattedLine(std::cout,"  Number of Frozen Virtual Orbitals:",  mopart.nFVirt);
    PERTURBFormattedLine(std::cout,"  Number of Correlated Orbitals:",      mopart.nCorrO);
    PERTURBFormattedLine(std::cout,"  Number of Correlated Electrons:",     mopart.nCorrE);
    PERTURBFormattedLine(std::cout,"  Number of Orbitals in RAS 1:",      mopart.nActOs[0]);
    PERTURBFormattedLine(std::cout,"  Number of Orbitals in RAS 2:",      mopart.nActOs[1]);
    PERTURBFormattedLine(std::cout,"  Number of Orbitals in RAS 3:",      mopart.nActOs[2]);
    PERTURBFormattedLine(std::cout,"  Maximum Number of Holes in RAS 1:",      mopart.mxHole);
    PERTURBFormattedLine(std::cout,"  Maximum Number of Electrons in RAS 3:",  mopart.mxElec);
    PERTURBFormattedLine(std::cout,"  Total Number of Determinants", this->NDet);

    std::cout << std::endl;
    PERTURBFormattedLine(std::cout,"* Parameters:");
    PERTURBFormattedLine(std::cout,"  Number of Roots Requested:",this->NStates);
    PERTURBFormattedLine(std::cout,"  Number of Determinants in external space:",   SDsize);
    PERTURBFormattedLine(std::cout,"  GMRES Maxmium Number of Iteration:", PTopts.maxIter);
    PERTURBFormattedLine(std::cout,"  GMRES Convergence Threshold:", PTopts.convCrit);

    std::cout << std::endl;
    std::cout << std::endl << bannerTop << std::endl << std::endl;

    this->mointsTF->printMORangesSummary();
    std::cout << std::endl << bannerTop << std::endl << std::endl;


  } // PERTURB::printPTHeader
    
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::printPTFooter() {

    std::cout << std::endl << "PERTURB Results:" << std::endl;
    std::cout << BannerTop << std::endl;

    for (auto i = 0ul; i < this->NStates; i++) {
      std::cout << std::fixed << std::right<< std::setprecision(10);
      std::cout.fill(' ');
      std::cout <<  "State:" << std::setw(4) << SoI_[i]
      << "  Energy (Hartree):" << std::setw(16) << this->StateEnergy[i] <<  std::endl;
    }
    std::cout << BannerTop << std::endl;

  } // PERTURB::printPTFooter


}; // namespace ChronusQ


