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

#include <newperturb.hpp>
#include <util/matout.hpp>
#include <cxxapi/output.hpp>
// #define _DEBUG_PTCATBUILD

namespace ChronusQ {

  /**
   * 
   *  /brief setup formatted line
   *
   */
  void DasPerturbFormattedLine(std::ostream &out, std::string s) {
    out << std::setw(45) << "  " + s << std::endl;
  }

  template <typename T>
  void DasPerturbFormattedLine(std::ostream &out, std::string s, T v) {
    out << std::setw(45) << "  " + s << v << std::endl;
  }

  template <typename T, typename U>
  void DasPerturbFormattedLine(std::ostream &out, std::string s, T v, U u) {
    out << std::setw(45) << "  " + s << v << u << std::endl;
  }

  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::printMRPTHeader() {

//    auto & ref = dynamic_cast<MCSCF<MatsT,IntsT>&>(*refMCwfn);
    auto& ptMOSpace = this->corrSpace; 
    std::cout << bannerTop << std::endl;
    std::cout << "*** MRPT Settings ***"  << std::endl;
    std::cout << std::left << std::setprecision(12) << std::scientific;
    std::cout << std::endl;
    DasPerturbFormattedLine(std::cout,"ENPT2: ", PTopts.ENPT ? "True" : "False");
    DasPerturbFormattedLine(std::cout,"GVVPT2: ", PTopts.GVVPT ? "True" : "False");
    DasPerturbFormattedLine(std::cout,"Virtual Orbitals Selected: ", PTopts.SELECTVIRTUAL);
    DasPerturbFormattedLine(std::cout,"Level shift applied: ", PTopts.LEVELSHIFT);
    DasPerturbFormattedLine(std::cout,"State-Averaged PT2 enable: ", PTopts.STATEAVERAGE ? "True":"False");
    if (PTopts.GVVPT)
      DasPerturbFormattedLine(std::cout,"Secondary States included: ", PTopts.SECONDARYROOTS);
#ifdef CQ_ENABLE_SPARSE
    DasPerturbFormattedLine(std::cout,"Sparse Threshold: ", PTopts.EPS);
#endif
    DasPerturbFormattedLine(std::cout,"External Space Dimension:",  PTFactory_->
      braCategoricalSpace()->nDeterminants());
    DasPerturbFormattedLine(std::cout,"Number of Roots Requested:",this->Target_States_.size());
    std::cout << BannerTop << std::endl;
    DasPerturbFormattedLine(std::cout,"*** PT Space Partition ***");
    DasPerturbFormattedLine(std::cout,"Number of Frozen Core Orbitals:",     ptMOSpace.nFCore);
    DasPerturbFormattedLine(std::cout,"Number of Frozen Virtual Orbitals:",  ptMOSpace.nFVirt);
    DasPerturbFormattedLine(std::cout,"Number of Correlated Orbitals:",      ptMOSpace.nCorrO);
    DasPerturbFormattedLine(std::cout,"Number of Correlated Electrons:",     ptMOSpace.nCorrE);
    std::cout << std::endl << bannerTop << std::endl;

    this->mointsTF->printMORangesSummary();
#ifdef _DEBUG_PTCATBUILD
    std::cout << std::endl << bannerTop << std::endl << std::endl;
    PTFactory_->braCategoricalSpace()->output(std::cout, "Categories Generated in Outer Space");
    std::cout << std::endl << bannerTop << std::endl << std::endl;
#endif

  } // DasPerturb::printMRPTHeader


  template <typename MatsT, typename IntsT>
  void DasPerturb<MatsT,IntsT>::printMRPTFooter() {

    std::cout << std::endl << bannerTop << std::endl << std::endl;
    std::cout << "MRPT2 Energies for Requested Electronic States: " << std::endl;
    std::cout << bannerTop << std::endl;
    for (size_t i = 0ul; i < Target_States_.size(); ++i) {
      std::cout << std::setw(10) << 
      std::fixed << std::setprecision(10) << "State [ " << 
      Target_States_[i]+1 << " ] : " << E2_[i] << std::endl;
      
    }

  } // DasPerturb::printMRPTFooter



}; // namespace ChronusQ
