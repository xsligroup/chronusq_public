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

#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::printCorrMOSpace() {
      
    auto & corrS = this->corrSpace;
  
    std::cout << std::left << std::endl;
    FormattedLine(std::cout,"* Correlated MO Space Partition:");
    FormattedLine(std::cout,"  Number of Electronic MOs:", corrS.nElecMO);
    if (ref_->nC == 4)
      FormattedLine(std::cout,"  Number of Negative MOs:", corrS.nNegMO);
    
    std::cout << std::endl;
    FormattedLine(std::cout,"  Number of Frozen Core Orbitals:",       corrS.nFCore);
    FormattedLine(std::cout,"  Number of Inactive Core Orbitals:",     corrS.nInact);
    FormattedLine(std::cout,"  Number of Correlated Orbitals:",        corrS.nCorrO);
    FormattedLine(std::cout,"  Number of Secondary Virtual Orbitals:", corrS.nSVirt);
    FormattedLine(std::cout,"  Number of Frozen Virtual Orbitals:",    corrS.nFVirt);
    
    std::cout << std::endl;
    FormattedLine(std::cout,"  Number of Electrons in Full Space:",  ref_->nO);
    FormattedLine(std::cout,"  Number of Correlated Electrons:",     corrS.nCorrE);
    if (ref_->nC == 1) {
      FormattedLine(std::cout,"  Number of Correlated Alpha Electrons:", corrS.nCorrEA);
      FormattedLine(std::cout,"  Number of Correlated Beta Electrons:",  corrS.nCorrEB);
    }
    
    std::cout << std::endl;
    std::cout << std::right << std::endl;
    this->mointsTF->printMORangesSummary();
    std::cout << std::left  << std::endl;
  
  }; // PostHartreeFock::printMOSpacePatition
  
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::print1RDMs() {

    if( printRDMs==0 ) return;

    std::cout << std::endl << "PostHartreeFock 1RDMs:" << std::endl;
    std::cout << bannerTop << std::endl;

    size_t nCorrO = this->corrSpace.nCorrO;
    size_t nS     = this->NStates;
    size_t nRDMEle= 0; // used for line breaker

    //Print full 1RDM
    if( printRDMs==1 ){
      std::cout << "Printing full 1-RDM" << std::endl;
      for (auto i = 0ul; i < nS; i++) {

        prettyPrintSmart(std::cout, "State " + std::to_string(i),
          this->oneRDM[i]->pointer(), nCorrO, nCorrO, nCorrO);

	if (this->SpinAnalysis)
          this->spinAnalysis(i);
      }
    } else if( printRDMs==2 ){
      std::cout.precision(2);
      std::cout << "Printing large (>" << std::fixed << rdmCut << ") diagonal elements of real 1-RDM" << std::endl;

      for (auto i = 0ul; i < nS; i++) {

        std::cout << std::endl << "State " << i+1 << ": ";
        nRDMEle = 0;

        for( auto ipp = 0ul; ipp < nCorrO; ipp++){

          if( std::real(this->oneRDM[i]->pointer()[ipp*nCorrO+ipp]) > rdmCut ){
            if( nRDMEle > 6 ){
              std::cout << std::endl << "         ";
              nRDMEle = 0;
            }
            std::cout << std::setw(3) << ipp+1 << "(" << std::real(this->oneRDM[i]->pointer()[ipp*nCorrO+ipp]) << ")" << " ";
            nRDMEle++;
          }
        }

        std::cout << std::endl;
	
	if (this->SpinAnalysis)
          this->spinAnalysis(i);

      }
    }

    std::cout << bannerTop << std::endl;

  }; //PostHartreeFock::print1RDMs

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::printMOInfo(std::ostream& out, 
                      size_t printMOLevel) {
    
    if (not printMOLevel)
      printMOLevel = printMOCoeffs;
    ref_->WaveFunction<MatsT,IntsT>::printMOInfo(out, printMOLevel);

  } //PostHartreeFock::printMOInfo



}; // namespace ChronusQ
