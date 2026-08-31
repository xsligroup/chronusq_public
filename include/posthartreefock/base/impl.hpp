 /* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
 *  
 *  This program is free software; you ca redistribute it and/or modify
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

#include <posthartreefock/base.hpp>

namespace ChronusQ {

  void PostHartreeFockBase::setupCorrelatedMOSpace(
    size_t nCorrO, size_t nCorrE, size_t nFCore, size_t nFVirt) {
  
    if (wfnRef_->nC == 4 and not this->FourCompNoPair)
      CErr("4C without no-pair approximation is not implemented.");

    if (wfnRef_->nC == 1) {
      corrSpace.nElecMO = wfnRef_->nAlphaOrbital();
    } else { 
      corrSpace.nElecMO = wfnRef_->nV + wfnRef_->nO;
    } 
    
    
    if (wfnRef_->nC == 4) {
      corrSpace.nMO = 2 * corrSpace.nElecMO;
      corrSpace.nNegMO = corrSpace.nElecMO;
    } else {
      corrSpace.nMO = corrSpace.nElecMO; 
    }

    if( (nCorrO + nFCore > corrSpace.nElecMO) or nCorrO == 0 
        or nCorrE + nFCore > wfnRef_->nO or nCorrE == 0)  
      CErr("Invalid correlated space, please modify the input!");
    
    if (wfnRef_->nC == 1) {
      if (wfnRef_->iCS) {
        corrSpace.nCorrEA = nCorrE / 2;
        corrSpace.nCorrEB = corrSpace.nCorrEA;  
      } else {
        int singleE   = wfnRef_->nOA - wfnRef_->nOB;
        corrSpace.nCorrEA = (nCorrE + singleE) / 2;
        corrSpace.nCorrEB = nCorrE - corrSpace.nCorrEA;  
      }

      if (nCorrO < corrSpace.nCorrEA) 
        CErr("The specified number of correlated orbitals can not hold all the correlated electrons"); 
    
      corrSpace.nInact = wfnRef_->nOA - corrSpace.nCorrEA - nFCore;
    } else { // other than 1C
  
      if (nCorrO < nCorrE)
        CErr("The specified number of correlated orbitals can not hold all the correlated electrons"); 
      
      corrSpace.nInact = wfnRef_->nO - nCorrE - nFCore;
    }    
    
    corrSpace.nFCore = nFCore;
    corrSpace.nFVirt = nFVirt;
    corrSpace.nCorrO = nCorrO;
    corrSpace.nCorrE = nCorrE;
    
    corrSpace.nSVirt = corrSpace.nElecMO - nFCore - corrSpace.nInact - nCorrO - nFVirt;
    // check if nSVirt overflows 
    if( corrSpace.nElecMO < (nCorrO + corrSpace.nInact + corrSpace.nSVirt + nFCore + nFVirt) )
      CErr("Invalid correlated space, please modify the input!");
    
    // initialize orbital index
    setCorrelatedSpaceAndReOrder();   
    
    // initialize indices in MOIntsTransformer
    setMORanges(); 

  }; // partition MOSpace
  

  void PostHartreeFockBase::turnOnStateAverage(const std::vector<double> & weight) {
    
    size_t NS = this->NStates;
    
    if( weight.size() != NS) 
      CErr("PostHartreeFockBase needs " + std::to_string(NS) + " weights for state average" );   

    this->StateAverage = true;
    this->SAWeight = std::vector<double>(NS);
    std::copy_n(weight.begin(), NS, this->SAWeight.begin());
  };
  
  /*
   *  Indicies Mapping
   *
   *  'N' -> Negative Energy MO, only for 4C 
   *  'C' -> Frozen Core orbitals 
   *  'I' -> Inactive Core 
   *  'A' -> Active Space for Correlation Calculations 
   *  'S' -> Secondary Virtual 
   *  'V' -> Frozen Virtual Orbitals 
   */
  void PostHartreeFockBase::setCorrelatedSpaceAndReOrder() {

    std::vector<size_t> nOrbs;
    std::vector<char> orbIdentifiers;
    
    size_t totalElecMO = this->corrSpace.nElecMO;
    
    // set up numbers
    if (wfnRef_->nC == 4) {
      nOrbs.push_back(corrSpace.nNegMO);
      orbIdentifiers.push_back('N');
    }

    nOrbs.push_back(corrSpace.nFCore);
    orbIdentifiers.push_back('C');

    nOrbs.push_back(corrSpace.nInact);
    orbIdentifiers.push_back('I');

    nOrbs.push_back(corrSpace.nCorrO);
    orbIdentifiers.push_back('A');

    nOrbs.push_back(corrSpace.nSVirt);
    orbIdentifiers.push_back('S');
    
    nOrbs.push_back(corrSpace.nFVirt);
    orbIdentifiers.push_back('V');
    
    auto & orbIndices = corrSpace.orbIndices;
    size_t scan_start = 0, scan_end = 0;

    if (orbIndices.size() == 0) {
      // initialize the MO indices as default
      for(auto i = 0; i < nOrbs.size(); i++) { 
        scan_end = scan_start + nOrbs[i];
        char & scan_id = orbIdentifiers[i];
        for(size_t j = scan_start; j < scan_end; j++) {
          orbIndices.push_back(scan_id);  
        }
        scan_start = scan_end;
      }
    } else if (orbIndices.size() == totalElecMO) {
      
      // generating swapping pairs if necessary 
      std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;
      moPairs.resize(2, {});
      
      for(auto i = 0; i < nOrbs.size(); i++) {
        scan_end = scan_start + nOrbs[i];
        char & scan_id = orbIdentifiers[i];
        
        for (size_t j = scan_start; j < scan_end; j++) {
          if (orbIndices[j] != scan_id) { 
            bool found = false;
            for (size_t k = scan_end; k < totalElecMO; k++) {
              if (orbIndices[k] == scan_id) {
                found = true;
                moPairs[0].push_back({j+1, k+1});
                orbIndices[k] = orbIndices[j];
                orbIndices[j] = scan_id;
                break;
              }
            }
            if (not found) 
              CErr("Inconsistency between correlated orbital partition and orbial indices");
          }
        }
        scan_start = scan_end;
      }

      this->swapMOs(moPairs,isAlpha);

    } else {
      CErr("Wrong size for corrSpace.orbIndices"); 
    }

  }; // PostHartreeFockBase::setActiveSpaceAndReOrder

}; // namespace ChronusQ


