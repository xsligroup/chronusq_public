/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
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

#include <mcwavefunction/base.hpp>
#include <detstringmanager.hpp>

namespace ChronusQ {

  void MCWaveFunctionBase::partitionMOSpace(std::vector<size_t> nActO, size_t nCorrE) {
  
    auto & mopart = this->MOPartition;
    auto & wfn = referenceWaveFunction();
    
    if (wfn.nC == 4 and not this->FourCompNoPair)
      CErr("4C without no-pair approximation is not implemented.");

    if (wfn.nC == 1) {
      mopart.nElecMO = wfn.nAlphaOrbital();
    } else { 
      mopart.nElecMO = wfn.nV + wfn.nO;
    } 
    
    
    if (wfn.nC == 4) {
      mopart.nMO = 2 * mopart.nElecMO;
      mopart.nNegMO = mopart.nElecMO;
    } else {
      mopart.nMO = mopart.nElecMO; 
    }

    size_t nCorrO = std::accumulate(nActO.begin(), nActO.end(), 0);
    
    if( nCorrO > mopart.nElecMO or nCorrO <=0 or (nCorrE > wfn.nO) or nCorrE <= 0)  
      CErr("This is not a valid active space, please modify the input!");
    
    if (wfn.nC == 1) {
      if (wfn.iCS) {
        mopart.nCorrEA = nCorrE / 2;
        mopart.nCorrEB = mopart.nCorrEA;  
      } else {
        int singleE   = wfn.nOA - wfn.nOB;
        mopart.nCorrEA = (nCorrE + singleE) / 2;
        mopart.nCorrEB = nCorrE - mopart.nCorrEA;  
      }

      if (nCorrO < mopart.nCorrEA) 
        CErr("The specified number of correlated orbitals can not hold all the correlated electrons"); 
    
      mopart.nInact = wfn.nOA - mopart.nCorrEA;
    } else { // other than 1C
  
      if (nCorrO < nCorrE)
        CErr("The specified number of correlated orbitals can not hold all the correlated electrons"); 
      
      mopart.nInact = wfn.nO - nCorrE;
    }    
    
    mopart.nFVirt = mopart.nElecMO - nCorrO - mopart.nInact;
    mopart.nCorrO = nCorrO;
    mopart.nCorrE = nCorrE;

    if( mopart.nElecMO < (nCorrO + mopart.nInact) )
      CErr("This is not a valid active space, please modify the input!");
    
  
    #define CONSTRUCT_CASSTRINGMANAGER(M_, N_) \
      std::dynamic_pointer_cast<DetStringManager>( \
        std::make_shared<CASStringManager>(M_, N_))

    // Initialize String Engine
    if (mopart.scheme == CAS) {
      if (wfn.nC == 1) {
        this->NDet = Comb(nCorrO, mopart.nCorrEA) *
                     Comb(nCorrO, mopart.nCorrEB);
        this->detStr = CONSTRUCT_CASSTRINGMANAGER(nCorrO, mopart.nCorrEA);
      
        this->detStrBeta = (wfn.iCS) ? detStr:
          CONSTRUCT_CASSTRINGMANAGER(nCorrO, mopart.nCorrEB);
      } else {
        this->NDet = Comb(nCorrO, nCorrE);
        this->detStr = CONSTRUCT_CASSTRINGMANAGER(nCorrO, nCorrE);
      }
    } else if (mopart.scheme == RAS) {
      size_t mhRas1 = std::min(nActO[0],mopart.mxHole);
      size_t meRas3 = std::min(nActO[2],mopart.mxElec);
      if(wfn.nC ==1) {
        if(mopart.mxHole > nActO[0]*2 or mopart.mxElec > nActO[2]*2)
          CErr("Invalid Input for RAS Specification");
        if((mopart.nCorrEA + mhRas1) > nCorrO)
          CErr("RAS3 is too small to generate full MxHole.");
        if(mopart.nCorrEA < meRas3)
          CErr("RAS1 is too small to generate full MxElec.");
        std::shared_ptr<RASStringManager> rasStr = std::make_shared<RASStringManager>(
            nActO, mopart.nCorrEA, mopart.mxHole, mopart.mxElec);
        std::shared_ptr<RASStringManager> rasStrBeta = (wfn.iCS) ? rasStr:
            std::make_shared<RASStringManager>(nActO, mopart.nCorrEB,
            mopart.mxHole, mopart.mxElec);
        mopart.fCat = genfCat(rasStr->LCategory(), rasStrBeta->LCategory(),
            rasStr->nCategory(), rasStrBeta->nCategory());
        this->NDet = rasStr->nString() * rasStrBeta->nString();
        this->detStr = std::dynamic_pointer_cast<DetStringManager>(rasStr);

        this->detStrBeta = (wfn.iCS) ? detStr:
                std::dynamic_pointer_cast<DetStringManager>(rasStrBeta);
      } else { // nC ==2 or 4
        if(mopart.mxHole > nActO[0] or mopart.mxElec > nActO[2])
          CErr("Invalid Input for RAS Specification");
        if((mopart.nCorrE + mhRas1) > nCorrO)
          CErr("RAS3 is too small to generate full MxHole.");
        if(mopart.nCorrE < meRas3)
          CErr("RAS1 is too small to generate full MxElec.");
        std::shared_ptr<RASStringManager> rasStr = std::make_shared<RASStringManager>(
               nActO, nCorrE, mopart.mxHole, mopart.mxElec);
        mopart.fCat = genfCat(rasStr->LCategory(), rasStr->nCategory());
        this->NDet = rasStr->nString();
        this->detStr = std::dynamic_pointer_cast<DetStringManager>(rasStr);
      }
      mopart.nActOs = nActO;  // Number of orbitals for each RAS space.
    } else
      CErr("DetStringManager other than CAS or RAS is not implemented"); 
    
    // initialize orbital index
    setActiveSpaceAndReOrder();   
    
    // initalize indices in MOIntsTransformer 
    setMORanges(); 

  }; // partiation MOSpace
  

  void MCWaveFunctionBase::turnOnStateAverage(const std::vector<double> & weight) {
    
    size_t NS = this->NStates;
    
    if( weight.size() != NS) 
      CErr("MCSCF needs "+std::to_string(NS)+" weights for state average" );   

    this->StateAverage = true;
    this->SAWeight = std::vector<double>(NS);
    std::copy_n(weight.begin(), NS, this->SAWeight.begin());
  };
  
  // Generate the category offset for nC=1: fCat[iCata+iCatb*nCat]
  std::vector<int> MCWaveFunctionBase::genfCat(std::vector<std::vector<size_t>> LCata,
    std::vector<std::vector<size_t>> LCatb, size_t nCata, size_t nCatb) {

    std::vector<int> fCat(nCata * nCatb, 0);
    size_t nDet = 0;

    size_t iCatb = 0;
    for (auto iEb = 0; iEb <= this->MOPartition.mxElec; iEb++)
    for (auto iHb = 0; iHb <= this->MOPartition.mxHole; iHb++, iCatb++) {
      size_t iCata = 0;
      for (auto iEa = 0; iEa <= this->MOPartition.mxElec; iEa++)
      for (auto iHa = 0; iHa <= this->MOPartition.mxHole; iHa++, iCata++) {
        if (LCata[iCata][3] == 0 or LCatb[iCatb][3] == 0) fCat[iCata + iCatb*nCata] = -2;
        else if ((iEb+iEa) > this->MOPartition.mxElec or (iHb+iHa) > this->MOPartition.mxHole)
          fCat[iCata + iCatb*nCata] = -1;
        else {
          fCat[iCata + iCatb*nCata] = nDet;
          nDet += LCata[iCata][3] * LCatb[iCatb][3];
        }
      }
    }

    return fCat;

  } // MCWaveFunctionBase::genfCat

  // Generate the category offset for nC>1: fCat[iCat]
  std::vector<int> MCWaveFunctionBase::genfCat(std::vector<std::vector<size_t>> LCat,
        size_t nCat) {

    std::vector<int> fCat(nCat, 0);
    size_t nDet = 0;

    for (auto iCat = 0; iCat < nCat; iCat++) {
      if (LCat[iCat][3] == 0) fCat[iCat] = -2;
      else {
        fCat[iCat] = nDet;
        nDet += LCat[iCat][3];
      }
    }

    return fCat;

  }; // MCWaveFunctionBase::genfCat
  

  /*
   *  'N' -> Negative Energy MO, only for 4C 
   *  'I' -> Inactive  / inactive
   *  'A' -> Active   "1, 2, 3" for RAS
   *  'S' -> Secondary / frozen virtual
   */
  void MCWaveFunctionBase::setActiveSpaceAndReOrder() {

    std::vector<size_t> nOrbs;
    std::vector<char> orbIdentifiers;
    
    if (referenceWaveFunction().nC == 4) {
      nOrbs.push_back(MOPartition.nNegMO);
      orbIdentifiers.push_back('N');
    }

    // set up numbers 
    nOrbs.push_back(MOPartition.nInact);
    orbIdentifiers.push_back('I');

    if (MOPartition.scheme != RAS) { 
      nOrbs.push_back(MOPartition.nCorrO);
      orbIdentifiers.push_back('A');
    } else {
      nOrbs.push_back(MOPartition.nActOs[0]);
      orbIdentifiers.push_back('1');
      nOrbs.push_back(MOPartition.nActOs[1]);
      orbIdentifiers.push_back('2');
      nOrbs.push_back(MOPartition.nActOs[2]);
      orbIdentifiers.push_back('3');
    }

    nOrbs.push_back(MOPartition.nFVirt);
    orbIdentifiers.push_back('S');
    
    auto & orbIndices = MOPartition.orbIndices;
    size_t scan_start = 0, scan_end = 0;

    if (orbIndices.size() == 0) {
      // initialize the MO indices
      for(auto i = 0; i < nOrbs.size(); i++) { 
        scan_end = scan_start + nOrbs[i];
        char & scan_id = orbIdentifiers[i];
        for(size_t j = scan_start; j < scan_end; j++) {
          orbIndices.push_back(scan_id);  
        }
        scan_start = scan_end;
      }
    } else if (orbIndices.size() == MOPartition.nMO) {
      
      // generating swapping pairs if necessary 
      std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;
      moPairs.resize(2, {});
      
      for(auto i = 0; i < nOrbs.size(); i++) {
        scan_end = scan_start + nOrbs[i];
        char & scan_id = orbIdentifiers[i];
        
        for (size_t j = scan_start; j < scan_end; j++) {
          if (orbIndices[j] != scan_id) { 
            bool found = false;
            for (size_t k = scan_end; k < MOPartition.nMO; k++) {
              if (orbIndices[k] == scan_id) {
                found = true;
                moPairs[0].push_back({j+1, k+1});
                orbIndices[k] = orbIndices[j];
                orbIndices[j] = scan_id;
                break;
              }
            }
            if (not found) CErr("inconsistency between number of orbitals and orbial indices");
          }
        }
        scan_start = scan_end;
      }

      this->swapMOs(moPairs,isAlpha);
    } else {
      CErr("wrong sizes set in MOPartition.orbIndices"); 
    }

  }; // MCWaveFunctionBase::setActiveSpaceAndReOrder

}; // namespace ChronusQ


