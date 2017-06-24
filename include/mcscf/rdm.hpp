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

#include <mcscf.hpp>

namespace ChronusQ {
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT, IntsT>::computeOneRDM(size_t i) {
    MCWaveFunction<MatsT, IntsT>::computeOneRDM(i);
    if (this->settings.doSCF) *oneRDMSOI = this->oneRDM[i];
  }; // MCSCF::computeOneRDM(i)
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::computeOneRDM() {

    MCWaveFunction<MatsT, IntsT>::computeOneRDM();
    
    if (this->StateAverage) {
      
      size_t nCorrO  = this->MOPartition.nCorrO;
      auto & weights = this->SAWeight;

      oneRDMSOI->clear();
      
      for (auto i = 0ul; i < this->NStates; i++)
        *oneRDMSOI += weights[i] * this->oneRDM[i];
    }    

  }; // MCSCF::computeOneRDM()

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT, IntsT>::computeTwoRDM(size_t i) {
    this->ciBuilder->computeTwoRDM(*this, this->CIVecs[i], *twoRDMSOI);
    
    size_t nCorrO  = this->MOPartition.nCorrO;
    if(this->detStr->scheme() == PRECOMPUTED_CONFIGURATION_DRIVEN_LIST) {
      auto & RDM2 = *twoRDMSOI;
      auto & RDM1 = *oneRDMSOI;
#pragma omp parallel for schedule(static) default(shared)       
      for (auto w = 0ul; w < nCorrO; w++)
      for (auto u = 0ul; u < nCorrO; u++)
      for (auto t = 0ul; t < nCorrO; t++)
        RDM2(t, u, u, w) -= RDM1(t, w);
    }
  }; // MCSCF::computeTwoRDM(i)
  
  // compute state-average two RDM
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::computeTwoRDM() {
    
    size_t nCorrO  = this->MOPartition.nCorrO;
    size_t nCorrO2 = nCorrO * nCorrO; 
    
    if(this->StateAverage) {
      
      auto twoRDMtmp = InCore4indexTPI<MatsT>(nCorrO); 
      auto & weights = this->SAWeight;
      twoRDMSOI->clear();

      for (auto i = 0ul; i < this->NStates; i++) {
        this->ciBuilder->computeTwoRDM(*this, this->CIVecs[i], twoRDMtmp);
        MatAdd('N', 'N', nCorrO2, nCorrO2, 
          MatsT(weights[i]), twoRDMtmp.pointer(), nCorrO2, 
          MatsT(1.), twoRDMSOI->pointer(), nCorrO2, twoRDMSOI->pointer(), nCorrO2);
      }  
      
      if(this->detStr->scheme() == PRECOMPUTED_CONFIGURATION_DRIVEN_LIST) {
        auto & RDM2 = *twoRDMSOI;
        auto & RDM1 = *oneRDMSOI;
#pragma omp parallel for schedule(static) default(shared)       
        for (auto w = 0ul; w < nCorrO; w++)
        for (auto u = 0ul; u < nCorrO; u++)
        for (auto t = 0ul; t < nCorrO; t++)
          RDM2(t, u, u, w) -= RDM1(t, w);
      }

    } else {
      CErr("Doesn't make sense to compute all two RDM for state specific case");
    }
    
  };  // MCSCF::computeTwoRDM     

}; // namespace ChronusQ
