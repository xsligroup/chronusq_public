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

#include <configinteraction.hpp>
#include <newcibuilder/impl.hpp>

namespace ChronusQ {
   
template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::computeTDM(
  size_t s1, size_t s2, std::shared_ptr<cqmatrix::Matrix<MatsT>> tdm) {
  ciBuilder->buildTDM(*CIVectors, *CIVectors, {s1, s2}, tdm);
} // ConfigInteraction::computeTDM 


template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::compute2TDM(
  size_t s1, size_t s2, std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM) {
  ciBuilder->buildTDM(*CIVectors, *CIVectors, {s1, s2}, nullptr, false, 1.0, twoTDM);
} // ConfigInteraction::compute2TDM


template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::computeRDMsForOrbitalRotations() {
  
  const auto& weights = this->SAWeight;
  oneRDMSOI->clear();
  twoRDMSOI->clear();
  
  size_t si = (this->StateAverage) ? 0ul: this->NStates - 1;
  
  for (; si < this->NStates; ++si) {
    const auto w = (this->StateAverage) ? weights[si] : 1.0;
#ifdef CQ_ENABLE_MPI 
    bool reduceRDM = (si == this->NStates - 1) ? true : false;    
#endif    
    ciBuilder->buildTDM(*CIVectors, *CIVectors, {si, si}, 
        oneRDMSOI, true, w, twoRDMSOI, true, w
#ifdef CQ_ENABLE_MPI 
        , reduceRDM, reduceRDM
#endif    
        );
  }

  size_t nCorrO = this->corrSpace.nCorrO;
  auto & RDM2 = *twoRDMSOI;
  auto & RDM1 = *oneRDMSOI;
/*
#pragma omp parallel for schedule(static) default(shared)       
  for (auto w = 0ul; w < nCorrO; w++)
  for (auto u = 0ul; u < nCorrO; u++)
  for (auto t = 0ul; t < nCorrO; t++) {
    RDM2(t, u, u, w) -= RDM1(t, w);
  }
*/  
} // ConfigInteraction::computeRDMsOfInterests
  
} // namespace ChronusQ
