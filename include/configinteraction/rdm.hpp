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
  
  if(!ciSettings.SparseDavidson) {
    std::shared_ptr<DistributedVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(CIVectors);
    ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {s1, s2}, tdm);
  }
  else {
#ifdef CQ_ENABLE_SPARSE
    std::shared_ptr<DistributedSparseVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedSparseVectors<MatsT>>(CIVectors);
    ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {s1, s2}, tdm);
#endif
  }
} // ConfigInteraction::computeTDM 


template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::compute2TDM(
  size_t s1, size_t s2, std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM) {
  
  if(!ciSettings.SparseDavidson) {
    std::shared_ptr<DistributedVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(CIVectors);
    ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {s1, s2}, nullptr, false, 1.0, twoTDM);
  }
  else {
#ifdef CQ_ENABLE_SPARSE
    std::shared_ptr<DistributedSparseVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedSparseVectors<MatsT>>(CIVectors);
    ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {s1, s2}, nullptr, false, 1.0, twoTDM);
#endif
  }
 
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
    if(!ciSettings.SparseDavidson) {
      std::shared_ptr<DistributedVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(CIVectors);
  
      ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {si, si}, 
        oneRDMSOI, true, w, twoRDMSOI, true, w
#ifdef CQ_ENABLE_MPI 
        , reduceRDM, reduceRDM
#endif    
        );
    }
    else {

#ifdef CQ_ENABLE_SPARSE
      std::shared_ptr<DistributedSparseVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedSparseVectors<MatsT>>(CIVectors);

      ciBuilder->buildTDM(*CIVectorsCast, *CIVectorsCast, {si, si},
        oneRDMSOI, true, w, twoRDMSOI, true, w
        , reduceRDM, reduceRDM
        );
#endif
    }
  }

  size_t nCorrO = this->corrSpace.nCorrO;
  auto & RDM2 = *twoRDMSOI;
  auto & RDM1 = *oneRDMSOI;

  //subtract one-electron terms from 2TDM.
  //These only occur in case 3b (two-electron excitation is all within the same DAS), so only subtract terms for which p,q==r,s are in the same DAS
  #pragma omp parallel for schedule(static) default(shared)
  for (auto q = 0ul; q < nCorrO; q++) {

    //find which DAS orbital q belongs to
    size_t orbitalRangeStart = 0;
    size_t orbitalRangeEnd = 0;
    for(const auto& activeSpace : ciSettings.activeSpaces) {
      orbitalRangeStart = orbitalRangeEnd;
      orbitalRangeEnd += activeSpace.nOrbitals;

      if(q >= orbitalRangeStart and q < orbitalRangeEnd) {
	break;
      }
    }

    for (auto p = 0ul; p < nCorrO; p++) {
      if(p >= orbitalRangeStart and p < orbitalRangeEnd) {
        for (auto s = 0ul; s < nCorrO; s++) {
	  if(s >= orbitalRangeStart and s < orbitalRangeEnd) {
            RDM2(p, q, q, s) -= RDM1(p, s);
	  }
        }
      }
    }  
  }
  
} // ConfigInteraction::computeRDMsOfInterests
  
} // namespace ChronusQ
