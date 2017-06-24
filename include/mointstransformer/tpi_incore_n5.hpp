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

#include <mointstransformer.hpp>
#include <util/timer.hpp>
#include <cqlinalg.hpp>
#include <matrix.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  std::shared_ptr<InCore4indexTPI<MatsT>> 
  MOIntsTransformer<MatsT,IntsT>::formAOTPIInCore(bool cacheAOTPI) {
      
      // cache AOTPI
      auto AOTPI = ints_cache_.getIntegral<InCore4indexTPI,MatsT>("AOTPI");

      if (not AOTPI) { 
        if(ss_.nC == 1) {
          AOTPI = std::make_shared<InCore4indexTPI<MatsT>>(
                    *std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(ss_.aoints_->TPI)); 
        } else if (ss_.nC == 2) {
          std::cout << "  * Using bare Coulomb Operator for 2e Integrals" << std::endl;
          AOTPI = std::make_shared<InCore4indexTPI<MatsT>>(
            std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(ss_.aoints_->TPI)
              ->template spatialToSpinBlock<MatsT>()); 
        } else if (ss_.nC == 4) {
          AOTPI = std::dynamic_pointer_cast<InCore4indexTPI<MatsT>>(
            std::make_shared<InCore4indexRelERI<MatsT>>(
              std::dynamic_pointer_cast<InCore4indexRelERI<IntsT>>(ss_.aoints_->TPI)
                ->template spatialToSpinBlock<MatsT>()));
        }

        if (cacheAOTPI) ints_cache_.addIntegral("AOTPI", AOTPI);
      } 

     return AOTPI; 
  }; // MOIntsTransformer::cacheAOTPIInCore
  
  /**
   *  \brief subset tranform TPI using incore  
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::subsetTransformTPIInCoreN5(
    const std::vector<std::pair<size_t,size_t>> & off_sizes, MatsT * MOTPI, 
    bool cacheIntermediates, TPI_TRANS_DELTA_TYPE delta) {
      
      auto AOTPI = formAOTPIInCore(cacheIntermediates); 
       
      size_t nAO = ss_.mo[0].dimension();
      auto MO = ss_.mo[0].pointer();
      
      MatsT * SCR = nullptr;

      if (delta != NO_KRONECKER_DELTA) {
        size_t np = off_sizes[0].second;
        size_t nq = off_sizes[1].second;
        size_t nr = off_sizes[2].second;
        size_t ns = off_sizes[3].second;
        SCR = CQMemManager::get().malloc<MatsT>(np*nq*nr*ns);
      } else { 
        SCR = MOTPI; 
      }
      
      if (ss_.nC != 4) {
        AOTPI->subsetTransform('N', MO, nAO, off_sizes, SCR);
      } else {
        auto AOTPI4C = std::dynamic_pointer_cast<InCore4indexRelERI<MatsT>>(AOTPI);
        AOTPI4C->subsetTransform('N', MO, nAO, off_sizes, SCR);
      }
      
      if (delta != NO_KRONECKER_DELTA) {
        size_t np = off_sizes[0].second;
        size_t nq = off_sizes[1].second;
        size_t nr = off_sizes[2].second;
        size_t ns = off_sizes[3].second;
        size_t npq = np*nq;
        size_t npqr = npq*nr;
        
        if (delta == KRONECKER_DELTA_PQ) {
          
          // MOTPI(p,r,s) = SCR(p,p,r,s)
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto r = 0ul; r < nr; r++) 
          for (auto s = 0ul; s < ns; s++) {
            size_t rsnr = r + s*nr;
            for (auto p = 0ul; p < np; p++) 
              MOTPI[p + np*rsnr] = SCR[p*(np+1) + npq*rsnr];  
          }
        
        } else if (delta == KRONECKER_DELTA_RS) {
          
          // MOTPI(p,q,r) = SCR(p,q,r,r)
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto p = 0ul; p < np; p++) 
          for (auto q = 0ul; q < nq; q++) {
            size_t pqnp = p + q*np;
            for (auto r = 0ul; r < nr; r++) 
              MOTPI[pqnp + npq*r] = SCR[npq*r*(nr+1) + pqnp];  
          }
        
        } else if (delta == KRONECKER_DELTA_PS) {
          
          // MOTPI(p,q,r) = SCR(p,q,r,p)
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto r = 0ul; r < nr; r++) 
          for (auto q = 0ul; q < nq; q++) {
            size_t qnprnpq = q*np + r*npq;
            for (auto p = 0ul; p < np; p++) 
              MOTPI[p + qnprnpq] = SCR[p*(npqr+1) + qnprnpq];  
          }
        
        } else if (delta == KRONECKER_DELTA_RQ) {

          // MOTPI(p,r,s) = SCR(p,r,r,s)
          size_t npr = np*nr;
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto p = 0ul; p < np; p++) 
          for (auto s = 0ul; s < ns; s++) {
            size_t psnpqr = p + s*npqr;
            size_t psnpr  = p + s*npr;
            for (auto r = 0ul; r < nr; r++) 
              MOTPI[r*np + psnpr] = SCR[r*(np+npq) + psnpqr];  
          }
        
        } else if (delta == KRONECKER_DELTA_PQ_RS) {
          
          // MOTPI(p,r) = SCR(p,p,r,r)
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto p = 0ul; p < np; p++) 
          for (auto r = 0ul; r < nr; r++) 
            MOTPI[p + r*np] = SCR[p*(1+np) + r*(1+nr)*npq];
          
        } else if (delta == KRONECKER_DELTA_PS_RQ) {
        
          // MOTPI(p,r) = SCR(p,r,r,p)
          #pragma omp parallel for schedule(static) collapse(2) default(shared)       
          for (auto p = 0ul; p < np; p++) 
          for (auto r = 0ul; r < nr; r++) 
            MOTPI[p + r*np] = SCR[p*(1+npqr) + r*(1+nr)*np];
        }
        
        CQMemManager::get().free(SCR);
      }

      return; 
  }; // MOIntsTransformer<MatsT,IntsT>::subsetTransformTPIInCoreN5
  
}; // namespace ChronusQ
