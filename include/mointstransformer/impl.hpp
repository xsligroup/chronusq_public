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
#include <mointstransformer/moranges.hpp>
#include <mointstransformer/hcore.hpp>
#include <mointstransformer/tpi_ssfock.hpp>
#include <mointstransformer/tpi_incore_n5.hpp>
#include <mointstransformer/tpi_full_direct.hpp>
#include <util/timer.hpp>

namespace ChronusQ {

/**
 *  \brief main interface to transform Core Hamiltonian:
 *      
 *     HCore(p,q) = \sum_{mu,nu} C(mu,p)^H * C(nu,q) * [h1e(mu,nu) 
 *               optional{ + gdFC * \sum_{lm,sg,i} C(lm,i)^H * C(sg,i) *
 *                              [(mu nu | lm sg) - xFC * (mu sg | lm nu)]}]  
 *
 *    gdFC ... two body part contribution factor
 *             1.0 for 2C/4C and 2.0 for 1C,
 *     xFC ... exchange contribution factor, 
 *             1.0 for 2C/4C and 0.5 for 1C
 *
 *  \param [in] pert      ... external EM pertubations
 *  \param [in] MOHCore   ... resulted HCore in MO basis
 *  \param [in] moType    ... mo indices for p and q, respectively 
 *  \param [in] deltaPQ   ... if index p and q are identical
 *  \param [in] coreIndex ... include core contributions if specified with core index char
 *                            otherwise MOHCore will just be H1e, default null char
 */
template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::transformHCore(EMPerturbation & pert,
  MatsT* MOHCore, const std::string & moType, bool deltaPQ, const char coreIndex) {
  auto off_sizes = parseMOType(moType);
  subsetTransformHCore(pert, off_sizes, MOHCore, deltaPQ, coreIndex);

}; // MOIntsTransformer::transformHCore

/**
 *  \brief main interface to transform G(D)
 *    
 *     G(D)(p, q) = gdFC * \sum_{mu,nu} C(mu,p)^H * C(nu,q) *  
 *                    \sum_{lm,sg} D(lm,sg) *
 *                       [(mu nu | lm sg) - xFC * (mu sg | lm nu)]}]  
 *
 *    gdFC ... two body part contribution factor
 *             1.0 for 2C/4C and 2.0 for 1C,
 *     xFC ... exchange contribution factor, 
 *             1.0 for 2C/4C and 0.5 for 1C
 *
 *  \param [in] pert      ... external EM pertubations
 *  \param [in] MOGD      ... resulted G(D) in MO basis
 *  \param [in] moType    ... mo indices for p and q, respectively 
 *  \param [in] deltaPQ   ... if index p and q are identical
 *  \param [in] cacheAOGD ... controls whether to cache AOGD, default false
 *  \param [in] cacheId   ... used to identify caches, default none
 */
template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::transformGD(EMPerturbation & pert,
  const cqmatrix::Matrix<MatsT> & Den, bool HerDen,
  MatsT* MOGD, const std::string & moType, bool deltaPQ,
  bool cacheAOGD, const std::string & cacheId) {
  auto off_sizes = parseMOType(moType);
  subsetTransformGD(pert, Den, HerDen, MOGD, off_sizes, deltaPQ, 
    cacheAOGD, cacheId);

} // MOIntsTransformer::transformGD

/**
 *  \brief main interface to transform G(InactD) 
 *    
 *     G(D)(p, q) = gdFC * \sum_{mu,nu} C(mu,p)^H * C(nu,q) *  
 *                    \sum_{lm,sg,i} C(lm,i)^H * C(sg,i) *
 *                       [(mu nu | lm sg) - xFC * (mu sg | lm nu)]}]  
 *
 *    gdFC ... two body part contribution factor
 *             1.0 for 2C/4C and 2.0 for 1C,
 *     xFC ... exchange contribution factor, 
 *             1.0 for 2C/4C and 0.5 for 1C
 *
 *  \param [in] pert      ... external EM pertubations
 *  \param [in] MOGD      ... resulted G(D) in MO basis
 *  \param [in] coreIndex ... core index for computing densities 
 *  \param [in] moType    ... mo indices for p and q, respectively 
 *  \param [in] deltaPQ   ... if index p and q are identical
 *  \param [in] cacheAOGD ... controls whether to cache AOGD, default false
 *  \param [in] cacheId   ... used to identify caches, default none
 */
template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::transformGD(EMPerturbation & pert,
  const char coreIndex, MatsT* MOGD, const std::string & moType,
  bool deltaPQ, bool cacheAOGD, const std::string & cacheId) {

  auto Den  = formInactDen(coreIndex);
  transformGD(pert, *Den, true, MOGD, moType, deltaPQ, cacheAOGD, cacheId);

} // MOIntsTransformer::transformHCore

/**
 *  \brief main interface to transform TPI
 *
 *    (p q | r s) = \sum_{mu,nu,lm,sg} C(mu,p)^H * C(lm,r)^H * 
 *                     [(mu nu | lm sg) - xFC * (mu sg | lm nu) ] * C(nu,q) * C(sg,s)
 *
 *     xFC ... exchange contribution factor, 
 *             1.0 for 2C/4C and 0.5 for 1C
 
 *  \param [in] pert               ... external EM pertubations
 *  \param [in] MOTPI              ... resulted TPI in MO basis
 *  \param [in] moType             ... mo indices for p, q, r and s, respectively 
 *  \param [in] cacheIntermediates ... controls whether to cache AO integrals for 
 *                                     InCore N5 and cache half transformed 
 *                                     integrals for InCore N6 and DIRECT N6
 *  \param [in] withExchange       ... whether compute Exchange part of the integrals
 *  \param [in] delta              ... whether two or more indice are indentical
 *                                     see options at include/mointstransformer.hpp 
 */
template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::transformTPI(EMPerturbation & pert, 
  MatsT* MOTPI, const std::string & moType, bool cacheIntermediates, 
  bool withExchange, TPI_TRANS_DELTA_TYPE delta) {
  
  // check delta
  if (delta != NO_KRONECKER_DELTA) {    
    if ((delta == KRONECKER_DELTA_PQ and getUniqueSymbol(moType[0]) != getUniqueSymbol(moType[1])) or
        (delta == KRONECKER_DELTA_RS and getUniqueSymbol(moType[2]) != getUniqueSymbol(moType[3])) or
        (delta == KRONECKER_DELTA_PS and getUniqueSymbol(moType[0]) != getUniqueSymbol(moType[3])) or
        (delta == KRONECKER_DELTA_RQ and getUniqueSymbol(moType[1]) != getUniqueSymbol(moType[2])) or
        (delta == KRONECKER_DELTA_PQ_RS and (getUniqueSymbol(moType[0]) != getUniqueSymbol(moType[1]) 
          or getUniqueSymbol(moType[2]) != getUniqueSymbol(moType[3]))) or
        (delta == KRONECKER_DELTA_PS_RQ and (getUniqueSymbol(moType[0]) != getUniqueSymbol(moType[3]) 
          or getUniqueSymbol(moType[1]) != getUniqueSymbol(moType[2]))))
      CErr("Specifiying KRONECKER_DELTA on different dimensions");
  }
  
  std::string C_moType = moType;
  std::string X_moType = ""; 
  TPI_TRANS_DELTA_TYPE C_delta = delta;
  TPI_TRANS_DELTA_TYPE X_delta = NO_KRONECKER_DELTA;
  bool swapCX = false;

  if (withExchange) {
    X_moType += moType[0];
    X_moType += moType[3];
    X_moType += moType[2];
    X_moType += moType[1];
    
    if      (delta == KRONECKER_DELTA_PQ)    X_delta = KRONECKER_DELTA_PS;
    else if (delta == KRONECKER_DELTA_RS)    X_delta = KRONECKER_DELTA_RQ;
    else if (delta == KRONECKER_DELTA_PS)    X_delta = KRONECKER_DELTA_PQ;
    else if (delta == KRONECKER_DELTA_RQ)    X_delta = KRONECKER_DELTA_RS;
    else if (delta == KRONECKER_DELTA_PQ_RS) X_delta = KRONECKER_DELTA_PS_RQ;
    else if (delta == KRONECKER_DELTA_PS_RQ) X_delta = KRONECKER_DELTA_PQ_RS;
    
    // KRONECKER_DELTA_PS    -> KRONECKER_DELTA_PQ
    // KRONECKER_DELTA_RQ    -> KRONECKER_DELTA_RS
    // KRONECKER_DELTA_PS_RQ -> KRONECKER_DELTA_PQ_RS
    if (delta == KRONECKER_DELTA_PS or 
        delta == KRONECKER_DELTA_RQ or
        delta == KRONECKER_DELTA_PS_RQ) { 
      std::swap(C_moType, X_moType);
      std::swap(C_delta,  X_delta);
      swapCX = true;
    }
    
    //std::cout << "C_moType = " << C_moType << ", C_delta =" << C_delta << std::endl;
    //std::cout << "X_moType = " << X_moType << ", X_delta =" << X_delta << std::endl;
 } 

  auto off_sizes = parseMOType(C_moType);

  // get the Coulomb part
  if (TPITransAlg_ == DIRECT_N6 or TPITransAlg_ == INCORE_N6) {
    subsetTransformTPISSFockN6(pert, off_sizes, MOTPI, C_moType, cacheIntermediates, C_delta);
  } else if (TPITransAlg_ == INCORE_N5) {
    subsetTransformTPIInCoreN5(off_sizes, MOTPI, cacheIntermediates, C_delta);
  } else {
    CErr("DIRECT_N5 NYI");
  }

  // get exchange part if needed
  // for 2c/4c: it's antisymmetrized integrals (pq|rs) - (ps|rq)
  // for 1c: it's scaled for spins, so exchange part will be scaled with 0.5. 
  //         So (pq|rs) - 0.5 * (ps|rq) will be computed
  if (withExchange) {
    
    size_t np = off_sizes[0].second;
    size_t nq = off_sizes[1].second;
    size_t nr = off_sizes[2].second;
    size_t ns = off_sizes[3].second;
    size_t npq = np * nq;
    size_t nps = np * ns;
    size_t npr = np * nr;
    size_t npqr = npq * nr;
    size_t nprs = nps * nr;
    
    MatsT fc = ss_.nC == 1 ? 0.5: 1.0;
    
    size_t qoff = off_sizes[1].first;
    size_t soff = off_sizes[3].first;
    bool qsSymm = nq == ns and qoff == soff; 
    
    if (C_delta == NO_KRONECKER_DELTA and qsSymm) {
      #pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto r = 0ul; r < nr; r++) 
      for (auto p = 0ul; p < np; p++) {
        MatsT tmp1, tmp2;
        size_t pqrs, psrq;
        size_t prnpq = p + r*npq;
        size_t prnps = p + r*nps;
        for (auto s = 0ul; s < ns; s++)
        for (auto q = 0ul; q <= s; q++) {
          pqrs = prnpq + q*np + s*npqr;
          psrq = prnps + s*np + q*nprs;
          tmp1 = MOTPI[pqrs];
          tmp2 = MOTPI[psrq];
          MOTPI[pqrs] = tmp1 - fc * tmp2;
          MOTPI[psrq] = tmp2 - fc * tmp1;
        }
      }
    } else {
      
      size_t SCRSize;
      if      (C_delta == NO_KRONECKER_DELTA) SCRSize = npqr*ns;
      else if (C_delta == KRONECKER_DELTA_RS) SCRSize = npqr;
      else if (C_delta == KRONECKER_DELTA_PQ) SCRSize = nprs;
      else if (C_delta == KRONECKER_DELTA_PQ_RS) SCRSize = npr;
      MatsT * SCR = CQMemManager::get().malloc<MatsT>(SCRSize);

      transformTPI(pert, SCR, X_moType, true, false, X_delta); 

      if (C_delta == NO_KRONECKER_DELTA) { 
        // (p,q,r,s) -= (p,s,q,r)
        #pragma omp parallel for schedule(static) collapse(2) default(shared)       
        for (auto s = 0ul; s < ns; s++)
        for (auto r = 0ul; r < nr; r++) 
        for (auto q = 0ul; q < nq; q++) 
        for (auto p = 0ul; p < np; p++) { 
          MOTPI[p + q*np + r*npq + s*npqr] -= fc * SCR[p + s*np + r*nps + q*nprs]; 
        }
      } else if (C_delta == KRONECKER_DELTA_PQ or C_delta == KRONECKER_DELTA_RS) {
        
        // for KRONECKER_DELTA_PQ: (p,r,s) -= (p,s,r)
        // for KRONECKER_DELTA_RS: (p,q,r) -= (p,r,q)
        size_t nt = KRONECKER_DELTA_PQ ? ns: nq;
        size_t npt = np*nt;
        MatsT fc_inv = swapCX ? 1./ fc: fc;
        
        #pragma omp parallel for schedule(static) collapse(2) default(shared)       
        for (auto t = 0ul; t < nt; t++) 
        for (auto r = 0ul; r < nr; r++) 
        for (auto p = 0ul; p < np; p++) { 
          MOTPI[p + r*np + t*npr] -= fc_inv * SCR[p + t*np + r*npt];
        }
        
        // (p,r,t) --> - (p,t,r) if swapCX 
        if (swapCX) {
          #pragma omp parallel for 
          for (auto p = 0ul; p < np; p++) {
            MatsT tmp;  
            size_t prt, ptr; 
            for (auto t = 0ul; t < nt; t++) 
            for (auto r = 0ul; r <= t; r++) {
              prt = p + r*np + t*npr;
              ptr = p + t*np + r*npt;
              tmp = MOTPI[prt];
              MOTPI[prt] = - MOTPI[ptr] * fc;
              MOTPI[ptr] = - tmp * fc;
            }
          } 
        }
      } else if (C_delta == KRONECKER_DELTA_PQ_RS) {
        // (p,p,r,r) -= (p,r,r,p)
        MatsT * C_term = MOTPI, *X_term = SCR; 
        if (swapCX) std::swap(X_term, C_term);

        MatAdd('N', 'N', np, nr, MatsT(1.), C_term, np, MatsT(-fc), X_term, np, MOTPI, np);
      
      }// delta

      CQMemManager::get().free(SCR);
    } 
  } // with Exchange  

} // MOIntsTransformer::transformTPI 

/**
 *  \brief main interface to direct transform TPI
 *
 *    (p q | r s) = \sum_{mu,nu,lm,sg} C(mu,p)^H * C(lm,r)^H * 
 *                     (mu nu | lm sg) * C(nu,q) * C(sg,s)
 *
 *     xFC ... exchange contribution factor, 
 *             1.0 for 2C/4C and 0.5 for 1C
 
 *  \param [in] pert               ... external EM pertubations
 *  \param [in] MOTPI              ... resulted TPI in MO basis
 *  \param [in] moType             ... mo indices for p, q, r and s, respectively 
 *  \param [in] delta              ... whether two or more indice are indentical
 *                                     see options at include/mointstransformer.hpp 
 */
template <typename MatsT, typename IntsT>
void MOIntsTransformer<MatsT,IntsT>::directTransformTPI(EMPerturbation & pert,
    MatsT* MOTPI, const std::string& moType) {
  
  auto off_sizes = parseMOType(moType);
  directTransformTPI(pert, MOTPI, off_sizes);

}


} // namespace ChronusQ
