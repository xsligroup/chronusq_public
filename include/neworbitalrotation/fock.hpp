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

#include <neworbitalrotation.hpp>
#include <mointstransformer/impl.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/ortho.hpp>

namespace ChronusQ {
  
/* 
 * form generalized Fock matrix 1
 *  
 * F1_pq = factor * (hCore_pq + \sum_tu 1PDM_{tu} [(pq|tu) - xFC * (pu|tq)])
 *
 *  for 1C,        factor = 2.0, xFC = 0.5
 *  for 2C and 4C, factor = 1.0, xFC = 1.0
 *
 */
template <typename MatsT, typename IntsT>
void NewOrbitalRotation<MatsT, IntsT>::formGeneralizedFock1(EMPerturbation & pert,
  cqmatrix::Matrix<MatsT> & oneRDM, MatsT * F1, const std::string & moType,
  bool deltaPQ) {
  
  // unpack pq_size
  const auto& corrS = postHF_.corrSpace;
  const auto pq_sizes = postHF_.mointsTF->parseMOType(moType);
  const size_t np   = pq_sizes[0].second;
  const size_t nq   = pq_sizes[1].second;
  const size_t npq  = np * nq;
  const size_t npqDim = deltaPQ ? np: npq;
  const size_t nCorrO  = corrS.nCorrO;
  const size_t nTOrb   = corrS.nMO;
  const size_t nInact  = corrS.nInact;
  const size_t nSVirt  = corrS.nSVirt;
  const size_t nCorrO2 = nCorrO * nCorrO;
  const size_t SCRDim  = std::max(nCorrO, std::max(nInact, nSVirt));
  if (not orbitalGradient_)
    orbitalGradient_ = std::make_shared<cqmatrix::Matrix<MatsT>>(corrS.nMO);

  
  // populate F1 as hCore_pq
  postHF_.mointsTF->transformHCore(pert, F1, moType, deltaPQ, 'i');
  
  const auto tsize  = postHF_.mointsTF->parseMOType("t");
  const size_t toff = tsize[0].first;
  const size_t nt   = tsize[0].second;
  const size_t nAO  = postHF_.reference()->nAlphaOrbital() * postHF_.reference()->nC;
  
  cqmatrix::Matrix<MatsT> Den(nAO);
  MatsT * SCR = CQMemManager::get().malloc<MatsT>(std::max(nAO*nt, npqDim)); 
  MatsT * tMO = postHF_.reference()->mo[0].pointer() + toff*nAO;

  // form densty D(nu,mu) = C(nu, u) 1PDM_{tu} C(mu, t)^*
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, 
    nAO, nt, nt, MatsT(1.), tMO, nAO, oneRDM.pointer(), nt, MatsT(0.), SCR, nAO);

  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
    nAO, nAO, nt, MatsT(1.), SCR, nAO, tMO, nAO, MatsT(0.), Den.pointer(), nAO);
  
  // get G(Dtu)_pq
  postHF_.mointsTF->transformGD(pert, Den, true, SCR, moType, deltaPQ, true, "tu");
  
  // add to F1
  if (postHF_.reference()->nC == 1) blas::scal(npqDim, MatsT(2.), F1, 1);
  blas::axpy(npqDim, MatsT(1.), SCR, 1, F1, 1);
  
  CQMemManager::get().free(SCR);
} // NewOrbitalRotation<MatsT>::formGeneralizedFock1

/* 
 * form generalized Fock matrix 2
 *  
 * F2_tq = \sum_u 1PDM_tu hCore_qu +  \sum_uvw 2PDM_tuvw (qu|vw)
 *
 */
template <typename MatsT, typename IntsT>
void NewOrbitalRotation<MatsT, IntsT>::formGeneralizedFock2(EMPerturbation & pert,
    cqmatrix::Matrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM, MatsT * F2, 
    const std::string & moType) {
  
  // unpack q_size
  const auto q_size = postHF_.mointsTF->parseMOType(std::string(moType));
  const size_t nq   = q_size[0].second;
  const size_t nCorrO = postHF_.corrSpace.nCorrO;
  const size_t nCorrO2 = nCorrO * nCorrO;
  const size_t nCorrO3 = nCorrO2 * nCorrO;
  
  MatsT * SCR = CQMemManager::get().malloc<MatsT>(nq * nCorrO3); 
  
  // populate SCR as hCore_qu 
  postHF_.mointsTF->transformHCore(pert, SCR, moType + "u", false, 'i');
  
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, 
    nCorrO, nq, nCorrO, MatsT(1.), oneRDM.pointer(), nCorrO,
    SCR, nq, MatsT(0.), F2, nCorrO);

  // populate ERI as (qu|vw)
  postHF_.mointsTF->transformTPI(pert, SCR, moType + "uvw", true, false);

  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, 
    nCorrO, nq, nCorrO3, MatsT(1.), twoRDM.pointer(), nCorrO,
    SCR, nq, MatsT(1.), F2, nCorrO);
  
  CQMemManager::get().free(SCR);
  
} // NewOrbitalRotation<MatsT>::formGeneralizedFock2

} // namespace ChronusQ
