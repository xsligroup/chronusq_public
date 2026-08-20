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

#include <posthartreefock.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

  /*
   * \brief Compute 1rdm
   */ 
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::computeOneRDM() {
    for (auto i = 0ul; i < NStates; i++) computeOneRDM(i);
  } // PostHartreeFock::computeOneRDM

  /*
   * \brief transform oneRDM(MO) to onePDM(AO) and put it back to ref_
   */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::rdm2pdm(cqmatrix::Matrix<MatsT> & rdm, double scale, bool isTDM) {

    // onePDM(AO)_{uv} = sum_{pq} C_{up} oneRDM(MO)_{pq} C^*_{qv}
    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;
    size_t fourCompOffset = (ref_->nC == 4) ? ref_->nAlphaOrbital() * 2: 0;
    size_t nCoreO = corrSpace.nFCore + corrSpace.nInact;
    size_t nCorrO = corrSpace.nCorrO;
    double fc1C = 0.0;
    if (isTDM == false) {
      fc1C = (ref_->nC == 1) ? 2.0 : 1.0;} 

    cqmatrix::Matrix<MatsT> tmpPDM(nAO);
    tmpPDM.clear();

    // Core
    for(auto i = fourCompOffset; i < nCoreO + fourCompOffset; i++) tmpPDM(i,i) = fc1C;

    // Active
    SetMat('R', nCorrO, nCorrO, scale, rdm.pointer(), nCorrO,
            tmpPDM.pointer() + (fourCompOffset+nCoreO)*(nAO+1), nAO);

    tmpPDM = tmpPDM.transform('C', ref_->mo[0].pointer(), nAO, nAO);
    
    // for 1C, only scalar part is rewritten
    if (ref_->nC == 1) ref_->onePDM->S() = tmpPDM;
    else *ref_->onePDM = tmpPDM.template spinScatter<MatsT>();
    
    ref_->ao2orthoDen();

  } //PostHartreeFock::rdm2pdm


  /* brief: Transforms RDM from the AO basis to the MO basis.
   * Formula: 1rdmMO = C^{\dagger} S 1rdmAO S C
   * Arguments: AO - RDM
   * Return: MO - RDM
   */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::pdm2rdm(cqmatrix::Matrix<MatsT> &rdmAO) {
  
    // oneRDM(MO) = C^{\dagger} S oneRDM(AO) S C
    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;
    size_t nMO = corrSpace.nMO;
    size_t nCorrO = corrSpace.nInact + corrSpace.nFCore;  

    //Obtain overlap: S
    cqmatrix::Matrix<MatsT> S(nAO);
    if(ref_->nC == 1){
      S = ref_->aoints_->overlap->matrix();
    } else if(ref_->nC == 2){
      std::fill_n(S.pointer(),nAO*nAO,MatsT(0.0));
      SetMat('N',nAO/2,nAO/2,MatsT(1.),ref_->aoints_->overlap->matrix().pointer(), nAO/2, S.pointer(),nAO);
      SetMat('N',nAO/2,nAO/2,MatsT(1.),ref_->aoints_->overlap->matrix().pointer(), nAO/2, S.pointer()+nAO*nAO/2+nAO/2,nAO);
    } else if(ref_->nC == 4){
      std::fill_n(S.pointer(),nAO*nAO,MatsT(0.));
      SetMat('N',nAO/4,nAO/4,MatsT(1.),ref_->aoints_->overlap->matrix().pointer(), nAO/4, S.pointer(),nAO);
      SetMat('N',nAO/4,nAO/4,MatsT(1./(2*SpeedOfLight()*SpeedOfLight())),ref_->aoints_->kinetic->matrix().pointer(), nAO/4, S.pointer()+nAO*nAO/4+nAO/4,nAO);
      SetMat('N',nAO/4,nAO/4,MatsT(1.),ref_->aoints_->overlap->matrix().pointer(), nAO/4, S.pointer()+nAO*nAO/2+nAO/2,nAO);
      SetMat('N',nAO/4,nAO/4,MatsT(1./(2*SpeedOfLight()*SpeedOfLight())),ref_->aoints_->kinetic->matrix().pointer(), nAO/4, S.pointer()+nAO*nAO*3/4+nAO*3/4,nAO);
    }

    //Create MO density matrix:
    cqmatrix::Matrix<MatsT> tmpdm1(nAO);
    cqmatrix::Matrix<MatsT> tmpdm2(nAO);
    
    // tmpdm1 = oneRDM(AO) S
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        nAO, nAO, nAO, MatsT(1.), rdmAO.pointer(), nAO, S.pointer(), nAO,
        MatsT(0.),tmpdm1.pointer(), nAO);
    // tmpdm2 = S tmpdm1
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        nAO, nAO, nAO, MatsT(1.), S.pointer(), nAO, tmpdm1.pointer(), nAO,
        MatsT(0.),tmpdm2.pointer(), nAO);
    // tmpdm1 = C^T tmpdm2 = C^T S oneRDM(AO) S
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
        nAO, nAO, nAO, MatsT(1.), ref_->mo[0].pointer(), nAO, tmpdm2.pointer(), nAO,
        MatsT(0.),tmpdm1.pointer(), nAO);
    // tmpdm2 = tmpdm1 C = C^T S oneRDM(SO) S C
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        nAO, nAO, nAO, MatsT(1.), tmpdm1.pointer(), nAO, ref_->mo[0].pointer(), nAO,
        MatsT(0.),tmpdm2.pointer(), nAO);

    if (ref_->nC == 1)  {
      ref_->onePDM->S() = tmpdm2;
    } else {
      *ref_->onePDM = tmpdm2.template spinScatter<MatsT>();
    }

  } //PostHartreeFock:: pdm2rdm

  /*
   * Compute the 2-RDM(MO basis) in the full space:
   * Contains all blocks
   * We will reduce this later for efficiency
   *
   */
  template <typename MatsT, typename IntsT>
  std::shared_ptr<InCore4indexTPI<MatsT>> PostHartreeFock<MatsT,IntsT>::computeFull2RDM(size_t i) {
    
    if (ref_->nC == 1) { CErr("NYI!!"); }

    size_t nMO    = corrSpace.nMO;
    size_t nCorrO = corrSpace.nCorrO;
    size_t nNegMO = corrSpace.nNegMO;
    size_t nFCore = corrSpace.nFCore;
    size_t nInact = corrSpace.nInact;

    // Occupied core = frozen-core + inactive (fully occupied, uncorrelated).
    // To exclude frozen-core electrons, move nFCore into coreStart and out of nCore.
    size_t coreStart = nNegMO;
    size_t nCore     = nFCore + nInact;
    size_t actOff    = nNegMO + nFCore + nInact;
    size_t occEnd    = actOff + nCorrO;
    
    // Compute the One and Two- RDMs for nCorrO
    computeOneRDM(i);
    auto twoRDM = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);
    twoRDM->clear();
    compute2RDM(i, i, twoRDM);

    // OneRDM Cache
    auto& RDM1 = *oneRDM[i];
    auto full2RDM = std::make_shared<InCore4indexTPI<MatsT>>(nMO);
    full2RDM->clear();
    
    // Full Space 1-RDM Look up:
    auto gamma = [&](size_t p, size_t q) -> MatsT {
      bool pC = (p >= coreStart and p < coreStart + nCore);
      bool qC = (q >= coreStart and q < coreStart + nCore);
      if (pC and qC) return (p == q) ? MatsT(1.0) : MatsT(0.0);
      bool pA = (p >= actOff and p < occEnd);
      bool qA = (q >= actOff and q < occEnd);
      if (pA and qA) return RDM1(p - actOff, q - actOff);
      return MatsT(0.0);
    };

    // Do the mixed blocks for 2RDM (core-active)
#pragma omp parallel for schedule(static) default(shared) collapse(2)
    for (size_t s = coreStart; s < occEnd; ++s)
    for (size_t r = coreStart; r < occEnd; ++r)
    for (size_t q = coreStart; q < occEnd; ++q)
    for (size_t p = coreStart; p < occEnd; ++p) {
      if (p >= actOff and q >= actOff and r >= actOff and s >= actOff) continue;
      (*full2RDM)(p,q,r,s) = gamma(p,q) * gamma(r,s) - gamma(p,s) * gamma(r,q);
    }
    
    // Do the correlated block of 2-RDM
#pragma omp parallel for schedule(static) default(shared) collapse(2)
    for (size_t s = 0; s < nCorrO; ++s)
    for (size_t r = 0; r < nCorrO; ++r)
      SetMat('N', nCorrO, nCorrO, MatsT(1.0),
          twoRDM->pointer() + (r + s * nCorrO) * (nCorrO * nCorrO), nCorrO,
          full2RDM->pointer() + actOff + actOff * nMO
            + (r + actOff) * (nMO * nMO) + (s + actOff) *(nMO * nMO * nMO), nMO);

    return full2RDM;

  } // PostHartreeFock::computeFull2RDM

}; // namespace ChronusQ


