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
  void PostHartreeFock<MatsT,IntsT>::rdm2pdm(cqmatrix::Matrix<MatsT> & rdm, double scale) {

    // onePDM(AO)_{uv} = sum_{pq} C_{up} oneRDM(MO)_{pq} C^*_{qv}
    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;
    size_t fourCompOffset = (ref_->nC == 4) ? ref_->nAlphaOrbital() * 2: 0;
    size_t nCoreO = corrSpace.nFCore + corrSpace.nInact;
    size_t nCorrO = corrSpace.nCorrO;

    double fc1C = (ref_->nC == 1) ? 2.0 : 1.0;

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
    else {
      *ref_->onePDM = tmpPDM.template spinScatter<MatsT>();
    }

    ref_->ao2orthoDen();

  } //PostHartreeFock::rdm2pdm

}; // namespace ChronusQ


