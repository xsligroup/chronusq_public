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

#include <mcwavefunction.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

  /*
   * \brief Compute 1rdm
   *         
   */ 
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::computeOneRDM(size_t i) {

    ciBuilder->computeOneRDM(*this, CIVecs[i], oneRDM[i]);

  } // MCWaveFunction::computeOneRDM

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::computeOneRDM() {

    for (auto i = 0ul; i < NStates; i++)
      MCWaveFunction<MatsT, IntsT>::computeOneRDM(i);

  } // MCWaveFunction::computeOneRDM

  /*
   * \brief Compute TDM
   *
   */
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::computeTDMs() {

    // allocate memory if not
    if (!TDMs.empty()) {
      TDMs.reserve(NStates);
      for (auto i = 0ul; i < NStates; i++) {
        TDMs.emplace_back(std::vector<cqmatrix::Matrix<MatsT>>(NStates,
                  cqmatrix::Matrix<MatsT>(MOPartition.nCorrO)));
        for (auto j = 0ul; j < NStates; j++) 
          ciBuilder->computeTDM(*this, CIVecs[i], CIVecs[j], TDMs[i][j]);
      }
    }

  } // MCWaveFunction::computeTDMs

  /*
   * \brief transform oneRDM(MO) to onePDM(AO)
   *
   */
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::rdm2pdm(cqmatrix::Matrix<MatsT> & rdm, double scale) {

    // onePDM(AO)_{uv} = sum_{pq} C_{up} oneRDM(MO)_{pq}^* C^*_{qv}
    size_t nAO = reference().nAlphaOrbital() * reference().nC;
    size_t fourCompOffset = (reference().nC == 4) ? reference().nAlphaOrbital() * 2: 0;
    size_t nI = MOPartition.nFCore + MOPartition.nInact;
    size_t nCorrO = MOPartition.nCorrO;

    double fc1C = (reference().nC == 1) ? 2.0 : 1.0;

    cqmatrix::Matrix<MatsT> tmpPDM(nAO);
    tmpPDM.clear();

    // Core
    for(auto i = fourCompOffset; i < nI+fourCompOffset; i++) tmpPDM(i,i) = fc1C;

    // Active
    SetMat('R', nCorrO, nCorrO, scale, rdm.pointer(), nCorrO,
            tmpPDM.pointer() + (fourCompOffset+nI)*(nAO+1), nAO);

    tmpPDM = tmpPDM.transform('C', reference().mo[0].pointer(), nAO, nAO);

    if (reference().nC == 1) reference().onePDM->S() = tmpPDM;
    // for 1C, only scalar part is rewritten
    else {
      *reference().onePDM = tmpPDM.template spinScatter<MatsT>();
    }

    reference().ao2orthoDen();

  } //MCWaveFunction::rdm2pdm



}; // namespace ChronusQ


