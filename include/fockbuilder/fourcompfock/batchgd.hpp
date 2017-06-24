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

#include <fockbuilder/fourcompfock.hpp>

namespace ChronusQ {

  /**   
   *  \brief Forms the 4C GD in Bathes
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formRawGDInBatches(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & onePDMs, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & coulombMatrices, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & exchangeMatrices,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & twoeHs) {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  };


  /**   
   *  \brief Forms the 4C Fock matrix using AO-direct
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formRawGDInBatchesDirect(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & onePDMs, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & coulombMatrices, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & exchangeMatrices,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & twoeHs) {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }; // FourCompFock<MatsT, IntsT>::formRawGDInBatchesDirect

  /*******************************************************************************/
  /*                                                                             */
  /* Compute memory requirement for build 4C GD in Batches                       */
  /* Returns:                                                                    */
  /*   size_t SCR size needed for one batch                                      */
  /*   IMPORTANT HERE: size are all in MatsT (dcomplex)                          */
  /*******************************************************************************/
  template <typename MatsT, typename IntsT>
  size_t FourCompFock<MatsT,IntsT>::formRawGDSCRSizePerBatch(SingleSlater<MatsT,IntsT> &ss,
    bool computeExchange, bool HerDen) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

    return 0;

  }; // FourCompFock<MatsT, IntsT>::formRawGDSCRSizePerBatch



}; // namespace ChronusQ
