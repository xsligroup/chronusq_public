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

#include <singleslater.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/distributedritpi.hpp>

namespace ChronusQ {

  /**
   *  \brief Contract the exact-exchange matrix directly from MO coefficients.
   *
   *  Used by both the regular Fock builder (FockBuilder::formRawGDInBatches)
   *  and the range-separated hybrid erfc-K path (KohnSham::formRangeSeparatedHybridExchange)

   */
  template <typename MatsT, typename IntsT>
  void contractExchangeKCoef(
      MPI_Comm comm, bool iCS, bool useCholeskyMOs,
      SingleSlater<MatsT,IntsT>& ss,
      const std::shared_ptr<TPIContractions<MatsT,IntsT>>& contraction,
      cqmatrix::PauliSpinorMatrices<MatsT>& Kout) {

    auto ritpi_incore = std::dynamic_pointer_cast<InCoreRITPIContraction<MatsT,IntsT>>(contraction);
    auto ritpi_dist   = std::dynamic_pointer_cast<DistributedRITPIContraction<MatsT,IntsT>>(contraction);
    const bool isRoot = MPIRank(comm) == 0;

    decltype(ss.getCholeskyMOs()) choleskyMOs;
    if (useCholeskyMOs and isRoot) choleskyMOs = ss.getCholeskyMOs();
    const size_t NB = ss.basisSet().nBasis;

    // Contract one spin block: handles MO selection, broadcast, dispatch.
    auto contractSpin = [&](int spin) -> cqmatrix::Matrix<MatsT> {
      size_t nO = 0;
      MatsT* mo = nullptr;
      if (isRoot) {
        if (useCholeskyMOs) { nO = choleskyMOs[spin]->nColumns(); mo = choleskyMOs[spin]->pointer(); }
        else { nO = (spin == 0) ? ss.nOA : ss.nOB; mo = ss.mo[spin].pointer(); }
      }

      cqmatrix::Matrix<MatsT> Kblock(NB);
      if (ritpi_incore) {
        if (isRoot and nO > 0) ritpi_incore->KCoefContract(comm, nO, mo, Kblock.pointer());
        else Kblock.clear();
      } else {
#ifdef CQ_ENABLE_MPI
        MPIBCast(&nO, 1, 0, comm);
#endif
        if (nO > 0) {
          auto* mo_buf = CQMemManager::get().malloc<MatsT>(NB * nO);
          if (isRoot) std::copy_n(mo, NB * nO, mo_buf);
#ifdef CQ_ENABLE_MPI
          MPIBCast(mo_buf, NB * nO, 0, comm);
#endif
          ritpi_dist->KCoefContract(comm, nO, mo_buf, Kblock.pointer());
          CQMemManager::get().free(mo_buf);
        } else {
          Kblock.clear();
        }
      }
      return Kblock;
    };

    auto AAblock = contractSpin(0);
    auto BBblock = iCS ? cqmatrix::Matrix<MatsT>(NB) : contractSpin(1);

    // Only root assembles spin-block exchange matrices.
    if (isRoot)
      Kout = iCS
        ? cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(AAblock)
        : cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(AAblock, BBblock);
  }

} // namespace ChronusQ
