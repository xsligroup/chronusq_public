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

#include <cqlinalg/blas3.hpp>
#include <cqlinalg/factorization.hpp>
#include <singleslater/hartreefock.hpp>

#include <response.hpp>

namespace ChronusQ {

template<typename MatsT, typename IntsT>
void HartreeFock<MatsT, IntsT>::computeFullNRStep(MatsT* orbRot) {

#if 1
  this->ao2orthoFock();
  this->MOFOCK();   // Update MOFock
  auto ptr = getPtr();
  PolarizationPropagator<HartreeFock<MatsT, IntsT>> resp(this->comm, FDR, ptr);
  resp.doNR              = true;
  resp.fdrSettings.bFreq = {0.001};

  std::cout << "  * STARTING HESSIAN INVERSE FOR NEWTON-RAPHSON\n";
  resp.run();
  std::cout << "\n\n";

  // Make sure that the Hessian gets deallocated
  MatsT* FM = resp.fullMatrix();
  CQMemManager::get().free(FM);

  MatsT* C = resp.fdrResults.SOL;

  size_t N = resp.getNSingleDim(false);
  std::copy_n(C, N, orbRot);
#else
  CErr();
#endif
};

template<typename MatsT, typename IntsT>
std::pair<double, MatsT*> HartreeFock<MatsT, IntsT>::getStab() {

  MatsT* CCPY = nullptr;
  double W    = 0;
#if 1
  this->MOFOCK();   // Update MOFock
  auto ptr = getPtr();
  PolarizationPropagator<HartreeFock<MatsT, IntsT>> resp(this->comm, RESIDUE, ptr);
  resp.doStab = true;

  std::cout << "  * STARTING HESSIAN DIAG FOR STABILITY CHECK\n";
  resp.run();
  std::cout << "\n\n";

  // Make sure that the Hessian gets deallocated
  MatsT* FM = resp.fullMatrix();
  CQMemManager::get().free(FM);

  MatsT* C = resp.resResults.VR;

  size_t N = resp.getNSingleDim(false);
  CCPY     = CQMemManager::get().malloc<MatsT>(N);
  std::copy_n(C, N, CCPY);

  W = resp.resResults.W[0];
#else
  CErr();
#endif

  return {W, CCPY};
};

};   // namespace ChronusQ
