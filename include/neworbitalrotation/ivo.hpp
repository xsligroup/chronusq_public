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
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/ortho.hpp>

// #define _DEBUG_OR_IVOS

namespace ChronusQ {
  
/* 
 * Generate improved virtual orbials 
 */
template <typename MatsT, typename IntsT>
void NewOrbitalRotation<MatsT, IntsT>::generateIVOs(EMPerturbation & pert,
  cqmatrix::Matrix<MatsT> & oneRDM) {
  
  const auto& mopart = postHF_.corrSpace;
  auto& mem    = postHF_.memManager;
  auto& ss     = *postHF_.reference();
  auto& mo     = ss.mo[0];
  
  size_t nC = ss.nC;
  
  if (nC == 1) CErr("1C IVO NYI");
  
  const size_t nAO     = mo.dimension();
  const size_t nTOrb   = mopart.nMO;
  const size_t nCorrO  = mopart.nCorrO;
  const size_t nCorrE  = mopart.nCorrE;
  const size_t nInact  = mopart.nInact;
  const size_t nSVirt  = mopart.nSVirt;
  const size_t coreOrbOff = nC == 4 ? nAO / 2 : 0;
  const size_t corrOrbOff = coreOrbOff + mopart.nFCore + nInact;
  const size_t virtOrbOff = corrOrbOff + nCorrO;
  
  /**************************************/
  /* Step 1: build scaled AO density    */
  /**************************************/ 
  cqmatrix::Matrix<MatsT> SCR(mem, nAO); 
  SCR.clear();
  
  // Cas contribution
  double scale = double(nCorrE - 1) / double(nCorrE);
  postHF_.rdm2pdm(oneRDM, scale);

  /****************************************************/
  /* Step 2: form and diagonalize virtual fock matrix */
  /****************************************************/
  
  // form fock through SingleSlater 
  ss.formFock(pert, false, 1.); 
  SCR = ss.fockMatrix->template spinGather<MatsT>(); 
  
  // transform to MO basis only with virtual block
  cqmatrix::Matrix<MatsT> virtualFock(mem, nSVirt);
  auto off_sizes = postHF_.mointsTF->parseMOType("ab");
  SCR.subsetTransform('N', mo.pointer(), nAO, off_sizes, virtualFock.pointer());  
  
  // diagonalize virtual fock Matrix
  dcomplex * EIVOs = mem.template malloc<dcomplex>(nSVirt);
  MatsT * U = mem.template malloc<MatsT>(nSVirt * nSVirt);
  MatsT * dummy = nullptr;
  GeneralEigen('N','V', nSVirt, virtualFock.pointer(), nSVirt, EIVOs, dummy, 1, U, nSVirt);

  /*****************************************************************/
  /* Step 3: transform virtal MO: C'(mu,b) = \sum_a C(mu,a) U(a,b) */
  /*****************************************************************/
  blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
    nAO, nSVirt, nSVirt, MatsT(1.), mo.pointer() + virtOrbOff * nAO, nAO,
    U, nSVirt, MatsT(0.), SCR.pointer(), nAO);
  SetMat('N', nAO, nSVirt, MatsT(1.), SCR.pointer(), nAO,
    mo.pointer() + virtOrbOff * nAO, nAO);
  
  // print out IVOs eignevalues
  std::fill_n(ss.eps1 + coreOrbOff, nTOrb - nSVirt, 0.);
  for(auto i = 0ul; i < nSVirt; i++) ss.eps1[i + virtOrbOff] = std::real(EIVOs[i]); 

  std::cout << std::endl << std::endl; 
  std::cout << " *-------------------------------------------------------*" << std::endl;      
  std::cout << " * Improved Virtual Orbitals (IVOs) have been generated. *" << std::endl;      
  std::cout << " *-------------------------------------------------------*" << std::endl;      
  std::cout << std::endl;

  ss.printEPS(std::cout);
  
  std::cout << std::endl << std::endl; 

  mem.free(EIVOs, U);

} // NewOrbitalRotation<MatsT>::generateIVOs

} // namespace ChronusQ
