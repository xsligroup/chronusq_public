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
#include <cqlinalg/ortho.hpp>

// #define DEBUG_NEWORBITALROTATION_GRADIENT

namespace ChronusQ {
  
/*
 * Build gradient block-wisely:
 *   IN -> inactive
 *   CO -> correlated
 *   FV -> frozen virtual
 *   NA -> negative MO for 4C only
 */
template <typename MatsT, typename IntsT>  
double NewOrbitalRotation<MatsT, IntsT>::computeOrbGradient(EMPerturbation & pert,
  cqmatrix::Matrix<MatsT> & oneRDM, InCore4indexTPI<MatsT> & twoRDM) {
  
  ProgramTimer::tick("Form Gradient");
  std::cout << "Compute Orbital Gradient" << std::endl;

  const auto& corrS = postHF_.corrSpace;
  auto& mem = postHF_.memManager;
  if (not orbitalGradient_) 
    orbitalGradient_ = std::make_shared<cqmatrix::Matrix<MatsT>>(mem, corrS.nMO);

  size_t nTOrb   = corrS.nMO;
  size_t nCorrO  = corrS.nCorrO;
  size_t nInact  = corrS.nInact;
  size_t nSVirt  = corrS.nSVirt;
  size_t nINCO   = nInact + nCorrO;
  size_t SCRDim  = std::max(nCorrO, std::max(nInact, nSVirt));
  
  // Allocate SCR
  MatsT * SCR  = mem.template malloc<MatsT>(SCRDim*SCRDim);
  MatsT * SCR2 = mem.template malloc<MatsT>(nCorrO*nInact);
  
  // zero out G
  std::fill_n(orbitalGradient_->pointer(), nTOrb * nTOrb, MatsT(0.));
  
  // point to start of the positive energy part 
  MatsT * G = orbitalGradient_->pointer() + corrS.nNegMO * (nTOrb + 1);
  
  // CO-CO block
  if (settings.rotate_within_correlated) {
    CErr("rotate MOs within correlated orbitals not implemented");
  } 
  
  // IN-CO block
  if (settings.rotate_inact_correlated) {
    // g_it = -2 * (F1_it - F2_ti)
    // populate SCR as F1_it and SCR2 as F2_ti
    this->formGeneralizedFock1(pert, oneRDM, SCR, "it");
    this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR2, "i");
    
    IMatCopy('T', nCorrO, nInact, MatsT(1.), SCR2, nCorrO, nInact);
    MatAdd('N', 'N', nInact, nCorrO, MatsT(1.), SCR, nInact,
      MatsT(-1.), SCR2, nInact, SCR, nInact);
    //  blas::scal(nInact * nCorrO, MatsT(2.), SCR, 1);
    
    SetMat('N', nInact, nCorrO, -MatsT(1.), SCR, nInact, G + nTOrb * nInact, nTOrb);
    SetMat('C', nInact, nCorrO,  MatsT(1.), SCR, nInact, G + nInact, nTOrb);
  } 
  

  // IN-FV block
  if (settings.rotate_inact_virtual) {
    // g_ia = -2 * F1_ia
    // pupulate SCR as F1_ia
    this->formGeneralizedFock1(pert, oneRDM, SCR, "ia");
    //  blas::scal(nInact * nSVirt, MatsT(2.), SCR, 1);
    SetMat('N', nInact, nSVirt, -MatsT(1.), SCR, nInact, G + nTOrb * nINCO, nTOrb);
    SetMat('C', nInact, nSVirt,  MatsT(1.), SCR, nInact, G + nINCO, nTOrb);
  } 

  // CO-FV block
  if (settings.rotate_correlated_virtual) { 
     // g_ta = -2 * F2_ta^*
     // pupulate SCR as F2_ta
    this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR, "a");
    // blas::scal(nCorrO * nSVirt, MatsT(2.), SCR, 1); 
    SetMat('R', nCorrO, nSVirt, -MatsT(1.), SCR, nCorrO, G + nTOrb * nINCO + nInact, nTOrb);
    SetMat('T', nCorrO, nSVirt,  MatsT(1.), SCR, nCorrO, G + nTOrb * nInact + nINCO, nTOrb);
  } 

  mem.free(SCR, SCR2);
  
  // repointing to the begining
  G = orbitalGradient_->pointer();
  
  // NA-PO block
  // which is the formualed similarly as IN-FV block and CO-FV block
  if (settings.rotate_negative_positive) {
    size_t nNegMO = corrS.nNegMO;
    SCR  = mem.template malloc<MatsT>(nNegMO * std::max(nInact, nCorrO)); 
    
    // IN-NA block
    if (nInact > 0) {
      this->formGeneralizedFock1(pert, oneRDM, SCR, "in");
      SetMat('N', nInact, nNegMO, -MatsT(1.), SCR, nInact, G + nNegMO, nTOrb);
      SetMat('C', nInact, nNegMO,  MatsT(1.), SCR, nInact, G + nNegMO * nTOrb, nTOrb);
    }

    // CO-NA Block
    this->formGeneralizedFock2(pert, oneRDM, twoRDM, SCR, "n");
    SetMat('R', nCorrO, nNegMO, -MatsT(1.), SCR, nCorrO, G + nNegMO + nInact, nTOrb);
    SetMat('T', nCorrO, nNegMO,  MatsT(1.), SCR, nCorrO, G + (nNegMO + nInact) * nTOrb, nTOrb);
    
    mem.free(SCR);
  }

  double orbitalGradientNorm =  lapack::lange(lapack::Norm::Fro, nTOrb, nTOrb, G, nTOrb); 

#ifdef DEBUG_NEWORBITALROTATION_GRADIENT
  prettyPrintSmart(std::cout, " new 1RDM ", oneRDM.pointer(), nCorrO, nCorrO, nCorrO);
  prettyPrintSmart(std::cout, " new 2RDM ", twoRDM.pointer(), nCorrO*nCorrO, nCorrO*nCorrO, nCorrO*nCorrO);
  std::cout << "new 2RDM norm = " << std::setprecision(16) << 
     lapack::lange(lapack::Norm::Fro, nCorrO*nCorrO, nCorrO*nCorrO, twoRDM.pointer(), nCorrO*nCorrO)
     << std::endl;
  prettyPrintSmart(std::cout, " OR Orbital Gradient ", G, nTOrb, nTOrb, nTOrb);
  std::cout << "OR Orbital Gradient Norm = " << orbitalGradientNorm << std::endl;
#endif  
  
  ProgramTimer::tock("Form Gradient");
  
  std::cout << "Compute Orbital Gradient End" << std::endl;
  
  return orbitalGradientNorm;
} // NewOrbitalRotation<MatsT>::computeOrbGradient

} // namespace ChronusQ
