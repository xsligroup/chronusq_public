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

#include <extrapolate.hpp>
#include <interpolate.hpp>

namespace ChronusQ {

/**
 *  \brief Control routine for DIIS and damping for SCF
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::modifyFock(EMPerturbation& pert) {

  ROOT_ONLY(this->comm);

  // Static Damping
  if( this->doingDamp ) fockDamping();

  // DIIS extrapolation
  if( this->scfControls.diisAlg == NONE ) return;

  // Prepare for DIIS
  vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
  vecShrdPtrMat<MatsT> den  = this->orbitalModifierDrivers.getOnePDM();
  size_t iDIIS = this->scfConv.nSCFIter % this->scfControls.nKeep;
  for( size_t i = 0; i < fock.size(); i++ ) {
    diisFock[iDIIS][i]   = *fock[i];
    diisOnePDM[iDIIS][i] = *den[i];
    diisError[iDIIS][i]  = orbGrad[i];
  }

  size_t nExtrap = std::min(this->scfConv.nSCFIter + 1, this->scfControls.nKeep);

  // Choose which algorithm
  if( this->scfControls.diisAlg == CDIIS ) {
    scfCDIIS(nExtrap, iDIIS);
  } else if( this->scfControls.diisAlg == EDIIS ) {
    scfEDIIS(nExtrap, iDIIS, pert);
  } else if( this->scfControls.diisAlg == CEDIIS ) {
    scfCEDIIS(nExtrap, iDIIS, pert);
  } else {
    CErr("Only CDIIS, EDIIS and CEDIIS is implemented so far", std::cout);
  }

};   // SingleSlater<T>::modifyFock

/**
 *  \brief Static Fock Damping routine
 *
 *  F^k = (1-dp)*F^k + dp*F^{k-1}
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::fockDamping() {

  // Don't damp the first iteration.  We don't want
  // to use the guess Fock and it's not saved anyway.
  if( this->scfConv.nSCFIter == 0 ) return;

  vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
  vecShrdPtrMat<MatsT> den  = this->orbitalModifierDrivers.getOnePDM();
  double dp                                              = this->scfControls.dampParam;
  for( size_t i = 0; i < fock.size(); i++ ) {

    // Damp the current orthonormal Fock matrix
    *fock[i] = (1 - dp) * *fock[i] + dp * prevFock[i];
    *den[i]  = (1 - dp) * *den[i] + dp * prevOnePDM[i];

    // Save a copy for next time
    prevFock[i]   = *fock[i];
    prevOnePDM[i] = *den[i];
  }
  ao2orthoFock(fock);
  ao2orthoDen(den);

};   // ConventionalSCF<RefType,MatsT,IntsT,IntsT>::fockDamping

/**
 *  \brief Commutator DIIS
 *
 *  Saves the AO fock and density matrices, evaluates the [F,D]
 *  commutator and uses this to extrapolate the Fock and density.
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::scfCDIIS(size_t nExtrap, size_t iDIIS) {

  //  Just save the Fock, density, and commutator for the first iteration
  if( this->scfConv.nSCFIter == 0 ) return;

  // Build the B matrix and return the coefficients for the extrapolation
  DIIS<MatsT> extrap(nExtrap, diisError);
  bool conv = extrap.extrapolate();
  if( conv ) {
    diisCombineMat(extrap.coeffs, nExtrap);
  } else {
    std::cout << "\n    *** WARNING: DIIS Inversion Failed -- "
              << " Defaulting to Fixed-Point step ***\n"
              << std::endl;
  }
};   // ConventionalSCF<RefType,MatsT,IntsT,IntsT>::scfDIIS

/**
 *  \brief Energy DIIS
 *
 *  Evaluates the EDIIS error metric and uses this to interpolate the Fock and density matrices.
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::scfEDIIS(size_t nExtrap, size_t iDIIS, EMPerturbation& pert) {

  // Compute Energy with current density matrix
  // The total energy from the last iteration was
  // computed by contracting the new density with the
  // old twoEH matrix. i.e. D_i * G[D_i-1]
  this->orbitalModifierDrivers.computeProperties(pert);
  diisEnergy[iDIIS] = this->orbitalModifierDrivers.getTotalEnergy();

  //  Just save the Fock, density, and coupling for the first iteration
  if( this->scfConv.nSCFIter == 0 ) return;

  // Evaluate Error metric/Interpolation matrix
  ediisErrorMetric(iDIIS, nExtrap);

  ENERGYDIIS<double> interp(nExtrap, diisBMat->pointer(), this->scfControls.nKeep, diisEnergy, iDIIS, 0.1 * this->scfControls.eneConvTol);
  bool conv = interp.interpolate();

  if( conv ) {
    // Convert double to MatsT and reorder indices with ediisMap
    std::vector<MatsT> c(nExtrap, MatsT(0.));
    for( size_t i = 0; i < nExtrap; i++ )
      c[i] = MatsT(interp.coeffs[i]);

    diisCombineMat(c, nExtrap);
  } else {
    std::cout << "\n    *** WARNING: EDIIS Optimization Failed -- "
              << " Defaulting to Fixed-Point step ***\n"
              << std::endl;
  }
};   // ConventionalSCF<RefType,MatsT,IntsT,IntsT>::scfEDIIS

/*
 *   Brief: Compute the Fock/Density Matrices using the Commutator/Energy DIIS method
 *
 */

template<typename MatsT>
void ConventionalSCF<MatsT>::scfCEDIIS(size_t nExtrap, size_t iDIIS, EMPerturbation& pert) {

  this->orbitalModifierDrivers.computeProperties(pert);
  diisEnergy[iDIIS] = this->orbitalModifierDrivers.getTotalEnergy();

  //  Just save the Fock, density, and commutator for the first iteration
  if( this->scfConv.nSCFIter == 0 ) return;

  // Compute CDIIS coeffs
  DIIS<MatsT> extrap(nExtrap, diisError);
  bool convCDIIS = extrap.extrapolate();

  // Compute EDIIS
  double diisSwitch = this->scfControls.cediisSwitch;
  ediisErrorMetric(iDIIS, nExtrap);
  bool convEDIIS  = true;
  double errorMax = computeFDCConv();
  ENERGYDIIS<double> interp(nExtrap, diisBMat->pointer(), this->scfControls.nKeep, diisEnergy, iDIIS, 0.1 * this->scfControls.eneConvTol);
  if( errorMax > 1E-3 * diisSwitch ) { convEDIIS = interp.interpolate(); }

  // If failed default to fixed step
  if( not convEDIIS ) { std::cout << "\n    *** WARNING: EDIIS Optimization Failed -- " << std::endl; }
  if( not convCDIIS ) {
    std::cout << "\n    *** WARNING: CDIIS Inversion Failed -- " << std::endl;
    for( auto& c : extrap.coeffs )
      c = MatsT(0.);
    extrap.coeffs[iDIIS] = MatsT(1.);
  }

  // Combine the CDIIS and EDIIS coeffs
  double scaleEDIIS = 1.;
  double scaleCDIIS = 0.;
  if( errorMax < diisSwitch and errorMax > 1E-3 * diisSwitch ) {
    // Smoothly transition from CDIIS to EDIIS
    scaleEDIIS = errorMax / diisSwitch;
    scaleCDIIS = 1. - scaleEDIIS;
  } else if( errorMax <= 1E-3 * diisSwitch ) {
    scaleEDIIS = 0.;
    scaleCDIIS = 1.;
  }
  for( size_t i = 0; i < nExtrap; i++ )
    extrap.coeffs[i] *= MatsT(scaleCDIIS);
  if( scaleEDIIS > 0. ) {
    for( size_t i = 0; i < nExtrap; i++ )
      extrap.coeffs[i] += MatsT(scaleEDIIS * interp.coeffs[i]);
  }
  diisCombineMat(extrap.coeffs, nExtrap);
}   // ConventionalSCF<RefType,MatsT,IntsT,IntsT> :: scfCEDIIS

/*
 *   Brief: Compute the error metric for energy DIIS which is
 *          (Fi-Fj)\cdot(Di-Dj)
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::ediisErrorMetric(size_t iDIIS, size_t nExtrap) {
  size_t N    = this->scfControls.nKeep;
  size_t nMat = diisOnePDM[iDIIS].size();

  diisBMat->clear();
  for( size_t a = 0; a < nMat; a++ ) {
    size_t NB = diisOnePDM[iDIIS][a].dimension();

    cqmatrix::Matrix<MatsT> dF(NB);
    cqmatrix::Matrix<MatsT> dD(NB);

    // Compute the coupling Matrix for EDIIS
    // (F_i - F_j)\cdot(D_i - D_j)
    size_t NB2 = NB * NB;
    for( size_t j = 0; j < nExtrap; j++ ) {
      dF.clear();
      dD.clear();

      dD = diisOnePDM[iDIIS][a] - diisOnePDM[j][a];
      dF = diisFock[iDIIS][a] - diisFock[j][a];

      diisBMat->pointer()[j + iDIIS * N] += 0.5 * std::real(blas::dot(NB2, dF.pointer(), 1, dD.pointer(), 1));
      diisBMat->pointer()[iDIIS + j * N] += diisBMat->pointer()[j + iDIIS * N];
    }
  }
}   // ConventionalSCF<RefType,MatsT,IntsT,IntsT> :: ediisErrorMetric

/*
 *   Brief: Take the coefficient vector from CDIIS, or CEDIIS and linearly combine
 *          DIIS fock and Density Matrices.
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::diisCombineMat(std::vector<MatsT> c, size_t nExtrap) {
  // Extrapolate Fock and density matrices using DIIS coefficients
  vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
  vecShrdPtrMat<MatsT> den  = this->orbitalModifierDrivers.getOnePDM();

  for( size_t a = 0; a < fock.size(); a++ ) {
    fock[a]->clear();
    den[a]->clear();
    for( auto j = 0; j < nExtrap; j++ ) {
      *fock[a] += c[j] * diisFock[j][a];
      *den[a] += c[j] * diisOnePDM[j][a];
    }
  }
  ao2orthoFock(fock);
  ao2orthoDen(den);
}   // ConventionalSCF<RefType,MatsT,IntsT,IntsT> :: diisCombineMat

/**
 *  \brief Allocates storage for different extrapolation approaches to SCF
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::allocExtrapStorage() {

  ROOT_ONLY(this->comm);

  diisFock.clear();
  diisFock.reserve(this->scfControls.nKeep);
  diisOnePDM.clear();
  diisOnePDM.reserve(this->scfControls.nKeep);
  diisError.clear();
  diisError.reserve(this->scfControls.nKeep);
  diisEnergy.clear();
  diisEnergy = std::vector<double>(this->scfControls.nKeep, 0.);
  diisBMat   = std::make_shared<cqmatrix::Matrix<double>>(this->scfControls.nKeep);
  diisBMat->clear();

  // Allocate memory to store previous orthonormal Focks and densities for DIIS
  vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
  if( this->scfControls.diisAlg != NONE ) {
    for( auto i = 0; i < this->scfControls.nKeep; i++ ) {
      std::vector<cqmatrix::Matrix<MatsT>> f;
      std::vector<cqmatrix::Matrix<MatsT>> d;
      std::vector<cqmatrix::Matrix<MatsT>> e;
      for( auto a = 0; a < fock.size(); a++ ) {
        f.emplace_back(fock[a]->dimension());
        d.emplace_back(fock[a]->dimension());
        e.emplace_back(fock[a]->dimension());
      }
      diisFock.push_back(f);
      diisOnePDM.push_back(d);
      diisError.push_back(e);
    }
  }

  // Allocate memory to store previous orthonormal Fock for damping
  if( this->scfControls.doDamp ) {
    prevFock.clear();
    prevOnePDM.clear();
    prevFock.reserve(fock.size());
    prevOnePDM.reserve(fock.size());
    for( size_t a = 0; a < fock.size(); a++ ) {
      prevFock.emplace_back(fock[a]->dimension());
      prevOnePDM.emplace_back(fock[a]->dimension());
    }
  }

};   // ConventionalSCF<RefType,MatsT,IntsT,IntsT>::allocExtrapStorage

template<typename MatsT>
void ConventionalSCF<MatsT>::computeOrbGradient(std::vector<cqmatrix::Matrix<MatsT>>& grad) {
  if( this->orbitalModifierDrivers.computeErrorVector ) {
    this->orbitalModifierDrivers.computeErrorVector(grad);
  } else {
    FDCommutator(grad);
  }
}

/**
 *  \brief Form the orthonormal [F,D] commutator and store in diisError
 *
 *  The commutator is formed using FD and its adjoint
 *              [F,D] = FD - Adj(FD)
 *
 *  TODO: This routine only works for R/UHF right now
 *
 */
template<typename MatsT>
void ConventionalSCF<MatsT>::FDCommutator(std::vector<cqmatrix::Matrix<MatsT>>& FDC) {

  ao2orthoDen();
  for( size_t a = 0; a < fockMatrixOrtho.size(); a++ ) {
    size_t NB = fockMatrixOrtho[a].dimension();
    cqmatrix::Matrix<MatsT> SCR(NB);
    FDC[a].clear();

    // Compute F*D
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), fockMatrixOrtho[a].pointer(), NB, onePDMOrtho[a].pointer(), NB, MatsT(0.), SCR.pointer(), NB);
    // Compute FD - DF
    MatAdd('N', 'C', NB, NB, MatsT(1.), SCR.pointer(), NB, MatsT(-1.), SCR.pointer(), NB, FDC[a].pointer(), NB);
  }
};   // ConventionaSCF<MatsT>::FDCommutator

template<typename MatsT>
double ConventionalSCF<MatsT>::computeFDCConv() {
  // Compute the Max element
  double maxGrad = 0.;
  for( size_t a = 0; a < orbGrad.size(); a++ ) {
    size_t NB = orbGrad[a].dimension();
    for( size_t i = 0; i < NB * NB; i++ )
      if( maxGrad < std::abs(orbGrad[a].pointer()[i]) ) maxGrad = std::abs(orbGrad[a].pointer()[i]);
  }
  return maxGrad;
};   // ConventionalSCF<MatsT> :: computeFDCConv

};   // namespace ChronusQ
