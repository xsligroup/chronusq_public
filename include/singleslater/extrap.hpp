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
#include <singleslater.hpp>
#include <util/matout.hpp>
#include <cqlinalg/blas1.hpp>

namespace ChronusQ {

 /**
   *  \brief Control routine for DIIS and damping for SCF
   *  
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::modifyFock() {

    ROOT_ONLY(comm); 

    // Static Damping
    scfControls.doDamp &= (std::dynamic_pointer_cast<ROFock<MatsT,IntsT>>(fockBuilder) == nullptr);
    if (scfControls.doDamp) fockDamping();

    // DIIS extrapolation
    if (scfControls.diisAlg == NONE) return;

    if (scfControls.diisAlg == CDIIS) {
      size_t nExtrap = std::min(scfConv.nSCFIter+1,scfControls.nKeep);
      scfDIIS(nExtrap);
    } else CErr("Only CDIIS is implemented so far",std::cout);

  }; // SingleSlater<T>::modifyFock



 /**
   *  \brief Static Fock Damping routine
   *
   *  F^k = (1-dp)*F^k + dp*F^{k-1}
   *  
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::fockDamping() {

    // Don't damp the first iteration.  We don't want 
    // to use the guess Fock and it's not saved anyway.
    if(scfConv.nSCFIter == 0) return;

    size_t NB = this->basisSet().nBasis;
    if( this->nC == 4 ) NB = 2 * NB;
    double dp = scfControls.dampParam;
   
    // Damp the current orthonormal Fock matrix 
    if( not savFile.exists() )
      *fockMatrixOrtho = (1-dp) * *fockMatrixOrtho + dp * *prevFock;
    else {

      cqmatrix::PauliSpinorMatrices<MatsT> FSCR(this->memManager, NB,
          fockMatrixOrtho->hasXY(), fockMatrixOrtho->hasZ());

      if (this->particle.charge < 0.)
        savFile.readData("/SCF/FOCK_ORTHO", FSCR);
      else
        savFile.readData("/PROT_SCF/FOCK_ORTHO", FSCR);

      *fockMatrixOrtho = (1-dp) * *fockMatrixOrtho + dp * FSCR;

    }

  }; // SingleSlater<T>::fockDamping


 /**
   *  \brief Commutator DIIS
   *
   *  Saves the AO fock and density matrices, evaluates the [F,D]
   *  commutator and uses this to extrapolate the Fock and density. 
   *
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::scfDIIS(size_t nExtrap) {

    // Save the current AO Fock and density matrices
    size_t NB    = this->basisSet().nBasis;
    if( this->nC == 4 ) NB = 2 * NB;
    size_t iDIIS = scfConv.nSCFIter % scfControls.nKeep;

    diisFock[iDIIS] = *this->fockMatrix;
    diisOnePDM[iDIIS] = *this->onePDM;

    // Evaluate orthonormal [F,D] and store in diisError
    FDCommutator(diisError[iDIIS]);

    scfConv.nrmFDC = 0.;
    for(auto E : diisError[iDIIS].SZYXPointers())
      scfConv.nrmFDC = std::max(scfConv.nrmFDC,blas::nrm2(NB*NB,E,1));

    // Just save the Fock, density, and commutator for the first iteration
    if (scfConv.nSCFIter == 0) return;
      
    // Build the B matrix and return the coefficients for the extrapolation
    DIIS<MatsT> extrap(nExtrap,diisError);


    if(extrap.extrapolate()) { 
      // Extrapolate Fock and density matrices using DIIS coefficients
      fockMatrix->clear();
      this->onePDM->clear();
      for(auto j = 0; j < nExtrap; j++) {
        *fockMatrix += extrap.coeffs[j] * diisFock[j];
        *this->onePDM += extrap.coeffs[j] * diisOnePDM[j];
      }
    } else {
      if( printLevel > 0 )
        std::cout << "\n    *** WARNING: DIIS Inversion Failed -- "
                  << " Defaulting to Fixed-Point step ***\n" << std::endl;
    }

    // Transform AO fock into the orthonormal basis
    ao2orthoFock();

  }; // SingleSlater<T>::scfDIIS



 /**
   *  \brief Allocates storage for different extrapolation approaches to SCF 
   *  
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::allocExtrapStorage() {

    ROOT_ONLY(comm);

    diisFock.clear();
    diisFock.reserve(scfControls.nKeep);
    diisOnePDM.clear();
    diisOnePDM.reserve(scfControls.nKeep);
    diisError.clear();
    diisError.reserve(scfControls.nKeep);

    // Allocate memory to store previous orthonormal Focks and densities for DIIS
    if (scfControls.diisAlg != NONE) {
      for(auto i = 0; i < scfControls.nKeep; i++) {
        diisFock.emplace_back(memManager, fockMatrix->dimension(),
                                  fockMatrix->hasXY(), fockMatrix->hasZ());
        diisFock.back().clear();
        diisOnePDM.emplace_back(memManager, fockMatrix->dimension(),
                                    fockMatrix->hasXY(), fockMatrix->hasZ());
        diisOnePDM.back().clear();
        diisError.emplace_back(memManager, fockMatrix->dimension(),
                                   fockMatrix->hasXY(), fockMatrix->hasZ());
        diisError.back().clear();
      }
    }

    // Allocate memory to store previous orthonormal Fock for damping 
    if( scfControls.doDamp and not savFile.exists() ) {
      SPIN_OPERATOR_ALLOC(this->basisSet().nBasis,prevFock);
    }

  }; // SingleSlater<T>::allocExtrapStorage



 /**
   *  \brief Deallocates storage for different extrapolation approaches to SCF 
   *  
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::deallocExtrapStorage() {

    ROOT_ONLY(comm);

    // Deallocate memory to store previous orthonormal Focks and densities for DIIS
    if (scfControls.diisAlg != NONE) {
      diisFock.clear();
      diisOnePDM.clear();
      diisError.clear();
    }

    // Deallocate memory to store previous orthonormal Fock for damping 
    prevFock = nullptr;

  }; // SingleSlater<T>::deallocExtrapStorage



 /**
   *  \brief Form the orthonormal [F,D] commutator and store in diisError 
   *  
   *  The commutator is formed using FD and its adjoint
   *              [F,D] = FD - Adj(FD)
   *
   *  TODO: This routine only works for R/UHF right now
   *
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::FDCommutator(cqmatrix::PauliSpinorMatrices<MatsT> &FDC) {

    size_t NB    = this->basisSet().nBasis;
    bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT,IntsT>>(fockBuilder) != nullptr);
    if( this->nC == 4 ) NB = 2 * NB;

    if(this->nC == 1) {
      cqmatrix::Matrix<MatsT> SCR(memManager, NB);

      // FD(S) = F(S)D(S)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), fockMatrixOrtho->S().pointer(), NB,
        onePDMOrtho->S().pointer(), NB, MatsT(0.), FDC.S().pointer(), NB);

      // FD(S) += F(z)D(z)
      if(nC == 2 or !iCS and !iRO) {
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), fockMatrixOrtho->Z().pointer(), NB,
          onePDMOrtho->Z().pointer(), NB, MatsT(0.), SCR.pointer(), NB);
        FDC.S() += SCR;
      }

      // Form {FD - DF}(S)
      SCR = FDC.S();
      MatAdd('N','C', NB, NB, MatsT(1.), FDC.S().pointer(), NB,
             MatsT(-1.), SCR.pointer(), NB, FDC.S().pointer(), NB);


      if(nC == 2 or !iCS and !iRO) {
        // FD(z) = F(S)D(z) + F(z)D(S)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), fockMatrixOrtho->S().pointer(), NB,
          onePDMOrtho->Z().pointer(), NB, MatsT(0.), FDC.Z().pointer(), NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), fockMatrixOrtho->Z().pointer(), NB,
          onePDMOrtho->S().pointer(), NB, MatsT(0.), SCR.pointer(), NB);
        FDC.Z() += SCR;

        // Form {FD - DF}(z)
        SCR = FDC.Z();
        MatAdd('N','C', NB, NB, MatsT(1.), FDC.Z().pointer(), NB,
               MatsT(-1.), SCR.pointer(), NB, FDC.Z().pointer(), NB);
      }
    } else {

      // Gather the orthonormal Fock and densities
      cqmatrix::Matrix<MatsT> FO(fockMatrixOrtho->template spinGather<MatsT>());
      cqmatrix::Matrix<MatsT> DO(onePDMOrtho->template spinGather<MatsT>());
      cqmatrix::Matrix<MatsT> SCR(memManager, 2*NB);

      // Compute FD product
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NB,2*NB,2*NB,MatsT(1.),FO.pointer(),2*NB,
           DO.pointer(),2*NB,MatsT(0.),SCR.pointer(),2*NB);
      
      // Compute FD - DF (Store in FO scratch)
      MatAdd('N','C',2*NB,2*NB,MatsT(1.),SCR.pointer(),2*NB,
             MatsT(-1.),SCR.pointer(),2*NB,FO.pointer(),2*NB);

      // Scatter Product into FDC
      FDC = FO.template spinScatter<MatsT>();

    }

  }; // SingleSlater<T>::FDCommutator





}; // namespace ChronusQ

