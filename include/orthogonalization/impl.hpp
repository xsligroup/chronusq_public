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

#include <orthogonalization.hpp>

namespace ChronusQ {

//#define _DEBUGORTHO

/**
 *  \brief Allocate, compute and store the orthonormalization matricies
 *  over the CGTO basis.
 *
 *  Computes either the Lowdin or Cholesky transformation matricies based
 *  on orthoType
 */
template<typename MatsT>
void Orthogonalization<MatsT>::computeOrtho() {
  if( not overlap ) CErr("Overlap has not been initialized");

  size_t NB  = overlap->dimension();
  size_t nSQ = NB * NB;

  forwardTrans->clear();
  backwardTrans->clear();

  // Allocate scratch
  MatsT* SCR1 = CQMemManager::get().template malloc<MatsT>(nSQ);
  std::fill_n(SCR1, nSQ, 0.);

  if( orthoType == LOWDIN ) {

    // Allocate more scratch
    MatsT* sE   = CQMemManager::get().template malloc<MatsT>(NB);
    MatsT* SCR2 = CQMemManager::get().template malloc<MatsT>(nSQ);

    // Diagonalize the overlap in scratch S = V * s * V**T
    std::copy_n(overlap->pointer(),nSQ,SCR1);
    HermetianEigen('V', 'U', NB, SCR1, NB, sE);

    if( std::abs(sE[0]) < 1e-10 ) CErr("Contracted Basis Set is Linearly Dependent!");

    // Compute X = V * s^{-1/2}
    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ )
        SCR2[i + j * NB] = SCR1[i + j * NB] / std::sqrt(sE[j]);

    // Compute O1 = X * V**T
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, MatsT(1.), SCR2, NB, SCR1, NB, MatsT(0.), forwardTrans->pointer(), NB);

    // Compute X = V * s^{1/2} in place (by multiplying by s)
    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ )
        SCR2[i + j * NB] = SCR2[i + j * NB] * sE[j];

    // Compute O2 = X * V**T
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, MatsT(1.), SCR2, NB, SCR1, NB, MatsT(0.), backwardTrans->pointer(), NB);

#ifdef _DEBUGORTHO
    // Debug code to validate the Lowdin orthogonalization
    prettyPrintSmart(std::cout,"Overlap", overlap->pointer(),NB,NB,NB);
    prettyPrintSmart(std::cout,"Operator Nonortho->Ortho Transformation", forwardTrans->pointer(),NB,NB,NB);
    prettyPrintSmart(std::cout,"Operator Ortho->Nonortho Transformation", backwardTrans->pointer(),NB,NB,NB);

    std::cerr << "Debugging Lowdin Orthogonalization" << std::endl;
    double maxDiff(-10000000);

    // Check that ortho1 and ortho2 are inverses of eachother
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), forwardTrans->pointer(), NB, backwardTrans->pointer(), NB, MatsT(0.), SCR1, NB);

    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ ) {

        if( i == j )
          maxDiff = std::max(maxDiff, std::abs(1. - SCR1[i + j * NB]));
        else
          maxDiff = std::max(maxDiff, std::abs(SCR1[i + j * NB]));
      }

    std::cerr << "  Ortho1 * Ortho2 = I: " << maxDiff << std::endl;

    // Check that ortho2 * ortho2 is the overlap
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), backwardTrans->pointer(), NB, backwardTrans->pointer(), NB, MatsT(0.), SCR1, NB);

    maxDiff = -100000;

    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ ) {

        maxDiff = std::max(maxDiff, std::abs(SCR1[i + j * NB] - overlap->pointer()[i + j * NB]));
      }

    std::cerr << "  Ortho2 * Ortho2 = S: " << maxDiff << std::endl;

    // Check that ortho1 * ortho1 is the inverse of the overlap
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), forwardTrans->pointer(), NB, forwardTrans->pointer(), NB, MatsT(0.), SCR1, NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), SCR1, NB, overlap->pointer(), NB, MatsT(0.), SCR2, NB);

    maxDiff = -10000;
    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ ) {

        if( i == j )
          maxDiff = std::max(maxDiff, std::abs(1. - SCR2[i + j * NB]));
        else
          maxDiff = std::max(maxDiff, std::abs(SCR2[i + j * NB]));
      }

    std::cerr << "  Ortho1 * Ortho1 * S = I: " << maxDiff << std::endl;

#endif

    // Free Scratch Space
    CQMemManager::get().free(sE, SCR2);

  } else if( orthoType == CHOLESKY ) {

    std::cout << "*** WARNING: Cholesky orthogonalization has not yet been confirmed ***" << std::endl;
    CErr("Cholesky orthogonalization should be checked thoroughly before using it");

    // Compute the Cholesky factorization of the overlap S = L * L**T
    lapack::potrf(lapack::Uplo::Lower,NB,SCR1,NB);

    // Copy the lower triangle to ortho2 (O2 = L)
    for( auto j = 0; j < NB; j++ )
      for( auto i = j; i < NB; i++ )
        backwardTrans->pointer()[i + j * NB] = SCR1[i + j * NB];

    // Compute the inverse of the overlap using the Cholesky factors
    lapack::potri(lapack::Uplo::Lower,NB,SCR1,NB);

    // O1 = O2**T * S^{-1}
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), backwardTrans->pointer(), NB, SCR1, NB, MatsT(0.), forwardTrans->pointer(), NB);

    // Remove upper triangle junk from O1
    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < j; i++ )
        forwardTrans->pointer()[i + j * NB] = 0.;

#ifdef _DEBUGORTHO
    // Debug code to validate the Lowdin orthogonalization

    std::cerr << "Debugging Cholesky Orthogonalization" << std::endl;

    // Debug code to validate the Cholesky orthogonalization
    MatsT* SCR2 = CQMemManager::get().template malloc<MatsT>(nSQ);

    double maxDiff = -1000;
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), forwardTrans->pointer(), NB, overlap->pointer(), NB, MatsT(0.), SCR1, NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), SCR1, NB, forwardTrans->pointer(), NB, MatsT(0.), SCR2, NB);

    for( auto j = 0; j < NB; j++ )
      for( auto i = 0; i < NB; i++ ) {

        if( i == j )
          maxDiff = std::max(maxDiff, std::abs(1. - SCR2[i + j * NB]));
        else
          maxDiff = std::max(maxDiff, std::abs(SCR2[i + j * NB]));
      }

    std::cerr << "Ortho1**T * S ** Ortho1 = I: " << maxDiff << std::endl;

    CQMemManager::get().free(SCR2);   // Free SCR2
#endif
  }

  CQMemManager::get().free(SCR1);
};

  //=======================
  // cqmatrix::Matrix Operators
  //=======================
  template<typename MatsT>
  cqmatrix::Matrix<MatsT> Orthogonalization<MatsT>::nonortho2ortho(cqmatrix::Matrix<MatsT> & mat) const {
    if( not overlap ) CErr("Overlap has not been initialized and computed");
    if( forwardTrans->dimension() != mat.dimension() ) CErr("Matrices are not the same dimension in nonortho2ortho");

    size_t NB = forwardTrans->dimension();
    return mat.transform('N', forwardTrans->pointer(), NB, NB);
  }

  template<typename MatsT>
  cqmatrix::Matrix<MatsT> Orthogonalization<MatsT>::ortho2nonortho(cqmatrix::Matrix<MatsT> & mat) const {
    if( not overlap ) CErr("Overlap has not been initialized and computed");
    if( backwardTrans->dimension() != mat.dimension() ) CErr("Matrices are not the same dimension in ortho2nonortho");

    size_t NB = backwardTrans->dimension();
    return mat.transform('N', backwardTrans->pointer(), NB, NB);
  }

  //======================
  // PauliSpinor Operators
  //======================
  template<typename MatsT>
  cqmatrix::PauliSpinorMatrices<MatsT> Orthogonalization<MatsT>::nonortho2ortho(cqmatrix::PauliSpinorMatrices<MatsT> & mat) const {
    if( not overlap ) CErr("Overlap has not been initialized and computed");
    if( forwardTrans->dimension() != mat.dimension() ) CErr("Matrices are not the same dimension in nonortho2ortho");

    size_t NB = forwardTrans->dimension();
    return mat.transform('N', forwardTrans->pointer(), NB, NB);
  }

  template<typename MatsT>
  cqmatrix::PauliSpinorMatrices<MatsT> Orthogonalization<MatsT>::ortho2nonortho(cqmatrix::PauliSpinorMatrices<MatsT> & mat) const {
    if( not overlap ) CErr("Overlap has not been initialized and computed");
    if( backwardTrans->dimension() != mat.dimension() ) CErr("Matrices are not the same dimension in ortho2nonortho");

    size_t NB = backwardTrans->dimension();
    return mat.transform('N', backwardTrans->pointer(), NB, NB);
  }

  //===================
  // Coefficients
  // NOTE: coefficients transform
  //       in the opposite direction
  //       to operators
  //===================
  template<typename MatsT>
  void Orthogonalization<MatsT>::nonortho2orthoCoeffs(cqmatrix::Matrix<MatsT>& mo) const {
    if( not overlap ) CErr("Overlap has not been initialized in nonortho2orthoCoeffs");

    size_t NB = mo.dimension();
    if( forwardTrans->dimension() != NB ) CErr("Matrices are not the same dimension in nonortho2orthoCoeffs");

    MatsT* SCR = CQMemManager::get().template malloc<MatsT>(NB * NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), backwardTrans->pointer(), NB, mo.pointer(), NB, MatsT(0.), SCR, NB);
    SetMat('N', NB, NB, MatsT(1.), SCR, NB, mo.pointer(), NB);
    CQMemManager::get().free(SCR);
  }

  template<typename MatsT>
  void Orthogonalization<MatsT>::ortho2nonorthoCoeffs(cqmatrix::Matrix<MatsT>& mo) const {
    if( not overlap ) CErr("Overlap has not been initialized in nonortho2orthoCoeffs");

    size_t NB = mo.dimension();
    if( forwardTrans->dimension() != NB ) CErr("Matrices are not the same dimension in ortho2nonorthoCoeffs");

    // Loop over components
    MatsT* SCR = CQMemManager::get().template malloc<MatsT>(NB * NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), forwardTrans->pointer(), NB, mo.pointer(), NB, MatsT(0.), SCR, NB);
    SetMat('N', NB, NB, MatsT(1.), SCR, NB, mo.pointer(), NB);
    CQMemManager::get().free(SCR);
  }

  template<typename MatsT>
  void Orthogonalization<MatsT>::nonortho2orthoCoeffs(std::vector<cqmatrix::Matrix<MatsT>>& mo) const {
    for( auto& m : mo )
      nonortho2orthoCoeffs(m);
  }

  template<typename MatsT>
  void Orthogonalization<MatsT>::ortho2nonorthoCoeffs(std::vector<cqmatrix::Matrix<MatsT>>& mo) const {
    for( auto& m : mo )
      ortho2nonorthoCoeffs(m);
  }

  template<typename MatsT>
  void Orthogonalization<MatsT>::nonortho2orthoCoeffs(std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>>&  mo) const {
    for( auto& m : mo )
      nonortho2orthoCoeffs(m);
  }

  template<typename MatsT>
  void Orthogonalization<MatsT>::ortho2nonorthoCoeffs(std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>>& mo) const {
    for( auto& m : mo )
      ortho2nonorthoCoeffs(m);
  }

  //===================
  // Orthogonalize N States
  //===================
  
  template<typename MatsT>
  void Orthogonalization<MatsT> :: orthogonalizeStates( cqmatrix::Matrix<MatsT>& mo, size_t nStates, size_t disp ) const {

    if( not overlap ) CErr("Overlap has not been initialized in orthogonalizeStates");
    if( nStates == 0 ) return;

    size_t NB = mo.dimension();
    if( overlap->dimension() != NB ) CErr("Matrices are not the same dimension in orthogonalizeStates");

    // Set pointer to start of states of interest
    MatsT* moPointer = mo.pointer() + disp*NB;
    
    // Compute Overlap Matrix for nStates
    cqmatrix::Matrix<MatsT> stateOverlap = overlap->transform( 'N', moPointer, nStates, NB);

#if 0
    prettyPrintSmart(std::cout, "MO Overlap", stateOverlap.pointer(),nStates,nStates,nStates);
#endif

    // Compute Orthogonalization Matrices
    Orthogonalization<MatsT> orthoMO(stateOverlap);

    // Transform MO's in place
    MatsT* SCR = CQMemManager::get().template malloc<MatsT>(NB * NB);
    MatsT* transPointer = orthoMO.backwardPointer()->pointer(); 
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, nStates, nStates, MatsT(1.), moPointer, NB, transPointer, nStates, MatsT(0.), SCR, NB);
    SetMat('N', NB, nStates, MatsT(1.), SCR, NB, moPointer, NB);
    CQMemManager::get().free(SCR);

#if 0
    // Test that they are orthogonalized
    // Compute Overlap Matrix for nStates
    cqmatrix::Matrix<MatsT> testOverlap = overlap->transform( 'N', moPointer, nStates, NB);

    // Compute Orthogonalization Matrices
    Orthogonalization<MatsT> testOrtho(stateOverlap);
    
    prettyPrintSmart(std::cout, "Test Overlap Matrix", testOrtho.overlapPointer()->pointer(),nStates,nStates,nStates);
#endif

  }

  // ====================================================================
  // Transform input operators with the gradient of the orthogonalization
  // ====================================================================
  template <typename MatsT>
  void Orthogonalization<MatsT> :: getOrthogonalizationGradients(
    std::vector<cqmatrix::Matrix<MatsT>>& gradOrtho,
    std::vector<cqmatrix::Matrix<MatsT>>& operators)
  {

    if( operators.size() == 0 )
      CErr("No operators to transform in getOrthogonalizationGradients");
    if( operators.size() != gradOrtho.size() )
      CErr("Different number of operators and output matrices in getOrthogonalizationGradients");

    size_t NB = operators[0].dimension();
    size_t nSQ = NB*NB;

    if( orthoType == LOWDIN ) {
      MatsT* sVecs   = CQMemManager::get().malloc<MatsT>(nSQ);
      MatsT* sE      = CQMemManager::get().malloc<MatsT>(NB);
      MatsT* weights = CQMemManager::get().malloc<MatsT>(nSQ);
      MatsT* SCR1    = CQMemManager::get().malloc<MatsT>(nSQ);
      MatsT* SCR2    = CQMemManager::get().malloc<MatsT>(nSQ);

      // Diagonalize the overlap in scratch S = V * s * V**T
      std::copy_n(overlap->pointer(),nSQ,sVecs);
      HermetianEigen('V','U',NB,sVecs,NB,sE);

      if( std::abs( sE[0] ) < 1e-10 )
        CErr("Contracted Basis Set is Linearly Dependent!");


      // Compute weights = (sqrt(si) + sqrt(sj))^-1
      for( auto i = 0; i < NB; i++ )
      for( auto j = 0; j < NB; j++ ) {
          weights[i*NB + j] = MatsT(1.) / (std::sqrt(sE[i]) + std::sqrt(sE[j]));
      }

      // Loop over gradient components
      for( auto iGrad = 0; iGrad < operators.size(); iGrad++ ) {

        //
        // dV/dR = V . (weights x V**T . O . V) . V**T
        //

        // O . V
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(1.),operators[iGrad].pointer(),NB,sVecs,NB,MatsT(0.),SCR1,NB);
        // V**T . O . V
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(1.),sVecs,NB,SCR1,NB,MatsT(0.),SCR2,NB);
        // weights x V**T . O . V
        std::transform(SCR2, SCR2+nSQ, weights, SCR2, std::multiplies<>());
        // (weights x V**T . O . V) . V**T
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
          NB,NB,NB,MatsT(1.),SCR2,NB,sVecs,NB,MatsT(0.),SCR1,NB);
        // dV/dR = V . (weights x V**T . O . V) . V**T
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(1.),sVecs,NB,SCR1,NB,MatsT(0.),gradOrtho[iGrad].pointer(),NB);

      } 

      CQMemManager::get().free(sVecs, sE, weights, SCR1, SCR2);

    }
    else if( orthoType == CHOLESKY ) {
      CErr("Cholesky orthogonalization gradients not yet implemented");
    }

  }
}
