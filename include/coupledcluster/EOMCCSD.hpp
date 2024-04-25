/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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

#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>
#include <util/math.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <functional>
#include <util/timer.hpp>
#include <coupledcluster/EOMCCSDVector.hpp>

#define EOMCC_CHECK_LAMBDA

namespace ChronusQ{

  template <typename MatsT>
  void swapVectorToFirst(size_t groundIndex, MatsT* M, size_t ldm) {
    MatsT* tmpVec = CQMemManager::get().malloc<MatsT>(ldm);

    SetMat('N', ldm, 1, 1.0, M + groundIndex * ldm, ldm, tmpVec, ldm);
    while (groundIndex > 0) {
      SetMat('N', ldm, 1, 1.0, M + (groundIndex - 1) * ldm, ldm, M + groundIndex * ldm, ldm);
      groundIndex--;
    }
    SetMat('N', ldm, 1, 1.0, tmpVec, ldm, M, ldm);

    CQMemManager::get().free(tmpVec);
  }

  template <typename MatsT>
  void findTrueGroundStateEOMCCEigen(size_t Hbar_dim_w0,
                                     MatsT* theta_w0, MatsT* VL_w0, MatsT* VR_w0, double e_conv) {

    size_t groundIndex = 0;
    double absGroundL0 = std::abs(VL_w0[0]);
    double secondLargestL0 = std::abs(VL_w0[Hbar_dim_w0]);
    if (secondLargestL0 > absGroundL0) {
      groundIndex = 1;
      std::swap(secondLargestL0, absGroundL0);
    }
    for (size_t i = 2; i < Hbar_dim_w0; i++) {
      if (std::abs(VL_w0[i * Hbar_dim_w0]) > absGroundL0) {
        groundIndex = i;
        absGroundL0 = std::abs(VL_w0[i * Hbar_dim_w0]);
      }
    }

    std::cout << "Found true ground state at index " << groundIndex
              << " with abs(L0) = " << absGroundL0
              <<", the next largest abs(L0) = " << secondLargestL0 << std::endl;

    if (groundIndex != 0) {

      std::cout << "Swap ground state to index 0 ..." << std::endl;

      swapVectorToFirst(groundIndex, theta_w0, 1);
      swapVectorToFirst(groundIndex, VL_w0, Hbar_dim_w0);
      swapVectorToFirst(groundIndex, VR_w0, Hbar_dim_w0);

      std::cout << "Swap finished" << std::endl;

    }
  }

  template <typename _F>
  void biOrthoNormalize(size_t N, size_t nR, RawVectors<_F> &VL, RawVectors<_F> &VR) {

    // Biorthonomalize VR_ and VL_
    // VL_^\dagger*VR_ = P*L*U
    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    // VR_ = VR_ * U^{-1}
    std::vector<int64_t> IPIV(nR);
    cqmatrix::Matrix<_F> LUMat(nR);

//    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
//               nR,nR,N,_F(1.),VL,N,VR,N,_F(0.),LUMat,nR);
    VL.dot_product(0, VR, 0, nR, nR, LUMat.pointer(), nR, false);
//    prettyPrintSmart(std::cout,"LUMat:  ",LUMat,nR,nR,nR);

    if (MPIRank() == 0) {
      lapack::getrf(nR, nR, LUMat.pointer(), nR, IPIV.data());
//    prettyPrintSmart(std::cout,"LUMat after LU:  ",LUMat,nR,nR,nR);


      // Compute the inverse of lower and upper triangular matrices in-place
      lapack::trtri(lapack::Uplo::Upper, lapack::Diag::NonUnit, nR, LUMat.pointer(), nR);
      lapack::trtri(lapack::Uplo::Lower, lapack::Diag::Unit, nR, LUMat.pointer(), nR);

      // VR_ = VR_ * U^{-1}
      blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper,
                 blas::Op::NoTrans, lapack::Diag::NonUnit, N, nR, _F(1.0), LUMat.pointer(), nR, VR.getPtr(), N);

      // Apply P^{T}
      //IPIV represents elementary permutation matrices, whose transpose are themselves.
      //P = P1 * P2 * ... * Pn
      // P^{T} = Pn^{T} * ... * P2^{T} * P1^{T}
      // = Pn * ... * P2 * P1
      Eigen::Map<
          Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      > VLMap(VL.getPtr(),N,nR);

      for (int i = 0; i < nR; i++){
        if (i != IPIV[i] - 1){
//          std::cout << i << "   "  << IPIV[i] - 1 << std::endl;
          VLMap.col(IPIV[i] - 1).swap(VLMap.col(i));
        }
      }


      // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
      blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Lower,
                 blas::Op::Trans, lapack::Diag::Unit, N, nR, _F(1.0), LUMat.pointer(), nR, VL.getPtr(), N);

//    prettyPrintSmart(std::cout,"New VL:  ",VL,N,nR,N);
//    prettyPrintSmart(std::cout,"New VR:  ",VR,N,nR,N);
    } // END ROOT_ONLY Section

    for (size_t i = 0; i < nR; i++) {
//      double norm = blas::nrm2(N, VR + i * N, 1);
//      blas::scal(N, 1.0/norm, VR + i * N, 1);
//      blas::scal(N, norm, VL + i * N, 1);
      double norm = VR.norm2F(i, 1);
      VR.scale(1.0/norm, i, 1);
      VL.scale(norm, i, 1);
    }

//    prettyPrintSmart(std::cout,"New VL after scale:  ",VL,N,nR,N);
//    prettyPrintSmart(std::cout,"New VR after scale:  ",VR,N,nR,N);

  }

  template <typename _F>
  void biOrthoNormalize(size_t nR, EOMCCSDVectorSet<_F> &VL, EOMCCSDVectorSet<_F> &VR) {

    // Biorthonomalize VR_ and VL_
    // VL_^\dagger*VR_ = P*L*U
    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    // VR_ = VR_ * U^{-1}
    std::vector<int64_t> IPIV(nR);
    cqmatrix::Matrix<_F> LUMat(nR);

    //    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
    //               nR,nR,N,_F(1.),VL,N,VR,N,_F(0.),LUMat,nR);
    VL.dot_product(0, VR, 0, nR, nR, LUMat.pointer(), nR, false);
    //    prettyPrintSmart(std::cout,"LUMat:  ",LUMat,nR,nR,nR);

    if (MPIRank() == 0) {
      lapack::getrf(nR, nR, LUMat.pointer(), nR, IPIV.data());
      //    prettyPrintSmart(std::cout,"LUMat after LU:  ",LUMat,nR,nR,nR);


      // Compute the inverse of lower and upper triangular matrices in-place
      lapack::trtri(lapack::Uplo::Upper, lapack::Diag::NonUnit, nR, LUMat.pointer(), nR);
      lapack::trtri(lapack::Uplo::Lower, lapack::Diag::Unit, nR, LUMat.pointer(), nR);
    }

    TA::get_default_world().gop.fence();
    LUMat.broadcast();
    MPIBCast(IPIV.data(), nR, 0, MPI_COMM_WORLD);

    cqmatrix::Matrix<_F> UMat(LUMat);
    UMat.setTriangle(blas::Uplo::Lower, 0.0, false);
    cqmatrix::Matrix<_F> LMat(LUMat);
    LMat.setTriangle(blas::Uplo::Upper, 0.0, true, 1.0);

    // VR_ = VR_ * U^{-1}
//    blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper,
//               blas::Op::NoTrans, lapack::Diag::NonUnit, N, nR, _F(1.0), LUMat, nR, VR.getPtr(), N);
    EOMCCSDVectorSet<_F> Vcopy(VR);
    Vcopy.multiply_matrix(0, blas::Op::NoTrans, nR, nR, _F(1.0), UMat.pointer(), nR, _F(0.0), VR, 0);

    // Apply P^{T}
    //IPIV represents elementary permutation matrices, whose transpose are themselves.
    //P = P1 * P2 * ... * Pn
    // P^{T} = Pn^{T} * ... * P2^{T} * P1^{T}
    // = Pn * ... * P2 * P1
    Eigen::Map<
        Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    > LMap(LMat.pointer(),nR,nR);

    for (int i = nR; i > 0; i--){
      if (i != IPIV[i - 1]){
//        std::cout << i - 1 << "   "  << IPIV[i - 1] - 1 << std::endl;
        LMap.col(IPIV[i - 1] - 1).swap(LMap.col(i - 1));
      }
    }

    LMat = LMat.T();

    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    Vcopy.set_data(0, nR, VL, 0);
    Vcopy.multiply_matrix(0, blas::Op::NoTrans, nR, nR, _F(1.0), LMat.pointer(), nR, _F(0.0), VL, 0);

    //    prettyPrintSmart(std::cout,"New VL:  ",VL,N,nR,N);
    //    prettyPrintSmart(std::cout,"New VR:  ",VR,N,nR,N);

    for (size_t i = 0; i < nR; i++) {
      //      double norm = blas::nrm2(N, VR + i * N, 1);
      //      blas::scal(N, 1.0/norm, VR + i * N, 1);
      //      blas::scal(N, norm, VL + i * N, 1);
      double norm = VR.norm2F(i, 1);
      VR.scale(1.0/norm, i, 1);
      VL.scale(norm, i, 1);
    }

    //    prettyPrintSmart(std::cout,"New VL after scale:  ",VL,N,nR,N);
    //    prettyPrintSmart(std::cout,"New VR after scale:  ",VR,N,nR,N);

  }


  template <typename MatsT, typename IntsT>
  inline double EOMCCSD<MatsT,IntsT>::signD(size_t &a, size_t &b, size_t &i, size_t &j) const {
    if (a == b or i == j) {
      a = 0;
      b = 0;
      i = 0;
      j = 0;
      return 0.0;
    }
    double sign = 1.0;
    if (a > b) {
      std::swap(a,b);
      sign *= -1.0;
    }
    if (i > j) {
      std::swap(i,j);
      sign *= -1.0;
    }
    return sign;
  }

  template <typename MatsT, typename IntsT>
  EOMCCSD<MatsT,IntsT>::EOMCCSD(const SafeFile &savFile,
                                CCIntermediates<MatsT> &intermediates,
                                const EOMSettings &eomSettings,
                                const CoupledClusterSettings &ccSettings):
      savFile_(savFile),
      intermediates_(intermediates), eomSettings(eomSettings),
      ccSettings_(ccSettings),
      vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
      T1_(intermediates.T->oneBody()), T2_(intermediates.T->twoBody()),
      Lg_(*intermediates.Lg),
      fockMatrix_ta(intermediates.fockMatrix),
      muMatrix(intermediates.muMatrix),
      antiSymMoints(intermediates.antiSymMoInts),
      tau(intermediates.tau),
      F_ae(intermediates.F_ae),
      F_mi(intermediates.F_mi),
      F_me(intermediates.F_me),
      W_mnij(intermediates.W_mnij),
      W_abef(intermediates.W_abef),
      W_mbej(intermediates.W_mbej),
      W_mnie(intermediates.W_mnie),
      W_amef(intermediates.W_amef),
      W_mbij(intermediates.W_mbij),
      W_abei(intermediates.W_abei),
      G_ae(intermediates.G_ae),
      G_mi(intermediates.G_mi),
      D_ai(intermediates.D_ai),
      D_abij(intermediates.D_abij) {

    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();
    nOVshift_ = nV * nO;
    nO2shift_ = nO * (nO - 1) / 2;
    nV2shift_ = nV * (nV - 1) / 2;

    Hbar_dim = nOVshift_ + nO2shift_ * nV2shift_;

    assignCVSIndices();

  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::run() {
    initilizeEOMCCSD();
    formEOMIntermediates();
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::full_diagonalization() {

    theta = CQMemManager::get().malloc<MatsT>(Hbar_dim);
    RawVectors<MatsT> VL(MPI_COMM_WORLD, Hbar_dim, Hbar_dim);
    RawVectors<MatsT> VR(MPI_COMM_WORLD, Hbar_dim, Hbar_dim);

    std::cout << " Start building the full matrix " << std::endl;

    auto beginBuild = tick();
    cqmatrix::Matrix<MatsT> fullMat(buildHbarCVS(false));
    std::cout << "buildHbar spend " << tock(beginBuild) << " s." << std::endl;

//    fullMat.output(std::cout, "Hbar", true);

//#define TEST_AGAINST_SIGMA
#ifdef TEST_AGAINST_SIGMA
    MatsT* diag = CQMemManager::get().malloc<MatsT>(Hbar_dim);
    beginBuild = tick();
    buildDiag(diag);
    std::cout << "buildDiag spend " << tock(beginBuild) << " s." << std::endl;

    MatsT* fullMat2 = CQMemManager::get().malloc<MatsT>(Hbar_dim * Hbar_dim);
    beginBuild = tick();
    buildHbar_sigma(fullMat2, false);
    std::cout << "buildHbar_sigma spend " << tock(beginBuild) << " s." << std::endl;

//    prettyPrintSmart(std::cout, "Hbar Ref", fullMat, Hbar_dim, Hbar_dim, Hbar_dim);

    std::cout << " Diagonal elements: " << std::endl;

    for (size_t i = 0; i < Hbar_dim; i++)
      std::cout << i << " : " << diag[i] << " ; " << fullMat2[i * Hbar_dim + i] << " diff " << diag[i] - fullMat(i,i) << std::endl;

    blas::axpy(Hbar_dim * Hbar_dim, -1.0, fullMat.pointer(), 1, fullMat2, 1);

    double errorNorm = blas::nrm2(Hbar_dim * Hbar_dim, fullMat2, 1);

//    prettyPrintSmart(std::cout, "Hbar Error", fullMat2, Hbar_dim, Hbar_dim, Hbar_dim);

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "Error Norm = " << errorNorm << std::endl;

    CQMemManager::get().free(diag, fullMat2);
#endif

    std::cout << " Finished building the full matrix " << std::endl;

    if (MPIRank() == 0) {
      SetLAThreads(GetNumThreads());

    std::cout << " Start Full diagonalization, Hbar_dim = "  << Hbar_dim << std::endl;

    beginBuild = tick();
    GeneralEigen('V', 'V', Hbar_dim, fullMat.pointer(), Hbar_dim, theta, VL.getPtr(), Hbar_dim, VR.getPtr(), Hbar_dim);
    std::cout << "GeneralEigen spend " << tock(beginBuild) << " s." << std::endl;

    std::cout << " Eigenvalues from the full matrix: " << std::endl;
    std::cout << std::fixed << std::setprecision(12);
    for ( auto i = 0; i < Hbar_dim; i++)
      std::cout << i << " EigV: " << theta[i] <<  std::endl;

//    prettyPrintSmart(std::cout, "VL", VL, Hbar_dim, Hbar_dim, Hbar_dim);
//    prettyPrintSmart(std::cout, "VR", VR, Hbar_dim, Hbar_dim, Hbar_dim);

//    MatsT* VLVR = CQMemManager::get().malloc<MatsT>(Hbar_dim * Hbar_dim);
//    VL.dot_product(0, VR, 0, Hbar_dim, Hbar_dim, VLVR, Hbar_dim);
//    prettyPrintSmart(std::cout, "VLVR", VLVR, Hbar_dim, Hbar_dim, Hbar_dim);
//    CQMemManager::get().free(VLVR);
      SetLAThreads(1);
    } // ROOT_ONLY section

    if (not eomSettings.oscillator_strength)
      return;




    // Building Hbar matrix including H0S H0D blocks
    size_t Hbar_dim_w0 = Hbar_dim + 1;
    MatsT* theta_w0 = nullptr;
    RawVectors<MatsT> VL_w0(MPI_COMM_WORLD, Hbar_dim_w0, Hbar_dim_w0);
    RawVectors<MatsT> VR_w0(MPI_COMM_WORLD, Hbar_dim_w0, Hbar_dim_w0);

    theta_w0 = CQMemManager::get().malloc<MatsT>(Hbar_dim_w0);

    std::cout << " Start building the full matrix " << std::endl;

    auto beginBuild_w0 = tick();
    cqmatrix::Matrix<MatsT> fullMat_w0(buildHbarCVS(true));
    std::cout << "buildHbar_w0 spend " << tock(beginBuild_w0) << " s." << std::endl;

//    fullMat_w0.output(std::cout, "Hbar_w0", true);

    std::cout << " Finished building the full matrix " << std::endl;

    if (eomSettings.save_hamiltonian and MPIRank() == 0) {
      savFile_.safeWriteData("/CC/HAMILTONIAN", fullMat_w0.pointer(), {Hbar_dim_w0,Hbar_dim_w0});
    }

    // Diagonalize the big Hbar matrix
    cqmatrix::Matrix<MatsT> fullMat_w0_copy(fullMat_w0);
    if (MPIRank() == 0) {
      SetLAThreads(GetNumThreads());
    std::cout << " Start Full diagonalization, Hbar_dim_w0 = "  << Hbar_dim_w0 << std::endl;

    beginBuild_w0 = tick();
    GeneralEigen('V', 'V', Hbar_dim_w0, fullMat_w0.pointer(), Hbar_dim_w0, theta_w0, VL_w0.getPtr(), Hbar_dim_w0, VR_w0.getPtr(), Hbar_dim_w0);
    findTrueGroundStateEOMCCEigen(Hbar_dim_w0, theta_w0, VL_w0.getPtr(), VR_w0.getPtr(), ccSettings_.eConv);
    std::cout << "GeneralEigen spend " << tock(beginBuild_w0) << " s." << std::endl;

    std::cout << " Eigenvalues from the full matrix: " << std::endl;
    std::cout << std::fixed << std::setprecision(12);
    for ( auto i = 0; i < Hbar_dim_w0; i++)
      std::cout << i << " EigV: " << theta_w0[i] <<  std::endl;
//    prettyPrintSmart(std::cout, "VL_w0", VL_w0, Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);
//    prettyPrintSmart(std::cout, "VR_w0", VR_w0, Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);
      SetLAThreads(1);
    }// END ROOT_ONLY section

#ifdef EOMCC_CHECK_LAMBDA
    if (not eomSettings.doCVS()) {

      // Begin Lambda check
      // Convert L1 and L2 amplitudes to a vector
      MatsT* L0p1 = CQMemManager::get().malloc<MatsT>(this->getHbarDim(true));
      Lg_.toRaw(L0p1, *this, true);

      if (MPIRank() == 0) {
      MatsT *L0raw = L0p1 + 1;

//      prettyPrintSmart(std::cout, "L0", L0raw, 1, Hbar_dim, 1);


      // Check if L0 is an eigenvector of Hbar with eigenvalue 0
      MatsT* L0p1Hbar = CQMemManager::get().malloc<MatsT>(Hbar_dim + 1);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 1, Hbar_dim_w0, Hbar_dim_w0, 1.0, L0p1, 1, fullMat_w0_copy.pointer(), Hbar_dim_w0, 0.0, L0p1Hbar, 1);
//      prettyPrintSmart(std::cout, "L0p1Hbar", L0p1Hbar, 1, Hbar_dim_w0, 1);
      std::cout << "L0 eig check error = " << std::scientific << std::setprecision(4)
                << blas::nrm2(Hbar_dim_w0, L0p1Hbar, 1) << std::endl;


      // Get L0 from VL_w0
      MatsT* L0_w0 = CQMemManager::get().malloc<MatsT>(Hbar_dim);
      SetMat('R', Hbar_dim, 1, 1.0/VL_w0.getPtr()[0], VL_w0.getPtr() + 1, Hbar_dim_w0, L0_w0, Hbar_dim);
//      prettyPrintSmart(std::cout, "L0_w0", L0_w0, 1, Hbar_dim, 1);

      blas::axpy(Hbar_dim, -1.0, L0raw, 1, L0_w0, 1);
      std::cout << "L0 check diff norm = " << std::scientific << std::setprecision(4)
                << blas::nrm2(Hbar_dim, L0_w0, 1) << std::endl;

      MatsT *L0VR = CQMemManager::get().malloc<MatsT>(Hbar_dim);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 1, Hbar_dim, Hbar_dim, 1.0, L0raw, 1, VR.getPtr(), Hbar_dim, 0.0, L0VR, 1);
//      prettyPrintSmart(std::cout, "L0VR", L0VR, 1, Hbar_dim, 1);

      // Get Hbar 0S and 0D blocks, this has to be done before diagonalize Hbar
      // because ZGEEV destroys Hbar
      MatsT *H0 = CQMemManager::get().malloc<MatsT>(Hbar_dim);
      SetMat('N', 1, Hbar_dim, 1.0, fullMat_w0_copy.pointer() + Hbar_dim_w0, Hbar_dim_w0, H0, 1);
//      prettyPrintSmart(std::cout, "H0", H0, 1, Hbar_dim, 1);

      // Solve for R0
      MatsT *H0VR = CQMemManager::get().malloc<MatsT>(Hbar_dim);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 1, Hbar_dim, Hbar_dim, 1.0, H0, 1, VR.getPtr(), Hbar_dim, 0.0, H0VR, 1);
//      prettyPrintSmart(std::cout, "H0VR", H0VR, 1, Hbar_dim, 1);
      for (size_t i = 0; i < Hbar_dim; i++) {
        H0VR[i] /= theta[i];
      }
//      prettyPrintSmart(std::cout, "H0VR/theta", H0VR, 1, Hbar_dim, 1);

      blas::axpy(Hbar_dim, 1.0, H0VR, 1, L0VR, 1);
      std::cout << "L0 Biorthonormality check error norm = " << std::scientific << std::setprecision(4)
                << blas::nrm2(Hbar_dim, L0VR, 1) << std::endl;

//      for (size_t i = 0; i < Hbar_dim; i++) {
//        H0VR[i] /= std::sqrt(1.0 + SmartConj(H0VR[i]) * H0VR[i]);
//      }
//      prettyPrintSmart(std::cout, "H0VR/theta after normalization", H0VR, 1, Hbar_dim, 1);

      CQMemManager::get().free(L0p1, L0p1Hbar, L0_w0, L0VR, H0, H0VR);
      } // END ROOT_ONLY section
      else {
        CQMemManager::get().free(L0p1);
      }
    }
#endif


    MatsT* VLVR_w0 = CQMemManager::get().malloc<MatsT>(Hbar_dim_w0 * Hbar_dim_w0);
    VL_w0.dot_product(0, VR_w0, 0, Hbar_dim_w0, Hbar_dim_w0, VLVR_w0, Hbar_dim_w0);
//    prettyPrintSmart(std::cout, "VLVR_w0", VLVR_w0, Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);

    MatsT* VRVR_w0 = CQMemManager::get().malloc<MatsT>(Hbar_dim_w0 * Hbar_dim_w0);
    VR_w0.dot_product(0, VR_w0, 0, Hbar_dim_w0, Hbar_dim_w0, VRVR_w0, Hbar_dim_w0);
//    prettyPrintSmart(std::cout, "VRVR_w0", VRVR_w0, Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);

    // biOrthoNormalize and check
    std::cout << "Start biorthonormalization ..." << std::endl;
    VL_w0.conjugate();
    biOrthoNormalize(Hbar_dim_w0, Hbar_dim_w0, VL_w0, VR_w0);
    VL_w0.conjugate();

//    VL_w0.print(std::cout, "New VL after biOrthoNormalize");
//    VR_w0.print(std::cout, "New VR after biOrthoNormalize");

    MatsT *SCR = CQMemManager::get().malloc<MatsT>(Hbar_dim_w0*Hbar_dim_w0);
    VL_w0.dot_product(0, VR_w0, 0, Hbar_dim_w0, Hbar_dim_w0, SCR, Hbar_dim_w0);
//    prettyPrintSmart(std::cout,"New VLVR",SCR,Hbar_dim_w0,Hbar_dim_w0,Hbar_dim_w0);
    for (size_t i = 0; i < Hbar_dim_w0; i++)
      SCR[(1+Hbar_dim_w0) * i] -= 1.0;
    std::cout << "Biorthonormality check error norm = " << std::scientific << std::setprecision(4)
              << blas::nrm2(Hbar_dim_w0*Hbar_dim_w0, SCR, 1) << std::endl;

    VR_w0.dot_product(0, VR_w0, 0, Hbar_dim_w0, Hbar_dim_w0, SCR, Hbar_dim_w0);
//    prettyPrintSmart(std::cout,"New VRVR",SCR,Hbar_dim_w0,Hbar_dim_w0,Hbar_dim_w0);
    for (size_t i = 0; i < Hbar_dim_w0; i++)
      SCR[(1+Hbar_dim_w0) * i] -= 1.0;
    std::cout << "R normality check error norm = " << std::scientific << std::setprecision(4)
              << blas::nrm2(Hbar_dim_w0, SCR, 1+Hbar_dim_w0) << std::endl;


    if (MPIRank() == 0) {
    for (size_t i = 0; i < Hbar_dim_w0; i++) {
      SetMat('N', Hbar_dim_w0, 1, theta_w0[i], VR_w0.getPtr(i), Hbar_dim_w0, SCR + i * Hbar_dim_w0, Hbar_dim_w0);
    }
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
               Hbar_dim_w0,Hbar_dim_w0,Hbar_dim_w0,
               MatsT(1.),fullMat_w0_copy.pointer(),Hbar_dim_w0,VR_w0.getPtr(),Hbar_dim_w0,MatsT(-1.),SCR,Hbar_dim_w0);
    std::cout << "Right eigenvector check error norm = " << std::scientific << std::setprecision(4)
              << blas::nrm2(Hbar_dim_w0*Hbar_dim_w0, SCR, 1) << std::endl;

    for (size_t i = 0; i < Hbar_dim_w0; i++) {
      SetMat('C', Hbar_dim_w0, 1, theta_w0[i], VL_w0.getPtr(i), Hbar_dim_w0, SCR + i, Hbar_dim_w0);
    }
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
               Hbar_dim_w0,Hbar_dim_w0,Hbar_dim_w0,
               MatsT(1.),VL_w0.getPtr(),Hbar_dim_w0,fullMat_w0_copy.pointer(),Hbar_dim_w0,MatsT(-1.),SCR,Hbar_dim_w0);
    std::cout << "Left eigenvector check error norm = " << std::scientific << std::setprecision(4)
              << blas::nrm2(Hbar_dim_w0*Hbar_dim_w0, SCR, 1) << std::endl;
    }



    // Begin evaluate density matrices
    if (MPIRank() == 0)
      IMatCopy('R', Hbar_dim_w0, Hbar_dim_w0, 1.0, VL_w0.getPtr(), Hbar_dim_w0, Hbar_dim_w0);

    // Convert Lg from the first vector of VL_w0
    MatsT* rawPtr = const_cast<MatsT*>(VL_w0.getPtr(0));
    if (MPIRank() != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(this->getHbarDim(true));
    TA::get_default_world().gop.template broadcast(rawPtr, this->getHbarDim(true), 0);
    Lg_.fromRaw(rawPtr, *this, true);
    if (MPIRank() != 0)
      CQMemManager::get().free(rawPtr);

    L_ = std::make_shared<EOMCCSDVectorSet<MatsT>>(vLabel_, oLabel_, eomSettings.nroots);
    std::dynamic_pointer_cast<EOMCCSDVectorSet<MatsT>>(L_)->fromRaw(MPI_COMM_WORLD, VL_w0, *this, true, 0, 1, eomSettings.nroots);
    R_ = std::make_shared<EOMCCSDVectorSet<MatsT>>(vLabel_, oLabel_, eomSettings.nroots);
    std::dynamic_pointer_cast<EOMCCSDVectorSet<MatsT>>(R_)->fromRaw(MPI_COMM_WORLD, VR_w0, *this, true, 0, 1, eomSettings.nroots);
    initilizeDensity();
    std::vector<double> oscStrength;
    std::vector<double> excitationE;
    std::cout << "----------------------------------------------" << std::endl;
    for (size_t i = 0; i < eomSettings.nroots ; i++){
      MatsT f = calcOscillatorStrength(i);
      std::cout << " Excited State: " << i << "  E = " << theta[i] << " Eh, "<< " f = " << f << std::endl;
      oscStrength.push_back(std::real(f));
      excitationE.push_back(std::real(theta[i]));
    }

    TA::get_default_world().gop.fence();
    // Write data to bin file
    if (savFile_.exists()) {
      savFile_.safeWriteData("/CC/EXCITATION_ENERGIES",excitationE.data(), {eomSettings.nroots});
      savFile_.safeWriteData("/CC/OSCILLATOR_STRENGTHS", oscStrength.data(), {eomSettings.nroots});
    }

    CQMemManager::get().free(theta_w0, VLVR_w0, VRVR_w0, SCR);

  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::initilizeEOMCCSD() {
    TAManager &TAmanager = TAManager::get();

    if (not W_mnie.is_initialized()){
      W_mnie = TAmanager.malloc<MatsT>("ooov");
    }

    if (not W_amef.is_initialized()){
      W_amef = TAmanager.malloc<MatsT>("vovv");
    }

    if (not W_mbij.is_initialized()){
      W_mbij = TAmanager.malloc<MatsT>("ovoo");
    }

    if (not W_abei.is_initialized()){
      W_abei = TAmanager.malloc<MatsT>("vvvo");
    }

  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formF_ae() {
    F_ae("a,e") -= 0.5 * this->T1_("a,m") * F_me("m,e");
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formF_mi() {
    F_mi("m,i") += 0.5 * this->T1_("e,i") * F_me("m,e");
  }

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_mnij() {
    W_mnij("m,n,i,j") += 0.25 * tau("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_abef() {
    W_abef("a,b,e,f") += 0.25 * tau("a,b,m,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_mbej() {
    W_mbej("m,b,e,j") -= 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  } 

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_mnie() {
    W_mnie("m,n,i,e") = - conj(this->antiSymMoints["vooo"]("e,i,m,n")) + this->T1_("f,i") * conj(this->antiSymMoints["vvoo"]("f,e,m,n"));
  } 

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_amef() {
    W_amef("a,m,e,f") = conj(this->antiSymMoints["vvvo"]("e,f,a,m")) - this->T1_("a,n") * conj(this->antiSymMoints["vvoo"]("e,f,n,m"));
  } 
  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_mbij() {
    W_mbij("m,b,i,j") = - this->antiSymMoints["vooo"]("b,m,i,j") - F_me("m,e") * this->T2_("b,e,i,j");
    W_mbij("m,b,i,j") += - this->T1_("b,n") * W_mnij("m,n,i,j");
    W_mbij("m,b,i,j") += - 0.5 * conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * tau("e,f,i,j");
    W_mbij("m,b,i,j") += - conj(this->antiSymMoints["vooo"]("e,i,m,n")) * this->T2_("b,e,j,n");
    W_mbij("m,b,i,j") += conj(this->antiSymMoints["vooo"]("e,j,m,n")) * this->T2_("b,e,i,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,j") = - this->antiSymMoints["vovo"]("b,m,e,j") - this->T2_("b,f,n,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    W_mbij("m,b,i,j") += this->T1_("e,i") * tmp("m,b,e,j");
    W_mbij("m,b,i,j") += - this->T1_("e,j") * tmp("m,b,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
  } 
  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formW_abei() {
    W_abei("a,b,e,i") = this->antiSymMoints["vvvo"]("a,b,e,i") - F_me("m,e") * this->T2_("a,b,m,i");
    W_abei("a,b,e,i") += this->T1_("f,i") * W_abef("a,b,e,f");
    W_abei("a,b,e,i") +=  0.5 * conj(this->antiSymMoints["vooo"]("e,i,m,n")) * tau("a,b,m,n");
    W_abei("a,b,e,i") += conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * this->T2_("a,f,m,i");
    W_abei("a,b,e,i") += -conj(this->antiSymMoints["vvvo"]("e,f,a,m")) * this->T2_("b,f,m,i");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,i") = - this->antiSymMoints["vovo"]("b,m,e,i") - this->T2_("b,f,n,i") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    W_abei("a,b,e,i") += - this->T1_("a,m") * tmp("m,b,e,i");
    W_abei("a,b,e,i") +=  this->T1_("b,m") * tmp("m,a,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
  }  

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formEOMIntermediates() {
    formF_ae();
    formF_mi();
    formW_mnij();
    formW_abef();
    formW_mbej();
    formW_mnie();
    formW_amef();
    formW_mbij();
    formW_abei();
  } 

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formR1_tilde(const TArray &R1, const TArray &R2, TArray &tildeR1) const {
    tildeR1("a,i") = F_ae("a,c") * R1("c,i");
    tildeR1("a,i") += - F_mi("k,i") * R1("a,k");
    tildeR1("a,i") += F_me("k,c") * R2("a,c,i,k");
    tildeR1("a,i") += W_mbej("k,a,c,i") * R1("c,k");
    tildeR1("a,i") += 0.5 * W_amef("a,k,c,d") * R2("c,d,i,k");
    tildeR1("a,i") += - 0.5 * W_mnie("k,l,i,d") * R2("a,d,k,l");
  }  

  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::formR2_tilde(const TArray &R1, const TArray &R2, TArray &tildeR2) const {
    tildeR2("a,b,i,j") = F_ae("b,e") * R2("a,e,i,j");
    tildeR2("a,b,i,j") += - F_ae("a,e") * R2("b,e,i,j");
    tildeR2("a,b,i,j") += - F_mi("k,j") * R2("a,b,i,k");
    tildeR2("a,b,i,j") += F_mi("k,i") * R2("a,b,j,k");
    tildeR2("a,b,i,j") += 0.5 * W_mnij("k,l,i,j") * R2("a,b,k,l");

    tildeR2("a,b,i,j") += 0.5 * W_abef("a,b,e,f") * R2("e,f,i,j");
    tildeR2("a,b,i,j") += W_mbej("k,b,c,j") * R2("a,c,i,k");
    tildeR2("a,b,i,j") += - W_mbej("k,a,c,j") * R2("b,c,i,k");
    tildeR2("a,b,i,j") += - W_mbej("k,b,c,i") * R2("a,c,j,k");
    tildeR2("a,b,i,j") += W_mbej("k,a,c,i") * R2("b,c,j,k");
    tildeR2("a,b,i,j") += W_abei("a,b,c,j") * R1("c,i"); // Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += - W_abei("a,b,c,i") * R1("c,j");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += -W_mbij("k,a,j,i") * R1("b,k");// Sign is different between Tianyuan's and literature
    tildeR2("a,b,i,j") += W_mbij("k,b,j,i") * R1("a,k");// Sign is different between Tianyuan's and literature

    TAManager &TAmanager = TAManager::get();
    //[Asthana:2019:4102] Break equations A8, A9 and A10 to save computational cost
    TArray tmp1 = TAmanager.malloc<MatsT>("vv");
    tmp1("b,d") = W_amef("b,k,d,c") * R1("c,k");
    tildeR2("a,b,i,j") += tmp1("b,d") * this->T2_("a,d,i,j");
    tildeR2("a,b,i,j") += - tmp1("a,d") * this->T2_("b,d,i,j");
    TArray tmp2 = TAmanager.malloc<MatsT>("oo");
    tmp2("l,j") = W_mnie("l,k,j,c") * R1("c,k");
    tildeR2("a,b,i,j") += - tmp2("l,j") * this->T2_("a,b,i,l");
    tildeR2("a,b,i,j") += tmp2("l,i") * this->T2_("a,b,j,l");

    //reuse tmp2 container
    tmp2("l,j") = - 0.5 * conj(antiSymMoints.at("vvoo")("d,c,k,l")) * R2("c,d,j,k");

    tildeR2("a,b,i,j") += tmp2("l,j") * this->T2_("a,b,i,l");
    tildeR2("a,b,i,j") += - tmp2("l,i") * this->T2_("a,b,j,l");
    //reuse tmp1 container
    tmp1("b,e") = 0.5 * conj(antiSymMoints.at("vvoo")("e,d,k,l")) * R2("b,d,k,l");

    tildeR2("a,b,i,j") += - tmp1("b,e") * this->T2_("a,e,i,j");

    tildeR2("a,b,i,j") += tmp1("a,e") * this->T2_("b,e,i,j");

    TAmanager.free("vv", std::move(tmp1));
    TAmanager.free("oo", std::move(tmp2));
  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::buildSigma(const TArray &V1, const TArray &V2,
                                        TArray &HV1, TArray &HV2, EOMCCEigenVecType vecType) const {
    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR1_tilde(V1, V2, HV1);
        formR2_tilde(V1, V2, HV2);
        break;
      case EOMCCEigenVecType::LEFT:
        TAManager &TAmanager = TAManager::get();
        TArray G_ae = TAmanager.malloc<MatsT>("vv");
        updateG_ae(V2, G_ae);
        TArray G_mi = TAmanager.malloc<MatsT>("oo");
        updateG_mi(V2, G_mi);
        formL1_tilde(V1, V2, G_ae, G_mi, HV1, false);
        formL2_tilde(V1, V2, G_ae, G_mi, HV2, false);
        TAmanager.free("vv", std::move(G_ae));
        TAmanager.free("oo", std::move(G_mi));
        break;
    }
  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::buildRightZeroBody(size_t nVec) {
    // Assumes EOMCCSDVectorSet R_ type

    std::shared_ptr<EOMCCSDVectorSet<MatsT>> VR = std::dynamic_pointer_cast<EOMCCSDVectorSet<dcomplex>>(R_);

    for (size_t i = 0; i < nVec; i++) {
      MatsT r0_1 = F_me("i,a").dot(VR->get(i).oneBody()("a,i")).get();
      MatsT r0_2 = conj(antiSymMoints["vvoo"]("a,b,i,j")).dot(VR->get(i).twoBody()("a,b,i,j")).get();
      TA::get_default_world().gop.fence();
      VR->get(i).zeroBody() = (r0_1 + 0.25 * r0_2) / theta[i];
    }
  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::buildHbar_sigma(MatsT * out, bool diagOnly) const {

    TA::get_default_world().gop.fence();

    std::vector<EOMCCSDVector<MatsT>> rs;

    TAManager &TAmanager = TAManager::get();
    TArray tmpR1 = TAmanager.malloc<MatsT>("vo");
    tmpR1("a,i") = 0.0 * tmpR1("a,i");

    TArray tmpR2 = TAmanager.malloc<MatsT>("vvoo");
    tmpR2("a,b,i,j") = 0.0 * tmpR2("a,b,i,j");

    rs.reserve(Hbar_dim);

    for (size_t i = 0; i < Hbar_dim; i++) {
      rs.emplace_back(vLabel_, oLabel_, tmpR1, tmpR2);
    }

    for (auto it = std::begin(tmpR1); it != std::end(tmpR1); ++it) {
      auto idx = it.index();
      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          rs[CVStoCompoundS(x[0],x[1])].oneBody().find(idx).get()[x] = 1.0;
        }
    }

    for (auto it = std::begin(tmpR2); it != std::end(tmpR2); ++it) {
      auto idx = it.index();
      auto &tile = it->get();
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              if (x[0] != x[1] and x[2] != x[3]){
                size_t a = x[0], b = x[1], i = x[2], j = x[3];
                double sign = 1.0;
                if (a > b) {
                  std::swap(a,b);
                  sign *= -1.0;
                }
                if (i > j) {
                  std::swap(i,j);
                  sign *= -1.0;
                }
                rs[CVStoCompoundD(a, b, i, j) + nOVshift_].twoBody().find(idx).get()[x] = sign;
              }
            }
    }

    TA::get_default_world().gop.fence();

    std::vector<EOMCCSDVector<MatsT>> hrs;
    hrs.reserve(Hbar_dim);

    for (size_t i = 0; i < Hbar_dim; i++) {
      hrs.emplace_back(vLabel_, oLabel_);
      const EOMCCSDVector<dcomplex> &Vi = rs[i];
      EOMCCSDVector<dcomplex> &AVi = hrs.back();
      buildSigma(Vi.oneBody(), Vi.twoBody(), AVi.oneBody(), AVi.twoBody(), EOMCCEigenVecType::RIGHT);
      AVi.enforceTwoBodySymmetry();
    }

    for (size_t i = 0; i < Hbar_dim; i++)
      if (diagOnly)
        out[i] = rs[i].dot(hrs[i]);
      else
        for (size_t j = 0; j < Hbar_dim; j++) {
          out[i + j * Hbar_dim] = rs[i].dot(hrs[j]);
        }

    TAmanager.free("vo", std::move(tmpR1));
    TAmanager.free("vvoo", std::move(tmpR2));
  }


  template <typename MatsT, typename IntsT>
  void EOMCCSD<MatsT,IntsT>::buildDiag(MatsT * diag) const {

//    for (size_t a = 0; a < NV; ++a)
//      for (size_t i = 0; i < NO; ++i) {
//        // F_ae[a][a] - F_mi[i][i] + W_mbej[i][a][a][i]
//      }
//
//    for (size_t a = 0; a < NV; ++a)
//      for (size_t b = 0; b < a; ++b)
//        for (size_t i = 0; i < NO; ++i)
//          for (size_t j = 0; j < i; ++j) {
//            //   F_ae[b][b] + F_ae[a][a] - F_mi[j][j] - F_mi[i][i]
//            // + W_mnij[i][j][i][j] + W_abef[a][b][a][b]
//            // + W_mbej[j][b][b][j] + W_mbej[j][a][a][j]
//            // + W_mbej[i][b][b][i] + W_mbej[i][a][a][i]
//            // - sum{l} T2[a][b][i][l] * conj(V[a][b][i][l])
//            // - sum{l} T2[a][b][l][j] * conj(V[a][b][l][j])
//            // - sum{e} T2[a][e][i][j] * conj(V[a][e][i][j])
//            // - sum{e} T2[e][b][i][j] * conj(V[e][b][i][j])
//          }
    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    std::fill_n(diag, Hbar_dim, MatsT(0.0));
    MatsT * diag2 = diag + nOVshift_;

    TA::foreach_inplace( F_ae, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t a = x[0], b = x[0];
        for (size_t i = 0; i < nO; ++i)
          diag[CVStoCompoundS(a,i)] += v;

        for (size_t b = a + 1; b < nV; ++b)
          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }

        for (size_t a = 0; a < b; ++a)
          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }
      }

    });

    TA::foreach_inplace( F_mi, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[1]) return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        x[1] = x[0];
        MatsT v = tile[x];

        size_t i = x[0], j = x[0];
        for (size_t a = 0; a < nV; ++a)
          diag[CVStoCompoundS(a,i)] -= v;

        for (size_t b = 0; b < nV; ++b)
          for (size_t a = 0; a < b; ++a)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[CVStoCompoundD(a, b, i, j)] -= v;
            }

        for (size_t b = 0; b < nV; ++b)
          for (size_t a = 0; a < b; ++a)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] -= v;
            }
      }

    });

    TA::foreach_inplace( W_mbej, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {
      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[3] or lobound[1] != lobound[2])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          x[3] = x[0];
          x[2] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], a = x[1], b = x[1], j = x[0];
          diag[CVStoCompoundS(a,i)] += v;

          for (size_t b = a + 1; b < nV; ++b)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }

          for (size_t b = a + 1; b < nV; ++b)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }

          for (size_t a = 0; a < b; ++a)
            for (size_t j = i + 1; j < nO; ++j) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }

          for (size_t a = 0; a < b; ++a)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_mnij, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
      return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t i = x[0], j = x[1];

          for (size_t b = 0; b < nV; ++b)
            for (size_t a = 0; a < b; ++a) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    TA::foreach_inplace( W_abef, [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] != lobound[2] or lobound[1] != lobound[3] or lobound[0] > lobound[1])
        return;

      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
        for(x[0] = lobound[0]; x[0] < std::min(x[1], static_cast<std::size_t>(upbound[0])); ++x[0]) {
          x[2] = x[0];
          x[3] = x[1];
          MatsT v = tile[x];

          size_t a = x[0], b = x[1];

          for (size_t j = 0; j < nO; ++j)
            for (size_t i = 0; i < j; ++i) {
              diag2[CVStoCompoundD(a, b, i, j)] += v;
            }
        }
      }

    });

    
    TA::foreach_inplace(T2_, antiSymMoints.at("vvoo"),
                        [&diag, &diag2, this, nV, nO](TA::Tensor<MatsT>& T2tile, const TA::Tensor<MatsT>& Vtile) {
      const auto& lobound = T2tile.range().lobound();
      const auto& upbound = T2tile.range().upbound();

      std::vector<std::size_t> x{0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              size_t a = x[0], b = x[1], i = x[2], j = x[3];
              MatsT v = SmartConj(T2tile[x]) * Vtile[x];

              if (i < j) {

                for (size_t a = 0; a < b; ++a)
                  diag2[CVStoCompoundD(a, b, i, j)] -= v;

                for (size_t b = a + 1; b < nV; ++b)
                  diag2[CVStoCompoundD(a, b, i, j)] -= v;

              }

              if (a < b) {

                for (size_t i = 0; i < j; ++i)
                  diag2[CVStoCompoundD(a, b, i, j)] -= v;

                for (size_t j = i + 1; j < nO; ++j)
                  diag2[CVStoCompoundD(a, b, i, j)] -= v;

              }
            }
    });

    TA::get_default_world().gop.fence();

    TA::get_default_world().gop.template reduce(diag, Hbar_dim, std::plus<MatsT>());

  }

  template <typename MatsT, typename IntsT>
  EOMCCSD<MatsT,IntsT>::~EOMCCSD() {
    if (theta) CQMemManager::get().free(theta);

    TAManager &TAmanager = TAManager::get();

    if (Rho_ij) TAmanager.free("oo", std::move(Rho_ij), true);
    if (Rho_ab) TAmanager.free("vv", std::move(Rho_ab), true);
    if (Rho_ia) TAmanager.free("ov", std::move(Rho_ia), true);
    if (Rho_ai) TAmanager.free("vo", std::move(Rho_ai), true);

  }

};
