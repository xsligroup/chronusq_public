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
#include <coupledcluster/MBExpansion.hpp>
#include <itersolver/davidson.hpp>
#include <itersolver.hpp>

#define EOMCC_CHECK_LAMBDA

namespace ChronusQ{
 
  template <typename MatsT>
  void swapVectorToFirst(size_t groundIndex, MatsT* M, size_t ldm);

  template <typename MatsT>
  void findTrueGroundStateEOMCCEigen(size_t Hbar_dim_w0,
                                     dcomplex* theta_w0, MatsT* VL_w0, MatsT* VR_w0, double e_conv);

  template <typename _F>
  void biOrthoNormalize(size_t N, size_t nR, RawVectors<_F> &VL, RawVectors<_F> &VR);

  template <typename _F>
  void biOrthoNormalize(size_t nR, MBExpansionSet<_F> &VL, MBExpansionSet<_F> &VR);

  template <typename _F>
  std::vector<size_t> getGuessIndices(size_t nGuess, size_t length, const EOMSettings& eomSettings,
                                      const _F *eomDiag, MPI_Comm comm); // getGuessIndices

  //template <typename MatsT>
  //inline double EOMCCBase<MatsT>::signT(size_t a, size_t b, size_t c, size_t i, size_t j, size_t k) const {
  //  if (a == b or a == c or b == c or i == j or i == k or j == k) {
  //    a = 0;
  //    b = 0;
  //    c = 0;
  //    i = 0;
  //    j = 0;
  //    k = 0;
  //    return 0.0;
  //  }
  //  double sign = 1.0;
  //  if (a > b) {
  //    std::swap(a,b);
  //    sign *= -1.0;
  //  }
  //  if (b > c) {
  //    std::swap(b,c);
  //    sign *= -1.0;
  //  }
  //  if (i > j) {
  //    std::swap(i,j);
  //    sign *= -1.0;
  //  }
  //  if (j > k) {
  //    std::swap(j,k);
  //    sign *= -1.0;
  //  }
  //  return sign;
  //}

  template <typename MatsT>
  EOMCCBase<MatsT>::EOMCCBase(const SafeFile &savFile,
                            CCIntermediates<MatsT> &intermediates,
                            const EOMSettings &eomSettings,
                            const CoupledClusterSettings &ccSettings):
     savFile_(savFile),
     intermediates_(intermediates), eomSettings(eomSettings),
     ccSettings_(ccSettings),
     fockMatrix_ta(intermediates.fockMatrix),
     muMatrix(intermediates.muMatrix),
     antiSymMoints(intermediates.antiSymMoInts),
     riMoints(intermediates.riMoInts){
     }

  template <typename MatsT>
  void EOMCCBase<MatsT>::prepEOMCC() {
    initializeEOMCC();
    formEOMIntermediates();
  }

  template <typename MatsT>
  void EOMCCBase<MatsT>::full_diagonalization() {

    theta = CQMemManager::get().malloc<dcomplex>(Hbar_dim);
    RawVectors<MatsT> VL(MPI_COMM_WORLD, Hbar_dim, Hbar_dim);
    RawVectors<MatsT> VR(MPI_COMM_WORLD, Hbar_dim, Hbar_dim);

    std::cout << " Start building the full matrix " << std::endl;

    auto beginBuild = tick();
    cqmatrix::Matrix<MatsT> fullMat = buildHbar(false);
    std::cout << "buildHbar spend " << tock(beginBuild) << " s." << std::endl;

    //fullMat.output(std::cout, "Hbar", true);

//#define TEST_AGAINST_SIGMA
#ifdef TEST_AGAINST_SIGMA
    MatsT* diag = CQMemManager::get().malloc<MatsT>(Hbar_dim);
    beginBuild = tick();
    buildDiag(diag, intermediates_.eps);
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

    //prettyPrintSmart(std::cout, "VL", VL.getPtr(), Hbar_dim, Hbar_dim, Hbar_dim);
    //prettyPrintSmart(std::cout, "VR", VR.getPtr(), Hbar_dim, Hbar_dim, Hbar_dim);

    MatsT* VLVR = CQMemManager::get().malloc<MatsT>(Hbar_dim * Hbar_dim);
    VL.dot_product(0, VR, 0, Hbar_dim, Hbar_dim, VLVR, Hbar_dim);
    //prettyPrintSmart(std::cout, "VLVR", VLVR, Hbar_dim, Hbar_dim, Hbar_dim);
    CQMemManager::get().free(VLVR);
      SetLAThreads(1);
    } // ROOT_ONLY section

    if (not eomSettings.oscillator_strength) {
      std::cout << " Eigenvalues from the full matrix: " << std::endl;
      std::cout << std::fixed << std::setprecision(12);
      for ( auto i = 0; i < Hbar_dim; i++)
        std::cout << i << " EigV: " << theta[i] <<  std::endl;

      return;
    }




    // Building Hbar matrix including H0S H0D blocks
    size_t Hbar_dim_w0 = Hbar_dim + 1;
    dcomplex* theta_w0 = nullptr;
    RawVectors<MatsT> VL_w0(MPI_COMM_WORLD, Hbar_dim_w0, Hbar_dim_w0);
    RawVectors<MatsT> VR_w0(MPI_COMM_WORLD, Hbar_dim_w0, Hbar_dim_w0);

    theta_w0 = CQMemManager::get().malloc<dcomplex>(Hbar_dim_w0);

    std::cout << " Start building the full matrix " << std::endl;

    auto beginBuild_w0 = tick();
    cqmatrix::Matrix<MatsT> fullMat_w0(buildHbar(true));
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

    //std::cout << " Eigenvalues from the full matrix: " << std::endl;
    //std::cout << std::fixed << std::setprecision(12);
    //for ( auto i = 0; i < Hbar_dim_w0; i++)
    //  std::cout << i << " EigV: " << theta_w0[i] <<  std::endl;
    //prettyPrintSmart(std::cout, "VL_w0", VL_w0.getPtr(), Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);
    //prettyPrintSmart(std::cout, "VR_w0", VR_w0.getPtr(), Hbar_dim_w0, Hbar_dim_w0, Hbar_dim_w0);
      SetLAThreads(1);
    }// END ROOT_ONLY section

#ifdef EOMCC_CHECK_LAMBDA
    if (not eomSettings.containActive()) {

      // Begin Lambda check
      // Convert L1 and L2 amplitudes to a vector
      MatsT* L0p1 = CQMemManager::get().malloc<MatsT>(this->getHbarDim(true));
      Lg_->toRaw(L0p1, true);

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
        H0VR[i] /= std::real(theta[i]);
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
      MatsT theta_i = 0.0;
      if constexpr (std::is_same_v<MatsT, dcomplex>) {
        theta_i = theta_w0[i];
      } else {
        theta_i = std::real(theta_w0[i]);
      }
      SetMat('N', Hbar_dim_w0, 1, theta_i, VR_w0.getPtr(i), Hbar_dim_w0, SCR + i * Hbar_dim_w0, Hbar_dim_w0);
    }
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
               Hbar_dim_w0,Hbar_dim_w0,Hbar_dim_w0,
               MatsT(1.),fullMat_w0_copy.pointer(),Hbar_dim_w0,VR_w0.getPtr(),Hbar_dim_w0,MatsT(-1.),SCR,Hbar_dim_w0);
    std::cout << "Right eigenvector check error norm = " << std::scientific << std::setprecision(4)
              << blas::nrm2(Hbar_dim_w0*Hbar_dim_w0, SCR, 1) << std::endl;

    for (size_t i = 0; i < Hbar_dim_w0; i++) {
      MatsT theta_i = 0.0;
      if constexpr (std::is_same_v<MatsT, dcomplex>) {
        theta_i = theta_w0[i];
      } else {
        theta_i = std::real(theta_w0[i]);
      }
      SetMat('C', Hbar_dim_w0, 1, theta_i, VL_w0.getPtr(i), Hbar_dim_w0, SCR + i, Hbar_dim_w0);
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
    TA::get_default_world().gop.broadcast(rawPtr, this->getHbarDim(true), 0);

    Lg_->fromRaw(rawPtr, true);
    if (MPIRank() != 0)
      CQMemManager::get().free(rawPtr);

    //VL_w0.print(std::cout, "right before from Raw");
    //VR_w0.print(std::cout, "right before from Raw");

    L_ = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, eomSettings.nroots, savFile_);
    std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(L_)->fromRaw(MPI_COMM_WORLD, VL_w0, *this, true, 0, 1, eomSettings.nroots);
    R_ = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, eomSettings.nroots, savFile_);
    std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(R_)->fromRaw(MPI_COMM_WORLD, VR_w0, *this, true, 0, 1, eomSettings.nroots);
    initializeDensity();
    std::vector<dcomplex> oscStrength;
    std::vector<dcomplex> excitationE;
    std::cout << "----------------------------------------------" << std::endl;
    for (size_t i = 0; i < eomSettings.nroots ; i++){
      dcomplex f = calcOscillatorStrength(i);
      std::cout << " Excited State: " << i << "  E = " << std::setprecision(12) << theta[i] << " Eh, "<< " f = " << f << std::endl;
      oscStrength.push_back(f);
      excitationE.push_back(theta[i]);
    }

    TA::get_default_world().gop.fence();
    // Write data to bin file
    if (savFile_.exists()) {
      savFile_.safeWriteData("/CC/EXCITATION_ENERGIES",excitationE.data(), {eomSettings.nroots});
      savFile_.safeWriteData("/CC/OSCILLATOR_STRENGTHS", oscStrength.data(), {eomSettings.nroots});
    }

    CQMemManager::get().free(theta_w0, VLVR_w0, VRVR_w0, SCR);
    if (theta) CQMemManager::get().free(theta);

  }

  template <typename MatsT>
  EOMCCBase<MatsT>::~EOMCCBase() {


    //TAManager &TAmanager = TAManager::get();

    //base if (Rho_ij) TAmanager.free("oo", std::move(Rho_ij), true);
    //base if (Rho_ab) TAmanager.free("vv", std::move(Rho_ab), true);
    //base if (Rho_ia) TAmanager.free("ov", std::move(Rho_ia), true);
    //base if (Rho_ai) TAmanager.free("vo", std::move(Rho_ai), true);

  }

  template <typename MatsT>
  void EOMCCBase<MatsT>:: davidsonSolve() { 

      std::shared_ptr<Davidson<MatsT>> davidson_r = nullptr;
      std::shared_ptr<Davidson<MatsT>> davidson_l = nullptr;
    // solving for R is required
      size_t Hbar_dim = getHbarDim();
      MatsT* eomDiag = CQMemManager::get().malloc<MatsT>(Hbar_dim);

      /// Functions for Davidson
      EOMCCEigenVecType eigenVecType = EOMCCEigenVecType::RIGHT;

      size_t nGuess = eomSettings.nroots * eomSettings.davidson_guess_multiplier;
      nGuess = std::min(nGuess, getHbarDim());
      dcomplex * curEig = CQMemManager::get().malloc<dcomplex>(nGuess);

      typename Davidson<MatsT>::VecsGen_t vecsGenerator =  this->EmptyDavidsonVectorBuilder(); //Davidson<MatsT>::VecsGen_t(); // Generator for new vector sets
      typename Davidson<MatsT>::LinearTrans_t sigmaBuilder = this->DavidsonResidualBuilder(eigenVecType); // Sigma vector builder
      typename Davidson<MatsT>::LinearTrans_t preConditioner = this->DavidsonPreconditionerBuilder(curEig, eomDiag); // Preconditioner

      // Clear cached TA objects
      TAManager::get().discard_cache();

      std::cout << std::endl << "Davidson-Liu algorithm for EOMCC:" << std::endl;
      std::cout << BannerMid << std::endl << std::endl;

      // Compute diagonal elements
      auto beginBuildDiag = tick();
      buildDiag(eomDiag, intermediates_.eps);
      std::cout << "  * Build diagonal elements for iterative solver spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginBuildDiag) << " s." << std::endl << std::endl;

    if (not eomSettings.skip_r) {
      std::cout << "Right eigensolver iterations:" << std::endl << std::endl;
      auto beginRightEig = tick();

      davidson_r = std::make_shared<Davidson<MatsT>> (MPI_COMM_WORLD, getHbarDim(),
                                  eomSettings.davidson_max_macro_iter,
                                  eomSettings.davidson_max_micro_iter,
                                  eomSettings.davidson_residual_conv, eomSettings.nroots,
                                  sigmaBuilder, preConditioner, vecsGenerator);

      davidson_r->setWhenSc(eomSettings.davidson_whenSc);
      davidson_r->setM(eomSettings.davidson_subspace_multiplier);
      davidson_r->setkG(eomSettings.davidson_guess_multiplier);
      davidson_r->setGramSchmidtRepeat(eomSettings.GramSchmidt_NRe);
      davidson_r->setGramSchmidtEps(eomSettings.GramSchmidt_eps);
      davidson_r->setResidueConvCheck(eomSettings.davidson_check_residual);
      davidson_r->setEigenVectorConvCheck(eomSettings.davidson_check_eigen_vector);
      davidson_r->setEigenValueConvCheck(eomSettings.davidson_check_eigen_value);
      davidson_r->setConvOnGramSchmidt(eomSettings.davidson_conv_on_GramSchmidt);
      davidson_r->setEigenVectorConvCriteria(eomSettings.davidson_eigen_vector_conv);
      davidson_r->setEigenValueConvCriteria(eomSettings.davidson_eigen_value_conv);
      davidson_r->setEigForT(curEig);
      davidson_r->setSaveOption(eomSettings.save_r);
      davidson_r->setSavePrefix("/CC/EOMRVECTOR");

      if (!eomSettings.davidson_Eref.empty())
        davidson_r->setEnergySpecific(eomSettings.davidson_Eref, eomSettings.davidson_ErefAbs);

      //if (eomSettings.davidson_Eref > 0.0) {
      //  davidson.useEnergySpecific(eomSettings.nroots, 0.0, eomSettings.davidson_Eref);
      //  if (eomSettings.davidson_sort_by_distance)
      //    davidson.setSortByDistance();
      //}

      //bool sortByDistance = eomSettings.davidson_sort_by_distance;
      std::vector<size_t> guessIndices = getGuessIndices(nGuess, getHbarDim(), eomSettings, eomDiag, MPI_COMM_WORLD);//, sortByDistance);
      // set guess for Davidson solver
      if (eomSettings.restart_r) {
        size_t nroots = eomSettings.nroots;
        MatsT * r_vector = CQMemManager::get().malloc<MatsT>(Hbar_dim);
        davidson_r->setGuess(nGuess, [&guessIndices, &r_vector, Hbar_dim, nroots, this] (size_t nGuess, SolverVectors<MatsT> &guessVec, size_t length) {
          guessVec.clear();
          if (auto derivedObj = dynamic_cast<MBExpansionSet<MatsT>*>(&guessVec)){
            int i = 0;
            for (i = 0; i< nroots; i++) {
              bool vecIexist = false;
              if (MPIRank() == 0) {
                vecIexist = savFile_.exists("/CC/EOMRVECTOR"+std::to_string(i));
                if (vecIexist) {
                  std::cout << "EOMCC Restart: right eigenvector of root " << i << " found in the restart file." << std::endl;
                  savFile_.readData("/CC/EOMRVECTOR"+std::to_string(i), r_vector);
                }
              }
              TA::get_default_world().gop.fence();
              MPIBCast(vecIexist, 0, MPI_COMM_WORLD);
              if (vecIexist) {
                MPIBCast(r_vector, Hbar_dim, 0, MPI_COMM_WORLD);
                TA::get_default_world().gop.fence();
                derivedObj->get(i).fromRaw(r_vector, false);
              } else {
                std::cout << "EOMCC Restart: right eigenvector of root " << i << " NOT found in the restart file, traditional guess appended thereafter." << std::endl;
                break; // stop if any of the root vector is missing
              }
            }
            for (; i < nGuess; i++) {
              derivedObj->set(guessIndices[i], i, 1.0);
            }
          }
        });
        MatsT CorrE;
        // if (MPIRank() == 0 and not this->ccSettings_.skipSCF) savFile_.readData("/CC/REFERENCE_ENERGY",&intermediates_.E_ref);
        if (MPIRank() == 0) savFile_.readData("/CC/CORRELATION_ENERGY",&CorrE);             
        intermediates_.E_cc = intermediates_.E_ref + std::real(CorrE);
        MPIBCast(&intermediates_.E_cc, 1, 0, MPI_COMM_WORLD);
        MPIBCast(&intermediates_.E_ref, 1, 0, MPI_COMM_WORLD);
        TA::get_default_world().gop.fence();
        CQMemManager::get().free(r_vector);
      }
      else { // not restart
        if (eomSettings.ccs_guess) { // better guess
          if (eomSettings.eom_implementation != EOM_IMPLEMENTATION::EOMCCSD and
              eomSettings.eom_implementation != EOM_IMPLEMENTATION::CVSEOMCCSD and
              eomSettings.eom_implementation != EOM_IMPLEMENTATION::EOMDIP_3h1p) {
            CErr("EOMCC CCS Guess is only implemented for EOMCCSD, CVSEOMCCSD, and EOMDIP_3h1p for now!");
          }
          size_t guess_dim = oneBodySize();
          nGuess = std::min(nGuess, guess_dim);
          MatsT * guess_vector = CQMemManager::get().malloc<MatsT>(guess_dim * nGuess);
          fillGuess(guess_vector, nGuess);
          davidson_r->setGuess(nGuess, [&guess_vector, guess_dim, this] (size_t nGuess, SolverVectors<MatsT> &guessVec, size_t length) {
            guessVec.clear();
            if (auto derivedObj = dynamic_cast<MBExpansionSet<MatsT>*>(&guessVec)){
              for (int i = 0; i< nGuess; i++) {
                derivedObj->get(i).partlyFromRaw(guess_vector + i * guess_dim, 1, false);
                }
              }
          });
          CQMemManager::get().free(guess_vector);

        }
        else { // traditional guess
          davidson_r->setGuess(nGuess, [&guessIndices] (size_t nGuess, SolverVectors<MatsT> &guessVec, size_t length) {
            guessVec.clear();
            for(size_t i = 0; i < nGuess; i++) guessVec.set(guessIndices[i], i, 1.0);
          });
        }
      }


      davidson_r->run();

      if (not davidson_r->hasConverged()) {
        CErr("EOMCC right Davidson iteration failed to converge!");
      }

      setR(davidson_r->VR());

        std::cout << "EOMCC energy results:" << std::endl << std::endl;

        std::cout << std::setw(18) << std::left <<  "  Excited states";
        std::cout << std::setw(34) << std::left << "Excitation Energy (Eh)";
        std::cout << std::setw(19) << std::left << "Total Energy (Eh)";
        std::cout << std::endl;
        std::cout << std::setw(18) << std::left <<  "  -------------";
        std::cout << std::setw(34) << std::left << "-----------------";
        std::cout << std::setw(18) << std::left << "-----------------";
        std::cout << std::endl;

        for (size_t i = 0; i < eomSettings.nroots ; i++){
          MatsT E_tot = intermediates_.E_cc;
          if constexpr (std::is_same_v<MatsT, double>) {
            E_tot += std::real(davidson_r->eigVal()[i]);
          } else {
            E_tot += davidson_r->eigVal()[i];
          }
          std::cout << std::setprecision(12) << std::fixed;
          std::cout << "      State "  << std::setw(6) << std::left << i+1;
          std::cout << std::setw(34) << std::left;
          if (std::abs(davidson_r->eigVal()[i]) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << davidson_r->eigVal()[i];
          std::cout << std::setw(34) << std::left;
          if (std::abs(E_tot) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << E_tot;
          std::cout << std::endl;
          
        }
      std::cout << "  * EOMCC right eigensolver spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginRightEig) << " s." << std::endl;
      if (MPIRank() == 0) savFile_.safeWriteData("/CC/EXCITATION_ENERGIES", davidson_r->eigVal(), {eomSettings.nroots});
      if (eomSettings.save_r) {
        MatsT * r_vector;
        r_vector = CQMemManager::get().malloc<MatsT>(Hbar_dim );
        for(size_t i = 0; i < eomSettings.nroots; i++) {
          if (auto VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(davidson_r->VR()))
            VR->get(i).toRaw(r_vector, false);
          else
            CErr("Invalid type for R_");
          TA::get_default_world().gop.fence();
          if (eomSettings.print_large_amplitude)
            print_largest_values_and_position(i, r_vector, Hbar_dim, eomSettings, intermediates_);
          if (MPIRank() == 0)
              savFile_.safeWriteData("/CC/EOMRVECTOR"+std::to_string(i), r_vector, {Hbar_dim});
        }
        TA::get_default_world().gop.fence();
        CQMemManager::get().free(r_vector);
      }

    }
    // solving for L is required
    if (eomSettings.oscillator_strength) {
        // Initialize density matrices for oscillator strength calculation
        initializeDensity();

        // Run CC lambda equations
        initializeGroundStateLambda();
        runLambda(); // has CVS version

        if (eomSettings.all_excited_dipole) {
          // transition dipole moments between excited states requested
          std::array<MatsT, 3> mu = calcGroundDipole();

          std::cout << "  Ground State Dipole Moment w/o Nuclear Contribution no Negative Sign (a.u.): " << std::endl
                    << "    X=" << std::setprecision(12) << std::fixed << mu[0]
                    << "    Y=" << std::setprecision(12) << std::fixed << mu[1]
                    << "    Z=" << std::setprecision(12) << std::fixed << mu[2]
                    << std::endl << std::endl;

          // Write data to bin file
          if (savFile_.exists())
            savFile_.safeWriteData("/CC/GROUND_STATE_DIPOLE", mu.data(), {3});
        }

        eigenVecType = EOMCCEigenVecType::LEFT;

        std::cout << BannerMid << std::endl << std::endl;
        std::cout << "Left eigensolver iterations:" << std::endl << std::endl;
        auto beginLeftEig = tick();

        davidson_l = std::make_shared<Davidson<MatsT>> (MPI_COMM_WORLD, getHbarDim(),
                                    eomSettings.davidson_max_macro_iter,
                                    eomSettings.davidson_max_micro_iter,
                                    eomSettings.davidson_residual_conv,
                                    eomSettings.nroots,
                                    sigmaBuilder, preConditioner, vecsGenerator);

        davidson_l->setWhenSc(eomSettings.davidson_whenSc);
        davidson_l->setM(eomSettings.davidson_subspace_multiplier);
        davidson_l->setkG(1);
        davidson_l->setGramSchmidtRepeat(eomSettings.GramSchmidt_NRe);
        davidson_l->setGramSchmidtEps(eomSettings.GramSchmidt_eps);
        davidson_l->setResidueConvCheck(eomSettings.davidson_check_residual);
        davidson_l->setEigenVectorConvCheck(eomSettings.davidson_check_eigen_vector);
        davidson_l->setEigenValueConvCheck(eomSettings.davidson_check_eigen_value);
        davidson_l->setConvOnGramSchmidt(eomSettings.davidson_conv_on_GramSchmidt);
        davidson_l->setEigenVectorConvCriteria(eomSettings.davidson_eigen_vector_conv);
        davidson_l->setEigenValueConvCriteria(eomSettings.davidson_eigen_value_conv);
        davidson_l->setEigForT(curEig);
        davidson_l->setSaveOption(eomSettings.save_l);
        davidson_l->setSavePrefix("/CC/EOMLVECTOR");

        // Reuse scratch spaces for right Davidson
        if (davidson_r) {
          davidson_l->setGuessScratch(davidson_r->getGuessScratch());
          davidson_l->setSubspaceScratch(davidson_r->getSubspaceScratch());
          davidson_l->setSigmaVecScratch(davidson_r->getSigmaVecScratch());
          davidson_l->setScratchR(davidson_r->getScratchR());
          davidson_l->setScratchS(davidson_r->getScratchS());

          davidson_r->clear_scratch();
        }

        if (!eomSettings.davidson_Eref.empty())
          davidson_l->setEnergySpecific(eomSettings.davidson_Eref, eomSettings.davidson_ErefAbs);

        //if (eomSettings.davidson_Eref > 0.0) {
        //  davidsonLeft.useEnergySpecific(eomSettings.nroots, 0.0, eomSettings.davidson_Eref);
        //  if (eomSettings.davidson_sort_by_distance)
        //    davidson_l->setSortByDistance();

        // guess for L
        if (eomSettings.skip_r or eomSettings.restart_l) {
          size_t nroots = eomSettings.nroots;
          MatsT * l_vector = CQMemManager::get().malloc<MatsT>(Hbar_dim);
          davidson_l->setGuess(nGuess, [&l_vector, Hbar_dim, nroots, eomDiag, this] (size_t nGuess, SolverVectors<MatsT> &guessVec, size_t length) {
            guessVec.clear();
            if (auto derivedObj = dynamic_cast<MBExpansionSet<MatsT>*>(&guessVec)){
              int i = 0;
              for (i = 0; i< nroots; i++) {
                bool vecIexist = false;
                if (MPIRank() == 0) {
                  vecIexist = savFile_.exists("/CC/EOMLVECTOR"+std::to_string(i));
                  if (vecIexist) {
                    std::cout << "EOMCC Restart: left eigenvector of root " << i << " found in the restart file." << std::endl;
                    savFile_.readData("/CC/EOMLVECTOR"+std::to_string(i), l_vector);
                  } else {
                    vecIexist = savFile_.exists("/CC/EOMRVECTOR"+std::to_string(i));
                    if (vecIexist) {
                      std::cout << "EOMCC Restart: left eigenvector of root " << i << " NOT found in the restart file. Reading the right eigenvector as guess..." << std::endl;
                      savFile_.readData("/CC/EOMRVECTOR"+std::to_string(i), l_vector);
                    }
                  }
                }
                TA::get_default_world().gop.fence();
                MPIBCast(vecIexist, 0, MPI_COMM_WORLD);
                if (vecIexist) {
                  MPIBCast(l_vector, Hbar_dim, 0, MPI_COMM_WORLD);
                  TA::get_default_world().gop.fence();
                  derivedObj->get(i).fromRaw(l_vector, false);
                } else {
                  std::cout << "EOMCC Restart: neither left nor right eigenvector of root " << i << " NOT found in the restart file, traditional guess appended thereafter." << std::endl;
                  break; // stop if any of the root vector is missing
                }
              }
              if (i == nGuess) return; // all guess vectors are set, return directly

              std::vector<size_t> guessIndices = getGuessIndices(nGuess, getHbarDim(), eomSettings, eomDiag, MPI_COMM_WORLD);//, sortByDistance);
              for (; i < nGuess; i++) {
                derivedObj->set(guessIndices[i], i, 1.0);
              }
            }
          });
          MatsT CorrE;
          // if (MPIRank() == 0 and not this->ccSettings_.skipSCF) savFile_.readData("/CC/REFERENCE_ENERGY",&intermediates_.E_ref);
          if (MPIRank() == 0) savFile_.readData("/CC/CORRELATION_ENERGY",&CorrE);             
          intermediates_.E_cc = intermediates_.E_ref + std::real(CorrE);
          TA::get_default_world().gop.fence();
          CQMemManager::get().free(l_vector);
        } else {
          davidson_l->setGuess(eomSettings.nroots, [davidson_r] (size_t nGuess, SolverVectors<MatsT> &guessVec, size_t length) {
            guessVec.clear();
            guessVec.set_data(0, nGuess, *davidson_r->VR(), 0);
            guessVec.conjugate(0, nGuess);
          });
        }

        davidson_l->run();
        if (not davidson_l->hasConverged()) {
          CErr("EOMCC left Davidson iteration failed to converge!");
        }
        davidson_l->clear_scratch();

        std::cout << "  * EOMCC left eigensolver spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginLeftEig) << " s." << std::endl;
        std::cout << BannerMid << std::endl;

        std::shared_ptr<MBExpansionSet<MatsT>> VL, VR;

        std::shared_ptr<MBExpansionSet<MatsT>> r_saved = nullptr;
        std::vector<dcomplex> excitation_energies(eomSettings.nroots);
        if (eomSettings.skip_r) {
          // read R from scratch file
          davidson_r = nullptr;
          size_t nroots = eomSettings.nroots;

          if (dynamic_cast<EOMRCCSD<MatsT>*>(this) != nullptr) // Restricted CC calculation, symmetry RCCSD
            r_saved = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, nroots, savFile_, MBTensorSymmetry::RCCSD);
          else
            r_saved = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, nroots, savFile_);
          r_saved->clear();
          MatsT * r_vector = CQMemManager::get().malloc<MatsT>(Hbar_dim);
          if (auto derivedObj = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(r_saved)){
            for (int i = 0; i< nroots; i++) {
              std::cout << "EOMCC Skip R: reading right eigenvector of root " << i << " from the restart file." << std::endl;
              if (MPIRank() == 0) savFile_.readData("/CC/EOMRVECTOR"+std::to_string(i), r_vector);
              TA::get_default_world().gop.fence();
              MPIBCast(r_vector, Hbar_dim, 0, MPI_COMM_WORLD);
              TA::get_default_world().gop.fence();
              r_saved->get(i).fromRaw(r_vector, false);
            }
          }
          if (MPIRank() == 0) savFile_.readData("/CC/EXCITATION_ENERGIES", excitation_energies.data());
          TA::get_default_world().gop.fence();
          MPIBCast(excitation_energies.data(), nroots, 0, MPI_COMM_WORLD);
          CQMemManager::get().free(r_vector);
        }
        switch (eomSettings.hbar_type) {
          case EOM_HBAR_TYPE::IMPLICIT:
            if (davidson_r)
                VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(davidson_r->VR());
            else
                VR = r_saved;
            VL = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(davidson_l->VR());
            break;
          case EOM_HBAR_TYPE::DEBUG:
            VR = std::make_shared<MBExpansionSet<MatsT>>(
                std::dynamic_pointer_cast<MBExpansionSetDebug<MatsT>>(davidson_r->VR())->getEOMCCSet());
            VL = std::make_shared<MBExpansionSet<MatsT>>(
                std::dynamic_pointer_cast<MBExpansionSetDebug<MatsT>>(davidson_l->VR())->getEOMCCSet());
            break;
          case EOM_HBAR_TYPE::EXPLICIT:
            VR = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, eomSettings.nroots, savFile_);
            VR->fromRaw(MPI_COMM_WORLD,
                        *std::dynamic_pointer_cast<RawVectors<MatsT>>(davidson_r->VR()),
                        *this, false, 0, 0, eomSettings.nroots);
            VL = std::make_shared<MBExpansionSet<MatsT>>(tensor_builder_, eomSettings.nroots, savFile_);
            VL->fromRaw(MPI_COMM_WORLD,
                        *std::dynamic_pointer_cast<RawVectors<MatsT>>(davidson_l->VR()),
                        *this, false, 0, 0, eomSettings.nroots);
            break;

        }

        setR(VR);
        setL(VL);
        if (davidson_r) {
          excitation_energies = std::vector<dcomplex>(davidson_r->eigVal(), davidson_r->eigVal() + eomSettings.nroots);
        }
        setTheta(excitation_energies.data(), eomSettings.nroots);

        buildRightZeroBody(eomSettings.nroots);

        if (eomSettings.davidson_biortho) {
          // BiOrthonormalization
          std::cout << "  *** BiOrthonormalize left and right eigenvectors ***" << std::endl;

          biOrthoNormalize(eomSettings.nroots, *VL, *VR);

        }

        std::cout << BannerMid << std::endl;

        std::cout << "EOMCC results:" << std::endl << std::endl;

        std::cout << std::setw(18) << std::left <<  "  Excited states";
        std::cout << std::setw(34) << std::left << "Excitation Energy (Eh)";
        std::cout << std::setw(19) << std::left << "Oscillator Strength";
        std::cout << std::endl;
        std::cout << std::setw(18) << std::left <<  "  -------------";
        std::cout << std::setw(34) << std::left << "-----------------";
        std::cout << std::setw(18) << std::left << "-----------------";
        std::cout << std::endl;

        auto beginOsc = tick();

        std::cout << "----------------------------------------------" << std::endl;
        std::vector<dcomplex> oscStrength;
        std::vector<dcomplex> excitationE;

        if (eomSettings.all_excited_dipole) {
          // transition dipole moments between excited states requested
          bool isRestricted = false;
          EOMCCSD<MatsT>* eomccsd_ptr = dynamic_cast<EOMCCSD<MatsT>*>(this);
          EOMRCCSD<MatsT>* eomrccsd_ptr = dynamic_cast<EOMRCCSD<MatsT>*>(this);
          if (eomccsd_ptr == nullptr) {
            if (eomrccsd_ptr != nullptr)
              isRestricted = true;
            else
              CErr("Transition dipole moments between excited states is only implemented for EOMCCSD now!");
          }

          cqmatrix::Matrix<MatsT> mu_g2x(eomSettings.nroots, 3), mu_x2g(eomSettings.nroots, 3);
          mu_g2x.clear();
          mu_x2g.clear();
          for (size_t i = 0; i < eomSettings.nroots ; i++){
            std::array<MatsT, 3> mu_g2x_i, mu_x2g_i;

            if (isRestricted) {
              mu_g2x_i = eomrccsd_ptr->calcGround2ExcitedTransitionDipole(i);
              mu_x2g_i = eomrccsd_ptr->calcExcited2GroundTransitionDipole(i);
            } else {
              mu_g2x_i = eomccsd_ptr->calcGround2ExcitedTransitionDipole(i);
              mu_x2g_i = eomccsd_ptr->calcExcited2GroundTransitionDipole(i);
            }

            for (size_t j = 0; j < 3; j++) {
              mu_g2x(i, j) = mu_g2x_i[j];
              mu_x2g(i, j) = mu_x2g_i[j];
            }

            MatsT DS = mu_g2x_i[0] * mu_x2g_i[0] + mu_g2x_i[1] * mu_x2g_i[1] + mu_g2x_i[2] * mu_x2g_i[2];

            dcomplex f = 2./3 * this->theta[i] * DS;

            std::cout << std::setprecision(12) << std::fixed;
            std::cout << "      State "  << std::setw(6) << std::left << i+1;
            std::cout << std::setw(34) << std::left;
            if (std::abs(davidson_l->eigVal()[i]) > 1e-6)
              std::cout << std::fixed << std::setprecision(12);
            else
              std::cout << std::scientific << std::setprecision(6);
            std::cout << davidson_l->eigVal()[i];
            std::cout << std::setw(34) << std::left;
            if (std::abs(f) > 1e-6)
              std::cout << std::fixed << std::setprecision(12);
            else
              std::cout << std::scientific << std::setprecision(6);
            std::cout << f;
            std::cout << std::endl;

            oscStrength.push_back(f);
            excitationE.push_back(excitation_energies[i]);
          }

          prettyPrintSmart(std::cout, "Ground to excited state transition electric dipole moments (a.u.):",
                           mu_g2x.pointer(), eomSettings.nroots, 3, eomSettings.nroots);
          prettyPrintSmart(std::cout, "Excited to ground state transition electric dipole moments (a.u.):",
                           mu_x2g.pointer(), eomSettings.nroots, 3, eomSettings.nroots);

          // Write data to bin file
          if (savFile_.exists()) {
            savFile_.safeWriteData("/CC/GROUND_TO_EXCITED_TRANSITION_DIPOLE", mu_g2x.pointer(), {3, eomSettings.nroots});
            savFile_.safeWriteData("/CC/EXCITED_TO_GROUND_TRANSITION_DIPOLE", mu_x2g.pointer(), {3, eomSettings.nroots});
          }

          // Implementation of transition dipole moments between excited states for EOMCCSD
          cqmatrix::NDArray<MatsT> mu_x2x({eomSettings.nroots, eomSettings.nroots, 3});
          mu_x2x.clear();
          for (size_t i = 0; i < eomSettings.nroots ; i++){
            for (size_t j = 0; j < eomSettings.nroots ; j++){
              std::array<MatsT, 3> mu_x2x_ij;
              if (isRestricted)
                mu_x2x_ij = eomrccsd_ptr->calcExcited2ExcitedTransitionDipole(i, j);
              else
                mu_x2x_ij = eomccsd_ptr->calcExcited2ExcitedTransitionDipole(i, j);
              for (size_t k = 0; k < 3; k++) {
                mu_x2x(i, j, k) = mu_x2x_ij[k];
              }
            }
          }

          prettyPrintSmart(std::cout, "Transition Dipole Moments between Excited States X Component (a.u.)",
                           &mu_x2x(0,0,0), eomSettings.nroots, eomSettings.nroots, eomSettings.nroots);
          prettyPrintSmart(std::cout, "Transition Dipole Moments between Excited States Y Component (a.u.)",
                           &mu_x2x(0,0,1), eomSettings.nroots, eomSettings.nroots, eomSettings.nroots);
          prettyPrintSmart(std::cout, "Transition Dipole Moments between Excited States Z Component (a.u.)",
                           &mu_x2x(0,0,2), eomSettings.nroots, eomSettings.nroots, eomSettings.nroots);

          // Write data to bin file
          if (savFile_.exists())
            savFile_.safeWriteData("/CC/EXCITED_TO_EXCITED_TRANSITION_DIPOLE", mu_x2x.pointer(), {3, eomSettings.nroots, eomSettings.nroots});

        } else // normal oscillator strength calculation
        for (size_t i = 0; i < eomSettings.nroots ; i++){
          dcomplex f = calcOscillatorStrength(i);

          std::cout << std::setprecision(12) << std::fixed;
          std::cout << "      State "  << std::setw(6) << std::left << i+1;
          std::cout << std::setw(34) << std::left;
          if (std::abs(davidson_l->eigVal()[i]) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << davidson_l->eigVal()[i];
          std::cout << std::setw(34) << std::left;
          if (std::abs(f) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << f;
          std::cout << std::endl;
          
          oscStrength.push_back(f);
          excitationE.push_back(excitation_energies[i]);
        }

        TA::get_default_world().gop.fence();
        // Write data to bin file
        if (savFile_.exists()) {
          savFile_.safeWriteData("/CC/EXCITATION_ENERGIES",excitation_energies.data(), {eomSettings.nroots});
          savFile_.safeWriteData("/CC/OSCILLATOR_STRENGTHS", oscStrength.data(), {eomSettings.nroots});
        }

        std::cout << std::endl << "  * Compute oscillator strength spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginOsc) << " s." << std::endl;

      }

      if (eomSettings.print_large_amplitude) {
        MatsT * r_vector;
        r_vector = CQMemManager::get().malloc<MatsT>(Hbar_dim );
        for(size_t i = 0; i < eomSettings.nroots; i++) {
          if (auto VR = std::dynamic_pointer_cast<MBExpansionSet<MatsT>>(R_))
            VR->get(i).toRaw(r_vector, false);
          else 
            CErr("Invalid type for R_");
          TA::get_default_world().gop.fence();
          if (eomSettings.print_large_amplitude)
            print_largest_values_and_position(i, r_vector, Hbar_dim, eomSettings, intermediates_);
        }
        TA::get_default_world().gop.fence();
        CQMemManager::get().free(r_vector);
      }


      if(davidson_r) davidson_r->clear_scratch();

      CQMemManager::get().free(eomDiag);
      if (curEig) CQMemManager::get().free(curEig);

    }


/// for full_diagonalization routine
  template <typename MatsT>
  void EOMCCBase<MatsT>::buildHbar_sigma(MatsT * out, bool diagOnly) const {

    TA::get_default_world().gop.fence();

    std::vector<MBExpansion<MatsT>> rs;
    rs.reserve(Hbar_dim);
    for (size_t i = 0; i < Hbar_dim; i++) {
      rs.emplace_back(tensor_builder_, true);
    }
    
    std::vector<MBExpansion<MatsT>> hrs;
    hrs.reserve(Hbar_dim);
    for (size_t i = 0; i < Hbar_dim; i++) {
      hrs.emplace_back(tensor_builder_, true);
    }
    

    for (size_t i = 0; i < Hbar_dim; i++) {
      MBExpansion<MatsT>& Vi = rs[i];
      MBExpansion<MatsT>& AVi = hrs[i];
      Vi.setElem(i, 1.0);
      buildSigma(Vi, AVi, EOMCCEigenVecType::RIGHT);
      AVi.enforceSymmetry();
    }

    for (size_t i = 0; i < Hbar_dim; i++)
      if (diagOnly)
        out[i] = rs[i].dot(hrs[i]);
      else
        for (size_t j = 0; j < Hbar_dim; j++) {
          out[i + j * Hbar_dim] = rs[i].dot(hrs[j]);
        }

  }

  template <typename MatsT>
  void EOMCCBase<MatsT>::print_largest_values_and_position(size_t i, MatsT* r_vector, size_t Hbar_dim, const EOMSettings& eomSettings, CCIntermediates<MatsT> & intermediates) {

    if (eomSettings.eom_type != EOM_TYPE::EE && not eomSettings.containActive()) {
        CErr("PrintLargeAmplitude only implemented with CVS-EOMCCSD");
    }
 
    // size evaluation
    size_t nO = TAManager::get().getRange(intermediates.oLabel).extent();
    size_t nV = TAManager::get().getRange(intermediates.vLabel).extent();
    size_t nCVSVActive = nV;
    size_t nCVSOActive = nO;
    size_t nCVSVContinuum =  nCVSVActive;
    size_t nCVSOCore      = eomSettings.cvs_core.size()    == 0? nCVSOActive : eomSettings.cvs_core.size();
    size_t nCVSoneBodySize = nCVSOCore * nCVSVContinuum;

    // Evaluate the three largest values and their position indices
    double maxMagnitude[3] = {0.0};
    size_t maxIndex[3] = {999};

    for (size_t l = 0; l < Hbar_dim; l++) {
        double magnitude = std::abs(r_vector[l]);

        // Check if the magnitude is larger than any of the current maximum magnitudes
        for (int j = 0; j < 3; j++) {
            if (magnitude > maxMagnitude[j]) {
                // Shift the previous maximum values
                for (int k = 2; k > j; k--) {
                    maxMagnitude[k] = maxMagnitude[k - 1];
                    maxIndex[k] = maxIndex[k - 1];
                }
                // Update the new maximum value
                maxMagnitude[j] = magnitude;
                maxIndex[j] = l;
                break;
            }
        }
    }
    // Print the three largest values and their position indices
    for (int l = 0; l < 3; l++) {
      std::cout << "Largest magnitude of state " << i+1 << ": " << std::fixed<<std::setprecision(4) <<maxMagnitude[l];
      if (maxIndex[l] < nCVSoneBodySize) {
        size_t ii, aa;
        aa = maxIndex[l] % nCVSVContinuum;
        ii = maxIndex[l] / nCVSVContinuum;
            std::cout << " i = " << ii << " th in core; a = " << aa << " th in continuum"<<std::endl;
      } else {
            std::cout << " Double excited state!" << std::endl;
      }
    }
    std::cout << std::endl;
  }



}
