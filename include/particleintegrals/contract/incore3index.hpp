/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

#include <integrals.hpp>
#include <util/matout.hpp>
#include <util/timer.hpp>
#include <util/threads.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/distributedritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/gradints/incore.hpp>

// For when ERI3J is stored with leading dimension of NB2
//#define LEADING_NB2

namespace ChronusQ {

  /**
   *  \brief Perform a Coulomb-type (34,12) RI-ERI contraction with
   *  a one-body operator.
   */   
  template <typename MatsT, typename IntsT>
  void InCoreRITPIContraction<MatsT, IntsT>::JContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    InCoreRITPI<IntsT> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NB2 = NB*NB;
    size_t NBRI = eri3j.nRIBasis();

    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    IntsT *AX = reinterpret_cast<IntsT*>(C.AX);

    // Extract the real part of X if X is Hermetian and if the ints are real
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;

    // Allocate scratch if IntsT and MatsT are different
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;

    if( extractRealPartX ) {

      X = CQMemManager::get().malloc<IntsT>(NB2);
      for(auto k = 0ul; k < NB2; k++) X[k] = std::real(C.X[k]);

    }

    if( allocAXScratch ) {

      AX = CQMemManager::get().malloc<IntsT>(NB2);
      std::fill_n(AX,NB2,0.);

    }


    auto Jtemp = CQMemManager::get().malloc<IntsT>(NBRI);
    std::fill_n(Jtemp, NBRI, IntsT(0.));
    // (ij|Q)S^{-1/2} -> ERI3J
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,1,NB2,IntsT(1.),eri3j.pointer(),NBRI,X,NB2,IntsT(0.),Jtemp,NBRI);
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB2,1,NBRI,IntsT(1.),eri3j.pointer(),NBRI,Jtemp,NBRI,IntsT(0.),AX,NB2);
#ifdef LEADING_NB2
    // Alternative code for when ERI3J is stored with leading dimension of NB2
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBRI,1,NB2,IntsT(1.),eri3j.pointer(),NB2,X,NB2,IntsT(0.),Jtemp,NBRI);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB2,1,NBRI,IntsT(1.),eri3j.pointer(),NB2,Jtemp,NBRI,IntsT(0.),AX,NB2);
#endif
    CQMemManager::get().free(Jtemp);

    // if Complex ints + Hermitian, conjugate
//    if( std::is_same<IntsT,dcomplex>::value and C.HER )
//      IMatCopy('R',NB,NB,IntsT(1.),AX,NB,NB);

    // If non-hermetian, transpose
//    if( not C.HER )  {

//      IMatCopy('T',NB,NB,IntsT(1.),AX,NB,NB);

//    }

    // Cleanup temporaries
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {

      std::copy_n(AX,NB2,C.AX);
      CQMemManager::get().free(AX);

    }

  }; // InCoreRITPIContraction::JContract



  template <typename MatsT, typename IntsT>
  void InCoreRITPIContraction<MatsT, IntsT>::KContract_real_impl(
      MPI_Comm, const double *Xr, double *AXr, double *Ktemp) const {

    InCoreRITPI<double> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<double>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NBRI = eri3j.nRIBasis();
    size_t NBNBRI = NB*NBRI;
    size_t NB2NBRI= NB*NBNBRI;
    auto mapCen2BfSt = this->mapCen2BfSt();
    size_t nAtoms = mapCen2BfSt.size();

    std::fill_n(Ktemp, NB2NBRI, double(0.));
    auto XrT = CQMemManager::get().malloc<double>(NB*NB);
    SetMat('T',NB,NB,double(1.),Xr,NB,XrT,NB);

    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      #pragma omp parallel for schedule(dynamic)
      for (size_t i = 0; i < nAtoms; ++i) {
        size_t bfStart = mapCen2BfSt[i];
        size_t bfEnd   = (i+1 < nAtoms) ? mapCen2BfSt[i+1] : NB;
        size_t nBf = bfEnd - bfStart;
        for(auto nu = bfStart; nu < bfEnd; nu++) 
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,nBf,nBf,double(1.),
                     eri3j.pointer()+nu*NBNBRI+bfStart*NBRI,NBRI,
                     XrT+bfStart*NB+bfStart,NB,double(0.),Ktemp+nu*NBNBRI+bfStart*NBRI,NBRI);
      }
    } else {
      #pragma omp parallel for
      for(auto nu = 0ul; nu < NB; nu++)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,NB,NB,double(1.),eri3j.pointer()+nu*NBNBRI,NBRI,XrT,NB,double(0.),Ktemp+nu*NBNBRI,NBRI);
    }
    CQMemManager::get().free(XrT);

    SetLAThreads(LAThreads);
    if (this->oneCenterK()) {
      for (size_t i = 0; i < nAtoms; ++i) {
        size_t bfStart = mapCen2BfSt[i];
        size_t bfEnd   = (i+1 < nAtoms) ? mapCen2BfSt[i+1] : NB;
        size_t nBf = bfEnd - bfStart;
        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,nBf,nBf,NBNBRI,double(1.),
                   eri3j.pointer()+bfStart*NBNBRI,NBNBRI,
                   Ktemp+bfStart*NBNBRI,NBNBRI,
                   double(0.),AXr+bfStart*NB+bfStart,NB);
      }
    } else {
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB,NB,NBNBRI,double(1.),eri3j.pointer(),NBNBRI,Ktemp,NBNBRI,double(0.),AXr,NB);
    } 

  }; // InCoreRITPIContraction::KContract_real_impl



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  a one-body operator.
   */   
  template <>
  void InCoreRITPIContraction<double, double>::KContract(
      MPI_Comm comm, TwoBodyContraction<double> &C) const {

    InCoreRITPI<double> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<double>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NBRI = eri3j.nRIBasis();

    auto Ktemp = CQMemManager::get().malloc<double>(NB*NB*NBRI);

    KContract_real_impl(comm, C.X, C.AX, Ktemp);

    CQMemManager::get().free(Ktemp);

  }; // InCoreRITPIContraction::KContract



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  a one-body operator.
   */   
  template <>
  void InCoreRITPIContraction<dcomplex, double>::KContract(
      MPI_Comm comm, TwoBodyContraction<dcomplex> &C) const {

    InCoreRITPI<double> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<double>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NBRI = eri3j.nRIBasis();
    size_t NB2 = NB*NB;
    dcomplex *X  = C.X;
    dcomplex *AX = C.AX;

    // Separate real and imaginary parts of density matrix
    double *Xr = CQMemManager::get().malloc<double>(NB2);
    double *Xi = CQMemManager::get().malloc<double>(NB2);
    #pragma omp parallel for
    for(auto k = 0ul; k < NB2; k++) {
      Xr[k] = std::real(X[k]);
      Xi[k] = std::imag(X[k]);
    }

    // Allocate result buffer and scratch
    double *AXr = CQMemManager::get().malloc<double>(NB2);  
    double *AXi = CQMemManager::get().malloc<double>(NB2);
    auto Ktemp = CQMemManager::get().malloc<double>(NB2*NBRI);

    KContract_real_impl(comm, Xr, AXr, Ktemp);
    KContract_real_impl(comm, Xi, AXi, Ktemp);

    // Combine real and imaginary parts
    for(auto k = 0ul; k < NB2; k++) {
      C.AX[k] = dcomplex(AXr[k], AXi[k]);
    }

    CQMemManager::get().free(Xr, Xi);
    CQMemManager::get().free(AXr, AXi);
    CQMemManager::get().free(Ktemp);

  }; // InCoreRITPIContraction<dcomplex, double>::KContract



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  a one-body operator.
   */   
  template <>
  void InCoreRITPIContraction<dcomplex, dcomplex>::KContract(
      MPI_Comm, TwoBodyContraction<dcomplex> &C) const {
    CErr("InCoreRITPIContraction<dcomplex, dcomplex>::KContract not implemented");
  }; // InCoreRITPIContraction<dcomplex, dcomplex>::KContract



  /**
   *  \brief Perform a Exchange-type (23,14) RI-ERI contraction with
   *  orbital coefficients.
   *  K_μν = Σ_{J,i} T^J_μi · T^J_νi*
   *  where T^J_μi = Σ_λ L^J_μλ · C_λi  (half-transformed integral)
   */
  template <typename MatsT, typename IntsT>
  void InCoreRITPIContraction<MatsT, IntsT>::KCoefContract(
      MPI_Comm comm, size_t NO, MatsT *C, MatsT *AX) const {
    ROOT_ONLY(comm);

    InCoreRITPI<IntsT> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NBRI = eri3j.nRIBasis();
    size_t NBNBRI = NB*NBRI;
    size_t NONBRI= NO*NBRI;

    std::fill_n(AX, NB*NB, MatsT(0.));
    MatsT *T = CQMemManager::get().malloc<MatsT>(NBNBRI*NO);

    // 1. T(i, J | nu) = C(lambda, i)^H @ L(J, lambda | nu)^T
    //    One-center: restrict λ to atom(ν)
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      std::fill_n(T, NONBRI * NB, MatsT(0.));
      #pragma omp parallel for schedule(dynamic)
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        for (size_t nu = bfStart; nu < bfEnd; nu++)
          blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::Trans,
                     NO, NBRI, nBf, MatsT(1.),
                     C + bfStart, NB,
                     eri3j.pointer() + nu * NBNBRI + bfStart * NBRI, NBRI,
                     MatsT(0.), T + nu * NONBRI, NO);
      }
    } else {
      #pragma omp parallel for
      for (auto nu = 0ul; nu < NB; nu++)
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::Trans,
                   NO, NBRI, NB, MatsT(1.), C, NB,
                   eri3j.pointer() + nu * NBNBRI, NBRI,
                   MatsT(0.), T + nu * NONBRI, NO);
    }
    SetLAThreads(LAThreads);

    // 2: K(mu, nu) = T(iJ, mu)^T @ T(iJ, nu) 
    //    One-center: per-atom syrk
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI, MatsT(1.),
                   T + bfStart * NONBRI, NONBRI,
                   MatsT(0.), AX + bfStart * NB + bfStart, NB);
        for (size_t q = 0; q < nBf; q++)
          for (size_t p = q + 1; p < nBf; p++)
            AX[(bfStart + q) + (bfStart + p) * NB] =
                AX[(bfStart + p) + (bfStart + q) * NB];
      }
    } else {
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI, MatsT(1.), T, NONBRI, MatsT(0.), AX, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = q + 1; p < NB; p++)
          AX[q + p * NB] = AX[p + q * NB];
    }

#ifdef LEADING_NB2
    // Alternative code for when ERI3J is stored with leading dimension of NB2
    size_t LAThreads = GetLAThreads();
    SetLAThreads(LAThreads);

    // K_{μν} = ∑_α ∑_i ∑_{λσ} L^α_{μλ} · C_{λi} · C_{σi}^* · L^α_{νσ}
    //        = ∑_α ∑_i (∑_{λ} L^α_{μλ} · C_{λi}) · (∑_{σ} C_{σi}^* · L^α_{νσ})*
    for (size_t alpha = 0; alpha < NBRI; ++alpha) {
      // Step 1: T_{μ i} = L_{μλ} ⋅ C_{λ i}
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NO,NB,MatsT(1.0),eri3j.pointer()+alpha*NB*NB,NB,C,NB,MatsT(0.0),Btemp1,NB);
      // Step 2: K += T ⋅ T^H
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NO,MatsT(1.0),Btemp1,NB,Btemp1,NB,MatsT(1.0),AX,NB);
    }
#endif
    CQMemManager::get().free(T);

  }; // InCoreRITPIContraction::KCoefContract



  template <>
  void InCoreRITPIContraction<dcomplex, double>::KCoefContract(
      MPI_Comm comm, size_t NO, dcomplex *C, dcomplex *AX) const {
    ROOT_ONLY(comm);

    InCoreRITPI<double> &eri3j = *std::dynamic_pointer_cast<InCoreRITPI<double>>(this->ints_);
    size_t NB = eri3j.nBasis();
    size_t NBRI = eri3j.nRIBasis();
    size_t NBNBRI = NB*NBRI;
    size_t NONBRI= NO*NBRI;

    std::fill_n(AX, NB*NB, dcomplex(0.));

    // Separate real and imaginary parts of MO coefficients
    double *Cr = CQMemManager::get().malloc<double>(NB*NO);
    double *Ci = CQMemManager::get().malloc<double>(NB*NO);
    #pragma omp parallel for
    for(auto k = 0ul; k < NB*NO; k++) {
      Cr[k] = std::real(C[k]);
      Ci[k] = std::imag(C[k]);
    }

    // Half-transformed buffers (real)
    double *T_re = CQMemManager::get().malloc<double>(NONBRI * NB);
    double *T_im = CQMemManager::get().malloc<double>(NONBRI * NB);

    // 1. T_re(i, J | ν) = Re(C)(λ, i)^T @ L(J, λ | ν)^T
    //    T_im(i, J | ν) = Im(C)(λ, i)^T @ L(J, λ | ν)^T
    //    One-center: restrict λ to atom(ν)
    size_t LAThreads = GetLAThreads();
    SetLAThreads(1);
    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();
      std::fill_n(T_re, NONBRI * NB, double(0.));
      std::fill_n(T_im, NONBRI * NB, double(0.));
      #pragma omp parallel for schedule(dynamic)
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        for (size_t nu = bfStart; nu < bfEnd; nu++) {
          blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                     NO, NBRI, nBf, double(1.),
                     Cr + bfStart, NB,
                     eri3j.pointer() + nu * NBNBRI + bfStart * NBRI, NBRI,
                     double(0.), T_re + nu * NONBRI, NO);
          blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                     NO, NBRI, nBf, double(1.),
                     Ci + bfStart, NB,
                     eri3j.pointer() + nu * NBNBRI + bfStart * NBRI, NBRI,
                     double(0.), T_im + nu * NONBRI, NO);
        }
      }
    } else {
      #pragma omp parallel for
      for (auto nu = 0ul; nu < NB; nu++) {
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                   NO, NBRI, NB, double(1.),
                   Cr, NB, eri3j.pointer() + nu * NBNBRI, NBRI,
                   double(0.), T_re + nu * NONBRI, NO);
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
                   NO, NBRI, NB, double(1.),
                   Ci, NB, eri3j.pointer() + nu * NBNBRI, NBRI,
                   double(0.), T_im + nu * NONBRI, NO);
      }
    }
    SetLAThreads(LAThreads);
    CQMemManager::get().free(Cr, Ci);
    
    double *K_temp = CQMemManager::get().malloc<double>(NB*NB);

    if (this->oneCenterK()) {
      auto mapCen2BfSt = this->mapCen2BfSt();
      size_t nAtoms = mapCen2BfSt.size();

      // 2a. Per-atom: Re(K_A) = T_re_A^T · T_re_A + T_im_A^T · T_im_A
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        size_t off = bfStart * NONBRI;
        size_t koff = bfStart * NB + bfStart;

        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI, double(1.),
                   T_re + off, NONBRI,
                   double(0.), K_temp + koff, NB);
        blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                   nBf, NONBRI, double(1.),
                   T_im + off, NONBRI,
                   double(1.), K_temp + koff, NB);

        for (size_t q = 0; q < nBf; q++)
          for (size_t p = 0; p <= q; p++) {
            double val = K_temp[(bfStart + q) + (bfStart + p) * NB];
            AX[(bfStart + p) + (bfStart + q) * NB] += val;
            if (p != q)
              AX[(bfStart + q) + (bfStart + p) * NB] += val;
          }
      }

      // 2b. Per-atom: Im(K_A) = T_re_A^T · T_im_A - (T_re_A^T · T_im_A)^T
      for (size_t iAtom = 0; iAtom < nAtoms; ++iAtom) {
        size_t bfStart = mapCen2BfSt[iAtom];
        size_t bfEnd = (iAtom + 1 < nAtoms) ? mapCen2BfSt[iAtom + 1] : NB;
        size_t nBf = bfEnd - bfStart;
        size_t off = bfStart * NONBRI;
        size_t koff = bfStart * NB + bfStart;

        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                   nBf, nBf, NONBRI, double(1.),
                   T_re + off, NONBRI,
                   T_im + off, NONBRI,
                   double(0.), K_temp + koff, NB);

        for (size_t q = 0; q < nBf; q++)
          for (size_t p = 0; p <= q; p++) {
            double imval = K_temp[(bfStart + p) + (bfStart + q) * NB]
                         - K_temp[(bfStart + q) + (bfStart + p) * NB];
            AX[(bfStart + p) + (bfStart + q) * NB] += dcomplex(0., imval);
            if (p != q)
              AX[(bfStart + q) + (bfStart + p) * NB] += dcomplex(0., -imval);
          }
      }

    } else {

      // 2a. Re(K) = T_re^T · T_re + T_im^T · T_im
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI, double(1.),
                 T_re, NONBRI,
                 double(0.), K_temp, NB);
      blas::syrk(blas::Layout::ColMajor, blas::Uplo::Lower, blas::Op::Trans,
                 NB, NONBRI, double(1.),
                 T_im, NONBRI,
                 double(1.), K_temp, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = 0; p <= q; p++) {
          double val = K_temp[q + p * NB];
          AX[p + q * NB] += val;
          if (p != q)
            AX[q + p * NB] += val;
        }

      // 2b. Im(K) = T_re^T · T_im - (T_re^T · T_im)^T
      blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                 NB, NB, NONBRI, double(1.),
                 T_re, NONBRI, T_im, NONBRI,
                 double(0.), K_temp, NB);
      for (size_t q = 0; q < NB; q++)
        for (size_t p = 0; p <= q; p++) {
          double imval = K_temp[p + q * NB] - K_temp[q + p * NB];
          AX[p + q * NB] += dcomplex(0., imval);
          if (p != q)
            AX[q + p * NB] += dcomplex(0., -imval);
        }
    }

    CQMemManager::get().free(T_re, T_im, K_temp);

  }; // InCoreRITPIContraction<dcomplex, double>::KCoefContract



  /**
   *  \brief Perform a Coulomb-type (34,12) RI-ERI contraction for Asymmtric ERi, using aux basis
   * 
   * Dimension: 
   *  - If we're calculating J^{ep} = \sum (ee|pp)P^{pp}, then density dimension is NB_Prot * NB_prot, 
   *    output dimension is NB_elec * NB_elec. No need to modify integrals
   * 
   *  - If we're calculating J^{pe} = \sum (pp|ee)P^{ee}, then density dimension is NB_elec * NB_elec,
   *    output dimension is NB_prot * NB_prot. Need to modify integrals to transpose. The flag contract second
   *    in base class TPIContration will be set to True.
   */   
  template <typename MatsT, typename IntsT>
  void InCoreAsymmRITPIContraction<MatsT, IntsT>::JContract(
      MPI_Comm, TwoBodyContraction<MatsT> &C) const {

    // Obtain info from original (ee|pp) ints
    InCoreAsymmRITPI<IntsT> &asymmInts = *std::dynamic_pointer_cast<InCoreAsymmRITPI<IntsT>>(this->ints_);
    size_t NB = asymmInts.nBasis();
    size_t snNB = asymmInts.snBasis();
    std::shared_ptr<InCoreRITPI<IntsT>> aux1 = asymmInts.getAux1();
    std::shared_ptr<InCoreRITPI<IntsT>> aux2 = asymmInts.getAux2();
    if(not aux1 and not aux2) CErr("No aux available in IncoreAsymmRITPIContration::JContract");

    // If contractSecond set to true, need to modify order to create (pp|ee) ints
    if(this->contractSecond){
      std::swap(NB, snNB);
      std::swap(aux1, aux2);
    }

    // If one aux basis is used, then aux dimension is in NBRI 
    // If both aux basis are used, then snNBRI stores the second aux dimension
    size_t NBRI, snNBRI = 0;

    // Define pointers for 3-index tensors to use in contraction
    IntsT* L3J = nullptr;
    IntsT* R3J = nullptr;
    IntsT* M2J = nullptr; 

    if(aux1){
      // if only use aux basis for 1st basis
      L3J = aux1->pointer();
      NBRI = aux1->nRIBasis();
      if (not aux2) R3J = asymmInts.pointer();
    }

    if(aux2){
      R3J = aux2->pointer();
      if(not aux1){
        // if only use aux basis for 2nd basis
        L3J = asymmInts.pointer();
        NBRI = aux2->nRIBasis();
      }else{
        // if use both aux basis  
        if (asymmInts.asymmCDalg() == ASYMM_CD_ALG::CONNECTOR) M2J = asymmInts.M2J();
        snNBRI = NBRI;
        NBRI = aux2->nRIBasis();

        if (M2J == nullptr and snNBRI != NBRI) {
          CErr("Missing middle matrix for LML' type asymmetric RI contraction.");
        }
      }
    }

    // X stores density
    IntsT *X  = reinterpret_cast<IntsT*>(C.X);
    // AX stores output matrix
    IntsT *AX = reinterpret_cast<IntsT*>(C.AX);

    // Extract the real part of X if X is Hermetian and if the ints are real
    const bool extractRealPartX = 
      C.HER and std::is_same<IntsT,double>::value and 
      std::is_same<MatsT,dcomplex>::value;
    
    // Allocate scratch if IntsT and MatsT are different
    const bool allocAXScratch = not std::is_same<IntsT,MatsT>::value;

    if( extractRealPartX ) {
      X = CQMemManager::get().malloc<IntsT>(snNB*snNB);
      for(auto k = 0ul; k < snNB*snNB; k++) X[k] = std::real(C.X[k]);
    }

    if( allocAXScratch ) {
      AX = CQMemManager::get().malloc<IntsT>(NB*NB);
    }
    std::fill_n(AX,NB*NB,0.);

    auto Jtemp = CQMemManager::get().malloc<IntsT>( NBRI );
    std::fill_n(Jtemp, NBRI, IntsT(0.));
    
    // R3J (NBRI by snNB^2) contracting with density (sbNB^2 by 1), generates a vector of length NBRI. 
    // auto gemm1Begin = tick();
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBRI,1,snNB*snNB,IntsT(1.),R3J,NBRI,X,snNB*snNB,IntsT(0.),Jtemp,NBRI);
    // double durGemm1 = tock(gemm1Begin);
    // std::cout << "  Asymm R3J X Density GEMM duration: " << durGemm1 << " s" << std::endl;

    // auto gemm2Begin = tick();
    if( !M2J ){
      // Left multiply by L3J.T (NB^2 by NBRI), to give output matrix (NB^2 by 1). 
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB*NB,1,NBRI,IntsT(1.),L3J,NBRI,Jtemp,NBRI,IntsT(0.),AX,NB*NB);
    } else{
      auto Jtemp1 = CQMemManager::get().malloc<IntsT>( snNBRI );
      std::fill_n(Jtemp1, snNBRI, IntsT(0.));
      // Left multiply by M2J (snNBRI by NBRI), to give temp vector (snNBRI by 1)
      if(this->contractSecond){
        blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,snNBRI,1,NBRI,IntsT(1.),M2J,NBRI,Jtemp,NBRI,IntsT(0.),Jtemp1,snNBRI);
      }else{
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,snNBRI,1,NBRI,IntsT(1.),M2J,snNBRI,Jtemp,NBRI,IntsT(0.),Jtemp1,snNBRI);
      }
      // Left multiply by L3J.T (NB^2 by snNBRI), to give output matrix (NB^2 by 1). 
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB*NB,1,snNBRI,IntsT(1.),L3J,snNBRI,Jtemp1,snNBRI,IntsT(0.),AX,NB*NB);
      CQMemManager::get().free(Jtemp1);
    }
    // double durGemm2 = tock(gemm2Begin);
    // std::cout << "  Asymm L3J X R3JXDensity GEMM duration: " << durGemm2 << " s" << std::endl;
    CQMemManager::get().free(Jtemp);

    // if Complex ints + Hermitian, conjugate
//    if( std::is_same<IntsT,dcomplex>::value and C.HER )
//      IMatCopy('R',NB,NB,IntsT(1.),AX,NB,NB);

    // If non-hermetian, transpose
//    if( not C.HER )  {

//      IMatCopy('T',NB,NB,IntsT(1.),AX,NB,NB);

//    }

    // Cleanup temporaries
    if( extractRealPartX ) CQMemManager::get().free(X);
    if( allocAXScratch ) {

      std::copy_n(AX,NB*NB,C.AX);
      CQMemManager::get().free(AX);

    }

  }; // InCoreAsymmRITPIContraction::JContract

}; // namespace ChronusQ
