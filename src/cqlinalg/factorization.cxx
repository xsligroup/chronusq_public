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
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/util.hpp>
#include <lapack.hh>
#include <lapack/fortran.h>

#include <cerr.hpp>
#include <util/matout.hpp>

namespace ChronusQ {


  // Real wraps DGETRF + DGETRI
  // Complex wraps ZGETRF + ZGETRI
  template <typename T>
  int LUInv(int N, T* A, int LDA) {

    int64_t *IPIV = CQMemManager::get().malloc<int64_t>(N);

    int64_t INFO = lapack::getrf(N,N,A,LDA,IPIV);

    if( INFO != 0 ) { CQMemManager::get().free(IPIV); return INFO; }

    lapack::getri(N,A,LDA,IPIV);
    
    CQMemManager::get().free(IPIV);

    return INFO;

  }; // LUInv

  template int LUInv<double>(int N, double* A, int LDA);
  template int LUInv<dcomplex>(int N, dcomplex* A, int LDA);


  template <typename T>
  int QR(int M, int N, T* A, int LDA, T* R, int LDR) {

    T *TAU = CQMemManager::get().malloc<T>(N);

    int INFO = lapack::geqrf(M,N,A,LDA,TAU);

    if( INFO != 0 ) { CQMemManager::get().free(TAU); return INFO; }

    int rCol = std::min(M,N);
    std::fill_n(R,rCol*LDR,0.);

    for(auto j = 0; j < N; j++) 
    for(auto i = 0; i <= j; i++)
      R[i + j*LDR] = A[i + j*LDA];


    INFO = lapack::ungqr(M,N,N,A,LDA,TAU);

    CQMemManager::get().free(TAU); 
    
    return INFO; 

  }
  template int QR<double>(int M, int N, double* A, int LDA, double* R, 
    int LDR);
  template int QR<dcomplex>(int M, int N, dcomplex* A, int LDA, dcomplex* R, 
    int LDR);





  template <typename _F>
  int TGEXC(bool WANTQ, bool WANTZ, int N, _F *A, int LDA, _F *B, int LDB,
    _F *Q, int LDQ, _F *Z, int LDZ, int IFST, int ILST);

  template<>
  int TGEXC(bool WANTQ, bool WANTZ, int N, double *A, int LDA, double *B, 
    int LDB, double *Q, int LDQ, double *Z, int LDZ, int IFST, int ILST) {

    int WANTQ_i = int(WANTQ);
    int WANTZ_i = int(WANTZ);

    int INFO;

    using namespace std::placeholders;
    auto gexc = std::bind(dtgexc_,&WANTQ_i,&WANTZ_i,&N,A,&LDA,B,&LDA,
      Q,&LDQ,Z,&LDZ,&IFST,&ILST,_1,_2,&INFO);

    int LWORK = getLWork<double>(gexc);
    double *WORK = CQMemManager::get().malloc<double>(LWORK);

    gexc(WORK,&LWORK);

    CQMemManager::get().free(WORK);

    return INFO;
  }

  template<>
  int TGEXC(bool WANTQ, bool WANTZ, int N, dcomplex *A, int LDA, dcomplex *B, 
    int LDB, dcomplex *Q, int LDQ, dcomplex *Z, int LDZ, int IFST, int ILST) {

    int WANTQ_i = int(WANTQ);
    int WANTZ_i = int(WANTZ);

    int INFO;

    ztgexc_(&WANTQ_i,&WANTZ_i,&N,A,&LDA,B,&LDA,Q,&LDQ,Z,&LDZ,
      &IFST,&ILST,&INFO);

    return INFO;
  }




  struct TGSEN_OUT {

    int INFO;
    int M;

    double PL;
    double PR;
    double DIF;

  };


  template <typename _F>
  TGSEN_OUT TGSEN(int IJOB, bool WANTQ, bool WANTZ, int * SELECT, int N, _F *A, 
    int LDA, _F *B, int LDB, dcomplex *ALPHA, _F *BETA, _F *Q, int LDQ, _F *Z, 
    int LDZ);


  template <>
  TGSEN_OUT TGSEN(int IJOB, bool WANTQ, bool WANTZ, int * SELECT, int N, 
    double *A, int LDA, double *B, int LDB, dcomplex *ALPHA, double *BETA, 
    double *Q, int LDQ, double *Z, int LDZ) {

    if( IJOB != 0 ) CErr("No proper logic for TGSEN.IJOB != 0");

    int WANTQ_i = int(WANTQ);
    int WANTZ_i = int(WANTZ);

    int INFO;

    int LIWORK = 1;
    int IWORK  = 0;

    TGSEN_OUT out;

    double *ALPHAR = CQMemManager::get().malloc<double>(N); 
    double *ALPHAI = CQMemManager::get().malloc<double>(N); 

    using namespace std::placeholders;
    auto gsen = std::bind(dtgsen_,&IJOB,&WANTQ_i,&WANTZ_i,SELECT,&N,A,&LDA,
      B,&LDB,ALPHAR,ALPHAI,BETA,Q,&LDQ,Z,&LDZ,&out.M,&out.PL,&out.PR,
      &out.DIF,_1,_2,&IWORK,&LIWORK,&INFO);

    int LWORK = getLWork<double>(gsen);
    double *WORK = CQMemManager::get().malloc<double>(LWORK);

    gsen(WORK,&LWORK);


    for(auto k = 0; k < N; k++)
      ALPHA[k] = dcomplex(ALPHAR[k],ALPHAI[k]);

    CQMemManager::get().free(ALPHAR,ALPHAI,WORK);

    return out;
  }



  template <>
  TGSEN_OUT TGSEN(int IJOB, bool WANTQ, bool WANTZ, int * SELECT, int N, 
    dcomplex *A, int LDA, dcomplex *B, int LDB, dcomplex *ALPHA, 
    dcomplex *BETA, dcomplex *Q, int LDQ, dcomplex *Z, int LDZ) {

    if( IJOB != 0 ) CErr("No proper logic for TGSEN.IJOB != 0");

    int WANTQ_i = int(WANTQ);
    int WANTZ_i = int(WANTZ);

    int INFO;

    int LIWORK = 1;
    int IWORK  = 0;

    TGSEN_OUT out;

    using namespace std::placeholders;
    auto gsen = std::bind(ztgsen_,&IJOB,&WANTQ_i,&WANTZ_i,SELECT,&N,A,&LDA,
      B,&LDB,ALPHA,BETA,Q,&LDQ,Z,&LDZ,&out.M,&out.PL,&out.PR,
      &out.DIF,_1,_2,&IWORK,&LIWORK,&INFO);

    int LWORK = getLWork<dcomplex>(gsen);
    dcomplex *WORK = CQMemManager::get().malloc<dcomplex>(LWORK);

    gsen(WORK,&LWORK);

    CQMemManager::get().free(WORK);

    return out;
  }

  template <typename _F>
  int OrdQZ2(char JOBVSL, char JOBVSR, int N, _F *A, int LDA, _F *B, 
    int LDB, dcomplex *ALPHA, _F *BETA, double hLim, double SIGMA, _F *VSL, int LDVSL, 
    _F *VSR, int LDVSR) {

    // Aux booleans for reordering functions
    bool WantQ = JOBVSL == 'V';
    bool WantZ = JOBVSR == 'V';

    double LARGE = 1.0e10;

    // Convert char to lapackpp friendly input
    lapack::Job JVL;
    lapack::Job JVR;

    if(JOBVSL == 'V')       JVL = lapack::Job::Vec;
    else if(JOBVSL == 'N')  JVL = lapack::Job::NoVec;
    else                    CErr("Invalid option for JOBVSL ( lapack::gges )");

    if(JOBVSR == 'V')       JVR = lapack::Job::Vec;
    else if(JOBVSR == 'N')  JVR = lapack::Job::NoVec;
    else                    CErr("Invalid option for JOBVSR ( lapack::gges )");


    // Get the initial QZ factorization
    int64_t SDIM = 0;
    int INFO = lapack::gges(JVL,JVR,lapack::Sort::NotSorted,nullptr,N,A,LDA,B,LDB,&SDIM,ALPHA,BETA,VSL,LDVSL,VSR,LDVSR);

    if( INFO != 0 ) CErr("QZ (lapack::gges) Failed in OrdQZ2");
    
    bool swap = true;

    int * SELECT = CQMemManager::get().malloc<int>(N);
    std::fill_n(SELECT,N,int(false));



    double betaTol = 1e-13;
    for(int i = 1; i < N; i++) 
    for(int j = i; j > 0; j--) {

        bool perfSwap = false;
        bool jm1Small = std::abs(BETA[j-1]) < betaTol;
        bool jSmall   = std::abs(BETA[j])   < betaTol;

        if( jSmall )        perfSwap = false; // if b[j] small, no swap
        else if( jm1Small ) perfSwap = true;  // swap if b[j-1] small
        else {

          // if (W[j-1] - SIGMA) > (W[j] - SIGMA), swap
          double Wjm1 = std::abs(ALPHA[j-1]/BETA[j-1] - SIGMA);
          double Wj   = std::abs(ALPHA[j]/BETA[j]     - SIGMA);
          if (std::real(ALPHA[j-1]/BETA[j-1]) < hLim) {
            Wjm1 += LARGE;
          }
          if (std::real(ALPHA[j]/BETA[j]) < hLim) {
            Wj += LARGE;
          }
          perfSwap = (Wjm1 > Wj) and (std::abs(Wjm1 - Wj) > 1e-10) ;

        }

        if( perfSwap ) {

          int ifst = j + 1;
          int ilst = j;

          TGEXC(WantQ,WantZ,N,A,LDA,B,LDA,VSL,LDVSL,VSR,LDVSR,ifst,
            ilst);

          // Recompute the Generalized eigenvalues from the permuted
          // Shur form
          TGSEN(0,WantQ,WantZ,SELECT,N,A,LDA,B,LDB,ALPHA,BETA,VSL,LDVSL,VSR,
            LDVSR);

        }

    }


    CQMemManager::get().free(SELECT);

    return INFO;

  }


  template <typename _F>
  int OrdQZ(char JOBVSL, char JOBVSR, int N, _F *A, int LDA, _F *B, 
    int LDB, dcomplex *ALPHA, _F *BETA, double SIGMA, _F *VSL, int LDVSL, 
    _F *VSR, int LDVSR) {

    // Aux booleans for reordering functions
    bool WantQ = JOBVSL == 'V';
    bool WantZ = JOBVSR == 'V';

    // Convert char to lapackpp friendly input
    lapack::Job JVL;
    lapack::Job JVR;

    if(JOBVSL == 'V')       JVL = lapack::Job::Vec;
    else if(JOBVSL == 'N')  JVL = lapack::Job::NoVec;
    else                    CErr("Invalid option for JOBVSL ( lapack::gges )");

    if(JOBVSR == 'V')       JVR = lapack::Job::Vec;
    else if(JOBVSR == 'N')  JVR = lapack::Job::NoVec;
    else                    CErr("Invalid option for JOBVSR ( lapack::gges )");


    // Get the initial QZ factorization
    int64_t SDIM = 0;
    int INFO = lapack::gges(JVL,JVR,lapack::Sort::NotSorted,nullptr,N,A,LDA,B,LDB,&SDIM,ALPHA,BETA,VSL,LDVSL,VSR,LDVSR);

    if( INFO != 0 ) CErr("QZ (lapack::gges) Failed in OrdQZ");
    
    bool swap = true;

    int * SELECT = CQMemManager::get().malloc<int>(N);
    std::fill_n(SELECT,N,int(false));

#if 0
    int nCount = 0;
    do {

      // Catch if something stupid has happened
      if(nCount++ > N) CErr("OrdQZ Failed!");

      if( swap ) {

        swap = false;
        double betaTol = 1e-13;
        for(auto j = 0; j < (N-1); j++) {

          bool perfSwap = false;
          bool jSmall  = std::abs(BETA[j])   < betaTol;
          bool j1Small = std::abs(BETA[j+1]) < betaTol;

          if( j1Small )     perfSwap = false; // if b[j+1] small, no swap
          else if( jSmall ) perfSwap = true;  // swap if b[j] small
          else {

            // if (W[j] - SIGMA) > (W[j+1] - SIGMA), swap
            double Wj  = std::abs(ALPHA[j]/BETA[j]     - SIGMA);
            double Wj1 = std::abs(ALPHA[j+1]/BETA[j+1] - SIGMA);
            perfSwap = (Wj > Wj1) and (std::abs(Wj - Wj1) > 1e-10) ;

          }


          if( perfSwap ) {

            swap = true;
            int ifst = j + 2;
            int ilst = j + 1;

            TGEXC(WantQ,WantZ,N,A,LDA,B,LDA,VSL,LDVSL,VSR,LDVSR,ifst,
              ilst);

          }
        }

        // Recompute the Generalized eigenvalues from the permuted
        // Shur form
        TGSEN(0,WantQ,WantZ,SELECT,N,A,LDA,B,LDB,ALPHA,BETA,VSL,LDVSL,VSR,
          LDVSR);


      }


    } while(swap);
#else


    double betaTol = 1e-13;
    for(int i = 1; i < N; i++) 
    for(int j = i; j > 0; j--) {

        bool perfSwap = false;
        bool jm1Small = std::abs(BETA[j-1]) < betaTol;
        bool jSmall   = std::abs(BETA[j])   < betaTol;

        if( jSmall )        perfSwap = false; // if b[j] small, no swap
        else if( jm1Small ) perfSwap = true;  // swap if b[j-1] small
        else {

          // if (W[j-1] - SIGMA) > (W[j] - SIGMA), swap
          double Wjm1 = std::abs(ALPHA[j-1]/BETA[j-1] - SIGMA);
          double Wj   = std::abs(ALPHA[j]/BETA[j]     - SIGMA);
          perfSwap = (Wjm1 > Wj) and (std::abs(Wjm1 - Wj) > 1e-10) ;

        }

        if( perfSwap ) {

          int ifst = j + 1;
          int ilst = j;

          TGEXC(WantQ,WantZ,N,A,LDA,B,LDA,VSL,LDVSL,VSR,LDVSR,ifst,
            ilst);

          // Recompute the Generalized eigenvalues from the permuted
          // Shur form
          TGSEN(0,WantQ,WantZ,SELECT,N,A,LDA,B,LDB,ALPHA,BETA,VSL,LDVSL,VSR,
            LDVSR);

        }

    }

#endif

    CQMemManager::get().free(SELECT);

    return INFO;

  }

  template 
  int OrdQZ2(char JOBVSL, char JOBVSR, int N, double *A, int LDA, double *B, 
    int LDB, dcomplex *ALPHA, double *BETA, double hLim, double SIMGA, double *VSL, 
    int LDVSL, double *VSR, int LDVSR);

  template 
  int OrdQZ2(char JOBVSL, char JOBVSR, int N, dcomplex *A, int LDA, dcomplex *B, 
    int LDB, dcomplex *ALPHA, dcomplex *BETA, double hLim, double SIMGA, dcomplex *VSL, 
    int LDVSL, dcomplex *VSR, int LDVSR);


  template 
  int OrdQZ(char JOBVSL, char JOBVSR, int N, double *A, int LDA, double *B, 
    int LDB, dcomplex *ALPHA, double *BETA, double SIMGA, double *VSL, 
    int LDVSL, double *VSR, int LDVSR);

  template 
  int OrdQZ(char JOBVSL, char JOBVSR, int N, dcomplex *A, int LDA, dcomplex *B, 
    int LDB, dcomplex *ALPHA, dcomplex *BETA, double SIMGA, dcomplex *VSL, 
    int LDVSL, dcomplex *VSR, int LDVSR);

  template<typename MatsT>
  void SVDInverse(const size_t N, MatsT* A, const size_t LDA, const double num) {
  
    // Compute SVD to determine which vectors are singular
    MatsT* U = CQMemManager::get().malloc<MatsT>(N*N);
    MatsT* VT = CQMemManager::get().malloc<MatsT>(N*N);
    double* sV = CQMemManager::get().malloc<double>(N);
  
    int info = lapack::gesvd(lapack::Job::AllVec, lapack::Job::AllVec, N, N, A, LDA, 
            sV, U, N, VT, N);
    if( info != 0 ) throw std::runtime_error("SVD Failed in SVD Inverse Function");
  
    // Zero out vectors that are singular or scale by inverse sing. value
    for( size_t i = 0; i < N; i++ ) {
      if( sV[i] > num ) {
        blas::scal(N, MatsT(1. / sV[i]), U + i * N, 1);
      } else {
        blas::scal(N, MatsT(0.),  U + i * N, 1);
        blas::scal(N, MatsT(0.), VT + i * N, 1);
      }
    }
  
    // Compute Matrix Inverse
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        N, N, N, 
        MatsT(1.), U, N, 
        VT, N, 
        MatsT(0.), A, LDA);
  
    CQMemManager::get().free(U,VT,sV);
  }

  template void SVDInverse(const size_t N, double* A, const size_t LDA, const double num);
  template void SVDInverse(const size_t N, dcomplex* A, const size_t LDA, const double num);

}; // namespace ChronusQ
