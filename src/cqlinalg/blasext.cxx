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
#include <cqlinalg/blasext.hpp>
#include <cerr.hpp>

#include <util/matout.hpp>
#include <util/math.hpp>

extern "C" {

#define mkl_domatcopy MKL_Domatcopy
void MKL_Domatcopy(
    char ordering, char trans,
    size_t rows, size_t cols,
    const double alpha,
    const double * A, size_t lda,
    double * B, size_t ldb);

#define mkl_zomatcopy MKL_Zomatcopy
void MKL_Zomatcopy(
    char ordering, char trans,
    size_t rows, size_t cols,
    const dcomplex alpha,
    const dcomplex * A, size_t lda,
    dcomplex * B, size_t ldb);

#define mkl_domatcopy2 MKL_Domatcopy2
void MKL_Domatcopy2(
    char ordering, char trans,
    size_t rows, size_t cols,
    const double alpha,
    const double * A, size_t lda, size_t stridea,
    double * B, size_t ldb, size_t strideb);

#define mkl_zomatcopy2 MKL_Zomatcopy2
void MKL_Zomatcopy2(
    char ordering, char trans,
    size_t rows, size_t cols,
    const dcomplex alpha,
    const dcomplex * A, size_t lda, size_t stridea,
    dcomplex * B, size_t ldb, size_t strideb);

#define mkl_domatadd MKL_Domatadd
void MKL_Domatadd(
    char ordering, char transa, char transb,
    size_t rows, size_t cols,
    const double alpha,
    const double * A, size_t lda,
    const double beta,
    const double * B, size_t ldb,
    double * C, size_t ldc);

#define mkl_zomatadd MKL_Zomatadd
void MKL_Zomatadd(
    char ordering, char transa, char transb,
    size_t rows, size_t cols,
    const dcomplex alpha,
    const dcomplex * A, size_t lda,
    const dcomplex beta,
    const dcomplex * B, size_t ldb,
    dcomplex * C, size_t ldc);

#define mkl_dimatcopy MKL_Dimatcopy
void MKL_Dimatcopy(
    const char ordering, const char trans,
    size_t rows, size_t cols,
    const double alpha,
    double * AB, size_t lda, size_t ldb);

#define mkl_zimatcopy MKL_Zimatcopy
void MKL_Zimatcopy(
    const char ordering, const char trans,
    size_t rows, size_t cols,
    const dcomplex alpha,
    dcomplex * AB, size_t lda, size_t ldb);
}


namespace ChronusQ {

  // MatAdd helper macros
  #define MATN

  #define MATR \
    .conjugate()

  #define MATT \
    .transpose()

  #define MATC \
    .adjoint()
  
  #define MAT_OP(OP, X) \
    X.block(0,0,M,N) MAT ## OP .template cast<_F3>()

  #define TRY_MAT_ADD(OPA, OPB) \
    if ( TRANSA == #OPA[0] and TRANSB == #OPB[0] ) \
       CMap.block(0,0,M,N).noalias() = \
         ALPHA * MAT_OP(OPA, AMap) + BETA * MAT_OP(OPB, BMap);

  // MatAdd generic template
  template <typename _F1, typename _F2, typename _F3, typename _FScale1, 
    typename _FScale2>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N,
              _FScale1 ALPHA, const _F1 *A, size_t LDA, _FScale2 BETA,
              const _F2 *B, size_t LDB, _F3 *C, size_t LDC){

    if( TRANSA == 'N' and TRANSB == 'N' ) {
      #pragma omp parallel
      {
        const _F1 *locA = A;
        const _F2 *locB = B;
        _F3 *locC = C;
        #pragma omp for
        for(int j = 0; j < N; j++) {
          locA = A + j*LDA;
          locB = B + j*LDB;
          locC = C + j*LDC;
          #pragma omp simd
          for(int i = 0; i < M; i++)
            locC[i] = ALPHA * locA[i] + BETA * locB[i];
        }
      }
    } else if (M == N) {

      if( TRANSA != 'N' ) 
        assert( reinterpret_cast<const void*>(A) != reinterpret_cast<void*>(C) );
      if( TRANSB != 'N' )
        assert( reinterpret_cast<const void*>(B) != reinterpret_cast<void*>(C) );

      Eigen::Map< const Eigen::Matrix<_F1,
          Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > AMap(A,LDA,N);
      Eigen::Map< const Eigen::Matrix<_F2,
          Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > BMap(B,LDB,N);
      Eigen::Map< Eigen::Matrix<_F3,
          Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > CMap(C,LDC,N);

      // Try all different combinations
      TRY_MAT_ADD(N, C)
      else TRY_MAT_ADD(C, N)
      else TRY_MAT_ADD(C, C)
      else TRY_MAT_ADD(N, T)
      else TRY_MAT_ADD(T, N)
      else TRY_MAT_ADD(T, T)
      else TRY_MAT_ADD(C, T)
      else TRY_MAT_ADD(T, C)
      else TRY_MAT_ADD(R, N)
      else TRY_MAT_ADD(N, R)
      else TRY_MAT_ADD(R, R)
      else TRY_MAT_ADD(R, C)
      else TRY_MAT_ADD(C, R)
      else TRY_MAT_ADD(R, T)
      else TRY_MAT_ADD(T, R)

    } else if ( TRANSA == 'N' and TRANSB == 'T') {
      #pragma omp parallel
      {
        const _F1 *locA = A;
        const _F2 *locB = B;
        _F3 *locC = C;
        #pragma omp for
        for(int j = 0; j < N; j++) {
          locA = A + j*LDA;
          locB = B + j;
          locC = C + j*LDC;
          for(int i = 0; i < M; i++)
            locC[i] = ALPHA * locA[i] + BETA * locB[i * LDB];
        }
      }
    } else if ( TRANSA == 'N' and TRANSB == 'C') {
      #pragma omp parallel
      {
        const _F1 *locA = A;
        const _F2 *locB = B;
        _F3 *locC = C;
        #pragma omp for
        for(int j = 0; j < N; j++) {
          locA = A + j*LDA;
          locB = B + j;
          locC = C + j*LDC;
          for(int i = 0; i < M; i++)
            locC[i] = ALPHA * locA[i] + BETA * SmartConj(locB[i * LDB]);
        }
      }
    } else {
      CErr("NYI in generic template");
    }

  }; // MatAdd generic template

#ifdef _CQ_MKL

  template<>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N, double ALPHA, 
    const double *A, size_t LDA, double BETA, const double *B, size_t LDB,
    double *C, size_t LDC) {

      mkl_domatadd('C',TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC);

  }; // MatAdd (real, real, real)

  template<>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N, dcomplex ALPHA, 
    const dcomplex *A, size_t LDA, dcomplex BETA, const dcomplex *B, size_t LDB,
    dcomplex *C, size_t LDC) {

      mkl_zomatadd('C',TRANSA,TRANSB,M,N,ALPHA,A,LDA,BETA,B,LDB,C,LDC);

  }; // MatAdd (complex, complex, complex)

#else

  template void MatAdd( char, char, size_t, size_t, dcomplex, const dcomplex*,
    size_t, dcomplex, const dcomplex*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const double*,
    size_t, double, const double*, size_t, double*, size_t);

#endif

  template<>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N, double ALPHA,
      const dcomplex *A, size_t LDA, double BETA, const dcomplex *B, size_t LDB,
      dcomplex *C, size_t LDC) {
    MatAdd(TRANSA,TRANSB,M,N,dcomplex(ALPHA),A,LDA,dcomplex(BETA),B,LDB,C,LDC);
  }

  template<>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N, dcomplex ALPHA,
      const dcomplex *A, size_t LDA, double BETA, const dcomplex *B, size_t LDB,
      dcomplex *C, size_t LDC) {
    MatAdd(TRANSA,TRANSB,M,N,ALPHA,A,LDA,dcomplex(BETA),B,LDB,C,LDC);
  }

  template<>
  void MatAdd(char TRANSA, char TRANSB, size_t M, size_t N, double ALPHA,
      const dcomplex *A, size_t LDA, dcomplex BETA, const dcomplex *B, size_t LDB,
      dcomplex *C, size_t LDC) {
    MatAdd(TRANSA,TRANSB,M,N,dcomplex(ALPHA),A,LDA,BETA,B,LDB,C,LDC);
  }

  template void MatAdd( char, char, size_t, size_t, double, const dcomplex*,
    size_t, double, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const dcomplex*,
    size_t, dcomplex, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const dcomplex*,
    size_t, double, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const dcomplex*,
    size_t, dcomplex, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const double*,
    size_t, double, const dcomplex*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const double*,
    size_t, double, const dcomplex*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const double*,
    size_t, dcomplex, const dcomplex*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const double*,
    size_t, dcomplex, const dcomplex*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const double*,
    size_t, double, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, double, const double*,
    size_t, dcomplex, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const double*,
    size_t, double, const double*, size_t, dcomplex*, size_t);

  template void MatAdd( char, char, size_t, size_t, dcomplex, const double*,
    size_t, dcomplex, const double*, size_t, dcomplex*, size_t);


  // Generic IMatCopy template
  template <typename _F, typename _FScale>
  void IMatCopy(char TRANS, size_t M, size_t N, _FScale ALPHA, _F *A, 
    size_t LDA, size_t LDB) {


      //assert( LDA == LDB );

      if( TRANS == 'N' and ALPHA == _FScale(1.) ) return;

      Eigen::Map< Eigen::Matrix<_F,
          Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > AMap(A,LDA,N);


      
      if( N == M  and LDA == LDB) {
        if( TRANS == 'T' ) AMap.block(0,0,M,N).transposeInPlace();
        if( TRANS == 'C' ) AMap.block(0,0,M,N).adjointInPlace();
      } else {

        if ( LDB < N ) CErr();

        Eigen::Matrix<_F,
            Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> ANew(LDB,M);

        if( TRANS == 'T' ) ANew = AMap.block(0,0,M,N).transpose();
        if( TRANS == 'C' ) ANew = AMap.block(0,0,M,N).adjoint();

        
        Eigen::Map< Eigen::Matrix<_F,
            Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > AMap1(A,LDB,M);
        
        if( TRANS == 'C' or TRANS == 'T' )
          AMap1 = ANew;

      }

      if( TRANS == 'R' ) 
        AMap.block(0,0,M,N).noalias() = AMap.block(0,0,M,N).conjugate();

      if( ALPHA != _FScale(1.) ) AMap.block(0,0,M,N) *= ALPHA;

  }; // Generic IMatCopy template


#ifdef _CQ_MKL

  template <>
  void IMatCopy(char TRANS, size_t M, size_t N, double ALPHA, double *A, 
    size_t LDA, size_t LDB) {

    mkl_dimatcopy('C',TRANS,M,N,ALPHA,A,LDA,LDB);

  }

  template <>
  void IMatCopy(char TRANS, size_t M, size_t N, dcomplex ALPHA, dcomplex *A, 
    size_t LDA, size_t LDB) {

    mkl_zimatcopy('C',TRANS,M,N,ALPHA,A,LDA,LDB);

  }

#else

  template void IMatCopy(char,size_t,size_t,double,double*,size_t,size_t);
  template void IMatCopy(char,size_t,size_t,dcomplex,dcomplex*,size_t,size_t);

/*
  template <>
  void IMatCopy(char TRANS, size_t M, size_t N, double ALPHA, double *A, 
    size_t LDA, size_t LDB) {

    char ORDER = 'c';

    int iM = M; int iN = N; int iLDA = LDA; int iLDB = LDB; 

    std::cerr << "IMATCOPY\n" << std::endl;

    CBLAS_TRANSPOSE iT = CblasNoTrans;
    if( TRANS == 'T' ) iT = CblasTrans;
    if( TRANS == 'C' ) iT = CblasConjTrans;
    if( TRANS == 'R' ) iT = CblasConjNoTrans;

//  dimatcopy_(&ORDER,&TRANS,&iM,&iN,&ALPHA,A,&iLDA,&iLDB);
    cblas_dimatcopy(CblasColMajor,iT,iM,iN,ALPHA,A,iLDA,iLDB);

  }

  template <>
  void IMatCopy(char TRANS, size_t M, size_t N, dcomplex ALPHA, dcomplex *A, 
    size_t LDA, size_t LDB) {

    char ORDER = 'C';

    int iM = M; int iN = N; int iLDA = LDA; int iLDB = LDB; 

    zimatcopy_(&ORDER,&TRANS,&iM,&iN,reinterpret_cast<double*>(&ALPHA),
      reinterpret_cast<double*>(A),&iLDA,&iLDB);

  }
*/
#endif
  
  template <>
  void IMatCopy(char TRANS, size_t M, size_t N, double ALPHA, dcomplex *A, 
    size_t LDA, size_t LDB) {

    IMatCopy(TRANS,M,N,dcomplex(ALPHA),A,LDA,LDB); 

  }


  // Generic HerMat template
  template <typename _F>
  void HerMat(char UPLO, size_t N, _F *A, size_t LDA) {

    assert( UPLO == 'L' or UPLO == 'U' );

    Eigen::Map< Eigen::Matrix<_F,
        Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > AMap(A,LDA,N);

    if( UPLO == 'L' )  
      AMap.block(0,0,N,N)= 
        AMap.block(0,0,N,N).template selfadjointView<Eigen::Lower>();
    else
      AMap.block(0,0,N,N) = 
        AMap.block(0,0,N,N).template selfadjointView<Eigen::Upper>();


    if(std::is_same<_F,dcomplex>::value)
    for(auto i = 0; i < N; i++)
      A[i + i*LDA] = std::real(A[i + i*LDA]);

  }; // Generic HerMat template

  template void HerMat(char, size_t, double*, size_t);
  template void HerMat(char, size_t, dcomplex*, size_t);

}; // namespace ChronusQ
