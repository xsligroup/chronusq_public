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
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blasext.hpp>
#include <cerr.hpp>

#include <util/matout.hpp>

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
}

namespace ChronusQ {

  template <typename T> T ComplexScale();
  template <> double ComplexScale(){ return -1; }
  template <> dcomplex ComplexScale(){ return dcomplex(0,1); }

  template <typename _F1, typename _F2>
  void SpinScatter(size_t M, size_t N, const _F1 *AA, size_t LDAA,
      const _F1 *AB, size_t LDAB, const _F1 *BA, size_t LDBA,
      const _F1 *BB, size_t LDBB, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA, bool BBeqAA) {

    _F2 YFACT = ComplexScale<_F2>();

    if (BBeqAA) {
      if(AS) SetMat('N',M,N,_F2(2.),AA,LDAA,AS,LDAS);
      if(AZ) SetMat('N',M,N,_F2(0.0),AA,LDAA,AZ,LDAZ);
    } else {
      if(AS) MatAdd('N','N',M,N,_F2(1.),AA,LDAA,_F2(1.) ,BB,LDBB,AS,LDAS);
      if(AZ) MatAdd('N','N',M,N,_F2(1.),AA,LDAA,_F2(-1.),BB,LDBB,AZ,LDAZ);
    }

    if (zeroABBA) {
      if(AY) SetMat('N',M,N,_F2(0.0),AA,LDAA,AY,LDAY);
      if(AX) SetMat('N',M,N,_F2(0.0),AA,LDAA,AX,LDAX);
    } else {
      if(AY) MatAdd('N','N',M,N, YFACT ,AB,LDAB, -YFACT ,BA,LDBA,AY,LDAY);
      if(AX) MatAdd('N','N',M,N,_F2(1.),AB,LDAB,_F2(1.) ,BA,LDBA,AX,LDAX);
    }

  };

  template <typename _F1, typename _F2>
  void SpinScatter(size_t M, size_t N, const _F1 *A, size_t LDA, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA, bool BBeqAA) {

/*
    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < M; i++) {
      AS[i + j*LDAS] = A[i + j*LDA] + A[(i+M) + (j+N)*LDA];
      AZ[i + j*LDAZ] = A[i + j*LDA] - A[(i+M) + (j+N)*LDA];
      AY[i + j*LDAY] = YFACT * (A[i + (j+N)*LDA] - A[(i+M) + j*LDA]);
      AX[i + j*LDAX] = A[i + (j+N)*LDA] + A[(i+M) + j*LDA];
    }
*/

    const _F1* A_AA = A;
    const _F1* A_AB = A_AA + N*LDA;
    const _F1* A_BA = A_AA + M;
    const _F1* A_BB = A_AB + M;

    SpinScatter(M,N,A_AA,LDA,A_AB,LDA,A_BA,LDA,A_BB,LDA,
                AS,LDAS,AZ,LDAZ,AY,LDAY,AX,LDAX,zeroABBA,BBeqAA);

  };

  template <typename _F1, typename _F2>
  void SpinScatter(size_t N, const _F1 *A, size_t LDA, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA, bool BBeqAA) {

   SpinScatter(N,N,A,LDA,AS,LDAS,AZ,LDAZ,AY,LDAY,AX,LDAX,zeroABBA,BBeqAA);

  };

  template <typename _F1, typename _F2>
  void SpinGather(size_t M, size_t N, _F1 *AA, size_t LDAA, _F1 *AB, size_t LDAB,
      _F1 *BA, size_t LDBA, _F1 *BB, size_t LDBB, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY, bool zeroZ) {

    _F2 YFACT = 0.5*ComplexScale<_F2>();

    if (zeroZ) {
      if(AA) SetMat('N',M,N,_F2(0.5),AS,LDAS,AA,LDAA);
      if(BB) SetMat('N',M,N,_F2(0.5),AS,LDAS,BB,LDBB);
    } else {
      if(AA) MatAdd('N','N',M,N,_F2(0.5),AS,LDAS,_F2(0.5) ,AZ,LDAZ,AA,LDAA);
      if(BB) MatAdd('N','N',M,N,_F2(0.5),AS,LDAS,_F2(-0.5),AZ,LDAZ,BB,LDBB);
    }

    if (zeroXY) {
      if(BA) SetMat('N',M,N,_F2(0.0),AS,LDAS,BA,LDBA);
      if(AB) SetMat('N',M,N,_F2(0.0),AS,LDAS,AB,LDAB);
    } else {
      if(BA) MatAdd('N','N',M,N,_F2(0.5),AX,LDAX, YFACT,AY,LDAY,BA,LDBA);
      if(AB) MatAdd('N','N',M,N,_F2(0.5),AX,LDAX,-YFACT,AY,LDAY,AB,LDAB);
    }

  };

  template <typename _F1, typename _F2>
  void SpinGather(size_t M, size_t N, _F1 *A, size_t LDA, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY, bool zeroZ) {

/*
    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < N; i++) {
      A[i + j*LDA]         = 0.5 * (AS[i + j*LDAS] + AZ[i + j*LDAZ]);
      A[(i+N) + (j+N)*LDA] = 0.5 * (AS[i + j*LDAS] - AZ[i + j*LDAZ]);
      A[(i+N) + j*LDA]     = 0.5 * (AX[i + j*LDAS] + YFACT * AY[i + j*LDAZ]);
      A[i + (j+N)*LDA]     = 0.5 * (AX[i + j*LDAS] - YFACT * AY[i + j*LDAZ]);
    }
*/
    _F1* A_AA = A;
    _F1* A_AB = A_AA + N*LDA;
    _F1* A_BA = A_AA + M;
    _F1* A_BB = A_AB + M;

    SpinGather(M,N,A_AA,LDA,A_AB,LDA,A_BA,LDA,A_BB,LDA,
                AS,LDAS,AZ,LDAZ,AY,LDAY,AX,LDAX,zeroXY,zeroZ);

  };

  template <typename _F1, typename _F2>
  void SpinGather(size_t N, _F1 *A, size_t LDA, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY, bool zeroZ) {

    SpinGather(N,N,A,LDA,AS,LDAS,AZ,LDAZ,AY,LDAY,AX,LDAX,zeroXY,zeroZ);

  };

  template
  void SpinScatter(size_t M, size_t N, const double *AA, size_t LDAA,
      const double *AB, size_t LDAB, const double *BA, size_t LDBA,
      const double *BB, size_t LDBB, double *AS, size_t LDAS,
      double *AZ, size_t LDAZ, double *AY, size_t LDAY, double *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t M, size_t N, const double *AA, size_t LDAA,
      const double *AB, size_t LDAB, const double *BA, size_t LDBA,
      const double *BB, size_t LDBB, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t M, size_t N, const dcomplex *AA, size_t LDAA,
      const dcomplex *AB, size_t LDAB, const dcomplex *BA, size_t LDBA,
      const dcomplex *BB, size_t LDBB, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t M, size_t N, const double *A, size_t LDA, double *AS, size_t LDAS,
      double *AZ, size_t LDAZ, double *AY, size_t LDAY, double *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t M, size_t N, const double *A, size_t LDA, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t M, size_t N, const dcomplex *A, size_t LDA, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t N, const double *A, size_t LDA, double *AS, size_t LDAS,
      double *AZ, size_t LDAZ, double *AY, size_t LDAY, double *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t N, const double *A, size_t LDA, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinScatter(size_t N, const dcomplex *A, size_t LDA, dcomplex *AS, size_t LDAS,
      dcomplex *AZ, size_t LDAZ, dcomplex *AY, size_t LDAY, dcomplex *AX, size_t LDAX,
      bool zeroABBA, bool zeroBB);

  template
  void SpinGather(size_t M, size_t N, double *AA, size_t LDAA, double *AB, size_t LDAB,
      double *BA, size_t LDBA, double *BB, size_t LDBB, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t M, size_t N, dcomplex *AA, size_t LDAA, dcomplex *AB, size_t LDAB,
      dcomplex *BA, size_t LDBA, dcomplex *BB, size_t LDBB, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t M, size_t N, dcomplex *AA, size_t LDAA, dcomplex *AB, size_t LDAB,
      dcomplex *BA, size_t LDBA, dcomplex *BB, size_t LDBB, const dcomplex *AS, size_t LDAS,
      const dcomplex *AZ, size_t LDAZ, const dcomplex *AY, size_t LDAY, const dcomplex *AX,
      size_t LDAX, bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t M, size_t N, double *A, size_t LDA, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t M, size_t N, dcomplex *A, size_t LDA, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t M, size_t N, dcomplex *A, size_t LDA, const dcomplex *AS, size_t LDAS,
      const dcomplex *AZ, size_t LDAZ, const dcomplex *AY, size_t LDAY, const dcomplex *AX,
      size_t LDAX, bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t N, double *A, size_t LDA, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t N, dcomplex *A, size_t LDA, const double *AS, size_t LDAS,
      const double *AZ, size_t LDAZ, const double *AY, size_t LDAY, const double *AX, size_t LDAX,
      bool zeroXY, bool zeroZ);

  template
  void SpinGather(size_t N, dcomplex *A, size_t LDA, const dcomplex *AS, size_t LDAS,
      const dcomplex *AZ, size_t LDAZ, const dcomplex *AY, size_t LDAY, const dcomplex *AX,
      size_t LDAX, bool zeroXY, bool zeroZ);

  /*
   *  Component Scatter a square matrix
   *
   *  A ==> [ A_LL, A_LS ] 
   *        [ A_SL, A_SS ]
   *
   */
  
  template <typename _F1, typename _F2>
  void ComponentScatter(size_t NL, size_t NS, const _F1 *A, size_t LDA,
      _F2 ScaleLL, _F2 *ALL, size_t LDALL, _F2 ScaleLS, _F2 *ALS, size_t LDALS,
      _F2 ScaleSL, _F2 *ASL, size_t LDASL, _F2 ScaleSS, _F2 *ASS, size_t LDASS,
      bool increment) { 

    size_t LS = NL*LDA;
    size_t SL = NL;
    size_t SS = LS + SL;
    
    if (increment) {
      if (ALL) MatAdd('N','N',NL,NL,ScaleLL,A   ,LDA,_F2(1.),ALL,LDALL,ALL,LDALL); 
      if (ALS) MatAdd('N','N',NL,NS,ScaleLS,A+LS,LDA,_F2(1.),ALS,LDALS,ALS,LDALS); 
      if (ASL) MatAdd('N','N',NS,NL,ScaleSL,A+SL,LDA,_F2(1.),ASL,LDASL,ASL,LDASL); 
      if (ASS) MatAdd('N','N',NS,NS,ScaleSS,A+SS,LDA,_F2(1.),ASS,LDASS,ASS,LDASS); 
    } else {
      if (ALL) SetMat('N',NL,NL,ScaleLL,A   ,LDA,ALL,LDALL); 
      if (ALS) SetMat('N',NL,NS,ScaleLS,A+LS,LDA,ALS,LDALS); 
      if (ASL) SetMat('N',NS,NL,ScaleSL,A+SL,LDA,ASL,LDASL); 
      if (ASS) SetMat('N',NS,NS,ScaleSS,A+SS,LDA,ASS,LDASS); 
    }
  };
  
  template <typename _F1, typename _F2>
  void ComponentScatter(size_t NL, size_t NS, const _F1 *A, size_t LDA,
      _F2 *ALL, size_t LDALL, _F2 *ALS, size_t LDALS,
      _F2 *ASL, size_t LDASL, _F2 *ASS, size_t LDASS,
      bool increment) { 

    ComponentScatter(NL,NS,A,LDA,_F2(1.),ALL,LDALL,_F2(1.),ALS,LDALS,
      _F2(1.),ASL,LDASL,_F2(1.),ASS,LDASS,increment);
  };

  /*
   *  Component Scatter a square matrix
   *
   *  [ A_LL, A_LS ]  ==>  A
   *  [ A_SL, A_SS ]       
   *
   */
  template <typename _F1, typename _F2>
  void ComponentGather(size_t NL, size_t NS, _F1 *A, size_t LDA,
    char TransLL, _F2 ScaleLL, const _F2 *ALL, size_t LDALL, 
    char TransLS, _F2 ScaleLS, const _F2 *ALS, size_t LDALS,
    char TransSL, _F2 ScaleSL, const _F2 *ASL, size_t LDASL, 
    char TransSS, _F2 ScaleSS, const _F2 *ASS, size_t LDASS,
    bool increment) {
  
    size_t LS = NL*LDA;
    size_t SL = NL;
    size_t SS = LS + SL;
  
    if (increment) {
      if (ALL) MatAdd(TransLL,'N',NL,NL,ScaleLL,ALL,LDALL,_F2(1.),A   ,LDA,A   ,LDA); 
      if (ALS) MatAdd(TransLS,'N',NL,NS,ScaleLS,ALS,LDALS,_F2(1.),A+LS,LDA,A+LS,LDA); 
      if (ASL) MatAdd(TransSL,'N',NS,NL,ScaleSL,ASL,LDASL,_F2(1.),A+SL,LDA,A+SL,LDA); 
      if (ASS) MatAdd(TransSS,'N',NS,NS,ScaleSS,ASS,LDASS,_F2(1.),A+SS,LDA,A+SS,LDA); 
    } else {
      if (ALL) SetMat(TransLL,NL,NL,ScaleLL,ALL,LDALL,A   ,LDA); 
      if (ALS) SetMat(TransLS,NL,NS,ScaleLS,ALS,LDALS,A+LS,LDA); 
      if (ASL) SetMat(TransSL,NS,NL,ScaleSL,ASL,LDASL,A+SL,LDA); 
      if (ASS) SetMat(TransSS,NS,NS,ScaleSS,ASS,LDASS,A+SS,LDA); 
    }
  };

  template <typename _F1, typename _F2>
  void ComponentGather(size_t NL, size_t NS, _F1 *A, size_t LDA,
    const _F2 *ALL, size_t LDALL, const _F2 *ALS, size_t LDALS,
    const _F2 *ASL, size_t LDASL, const _F2 *ASS, size_t LDASS,
    bool increment) {
    
    ComponentGather(NL,NS,A,LDA,'N',_F2(1.),ALL,LDALL,'N',_F2(1.),ALS,LDALS,
      'N',_F2(1.),ASL,LDASL,'N',_F2(1.),ASS,LDASS,increment);

  };

  template 
  void ComponentScatter(size_t NL, size_t NS, const double *A, size_t LDA, 
    double ScaleLL, double *ALL, size_t LDALL, double ScaleLS, double *ALS, size_t LDALS, 
    double ScaleSL, double *ASL, size_t LDASL, double ScaleSS, double *ASS, size_t LDASS, 
    bool increment);
  
  template 
  void ComponentScatter(size_t NL, size_t NS, const double *A, size_t LDA, 
    dcomplex ScaleLL, dcomplex *ALL, size_t LDALL, dcomplex ScaleLS, dcomplex *ALS, size_t LDALS, 
    dcomplex ScaleSL, dcomplex *ASL, size_t LDASL, dcomplex ScaleSS, dcomplex *ASS, size_t LDASS, 
    bool increment);
  
  template 
  void ComponentScatter(size_t NL, size_t NS, const dcomplex *A, size_t LDA, 
    dcomplex ScaleLL, dcomplex *ALL, size_t LDALL, dcomplex ScaleLS, dcomplex *ALS, size_t LDALS, 
    dcomplex ScaleSL, dcomplex *ASL, size_t LDASL, dcomplex ScaleSS, dcomplex *ASS, size_t LDASS, 
    bool increment);
  
  template 
  void ComponentScatter(size_t NL, size_t NS, const double *A, size_t LDA, 
    double *ALL, size_t LDALL, double *ALS, size_t LDALS, 
    double *ASL, size_t LDASL, double *ASS, size_t LDASS, bool increment);
    
  template 
  void ComponentScatter(size_t NL, size_t NS, const double *A, size_t LDA, 
    dcomplex *ALL, size_t LDALL, dcomplex *ALS, size_t LDALS, 
    dcomplex *ASL, size_t LDASL, dcomplex *ASS, size_t LDASS, bool increment);

  template 
  void ComponentScatter(size_t NL, size_t NS, const dcomplex *A, size_t LDA, 
    dcomplex *ALL, size_t LDALL, dcomplex *ALS, size_t LDALS, 
    dcomplex *ASL, size_t LDASL, dcomplex *ASS, size_t LDASS, bool increment);

  template 
  void ComponentGather(size_t NL, size_t NS, double *A, size_t LDA, 
    char TransLL, double ScaleLL, const double *ALL, size_t LDALL, 
    char TransLS, double ScaleLS, const double *ALS, size_t LDALS, 
    char TransSL, double ScaleSL, const double *ASL, size_t LDASL, 
    char TransSS, double ScaleSS, const double *ASS, size_t LDASS, 
    bool increment);

  template 
  void ComponentGather(size_t NL, size_t NS, dcomplex *A, size_t LDA, 
    char TransLL, double ScaleLL, const double *ALL, size_t LDALL, 
    char TransLS, double ScaleLS, const double *ALS, size_t LDALS, 
    char TransSL, double ScaleSL, const double *ASL, size_t LDASL, 
    char TransSS, double ScaleSS, const double *ASS, size_t LDASS, 
    bool increment);

  template 
  void ComponentGather(size_t NL, size_t NS, dcomplex *A, size_t LDA, 
    char TransLL, dcomplex ScaleLL, const dcomplex *ALL, size_t LDALL, 
    char TransLS, dcomplex ScaleLS, const dcomplex *ALS, size_t LDALS, 
    char TransSL, dcomplex ScaleSL, const dcomplex *ASL, size_t LDASL, 
    char TransSS, dcomplex ScaleSS, const dcomplex *ASS, size_t LDASS, 
    bool increment);

  template 
  void ComponentGather(size_t NL, size_t NS, double *A, size_t LDA, 
    const double *ALL, size_t LDALL, const double *ALS, size_t LDALS, 
    const double *ASL, size_t LDASL, const double *ASS, size_t LDASS, bool increment);
  
  template 
  void ComponentGather(size_t NL, size_t NS, dcomplex *A, size_t LDA, 
    const double *ALL, size_t LDALL, const double *ALS, size_t LDALS, 
    const double *ASL, size_t LDASL, const double *ASS, size_t LDASS, bool increment);
  
  template 
  void ComponentGather(size_t NL, size_t NS, dcomplex *A, size_t LDA, 
    const dcomplex *ALL, size_t LDALL, const dcomplex *ALS, size_t LDALS, 
    const dcomplex *ASL, size_t LDASL, const dcomplex *ASS, size_t LDASS, bool increment);


  template <typename _F1, typename _F2, typename _FScale>
  void SetMat(char TRANS, size_t M, size_t N, _FScale ALPHA, const _F1 *A, size_t LDA,
    size_t SA, _F2 *B, size_t LDB, size_t SB) {

    assert( TRANS == 'N' or TRANS == 'R' or TRANS == 'T' or TRANS == 'C');

    using namespace Eigen;

    typedef const Matrix<_F1,Dynamic,Dynamic,ColMajor> F1Mat;
    typedef Matrix<_F2,Dynamic,Dynamic,ColMajor> F2Mat;
    typedef Stride<Dynamic,Dynamic> DynamicStride; 

    typedef Map<F1Mat,0,DynamicStride> F1Map;
    typedef Map<F2Mat,0,DynamicStride> F2Map;


    F1Map AMap(A,M,N, DynamicStride(LDA,SA));

    if      ( TRANS == 'N' ) {
      F2Map BMap(B,M,N, DynamicStride(LDB,SB));
      BMap = ALPHA * AMap;
    } else if ( TRANS == 'R' ) {
      F2Map BMap(B,M,N, DynamicStride(LDB,SB));
      BMap = ALPHA * AMap.conjugate();
    } else if ( TRANS == 'T' ) {
      F2Map BMap(B,N,M, DynamicStride(LDB,SB));
      BMap = ALPHA * AMap.transpose();
    } else if ( TRANS == 'C' ) {
      F2Map BMap(B,N,M, DynamicStride(LDB,SB));
      BMap = ALPHA * AMap.adjoint();
    }

  }

#ifdef _CQ_MKL

  template <>
  void SetMat(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, size_t SA, double *B, size_t LDB, size_t SB) {

    if( SA != 1 or SB != 1)
      mkl_domatcopy2('C',TRANS,M,N,ALPHA,A,LDA,SA,B,LDB,SB);
    else
      mkl_domatcopy('C',TRANS,M,N,ALPHA,A,LDA,B,LDB);

  };

  template <>
  void SetMat(char TRANS, size_t M, size_t N, dcomplex ALPHA, const dcomplex *A,
    size_t LDA, size_t SA, dcomplex *B, size_t LDB, size_t SB) {

    if( SA != 1 or SB != 1)
      mkl_zomatcopy2('C',TRANS,M,N,ALPHA,A,LDA,SA,B,LDB,SB);
    else
      mkl_zomatcopy('C',TRANS,M,N,ALPHA,A,LDA,B,LDB);

  };

#else

  template
  void SetMat(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, size_t SA, double *B, size_t LDB, size_t SB);

  template
  void SetMat(char TRANS, size_t M, size_t N, dcomplex ALPHA, const dcomplex *A,
    size_t LDA, size_t SA, dcomplex *B, size_t LDB, size_t SB);


#endif

  template
  void SetMat(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, size_t SA, dcomplex *B, size_t LDB, size_t SB);

  template
  void SetMat(char TRANS, size_t M, size_t N, double ALPHA, const dcomplex *A,
    size_t LDA, size_t SA, dcomplex *B, size_t LDB, size_t SB);

  template
  void SetMat(char TRANS, size_t M, size_t N, dcomplex ALPHA, const double *A,
    size_t LDA, size_t SA, dcomplex *B, size_t LDB, size_t SB);


  /// \brief A2c = [ A  0 ]
  ///              [ 0  A ]
  template <typename _F1, typename _F2>
  void SetMatDiag(size_t M, size_t N, const _F1 *A, size_t LDA, _F2 *A2c, size_t LD2c) {
    SetMat('N',M,N,1.,A,LDA,A2c,LD2c);
    SetMat('N',M,N,1.,A,LDA,A2c + M + N*LD2c,LD2c);
    SetMat('N',M,N,0.,A,LDA,A2c + M,LD2c);
    SetMat('N',M,N,0.,A,LDA,A2c + N*LD2c,LD2c);
  }

  template void SetMatDiag(size_t, size_t, const double*, size_t, double*, size_t);
  template void SetMatDiag(size_t, size_t, const double*, size_t, dcomplex*, size_t);
  template void SetMatDiag(size_t, size_t, const dcomplex*, size_t, dcomplex*, size_t);




  template<>
  void SetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, dcomplex *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,A,LDA,1,reinterpret_cast<double*>(B),2*LDB,2);

  }; // SetMatRE (complex)

  template<>
  void SetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, double *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,A,LDA,1,B,LDB,1);
    

  }; // SetMatRE (real)


  template<>
  void SetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, dcomplex *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,A,LDA,1,reinterpret_cast<double*>(B)+1,2*LDB,2);

  }; // SetMatIM (complex)

  template<>
  void SetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, double *B, size_t LDB) {

    assert(false);

  }; // SetMatRM (real)

  template<>
  void GetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const dcomplex *A,
    size_t LDA, double *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,reinterpret_cast<const double*>(A),2*LDA,2,B,LDA,1);

  }; // GetMatRE (complex)

  template<>
  void GetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, double *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,A,LDA,1,B,LDA,1);

  }; // GetMatRE (real)

  template<>
  void GetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const dcomplex *A,
    size_t LDA, double *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,reinterpret_cast<const double*>(A)+1,2*LDA,2,B,LDA,1);

  }; // GetMatIM (complex)

  template<>
  void GetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, double *B, size_t LDB) {

    SetMat(TRANS,M,N,0.,A,LDA,1,B,LDA,1);

  }; // GetMatIM (real)









  // Non-contiguous sub matrix operations

  template <typename _F1, typename _F2>
  void SubMatSet(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut) {

    
    Eigen::Map<
      Eigen::Matrix<_F1,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ABigMap(ABig,LDAB,N);

    Eigen::Map<
      Eigen::Matrix<_F2,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ASmallMap(ASmall,LDAS,NSub);

    size_t i(0);
    for( auto& iCut : SubMatCut ) {
      size_t deltaI = iCut.second - iCut.first;
      size_t j(0);
    for( auto& jCut : SubMatCut ) {
      size_t deltaJ = jCut.second - jCut.first;
    
      ASmallMap.block(i,j,deltaI,deltaJ).noalias() =
        ABigMap.block(iCut.first,jCut.first,deltaI,deltaJ);
    
      j += deltaJ;
    }
      i += deltaI;
    }
  };

  template <typename _F1, typename _F2>
  void SubMatGet(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut) {

    
    Eigen::Map<
      Eigen::Matrix<_F1,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ABigMap(ABig,LDAB,N);

    Eigen::Map<
      Eigen::Matrix<_F2,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ASmallMap(ASmall,LDAS,NSub);

    size_t i(0);
    for( auto& iCut : SubMatCut ) {
      size_t deltaI = iCut.second - iCut.first;
      size_t j(0);
    for( auto& jCut : SubMatCut ) {
      size_t deltaJ = jCut.second - jCut.first;
    
      ABigMap.block(iCut.first,jCut.first,deltaI,deltaJ).noalias() =
        ASmallMap.block(i,j,deltaI,deltaJ);
    
      j += deltaJ;
    }
      i += deltaI;
    }
  };

  template <typename _F1, typename _F2>
  void SubMatInc(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut) {

    
    Eigen::Map<
      Eigen::Matrix<_F1,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ABigMap(ABig,LDAB,N);

    Eigen::Map<
      Eigen::Matrix<_F2,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ASmallMap(ASmall,LDAS,NSub);

    size_t i(0);
    for( auto& iCut : SubMatCut ) {
      size_t deltaI = iCut.second - iCut.first;
      size_t j(0);
    for( auto& jCut : SubMatCut ) {
      size_t deltaJ = jCut.second - jCut.first;
    
      ASmallMap.block(i,j,deltaI,deltaJ).noalias() +=
        ABigMap.block(iCut.first,jCut.first,deltaI,deltaJ);
    
      j += deltaJ;
    }
      i += deltaI;
    }
  };

  template <typename _F1, typename _F2>
  void IncBySubMat(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut) {

    
    Eigen::Map<
      Eigen::Matrix<_F1,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ABigMap(ABig,LDAB,N);

    Eigen::Map<
      Eigen::Matrix<_F2,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
        ASmallMap(ASmall,LDAS,NSub);

    size_t i(0);
    for( auto& iCut : SubMatCut ) {
      size_t deltaI = iCut.second - iCut.first;
      size_t j(0);
    for( auto& jCut : SubMatCut ) {
      size_t deltaJ = jCut.second - jCut.first;
    
      ABigMap.block(iCut.first,jCut.first,deltaI,deltaJ).noalias() +=
        ASmallMap.block(i,j,deltaI,deltaJ);
    
      j += deltaJ;
    }
      i += deltaI;
    }
  };

  // Instantiate functions

  template 
  void SubMatSet(size_t M, size_t N, size_t MSub, size_t NSub, double *ABig, 
    size_t LDAB, double *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  template 
  void SubMatSet(size_t M, size_t N, size_t MSub, size_t NSub, dcomplex *ABig, 
    size_t LDAB, dcomplex *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  template 
  void SubMatGet(size_t M, size_t N, size_t MSub, size_t NSub, double *ABig, 
    size_t LDAB, double *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  template 
  void SubMatInc(size_t M, size_t N, size_t MSub, size_t NSub, double *ABig, 
    size_t LDAB, double *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  template 
  void IncBySubMat(size_t M, size_t N, size_t MSub, size_t NSub, double *ABig, 
    size_t LDAB, double *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  template 
  void IncBySubMat(size_t M, size_t N, size_t MSub, size_t NSub, dcomplex *ABig, 
    size_t LDAB, dcomplex *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);


  template <typename Apha, typename ATyp, typename BTyp, typename CTyp>
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, Apha alpha,
    ATyp* A, size_t LDA, std::vector<BTyp*> V, size_t LDB, CTyp* SCR,
    std::vector<CTyp*> U, size_t LDC) {
 
    size_t X;
    if( KB % KA == 0 )
      X = KB / KA;
    else
      CErr("KB is not a multiple of KA; No implicit identity tensor possible");
 
    if( V.size() != U.size() )
      CErr("Different number of input and output matrices in TransformLeft");
 
    for(auto iMat = 0; iMat < V.size(); iMat++) {
 
      BTyp* B = V[iMat];
      CTyp* C = U[iMat];
      bool inPlace = ((void*)B != (void*)C);
      CTyp* GemmOut = inPlace ? C : SCR;
      size_t GemmLD = inPlace ? LDC : M;

      for(auto iBlock = 0; iBlock < X; iBlock++) {
 
        blas::gemm(
          blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
          M, N, KA,
          alpha,
          A, LDA,
          B + iBlock*KA, LDB,
          CTyp(0.),
          GemmOut, GemmLD);

        if( !inPlace )
          SetMat('N', M, N, CTyp(1.), GemmOut, M, C + iBlock*KA, LDC);
        else
          GemmOut += KA;
      }
 
    }
  }
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, double alpha,
    double* A, size_t LDA, std::vector<double*> V, size_t LDB, double* SCR,
    std::vector<double*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, dcomplex alpha,
    double* A, size_t LDA, std::vector<double*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, double alpha,
    dcomplex* A, size_t LDA, std::vector<double*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, double alpha,
    double* A, size_t LDA, std::vector<dcomplex*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, double alpha,
    dcomplex* A, size_t LDA, std::vector<dcomplex*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, dcomplex alpha,
    double* A, size_t LDA, std::vector<dcomplex*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, dcomplex alpha,
    dcomplex* A, size_t LDA, std::vector<double*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);
 
  template
  void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, dcomplex alpha,
    dcomplex* A, size_t LDA, std::vector<dcomplex*> V, size_t LDB, dcomplex* SCR,
    std::vector<dcomplex*> U, size_t LDC);


  template <typename ATyp, typename TTyp, typename BTyp>
  void PairTransformation(char TRANST, const TTyp * T1, size_t LDT1, size_t OffK, 
    const TTyp * T2, size_t LDT2, size_t OffL,
    char TRANSA, const ATyp * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, BTyp * B, size_t NK, size_t NL, ATyp * ASCR, BTyp * SCR, bool increment) {

    blas::Op OP_TRANS;
    if (TRANST == 'T') {
      OP_TRANS = blas::Op::Trans;
    } else if (TRANST == 'C') {
      OP_TRANS = blas::Op::ConjTrans;
    } else if (TRANST == 'N') {
      OP_TRANS = blas::Op::NoTrans;
      OffK *= LDT1;
      OffL *= LDT2;
    } else CErr("Wrong TRANST in PairTransformation");
    
    const ATyp * AP = nullptr; 
    
    // Pre-processing A 
    if(TRANSA == 'N') AP = A;
    else if(TRANSA == 'R') {
      SetMat('R', NI*NJ, NM, ATyp(1.), A, NI*NJ, ASCR, NI*NJ);
      AP = ASCR; 
    } else if (TRANSA == 'C' or TRANSA == 'T') {
      SetMat(TRANSA, NM, NI*NJ, ATyp(1.), A, NM, ASCR, NI*NJ);
      AP = ASCR; 
    } else CErr("Wrong TRANSA in PairTransformation");

    // first transformation
    // SCR(J M, K) = A(I, J M)^H @ T(I, K)
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, OP_TRANS, 
      NJ*NM, NK, NI, ATyp(1.), AP, NI, T1+OffK, LDT1, BTyp(0.), SCR, NJ*NM);

    // second transformation
    BTyp outFactor = increment ? 1.0 : 0.0;
    // B(M K, L) = SCR(J, M K)^H @ T(J, L) 
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, OP_TRANS, 
      NM*NK, NL, NJ, BTyp(1.), SCR, NJ, T2+OffL, LDT2, outFactor, B, NM*NK);
    
    // Post-processing B
    if (TRANSB == 'N')      IMatCopy('T', NM, NK*NL, BTyp(1.), B, NM, NK*NL); 
    else if (TRANSB == 'R') IMatCopy('C', NM, NK*NL, BTyp(1.), B, NM, NK*NL); 
    else if (TRANSB == 'C') IMatCopy('R', NM, NK*NL, BTyp(1.), B, NM, NM); 
    else if (TRANSB == 'T') ; // do nothing
    else CErr("Wrong TRANSB in PairTransformation");
  
  } // PairTransformation

  template 
  void PairTransformation(char TRANST, const double * T1, size_t LDT1, size_t OffK, 
    const double * T2, size_t LDT2, size_t OffL,
    char TRANSA, const double * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, double * B, size_t NK, size_t NL, double * ASCR, double * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const double * T1, size_t LDT1, size_t OffK, 
    const double * T2, size_t LDT2, size_t OffL,
    char TRANSA, const dcomplex * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, dcomplex * ASCR, dcomplex * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const dcomplex * T1, size_t LDT1, size_t OffK, 
     const dcomplex * T2, size_t LDT2, size_t OffL,
    char TRANSA, const double * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, double * ASCR, dcomplex * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const dcomplex * T1, size_t LDT1, size_t OffK, 
    const dcomplex * T2, size_t LDT2, size_t OffL,
    char TRANSA, const dcomplex * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, dcomplex * ASCR, dcomplex * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const double * T, size_t LDT, size_t OffK, size_t OffL,
    char TRANSA, const double * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, double * B, size_t NK, size_t NL, double * ASCR, double * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const double * T, size_t LDT, size_t OffK, size_t OffL,
    char TRANSA, const dcomplex * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, dcomplex * ASCR, dcomplex * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const dcomplex * T, size_t LDT, size_t OffK, size_t OffL,
    char TRANSA, const double * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, double * ASCR, dcomplex * SCR, bool increment);

  template 
  void PairTransformation(char TRANST, const dcomplex * T, size_t LDT, size_t OffK, size_t OffL,
    char TRANSA, const dcomplex * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, dcomplex * B, size_t NK, size_t NL, dcomplex * ASCR, dcomplex * SCR, bool increment);

}; // namespace ChronusQ


