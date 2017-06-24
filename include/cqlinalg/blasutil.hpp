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

#include <cqlinalg/cqlinalg_config.hpp>

namespace ChronusQ {



  // Spin Scatter a matrix
  template <typename _F1, typename _F2>
  void SpinScatter(size_t N, const _F1 *A, size_t LDA, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA = false, bool BBeqAA = false);

  // Spin Scatter a matrix
  template <typename _F1, typename _F2>
  void SpinScatter(size_t M, size_t N, const _F1 *A, size_t LDA, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA = false, bool BBeqAA = false);

  // Spin Scatter matrix blocks
  template <typename _F1, typename _F2>
  void SpinScatter(size_t M, size_t N, const _F1 *AA, size_t LDAA, const _F1 *AB, size_t LDAB,
      const _F1 *BA, size_t LDBA, const _F1 *BB, size_t LDBB, _F2 *AS, size_t LDAS,
      _F2 *AZ, size_t LDAZ, _F2 *AY, size_t LDAY, _F2 *AX, size_t LDAX,
      bool zeroABBA = false, bool BBeqAA = false);

  // Spin Gather a matrix
  template <typename _F1, typename _F2>
  void SpinGather(size_t N, _F1 *A, size_t LDA, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY = false, bool zeroZ = false);

  // Spin Gather a matrix
  template <typename _F1, typename _F2>
  void SpinGather(size_t M, size_t N, _F1 *A, size_t LDA, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY = false, bool zeroZ = false);

  // Spin Gather a matrix
  template <typename _F1, typename _F2>
  void SpinGather(size_t M, size_t N, _F1 *AA, size_t LDAA, _F1 *AB, size_t LDAB,
      _F1 *BA, size_t LDBA, _F1 *BB, size_t LDBB, const _F2 *AS, size_t LDAS,
      const _F2 *AZ, size_t LDAZ, const _F2 *AY, size_t LDAY, const _F2 *AX, size_t LDAX,
      bool zeroXY = false, bool zeroZ = false);

  // Component Scatter a matrix 
  template <typename _F1, typename _F2>
  void ComponentScatter(size_t NL, size_t NS, const _F1 *A, size_t LDA, 
      _F2 ScaleLL, _F2 *ALL, size_t LDALL, _F2 ScaleLS, _F2 *ALS, size_t LDALS, 
      _F2 ScaleSL, _F2 *ASL, size_t LDASL, _F2 ScaleSS, _F2 *ASS, size_t LDASS, 
      bool increment = false);
  
  // Component Scatter a matrix 
  template <typename _F1, typename _F2>
  void ComponentScatter(size_t NL, size_t NS, const _F1 *A, size_t LDA, _F2 *ALL, size_t LDALL,
      _F2 *ALS, size_t LDALS, _F2 *ASL, size_t LDASL, _F2 *ASS, size_t LDASS, bool increment = false);

  // Component Gatter a matrix 
  template <typename _F1, typename _F2>
  void ComponentGather(size_t NL, size_t NS, _F1 *A, size_t LDA, 
      char TransLL, _F2 ScaleLL, const _F2 *ALL, size_t LDALL, char TransLS, const _F2 ScaleLS, _F2 *ALS, size_t LDALS, 
      char TransSL, _F2 ScaleSL, const _F2 *ASL, size_t LDASL, char TransSS, const _F2 ScaleSS, _F2 *ASS, size_t LDASS, 
      bool increment = false);

  // Component Gatter a matrix 
  template <typename _F1, typename _F2>
  void ComponentGather(size_t NL, size_t NS, _F1 *A, size_t LDA, const _F2 *ALL, size_t LDALL,
      const _F2 *ALS, size_t LDALS, const _F2 *ASL, size_t LDASL, const _F2 *ASS, size_t LDASS, bool increment = false);


  // B = ALPHA * OP(A)
  template <typename _F1, typename _F2, typename _FScale>
  void SetMat(char TRANS, size_t M, size_t N, _FScale ALPHA, const _F1 *A, size_t LDA,
    size_t SA, _F2 *B, size_t LDB, size_t SB);

  template <typename _F1, typename _F2, typename _FScale>
  void SetMat(char TRANS, size_t M, size_t N, _FScale ALPHA, const _F1 *A, size_t LDA,
    _F2 *B, size_t LDB) {

    SetMat(TRANS,M,N,ALPHA,A,LDA,1,B,LDB,1);

  }


  // RE(B) = ALPHA * OP(A)
  template <typename _F>
  void SetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, _F *B, size_t LDB);

  // IM(B) = ALPHA * OP(A)
  template <typename _F>
  void SetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const double *A,
    size_t LDA, _F *B, size_t LDB);

  // B = ALPHA * OP(RE(A))
  template <typename _F>
  void GetMatRE(char TRANS, size_t M, size_t N, double ALPHA, const _F *A,
    size_t LDA, double *B, size_t LDB);

  // B = ALPHA * OP(IM(A))
  template <typename _F>
  void GetMatIM(char TRANS, size_t M, size_t N, double ALPHA, const _F *A,
    size_t LDA, double *B, size_t LDB);

  // A2c = [ A  0 ]
  //       [ 0  A ]
  template <typename _F1, typename _F2>
  void SetMatDiag(size_t M, size_t N, const _F1 *A, size_t LDA, _F2 *A2c, size_t LD2c);
















  // Sub-matrix utilities for multiple, non-contiguous blocks


  /**
   *  \brief Set a packed submatrix to the corresponding blocks of
   *  the full matrix using sub matrix cuts
   *
   *  ASmall = ABig submatrix
   *
   *  \param [in] M       Number of rows in ABig
   *  \param [in] N       Number of columns in ABig
   *  \param [in] MSub    Number of rows in ASmall
   *  \param [in] NSub    Number of columns in ASmall
   *  \param [in] ABig    Raw storage of supermatrix
   *  \param [in] LDAB    Leading dimension of ABig
   *  \param [in/out] ASmall  Raw storage of submatrix
   *  \param [in] LDAS    Leading dimension of ASmall
   *  \param [in] SubMatCut Vector of pairs specifing the blocks of
   *                        the super matrix to be used
   *  
   */ 
  template <typename _F1, typename _F2>
  void SubMatSet(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);




  /**
   *  \brief Pack subblock of a matrix into a contiguous array
   *  using sub matrix cuts
   *
   *  ABig submatrix = ASmall
   *
   *  \param [in] M       Number of rows in ABig
   *  \param [in] N       Number of columns in ABig
   *  \param [in] MSub    Number of rows in ASmall
   *  \param [in] NSub    Number of columns in ASmall
   *  \param [in/out] ABig    Raw storage of supermatrix
   *  \param [in] LDAB    Leading dimension of ABig
   *  \param [in] ASmall  Raw storage of submatrix
   *  \param [in] LDAS    Leading dimension of ASmall
   *  \param [in] SubMatCut Vector of pairs specifing the blocks of
   *                        the super matrix to be used
   *  
   */ 
  template <typename _F1, typename _F2>
  void SubMatGet(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);




  /**
   *  \brief Increment a packed submatrix to the corresponding blocks of
   *  the full matrix using sub matrix cuts
   *
   *  ASmall += ABig submatrix
   *
   *  \param [in] M       Number of rows in ABig
   *  \param [in] N       Number of columns in ABig
   *  \param [in] MSub    Number of rows in ASmall
   *  \param [in] NSub    Number of columns in ASmall
   *  \param [in] ABig    Raw storage of supermatrix
   *  \param [in] LDAB    Leading dimension of ABig
   *  \param [in/out] ASmall  Raw storage of submatrix
   *  \param [in] LDAS    Leading dimension of ASmall
   *  \param [in] SubMatCut Vector of pairs specifing the blocks of
   *                        the super matrix to be used
   *  
   */ 
  template <typename _F1, typename _F2>
  void SubMatInc(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);




  /**
   *  \brief Increment a subblock of a matrix by a packed array
   *  using sub matrix cuts
   *
   *  ABig submatrix += ASmall
   *
   *  \param [in] M       Number of rows in ABig
   *  \param [in] N       Number of columns in ABig
   *  \param [in] MSub    Number of rows in ASmall
   *  \param [in] NSub    Number of columns in ASmall
   *  \param [in/out] ABig    Raw storage of supermatrix
   *  \param [in] LDAB    Leading dimension of ABig
   *  \param [in] ASmall  Raw storage of submatrix
   *  \param [in] LDAS    Leading dimension of ASmall
   *  \param [in] SubMatCut Vector of pairs specifing the blocks of
   *                        the super matrix to be used
   *  
   */ 
  template <typename _F1, typename _F2>
  void IncBySubMat(size_t M, size_t N, size_t MSub, size_t NSub, _F1 *ABig, 
    size_t LDAB, _F2 *ASmall, size_t LDAS, 
    std::vector<std::pair<size_t,size_t>> &SubMatCut);

  /**
   * \brief Multiply a series of matrices on the right with a given matrix
   *
   * C_i = \alpha AB_i for all B_i in V
   *
   * If the matrices are of different dimension, the A matrix is assumed to
   * be tensored with the appropriate identity matrix to give the dimension of
   * the A matrix.
   *
   * e.g.
   *
   * B = | a b |  A = | e |  -->  C = | ea eb |
   *     | c d |                      | ec ed |
   *
   * \param[in] M      Number of rows in A
   * \param[in] N      Numer of columns in B
   * \param[in] KA     Number of columns in A
   * \param[in] KB     Number of rows in B. Must be a multiple of KA.
   * \param[in] alpha  Scalar in multiplication
   * \param[in] A      Pointer to A
   * \param[in] LDA    Leading dimension of A
   * \param[in] V      Vector of pointers to B matrices
   * \param[in] LDB    Leading dimension of B
   * \param[in] SCR    Scratch space, at least M*N. Only referenced if any B_i
   *                   is the same as any C_i
   * \param[out] U     Vector of pointers to C matrices. Can be the same as V.
   * \param[out] LDC   Leading dimension of C
   */
   template <typename Apha, typename ATyp, typename BTyp, typename CTyp>
   void TransformLeft(size_t M, size_t N, size_t KA, size_t KB, Apha alpha,
     ATyp* A, size_t LDA, std::vector<BTyp*> V, size_t LDB, CTyp* SCR,
     std::vector<CTyp*> U, size_t LDC);

  /**
   * \brief Transform tensor A with transformation matrix T1 and T2 in pair 
   *        (first two dim of A)
   *
   * B(K, L, M) = T1(I, K+OffK)^H A(I, J, M) @ T2(J, L+OffL) 
   *
   * \param [in]  TRANSA    Whether transpose/adjoint A as A(IJ, M) before transformation
   * \param [in]  DI        Dimension of index I in A
   * \param [in]  DJ        Dimension of index J in A
   * \param [in]  DM        Dimension of index M in A
   * \param [in]  A         Tensor before transformation
   *
   * \param [in]  TRANST    Whether transpose/adjoint T1 and T2 
   * \param [in]  T1        Transformation matrix T1
   * \param [in]  LDT1      Leading dimension of T1
   * \param [in]  OffK      offset in index K
   * \param [in]  T2        Transformation matrix T2
   * \param [in]  LDT2      Leading dimension of T2
   * \param [in]  OffL      offset in index L 
   *
   * \param [in]  TRANSB    Whether transpose/adjoint B as B(KL, M) after transformation
   * \param [in]  DK        Dimension of index K in B
   * \param [in]  DL        Dimension of index L in B
   * \param [out] B         Tensor after transformation, can be same as A
   *
   * \param[in]   ASCR      Scratch space for A at least NI*NJ*NM of ATyp. 
   *                        Only reference when TRANSA != 'N'
   * \param[in]   SCR       Scratch space for intermediates, at least NJ*NM*NK of BTyp. 
   *
   * \param [in]  increment  Perform B += result if true
   *
   */

  template <typename ATyp, typename TTyp, typename BTyp>
  void PairTransformation(char TRANST, const TTyp * T1, size_t LDT1, size_t OffK, 
    const TTyp * T2, size_t LDT2, size_t OffL,
    char TRANSA, const ATyp * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, BTyp * B, size_t NK, size_t NL, ATyp * ASCR, BTyp * SCR, bool increment);
  
  template <typename ATyp, typename TTyp, typename BTyp>
  void PairTransformation(char TRANST, const TTyp * T, size_t LDT, size_t OffK, size_t OffL,
    char TRANSA, const ATyp * A, size_t NI, size_t NJ, size_t NM, 
    char TRANSB, BTyp * B, size_t NK, size_t NL, ATyp * ASCR, BTyp * SCR, bool increment) {
      
    PairTransformation(TRANST, T, LDT, OffK, T, LDT, OffL, TRANSA, 
      A, NI, NJ, NM, TRANSB, B, NK, NL, ASCR, SCR, increment);
    
  };

}; // namespace ChronusQ


