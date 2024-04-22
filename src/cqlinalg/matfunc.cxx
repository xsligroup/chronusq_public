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
#include <cqlinalg/matfunc.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/eig.hpp>

namespace ChronusQ {


  template <typename F, typename _F1, typename _F2>
  void MatDiagFunc(const F &func, size_t N, _F1 *A, size_t LDA, _F2 *B, 
    size_t LDB) {

    // Allocate space for eigenvalues and scratch
    double *W = CQMemManager::get().malloc<double>(N);
    _F1* SCR  = CQMemManager::get().malloc<_F1>(N*N);
    _F2* SCR2  = CQMemManager::get().malloc<_F2>(N*N);

    std::copy_n(A,N*N,SCR); // Copy A to SCR

    // A = V * a * V**H
    HermetianEigen('V','U',N,SCR,N,W);

    // Compute X = V * func(a)
    for(auto j = 0; j < N; j++)
    for(auto i = 0; i < N; i++)
      SCR2[i + j*N] = SCR[i + j*N] * func(W[j]);
    
    // Compute B**H = V * X**H
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,N,N,N,_F2(1.),SCR,N,SCR2,N,_F2(0.),B,LDB);

    // FIXME: Use MKL for this transpose when direct is merged in
    Eigen::Map<Eigen::Matrix<_F2,Eigen::Dynamic,Eigen::Dynamic,
      Eigen::ColMajor>> BMap(B,LDB,N);

    BMap.adjointInPlace(); 

    CQMemManager::get().free(SCR,SCR2,W);

  };

  template void MatDiagFunc(const std::function<double(double)> &func, size_t N, dcomplex *A, size_t LDA,
                            dcomplex *B, size_t LDB);
  template void MatDiagFunc(const std::function<double(double)> &func, size_t N, double *A, size_t LDA,
                            double *B, size_t LDB);


  template <typename _FExp, typename _F1, typename _F2>
  void MatExp(char ALG, size_t N, _FExp ALPHA, _F1 *A, size_t LDA, 
    _F2 *ExpA, size_t LDEXPA) {

    assert(ALG == 'D');
    assert(std::real(ALPHA) < 1e-14);

    double AIM = std::is_same<_FExp,dcomplex>::value ? std::imag(ALPHA) : 0.;

    MatDiagFunc([&](double x) -> _F2 { 
        return dcomplex(std::cos(AIM*x),std::sin(AIM*x)); 
      }, N,A,LDA,ExpA,LDEXPA);

  };


//template void MatExp(char,size_t,double,double*,size_t,double*,size_t);

  template void MatExp(char,size_t,dcomplex,dcomplex*,size_t,dcomplex*,size_t);


  template <typename MatsU>
  void MatExp(size_t N, MatsU *A, size_t LDA,
    MatsU *ExpA, size_t LDEXPA) {

    // allocate memory
    size_t N2 = N*N;
    MatsU * OddTerms  = CQMemManager::get().template malloc<MatsU>(N2);
    MatsU * EvenTerms = CQMemManager::get().template malloc<MatsU>(N2);

    std::fill_n(ExpA, N2, MatsU(0.));
    // form zero order term
    for (auto i = 0ul; i < N; i++) ExpA[i + i*N] = MatsU(1.);

    // form 1st order term
    std::copy_n(A, N2, OddTerms);

    double residue;
    double small_number = std::numeric_limits<double>::epsilon();
    bool converged = false;
    size_t maxIter = 200;
    for (auto iter = 0ul; iter < maxIter; iter+=2) {
      // add the odd term
      MatAdd('N', 'N', N, N, MatsU(1.), OddTerms, N, MatsU(1.), ExpA, N, ExpA, N);

      residue = lapack::lange(lapack::Norm::Fro, N, N, OddTerms, N);
      if (residue <= small_number) {
        converged = true;
        break;
      }

      // form and add next even term
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        N, N, N, MatsU(1./(iter+2)), A, N, OddTerms, N,
        MatsU(0.), EvenTerms, N);
      MatAdd('N', 'N', N, N, MatsU(1.), EvenTerms, N, MatsU(1.), ExpA, N, ExpA, N);

      residue = lapack::lange(lapack::Norm::Fro, N, N, EvenTerms, N);
      if (residue <= small_number) {
        converged = true;
        break;
      }

      // form next odd term
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        N, N, N, MatsU(1./(iter+3)), A, N, EvenTerms, N,
        MatsU(0.), OddTerms, N);

    }

    CQMemManager::get().free(OddTerms, EvenTerms);

    if(not converged) throw std::runtime_error("Matrix Exponential calculation failed to converge");
  }; // MCSCF::MatExpT

  template void MatExp(size_t N, double *A, size_t LDA, double *ExpA, size_t LDEXPA);
  template void MatExp(size_t N, dcomplex *A, size_t LDA, dcomplex *ExpA, size_t LDEXPA);

}; // namespace ChronusQ
