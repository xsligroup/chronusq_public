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

#include <chronusq_sys.hpp>
#include <cxxapi/output.hpp>
#include <Eigen/Core>

namespace ChronusQ {

// Smallest value to print
constexpr long double PRINT_SMALL = 1e-10;

template <typename T, 
          typename std::enable_if<
                     std::is_same<T,double>::value, int>::type = 0>
void printValWithCheck(std::ostream& out, T val, const size_t printWidth) {
  if(std::abs(val) > PRINT_SMALL) out << std::setw(printWidth) << val; 
  else if(std::isnan(val))        out << std::setw(printWidth) << "NAN";
  else if(std::isinf(val))        out << std::setw(printWidth) << "INF";
  else                            out << std::setw(printWidth) << 0.;
} // check for double

template <typename T, 
          typename std::enable_if<
                     std::is_integral<T>::value, int>::type = 0>
void printValWithCheck(std::ostream& out, T val, const size_t printWidth) {
  out << std::setw(printWidth) << +val; 
} // No check for int types

/**
 *  \brief Base routine to print out a matrix given raw storage and dimension
 *  parameters in a standard (pretty) way.
 *
 *  \warning T should be a "real" type (float, double, etc)
 *
 *  \param [in/out] out        Output device
 *  \param [in]     A          Raw storage of the matrix A
 *  \param [in]     M          Number of rows of the matrix A
 *  \param [in]     N          Number of cols of the matrix A
 *  \param [in]     LDA        Leading dimension of the matrix A
 *  \param [in]     colStride  Column stride of A
 *  \param [in]     list       Number of columns to print before page break
 *  \param [in]     printWidth Field with of a matrix column
 */
template <typename T>
void prettyPrintSmartBase(std::ostream& out, const T* A, const size_t M,
  const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  int end,endLT;
  out << bannerTop;

  out << std::scientific << std::left << std::setprecision(8);
  for(size_t i = 0; i < N; i += list) {
    out << std::endl;
    end = list;
    out << std::setw(5) << " ";
    if((i + list) >= N) end = N - i;
    out << std::right;
    for(size_t k = i; k < i+end; k++) out << std::setw(printWidth) << k+1;
    out << std::endl;
    for(size_t j = 0; j < M; j++) {
      out << std::setw(5) << std::left << j+1;
      out << std::right;
      for(size_t n = i; n < i+end; n++) {
        printValWithCheck(out, A[j*colStride + n*LDA], printWidth);
      }
      out << std::endl;
    };
  };
  out << bannerEnd << std::endl;

}; // prettyPrintSmartBase

/**
 *  \brief Print a real matrix in a standard way given raw storage and
 *  dimension parameters
 *
 *  \param [in/out] out        Output device
 *  \param [in]     A          Raw storage of the matrix A
 *  \param [in]     M          Number of rows of the matrix A
 *  \param [in]     N          Number of cols of the matrix A
 *  \param [in]     LDA        Leading dimension of the matrix A
 *  \param [in]     colStride  Column stride of A
 *  \param [in]     list       Number of columns to print before page break
 *  \param [in]     printWidth Field with of a matrix column
 */
template <typename T, 
          typename std::enable_if<
                     std::is_same<T,double>::value,int>::type = 0>
void prettyPrintSmart(std::ostream& out, std::string str, const T* A,
  const size_t M, const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  out.precision(10);
  out.fill(' ');
  out << std::endl << str + ": " << std::endl;

  prettyPrintSmartBase(out,A,M,N,LDA,colStride,list,printWidth);

}; // prettyPrintSmart (T = double)

/**
 *  \brief Print a complex matrix in a standard way given raw storage and
 *  dimension parameters.
 *
 *  Prints real and imaginary parts as separate matricies
 *
 *  \param [in/out] out        Output device
 *  \param [in]     A          Raw storage of the matrix A
 *  \param [in]     M          Number of rows of the matrix A
 *  \param [in]     N          Number of cols of the matrix A
 *  \param [in]     LDA        Leading dimension of the matrix A
 *  \param [in]     colStride  Column stride of A
 *  \param [in]     list       Number of columns to print before page break
 *  \param [in]     printWidth Field with of a matrix column
 */
template <typename T, 
          typename std::enable_if<
                     std::is_same<T,dcomplex>::value,int>::type = 0>
void prettyPrintSmart(std::ostream& out, std::string str, const T* A,
  const size_t M, const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 5, const size_t printWidth = 16) {

  out.precision(10);
  out.fill(' ');

  out << std::endl << "Re[" << str + "]: " << std::endl;
  prettyPrintSmartBase(out,reinterpret_cast<const double*>(A),
                       M,N,2*LDA,2*colStride,list,printWidth);

  out << std::endl << "Im[" << str + "]: " << std::endl;
  prettyPrintSmartBase(out,reinterpret_cast<const double*>(A)+1,
                       M,N,2*LDA,2*colStride,list,printWidth);

}; // prettyPrintSmart (T = dcomplex)

/**
 *  \brief Print a int matrix in a standard way given raw storage and
 *  dimension parameters
 *
 *  \param [in/out] out        Output device
 *  \param [in]     A          Raw storage of the matrix A
 *  \param [in]     M          Number of rows of the matrix A
 *  \param [in]     N          Number of cols of the matrix A
 *  \param [in]     LDA        Leading dimension of the matrix A
 *  \param [in]     colStride  Column stride of A
 *  \param [in]     list       Number of columns to print before page break
 *  \param [in]     printWidth Field with of a matrix column
 */
template <typename T, 
          typename std::enable_if<std::is_integral<T>::value,int>::type = 0>
void prettyPrintSmart(std::ostream& out, std::string str, const T* A,
  const size_t M, const size_t N, const size_t LDA, const size_t colStride = 1, 
  const size_t list = 10, const size_t printWidth = 8) {

  out.fill(' ');
  out << std::endl << str + ": " << std::endl;

  prettyPrintSmartBase(out,A,M,N,LDA,colStride,list,printWidth);

}; // prettyPrintSmart (T = int)

template <typename T, 
          typename std::enable_if<
                     std::is_same<T,double>::value,int>::type = 0>
void mathematicaPrint(std::ostream& out, std::string str, const T* A,
  const size_t M, const size_t N, const size_t LDA, 
  const size_t colStride = 1) { 

  out << str << ":   \n";

  out << "{ \n";
  for(auto i = 0; i < M; i++){
    if( N != 1) out << "{  ";
  for(auto j = 0; j < N; j++){
    out << std::setprecision(12) << A[colStride*i + j*LDA];
    if(j != N-1) out << ", ";
  }
    if(N != 1) out << "}";
    if( i != M-1 ) out << ",";
  }
  out << "\n}\n";
  

};

template <typename T, 
          typename std::enable_if<
                     std::is_same<T,dcomplex>::value,int>::type = 0>
void mathematicaPrint(std::ostream& out, std::string str, const T* A,
  const size_t M, const size_t N, const size_t LDA,
  const size_t colStride = 1) { 

  out << str << ":   \n";

  out << "{ \n";
  for(auto i = 0; i < M; i++){
    if( N != 1 ) out << "{  ";
  for(auto j = 0; j < N; j++){
    out << "Complex[";
    out << std::setprecision(12) << std::real(A[colStride*i + j*LDA]);
    out << ",";
    out << std::setprecision(12) << std::imag(A[colStride*i + j*LDA]);
    out << "]";
    if(j != N-1) out << ", ";
  }
    if( N != 1 ) out << "}";
    if( i != M-1 ) out << ",";
  }
  out << "\n}\n";
  

};

}; // namespace ChronusQ

