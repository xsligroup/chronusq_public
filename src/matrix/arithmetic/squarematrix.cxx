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

#include <matrix.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  SquareMatrix<MatsT>::SquareMatrix( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ):
      SquareMatrix(scaled.matrix().memManager(), scaled.matrix().dimension()) {
    if (scaled.isPauli() and dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(scaled.matrix()).hasZ())
      CErr("Cannot create a SquareMatrix from a PauliSpinorSquareMatrices with XYZ components.");
    SetMat('N',dimension(),dimension(),scaled.scalar(),scaled.matrix().pointer(),dimension(),pointer(),dimension());
  }
  template SquareMatrix<double>::SquareMatrix( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>::SquareMatrix( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>::SquareMatrix( const ScaledSquareMatrix<dcomplex, double>& );
  template SquareMatrix<dcomplex>::SquareMatrix( const ScaledSquareMatrix<double, dcomplex>& );
  template SquareMatrix<dcomplex>::SquareMatrix( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  SquareMatrix<MatsT>&
  SquareMatrix<MatsT>::operator=( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) {
    if (dimension() != scaled.matrix().dimension())
      CErr("Cannot assign SquareMatrix of different size.");
    if (std::is_same<MatsT, double>::value and
        std::is_same<MatsU, dcomplex>::value)
      CErr("Cannot assign a complex SquareMatrix object to a real one.");
    if (scaled.isPauli() and dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(scaled.matrix()).hasZ())
      CErr("Cannot assign a SquareMatrix from a PauliSpinorSquareMatrices with XYZ components.");
    SetMat('N',dimension(),dimension(),scaled.scalar(),scaled.matrix().pointer(),dimension(),pointer(),dimension());
    return *this;
  }
  template SquareMatrix<double>&
  SquareMatrix<double>::operator=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>&
  SquareMatrix<dcomplex>::operator=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>&
  SquareMatrix<dcomplex>::operator=( const ScaledSquareMatrix<dcomplex, double>& );
  template SquareMatrix<dcomplex>&
  SquareMatrix<dcomplex>::operator=( const ScaledSquareMatrix<double, dcomplex>& );
  template SquareMatrix<dcomplex>&
  SquareMatrix<dcomplex>::operator=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator=( const SquareMatrix<MatsT> &other ) {
    if (this != &other) { // self-assignment check expected
      if (N_ != other.N_)
        CErr("Cannot assign SquareMatrix of different size.");
      std::copy_n(other.ptr_, N_*N_, ptr_);
    }
    return *this;
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator=( const SquareMatrix<double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator=( const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator=( SquareMatrix<MatsT> &&other ) {
    if (this != &other) { // self-assignment check expected
      if (N_ != other.N_)
        CErr("Cannot assign SquareMatrix of different size.");
      memManager_.free(ptr_);
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator=( SquareMatrix<double>&& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator=( SquareMatrix<dcomplex>&& );

  template <typename MatsT>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator*=( MatsT scalar ) {
    blas::scal(dimension()*dimension(), scalar, pointer(), 1);
    return *this;
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator*=( double scalar );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator*=( dcomplex scalar );

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator+=( const SquareMatrix<MatsU> &other ) {
    if (dimension() != other.dimension())
      CErr("Cannot add two SquareMatrix of different size.");
    if (std::is_same<MatsT, double>::value and
        std::is_same<MatsU, dcomplex>::value)
      CErr("Cannot assign a complex SquareMatrix object to a real one.");
    blas::axpy(dimension()*dimension(), 1.0, other.pointer(), 1, pointer(), 1);
    return *this;
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator+=( const SquareMatrix<double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const SquareMatrix<double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> SquareMatrix<MatsT>::operator+( const SquareMatrix<MatsU> &other ) const {
    SquareMatrix<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> sum(*this);
    sum += other;
    return sum;
  }
  template SquareMatrix<double> SquareMatrix<double>::operator+( const SquareMatrix<double>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::operator+( const SquareMatrix<double>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<double>::operator+( const SquareMatrix<dcomplex>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::operator+( const SquareMatrix<dcomplex>& ) const;

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> SquareMatrix<MatsT>::operator-( const SquareMatrix<MatsU> &other ) const {
    return *this + (-other);
  }
  template SquareMatrix<double> SquareMatrix<double>::operator-( const SquareMatrix<double>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::operator-( const SquareMatrix<double>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<double>::operator-( const SquareMatrix<dcomplex>& ) const;
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::operator-( const SquareMatrix<dcomplex>& ) const;

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator+=( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) {
    if (dimension() != scaled.matrix().dimension())
      CErr("Cannot add two SquareMatrix of different size.");
    if (scaled.isPauli() and dynamic_cast<const PauliSpinorSquareMatrices<MatsU>&>(scaled.matrix()).hasZ())
      CErr("Cannot assign a SquareMatrix from a PauliSpinorSquareMatrices with XYZ components.");
    blas::axpy(dimension()*dimension(), scaled.scalar(), scaled.matrix().pointer(), 1, pointer(), 1);
    return *this;
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator+=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const ScaledSquareMatrix<dcomplex, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const ScaledSquareMatrix<double, dcomplex>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator+=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator-=( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) {
    return operator+=(-scaled);
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator-=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const ScaledSquareMatrix<double, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const ScaledSquareMatrix<dcomplex, double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const ScaledSquareMatrix<double, dcomplex>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const ScaledSquareMatrix<dcomplex, dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  SquareMatrix<MatsT>::operator+( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) const {
    PauliSpinorSquareMatrices<typename std::conditional<
        (std::is_same<MatsT, dcomplex>::value or
         std::is_same<MatsU, dcomplex>::value or
         std::is_same<ScalarT, dcomplex>::value),
        dcomplex, double>::type> sum(*this);
    sum += scaled;
    return sum;
  }
  template PauliSpinorSquareMatrices<double> SquareMatrix<double>::operator+( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator+( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator+( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator+( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator+( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator+( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator+( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator+( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;

  template PauliSpinorSquareMatrices<double>
  operator+( const ScaledSquareMatrix<double, double>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, double>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, dcomplex>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, dcomplex>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, double>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, double>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<double, dcomplex>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator+( const ScaledSquareMatrix<dcomplex, dcomplex>&, const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  template <typename ScalarT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type>
  SquareMatrix<MatsT>::operator-( const ScaledSquareMatrix<ScalarT, MatsU> &scaled ) const {
    return *this + (-scaled);
  }
  template PauliSpinorSquareMatrices<double> SquareMatrix<double>::operator-( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator-( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator-( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<double>::operator-( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator-( const ScaledSquareMatrix<double, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator-( const ScaledSquareMatrix<dcomplex, double>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator-( const ScaledSquareMatrix<double, dcomplex>& ) const;
  template PauliSpinorSquareMatrices<dcomplex> SquareMatrix<dcomplex>::operator-( const ScaledSquareMatrix<dcomplex, dcomplex>& ) const;

  template PauliSpinorSquareMatrices<double>
  operator-( const ScaledSquareMatrix<double, double>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, double>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, dcomplex>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, dcomplex>&, const SquareMatrix<double>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, double>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, double>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<double, dcomplex>&, const SquareMatrix<dcomplex>& );
  template PauliSpinorSquareMatrices<dcomplex>
  operator-( const ScaledSquareMatrix<dcomplex, dcomplex>&, const SquareMatrix<dcomplex>& );

  template <typename MatsT>
  template <typename MatsU>
  SquareMatrix<MatsT>& SquareMatrix<MatsT>::operator-=( const SquareMatrix<MatsU> &other ) {
    return operator+=(-other);
  }
  template SquareMatrix<double>& SquareMatrix<double>::operator-=( const SquareMatrix<double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const SquareMatrix<double>& );
  template SquareMatrix<dcomplex>& SquareMatrix<dcomplex>::operator-=( const SquareMatrix<dcomplex>& );
  
  // default TRANS is 'T'
  template <typename MatsT>
  SquareMatrix<MatsT> SquareMatrix<MatsT>::T(char TRANS) {
    SquareMatrix<MatsT> out(memManager(), dimension());
    SetMat(TRANS,dimension(),dimension(),MatsT(1.0),pointer(),dimension(),out.pointer(),dimension());
    return out;
  }
  template SquareMatrix<double>   SquareMatrix<double>::T(char TRANS);
  template SquareMatrix<dcomplex> SquareMatrix<dcomplex>::T(char TRANS);


  template <typename MatsT>
  void SquareMatrix<MatsT>::setTriangle(blas::Uplo upLo, MatsT value, bool setDiag, MatsT diagValue) {
    switch (upLo) {
    case blas::Uplo::Upper:
      for (size_t j = 1; j < N_; j++)
        for (size_t i = 0; i < j; i++)
          operator()(i,j) = value;
      break;
    case blas::Uplo::Lower:
      for (size_t j = 0; j < N_; j++)
        for (size_t i = j + 1; i < N_; i++)
          operator()(i,j) = value;
      break;
    }

    if (setDiag)
      for (size_t i = 0; i < N_; i++)
        operator()(i,i) = diagValue;
  }
  template void SquareMatrix<double>::setTriangle(blas::Uplo upLo, double value, bool setDiag, double diagValue);
  template void SquareMatrix<dcomplex>::setTriangle(blas::Uplo upLo, dcomplex value, bool setDiag, dcomplex diagValue);


}; // namespace ChronusQ
