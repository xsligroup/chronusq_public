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
#include <memmanager.hpp>
#include <cerr.hpp>
#include <util/matout.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  template <typename ScalarT, typename MatsT>
  class ScaledSquareMatrix;

  template <typename MatsT>
  class PauliSpinorSquareMatrices;

  template <typename MatsT>
  class SquareMatrix {

    template <typename MatsU>
    friend class SquareMatrix;

  protected:
    size_t N_;
    CQMemManager &memManager_; ///< CQMemManager to allocate matricies
    MatsT *ptr_ = nullptr;     ///< Raw matrix storage (2 index)

  public:

    // Constructor
    SquareMatrix() = delete;
    SquareMatrix(CQMemManager &mem, size_t n):
        N_(n), memManager_(mem) {
      malloc();
    }
    SquareMatrix( const SquareMatrix &other ):
        SquareMatrix(other.memManager_, other.N_) {
      std::copy_n(other.ptr_, N_*N_, ptr_);
    }
    template <typename MatsU>
    SquareMatrix( const SquareMatrix<MatsU> &other, int = 0 ):
        SquareMatrix(other.memManager_, other.N_) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real SquareMatrix from a Complex one.");
      std::copy_n(other.ptr_, N_*N_, ptr_);
    }
    SquareMatrix( SquareMatrix &&other ):
        N_(other.N_), memManager_(other.memManager_),
        ptr_(other.ptr_) { other.ptr_ = nullptr; }
    template <typename MatsU>
    SquareMatrix( const PauliSpinorSquareMatrices<MatsU> &other ):
        SquareMatrix(other.memManager_, other.N_) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real SquareMatrix from a Complex one.");
      if (other.hasZ())
        CErr("Cannot create a SquareMatrix from a PauliSpinorSquareMatrices"
             " with XYZ components.");
      std::copy_n(other.pointer(), N_*N_, ptr_);
    }
    SquareMatrix( PauliSpinorSquareMatrices<MatsT> &&other ):
        N_(other.N_), memManager_(other.memManager_),
        ptr_(other.ptr_) {
      if (other.hasZ())
        CErr("Cannot create a SquareMatrix from a PauliSpinorSquareMatrices"
             " with XYZ components.");
      other.ptr_ = nullptr;
    }
    template <typename ScalarT, typename MatsU>
    SquareMatrix( const ScaledSquareMatrix<ScalarT, MatsU>& );

    SquareMatrix& operator=( const SquareMatrix &other );
    SquareMatrix& operator=( SquareMatrix &&other );

    template <typename ScalarT, typename MatsU>
    SquareMatrix& operator=( const ScaledSquareMatrix<ScalarT, MatsU>& );

    CQMemManager& memManager() const { return memManager_; }
    size_t dimension() const{ return N_; }

    SquareMatrix& operator*=( MatsT );
    ScaledSquareMatrix<double, MatsT> operator-() const {
      return ScaledSquareMatrix<double, MatsT>(-1.0, *this);
    }

    template <typename MatsU>
    SquareMatrix& operator+=( const SquareMatrix<MatsU>& );
    template <typename MatsU>
    SquareMatrix& operator-=( const SquareMatrix<MatsU>& );
    template <typename MatsU>
    SquareMatrix<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator+( const SquareMatrix<MatsU>& ) const;
    template <typename MatsU>
    SquareMatrix<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator-( const SquareMatrix<MatsU>& ) const;

    template <typename ScalarT, typename MatsU>
    SquareMatrix& operator+=( const ScaledSquareMatrix<ScalarT, MatsU>& );
    template <typename ScalarT, typename MatsU>
    SquareMatrix& operator-=( const ScaledSquareMatrix<ScalarT, MatsU>& );
    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value or
     std::is_same<ScalarT, dcomplex>::value),
    dcomplex, double>::type>
    operator+( const ScaledSquareMatrix<ScalarT, MatsU>& ) const;
    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value or
     std::is_same<ScalarT, dcomplex>::value),
    dcomplex, double>::type>
    operator-( const ScaledSquareMatrix<ScalarT, MatsU>& ) const;

    MatsT& operator()(size_t p, size_t q) {
      return ptr_[p + q*N_];
    }
    MatsT operator()(size_t p, size_t q) const {
      return ptr_[p + q*N_];
    }

    // Matrix direct access
    MatsT* pointer() { return ptr_; }
    const MatsT* pointer() const { return ptr_; }

    virtual std::vector<MatsT*> SZYXPointers() {
      return { pointer() };
    }

    SquareMatrix<double> real_part() {
      SquareMatrix<double> realMat(memManager_, N_);
      GetMatRE('N',N_,N_,1.,pointer(),N_,realMat.pointer(),N_);
      return realMat;
    }
    
    // transform and return the transformed matrix
    SquareMatrix<MatsT> T(char TRANS = 'T');
    
    void clear() {
      std::fill_n(ptr_,N_*N_,MatsT(0.));
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      std::string matStr;
      if (s == "")
        matStr = "Square Matrix";
      else
        matStr = "Square Matrix[" + s + "]";
      if (printFull)
        prettyPrintSmart(out, matStr, pointer(), N_, N_, N_);
      else {
        out << matStr << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

#ifdef CQ_ENABLE_MPI
      // BCast matrix to all MPI processes
      if( MPISize(comm) > 1 ) {
        std::cerr  << "  *** Scattering a matrix ***\n";
        size_t N_bcast = N_;
        MPIBCast(N_bcast,root,comm);

        if (N_bcast != N_) {
          N_ = N_bcast;
          malloc();
        }

        MPIBCast(ptr_,N_*N_,root,comm);
      }
#endif

    }

    template <typename MatsU>
    PauliSpinorSquareMatrices<MatsU> spinScatter(
        bool hasXY = true, bool hasZ = true) const;

    template <typename MatsU>
    SquareMatrix<MatsU> spatialToSpinBlock() const;
    
    template <typename MatsU>
    void componentScatter(SquareMatrix<MatsU> & LL,
                          SquareMatrix<MatsU> & LS,
                          SquareMatrix<MatsU> & SL,
                          SquareMatrix<MatsU> & SS,
                          bool increment = false) const;
     
    template <typename MatsU>
    void componentGather(const SquareMatrix<MatsU> & LL,
                         const SquareMatrix<MatsU> & LS,
                         const SquareMatrix<MatsU> & SL,
                         const SquareMatrix<MatsU> & SS,
                         bool increment = false);
     
    template <typename MatsU>
    static SquareMatrix<MatsT>
    componentGatherBuild(const SquareMatrix<MatsU> & LL,
                         const SquareMatrix<MatsU> & LS,
                         const SquareMatrix<MatsU> & SL,
                         const SquareMatrix<MatsU> & SS);
    
    
    
    template <typename TransT>
    SquareMatrix<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const;

    template <typename TransT, typename OutT>
    void subsetTransform(
        char TRANS, const TransT* T, int LDT,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        OutT* out, bool increment = false) const;

    void malloc() {
      if (ptr_) memManager_.free(ptr_);
      size_t N2 = N_*N_;
      if (N2 != 0) {
        try { ptr_ = memManager_.malloc<MatsT>(N2); }
        catch(...) {
          std::cout << std::fixed;
          std::cout << "Insufficient memory for the full INTS matrix ("
                    << (N2/1e9) * sizeof(double) << " GB)" << std::endl;
          std::cout << std::endl << memManager_ << std::endl;
          throw std::bad_alloc();
        }
      }
    }

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<SquareMatrix<MatsU>>
    convert(const std::shared_ptr<SquareMatrix<MatsT>>&);

    ~SquareMatrix() {
      if(ptr_) memManager_.free(ptr_);
    }

  }; // class SquareMatrix

  template <typename ScalarT, typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type> operator+(
      const ScaledSquareMatrix<ScalarT, MatsT> &lhs, const SquareMatrix<MatsU> &rhs ) {
    return rhs + lhs;
  }

  template <typename ScalarT, typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type> operator-(
      const ScaledSquareMatrix<ScalarT, MatsT> &lhs, const SquareMatrix<MatsU> &rhs ) {
    return lhs + (-rhs);
  }

  template <typename MatsT>
  std::ostream& operator<<(std::ostream&, const SquareMatrix<MatsT>&);

}; // namespace ChronusQ
