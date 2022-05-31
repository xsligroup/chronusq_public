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
#include <matrix.hpp>

namespace ChronusQ {

  enum class PAULI_SPINOR_COMPS : size_t {
    S = 0, Z = 1, Y = 2, X = 3
  };
  
  /**
   *  \brief Templated class to handle the storage of
   *  one-electron two-spinor integral matrices
   *  Scalar, X, Y, and Z in Pauli representation.
   *
   */
  template <typename MatsT>
  class PauliSpinorSquareMatrices : public SquareMatrix<MatsT> {

    template <typename MatsU>
    friend class PauliSpinorSquareMatrices;
    template <typename MatsU>
    friend class SquareMatrix;

  protected:
    std::vector<SquareMatrix<MatsT>> components_;

  public:

    // Constructor
    PauliSpinorSquareMatrices() = delete;
    PauliSpinorSquareMatrices( const PauliSpinorSquareMatrices & ) = default;
    PauliSpinorSquareMatrices( PauliSpinorSquareMatrices && ) = default;
    PauliSpinorSquareMatrices(CQMemManager &mem, size_t n, bool hasXY = true, bool hasZ = true):
        SquareMatrix<MatsT>(mem, n) {
      size_t nZYX = hasXY ? 3 : hasZ ? 1 : 0;
      components_.reserve(nZYX);
      for (size_t i = 0; i < nZYX; i++)
        components_.emplace_back(mem, n);
    }
    template <typename MatsU>
    PauliSpinorSquareMatrices( const SquareMatrix<MatsU> &other,
                               bool addXY = false, bool addZ = false ):
        SquareMatrix<MatsT>(other, 0) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real PauliSpinorSquareMatrices from a Complex one.");
      size_t nZYX = addXY ? 3 : addZ ? 1 : 0;
      components_.reserve(nZYX);
      while (nZYX > 0) {
        components_.emplace_back(this->memManager(), this->dimension());
        components_.back().clear();
        nZYX--;
      }
    }
    template <typename MatsU>
    PauliSpinorSquareMatrices( const PauliSpinorSquareMatrices<MatsU> &other,
                               bool addXY = false, bool addZ = false ):
        SquareMatrix<MatsT>(other, 0) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real PauliSpinorSquareMatrices from a Complex one.");
      size_t nZYX = std::max(addXY ? 3ul : addZ ? 1ul : 0ul, other.components_.size());
      components_.reserve(nZYX);
      for (auto &p : other.components_)
        components_.emplace_back(p);
      while (nZYX > other.components_.size()) {
        components_.emplace_back(this->memManager(), this->dimension());
        components_.back().clear();
        nZYX--;
      }
    }
    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices( const ScaledSquareMatrix<ScalarT, MatsU>&,
                               bool addXY = false, bool addZ = false );

    PauliSpinorSquareMatrices& operator=( const PauliSpinorSquareMatrices &other );
    PauliSpinorSquareMatrices& operator=( PauliSpinorSquareMatrices &&other );
    template <typename MatsU>
    PauliSpinorSquareMatrices& operator=( const SquareMatrix<MatsU> &other ) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real PauliSpinorSquareMatrices from a Complex SquareMatrix.");
      if (this->N_ != other.dimension())
        CErr("Cannot assign a SquareMatrix of different size to PauliSpinorSquareMatrices.");
      SquareMatrix<MatsT>::operator=(other);
      for (auto &comp : components_) comp.clear();
      return *this;
    }
    template <typename MatsU>
    PauliSpinorSquareMatrices& operator=( SquareMatrix<MatsU> &&other ) {
      if (std::is_same<MatsU, dcomplex>::value
          and std::is_same<MatsT, double>::value)
        CErr("Cannot create a Real PauliSpinorSquareMatrices from a Complex SquareMatrix.");
      if (this->N_ != other.dimension())
        CErr("Cannot assign a SquareMatrix of different size to PauliSpinorSquareMatrices.");
      SquareMatrix<MatsT>::operator=(std::move(other));
      for (auto &comp : components_) comp.clear();
      return *this;
    }
    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices& operator=( const ScaledSquareMatrix<ScalarT, MatsU>& );

    PauliSpinorSquareMatrices& operator*=( MatsT );

    ScaledSquareMatrix<double, MatsT> operator-() const {
      return ScaledSquareMatrix<double, MatsT>(-1.0, *this);
    }

    template <typename MatsU>
    PauliSpinorSquareMatrices& operator+=( const PauliSpinorSquareMatrices<MatsU>& );
    template <typename MatsU>
    PauliSpinorSquareMatrices& operator-=( const PauliSpinorSquareMatrices<MatsU>& );
    template <typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator+( const PauliSpinorSquareMatrices<MatsU>& ) const;
    template <typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator-( const PauliSpinorSquareMatrices<MatsU>& ) const;

    template <typename MatsU>
    PauliSpinorSquareMatrices& operator+=( const SquareMatrix<MatsU>& );
    template <typename MatsU>
    PauliSpinorSquareMatrices& operator-=( const SquareMatrix<MatsU>& );
    template <typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator+( const SquareMatrix<MatsU>& ) const;
    template <typename MatsU>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<MatsU, dcomplex>::value),
    dcomplex, double>::type> operator-( const SquareMatrix<MatsU>& ) const;

    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices& operator+=( const ScaledSquareMatrix<ScalarT, MatsU>& );
    template <typename ScalarT, typename MatsU>
    PauliSpinorSquareMatrices& operator-=( const ScaledSquareMatrix<ScalarT, MatsU>& );
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

    bool hasXY() const { return components_.size() == 3; }
    bool hasZ() const { return components_.size() >= 1; }
    size_t nComponent() const { return components_.size() + 1; }

    SquareMatrix<MatsT>& operator[](PAULI_SPINOR_COMPS comp) {
      if (comp == PAULI_SPINOR_COMPS::S)
        return *this;
      size_t i = static_cast<size_t>(comp) - 1;
      if (i >= components_.size())
        CErr("Requested component is NOT in this PauliSpinorSquareMatrices object.");
      return components_[i];
    }

    const SquareMatrix<MatsT>& operator[](PAULI_SPINOR_COMPS comp) const {
      if (comp == PAULI_SPINOR_COMPS::S)
        return *this;
      size_t i = static_cast<size_t>(comp) - 1;
      if (i >= components_.size())
        CErr("Requested component is NOT in this PauliSpinorSquareMatrices object.");
      return components_[i];
    }

    SquareMatrix<MatsT>& S() { return operator[](PAULI_SPINOR_COMPS::S); }
    const SquareMatrix<MatsT>& S() const { return operator[](PAULI_SPINOR_COMPS::S); }
    SquareMatrix<MatsT>& X() { return operator[](PAULI_SPINOR_COMPS::X); }
    const SquareMatrix<MatsT>& X() const { return operator[](PAULI_SPINOR_COMPS::X); }
    SquareMatrix<MatsT>& Y() { return operator[](PAULI_SPINOR_COMPS::Y); }
    const SquareMatrix<MatsT>& Y() const { return operator[](PAULI_SPINOR_COMPS::Y); }
    SquareMatrix<MatsT>& Z() { return operator[](PAULI_SPINOR_COMPS::Z); }
    const SquareMatrix<MatsT>& Z() const { return operator[](PAULI_SPINOR_COMPS::Z); }

    virtual std::vector<MatsT*> SZYXPointers() {
      if (hasXY())
        return { S().pointer(), Z().pointer(), Y().pointer(), X().pointer() };
      if (hasZ())
        return { S().pointer(), Z().pointer() };
      return { S().pointer() };
    }
    
    PauliSpinorSquareMatrices<double> real_part() {
      PauliSpinorSquareMatrices<double> realMats(S().real_part(), false, false);
      realMats.components_.reserve(components_.size());
      for (SquareMatrix<MatsT> mat : components_)
        realMats.components_.emplace_back(mat.real_part());
      return realMats;
    }

    void clear() {
      SquareMatrix<MatsT>::clear();
      for (SquareMatrix<MatsT>& c : components_)
        c.clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "PauliSpinorOEI";
        else
          oeiStr = "PauliSpinorOEI[" + s + "]";
        prettyPrintSmart(out, oeiStr+".S", S().pointer(),
            this->dimension(), this->dimension(), this->dimension());
        if(hasZ())
          prettyPrintSmart(out, oeiStr+".Z", Z().pointer(),
              this->dimension(), this->dimension(), this->dimension());
        if(hasXY()) {
          prettyPrintSmart(out, oeiStr+".Y", Y().pointer(),
              this->dimension(), this->dimension(), this->dimension());
          prettyPrintSmart(out, oeiStr+".X", X().pointer(),
              this->dimension(), this->dimension(), this->dimension());
        }
      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Pauli spinor one-electron integrals";
        else
          oeiStr = "PauliSpinorOEI[" + s + "]";
        out << oeiStr;
        if(not hasZ())
          out << " without XYZ components";
        else if(not hasXY())
          out << " without XY components";
        out << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      SquareMatrix<MatsT>::broadcast(comm, root);
      for (SquareMatrix<MatsT>& comp : components_)
        comp.broadcast(comm, root);
    }

    template <typename MatsU>
    SquareMatrix<MatsU> spinGather() const;

    template <typename MatsU>
    std::vector<SquareMatrix<MatsU>> spinGatherToBlocks(
        bool genABBA = true, bool genBB = true) const;

    template <typename MatsU>
    static PauliSpinorSquareMatrices<MatsT>
    spinBlockScatterBuild(const SquareMatrix<MatsU> &AA, bool hasXY = false, bool hasZ = false);

    template <typename MatsU>
    static PauliSpinorSquareMatrices<MatsT>
    spinBlockScatterBuild(const SquareMatrix<MatsU> &AA, const SquareMatrix<MatsU> &BB,
                          bool hasXY = false, bool hasZ = true);

    template <typename MatsU>
    static PauliSpinorSquareMatrices<MatsT>
    spinBlockScatterBuild(const SquareMatrix<MatsU> &AA, const SquareMatrix<MatsU> &AB,
                          const SquareMatrix<MatsU> &BA, const SquareMatrix<MatsU> &BB,
                          bool hasXY = true, bool hasZ = true);
    
    template <typename MatsU>
    void componentScatter(PauliSpinorSquareMatrices<MatsU> & LL, 
                          PauliSpinorSquareMatrices<MatsU> & LS,
                          PauliSpinorSquareMatrices<MatsU> & SL,
                          PauliSpinorSquareMatrices<MatsU> & SS,
                          bool increment = false) const; 
    
    template <typename MatsU>
    void componentGather(const PauliSpinorSquareMatrices<MatsU> & LL, 
                         const PauliSpinorSquareMatrices<MatsU> & LS,
                         const PauliSpinorSquareMatrices<MatsU> & SL,
                         const PauliSpinorSquareMatrices<MatsU> & SS,
                         bool increment = false); 

    template <typename MatsU>
    static PauliSpinorSquareMatrices<MatsT> 
    componentGatherBuild(const PauliSpinorSquareMatrices<MatsU> & LL, 
                         const PauliSpinorSquareMatrices<MatsU> & LS,
                         const PauliSpinorSquareMatrices<MatsU> & SL,
                         const PauliSpinorSquareMatrices<MatsU> & SS); 
    
    template <typename MatsU>
    void componentAdd(const char TRANS, MatsU scale, const std::string & comp,
      const PauliSpinorSquareMatrices<MatsU> & pauli);

    void symmetrizeLSSL(char TRANS, bool get_SL_from_LS = true);
    
    template <typename TransT>
    PauliSpinorSquareMatrices<typename std::conditional<
    (std::is_same<MatsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      PauliSpinorSquareMatrices<typename std::conditional<
      (std::is_same<MatsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transMats(this->memManager(), NT, hasXY(), hasZ());
      transMats.S() = S().transform(TRANS, T, NT, LDT);
      if (hasZ())
        transMats.Z() = Z().transform(TRANS, T, NT, LDT);
      if (hasXY()) {
        transMats.Y() = Y().transform(TRANS, T, NT, LDT);
        transMats.X() = X().transform(TRANS, T, NT, LDT);
      }
      return transMats;
    }

    virtual ~PauliSpinorSquareMatrices() {}

  }; // class PauliSpinorSquareMatrices

  template <typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator+(
      const SquareMatrix<MatsT> &lhs, const PauliSpinorSquareMatrices<MatsU> &rhs ) {
    return rhs + lhs;
  }

  template <typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value),
  dcomplex, double>::type> operator-(
      const SquareMatrix<MatsT> &lhs, const PauliSpinorSquareMatrices<MatsU> &rhs ) {
    return lhs + (-rhs);
  }

  template <typename ScalarT, typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type> operator+(
      const ScaledSquareMatrix<ScalarT, MatsT> &lhs, const PauliSpinorSquareMatrices<MatsU> &rhs ) {
    return rhs + lhs;
  }

  template <typename ScalarT, typename MatsT, typename MatsU>
  PauliSpinorSquareMatrices<typename std::conditional<
  (std::is_same<MatsT, dcomplex>::value or
   std::is_same<MatsU, dcomplex>::value or
   std::is_same<ScalarT, dcomplex>::value),
  dcomplex, double>::type> operator-(
      const ScaledSquareMatrix<ScalarT, MatsT> &lhs, const PauliSpinorSquareMatrices<MatsU> &rhs ) {
    return lhs + (-rhs);
  }

}; // namespace ChronusQ
