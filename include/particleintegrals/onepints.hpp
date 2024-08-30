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
#include <particleintegrals.hpp>
#include <matrix.hpp>
#include <libint2/shell.h>

namespace ChronusQ {

  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  one electron integral matrix in a finite basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class OnePInts : public ParticleIntegrals {

    template <typename IntsU>
    friend class OnePInts;

  protected:
    cqmatrix::Matrix<IntsT> mat_; ///< One Particle integrals (2 index)

  public:

    // Constructor
    OnePInts() = delete;
    OnePInts(size_t nb):
        ParticleIntegrals(nb), mat_(nb) {}
    OnePInts( const OnePInts &other ) = default;
    template <typename IntsU>
    OnePInts( const OnePInts<IntsU> &other, int = 0 ):
        ParticleIntegrals(other), mat_(other.mat_) {}
    OnePInts( OnePInts &&other ) = default;
    OnePInts( const cqmatrix::Matrix<IntsT> &other ):
        ParticleIntegrals(other.dimension()),
        mat_(other) {}
    template <typename IntsU>
    OnePInts( const cqmatrix::Matrix<IntsU> &other, int = 0 ):
        ParticleIntegrals(other.dimension()),
        mat_(other) {}
    OnePInts( cqmatrix::Matrix<IntsT> &&other ):
        ParticleIntegrals(other.dimension()),
        mat_(std::move(other)) {}

    OnePInts& operator=( const OnePInts &other ) {
      if (this != &other) {
        NB = other.nBasis();
        mat_ = other.mat_;
      }
      return *this;
    }
    OnePInts& operator=( OnePInts &&other ) {
      if (this != &other) {
        NB = other.nBasis();
        mat_ = std::move(other.mat_);
      }
      return *this;
    }
    template <typename IntsU>
    OnePInts& operator=( const cqmatrix::Matrix<IntsU> &other ) {
      NB = other.dimension();
      mat_ = other;
      return *this;
    }
    template <typename IntsU>
    OnePInts& operator=( cqmatrix::Matrix<IntsU> &&other ) {
      NB = other.dimension();
      mat_ = std::move(other);
      return *this;
    }

    IntsT& operator()(size_t p, size_t q) {
      return mat_(p,q);
    }
    IntsT operator()(size_t p, size_t q) const {
      return mat_(p,q);
    }

    // Matrix direct access
    cqmatrix::Matrix<IntsT>& matrix() { return mat_; }
    const cqmatrix::Matrix<IntsT>& matrix() const { return mat_; }
    IntsT* pointer() { return mat_.pointer(); }
    const IntsT* pointer() const { return mat_.pointer(); }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override;

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override {
      
      CErr("One-Particle integral evaluation using two different basis is not implemented"); 

    };

    virtual void clear() override { mat_.clear(); }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      std::string opiStr;
      if (s == "")
        opiStr = "One-Particle integral";
      else
        opiStr = "OPI[" + s + "]";
      if (printFull)
        prettyPrintSmart(out, opiStr, pointer(), this->nBasis(),
                         this->nBasis(), this->nBasis());
      else {
        out << opiStr << std::endl;
      }
    }

    void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);
      mat_.broadcast(comm, root);
    }

    template <typename IntsU>
    cqmatrix::PauliSpinorMatrices<IntsU> spinScatter() const {
      return mat_.template spinScatter<IntsU>();
    }

    template <typename IntsU>
    OnePInts<IntsU> spatialToSpinBlock() const {
      return mat_.template spatialToSpinBlock<IntsU>();
    }

    template <typename TransT>
    OnePInts<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      return mat_.transform(TRANS, T, NT, LDT);
    }

    template <typename TransT, typename OutT>
    void subsetTransform(
        char TRANS, const TransT* T, int LDT,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        OutT* out, bool increment = false) const {
      return mat_.subsetTransform(TRANS, T, LDT, off_size, out, increment);
    }

    static void OnePDriverLibint(libint2::Operator, Molecule&,
        BasisSet&, std::vector<IntsT*>, Particle p, size_t deriv=0, size_t S0a=0);
    void OnePDriverLibcint(OPERATOR, const Molecule&,
        const BasisSet&, const HamiltonianOptions&);
    template <size_t NOPER, bool SYMM, typename F>
    static void OnePDriverLocal(const F&,
        std::vector<libint2::Shell>&, std::vector<IntsT*>);

    // Pointer convertor
    template <typename IntsU>
    static std::shared_ptr<OnePInts<IntsU>>
    convert(const std::shared_ptr<OnePInts<IntsT>>&);

    ~OnePInts() {}

  }; // class OnePInts

}; // namespace ChronusQ
