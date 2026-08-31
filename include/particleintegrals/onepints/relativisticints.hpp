/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <particleintegrals/onepints.hpp>
#include <matrix.hpp>

namespace ChronusQ {
  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  one electron integral matrices O, pOp, and three pOxp in
   *  a finite basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class OnePRelInts : public OnePInts<IntsT> {

    template <typename IntsU>
    friend class OnePRelInts;

  protected:
    cqmatrix::PauliSpinorMatrices<IntsT> smallComponent_;

  public:

    // Constructor
    OnePRelInts() = delete;
    OnePRelInts( const OnePRelInts & ) = default;
    OnePRelInts( OnePRelInts && ) = default;
    OnePRelInts(size_t nb, bool SORelativistic):
        OnePInts<IntsT>(nb),
        smallComponent_(nb, SORelativistic, SORelativistic) {}

    template <typename IntsU>
    OnePRelInts( const OnePRelInts<IntsU> &other, int = 0 ):
        OnePInts<IntsT>(other.nBasis()),
        smallComponent_(other.smallComponent_) {}

    OnePRelInts( const cqmatrix::PauliSpinorMatrices<IntsT> &other ):
    OnePInts<IntsT>(other.nRows()),
    smallComponent_(other) {}

    template <typename IntsU>
    OnePRelInts( const cqmatrix::PauliSpinorMatrices<IntsU> &other, int = 0 ):
        OnePInts<IntsT>(other.nRows()),
        smallComponent_(other) {}

    OnePRelInts( cqmatrix::PauliSpinorMatrices<IntsT> &&other ):
        OnePInts<IntsT>(other.nRows()),
        smallComponent_(std::move(other)) {}

    bool hasSpinOrbit() const { return smallComponent_.hasXY() and smallComponent_.hasZ(); }

    cqmatrix::Matrix<IntsT>& scalar() { return smallComponent_.S(); }
    const cqmatrix::Matrix<IntsT>& scalar() const { return smallComponent_.S(); }
    cqmatrix::Matrix<IntsT>& SOX() { return smallComponent_.X(); }
    const cqmatrix::Matrix<IntsT>& SOX() const { return smallComponent_.X(); }
    cqmatrix::Matrix<IntsT>& SOY() { return smallComponent_.Y(); }
    const cqmatrix::Matrix<IntsT>& SOY() const { return smallComponent_.Y(); }
    cqmatrix::Matrix<IntsT>& SOZ() { return smallComponent_.Z(); }
    const cqmatrix::Matrix<IntsT>& SOZ() const { return smallComponent_.Z(); }

    cqmatrix::PauliSpinorMatrices<IntsT>& SZYX() { return smallComponent_; }
    const cqmatrix::PauliSpinorMatrices<IntsT>& SZYX() const { return smallComponent_; }

    std::vector<IntsT*> SOXYZPointers() {
      if (!hasSpinOrbit())
        return std::vector<IntsT*>();
      return { SOX().pointer(), SOY().pointer(), SOZ().pointer() };
    }

    template <typename IntsU>
    cqmatrix::Matrix<IntsU> formW() const;

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override;

    void OnePRelDriverLibcint(const Molecule&,
        const BasisSet&, const HamiltonianOptions &options);

    virtual void clear() override {
      OnePInts<IntsT>::clear();
      smallComponent_.clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "RelOPI";
        else
          oeiStr = "RelOPI[" + s + "]";
        prettyPrintSmart(out, oeiStr+".LL", this->pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        prettyPrintSmart(out, oeiStr+".SS.S", scalar().pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        if(this->hasSpinOrbit()) {
          prettyPrintSmart(out, oeiStr+".SS.X", SOX().pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
          prettyPrintSmart(out, oeiStr+".SS.Y", SOY().pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
          prettyPrintSmart(out, oeiStr+".SS.Z", SOZ().pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
        }
      } else {
        std::string opiStr;
        if (s == "")
          opiStr = "Relativistic one-particle integral";
        else
          opiStr = "RelOPI[" + s + "]";
        out << opiStr;
        if(this->hasSpinOrbit())
          out << " with";
        else
          out << " without";
        out << " spin-orbit integrals" << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      OnePInts<IntsT>::broadcast(comm, root);
      smallComponent_.broadcast();
    }

    template <typename TransT>
    OnePRelInts<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
      char TRANS, const TransT* T, int NT, int LDT) const {
        OnePRelInts<typename std::conditional<
        (std::is_same<IntsT, dcomplex>::value or
        std::is_same<TransT, dcomplex>::value),
        dcomplex, double>::type> transInts(NT, hasSpinOrbit());

        transInts.matrix() = this->matrix().transform(TRANS, T, NT, LDT);
        transInts.smallComponent_ = smallComponent_.transform(TRANS, T, NT, LDT);
        return transInts;
      }

    virtual ~OnePRelInts() {}

  }; // class OneERelInts

}; // namespace ChronusQ
