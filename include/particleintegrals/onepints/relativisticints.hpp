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
#include <particleintegrals/onepints.hpp>

namespace ChronusQ {

  enum class REL_INTS_COMPS : size_t {
    SCALAR = 0, SOZ = 1, SOY = 2, SOX = 3, O = 4
  };

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
    std::vector<OnePInts<IntsT>> components_;

  public:

    // Constructor
    OnePRelInts() = delete;
    OnePRelInts( const OnePRelInts & ) = default;
    OnePRelInts( OnePRelInts && ) = default;
    OnePRelInts(size_t nb, bool SORelativistic):
        OnePInts<IntsT>(nb) {
      size_t nRel = SORelativistic ? 4 : 1;
      components_.reserve(nRel);
      for (size_t i = 0; i < nRel; i++)
        components_.emplace_back(nb);
    }

    template <typename IntsU>
    OnePRelInts( const OnePRelInts<IntsU> &other, int = 0 ):
        OnePInts<IntsT>(other, 0) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real OnePRelInts from a Complex one.");
      components_.reserve(other.components_.size());
      for (auto &p : other.components_)
        components_.emplace_back(p);
    }

    bool hasSpinOrbit() const { return components_.size() == 4; }

    OnePInts<IntsT>& operator[](REL_INTS_COMPS comp) {
      if (comp == REL_INTS_COMPS::O)
        return *this;
      size_t i = static_cast<size_t>(comp);
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return components_[i];
    }

    const OnePInts<IntsT>& operator[](REL_INTS_COMPS comp) const {
      if (comp == REL_INTS_COMPS::O)
        return *this;
      size_t i = static_cast<size_t>(comp);
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return components_[i];
    }

    OnePInts<IntsT>& scalar() { return operator[](REL_INTS_COMPS::SCALAR); }
    const OnePInts<IntsT>& scalar() const { return operator[](REL_INTS_COMPS::SCALAR); }
    OnePInts<IntsT>& SOX() { return operator[](REL_INTS_COMPS::SOX); }
    const OnePInts<IntsT>& SOX() const { return operator[](REL_INTS_COMPS::SOX); }
    OnePInts<IntsT>& SOY() { return operator[](REL_INTS_COMPS::SOY); }
    const OnePInts<IntsT>& SOY() const { return operator[](REL_INTS_COMPS::SOY); }
    OnePInts<IntsT>& SOZ() { return operator[](REL_INTS_COMPS::SOZ); }
    const OnePInts<IntsT>& SOZ() const { return operator[](REL_INTS_COMPS::SOZ); }

    std::vector<OnePInts<IntsT>>& SZYX() {
      return components_;
    }
    const std::vector<OnePInts<IntsT>>& SZYX() const {
      return components_;
    }
    std::vector<IntsT*> SOXYZPointers() {
      if (this->components_.size() <= 1)
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
      for (OnePInts<IntsT>& c : components_)
        c.clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "RelOPI";
        else
          oeiStr = "RelOPI[" + s + "]";
        prettyPrintSmart(out, oeiStr+".O", this->pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        prettyPrintSmart(out, oeiStr+".Scalar", scalar().pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        if(this->hasSpinOrbit()) {
          prettyPrintSmart(out, oeiStr+".SOX", SOX().pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
          prettyPrintSmart(out, oeiStr+".SOY", SOY().pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
          prettyPrintSmart(out, oeiStr+".SOZ", SOZ().pointer(),
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

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t nRel = components_.size();
        MPIBCast(nRel,root,comm);

        if (components_.size() != nRel) {
          components_.clear();
          components_.reserve(nRel);
          for (size_t i = 0; i < nRel; i++) {
            components_.emplace_back(this->NB);
          }
        }

        for (OnePInts<IntsT>& comp : components_)
          comp.broadcast(comm, root);
      }
#endif
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
      transInts[REL_INTS_COMPS::O] =
          (*this)[REL_INTS_COMPS::O].transform(TRANS, T, NT, LDT);
      transInts.scalar() = scalar().transform(TRANS, T, NT, LDT);
      if (hasSpinOrbit()) {
        transInts.SOX() = SOX().transform(TRANS, T, NT, LDT);
        transInts.SOY() = SOY().transform(TRANS, T, NT, LDT);
        transInts.SOZ() = SOZ().transform(TRANS, T, NT, LDT);
      }
      return transInts;
    }

    virtual ~OnePRelInts() {}

  }; // class OneERelInts

}; // namespace ChronusQ
