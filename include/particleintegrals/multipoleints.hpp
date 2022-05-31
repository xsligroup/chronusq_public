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
#include <particleintegrals/onepints.hpp>
#include <particleintegrals/vectorints.hpp>

namespace ChronusQ {

  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  dipole, quadrupole, and octupole integral matrices in a finite
   *  basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class MultipoleInts : public ParticleIntegrals {

    template <typename IntsU>
    friend class MultipoleInts;

  protected:
    size_t highOrder_ = 0;
    bool symmetric_ = false;
    std::vector<VectorInts<IntsT>> components_;

    void orderCheck(size_t order) const {
      if (order > highOrder())
        CErr("Order of XYZ components exceeds highest order in this MultipoleInts.");
      if (order <= 0)
        CErr("No zeroth-order components in MultipoleInts.");
    }

    size_t cumeComponents(size_t order) const {
      size_t count = 0;
      for (size_t i = 1; i <= order; i++)
        count += VectorInts<IntsT>::nComponentsOnOrder(i, symmetric());
      return count;
    }

    std::pair<size_t, size_t> index(std::vector<size_t> comps) const {
      size_t order = comps.size();
      orderCheck(order--);
      return std::make_pair(order, components_[order].index(comps));
    }

    std::pair<size_t, size_t> index(size_t i) const {
      size_t order = 1, orderNComp = 3;
      while (i >= orderNComp) {
        i -= orderNComp;
        orderNComp = VectorInts<IntsT>::nComponentsOnOrder(++order, symmetric());
      }
      return std::make_pair(order-1, i);
    }

    std::string label(size_t i) const {
      std::pair<size_t, size_t> indices = index(i);
      return VectorInts<IntsT>::indexToLabel(
            indices.second, indices.first+1, symmetric());
    }

  public:

    // Constructor
    MultipoleInts() = delete;
    MultipoleInts( const MultipoleInts & ) = default;
    MultipoleInts( MultipoleInts && ) = default;
    MultipoleInts(CQMemManager &mem, size_t nb, size_t order, bool symm):
        ParticleIntegrals(mem, nb), highOrder_(order), symmetric_(symm) {
      if (order == 0)
        CErr("MultipoleInts order must be at least 1.");
      components_.reserve(order);
      for (size_t i = 1; i <= order; i++) {
        components_.emplace_back(mem, nb, i, symm);
      }
    }

    template <typename IntsU>
    MultipoleInts( const MultipoleInts<IntsU> &other, int = 0 ):
        ParticleIntegrals(other),
        highOrder_(other.highOrder_), symmetric_(other.symmetric_) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real MultipoleInts from a Complex one.");
      components_.reserve(other.components_.size());
      for (auto &p : other.components_)
        components_.emplace_back(p);
    }

    MultipoleInts& operator=( const MultipoleInts &other ) {
      if (this != &other) {
        highOrder_ = other.highOrder_;
        symmetric_ = other.symmetric_;
        components_.clear();
        for (const auto & comp : other.components_)
          components_.push_back(comp);
      }
      return *this;
    }
    MultipoleInts& operator=( MultipoleInts &&other ) {
      if (this != &other) {
        highOrder_ = other.highOrder_;
        symmetric_ = other.symmetric_;
        components_.clear();
        for (const auto & comp : other.components_)
          components_.push_back(comp);
      }
      return *this;
    }

    size_t size() const { return cumeComponents(highOrder()); }
    bool symmetric() const { return symmetric_; }
    size_t highOrder() const { return highOrder_; }

    OnePInts<IntsT>& operator[](size_t i) {
      std::pair<size_t, size_t> indices = index(i);
      return components_[indices.first][indices.second];
    }

    const OnePInts<IntsT>& operator[](size_t i) const {
      std::pair<size_t, size_t> indices = index(i);
      return components_[indices.first][indices.second];
    }

    OnePInts<IntsT>& operator[](std::string s) {
      orderCheck(s.size());
      return components_[s.size()-1][s];
    }

    const OnePInts<IntsT>& operator[](std::string s) const {
      orderCheck(s.size());
      return components_[s.size()-1][s];
    }

    VectorInts<IntsT>& getByOrder(size_t order) {
      orderCheck(order--);
      return components_[order];
    }

    const VectorInts<IntsT>& getByOrder(size_t order) const {
      orderCheck(order--);
      return components_[order];
    }

    std::vector<IntsT*> pointersByOrder(size_t order) {
      return getByOrder(order).pointers();
    }

    std::vector<IntsT*> dipolePointers() {
      return pointersByOrder(1);
    }

    std::vector<IntsT*> quadrupolePointers() {
      return pointersByOrder(2);
    }

    std::vector<IntsT*> octupolePointers() {
      return pointersByOrder(3);
    }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) 
    {
      CErr("computeAOInts with two basis sets is NYI for MultipoleInts class");
    }

    virtual void clear() {
      for (VectorInts<IntsT>& c : components_)
        c.clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "MultipoleInts.";
        else
          oeiStr = "MultipoleInts[" + s + "].";
        for (size_t i = 0; i < size(); i++) {
          std::string label_i = label(i);
          prettyPrintSmart(out, oeiStr+label_i, operator[](label_i).pointer(),
                           this->nBasis(), this->nBasis(), this->nBasis());
        }
      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Multipole integral";
        else
          oeiStr = "MultipoleInts[" + s + "]";
        out << oeiStr << " up to ";
        switch (highOrder()) {
        case 1:
          out << "Dipole";
          break;
        case 2:
          out << "Quadrupole";
          break;
        case 3:
          out << "Octupole";
          break;
        default:
          out << highOrder() << "-th order";
          break;
        }
        out << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t highOrder_bcast = highOrder_;
        bool symmetric_bcast = symmetric_;
        MPIBCast(highOrder_bcast,root,comm);
        MPIBCast(symmetric_bcast,root,comm);

        if (highOrder_bcast != highOrder_ or symmetric_bcast != symmetric_) {
          highOrder_ = highOrder_bcast;
          symmetric_ = symmetric_bcast;
          components_.clear();
          components_.reserve(highOrder_);
          for (size_t i = 1; i <= highOrder_; i++) {
            components_.emplace_back(memManager_, NB, i, symmetric_);
          }
        }

        for (VectorInts<IntsT>& comp : components_)
          comp.broadcast(comm, root);
      }
#endif
    }

    template <typename IntsU>
    MultipoleInts<IntsU> spatialToSpinBlock() const {
      MultipoleInts<IntsU> spinBlockInts(memManager_, NB * 2, highOrder_, symmetric_);
      size_t size = components_.size();
      for (size_t i = 0; i < size; i++) {
        spinBlockInts.components_[i] = components_[i].template spatialToSpinBlock<IntsU>();
      }
      return spinBlockInts;
    }

    template <typename TransT>
    MultipoleInts<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      MultipoleInts<typename std::conditional<
      (std::is_same<IntsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(
          memManager(), NT, highOrder(), symmetric());
      for (size_t i = 0; i < size(); i++)
        transInts[i] = (*this)[i].transform(TRANS, T, NT, LDT);
      return transInts;
    }

    ~MultipoleInts() {}

  }; // class MultipoleInts

}; // namespace ChronusQ
