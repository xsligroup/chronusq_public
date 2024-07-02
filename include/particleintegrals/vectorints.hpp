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
#include <particleintegrals/onepints/relativisticints.hpp>

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
  class VectorInts : public ParticleIntegrals {

    template <typename IntsU>
    friend class VectorInts;

  public:

    static size_t nComponentsOnOrder(size_t order, bool symmetric) {
      if (symmetric) {
        return (order + 2) * (order + 1) / 2;
      } else {
        size_t pow = 3;
        for (size_t i = 2; i <= order; i++) {
          pow *= 3;
        }
        return pow;
      }
    }

    static std::string indexToLabel(size_t i, size_t order, bool symmetric) {
      std::string s;
      if (symmetric) {
        size_t yz = static_cast<size_t>(std::sqrt(2*i+0.25)-0.5);
        std::array<size_t, 3> ncomps{0,0,0};
        ncomps[0] = order - yz;
        ncomps[2] = i - yz*(yz+1)/2;
        ncomps[1] = yz - ncomps[2];
        for (size_t j = 0; j < 3; j++)
          for (size_t k = 0; k < ncomps[j]; k++)
            s += static_cast<char>('X'+j);
      } else
        for (size_t j = 0; j < order; j++) {
          s = static_cast<char>('X'+i%3) + s;
          i /= 3;
        }
      return s;
    }

  protected:
    size_t order_ = 0;
    bool symmetric_ = false;
    std::vector<std::shared_ptr<OnePInts<IntsT>>> components_;

    void orderCheck(size_t order) const {
      if (order != order_)
        CErr("Order of XYZ components does not match order in this VectorInts.");
    }

    size_t nComponents() const {
      return nComponentsOnOrder(order(), symmetric());
    }

    static std::vector<size_t> indices(std::string s) {
      std::vector<size_t> inds(s.size());
      std::transform(s.begin(), s.end(), inds.begin(),
          [](unsigned char c){ return std::toupper(c) - 'X';});
      return inds;
    }

    size_t index(std::vector<size_t> comps) const {
      size_t order = comps.size();
      orderCheck(order);
      size_t count = 0;
      if (symmetric()) {
        std::array<size_t, 3> ncomps{0,0,0};
        for (size_t i = 0; i < order; i++) {
          if (comps[i] > 2)
            CErr("VectorInts Component label must be"
                 " combinations of 'X', 'Y', and 'Z'");
          ncomps[comps[i]]++;
        }
        size_t yz = ncomps[1] + ncomps[2];
        count = yz * (yz+1) / 2 + ncomps[2];
      } else {
        size_t pow = 1;
        for (int i = order - 1; i >= 0; i--) {
          if (comps[i] > 2)
            CErr("VectorInts Component label must be"
                 " combinations of 'X', 'Y', and 'Z'");
          count += comps[i] * pow;
          pow *= 3;
        }
      }
      return count;
    }

  public:

    // Constructor
    VectorInts() = delete;
    VectorInts( const VectorInts & ) = default;
    VectorInts( VectorInts && ) = default;
    VectorInts(size_t nb, size_t order, bool symm):
        ParticleIntegrals(nb), order_(order), symmetric_(symm) {
      if (order == 0)
        CErr("VectorInts order must be at least 1.");
      size_t size = nComponents();
      components_.reserve(size);
      for (size_t i = 0; i < size; i++) {
        components_.emplace_back(std::make_shared<OnePInts<IntsT>>(nb));
      }
    }

    template <typename IntsU>
    VectorInts( const VectorInts<IntsU> &other, int = 0 ):
        ParticleIntegrals(other),
        order_(other.order_), symmetric_(other.symmetric_) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real VectorInts from a Complex one.");
      components_.reserve(other.components_.size());
      for (auto &p : other.components_)
        components_.emplace_back(std::make_shared<OnePInts<IntsT>>(*p));
    }

    VectorInts& operator=( const VectorInts &other ) {
      if (this != &other) {
        order_ = other.order_;
        symmetric_ = other.symmetric_;
        components_.clear();
        for (const auto & comp : other.components_)
          components_.push_back(comp);
      }
      return *this;
    }
    VectorInts& operator=( VectorInts &&other ) {
      if (this != &other) {
        order_ = other.order_;
        symmetric_ = other.symmetric_;
        components_.clear();
        for (const auto & comp : other.components_)
          components_.push_back(comp);
      }
      return *this;
    }

    size_t size() const { return components_.size(); }
    bool symmetric() const { return symmetric_; }
    size_t order() const { return order_; }

    std::shared_ptr<OnePInts<IntsT>>& operator[](size_t i) {
      return components_[i];
    }

    const std::shared_ptr<OnePInts<IntsT>>& operator[](size_t i) const {
      return components_[i];
    }

    std::shared_ptr<OnePInts<IntsT>>& operator[](std::string s) {
      return components_[index(indices(s))];
    }

    const std::shared_ptr<OnePInts<IntsT>>& operator[](std::string s) const {
      return components_[index(indices(s))];
    }

    std::vector<IntsT*> pointers() {
      std::vector<IntsT*> ps(nComponents());
      std::transform(components_.begin(),
          components_.end(), ps.begin(),
          [](std::shared_ptr<OnePInts<IntsT>>& opi){ return opi->pointer(); });
      return ps;
    }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
                               OPERATOR, const HamiltonianOptions&) override;

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override {
      
      CErr("Vector integral evaluation using two different basis is not implemented"); 

    };

    virtual void clear() override {
      for (std::shared_ptr<OnePInts<IntsT>>& c : components_)
        c->clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (printFull) {
        std::string opiStr;
        if (s == "")
          opiStr = "VectorInts.";
        else
          opiStr = "VectorInts[" + s + "].";
        for (size_t i = 0; i < size(); i++)
          prettyPrintSmart(out, opiStr+indexToLabel(i, order(), symmetric()),
                           operator[](i)->pointer(), this->nBasis(),
                           this->nBasis(), this->nBasis());
      } else {
        switch (order()) {
        case 1:
          out << "1st order ";
          break;
        case 2:
          out << "2nd order ";
          break;
        case 3:
          out << "3rd order ";
          break;
        default:
          out << order() << "-th order ";
          break;
        }
        if (s == "")
          out << "vector integral";
        else
          out << "VectorInts[" + s + "]";
        out << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t order_bcast = order_;
        bool symmetric_bcast = symmetric_;
        MPIBCast(order_bcast,root,comm);
        MPIBCast(symmetric_bcast,root,comm);

        if (order_bcast != order_ or symmetric_bcast != symmetric_) {
          order_ = order_bcast;
          symmetric_ = symmetric_bcast;
          components_.clear();
          size_t size = nComponents();
          components_.reserve(size);
          for (size_t i = 0; i < size; i++) {
            components_.emplace_back(std::make_shared<OnePInts<IntsT>>(NB));
          }
        }

        for (std::shared_ptr<OnePInts<IntsT>>& comp : components_)
          comp->broadcast(comm, root);
      }
#endif
    }

    template <typename IntsU>
    VectorInts<IntsU> spatialToSpinBlock() const {
      VectorInts<IntsU> spinBlockInts(NB * 2, order_, symmetric_);
      size_t size = nComponents();
      for (size_t i = 0; i < size; i++) {
        *spinBlockInts.components_[i] = components_[i]->template spatialToSpinBlock<IntsU>();
      }
      return spinBlockInts;
    }

    template <typename TransT>
    VectorInts<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      VectorInts<typename std::conditional<
      (std::is_same<IntsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(NT, order(), symmetric());
      for (size_t i = 0; i < size(); i++){
        std::shared_ptr<OnePInts<IntsT>> comp = (*this)[i];
        
        if (std::shared_ptr<OnePRelInts<IntsT>> relInt
            = std::dynamic_pointer_cast<OnePRelInts<IntsT>>(comp)) {
          transInts[i] = std::make_shared<OnePRelInts<typename std::conditional<
            (std::is_same<IntsT, dcomplex>::value or std::is_same<TransT, dcomplex>::value),
            dcomplex, double>::type>>( relInt->transform(TRANS, T, NT, LDT), 0 );
            continue;
        }
        transInts[i] = std::make_shared<OnePInts<typename std::conditional<
            (std::is_same<IntsT, dcomplex>::value or std::is_same<TransT, dcomplex>::value),
            dcomplex, double>::type>>( comp->transform(TRANS, T, NT, LDT), 0 );
      }
      return transInts;
    }

    void convert2OnePRelInts(size_t nb, bool SORelativistic){
      size_t originalSize = this->size();
      components_.clear();
      components_.reserve(originalSize);
      for (size_t i = 0; i < originalSize; i++) 
        components_.emplace_back(std::make_shared<OnePRelInts<IntsT>>(nb, SORelativistic));
    }

    ~VectorInts() {}

  }; // class VectorInts

}; // namespace ChronusQ
