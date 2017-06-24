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
#include <particleintegrals/twopints.hpp>

namespace ChronusQ {

  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  gradient integral matrices in a finite basis set.
   *
   *  Templated over IntClass to be general to the underlying storage
   */
  template <template <typename> class IntClass, typename IntsT>
  class GradInts : public ParticleIntegrals {

    template <template <typename> class, typename>
    friend class GradInts;

  protected:
    typedef std::shared_ptr<IntClass<IntsT>> IntPtr;
    typedef std::shared_ptr<const IntClass<IntsT>> ConstIntPtr;

    std::vector<IntPtr> components_;
    size_t nAtoms_;

  public:

    //
    // Constructors
    //
    
    GradInts() = delete;
    GradInts( const GradInts & ) = default;
    GradInts( GradInts && ) = default;

    // Main constructor
    GradInts(size_t nBasis, size_t nAtoms):
        ParticleIntegrals(nBasis), nAtoms_(nAtoms) {

      components_.reserve(3*nAtoms_);

      for (size_t i = 0; i < 3*nAtoms_; i++) {
        components_.emplace_back(
          std::make_shared<IntClass<IntsT>>(nBasis)
        );
      }

    };

    // Two basis constructor
    GradInts(size_t nBasis, size_t snBasis, size_t nAtoms):
        ParticleIntegrals(nBasis), nAtoms_(nAtoms) {

      components_.reserve(3*nAtoms_);

      for (size_t i = 0; i < 3*nAtoms_; i++) {
        components_.emplace_back(
          std::make_shared<IntClass<IntsT>>(nBasis, snBasis)
        );
      }

    };

    // Constructor from a vector of pointers to integrals
    template <template <typename> class IntSubClass>
    GradInts(size_t nBasis, size_t nAtoms,
      std::vector<std::shared_ptr<IntSubClass<IntsT>>> integrals) :
      ParticleIntegrals(nBasis), nAtoms_(nAtoms)
    {

      components_.reserve(3*nAtoms_);

      assert(integrals.size() == 3*nAtoms_);

      for (size_t i = 0; i < 3*nAtoms_; i++) {
        components_.emplace_back(
          std::dynamic_pointer_cast<IntClass<IntsT>>(integrals[i])
        );
      }

      if ( std::any_of(components_.begin(),
                       components_.end(),
                       [](auto& p) { return p == nullptr; }) ) {
        CErr("Couldn't convert to parent class in GradInts constructor!");
      }

    }

    // Integral matrix type converter constructor
    template <typename IntsU>
    GradInts( const GradInts<IntClass,IntsU> &other, int = 0 ):
        ParticleIntegrals(other), nAtoms_(other.nAtoms_) {

      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real GradInts from a Complex one.");

      components_.reserve(other.components_.size());

      for (auto &p : other.components_)
        components_.emplace_back(std::make_shared(*p));
    }

    // TODO
    // // Integral class converter constructor
    // template <template <typename> class OtherIntClass>
    // GradInts( const GradInts<OtherIntClass,IntsT> &other, char = ' ' ) :
    //   ParticleIntegrals(other), nAtoms_(other.nAtoms_) {

    // }


    //
    // Access
    //

    // Element access by internal storage
    IntPtr& operator[](size_t i) {
      return components_[i];
    }
    const ConstIntPtr& operator[](size_t i) const {
      return components_[i];
    }

    // Element access by atom index and cartesian index
    IntPtr& integralByAtomCart(size_t atom, size_t xyz) {
      return components_[3*atom + xyz];
    }
    const ConstIntPtr& integralByAtomCart(size_t atom, size_t xyz) const {
      return components_[3*atom + xyz];
    }

    // Size access
    size_t nAtoms() const { return nAtoms_; }
    size_t size() const { return components_.size(); }


    //
    // Interface methods
    //

    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);

    virtual void clear() {
      for (IntPtr& c : components_)
        c->clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {

      if (printFull) {

        std::string oeiStr;
        if (s == "")
          oeiStr = "GradInts.";
        else
          oeiStr = "GradInts[" + s + "].";

        std::cout << oeiStr << std::endl;

        for (size_t i = 0; i < components_.size(); i++) {
          oeiStr = "Grad " + s + " " + std::to_string(i);
          components_[i]->output(out, oeiStr, printFull);
        }

      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Gradient Integral";
        else
          oeiStr = "GradInts[" + s + "]";
        out << oeiStr << " with " << nAtoms_ << " atoms.";
        out << std::endl;
      }

    }

    // Conversion to other type
    template <typename TransT>
    GradInts<IntClass, typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {

      GradInts<IntClass, typename std::conditional<
      (std::is_same<IntsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(NT, nAtoms_);

      for (size_t i = 0; i < components_.size(); i++)
        transInts[i] = (*this)[i].transform(TRANS, T, NT, LDT);

      return transInts;
    }

    ~GradInts() {}

  }; // class GradInts



  /**
   *  \brief Class to handle the contraction of the gradient two-body potential
   */
  template <typename MatsT, typename IntsT>
  class GradContractions {

    template <typename MatsU, typename IntsU>
    friend class GradContractions;

  protected:
    GradInts<TwoPInts,IntsT>& grad_;

  public:

    // Constructors
    GradContractions() = delete;
    GradContractions(GradInts<TwoPInts,IntsT>& grads): grad_(grads) { }

    template <typename MatsU>
    GradContractions( const GradContractions<MatsU,IntsT> &other, int dummy = 0 ):
      GradContractions(other.grad_) {}
    template <typename MatsU>
    GradContractions( GradContractions<MatsU,IntsT> &&other, int dummy = 0 ):
      GradContractions(other.grad_) {}

    GradContractions( const GradContractions &other ):
      GradContractions(other, 0) {}
    GradContractions( GradContractions &&other ):
      GradContractions(std::move(other), 0) {}

    GradInts<TwoPInts,IntsT>& ints() { return grad_; }
    const GradInts<TwoPInts,IntsT>& ints() const { return grad_; }


    /**
     *  Contract the gradient of the two-body potential with one-body operators
     *
     *  The outer vector of contractions is over gradient components (i.e. it
     *  should be 3*NAtoms long) and the inner vector is over contractions 
     *  (similar to in the non-gradient contraction case)
     *
     *  \param [in] comm MPI communicator for this parallel operation (usually the MPI world)
     *  \param [in] screen Whether or not to screen the integrals
     *  \param [in/out] contList List of one body operators for contraction.
     *  \param [in] pert Perturbation used in GIAO calculations, but GIAOs are NYI in this method
     */
    virtual void gradTwoBodyContract(
      MPI_Comm,
      const bool,
      std::vector<std::vector<TwoBodyContraction<MatsT>>>&,
      EMPerturbation&) const = 0;

    // Contract all into the same storage
    inline void gradTwoBodyContract(
      MPI_Comm comm,
      const bool screen,
      std::vector<TwoBodyContraction<MatsT>>& list,
      EMPerturbation& pert) const {
      
      // Copy into individual contractions (using the same storage)
      size_t N = this->grad_.size();
      std::vector<std::vector<TwoBodyContraction<MatsT>>> longList(N, list);

      gradTwoBodyContract(comm, screen, longList, pert);
    }

    // Other defaults
    inline void gradTwoBodyContract(
      MPI_Comm comm,
      const bool screen,
      std::vector<std::vector<TwoBodyContraction<MatsT>>>& contList) const {
      
      EMPerturbation pert;
      gradTwoBodyContract(comm,screen,contList,pert);
    }

    inline void gradTwoBodyContract(
      MPI_Comm comm,
      const bool screen,
      std::vector<TwoBodyContraction<MatsT>>& contList) const {
      
      EMPerturbation pert;
      gradTwoBodyContract(comm,screen,contList,pert);
    }
    
    inline void gradTwoBodyContract(
        MPI_Comm comm,
        std::vector<std::vector<TwoBodyContraction<MatsT>>> &contList,
        EMPerturbation &pert) const {
      gradTwoBodyContract(comm,true,contList,pert);
    }

    inline void gradTwoBodyContract(
        MPI_Comm comm,
        std::vector<TwoBodyContraction<MatsT>> &contList,
        EMPerturbation &pert) const {
      gradTwoBodyContract(comm,true,contList,pert);
    }

    inline void gradTwoBodyContract(
        MPI_Comm comm,
        std::vector<std::vector<TwoBodyContraction<MatsT>>> &contList) const {
      gradTwoBodyContract(comm,true,contList);
    }

    inline void gradTwoBodyContract(
        MPI_Comm comm,
        std::vector<TwoBodyContraction<MatsT>> &contList) const {
      gradTwoBodyContract(comm,true,contList);
    }

    // Destructor
    virtual ~GradContractions() {}

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<GradContractions<MatsU,IntsT>>
    convert(const std::shared_ptr<GradContractions<MatsT,IntsT>>&);

    bool contractSecond = false;

  }; // class GradientContractions

}; // namespace ChronusQ

