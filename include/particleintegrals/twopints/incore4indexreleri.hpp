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
#include <particleintegrals/onepints/relativisticints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

namespace ChronusQ {

  // define a macro to handle ERI and 4C ERI subsetTransfrom Calls
  #define INCORE4CINDEXERI_SUBSETTRANSFROM(ERI, ERI4C, ...) \
    if (ERI4C) ERI4C->subsetTransform(__VA_ARGS__); \
    else ERI->subsetTransform(__VA_ARGS__);
  
  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  electron repulsion integral matrices V, pVp, and three pVxp in
   *  a finite basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class InCore4indexRelERI : public InCore4indexTPI<IntsT> {

    template <typename IntsU>
    friend class InCore4indexRelERI;

  protected:
    std::vector<InCore4indexTPI<IntsT>> components_;

  public:

    // Constructor
    InCore4indexRelERI() = delete;
    InCore4indexRelERI( const InCore4indexRelERI & ) = default;
    InCore4indexRelERI( InCore4indexRelERI && ) = default;
    InCore4indexRelERI(size_t nb, size_t nRel):
        InCore4indexTPI<IntsT>(nb) {
      components_.reserve(nRel);
      for (size_t i = 0; i < nRel; i++)
        components_.emplace_back(nb);
    }

    template <typename IntsU>
    InCore4indexRelERI( const InCore4indexRelERI<IntsU> &other, int = 0 ):
        InCore4indexTPI<IntsT>(other, 0) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real InCore4indexRelERI from a Complex one.");
      components_.reserve(other.components_.size());
      for (auto &p : other.components_)
        components_.emplace_back(p);
    }

    size_t nRelComp() const { return components_.size(); }

    InCore4indexTPI<IntsT>& operator[](size_t i) {
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return components_[i];
    }

    InCore4indexTPI<IntsT>& operator[](REL_INTS_COMPS comp) {
      if (comp == REL_INTS_COMPS::O)
        return *this;
      return operator[](static_cast<size_t>(comp));
    }

    const InCore4indexTPI<IntsT>& operator[](size_t i) const {
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return components_[i];
    }

    const InCore4indexTPI<IntsT>& operator[](REL_INTS_COMPS comp) const {
      if (comp == REL_INTS_COMPS::O)
        return *this;
      return operator[](static_cast<size_t>(comp));
    }

    std::vector<InCore4indexTPI<IntsT>>& relComps() {
      return components_;
    }
    const std::vector<InCore4indexTPI<IntsT>>& relComps() const {
      return components_;
    }

    
    // Computation interfaces
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
        EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &hamiltonianOptions) override {

      if(hamiltonianOptions.Libcint) {

        // Use Libcint to compute nonrelativistic and DCB integrals
        computeERICINT(basisSet, mol, emPert, op, hamiltonianOptions);

        //if( hamiltonianOptions.Gauge )
          //computeERIGauge(basisSet, mol, emPert, op, hamiltonianOptions);

      } else {

        // Use Libint to compute nonrelativistic
        InCore4indexTPI<IntsT>::computeAOInts(basisSet, mol, emPert, op, hamiltonianOptions);

        // Use Libint to compute DCB integrals
        if (hamiltonianOptions.DiracCoulomb or hamiltonianOptions.Gaunt)
          computeERIDCB(basisSet, mol, emPert, op, hamiltonianOptions);

        // use in house code to compute gauge integral
        if (hamiltonianOptions.Gauge)
          computeERIGauge(basisSet, mol, emPert, op, hamiltonianOptions);
      }

    }

    /// Evaluate Spin-Own-Orbit ERIs in the CGTO basis
    void computeERIDCB(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);
    void computeERIGauge(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);
    void computeERICINT(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);

    virtual void clear() override {
      InCore4indexTPI<IntsT>::clear();
      for (InCore4indexTPI<IntsT>& c : components_)
        c.clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "RelIncore4IndexERI";
        else
          oeiStr = "RelIncore4IndexERI[" + s + "]";
        prettyPrintSmart(out, oeiStr+".NR", this->pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        for (size_t i = 0; i < nRelComp(); i++) {
          prettyPrintSmart(out, oeiStr+"."+std::to_string(i), (*this)[i].pointer(),
              this->nBasis(), this->nBasis(), this->nBasis());
        }
      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Relativistic incore 4-index electron repulsion integral";
        else
          oeiStr = "RelIncore4IndexERI[" + s + "]";
        out << oeiStr;
        out << " with " << nRelComp();
        out << " DCB integrals" << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      InCore4indexTPI<IntsT>::broadcast(comm, root);

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

        for (InCore4indexTPI<IntsT>& comp : components_)
          comp.broadcast(comm, root);
      }
#endif
    }

    // Pauil Matrice representation to spinor representation
    template <typename IntsU>
    InCore4indexRelERI<IntsU> spatialToSpinBlock() const;
    
    template <typename TransT>
    InCore4indexRelERI<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      InCore4indexRelERI<typename std::conditional<
      (std::is_same<IntsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(NT, nRelComp());
      transInts[REL_INTS_COMPS::O] =
          (*this)[REL_INTS_COMPS::O].transform(TRANS, T, NT, LDT);
      for (size_t i = 0; i < nRelComp(); i++)
        transInts[i] = (*this)[i].transform(TRANS, T, NT, LDT);
      return transInts;
    }

    template <typename TransT, typename OutT>
    void subsetTransform(
        char TRANS, const TransT* T, int LDT,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        OutT* out, bool increment = false) const;
    
    template <typename TransT, typename OutT>
    void subsetTransformWithLSComps(
        const std::string & LSComps, char TRANS, 
        const TransT* TL, int LDTL, const TransT* TS, int LDTS,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        const IntsT * in, OutT* out, bool increment = false) const;
    
    virtual ~InCore4indexRelERI() {}

  }; // class InCore4indexRelERI

  template <typename MatsT, typename IntsT>
  class InCore4indexRelERIContraction : public InCore4indexTPIContraction<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCore4indexRelERIContraction;

  public:

    // Constructors

    InCore4indexRelERIContraction() = delete;
    InCore4indexRelERIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      InCore4indexTPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCore4indexRelERIContraction(
        const InCore4indexRelERIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCore4indexRelERIContraction(other.ints_) {}
    template <typename MatsU>
    InCore4indexRelERIContraction(
        InCore4indexRelERIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCore4indexRelERIContraction(other.ints_) {}

    InCore4indexRelERIContraction( const InCore4indexRelERIContraction &other ):
      InCore4indexRelERIContraction(other, 0) {}
    InCore4indexRelERIContraction( InCore4indexRelERIContraction &&other ):
      InCore4indexRelERIContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual ~InCore4indexRelERIContraction() {}

  }; // class InCore4indexRelERIContraction

}; // namespace ChronusQ
