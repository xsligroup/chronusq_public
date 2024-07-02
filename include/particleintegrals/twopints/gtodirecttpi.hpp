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

#include <particleintegrals/twopints.hpp>

namespace ChronusQ {

  template <typename IntsT>
  class DirectTPI : public TwoPInts<IntsT> {

    template <typename IntsU>
    friend class DirectTPI;

  protected:
    BasisSet &basisSet_;
    BasisSet &basisSet2_;
    double threshSchwarz_ = 1e-12; ///< Schwarz screening threshold
    double* schwarz_ = nullptr;   ///< Schwarz bounds for the TPIs
    double* schwarz2_ = nullptr;  ///< second Schwarz bounds for the TPIs
    Molecule &molecule_;

  public:

    // Constructor
    DirectTPI() = delete;
    DirectTPI(BasisSet &basis, BasisSet &basis2, Molecule &mol, double threshSchwarz):
        TwoPInts<IntsT>(basis.nBasis, basis2.nBasis), 
        basisSet_(basis), basisSet2_(basis2),molecule_(mol),
        threshSchwarz_(threshSchwarz) {}
    DirectTPI( const DirectTPI &other ):
        DirectTPI(other.basisSet(), other.basisSet2(), other.molecule(), other.threshSchwarz_) {
      std::copy_n(other.schwarz_, basisSet_.nShell*basisSet_.nShell, schwarz_);
      std::copy_n(other.schwarz2_, basisSet2_.nShell*basisSet2_.nShell, schwarz2_);
    }
    template <typename IntsU>
    DirectTPI( const DirectTPI<IntsU> &other, int = 0 ):
        DirectTPI(other.basisSet_, other.basisSet2_, other.molecule_, other.threshSchwarz_) {
      if (other.schwarz_) {
        size_t NS = basisSet().nShell;
        schwarz_ = CQMemManager::get().malloc<double>(NS*NS);
        std::copy_n(other.schwarz_, NS*NS, schwarz_);
      }
      if (other.schwarz2_) {
        size_t NS = basisSet2().nShell;
        schwarz2_ = CQMemManager::get().malloc<double>(NS*NS);
        std::copy_n(other.schwarz2_, NS*NS, schwarz2_);
      }
    }
    DirectTPI( DirectTPI &&other ): TwoPInts<IntsT>(std::move(other)),
      basisSet_(other.basisSet_), basisSet2_(other.basisSet2_), 
      threshSchwarz_(other.threshSchwarz_),
      molecule_(other.molecule_),
      schwarz_(other.schwarz_), schwarz2_(other.schwarz2_) {
      other.schwarz_ = nullptr; 
      other.schwarz2_ = nullptr;
    }

    BasisSet& basisSet() { return basisSet_; }
    BasisSet& basisSet2() { return basisSet2_; }
    Molecule& molecule() { return molecule_; }
    const BasisSet& basisSet() const { return basisSet_; }
    const BasisSet& basisSet2() const { return basisSet2_; }
    const Molecule& molecule() const { return molecule_; }
    double threshSchwarz() const { return threshSchwarz_; }
    double*& schwarz()  { return schwarz_; }
    double*& schwarz2() { return schwarz2_; }

    // Single element interfaces
    virtual IntsT operator()(size_t p, size_t q, size_t r, size_t s) const override {
      CErr("NYI");
      return 0;
    }
    virtual IntsT operator()(size_t pq, size_t rs) const override {
      CErr("NYI");
      return 0;
    }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override {}

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override {}

    virtual void clear() override {}

    void computeSchwarz();

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (s == "")
        out << "  Two particle integral:" << std::endl;
      else
        out << "  TPI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "DIRECT";
      out << std::endl;

      out << "    * Schwarz Screening Threshold = "
          << threshSchwarz() << std::endl;

      if (printFull) {
        CErr("Printing Full ERI tensor for Direct contraction NYI.", out);
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      TwoPInts<IntsT>::broadcast(comm, root);
#ifdef CQ_ENABLE_MPI
      CErr("DirectTPI::broadcast() NYI");
#endif
    }

    virtual ~DirectTPI() {
      if(schwarz_)  CQMemManager::get().free(schwarz_);
      if(schwarz2_) CQMemManager::get().free(schwarz2_);
    }

  }; // class DirectTPI

  template <typename MatsT, typename IntsT>
  class GTODirectTPIContraction : public TPIContractions<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class GTODirectTPIContraction;

  public:

    // Constructors

    GTODirectTPIContraction() = delete;
    GTODirectTPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      TPIContractions<MatsT,IntsT>(tpi) {

      if (typeid(*tpi) != typeid(DirectTPI<IntsT>))
        CErr("GTODirectTPIContraction expect a DirectTPI reference.");

    }

    template <typename MatsU>
    GTODirectTPIContraction(
      const GTODirectTPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      GTODirectTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }
    template <typename MatsU>
    GTODirectTPIContraction(
      GTODirectTPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      GTODirectTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }

    GTODirectTPIContraction( const GTODirectTPIContraction &other ):
      GTODirectTPIContraction(other, 0) {
      this->contractSecond = other.contractSecond;
    }
    GTODirectTPIContraction( GTODirectTPIContraction &&other ):
      GTODirectTPIContraction(std::move(other), 0) {
      this->contractSecond = other.contractSecond;
    }

    /**
     *  \brief Perform various tensor contractions of the ERI tensor
     *  directly. Wraps other helper functions and provides
     *  loop structure
     *
     *  Currently supports
     *    - Coulomb-type (34,12) contractions
     *    - Exchange-type (23,12) contractions
     *
     *  Works with both real and complex matricies
     *
     *  \param [in/out] list Contains information pertinent to the
     *    matricies to be contracted with. See TwoBodyContraction
     *    for details
     */
    virtual void twoBodyContract(
        MPI_Comm c,
        const bool screen,
        std::vector<TwoBodyContraction<MatsT>> &list,
        EMPerturbation &pert) const {

        directScaffoldNew(c, screen, list);
    }

    void directScaffoldNew(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&) const;
    
    size_t directScaffoldNewSCRSize() const;

    void direct4CScaffold(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&) const;

    void directScaffold(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        EMPerturbation&) const;

    virtual ~GTODirectTPIContraction() {}

  }; // class GTODirectTPIContraction

}; // namespace ChronusQ
