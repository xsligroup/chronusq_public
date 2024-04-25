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
//#include <util/time.hpp>
#include <cxxapi/output.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  template <typename IntsT>
  class InCore4indexTPI : public TwoPInts<IntsT> {

    template <typename IntsU>
    friend class InCore4indexTPI;

  protected:
    size_t NB2, sNB2, NB3;
    IntsT* TPI = nullptr; ///< Two particle integrals (4 index)

  public:

    // Constructor
    InCore4indexTPI() = delete;
    InCore4indexTPI(size_t nb, size_t snb = 0):
        TwoPInts<IntsT>(nb, snb) {
      NB2  = this->nBasis()*this->nBasis();
      sNB2 = this->snBasis()*this->snBasis();
      NB3  = NB2 * this->snBasis();
      malloc();
    }
    InCore4indexTPI( const InCore4indexTPI &other ):
        InCore4indexTPI(other.nBasis(), other.snBasis()) {
      std::copy_n(other.TPI, NB2*sNB2, TPI);
    }
    template <typename IntsU>
    InCore4indexTPI( const InCore4indexTPI<IntsU> &other, int = 0 ):
        InCore4indexTPI(other.nBasis(), other.snBasis()) {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real InCore4indexTPI from a Complex one.");
      std::copy_n(other.TPI, NB2*sNB2, TPI);
    }
    InCore4indexTPI( InCore4indexTPI &&other ): TwoPInts<IntsT>(std::move(other)),
      NB2(other.NB2), sNB2(other.sNB2), NB3(other.NB3), TPI(other.TPI) {
      other.TPI = nullptr;
    }

    InCore4indexTPI& operator=( const InCore4indexTPI &other ) {
      if (this != &other) { // self-assignment check expected
        if (this->nBasis() != other.nBasis() or this->snBasis() != other.snBasis()) {
          this->NB = other.NB;
          this->sNB = other.sNB;
          NB2  = other.NB2;
          sNB2 = other.sNB2;
          NB3  = other.NB3;
          malloc(); // reallocate memory
        }
        std::copy_n(other.TPI, NB2*sNB2, TPI);
      }
      return *this;
    }
    InCore4indexTPI& operator=( InCore4indexTPI &&other ) {
      if (this != &other) { // self-assignment check expected
        CQMemManager::get().free(TPI);
        this->NB  = other.NB;
        this->sNB = other.sNB;
        NB2  = other.NB2;
        sNB2 = other.sNB2;
        NB3 = other.NB3;
        TPI = other.TPI;
        other.TPI = nullptr;
      }
      return *this;
    }

    // Single element interfaces
    virtual IntsT operator()(size_t p, size_t q, size_t r, size_t s) const {
      return TPI[p + q*this->nBasis() + r*NB2 + s*NB3];
    }
    IntsT& operator()(size_t p, size_t q, size_t r, size_t s) {
      return TPI[p + q*this->nBasis() + r*NB2 + s*NB3];
    }
    virtual IntsT operator()(size_t pq, size_t rs) const {
      return TPI[pq + rs*NB2];
    }
    IntsT& operator()(size_t pq, size_t rs) {
      return TPI[pq + rs*NB2];
    }

    // Tensor direct access
    IntsT* pointer() { return TPI; }
    const IntsT* pointer() const { return TPI; }

    // Computation interfaces
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
        EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &hamiltonianOptions) {

      // Use Libcint to compute nonrelativistic integrals
      if ( hamiltonianOptions.Libcint ) computeERIGCCINT(basisSet, mol, emPert, op, hamiltonianOptions);
      // Use Libint to compute nonrelativistic integrals
      else if (hamiltonianOptions.basisType == COMPLEX_GIAO)
        computeERINR(basisSet, basisSet, mol, emPert, op, hamiltonianOptions);
      else
        computeERIGCNR(basisSet, mol, emPert, op, hamiltonianOptions);

    }

    /// Evaluate nonrealtivistic ERIs in the CGTO basis
    void computeERINR(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);
    void computeERIGCNR(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);
    void computeERINRCINT(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);
    void computeERIGCCINT(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&);


    virtual void computeAOInts(BasisSet &basisSet, BasisSet &basisSet2, 
      Molecule &mol, EMPerturbation &emPert, OPERATOR op, 
      const HamiltonianOptions &hamiltonianOptions) { 

      // Use Libint to compute nonrelativistic integrals
      computeERINR(basisSet, basisSet2, mol, emPert, op, hamiltonianOptions);     

    };

    virtual void clear() {
      std::fill_n(TPI, NB2*sNB2, IntsT(0.));
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (s == "")
        out << "  Two particle integral:" << std::endl;
      else
        out << "  TPI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "INCORE (Gemm)";
      out << std::endl;
      if (printFull) {
        out << bannerTop << std::endl;
        size_t NB = this->nBasis();
        size_t sNB = this->snBasis();
        out << std::scientific << std::left << std::setprecision(8);
        for(auto i = 0ul; i < NB; i++)
        for(auto j = 0ul; j < NB; j++)
        for(auto k = 0ul; k < sNB; k++)
        for(auto l = 0ul; l < sNB; l++){
          if (std::abs(operator()(i,j,k,l)) > PRINT_SMALL) {
            out << "    (" << i << "," << j << "|" << k << "," << l << ")  ";
            out << operator()(i,j,k,l) << std::endl;
          }
        };
        out << bannerEnd << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {
      TwoPInts<IntsT>::broadcast(comm, root);

#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t NB2_bcast = NB2;
        size_t sNB2_bcast = sNB2;
        MPIBCast(NB2_bcast,root,comm);
        MPIBCast(sNB2_bcast,root,comm);

        if (NB2_bcast != NB2 or sNB2_bcast != sNB2) {
          NB2  = this->nBasis()*this->nBasis();
          sNB2 = this->snBasis()*this->snBasis();
          NB3  = NB2 * this->snBasis();
          malloc();
        }

        MPIBCast(TPI,NB2*sNB2,root,comm);
      }
#endif

    }

    void malloc() {
      if(TPI) CQMemManager::get().free(TPI);
      size_t NB4 = NB2*sNB2;
      try { TPI = CQMemManager::get().template malloc<IntsT>(NB4); }
      catch(...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for the full TPI tensor ("
                  << (NB4/1e9) * sizeof(double) << " GB)" << std::endl;
        std::cout << std::endl << CQMemManager::get() << std::endl;
        CErr();
      }
    }

    template <typename IntsU>
    InCore4indexTPI<IntsU> spatialToSpinBlock(char TRANS1 = 'I', char TRANS2 = 'I') const;

    template <typename TransT>
    InCore4indexTPI<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const;

    template <typename TransT, typename OutT>
    void subsetTransform(
        char TRANS, const TransT* T, int LDT,
        const std::vector<std::pair<size_t,size_t>> &off_size,
        OutT* out, bool increment = false) const;

    virtual ~InCore4indexTPI() {
      if(TPI) CQMemManager::get().free(TPI);
    }

  }; // class InCore4indexTPI

  template <typename MatsT, typename IntsT>
  class InCore4indexTPIContraction : public TPIContractions<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCore4indexTPIContraction;

  public:

    // Constructors

    InCore4indexTPIContraction() = delete;
    InCore4indexTPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      TPIContractions<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCore4indexTPIContraction(
      const InCore4indexTPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCore4indexTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }
    template <typename MatsU>
    InCore4indexTPIContraction(
      InCore4indexTPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCore4indexTPIContraction(other.ints_) {
      this->contractSecond = other.contractSecond;
    }

    InCore4indexTPIContraction( const InCore4indexTPIContraction &other ):
      InCore4indexTPIContraction(other, 0) {
      this->contractSecond = other.contractSecond;
    }
    InCore4indexTPIContraction( InCore4indexTPIContraction &&other ):
      InCore4indexTPIContraction(std::move(other), 0) {
      this->contractSecond = other.contractSecond;
    }

    // Computation interfaces
    virtual void twoBodyContract(
        MPI_Comm comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        EMPerturbation&) const;

    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual ~InCore4indexTPIContraction() {}

  }; // class InCore4indexTPIContraction

}; // namespace ChronusQ
