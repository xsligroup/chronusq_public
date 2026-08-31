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
#include <integrals.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>

namespace ChronusQ {

  enum class REL_INTS_COMPS : size_t {
    SCALAR = 0, SOZ = 1, SOY = 2, SOX = 3, O = 4
  };

  template <typename IntsT>
  struct TPIContractionPointers {
    std::vector<IntsT*> pointers;
    std::vector<size_t> aux_dims;

    TPIContractionPointers() {
      pointers.reserve(2);
      aux_dims.reserve(1);
    }

    TPIContractionPointers(IntsT *tpi): pointers({tpi}) {}

    TPIContractionPointers(IntsT *tpi_1, IntsT *tpi_2, size_t aux_dim)
    : pointers({tpi_1, tpi_2}), aux_dims({aux_dim}) {}

    bool isRI() const {
      return pointers.size() > 1;
    }
  };

  template <typename IntsT>
  class IncoreTPIList: public ParticleIntegrals {
  public:

    IncoreTPIList() = delete;
    IncoreTPIList(size_t nb):
      ParticleIntegrals(nb) {}

    /**
     * Return pointers for two-body contraction
     * @param index the index of requested integral
     * @return pointers for two-body contraction
     */
    virtual TPIContractionPointers<IntsT> getPointers(size_t index) = 0;

    /// Evaluate AO Integrals according to one basis sets
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) = 0;
    /// Evaluate AO Integrals according to two basis sets
    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
                               OPERATOR, const HamiltonianOptions&) override {

      CErr("Relativistic TPI with two sets of basis NYI");

    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      ParticleIntegrals::broadcast(comm, root);
    }
  };

  template <typename IntsT>
  class Incore4indexTPIList : public IncoreTPIList<IntsT> {
  protected:
    std::vector<std::shared_ptr<InCore4indexTPI<IntsT>>> components_;

  public:
    Incore4indexTPIList(size_t nb, size_t nComp):
      IncoreTPIList<IntsT>(nb) {
      components_.reserve(nComp);
      for (size_t i = 0; i < nComp; i++)
        components_.emplace_back(std::make_shared<InCore4indexTPI<IntsT>>(nb));
    }

    virtual TPIContractionPointers<IntsT> getPointers(size_t index) override {
      return components_[index]->pointer();
    }

    InCore4indexTPI<IntsT>& operator[](size_t i) {
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return *components_[i];
    }

    const InCore4indexTPI<IntsT>& operator[](size_t i) const {
      if (i >= components_.size())
        CErr("Requested component is NOT in this RelativisticInts object.");
      return *components_[i];
    }

    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) override {

      if(hamiltonianOptions.Libcint) {

        // Use Libcint to compute DCB integrals
        computeERICINT(basisSet, mol, emPert, op, hamiltonianOptions);

      } else {

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
      for (auto c : components_)
        c->clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Incore4indexTPIList";
        else
          oeiStr = s;
        for (size_t i = 0; i < components_.size(); i++) {
          prettyPrintSmart(out, oeiStr+"["+std::to_string(i)+"]", (*this)[i].pointer(),
                           this->nBasis(), this->nBasis(), this->nBasis());
        }
      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "List of incore 4-index two particle integral";
        else
          oeiStr = "Incore4indexTPIList[" + s + "]";
        out << oeiStr;
        out << " with " << components_.size();
        out << " DCB integrals" << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      IncoreTPIList<IntsT>::broadcast(comm, root);
#ifdef CQ_ENABLE_MPI
      if( MPISize(comm) > 1 ) {
        size_t nRel = components_.size();
        MPIBCast(nRel,root,comm);

        if (components_.size() != nRel) {
          components_.clear();
          components_.reserve(nRel);
          for (size_t i = 0; i < nRel; i++) {
            components_.emplace_back(std::make_shared<InCore4indexTPI<IntsT>>(this->NB));
          }
        }

        for (auto comp : components_)
          comp->broadcast(comm, root);
      }
#endif
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
  };

  template <typename IntsT>
  class IncoreDCRITPIList : public IncoreTPIList<IntsT> {
  protected:
    std::shared_ptr<InCoreCholeskyRIERI<IntsT>> LLLL_;
    std::vector<std::shared_ptr<InCoreAsymmRITPI<IntsT>>> components_;
  public:

    IncoreDCRITPIList() = delete;
    IncoreDCRITPIList(size_t nb, std::shared_ptr<InCoreCholeskyRIERI<IntsT>> LLLL):
      IncoreTPIList<IntsT>(nb), LLLL_(LLLL), components_(4, nullptr) {}

    std::shared_ptr<InCoreCholeskyRIERI<IntsT>> LLLL_term() const {
      return LLLL_;
    }

    void set_LLLL_term(std::shared_ptr<InCoreCholeskyRIERI<IntsT>> LLLL) {
      LLLL_ = LLLL;
    }

    std::shared_ptr<InCoreAsymmRITPI<IntsT>> asymm_term(size_t index) const {
      return components_[index];
    }

    void set_asymm_term(size_t index, std::shared_ptr<InCoreAsymmRITPI<IntsT>> asymm_ints) {
      components_[index] = asymm_ints;
    }

    size_t nRIBasis() const {
      return LLLL_->nRIBasis();
    }

    virtual TPIContractionPointers<IntsT> getPointers(size_t index) override {
      return TPIContractionPointers<IntsT>(components_[index]->pointer(), LLLL_->pointer(), LLLL_->nRIBasis());
    }

    /// Evaluate AO Integrals according to one basis sets
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) override {

      if (not hamiltonianOptions.DiracCoulomb and not hamiltonianOptions.DiracCoulombSSSS)
        CErr("DC terms not requested in hamiltonianOptions");

      if (hamiltonianOptions.Libcint) {
        if (LLLL_->nRIBasis() == 0) // Missing LLLL riERI for DC riERI terms
          LLLL_->computeAOInts(basisSet, mol, emPert, op, hamiltonianOptions);

        computeCholeskyDCERI_CINT(basisSet, mol);
      } else
        CErr("Libint version of CholeskyDCERI NYI");
    }

    void computeCholeskyDCERI_CINT(BasisSet &originalBasisSet, Molecule &mol);

    virtual void clear() override {
      LLLL_ = nullptr;
      for (auto &c : components_)
        c = nullptr;
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      CErr("NYI");
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      IncoreTPIList<IntsT>::broadcast(comm, root);
#ifdef CQ_ENABLE_MPI
      CErr("IncoreDCRITPIList::broadcast NYI");
      if( MPISize(comm) > 1 ) {
        LLLL_->broadcast(comm, root);

        size_t nRel = components_.size();
        MPIBCast(nRel,root,comm);

        if (components_.size() != nRel) {
          components_.clear();
          components_.reserve(nRel);
          for (size_t i = 0; i < nRel; i++) {
//            components_.emplace_back(std::make_shared<InCore4indexTPI<IntsT>>(this->memManager_, this->NB));
          }
        }

        for (auto comp : components_)
          comp->broadcast(comm, root);
      }
#endif
    }
  };

  template <typename IntsT>
  class IncoreSSSSRITPIList : public IncoreTPIList<IntsT> {
  protected:
    std::shared_ptr<IncoreDCRITPIList<IntsT>> DC_;
  public:

    IncoreSSSSRITPIList() = delete;
    IncoreSSSSRITPIList(std::shared_ptr<IncoreDCRITPIList<IntsT>> DC):
        IncoreTPIList<IntsT>(DC->nBasis()), DC_(DC) {}

    std::shared_ptr<IncoreDCRITPIList<IntsT>> DC_term() const {
      return DC_;
    }

    void set_DC_term(std::shared_ptr<IncoreDCRITPIList<IntsT>> DC) {
      DC_ = DC;
    }

    size_t nRIBasis() const {
      return DC_->nRIBasis();
    }

    virtual TPIContractionPointers<IntsT> getPointers(size_t index) override {
      IntsT *leftPtr = nullptr, *rightPtr = nullptr;

      if (index < 4) {
        leftPtr = DC_->asymm_term(index)->pointer();
        rightPtr = DC_->asymm_term(0)->pointer();

      } else if (index < 16) {
        index -= 4;
        leftPtr = DC_->asymm_term(index / 3)->pointer();
        rightPtr = DC_->asymm_term(index % 3 + 1)->pointer();

      } else
        CErr("Invalid index " + std::to_string(index) + " for IncoreSSSSRITPIList");

      return TPIContractionPointers<IntsT>(leftPtr, rightPtr, DC_->nRIBasis());
    }

    /// Evaluate AO Integrals according to one basis sets
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) override {

      if (not hamiltonianOptions.DiracCoulombSSSS)
        CErr("DC-SSSS terms not requested in hamiltonianOptions");

      if (DC_->asymm_term(0) == nullptr) // ri-SSLL term not built yet
        DC_->computeAOInts(basisSet, mol, emPert, op, hamiltonianOptions);

    }

    virtual void clear() override {
      DC_ = nullptr;
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      CErr("NYI");
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      IncoreTPIList<IntsT>::broadcast(comm, root);
      DC_->broadcast(comm, root);
    }

  };

  template <typename IntsT>
  class IncoreGauntRITPIList : public IncoreTPIList<IntsT> {
  protected:
    std::shared_ptr<InCoreRITPI<IntsT>> LLLL_;
    std::vector<std::shared_ptr<InCoreAsymmRITPI<IntsT>>> components_;
  public:

    IncoreGauntRITPIList() = delete;
    IncoreGauntRITPIList(size_t nb, std::shared_ptr<InCoreRITPI<IntsT>> LLLL):
        IncoreTPIList<IntsT>(nb), LLLL_(LLLL), components_(4, nullptr) {}

    std::shared_ptr<InCoreRITPI<IntsT>> LLLL_term() const {
      return LLLL_;
    }

    void set_LLLL_term(std::shared_ptr<InCoreRITPI<IntsT>> LLLL) {
      LLLL_ = LLLL;
    }

    std::shared_ptr<InCoreAsymmRITPI<IntsT>> asymm_term(size_t index) const {
      return components_[index];
    }

    void set_asymm_term(size_t index, std::shared_ptr<InCoreAsymmRITPI<IntsT>> asymm_ints) {
      components_[index] = asymm_ints;
    }

    size_t nRIBasis() const {
      return LLLL_->nRIBasis();
    }

    virtual TPIContractionPointers<IntsT> getPointers(size_t index) override {
      CErr("IncoreGauntRITPIList::getPointers NYI");
      return TPIContractionPointers<IntsT>();
    }

    /// Evaluate AO Integrals according to one basis sets
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
                               EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) override {

      if (not hamiltonianOptions.Gaunt)
        CErr("Gaunt terms not requested in hamiltonianOptions");

      if (hamiltonianOptions.Libcint)
        computeCholeskyGauntERI_CINT(basisSet, mol);
      else
        CErr("Libint version of CholeskyDCERI NYI");
    }

    void computeCholeskyGauntERI_CINT(BasisSet &originalBasisSet, Molecule &mol);

    virtual void clear() override {
      for (auto &c : components_)
        c = nullptr;
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const override {
      CErr("NYI");
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      IncoreTPIList<IntsT>::broadcast(comm, root);
#ifdef CQ_ENABLE_MPI
      CErr("IncoreGauntRITPIList::broadcast NYI");
      if( MPISize(comm) > 1 ) {
        LLLL_->broadcast(comm, root);

        size_t nRel = components_.size();
        MPIBCast(nRel,root,comm);

        if (components_.size() != nRel) {
          components_.clear();
          components_.reserve(nRel);
          for (size_t i = 0; i < nRel; i++) {
//            components_.emplace_back(std::make_shared<InCore4indexTPI<IntsT>>(this->memManager_, this->NB));
          }
        }

        for (auto comp : components_)
          comp->broadcast(comm, root);
      }
#endif
    }
  };

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
  class InCoreRelERI : public InCoreTPI<IntsT> {

    template <typename IntsU>
    friend class InCoreRelERI;

  protected:
    std::shared_ptr<InCoreTPI<IntsT>> LLLL_ = nullptr;
    std::shared_ptr<IncoreTPIList<IntsT>> DC_ = nullptr; // 4
    std::shared_ptr<IncoreTPIList<IntsT>> gaunt_ = nullptr; // 19
    std::shared_ptr<IncoreTPIList<IntsT>> SSSS_ = nullptr; // 16
    std::shared_ptr<IncoreTPIList<IntsT>> gauge_ = nullptr; // 26

  public:

    // Constructor
    InCoreRelERI() = delete;
    InCoreRelERI( const InCoreRelERI & ) = default;
    InCoreRelERI( InCoreRelERI && ) = default;

    // The constructor only initialize super class
    InCoreRelERI(size_t nb):
      InCoreTPI<IntsT>(nb) {}

    // The constructor for incore 4-index version
    InCoreRelERI(size_t nb, bool DC, bool gaunt, bool SSSS, bool gauge):
      InCoreTPI<IntsT>(nb),
      LLLL_(std::make_shared<InCore4indexTPI<IntsT>>(nb)),
      DC_(DC? std::make_shared<Incore4indexTPIList<IntsT>>(nb, 4) : nullptr),
      gaunt_(gaunt? std::make_shared<Incore4indexTPIList<IntsT>>(nb, 19) : nullptr),
      SSSS_(SSSS? std::make_shared<Incore4indexTPIList<IntsT>>(nb, 16) : nullptr),
      gauge_(gauge? std::make_shared<Incore4indexTPIList<IntsT>>(nb, 26) : nullptr) {}

    // The constructor for incore RI version
    InCoreRelERI(std::shared_ptr<InCoreCholeskyRIERI<IntsT>> LLLL_riERI,
                       const HamiltonianOptions &hamiltonianOptions,
                       const CDRIIntsOptions& riOptions):
        InCoreTPI<IntsT>(LLLL_riERI->nBasis()), LLLL_(LLLL_riERI) {

      size_t NB = LLLL_riERI->nBasis();
      
      if (riOptions.CDRI_LLLL) {
        set_LLLL_term(LLLL_riERI);
      } else {
        set_LLLL_term(std::make_shared<InCore4indexTPI<IntsT>>(NB));
      }

      if (riOptions.CDRI_SSLL or riOptions.CDRI_SSSS) {
        std::shared_ptr<IncoreDCRITPIList<IntsT>> SSLL_riERI = std::make_shared<IncoreDCRITPIList<IntsT>>(NB, LLLL_riERI);

        if (hamiltonianOptions.DiracCoulomb)
          if (riOptions.CDRI_SSLL) {
            set_DC_terms(SSLL_riERI);
          } else {
            set_DC_terms(std::make_shared<Incore4indexTPIList<IntsT>>(NB, 4));
          }

        if (hamiltonianOptions.DiracCoulombSSSS)
          if (riOptions.CDRI_SSSS){
            set_SSSS_terms(std::make_shared<IncoreSSSSRITPIList<IntsT>>(SSLL_riERI));
          } else {
            set_SSSS_terms(std::make_shared<Incore4indexTPIList<IntsT>>(NB, 16));
          }
      }

      if (hamiltonianOptions.Gaunt)
        set_gaunt_terms(std::make_shared<Incore4indexTPIList<IntsT>>(NB, 19));

      if (hamiltonianOptions.Gauge)
        set_gauge_terms(std::make_shared<Incore4indexTPIList<IntsT>>(NB, 26));
    }

    size_t nRelComp() const { return (DC_ ? 4 : 0) + (gaunt_ ? 19 : 0)
                                + (SSSS_ ? 16 : 0) + (gauge_ ? 26 : 0); }

    std::shared_ptr<InCoreTPI<IntsT>> LLLL_term() const {
      return LLLL_;
    }

    void set_LLLL_term(std::shared_ptr<InCoreTPI<IntsT>> LLLL) {
      LLLL_ = LLLL;
    }

    std::shared_ptr<IncoreTPIList<IntsT>> DC_terms() const {
      return DC_;
    }

    void set_DC_terms(std::shared_ptr<IncoreTPIList<IntsT>> DC) {
      DC_ = DC;
    }

    std::shared_ptr<IncoreTPIList<IntsT>> gaunt_terms() const {
      return gaunt_;
    }

    void set_gaunt_terms(std::shared_ptr<IncoreTPIList<IntsT>> gaunt) {
      gaunt_ = gaunt;
    }

    std::shared_ptr<IncoreTPIList<IntsT>> SSSS_terms() const {
      return SSSS_;
    }

    void set_SSSS_terms(std::shared_ptr<IncoreTPIList<IntsT>> SSSS) {
      SSSS_ = SSSS;
    }

    std::shared_ptr<IncoreTPIList<IntsT>> gauge_terms() const {
      return gauge_;
    }

    void set_gauge_terms(std::shared_ptr<IncoreTPIList<IntsT>> gauge) {
      gauge_ = gauge;
    }

    // Single element interfaces
    virtual IntsT operator()(size_t p, size_t q, size_t r, size_t s) const {
      return LLLL_->operator()(p,q,r,s);
    }
    virtual IntsT operator()(size_t pq, size_t rs) const {
      return LLLL_->operator()(pq,rs);
    }

    // Tensor direct access
    virtual IntsT* pointer() override { return LLLL_->pointer(); }
    virtual const IntsT* pointer() const override { return LLLL_->pointer(); }

    // Two-body contraction helper
    TPIContractionPointers<IntsT> getPointers(TWOBODY_CONTRACTION_TYPE contType, size_t index) {
      TPIContractionPointers<IntsT> ERI4s;

      switch (contType) {
        case COULOMB:
        case EXCHANGE:
          if (std::shared_ptr<InCore4indexTPI<IntsT>> LLLL4I
              = std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(LLLL_)) {
            ERI4s.pointers.push_back(LLLL4I->pointer());
          } else if (std::shared_ptr<InCoreRITPI<IntsT>> LLLLRI
              = std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(LLLL_)) {
            ERI4s.pointers.push_back(LLLLRI->pointer());
            ERI4s.pointers.push_back(LLLLRI->pointer());
            ERI4s.aux_dims.push_back(LLLLRI->nRIBasis());
          }
          break;
        case DC_COULOMB:
        case DC_EXCHANGE:
          ERI4s = DC_->getPointers(index);
          break;
        case SSSS_COULOMB:
        case SSSS_EXCHANGE:
          ERI4s = SSSS_->getPointers(index);
          break;
        case GAUNT_COULOMB:
        case GAUNT_EXCHANGE:
          ERI4s = gaunt_->getPointers(index);
          break;
        case GAUGE_COULOMB:
        case GAUGE_EXCHANGE:
          ERI4s = gauge_->getPointers(index);
          break;
        default:
          CErr("Unsupported two-body contraction type for InCoreRelERI::getPointers.");
          break;
      }

      return ERI4s;
    }
    
    // Computation interfaces
    virtual void computeAOInts(BasisSet &basisSet, Molecule &mol,
        EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &hamiltonianOptions) override {

      // Use Libcint or Libint to compute nonrelativistic
      LLLL_->computeAOInts(basisSet, mol, emPert, op, hamiltonianOptions);

      HamiltonianOptions tmpOptions(hamiltonianOptions);
      tmpOptions.BareCoulomb = false; // Do bare Coulomb only in Restricted-Kinetic balance (RKB)
      tmpOptions.DiracCoulomb = false; // Dirac-Coulomb without SSSS
      tmpOptions.DiracCoulombSSSS = false; // SSSS to Dirac-Coulomb
      tmpOptions.Gaunt = false; // Gaunt
      tmpOptions.Gauge = false; // Gauge

      if (DC_) {
        tmpOptions.DiracCoulomb = true;
        DC_->computeAOInts(basisSet, mol, emPert, op, tmpOptions);
        tmpOptions.DiracCoulomb = false;
      }

      if (gaunt_) {
        tmpOptions.Gaunt = true;
        gaunt_->computeAOInts(basisSet, mol, emPert, op, tmpOptions);
        tmpOptions.Gaunt = false;
      }

      if (SSSS_) {
        tmpOptions.DiracCoulombSSSS = true;
        SSSS_->computeAOInts(basisSet, mol, emPert, op, tmpOptions);
        tmpOptions.DiracCoulombSSSS = false;
      }

      if (gauge_) {
        tmpOptions.Gauge = true;
        gauge_->computeAOInts(basisSet, mol, emPert, op, tmpOptions);
        tmpOptions.Gauge = false;
      }

    }
    virtual void computeAOInts(BasisSet &basisSet, BasisSet &basisSet2,
                               Molecule &mol, EMPerturbation &emPert, OPERATOR op,
                               const HamiltonianOptions &hamiltonianOptions) {

      CErr("Relativistic TPI with two sets of basis NYI");

    };

    virtual void clear() {
      LLLL_->clear();
      if (DC_) DC_->clear();
      if (gaunt_) gaunt_->clear();
      if (SSSS_) SSSS_->clear();
      if (gauge_) gauge_->clear();
    }

    virtual void output(std::ostream &out, const std::string &s = "",
                        bool printFull = false) const {
      if (printFull) {
        std::string oeiStr;
        if (s == "")
          oeiStr = "RelIncoreERI";
        else
          oeiStr = "RelIncoreERI[" + s + "]";
        prettyPrintSmart(out, oeiStr+".LLLL", this->pointer(),
            this->nBasis(), this->nBasis(), this->nBasis());
        if (DC_) DC_->output(out, oeiStr+".DC", printFull);
        if (gaunt_) gaunt_->output(out, oeiStr+".Gaunt", printFull);
        if (SSSS_) SSSS_->output(out, oeiStr+".SSSS", printFull);
        if (gauge_) gauge_->output(out, oeiStr+".Gauge", printFull);
      } else {
        std::string oeiStr;
        if (s == "")
          oeiStr = "Relativistic incore electron repulsion integral";
        else
          oeiStr = "RelIncoreERI[" + s + "]";
        out << oeiStr;
        out << " with " << nRelComp();
        out << " DCB integrals" << std::endl;
      }
    }

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) override {
      InCoreTPI<IntsT>::broadcast(comm, root);

      LLLL_->broadcast(comm, root);
      if (DC_) DC_->broadcast(comm, root);
      if (gaunt_) gaunt_->broadcast(comm, root);
      if (SSSS_) SSSS_->broadcast(comm, root);
      if (gauge_) gauge_->broadcast(comm, root);
    }

    // Pauil Matrice representation to spinor representation
    template <typename IntsU>
    Incore4indexTPIList<IntsU> spatialToSpinBlock() const;
    
    template <typename TransT>
    InCoreRelERI<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<TransT, dcomplex>::value),
    dcomplex, double>::type> transform(
        char TRANS, const TransT* T, int NT, int LDT) const {
      InCoreRelERI<typename std::conditional<
      (std::is_same<IntsT, dcomplex>::value or
       std::is_same<TransT, dcomplex>::value),
      dcomplex, double>::type> transInts(NT, nRelComp());
      transInts[REL_INTS_COMPS::O] =
          (*this)[REL_INTS_COMPS::O].transform(TRANS, T, NT, LDT);
      for (size_t i = 0; i < nRelComp(); i++)
        transInts[i] = (*this)[i].transform(TRANS, T, NT, LDT);
      return transInts;
    }

    void RI_direct_error(BasisSet &basisSet, Molecule &mol,
                         EMPerturbation &emPert, OPERATOR op,
                         const HamiltonianOptions &hamiltonianOptions,
                         double largeValueThreshold) const;
    
    virtual ~InCoreRelERI() {}

  }; // class InCoreRelERI

  template <typename MatsT, typename IntsT>
  class InCoreRelERIContraction : public InCoreTPIContraction<MatsT,IntsT> { // TODO: rename to InCoreRelERIContraction

    template <typename MatsU, typename IntsU>
    friend class InCoreRelERIContraction;

  public:

    // Constructors

    InCoreRelERIContraction() = delete;
    InCoreRelERIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      InCoreTPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCoreRelERIContraction(
        const InCoreRelERIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCoreRelERIContraction(other.ints_) {}
    template <typename MatsU>
    InCoreRelERIContraction(
        InCoreRelERIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCoreRelERIContraction(other.ints_) {}

    InCoreRelERIContraction( const InCoreRelERIContraction &other ):
      InCoreRelERIContraction(other, 0) {}
    InCoreRelERIContraction( InCoreRelERIContraction &&other ):
      InCoreRelERIContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual ~InCoreRelERIContraction() {}

  }; // class InCoreRelERIContraction

}; // namespace ChronusQ
