/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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
#ifdef CQ_HAS_TA
#include <itersolver/solvervectors.hpp>
#include <coupledcluster/TAManager.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  class EOMCCSD;

  template <typename MatsT>
  class EOMCCSDVector {
  protected:
    using TArray = TA::TArray<MatsT>;

    char vLabel_ = 'v';
    char oLabel_ = 'o';
    MatsT V0_ = 0.0;
    TArray V1_;
    TArray V2_;

  public:
    EOMCCSDVector(char vLabel, char oLabel);
    EOMCCSDVector(char vLabel, char oLabel, const TArray &V1, const TArray &V2);
    EOMCCSDVector(char vLabel, char oLabel, MatsT V0, const TArray &V1, const TArray &V2);
    EOMCCSDVector(const EOMCCSDVector<MatsT>&);
    EOMCCSDVector(EOMCCSDVector<MatsT>&&);
    ~EOMCCSDVector();

    size_t length(bool includeZeroBody = false) const {
      TAManager &TAmanager = TAManager::get();
      size_t nO = TAmanager.getRange(oLabel_).extent();
      size_t nV = TAmanager.getRange(vLabel_).extent();
      return nO * nV + nO*(nO - 1)*nV*(nV-1)/4 + (includeZeroBody ? 1 : 0);
    }

    EOMCCSDVector<MatsT>& operator=(const EOMCCSDVector<MatsT>&);

    void swap(EOMCCSDVector<MatsT> &other);

    void assign(MatsT V0, const TArray &V1, const TArray &V2);
    void assign(const TArray &V1, const TArray &V2);

    void initialize();

    MatsT dot(const EOMCCSDVector<MatsT> &other, bool conjA = true) const;

    void axpy(MatsT alpha, const EOMCCSDVector<MatsT> &X);

    double norm() const;

    double absmax() const;

    void scale(MatsT factor);

    void conjugate();

    void normalize();

    void enforceTwoBodySymmetry();

    void projectOut(const EOMCCSDVector<MatsT> &other, bool normalized = false);

    MatsT& zeroBody() { return V0_; }
    const MatsT& zeroBody() const { return V0_; }
    TArray& oneBody() { return V1_; }
    const TArray& oneBody() const { return V1_; }
    TArray& twoBody() { return V2_; }
    const TArray& twoBody() const { return V2_; }

    void setElem(size_t idx, MatsT elem);
    void setZeroBodyElem(MatsT elem);
    void setOneBodyElem(size_t a, size_t i, MatsT elem);
    void setTwoBodyElem(size_t a, size_t b, size_t i, size_t j, MatsT elem);

    void print(std::ostream& out, std::string str) const {
      out << str << "[zero body]: " << V0_ << std::endl;
      out << str << "[one body]:" << std::endl << V1_ << std::endl;
      out << str << "[two body]:" << std::endl << V2_ << std::endl;
    }

    void toRaw(MatsT *raw, bool includeZeroBody = false) const;

    template <typename IntsT>
    void toRaw(MatsT *raw, const EOMCCSD<MatsT,IntsT> &eom, bool includeZeroBody = false) const;

    void fromRaw(const MatsT *raw, bool hasZeroBody = false);

    template <typename IntsT>
    void fromRaw(const MatsT *raw, const EOMCCSD<MatsT,IntsT> &eom, bool hasZeroBody = false);

  };

  template <typename MatsT>
  class EOMCCSDVectorSet : public SolverVectors<MatsT> {
  protected:
    using TArray = TA::TArray<MatsT>;

    char vLabel_ = 'v';
    char oLabel_ = 'o';
    
    std::vector<EOMCCSDVector<MatsT>> vecs_;

    void initialize(size_t nVec);
    
  public:
    EOMCCSDVectorSet(char vLabel, char oLabel, size_t nVec): vLabel_(vLabel), oLabel_(oLabel) {
      initialize(nVec);
    }

    virtual size_t length() const override {
      TAManager &TAmanager = TAManager::get();
      size_t nV = TAmanager.getRange(vLabel_).extent();
      size_t nO = TAmanager.getRange(oLabel_).extent();
      return nO * nV + nO*(nO - 1)*nV*(nV-1)/4;
    }
    size_t length(bool includeZeroBody) const {
      return length() + (includeZeroBody ? 1 : 0);
    }
    template <typename IntsT>
    size_t length(const EOMCCSD<MatsT,IntsT> &eom, bool includeZeroBody) const {
      return eom.getHbarDim() + (includeZeroBody ? 1 : 0);
    }

    virtual size_t size() const override { return vecs_.size(); }

    // Get element
    virtual MatsT get(size_t i, size_t j) const override {
      CErr("Get element in EOMCCSDVectorSet object is invalid.");
      return 0.0;
    }
    // Set element
    virtual void set(size_t i, size_t j, MatsT value) override {
      vecs_[j].setElem(i, value);
    }

    // Get EOMCCSDVector
    EOMCCSDVector<MatsT>& get(size_t i) {
      this->sizeCheck(i, "in EOMCCSDVectorSet<MatsT>::get");
      return vecs_[i];
    }
    const EOMCCSDVector<MatsT>& get(size_t i) const {
      this->sizeCheck(i, "in EOMCCSDVectorSet<MatsT>::get");
      return vecs_[i];
    }
    // Set EOMCCSDVector
    virtual void set(size_t i, const EOMCCSDVector<MatsT>& other) {
      vecs_[i] = other;
    }

    using SolverVectors<MatsT>::clear;
    void clear(size_t shift, size_t nVec) override {
      if (nVec == 0) return;
      this->sizeCheck(shift + nVec, "in EOMCCSDVectorSet<MatsT>::clear");
      for (size_t i = shift; i < vecs_.size(); i++)
        vecs_[i].scale(0.0);
    }

    using SolverVectors<MatsT>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      if (nVec == 0) {
        out << std::endl << str + ": " << std::endl;
        return;
      }
      this->sizeCheck(shift + nVec, "in EOMCCSDVectorSet<MatsT>::print");
      for (size_t i = 0; i < nVec; i++) {
        get(i + shift).print(out, str + "(" + std::to_string(i) + ")");
      }
    }

    RawVectors<MatsT> toRaw(MPI_Comm c,
                            bool includeZeroBody = false, size_t shift = 0,
                            size_t nVec = std::numeric_limits<size_t>::max()) const;

    template <typename IntsT>
    RawVectors<MatsT> toRaw(MPI_Comm c, const EOMCCSD<MatsT,IntsT> &eom,
                            bool includeZeroBody = false, size_t shift = 0,
                            size_t nVec = std::numeric_limits<size_t>::max()) const;

    void fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw, bool hasZeroBody = false,
                 size_t shiftThis = 0, size_t shiftRaw = 0, size_t nVec = std::numeric_limits<size_t>::max());

    template <typename IntsT>
    void fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw,
                 const EOMCCSD<MatsT,IntsT> &eom, bool hasZeroBody = false,
                 size_t shiftThis = 0, size_t shiftRaw = 0, size_t nVec = std::numeric_limits<size_t>::max());

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 MatsT alpha, MatsT const *B, int64_t ldb,
                                 MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                             int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) override;

    using SolverVectors<MatsT>::scale;
    virtual void scale(MatsT scalar, size_t shift, size_t nVec) override;

    using SolverVectors<MatsT>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) override;

//    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
//                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, MatsT *R = nullptr, int LDR = 0) override;

    using SolverVectors<MatsT>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<MatsT>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;
    
    virtual ~EOMCCSDVectorSet() {}
    
  };



  template <typename MatsT>
  class EOMCCSDVectorSetDebug : public SolverVectors<MatsT> {

  protected:
    using TArray = TA::TArray<MatsT>;
    EOMCCSDVectorSet<MatsT> eomccSet_;
    RawVectors<MatsT> rawSet_;

  public:
    EOMCCSDVectorSetDebug(char vLabel, char oLabel, size_t nVec,
                          MPI_Comm c):
                          eomccSet_(vLabel, oLabel, nVec),
                          rawSet_(c, eomccSet_.length(), nVec) {}

    double compareDebug(size_t shift = 0, size_t nVec = std::numeric_limits<size_t>::max());

    EOMCCSDVectorSet<MatsT>& getEOMCCSet() {
      return eomccSet_;
    }

    RawVectors<MatsT>& getRawSet() {
      return rawSet_;
    }

    virtual size_t length() const override { return eomccSet_.length(); }
    virtual size_t size() const override { return eomccSet_.size(); }

    // Get element
    virtual MatsT get(size_t i, size_t j) const override {
      return eomccSet_.get(i, j);
    }
    // Set element
    virtual void set(size_t i, size_t j, MatsT value) override {
      eomccSet_.set(i, j, value);
      rawSet_.set(i, j, value);
    }

    using SolverVectors<MatsT>::clear;
    void clear(size_t shift, size_t nVec) override {
      eomccSet_.clear(shift, nVec);
      rawSet_.clear(shift, nVec);
    }

    using SolverVectors<MatsT>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      eomccSet_.print(out, str, shift, nVec);
      rawSet_.print(out, str, shift, nVec);
    }

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 MatsT alpha, MatsT const *B, int64_t ldb,
                                 MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                             int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) override;

    using SolverVectors<MatsT>::scale;
    virtual void scale(MatsT scalar, size_t shift, size_t nVec) override;

    using SolverVectors<MatsT>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) override;

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, MatsT *R = nullptr, int LDR = 0) override;

    using SolverVectors<MatsT>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<MatsT>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual ~EOMCCSDVectorSetDebug() {}
  };

}; // namespace ChronusQ
#endif
