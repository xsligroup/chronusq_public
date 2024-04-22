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
#include <util/matout.hpp>

namespace ChronusQ {

  /**
   *  \brief Templated class to handle the storage of
   *  one electron integral matrix within arbitrary active spaces
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets
   *  (e.g., GTO and GIAO)
   */
  template <typename IntsT>
  class DASOnePInts : public ParticleIntegrals {

  protected:

    cqmatrix::Matrix<IntsT> mat_; ///< DAS One Particle Ints (2 index)
    size_t nb1_  = 0ul;
    size_t nb2_  = 0ul;

  public:

    // Constructor
    DASOnePInts() = delete;
    DASOnePInts(size_t nb1, size_t nb2):
        ParticleIntegrals(0ul), mat_(nb1, nb2), nb1_(nb1), nb2_(nb2){}
    DASOnePInts(const DASOnePInts &other ) = default;
    DASOnePInts(DASOnePInts &&other ) = default;
    
    size_t nBasis1() const { return nb1_; } 
    size_t nBasis2() const { return nb2_; } 

    DASOnePInts& operator=(const DASOnePInts &other ) = default;

    DASOnePInts& operator=(DASOnePInts &&other ) = default;

    IntsT& operator()(size_t p, size_t q) {
      return mat_(p, q);
    }
    IntsT operator()(size_t p, size_t q) const {
      return mat_(p, q);
    }

    // Matrix direct access
    IntsT* pointer() { return mat_.pointer(); }
    const IntsT* pointer() const { return mat_.pointer(); }

    // disable Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) { 
      CErr("computeAOInts is not supported in GASOnePInt"); 
    };
    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) { 
      CErr("computeAOInts is not supported in GASOnePInt"); 
    };

    void clear() {
      mat_.clear();
    }

    void output(std::ostream &out, const std::string &s = "",
                   bool printFull = false) const {
      std::string opiStr;
      if (s == "")
        opiStr = "One-Particle integral";
      else
        opiStr = "OPI[" + s + "]";
      
      opiStr += ", Dimension: (" + std::to_string(nb1_) + "," 
               + std::to_string(nb2_) + ")"; 
      
      if (printFull)
        prettyPrintSmart(out, opiStr, pointer(), this->nBasis1(),
                         this->nBasis2(), this->nBasis1());
      else {
        out << opiStr << std::endl;
      }
    }

    ~DASOnePInts() {}

  }; // class DASOnePInts

  template <typename IntsT>
  class GASTwoPInts : public ParticleIntegrals {

  protected:

    IntsT * ptr_; ///< GAS Two Particle Ints (4 index)
    size_t nb1_;
    size_t nb2_;
    size_t nb3_;
    size_t nb4_;
    size_t nb12_;
    size_t nb123_;
    size_t nb1234_;

  public:

    // Constructor
    GASTwoPInts() = delete;
    GASTwoPInts(size_t nb1, size_t nb2,
        size_t nb3, size_t nb4):
        ParticleIntegrals(0ul), nb1_(nb1), nb2_(nb2),
        nb3_(nb3), nb4_(nb4) { 
      nb12_ = nb1_ * nb2_;
      nb123_ = nb12_ * nb3_;
      nb1234_ = nb123_ * nb4_;
      ptr_ = CQMemManager::get().malloc<IntsT>(nb1234_);
    }
    GASTwoPInts( const GASTwoPInts &other ): 
        GASTwoPInts(other.nb1_, other.nb2_,
        other.nb3_, other.nb4_) {
      std::copy_n(other.ptr_, nb1234_, ptr_);
    };
    GASTwoPInts( GASTwoPInts &&other ): 
        ParticleIntegrals(0ul), 
        nb1_(other.nb1_), nb2_(other.nb2_),
        nb3_(other.nb3_), nb4_(other.nb4_),
        ptr_(other.ptr_) {
      nb12_ = nb1_ * nb2_;
      nb123_ = nb12_ * nb3_;
      nb1234_ = nb123_ * nb4_;
      other.ptr_ = nullptr;
    };
    
    size_t nBasis1() const { return nb1_; } 
    size_t nBasis2() const { return nb2_; } 
    size_t nBasis3() const { return nb3_; } 
    size_t nBasis4() const { return nb4_; } 

    GASTwoPInts& operator=( const GASTwoPInts &other ) {
      if (nb1234_ != other.nb1234_) {
        CQMemManager::get().free(ptr_);
        ptr_ = CQMemManager::get().malloc<IntsT>(nb1234_);
        nb1_  = other.nb1_;
        nb2_  = other.nb2_;
        nb3_  = other.nb3_;
        nb4_  = other.nb4_;
      }
      std::copy_n(other.ptr_, nb1234_, ptr_);
      return *this; 
    }
    GASTwoPInts& operator=( GASTwoPInts &&other ) {
      return GASTwoPInts(other);
    }

    IntsT& operator()(size_t p, size_t q, size_t r, size_t s) {
      return ptr_[p + q * nb1_ + r * nb12_ + s * nb123_];
    }
    IntsT operator()(size_t p, size_t q, size_t r, size_t s) const {
      return ptr_[p + q * nb1_ + r * nb12_ + s * nb123_];
    }

    // Matrix direct access
    IntsT* pointer() { return ptr_; }
    const IntsT* pointer() const { return ptr_; }
    std::vector<size_t> offsets() const { return {1, nb1_, nb12_, nb123_}; }   

    // disable Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) { 
      CErr("computeAOInts is not supported in GASTwoPInt"); 
    };
    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) { 
      CErr("computeAOInts is not supported in GASTwoPInt"); 
    };

    void clear() {
      std::fill_n(ptr_, nb1234_, IntsT(0.));
    }

    void output(std::ostream &out, const std::string &s = "",
                   bool printFull = false) const {
      std::string opiStr;
      if (s == "")
        opiStr = "Two-Particle integral";
      else
        opiStr = "TPI[" + s + "]";
      
      opiStr += ", Dimension: (" + std::to_string(nb1_) + "," 
               + std::to_string(nb2_) + "," 
               + std::to_string(nb3_) + "," 
               + std::to_string(nb4_) + ")"; 
      
      if (printFull)
        prettyPrintSmart(out, opiStr, pointer(), nb12_,
                         nb3_*nb4_, nb12_);
      else {
        out << opiStr << std::endl;
      }
    }

    ~GASTwoPInts() { CQMemManager::get().free(ptr_); }

  }; // class GASTwoPInts

}; // namespace ChronusQ
