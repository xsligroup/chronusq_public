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

#include <fockbuilder.hpp>


namespace ChronusQ {

  /**
   * \brief The FourCompFock class
   */
  template <typename MatsT, typename IntsT>
  class FourCompFock : public FockBuilder<MatsT,IntsT> {
  public:

    // Constructors
    FourCompFock() = delete;
    FourCompFock(HamiltonianOptions hamiltonianOptions):
        FockBuilder<MatsT,IntsT>(hamiltonianOptions) {}

    // Different type
    template <typename MatsU>
    FourCompFock(const FourCompFock<MatsU,IntsT> &other):
        FockBuilder<MatsT,IntsT>(other){}
    template <typename MatsU>
    FourCompFock(FourCompFock<MatsU,IntsT> &&other):
        FockBuilder<MatsT,IntsT>(other){}

    // Virtual destructor
    virtual ~FourCompFock() {}


    // Public member functions

    // Specialized formGD function for the 4C Hamiltonian 
    // (see include/fockbuilder/fourcompfock/impl.hpp)
    void formGD(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1., bool HerDen = true);
    void formGDInCore(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1., bool HerDen = true);
    void formGDDirect(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1., bool HerDen = true);
    void formGD3Index(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1., bool HerDen = true);
    
    void formRawGDInBatches(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool, double, bool, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &,
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &);

    void formRawGDInBatchesDirect(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool, double, bool, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &,
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> &);
    
    size_t formRawGDSCRSizePerBatch(SingleSlater<MatsT,IntsT> &, bool, bool) const;

    // Form a fock matrix (see include/fockbuilder/impl.hpp for docs)
    virtual void formFock(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1.);

    // Compute the gradient
    virtual void getGrad() {
      CErr("4CHF Fock gradient NYI",std::cout);
    }


  };  
  
  
}

