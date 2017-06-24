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
   * \brief The ROFock class
   */
  template <typename MatsT, typename IntsT>
  class ROFock : public FockBuilder<MatsT,IntsT> {
  public:

    // Constructors
    ROFock() = delete;
    ROFock(HamiltonianOptions hamiltonianOptions):
        FockBuilder<MatsT,IntsT>(hamiltonianOptions) {}

    // Different type
    template <typename MatsU>
    ROFock(const ROFock<MatsU,IntsT> &other):
        FockBuilder<MatsT,IntsT>(other){}
    template <typename MatsU>
    ROFock(ROFock<MatsU,IntsT> &&other):
        FockBuilder<MatsT,IntsT>(other){}

    // Virtual destructor
    virtual ~ROFock() {}


    // Public member functions

    // Form an Roothaan's effective fock for ROHF (see include/fockbuilder/ROFock/impl.hpp for docs)
    void rohfFock(SingleSlater<MatsT,IntsT> &);

    // Form a fock matrix (see include/fockbuilder/impl.hpp for docs)
    virtual void formFock(SingleSlater<MatsT,IntsT> &, EMPerturbation &, bool increment = false, double xHFX = 1.);

    // Compute the gradient
    virtual void getGrad() {
      CErr("ROHF Fock gradient NYI",std::cout);
    }


  };

}
