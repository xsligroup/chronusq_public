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
#include <particleintegrals/twopints/incore4indextpi.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  class InCore4indexGradContraction : public GradContractions<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCore4indexGradContraction;

  public:

    InCore4indexGradContraction() = delete;
    InCore4indexGradContraction(GradInts<TwoPInts,IntsT> &grad):
      GradContractions<MatsT,IntsT>(grad) {}

    template <typename MatsU>
    InCore4indexGradContraction(
        const InCore4indexGradContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCore4indexGradContraction(other.grad_) {}
    template <typename MatsU>
    InCore4indexGradContraction(
        InCore4indexGradContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCore4indexGradContraction(other.grad_) {}

    InCore4indexGradContraction( const InCore4indexGradContraction &other ):
      InCore4indexGradContraction(other, 0) {}
    InCore4indexGradContraction( InCore4indexGradContraction &&other ):
      InCore4indexGradContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void gradTwoBodyContract(
        MPI_Comm comm,
        const bool,
        std::vector<std::vector<TwoBodyContraction<MatsT>>> &list,
        EMPerturbation&) const;

    virtual ~InCore4indexGradContraction() {}

  }; // Class InCore4indexGradContraction

}; // namespace ChronusQ
