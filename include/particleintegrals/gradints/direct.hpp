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

#include <particleintegrals/twopints.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  class DirectGradContraction : public GradContractions<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class DirectGradContraction;

  public:

    DirectGradContraction() = delete;
    DirectGradContraction(GradInts<TwoPInts,IntsT> &grad):
      GradContractions<MatsT,IntsT>(grad) {}

    template <typename MatsU>
    DirectGradContraction(
        const DirectGradContraction<MatsU,IntsT> &other, int dummy = 0 ):
      DirectGradContraction(other.grad_) {}
    template <typename MatsU>
    DirectGradContraction(
        DirectGradContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      DirectGradContraction(other.grad_) {}

    DirectGradContraction( const DirectGradContraction &other ):
      DirectGradContraction(other, 0) {}
    DirectGradContraction( DirectGradContraction &&other ):
      DirectGradContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void gradTwoBodyContract(
        MPI_Comm comm, 
        const bool screen,
        std::vector<std::vector<TwoBodyContraction<MatsT>>>& cList,
        EMPerturbation&) const {

      if (std::is_same<IntsT,dcomplex>::value)
        CErr("GIAO gradients NYI!");

      // Hackily turns off screeing if threshSchwarz_ is 0.0
      DirectTPI<IntsT> &tpi = dynamic_cast<DirectTPI<IntsT>&>(*this->grad_[0]);
      if (tpi.threshSchwarz() == 0.0) {
        directScaffoldGrad(comm, false, cList);
      } else {
        directScaffoldGrad(comm, screen, cList);
      }
    };

    virtual void gradTwoBodyTraceContract(
        MPI_Comm comm, bool screen,
        std::vector<TwoBodyContraction<MatsT>>& twoBodyContraction,
        const std::vector<const MatsT*>& traceDensities,
        const std::vector<double>& traceCoeffs,
        std::vector<double>& gradientOut,
        EMPerturbation&) const{

      if (std::is_same<IntsT,dcomplex>::value)
        CErr("trace-mode gradients NYI for complex integrals");

      // Hackily turns off screening if threshSchwarz_ is 0.0
      DirectTPI<IntsT> &tpi = dynamic_cast<DirectTPI<IntsT>&>(*this->grad_[0]);
      const bool doScreen = (tpi.threshSchwarz() == 0.0) ? false : screen;

      std::vector<std::vector<TwoBodyContraction<MatsT>>> cList(this->grad_.size(), twoBodyContraction);

      directScaffoldGradImpl(comm, doScreen, cList, traceDensities, traceCoeffs, &gradientOut);
    };

    void directScaffoldGrad(
        MPI_Comm,
        const bool,
        std::vector<std::vector<TwoBodyContraction<MatsT>>>&) const;

    // Shared scaffold for both output modes. A null gradientOut writes
    // derivative Fock matrices into cList[..].AX; a non-null gradientOut
    // accumulates scalar gradient contributions directly, without ever
    // forming those matrices. The mode is resolved once, outside the
    // quartet loops.
    void directScaffoldGradImpl(
        MPI_Comm,
        const bool,
        std::vector<std::vector<TwoBodyContraction<MatsT>>>&,
        const std::vector<const MatsT*>&,
        const std::vector<double>&,
        std::vector<double>*) const;

    virtual ~DirectGradContraction() {}

  }; // Class DirectGradContraction

}; // namespace ChronusQ
