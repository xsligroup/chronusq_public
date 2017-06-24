/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of WashinGIAOn)
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

#include <particleintegrals/twopints/gtodirecttpi.hpp>

namespace ChronusQ {

  class GIAODirectERIContraction : public GTODirectTPIContraction<dcomplex,dcomplex> {

    // typeid for object pointer
    template <typename T> const std::type_info& pointer_to_typeid(const T& v) { return typeid(v); }

  public:
    typedef dcomplex ResultsT;

    // Constructors

    GIAODirectERIContraction() = delete;
    GIAODirectERIContraction(std::shared_ptr<TwoPInts<dcomplex>> eri):
      GTODirectTPIContraction<dcomplex,dcomplex>(eri) {

      const std::type_info& ti_eri = pointer_to_typeid(*eri);
      const std::type_info& ti_direct_tpi = typeid(DirectTPI<dcomplex>);
      if (ti_eri != ti_direct_tpi)
        CErr("GIAODirectERIContraction expect a DirectERI<dcomplex> reference.");

    }

    GIAODirectERIContraction( const GIAODirectERIContraction &other ):
      GIAODirectERIContraction(other.ints_) {}
    GIAODirectERIContraction( GIAODirectERIContraction &&other ):
      GIAODirectERIContraction(other.ints_) {}

    virtual void twoBodyContract(
        MPI_Comm c,
        const bool screen,
        std::vector<TwoBodyContraction<dcomplex>> &list,
        EMPerturbation &pert) const;

    virtual ~GIAODirectERIContraction() {}

  }; // class GIAODirectERIContraction

}; // namespace ChronusQ
