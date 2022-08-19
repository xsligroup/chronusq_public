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

  template <typename MatsT, typename IntsT>
  class GTODirectRelERIContraction : public GTODirectTPIContraction<MatsT,IntsT> {
  
  protected:
    
    size_t libcintCacheSize(const TWOBODY_CONTRACTION_TYPE &, int *, const int, 
      int *, const int, double *) const; 
  
  public:

    IntsT *ERI4DCB = nullptr; // Storage for 3-index DCB ERI intermediates

    // Constructors

    GTODirectRelERIContraction() = delete;
    GTODirectRelERIContraction(TwoPInts<IntsT> &tpi):
      GTODirectTPIContraction<MatsT,IntsT>(tpi) {

      if (typeid(tpi) != typeid(DirectTPI<IntsT>))
        CErr("GTODirectRelERIContraction expect a DirectTPI<IntsT> reference.");

    }

    GTODirectRelERIContraction( const GTODirectRelERIContraction &other ):
      GTODirectRelERIContraction(other.ints_) {}
    GTODirectRelERIContraction( GTODirectRelERIContraction &&other ):
      GTODirectRelERIContraction(other.ints_) {}

    // Computation interfaces
    virtual void twoBodyContract(
        MPI_Comm comm,
        const bool screen,
        std::vector<TwoBodyContraction<MatsT>> &list,
        EMPerturbation&,
        const bool computeExchange = true, 
        const TYPE_4C type4C = TYPE_4C::All,
        const APPROXIMATION_TYPE_4C approximate4C = APPROXIMATION_TYPE_4C::None) 
        const {
      
      if (type4C==TYPE_4C::All) 
        directScaffoldLibcint(comm, screen, list, computeExchange, approximate4C); 
      else if (type4C==TYPE_4C::SpinFree) 
        directScaffoldLibcintSpinFree(comm, screen, list, computeExchange, approximate4C); 
      else if (type4C==TYPE_4C::SpinDependent) 
        directScaffoldLibcintSpinDependent(comm, screen, list, computeExchange, approximate4C); 
      //else CErr("twoBodyContract with Coulomb Only is deprecated "); 
     // directScaffold(comm, screen, list);
//      twoBodyContract3Index(comm, list);
    }

    virtual void twoBodyRelContract(
        MPI_Comm comm,
        const bool screen,
        std::vector<TwoBodyRelContraction<MatsT>> &list,
        EMPerturbation&,
        const bool computeExchange = true) const {
      
      if (computeExchange) CErr("Exchange Term NYI in twoBodyRelContract"); 
      else directRelScaffoldLibcintCoulombOnly(comm, screen, list);
    
    }
    
    void directScaffold(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&) const;

    void directScaffoldLibcint(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        const bool,
        const APPROXIMATION_TYPE_4C) const;

    void directScaffoldLibcintSpinFree(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        const bool,
        const APPROXIMATION_TYPE_4C) const;

    void directScaffoldLibcintSpinDependent(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<MatsT>>&,
        const bool,
        const APPROXIMATION_TYPE_4C) const;

    void directRelScaffoldLibcintCoulombOnly(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyRelContraction<MatsT>>&) const;
   
    size_t directRelScaffoldLibcintSCRSize(
      const TWOBODY_CONTRACTION_TYPE &,
      const bool) const;
    
    void twoBodyContract3Index(
        MPI_Comm,
        std::vector<TwoBodyContraction<MatsT>>&) const;

    void JContract3Index(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    void KContract3Index(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    void computeERI3Index(size_t); // Evaluate Spin-Own-Orbit ERIs in the CGTO basis

    virtual ~GTODirectRelERIContraction() {}

  }; // class GTODirectRelERIContraction

}; // namespace ChronusQ
