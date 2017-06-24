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

#include <integrals.hpp>
#include <particleintegrals/contract/incore.hpp>
#include <particleintegrals/contract/direct.hpp>
#include <particleintegrals/contract/direct4C.hpp>
#include <particleintegrals/onepints/relativisticints.hpp>


namespace ChronusQ {

  template <>
  void InCore4indexTPIContraction<double,dcomplex>::JContract(
      MPI_Comm, TwoBodyContraction<double>&) const { CErr("NYI"); }
  template <>
  void InCore4indexTPIContraction<double,dcomplex>::KContract(
      MPI_Comm, TwoBodyContraction<double>&) const { CErr("NYI"); }

  template <>
  void InCoreRITPIContraction<double,dcomplex>::JContract(
      MPI_Comm, TwoBodyContraction<double>&) const { CErr("NYI"); }
  template <>
  void InCoreRITPIContraction<double,dcomplex>::KContract(
      MPI_Comm, TwoBodyContraction<double>&) const { CErr("NYI"); }

  template class OnePInts<double>;
  template class OnePInts<dcomplex>;

  template class OnePRelInts<double>;
  template class OnePRelInts<dcomplex>;

  template OnePInts<dcomplex>::OnePInts(const OnePInts<double>&, int);
  template OnePRelInts<dcomplex>::OnePRelInts(const OnePRelInts<double>&, int);
  template MultipoleInts<dcomplex>::MultipoleInts(const MultipoleInts<double>&, int);
  template DirectTPI<dcomplex>::DirectTPI(const DirectTPI<double>&, int);
  template InCore4indexTPI<dcomplex>::InCore4indexTPI(const InCore4indexTPI<double>&, int);
  template InCoreRITPI<dcomplex>::InCoreRITPI(const InCoreRITPI<double>&, int);
  template InCoreAsymmRITPI<dcomplex>::InCoreAsymmRITPI(const InCoreAsymmRITPI<double>&, int);
  template InCoreAuxBasisRIERI<dcomplex>::InCoreAuxBasisRIERI(const InCoreAuxBasisRIERI<double>&, int);

  template class InCore4indexTPIContraction<double, double>;
  template class InCore4indexTPIContraction<dcomplex, double>;
  template class InCore4indexTPIContraction<dcomplex, dcomplex>;

  template class InCoreRITPIContraction<double, double>;
  template class InCoreRITPIContraction<dcomplex, double>;
  template class InCoreRITPIContraction<dcomplex, dcomplex>;

  template class InCoreAsymmRITPIContraction<double, double>;
  template class InCoreAsymmRITPIContraction<dcomplex, double>;
  template class InCoreAsymmRITPIContraction<dcomplex, dcomplex>;

  template class GTODirectTPIContraction<double, double>;
  template class GTODirectTPIContraction<dcomplex, double>;
  template class GTODirectTPIContraction<dcomplex, dcomplex>;

  template class InCore4indexRelERIContraction<double, double>;
  template class InCore4indexRelERIContraction<dcomplex, double>;
  template class InCore4indexRelERIContraction<dcomplex, dcomplex>;

  template class GTODirectRelERIContraction<double, double>;
  template class GTODirectRelERIContraction<dcomplex, double>;
  template class GTODirectRelERIContraction<dcomplex, dcomplex>;

  template class Integrals<double>;
  template class Integrals<dcomplex>;

  template class GradContractions<double,double>;
  template class GradContractions<dcomplex,double>;
  template class GradContractions<dcomplex,dcomplex>;

  template class InCore4indexGradContraction<double, double>;
  template class InCore4indexGradContraction<dcomplex, double>;
  template class InCore4indexGradContraction<dcomplex, dcomplex>;

  template class DirectGradContraction<double, double>;
  template class DirectGradContraction<dcomplex, double>;
  template class DirectGradContraction<dcomplex, dcomplex>;


}; // namespace ChronusQ
