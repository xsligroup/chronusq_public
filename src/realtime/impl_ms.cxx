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

#include <realtime/realtimemultislater/impl.hpp>

namespace ChronusQ {


  template <>
  void RTMS::scal(double* source, size_t vecSize_, dcomplex actor) {
    CErr("Invalid type combination for scal");
  };

  template <>
  void RTMS::dot(dcomplex* source_1, dcomplex* source_2, size_t vecSize_, double& result) {
    CErr("Invalid type combination for dot");
  };
  
  template <>
  void RealTimeMultiSlater<dcomplex, dcomplex>::propagateWFN_SSO(bool, bool) {
    CErr("Invalid");
  };

  template <>
  void RealTimeMultiSlater<dcomplex, double>::propagateWFN_SSO(bool, bool) {
    CErr("Invalid");
  };


  template class RealTimeMultiSlater<double, double >; 
  template class RealTimeMultiSlater<dcomplex, double>; 
}; // namespace ChronusQ
