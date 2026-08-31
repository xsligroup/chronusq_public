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

#include <realtime/realtimemultislater/impl.hpp>
#include <realtime/realtimeci/impl.hpp>

namespace ChronusQ {

  template <>
  void RTMS::scal(std::shared_ptr<SolverVectors<double>> source, dcomplex factor) {
    CErr("Invalid type combination for scal");
  };

  template <>
  void RTMS::normalize(std::shared_ptr<SolverVectors<double>> source, dcomplex &result) {
    CErr("Invalid type combination for norm");
  };
  /*
  template <>
  void RTMS::dot(std::shared_ptr<SolverVectors<dcomplex>> source_1, std::shared_ptr<SolverVectors<dcomplex>> source_2, double &result) {
    CErr("Invalid type combination for dot");
  };
  */
  
  template <>
  void RealTimeMultiSlaterBase<dcomplex, dcomplex>::propagateWFN_SSO(bool, bool) {
    CErr("Invalid");
  };

  template <>
  void RealTimeMultiSlaterBase<dcomplex, double>::propagateWFN_SSO(bool, bool) {
    CErr("Invalid");
  };
  template <>
  void RealTimeMultiSlaterBase<double, double>::propagateWFN_RK4(bool, bool) {
    CErr("Invalid");
  };

  template class RealTimeMultiSlaterBase<double, double >; 
  template class RealTimeMultiSlaterBase<dcomplex, double>; 
  template class RealTimeCI<double, double >;
  template class RealTimeCI<dcomplex, double>;
}; // namespace ChronusQ
