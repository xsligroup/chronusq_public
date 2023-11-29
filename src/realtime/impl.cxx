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

#include <realtime/impl.hpp>
#include <orbitalmodifiernew/realtimeSCF/impl.hpp>

namespace ChronusQ {

  template class RealTime<HartreeFock,double>;
  template class RealTime<HartreeFock,dcomplex>;
  template class RealTime<KohnSham,double>;
  template class RealTime<KohnSham,dcomplex>;
  template class RealTime<NEOSS,double>;
  template class RealTime<NEOSS,dcomplex>;

  template void RealTime<HartreeFock,double>::alloc<double>();
  template void RealTime<HartreeFock,double>::alloc<dcomplex>();
  template void RealTime<HartreeFock,dcomplex>::alloc<dcomplex>();
  template void RealTime<KohnSham,double>::alloc<double>();
  template void RealTime<KohnSham,double>::alloc<dcomplex>();
  template void RealTime<KohnSham,dcomplex>::alloc<dcomplex>();

  template void RealTime<NEOSS,double>::alloc<double>();
  template void RealTime<NEOSS,double>::alloc<dcomplex>();
  template void RealTime<NEOSS,dcomplex>::alloc<dcomplex>();

  //template class RealTimeSCF<HartreeFock,double,double>;
  template class RealTimeSCF<HartreeFock,dcomplex,double>;
  template class RealTimeSCF<HartreeFock,dcomplex,dcomplex>;
  //template class RealTimeSCF<KohnSham,double,double>;
  template class RealTimeSCF<KohnSham,dcomplex,double>;
  template class RealTimeSCF<KohnSham,dcomplex,dcomplex>;
  //template class RealTimeSCF<NEOSS,double,double>;
  template class RealTimeSCF<NEOSS,dcomplex,double>;
  template class RealTimeSCF<NEOSS,dcomplex,dcomplex>;

}; // namespace ChronusQ
