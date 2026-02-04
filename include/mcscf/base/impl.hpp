/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you ca redistribute it and/or modify
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

#include <mcwavefunction/base.hpp>
#include <detstringmanager.hpp>

namespace ChronusQ {

  void MCSCFBase::turnOnStateAverage(const std::vector<double> & weight) {
    
    size_t NS = this->NStates;
    
    if( weight.size() != NS) 
      CErr("MCSCF needs "+std::to_string(NS)+" weights for state average" );   

    this->StateAverage = true;
    this->SAWeight = std::vector<double>(NS);
    std::copy_n(weight.begin(), NS, this->SAWeight.begin());
  };
  
}; // namespace ChronusQ
