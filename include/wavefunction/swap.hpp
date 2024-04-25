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

#include <wavefunction/base.hpp>
#include <matrix.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::swapMOs(
      std::vector<std::vector<std::pair<size_t, size_t>>> & moP, SpinType spin){
 
    if (moP[spin].empty()) return;
    
    auto MO = mo[spin].pointer();
    size_t LDMO = mo[spin].dimension();
    MatsT * SCR = CQMemManager::get().malloc<MatsT>(LDMO);

    if( spin==0 ) std::cout << "  * the following MOs are swapped" << std::endl;
    
    else if( spin==1 ) std::cout << "  * the following beta MOs are swapped" << std::endl;

    for( auto & pair: moP[spin] ){

      std::cout << "    " << std::setw(5) << pair.first << " <--> " 
                << std::setw(5) << pair.second << std::endl;

      MatsT * first_p = MO + (pair.first-1) * LDMO;
      MatsT * second_p = MO + (pair.second-1) * LDMO;

      std::copy_n(first_p,  LDMO, SCR);
      std::copy_n(second_p, LDMO, first_p);
      std::copy_n(SCR,      LDMO, second_p);

    }
    
    CQMemManager::get().free(SCR);

  }; // MCWaveFunction:::swapMOs

}; // namespace ChronusQ
