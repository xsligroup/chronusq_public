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

#include <mcscf.hpp>
#include <mcwavefunction.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MultiParticleMCWaveFunction<MatsT,IntsT>::computeOneRDM(size_t i)
  {
    // Build the vector of reference wrappers for the RDMs
    std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> rdmRefs;
    ApplyToEach([&rdmRefs,i](SubMCWfnPtr & mcwfn){rdmRefs.emplace_back(mcwfn->oneRDM[i]);});

    multiparticleciBuilder->computeOneRDM(*this,this->CIVecs[i],rdmRefs);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleMCWaveFunction<MatsT,IntsT>::computeOneRDM()
  {
    for(size_t i = 0; i < this->NStates; i++)
    {
        computeOneRDM(i);
    }
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleMCWaveFunction<MatsT,IntsT>::computeTDMs()
  {
    CErr("TDMs for MultiParticleMCWaveFunction NYI!");    
  }

}; // namespace ChronusQ
