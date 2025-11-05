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

#include <realtime.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void RealTimeNEOCI<MatsT, IntsT>::genCubes()
    {
      // Current iteration number
      std::string curiter = std::to_string(this->curState.iStep);

      // Because this is called AFTER the dipole is calculated
      // the onePDM is already populated in the singleslater reference
      // and we can simply evaluate the cube normally
      auto ecube = this->intScheme.rtcubes[PAR_TYPE::ELECTRONIC];
      auto pcube = this->intScheme.rtcubes[PAR_TYPE::PROTONIC];
      std::string ecube_name;
      std::string pcube_name;
      if(this->intScheme.cubeOptsRTMS.cubeFileName.empty())
      {
        ecube_name = "RTCI_ELEC";
        pcube_name = "RTCI_PROT";
      }
      else
      {
        ecube_name = this->intScheme.cubeOptsRTMS.cubeFileName + "_RTCI_ELEC";
        pcube_name = this->intScheme.cubeOptsRTMS.cubeFileName + "_RTCI_PROT";
      }

      SingleSlater<MatsT,IntsT> * ess_ptr = &neomcref_->ewfn_->reference();
      SingleSlater<MatsT,IntsT> * pss_ptr = &neomcref_->pwfn_->reference();
        
      ecube_name += "_"+std::to_string(int(this->curState.iStep)); 
      pcube_name += "_"+std::to_string(int(this->curState.iStep)); 

      ecube->evalDenCube(ecube_name,ess_ptr->onePDM);
      pcube->evalDenCube(pcube_name,pss_ptr->onePDM);

    }

}; // namespace ChronusQ
