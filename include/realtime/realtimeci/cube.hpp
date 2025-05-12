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
    void RealTimeCI<MatsT, IntsT>::genCubes()
    {
        // Current iteration number
        std::string curiter = std::to_string(this->curState.iStep);

        // Because this is called AFTER the dipole is calculated
        // the onePDM is already populated in the singleslater reference
        // and we can simply evaluate the cube normally
        auto cube = this->intScheme.rtcubes[PAR_TYPE::ELECTRONIC];

        std::string cube_name;
        if(this->intScheme.cubeOptsRTMS.cubeFileName.empty())
        {
            cube_name = "RTCI";
        }
        else
        {
            cube_name = this->intScheme.cubeOptsRTMS.cubeFileName + "_RTCI";
        }

        SingleSlater<MatsT,IntsT> * ss_ptr = &reference_->reference();
        cube_name += "_"+std::to_string(int(this->curState.iStep)); 
        
        cube->evalDenCube(cube_name,ss_ptr->onePDM);

    }

}; // namespace ChronusQ
