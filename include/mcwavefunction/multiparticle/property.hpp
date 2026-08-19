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

#include <mcwavefunction.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::computeMultipole()
    {
        for(size_t i = 0; i < this->NStates; i++)
            this->computeMultipole(i);
    }

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::computeMultipole(size_t i)
    {
        ApplyToEach([i](SubMCWfnPtr mcwfn){mcwfn->rdm2pdm(i);});
        EMPerturbation emPert;
        mcSSref_->computeMultipole(emPert);

        std::cout << std::endl;
        std::cout << " *---------------------------------------------*" << std::endl;
        std::cout << " *          Multipole for State " + std::to_string(i+1) + "             *" << std::endl;
        std::cout << " *---------------------------------------------*" << std::endl;
        mcSSref_->printMultipoles(std::cout);

        this->elecDipoles.push_back(mcSSref_->elecDipole);
        this->elecQuadrupoles.push_back(mcSSref_->elecQuadrupole);
        this->elecOctupoles.push_back(mcSSref_->elecOctupole);

        return;
    }

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::formNaturalOrbs(size_t root)
    {
        ApplyToEach([root](SubMCWfnPtr mcwfn){mcwfn->formNaturalOrbs(root);});
    }

}; // namespace ChronusQ
