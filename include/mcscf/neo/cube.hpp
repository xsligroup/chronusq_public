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

#include <mcscf.hpp>
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>
#include <cubegen.hpp>


namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::runCube(std::vector<std::shared_ptr<CubeGen>> cubes)
    {
        SingleSlater<MatsT,IntsT> * ess = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->ewfn_->reference());
        SingleSlater<MatsT,IntsT> * pss = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());

        // cubes[0] is the electronic cube
        // cubes[1] is the protonic
        std::shared_ptr<CubeGen> ecube = cubes[PAR_TYPE::ELECTRONIC];
        std::shared_ptr<CubeGen> pcube = cubes[PAR_TYPE::PROTONIC];

        // Generate density cubes
        auto ePDMs = this->getOnePDM();
        auto pPDMs = this->getPOnePDM();

        std::string cube_name;
        if(this->cubeOptsMC.cubeFileName.empty()) {
            cube_name = "NEOMCSCF";
        } else {
            cube_name = this->cubeOptsMC.cubeFileName;
            cube_name = cube_name + "_MCWFN";
        }

        std::function<double(MatsT)> ReOrMag = [](MatsT x){return std::real(x);};
        // density cube
        if(this->cubeOptsMC.denCube)
        {

            // Handle which roots to generate
            size_t NRoots = this->NStates;
            std::vector<size_t> RootsToCube;
            if(this->cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::GS)
            {
                RootsToCube.push_back(0);
            }
            else if(this->cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::ALL)
            {
                for(size_t i = 0; i < NRoots; i++)
                    RootsToCube.push_back(i);
            }
            else if(this->cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::AVERAGE)
            {
                CErr("State Average Density cube NYI for NEOMCSCF");
            }
            else
            {
                RootsToCube=this->cubeOptsMC.custom_root_request;
            }

            for(const auto Root : RootsToCube)
            {
                std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> etemp_pdm = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(ePDMs[Root]->nRows(),false,false);
                etemp_pdm->S() = *ePDMs[Root];
                ecube->evalDenCube(cube_name+"_ELEC_"+std::to_string(Root+1),etemp_pdm);
                std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> ptemp_pdm = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(pPDMs[Root]->nRows(),false,false);
                ptemp_pdm->S() = *pPDMs[Root];
                pcube->evalDenCube(cube_name+"_PROT_"+std::to_string(Root+1),ptemp_pdm,1.0);

            }
        }
        // Orbital cubes
        if (this->cubeOptsMC.orbCube) 
        {
            // Electronic orbitals 
            {
                size_t NB = ess->nAlphaOrbital();
                size_t NOrb = NB * ess->nC;

                // Handle which orbitals to generate
                std::vector<size_t> OrbsToCube;
                if(this->cubeOptsMC.whichMO == MO_CLASSES::ALL)
                {
                    for(size_t i = 0; i < NOrb; i++)
                        OrbsToCube.push_back(i);
                }
                // Base case is custom vector
                else
                {
                    OrbsToCube = this->cubeOptsMC.custom_orb_request;
                }

                // SMG 04/23/25
                // Assuming 1C Real orbitals for now
                std::string nextCubes = cube_name;
                ecube->evalOrbCube(nextCubes,ess->mo[0].pointer(),NOrb,OrbsToCube,ReOrMag);

            }

            // Protonic orbitals 
            {
                size_t NB = pss->nAlphaOrbital();
                size_t NOrb = NB * ess->nC;

                // Handle which orbitals to generate
                std::vector<size_t> OrbsToCube;
                if(this->cubeOptsMC.whichMO == MO_CLASSES::ALL)
                {
                    for(size_t i = 0; i < NOrb; i++)
                        OrbsToCube.push_back(i);
                }
                // Base case is custom vector
                else
                {
                    OrbsToCube = this->cubeOptsMC.custom_orb_request;
                }

                // SMG 04/23/25
                // Assuming 1C Real orbitals for now
                std::string nextCubes = cube_name+"_PROT";
                pcube->evalOrbCube(nextCubes,pss->mo[0].pointer(),NOrb,OrbsToCube,ReOrMag);

            }

        }

    }

}; // namespace ChronusQ

