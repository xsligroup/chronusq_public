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

#include <singleslater.hpp>
#include <singleslater/base.hpp>
#include <singleslater/neoss.hpp>

namespace ChronusQ {

  /*
   * Brief: Make cube files
   *
   */
  template<typename MatsT,typename IntsT>
  void SingleSlater<MatsT,IntsT> :: runCube(std::vector<std::shared_ptr<CubeGen>> cubes, EMPerturbation &emPert, std::string prefix, std::shared_ptr<Molecule> mol) {

      // Electron only for SingleSlater 
      std::shared_ptr<CubeGen> cube;
      if(this->particle.charge == 1.0)
      {
        cube = cubes[PAR_TYPE::PROTONIC];
      }
      else
      {
        cube = cubes[PAR_TYPE::ELECTRONIC]; 
      }

      if(mol){
        cube->updateMolAndBasis(mol);
      }

      std::string cube_name;
      if(cubeOptsSS.cubeFileName.empty() and prefix.empty()) {
        cube_name = "SCF";
      } else {
        cube_name = cubeOptsSS.cubeFileName;
        if(prefix.empty())
          cube_name = cube_name + "_SCF";
        else
          cube_name = cube_name + prefix;
      }

      if(this->particle.charge == 1.0)
      {
        cube_name += "_PROT";
      }

      // density cube 
      if (cubeOptsSS.denCube) {

        cube->evalDenCube(cube_name,this->onePDM);

      }

      // Orbital Cubes
      if (cubeOptsSS.orbCube) {

        size_t NB = this->nAlphaOrbital();
        size_t NOrb = NB * this->nC;

        // Handle which orbitals to generate
        std::vector<size_t> OrbsToCube;
        if(cubeOptsSS.whichMO == MO_CLASSES::ALL)
        {
          for(size_t i = 0; i < NOrb; i++)
            OrbsToCube.push_back(i);
        }
        // Base case is custom vector
        else
        {
          OrbsToCube = cubeOptsSS.custom_orb_request;
        }

        // Iterate through the different possibilities of 
        // orbital type (real/complex, alpha/beta, large/small) 
        // Control flow generates the same cubes as would be
        // printed out (with same naming scheme) as
        // wavefunction/print.hpp

        // Check if complex (for naming files)
        bool is_complex = std::is_same<MatsT,dcomplex>::value;

        // Functions passed to CubeGen
        std::function<double(MatsT)> ReOrMag;
        std::function<double(MatsT)> ImOrPhase;
        // Check if user requested Magnitude and phase rather than
        // real and imaginary
        bool MagPhase = cubeOptsSS.MagnitudeAndPhase;
        if(MagPhase)
        {
          ReOrMag = [](MatsT x){return std::abs(x);};
          ImOrPhase = [](MatsT x){return std::arg(x);};
        }
        else
        {
          ReOrMag = [](MatsT x){return std::real(x);};
          ImOrPhase = [](MatsT x){return std::imag(x);};
        }

        // Large Alpha real/mag
        // Always evaluated
        {
          std::string nextCubes = cube_name;
          if(this->nC >= 2 || ! this->iCS)
            nextCubes += "_ALPHA";
          if(this->nC == 4)
            nextCubes += "_LARGE";
          if(is_complex)
            nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL"; 
        
          cube->evalOrbCube(nextCubes,this->mo[0].pointer(),NOrb,OrbsToCube,ReOrMag);
        }

        // Large Alpha imag/phase
        if(is_complex)
        {
          std::string nextCubes = cube_name;
          if(this->nC >= 2 || ! this->iCS)
            nextCubes += "_ALPHA";
          if(this->nC == 4)
            nextCubes += "_LARGE";
          if(is_complex)
            nextCubes += MagPhase ? "_PHASE" : "_IMAG"; 

          cube->evalOrbCube(nextCubes,this->mo[0].pointer(),NOrb,OrbsToCube,ImOrPhase);
        }

        // Beta pieces
        // Also prevents Proton NEO cubes from generating beta files!
        if((this->nC >= 2 || ! this->iCS) && this->particle.charge != 1.0)
        {
          // Large Beta real/mag
          if(this->nC == 1)
          {
            std::string nextCubes = cube_name + "_BETA";
            if(is_complex)
              nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL"; 
            cube->evalOrbCube(nextCubes,this->mo[1].pointer(),NOrb,OrbsToCube,ReOrMag); 
          }
          else
          {
            std::string nextCubes = cube_name + "_BETA";
            if(this->nC == 4)
              nextCubes += "_LARGE";
            if(is_complex)
              nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL"; 
            cube->evalOrbCube(nextCubes,this->mo[0].pointer()+(this->nC/2)*NB,NOrb,OrbsToCube,ReOrMag);
          }

          if(is_complex)
          {
            if(this->nC == 1)
            {
              std::string nextCubes = cube_name + "_BETA";
              nextCubes += MagPhase ? "_PHASE" : "_IMAG"; 
              cube->evalOrbCube(nextCubes,this->mo[1].pointer(),NOrb,OrbsToCube,ImOrPhase); 
            }
            else
            {
              std::string nextCubes = cube_name + "_BETA";
              if(this->nC == 4)
                nextCubes += "_LARGE";
              if(is_complex)
                nextCubes += MagPhase ? "_PHASE" : "_IMAG"; 
              cube->evalOrbCube(nextCubes,this->mo[0].pointer()+(this->nC/2)*NB,NOrb,OrbsToCube,ImOrPhase);
            }


          }
          
        }
      

        // Small pieces go here should they be implemented
        //if(this->nC == 4)
        // Small Alpha real/mag

        // Small Alpha imag/phase

        // Small Beta real/mag

        // Small Beta imag/phase

      }

  }

  template<typename MatsT,typename IntsT>
  void NEOSS<MatsT,IntsT> :: runCube(std::vector<std::shared_ptr<CubeGen>> cubes, EMPerturbation &emPert, std::string prefix, std::shared_ptr<Molecule> mol) {

      applyToEach([&](SubSSPtr & ss){
        auto SS = std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ss);
        // Copy the SingleSlater options
        SS->cubeOptsSS = this->cubeOptsSS;
        SS->runCube(cubes,emPert,prefix,mol);});
      return;
  }


}; // namespace ChronusQ