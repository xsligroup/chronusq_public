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

#include <mcwavefunction.hpp>

namespace ChronusQ {

  /*
   * Brief: Make cube files
   *
   */
  template<typename MatsT,typename IntsT>
  void MCWaveFunction<MatsT,IntsT> :: runCube(std::vector<std::shared_ptr<CubeGen>> cubes) {

      SingleSlater<MatsT,IntsT> * ss_ptr = &reference();

      // Currently no NEO-CI: the first cube is always the electronic cube
      auto cube = cubes[0]; 

      std::string cube_name;
      if(cubeOptsMC.cubeFileName.empty()) {
        cube_name = "MCWFN";
      } else {
        cube_name = cubeOptsMC.cubeFileName;
        cube_name = cube_name + "_MCWFN";
      }

      // density cube 
      if (cubeOptsMC.denCube) {
        
        // Handle which roots to generate
        size_t NRoots = this->NStates;
        std::vector<size_t> RootsToCube;
        if(cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::GS)
        {
          RootsToCube.push_back(0);
        }
        else if(cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::ALL)
        {
          for(size_t i = 0; i < NRoots; i++)
            RootsToCube.push_back(i);
        }
        else if(cubeOptsMC.whichCIRoots == CI_CUBE_ROOT_CLASSES::AVERAGE)
        {
          std::shared_ptr<cqmatrix::Matrix<MatsT>> SARDM = std::make_shared<cqmatrix::Matrix<MatsT>>(this->oneRDM[0].nRows());
          SARDM->clear();
          // Build the state averaged density
          for(size_t i = 0; i < NRoots; i++)
          {
            *SARDM += this->SAWeight[i] * this->oneRDM[i];
          }
          rdm2pdm(*SARDM);
          std::string sa_cube_name = cube_name+"_SA";
          cube->evalDenCube(sa_cube_name,ss_ptr->onePDM);
        }
        else
        {
          RootsToCube=cubeOptsMC.custom_root_request;
        }
        
        // Update 1PDM
        if(cubeOptsMC.whichCIRoots != CI_CUBE_ROOT_CLASSES::AVERAGE)
        {
          for(const auto Root : RootsToCube)
          {
            if(Root > NRoots)
              CErr("Requesting cube for Root #" + std::to_string(Root+1) + " which is greater than CI Roots available!");
            rdm2pdm(this->oneRDM[Root]);
            std::string cube_name_with_root = cube_name+"_ROOT_"+std::to_string(Root+1);
            cube->evalDenCube(cube_name_with_root,ss_ptr->onePDM);
          }
        }

      }

      // Orbital Cubes
      if (cubeOptsMC.orbCube) {

        size_t NB = ss_ptr->nAlphaOrbital();
        size_t NOrb = NB * ss_ptr->nC;

        // Handle which orbitals to generate
        std::vector<size_t> OrbsToCube;
        if(cubeOptsMC.whichMO == MO_CLASSES::ALL)
        {
          for(size_t i = 0; i < NOrb; i++)
            OrbsToCube.push_back(i);
        }
        // Base case is custom vector
        else
        {
          OrbsToCube = cubeOptsMC.custom_orb_request;
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
        bool MagPhase = cubeOptsMC.MagnitudeAndPhase;
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
          if(ss_ptr->nC >= 2 || ! ss_ptr->iCS)
            nextCubes += "_ALPHA";
          if(ss_ptr->nC == 4)
            nextCubes += "_LARGE";
          if(is_complex)
            nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL";

          cube->evalOrbCube(nextCubes,ss_ptr->mo[0].pointer(),NOrb,OrbsToCube,ReOrMag);
        }

        // Large Alpha imag/phase
        if(is_complex)
        {
          std::string nextCubes = cube_name;
          if(ss_ptr->nC >= 2 || ! ss_ptr->iCS)
            nextCubes += "_ALPHA";
          if(ss_ptr->nC == 4)
            nextCubes += "_LARGE";
          if(is_complex)
            nextCubes += MagPhase ? "_PHASE" : "_IMAG";

          cube->evalOrbCube(nextCubes,ss_ptr->mo[0].pointer(),NOrb,OrbsToCube,ImOrPhase);
        }

        // Beta pieces
        if(ss_ptr->nC >= 2 || ! ss_ptr->iCS)
        {
          // Large Beta real/mag
          if(ss_ptr->nC == 1)
          {
            std::string nextCubes = cube_name + "_BETA";
            if(is_complex)
              nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL";
            cube->evalOrbCube(nextCubes,ss_ptr->mo[1].pointer(),NOrb,OrbsToCube,ReOrMag);
          }
          else
          {
            std::string nextCubes = cube_name + "_BETA";
            if(ss_ptr->nC == 4)
              nextCubes += "_LARGE";
            if(is_complex)
              nextCubes += MagPhase ? "_MAGNITUDE" : "_REAL";
            cube->evalOrbCube(nextCubes,ss_ptr->mo[0].pointer()+(ss_ptr->nC/2)*NB,NOrb,OrbsToCube,ReOrMag);
          }

          if(is_complex)
          {
            if(ss_ptr->nC == 1)
            {
              std::string nextCubes = cube_name + "_BETA";
              nextCubes += MagPhase ? "_PHASE" : "_IMAG";
              cube->evalOrbCube(nextCubes,ss_ptr->mo[1].pointer(),NOrb,OrbsToCube,ImOrPhase);
            }
            else
            {
              std::string nextCubes = cube_name + "_BETA";
              if(ss_ptr->nC == 4)
                nextCubes += "_LARGE";
              if(is_complex)
                nextCubes += MagPhase ? "_PHASE" : "_IMAG";
              cube->evalOrbCube(nextCubes,ss_ptr->mo[0].pointer()+(ss_ptr->nC/2)*NB,NOrb,OrbsToCube,ImOrPhase);
            }

          }

        }

        // Small pieces go here should they be implemented
        //if(ss_ptr->nC == 4)
        // Small Alpha real/mag

        // Small Alpha imag/phase

        // Small Beta real/mag

        // Small Beta imag/phase
      }
  }

}; // namespace ChronusQ

