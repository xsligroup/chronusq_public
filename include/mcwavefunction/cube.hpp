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

  /*
   * Brief: Make cube files
   *
   */
  template<typename MatsT,typename IntsT>
  void MCWaveFunction<MatsT,IntsT> :: runCube(std::vector<std::shared_ptr<CubeGen>> cubes, EMPerturbation &emPert) {

      SingleSlater<MatsT,IntsT> * ss_ptr = &reference();

      // Currently no NEO-CI
      auto cube = cubes[PAR_TYPE::ELECTRONIC]; 

      std::string cube_name;
      if(cubeOptsMC.cubeFileName.empty()) {
        cube_name = "MCWFN";
      } else {
        cube_name = cubeOptsMC.cubeFileName;
        cube_name = cube_name + "_MCWFN";
      }

      // density cube 
      if (cubeOptsMC.denCube) {

        // Update 1PDM
        // TODO: Add excited states
        rdm2pdm(this->oneRDM[0]);

        cube->evalDenCube(cube_name,ss_ptr->onePDM);

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

