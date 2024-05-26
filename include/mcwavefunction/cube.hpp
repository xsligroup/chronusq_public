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
  void MCWaveFunction<MatsT,IntsT> :: runCube(std::shared_ptr<CubeGen> cube) {

      SingleSlater<MatsT,IntsT> * ss_ptr = &reference();

      std::string cube_name;
      if(cube->getCubeFileName().empty()) {
        cube_name = "MCWFN";
      } else {
        cube_name = cube->getCubeFileName();
        cube_name = cube_name + "_MCWFN";
      }

      // density cube 
      if (cube->getDenEval()) {

        // Update 1PDM
        // TODO: Add excited states
        rdm2pdm(this->oneRDM[0]);

        cube->evalDenCube(cube_name,ss_ptr->onePDM);

      }

  }

}; // namespace ChronusQ

