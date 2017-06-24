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

#include <chronusq_sys.hpp>
#include <cubegen.hpp>

namespace ChronusQ {

    /**
     * @brief Generates the cubefile for a density component
     * 
     * This is called from evalCube.
    */
    template <typename LocMatsT, typename ValManipOp>
    void CubeGen::evalOrbCompCube(LocMatsT * MO, size_t LDMO, size_t MOIndex, ValManipOp op)
    {
      // gets number of basis sets
      size_t NB = basis_->nBasis;

      double * BASIS;

      LocMatsT * ThisMO = MO + MOIndex * LDMO;


      // SMG Debug printing
      //std::cout << "MOIndex: " << MOIndex << std::endl;
      //std::cout << "Re, Im, Mag, Phase" << std::endl;
      //for(size_t i = 0; i < NB; i++)
      //  std::cout << std::right << std::setprecision(6) << std::setw(14) << 
      //   std::real(ThisMO[i]) <<  "   " <<
      //   std::imag(ThisMO[i]) <<  "   " << 
      //   std::abs(ThisMO[i]) <<  "   " << 
      //   std::arg(ThisMO[i]) << std::endl;

      for(auto ix = 0l; ix < voxelGrid_[0]; ix++) {
        for(auto iy = 0l; iy < voxelGrid_[1]; iy++) {
          for(auto iz = 0l; iz < voxelGrid_[2]; iz++) {

            std::vector<LocMatsT> SCR(NB,0.);

            LocMatsT val = 0;

            BASIS = EvalShellSetAtPoint(ix,iy,iz);

            val = blas::dot(NB,ThisMO,1,BASIS,1);

            *cubeFile_ << std::right << std::setw(15) << std::setprecision(5)
            << std::scientific << std::uppercase << op(val);

            if( iz % 6 == 5 ) *cubeFile_ << "\n";

          }
          *cubeFile_ << "\n";
        }
      }
    };

}
