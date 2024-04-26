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
     * @brief Generates the cubefile for charge density
     * 
     * This is called from evalCube.
    */
    template <typename LocMatsT>
    void CubeGen::evalCDCube(std::shared_ptr<cqmatrix::PauliSpinorMatrices<LocMatsT>> oPDM_)
    {
      // gets number of basis sets
      size_t NB = ref_->nAlphaOrbital();

      std::vector<libint2::Shell> &shells = ref_->basisSet().shells;

      for(auto ix = 0l; ix < voxelGrid_[0]; ix++) {
        for(auto iy = 0l; iy < voxelGrid_[1]; iy++) {
          for(auto iz = 0l; iz < voxelGrid_[2]; iz++) {

            // voxel calculation
            std::array<double,3> pt = {
            (ix-(int)voxelGrid_[0]/2) * voxelUnits_[0],
            (iy-(int)voxelGrid_[1]/2) * voxelUnits_[1],
            (iz-(int)voxelGrid_[2]/2) * voxelUnits_[2]
            };

            std::vector<double> BASIS(NB,0.);
            std::vector<LocMatsT> SCR(NB,0.);

            LocMatsT val = 0;

            evalShellSet(NOGRAD,ref_->basisSet().shells,&pt[0],1,&BASIS[0],false);

            blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, 1, NB, NB,1.,
            &BASIS[0], NB, oPDM_->S().pointer(), NB,0., &SCR[0], 1);

            val = blas::dotu(NB,&SCR[0],1,&BASIS[0],1);

            *cubeFile_ << std::right << std::setw(15) << std::setprecision(5)
            << std::scientific << std::uppercase << std::real(ref_->particle.charge *  val);

            if( iz % 6 == 5 ) *cubeFile_ << "\n";

          }
          *cubeFile_ << "\n";
        }
      }
    };

}
