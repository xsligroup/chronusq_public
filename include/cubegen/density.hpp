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
    template <typename LocMatsT>
    void CubeGen::evalDenCompCube(LocMatsT *oPDM_)
    {
      // gets number of basis sets
      size_t NB = basis_->nBasis;

      double * BASIS;

      for(auto ix = 0l; ix < voxelGrid_[0]; ix++) {
        for(auto iy = 0l; iy < voxelGrid_[1]; iy++) {
          for(auto iz = 0l; iz < voxelGrid_[2]; iz++) {

            std::vector<LocMatsT> SCR(NB,0.);

            LocMatsT val = 0;

            BASIS = EvalShellSetAtPoint(ix,iy,iz);

            blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, 1, NB, NB,1.,
            BASIS, NB, oPDM_, NB,0., &SCR[0], 1);

            val = blas::dotu(NB,&SCR[0],1,&BASIS[0],1);

            *cubeFile_ << std::right << std::setw(20) << std::setprecision(12)
            << std::scientific << std::uppercase << std::real(particleCharge_ *  val);

            if( iz % 6 == 5 ) *cubeFile_ << "\n";

          }
          *cubeFile_ << "\n";
        }
      }
    };

    template <>
    void CubeGen::evalDenCompCube(dcomplex *oPDM_)
    {
      // gets number of basis sets
      size_t NB = basis_->nBasis;

      double * BASIS;
      dcomplex * GIAO_BASIS;

      // debug only
      //double test_value = 0;

    // GIAO
    if( basis_->basisType == COMPLEX_GIAO ) {

      std::cout << "Generating Denfile for GIAO" << std::endl;

      for(auto ix = 0l; ix < voxelGrid_[0]; ix++) {
        for(auto iy = 0l; iy < voxelGrid_[1]; iy++) {
          for(auto iz = 0l; iz < voxelGrid_[2]; iz++) {

            std::vector<dcomplex> SCR(NB,0.);

            dcomplex val = 0;

            GIAO_BASIS = EvalShellSetAtPointGIAO(ix,iy,iz);

            // ShiChao's Thesis, D.9
            // \rho = \Sigma_\munu P_\munu \Xi_\mu \Xi_\nu*

            // First Do (Chi^T*P)
            blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, 1, NB, NB, dcomplex(1.),
            GIAO_BASIS, NB, oPDM_, NB, dcomplex(0.), &SCR[0], 1);

            // Next Do (Chi^T*P)* Chi
            //val = std::real(blas::dotu(NB,&SCR[0],1,&BASIS[0],1));
            for (size_t nu = 0; nu < NB; nu++) 
              val += SCR[nu] * std::conj(GIAO_BASIS[nu]);

            //test_value += particleCharge_*std::real(val)*voxelUnits_[0]*voxelUnits_[1]*voxelUnits_[2]; 

            *cubeFile_ << std::right << std::setw(20) << std::setprecision(12)
            << std::scientific << std::uppercase << std::real(particleCharge_ *  std::real(val));

            if( iz % 6 == 5 ) *cubeFile_ << "\n";

          }
          *cubeFile_ << "\n";
        }
      }
      //*cubeFile_ << std::right << std::setw(20) << std::setprecision(12) << std::scientific << std::uppercase << test_value;

    // Complex GTO
    } else {

      for(auto ix = 0l; ix < voxelGrid_[0]; ix++) {
        for(auto iy = 0l; iy < voxelGrid_[1]; iy++) {
          for(auto iz = 0l; iz < voxelGrid_[2]; iz++) {

            std::vector<dcomplex> SCR(NB,0.);

            dcomplex val = 0;

            BASIS = EvalShellSetAtPoint(ix,iy,iz);

            blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans, 1, NB, NB,1.,
            BASIS, NB, oPDM_, NB,0., &SCR[0], 1);

            val = blas::dotu(NB,&SCR[0],1,&BASIS[0],1);

            *cubeFile_ << std::right << std::setw(20) << std::setprecision(12)
            << std::scientific << std::uppercase << std::real(particleCharge_ * std::real(val));

            if( iz % 6 == 5 ) *cubeFile_ << "\n";

          }
          *cubeFile_ << "\n";
        }
      }

    }
    };

}
