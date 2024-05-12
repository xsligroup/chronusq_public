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
#include <cubegen.hpp>

namespace ChronusQ {

    /**
     * @brief Generates the dimensions of voxel grid
     * 
     * This uses the user input of coarse, medium, or fine
     * in the RES option, and uses that to calculate
     * the dimensions of the voxel grid, creating correct spacing
     * and quality of the output
    */
    void CubeGen::calculateVoxelDimensions() {

      std::array<double,3> maxDimensions = {0.0,0.0,0.0};
      std::map<RES_TYPE,double> resolution = {
          {RES_TYPE::COARSE, 1.0/3},
          {RES_TYPE::MEDIUM, 1.0/6},
          {RES_TYPE::FINE, 1.0/12}
      };
      // set steps based on resolution
      voxelUnits_ = {resolution[res_], resolution[res_], resolution[res_]};
      
      // get array of the greatest dimensions to create
      // total grid lengths
      for( auto &atom : mol_->atoms ) {
        for(int i = 0; i < std::size(maxDimensions); i++) {
          maxDimensions[i] = std::max(maxDimensions[i], std::abs(atom.coord[i]));
        }
      }

      std::array<double,3> maxGrid = {0,0,0};
      // create the correct amount of space for molecule
      for(int i = 0; i < std::size(maxGrid); i++) {
        // add padding
          maxGrid[i] = std::abs(maxDimensions[i]) * 2 + cubePadding_;
      }
      // using largest dimensions, convert space to the voxel spaces
      for(int i = 0; i < std::size(maxDimensions); i++) {
          voxelGrid_[i] = std::ceil(maxGrid[i] * 2 / voxelUnits_[i]);
      }

    }
    
    /**
     * @brief Calculates the 3d coordinates of center of system
     * 
     * This uses the inputs of particles to calculate the center
     * point of the system. This is not currently used, but could
     * have useful future application
     * @returns a vector of the center of system
    */
    std::vector<double> CubeGen::calcCenter() {
      std::vector<double> maxDimensions = {0.0,0.0,0.0};

      for( auto &atom : mol_->atoms ) {
        for(int i = 0; i < 3; i++) {
          maxDimensions[i] += atom.coord[i];
        }
      }

      size_t num_atoms = mol_->atoms.size();

      std::vector<double> centerPoint = {maxDimensions[0] / num_atoms,
                                           maxDimensions[1] / num_atoms,
                                           maxDimensions[2] / num_atoms};

      return centerPoint;
    }

    double * CubeGen::EvalShellSetAtPoint(int ix, int iy, int iz)
    {
        size_t NB = basis_->nBasis;

        // If we precomputed the basis we can return the cached basis
        if(evaluated_basis.size()==NB * voxelGrid_[0] * voxelGrid_[1] * voxelGrid_[2])
        {
            return &(evaluated_basis.data()[NB*(ix + iy * voxelGrid_[0] + iz * voxelGrid_[0] * voxelGrid_[1])]);
        }

        std::array<double,3> pt = {
            (ix-(int)voxelGrid_[0]/2) * voxelUnits_[0],
            (iy-(int)voxelGrid_[1]/2) * voxelUnits_[1],
            (iz-(int)voxelGrid_[2]/2) * voxelUnits_[2]
            };
        // Possibly add a parameter here for treading of the basis evaluation
        // with default value 0
        evalShellSet(NOGRAD,basis_->shells,&pt[0],1,&evaluated_basis[0],false);
        return &evaluated_basis[0];
    }

    void CubeGen::ComputeBasis()
    {
      // Guard against multiply trying to call this function
      if(evaluated_basis.size())
        return;
      
      size_t NB = basis_->nBasis;
      size_t TotalBasisSize = NB * voxelGrid_[0] * voxelGrid_[1] * voxelGrid_[2];

      try
      {
        std::cout << "Requested Total Cube Basis Size: " << TotalBasisSize << std::endl;
        evaluated_basis.resize(TotalBasisSize);
      }
      catch(...)
      {
        std::cout << " *** Attempted to precompute basis evaluation for CubeGen"  << std::endl
                  << "     but didn't have sufficient memory allocated.  Basis"   << std::endl
                  << "     will be re-computed for every cube generated!    ***"  << std::endl;
        // Still need maintain capacity for a single basis evaluation so
        // a double * can always be returned to the start of the memory
        // In the future can multiply by number of threads so in the 
        // future the cube evaluation can still be omp parallelized
        evaluated_basis.resize(NB);// * GetNumThreads());
        return;
      }

      // Calculate the amount of scratch space needed for evalShellSet
      size_t SCRSize = 0;
      // Might be able to batch points in the future, for now point by point over the grid evaluation
      size_t npts = 1;
      size_t nShSize = basis_->shells.size();
      // r contribution
      SCRSize += 3 * npts * nShSize;
      // rSq contribution
      SCRSize += npts * nShSize;

      int LMax = 0;
      for(auto iSh = 0; iSh < nShSize; iSh++)
      {
        LMax=std::max(basis_->shells[iSh].contr[0].l,LMax);
      }
      size_t shSizeCar = ((LMax+1)*(LMax+2))/2;

      // SCR_Car Contribution
      SCRSize +=  shSizeCar;

      // Allocate enough memory 
      size_t nThreads = GetNumThreads();

      double * SCR =  CQMemManager::get().malloc<double>(nThreads*SCRSize);
      double * tSCR;

      int ix,iy,iz;

#pragma omp parallel default(shared) private(tSCR, ix, iy, iz)
{
      auto iThread = GetThreadID();
      tSCR = SCR + iThread * SCRSize;
      for(iz = iThread; iz < voxelGrid_[2]; iz+=nThreads){
        for(iy = 0l; iy < voxelGrid_[1]; iy++){
          for(ix = 0l; ix < voxelGrid_[0]; ix++)
      {
        std::array<double,3> pt = {
            (ix-(int)voxelGrid_[0]/2) * voxelUnits_[0],
            (iy-(int)voxelGrid_[1]/2) * voxelUnits_[1],
            (iz-(int)voxelGrid_[2]/2) * voxelUnits_[2]
            };
        size_t offset = NB * (ix + iy * voxelGrid_[0] + iz * voxelGrid_[0] * voxelGrid_[1]);

        evalShellSet(NOGRAD,basis_->shells,&pt[0],1,&(evaluated_basis.data()[offset]),false,tSCR);
      } // iz
      } // iy
      } // ix
}

      // Free Scratch Memory
      CQMemManager::get().free(SCR);
      
      return;
    }

}
