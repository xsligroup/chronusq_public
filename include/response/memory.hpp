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

#include <response/tbase.hpp>

namespace ChronusQ {

  template <typename T>
  void ResponseTBase<T>::allocResidueResults() {

    bool isRoot = MPIRank(comm_) == 0;

    if( isRoot and genSettings.printLevel > 0 )
      std::cout << "  * ALLOCATING RESIDUE RESPONSE RESULTS\n";

    nSingleDim_ = getNSingleDim(genSettings.doTDA);
    size_t N = nSingleDim_;

    size_t nVec = genSettings.doFull ? nSingleDim_ : resSettings.nRoots;

    //std::cerr << " NV " << nVec << std::endl;

    resResults.W  = CQMemManager::get().malloc<double>(nVec);

    if( isRoot ) { // Don't allocate the vector results on non-root process
      resResults.VR = resSettings.needVR ? 
        CQMemManager::get().malloc<T>(nVec*N) : nullptr;
      resResults.VL = resSettings.needVL ? 
        CQMemManager::get().malloc<T>(nVec*N) : nullptr;
    }


  };

  template <typename T>
  void ResponseTBase<T>::allocFDRResults() {

    bool isRoot = MPIRank(comm_) == 0;
    if( isRoot and genSettings.printLevel > 0 )
      std::cout << "  * ALLOCATING FDR RESPONSE RESULTS\n";

    std::vector<ResponseOperator> ops = genSettings.bOps;

    fdrSettings.nRHS = 0;
    for(auto &op : ops) fdrSettings.nRHS += OperatorSize[op];

    nSingleDim_ = getNSingleDim(genSettings.doTDA);

    size_t nOmega = fdrSettings.bFreq.size();
    size_t nRHSN  = fdrSettings.nRHS * nSingleDim_;

    if( not isRoot ) return; // Don't allocate the results on non-root process

    // RHS allocated delagated to formRHS
    if( fdrSettings.dampFactor == 0. and not fdrSettings.forceDamp ) {

      fdrResults.SOL = CQMemManager::get().malloc<T>(nRHSN*nOmega);

    } else {

      dfdrResults.SOL = 
        CQMemManager::get().malloc<dcomplex>(nRHSN*nOmega);

    }

  };
};


