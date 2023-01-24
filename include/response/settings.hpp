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

#include <cqlinalg/cqlinalg_config.hpp>

#include <fields.hpp>
#include <response/enums.hpp>

namespace ChronusQ {

  struct ResponseSettings {

    ResponseType jobType = RESIDUE;

    // Convergence
    double convCrit = 1e-7;
    int64_t maxIter = 500;


    // Matrix Handle
    bool doFull        = true; ///< Solve Problem by (Sca)LAPACK
    bool formFullMat   = true; ///< Form The Full Matrix
    bool distMatFromRoot = false; ///< Form mat on root process and distribute
    bool formMatDist     = false; ///< Form matrix distributed

#ifdef CQ_ENABLE_MPI
    bool isDist(){ return distMatFromRoot or formMatDist; } 
#else
    bool isDist(){ return false; } 
#endif

#ifdef CQ_ENABLE_MPI
    int64_t MB = 2; ///< BLACS distribution factor
#endif

    // Matrix Properties
    bool doTDA    = false;
    bool doSA     = false;
    bool matIsHer = false;

    int  printLevel = 1;
    bool evalProp   = true;

    std::vector<ResponseOperator> aOps  = AllOps;
    std::vector<ResponseOperator> bOps  = { LenElectricDipole };

  };

  struct FDResponseSettings {

    size_t              order      = 1;
    double              dampFactor = 0.;
    bool                forceDamp  = false;
    std::vector<double> bFreq      = { 0. };

    size_t              nRHS       = 0;

    bool                needP      = false;
    bool                needQ      = false;
  };

  struct ResidueResponseSettings {

    size_t nRoots = 3;
    size_t needVR = true;
    size_t needVL = false;


    // Direct settings
    double deMin = 0.;


    // GPLHR specific settings
    size_t gplhr_m     = 3;
    double gplhr_sigma = 0.;
    bool useGDiag = true;
  };


}; // namespace ChronusQ

