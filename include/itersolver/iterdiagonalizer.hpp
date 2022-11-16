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

#include <itersolver.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  template <typename _F>
  void IterDiagonalizer<_F>::run() {

    if( MPISize(this->comm_) > 1 ) {

      MPIBCast(nRoots_,0,this->comm_);
      MPIBCast(this->N_,0,this->comm_);
      MPIBCast(this->mSS_,0,this->comm_);

    }

    bool isRoot = MPIRank(this->comm_) == 0;
    this->mSS_ = std::min(this->mSS_,this->N_);
    this->nGuess_ = std::min(this->nGuess_,this->N_);
    
    if( isRoot ) {
      std::cout << "\n  * IterDiagonalizer will solve for " << nRoots_ 
        << " eigenroots\n\n";

      std::cout << std::left;
      std::cout << std::setw(30) << "    * Problem Dimension        = " 
                                 << this->N_ << "\n";
      std::cout << std::setw(30) << "    * Initial of Guess Vectors = " 
                                 << this->nGuess_ << "\n";
      std::cout << std::setw(30) << "    * Maximum Subspace         = " 
                                 << this->mSS_ << "\n";
      std::cout << std::setw(30) << "    * Maximum Micro Iterations = " 
                                 << this->maxMicroIter_ << "\n";
      std::cout << std::setw(30) << "    * Maximum Macro Iterations = " 
                                 << this->maxMacroIter_ << "\n";
      std::cout << std::setw(30) << "    * Residual Conv Crit       = " 
                                 << std::scientific << std::setprecision(4)
                                 << this->convCrit_ << "\n";

      std::cout << "\n\n" << std::endl;
    }

    alloc(); // Allocate Scratch space

    // Macro iterations UNFINISHED
    for(auto iMacro = 0; iMacro < this->maxMacroIter_; iMacro++) {

      this->converged_ = runMicro();
      if( this->converged_ ) break;

      restart();
    }

  }; // iterDiagonalizer::run 

}; // namespace ChronusQ

