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

#include <cqlinalg/eig.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/ortho.hpp>

#include <itersolver.hpp>

#include <util/matout.hpp>
#include <util/timer.hpp>

namespace ChronusQ {

  template <typename T>
  void ResponseTBase<T>::runFullResidue() {

    ROOT_ONLY(comm_);

    if( this->genSettings.printLevel > 0 )
      std::cout << "PERFORMING THE FULL DIAGONALIZATION\n";

    T* fMatUse = formFullFromMemory();


    //std::cerr << "Top of runFullResidue\n";
    char JOBVR = resSettings.needVR ? 'V' : 'N';
    char JOBVL = resSettings.needVL ? 'V' : 'N';

    ProgramTimer::tick("Full Diagonalize");

    if( genSettings.matIsHer ) {

      //std::cerr << "Hermetian Diag\n";
      //std::cerr << " FM " << fMatUse << std::endl;

      T * mat = fMatUse;
      if( resSettings.needVR )  {
        std::copy_n(fMatUse,nSingleDim_*nSingleDim_,resResults.VR);
        mat = resResults.VR;
      }

      HermetianEigen(JOBVR,'U',nSingleDim_,mat,nSingleDim_,resResults.W);

    } else {

      //std::cerr << "NonHermetian Diag\n";
      //std::cerr << " FM " << fMatUse << std::endl;

      dcomplex *W = 
        CQMemManager::get().malloc<dcomplex>(nSingleDim_);

      GeneralEigen(JOBVL, JOBVR, nSingleDim_, fMatUse, nSingleDim_,
                   W, resResults.VL, nSingleDim_, resResults.VR, nSingleDim_);
      
      for(auto k = 0; k < nSingleDim_; k++)
        resResults.W[k] = std::real(W[k]);

      CQMemManager::get().free(W);

      //CErr();

    }

    ProgramTimer::tock("Full Diagonalize");


    /*
    // Ensure proper orthogonality of eigenvectors
    eigVecNorm();


    if( fMatUse != fullMatrix_ ) CQMemManager::get().free(fMatUse);

    // Write the eigensystem to disk
    if( savFile.exists() ) {

      savFile.safeWriteData("/RESP/RESIDUE/EIGENVALUES",resResults.W,
        {resSettings.nRoots});
      
      if( resSettings.needVR )
      savFile.safeWriteData("/RESP/RESIDUE/EIGENVECTORS",resResults.VR,
        {resSettings.nRoots,nSingleDim_});

    }

    */

    if( fMatUse != fullMatrix_ ) CQMemManager::get().free(fMatUse);

  };

  template <typename T>
  void ResponseTBase<T>::runIterResidue() {

    if((not genSettings.isDist()) and genSettings.formFullMat) 
      ROOT_ONLY(comm_);

    bool isRoot = MPIRank(comm_) == 0;
    bool isDist = this->genSettings.isDist();


    typename GPLHR<T>::LinearTrans_t lt = [&](size_t nVec, SolverVectors<T> &V, SolverVectors<T> &AV) {

      iterLinearTrans(nVec,V,AV);

    };

    ProgramTimer::tick("Iter Diagonalize");


    typename GPLHR<T>::Shift_t pc =
      bool(PC_) ? PC_ :
      [&](size_t nVec, T shift, SolverVectors<T> &V, SolverVectors<T> &AV) {

      //if( not this->fullMatrix_ ) CErr();
      std::cout << "in PC of GPLHR" << std::endl;
      AV.set_data(0, nVec, V, 0);
    };


    MPI_Comm gplhrComm = (isDist or not genSettings.formFullMat) 
      ? comm_ : rcomm_;
    
    GPLHR<T> gplhr(gplhrComm,nSingleDim_,
      genSettings.maxIter,genSettings.convCrit,
      resSettings.nRoots,lt,pc);

    gplhr.setM(resSettings.gplhr_m);
    gplhr.sigma   = resSettings.gplhr_sigma;
    gplhr.hardLim = resSettings.deMin;

    if( hasResGuess_ )
      gplhr.setGuess(resSettings.nRoots,
          [&](size_t nG, SolverVectors<T> & G, size_t LDG){ this->resGuess(nG,G,LDG); });



    gplhr.run();

    if( isRoot ) {

      // Extract data from GPLHR storage
      for(auto k = 0; k < resSettings.nRoots; k++)
        resResults.W[k] = std::real(gplhr.eigVal()[k]);

      std::copy_n(tryGetRawVectorsPointer(*gplhr.VR()),
                  this->nSingleDim_ * resSettings.nRoots,
                  resResults.VR);

    }

    ProgramTimer::tock("Iter Diagonalize");

  };

}; // namespace ChronusQ

