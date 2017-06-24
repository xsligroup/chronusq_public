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

#include <singleslater.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::populationAnalysis() {

    const size_t NB = this->basisSet().nBasis;
    MatsT* SCR  = CQMemManager::get().malloc<MatsT>(NB*NB);
    std::fill_n(SCR,NB*NB,MatsT(0.));

    // Molecule object to use
    Molecule inputMol = this->molecule();

    // For protonic SS, charge analysis are done for only proton atoms 
    if (this->particle.charge > 0)  inputMol = inputMol.retainQNuc();

    // Mulliken population analysis
    mullikenCharges.clear();
    if (nC != 4) {
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->S().pointer(),NB,MatsT(0.),SCR,NB);
    }
    else {
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->kinetic->pointer(),NB,
           this->onePDM->S().pointer()+2*NB*NB+NB,2*NB,MatsT(0.),SCR,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->S().pointer(),2*NB,MatsT(1./(2.*SpeedOfLight*SpeedOfLight)),SCR,NB);
    }

    for(auto iAtm = 0; iAtm < inputMol.nAtoms; iAtm++) {

      size_t iEnd;
      if( iAtm == inputMol.nAtoms-1 )
        iEnd = NB;
      else
        iEnd = this->basisSet().mapCen2BfSt[iAtm+1];

      size_t iSt = this->basisSet().mapCen2BfSt[iAtm];

      mullikenCharges.emplace_back(inputMol.atoms[iAtm].nucCharge);
      for(auto i = iSt; i < iEnd; i++)
        mullikenCharges.back() -= std::real(SCR[i*(NB+1)]);
    } 


    // Lowdin population analysis
    lowdinCharges.clear();

    for(auto iAtm = 0; iAtm < inputMol.nAtoms; iAtm++) {

      size_t iEnd;
      if( iAtm == inputMol.nAtoms-1 ){
        if( nC == 4 )
          iEnd = 2*NB;
	else
          iEnd = NB;
      }
      else
        iEnd = this->basisSet().mapCen2BfSt[iAtm+1];

      size_t iSt = this->basisSet().mapCen2BfSt[iAtm];

      lowdinCharges.emplace_back(inputMol.atoms[iAtm].nucCharge);
      for(auto i = iSt; i < iEnd; i++)
        lowdinCharges.back() -= std::real(this->onePDMOrtho->S()(i,i));
    } 


    CQMemManager::get().free(SCR);


  };

}; // namespace ChronusQ


