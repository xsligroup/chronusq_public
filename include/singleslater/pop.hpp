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
    const size_t nDen = this->onePDM->nComponent();
    MatsT* SCR  = CQMemManager::get().malloc<MatsT>(nDen*NB*NB);
    std::fill_n(SCR,nDen*NB*NB,MatsT(0.));

    // Molecule object to use
    Molecule inputMol = this->molecule();

    // For protonic SS, charge analysis are done for only proton atoms 
    if (this->particle.charge > 0)  inputMol = inputMol.retainQNuc();

    // Mulliken population analysis
    mullikenCharges.clear();
    if (nC != 4) {
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->S().pointer(),NB,MatsT(0.),SCR,NB);

      if( this->onePDM->hasZ() )
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Z().pointer(),NB,MatsT(0.),SCR+DENSITY_TYPE::MZ*NB*NB,NB);
      if( this->onePDM->hasXY() ){
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Y().pointer(),NB,MatsT(0.),SCR+DENSITY_TYPE::MY*NB*NB,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->X().pointer(),NB,MatsT(0.),SCR+DENSITY_TYPE::MX*NB*NB,NB);
      }

    }
    else {

      // TODO: 4c values need verified.
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->kinetic->pointer(),NB,
           this->onePDM->S().pointer()+2*NB*NB+NB,2*NB,MatsT(0.),SCR,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->S().pointer(),2*NB,MatsT(1./(2.*SpeedOfLight*SpeedOfLight)),SCR,NB);

      if( this->onePDM->hasZ() ){
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->kinetic->pointer(),NB,
           this->onePDM->Z().pointer()+2*NB*NB+NB,2*NB,MatsT(0.),SCR+DENSITY_TYPE::MZ*NB*NB,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Z().pointer(),2*NB,MatsT(1./(2.*SpeedOfLight*SpeedOfLight)),SCR+DENSITY_TYPE::MZ*NB*NB,NB);
      }

      if( this->onePDM->hasXY() ){
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->kinetic->pointer(),NB,
           this->onePDM->Y().pointer()+2*NB*NB+NB,2*NB,MatsT(0.),SCR+DENSITY_TYPE::MY*NB*NB,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Y().pointer(),2*NB,MatsT(1./(2.*SpeedOfLight*SpeedOfLight)),SCR+DENSITY_TYPE::MY*NB*NB,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->kinetic->pointer(),NB,
           this->onePDM->X().pointer()+2*NB*NB+NB,2*NB,MatsT(0.),SCR+DENSITY_TYPE::MX*NB*NB,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->X().pointer(),2*NB,MatsT(1./(2.*SpeedOfLight*SpeedOfLight)),SCR+DENSITY_TYPE::MX*NB*NB,NB);
      }

    }

    // Building lookup for basis function angular momentum
    std::vector<size_t> bfAngMom(NB);
    size_t nShell = this->basisSet().shells.size();

    for (size_t iShell = 0; iShell < nShell; iShell++) {

      size_t bfSt = this->basisSet().mapSh2Bf[iShell];
      size_t shellSize = this->basisSet().shells[iShell].size();
      size_t L = this->basisSet().shells[iShell].contr[0].l;

      for (size_t k = 0; k < shellSize; ++k) {
        bfAngMom[bfSt + k] = L;
      }

    }

    // Loop over atoms
    for(auto iAtm = 0; iAtm < inputMol.nAtoms; iAtm++) {

      // Initializing everything needed for components and l
      mullikenCharges.emplace_back();
      auto &mullikenAtom = mullikenCharges.back();
      mullikenAtom.init(this->basisSet().maxL+1);

      size_t iEnd;
      if( iAtm == inputMol.nAtoms-1 )
        iEnd = NB;
      else
        iEnd = this->basisSet().mapCen2BfSt[iAtm+1];

      size_t iSt = this->basisSet().mapCen2BfSt[iAtm];

      if (this->particle.charge < 0)
        mullikenAtom.add(DENSITY_TYPE::SCALAR,inputMol.atoms[iAtm].nucCharge);
      else
        mullikenAtom.add(DENSITY_TYPE::SCALAR,std::real(0.));
      for(auto i = iSt; i < iEnd; i++){

        size_t bfL = bfAngMom[i];

        // The scalar total density is scaled by -1 because output is Z-q
        mullikenAtom.addL(DENSITY_TYPE::SCALAR,bfL,-1.0,1.0,(-1.0 * this->particle.charge) * std::real(SCR[i*(NB+1)]));
        if( this->onePDM->hasZ() )
          mullikenAtom.addL(DENSITY_TYPE::MZ,bfL,1.0,1.0,(-1.0 * this->particle.charge) * std::real(SCR[DENSITY_TYPE::MZ*NB*NB+i*(NB+1)]));
        if( this->onePDM->hasXY() ){
          mullikenAtom.addL(DENSITY_TYPE::MY,bfL,1.0,1.0,(-1.0 * this->particle.charge) * std::real(SCR[DENSITY_TYPE::MY*NB*NB+i*(NB+1)]));
          mullikenAtom.addL(DENSITY_TYPE::MX,bfL,1.0,1.0,(-1.0 * this->particle.charge) * std::real(SCR[DENSITY_TYPE::MX*NB*NB+i*(NB+1)]));
        }
      }
    } 

    // Lowdin population analysis
    lowdinCharges.clear();

    for(auto iAtm = 0; iAtm < inputMol.nAtoms; iAtm++) {

      // Initializing everything needed for components and l
      lowdinCharges.emplace_back();
      auto &lowdinAtom = lowdinCharges.back();
      // L decomposition NYI
      lowdinAtom.init(0);

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

      if (this->particle.charge < 0)
        lowdinAtom.add(DENSITY_TYPE::SCALAR,inputMol.atoms[iAtm].nucCharge);
      else
        lowdinAtom.add(DENSITY_TYPE::SCALAR,std::real(0.));
      for(auto i = iSt; i < iEnd; i++)
        lowdinAtom.add(DENSITY_TYPE::SCALAR,-1.0*(-1.0 * this->particle.charge) * std::real(this->onePDMOrtho->S()(i,i)));
    } 

    CQMemManager::get().free(SCR);

  };

}; // namespace ChronusQ


