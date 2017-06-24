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

#include <orbitalmodifiernew.hpp>
#include <orbitalmodifiernew/orbitaloptimizernew/impl.hpp>
#include <orbitalmodifiernew/conventionalSCFnew/impl.hpp>
//#include <orbitalmodifiernew/newtonRaphsonSCFnew/impl.hpp>
#include <orbitalmodifiernew/realtimeSCF/impl.hpp>

namespace ChronusQ {

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
  void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::ao2orthoFock(std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fockSquareAO) {

  ROOT_ONLY(this->mpiComm);
  if(fockSquareAO.empty()) fockSquareAO = this->singleSlaterSystem.getFock();
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho = this->singleSlaterSystem.getOrtho();

  for( size_t i = 0; i < fockSquareAO.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) {
      fockSquareOrtho[i] = ortho[i]->nonortho2ortho(*fockSquareAO[i]);
    } else {
      fockSquareOrtho[i] = *(fockSquareAO[i]);
    }
  }
};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
  void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::ao2orthoDen(std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> tempOnePDMSquareAO) {

  ROOT_ONLY(this->mpiComm);
  if(tempOnePDMSquareAO.empty()) tempOnePDMSquareAO = this->singleSlaterSystem.getOnePDM();
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho = this->singleSlaterSystem.getOrtho();

  for( size_t i = 0; i < tempOnePDMSquareAO.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) {
      onePDMSquareOrtho[i] = ortho[i]->ortho2nonortho(*tempOnePDMSquareAO[i]);
    } else {
      onePDMSquareOrtho[i] = *(tempOnePDMSquareAO[i]);
    }
  }
};


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::ortho2aoDen(std::vector<cqmatrix::Matrix<MatsT>> tempOnePDMSquareOrtho) {

  ROOT_ONLY(this->mpiComm);
  //if(tempOnePDMSquareOrtho.empty()) tempOnePDMSquareOrtho = this->singleSlaterSystem.getOnePDMOrtho();
  std::vector<std::shared_ptr<Orthogonalization<MatsT>>> ortho = this->singleSlaterSystem.getOrtho();

  for( size_t i = 0; i < tempOnePDMSquareOrtho.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) {
      onePDMSquareAO[i] = ortho[i]->nonortho2ortho(tempOnePDMSquareOrtho[i]);
    } else {
      onePDMSquareAO[i] = tempOnePDMSquareOrtho[i];
    }
  }
};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::diagOrthoFock() {
  ROOT_ONLY(this->mpiComm);
  for( size_t i = 0; i < this->singleSlaterSystem.moCoefficients.size(); i++ ) {
    size_t NB = this->singleSlaterSystem.moCoefficients[i].get().dimension();
    if(NB != fockSquareOrtho[i].dimension() ) CErr("Ortho Fock and MO dimensions do not match");
    std::copy_n(fockSquareOrtho[i].pointer(), NB * NB, this->singleSlaterSystem.moCoefficients[i].get().pointer());
    int INFO  = HermetianEigen('V', 'L', NB, this->singleSlaterSystem.moCoefficients[i].get().pointer(), NB, this->singleSlaterSystem.moEigenvalues[i]);
    if( INFO != 0 ) {
      std::cout << "Attempted to diagonalize " << i << "the Fock Matrix" << std::endl;
      CErr("HermetianEigen failed in Fock", std::cout);
    }
  }
};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::ortho2aoMOs() {
  ROOT_ONLY(this->mpiComm);
  vecShrdPtrOrtho<MatsT> ortho = this->singleSlaterSystem.getOrtho();

  for( size_t i = 0; i < this->singleSlaterSystem.moCoefficients.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) { ortho[i]->ortho2nonorthoCoeffs(this->singleSlaterSystem.moCoefficients[i]); }
  }
};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::ao2orthoMOs() {
  ROOT_ONLY(this->mpiComm);
  vecShrdPtrOrtho<MatsT> ortho = this->singleSlaterSystem.getOrtho();

  for( size_t i = 0; i < this->singleSlaterSystem.moCoefficients.size(); i++ ) {
    if( ortho[i]->hasOverlap() ) { ortho[i]->nonortho2orthoCoeffs(this->singleSlaterSystem.moCoefficients[i]); }
  }
};

};   // Namespace ChronusQ
