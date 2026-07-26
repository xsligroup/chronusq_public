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

#include <fockbuilder/interparticlefock.hpp>

#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/gradints/direct.hpp>
#include <particleintegrals/gradints/incore.hpp>

#include <dft.hpp>
#include <dft/util.hpp>
#include <grid/integrator.hpp>
#include <util/threads.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void InterParticleFockBuilder<MatsT,IntsT>::formInterParticleCoulomb(
    SingleSlater<MatsT,IntsT>& ss, bool increment)
  {
    // Validate internal pointers
    if( this->aux_ss == nullptr )
      CErr("aux_ss uninitialized in formInterParticleCoulomb!");
    if( this->outMat == nullptr )
      CErr("outMat uninitialized in formInterParticleCoulomb!");
    if( contraction == nullptr )
      CErr("contraction uninitialized in formInterParticleCoulomb!");

    // Decide onePDM to use
    cqmatrix::PauliSpinorMatrices<MatsT>& contract1PDM
        = increment ? *this->aux_ss->deltaOnePDM : *this->aux_ss->onePDM;

    if(not increment)
      this->outMat->clear();

    std::vector<TwoBodyContraction<MatsT>> contract =
      { {contract1PDM.S().pointer(), this->outMat->pointer(), true, COULOMB} };

    EMPerturbation pert;

    auto beginContract = tick();
    contraction->twoBodyContract(ss.comm, contract, pert);

    if(contraction->printContractionTiming)
      std::cout << "        " << std::left << std::setw(38)
                << "InterParticle-Coulomb-Contraction duration = "
                << tock(beginContract) << " s\n" << std::endl;
  }

  template <typename MatsT, typename IntsT>
  std::vector<double> InterParticleFockBuilder<MatsT,IntsT>::formInterParticleCoulombGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& pert, double xHFX) {

    if(this->aux_ss == nullptr)
      CErr("aux_ss uninitialized in formInterParticleCoulombGrad!");
    if(gradTPI == nullptr)
      CErr("gradTPI uninitialized in formInterParticleCoulombGrad!");

    const size_t NB = ss.basisSet().nBasis;
    const size_t nGrad = 3*ss.molecule().nAtoms;

    // Form contraction
    std::unique_ptr<GradContractions<MatsT,IntsT>> contract;
    if(std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>((*gradTPI)[0]))
      contract = std::make_unique<InCore4indexGradContraction<MatsT,IntsT>>(*gradTPI);
    else if(std::dynamic_pointer_cast<DirectTPI<IntsT>>((*gradTPI)[0]))
      contract = std::make_unique<DirectGradContraction<MatsT,IntsT>>(*gradTPI);
    else
      CErr("Gradients of RI inter-particle integrals NYI!");

    // Assume that the order of the TPI is the same between the regular and
    //   gradient integrals
    contract->contractSecond = contraction->contractSecond;
    contract->isCross = true;            // Important for GIAO JContract
    contract->traceDensity = ss.onePDM;  // Determines the density to screen in Direct AO two-e gradient formation 

    // Create contraction list
    std::vector<std::vector<TwoBodyContraction<MatsT>>> cList;

    std::vector<cqmatrix::Matrix<MatsT>> JList;
    JList.reserve(nGrad);

    for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {
      std::vector<TwoBodyContraction<MatsT>> tempCont;

      // Coulomb
      JList.emplace_back(NB);
      JList.back().clear();
      tempCont.push_back(
         {this->aux_ss->onePDM->S().pointer(), JList.back().pointer(), true, COULOMB}
      );

      cList.push_back(tempCont);
    }

    // Contract to J/K
    contract->gradTwoBodyContract(MPI_COMM_WORLD, true, cList, pert);

    // Contract to gradient
    std::vector<double> gradient;
    cqmatrix::PauliSpinorMatrices<MatsT> twoEGrad(NB, false, false);
    const double prefactor = 2. * ss.particle.charge * this->aux_ss->particle.charge;
    for(size_t iGrad = 0; iGrad < nGrad; iGrad++) {
      twoEGrad.S() = prefactor * JList[iGrad];
      double gradVal = ss.template computeOBProperty<SCALAR>(
        twoEGrad.S().pointer()
      );
      gradient.push_back(0.25*gradVal);
    }

    return gradient;

  }

  template <typename MatsT, typename IntsT>
  void InterParticleFockBuilder<MatsT,IntsT>::formFock(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& empert, bool increment,
    double xHFX)
  {
    if( this->upstream == nullptr )
      CErr("Upstream FockBuilder uninitialized in formFock!");
    if( this->aux_ss == nullptr )
      CErr("aux_ss uninitialized in formFock!");

    // Call all upstream FockBuilders (intra-particle + any other interactions)
    this->upstream->formFock(ss, empert, increment, xHFX);

    formInterParticleCoulomb(ss, increment);

    ROOT_ONLY(ss.comm);

    // Compute correct sign for interparticle Coulomb interaction
    // Factor of 2 comes from our convention, it will be scaled in SingleSlater::computeEnergy
    double prefactor = 2. * ss.particle.charge * this->aux_ss->particle.charge;
    *ss.twoeH      += prefactor * *this->outMat;
    *ss.fockMatrix += prefactor * *this->outMat;
  }

  template <typename MatsT, typename IntsT>
  std::vector<double> InterParticleFockBuilder<MatsT,IntsT>::getGDGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& pert, double xHFX) {
    if(this->upstream == nullptr)
      CErr("Upstream FockBuilder uninitialized in getGDGrad!");

    auto gradient = this->upstream->getGDGrad(ss,pert,xHFX);
    auto pairGradient = formInterParticleCoulombGrad(ss,pert,xHFX);
    std::transform(gradient.begin(),gradient.end(),pairGradient.begin(),
                   gradient.begin(),std::plus<double>());
    return gradient;
  };

}
