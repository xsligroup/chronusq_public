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

#include <mcwavefunction.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

 /*
  * \brief Perform mulliken analysis for each state
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. SingleSlater->populationAnalysis()
  *
  */
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::populationAnalysis(size_t i) {

    std::cout << std::endl << "Population Analysis for State " << i+1 << ": ";
    SingleSlater<MatsT,IntsT> * ss_ptr = &reference();

    // transform oneRDM to AO basis
    rdm2pdm(this->oneRDM[i]);

    ss_ptr->populationAnalysis();
    ss_ptr->printMiscProperties(std::cout);

  }; // MCWaveFunction::populationAnalysis

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::populationAnalysis() {

    for (auto i = 0ul; i < this->NStates; i++) {

      MCWaveFunction::populationAnalysis(i);

    }

  }; // MCWaveFunction::populationAnalysis


 /*
  * \brief Compute oscillator strength for MC wavefunction
  *         using AO dipole and MO coefficients and MO TDM
  *         Only for 1C and 2C
  *         s1: initial state
  *         s2: final state
  */ 
  template <typename MatsT, typename IntsT>
  double MCWaveFunction<MatsT,IntsT>::oscillator_strength(size_t s2, size_t s1) {

    if (referenceWaveFunction().nC == 4) CErr("4C has no dipole.");

    auto  &mem = memManager;
    size_t nAO = reference().nAlphaOrbital() * reference().nC;;
    size_t nCorrO = MOPartition.nCorrO;
    size_t nInact = MOPartition.nInact;

    // compute transition density matrix for specific state
    SquareMatrix<MatsT> tmpTDM1(mem,nCorrO);
    SquareMatrix<MatsT> tmpTDM2(mem,nCorrO);
    ciBuilder->computeTDM(*this, CIVecs[s1], CIVecs[s2], tmpTDM1);
    ciBuilder->computeTDM(*this, CIVecs[s2], CIVecs[s1], tmpTDM2);

    MatsT D = MatsT(0.);

    // dipole AO -> MO transformation
    auto MOdipole = moints.getIntegral<VectorInts,MatsT>("MOdipole");

    if (not MOdipole) {
      std::shared_ptr<VectorInts<IntsT>> AOdipole = 
                std::make_shared<VectorInts<IntsT>>(mem, nAO, 1, true);
      std::shared_ptr<VectorInts<MatsT>> MOdipole_scr =
                std::make_shared<VectorInts<MatsT>>(mem, nCorrO, 1, true);

      std::vector<std::pair<size_t, size_t>> active(2, {MOPartition.nFCore+nInact, nCorrO});

      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        if (referenceWaveFunction().nC == 1)
                 (*AOdipole)[iXYZ] = (*reference().aoints.lenElectric)[iXYZ];
        else if (referenceWaveFunction().nC == 2)
          (*AOdipole)[iXYZ] = (*reference().aoints.lenElectric)[iXYZ]
                                .template spatialToSpinBlock<IntsT>();
        (*AOdipole)[iXYZ].subsetTransform('N',reference().mo[0].pointer(),
                nAO, active, (*MOdipole_scr)[iXYZ].pointer(), false);
      }

      moints.addIntegral("MOdipole", MOdipole_scr);
    }

    MOdipole = moints.getIntegral<VectorInts,MatsT>("MOdipole");

    // dipole strength D = Tr(TDM \dot MOdiple) Tr(TDM^* \dot MOdipole)
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      D += blas::dot(nCorrO*nCorrO,tmpTDM1.pointer(),1,(*MOdipole)[iXYZ].pointer(),1)
          *blas::dot(nCorrO*nCorrO,tmpTDM2.pointer(),1,(*MOdipole)[iXYZ].pointer(),1);
    }

    // oscillator strength f = 2/3 (E2 - E1) D.
    double f = (2./3.) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(D);

    // output
    std::cout << "Excited State: " << std::setw(3) << std::right << s2+1
              << " to state: " << std::setw(3) << std::right << s1+1 << ":";
    std::cout << std::setw(15) << std::right << "E(Eh) = "
              << std::setprecision(8) << std::fixed << (StateEnergy[s2] - StateEnergy[s1]);
    std::cout << std::setw(15) << std::right << "f = "
              << std::setprecision(6) << std::fixed << f << std::endl;

    tmpTDM1.clear();
    tmpTDM2.clear();

    return f;

  } // MCWaveFunction::oscillator_strength








}; // namespace ChronusQ


