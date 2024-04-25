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

#include <posthartreefock.hpp>
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
  void PostHartreeFock<MatsT,IntsT>::populationAnalysis(size_t i) {
    
    ROOT_ONLY(this->comm);

    std::cout << std::endl << "Population Analysis for State " << i+1 << ": ";

    // transform oneRDM to AO basis
    rdm2pdm(*this->oneRDM[i]);

    ref_->populationAnalysis();
    ref_->printMiscProperties(std::cout);

  }; // PostHartreeFock::populationAnalysis

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::populationAnalysis() {
    
    ROOT_ONLY(this->comm);

    for (auto i = 0ul; i < this->NStates; i++) {
      populationAnalysis(i);
    }

  }; // PostHartreeFock::populationAnalysis

 /*
  * \brief Compute oscillator strength for MC wavefunction
  *         using AO dipole and MO coefficients and MO TDM
  *         Only for 1C and 2C
  *         s1: initial state
  *         s2: final state
  */ 
  template <typename MatsT, typename IntsT>
  double PostHartreeFock<MatsT,IntsT>::oscillator_strength(size_t s2, size_t s1) {

    if (ref_->nC == 4) CErr("4C has no dipole.");

    size_t nAO = ref_->nAlphaOrbital() * ref_->nC;;
    size_t nCorrO = corrSpace.nCorrO;
    size_t nCoreO = corrSpace.nInact + corrSpace.nFCore;

    // compute transition density matrix for specific state
    auto tmpTDM1 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    auto tmpTDM2 = std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO);
    computeTDM(s1, s2, tmpTDM1);
    computeTDM(s2, s1, tmpTDM2);
     
    double f;
    if (MPIRank(this->comm) == 0) {
     
      // dipole AO -> MO transformation
      auto MOdipole = moints.getIntegral<VectorInts, MatsT>("MOdipole");

      if (not MOdipole) {
        VectorInts<IntsT> AOdipole(nAO, 1, true);
        std::shared_ptr<VectorInts<MatsT>> MOdipole_scr =
                  std::make_shared<VectorInts<MatsT>>(nCorrO, 1, true);
        
        auto corrOffs = mointsTF->parseMOType("tu");

        for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
          if (ref_->nC == 1)
             AOdipole[iXYZ] = (*ref_->aoints_->lenElectric)[iXYZ];
          else if (ref_->nC == 2)
             AOdipole[iXYZ] = std::make_shared<OnePInts<IntsT>>((*ref_->aoints_->lenElectric)[iXYZ]
                                  ->template spatialToSpinBlock<IntsT>());
          AOdipole[iXYZ]->subsetTransform('N',ref_->mo[0].pointer(),
                  nAO, corrOffs, (*MOdipole_scr)[iXYZ]->pointer(), false);
        }

        moints.addIntegral("MOdipole", MOdipole_scr);
      }

      MOdipole = moints.getIntegral<VectorInts,MatsT>("MOdipole");

      MatsT D = MatsT(0.);
      // dipole strength D = Tr(TDM \dot MOdiple) Tr(TDM^* \dot MOdipole)
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
        D += blas::dotu(nCorrO*nCorrO,tmpTDM1->pointer(),1,(*MOdipole)[iXYZ]->pointer(),1)
            *blas::dotu(nCorrO*nCorrO,tmpTDM2->pointer(),1,(*MOdipole)[iXYZ]->pointer(),1);
      }

      // oscillator strength f = 2/3 (E2 - E1) D.
      f = (2./3.) * (StateEnergy[s2] - StateEnergy[s1]) * std::real(D);

      // output
      std::cout << "\nExcited State: " << std::setw(3) << std::right << s2+1
                << " to state: " << std::setw(3) << std::right << s1+1 << ":";
      std::cout << std::setw(15) << std::right << "E(Eh) = "
                << std::setprecision(8) << std::fixed << (StateEnergy[s2] - StateEnergy[s1]);
      std::cout << std::setw(15) << std::right << "f = "
                << std::setprecision(12) << std::fixed << f << std::endl;
    }
    MPIBCast(f, 0, this->comm);

    return f;

  } // PostHartreeFock::oscillator_strength

}; // namespace ChronusQ


