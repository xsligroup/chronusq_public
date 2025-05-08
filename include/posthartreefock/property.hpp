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

#include <configinteraction.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

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
  * \brief Store the AO representation of the PDM in bin
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. save SS->onePDM
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::saveOnePDMs(size_t i) {
    
    ROOT_ONLY(this->comm);

    std::cout << "Saving PDM for State " << i+1 << " as requested." << std::endl;

    // transform oneRDM to AO basis
    rdm2pdm(*this->oneRDM[i]);

    // Convert to AO basis and update PDM in ref
    std::string rdmStr = "POSTHF/RDM-"+std::to_string(i+1);
    savFile.safeWriteData(rdmStr, *ref_->onePDM);

  }; // PostHartreeFock::saveOnePDMs(i)

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::saveOnePDMs() {
   
    for (auto i : this->saveOnePDM_states) {

      PostHartreeFock::saveOnePDMs(i);

    }

  }; // PostHartreeFock::saveOnePDMs()

  /* Active space spin overlap matrices
 
     S_kl = C_k^\dagger S C_l

     stored as

     Saa, Sab, Sba, Sbb

     where 
     k,l are spin components;
     C_k is the k-spin component MO coefficients
     S_kl is the kl-spin overlap matrix
     S is the overlap matrix

  */
  template <typename MatsT, typename IntsT>
  std::vector<cqmatrix::Matrix<MatsT>> PostHartreeFock<MatsT,IntsT>::spinOverlap() {

    size_t nC = ref_->nC;
    size_t NBasis = ref_->nAlphaOrbital();
    size_t nCorrO = this->corrSpace.nCorrO;
    size_t nInact = this->corrSpace.nInact;

    std::vector<cqmatrix::Matrix<MatsT>> spin_overlap;

    cqmatrix::Matrix<MatsT> S = ref_->aoints_->overlap->matrix();
    cqmatrix::Matrix<MatsT> CMOa(NBasis, nCorrO);
    cqmatrix::Matrix<MatsT> CMOb(NBasis, nCorrO);
    cqmatrix::Matrix<MatsT> scratch(nCorrO, NBasis);
   
    if(nC == 1) {
      SetMat('N', NBasis, NBasis, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nInact, NBasis, CMOa.pointer(), NBasis);
      SetMat('N', NBasis, NBasis, MatsT(1.0), ref_->moCoefficients[1].get().pointer()+NBasis*nInact, NBasis, CMOb.pointer(), NBasis);
    }

    else if(nC == 2) {
      SetMat('N', NBasis*nC, NBasis*nC, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nC*nInact, NBasis*nC, CMOa.pointer(), NBasis);
      SetMat('N', NBasis*nC, NBasis*nC, MatsT(1.0), ref_->moCoefficients[0].get().pointer()+NBasis*nC*nInact+NBasis, NBasis*nC, CMOb.pointer(), NBasis);
    }

    else {
      CErr("4-Component spin expectation values NYI - Instance 1");
    }

    if (nC == 1 || nC == 2) {
      spin_overlap = std::vector<cqmatrix::Matrix<MatsT>>(4, cqmatrix::Matrix<MatsT>(nCorrO));

      // Computing aa/ab overlaps
      // Ca^\dagger * S
       blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOa.pointer(), NBasis, 
                 S.pointer(), NBasis, 
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);

      // Saa = Ca^\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[0].pointer(), nCorrO);

      // Sab = Ca^\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[1].pointer(), nCorrO);

      // Computing bb/ba overlaps
      // Cb^\dagger * S
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
                 nCorrO, NBasis, NBasis, 
                 MatsT(1.0), 
                 CMOb.pointer(), NBasis, 
                 S.pointer(), NBasis, 
                 MatsT(0.0), 
                 scratch.pointer(), nCorrO);

      // Sba = Cb\dagger S * Ca
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOa.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[2].pointer(), nCorrO);

      // Sbb = Cb\dagger S * Cb
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
                 nCorrO, nCorrO, NBasis, 
                 MatsT(1.0), 
                 scratch.pointer(), nCorrO, 
                 CMOb.pointer(), NBasis, 
                 MatsT(0.0), 
                 spin_overlap[3].pointer(), nCorrO);
    }

    else {
      CErr("4-Component spin expectation values NYI - Instance 2");
    }

    return spin_overlap;

  };



  /*
  * \brief Spin analysis for each state
  *         1. Transform 1RDM back to AO basis, copy to SS->onePDM
  *         2. SingleSlater->populationAnalysis()
  *
  */
  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAnalysis(size_t i, std::vector<cqmatrix::Matrix<MatsT>>* spin_overlap) {

    MatsT ssq;
    MatsT sz = MatsT(0.0);
    MatsT ssq0 = MatsT(0.0);
    MatsT ssqx = MatsT(0.0);
    MatsT ssqy = MatsT(0.0);
    MatsT ssqz = MatsT(0.0);
    size_t nCorrO = this->corrSpace.nCorrO;   
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);

    twoRDMSOI->clear(); 
    
    std::cout << std::endl << "Spin Analysis for State " << i+1 << ": " << std::endl;

    ConfigurationInteraction<MatsT,IntsT>* ci = dynamic_cast<ConfigurationInteraction<MatsT,IntsT>*>(this);
    if (ci == nullptr) {
      CErr("Dyanmic cast from PostHartreeFock to ConfigurationInteraction failed.");
    }

    ci->compute2TDM(i, i, twoRDMSOI);

    #pragma omp declare reduction (+ : std::complex<double> : omp_out += omp_in) initializer(omp_priv = std::complex<double>(0, 0))
    #pragma omp parallel for collapse(4) schedule(static) default(shared) reduction(+:ssq0, sz, ssqx, ssqy, ssqz)
    for (auto p = 0ul; p < nCorrO; p++) {
      for (auto q = 0ul; q < nCorrO; q++) {
        for (auto r = 0ul; r < nCorrO; r++) {
          for (auto s = 0ul; s < nCorrO; s++) {
          
            MatsT actual2RDM;

            if (r == 0 and s == 0) {
              // sum_k <S_k^2> = (3/4) sum_pq ( gamma_pq (Saa+Sbb)_pq )
              ssq0 += (*this->oneRDM[i])(p,q) * (((*spin_overlap)[0])(p,q) + ((*spin_overlap)[3])(p,q));
   
              // <Sz> = (1/2) sum_pq ( gamma_pq (Saa-Sbb)_pq )
              sz += (*this->oneRDM[i])(p,q) * (((*spin_overlap)[0])(p,q) - ((*spin_overlap)[3])(p,q));
            }

            if (q == r) {
              actual2RDM = (*twoRDMSOI)(p,q,r,s) - (*this->oneRDM[i])(p,s);
            }

            else {
              actual2RDM = (*twoRDMSOI)(p,q,r,s);
            }

            // <S_k(1)S_k(2)> = sgn(k) (1/4) sum_pqrs ( Gamma_pqrs S_pq S_rs )
            ssqx += (((*spin_overlap)[1])(p,q) + ((*spin_overlap)[2])(p,q)) * actual2RDM * (((*spin_overlap)[1])(r,s) + ((*spin_overlap)[2])(r,s));
            ssqy += (((*spin_overlap)[1])(p,q) - ((*spin_overlap)[2])(p,q)) * actual2RDM * (((*spin_overlap)[1])(r,s) - ((*spin_overlap)[2])(r,s));
            ssqz += (((*spin_overlap)[0])(p,q) - ((*spin_overlap)[3])(p,q)) * actual2RDM * (((*spin_overlap)[0])(r,s) - ((*spin_overlap)[3])(r,s));

          }
        }
      }
    }

    ssq = MatsT(0.25)*ssqx + MatsT(0.25)*ssqz - MatsT(0.25)*ssqy + MatsT(0.75)*ssq0;

    std::cout << "<Sz> : " << MatsT(0.5)*sz << std::endl;
    std::cout << "<S^2> : " << ssq << std::endl << std::endl; 

    // for debugging
    #if 0 
    std::cout << "State : " << i+1 << "   Sx^2 : " << MatsT(0.25)*ssqx << std::endl;
    std::cout << "State : " << i+1 << "   Sy^2 : " << MatsT(-0.25)*ssqy << std::endl;
    std::cout << "State : " << i+1 << "   Sz^2 : " << MatsT(0.25)*ssqz << std::endl;
    std::cout << "State : " << i+1 << "   S0^2 : " << MatsT(0.75)*ssq0 << std::endl;
    #endif



  }; // PostHartreeFock::spinAnalysis

  template <typename MatsT, typename IntsT>
  void PostHartreeFock<MatsT,IntsT>::spinAnalysis() {

    spin_overlap = spinOverlap();

    for (auto i = 0ul; i < this->NStates; i++) {

      PostHartreeFock::spinAnalysis(i, &spin_overlap);

    }

  }; 




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
      auto MOdipole = moints->getIntegral<VectorInts, MatsT>("MOdipole");

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

        moints->addIntegral("MOdipole", MOdipole_scr);
      }

      MOdipole = moints->getIntegral<VectorInts,MatsT>("MOdipole");

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


