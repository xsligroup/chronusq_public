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
#define __DEBUGTPB__

#pragma once

#include <orbitalmodifiernew/realtimeSCF.hpp>
#include <cxxapi/output.hpp>
#include <physcon.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <matrix.hpp>

#include <util/matout.hpp>
#include <util/timer.hpp>
#include <unsupported/Eigen/MatrixFunctions>


namespace ChronusQ {

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::computeTau() {

  ROOT_ONLY(this->mpiComm);
  
  // Handle type conversion
  if (typeid(this->singleSlaterSystem) != typeid(NEOSS<MatsT,IntsT>)) 
    CErr("Tau term is only implemented for NEO protons");
  NEOSS<MatsT,IntsT>& neoss = dynamic_cast<NEOSS<MatsT,IntsT>&>(this->singleSlaterSystem);

  // call function only for protonic subsystem
  SingleSlater<MatsT, IntsT>& pss = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(neoss.getSubSSBase(std::string("Protonic")))));
  pss.computeTau();

}


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::addTauToFock() {

  ROOT_ONLY(this->mpiComm);
  
  // Handle type conversion
  if (typeid(this->singleSlaterSystem) != typeid(NEOSS<MatsT,IntsT>)) 
    CErr("Tau term is only implemented for NEO protons");
  NEOSS<MatsT,IntsT>& neoss = dynamic_cast<NEOSS<MatsT,IntsT>&>(this->singleSlaterSystem);
  
  // call function only for protonic subsystem
  SingleSlater<MatsT, IntsT>& pss = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(neoss.getSubSSBase(std::string("Protonic")))));
  pss.addTauToFock();

}


/**
 *  \brief Form the adjoint of the unitary propagator
 *
 *  \f[
 *    U = \exp\left( -i \delta t F \right)
 *      = \exp\left( -\frac{i\delta t}{2}
 *                    \left(F^S \otimes I_2 + F^k \sigma_k\right) \right)
 *      = \frac{1}{2}U^S \otimes I_2 + \frac{1}{2} U^k \otimes \sigma_k
 *  \f]
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::formPropagatorForAll(std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fockSquareAO) {

  ROOT_ONLY(this->mpiComm);

  ProgramTimer::tick("Propagator Formation");
  this->ao2orthoFock(fockSquareAO);
  for( size_t i = 0; i < this->fockSquareOrtho.size(); i++ ) {
    size_t NB = this->fockSquareOrtho[i].dimension();
    MatExp('D',NB,dcomplex(0.,-integrationProgress.currentDeltaT),
           this->fockSquareOrtho[i].pointer(),NB,unitarySquareOrtho[i].pointer(),NB);
  }
  ProgramTimer::tock("Propagator Formation");
#if 0
  prettyPrintSmart(std::cout,"U",unitarySquareOrtho[i].pointer(),NB,NB,NB);
#endif
}; // RealTime::formPropagatorForAll

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::propagateDenForAll() {

  ROOT_ONLY(this->mpiComm);

  ProgramTimer::tick("Propagate Density");

  for( size_t i = 0; i < unitarySquareOrtho.size(); i++ ) {
    size_t NB = this->fockSquareOrtho[i].dimension();
    cqmatrix::Matrix<MatsT> SCR(NB);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, dcomplex(1.),
               unitarySquareOrtho[i].pointer(), NB,
               this->previousOnePDMSquareOrtho[i].pointer(), NB, dcomplex(0.), SCR.pointer(), NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, dcomplex(1.), SCR.pointer(), NB,
               unitarySquareOrtho[i].pointer(), NB, dcomplex(0.), this->previousOnePDMSquareOrtho[i].pointer(), NB);

  }

  ProgramTimer::tock("Propagate Density");

}; // RealTime::propagateDenForAll


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::unitaryProgatationForAll(std::vector<cqmatrix::Matrix<MatsT>>& onePDMSquareOrthoSave, bool startMMUTStep, bool finalMMUTStep) {

  if(tdSCFOptions.includeTau) addTauToFock();

  // Form F, U, and propagate P to next step using U (in orthonormal basis)
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fock_k = this->singleSlaterSystem.getFock();
  formPropagatorForAll(fock_k);
  propagateDenForAll();
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Perform Explicit Magnus 2 when
  //   - Start or Finish MMUT and the restart algorithm is set to be Magnus2
  //   - The integration algorithm for each step is set to be Magnus2
  if( (tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTModifiedMidpoint and (finalMMUTStep or startMMUTStep) and tdSCFOptions.restartAlgorithm == RestartAlgorithm::ExplicitMagnus2) 
      or tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2) {
    if(MPIRank(this->mpiComm) == 0){
      for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
        this->onePDMSquareOrtho[i] = this->previousOnePDMSquareOrtho[i];
        this->previousOnePDMSquareOrtho[i] = onePDMSquareOrthoSave[i];
      }
    }
    this->formFock(false, integrationProgress.currentTime + tdSCFOptions.deltaT); // F(k+1)
    // For traveling basis, ddding in tau term into protonic fock matrix
    if(tdSCFOptions.includeTau) addTauToFock();
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fock_k1 = this->singleSlaterSystem.getFock();
    for( size_t i = 0; i < fock_k.size(); i++ )
      *fock_k[i] = 0.5 * (*fock_k[i] + *fock_k1[i]); // compute 0.5 * (F(k) + F(k+1))
    formPropagatorForAll(fock_k);
    propagateDenForAll();
    if(MPIRank(this->mpiComm) == 0)
      for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) 
        this->onePDMSquareOrtho[i] = onePDMSquareOrthoSave[i];

    MPI_Barrier(MPI_COMM_WORLD);
  }  // End 2nd order magnus
} // RealTimeSCF<singleSlaterT,MatsT,IntsT>::unitaryProgatationForAll

  /**
   *  Performs density propagation using the specified algorithm(s).
   *  For NEO calculations, electronic/protonic subsystem can be propagated using different algorithms
   *
   *  \param [in] onePDMSquareOrthoSave copy of current orthonormal density that will be used in Magnus2 propagation 
   *  \param [in] startMMUTStep         whether this is the first MMUT step, in which case we do Magnus2 propagation
   *  \param [in] finalMMUTStep         whether this is the final MMUT step, in which case we do Magnus2 propagation
   * 
   *  Upon entry: previousOnePDMSquareOrtho stores the density matrix/matrices that need to be propagated
   *  Upon exit:  previousOnePDMSquareOrtho stores the density matrix/matrices that have been propagated to the next timestep
   */ 
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::doPropagation(std::vector<cqmatrix::Matrix<MatsT>>& onePDMSquareOrthoSave, bool startMMUTStep, bool finalMMUTStep) {

  EMPerturbation pert_tp5 = tdEMPerturbation.getPert(integrationProgress.currentTime);
  EMPerturbation pert_t1 = tdEMPerturbation.getPert(integrationProgress.currentTime);

  // Propagtion for non-NEO calculations:
  if (typeid(this->singleSlaterSystem) != typeid(NEOSS<MatsT,IntsT>)) {

    // RK4 Propagation
    if ( tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTRungeKuttaOrderFour ){
      this->singleSlaterSystem.setOnePDMOrtho(previousOnePDMSquareOrtho.data());
      this->singleSlaterSystem.RK4Propagation(false, integrationProgress.currentDeltaT, false, pert_tp5, pert_t1);
      previousOnePDMSquareOrtho = this->singleSlaterSystem.getOnePDMOrtho();
    } else {  
    // Unitary Propagation (MMUT or Magnus2)
      unitaryProgatationForAll(onePDMSquareOrthoSave, startMMUTStep, finalMMUTStep);    
    }

  } else {
  // Propagation for NEO calculations:
 
    // Unitary Propagation For All Density (when neither of subsystems uses RK4 )
    if ( tdSCFOptions.integrationAlgorithm == tdSCFOptions.protIntegrationAlgorithm 
        and (tdSCFOptions.integrationAlgorithm != RealTimeAlgorithm::RTRungeKuttaOrderFour) ){
      //std::cout << "NEO unitaryProgatationForAll" << std::endl;
      unitaryProgatationForAll(onePDMSquareOrthoSave, startMMUTStep, finalMMUTStep);    
    } else {
    // Separate Propagation For Each Density 

      // Transfer densities we want to propagate from RealTimeSCF object to each SingleSlater object
      NEOSS<MatsT,IntsT>& neoss = dynamic_cast<NEOSS<MatsT,IntsT>&>(this->singleSlaterSystem);
      SingleSlater<MatsT, IntsT>& ess = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(neoss.getSubSSBase(std::string("Electronic")))));
      SingleSlater<MatsT, IntsT>& pss = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(neoss.getSubSSBase(std::string("Protonic")))));
      neoss.setOnePDMOrtho(previousOnePDMSquareOrtho.data());
      
      // Propagate protonic subsystem density
      if(tdSCFOptions.protIntegrationAlgorithm == RealTimeAlgorithm::RTForwardEuler or tdSCFOptions.protIntegrationAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2){
        bool doMagnus2 = (tdSCFOptions.protIntegrationAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2);
        pss.unitaryPropagation(tdSCFOptions.includeTau, integrationProgress.currentDeltaT, doMagnus2, pert_t1);
      } else {
        pss.RK4Propagation(tdSCFOptions.includeTau, integrationProgress.currentDeltaT, false, pert_tp5, pert_t1);
      }

      // Propagate electronic subsystem density
      if(tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTForwardEuler or tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2){
        bool doMagnus2 = (tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2);
        ess.unitaryPropagation(false, integrationProgress.currentDeltaT, doMagnus2, pert_t1);
      } else {
        ess.RK4Propagation(false, integrationProgress.currentDeltaT, false, pert_tp5, pert_t1);
      }
      
      // Put updated pssOnePDMOrtho and essOnePDMOrtho back into previousOnePDMSquareOrtho
      previousOnePDMSquareOrtho = neoss.getOnePDMOrtho();
    }
  }
} // RealTimeSCF<singleSlaterT,MatsT,IntsT>::doPropagation



}; // namespace ChronusQ

