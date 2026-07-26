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
#include <algorithm>


namespace ChronusQ {

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::computeTau() {

  // NOTE: per-subsystem computeTau (rt.hpp) computes on root and broadcasts tau
  //   to all ranks, so it must be invoked on every rank (no ROOT_ONLY here).
  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    for(const auto& label : this->singleSlaterSystem.getOrder())
      if(label != "E") this->singleSlaterSystem.getSubSS(label)->computeTau();
  } else
    CErr("Tau term is only implemented for non-electronic quantum subsystems");

}


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::addTauToFock() {

  ROOT_ONLY(this->mpiComm);

  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    for(const auto& label : this->singleSlaterSystem.getOrder())
      if(label != "E") this->singleSlaterSystem.getSubSS(label)->addTauToFock();
  } else
    CErr("Tau term is only implemented for non-electronic quantum subsystems");

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
    size_t NB = this->fockSquareOrtho[i].nRows();
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
    size_t NB = this->fockSquareOrtho[i].nRows();
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
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::unitaryProgatationForAll(
  std::vector<cqmatrix::Matrix<MatsT>>& onePDMSquareOrthoSave, RealTimeAlgorithm algorithm, bool startMMUTStep, bool finalMMUTStep) {

  if(tdSCFOptions.includeTau) addTauToFock();

  // Form F, U, and propagate P to next step using U (in orthonormal basis)
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fock_k = this->singleSlaterSystem.getFock();
  formPropagatorForAll(fock_k);
  propagateDenForAll();
  
  MPI_Barrier(MPI_COMM_WORLD);

  // Perform Explicit Magnus 2 when
  //   - Start or Finish MMUT and the restart algorithm is set to be Magnus2
  //   - The integration algorithm for each step is set to be Magnus2
  if( (algorithm == RealTimeAlgorithm::RTModifiedMidpoint and (finalMMUTStep or startMMUTStep) and tdSCFOptions.restartAlgorithm == RestartAlgorithm::ExplicitMagnus2)
      or algorithm == RealTimeAlgorithm::RTExplicitMagnus2) {
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
   *  MultiParticleSS subsystems can use different propagation algorithms.
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
  if constexpr (not std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {

    // RK4 Propagation
    if ( tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTRungeKuttaOrderFour ){
      this->singleSlaterSystem.setOnePDMOrtho(previousOnePDMSquareOrtho.data());
      this->singleSlaterSystem.RK4Propagation(false, integrationProgress.currentDeltaT, false, pert_tp5, pert_t1);
      previousOnePDMSquareOrtho = this->singleSlaterSystem.getOnePDMOrtho();
    } 
    else if (tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::ElectronicBornOppenheimer){
      CErr("Electronic Born-Oppenheimer calcualtion only valid in a NEO context!");
    }
    else if (isUnitaryRTAlgorithm(tdSCFOptions.integrationAlgorithm)) {
    // Unitary Propagation (ForwardEuler, MMUT or Magnus2)
      unitaryProgatationForAll(onePDMSquareOrthoSave, tdSCFOptions.integrationAlgorithm, startMMUTStep, finalMMUTStep);
    }
    else {
      CErr("Using an improper RT propagation algorithm");
    }

  } 
  else {
  // Propagation for NEO calculations:
    auto& multiSS = this->singleSlaterSystem;
    auto labels = multiSS.getOrder();
    if(labels.empty()) CErr("No quantum subsystems available for RT propagation");

    RealTimeAlgorithm firstAlgorithm = tdSCFOptions.subsystemIntegrationAlgorithms.at(labels.front());
    bool sameUnitaryAlgorithm = isUnitaryRTAlgorithm(firstAlgorithm) and
      std::all_of(labels.begin(), labels.end(), [&](const std::string& label) {
        return tdSCFOptions.subsystemIntegrationAlgorithms.at(label) == firstAlgorithm;
      });

    if(sameUnitaryAlgorithm) {
      unitaryProgatationForAll(onePDMSquareOrthoSave, firstAlgorithm, startMMUTStep, finalMMUTStep);
      return;
    }

    // Mixed algorithms propagate each subsystem separately. 
    // Place electronic density to last such that BORT can converge with other densities at (t+dt) instead of (t)
    std::stable_partition(labels.begin(), labels.end(), [&](const std::string& label) {
      return tdSCFOptions.subsystemIntegrationAlgorithms.at(label) != RealTimeAlgorithm::ElectronicBornOppenheimer;
    });

    multiSS.setOnePDMOrtho(previousOnePDMSquareOrtho.data());

    for(const auto& label : labels) {
      auto ss = multiSS.getSubSS(label);
      RealTimeAlgorithm algorithm = tdSCFOptions.subsystemIntegrationAlgorithms.at(label);
      double deltaT = getPropagationTimeStep(algorithm, startMMUTStep, finalMMUTStep);
      bool includeTau = tdSCFOptions.includeTau and label != "E";

      RTFockFormation formTargetFock = [&multiSS, label](EMPerturbation& pert, bool increment) {
        multiSS.formFockForTargets(pert, {label}, increment);
      };

      if(algorithm == RealTimeAlgorithm::RTForwardEuler)
        ss->unitaryPropagation(includeTau, deltaT, false, pert_t1, formTargetFock);
      else if(algorithm == RealTimeAlgorithm::RTExplicitMagnus2)
        ss->unitaryPropagation(includeTau, deltaT, true, pert_t1, formTargetFock);
      else if(algorithm == RealTimeAlgorithm::RTModifiedMidpoint) {
        bool doMagnus2 = (startMMUTStep or finalMMUTStep) and tdSCFOptions.restartAlgorithm == RestartAlgorithm::ExplicitMagnus2;
        ss->unitaryPropagation(includeTau, deltaT, doMagnus2, pert_t1, formTargetFock);
      } else if(algorithm == RealTimeAlgorithm::RTRungeKuttaOrderFour)
        ss->RK4Propagation(includeTau, deltaT, false, pert_tp5, pert_t1, formTargetFock);
      else if(algorithm == RealTimeAlgorithm::ElectronicBornOppenheimer) {
        if(label != "E") CErr("Electronic Born-Oppenheimer propagation is only valid for subsystem E");
        electronicBornOppenheimer();
      } else
        CErr("Using an improper RT propagation algorithm for subsystem " + label);
    }

    previousOnePDMSquareOrtho = multiSS.getOnePDMOrtho();
  }
} // RealTimeSCF<singleSlaterT,MatsT,IntsT>::doPropagation

  /**
   *  Performs an SCF calculation on the electronic wavefunction in the Electronic Born-Oppenheimer approximation
   */ 
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::electronicBornOppenheimer()
{
  if constexpr (not std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    CErr("Electronic Born-Oppenheimer propagation requires MultiParticleSS");
  } else {
    auto& multiSS = this->singleSlaterSystem;
    SCFControls originalControls = multiSS.scfControls;

    // Tell the MultiParticleSS object to only optimize the electronic wavefunction
    multiSS.scfControls.NEOSubSystemOpt = {"E"};
    multiSS.scfControls.NEOStepwiseOpt = false;

    // For restart jobs, energyOnly is set to True to avoid doing an unnecessary SCF calculation at the start
    // We need to overwrite that here to make sure the SCF runs
    multiSS.scfControls.scfAlg = _CONVENTIONAL_SCF;
    multiSS.scfControls.energyOnly = false;

    // Overwrite the convergence thresholds
    multiSS.scfControls.rmsdPConvTol = tdSCFOptions.BORTAccuracy;
    multiSS.scfControls.maxdPConvTol = tdSCFOptions.BORTAccuracy*100;
    multiSS.scfControls.eneConvTol   = tdSCFOptions.BORTAccuracy*100;

    // Set printing for the BORT SCF calculation
    multiSS.scfControls.printLevel = tdSCFOptions.BORTPrintLevel;

    ConventionalSCFNew<MultiParticleSS,MatsT,IntsT> conventionalSCF(multiSS.scfControls, multiSS, this->mpiComm);

    // Call SCF
    multiSS.initializeSCF();

    // Use the current field
    EMPerturbation emPert = tdEMPerturbation.getPert(integrationProgress.currentTime);
    if(tdSCFOptions.includeSCFField)
      for(auto& field : staticEMPerturbation.fields) emPert.addField(field);
    conventionalSCF.run(emPert);

    // Reset the SCF controls
    multiSS.scfControls = originalControls;

    // Print the converged SCF energy
    std::cout << "    Electronic Born-Oppenheimer Converged SCF Energy: "
              << std::fixed << std::setprecision(10) << std::right
              << multiSS.getTotalEnergy() << std::endl;
  }
}



}; // namespace ChronusQ

