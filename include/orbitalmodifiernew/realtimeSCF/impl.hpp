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

#include <orbitalmodifiernew/realtimeSCF.hpp>
#include <orbitalmodifiernew/realtimeSCF/propagation.hpp>
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
std::vector<RealTimeAlgorithm>
RealTimeSCF<singleSlaterT,MatsT,IntsT>::getDensityIntegrationAlgorithms() {

  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    std::vector<RealTimeAlgorithm> algorithms;

    // Match each flattened density block with its subsystem algorithm.
    for(const auto& label : this->singleSlaterSystem.getOrder()) {
      auto iter = tdSCFOptions.subsystemIntegrationAlgorithms.find(label);
      if(iter == tdSCFOptions.subsystemIntegrationAlgorithms.end())
        CErr("Missing RT integration algorithm for quantum subsystem " + label);
      size_t nDensities = this->singleSlaterSystem.getSubSS(label)->getOnePDM().size();
      algorithms.insert(algorithms.end(), nDensities, iter->second);
    }

    if(algorithms.size() != this->onePDMSquareOrtho.size())
      CErr("Mismatch between RT density blocks and subsystem integration algorithms");

    return algorithms;
  } else
    return std::vector<RealTimeAlgorithm>(this->onePDMSquareOrtho.size(), tdSCFOptions.integrationAlgorithm);
}


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
double RealTimeSCF<singleSlaterT,MatsT,IntsT>::getPropagationTimeStep(
  RealTimeAlgorithm algorithm, bool startMMUTStep, bool finalMMUTStep) const {

  if(algorithm == RealTimeAlgorithm::RTModifiedMidpoint and not (startMMUTStep or finalMMUTStep))
    return 2. * tdSCFOptions.deltaT;

  return tdSCFOptions.deltaT;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::run(EMPerturbation &perturbation) {

  staticEMPerturbation = perturbation;

  ProgramTimer::tick("Real Time Total");

  printRunHeader(staticEMPerturbation);

  bool startMMUTStep(false); //   startStep for the MMUT iterations (propagate delta t)
  bool finalMMUTStep(false); //   finalStep for the MMUT iterations (propagate delta t)
  bool normalMMUTStep(true); // normal step for the MMUT iterations (propagate 2 * delta t)

  // Disable RI-K-Coefficient Contraction during RT as coefficients are not updated.
  // Use density contration during RT instead
  this->singleSlaterSystem.setDenEqCoeff(false);
    
  // Evaluate D3 correction energy for current geometry
#ifdef CQ_HAS_D3
    if (this->singleSlaterSystem.d3Utils) {
      this->singleSlaterSystem.d3Utils->evaluate(this->singleSlaterSystem.molecule());
    }
#endif

  // Initialize RT or Restore 1PDM Ortho and integration process (on root process)
  // If first step,form density in AO and Ortho basis
  if ( tdSCFOptions.restoreFromStep != 0 ) {
    this->restoreState();
  } else {
    //this->singleSlaterSystem.formDensity();
    this->ao2orthoDen();  // transform and save the orthonormal density in onePDMSquareOrtho (on root process)
    this->ao2orthoFock(); // transform and save the orthonormal Fock matrix in fockSquareOrtho (on root process)
  }

  // Recompute one-electron integrals (including the metric) (on root process).
  this->singleSlaterSystem.formCoreH(perturbation,false);

  if( tdSCFOptions.includeTau ) 
    computeTau();

  // Putting current orthonormal density in previousOnePDMSquareOrtho for propagation 
  if ( MPIRank(this->mpiComm) == 0 )
    for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) 
      previousOnePDMSquareOrtho[i] = this->onePDMSquareOrtho[i];

  const auto densityIntegrationAlgorithms = getDensityIntegrationAlgorithms();
  if(densityIntegrationAlgorithms.empty())
    CErr("No density matrices available for real-time propagation");

  const bool hasMMUT = std::find(densityIntegrationAlgorithms.begin(),
    densityIntegrationAlgorithms.end(), RealTimeAlgorithm::RTModifiedMidpoint) !=
    densityIntegrationAlgorithms.end();
  const RealTimeAlgorithm uniformAlgorithm = densityIntegrationAlgorithms.front();
  const bool useUniformUnitary = isUnitaryRTAlgorithm(uniformAlgorithm) and
    std::all_of(densityIntegrationAlgorithms.begin(), densityIntegrationAlgorithms.end(),
      [&](RealTimeAlgorithm algorithm) { return algorithm == uniformAlgorithm; });

  // Recover current time for the step information
  integrationProgress.currentTime = tdSCFOptions.restoreFromStep * tdSCFOptions.deltaT;

  for(integrationProgress.currentStep = tdSCFOptions.restoreFromStep;
      integrationProgress.currentStep < tdSCFOptions.maxSteps;
      integrationProgress.currentTime += tdSCFOptions.deltaT, integrationProgress.currentStep++) {

    ProgramTimer::tick("Real Time Iter");

    // Perturbation for the current time
    EMPerturbation currentPerturbation = tdEMPerturbation.getPert(integrationProgress.currentTime);

    // "Start" the MMUT if the current step is the first step or a restart step  or if field is discontinuous
    // "Finish" the MMUT if the current step is the last step or a irestart step or if field is discontinuous
    if(hasMMUT) {
      startMMUTStep = finalMMUTStep or (integrationProgress.currentStep == tdSCFOptions.restoreFromStep );
      finalMMUTStep = (integrationProgress.currentStep == tdSCFOptions.maxSteps ) or tdEMPerturbation.isFieldDiscontinuous(integrationProgress.currentTime, tdSCFOptions.deltaT);
      if(tdSCFOptions.iRestart > 0 and (integrationProgress.currentStep + 1) % tdSCFOptions.iRestart == 0) finalMMUTStep = true;
    };

    // The global time step is used when all densities are propagated together.
    integrationProgress.currentDeltaT = getPropagationTimeStep(
      useUniformUnitary ? uniformAlgorithm : tdSCFOptions.integrationAlgorithm,
      startMMUTStep, finalMMUTStep);

    // Set up the densities used for Fock formation and propagation
    normalMMUTStep = not (finalMMUTStep or startMMUTStep);
    std::vector<cqmatrix::Matrix<MatsT>> onePDMSquareOrthoSave;
    if(MPIRank(this->mpiComm) == 0){
      if(printLevel > 0 and startMMUTStep) std::cout << "  *** Starting MMUT ***\n";
      if(printLevel > 0 and finalMMUTStep) std::cout << "  *** Finishing MMUT ***\n";

      // When every subsystem uses the same unitary algorithm, set up their magnus2 propagation:
      //   Magnus2 uses two propagations from the same starting density:
      //   1. Save P(t). F(t) is formed below from onePDMSquareOrtho.
      //   2. Propagate P(t) with F(t) to obtain a trial P(t+dt).
      //   3. Form F(t+dt) from the trial density.
      //   4. Restore P(t), average F(t) and F(t+dt), and propagate again.
      //   Here we use onePDMSquareOrthoSave preserves P(t)
      if(useUniformUnitary and
         (uniformAlgorithm == RealTimeAlgorithm::RTExplicitMagnus2 or
         (uniformAlgorithm == RealTimeAlgorithm::RTModifiedMidpoint and
          not normalMMUTStep and tdSCFOptions.restartAlgorithm == RestartAlgorithm::ExplicitMagnus2)))
        onePDMSquareOrthoSave = this->previousOnePDMSquareOrtho;

      // Set up the density containers for each subsystem's propagation:
      //   A normal MMUT step performs the two-step update
      //     P(t+dt) = U P(t-dt) U^*, where U = exp[-i 2 dt F(t)].
      //   1. previousOnePDMSquareOrtho initially stores the newest P(t), while
      //      onePDMSquareOrtho retains P(t-dt) from the preceding MMUT step.
      //   2. Swap them so onePDMSquareOrtho stores P(t), which is used below
      //      to form F(t) and the propagator U.
      //   3. The swap leaves previousOnePDMSquareOrtho storing P(t-dt), which
      //      is propagated with U to produce P(t+dt).
      //
      //   Every other algorithm, including an MMUT start or finish step, uses
      //   a one-step update starting from P(t):
      //   1. previousOnePDMSquareOrtho stores P(t).
      //   2. Copy it into onePDMSquareOrtho to form F(t).
      //   3. Leave previousOnePDMSquareOrtho at P(t) as the propagation source.
      for(size_t i = 0; i < this->onePDMSquareOrtho.size(); i++) {
        if(densityIntegrationAlgorithms[i] == RealTimeAlgorithm::RTModifiedMidpoint and normalMMUTStep)
          this->onePDMSquareOrtho[i].swap(this->previousOnePDMSquareOrtho[i]);
        else
          this->onePDMSquareOrtho[i] = this->previousOnePDMSquareOrtho[i];
      }
    }

    this->formFock(false, integrationProgress.currentTime);
    this->singleSlaterSystem.computeEnergy(currentPerturbation);
    // Add D3 correction to the total energy
#ifdef CQ_HAS_D3
    if (this->singleSlaterSystem.d3Utils) {
      this->singleSlaterSystem.totalEnergy += this->singleSlaterSystem.d3Utils->result().energy;
    }
#endif
    this->singleSlaterSystem.computeProperties(currentPerturbation, {ELECTRIC_DIPOLE});
    if( printLevel > 0 and (MPIRank(this->mpiComm) == 0) ) printIteration();

    doPropagation(onePDMSquareOrthoSave, startMMUTStep, finalMMUTStep);

    this->saveState(currentPerturbation); // Save the current energy, dipole, and propagated density every iSave steps
    this->saveCube(); // Output Cube
    ProgramTimer::tock("Real Time Iter");

  } // Time loop

  // After RT, Density is propagated to maxstep+1. 
  // Transfer maxstep+1 density from RealtimeSCF object back to SingleSlater object
  this->ortho2aoDen(this->previousOnePDMSquareOrtho);                // Transfrom density to AO 
  this->singleSlaterSystem.setOnePDMAO(this->onePDMSquareAO.data()); // Scatter to spin blocks and copy back to SingleSlater object 
  this->singleSlaterSystem.ao2orthoDen();                             

  // Allow next RT in Ehrenfest to start from the last saved density (which is propagated to maxstep+1)
  if(tdSCFOptions.doMD)  tdSCFOptions.restoreFromStep = -1; 

  ProgramTimer::tock("Real Time Total");

}; // RealTime::doPropagation


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT, IntsT>::formFock(bool increment, double time) {

  //this->singleSlaterSystem.setOnePDMOrtho(this->onePDMSquareOrtho.data());
  //this->singleSlaterSystem.checkIdempotency();

  //this->singleSlaterSystem.ortho2aoDen();

  // Transfrom density to AO (on root process)
  this->ortho2aoDen(this->onePDMSquareOrtho);
  // Scatter to spin blocks (on root process) and broacast if doing MPI
  this->singleSlaterSystem.setOnePDMAO(this->onePDMSquareAO.data());

  ProgramTimer::timeOp("Form Fock", [&]() {
      // Get perturbation for the current time and build a Fock matrix
      EMPerturbation pert_t = tdEMPerturbation.getPert(time);

      // check whether to update gaunt and gauge terms or re-use the saved matrices
      if (this->singleSlaterSystem.nC == 4){
        if ((integrationProgress.currentStep % tdSCFOptions.rtGaunt)==0){
          this->singleSlaterSystem.fockBuilder->hamiltonianOptions_.updateGaunt = true;
        } 
        else{
          this->singleSlaterSystem.fockBuilder->hamiltonianOptions_.updateGaunt = false;
        }

        if ((integrationProgress.currentStep % tdSCFOptions.rtGauge)==0) {
          this->singleSlaterSystem.fockBuilder->hamiltonianOptions_.updateGauge = true;
        } 
        else{
          this->singleSlaterSystem.fockBuilder->hamiltonianOptions_.updateGauge = false;
        } 
      } 

      // Add the SCF Perturbation
      if ( tdSCFOptions.includeSCFField ) for( auto& field : staticEMPerturbation.fields ) pert_t.addField(field );
      this->singleSlaterSystem.formFock(pert_t,increment);
  });

}; // RealTime::formFock


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::createRTDataSets(size_t maxPoints) {

  ROOT_ONLY(this->mpiComm);
  
  if(tdSCFOptions.doMD){
    //integrationProgress.maxSavePoints = (std::max(tdSCFOptions.rtMaxStepsPerMDStep/tdSCFOptions.iSave, static_cast<size_t>(1)) + 1) * tdSCFOptions.totalMDSteps;
    integrationProgress.maxSavePoints = std::max(tdSCFOptions.rtMaxStepsPerMDStep/tdSCFOptions.iSave, static_cast<size_t>(1)) * tdSCFOptions.totalMDSteps +1;
  } else{
    if( maxPoints == 0 ) integrationProgress.maxSavePoints = tdSCFOptions.maxSteps/tdSCFOptions.iSave + 2;
    else integrationProgress.maxSavePoints = maxPoints;
  }

  if(tdSCFOptions.restoreFromStep != 0 ) return;

  savFile.createGroup("RTNEW");

  savFile.createDataSet<size_t>("RTNEW/ISAVE", {1});
  savFile.createDataSet<size_t>("RTNEW/ICUBE", {1});
  savFile.createDataSet<size_t>("RTNEW/LASTSAVEPOINT", {1});
  savFile.createDataSet<size_t>("RTNEW/MAXSAVEPOINTS", {1});
  savFile.createDataSet<size_t>("RTNEW/MAXSTEPS", {1});

  size_t maxDim = tdSCFOptions.doMD ? tdSCFOptions.rtMaxStepsPerMDStep*tdSCFOptions.totalMDSteps : tdSCFOptions.maxSteps;
  savFile.createDataSet<size_t>("RTNEW/SAVESTEP", {integrationProgress.maxSavePoints});
  savFile.createDataSet<size_t>("RTNEW/STEP",   {maxDim});
  savFile.createDataSet<double>("RTNEW/TIME",   {maxDim});
  savFile.createDataSet<double>("RTNEW/ENERGY", {maxDim});
  savFile.createDataSet<double>("RTNEW/LEN_ELEC_DIPOLE",       {maxDim*3});
  savFile.createDataSet<double>("RTNEW/LEN_ELEC_DIPOLE_FIELD", {maxDim*3});

  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    for(const auto& label : this->singleSlaterSystem.getOrder()) {
      const std::string prefix = "RTNEW/" + label;
      savFile.createGroup(prefix);
      savFile.createDataSet<double>(prefix + "/LEN_DIPOLE", {maxDim*3});
    }
  }

  for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
    size_t nBasis = this->onePDMSquareOrtho[i].nRows();
    savFile.createDataSet<dcomplex>("RTNEW/TD_1PDM_ORTHO"+std::to_string(i), {integrationProgress.maxSavePoints*nBasis*nBasis});
    savFile.createDataSet<double>("RTNEW/ORBITALPOPULATION"+std::to_string(i), {integrationProgress.maxSavePoints*nBasis});
    if(tdSCFOptions.saveOnePDM)
      savFile.createDataSet<dcomplex>("RTNEW/TD_1PDM"+std::to_string(i), {integrationProgress.maxSavePoints*nBasis*nBasis});
  }

}; // RealTime::createRTDataSets

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::saveState(EMPerturbation& currentPerturbation) {

  ROOT_ONLY(this->mpiComm);

  savFile.safeWriteData("RTNEW/ISAVE", &tdSCFOptions.iSave, {1});
  savFile.safeWriteData("RTNEW/ICUBE", &tdSCFOptions.iCube, {1});
  savFile.safeWriteData("RTNEW/MAXSAVEPOINTS", &integrationProgress.maxSavePoints, {1});
  savFile.safeWriteData("RTNEW/MAXSTEPS", &tdSCFOptions.maxSteps, {1});

  integrationProgress.time.push_back(integrationProgress.currentTime);
  integrationProgress.energy.push_back(this->singleSlaterSystem.totalEnergy);
  integrationProgress.electricDipole.push_back(this->singleSlaterSystem.elecDipole);
  if( currentPerturbation.fields.size() > 0 ) integrationProgress.electricDipoleField.push_back(currentPerturbation.getDipoleAmp(Electric) );
  
  savFile.partialWriteData("RTNEW/TIME", &integrationProgress.currentTime, {integrationProgress.currentStep},{1},{0},{1});
  savFile.partialWriteData("RTNEW/STEP", &integrationProgress.currentStep, {integrationProgress.currentStep},{1},{0},{1});
  savFile.partialWriteData("RTNEW/ENERGY", &this->singleSlaterSystem.totalEnergy, {integrationProgress.currentStep},{1},{0},{1});
  savFile.partialWriteData("RTNEW/LEN_ELEC_DIPOLE", &this->singleSlaterSystem.elecDipole[0],{integrationProgress.currentStep*3}, {3},{0},{3});

  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    const auto& subsystemDipoles = this->singleSlaterSystem.getSubsystemDipoles();
    for(const auto& label : this->singleSlaterSystem.getOrder())
      savFile.partialWriteData("RTNEW/" + label + "/LEN_DIPOLE",
        &subsystemDipoles.at(label)[0], {integrationProgress.currentStep*3},
        {3}, {0}, {3});
  }

  std::array<double,3> elecDipoleField = currentPerturbation.getDipoleAmp(Electric);
  if (integrationProgress.electricDipoleField.size() > 0)
    savFile.partialWriteData("RTNEW/LEN_ELEC_DIPOLE_FIELD",&elecDipoleField[0], {integrationProgress.currentStep*3}, {3},{0},{3});

  if (integrationProgress.currentStep == tdSCFOptions.restoreFromStep   // Save on entry
      or integrationProgress.currentStep == tdSCFOptions.maxSteps-1     // Save on exit
      or integrationProgress.currentStep % tdSCFOptions.iSave == 0)     // Save every iSave step
  {

    if(integrationProgress.currentStep == tdSCFOptions.restoreFromStep
        and tdSCFOptions.restoreFromStep != 0) {
      integrationProgress.lastSavePoint++;
      return; // For the first step after restart, disable saving
    }
    
    std::cout << "  *** Saving step #"<<integrationProgress.currentStep<<"( t = "<<integrationProgress.currentTime<<" au) to binary file ***" << std::endl;
    
    savFile.safeWriteData("RTNEW/LASTSAVEPOINT", &(integrationProgress.lastSavePoint), {1});
    savFile.partialWriteData("RTNEW/SAVESTEP", &integrationProgress.currentStep, {integrationProgress.lastSavePoint},{1},{0},{1});

    for (size_t i = 0; i < this->onePDMSquareOrtho.size(); i++) {
      size_t nBasis = this->onePDMSquareOrtho[i].nRows();
      // TODO: Determine which density to save. 
      // For Ehrenfest dynamics with many RT restarts, we need to want to save the density AFTER Propagation to allow appropriate next restart
      // the current density (time k,   before propagation) is saved in onePDMSquareOrtho
      // the next    density (time k+1, after  propagation) is saved in previousOnePDMSquareOrtho
      savFile.partialWriteData("RTNEW/TD_1PDM_ORTHO" + std::to_string(i), this->previousOnePDMSquareOrtho[i].pointer(),{integrationProgress.lastSavePoint*nBasis*nBasis}, {nBasis * nBasis}, {0}, {nBasis * nBasis});
      if(tdSCFOptions.saveOnePDM){
        std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> tempOnePDM = this->singleSlaterSystem.getOnePDM();
        savFile.partialWriteData("RTNEW/TD_1PDM" + std::to_string(i), tempOnePDM[i]->pointer(),{integrationProgress.lastSavePoint*nBasis*nBasis}, {nBasis * nBasis}, {0}, {nBasis * nBasis});
      }
    }

    integrationProgress.lastSavePoint++;
  }

}; // RealTime::saveState

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::restoreState() {
   
  // Checking if savfile exists (on root process)
  int flag = 0;  // 0 means everything is okay, 1 means error out
  if (MPIRank(this->mpiComm) == 0) flag = this->savFile.exists() ?  0 : 1;
  #ifdef CQ_ENABLE_MPI
  if( MPISize(this->mpiComm) > 1 ) MPIBCast(&flag, 1, 0, this->mpiComm);
  if (flag != 0)                   MPI_Abort(this->mpiComm, 1); 
  #endif
  if (flag != 0)  CErr("SavFile not found!");

  // Restore integration progress and 1PDM ortho (on root process)
  if( MPIRank(this->mpiComm) == 0 ) {

    size_t maxSavePoints, lastSavePoint;
    savFile.readData("RTNEW/MAXSAVEPOINTS", &maxSavePoints);
    if ( maxSavePoints != integrationProgress.maxSavePoints ) CErr("Mismatched requested and saved propagation length!");
    savFile.readData("RTNEW/LASTSAVEPOINT", &lastSavePoint);

    size_t iSave;
    savFile.readData("RTNEW/ISAVE", &iSave);

    size_t restartStep;
    if (tdSCFOptions.restoreFromStep < 0) {
      integrationProgress.lastSavePoint = lastSavePoint;
      if (tdSCFOptions.doMD) {
        size_t scaleFactor = std::min(tdSCFOptions.rtMaxStepsPerMDStep, iSave);
        restartStep = integrationProgress.lastSavePoint*scaleFactor-1;
      } else {
        restartStep = integrationProgress.lastSavePoint*iSave;
      }
    } else {
      if (tdSCFOptions.doMD){
        integrationProgress.lastSavePoint = tdSCFOptions.restoreFromStep;
        size_t scaleFactor = std::min(tdSCFOptions.rtMaxStepsPerMDStep, iSave);
        restartStep = integrationProgress.lastSavePoint*scaleFactor-1;
      } else {
        integrationProgress.lastSavePoint = tdSCFOptions.restoreFromStep/iSave;
        if(integrationProgress.lastSavePoint > lastSavePoint) integrationProgress.lastSavePoint = lastSavePoint;
        restartStep = integrationProgress.lastSavePoint*iSave;
      }
    }

    // Restore time dependent density
    try {
      savFile.partialReadData("RTNEW/STEP", &integrationProgress.currentStep, {restartStep}, {1}, {0}, {1});
      savFile.partialReadData("RTNEW/TIME", &integrationProgress.currentTime, {restartStep}, {1}, {0}, {1});

      for (size_t i = 0; i < this->onePDMSquareOrtho.size(); i++) {
        size_t nBasis = this->onePDMSquareOrtho[i].nRows();
        // NOTE: Reading from previous RT maxstep+1 density (density after the last propagation)
        if(savFile.getDims("RTNEW/TD_1PDM_ORTHO" + std::to_string(i)) != std::vector<size_t>{integrationProgress.maxSavePoints*nBasis*nBasis})
          CErr("Mismatched requested and saved propagation length!");
        savFile.partialReadData("RTNEW/TD_1PDM_ORTHO" + std::to_string(i), this->onePDMSquareOrtho[i].pointer(), {integrationProgress.lastSavePoint*nBasis*nBasis}, {nBasis * nBasis}, {0}, {nBasis * nBasis});
      }
    } catch(...) { }

    if( printLevel > 0 ) {
      std::cout << "  *** Restoring from step " << integrationProgress.currentStep << " (";
      std::cout << std::setprecision(4) << integrationProgress.currentTime << " au) ***" << std::endl;
    }
    
    // Note: since density is already previous maxstep+1, we need to add one to current step
    integrationProgress.currentStep++;
    tdSCFOptions.restoreFromStep = integrationProgress.currentStep;
  }

  // Broadcast the integration process to all MPI processes
  #ifdef CQ_ENABLE_MPI
  if( MPISize(this->mpiComm) > 1 ){
    std::cout  << "  *** Scattering the integrationProgress ***\n";
    MPIBCast(&(integrationProgress.lastSavePoint), 1, 0, this->mpiComm);
    MPIBCast(&(integrationProgress.currentStep),   1, 0, this->mpiComm);
    MPIBCast(&(integrationProgress.currentTime),   1, 0, this->mpiComm);
    MPIBCast(&(tdSCFOptions.restoreFromStep),      1, 0, this->mpiComm);
  }
  #endif

}; // RealTime::restoreState


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::saveCube() {

  ROOT_ONLY(this->mpiComm);

  if (tdSCFOptions.iCube > 0) {
  if (integrationProgress.currentStep % tdSCFOptions.iCube == 0) {

    std::cout << tdSCFOptions.iCube << std::endl;
    std::cout << "  *** Saving density of step #"<<integrationProgress.currentStep<<"( t = "<<integrationProgress.currentTime<<" au) in cube ***" << std::endl;

    // NEO
    if (std::is_same<NEOSS<MatsT,IntsT>,singleSlaterT<MatsT,IntsT>>::value) {

      auto ss_all = dynamic_cast<NEOSS<MatsT,IntsT>*>(&this->singleSlaterSystem);
      SingleSlater<MatsT, IntsT>& ess = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(ss_all->getSubSSBase(std::string("Electronic")))));
      SingleSlater<MatsT, IntsT>& pss = dynamic_cast<SingleSlater<MatsT, IntsT>&>((*(ss_all->getSubSSBase(std::string("Protonic")))));
      auto cube = tdSCFOptions.rtcubes[PAR_TYPE::ELECTRONIC];
      auto pcube = tdSCFOptions.rtcubes[PAR_TYPE::PROTONIC];
      std::string cube_name,pcube_name;
      if(tdSCFOptions.cubeOptsRT.cubeFileName.empty()) {
        pcube_name = "RT_P_" + std::to_string(integrationProgress.currentStep);
        cube_name  = "RT_"   + std::to_string(integrationProgress.currentStep);
      } else {
        cube_name = tdSCFOptions.cubeOptsRT.cubeFileName;
        pcube_name = cube_name + "_RT_P_" + std::to_string(integrationProgress.currentStep);
        cube_name  = cube_name + "_RT_"   + std::to_string(integrationProgress.currentStep);
      }

      // density cube 
      if (tdSCFOptions.cubeOptsRT.denCube) {
        cube->evalDenCube(cube_name,ess.onePDM,true);
        pcube->evalDenCube(pcube_name,pss.onePDM,true);
      }

    } else {

      auto ss = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->singleSlaterSystem);
      auto cube = tdSCFOptions.rtcubes[PAR_TYPE::ELECTRONIC];
      std::string cube_name;
      if(tdSCFOptions.cubeOptsRT.cubeFileName.empty()) {
        cube_name = "RT_" + std::to_string(integrationProgress.currentStep);
      } else {
        cube_name = tdSCFOptions.cubeOptsRT.cubeFileName;
        cube_name = cube_name + "_RT_" + std::to_string(integrationProgress.currentStep);
      }

      // density cube 
      if (tdSCFOptions.cubeOptsRT.denCube) {


        cube->evalDenCube(cube_name,ss->onePDM,true);
      }

    }

  }
  }
}; // RealTime::saveCube

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s) const {
  out << std::setw(38) << "  " + s << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, double v) const {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, size_t v) const {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, std::string v) const {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, double v, std::string u) const {
  out << std::setw(38) << "  " + s << v << u << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::printRunHeader(EMPerturbation& perturbation) const {

  // No printing if silent
  if( this->printLevel == 0 ) return;

  std::cout << BannerTop << std::endl;
  std::cout << "Real-Time Propagation Settings:" << std::endl << std::endl;

  std::cout << std::left << std::setprecision(7);
  std::string AUTime = " \u0127 / Eh";

  RTFormattedLineNew(std::cout,"* Simulation Parameters:");

  //int nSteps = tdSCFOptions.tMax / tdSCFOptions.deltaT;
  RTFormattedLineNew(std::cout, "Simulation Time:", tdSCFOptions.tMax, AUTime);
  RTFormattedLineNew(std::cout, " ", tdSCFOptions.tMax * FSPerAUTime(), " fs");
  RTFormattedLineNew(std::cout, "Number of Steps:", tdSCFOptions.maxSteps);
  RTFormattedLineNew(std::cout, "Step Size:", tdSCFOptions.deltaT, AUTime);
  RTFormattedLineNew(std::cout, " ", tdSCFOptions.deltaT * FSPerAUTime() , " fs");

  std::cout << std::endl;
  RTFormattedLineNew(std::cout,"* Integration Parameters:");

  bool hasMMUT = tdSCFOptions.integrationAlgorithm == RealTimeAlgorithm::RTModifiedMidpoint;
  if constexpr (std::is_same_v<singleSlaterT<MatsT,IntsT>, MultiParticleSS<MatsT,IntsT>>) {
    hasMMUT = false;
    for(const auto& label : this->singleSlaterSystem.getOrder()) {
      RealTimeAlgorithm algorithm = tdSCFOptions.subsystemIntegrationAlgorithms.at(label);
      RTFormattedLineNew(std::cout, label + " Integration:", realTimeAlgorithmName(algorithm));
      hasMMUT = hasMMUT or algorithm == RealTimeAlgorithm::RTModifiedMidpoint;
    }
  } else
    RTFormattedLineNew(std::cout,"Electronic Integration:",
      realTimeAlgorithmName(tdSCFOptions.integrationAlgorithm));

  if(hasMMUT) {
    std::string rstString;
    if(tdSCFOptions.restartAlgorithm == RestartAlgorithm::ForwardEuler )
      rstString = "Forward Euler";
    else if(tdSCFOptions.restartAlgorithm == RestartAlgorithm::ExplicitMagnus2 )
      rstString = "Explicit 2nd Order Magnus";
    RTFormattedLineNew(std::cout, "Restarting MMUT every ", tdSCFOptions.iRestart, " steps with a(n) " + rstString + " step");
  }

  if(tdEMPerturbation.fields.size() > 0 ) {
    std::cout << std::endl;
    RTFormattedLineNew(std::cout,"* Perturbation:\n");

    for(auto &field : tdEMPerturbation.fields) {
      std::cout << std::setw(4) << " ";
      std::cout << "Field " << std::distance(&field,&tdEMPerturbation.fields[0]) + 1<< ":  ";

      auto amp = field->getAmp(0);
      if( dynamic_cast<TDDipoleField&>(*field).emFieldTyp == Electric )
        std::cout << "Electric";
      else
        std::cout << "Magnetic";
      std::cout << " ";

      if( amp.size() == 3 ) std::cout << "Dipole";
      std::cout << " Field\n";
      std::cout << std::setw(4) << " ";
      std::cout << std::setw(20) << " * Amplitude (AU)" << "{ ";
      for(auto i = 0; i < amp.size(); i++) {
        std::cout << amp[i]; if(i != amp.size() - 1) std::cout << ", ";
      }
      std::cout << " }\n";

      std::cout << std::setw(4) << " ";
      try {
        StepField &env = dynamic_cast<StepField&>(*field->envelope);
        std::cout << std::setw(20) << " * Step Field";
        std::cout << std::setw(9) << "TON = "  << std::setw(10) << env.tOn;
        std::cout << "   ";
        std::cout << std::setw(9) << "TOFF = " << std::setw(10) << env.tOff;
        std::cout << std::endl;
      } catch(...) { }
    }
  }

  std::cout << std::endl;
  RTFormattedLineNew(std::cout,"* Misc Parameters:");

  std::string expString;
  if(tdSCFOptions.propagatorAlgorithm == PropagatorAlgorithm::Diagonalization )
    expString = "Eigen Decomposition";
  else if(tdSCFOptions.propagatorAlgorithm == PropagatorAlgorithm::TaylorExpansion )
    expString = "Taylor Expansion";

  RTFormattedLineNew(std::cout,"Matrix Exponential Method:",expString);
  std::cout << std::endl << BannerTop << std::endl;
  std::cout << std::endl << std::fixed << std::right;

  if( this->printLevel == 1 ) {
    std::cout << std::setprecision(4);
    std::cout << std::setw(11) << "Time (a.u.)" << " ";

    std::cout << std::setprecision(10);
    std::cout << std::setw(16) << "Energy (Eh)" << " ";

    std::cout << std::setprecision(8);
    std::cout << std::setw(16) << "Dipole (X)" << " ";
    std::cout << std::setw(16) << "Dipole (Y)" << " ";
    std::cout << std::setw(16) << "Dipole (Z)" << " ";

    std::cout << std::endl << bannerTop << std::endl << std::endl;
  }

  if( this->printLevel == -1 ) this->printLevel = 0;

};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::printIteration(bool printDiff) {
  if(this->printLevel == 1) printStepSummary();
  else if(this->printLevel > 1) printStepDetail();
  if( printDen ) this->singleSlaterSystem.onePDM->output(std::cout, "OnePDM at t=" + std::to_string(integrationProgress.currentTime), true);
  if (tdSCFOptions.Rtprintden != 0) {
    int Rtprintdenstep = 0;
    Rtprintdenstep = integrationProgress.currentStep % tdSCFOptions.Rtprintden;
    if (Rtprintdenstep ==0) {
      this->singleSlaterSystem.onePDM->output(std::cout, "OnePDM at t=" + std::to_string(integrationProgress.currentTime), true);
    } 
  }
  if( tdSCFOptions.orbitalPopFreq != 0 && integrationProgress.currentStep % tdSCFOptions.orbitalPopFreq == 0) 
      this->singleSlaterSystem.printOrbitalPopulation(std::cout);
};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::printStepSummary() {
  std::cout << std::fixed << std::right;

  std::cout << std::setprecision(4);
  std::cout << std::setw(11) << integrationProgress.currentTime << " ";

  std::cout << std::setprecision(10);
  std::cout << std::setw(16) << this->singleSlaterSystem.getTotalEnergy() << " ";

  std::cout << std::setprecision(8);
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[0] << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[1] << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[2] << " ";

  std::cout << std::endl;

};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::printStepDetail() {
  std::cout << bannerTop << "\n\n";
  std::cout << std::fixed << std::right;
  std::cout << "Step: " << std::setw(7) << integrationProgress.currentStep << '\n';

  std::cout << std::setprecision(5) << "Time: ";
  std::cout << std::setw(11) << integrationProgress.currentTime << " (au) | ";
  std::cout << std::setw(11) << integrationProgress.currentTime * FSPerAUTime() << " (fs)\n";

  std::cout << std::setprecision(12) << "Energy: ";
  std::cout << std::setw(24) << this->singleSlaterSystem.totalEnergy << " (Hartree)\n";

  std::cout << std::setprecision(8) << "Dipole: ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[0] / EBohrPerDebye() << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[1] / EBohrPerDebye() << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[2] / EBohrPerDebye() << " (Debye)";
  std::cout << std::endl;
};

}; // namespace ChronusQ

#if 0
// set the density matrix to use the current orthonormalized density matrix
  if( std::is_same<NEOSS<MatsT,IntsT>,singleSlaterT<MatsT,IntsT>>::value ) {
    auto neoSingleSlaterSystem = dynamic_cast<NEOSS<MatsT,IntsT>*>(&this->singleSlaterSystem);
    auto neoMap = neoSingleSlaterSystem->getSubsystemMap();
    auto neoSubsystemOrder = neoSingleSlaterSystem->getOrder();
    assert( !neoMap.empty() );
    // Loop over all subsystems
    size_t i = 0;
    for( auto& neoSubsystemLabel: neoSubsystemOrder ) {
      auto &neoSubsystem = neoMap[neoSubsystemLabel];
      if(neoSubsystem.get()->nC == 1) {
        if(neoSubsystem.get()->iCS) {
          *neoSubsystem.get()->onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i]);
        } else {
          *neoSubsystem.get()->onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i],this->onePDMSquareOrtho[i+1]);
          i++;
        }
      } else {
        *neoSubsystem.get()->onePDMOrtho = this->onePDMSquareOrtho[i].template spinScatter<MatsT>();
      }
      // XSLI TODO: move ortho2aoDen to OribtalModifier
      neoSubsystem.get()->ortho2aoDen();
      i++;
    }
  } else for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
    if(this->singleSlaterSystem.nC == 1) {
      if(this->singleSlaterSystem.iCS) {
        *this->singleSlaterSystem.onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i]);
      } else {
        *this->singleSlaterSystem.onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i],this->onePDMSquareOrtho[i+1]);
        i++;
      }
    } else {
      *this->singleSlaterSystem.onePDMOrtho = this->onePDMSquareOrtho[i].template spinScatter<MatsT>();
    }
    // XSLI TODO: move ortho2aoDen to OribtalModifier
    this->singleSlaterSystem.ortho2aoDen();
  }
#endif
