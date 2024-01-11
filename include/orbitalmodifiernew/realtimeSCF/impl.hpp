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
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::run(EMPerturbation &perturbation) {

  staticEMPerturbation = perturbation;

  ProgramTimer::tick("Real Time Total");

  printRunHeader(staticEMPerturbation);

  bool startStep(false); // startStep the MMUT iterations
  bool finalStep(false); // Wrap up the MMUT iterations
  bool normalStep(true); // Normal step for the MMUT iterations

  // Initialize the RT iterations (on root process).
  this->singleSlaterSystem.formDensity();
  this->ao2orthoDen(); // transform and save the orthonormal density in onePDMSquareOrtho (on root process)
  this->ao2orthoFock(); // transform and save the orthonormal Fock matrix in fockSquareOrtho (on root process)
  // Disable RI-K-Coefficient Contraction during RT as coefficients are not updated.
  // Use density contration during RT instead
  this->singleSlaterSystem.setDenEqCoeff(false);

  // Restore 1PDM Ortho and integration process (on root process)
  if ( tdSCFOptions.restoreFromStep != 0 ) this->restoreState();
  if ( MPIRank(this->mpiComm) == 0)
    for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) 
      previousOnePDMSquareOrtho[i] = this->onePDMSquareOrtho[i];

  //integrationProgress.currentTime = tdSCFOptions.restoreFromStep * tdSCFOptions.deltaT;

  for(integrationProgress.currentStep = tdSCFOptions.restoreFromStep;
      integrationProgress.currentStep <= tdSCFOptions.maxSteps;
      integrationProgress.currentTime += tdSCFOptions.deltaT, integrationProgress.currentStep++) {

    ProgramTimer::tick("Real Time Iter");

    // Perturbation for the current time
    EMPerturbation currentPerturbation = tdEMPerturbation.getPert(integrationProgress.currentTime);

    // "Start" the MMUT if the current step is the first step or a restart step or field discontinuous
    // "Finish" the MMUT if the current step is the last step or a restart step or field discontinuous
    if(tdSCFOptions.integrationAlgorithm == RTModifiedMidpoint ) {
      startStep = finalStep or (integrationProgress.currentStep == tdSCFOptions.restoreFromStep );
      finalStep = (integrationProgress.currentStep == tdSCFOptions.maxSteps ) or tdEMPerturbation.isFieldDiscontinuous(integrationProgress.currentTime, tdSCFOptions.deltaT);
      if(tdSCFOptions.iRestart > 0 and (integrationProgress.currentStep + 1) % tdSCFOptions.iRestart == 0) finalStep = true;
    };

    // Determine the step type, half or full step, for the current integration step for MMUT
    if(tdSCFOptions.integrationAlgorithm == RTModifiedMidpoint ) {
      if(startStep or finalStep) integrationProgress.currentDeltaT = tdSCFOptions.deltaT;
      else integrationProgress.currentDeltaT = 2. * tdSCFOptions.deltaT;
    } else {
      integrationProgress.currentDeltaT = tdSCFOptions.deltaT;
    }

    normalStep = true;
    if(finalStep or startStep) normalStep = false;
    std::vector<SquareMatrix<MatsT>> onePDMSquareOrthoSave;
    if(MPIRank(this->mpiComm) == 0){
      if(normalStep) {
        std::swap(this->onePDMSquareOrtho,this->previousOnePDMSquareOrtho);
      } else if(startStep or finalStep) {
        if(printLevel > 0 and startStep) std::cout << "  *** Starting MMUT ***\n";
        if(printLevel > 0 and finalStep) std::cout << "  *** Finishing MMUT ***\n";
        for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
          this->onePDMSquareOrtho[i] = this->previousOnePDMSquareOrtho[i];
          if(tdSCFOptions.restartAlgorithm == ExplicitMagnus2) onePDMSquareOrthoSave.emplace_back(this->previousOnePDMSquareOrtho[i]);
        }
      }
    }

    this->formFock(false, integrationProgress.currentTime);
    this->singleSlaterSystem.computeEnergy(currentPerturbation);
    this->singleSlaterSystem.computeProperties(currentPerturbation);
    this->saveState(currentPerturbation); // Save the current state for every iSave steps
    std::vector<std::shared_ptr<SquareMatrix<MatsT>>> fock_k = this->singleSlaterSystem.getFock();
    formPropagator(fock_k);
    doPropagation();
    if( printLevel > 0 and (MPIRank(this->mpiComm) == 0) ) printIteration();

    MPI_Barrier(MPI_COMM_WORLD);

    // Explicit Magnus 2
    if ((finalStep or startStep) and tdSCFOptions.restartAlgorithm == ExplicitMagnus2 ) {
      if(MPIRank(this->mpiComm) == 0){
      	for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
      	  this->onePDMSquareOrtho[i] = this->previousOnePDMSquareOrtho[i];
      	  this->previousOnePDMSquareOrtho[i] = onePDMSquareOrthoSave[i];
      	}
      }
      this->formFock(false, integrationProgress.currentTime + tdSCFOptions.deltaT); // F(k+1)
      std::vector<std::shared_ptr<SquareMatrix<MatsT>>> fock_k1 = this->singleSlaterSystem.getFock();
      for( size_t i = 0; i < fock_k.size(); i++ )
        *fock_k[i] = 0.5 * (*fock_k[i] + *fock_k1[i]); // compute 0.5 * (F(k) + F(k+1))
      formPropagator(fock_k);
      doPropagation();
      if(MPIRank(this->mpiComm) == 0)
        for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) 
          this->onePDMSquareOrtho[i] = onePDMSquareOrthoSave[i];

      MPI_Barrier(MPI_COMM_WORLD);
    }  // End 2nd order magnus

    ProgramTimer::tock("Real Time Iter");

  } // Time loop

  ProgramTimer::tock("Real Time Total");

}; // RealTime::doPropagation


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT, IntsT>::formFock(bool increment, double time) {

  //this->singleSlaterSystem.setOnePDMOrtho(this->onePDMSquareOrtho.data());
  //this->singleSlaterSystem.ortho2aoDen();

  // Transfrom density to AO (on root process)
  this->ortho2aoDen(this->onePDMSquareOrtho);
  // Scatter to spin blocks (on root process) and broacast if doing MPI
  this->singleSlaterSystem.setOnePDMAO(this->onePDMSquareAO.data());

  ProgramTimer::timeOp("Form Fock", [&]() {
      // Get perturbation for the current time and build a Fock matrix
      EMPerturbation pert_t = tdEMPerturbation.getPert(time);
      // Add the SCF Perturbation
      if ( tdSCFOptions.includeSCFField ) for( auto& field : staticEMPerturbation.fields ) pert_t.addField(field );
      this->singleSlaterSystem.formFock(pert_t,increment);
  });

}; // RealTime::formFock

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
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::formPropagator(std::vector<std::shared_ptr<SquareMatrix<MatsT>>> fockSquareAO) {

  ROOT_ONLY(this->mpiComm);

  ProgramTimer::tick("Propagator Formation");
  this->ao2orthoFock(fockSquareAO);
  for( size_t i = 0; i < this->fockSquareOrtho.size(); i++ ) {
    size_t NB = this->fockSquareOrtho[i].dimension();
    MatExp('D',NB,dcomplex(0.,-integrationProgress.currentDeltaT),
           this->fockSquareOrtho[i].pointer(),NB,unitarySquareOrtho[i].pointer(),NB, this->memManager);
  }
  ProgramTimer::tock("Propagator Formation");
#if 0
  prettyPrintSmart(std::cout,"U",unitarySquareOrtho[i].pointer(),NB,NB,NB);
#endif
}; // RealTime::formPropagator

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::doPropagation() {

  ROOT_ONLY(this->mpiComm);

  ProgramTimer::tick("Propagate Density");

  for( size_t i = 0; i < unitarySquareOrtho.size(); i++ ) {
    size_t NB = this->fockSquareOrtho[i].dimension();
    SquareMatrix<MatsT> SCR(this->memManager,NB);

    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, dcomplex(1.),
               unitarySquareOrtho[i].pointer(), NB,
               this->previousOnePDMSquareOrtho[i].pointer(), NB, dcomplex(0.), SCR.pointer(), NB);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, dcomplex(1.), SCR.pointer(), NB,
               unitarySquareOrtho[i].pointer(), NB, dcomplex(0.), this->previousOnePDMSquareOrtho[i].pointer(), NB);

  }

  ProgramTimer::tock("Propagate Density");

}; // RealTime::doPropagation


template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::createRTDataSets(size_t maxPoints) {

  ROOT_ONLY(this->mpiComm);

  if( maxPoints == 0 ) integrationProgress.maxSavePoints = tdSCFOptions.maxSteps/tdSCFOptions.iSave + 2;
  else integrationProgress.maxSavePoints = maxPoints;

  if(tdSCFOptions.restoreFromStep != 0 ) return;

  savFile.createGroup("RTNEW");

  savFile.createDataSet<size_t>("RTNEW/ISAVE", {1});
  savFile.createDataSet<size_t>("RTNEW/LASTSAVEPOINT", {1});
  savFile.createDataSet<size_t>("RTNEW/MAXSAVEPOINTS", {1});

  savFile.createDataSet<size_t>("RTNEW/STEP", {integrationProgress.maxSavePoints});
  savFile.createDataSet<double>("RTNEW/TIME", {integrationProgress.maxSavePoints});
  savFile.createDataSet<double>("RTNEW/ENERGY", {integrationProgress.maxSavePoints});
  savFile.createDataSet<double>("RTNEW/LEN_ELEC_DIPOLE", {integrationProgress.maxSavePoints*3});
  savFile.createDataSet<double>("RTNEW/LEN_ELEC_DIPOLE_FIELD", {integrationProgress.maxSavePoints*3});

  for( size_t i = 0; i < this->onePDMSquareOrtho.size(); i++ ) {
    size_t nBasis = this->onePDMSquareOrtho[i].dimension();
    savFile.createDataSet<dcomplex>("RTNEW/TD_1PDM_ORTHO"+std::to_string(i), {integrationProgress.maxSavePoints*nBasis*nBasis});
    savFile.createDataSet<double>("RTNEW/ORBITALPOPULATION"+std::to_string(i), {integrationProgress.maxSavePoints*nBasis});
  }

}; // RealTime::createRTDataSets

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::saveState(EMPerturbation& currentPerturbation) {

  ROOT_ONLY(this->mpiComm);

  savFile.safeWriteData("RTNEW/ISAVE", &tdSCFOptions.iSave, {1});
  savFile.safeWriteData("RTNEW/MAXSAVEPOINTS", &integrationProgress.maxSavePoints, {1});

  integrationProgress.time.push_back(integrationProgress.currentTime);
  integrationProgress.energy.push_back(this->singleSlaterSystem.totalEnergy);
  integrationProgress.electricDipole.push_back(this->singleSlaterSystem.elecDipole);
  if( currentPerturbation.fields.size() > 0 ) integrationProgress.electricDipoleField.push_back(currentPerturbation.getDipoleAmp(Electric) );

  if (integrationProgress.currentStep == tdSCFOptions.restoreFromStep
      or integrationProgress.currentStep == tdSCFOptions.maxSteps
      or integrationProgress.currentStep % tdSCFOptions.iSave == 0) {

    std::cout << "  *** Saving step #"<<integrationProgress.currentStep<<"( t = "<<integrationProgress.currentTime<<" au) to binary file ***" << std::endl;
    savFile.safeWriteData("RTNEW/LASTSAVEPOINT", &(integrationProgress.lastSavePoint), {1});
    savFile.partialWriteData("RTNEW/TIME", &integrationProgress.currentTime, {integrationProgress.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("RTNEW/STEP", &integrationProgress.currentStep, {integrationProgress.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("RTNEW/ENERGY", &this->singleSlaterSystem.totalEnergy, {integrationProgress.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("RTNEW/LEN_ELEC_DIPOLE", &this->singleSlaterSystem.elecDipole[0],{integrationProgress.lastSavePoint*3}, {3},{0},{3});
    std::array<double,3> elecDipoleField = currentPerturbation.getDipoleAmp(Electric);
    if (integrationProgress.electricDipoleField.size() > 0)
      savFile.partialWriteData("RTNEW/LEN_ELEC_DIPOLE_FIELD",&elecDipoleField[0], {integrationProgress.lastSavePoint*3}, {3},{0},{3});

    for (size_t i = 0; i < this->onePDMSquareOrtho.size(); i++) {
      size_t nBasis = this->onePDMSquareOrtho[i].dimension();
      savFile.partialWriteData("RTNEW/TD_1PDM_ORTHO" + std::to_string(i), this->onePDMSquareOrtho[i].pointer(),{integrationProgress.lastSavePoint*nBasis*nBasis}, {nBasis * nBasis}, {0}, {nBasis * nBasis});
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

    hsize_t maxSavePoints, lastSavePoint;
    savFile.readData("RTNEW/MAXSAVEPOINTS", &maxSavePoints);
    if ( maxSavePoints != integrationProgress.maxSavePoints ) CErr("Mismatched requested and saved propagation length!");
    savFile.readData("RTNEW/LASTSAVEPOINT", &lastSavePoint);

    size_t iSave;
    savFile.readData("RTNEW/ISAVE", &iSave);

    if(tdSCFOptions.restoreFromStep < 0) integrationProgress.lastSavePoint = lastSavePoint;
    else {
      integrationProgress.lastSavePoint = tdSCFOptions.restoreFromStep/iSave;
      if(integrationProgress.lastSavePoint> lastSavePoint) integrationProgress.lastSavePoint = lastSavePoint;
    }

    // Restore time dependent density
    try {
      savFile.partialReadData("RTNEW/STEP", &integrationProgress.currentStep, {integrationProgress.lastSavePoint}, {1}, {0}, {1});
      savFile.partialReadData("RTNEW/TIME", &integrationProgress.currentTime, {integrationProgress.lastSavePoint}, {1}, {0}, {1});

      for (size_t i = 0; i < this->onePDMSquareOrtho.size(); i++) {
        size_t nBasis = this->onePDMSquareOrtho[i].dimension();
        if(savFile.getDims("RTNEW/TD_1PDM_ORTHO" + std::to_string(i)) != std::vector<hsize_t>{integrationProgress.maxSavePoints*nBasis*nBasis})
          CErr("Mismatched requested and saved propagation length!");
        savFile.partialReadData("RTNEW/TD_1PDM_ORTHO" + std::to_string(i), this->onePDMSquareOrtho[i].pointer(), {integrationProgress.lastSavePoint*nBasis*nBasis}, {nBasis * nBasis}, {0}, {nBasis * nBasis});
      }
    } catch(...) { }

    if( printLevel > 0 ) {
      std::cout << "  *** Restoring from step " << integrationProgress.currentStep << " (";
      std::cout << std::setprecision(4) << integrationProgress.currentTime << " au) ***" << std::endl;
    }

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
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::orbitalPop() {

  // Spin-Gather Ortho Density
  std::vector<SquareMatrix<MatsT>> orthoDen;
  if(this->singleSlaterSystem.nC == 1 ){
    orthoDen = this->singleSlaterSystem.onePDMOrtho->template spinGatherToBlocks<MatsT>(false);
  } else {
    orthoDen.push_back(this->singleSlaterSystem.onePDMOrtho->template spinGather<MatsT>());
  }

  // Transform a copy of the MOs because
  std::vector<SquareMatrix<MatsT>> orthoMO = this->singleSlaterSystem.mo;

  this->singleSlaterSystem.orthoAB->nonortho2orthoCoeffs(orthoMO);

  // Transform alpha Density and compute populations
  size_t NB = orthoMO[0].dimension();
  std::vector<double> population;
  std::vector<SquareMatrix<MatsT>> moDen;
  moDen.push_back( orthoDen[0].transform('N',orthoMO[0].pointer(),NB,NB) );
  for( size_t i=0; i<NB; ++i)
    population.push_back( std::real(moDen[0](i,i)) );

  // UHF Beta populations
  if(this->singleSlaterSystem.nC == 1 and not this->singleSlaterSystem.iCS ){
    moDen.push_back( orthoDen[1].transform('N',orthoMO[1].pointer(),NB,NB) );
    for( size_t i=0; i<NB; ++i)
      population.push_back( std::real(moDen[1](i,i)) );
  }

  // Printing
  if( this->printLevel > 1 ) {

    size_t orbPerRow = 5;
    auto printBlock = [&](std::string header, size_t& start, size_t n){
      std::cout << header << std::endl;
      std::cout << std::fixed << std::setprecision(11);

      for(auto idx = 0; idx < n; idx += orbPerRow) {

        size_t end = idx + orbPerRow < n ? orbPerRow : n - idx;
        for(auto idummy = idx; idummy < idx+end; idummy++) {
          std::cout << std::setw(15) << population[start+idummy];
        }
        std::cout << '\n';
      }
      start += n;
    };

    if( this->printLevel > 3 ){
      moDen[0].output(std::cout, "MO Density Matrix", true);
      if(this->singleSlaterSystem.nC == 1 and not this->singleSlaterSystem.iCS )
        moDen[1].output(std::cout, "MO Beta Density Matrix", true);
    }

    size_t start = 0;
    if(this->singleSlaterSystem.nC == 1 ) {
      printBlock("Alpha occupied orbitals", start, this->singleSlaterSystem.nOA);
      printBlock("Alpha virtual orbitals", start, this->singleSlaterSystem.nVA);
      if( not this->singleSlaterSystem.iCS ){
        printBlock("Beta occupied orbitals", start, this->singleSlaterSystem.nOB);
        printBlock("Beta virtual orbitals", start, this->singleSlaterSystem.nVB);
      }
    }
    else if(this->singleSlaterSystem.nC == 2 ){
      printBlock("Occupied orbitals", start, this->singleSlaterSystem.nO);
      printBlock("Virtual orbitals", start, this->singleSlaterSystem.nV);
    } else if(this->singleSlaterSystem.nC == 4 ){
      start += NB/2;
      printBlock("Positive Energy Occupied orbitals", start, this->singleSlaterSystem.nO);
      printBlock("Positive Energy Virtual orbitals", start, this->singleSlaterSystem.nV);
    }
    std::cout << std::flush;
  }

}; // RealTime :: orbitalPop

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s) {
  out << std::setw(38) << "  " + s << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, double v) {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, size_t v) {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, std::string v) {
  out << std::setw(38) << "  " + s << v << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::RTFormattedLineNew(std::ostream &out, std::string s, double v, std::string u) {
  out << std::setw(38) << "  " + s << v << u << std::endl;
}

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::printRunHeader(EMPerturbation& perturbation) {

  // No printing if silent
  if( this->printLevel == 0 ) return;

  std::cout << BannerTop << std::endl;
  std::cout << "Real-Time Propagation Settings:" << std::endl << std::endl;

  std::cout << std::left << std::setprecision(7);
  std::string AUTime = " \u0127 / Eh";

  RTFormattedLineNew(std::cout,"* Simulation Parameters:");

  //int nSteps = tdSCFOptions.tMax / tdSCFOptions.deltaT;
  RTFormattedLineNew(std::cout, "Simulation Time:", tdSCFOptions.tMax, AUTime);
  RTFormattedLineNew(std::cout, " ", tdSCFOptions.tMax * FSPerAUTime, " fs");
  RTFormattedLineNew(std::cout, "Number of Steps:", tdSCFOptions.maxSteps);
  RTFormattedLineNew(std::cout, "Step Size:", tdSCFOptions.deltaT, AUTime);
  RTFormattedLineNew(std::cout, " ", tdSCFOptions.deltaT * FSPerAUTime , " fs");

  std::cout << std::endl;
  RTFormattedLineNew(std::cout,"* Integration Parameters:");

  std::string methString;
  if(tdSCFOptions.integrationAlgorithm == RTModifiedMidpoint )
    methString = "Modified Midpoint Unitary Transformation (MMUT)";
  else if(tdSCFOptions.integrationAlgorithm == RTExplicitMagnus2)
    methString = "Explicit 2nd Order Magnus";

  RTFormattedLineNew(std::cout,"Electronic Integration:",methString);

  if(tdSCFOptions.integrationAlgorithm == RTModifiedMidpoint ) {
    std::string rstString;
    if(tdSCFOptions.restartAlgorithm == RTForwardEuler )
      rstString = "Forward Euler";
    else if(tdSCFOptions.restartAlgorithm == RTExplicitMagnus2 )
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
  if(tdSCFOptions.propagatorAlgorithm == Diagonalization )
    expString = "Eigen Decomposition";
  else if(tdSCFOptions.propagatorAlgorithm == TaylorExpansion )
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
        } else{
        }
      }
  
  //if(tdSCFOptions.iPrint != 0 && integrationProgress.currentStep % tdSCFOptions.iPrint == 0) orbitalPop();
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
  std::cout << std::setw(11) << integrationProgress.currentTime * FSPerAUTime << " (fs)\n";

  std::cout << std::setprecision(12) << "Energy: ";
  std::cout << std::setw(24) << this->singleSlaterSystem.totalEnergy << " (Hartree)\n";

  std::cout << std::setprecision(8) << "Dipole: ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[0] * EBohrPerDebye << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[1] * EBohrPerDebye << " ";
  std::cout << std::setw(16) << this->singleSlaterSystem.elecDipole[2] * EBohrPerDebye << " (Debye)";
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
          *neoSubsystem.get()->onePDMOrtho = PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i]);
        } else {
          *neoSubsystem.get()->onePDMOrtho = PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i],this->onePDMSquareOrtho[i+1]);
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
        *this->singleSlaterSystem.onePDMOrtho = PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i]);
      } else {
        *this->singleSlaterSystem.onePDMOrtho = PauliSpinorSquareMatrices<MatsT>::spinBlockScatterBuild(this->onePDMSquareOrtho[i],this->onePDMSquareOrtho[i+1]);
        i++;
      }
    } else {
      *this->singleSlaterSystem.onePDMOrtho = this->onePDMSquareOrtho[i].template spinScatter<MatsT>();
    }
    // XSLI TODO: move ortho2aoDen to OribtalModifier
    this->singleSlaterSystem.ortho2aoDen();
  }
#endif
