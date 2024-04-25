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

#include <realtime.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <matrix.hpp>

#include <util/matout.hpp>
#include <util/timer.hpp>
#include <unsupported/Eigen/MatrixFunctions>

template <size_t N, typename T>
std::array<T,N> valarray2array(const std::valarray<T> &x) {
  assert( x.size() == N );

  std::array<T,N> arr;
  std::copy_n(&x[0],N,arr.begin());
  return arr;
 
};

namespace ChronusQ {

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::doPropagation() {

    ProgramTimer::tick("Real Time Total");

    printRTHeader();

    if ( savFile.exists() )
      if ( restart )
        restoreState();

    // Upon entry to RT, assume only the orthonormal density is valid
    for(auto idx = 0; idx < systems_.size(); idx++) {

      // Allow the printing of contraction timing during propagation
      systems_[idx]->TPI->printContractionTiming = this->printContractionTiming? true : false;
      if(auto neofock = std::dynamic_pointer_cast<NEOFockBuilder<dcomplex,double>>(systems_[idx]->fockBuilder)){
        neofock->setPrintContractionTiming(this->printContractionTiming);
      } else if(auto neoks = std::dynamic_pointer_cast<NEOKohnShamBuilder<dcomplex,double>>(systems_[idx]->fockBuilder)){
        if(auto neofock = dynamic_cast<NEOFockBuilder<dcomplex,IntsT>*>(neoks->getUpstream())){
          neofock->setPrintContractionTiming(this->printContractionTiming);
        }else{
          CErr("The upstream of NEOKohnShamBuilder is set up incorrectly! Expecting a NEOFockBuilder ptr");
        }
      }

      // For propagation, all RI-K contractions are done with density instead of coefficients
      systems_[idx]->setDenEqCoeff(false);
      systems_[idx]->computeOrtho();
      systems_[idx]->ortho2aoDen();
      systems_[idx]->ortho2aoMOs();
    }

    bool Start(false); // Start the MMUT iterations
    bool FinMM(false); // Wrap up the MMUT iterations

    size_t maxStep = (size_t)((intScheme.tMax + intScheme.deltaT/4)/intScheme.deltaT);

    curState.xTime = intScheme.restoreStep * intScheme.deltaT;

    for( curState.iStep = intScheme.restoreStep; 
         curState.iStep <= maxStep;
         curState.xTime += intScheme.deltaT, curState.iStep++) {

      ProgramTimer::tick("Real Time Iter");

      // Perturbation for the current time
      EMPerturbation pert_t = pert.getPert(curState.xTime);



    // Determine the step type for the current integration step 
    if( intScheme.intAlg == MMUT ) {

      // "Start" the MMUT if this is the first step or we just
      // "Finished" the MMUT segment
      Start = ( curState.iStep == intScheme.restoreStep ) or FinMM;

      // "Start" the MMUT if the current step index is a restart
      // step
      if( intScheme.iRstrt > 0 ) 
        Start = Start or ( curState.iStep % intScheme.iRstrt == 0 );

      // "Finish" the MMUT if this is the last step
      // NOTE: To compare to Gaussian, do NOT do this restart for Ehrenfest
      FinMM = ( curState.iStep == maxStep );

      // TODO: "Finish" the MMUT if the field turns on or off
      FinMM = FinMM or pert.isFieldDiscontinuous(curState.xTime, intScheme.deltaT);
        
      // "Finish" the MMUT if the next step will be a restart step
      if( intScheme.iRstrt > 0 ) 
        FinMM = FinMM or ( (curState.iStep + 1) % intScheme.iRstrt == 0 );

      // If "Starting" or "Finishing" the MMUT, the step type is
      // the specified restart step type, else it is the MMUT step
      if( Start or FinMM ) curState.curStep = intScheme.rstStep;
      else                 curState.curStep = ModifiedMidpoint;

      if( (Start or FinMM) && printLevel > 0 )
        std::cout <<"  *** Restarting MMUT ***\n";

    // For non leapfrog scheme, the step type is constant
    } else if ( intScheme.intAlg == ExpMagnus2 )
      curState.curStep = ExplicitMagnus2;





      for(auto idx = 0; idx < systems_.size(); idx++) {
        // Handle density copies / swaps for the current step
        //  + Determine the step size
          
        if( curState.curStep == ModifiedMidpoint ) {
          // Swap the saved density with the SingleSlater density
            
          // DOSav(k) = DO(k)
          // DO(k)    = DO(k-1)
          std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> tmp = DOSav[idx];
          DOSav[idx] = systems_[idx]->onePDMOrtho;
          systems_[idx]->onePDMOrtho = tmp;


          curState.stepSize = 2. * intScheme.deltaT;

        } else {
          // Save a copy of the SingleSlater density in the saved density
          // storage

          // DOSav(k) = DO(k)
          *DOSav[idx] = *systems_[idx]->onePDMOrtho;
     
          curState.stepSize = intScheme.deltaT;

        }
      }

      // Form the Fock matrix at the current time
      for(auto idx = 0; idx < systems_.size(); idx++) {
        this->formFock(false,curState.xTime,idx);
      }

      // Compute properties for D(k)
      propagator_.computeEnergy(pert_t);
      propagator_.computeProperties(pert_t);

      // Save data
      // TODO: Fix this when we have a stable definition of MD + electronic steps
      saveState(pert_t);

      // Save D(k) if doing Magnus 2
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> den_k;
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> denOrtho_k;
      for(auto idx = 0; idx < systems_.size(); idx++) {
        if ( curState.curStep == ExplicitMagnus2 ) {
          den_k.push_back(std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(*systems_[idx]->onePDM));
          denOrtho_k.push_back(std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(*systems_[idx]->onePDMOrtho));
        }
      }

      // Printing out real-time density
      if( printDen ) propagator_.onePDM->output(std::cout, "OnePDM at t=" + std::to_string(curState.xTime), true);

      // Print progress line in the output file
      printRTStep();
      if( this->orbitalPopFreq != 0 &&
          curState.iStep % this->orbitalPopFreq == 0) {
        orbitalPop();
      }




      for(auto idx = 0; idx < systems_.size(); idx++) {
        // Orthonormalize the AO Fock matrix
        // F(k) -> FO(k)
        systems_[idx]->ao2orthoFock();


        // Form the propagator from the orthonormal Fock matrix
        // FO(k) -> U**H(k) = exp(- i * dt * FO(k) )
        formPropagator(idx);

        // Propagator the orthonormal density matrix
        // DO (in propagator_) will now store DO(k+1)
        //
        // DO(k+1) = U**H(k) * DO * U(k)
        // - Where DO is what is currently stored in propagator_
        //
        // ***
        // This function also transforms DO(k+1) to the AO
        // basis ( DO(k+1) -> D(k+1) in propagator_ ) and
        // computes the change in density from the previous 
        // AO density ( delD = D(k+1) - D(k) ) 
        // ***
        propagateWFN(idx);
      }

      //
      // Second order magnus
      //
      if ( curState.curStep == ExplicitMagnus2 ) {

        for(auto idx = 0; idx < systems_.size(); idx++) {
          // F(k)
          cqmatrix::PauliSpinorMatrices<dcomplex> fock_k(*systems_[idx]->fockMatrix);
          
          // F(k + 1)
          formFock(false, curState.xTime + intScheme.deltaT, idx);

          // Store 0.5 * ( F(k) + F(k+1) ) in propagator_.fockMatrix
          *systems_[idx]->fockMatrix = 0.5 * (fock_k + *systems_[idx]->fockMatrix);
        }


        for(auto idx = 0; idx < systems_.size(); idx++) {
          // Restore old densities
          *systems_[idx]->onePDM = *den_k[idx];
          *systems_[idx]->onePDMOrtho = *denOrtho_k[idx];
        }

        // Repeat formation of propagator and propagation
        for(auto idx = 0; idx < systems_.size(); idx++) {
          systems_[idx]->ao2orthoFock();
          formPropagator(idx);
          propagateWFN(idx);
        }

      }  // End 2nd order magnus

      ProgramTimer::tock("Real Time Iter");

    } // Time loop


    ProgramTimer::tock("Real Time Total");

  //mathematicaPrint(std::cerr,"Dipole-X",&data.ElecDipole[0][0],
  //  curState.iStep,1,curState.iStep,3);

  }; // RealTime::doPropagation


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
  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::formPropagator(size_t idx) {

    ProgramTimer::tick("Propagator Formation");

    size_t NB = systems_[idx]->nAlphaOrbital();

    // Form U

    // Restricted
    if( not UH[idx]->hasZ() ) {
      // See docs for factor of 2
      MatExp('D',NB,dcomplex(0.,-curState.stepSize/2.),
        systems_[idx]->fockMatrixOrtho->S().pointer(),NB,UH[idx]->S().pointer(),NB);

      blas::scal(NB*NB,dcomplex(2.),UH[idx]->S().pointer(),1);

    // Unrestricted
    } else if( not UH[idx]->hasXY() ) {

      // Transform SCALAR / MZ -> ALPHA / BETA
      std::vector<cqmatrix::Matrix<dcomplex>> Fblocks =
          systems_[idx]->fockMatrixOrtho->template spinGatherToBlocks<dcomplex>(false);
      std::vector<cqmatrix::Matrix<dcomplex>> UHblocks;
      UHblocks.reserve(2);
      UHblocks.emplace_back(NB);
      UHblocks.emplace_back(NB);

      MatExp('D',NB,dcomplex(0.,-curState.stepSize),
        Fblocks[0].pointer(),NB,UHblocks[0].pointer(),NB);
      MatExp('D',NB,dcomplex(0.,-curState.stepSize),
        Fblocks[1].pointer(),NB,UHblocks[1].pointer(),NB);

      // Transform ALPHA / BETA -> SCALAR / MZ
      *UH[idx] = cqmatrix::PauliSpinorMatrices<dcomplex>::
          spinBlockScatterBuild<dcomplex>(UHblocks[0],UHblocks[1]);

    // Generalized (2C) or 4C
    } else {

      size_t nC = systems_[idx]->nC;
      cqmatrix::Matrix<dcomplex> FnC(systems_[idx]->fockMatrixOrtho->template spinGather<dcomplex>());
      cqmatrix::Matrix<dcomplex> UHnC(nC*NB);

      MatExp('D',nC*NB,dcomplex(0.,-curState.stepSize),
             FnC.pointer(),nC*NB,UHnC.pointer(),nC*NB);

      *UH[idx] = UHnC.template spinScatter<dcomplex>();

    }

    ProgramTimer::tock("Propagator Formation");
#if 0

    prettyPrintSmart(std::cout,"UH Scalar",UH[idx]->S().pointer(),NB,NB,NB);
    if( UH[idx]->hasZ() )
      prettyPrintSmart(std::cout,"UH MZ",UH[idx]->Z().pointer(),NB,NB,NB);

#endif
    
  }; // RealTime::formPropagator



  
  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::propagateWFN(size_t idx) {

    ProgramTimer::tick("Propagate WFN");

    size_t NB = systems_[idx]->nAlphaOrbital();
    size_t NC = systems_[idx]->nC;
    dcomplex *SCR  = CQMemManager::get().malloc<dcomplex>(NC*NC*NB*NB);
    dcomplex *SCR1 = CQMemManager::get().malloc<dcomplex>(NC*NC*NB*NB);

    if( not UH[idx]->hasXY() ) {

      // Create X(S) = (U**H * DO)(S) in SCR 

      // SCR = 0.5 * U(S)**H * DO(S)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,dcomplex(0.5),UH[idx]->S().pointer(),NB,
          systems_[idx]->onePDMOrtho->S().pointer(),NB,dcomplex(0.),SCR,NB);

      // SCR += 0.5 * U(Z)**H * DO(Z)
      if( UH[idx]->hasZ() )
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,dcomplex(0.5),UH[idx]->Z().pointer(),NB,
          systems_[idx]->onePDMOrtho->Z().pointer(),NB,dcomplex(1.),SCR,NB);





      // Create X(Z) = (U**H * DO)(Z) in SCR1
        
      if( systems_[idx]->onePDMOrtho->hasZ() ) {

        // SCR1 = 0.5 * U(S)**H * DO(Z)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,dcomplex(0.5),UH[idx]->S().pointer(),NB,
          systems_[idx]->onePDMOrtho->Z().pointer(),NB,dcomplex(0.),SCR1,NB);

        // SCR1 += 0.5 * U(Z)**H * DO(S)
        if( UH[idx]->hasZ() )
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,dcomplex(0.5),UH[idx]->Z().pointer(),NB,
            systems_[idx]->onePDMOrtho->S().pointer(),NB,dcomplex(1.),SCR1,NB);

      }






      // DO(S) = 0.5 * ( X(S) * U(S) + X(Z) * U(Z) )
      //       = 0.5 * ( SCR  * U(S) + SCR1 * U(Z) )

      // DO(S) = 0.5 * SCR * U(S)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,dcomplex(0.5),SCR,NB,UH[idx]->S().pointer(),NB,
           dcomplex(0.),systems_[idx]->onePDMOrtho->S().pointer(),NB);
 
      // DO(S) += 0.5 * SCR1 * U(Z)
      if( UH[idx]->hasZ() )
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,dcomplex(0.5),SCR1,NB,UH[idx]->Z().pointer(),NB,
             dcomplex(1.),systems_[idx]->onePDMOrtho->S().pointer(),NB);


      if( UH[idx]->hasZ() ) {

        // DO(Z) = 0.5 * ( X(S) * U(Z) + X(Z) * U(S) )
        //       = 0.5 * ( SCR  * U(Z) + SCR1 * U(S) )
          
        // DO(Z) = 0.5 * SCR * U(Z)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,dcomplex(0.5),SCR,NB,UH[idx]->Z().pointer(),NB,
             dcomplex(0.),systems_[idx]->onePDMOrtho->Z().pointer(),NB);
 
        // DO(Z) += 0.5 * SCR1 * U(S)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,dcomplex(0.5),SCR1,NB,UH[idx]->S().pointer(),NB,
             dcomplex(1.),systems_[idx]->onePDMOrtho->Z().pointer(),NB);

      }

    } else {

      size_t nC = systems_[idx]->nC;

      // Gather DO
      cqmatrix::Matrix<dcomplex> DO(systems_[idx]->onePDMOrtho->template spinGather<dcomplex>());

      // Gather UH into SCR
      cqmatrix::Matrix<dcomplex> UHblockForm(UH[idx]->template spinGather<dcomplex>());

      // SCR1 = U**H * DO
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nC*NB,nC*NB,nC*NB,dcomplex(1.),UHblockForm.pointer(),nC*NB,
           DO.pointer(),nC*NB,dcomplex(0.),SCR1,nC*NB);

      // DO = SCR1 * U
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,nC*NB,nC*NB,nC*NB,dcomplex(1.),SCR1,nC*NB,
           UHblockForm.pointer(),nC*NB,dcomplex(0.),DO.pointer(),nC*NB);

      // Scatter DO
      *systems_[idx]->onePDMOrtho = DO.template spinScatter<dcomplex>();

    }


    systems_[idx]->ortho2aoDen();

    CQMemManager::get().free(SCR,SCR1);

    ProgramTimer::tock("Propagate WFN");

  }; // RealTime::propagatorWFN


  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::createRTDataSets(size_t maxPoints) {

    if( restart ) return;

    if( maxPoints == 0 )
      maxPoints = intScheme.tMax / intScheme.deltaT + 1;
    
    savFile.createGroup("RT");

    savFile.createDataSet<double>("RT/TIME", {maxPoints});
    savFile.createDataSet<double>("RT/ENERGY", {maxPoints});
    savFile.createDataSet<double>("RT/LEN_ELEC_DIPOLE", {maxPoints,3});
    savFile.createDataSet<double>("RT/LEN_ELEC_DIPOLE_FIELD", {maxPoints,3});

    if( this->orbitalPopFreq != 0 ) {
      hsize_t nPop = (maxPoints / this->orbitalPopFreq);
      if( this->orbitalPopFreq != 1 &&
          (maxPoints-1) % this->orbitalPopFreq == 0 )
        nPop += 1;
      if( maxPoints == 1 )
        nPop = 1;
      hsize_t nOrbs = propagator_.nOrbital();
      savFile.createDataSet<double>("RT/ORBITALPOPULATION", {nPop, nOrbs});
    }
  }; // RealTime::createRTDataSets


  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::restoreState() {

    hsize_t maxPoints = intScheme.tMax / intScheme.deltaT + 1;

    if ( savFile.getDims("RT/TIME")[0] != maxPoints )
      CErr("Mismatched requested and saved propagation length!");

    // Restore time dependent density
    try {
      savFile.readData("RT/TD_1PDM", *propagator_.onePDM);
      savFile.readData("RT/TD_1PDM_ORTHO", *propagator_.onePDMOrtho);
    } catch(...) { }

    // Find last time step that was checkpointed
    double* timeData = CQMemManager::get().malloc<double>(maxPoints);
    savFile.readData("RT/TIME", timeData);
    int offset = *timeData < 1e-10 ? -1 : 0;
    size_t restoreStep = offset + std::distance( timeData, 
      std::find_if( timeData+1, timeData+maxPoints,
        [](double x){ return x < 1e-10; }
      )
    );
    CQMemManager::get().free(timeData);

    if( printLevel > 0 ) {
      std::cout << "  *** Restoring from step " << restoreStep << " (";
      std::cout << std::setprecision(4) << restoreStep * intScheme.deltaT;
      std::cout << " AU) ***" << std::endl;
    }

    intScheme.restoreStep = restoreStep;

  }; // RealTime::restoreState


  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::saveState(EMPerturbation& pert_t) {
    

    data.Time.push_back(curState.xTime);
    data.Energy.push_back(propagator_.totalEnergy);
    data.ElecDipole.push_back(propagator_.elecDipole);
    if( pert_t.fields.size() > 0 )
      data.ElecDipoleField.push_back( pert_t.getDipoleAmp(Electric) );

    // Write to file
    if( savFile.exists() ) {

      hsize_t nSteps = 0;
      size_t maxStep = (size_t)((intScheme.tMax + intScheme.deltaT/4)/intScheme.deltaT);

      if (curState.iStep % intScheme.iSave == 0 
          and curState.iStep != intScheme.restoreStep)
        nSteps = intScheme.iSave;
      else if( curState.iStep == maxStep )
        nSteps = (curState.iStep - intScheme.restoreStep) % intScheme.iSave + 1;

      hsize_t lastPos = curState.iStep - nSteps + 1;
      hsize_t memLastPos = data.Time.size() - nSteps;

      if (nSteps != 0) {
        if( printLevel > 0 )
          std::cout << "  *** Saving data to binary file ***" << std::endl;
        savFile.partialWriteData("RT/TIME", &data.Time[0], {lastPos},
            {nSteps}, {memLastPos}, {data.Time.size()});
        savFile.partialWriteData("RT/ENERGY", &data.Energy[0], {lastPos},
            {nSteps}, {memLastPos}, {data.Time.size()});
        savFile.partialWriteData("RT/LEN_ELEC_DIPOLE", &data.ElecDipole[0][0],
          {lastPos, 0}, {nSteps, 3}, {memLastPos, 0}, {data.Time.size(), 3});

        if( data.ElecDipoleField.size() > 0 )
          savFile.partialWriteData("RT/LEN_ELEC_DIPOLE_FIELD",
            &data.ElecDipoleField[0][0], {lastPos,0}, {nSteps,3},
            {memLastPos, 0}, {data.Time.size(), 3});

        savFile.safeWriteData("RT/TD_1PDM", *propagator_.onePDM);
        if ( curState.curStep == ModifiedMidpoint )
          savFile.safeWriteData("RT/TD_1PDM_ORTHO",*DOSav[0]);
        else
          savFile.safeWriteData("RT/TD_1PDM_ORTHO",*propagator_.onePDMOrtho);
      }
    }
  }; // RealTime::saveState

  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::orbitalPop() {

    // Spin-Gather Ortho Density
    std::vector<cqmatrix::Matrix<dcomplex>> orthoDen;
    if( propagator_.nC == 1 ){
       orthoDen = propagator_.onePDMOrtho->template spinGatherToBlocks<dcomplex>(false);
    } else {
       orthoDen.push_back(propagator_.onePDMOrtho->template spinGather<dcomplex>());
    }

    // Transform a copy of the MOs because
    std::vector<cqmatrix::Matrix<dcomplex>> orthoMO = propagator_.mo;

    propagator_.orthoAB->nonortho2orthoCoeffs(orthoMO);

    // Transform alpha Density and compute populations
    size_t NB = orthoMO[0].dimension();
    std::vector<double> population;
    std::vector<cqmatrix::Matrix<dcomplex>> moDen;
    moDen.push_back( orthoDen[0].transform('N',orthoMO[0].pointer(),NB,NB) );
    for( size_t i=0; i<NB; ++i)
        population.push_back( std::real(moDen[0](i,i)) );

    // UHF Beta populations
    if( propagator_.nC == 1 and not propagator_.iCS ){
      moDen.push_back( orthoDen[1].transform('N',orthoMO[1].pointer(),NB,NB) );
      for( size_t i=0; i<NB; ++i)
          population.push_back( std::real(moDen[1](i,i)) );
    }
          
    if( savFile.exists() ) {
      size_t fullDim  = population.size();
      hsize_t location = curState.iStep / this->orbitalPopFreq;
      savFile.partialWriteData("RT/ORBITALPOPULATION", population.data(),
        {location, 0}, {1, fullDim}, {0, 0}, {1, fullDim});
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
        if( propagator_.nC == 1 and not propagator_.iCS ) 
            moDen[1].output(std::cout, "MO Beta Density Matrix", true);
      }

      size_t start = 0;
      if( propagator_.nC == 1 ) {
        printBlock("Alpha occupied orbitals", start, propagator_.nOA);
        printBlock("Alpha virtual orbitals", start, propagator_.nVA);
        if( not propagator_.iCS ){
          printBlock("Beta occupied orbitals", start, propagator_.nOB);
          printBlock("Beta virtual orbitals", start, propagator_.nVB);
        }
      }
      else if( propagator_.nC == 2 ){
        printBlock("Occupied orbitals", start, propagator_.nO);
        printBlock("Virtual orbitals", start, propagator_.nV);
      } else if( propagator_.nC == 4 ){
        start += NB/2;
        printBlock("Positive Energy Occupied orbitals", start, propagator_.nO);
        printBlock("Positive Energy Virtual orbitals", start, propagator_.nV);
      }
      std::cout << std::flush;
    }

  }; // RealTime :: orbitalPop

  /*
  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::orbitalPop() { 

    bool unrestricted = (propagator_.mo.size() > 1);
    bool twocomp = (propagator_.nC == 2);

    // Gather the orthonormal density into the full spinor form
    auto fullDen = propagator_.onePDMOrtho->template spinGather<dcomplex>();

    if( propagator_.nC > 2 ) {
      CErr("Real time orbital population not implemented for > 2c");
    }

    //
    // Form the ortho to MO transform in the full spinor form
    //
    
    // Allocation
    size_t fullDim = fullDen.dimension();
    cqmatrix::Matrix<dcomplex> orthoTrans(fullDim); 
    orthoTrans.clear();
    auto& mos = propagator_.mo;

    size_t moDim = mos[0].dimension();
    auto bbOffset = 2*moDim*moDim + moDim;

    // Transform AO MO to ortho MO
    dcomplex* ortho = propagator_.ortho[1].pointer();
    size_t orthoDim = propagator_.ortho[1].dimension();
    std::vector<dcomplex*> moPointers;
    std::vector<dcomplex*> outPointers;
    for(auto i = 0; i < mos.size(); i++) {
      moPointers.push_back(mos[i].pointer());
      outPointers.push_back(orthoTrans.pointer() + i*bbOffset);
    }
    if(propagator_.iCS) {
      moPointers.push_back(mos[0].pointer());
      outPointers.push_back(orthoTrans.pointer() + bbOffset);
    }

    TransformLeft(orthoDim, moDim, orthoDim, moDim, dcomplex(1.), ortho, orthoDim,
      moPointers, moDim, (dcomplex*)nullptr, outPointers, fullDim);

    // Transform the orthonormal density into the MO basis
    cqmatrix::Matrix<dcomplex> moDen = fullDen.transform('N', orthoTrans.pointer(),
      fullDim, fullDim);

    std::vector<double> population;
    for(auto i = 0; i < fullDim; i++) {
      population.push_back(std::real(moDen(i,i)));
    }

    if( savFile.exists() ) {
      hsize_t location = curState.iStep / this->orbitalPopFreq;
      savFile.partialWriteData("RT/ORBITALPOPULATION", population.data(),
        {location, 0}, {1, fullDim}, {0, 0}, {1, fullDim});
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

      if( this->printLevel > 3 )
        moDen.output(std::cout, "MO Density Matrix", true);

      size_t start = 0;
      if( not twocomp ) {
        printBlock("Alpha occupied orbitals", start, propagator_.nOA);
        printBlock("Alpha virtual orbitals", start, propagator_.nVA);
        printBlock("Beta occupied orbitals", start, propagator_.nOB);
        printBlock("Beta virtual orbitals", start, propagator_.nVB);
      }
      else {
        printBlock("Occupied orbitals", start, propagator_.nO);
        printBlock("Virtual orbitals", start, propagator_.nV);
      }
      std::cout << std::flush;
    }

  };
*/


  template <template <typename, typename> class _SSTyp, typename IntsT>
  void RealTime<_SSTyp,IntsT>::updateAOProperties(double currentTime) {
    // Form AO density
    //propagator_.computeOrtho();
    //propagator_.ortho2aoDen();
    for(auto idx = 0; idx < systems_.size(); idx++) {
      systems_[idx]->computeOrtho();
      systems_[idx]->ortho2aoDen();
      systems_[idx]->ortho2aoMOs();
    }


    // Form fock matrix
    for(auto idx = 0; idx < systems_.size(); idx++) {
      this->formFock(false,currentTime,idx);
    }
    // Compute properties
    EMPerturbation pert_t = pert.getPert(currentTime);
    propagator_.computeProperties(pert_t);
    // Print progress line in the output file
    // printRTStep();


  }; // RealTime::updateAOProperties

}; // namespace ChronusQ

