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
#include <mcscf.hpp>
#include <mcscf/print.hpp>
#include <mcscf/rdm.hpp>
#include <mcscf/cisolver.hpp>
#include <util/matout.hpp>
#include <orbitalrotation.hpp>

// #define DEBUG_MCSCF_IMPL

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::run(EMPerturbation & externalPert) {
    
    ProgramTimer::tick("MCSCF Total");
    // Create combined perturbation
    EMPerturbation pert;
    // Add on the MCSCF Perturbation if present
    for( auto& field : this->mcscfPert.fields )
      pert.addField( field );
    // Finally add any additional Perturbations
    for( auto& field : externalPert.fields )
      pert.addField( field );

    // allocating memeory
    this->alloc(); 
    
    // Initial Printing
    this->printMCSCFHeader(pert);
    
    // MCSCF Intial CI solution
    std::cout << "Cycle 0:\n" << std::endl;
    FormattedLine(std::cout, "AO to MO Intergral Transformation ...");
    
    ProgramTimer::tick("Solve CI");
    
    ProgramTimer::tick("Integral Trans");
    if (!this->readCI or this->settings.doSCF)
      MCWaveFunction<MatsT,IntsT>::transformInts(pert);
    ProgramTimer::tock("Integral Trans");
    
    std::cout << std::left << std::setprecision(10); 
    FormattedLine(std::cout, "Inactive Energy:", this->InactEnergy);

    // Precompute the field - nuclear diagonal contributions
    precompute_NucEField(pert);

    ProgramTimer::tick("Diagonalization");
    if (!this->readCI)
      this->ciSolver->solveCI(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), pert);
    ProgramTimer::tock("Diagonalization");
    
    ProgramTimer::tock("Solve CI");

    // Initial 1RDM construction
    MCWaveFunction<MatsT, IntsT>::computeOneRDM();
     
    // MCSCF Cycles 
    
    if(this->settings.doSCF) {

      this->printStateEnergy();
      
      std::vector<double> EPrev = std::vector<double>(this->NStates);
      std::fill_n(EPrev.begin(), this->NStates, 0.);
      double EDiff      = 0.;
      bool   converged  = false;
      
      for (auto iter = 0ul; iter < settings.maxSCFIter; iter++) {
        
        ProgramTimer::tick("Orbital Rotation");
        // Exam Energy
        if(this->StateAverage) {
          EDiff = 0.;
          std::vector<double> EDiff2 = std::vector<double>(this->NStates);
          for (auto i = 0ul; i < this->NStates; i++)
            EDiff2[i] = this->StateEnergy[i] - EPrev[i];
          
          EDiff = *std::max_element(EDiff2.begin(), EDiff2.end(), 
            [&] (double a, double b) { return std::abs(a) < std::abs(b); });

        } else {
          EDiff = this->StateEnergy.back() - EPrev.back(); 
        }
        
        std::cout << "  (Maximum) Energy Difference = " << std::setw(18) 
                  <<  std::right << EDiff << std::left << std::endl; 
        if(std::abs(EDiff) <= settings.scfEnergyConv) converged = true; 
        
        // compute RDMs
        if (this->StateAverage) {
          this->computeOneRDM(); 
          this->computeTwoRDM(); 
        } else {
          this->computeOneRDM(this->NStates - 1); 
          this->computeTwoRDM(this->NStates - 1);
        }
        
        // this->print1RDMs();
        
        // compute gradient 
        double orbitalGradientNorm = 
          moRotator->computeOrbGradient(pert, *oneRDMSOI, *twoRDMSOI);
        
        // Gradient Convergence exam
        std::cout << "  Orbital Gradient Residue   = " << std::setw(18) 
                  << std::right << orbitalGradientNorm << std::left << std::endl;
        
        if(converged and orbitalGradientNorm < settings.scfGradientConv) break; 
          
        converged = false;
        
        // start of the new cycle 
        std::cout << "\n\nCycle " << iter+1 << ":\n" << std::endl;
        
        // compute hession diagonal and rotate orbitals
        FormattedLine(std::cout, "Performing Orbital Rotation ...");
        
        moRotator->rotateMO(pert, *oneRDMSOI, *twoRDMSOI);
        
        ProgramTimer::tock("Orbital Rotation");
        
        this->mointsTF->clearAllCache();

        // print energy and update EPrev
        std::copy_n(this->StateEnergy.begin(), this->NStates, EPrev.begin());

        ProgramTimer::tick("Solve CI");
        
        // Re-transform intgrals and solve new CI
        FormattedLine(std::cout, "Redo AO to MO Intergral Transformation ...");
        ProgramTimer::tick("Integral Trans");
        MCWaveFunction<MatsT,IntsT>::transformInts(pert);
        // Don't need to call precompute_NucEField here since pert cannot
        // change between above and here
        ProgramTimer::tock("Integral Trans");

        std::cout << std::left << std::setprecision(10); 
        FormattedLine(std::cout, "Inactive Energy:", this->InactEnergy);
        ProgramTimer::tick("Diagonalization");
        this->ciSolver->solveCI(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), pert);
        ProgramTimer::tock("Diagonalization");
        
        ProgramTimer::tock("Solve CI");
      
        this->printStateEnergy();
      
        saveCurrentStates();
      
      } // SCF Iteration

      if(not converged) 
        CErr("\n MCSCF failed to converged in " + std::to_string(settings.maxSCFIter) + " cycles !");
    
      // compute 1RDMs
      if (this->StateAverage) this->computeOneRDM(); 
      else this->computeOneRDM(this->NStates - 1); 
      
      // generate IVOs as needed
      ProgramTimer::tick("Gen IVOs");
      if (this->settings.doIVOs) moRotator->generateIVOs(pert, *oneRDMSOI);   
      ProgramTimer::tock("Gen IVOs");

    } // doSCF

    // Final printing and save states for restart
    std::cout << "\n\nMCSCF Complete!" << std::endl;
    std::cout << bannerEnd << std::endl;

    this->printMCSCFFooter();

    // property calculation
    ProgramTimer::tick("Property Eval");

    // dipole moment
    if( this->multipoleMoment )
      MCWaveFunction<MatsT,IntsT>::computeMultipole();

    // mulliken analysis
    if (this->PopulationAnalysis) {
      std::cout<<"\n\nPopulation analysis in mcscf."<<std::endl;
      MCWaveFunction<MatsT,IntsT>::populationAnalysis();
    }

    // oscillator strength
    if (this->NosS1) {


      this->osc_str = this->memManager.template malloc<double>(this->NosS1*this->NStates);
      for (size_t s1 = 0ul; s1 < this->NosS1; s1++)
      for (size_t s2 = 0ul; s2 < this->NStates; s2++){
//        if (s2 < this->NosS1) this->osc_str[s2+s1*this->NStates] = 0.;
        if (s2 <= s1) this->osc_str[s2+s1*this->NStates] = 0.;
        else this->osc_str[s2+s1*this->NStates] = 
                MCWaveFunction<MatsT,IntsT>::oscillator_strength(s2,s1);
      }

    }

    ProgramTimer::tock("Property Eval");

    saveCurrentStates(true);

    ProgramTimer::tock("MCSCF Total");
 
  }; //MCSCF::run
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::saveCurrentStates( bool saveProp ) {
    
    MCWaveFunction<MatsT, IntsT>::saveCurrentStates(saveProp);
    
    // only save MO when doing orbital rotation
    if (settings.doSCF and this->savFile.exists()) {
      auto mo_dim = this->reference().mo[0].dimension();
      this->savFile.safeWriteData("SCF/MO1", this->reference().mo[0].pointer(), {mo_dim, mo_dim});
    }

  }; // MCSCF::saveCurrentStates


  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::alloc() {
    
    MCWaveFunction<MatsT,IntsT>::alloc();
    
    ciSolver = std::make_shared<CISolver<MatsT,IntsT>>(settings.ciAlg, 
      settings.maxCIIter, settings.ciVectorConv,
      settings.maxDavidsonSpace, settings.nDavidsonGuess,
      settings.energyRefs);
    
    if (this->settings.doSCF) {
      
      this->cacheHalfTransTPI_ = true;
      
      size_t nCorrO = this->MOPartition.nCorrO;
      oneRDMSOI = std::make_shared<SquareMatrix<MatsT>>(this->memManager,nCorrO);
      twoRDMSOI = std::make_shared<InCore4indexTPI<MatsT>>(this->memManager, nCorrO);     
      
      moRotator = std::make_shared<OrbitalRotation<MatsT, IntsT>>(
        dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), settings.ORSettings);
    }
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::dealloc() {
    oneRDMSOI = nullptr;
    twoRDMSOI = nullptr;
    ciSolver  = nullptr;
    moRotator = nullptr;
  }

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::precompute_NucEField(EMPerturbation & pert)
  {
    // Zero out in case this has already been calculated & stored
    // (for example from RT-CI)
    this->EFieldNuc = 0.0;

    if(!pert_has_type(pert,Electric))
      return;

    std::array<double,3> nucmoment = {0.0,0.0,0.0};
    for(auto & atom : this->reference().molecule().atoms)
    {
      if(atom.quantum) continue;
      MatAdd('N','N',3,1,1.,&nucmoment[0],3,atom.nucCharge,&atom.coord[0],3,&nucmoment[0],3);
    }
      
    auto elecDipoleField = pert.getDipoleAmp(Electric);
    this->EFieldNuc+=blas::dot(3,&nucmoment[0],1,&elecDipoleField[0],1);
  
  // The code below takes an alternate approach:  Instead of folding the field
  // contributions into the one electron integrals, it calculates the dipole 
  // moment for each individual basis function (in this case, slater determinants)
  // as well as the cross terms between them (ie, E\dot<SD_i|x,y,z|SD_j>).
  // This is far less computationally efficient than folding the field contribution
  // into the integrals, but might have future use if say the dipole moment of a 
  // single non-aufbau slater determinant is required for some reason.

  // To use the below code, you must also make corresponding changes in the 
  // hcore portion of mointstransformer in order to prevent the field being
  // double counted by both folding into the integrals and calculating each 
  // slater determinants dipole moment

//    if(this->reference().nC!=1)
//      CErr("CI Dipole Diagonal EField for nC!=1 NYI");
//
//    SingleSlater<MatsT,IntsT> * ss_ptr = &(this->reference());
//
//    size_t NDet = this->NDet;
//    size_t nAO = this->reference().nAlphaOrbital();
//    size_t nI = this->MOPartition.nFCore + this->MOPartition.nInact;
//    size_t nCorrO = this->MOPartition.nCorrO;
//
//    this->EFieldDiag = this->memManager.template malloc<MatsT>(NDet*NDet);
//    std::fill_n(this->EFieldDiag,NDet*NDet,MatsT(0.0));
//
//    MatsT * dummyCIi = this->memManager.template malloc<MatsT>(NDet);
//    MatsT * dummyCIj = this->memManager.template malloc<MatsT>(NDet);
//    SquareMatrix<MatsT> dummyRDM(this->memManager,NDet);
//    SquareMatrix<MatsT> dummyRDM2(this->memManager,NDet);
//    SquareMatrix<MatsT> dummyPDM(this->memManager,nAO);
//
//    auto elecDipoleField = pert.getDipoleAmp(Electric);
//    
//    // Loop over determinants
//    for(size_t i = 0; i < NDet; i++)
//    {
//      // Zero out the CI vector and RDM
//      std::fill_n(dummyCIi,NDet,MatsT(0.0));
//      dummyCIi[i]=1.0;
//      for(size_t j = 0; j < NDet; j++)
//      {
//        if(i==j)
//        {
//          std::fill_n(dummyCIj,NDet,MatsT(0.0));
//          dummyCIj[j]=1.0;
//          dummyRDM.clear();
//          dummyPDM.clear();
//          for(size_t k = 0; k < 3; k++)
//            ss_ptr->elecDipole[k] = 0.0;
//
//          // Fill in the dummy CI Vector and turn it into the RDM
//          this->ciBuilder->computeTDM(*this,dummyCIi,dummyCIj,dummyRDM);
//
//          // Convert the RDM to a PDM
//          // Fold in core orbitals
//          for(size_t k = 0; k < nI; k++)
//            dummyPDM(k,k) = 2.0;
//          // Add RDM contributions
//          SetMat('R',nCorrO,nCorrO,1.0,dummyRDM.pointer(),nCorrO,dummyPDM.pointer()+nI*(nAO+1),nAO);
//          // Transform
//          dummyPDM = dummyPDM.transform('C',this->reference().mo[0].pointer(),nAO,nAO);
//          // Set the reference PDM to this PDM
//          this->reference().onePDM->S() = dummyPDM;
//
//          ss_ptr->computeMultipole(pert);
//
//          this->EFieldDiag[i+j*NDet] = ss_ptr->elecDipole[0]*elecDipoleField[0]+
//                                      ss_ptr->elecDipole[1]*elecDipoleField[1]+
//                                      ss_ptr->elecDipole[2]*elecDipoleField[2];
//        }
//        else
//        {
//          std::fill_n(dummyCIj,NDet,MatsT(0.0));
//          dummyCIj[j]=1.0;
//          dummyRDM.clear();
//          dummyRDM2.clear();
//          dummyPDM.clear();
//          // Fill in the dummy CI Vector and turn it into the RDM
//          this->ciBuilder->computeTDM(*this,dummyCIi,dummyCIj,dummyRDM);
//          this->ciBuilder->computeTDM(*this,dummyCIj,dummyCIi,dummyRDM2);
//
//          auto MOdipole = this->moints.template getIntegral<VectorInts,MatsT>("MOdipole");
//          if (!MOdipole) {
//            std::shared_ptr<VectorInts<IntsT>> AOdipole =
//                      std::make_shared<VectorInts<IntsT>>(this->memManager, nAO, 1, true);
//            std::shared_ptr<VectorInts<MatsT>> MOdipole_scr =
//                      std::make_shared<VectorInts<MatsT>>(this->memManager, nCorrO, 1, true);
//
//            std::vector<std::pair<size_t, size_t>> active(2, {this->MOPartition.nFCore+this->MOPartition.nInact, nCorrO});
//
//            for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
//              (*AOdipole)[iXYZ] = std::make_shared<OnePInts<IntsT>>( *((*this->reference().aoints_->lenElectric)[iXYZ]) );
//              
//              (*AOdipole)[iXYZ]->subsetTransform('N',this->reference().mo[0].pointer(),
//                  nAO, active, (*MOdipole_scr)[iXYZ]->pointer(), false);
//            }
//
//            this->moints.addIntegral("MOdipole", MOdipole_scr);
//          }
//
//          MOdipole = this->moints.template getIntegral<VectorInts,MatsT>("MOdipole");
//
//          MatsT Etemp = 0.0;
//          for(size_t ixyz = 0; ixyz < 3; ixyz++)
//          {
//            Etemp += elecDipoleField[ixyz]*blas::dotu(nCorrO*nCorrO,dummyRDM.pointer(),1,(*MOdipole)[ixyz]->pointer(),1);
//          }
//          this->EFieldDiag[i+j*NDet]=-Etemp;
//          this->EFieldDiag[j+i*NDet]=-Etemp;
// 
//        }
//      }
//    }
//    // prettyPrintSmart(std::cout,"Matrix-Dipole additions",this->EFieldDiag,NDet,NDet,NDet);
//
  }
}; // namespace ChronusQ

