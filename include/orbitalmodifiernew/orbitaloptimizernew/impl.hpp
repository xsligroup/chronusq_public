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

namespace ChronusQ {

/**
 *  \brief Performs the self-consistant-field procedure given a set of
 *  orbitals.
 *
 *  \warning SCF procedure assumes that the 1PDM and orbital (mo1/2) storage
 *  has been populated in some way.
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::run(EMPerturbation& pert) {

  // Perform Sanity Check of MO's and Eigenvalues
  if( this->singleSlaterSystem.moCoefficients.size() != this->singleSlaterSystem.moEigenvalues.size() ) CErr("Number of MO's and number of Eigenvalues are not equal!");

  ProgramTimer::tick("SCF Total");

  // Initialize type independent parameters
  bool isConverged = false;

  if( scfControls.printLevel > 0 ) {
    printRunHeader(pert);
    printHeaderFinal();
    //printIteration(std::cout, false);
  }

  // Copy Guess onePDM
  // vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
  // for( size_t a = 0; a < onePDM.size(); a++) prevOnePDM[a] = *onePDM[a];

  doingDamp = scfControls.doDamp;

  // If guess is SAD or READDEN, then the first formDensity call can be skipped
  bool skipFormingDensity = (scfControls.guess = SAD) or (scfControls.guess = READDEN);

  // SCF procedure
  for( this->scfConv.nSCFIter = 0; this->scfConv.nSCFIter < scfControls.maxSCFIter; this->scfConv.nSCFIter++ ) {

    ProgramTimer::tick("SCF Iter");

    // If we have a guess density in the first step, then we skip formDensity
    // Otherwise, we form density from MO coefficients
    if( not(this->scfConv.nSCFIter == 0 and skipFormingDensity) ) {
      this->singleSlaterSystem.formDensity();
      // Coefficients and Density represent the same wavefunctions
      this->singleSlaterSystem.setDenEqCoeff(true);
    }

    ProgramTimer::timeOp("Form Fock", [&]() { this->singleSlaterSystem.formFock(pert, false); });

    // Evaluate convergence
    isConverged = evaluateProgress(pert);

    // For energyonly calculations, the orbital energies need be populated.
    // The guess MOs are used to compute the orbital energies
    if(scfControls.energyOnly) this->computeEigenvalues(pert);

    // Save current state of the wave function (method specific)
    this->singleSlaterSystem.saveCurrentState(not scfControls.energyOnly);

    // Print out iteration information
    if( scfControls.printLevel > 0 and (MPIRank(this->mpiComm) == 0) )
      printIteration(this->scfConv.nSCFIter != 0);

    // Exit loop on convergence
    if( isConverged ) break;

    // Get new orbitals and densities from current state:
    //   C/D(k) -> C/D(k + 1)
    this->getNewOrbitals(pert);
    // Coefficients and Density represent different wavefunctions
    this->singleSlaterSystem.setDenEqCoeff(false);

    ProgramTimer::tock("SCF Iter");

  };   // Iteration loop

#ifdef CQ_ENABLE_MPI

      // Broadcast the updated MO's/Eigenvalues to all MPI processes
    if( MPISize(this->mpiComm) > 1 ) {
      for( auto& m : this->singleSlaterSystem.moCoefficients ){
        std::cerr  << "  *** Scattering the MOs ***\n";
        size_t Nmo = m.get().dimension();
        MPIBCast(m.get().pointer(),Nmo*Nmo,0,this->mpiComm);
      }

      std::cerr  << "  *** Scattering EPS ***\n";
      size_t cnt = 0;
      for( auto* e : this->singleSlaterSystem.moEigenvalues ){
        size_t Nmo = this->singleSlaterSystem.moCoefficients[cnt].get().dimension();
        MPIBCast(e,Nmo,0,this->mpiComm);
        ++cnt;
      }
    }
#endif

  // Populate moFock after SCF
  this->singleSlaterSystem.MOFOCK();
  
  // Compute final properties
  this->singleSlaterSystem.computeProperties(pert);

  // printSCFFooter(isConverged);
  if( not isConverged )
    CErr(std::string("SCF Failed to converged within ") + std::to_string(scfControls.maxSCFIter) + std::string(" iterations"));
  else if( scfControls.printLevel > 0 ) {
    std::cout << std::endl
              << "SCF Completed: E(" << scfControls.refShortName_ << ") = " << std::fixed << std::setprecision(10) << this->singleSlaterSystem.getTotalEnergy() << " Eh after "
              << this->scfConv.nSCFIter << " SCF Iterations" << std::endl;
  }

  if( scfControls.printLevel > 0 ) std::cout << BannerEnd << std::endl;

  if( scfControls.printLevel > 1 ) this->singleSlaterSystem.printProperties();

  // Print MO occupations 
  //this->singleSlaterSystem.printOrbitalPopulation(std::cout);

  ProgramTimer::tock("SCF Total");

  // Save final results to bin file
  this->singleSlaterSystem.saveCurrentState(not scfControls.energyOnly);

};   // OrbitalOptimizer<MatsT,IntsT>::SCF()

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::printRunHeader(EMPerturbation& pert) {

  std::cout << BannerTop << std::endl;
  std::cout << "Self Consistent Field (SCF) Settings:" << std::endl << std::endl;

  std::cout << std::setw(38) << std::left << "  Reference:" << scfControls.refLongName_ << std::endl;

  std::cout << std::setprecision(6) << std::scientific;
  // out << std::setw(38) << std::left << "  Density Convergence Tolerance:" << scfControls.denConvTol << std::endl;

  std::cout << std::setw(38) << std::left << "  Energy Difference:" << scfControls.eneConvTol << std::endl;

  std::cout << std::setw(38) << std::left << "  RMS Density Difference:" << scfControls.rmsdPConvTol << std::endl;

  std::cout << std::setw(38) << std::left << "  Max Density Difference:" << scfControls.maxdPConvTol << std::endl;

  std::cout << std::setw(38) << std::left << "  Maximum Number of SCF Cycles:" << scfControls.maxSCFIter << std::endl;

  std::cout << std::setw(38) << std::left << "  SCF Algorithm:";
  if( scfControls.scfAlg == _CONVENTIONAL_SCF )
    std::cout << "Conventional SCF";
  else if( scfControls.scfAlg == _NEWTON_RAPHSON_SCF )
    std::cout << "Newton-Raphson (2nd Order)";
  else if( scfControls.scfAlg == _SKIP_SCF )
    std::cout << "Skip SCF";
  std::cout << "\n";

  // Field print
  if( pert.fields.size() != 0 ) {

    std::cout << "\n\n  * SCF will be performed in the presence of an EM "
        << "perturbation:\n\n";

    for( auto& field : pert.fields ) {

      auto amp = field->getAmp();

      std::cout << "     * ";
      if( field->emFieldTyp == Electric )
        std::cout << "Electric";
      else
        std::cout << "Magnetic";

      std::cout << " ";

      if( field->size == 3 )
        std::cout << "Dipole";
      else if( field->size == 6 )
        std::cout << "Quadrupole";
      else if( field->size == 10 )
        std::cout << "Octupole";

      std::cout << " Field: ";
      std::cout << "{ ";
      for( auto i = 0; i < amp.size(); i++ ) {
        std::cout << amp[i];
        if( i != amp.size() - 1 ) std::cout << ", ";
      }
      std::cout << " }\n";
    }
  }

};   // OrbitalOptimizer<MatsT>::printSCFHeader

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>:: printHeaderFinal() const {
  std::cout << std::endl << BannerMid << std::endl << std::endl;
  std::cout << std::setw(16) << std::right  << "SCF Iteration";
  std::cout << std::setw(18) << std::right  << "Energy (Eh) ";
  std::cout << std::setw(18) << std::right  << " \u0394E (Eh)   ";  // Increase width by 1 to account for unicode characters
  std::cout << std::setw(17) << std::right  << "   RMS|\u0394D| ";
  std::cout << std::setw(18) << std::right  << " Max|\u0394D|";
  std::cout << std::endl;

  std::cout << std::setw(16) << std::right  << "-------------";
  std::cout << std::setw(18) << std::right  << "----------- ";
  std::cout << std::setw(16) << std::right  << "-------  ";
  std::cout << std::setw(16) << std::right  << "-------";
  std::cout << std::setw(18) << std::right  << "-------";
  std::cout << std::endl;
}

/**
 *  \brief Print the current convergence information of the SCF
 *  procedure
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::printIteration(bool printDiff) {

  // SCF Iteration
  std::cout << std::setw(10) << std::right << "SCFIt:";

  if( printDiff )
    std::cout << std::setw(6) << std::right << this->scfConv.nSCFIter;
  else
    std::cout << std::setw(6) << std::right << 0;

  // Current Total Energy
  std::cout << std::setw(18) << std::fixed << std::setprecision(10) << std::right << this->singleSlaterSystem.getTotalEnergy();

  if( printDiff ) {
    std::cout << std::scientific << std::setprecision(7);
    // Current Change in Energy
    std::cout << std::setw(18) << std::right << this->scfConv.deltaEnergy;
    std::cout << std::setw(18) << std::right << this->scfConv.rmsdP;
    std::cout << std::setw(18) << std::right << this->scfConv.maxdP;
  }

  std::cout << std::endl;
};   // OrbitalOptimizer<MatsT>::printSCFProg

/**
 *  \brief Evaluate SCF convergence based on various criteria.
 *
 *  Checks the norm of [F,D], if converged -> SCF converged.
 *
 *  Checks change in energy and density between SCF iterations,
 *    if *both* converged -> SCF converged.
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
bool OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::evaluateProgress(EMPerturbation& pert) {

  bool isConverged = true;

  // Compute all SCF convergence information on root process
  if( MPIRank(this->mpiComm) == 0 ) {

    // Compute new energy (with new Density)
    this->singleSlaterSystem.computeEnergy(pert);
    scfConv.currentEnergy = this->singleSlaterSystem.getTotalEnergy();

    if(not scfControls.energyOnly) {

      // Check energy convergence
      bool energyConv = false;
      if(this->scfConv.nSCFIter == 0) {
        scfConv.previousEnergy = scfConv.currentEnergy;
      } else {
        this->scfConv.deltaEnergy = scfConv.currentEnergy - scfConv.previousEnergy;
        energyConv = this->scfConv.deltaEnergy < scfControls.smallEnergy;
        energyConv = energyConv and (std::abs(this->scfConv.deltaEnergy) < scfControls.eneConvTol);
        scfConv.previousEnergy = scfConv.currentEnergy;
      }

      // Check density convergence
      bool denConv = false;
      if(this->scfConv.nSCFIter == 0) {
        // Allocate prevOnePDM
        vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
        for( size_t a = 0; a < onePDM.size(); a++ ) {
          prevOnePDM.emplace_back(onePDM[a]->dimension());
          prevOnePDM[a] = *onePDM[a];
        }
      } else {
        vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
        this->scfConv.rmsdP = 0.;
        this->scfConv.maxdP = 0.;
        // Here, onePDM is in the full spin-block form
        for( size_t a = 0; a < onePDM.size(); a++ ) {
          cqmatrix::Matrix<MatsT> dDen = *onePDM[a] - prevOnePDM[a];
          prevOnePDM[a] = *onePDM[a];
          size_t NB = onePDM[a]->dimension();
          this->scfConv.rmsdP += blas::nrm2(NB*NB,dDen.pointer(),1) / NB;
          // std::abs returns the norm of a complex number
          for( size_t b = 0; b < NB*NB; b++) {
            if(this->scfConv.maxdP < std::abs(*(dDen.pointer()+b)))
              this->scfConv.maxdP = std::abs(*(dDen.pointer()+b));
          };
        }
        denConv = this->scfConv.rmsdP < scfControls.rmsdPConvTol;
        denConv = denConv and (this->scfConv.maxdP < scfControls.maxdPConvTol);
      }

      isConverged = energyConv and denConv;

      // Toggle damping based on energy difference
      if(scfControls.doDamp ) {
        bool largeEDiff = std::abs(this->scfConv.deltaEnergy) > scfControls.dampError;

        if( doingDamp and not largeEDiff and scfControls.dampParam > 0. ) {

          if( scfControls.printLevel > 0 )
            std::cout << "    *** Damping Disabled - Energy Difference Fell Below " << scfControls.dampError << " ***" << std::endl;

          doingDamp = false;

        } else if( not doingDamp and largeEDiff and scfControls.dampParam <= 0. ) {

          if( scfControls.printLevel > 0 )
            std::cout << "    *** Damping Enabled Due to " << scfControls.dampError << " Oscillation in Energy ***" << std::endl;

          doingDamp = true;
        }
      }
    }
  }

#ifdef CQ_ENABLE_MPI
  // Broadcast whether or not we're converged to ensure that all
  // MPI processes exit the SCF simultaneously
  if( MPISize(this->mpiComm) > 1 ) MPIBCast(isConverged, 0, this->mpiComm);
#endif

  return isConverged;

};   // OrbitalOptimizer<MatsT>::evalConver

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
double OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::computeDensityConv() {
    vecShrdPtrMat<MatsT> currDen = this->singleSlaterSystem.getOnePDM();

    // Compute RMS change in Density
    double rmsDen = 0.;
    for( size_t a=0; a<currDen.size(); a++ ) {
        cqmatrix::Matrix<MatsT> dDen = *currDen[a] - prevOnePDM[a];
        size_t NB = currDen[a]->dimension();
        rmsDen += blas::nrm2(NB*NB,dDen.pointer(),1) / NB;
        prevOnePDM[a] = *currDen[a];
    }
    return rmsDen;
};  // OrbitalOptimizer<MatsT> :: computeDensityConv


/*
 *  Brief: Computes the orbital energies from a set of MO's
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>::computeEigenvalues(EMPerturbation& pert) {

  vecShrdPtrMat<MatsT> fock = this->singleSlaterSystem.getFock();

  for( size_t i = 0; i < this->singleSlaterSystem.moCoefficients.size(); i++ ) {
    size_t NB = fock[i]->dimension();

    cqmatrix::Matrix<MatsT> moFock = fock[i]->transform('N', this->singleSlaterSystem.moCoefficients[i].get().pointer(), NB, NB);
    for( size_t a = 0; a < NB; a++ )
      this->singleSlaterSystem.moEigenvalues[i][a] = std::real(moFock(a, a));
  }
};


};   // namespace ChronusQ
