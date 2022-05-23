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

#include <cxxapi/input.hpp>
#include <cxxapi/output.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>
#include <cxxapi/procedural.hpp>

#include <util/files.hpp>
#include <util/matout.hpp>
#include <util/mpi.hpp>
#include <util/threads.hpp>
#include <util/timer.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <integrals.hpp>
#include <singleslater.hpp>
#include <response.hpp>
#include <realtime.hpp>
#include <itersolver.hpp>
#include <coupledcluster.hpp>
#include <mcwavefunction.hpp>
#include <mcscf.hpp>

#include <findiff/geomgrad.hpp>
#include <particleintegrals/gradints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/gradints/incore.hpp>
#include <particleintegrals/gradints/direct.hpp>

#include <cqlinalg/blasext.hpp>
#include <cqlinalg/eig.hpp>

#include <geometrymodifier.hpp>
#include <geometrymodifier/moleculardynamics.hpp>
#include <geometrymodifier/singlepoint.hpp>
#include <physcon.hpp>

#include <corehbuilder/x2c.hpp>
#include <corehbuilder/nonrel.hpp>
#include <fockbuilder/matrixfock.hpp>

#include <fockbuilder/neofock.hpp>

//#include <cubegen.hpp>

namespace ChronusQ {

  template class NEOKohnShamBuilder<double,double>;
  template class NEOKohnShamBuilder<dcomplex,double>;
  template class NEOKohnShamBuilder<dcomplex,dcomplex>;

#ifdef ENABLE_BCAST_COUNTER
  int bcastCounter = 0;
#endif

  void RunChronusQ(std::string inFileName,
    std::string outFileName, std::string rstFileName,
    std::string scrFileName) {

    // Check to make sure input and output file name are different.
    if( inFileName == outFileName )
      CErr("Input file name and output file name cannot be identical.");

    int rank = MPIRank();
    int size = MPISize();

    // Redirect output to output file if not STDOUT
    std::shared_ptr<std::ofstream> outfile;
    std::streambuf *coutbuf = std::cout.rdbuf();

    if( outFileName.compare("STDOUT") and (rank == 0) ) {

      outfile = std::make_shared<std::ofstream>(outFileName);
      std::cout.rdbuf(outfile->rdbuf());

    }


    // Setup MPI rank files

    std::shared_ptr<std::ofstream> rankfile;
    std::streambuf *cerrbuf = std::cerr.rdbuf();

    if( size > 1 ) {
      
      std::string rankFileName = outFileName + ".mpi." + 
        std::to_string(rank);

      rankfile = std::make_shared<std::ofstream>(rankFileName);
      std::cerr.rdbuf(rankfile->rdbuf());

      std::cerr << "Hello from RANK = " << rank << " / SIZE = " << size 
                << "\n\n";

    }


    std::ostream &output = (rank == 0) ? std::cout : std::cerr;


    // Output CQ header
    CQOutputHeader(output);
    if(rankfile and rank == 0) CQOutputHeader(std::cerr);


    // Parse Input File
    CQInputFile input(inFileName);


    CQINPUT_VALID(output,input);

    // Dump contents of input file into output file
    if( rank == 0 ) {
      std::cout << "\n\n\n";
      std::cout << "Input File:\n" << BannerTop << std::endl;
      std::ifstream inStream(inFileName);
      std::istreambuf_iterator<char> begin_src(inStream);
      std::istreambuf_iterator<char> end_src;
      std::ostreambuf_iterator<char> begin_dest(std::cout);
      std::copy(begin_src,end_src,begin_dest);
      inStream.close();
      std::cout << BannerEnd << "\n\n\n" << std::endl;
    }


    // TEMPORARY
    bool doTemp = true;

    // Determine JOB type
    JobType jobType;
    
    try {
      jobType = parseJob(input.getData<std::string>("QM.JOB"));
    } catch (...) {
      CErr("Must Specify QM.JOB",output);
    }

    // Break into sequence of individual jobs
    std::vector<JobType> jobs;
    if( jobType != SCF ) {
      jobs.push_back(SCF);
    }
    jobs.push_back(jobType);

    // Check if we're doing NEO
    bool doNEO = false;
    if ( input.containsSection("SCF") ) {
      try {
        doNEO = input.getData<bool>("SCF.NEO");
      } catch(...) { ; }
    }

    auto memManager = CQMiscOptions(output,input);


    // Create Molecule and BasisSet objects
    Molecule mol(std::move(CQMoleculeOptions(output,input,scrFileName)));

    // Create BasisSet objects
    std::shared_ptr<BasisSet> basis = CQBasisSetOptions(output,input,mol,"BASIS");
    std::shared_ptr<BasisSet> dfbasis = CQBasisSetOptions(output,input,mol,"DFBASIS");

    auto aoints = CQIntsOptions(output,input,*memManager,mol,basis,dfbasis,nullptr);

    // Create BasisSet and integral objects for nuclear orbitals 
    std::shared_ptr<BasisSet> prot_basis = 
      doNEO ? CQBasisSetOptions(output,input,mol,"PBASIS") : nullptr;
    auto prot_aoints = 
      doNEO ? CQIntsOptions(output,input,*memManager,mol,prot_basis,dfbasis,nullptr,"PINTS"): nullptr;
    auto ep_aoints   = 
      doNEO? CQIntsOptions(output,input,*memManager,mol,basis,dfbasis,prot_basis,"EPINTS") : nullptr;

    std::shared_ptr<SingleSlaterBase> ss  = nullptr;

    SingleSlaterOptions ssOptions;

    // EM Perturbation for SCF
    EMPerturbation emPert;

    // SCF options
    SCFControls scfControls = CQSCFOptions(output,input,emPert);

    // Create the SingleSlater object
    if (doNEO) {
      std::tie(ss, ssOptions) = CQNEOSSOptions(output,input,*memManager,mol,
                                              *basis,*prot_basis,
                                               aoints, prot_aoints,
                                               ep_aoints);
      ss->scfControls = scfControls;
      // For NEO only one ModifyOrbitals needs to be made since it is a
      // driver for both the NEOSingleSlater and the aux_neoss
      ss->buildModifyOrbitals();
    } else {
      ssOptions = CQSingleSlaterOptions(output,input,mol,*basis,aoints);
      ssOptions.scfControls = scfControls;
      ss = ssOptions.buildSingleSlater(output,*memManager,mol,*basis,aoints);
      ss->buildModifyOrbitals();

      // MO swapping
      HandleOrbitalSwaps(output, input, *ss);
    }

    bool rstExists = false;
    if( (ss->scfControls.guess == READMO or 
         ss->scfControls.guess == READDEN or
         ss->scfControls.prot_guess == READMO or
         ss->scfControls.prot_guess == READDEN)
        and scrFileName.empty())
      rstExists = true;
    else if( (ss->scfControls.guess == READDEN or
              ss->scfControls.prot_guess == READDEN)
             and not scrFileName.empty() )
      ss->scrBinFileName = scrFileName;
    else if( ss->scfControls.guess == FCHKMO or
             ss->scfControls.prot_guess == FCHKMO )
      ss->fchkFileName = scrFileName;

    // Create the restart and scratch files
    SafeFile rstFile(rstFileName, rstExists);

    if( not rstExists and rank == 0 ) rstFile.createFile();


    if( rank == 0 ) {
      ss->savFile     = rstFile;
      aoints->savFile = rstFile;
      if (doNEO) { 
        prot_aoints->savFile = rstFile;
        ep_aoints->savFile   = rstFile;
      }
    }

    // If doing NEO, propagate setup to subsystems
    if(auto neoss = std::dynamic_pointer_cast<NEOBase>(ss)) {
      neoss->setSubSetup();
    }

    // Done setting up
    //
    // START OF REAL PROCEDURAL SECTION

    for( auto& job: jobs ) {

      bool firstStep = true;
      std::shared_ptr<RealTimeBase> rt = nullptr;

//      if (ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::OFF) {
//        compute_X2C_CoreH_Fock(*memManager, mol, *basis, aoints, emPert, ss, ssOptions);
//      }

      JobType elecJob = CQGeometryOptions(output, input, job, mol, ss, rt,
        ep_aoints, emPert);

      // Loop over various structures
      while( mol.geometryModifier->hasNext() ) {

        // Update geometry
        mol.geometryModifier->electronicPotentialEnergy = ss->totalEnergy;
        mol.geometryModifier->update(true, mol, firstStep);
        // Update basis to the new geometry
        basis->updateNuclearCoordinates(mol);
        if( dfbasis != nullptr ) dfbasis->updateNuclearCoordinates(mol);

        if( doNEO ) {
          prot_basis->updateNuclearCoordinates(mol);
        }

        // Update integrals
        // TODO: Time dependent field?
        aoints->computeAOTwoE(*basis, mol, emPert);

        if (doNEO) { 
          if(auto p = std::dynamic_pointer_cast<Integrals<double>>(prot_aoints))
            prot_aoints->computeAOTwoE(*prot_basis, mol, emPert);
          else
            CErr("NEO with complex integrals NYI!",output);
          if(auto p = std::dynamic_pointer_cast<Integrals<double>>(ep_aoints))
            ep_aoints->computeAOTwoE(*basis, *prot_basis, mol, emPert); 
        }

        // Note, these guessSSOptions does not apply to NEO guess
        SingleSlaterOptions guessSSOptions(ssOptions);
        guessSSOptions.refOptions.isKSRef = false;
        guessSSOptions.refOptions.nC = 1;
        guessSSOptions.hamiltonianOptions.OneEScalarRelativity = false;
        guessSSOptions.hamiltonianOptions.OneESpinOrbit = false;

        // Run SCF job
        if( elecJob == SCF ) {

          if (ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::OFF) {
            compute_X2C_CoreH_Fock(*memManager, mol, *basis, aoints, emPert, ss, ssOptions);
          }

          ss->formCoreH(emPert, true);
          if(firstStep) {
            ss->formGuess(guessSSOptions);
            ss->formFock(emPert, false);
          }
          ss->runModifyOrbitals(emPert);
        }

        // Run RT job
        if( elecJob == RT ) {

          // Initialize core hamiltonian
          rt->formCoreH(emPert);

          // Get correct time length
          if( !firstStep ) {
            rt->intScheme.restoreStep = rt->curState.iStep;
            rt->intScheme.tMax = rt->intScheme.tMax + rt->intScheme.nSteps*rt->intScheme.deltaT;
          }

          if( MPISize() > 1 ) CErr("RT + MPI NYI!",output);

          rt->doPropagation();

        }


        if( elecJob == RESP ) {


          if( ss->scfControls.scfAlg == _SKIP_SCF and ss->scfControls.guess == READDEN )
            CErr("READDEN + SKIP + RESPONSE disabled. Use READMO instead.");

          auto resp = CQResponseOptions(output,input,ss);
          resp->savFile = rstFile;
          resp->run();

          if( MPIRank(MPI_COMM_WORLD) == 0 ) resp->printResults(output);
          MPI_Barrier(MPI_COMM_WORLD);

        }

        
        if( elecJob == CC ){

          // FIXME: Need to implement NEO-CC
          if (doNEO)
            CErr("NEO-CC NYI!",output);

          #ifdef CQ_HAS_TA
            auto cc = CQCCOptions(output, input, ss);
            cc->run(); 
          #else
            CErr("TiledArray must be compiled to use Coupled-Cluster code!");
          #endif
        }

        if ( elecJob == MR ) {

          if (doNEO)
            CErr("NEO-MCSCF NYI!",output);

          auto mcscf = CQMCSCFOptions(output,input,ss);
          mcscf->savFile = rstFile;
          mcscf->run(emPert);
        }

        firstStep = false;

      } // Loop over geometries
    } // Loop over different jobs

    
    ProgramTimer::tock("Chronus Quantum");
    printTimerSummary(std::cout);
     
    // Output CQ footer
    CQOutputFooter(output);

    // Reset std::cout and std::cerr
    if(outfile)  std::cout.rdbuf(coutbuf);
    if(rankfile) std::cerr.rdbuf(cerrbuf);

  }; // RunChronusQ


 }; // namespace ChronusQ
