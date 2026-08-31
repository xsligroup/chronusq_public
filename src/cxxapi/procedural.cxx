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
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>
#include <cxxapi/procedural.hpp>

#include <filesystem>
#include <algorithm>
#include <util/files.hpp>
#include <util/mpi.hpp>
#include <util/threads.hpp>
#include <util/timer.hpp>

#include <cubegen.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <basisset/remove_linear_dep_shells.hpp>
#include <integrals.hpp>
#include <singleslater.hpp>
#include <singleslater/multiparticless.hpp>
#include <coupledcluster.hpp>
#include <mcwavefunction.hpp>
#include <mcscf.hpp>
#include <newperturb.hpp>
#include <perturb.hpp>
#include <mp.hpp>

#include <findiff/geomgrad.hpp>
#include <particleintegrals/gradints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/gradints/incore.hpp>
#include <particleintegrals/gradints/direct.hpp>

#include <cqlinalg/blasext.hpp>

#include <geometrymodifier.hpp>
#include <geometrymodifier/moleculardynamics.hpp>
#include <geometrymodifier/singlepoint.hpp>
#include <physcon.hpp>

#include <corehbuilder/x2c.hpp>
#include <corehbuilder/nonrel.hpp>
#include <fockbuilder/matrixfock.hpp>

#include <fockbuilder/neofock.hpp>
#include <fockbuilder/interparticlefock.hpp>
#include <itersolver.hpp>

#include <unistd.h>
#include <limits.h>

#include <coupledcluster/TAManager.hpp>
#include <orbitalmodifiernew.hpp>
#include <gauxcutils.hpp>
#include <d3utils.hpp>
//#include <TiledArray/util/bug.h>

#include <intermediates.hpp>

#include <job.hpp>

namespace ChronusQ {

  template class NEOKohnShamBuilder<double,double>;
  template class NEOKohnShamBuilder<dcomplex,double>;
  template class NEOKohnShamBuilder<dcomplex,dcomplex>;

#ifdef ENABLE_BCAST_COUNTER
  int bcastCounter = 0;
#endif

  void RunChronusQ(std::string inFileName,
    std::string outFileName, std::string rstFileName,
    std::string scrFileName, bool rstExists) {

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
      std::string rankFileName = outFileName + ".mpi." + std::to_string(rank);
      rankfile = std::make_shared<std::ofstream>(rankFileName);
      std::cerr.rdbuf(rankfile->rdbuf());
      std::cerr << "Hello from RANK = " << rank << " / SIZE = " << size << std::endl;

#ifndef HOST_NAME_MAX // not defined on MacOS or other BSDs
#define HOST_NAME_MAX 1024
#endif
      char hostname[HOST_NAME_MAX];
      gethostname(hostname, HOST_NAME_MAX);
      std::cerr << "HostName = " << hostname << std::endl << std::endl;

      if (rank != 0) {
        std::cout.rdbuf(rankfile->rdbuf());
      }
      int i = 0;

#ifndef NDEBUG
//      while (i == 0) {
//        sleep(10);
//      }
#endif
//      TA::launch_lldb_xterm();
      MPI_Barrier(MPI_COMM_WORLD);
    }

    std::ostream &output = (rank == 0) ? std::cout : std::cerr;

    // Output CQ header
    CQOutputHeader(output);
    if(rankfile and rank == 0) CQOutputHeader(std::cerr);

    CQIntermediates::getInstance().clear();

    // Parse Input File
    CQInputFile input(inFileName);
    //SCFOptions scfOptions;
    TDSCFOptions tdSCFOptions;
    SingleSlaterGuessOptions ssGuessOptions;
    input.parse();

    // Misc options initializes CQMemManager
    CQMiscOptions(output,input);

    if (input.containsSection("SCF")) {
      ssGuessOptions.parseSection(input.getSection("SCF"));
    }

    // DeltaT and TMax will be overwritten if Dynamics Section is set
    if (input.containsSection("RT")) {
      tdSCFOptions.parseSection(input.getSection("RT"));
    }


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


      std::cout << "Parsed Input File:\n" << BannerTop << std::endl;
      std::cout << input << std::endl;
      std::cout << BannerEnd << "\n\n\n" << std::endl;
    }

    CQINPUT_VALID(output,input);

    CQPhysConSetOptions(output, input);

    // TEMPORARY
    bool doTemp = true;

    // Determine JOB type
    JobType jobType;
    
    try {
      jobType = parseJob(input.getData<std::string>("QM/JOB"));
    } catch (...) {
      CErr("Must Specify QM/JOB",output);
    }

    // Break into sequence of individual jobs
    std::vector<CQJob> jobs;
    if( jobType != JobType::SCF
        and not ((jobType == JobType::CC or jobType == JobType::EOMCC)
                  and ((input.containsData("CC/SKIPSCF") and input.getData<bool>("CC/SKIPSCF")) 
                    or (input.containsData("CC/SKIPCC") and input.getData<bool>("CC/SKIPCC"))))) {
      jobs.push_back(JobType::SCF);
    }
    // if RT MR propagation add MR calculation
    if( jobType == JobType::RT ) {
        if( input.containsSection("MCSCF")) {
            jobs.push_back(JobType::CI);
        }
    }
    jobs.push_back(jobType);

    // Check if we're doing NEO
    bool doNEO = false;
    if ( input.containsSection("SCF") ) {
      try {
        doNEO = input.getData<bool>("SCF/NEO");
      } catch(...) { ; }
    }

    Molecule mol(std::move(CQMoleculeOptions(output,input,scrFileName))); // Create Molecule object

    std::shared_ptr<BasisSet> basis = CQBasisSetOptions(output,input,mol,"BASIS"); // Create BasisSet object
    read_option_and_remove_linear_dependency(input, *basis, mol, output); // Remove linear dependent basis functions if requested

    std::shared_ptr<BasisSet> dfbasis = CQBasisSetOptions(output,input,mol,"DFBASIS"); // Create BasisSet object for DFBasis if defined
    std::shared_ptr<BasisSet> guessbasis = input.containsSection("GUESSBASIS") ? CQBasisSetOptions(output,input,mol,"GUESSBASIS") : nullptr; // Guess basis set for density projection
    IntegralOptions aoints_options = getIntegralOptions(output,input,basis,dfbasis,nullptr,"INTS");

    std::vector<QuantumSubsystem> quantumSubsystems;

    // Build electronic subsystem
    quantumSubsystems.push_back({"E","E","QM","BASIS","DFBASIS","GUESSBASIS","INTS",basis,dfbasis,guessbasis,aoints_options});

    // Build other quantum subsystems
    for(const auto& label : mol.getQuantumSystemLabels()) {
      auto subs = buildQuantumSubsystems(output,input,mol,label);
      quantumSubsystems.insert(quantumSubsystems.end(),
        std::make_move_iterator(subs.begin()),std::make_move_iterator(subs.end()));
    }
    const bool doRT = jobType == JobType::RT or jobType == JobType::EHRENFEST;
    if(doRT) CQRTExpandSubsystemAlgorithms(tdSCFOptions, quantumSubsystems);

    std::unordered_map<std::string, std::shared_ptr<BasisSet>> subsystemBasis;
    for(auto& sys : quantumSubsystems) subsystemBasis[sys.label] = sys.basis;
    
    // Build pair interactions between quantum subsystems
    std::vector<QuantumPairInteraction> quantumPairInteractions;
    for(size_t i = 0; i < quantumSubsystems.size(); ++i) {
      for(size_t j = i + 1; j < quantumSubsystems.size(); ++j) {
        quantumPairInteractions.push_back(buildQuantumPairInteraction(output,input,mol,quantumSubsystems[i],quantumSubsystems[j]));
      }
    }

    // Print Quantum Subsystem setup
    printQuantumSetup(output, input, quantumSubsystems, quantumPairInteractions);

    // Build Two-body integrals for all quantum subsystems and pair interactions
    buildTwoBodyIntegrals(output, mol, quantumSubsystems, quantumPairInteractions);

    // Electronic two-body integrals
    auto aoints = quantumSubsystems[0].integrals;

    std::shared_ptr<SingleSlaterBase> ss  = nullptr;

    SingleSlaterOptions ssOptions;

    // EM Perturbation for SCF
    EMPerturbation emPert;

    // SCF options
    SCFControls scfControls = CQSCFOptions(output,input,emPert);

    // Cubes for all quantum particles
    std::vector<std::shared_ptr<CubeGen>> cubes;
    // The first cube is always the electronic cube
    cubes.push_back(CQCUBEOptions(output,input,std::make_shared<Molecule>(mol),basis,emPert,-1.0));
    auto cube = cubes[0];

    // Create the SingleSlater object
    if (doNEO) {
      // Get Particle properties for each quantum subsystem
      for(auto& sys : quantumSubsystems)
        sys.ssOptions = CQSingleSlaterOptions(output,input,mol,*sys.basis);

      // Build MultiParticleSS object
      ss = CQMultiParticleSSOptions(output,input,mol,quantumSubsystems,quantumPairInteractions,scfControls);
      ss->scfControls = scfControls;
      ssOptions = quantumSubsystems.front().ssOptions;

      
      if( auto multiParticleSS = std::dynamic_pointer_cast<MultiParticleSSBase>(ss)) {

        // MO swapping applies to the electronic subsystem
        HandleOrbitalSwaps(output, input, *multiParticleSS->getSubSSBase("E"), "");

        ParseCubeSubsection(output, input, "SCF", ss->cubeOptsSS, cube);
        for(const auto& label : multiParticleSS->getLabels()) {
          if(label == "E") continue; // Electronic is already built 
          cubes.push_back(CQCUBEOptions(output, input, std::make_shared<Molecule>(mol),
                          subsystemBasis.at(label), emPert, multiParticleSS->getSubSSBase(label)->particle.charge));
        }
      }
    } else {
      if(not scfControls.subsystemGuesses.empty())
        CErr("Per-subsystem SCF guess controls require SCF/NEO = TRUE", output);
      ssOptions = CQSingleSlaterOptions(output,input,mol,*basis);
      ssOptions.scfControls = scfControls;
      ss = ssOptions.buildSingleSlater(output,mol,*basis, aoints, aoints_options);

      // MO swapping
      HandleOrbitalSwaps(output, input, *ss, "");

      ParseOrbitalPropSubsection(output, input, ss);
      ParseCubeSubsection(output, input, "SCF", ss->cubeOptsSS, cube);
    }

    // GAUXC
    std::shared_ptr<GauXCOptions> gauxcOptions;
    if (ssOptions.refOptions.isKSRef and ssOptions.intParam.useGauXC) {
      gauxcOptions = std::make_shared<GauXCOptions>(CQGauXCOptions(output, input, ssOptions,
        doNEO ? &quantumSubsystems : nullptr));
      if(doNEO)
        ss->gauxcUtils = gauxcOptions->buildGauXCUtils(quantumSubsystems, ss->molecule(), MPI_COMM_WORLD);
      else
        ss->gauxcUtils = gauxcOptions->buildGauXCUtils(basis, ss->molecule(), MPI_COMM_WORLD);
      ss->setupRangeSeparatedHybridExchange();
    }
    
    // Dispersion correction
#ifdef CQ_HAS_D3
    if(ssOptions.refOptions.useD3) {
      ss->d3Utils = std::make_shared<ChronusQ::D3Utils>(ssOptions.refOptions.d3RefString, ssOptions.refOptions.d3ModelString, true, false, true);
    }
#endif

    if(not doNEO) {
      if((ss->scfControls.guess == READMO or
          ss->scfControls.guess == READDEN) and not scrFileName.empty())
        ss->scrBinFileName = scrFileName;
      else if(ss->scfControls.guess == FCHKMO)
        ss->fchkFileName = scrFileName;
    } else {
      const bool usesBinGuess = std::any_of(
        quantumSubsystems.begin(), quantumSubsystems.end(), [](const auto& sys) {
          return sys.ssOptions.scfControls.guess == READMO or sys.ssOptions.scfControls.guess == READDEN;
        });
      const bool usesFchkGuess = std::any_of(
        quantumSubsystems.begin(), quantumSubsystems.end(), [](const auto& sys) {
          return sys.ssOptions.scfControls.guess == FCHKMO;
        });
      if(usesBinGuess and not scrFileName.empty())
        ss->scrBinFileName = scrFileName;
      else if(usesFchkGuess)
        ss->fchkFileName = scrFileName;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Create the restart and scratch files
    if( not rstExists and rank == 0 ) {
      SafeFile rstFile(rstFileName, rstExists);
      rstFile.createFile();
    }

    SafeFile rstFile(rstFileName, true);

    if( rank == 0 ) {

      // Save mol and basis data to bin
      mol.save(rstFile);
      basis->save(rstFile, mol);

      ss->savFile     = rstFile;
      // Attach savFile to every subsystem's SingleSlater and its integrals
      for(auto& sys : quantumSubsystems) {
        // the subsystem's own integrals object
        if(sys.integrals) sys.integrals->savFile = rstFile;
      }

      // Attach savFile to every pair interaction's integrals
      for(auto& pair : quantumPairInteractions) {
        if(pair.integrals) pair.integrals->savFile = rstFile;
      }
    }

    // Save reference info to bin file
    saveRefs( ssOptions, ss );

    // If doing NEO, propagate setup to subsystems
    if(auto multiParticleSS = std::dynamic_pointer_cast<MultiParticleSSBase>(ss)) {
      multiParticleSS->setSubSetup();
      multiParticleSS->saveSubsystemReferenceTypes();
    }

    // If we are doing RTCI we need a pointer to an mcscf object that is in this scope
    std::shared_ptr<MCWaveFunctionBase> mcwfn(nullptr);
    std::shared_ptr<TDEMPerturbation> tdPert = std::make_shared<TDEMPerturbation>();
    std::shared_ptr<RealTimeBase> rt;

    // Done setting up
    //
    // START OF REAL PROCEDURAL SECTION

    for( auto& job: jobs ) {

      bool firstStep = true;
//      if (ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::OFF) {
//        compute_X2C_CoreH_Fock( mol, *basis, aoints, emPert, ss, ssOptions);
//      }

      JobType elecJob = CQGeometryOptions(output, input, rstFile, job.jobType, mol, ss, mcwfn, rt, tdPert,
        quantumSubsystems, quantumPairInteractions, emPert, tdSCFOptions, gauxcOptions);

      // Loop over various structures
      while( mol.geometryModifier->hasNext() ) {

        // Update geometry
        mol.geometryModifier->electronicPotentialEnergy = ss->totalEnergy;
        mol.geometryModifier->update(true, mol, firstStep, tdSCFOptions, ss, emPert, cubes);
        // Update basis to the new geometry
        for(auto& sys : quantumSubsystems) {
          sys.basis->updateNuclearCoordinates(mol);
        }

        // Calculate integrals 
        // TODO: Time dependent field?
        if (elecJob == JobType::RT and !tdSCFOptions.doMD) {
          // For Real-time jobs, since basis functions are frozen, integrals do not need to be re-calculated.
          //                     assume we can re-use the same integrals from SCF job
          std::cout << "Skipping integral calculations for RT job. Assuming it's pre-computed." << std::endl;
        } else if ((elecJob == JobType::CC or elecJob == JobType::EOMCC)
                   and not (input.containsData("CC/SKIPSCF") and input.getData<bool>("CC/SKIPSCF"))) {
          std::cout << "Skipping integral calculations for CC. Assuming it's pre-computed." << std::endl;
        } else {
          // Symmetric (intra-subsystem) two-body integrals
          for(auto& sys : quantumSubsystems)
            sys.integrals->computeAOTwoE(*sys.basis, mol, emPert);
          // Asymmetric (cross/pair) two-body integrals
          for(auto& pair : quantumPairInteractions) {
            auto& basisA = subsystemBasis.at(pair.labelA);
            auto& basisB = subsystemBasis.at(pair.labelB);
            pair.integrals->computeAOTwoE(*basisA, *basisB, mol, emPert);
          }
          // Clean up raw ERI that might have been used during the asymmetric RI integral computation
          for(auto& sys : quantumSubsystems) {
            if(auto ints = std::dynamic_pointer_cast<Integrals<double>>(sys.integrals)) {
              if(auto ri = std::dynamic_pointer_cast<InCoreRITPI<double>>(ints->TPI))
                ri->clearRawERI();
            }
          }
        }

        // Rebuild the GauXC molecular grid with new geometry and new basis
        if ((job.jobType == JobType::BOMD or job.jobType == JobType::EHRENFEST) and (ss->gauxcUtils and gauxcOptions)) {
          if(doNEO)
            ss->gauxcUtils = gauxcOptions->buildGauXCUtils(quantumSubsystems, mol, MPI_COMM_WORLD);
          else
            ss->gauxcUtils = gauxcOptions->buildGauXCUtils(basis, mol, MPI_COMM_WORLD);
          ss->setupRangeSeparatedHybridExchange();
        }

        // Note, these guessSSOptions does not apply to NEO guess
        SingleSlaterOptions guessSSOptions(ssOptions);
        guessSSOptions.scfControls.guessBasis = guessbasis;
        guessSSOptions.scfControls.scfGuessOutFile = rstFileName;
        //guessSSOptions.refOptions.isKSRef = false;
        //guessSSOptions.refOptions.nC = 1;
        //guessSSOptions.hamiltonianOptions.OneEScalarRelativity = false;
        //guessSSOptions.hamiltonianOptions.OneESpinOrbit = false;

        // Run SCF job
        if( elecJob == JobType::SCF ) {

          if (ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::OFF) {
            if (doNEO) {
              auto multiBase = std::dynamic_pointer_cast<MultiParticleSSBase>(ss);
              if(!multiBase) CErr("Expected MultiParticleSS for NEO X2C", output);
              auto essbase = multiBase->getSubSSBase("E");
              compute_X2C_CoreH_Fock(mol, *basis, aoints, emPert, essbase, ssOptions);
            }  else
              compute_X2C_CoreH_Fock( mol, *basis, aoints, emPert, ss, ssOptions);
            if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::FOCK) {
              // If X2C.FOCK, the 2-component computation should be energy-only
              ss->scfControls.energyOnly = true;
              if (doNEO) {
                CErr("NEO with mmfX2C never tested!", output);
              }
            }
          }
          ss->formCoreH(emPert, true);
          //if(firstStep) ss->formGuess(guessSSOptions);

          auto conventionalSCF = buildConventionalSCF(ss->scfControls, *ss);

          if(conventionalSCF!=nullptr){
            ss->formGuess(emPert, guessSSOptions);
            ss->initializeSCF();
            conventionalSCF->run(emPert);
            if(cube) ss->runCube(cubes);
          }
        }

        if ( elecJob == JobType::MP2 ) {
          if( doNEO )
            CErr("NEO/multicomponent post-Hartree-Fock (MP2) is not yet implemented.");
          auto mp2 = CQMP2Options(output,input,ss);
          mp2->runMP2(emPert);
        }


        // Run RT job
        if( elecJob == JobType::RT ) {

          // Guard for 4C RT as this is untested code
          if(ssOptions.refOptions.nC == 4) CErr("Four Component Real Time NYI!");

          if (mcwfn){
            // Multi-slater (MC) RT still uses the legacy RealTimeCI machinery
            rt->intScheme.cubeOptsRTMS = mcwfn->cubeOptsMC;
            rt->intScheme.rtcubes = cubes;
            rt->run(firstStep, emPert);
          } else {
            // Single-slater RT is handled by the RealTimeSCF driver

            // Handle RT cube files
            if (cube) {
              tdSCFOptions.cubeOptsRT = ss->cubeOptsSS; 
              tdSCFOptions.rtcubes = cubes;
            }

            auto realtimeSCF = buildRealTimeSCF(tdSCFOptions, *tdPert, *ss);
            if(realtimeSCF) {
              realtimeSCF->initialize(0);
              realtimeSCF->run(emPert);
            }
          }
        }


        if( elecJob == JobType::RESP ) {

          std::cout << "Warning: Response Theory is experimental. Use at your own risk!" << std::endl;
          if( ss->scfControls.scfAlg == _SKIP_SCF and ss->scfControls.guess == READDEN )
            CErr("READDEN + SKIP + RESPONSE disabled. Use READMO instead.");

          auto resp = CQResponseOptions(output,input,ss,emPert);
          resp->savFile = rstFile;
          resp->run();

          if( MPIRank(MPI_COMM_WORLD) == 0 ) resp->printResults(output);
          MPI_Barrier(MPI_COMM_WORLD);

        }


        if( elecJob == JobType::CC or elecJob == JobType::EOMCC ){

          // FIXME: Need to implement NEO-CC
          if (doNEO)
            CErr("NEO-CC NYI!",output);

#ifdef CQ_HAS_TA

          if (input.containsData("CC/SKIPSCF") and input.getData<bool>("CC/SKIPSCF")) {
            // compute the necessary 1e ints
            std::vector<std::pair<OPERATOR,size_t>> ops{{LEN_ELECTRIC_MULTIPOLE,1}};
            aoints->computeAOOneP(mol, *basis, emPert, ops, ssOptions.hamiltonianOptions);
          }
          if (std::shared_ptr<SingleSlater<dcomplex,double>> ccref = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {
            runCoupledCluster(jobType, mol, ccref, aoints,  rstFile, input, output);
          } else if (std::shared_ptr<SingleSlater<double,double>> ccref = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {
            runCoupledCluster(jobType, mol, ccref, aoints,  rstFile, input, output);
          } else {
            CErr("Unsupported reference type for coupled cluster!", output);
          }
          TAManager::get().discard_cache();
          std::cout << TAManager::get() << std::endl;

#else
          CErr("TiledArray must be compiled to use Coupled-Cluster code!");
#endif
        }

        if ( elecJob == JobType::CI or elecJob == JobType::PT ) {

          EMPerturbation additionalPert; // in other places we might have an additional perturbation to mcscf 

          if (input.containsSection("MCSCF") and input.containsSection("CI")){
            CErr("Sections for the new and old CI codes are specified. Please specify either [MCSCF] or [CI]!");
          }

          if (input.containsSection("PERTURB") and input.containsSection("CI"))
            CErr("[PERTURB] section does not work with the new [CI] section, please use the [MCSCF] section.");
          if (input.containsSection("MRPT") and input.containsSection("MCSCF"))
            CErr("[MRPT] section does not work with the old [MCSCF] section, please use the [CI] section.");
          if (input.containsSection("PERTURB") and !input.containsSection("MCSCF"))
            CErr("Only [PERTURB] section was specified, please also specify the [MCSCF] section.");
          if (input.containsSection("MRPT") and !input.containsSection("CI"))
            CErr("Only [MRPT] section was specified, please also specify the [CI] section.");
            
          if (input.containsSection("MCSCF")) {
            std::shared_ptr<MCSCFBase> mcscf= CQMCSCFOptions(output,input,ss,mcwfn,emPert,cube,doNEO,quantumSubsystems,quantumPairInteractions);
            mcscf->savFile = rstFile;
            mcscf->run(additionalPert);
            ParseCubeSubsection(output,input,"MCSCF",mcwfn->cubeOptsMC,cube);
            if(cube) mcscf->runCube(cubes);
           
            if (input.containsSection("PERTURB")) {
              auto perturb = CQPerturbOptions(output,input,mcwfn);
              perturb->savFile = rstFile;
              perturb->run(emPert);  /// LXL: need to change the fields accordingly
            }
          }

          if (input.containsSection("CI")) {
            if (ssOptions.refOptions.nC == 1) {
              ss = ss->convert1CSSToGHFSS(ssOptions, output, input, emPert);
            } else {
              if (input.containsData("CI/NACTELECA") or 
                  input.containsData("CI/NACTELECB") or
                  input.containsData("CI/NACTORBA") or
                  input.containsData("CI/NACTORBB"))
                  CErr("Cannot specify CI/NActOA, CI/NActOB, CI/NActEA or CI/NActEB for unrestricted DAS.");
            }
            auto ci = CQCIOptions(output,input,ss,emPert,cube);
            ci->savFile = rstFile;
            ci->run(additionalPert);
            if(cube) ci->runCube(cubes);
            
            if (input.containsSection("MRPT")) {
              auto mrpt = CQMRPTSettings(output,input,ci);
              mrpt->savFile = rstFile;
              mrpt->run(additionalPert);  
            }
          }
        }

        firstStep = false;

      } // Loop over geometries
    } // Loop over different jobs


    CQMemManager::get().printHighWaterMark(output);

    ProgramTimer::tock("Chronus Quantum");
    printTimerSummary(std::cout);
     
    // Output CQ footer
    CQOutputFooter(output);

    CQIntermediates::getInstance().clear();

    // Reset std::cout and std::cerr
    if (rank == 0) {
      if (outfile) std::cout.rdbuf(coutbuf);
    } else {
      if (rankfile) std::cout.rdbuf(coutbuf);
    }
    if (rankfile) std::cerr.rdbuf(cerrbuf);

  }; // RunChronusQ

  void CQParser(std::string inFileName, std::string outFileName) {

  };

}; // namespace ChronusQ
