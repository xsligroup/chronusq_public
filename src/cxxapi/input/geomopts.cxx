/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <physcon.hpp>
#include <chronusq_sys.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/procedural.hpp>
#include <cerr.hpp>

#include <geometrymodifier.hpp>
#include <geometrymodifier/moleculardynamics.hpp>
#include <geometrymodifier/singlepoint.hpp>
#include <singleslater/multiparticless.hpp>

namespace ChronusQ {

  std::set<std::string> CQDYNAMICS_VALID(const std::map<std::string, std::string>& inputSection) {
    std::set<std::string> allowedKeywords = {
      "NNUCPGRAD",
      "NELECPNUC",
      "TMAX",
      "DELTAT",
      "TPB",
      "RESTART",
      "INIT_PERT",
      "PERT_VALUE_X",
      "PERT_VALUE_Y",
      "PERT_VALUE_Z",
      "SAVEALLGEOMETRY",
      "PRINTPROPERTY",
      "PROJECT_ORTHO_DEN",
      "ONLYMOVEH",
      "VELOCITY"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  JobType CQGeometryOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<MCWaveFunctionBase> mcscf,
    std::shared_ptr<RealTimeBase>& rt,
    std::shared_ptr<TDEMPerturbation>& tdPert,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions,
    const std::shared_ptr<GauXCOptions>& gauxcOptions)
  {

    JobType elecJob = job;
    if( job == JobType::BOMD or job == JobType::EHRENFEST or job == JobType::RT ) {
      elecJob = CQDynamicsOptions(out, input, rstFile, job, mol, ss, mcscf, rt, tdPert,
        quantumSubsystems, quantumPairInteractions, emPert, tdSCFOptions, gauxcOptions);
    }
    // add else if job == OPT
    else {
      // Single point job
      mol.geometryModifier = std::make_shared<SinglePoint>();
    }

    return elecJob;
  }

  void createGradientIntegrals(CQInputFile& input, Molecule& mol,
    std::shared_ptr<SingleSlaterBase> ss,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions)
  {

#define ADD_GRAD_INCORE(T) \
  std::vector<std::shared_ptr<InCore4indexTPI<T>>> gints;\
  for ( auto i = 0; i < mol.atoms.size() * 3; i++ ) {\
    auto newg = ss2 ? \
      std::make_shared<InCore4indexTPI<T>>( \
        ss1->basisSet().nBasis, \
        ss2->basisSet().nBasis) : \
      std::make_shared<InCore4indexTPI<T>>( \
        ss1->basisSet().nBasis); \
    gints.push_back(newg); \
  } \
    \
  auto casted = dynamic_cast<Integrals<T>*>(ints); \
  casted->gradERI = std::make_shared<GradInts<TwoPInts,T>>( \
    ss1->basisSet().nBasis, mol.atoms.size(), gints \
  );

#define ADD_GRAD_DIRECT(T, THRESH) \
  std::vector<std::shared_ptr<DirectTPI<T>>> gints;\
  for ( auto i = 0; i < mol.atoms.size() * 3; i++ ) {\
    auto newg = ss2 ? \
      std::make_shared<DirectTPI<T>>( \
        ss1->basisSet(), ss2->basisSet(), mol, THRESH) :  \
      std::make_shared<DirectTPI<T>>( \
        ss1->basisSet(), ss1->basisSet(), mol, THRESH); \
    gints.push_back(newg); \
  } \
    \
  auto casted = dynamic_cast<Integrals<T>*>(ints); \
  casted->gradERI = std::make_shared<GradInts<TwoPInts,T>>( \
    ss1->basisSet().nBasis, mol.atoms.size(), gints \
  );


    // Function to create a new gradient integral TPI and add it to ints
    auto createGradInt = [&](IntegralsBase* ints,
      std::shared_ptr<SingleSlaterBase> ss1,
      std::shared_ptr<SingleSlaterBase> ss2,
      std::string section) {

      std::string GRAD_ALG = "DIRECT";
      OPTOPT( GRAD_ALG = input.getData<std::string>(section +"/GRADALG"););
      double GRAD_SCHWARZ = 1e-12;
      OPTOPT( GRAD_SCHWARZ = input.getData<double>(section +"/GRADSCHWARZ"););

      bool cmplx_ints = dynamic_cast<Integrals<dcomplex>*>(ints);
      if ( cmplx_ints && GRAD_ALG == "DIRECT") {
        std::cout << "DIRECT is not available for GIAO gradients. Changing to INCORE" << std::endl;
#ifdef CQ_ENABLE_MPI
        if( MPISize(ss1->comm) > 1 ) {
          CErr("DIRECT GIAO gradients are not available for MPI",std::cout);
        }
#endif
        GRAD_ALG = "INCORE";
      }

      if( GRAD_ALG == "INCORE" ) {

        if( cmplx_ints ) {
          ADD_GRAD_INCORE(dcomplex);
        }
        else {
          ADD_GRAD_INCORE(double);
        }
      }
      else {
        if( cmplx_ints ) {
          ADD_GRAD_DIRECT(dcomplex, GRAD_SCHWARZ);
        }
        else {
          ADD_GRAD_DIRECT(double, GRAD_SCHWARZ);
        }
      }

    };

    auto multiBase = std::dynamic_pointer_cast<MultiParticleSSBase>(ss);
    std::unordered_map<std::string,std::shared_ptr<SingleSlaterBase>> subSS;

    for(auto& sys : quantumSubsystems) {
      auto subsystem = multiBase ? multiBase->getSubSSBase(sys.label) : ss;
      if(!subsystem or !sys.integrals)
        CErr("Missing subsystem data while building gradient integrals for " + sys.label);
      subSS.emplace(sys.label, subsystem);
      createGradInt(sys.integrals.get(), subsystem, nullptr, sys.intsSection);
    }

    for(auto& pair : quantumPairInteractions) {
      if(!multiBase or !pair.integrals)
        CErr("Missing multiparticle data while building gradient integrals for " +
          pair.labelA + "-" + pair.labelB);

      auto ssA = subSS.at(pair.labelA);
      auto ssB = subSS.at(pair.labelB);
      createGradInt(pair.integrals.get(), ssA, ssB, pair.intsSection);

      if(auto multiSS = std::dynamic_pointer_cast<MultiParticleSS<double,double>>(ss)) {
        auto ints = std::dynamic_pointer_cast<Integrals<double>>(pair.integrals);
        if(!ints) CErr("Invalid gradient integral type for " + pair.labelA + "-" + pair.labelB);
        bool contractSecond = multiSS->getCrossTPIs(pair.labelA,pair.labelB).first;
        multiSS->addGradientIntegrals(pair.labelA,pair.labelB,ints->gradERI,contractSecond);
      }
      else if(auto multiSS = std::dynamic_pointer_cast<MultiParticleSS<dcomplex,double>>(ss)) {
        auto ints = std::dynamic_pointer_cast<Integrals<double>>(pair.integrals);
        if(!ints) CErr("Invalid gradient integral type for " + pair.labelA + "-" + pair.labelB);
        bool contractSecond = multiSS->getCrossTPIs(pair.labelA,pair.labelB).first;
        multiSS->addGradientIntegrals(pair.labelA,pair.labelB,ints->gradERI,contractSecond);
      }
      else if(auto multiSS = std::dynamic_pointer_cast<MultiParticleSS<dcomplex,dcomplex>>(ss)) {
        auto ints = std::dynamic_pointer_cast<Integrals<dcomplex>>(pair.integrals);
        if(!ints) CErr("Invalid gradient integral type for " + pair.labelA + "-" + pair.labelB);
        bool contractSecond = multiSS->getCrossTPIs(pair.labelA,pair.labelB).first;
        multiSS->addGradientIntegrals(pair.labelA,pair.labelB,ints->gradERI,contractSecond);
      }
      else {
        CErr("Pair gradient integrals require a MultiParticleSS object");
      }
    }

  }

  JobType CQDynamicsOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<MCWaveFunctionBase> mcwfn,
    std::shared_ptr<RealTimeBase>& rt,
    std::shared_ptr<TDEMPerturbation>& tdPert,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions,
    const std::shared_ptr<GauXCOptions>& gauxcOptions)
  {

    JobType elecJob;


    if( job == JobType::BOMD or job == JobType::EHRENFEST ) {

      if( not input.containsSection("DYNAMICS") )
        CErr("Dynamics Section must be specified for BOMD/EHRENFEST/RT job",out);

      double tMax, deltaT;
      try {
        tMax = input.getData<double>("DYNAMICS/TMAX");
      } catch(...) {
        CErr("Must specify DYNAMICS/TMAX for simulation length");
      }

      try {
        deltaT = input.getData<double>("DYNAMICS/DELTAT");
      } catch(...) {
        CErr("Must specify DYNAMICS/DELTAT for integration time step");
      }

      // Create geometry updater
      MDOptions mdOpt(tMax, deltaT);

      OPTOPT( mdOpt.nMidpointFockSteps = input.getData<size_t>("DYNAMICS/NNUCPGRAD"); )
      OPTOPT( mdOpt.nElectronicSteps = input.getData<size_t>("DYNAMICS/NELECPNUC"); )

      OPTOPT( mdOpt.saveAllGeometry = input.getData<bool>("DYNAMICS/SAVEALLGEOMETRY");)
      OPTOPT( mdOpt.printProperty = input.getData<bool>("DYNAMICS/PRINTPROPERTY");)
      OPTOPT( mdOpt.onlyMoveH = input.getData<bool>("DYNAMICS/ONLYMOVEH");)

      // Parsing restart options
      std::string restart = "FALSE";
      OPTOPT( restart = input.getData<std::string>("DYNAMICS/RESTART");)
      trim(restart);
      if (not restart.compare("TRUE")) {
        mdOpt.restoreFromNuclearStep = -1;
        std::cout << "Restart Option Found!" << std::endl;
        std::cout << "Restart MD from the last saved point" << std::endl;
      } else if (not restart.compare("FALSE")) {
        mdOpt.restoreFromNuclearStep = 0; // Default value
      } else {
        std::istringstream iss(restart);
        double inputTime;
        if (iss >> inputTime && iss.eof()) { // Checks for valid double and consumes entire input
          if (inputTime == -1) {
            // Special case for explicit -1 input
            mdOpt.restoreFromNuclearStep = -1;
          } else if (inputTime >= 0) {
            // Valid positive double handling
            mdOpt.restoreFromNuclearStep = static_cast<long int>(inputTime / deltaT);
            double restartTime = mdOpt.restoreFromNuclearStep * deltaT;
            if (restartTime > 0 && restartTime < tMax) {
              std::cout << "Restart MD from time=" << restartTime << " AU" << std::endl;
            } else if (restartTime >= 0 && restartTime < deltaT){
              CErr("Invalid Restart Step! Restart time needs to be larger than deltaT.");
            } else {
              CErr("Invalid Restart Step! Time out of bounds (Need to be 0~TMax).");
            }
          } else {
            CErr("Input must be a non-negative double or -1.");
          }
        } else {
          CErr("Invalid input for DYNAMICS/RESTART");
        }
      }



      // TODO: we need to have a separate GUESS section for MD
      if( job == JobType::BOMD )
        mdOpt.nMidpointFockSteps = 0;

      auto md = std::make_shared<MolecularDynamics>(mdOpt, mol, rstFile, MPI_COMM_WORLD);
      
      // If doNEO, Choose how to move quantum proton basis function centers during dynamics simulations
      // Default is 'fixed'
      if(mol.atomsQ.size() > 0) {
        bool useTPB = false;
        OPTOPT( useTPB = input.getData<bool>("DYNAMICS/TPB");)
        if (useTPB) {
          md->NEODynamicsOpts.tpb = true;
          md->NEODynamicsOpts.includeQProtKE = (job == JobType::BOMD) ? true : false ;
        } else {
          std::cout << "Quantum Proton will be fixed during dynamics" << std::endl;
        }
        if (md->mdOptions.onlyMoveH) CErr("Only Move Hydrogen is not supported for NEO calculations");
      }

      mol.geometryModifier = md;

      // TODO: operator overload << such that we can print out information of md to output
      // Temporary sketchy workaround:
      std::cout << std::endl;
      std::cout<< "================================================================================" << std::endl;
      std::cout << "Molecular Dynamics Information " << std::endl;
      std::cout << "TMAX:                    " << tMax << std::endl;
      std::cout << "DeltaT:                  " << deltaT << std::endl;
      std::cout << "JobType:                 " <<  (job==JobType::BOMD? "BOMD" : "Ehrenfest") << std::endl;
      std::cout << "DoNEO:                   " << (mol.atomsQ.size()>0? "True" : "False") << std::endl;
      if(mol.atomsQ.size()>0) {
        std::cout << "Traveling Proton Basis:  " << (md->NEODynamicsOpts.tpb? "True" : "False") << std::endl;
        std::cout << "QProt KE Included:       " << (md->NEODynamicsOpts.includeQProtKE? "True" : "False") << std::endl;
      }
      std::cout << "Only Move Hydrogen:     " << (md->mdOptions.onlyMoveH? "True" : "False") << std::endl;
      std::cout<< "================================================================================" << std::endl;
      std::cout << std::endl;

      createGradientIntegrals(input, mol, ss, quantumSubsystems, quantumPairInteractions);

      // Provide definition for gradient calculations
      md->gradientGetter = [&, ss](){ return ss->getGrad(emPert,false,false); };

      std::unordered_map<std::string,std::shared_ptr<BasisSet>> subsystemBasis;
      for(auto& sys : quantumSubsystems) {
        subsystemBasis.emplace(sys.label,sys.basis);
        if(md->NEODynamicsOpts.tpb and sys.label != "E")
          sys.integrals->options_.includeTau = true;
      }

      md->updateBasisIntsHamiltonian = [&, ss, subsystemBasis, gauxcOptions](){
        // Update every subsystem basis and its intra-particle integrals.
        for(auto& sys : quantumSubsystems) {
          sys.basis->updateNuclearCoordinates(mol);
          sys.integrals->computeAOTwoE(*sys.basis,mol,emPert);
        }

        // Update every pair interaction with the same subsystem ordering used at setup.
        for(auto& pair : quantumPairInteractions)
          pair.integrals->computeAOTwoE(*subsystemBasis.at(pair.labelA),
            *subsystemBasis.at(pair.labelB),mol,emPert);

        // Rebuild the GauXC molecular grid with new geometry and new basis
        if(ss->gauxcUtils and gauxcOptions) {
          if(quantumSubsystems.size() > 1)
            ss->gauxcUtils = gauxcOptions->buildGauXCUtils(quantumSubsystems, mol, MPI_COMM_WORLD);
          else
            ss->gauxcUtils = gauxcOptions->buildGauXCUtils(quantumSubsystems.front().basis, mol, MPI_COMM_WORLD);
          ss->setupRangeSeparatedHybridExchange();
        }

        // Update one-electron integrals, metric transformations and Fock matrices.
        ss->formCoreH(emPert,false);
        ss->formFock(emPert,false);
      };

      if(job == JobType::EHRENFEST) {
        md->finalMidpointFock = [&, ss, md](){

          // Update basis, integrals, and hamiltonian
          md->updateBasisIntsHamiltonian();

          // Transform the propagated density with the metric at the new geometry.
          ss->ortho2aoDen();

          // Recompute the Fock matrix and energy from the transformed density.
          ss->formFock(emPert,false);
          ss->computeEnergy(emPert);
#ifdef CQ_HAS_D3
          if(ss->d3Utils) {
            ss->d3Utils->evaluate(ss->molecule());
            ss->totalEnergy += ss->d3Utils->result().energy;
          }
#endif
          return ss->totalEnergy;
        };
      }

      md->pertFirstAtom = [&mol, md](){
        // Apply perturbation for first atom
        mol.atoms[0].coord[0] += md->mdOptions.pert_val_x;
        mol.atoms[0].coord[1] += md->mdOptions.pert_val_y;
        mol.atoms[0].coord[2] += md->mdOptions.pert_val_z;
        mol.update();
      };

      // Parse initial velocity

      std::string velocityStr;
      OPTOPT( velocityStr = input.getData<std::string>("DYNAMICS/VELOCITY");)
      if ( not velocityStr.empty() ) {
        if (mdOpt.restoreFromNuclearStep != 0 )
          CErr("Restart with a newly specified velocity NYI!");
        md->parseVelocityFromInput(mol, velocityStr, out);
      }


      // Whether to perturb the first atom's geometry
      OPTOPT( md->mdOptions.pertFirstAtom = input.getData<bool>("DYNAMICS/INIT_PERT");)
      OPTOPT( md->mdOptions.pert_val_x = input.getData<double>("DYNAMICS/PERT_VALUE_X");)
      OPTOPT( md->mdOptions.pert_val_y = input.getData<double>("DYNAMICS/PERT_VALUE_Y");)
      OPTOPT( md->mdOptions.pert_val_z = input.getData<double>("DYNAMICS/PERT_VALUE_Z");)
      
      OPTOPT( md->mdOptions.projectOrthoDen = input.getData<bool>("DYNAMICS/PROJECT_ORTHO_DEN");)
      

      // Set up electronic jobs for each MD type
      if( job == JobType::BOMD ) {
        elecJob = JobType::SCF;
      } else if( job == JobType::EHRENFEST ) {

        // Determint deltaT in RT by # of Midpoint and RT steps specified in dynamics section (settings in RT section is disabled)
        // TODO: Error out when both RT and Dynamics Section have conflicting TMax and DeltaT for RT job
        tdSCFOptions.deltaT =  mdOpt.timeStepAU/(mdOpt.nMidpointFockSteps*mdOpt.nElectronicSteps);
        tdSCFOptions.totalMDSteps = mdOpt.nMidpointFockSteps*mdOpt.nNuclearSteps;
        tdSCFOptions.rtMaxStepsPerMDStep = mdOpt.nElectronicSteps;
        tdSCFOptions.doMD = true;
        tdSCFOptions.includeTau = md->NEODynamicsOpts.tpb;

        elecJob = JobType::RT;
      }

    }
    else if( job == JobType::RT ) {
      // Single point job
      mol.geometryModifier = std::make_shared<SinglePoint>(MPI_COMM_WORLD);
      elecJob = JobType::RT;
      // Handle field specification
      try {
        // Get raw string from input
        std::string fieldSpec = input.getData<std::string>("RT/FIELD");
        std::istringstream fieldStream(fieldSpec);
        // Loop over field specification lines
        for(std::string fieldStr; std::getline(fieldStream, fieldStr); ) {
        //parseRTField(fieldStr, out, tdPert);
        //tdPert->addField(parseRTField(fieldStr, out));
        auto parsedfield = parseRTField(fieldStr, out);
        if (parsedfield)
          tdPert->addField(parsedfield);
        }
      } catch( std::runtime_error &e ) {
        throw;
      } catch(...) { 
        out << "  *** No TD Field Defaulting to Trivial Propagation ***\n";
      }
      if(mcwfn){
        rt = CQRealTimeOptions(out,input,ss,mcwfn,tdPert, emPert);
        rt->setTDPerturbation(*tdPert);
        rt->savFile = ss->savFile;
        rt->createRTDataSets(0);
      }

    }

    return elecJob;
  }

#undef ADD_GRAD_INCORE
#undef ADD_GRAD_DIRECT

}
