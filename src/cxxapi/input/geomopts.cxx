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
#include <physcon.hpp>
#include <chronusq_sys.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/procedural.hpp>
#include <cerr.hpp>

#include <geometrymodifier.hpp>
#include <geometrymodifier/moleculardynamics.hpp>
#include <geometrymodifier/singlepoint.hpp>

namespace ChronusQ {

  void CQDYNAMICS_VALID( std::ostream& out, CQInputFile& input ) {
    std::vector<std::string> allowedKeywords = {
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
      "SAVEALLGEOMETRY"
    };
  }

  JobType CQGeometryOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<MCWaveFunctionBase> mcscf,
    std::shared_ptr<RealTimeBase>& rt, 
    std::shared_ptr<TDEMPerturbation>& tdPert, std::shared_ptr<IntegralsBase> epints,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions)
  {

    JobType elecJob = job;
    if( job == JobType::BOMD or job == JobType::EHRENFEST or job == JobType::RT ) {
      elecJob = CQDynamicsOptions(out, input, rstFile, job, mol, ss, mcscf, rt, tdPert, epints, emPert, tdSCFOptions);
    }
    // add else if job == OPT
    else {
      // Single point job
      mol.geometryModifier = std::make_shared<SinglePoint>();
    }

    return elecJob;
  }

  // Function to get an IntegralsBase pointer out
  IntegralsBase* extractIntPtr(std::shared_ptr<SingleSlaterBase> ss) {
    IntegralsBase* ints = nullptr;

    if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) ) {
      ints = ss_t->aoints_.get();
    }
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) ) {
      ints = ss_t->aoints_.get();
    }
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) ) {
      ints = ss_t->aoints_.get();
    }

    return ints;
  }

  void createGradientIntegrals(CQInputFile& input, Molecule& mol,
    std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<IntegralsBase> epints)
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

#define ADD_GRAD_DIRECT(T) \
  std::vector<std::shared_ptr<DirectTPI<T>>> gints;\
  for ( auto i = 0; i < mol.atoms.size() * 3; i++ ) {\
    auto newg = ss2 ? \
      std::make_shared<DirectTPI<T>>( \
        ss1->basisSet(), ss2->basisSet(), mol, 1e-12) :  \
      std::make_shared<DirectTPI<T>>( \
        ss1->basisSet(), ss1->basisSet(), mol, 1e-12); \
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
      OPTOPT( GRAD_ALG = input.getData<std::string>(section +".GRADALG"););

      bool cmplx_ints = dynamic_cast<Integrals<dcomplex>*>(ints);


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
          ADD_GRAD_DIRECT(dcomplex);
        }
        else {
          ADD_GRAD_DIRECT(double);
        }
      }

    };

    // Begin actual work!

    // NEO
    if( auto neoss = std::dynamic_pointer_cast<NEOBase>(ss) ) {
      auto labels = neoss->getLabels();
      for( auto ilbl = 0; ilbl < labels.size(); ilbl++ ) {
        auto ss1 = neoss->getSubSSBase(labels[ilbl]);
        IntegralsBase* subints = extractIntPtr(ss1);
        // TODO: Generalize this
        std::string section = labels[ilbl] == "Protonic" ? "PINTS" : "INTS";
        createGradInt(subints, ss1, nullptr, section);

        auto intcast = dynamic_cast<Integrals<double>*>(subints);


        for(auto jlbl = ilbl + 1; jlbl < labels.size(); jlbl++) {
          auto ss2 = neoss->getSubSSBase(labels[jlbl]);
          createGradInt(epints.get(), ss2, ss1, "EPINTS");
          bool contractSecond = labels[ilbl] == "Protonic";

          if(auto neoss_t = std::dynamic_pointer_cast<NEOSS<double,double>>(neoss)) {
            auto epcasted = std::dynamic_pointer_cast<Integrals<double>>(epints);
            neoss_t->addGradientIntegrals(labels[ilbl], labels[jlbl], epcasted->gradERI, contractSecond);
          }
          else if(auto neoss_t = std::dynamic_pointer_cast<NEOSS<dcomplex,double>>(neoss)) {
            auto epcasted = std::dynamic_pointer_cast<Integrals<double>>(epints);
            neoss_t->addGradientIntegrals(labels[ilbl], labels[jlbl], epcasted->gradERI, contractSecond);
          }
          else if(auto neoss_t = std::dynamic_pointer_cast<NEOSS<dcomplex,dcomplex>>(neoss)) {
            auto epcasted = std::dynamic_pointer_cast<Integrals<dcomplex>>(epints);
            neoss_t->addGradientIntegrals(labels[ilbl], labels[jlbl], epcasted->gradERI, contractSecond);
          }
          else {
            CErr("No successful NEOSS cast. This should never happen!");
          }
        }
      }
    }
    // Conventional
    else {
      IntegralsBase* eints = extractIntPtr(ss);
      createGradInt(eints, ss, nullptr, "INTS");
    }

  }

  JobType CQDynamicsOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<MCWaveFunctionBase> mcscf,
    std::shared_ptr<RealTimeBase>& rt, 
    std::shared_ptr<TDEMPerturbation>& tdPert, std::shared_ptr<IntegralsBase> epints,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions)
  {

    JobType elecJob;


    if( job == JobType::BOMD or job == JobType::EHRENFEST ) {

      if( not input.containsSection("DYNAMICS") )
        CErr("Dynamics Section must be specified for BOMD/EHRENFEST/RT job",out);

      double tMax, deltaT;
      try {
        tMax = input.getData<double>("DYNAMICS.TMAX");
      } catch(...) {
        CErr("Must specify DYNAMICS.TMAX for simulation length");
      }

      try {
        deltaT = input.getData<double>("DYNAMICS.DELTAT");
      } catch(...) {
        CErr("Must specify DYNAMICS.DELTAT for integration time step");
      }

      // Create geometry updater
      MDOptions mdOpt(tMax, deltaT);

      OPTOPT( mdOpt.nMidpointFockSteps = input.getData<size_t>("DYNAMICS.NNUCPGRAD"); )
      OPTOPT( mdOpt.nElectronicSteps = input.getData<size_t>("DYNAMICS.NELECPNUC"); )

      OPTOPT( mdOpt.saveAllGeometry = input.getData<bool>("DYNAMICS.SAVEALLGEOMETRY");)

      // Parsing restart options
      std::string restart = "FALSE";
      OPTOPT( restart = input.getData<std::string>("DYNAMICS.RESTART");)
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
          CErr("Invalid input for DYNAMICS.RESTART");
        }
      }



      // TODO: we need to have a separate GUESS section for MD
      if( job == JobType::BOMD )
        mdOpt.nMidpointFockSteps = 0;

      auto md = std::make_shared<MolecularDynamics>(mdOpt, mol, rstFile);
      
      // If doNEO, Choose how to move quantum proton basis function centers during dynamics simulations
      // Default is 'fixed'
      if(mol.atomsQ.size() > 0) {
        bool useTPB = false;
        OPTOPT( useTPB = input.getData<bool>("DYNAMICS.TPB");)
        if (useTPB) {
          md->NEODynamicsOpts.tpb = true;
          md->NEODynamicsOpts.includeQProtKE = (job == JobType::BOMD) ? true : false ;
        } else {
          std::cout << "Quantum Proton will be fixed during dynamics" << std::endl;
        }
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
      std::cout<< "================================================================================" << std::endl;
      std::cout << std::endl;

      createGradientIntegrals(input, mol, ss, epints);



      // Provide definition for gradient calculations
      md->gradientGetter = [&, ss](){ return ss->getGrad(emPert,false,false); };

      // Provide definitions for 
      //     - std::function<void()> updateBasisIntsHamiltonian; (for Ehrenfest and BOMD)
      //     - std::function<double()> finalMidpointFock; (only for Ehrenfest)
      // For NEOSS and regular singleslater, the definitions are different
      if( auto neoss = std::dynamic_pointer_cast<NEOBase>(ss) ) {
        // Obtain NEO integrals and basis as a vector
        std::vector<IntegralsBase*> ints;
        std::vector<BasisSet*> bases;
        BasisSet* ebasis = nullptr;
        BasisSet* pbasis = nullptr;
        IntegralsBase* pint = nullptr;
        auto labels = neoss->getLabels();
        for( auto label: labels ) {
          auto subss = neoss->getSubSSBase(label);
          ints.push_back(extractIntPtr(subss));
          bases.push_back(&subss->basisSet());
          if( label == "Electronic" ) {
            ebasis = bases.back();
          } else if( label == "Protonic" ) {
            pbasis = bases.back();
            pint = ints.back();
          }
        }

        if(md->NEODynamicsOpts.tpb)
          pint->options_.includeTau = true;


        md->updateBasisIntsHamiltonian = [=, &mol, &emPert](){
          // Update basis and two-e integrals at new geometry
          for( auto isub = 0; isub < ints.size(); isub++ ) {
            bases[isub]->updateNuclearCoordinates(mol);
            ints[isub]->computeAOTwoE(*bases[isub], mol, emPert);
          }
          epints->computeAOTwoE(*ebasis, *pbasis, mol, emPert);
          // Update 1-e integrals, including S metric and transformation matrix
          ss->formCoreH(emPert,false);
          ss->formFock(emPert,false);
        };

        if (job == JobType::EHRENFEST) {
          md->finalMidpointFock = [=, &mol, &emPert](){
            // Update basis, integrals, and hamiltonian
            md->updateBasisIntsHamiltonian();

            // Transform ortho density with new metric for property and gradient evaluation
            if( auto ss_t = std::dynamic_pointer_cast<NEOSS<double,double>>(ss) )           ss_t->ortho2aoDen();
            else if( auto ss_t = std::dynamic_pointer_cast<NEOSS<dcomplex,double>>(ss) )    ss_t->ortho2aoDen();
            else if( auto ss_t = std::dynamic_pointer_cast<NEOSS<dcomplex,dcomplex>>(ss) )  ss_t->ortho2aoDen();
            else CErr("Unsuccessful Cast!");

            // Recompute fock matrix and get updated energy
            ss->formFock(emPert,false);
            ss->computeEnergy(emPert);
            return ss->totalEnergy;
          };
        }

        md->pertFirstAtom = [=, &mol, &emPert](){
          // Apply perturbation for first atom
          mol.atoms[0].coord[0] += md->mdOptions.pert_val_x;
          mol.atoms[0].coord[1] += md->mdOptions.pert_val_y;
          mol.atoms[0].coord[2] += md->mdOptions.pert_val_z;
          mol.update();
        };

      } // End definitions for updateBasisIntsHamiltonian and finalMidpointFock for when ss is NEOSS
      else {
        auto aoints = extractIntPtr(ss);
        BasisSet* basis = &ss->basisSet();

        md->updateBasisIntsHamiltonian = [=, &mol, &emPert]() {
          // Update basis and two-e integrals at new geometry
          basis->updateNuclearCoordinates(mol);
          aoints->computeAOTwoE(*basis, mol, emPert);
          // Update 1-e integrals, including S metric and transformation matrix
          ss->formCoreH(emPert, false);
          ss->formFock(emPert,false);
        };

        if (job == JobType::EHRENFEST) {
          md->finalMidpointFock = [=, &mol, &emPert](){
            // Update basis, integrals, and hamiltonian
            md->updateBasisIntsHamiltonian();

            // Transform ortho density with new metric for property and gradient evaluation
            if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) )           ss_t->ortho2aoDen();
            else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) )    ss_t->ortho2aoDen();
            else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) )  ss_t->ortho2aoDen();
            else CErr("Unsuccessful Cast!");

            // Recompute fock matrix and get updated energy
            ss->formFock(emPert,false);
            ss->computeEnergy(emPert);
            return ss->totalEnergy;
          };
        }

        md->pertFirstAtom = [&, aoints, basis, ss](){
          // Apply perturbation for first atom
          mol.atoms[0].coord[0] += md->mdOptions.pert_val_x;
          mol.atoms[0].coord[1] += md->mdOptions.pert_val_y;
          mol.atoms[0].coord[2] += md->mdOptions.pert_val_z;
          mol.update();
        };

      }  // End definitions for updateBasisIntsHamiltonian and finalMidpointFock for when ss is regular SingleSlater

      // Parse initial velocity

      std::string velocityStr;
      OPTOPT( velocityStr = input.getData<std::string>("DYNAMICS.VELOCITY");)
      if ( not velocityStr.empty() ) {
        if (mdOpt.restoreFromNuclearStep != 0 )
          CErr("Restart with a newly specified velocity NYI!");
        md->parseVelocityFromInput(mol, velocityStr, out);
      }


      // Whether to perturb the first atom's geometry
      OPTOPT( md->mdOptions.pertFirstAtom = input.getData<bool>("DYNAMICS.INIT_PERT");)
      OPTOPT( md->mdOptions.pert_val_x = input.getData<double>("DYNAMICS.PERT_VALUE_X");)
      OPTOPT( md->mdOptions.pert_val_y = input.getData<double>("DYNAMICS.PERT_VALUE_Y");)
      OPTOPT( md->mdOptions.pert_val_z = input.getData<double>("DYNAMICS.PERT_VALUE_Z");)
      

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
      mol.geometryModifier = std::make_shared<SinglePoint>();
      elecJob = JobType::RT;
      // Handle field specification
      try {
        // Get raw string from input
        std::string fieldSpec = input.getData<std::string>("RT.FIELD");
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
      if(mcscf){
        rt = CQRealTimeOptions(out,input,ss,mcscf,tdPert, emPert);
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
