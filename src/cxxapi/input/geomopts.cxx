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
    };
  }

  JobType CQGeometryOptions(std::ostream& out, CQInputFile& input, 
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss,
    std::shared_ptr<RealTimeBase>& rt, std::shared_ptr<IntegralsBase> epints,
    EMPerturbation& emPert)
  {

    JobType elecJob = job;
    if( job == JobType::BOMD or job == JobType::EHRENFEST or job == JobType::RT ) {
      elecJob = CQDynamicsOptions(out, input, job, mol, ss, rt, epints, emPert);
    }
    // add else if job == OPT
    else {
      // Single point job
      MolecularOptions molOpt(0.0, 0.0);
      mol.geometryModifier = std::make_shared<SinglePoint>(molOpt);
    }

    return elecJob;
  }

  // Function to get an IntegralsBase pointer out
  IntegralsBase* extractIntPtr(std::shared_ptr<SingleSlaterBase> ss) {
    IntegralsBase* ints = nullptr;

    if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) ) {
      ints = &ss_t->aoints;
    }
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) ) {
      ints = &ss_t->aoints;
    }
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) ) {
      ints = &ss_t->aoints;
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
        ss1->memManager, ss1->basisSet().nBasis, \
        ss2->basisSet().nBasis) : \
      std::make_shared<InCore4indexTPI<T>>( \
        ss1->memManager, ss1->basisSet().nBasis); \
    gints.push_back(newg); \
  } \
    \
  auto casted = dynamic_cast<Integrals<T>*>(ints); \
  casted->gradERI = std::make_shared<GradInts<TwoPInts,T>>( \
    ss1->memManager, ss1->basisSet().nBasis, mol.atoms.size(), gints \
  );

#define ADD_GRAD_DIRECT(T) \
  std::vector<std::shared_ptr<DirectTPI<T>>> gints;\
  for ( auto i = 0; i < mol.atoms.size() * 3; i++ ) {\
    auto newg = ss2 ? \
      std::make_shared<DirectTPI<T>>( \
        ss1->memManager, ss1->basisSet(), ss2->basisSet(), mol, 1e-12) :  \
      std::make_shared<DirectTPI<T>>( \
        ss1->memManager, ss1->basisSet(), ss1->basisSet(), mol, 1e-12); \
    gints.push_back(newg); \
  } \
    \
  auto casted = dynamic_cast<Integrals<T>*>(ints); \
  casted->gradERI = std::make_shared<GradInts<TwoPInts,T>>( \
    ss1->memManager, ss1->basisSet().nBasis, mol.atoms.size(), gints \
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
          createGradInt(epints.get(), ss1, ss2, "EPINTS");
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

  JobType CQDynamicsOptions(std::ostream& out, CQInputFile& input, 
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss,
    std::shared_ptr<RealTimeBase>& rt, std::shared_ptr<IntegralsBase> epints,
    EMPerturbation& emPert)
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
      MolecularOptions molOpt(tMax, deltaT);

      OPTOPT( molOpt.nMidpointFockSteps = input.getData<size_t>("DYNAMICS.NNUCPGRAD"); )
      OPTOPT( molOpt.nElectronicSteps = input.getData<size_t>("DYNAMICS.NELECPNUC"); )

      // TODO: we need to have a separate GUESS section for MD
      if( job == JobType::BOMD )
        molOpt.nMidpointFockSteps = 0;

      auto md = std::make_shared<MolecularDynamics>(molOpt, mol);
      mol.geometryModifier = md;

      // Handle gradient integrals
      createGradientIntegrals(input, mol, ss, epints);

      // Set gradient computation methods
      if( job == JobType::BOMD ) {
        md->gradientGetter = [&, ss](){ return ss->getGrad(emPert,false,false); };
        elecJob = JobType::SCF;
      }
      else if( job == JobType::EHRENFEST ) {

        rt = CQRealTimeOptions(out,input,ss,emPert);
        rt->savFile = ss->savFile;
        rt->intScheme.deltaT = molOpt.timeStepAU/
                               (molOpt.nMidpointFockSteps*molOpt.nElectronicSteps);
        rt->createRTDataSets(molOpt.nElectronicSteps*molOpt.nMidpointFockSteps*molOpt.nNuclearSteps+1);
        rt->intScheme.nSteps = molOpt.nElectronicSteps;
        rt->intScheme.tMax = rt->intScheme.nSteps * rt->intScheme.deltaT;

        int printLevel = -1;
        try {
          printLevel = input.getData<int>("RT.PRINTLEVEL");
        } catch(...) { }

        md->gradientGetter = [&, printLevel, rt](){
          rt->printLevel = printLevel;
          return rt->getGrad(emPert);
        };

        if( auto neoss = std::dynamic_pointer_cast<NEOBase>(ss) ) {
          // FIXME: Generalize this to account for more than two subsystems
          std::vector<IntegralsBase*> ints;
          std::vector<BasisSet*> bases;
          BasisSet* ebasis = nullptr;
          BasisSet* pbasis = nullptr;
          auto labels = neoss->getLabels();
          for( auto label: labels ) {
            auto subss = neoss->getSubSSBase(label);
            ints.push_back(extractIntPtr(subss));
            bases.push_back(&subss->basisSet());
            if( label == "Electronic" )
              ebasis = bases.back();
            else if( label == "Protonic" )
              pbasis = bases.back();
          }

          md->finalMidpointFock = [=, &mol, &emPert](double t){
            for( auto isub = 0; isub < ints.size(); isub++ ) {
              bases[isub]->updateNuclearCoordinates(mol);
              ints[isub]->computeAOTwoE(*bases[isub], mol, emPert);
            }
            epints->computeAOTwoE(*ebasis, *pbasis, mol, emPert);
            rt->formCoreH(emPert);
            rt->updateAOProperties(t);
            return rt->totalEnergy();
          };
        }
        else {
          auto aoints = extractIntPtr(ss);
          BasisSet* basis = &ss->basisSet();
          md->finalMidpointFock = [&, aoints, basis, rt](double t){
            basis->updateNuclearCoordinates(mol);
            aoints->computeAOTwoE(*basis, mol, emPert);
            rt->formCoreH(emPert);
            rt->updateAOProperties(t);
            return rt->totalEnergy();
          };
        }

        elecJob = JobType::RT;
      }

    }
    else if( job == JobType::RT ) {

      rt = CQRealTimeOptions(out,input,ss,emPert);

      rt->savFile = ss->savFile;
      rt->createRTDataSets(0);
      // Single point job
      MolecularOptions molOpt(0.0, 0.0);
      mol.geometryModifier = std::make_shared<SinglePoint>(molOpt);
      elecJob = JobType::RT;

    }

    return elecJob;
  }

#undef ADD_GRAD_INCORE
#undef ADD_GRAD_DIRECT

}
