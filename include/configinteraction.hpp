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

#include <chronusq_sys.hpp>
#include <cerr.hpp>
#include <singleslater.hpp>
#include <posthartreefock/base/space.hpp>
#include <mcscf.hpp>
#include <detfactory.hpp>
#include <newcibuilder.hpp>
#include <neworbitalrotation.hpp>

namespace ChronusQ {
  
struct CISettings {
   
   CIDiagonalizationAlgorithm ciAlg = CI_DEFAULT; 
   std::string ciSigma2eContAlg = "SIMPLE";

   // This stores the input from user
   // might be different from what's being used
   // in DetFactory
   std::vector<ActiveSpaceParameters> activeSpaces;
   std::vector<std::vector<size_t>> refOcc;
   int maxInterSpaceEx = -1;

   // for davidson and gplhr
   size_t maxCIIter        = 128;        
   double ciVectorConv     = 1.0e-6;
   bool checkEigenValue = true;
   bool checkEigenVector = true;
   bool checkResidue = false;
   size_t maxDavidsonSpace = 20;
   size_t nDavidsonGuess   = 3;
   std::vector<std::pair<double, size_t>> energyRefs;
   
   //Sparse Davidson
   bool SparseDavidson = false;
   double SparseDavidsonEps = 1e-12;

   //Threads
   size_t nThreads_;

   // SCF Settings 
   bool doSCF           = false;
   bool doIVOs          = false;
   
   size_t maxSCFIter         = 0;
   double scfEnergyConv      = 1.0e-8; 
   double scfGradientConv    = 1.0e-4;
   OrbitalRotationSettings ORSettings;

   CISettings() {
     ORSettings.rotate_within_correlated = false;
   }
   
   CISettings(const CISettings &) = default;
   CISettings(CISettings &&)      = default;

   void print(bool);
}; // struct CISettings 

template <typename MatsT, typename IntsT>
class ConfigurationInteraction: public PostHartreeFock<MatsT, IntsT> {

public:
  
  // potentially to be expanded to sparse vectors
  //std::shared_ptr<DistributedVectors<MatsT>> CIVectors = nullptr;
  std::shared_ptr<SolverVectors<MatsT>> CIVectors = nullptr;

  CISettings ciSettings;
  std::shared_ptr<DeterminantFactory> detFactory = nullptr;
  std::shared_ptr<NewCIBuilder<MatsT>> ciBuilder  = nullptr;
  std::shared_ptr<NewOrbitalRotation<MatsT,IntsT>> moRotator = nullptr;  

  // Reduced density Matrices (RDMs) only span over correlated space
  // SOI: state of interest, for orbital rotation
  // it's either state specific or state averaged RDM
  std::shared_ptr<cqmatrix::Matrix<MatsT>>    oneRDMSOI = nullptr;
  std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = nullptr;
  
  // Disable default, copy and move constructors
  ConfigurationInteraction()                          = delete;
  ConfigurationInteraction(const ConfigurationInteraction &) = delete;
  ConfigurationInteraction(ConfigurationInteraction &&)      = delete;

  /**
   *  \brief ConfigInteraction Constructor.
   *
   *  Stores references to a "reference" SingleSlater object and
   *  makes a copy of the reference into a complex
   *  SingleSlater object for the propagation.
   */ 
  template <typename MatsU>
  ConfigurationInteraction(std::shared_ptr<SingleSlater<MatsU,IntsT>> ref, size_t NS) :
    PostHartreeFock<MatsT,IntsT>(ref, NS) { };  // ConfigInteraction constructor

  ~ConfigurationInteraction(){ dealloc(); }

  void initialization();
  size_t nDeterminants() const { return  detFactory->ketCategoricalSpace()->nDeterminants(); }
  void run(EMPerturbation &) override;
  void computeTDM(size_t s1, size_t s2, std::shared_ptr<cqmatrix::Matrix<MatsT>> tdm) override;
  void compute2TDM(size_t s1, size_t s2, std::shared_ptr<InCore4indexTPI<MatsT>> twoTDM) override;
  void compute2RDM(size_t s1, size_t s2, std::shared_ptr<InCore4indexTPI<MatsT>> twoRDM) override;
 
  void solveCI();
  
  void computeRDMsForOrbitalRotations(); 

  void saveCurrentStates() override;

  void printStateEnergy();
  void printCIHeader();
  void printCIFooter();

  // Memory functions
  void alloc();
  void dealloc();

}; // class ConfigInteraction
   
} // namespace ChronusQ
