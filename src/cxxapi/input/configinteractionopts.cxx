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
#include <cxxapi/options.hpp>
#include <cerr.hpp>
#include <configinteraction.hpp>
#include <fockbuilder/rofock.hpp>
#include <regex>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQCI_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "JOBTYPE",
      "NROOTS",
      "NACTORB",
      "NACTELEC",
      "OCCRESTRICTIONS",
      "REFERENCEOCC",
      "MAXINTERSPACEEX",
      "CIDIAGALG",
      "CICONV",
      "CISIGMA2EALG",
      "MAXCIITER",
      "MAXSCFITER",
      "STATEAVERAGE",
      "SCFENECONV",
      "SCFGRADCONV",
      "SCFALG",
      "ROTATENEGORBS",
      "HESSDIAGSCALE",
      "ACTIVEORBITAL",
      "SWAPMO",
      // "INORBITAL",
      // "FVORBITAL",
      "POPULATION",
      "OSCISTREN",
      "GENIVO",
      "PRINTMOS",
      "PRINTRDMS",
      "MAXDAVIDSONSPACE",
      "NDAVIDSONGUESS",
      // Parse in the future
      // "RAS1MAXHOLE",
      // "RAS3MAXELEC",
      // "CASORBITAL",
      // "RAS1ORBITAL",
      // "RAS2ORBITAL",
      // "RAS3ORBITAL",

      "DAS"
    };

    // Specified keywords
    std::vector<std::string> ciKeywords = input.getDataInSection("CI");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : ciKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword CI." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  \brief Construct a CI object using the input file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] ss     SingleSlater reference
   *                     
   *
   *  \returns shared_ptr to a RealTimeBase object
   *    constructed from the input options.
   *
   */ 
  std::shared_ptr<PostHartreeFockBase> CQCIOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> & ss, EMPerturbation& scfPert ) {

    if( not input.containsSection("CI") )
      CErr("CI section must be specified for CI job",std::cout);
    
    std::string jobType;
    
    try {
      jobType = input.getData<std::string>("CI.JOBTYPE");
    } catch(...) {
      CErr("A specific job Type is needed for CI job");
    }
    
    //trim spaces
    trim(jobType);

    // MC Methods Keywords
    std::vector<std::string> MCMethods {
      "CI",
      "SCF"
    };
    
    // Construct valid job types 
    std::vector<std::string> CASJobs, RASJobs, GASJobs, SelectedJobs;
    for (auto &m: MCMethods) {
      CASJobs.emplace_back("CAS" + m);
      RASJobs.emplace_back("RAS" + m);
      GASJobs.emplace_back("GAS" + m);
      SelectedJobs.emplace_back("SELECTED" + m);
    }
     
    // Determine Scheme    
    bool isCASJob = 
      std::find(CASJobs.begin(),CASJobs.end(),jobType) != CASJobs.end();
    bool isRASJob = 
      std::find(RASJobs.begin(),RASJobs.end(),jobType) != RASJobs.end();
    bool isGASJob =
      std::find(GASJobs.begin(),GASJobs.end(),jobType) != GASJobs.end();
    bool isSelectedJob = 
      std::find(SelectedJobs.begin(), SelectedJobs.end(), jobType) != SelectedJobs.end();
    
    if(not isCASJob and not isRASJob and not isGASJob and not isSelectedJob) 
      CErr(jobType + " is not a valid CI.JOBTYPE",out);

    if(isSelectedJob)
      CErr(jobType + " NYI");
    
    // erase scheme and get methods
    if(isSelectedJob) jobType.erase(0,8);
    else jobType.erase(0,3);
    
    bool isCI  = not jobType.compare("CI");
    bool isSCF = not jobType.compare("SCF");

    // RASSCF NYI
    if( isRASJob and isSCF )
      CErr("RASSCF not yet implemented",out);

    // 1c+RAS NYI
    if( ss->nC==1 and isRASJob )
      CErr("1c + RAS not yet implemented",out);

    // See if reference is RO for USCF check
    #define IsRO(MT,IT) \
    std::dynamic_pointer_cast<SingleSlater<MT,IT>>(ss) ? (std::dynamic_pointer_cast<ROFock<MT,IT>>(std::dynamic_pointer_cast<SingleSlater<MT,IT>>(ss)->fockBuilder) != nullptr) : false
    bool isRO = IsRO(double,double) || IsRO(double,dcomplex) || IsRO(dcomplex,dcomplex);

    // USCF MOs NYI
    if( ss->nC==1 and not ss->iCS and not isRO ){
      CErr("Unrestricted MOs not yet implemented",out);
    }

    // parse number of roots
    size_t nR; 
	  try {
      nR = input.getData<int>("CI.NROOTS"); 
    } catch (...) {
      nR = 1;
    }
    
    // parse space partition
    std::string sActO;
    std::vector<size_t> nActOs;
	  size_t nActE;

	  try {
      sActO = input.getData<std::string>("CI.NACTORB");
    } catch(...) {
      CErr("Must specify CI.NActO for # active orbitals");
    }
    try {
      nActE = input.getData<int>("CI.NACTELEC");
    } catch(...) {
      CErr("Must specify CI.NActE for # active electrons");
    }

    // XSLI: obsolete old code that defines the number of orbitals in ea active space
    // In the new parser, we just need to know the total number of active orbitals
    std::vector<std::string> nactoTokens;
    split(nactoTokens, sActO, " ,;");
    for (auto & nacto_i: nactoTokens) nActOs.push_back(std::stoul(nacto_i));
    // construct object
    std::shared_ptr<PostHartreeFockBase> ci = nullptr;
    CISettings * ciSettings;   
    
    // Construct CI object
    #define CONSTRUCT_CI_OBJ(_MT,_IT)             \
    if( not found ) try { \
	auto derived_ss = std::dynamic_pointer_cast<SingleSlater<_MT,_IT>>(ss); \
        if (derived_ss) { \
          auto ciObj = std::make_shared<ConfigurationInteraction<_MT,_IT>>( derived_ss, nR); \
          ci = std::dynamic_pointer_cast<PostHartreeFockBase>(ciObj); \
          ciSettings = &(ciObj->ciSettings); \
          found = true;  \
        } \
      } catch(...) { }

    bool found = false;
    CONSTRUCT_CI_OBJ(double,double);   
    CONSTRUCT_CI_OBJ(dcomplex,double); 
    CONSTRUCT_CI_OBJ(dcomplex,dcomplex);

    /*
     * TODO: to have abilities of splitting the space
     */  
    //size_t maxNActO = *std::max_element(nActOs.begin(), nActOs.end());

    // OPTOPT( ci.FourCompNoPair = input.getData<bool>("CI.FOURCOMPNOPAIR"));
    
    // set up scheme
    // if      (isCASJob) ci->MOPartition.scheme = CAS;
    // else if (isRASJob) ci->MOPartition.scheme = RAS;
    

    // if(not isRASJob) {
    //   if (nActOs.size() != 1) CErr("Wrong input of CI.NACTO for" + jobType);
    // } else if (isRASJob) {
    //   if (nActOs.size() != 3) CErr("Wrong input of CI.NACTO for" + jobType);
    //   try {
    //     ci->MOPartition.mxHole = input.getData<int>("CI.RAS1MAXHOLE");
    //   } catch(...) {
    //     CErr("Must specify CI.RAS1MAXHOLE for a RAS job");
    //   }
    //   try {
    //     ci->MOPartition.mxElec = input.getData<int>("CI.RAS3MAXELEC");
    //   } catch(...) {
    //     CErr("Must specify CI.RAS3MAXELEC for a RAS Job");
    //   }
    // }   
    
    std::cout << std::endl << std::endl << std::endl << std::endl;

    std::cout << "           ************************************" << std::endl;      
    std::cout << "           *                                  *" << std::endl;      
    std::cout << "           *  Configuration Interaction (CI)  *" << std::endl;      
    std::cout << "           *                                  *" << std::endl;      
    std::cout << "           ************************************" << std::endl;      
    
    std::cout << std::endl <<BannerTop << std::endl;
	
    ci->setupCorrelatedMOSpace(std::accumulate(nActOs.begin(), nActOs.end(), 0), nActE, 0, 0);
    size_t corrOOffset = ci->corrSpace.nNegMO + ci->corrSpace.nFCore + ci->corrSpace.nInact;
    OPTOPT( ciSettings->maxInterSpaceEx = input.getData<int>("CI.MAXINTERSPACEEX"); )
    ConstructActiveSpaces(out, input, nActOs, nActE, corrOOffset, ciSettings->maxInterSpaceEx,
                          ciSettings->activeSpaces, ciSettings->refOcc, "CI");
    //ReadReferenceOcc(out, input, ciSettings->refOcc, "CI");

#if 0   
    // MO Selections or Swaps
    std::string casMOStrings, icMOStrings, svMOStrings;
    std::vector<std::string> rasMOStrings(3);
    OPTOPT( casMOStrings    = input.getData<std::string>("CI.CASORBITAL"));
    // OPTOPT( rasMOStrings[0] = input.getData<std::string>("CI.RAS1ORBITAL"));
    // OPTOPT( rasMOStrings[1] = input.getData<std::string>("CI.RAS2ORBITAL"));
    // OPTOPT( rasMOStrings[2] = input.getData<std::string>("CI.RAS3ORBITAL"));
    OPTOPT( icMOStrings     = input.getData<std::string>("CI.INORBITAL"));
    OPTOPT( svMOStrings     = input.getData<std::string>("CI.FVORBITAL"));

    #define SET_ORBITAL_INDEX(ORBINDEX, INPUTSTRING, C) \
      if (not ORBINDEX.empty()) { \
        std::vector<std::string> moTokens; \
        split(moTokens, INPUTSTRING, ", "); \
        for (auto & mo: moTokens)  { \
          std::vector<std::string> mo2; \
          split(mo2, mo, "-"); \
          if (mo2.size() == 1) { \
            ORBINDEX[std::stoul(mo2[0])-1] = C; \
          } else if (mo2.size() == 2) { \
            for (auto i = std::stoul(mo2[0]); i <= std::stoul(mo2[1]); i++) \
              ORBINDEX[i-1] = C; \
          } else CErr("Unrecogonized pattern in orbital selection"); }}
    
    #define FILL_DEFAULT_INDEX(ORBINDEX, ITER, C, N) \
      { for (auto i = 0ul; i < N; i++) { \
          while (ITER < ORBINDEX.size()) { \
            if (ORBINDEX[ITER] == 'N') break; \
            ITER++; \
          } \
          ORBINDEX[ITER] = C; }}
    
    // move the whole thing to HandlePostHFSelectMO
    bool selectMO = not icMOStrings.empty() or not svMOStrings.empty();
     
    if (isCASJob) 
      selectMO = selectMO or not casMOStrings.empty();
    else if (isRASJob)
      selectMO = selectMO or not rasMOStrings[0].empty() or 
                 not rasMOStrings[1].empty() or not rasMOStrings[2].empty();

    if (selectMO) {

      std::cout << "  * Selecting Active Space Explicitly:" << std::endl;
      
      // accomondate cases for no no-pair approximation
      size_t fourCOffSet = ci->MOPartition.nNegMO;  
      
      std::vector<char> inputOrbIndices(ci->MOPartition.nMO, 'N');
      
      // parse input
      SET_ORBITAL_INDEX(inputOrbIndices, icMOStrings, 'I');
      SET_ORBITAL_INDEX(inputOrbIndices, svMOStrings, 'S');
      if (isCASJob) {
        SET_ORBITAL_INDEX(inputOrbIndices, casMOStrings, 'A');
      } else if (isRASJob) {
        for (auto i = 0; i < 3; i++) {
          char i_char = '1' + i;
          SET_ORBITAL_INDEX(inputOrbIndices, rasMOStrings[i], i_char);
        }
      }
      
      // fill default index for those undefined ones
      size_t mo_iter = fourCOffSet;
      if (icMOStrings.empty()) {
        size_t n_char = ci->MOPartition.nInact;
        FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'I', n_char);
      }
      
      if (isCASJob) {
        if (casMOStrings.empty()) {
          size_t n_char = ci->MOPartition.nCorrO;
          FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'A', n_char);
        }
      } else if (isRASJob) {
        for (auto i = 0; i < 3; i++) { 
          char i_char = 'A';
          size_t n_char = ci->MOPartition.nActOs[i];
          if (rasMOStrings[i].empty())
            FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, i_char, n_char);
        }
      }
      if (svMOStrings.empty()) 
        FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'S', ci->MOPartition.nFVirt);
      
      std::cout << std::endl;

      std::cout << "    Construct Orbital Indices as:" << std::endl;
      
      // print 10 per line
      for (auto i = 0ul, sPerLine = 10ul; i < inputOrbIndices.size(); i++) {
        
        if (i % sPerLine == 0) 
          std::cout << "      MO " << std::setw(5) << i + 1 << " ~ " 
                    << std::setw(5) 
                    << std::min(i + sPerLine, inputOrbIndices.size()) << ":    ";
        
        std::cout << inputOrbIndices[i] << "  ";

        if ( (i+1) % sPerLine == 0) std::cout << std::endl;
      } 
      
      std::cout << std::endl << std::endl;
      
      ci->MOPartition.orbIndices = inputOrbIndices;
      ci->setActiveSpaceAndReOrder();
    }    
#endif



    // Parse CI Options

    // Change default based on # determinants
    std::string ciALG;
    OPTOPT( ciALG = input.getData<std::string>("CI.CIDIAGALG");)
    trim(ciALG);

    if( not ciALG.compare("FULLMATRIX") ) {
      ciSettings->ciAlg = CIDiagonalizationAlgorithm::CI_FULL_MATRIX;
    } else if( not ciALG.compare("DAVIDSON") ) {
      ciSettings->ciAlg = CIDiagonalizationAlgorithm::CI_DAVIDSON;
      OPTOPT( ciSettings->maxCIIter = 
                input.getData<size_t>("CI.MAXCIITER"); )
      OPTOPT( ciSettings->ciVectorConv = 
                input.getData<double>("CI.CICONV"); )
      OPTOPT( ciSettings->maxDavidsonSpace = 
                input.getData<size_t>("CI.MAXDAVIDSONSPACE");)
      OPTOPT( ciSettings->nDavidsonGuess = 
                input.getData<size_t>("CI.NDAVIDSONGUESS");)
    } else if(not ciALG.empty())
      CErr(ciALG + "is not a valid CI.CIDIAGALG",out);
    
    std::string ciSigma2eALG = "NAIVE";
    OPTOPT( ciSigma2eALG = input.getData<std::string>("CI.CISIGMA2EALG");)
    trim(ciSigma2eALG);
    ciSettings->ciSigma2eContAlg = ciSigma2eALG;
    
    // Parse Orbital Rotation Options
    if (isSCF) {
      
      ciSettings->doSCF = true;
      
      if (ss->nC == 4) {
        // default as true
        ciSettings->ORSettings.rotate_negative_positive = true;
        OPTOPT(ciSettings->ORSettings.rotate_negative_positive
          = input.getData<bool>("CI.ROTATENEGORBS"); )
      }

      OPTOPT(ciSettings->doIVOs = input.getData<bool>("CI.GENIVO"); )

      bool StateAverage = false;
      OPTOPT( StateAverage = input.getData<bool>("CI.STATEAVERAGE");)
      if(StateAverage) {
        //TODO: make as input in the future
        std::vector<double> SAWeights = std::vector<double>(nR, 1./nR);
        ci->turnOnStateAverage(SAWeights);
      }
      
      size_t maxSCFIter = 128; // default
      OPTOPT( maxSCFIter = input.getData<size_t>("CI.MAXSCFITER"); )
      ciSettings->maxSCFIter = maxSCFIter;
       
      auto & ORSettings = ciSettings->ORSettings;
      
      std::string scfALG = "AQ2nd";
      
      OPTOPT( scfALG = input.getData<std::string>("CI.SCFALG");)
       
      if( not scfALG.compare("AQ2nd") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_APPROX_QUASI_2ND_ORDER;
      } else if( not scfALG.compare("Q2nd") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_QUASI_2ND_ORDER;
      } else if( not scfALG.compare("2nd") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_2ND_ORDER;
        CErr("Second Order method is not implemented yet");
      } else {
        CErr(scfALG + " is not a valid CI.SCFALG",out);
      }

      OPTOPT( ciSettings->scfEnergyConv = 
                input.getData<double>("CI.SCFENECONV"); )
      
      OPTOPT( ciSettings->scfGradientConv = 
                input.getData<double>("CI.SCFGRADCONV"); )
      
      OPTOPT( ORSettings.hessianDiagScale = 
              input.getData<double>("CI.HESSDIAGSCALE"); )
     
   } // SCF Options

   // Properties 
   HandlePostHFProperties(out, input, ci, "CI");

   // RDMs
   HandlePostHFRDMPrinting(out, input, ci, "CI");

   // MO swapping
   // Should occur after active orbital selection
   HandlePostHFOrbitalSwaps(out, input, ss, ci, "CI");
   
   return ci;

  }; // CQCIOptions

}; // namespace ChronusQ

