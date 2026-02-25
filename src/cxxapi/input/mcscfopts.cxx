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
#include <physcon.hpp>
#include <mcwavefunction/base.hpp>
#include <mcscf.hpp>
#include <fockbuilder/rofock.hpp>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQMCSCF_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "JOBTYPE",
      "NROOTS",
      "NACTO",
      "NACTE",
      "NACTP",
      "NATORB",
      "NATORBREDIAG",
      "RAS1MAXHOLE",
      "RAS3MAXELEC",
      "READCI",
      "CIDIAGALG",
      "CICONV",
      "MAXCIITER",
      "MAXSCFITER",
      "STATEAVERAGE",
      "SAWEIGHTS",
      "SCFENECONV",
      "SCFGRADCONV",
      "SCFALG",
      "ROTATENEGORBS",
      "HESSDIAGSCALE",
      "CASORBITAL",
      "RAS1ORBITAL",
      "RAS2ORBITAL",
      "RAS3ORBITAL",
      "SWAPMO",
      "INORBITAL",
      "FVORBITAL",
      "POPULATION",
      "OSCISTREN",
      "GENIVO",
      "FIELD",
      "PRINTMOS",
      "PRINTRDMS",
      "PRINTSPIN",
      "MAXDAVIDSONSPACE",
      "NDAVIDSONGUESS",
      "PRINTMULT",
      "CUBE",
      "NDETPRINT",
    };

    // Specified keywords
    std::vector<std::string> mcscfKeywords = input.getDataInSection("MCSCF");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : mcscfKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword MCSCF." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  \brief Handle NRoots input parsing for multiple-root calculation.
   *  \details Parses the NRoots input string to extract energy references and number of roots.
   *         Supports both low-energy roots and energy-specific roots with specified thresholds.
   *         The format of the input string include lines with either a single number (for low-energy roots)
   *         or pairs of numbers (energy threshold and number of roots) separated by spaces, tabs, or commas.
   *         e.g.:
   *         NROOTS = "5
   *         0.3  3
   *         2.5  2"
   *         This will return 10 for a total of 10 roots, and
   *         energyRefs will contain: {{-inf, 5}, {0.3, 3}, {2.5, 2}}
   *         for 5 low-energy roots, 3 roots above 0.3, and 2 roots above 2.5.
   *
   *  \param nroots String input for NRoots.
   *  \param energyRefs Vector of pairs to store energy references and number of roots.
   *         Will be cleared at the beginning and filled during parsing as return.
   *         Each pair consists of (energy reference, number of roots).
   *
   *  \return Total number of roots requested.
   *
   */
  size_t HandleNRootsInput(std::string nroots,
                    std::vector<std::pair<double, size_t>> &energyRefs) {

    energyRefs.clear();

    size_t nRoots = 0;
    size_t nlowRoots = 0;
    bool lowRoots = false;
    std::vector<std::string> nRTokens;
    std::istringstream nRStream(nroots);
    for( std::string line; std::getline(nRStream, line); ) {
      //std::cout<<line<<std::endl;
      split(nRTokens, line, " \t,");

      if( nRTokens.size() == 0) continue;
      else if ( nRTokens.size() == 1) {
        // only allow one entry of the number of low energy roots
        if (!lowRoots) {
          nRoots += std::stoul(nRTokens[0]);
          nlowRoots += std::stoul(nRTokens[0]);
          lowRoots = true;
          continue;
        }
        else CErr("Too many input for number of low energy roots.");
      }
      else if ( nRTokens.size() > 2 ) CErr("Too many parameters for Energy Specific Settings in one line.");

      // default energy unit is a.u.
      std::regex EThreshold("([+-]?[0-9]+([.][0-9]*)?|[.][0-9]+)");
      std::smatch ethres;
      std::regex_search(nRTokens[0], ethres, EThreshold);
      double Ethres;
      if (ethres.str(1).size()>0) {
        Ethres = std::stod(ethres.str(1));
      }

      auto const regexAU = std::regex("a.u.|au|Hartree",std::regex_constants::icase);
      auto const regexEV = std::regex("ev|electronvolt",std::regex_constants::icase);
      auto const regexNM = std::regex("nm|nanometer",std::regex_constants::icase);
      double unit = 1.;
      if ( std::regex_search(nRTokens[0], regexAU) ) unit = 1.;
      else if ( std::regex_search(nRTokens[0], regexEV) ) unit = EVPerHartree;
      else if ( std::regex_search(nRTokens[0], regexNM) ) unit = NMPerHartree;

      energyRefs.emplace_back(Ethres / unit, std::stoul(nRTokens[1]));
      nRoots += std::stoul(nRTokens[1]);
    }
    
    if (lowRoots){
      // lowRoots is True so we did something like 
    // nRoots:
    //     NROOTS:
    // 5
    // 0.   5
    // 100. 10
    // This should yield 10 roots at 0.0 and 10 roots at 100.0
      bool found = false;
      for (auto& [energy, energy_roots] : energyRefs)
      {
        if (energy == 0.0 and not found){
          energy_roots += nlowRoots;
          found = true;
        }
      }
      if (not found) {
        energyRefs.emplace_back(-std::numeric_limits<double>::infinity(), nlowRoots);
      }
    }

    return nRoots;

  } // HandleNRootsInput

  std::vector<double> HandleSAWeightsInput(CQInputFile &input,
    size_t nR) {

    std::string sSAWeights;
    std::vector<double> SAWeights = std::vector<double>(nR, 1./nR);

    OPTOPT( sSAWeights = input.getData<std::string>("MCSCF.SAWEIGHTS"); )
    if (!sSAWeights.empty()) {
      std::cout << "  * Manual State Average Weights Detected: " << std::endl;
      //Parse SAWeights string
      std::vector<std::string> saTokens;
      split(saTokens, sSAWeights, " ,;");
      if (saTokens.size() != nR) CErr("Number of Input State Average Weights does not match number of roots.");
      for (auto i=0; i < nR; i++)
        SAWeights[i] = std::stod(saTokens[i]);
    }

    // post process
    double sum = std::accumulate(SAWeights.begin(), SAWeights.end(), 0.);
    if (sum != 1) {
      std::cout << "   Rescale State Average Weights Sum to 1." << std::endl;
      for (auto i = 0ul; i < SAWeights.size(); i++)
        SAWeights[i] = SAWeights[i]/sum;
    }
    return SAWeights;

  } // HandleSAWeightsInput

  void HandleRDMPrinting(std::ostream &out, CQInputFile &input,
    std::shared_ptr<MCWaveFunctionBase> &mcwf) {

    // Parse RDM printing
    std::string printRDMString;
    OPTOPT( printRDMString = input.getData<std::string>("MCSCF.PRINTRDMS"));
    if ( not printRDMString.empty() ) {
      std::cout << "  * Printing RDM detected: " << std::endl;

      std::vector<std::string> rdmTokens;
      split(rdmTokens, printRDMString, " \t,");

      if( rdmTokens.size() != 1 and rdmTokens.size() != 2 ) CErr("Need 1 or 2 entries in single line for RDM printing");

      // Parse rdmCut if present
      if( rdmTokens.size() == 2 ) mcwf->rdmCut=std::stod(trim(rdmTokens[1]));

      try { mcwf->printRDMs = std::stoi(trim(rdmTokens[0])); }
      catch(...) {
        CErr("Invalid PRINTRDMS input. Please use number 0 ~ 2.");
      }
      if (mcwf->printRDMs >= 3 ) CErr("MCSCF print RDM level is not valid!");

    }

  }

  void HandleDetPrinting(std::ostream & out, CQInputFile & input,
    std::shared_ptr<MCWaveFunctionBase> & mcwf)
  {
    std::string DetPrintString;
    OPTOPT(DetPrintString = input.getData<std::string>("MCSCF.NDETPRINT"));

    // If not specified return
    if(DetPrintString.empty())
    {return;}
    if(DetPrintString=="ALL")
    {
      std::cout << "Printing the entire CI Vector...Beware!" << std::endl;
      mcwf->NDetPrint=DetPrint::ALLDET;
      return;
    }
    try{mcwf->NDetPrint=std::stoi(DetPrintString);}
    catch(...)
    {
      CErr("Unrecognized options for NDETPRINT");
    }
    return;
  }


  std::unordered_map<std::string,int> MCSCFSpinMap = {
    { "A" , 0  },
    { "B" , 1  }
  };

  void HandleMCSCFOrbitalSwaps(std::ostream &out, CQInputFile &input,
    std::shared_ptr<MCWaveFunctionBase> &mcwf) {

    // MO swapping
    std::string swapMOStrings;
    OPTOPT( swapMOStrings = input.getData<std::string>("MCSCF.SWAPMO"));
    if ( not swapMOStrings.empty() ) {
      std::cout << "  * Manually MO Swapping Detected: " << std::endl;

      // Pair function for [MCSCF] MO swap
      std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;
      moPairs.resize(2, {});

      std::vector<std::string> moTokens;
      //Loop over lines of mo swapping
      std::istringstream moStream(swapMOStrings);

      for( std::string line; std::getline(moStream, line); ) {
        split(moTokens, line, " \t,");

        if( moTokens.size() == 0 ) continue;
        else if( moTokens.size() != 2 and moTokens.size() != 3 ) CErr("Need 2 or 3 entries in single line for swapping");

        // Parse spin if present
        std::string spinDir("A");
        if( moTokens.size() == 3 ) spinDir=moTokens[2];
        trim(spinDir);

        // mo[1] NYI for MCSCF
        if( spinDir == "B" ) CErr("Swapping of beta MOs NYI for MCSCF");

        moPairs[MCSCFSpinMap[spinDir]].emplace_back(std::stoul(moTokens[0]), std::stoul(moTokens[1]));
      }

      mcwf->swapMOs(moPairs,isAlpha);

    }

  }

  void set_orbital_index(std::vector<char> &orbindex, std::string inputstring, char C) {

    if (not orbindex.empty()) {
      std::vector<std::string> moTokens;
      split(moTokens, inputstring, ", ");
      for (auto & mo: moTokens) {
        std::vector<std::string> mo2;
        split(mo2, mo, "-");
        if (mo2.size() == 1) {
          orbindex[std::stoul(mo2[0])-1] = C;
        } else if (mo2.size() == 2) {
          for (auto i = std::stoul(mo2[0]); i <= std::stoul(mo2[1]); i++)
            orbindex[i-1] = C;
        } else CErr("Unrecogonized pattern in orbital selection"); 
      }
    }
  }

  void fill_default_index(std::vector<char> &orbindex, size_t iter, char C, size_t N) {

    for (auto i = 0ul; i < N; i++) {
      while (iter < orbindex.size()) {
        if (orbindex[iter] == 'N') break;
        iter++;
      }
      orbindex[iter] = C;
    }
  }

  void print_orbIndices(std::vector<char> orbindex) {

    std::cout << std::endl;
    std::cout << "    Construct Orbital Indices as:" << std::endl;

    // print 10 per line
    for (auto i = 0ul, sPerLine = 10ul; i < orbindex.size(); i++) {
      if (i % sPerLine == 0)
        std::cout << "      MO " << std::setw(5) << i + 1 << " ~ "
                  << std::setw(5)
                  << std::min(i + sPerLine, orbindex.size()) << ":    ";

        std::cout << orbindex[i] << "  ";

        if ( (i+1) % sPerLine == 0) std::cout << std::endl;
      }

      std::cout << std::endl << std::endl;

  }


  /**
   *
   * \brief Parse an orbital selection string of comma-separated indices or ranges
   * \example An input "1, 3-5, 7" will return a vector of [0,2,3,4,6]
   *
   * \param [in] input_str Input string to be parsed
   *
   * \return A vector of 0-based orbital indices
   *
   */
  std::vector<size_t> parseOrbitalSelectionInput(std::string input_str) {
    std::vector<size_t> orbitals;
    std::vector<std::string> moTokens;
    split(moTokens, input_str, ", ");
    for (auto & mo: moTokens)  {
      std::vector<std::string> mo2;
      split(mo2, mo, "-");
      if (mo2.size() == 1) {
        orbitals.push_back(std::stoul(mo2[0])-1);
      } else if (mo2.size() == 2) {
        for (size_t i = std::stoul(mo2[0]), iEnd = std::stoul(mo2[1]); i <= iEnd; i++)
          orbitals.push_back(i-1);
      } else CErr("Unrecognized pattern in orbital selection");
    }
    std::sort(orbitals.begin(), orbitals.end());
    orbitals.erase(std::unique(orbitals.begin(), orbitals.end()), orbitals.end());
    return orbitals;
  }

  /**
   *
   * \brief Parse an orbital selection vector to a string of comma-separated indices or ranges
   * \example An input vector of [0,2,3,4,6] will return a "1, 3-5, 7" string
   *
   * \param [in] orbitals Input string to be parsed
   *
   * \return A string of 1-based comma-separated indices or ranges
   *
   */
  std::string orbitalSelectionToString(std::vector<size_t> orbitals) {
    if (orbitals.empty()) return "";

    size_t rangeStart = orbitals[0];
    std::string s = std::to_string(rangeStart + 1);

    for (size_t i = 1; i < orbitals.size(); i++)  {
      if (orbitals[i] - 1 != orbitals[i-1]) {
        if (orbitals[i-1] != rangeStart)
          s += "-" + std::to_string(orbitals[i-1] + 1);
        rangeStart = orbitals[i];
        s += "," + std::to_string(rangeStart + 1);
      }
    }

    if (orbitals.back() != rangeStart)
      s += "-" + std::to_string(orbitals.back() + 1);

    return s;
  }


  /**
   *
   * \brief Primary interface for which an MCSCF job is built
   *
   * \return A MCSCFBase pointer from which MCSCF calculations are run
   *
   */

  std::shared_ptr<MCSCFBase> CQMCSCFOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> &ss, std::shared_ptr<MCWaveFunctionBase> & mcwfn, EMPerturbation& scfPert, std::shared_ptr<CubeGen> cube, bool doNEO ) 
  {

    std::cout << std::endl << std::endl << std::endl << std::endl;

    std::cout << "           *********************************************************" << std::endl;      
    std::cout << "           *                                                       *" << std::endl;      
    std::cout << "           *  Multi Configurational Self Consistent Field (MCSCF)  *" << std::endl;      
    std::cout << "           *                                                       *" << std::endl;      
    std::cout << "           *********************************************************" << std::endl;      
    
    std::cout << std::endl << BannerTop << std::endl;
 

    std::string prefix = "";
    std::shared_ptr<MCWaveFunctionBase> auxmcwfn;
    std::shared_ptr<MCSCFJobType> mcscfjob = nullptr;
    if(!doNEO)
    {
      mcwfn = CQBuildMCWaveFunction(out,input,ss,scfPert,cube,prefix,mcscfjob);
    }
    else
    {
      std::shared_ptr<NEOBase> neobase;
      std::shared_ptr<SingleSlaterBase> tempss;

      #define GETNEOBASE(_MT,_IT)             \
      if( not found ) try {                          \
        neobase = std::dynamic_pointer_cast<NEOBase>( \
            std::dynamic_pointer_cast<NEOSS<_MT,_IT>>(ss)); \
        found = true;                 \
      } catch(...) { }

      bool found = false;
      GETNEOBASE(double,double);
      GETNEOBASE(dcomplex,double);
      GETNEOBASE(dcomplex,dcomplex);

      tempss = neobase->getSubSSBase("Electronic");
      auxmcwfn = CQBuildMCWaveFunction(out,input,tempss,scfPert,cube,prefix,mcscfjob);

      #define MAKENEOMCWFN(_MT,_IT) \
      if(not found) try { \
       mcwfn = std::dynamic_pointer_cast<MCWaveFunctionBase>( \
              std::make_shared<NEOMCWaveFunction<_MT,_IT>>( \
                dynamic_cast<NEOSS<_MT,_IT>&>(*ss),auxmcwfn->NStates)); \
       found = true; \
       } catch(...) { }

      found = false;
      MAKENEOMCWFN(double,double);
      MAKENEOMCWFN(dcomplex,double);
      MAKENEOMCWFN(dcomplex,dcomplex);

      mcwfn->addMCWaveFunction(auxmcwfn,"Electronic");

      prefix = "PROT";
      tempss = neobase->getSubSSBase("Protonic");
      auxmcwfn = CQBuildMCWaveFunction(out,input,tempss,scfPert,cube,prefix,mcscfjob);
      mcwfn->addMCWaveFunction(auxmcwfn,"Protonic");
      prefix="";

    }

    return CQBuildMCSCFOptions(out,input,mcwfn,scfPert,cube,prefix,mcscfjob);
  }


  // A struct for convenience for parsing the exact job type requested by the user
  struct MCSCFJobType
  {
    // Data members 
    std::string jobType;
    bool isCASJob;
    bool isRASJob;
    bool isDMRGJob;

    bool isCI;
    bool isSCF;
    bool isPT2;
    bool isPDFT;

    // MC Methods Keywords
    std::vector<std::string> MCMethods {
      "CI",
      "SCF"
      // "PT2"
      // "PDFT",
      // "SELECTIVE",
    };

    void parseMCSCFJobType(std::ostream&out,CQInputFile&input,std::string&prefix)
    {
      try {
        jobType = input.getData<std::string>(prefix+"MCSCF.JOBTYPE");
      } catch(...) {
        CErr("A specific job Type is needed for MCSCF job");
      }
      //trim spaces
      trim(jobType);

      // Construct valid job types 
      std::vector<std::string> CASJobs, RASJobs, DMRGJobs;
      for (auto &m: MCMethods) {
        CASJobs.emplace_back("CAS" + m);
        RASJobs.emplace_back("RAS" + m);
        DMRGJobs.emplace_back("DMRG" + m);
      }

      // Determine Scheme    
      isCASJob = std::find(CASJobs.begin(),CASJobs.end(),jobType) != CASJobs.end();
      isRASJob = std::find(RASJobs.begin(),RASJobs.end(),jobType) != RASJobs.end();
      isDMRGJob = std::find(DMRGJobs.begin(),DMRGJobs.end(),jobType) != DMRGJobs.end();

      if(!isCASJob and !isRASJob and !isDMRGJob) 
        CErr(jobType + " is not a valid MCSCF.JOBTYPE",out);
      
        // erase scheme and get methods
      if(isDMRGJob) jobType.erase(0,4);
      else jobType.erase(0,3);
    
      isCI  = not jobType.compare("CI");
      isSCF = not jobType.compare("SCF");
      isPT2 = not jobType.compare("PT2");
      isPDFT = not jobType.compare("PDFT");

      // Disabling functionality NYI
      // RASSCF NYI
      if( isRASJob and isSCF )
        CErr("RASSCF not yet implemented");

    };
  };

  // Builds the MCSCF options
  // Takes a MCWaveFunctionBase object, which may be of type MCWaveFunction or MultiComponentMCWaveFunction
  std::shared_ptr<MCSCFBase> CQBuildMCSCFOptions(std::ostream &out,
    CQInputFile &input, std::shared_ptr<MCWaveFunctionBase> &mcwfn,  EMPerturbation& scfPert, std::shared_ptr<CubeGen> cube, std::string & prefix,std::shared_ptr<MCSCFJobType>& mcscfjob) {

    if( not input.containsSection(prefix+"MCSCF") )
      CErr("MCSCF section must be specified for MCSCF job",out);

    if(!mcscfjob)
    {
      mcscfjob = std::make_shared<MCSCFJobType>();
      mcscfjob->parseMCSCFJobType(out,input,prefix);
    }

    // construct object
    std::shared_ptr<MCSCFBase> mcscf;
    MCSCFSettings * mcscfSettings;   
    bool found = false;    
    #define CONSTRUCT_MCSCF(_MT,_IT)             \
    if( not found ) try {                          \
      if(std::dynamic_pointer_cast<MCWaveFunction<_MT,_IT>>(mcwfn)){ \
      mcscf = std::dynamic_pointer_cast<MCSCFBase>( \
          std::make_shared<MCSCF<_MT,_IT>>( \
            std::dynamic_pointer_cast<MCWaveFunction<_MT,_IT>>(mcwfn))); \
      found = true;}                \
  	} catch(...) { }

    CONSTRUCT_MCSCF(double,double);
    CONSTRUCT_MCSCF(dcomplex,double);
    CONSTRUCT_MCSCF(dcomplex,dcomplex);

    mcscfSettings = &mcscf->settings;

    std::string nRoots;
    size_t nR;
    std::vector<std::pair<double, size_t>> EnergyRefs;
    OPTOPT(nRoots = input.getData<std::string>(prefix+"MCSCF.NROOTS");)
    if ( not nRoots.empty() ) {
      nR = HandleNRootsInput(nRoots, EnergyRefs);
    }

    if( mcwfn->readCI and nR > 1 ) CErr("READCI not implemented for more than 1 state");
 
    // Parse CI Options
    if (mcscfjob->isCI or mcscfjob->isSCF) {

      // Check if number of roots is valid
      if( nR > mcscf->NDet ) CErr("# roots > # determinants");

      // Change default based on # determinants
      std::string ciALG;
      if( mcscf->NDet<750 ) ciALG = "FULLMATRIX";
      else ciALG = "DAVIDSON";
      OPTOPT( ciALG = input.getData<std::string>(prefix+"MCSCF.CIDIAGALG");)
      trim(ciALG);

      if( not ciALG.compare("FULLMATRIX") ) {
        mcscfSettings->ciAlg = CIDiagonalizationAlgorithm::CI_FULL_MATRIX;
      } else if( not ciALG.compare("DAVIDSON") ) {
        mcscfSettings->ciAlg = CIDiagonalizationAlgorithm::CI_DAVIDSON;
        OPTOPT( mcscfSettings->maxCIIter = 
                  input.getData<size_t>(prefix+"MCSCF.MAXCIITER"); )
        OPTOPT( mcscfSettings->ciVectorConv = 
                  input.getData<double>(prefix+"MCSCF.CICONV"); )
        OPTOPT( mcscfSettings->maxDavidsonSpace = 
                  input.getData<size_t>(prefix+"MCSCF.MAXDAVIDSONSPACE");)
        OPTOPT( mcscfSettings->nDavidsonGuess = 
                  input.getData<size_t>(prefix+"MCSCF.NDAVIDSONGUESS");)
        if (!EnergyRefs.empty()) mcscfSettings->energyRefs = EnergyRefs;
      } else if( not ciALG.compare("SKIP") )
      {
        mcscfSettings->ciAlg = CIDiagonalizationAlgorithm::SKIP;
      } 
      else CErr(ciALG + "is not a valid MCSCF.CIDIAGALG",out);
      
    } // CI Options
    
    // Parse Orbital Rotation Options
    if (mcscfjob->isSCF) {
      
      mcscfSettings->doSCF = true;
      
      if (mcwfn->getnC() == 4) {
        // default as true
        mcscfSettings->ORSettings.rotate_negative_positive = true;
        OPTOPT(mcscfSettings->ORSettings.rotate_negative_positive
          = input.getData<bool>(prefix+"MCSCF.ROTATENEGORBS"); )
      }

      OPTOPT(mcscfSettings->doIVOs = input.getData<bool>(prefix+"MCSCF.GENIVO"); )

      bool StateAverage = (nR > 1);
      OPTOPT( StateAverage = input.getData<bool>(prefix+"MCSCF.STATEAVERAGE");)
      if(StateAverage) {
        std::vector<double> SAWeights = HandleSAWeightsInput(input, nR);
        mcscf->turnOnStateAverage(SAWeights);
      }
      
      size_t maxSCFIter = 128; // default
      OPTOPT( maxSCFIter = input.getData<size_t>(prefix+"MCSCF.MAXSCFITER"); )
      mcscfSettings->maxSCFIter = maxSCFIter;
       
      auto & ORSettings = mcscfSettings->ORSettings;
      
      std::string scfALG = "AQ2ND";
      
      OPTOPT( scfALG = input.getData<std::string>(prefix+"MCSCF.SCFALG");)
       
      if( not scfALG.compare("AQ2ND") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_APPROX_QUASI_2ND_ORDER;
      } else if( not scfALG.compare("Q2ND") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_QUASI_2ND_ORDER;
        CErr("Quasi Second Order method is not implemented yet");
      } else if( not scfALG.compare("2ND") ) {
        ORSettings.alg = OrbitalRotationAlgorithm::ORB_ROT_2ND_ORDER;
        CErr("Second Order method is not implemented yet");
      } else {
        CErr(scfALG + " is not a valid MCSCF.SCFALG",out);
      }

      OPTOPT( mcscfSettings->scfEnergyConv = 
                input.getData<double>(prefix+"MCSCF.SCFENECONV"); )
      
      OPTOPT( mcscfSettings->scfGradientConv = 
                input.getData<double>(prefix+"MCSCF.SCFGRADCONV"); )
      
      OPTOPT( ORSettings.hessianDiagScale = 
              input.getData<double>(prefix+"MCSCF.HESSDIAGSCALE"); )
     
   } // SCF Options

   // Natural orbitals
   OPTOPT( mcscfSettings->NatOrbs = input.getData<int>(prefix+"MCSCF.NATORB"); )
   OPTOPT( mcscfSettings->NatOrbRediag = input.getData<bool>(prefix+"MCSCF.NATORBREDIAG"); )

   // Mulliken charge analysis
   OPTOPT( mcscfSettings->PopulationAnalysis = input.getData<bool>(prefix+"MCSCF.POPULATION"); )

   // Mulliken charge analysis
   OPTOPT( mcscfSettings->SpinAnalysis = input.getData<bool>(prefix+"MCSCF.PRINTSPIN"); )

   // Oscillator strength
   OPTOPT( mcscfSettings->NosS1 = input.getData<size_t>(prefix+"MCSCF.OSCISTREN"); )

   // Multipole moments
   OPTOPT( mcscfSettings->multipoleMoment = input.getData<bool>(prefix+"MCSCF.PRINTMULT"); )

   // MCSCF Field
   std::string fieldStr;
   OPTOPT( fieldStr = input.getData<std::string>(prefix+"MCSCF.FIELD");)
   EMPerturbation parsedField;
   handleField(fieldStr, parsedField, scfPert);
   mcscf->mcscfPert.addField(parsedField);

   if( pert_has_type(mcscf->mcscfPert,Magnetic) ) {
      if (mcwfn->getBasisType() != COMPLEX_GIAO) 
        CErr("NYI - MCSCF with magnetic field only works with GIAO!");
   } else if (pert_has_type(mcscf->mcscfPert, Electric) && mcwfn->getnC() == 4) {
       CErr("NYI - 4C MCSCF with electric field might work, but unverified.");
   }

  return mcscf;

  }; // CQCIOptions


  /**
   *  \brief Construct a MCWaveFunction object using the input 
   *  file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] ss     SingleSlater reference
   *                     
   *
   *  \returns shared_ptr to a MCWaveFunction object
   *    constructed from the input options.
   *
   */ 
  std::shared_ptr<MCWaveFunctionBase> CQBuildMCWaveFunction(std::ostream & out,
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> & ss, EMPerturbation & scfPert, std::shared_ptr<CubeGen> cube, std::string & prefix,std::shared_ptr<MCSCFJobType>& mcscfjob)
  {
    // construct object
    std::shared_ptr<MCWaveFunctionBase> mcwfn;
     
    // parse number of roots
    size_t nR = 1;
    std::string nRoots;
    std::vector<std::pair<double, size_t>> EnergyRefs;
    OPTOPT(nRoots = input.getData<std::string>(prefix+"MCSCF.NROOTS");)
    if ( not nRoots.empty() ) {
      nR = HandleNRootsInput(nRoots, EnergyRefs);
    }

    // Construct the MCWaveFuction object
    #define CONSTRUCT_MCWFN(_MT,_IT)             \
    if( not found ) try {                          \
      mcwfn = std::dynamic_pointer_cast<MCWaveFunctionBase>( \
          std::make_shared<MCWaveFunction<_MT,_IT>>( \
            dynamic_cast<SingleSlater<_MT,_IT>& >(*ss), nR)); \
      found = true;                 \
    } catch(...) { }

    bool found = false;
    CONSTRUCT_MCWFN( double,   double   );
    CONSTRUCT_MCWFN( dcomplex, double   );
    CONSTRUCT_MCWFN( dcomplex, dcomplex );

    if(!mcscfjob)
    {
      mcscfjob = std::make_shared<MCSCFJobType>();
      mcscfjob->parseMCSCFJobType(out,input,prefix);
    }

   // 1c+RAS NYI
    if( ss->nC==1 and mcscfjob->isRASJob )
      CErr("1c + RAS not yet implemented",out);

    // See if reference is RO for USCF check
    #define IsRO(MT,IT) \
    (std::dynamic_pointer_cast<SingleSlater<MT,IT>>(ss) ? ([](auto const& ss_sp) { \
        auto fb = ss_sp->fockBuilder; \
        if (auto neo = std::dynamic_pointer_cast<NEOFockBuilder<MT,IT>>(fb)) { \
          /* getNonNEOUpstream() returns raw pointer; do not re-own it */ \
          auto* up = neo->getNonNEOUpstream(); \
          return dynamic_cast<ROFock<MT,IT>*>(up) != nullptr; \
        } \
        return std::dynamic_pointer_cast<ROFock<MT,IT>>(fb) != nullptr; \
      })(std::dynamic_pointer_cast<SingleSlater<MT,IT>>(ss)) \
    : false)
    bool isRO = IsRO(double,double) || IsRO(double,dcomplex) || IsRO(dcomplex,dcomplex);

    // USCF MOs NYI
    if( ss->nC==1 and not ss->iCS and not isRO and ss->particle.charge < 0){
      CErr("Unrestricted MOs not yet implemented",out);
    }

    // Build the MCWaveFunction

    // set up scheme
    if      (mcscfjob->isCASJob) mcwfn->MOPartition.scheme = CAS;
    else if (mcscfjob->isRASJob) mcwfn->MOPartition.scheme = RAS;

    // Parse space partition
    std::string sActO;
    std::vector<size_t> nActO;
  	size_t nActE = 0;
    size_t nActP = 0;

    try {
      sActO = input.getData<std::string>(prefix+"MCSCF.NACTO");
    } catch(...) {
      CErr("Must specify MCSCF.NActO for # active orbitals");
    }
    OPTOPT(nActE = input.getData<int>(prefix+"MCSCF.NACTE"));
    OPTOPT(nActP = input.getData<int>(prefix+"MCSCF.NACTP"));
    if(!nActE && !nActP)
      CErr("Must specify active number of particles in MCSCF with NACTE or NACTP");
   
    std::vector<std::string> nactoTokens;
    split(nactoTokens, sActO, " ,;");
    for (auto & nacto_i: nactoTokens)
      nActO.push_back(std::stoul(nacto_i));

    if(not mcscfjob->isRASJob) {
      if (nActO.size() != 1) CErr("Wrong input of MCSCF.NACTO for" + mcscfjob->jobType);
    } else if (mcscfjob->isRASJob) {
      if (nActO.size() != 3) CErr("Wrong input of MCSCF.NACTO for" + mcscfjob->jobType);
      try {
        mcwfn->MOPartition.mxHole = input.getData<int>(prefix+"MCSCF.RAS1MAXHOLE");
      } catch(...) {
        CErr("Must specify MCSCF.RAS1MAXHOLE for a RAS job");
      }
      try {
        mcwfn->MOPartition.mxElec = input.getData<int>(prefix+"MCSCF.RAS3MAXELEC");
      } catch(...) {
        CErr("Must specify MCSCF.RAS3MAXELEC for a RAS Job");
      }
    }   

    mcwfn->partitionMOSpace(nActO,ss->particle.charge < 0 ? nActE : nActP);

    // MO Selections or Swaps
    std::string casMOStrings, fcMOStrings, fvMOStrings;
    std::vector<std::string> rasMOStrings(3);
    OPTOPT( casMOStrings    = input.getData<std::string>(prefix+"MCSCF.CASORBITAL"));
    OPTOPT( rasMOStrings[0] = input.getData<std::string>(prefix+"MCSCF.RAS1ORBITAL"));
    OPTOPT( rasMOStrings[1] = input.getData<std::string>(prefix+"MCSCF.RAS2ORBITAL"));
    OPTOPT( rasMOStrings[2] = input.getData<std::string>(prefix+"MCSCF.RAS3ORBITAL"));
    OPTOPT( fcMOStrings     = input.getData<std::string>(prefix+"MCSCF.INORBITAL"));
    OPTOPT( fvMOStrings     = input.getData<std::string>(prefix+"MCSCF.FVORBITAL"));

    bool selectMO = not fcMOStrings.empty() or not fvMOStrings.empty();
    
    if (mcscfjob->isCASJob or mcscfjob->isDMRGJob) 
      selectMO = selectMO or not casMOStrings.empty();
    else if (mcscfjob->isRASJob)
      selectMO = selectMO or not rasMOStrings[0].empty() or 
                 not rasMOStrings[1].empty() or not rasMOStrings[2].empty();

    if (selectMO) {

      std::cout << "  * Selecting Active Space Explicitly:" << std::endl;
      
      // accommodate cases for no no-pair approximation
      size_t fourCOffSet = mcwfn->MOPartition.nNegMO;  
      
      std::vector<char> inputOrbIndices(mcwfn->MOPartition.nMO, 'N');
      
      // parse input
      set_orbital_index(inputOrbIndices, fcMOStrings, 'I');
      set_orbital_index(inputOrbIndices, fvMOStrings, 'S');
      if (mcscfjob->isCASJob or mcscfjob->isDMRGJob) {
        set_orbital_index(inputOrbIndices, casMOStrings, 'A');
      } else if (mcscfjob->isRASJob) {
        for (auto i = 0; i < 3; i++) {
          char i_char = '1' + i;
          set_orbital_index(inputOrbIndices, rasMOStrings[i], i_char);
        }
      }
      
      // fill default index for those undefined ones
      size_t mo_iter = fourCOffSet;
      if (fcMOStrings.empty()) {
        size_t n_char = mcwfn->MOPartition.nInact;
        fill_default_index(inputOrbIndices, mo_iter, 'I', n_char);
      }
      
      if (mcscfjob->isCASJob or mcscfjob->isDMRGJob) {
        if (casMOStrings.empty()) {
          size_t n_char = mcwfn->MOPartition.nCorrO;
          fill_default_index(inputOrbIndices, mo_iter, 'A', n_char);
        }
      } else if (mcscfjob->isRASJob) {
        for (auto i = 0; i < 3; i++) { 
          char i_char = '1' + i;
          size_t n_char = mcwfn->MOPartition.nActOs[i];
          if (rasMOStrings[i].empty())
            fill_default_index(inputOrbIndices, mo_iter, i_char, n_char);
        }
      }
      if (fvMOStrings.empty()) 
        fill_default_index(inputOrbIndices, mo_iter, 'S', mcwfn->MOPartition.nFVirt);

      print_orbIndices(inputOrbIndices);
     
      mcwfn->MOPartition.orbIndices = inputOrbIndices;
      mcwfn->setActiveSpaceAndReOrder();
    }    

    // MO swapping
    // Should occur after active orbital selection
    HandleMCSCFOrbitalSwaps(out, input, mcwfn);

    // Handling of printing determinants
    HandleDetPrinting(std::cout,input,mcwfn);

    // Printing Options
    // MOs
    if ( input.containsData("MCSCF.PRINTMOS") ) {
      try { mcwfn->printMOCoeffs = input.getData<size_t>("MCSCF.PRINTMOS"); }
      catch(...) {
        CErr("Invalid PRINTMOS input. Please use number 0 ~ 9.");
      }
    }
    if (mcwfn->printMOCoeffs >= 10 ) CErr("MCSCF print level is not valid!");

    // RDMs
    HandleRDMPrinting(out, input, mcwfn);


    // ReadCI
    OPTOPT( mcwfn->readCI = input.getData<bool>("MCSCF.READCI");)

   if( cube ){
     auto &cubeOptions = cube->getCubeOptions();
     mcwfn->cubeOptsMC = cubeOptions;
   }

   // Check for CubeGen subsection
   if( input.containsSection(prefix+"MCSCF.CUBE") ){
     std::cout << " Found [MCSCF.CUBE] section" << std::endl;
     CQCUBE_VALID(out,input,prefix+"MCSCF.");
     CQCUBEOptionalKeywords(out,input,mcwfn->cubeOptsMC,prefix+"MCSCF.");
   }
 
    return mcwfn;
  }

}; // namespace ChronusQ

