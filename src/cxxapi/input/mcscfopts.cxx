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
      "CUBE"
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

  size_t HandleNRootsInput(std::string nroots,
                    std::vector<std::pair<double, size_t>> &energyRefs) {

    size_t nRoots = 0;
    bool lowRoots = false;
    std::vector<std::string> nRTokens;
    std::istringstream nRStream(nroots);
    for( std::string line; std::getline(nRStream, line); ) {
      std::cout<<line<<std::endl;
      split(nRTokens, line, " \t,");

      if( nRTokens.size() == 0) continue;
      else if ( nRTokens.size() == 1) {
        // only allow one entry of the number of low energy roots
        if (!lowRoots) {
          nRoots += std::stoul(nRTokens[0]);
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

  std::unordered_map<std::string,int> MCSCFSpinMap = {
    { "A" , 0  },
    { "B" , 1  }
  };

  void HandleMCSCFOrbitalSwaps(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> &ss, std::shared_ptr<MCWaveFunctionBase> &mcwf) {

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
   *  \brief Construct a MCSCF object using the input 
   *  file.
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
  std::shared_ptr<MCWaveFunctionBase> CQMCSCFOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> &ss, EMPerturbation& scfPert, std::shared_ptr<CubeGen> cube ) {

    if( not input.containsSection("MCSCF") )
      CErr("MCSCF section must be specified for MCSCF job",out);

    std::string jobType;
    
    try {
      jobType = input.getData<std::string>("MCSCF.JOBTYPE");
    } catch(...) {
      CErr("A specific job Type is needed for MCSCF job");
    }
    
    //trim spaces
    trim(jobType);

    // MC Methods Keywords
    std::vector<std::string> MCMethods {
      "CI",
      "SCF"
      // "PT2"
      // "PDFT",
      // "SELECTIVE",
    };
    
    // Construct valid job types 
    std::vector<std::string> CASJobs, RASJobs, DMRGJobs;
    for (auto &m: MCMethods) {
      CASJobs.emplace_back("CAS" + m);
      RASJobs.emplace_back("RAS" + m);
      DMRGJobs.emplace_back("DMRG" + m);
    }
     
    // Determine Scheme    
    bool isCASJob = 
      std::find(CASJobs.begin(),CASJobs.end(),jobType) != CASJobs.end();
    bool isRASJob = 
      std::find(RASJobs.begin(),RASJobs.end(),jobType) != RASJobs.end();
    bool isDMRGJob = 
      std::find(DMRGJobs.begin(),DMRGJobs.end(),jobType) != DMRGJobs.end();
    
    if(not isCASJob and not isRASJob and not isDMRGJob) 
      CErr(jobType + " is not a valid MCSCF.JOBTYPE",out);
    
    // erase scheme and get methods
    if(isDMRGJob) jobType.erase(0,4);
    else jobType.erase(0,3);
    
    bool isCI  = not jobType.compare("CI");
    bool isSCF = not jobType.compare("SCF");
    bool isPT2 = not jobType.compare("PT2");
    bool isPDFT = not jobType.compare("PDFT");

    // Disabling functionality NYI
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

    // construct object
    std::shared_ptr<MCWaveFunctionBase> mcscf;
    MCSCFSettings * mcscfSettings;   
    
    // parse number of roots
    size_t nR = 1;
    std::string nRoots;
    std::vector<std::pair<double, size_t>> EnergyRefs;
    OPTOPT(nRoots = input.getData<std::string>("MCSCF.NROOTS");)
    if ( not nRoots.empty() ) {
      nR = HandleNRootsInput(nRoots, EnergyRefs);
    }
    
    // Construct MCSCF object
    #define CONSTRUCT_MCSCF(_MT,_IT)             \
    if( not found ) try {                          \
      mcscf = std::dynamic_pointer_cast<MCWaveFunctionBase>( \
          std::make_shared<MCSCF<_MT,_IT>>( \
            dynamic_cast<SingleSlater<_MT,_IT>& >(*ss), nR)); \
      mcscfSettings = &(std::dynamic_pointer_cast<MCSCF<_MT,_IT>>(mcscf)->settings); \
      found = true;                 \
	} catch(...) { }

    
    bool found = false;
    if (isCI or isSCF) {
      CONSTRUCT_MCSCF( double,   double   );
      CONSTRUCT_MCSCF( dcomplex, double   );
      CONSTRUCT_MCSCF( dcomplex, dcomplex );
    } else {
      CErr("Not Implemented yet");
    }

    // OPTOPT( mcscf.FourCompNoPair = input.getData<bool>("MCSCF.FOURCOMPNOPAIR"));
    
    // set up scheme
    if      (isCASJob) mcscf->MOPartition.scheme = CAS;
    else if (isRASJob) mcscf->MOPartition.scheme = RAS;
    
    // Parse space partition
    std::string sActO;
    std::vector<size_t> nActO;
	size_t nActE;

	try {
      sActO = input.getData<std::string>("MCSCF.NACTO");
    } catch(...) {
      CErr("Must specify MCSCF.NActO for # active orbitals");
    }
    try {
      nActE = input.getData<int>("MCSCF.NACTE");
    } catch(...) {
      CErr("Must specify MCSCF.NActE for # active electrons");
    }
    
    std::vector<std::string> nactoTokens;
    split(nactoTokens, sActO, " ,;");
    for (auto & nacto_i: nactoTokens)
      nActO.push_back(std::stoul(nacto_i));

    if(not isRASJob) {
      if (nActO.size() != 1) CErr("Wrong input of MCSCF.NACTO for" + jobType);
    } else if (isRASJob) {
      if (nActO.size() != 3) CErr("Wrong input of MCSCF.NACTO for" + jobType);
      try {
        mcscf->MOPartition.mxHole = input.getData<int>("MCSCF.RAS1MAXHOLE");
      } catch(...) {
        CErr("Must specify MCSCF.RAS1MAXHOLE for a RAS job");
      }
      try {
        mcscf->MOPartition.mxElec = input.getData<int>("MCSCF.RAS3MAXELEC");
      } catch(...) {
        CErr("Must specify MCSCF.RAS3MAXELEC for a RAS Job");
      }
    }   
    
    std::cout << std::endl << std::endl << std::endl << std::endl;

    std::cout << "           *********************************************************" << std::endl;      
    std::cout << "           *                                                       *" << std::endl;      
    std::cout << "           *  Multi Configurational Self Consistent Field (MCSCF)  *" << std::endl;      
    std::cout << "           *                                                       *" << std::endl;      
    std::cout << "           *********************************************************" << std::endl;      
    
    std::cout << std::endl <<BannerTop << std::endl;
	
    mcscf->partitionMOSpace(nActO, nActE);
    
    // MO Selections or Swaps
    std::string casMOStrings, fcMOStrings, fvMOStrings;
    std::vector<std::string> rasMOStrings(3);
    OPTOPT( casMOStrings    = input.getData<std::string>("MCSCF.CASORBITAL"));
    OPTOPT( rasMOStrings[0] = input.getData<std::string>("MCSCF.RAS1ORBITAL"));
    OPTOPT( rasMOStrings[1] = input.getData<std::string>("MCSCF.RAS2ORBITAL"));
    OPTOPT( rasMOStrings[2] = input.getData<std::string>("MCSCF.RAS3ORBITAL"));
    OPTOPT( fcMOStrings     = input.getData<std::string>("MCSCF.INORBITAL"));
    OPTOPT( fvMOStrings     = input.getData<std::string>("MCSCF.FVORBITAL"));

//    #define SET_ORBITAL_INDEX(ORBINDEX, INPUTSTRING, C) \
//      if (not ORBINDEX.empty()) { \
//        std::vector<std::string> moTokens; \
//        split(moTokens, INPUTSTRING, ", "); \
//        for (auto & mo: moTokens)  { \
//          std::vector<std::string> mo2; \
//          split(mo2, mo, "-"); \
//          if (mo2.size() == 1) { \
//            ORBINDEX[std::stoul(mo2[0])-1] = C; \
//          } else if (mo2.size() == 2) { \
//            for (auto i = std::stoul(mo2[0]); i <= std::stoul(mo2[1]); i++) \
//              ORBINDEX[i-1] = C; \
//          } else CErr("Unrecogonized pattern in orbital selection"); }}
//    
//    #define FILL_DEFAULT_INDEX(ORBINDEX, ITER, C, N) \
//      { for (auto i = 0ul; i < N; i++) { \
//          while (ITER < ORBINDEX.size()) { \
//            if (ORBINDEX[ITER] == 'N') break; \
//            ITER++; \
//          } \
//          ORBINDEX[ITER] = C; }}

    bool selectMO = not fcMOStrings.empty() or not fvMOStrings.empty();
    
    if (isCASJob or isDMRGJob) 
      selectMO = selectMO or not casMOStrings.empty();
    else if (isRASJob)
      selectMO = selectMO or not rasMOStrings[0].empty() or 
                 not rasMOStrings[1].empty() or not rasMOStrings[2].empty();

    if (selectMO) {

      std::cout << "  * Selecting Active Space Explicitly:" << std::endl;
      
      // accommodate cases for no no-pair approximation
      size_t fourCOffSet = mcscf->MOPartition.nNegMO;  
      
      std::vector<char> inputOrbIndices(mcscf->MOPartition.nMO, 'N');
      
      // parse input
//      SET_ORBITAL_INDEX(inputOrbIndices, fcMOStrings, 'I');
//      SET_ORBITAL_INDEX(inputOrbIndices, fvMOStrings, 'S');
      set_orbital_index(inputOrbIndices, fcMOStrings, 'I');
      set_orbital_index(inputOrbIndices, fvMOStrings, 'S');
      if (isCASJob or isDMRGJob) {
//        SET_ORBITAL_INDEX(inputOrbIndices, casMOStrings, 'A');
        set_orbital_index(inputOrbIndices, casMOStrings, 'A');
      } else if (isRASJob) {
        for (auto i = 0; i < 3; i++) {
          char i_char = '1' + i;
//          SET_ORBITAL_INDEX(inputOrbIndices, rasMOStrings[i], i_char);
          set_orbital_index(inputOrbIndices, rasMOStrings[i], i_char);
        }
      }
      
      // fill default index for those undefined ones
      size_t mo_iter = fourCOffSet;
      if (fcMOStrings.empty()) {
        size_t n_char = mcscf->MOPartition.nInact;
//        FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'I', n_char);
        fill_default_index(inputOrbIndices, mo_iter, 'I', n_char);
      }
      
      if (isCASJob or isDMRGJob) {
        if (casMOStrings.empty()) {
          size_t n_char = mcscf->MOPartition.nCorrO;
//          FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'A', n_char);
          fill_default_index(inputOrbIndices, mo_iter, 'A', n_char);
        }
      } else if (isRASJob) {
        for (auto i = 0; i < 3; i++) { 
          char i_char = '1' + i;
          size_t n_char = mcscf->MOPartition.nActOs[i];
          if (rasMOStrings[i].empty())
//            FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, i_char, n_char);
            fill_default_index(inputOrbIndices, mo_iter, i_char, n_char);
        }
      }
      if (fvMOStrings.empty()) 
//        FILL_DEFAULT_INDEX(inputOrbIndices, mo_iter, 'S', mcscf->MOPartition.nFVirt);
        fill_default_index(inputOrbIndices, mo_iter, 'S', mcscf->MOPartition.nFVirt);

      print_orbIndices(inputOrbIndices);
//      std::cout << std::endl;
//
//      std::cout << "    Construct Orbital Indices as:" << std::endl;
//      
//      // print 10 per line
//      for (auto i = 0ul, sPerLine = 10ul; i < inputOrbIndices.size(); i++) {
//        
//        if (i % sPerLine == 0) 
//          std::cout << "      MO " << std::setw(5) << i + 1 << " ~ " 
//                    << std::setw(5) 
//                    << std::min(i + sPerLine, inputOrbIndices.size()) << ":    ";
//        
//        std::cout << inputOrbIndices[i] << "  ";
//
//        if ( (i+1) % sPerLine == 0) std::cout << std::endl;
//      } 
//      
//      std::cout << std::endl << std::endl;
      
      mcscf->MOPartition.orbIndices = inputOrbIndices;
      mcscf->setActiveSpaceAndReOrder();
    }    

    OPTOPT( mcscf->readCI = input.getData<bool>("MCSCF.READCI");)

    if( mcscf->readCI and nR>1 ) CErr("READCI not implemented for more than 1 state");
 
    // Parse CI Options
    if (isCI or isSCF) {

      // Change default based on # determinants
      std::string ciALG;
      if( mcscf->NDet<750 ) ciALG = "FULLMATRIX";
      else ciALG = "DAVIDSON";
      OPTOPT( ciALG = input.getData<std::string>("MCSCF.CIDIAGALG");)
      trim(ciALG);

      if( not ciALG.compare("FULLMATRIX") ) {
        mcscfSettings->ciAlg = CIDiagonalizationAlgorithm::CI_FULL_MATRIX;
      } else if( not ciALG.compare("DAVIDSON") ) {
        mcscfSettings->ciAlg = CIDiagonalizationAlgorithm::CI_DAVIDSON;
        OPTOPT( mcscfSettings->maxCIIter = 
                  input.getData<size_t>("MCSCF.MAXCIITER"); )
        OPTOPT( mcscfSettings->ciVectorConv = 
                  input.getData<double>("MCSCF.CICONV"); )
        OPTOPT( mcscfSettings->maxDavidsonSpace = 
                  input.getData<size_t>("MCSCF.MAXDAVIDSONSPACE");)
        OPTOPT( mcscfSettings->nDavidsonGuess = 
                  input.getData<size_t>("MCSCF.NDAVIDSONGUESS");)
        if (!EnergyRefs.empty()) mcscfSettings->energyRefs = EnergyRefs;
      } else CErr(ciALG + "is not a valid MCSCF.CIDIAGALG",out);
      
    } // CI Options
    
    // Parse Orbital Rotation Options
    if (isSCF) {
      
      mcscfSettings->doSCF = true;
      
      if (ss->nC == 4) {
        // default as true
        mcscfSettings->ORSettings.rotate_negative_positive = true;
        OPTOPT(mcscfSettings->ORSettings.rotate_negative_positive
          = input.getData<bool>("MCSCF.ROTATENEGORBS"); )
      }

      OPTOPT(mcscfSettings->doIVOs = input.getData<bool>("MCSCF.GENIVO"); )

      bool StateAverage = false;
      OPTOPT( StateAverage = input.getData<bool>("MCSCF.STATEAVERAGE");)
      if(StateAverage) {
        std::vector<double> SAWeights = HandleSAWeightsInput(input, nR);
        mcscf->turnOnStateAverage(SAWeights);
      }
      
      size_t maxSCFIter = 128; // default
      OPTOPT( maxSCFIter = input.getData<size_t>("MCSCF.MAXSCFITER"); )
      mcscfSettings->maxSCFIter = maxSCFIter;
       
      auto & ORSettings = mcscfSettings->ORSettings;
      
      std::string scfALG = "AQ2ND";
      
      OPTOPT( scfALG = input.getData<std::string>("MCSCF.SCFALG");)
       
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
                input.getData<double>("MCSCF.SCFENECONV"); )
      
      OPTOPT( mcscfSettings->scfGradientConv = 
                input.getData<double>("MCSCF.SCFGRADCONV"); )
      
      OPTOPT( ORSettings.hessianDiagScale = 
              input.getData<double>("MCSCF.HESSDIAGSCALE"); )
     
   } // SCF Options

   // Mulliken charge analysis
   OPTOPT( mcscf->PopulationAnalysis = input.getData<bool>("MCSCF.POPULATION"); )

   // Mulliken charge analysis
   OPTOPT( mcscf->SpinAnalysis = input.getData<bool>("MCSCF.PRINTSPIN"); )

   // Oscillator strength
   OPTOPT( mcscf->NosS1 = input.getData<size_t>("MCSCF.OSCISTREN"); )

   // Multipole moments
   OPTOPT( mcscf->multipoleMoment = input.getData<bool>("MCSCF.PRINTMULT"); )

   // MCSCF Field
   std::string fieldStr;
   OPTOPT(
      fieldStr = input.getData<std::string>("MCSCF.FIELD");
   )
   EMPerturbation parsedField;
   handleField(fieldStr, parsedField, scfPert);
   mcscf->mcscfPert.addField(parsedField);

   if( pert_has_type(mcscf->mcscfPert,Magnetic) ) {
      if (ss->basisSet_.basisType != COMPLEX_GIAO) 
        CErr("NYI - MCSCF with magnetic field only works with GIAO!");
   } else if (pert_has_type(mcscf->mcscfPert, Electric) && ss->nC == 4) {
       CErr("NYI - 4C MCSCF with electric field might work, but unverified.");
   }

   // Printing Options
   // MOs
   if ( input.containsData("MCSCF.PRINTMOS") ) {
     try { mcscf->printMOCoeffs = input.getData<size_t>("MCSCF.PRINTMOS"); }
     catch(...) {
       CErr("Invalid PRINTMOS input. Please use number 0 ~ 9.");
     }
   }
   if (mcscf->printMOCoeffs >= 10 ) CErr("MCSCF print level is not valid!");

   // RDMs
   HandleRDMPrinting(out, input,  mcscf);

   // MO swapping
   // Should occur after active orbital selection
   HandleMCSCFOrbitalSwaps(out, input, ss, mcscf);

   if( cube ){

     auto &cubeOptions = cube->getCubeOptions();
     mcscf->cubeOptsMC = cubeOptions;

   }

   // Check for CubeGen subsection
   if( input.containsSection("MCSCF.CUBE") ){
     std::cout << " Found [MCSCF.CUBE] section" << std::endl;
     CQCUBE_VALID(out,input,"MCSCF.");
     CQCUBEOptionalKeywords(out,input,mcscf->cubeOptsMC,"MCSCF.");
   }
   
   return mcscf;

  }; // CQCIOptions

}; // namespace ChronusQ

