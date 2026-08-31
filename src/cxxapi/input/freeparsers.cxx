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
#include <cerr.hpp>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <regex>
#include <string>
#include <orbitalmodifieroptions.hpp>
#include <ksrefs.hpp>

namespace ChronusQ {

  namespace {

    std::string canonicalRealTimeAlgorithm(std::string value) {
      std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return std::toupper(c); });

      if(value == "MMUT" or value == "MODIFIEDMIDPOINT") return "MMUT";
      if(value == "FORWARDEULER" or value == "EULER") return "FORWARDEULER";
      if(value == "EXPLICITMAGNUS2" or value == "EXPLICITMAGNUSTWO" or
         value == "MAGNUS2" or value == "MAGNUSTWO") return "MAGNUS2";
      if(value == "RK4" or value == "RUNGEKUTTAFOURTHORDER") return "RK4";
      if(value == "SSO") return "SSO";
      if(value == "BORT") return "BORT";

      CErr("Unrecognized RT integration algorithm: " + value);
      return "";
    }

    RealTimeAlgorithm parseRealTimeAlgorithm(const std::string& value) {
      const auto algorithm = canonicalRealTimeAlgorithm(value);
      if(algorithm == "MMUT") return RealTimeAlgorithm::RTModifiedMidpoint;
      if(algorithm == "FORWARDEULER") return RealTimeAlgorithm::RTForwardEuler;
      if(algorithm == "MAGNUS2") return RealTimeAlgorithm::RTExplicitMagnus2;
      if(algorithm == "RK4") return RealTimeAlgorithm::RTRungeKuttaOrderFour;
      if(algorithm == "SSO") return RealTimeAlgorithm::RTSymplecticSplitOperator;
      if(algorithm == "BORT") return RealTimeAlgorithm::ElectronicBornOppenheimer;
      return RealTimeAlgorithm::Uninitialized;
    }

  }

  std::string doubleToString(double value) {
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
  }

  std::string escapeRegEx(const std::string &s) {                                 
    static const std::regex metacharacters(R"([\.\^\$\+\(\)\[\]\{\}\|\?\*])");       
    return std::regex_replace(s, metacharacters, "\\$&");                         
  }
  
  std::string expandSciNote(std::string value) {
    std::regex del("D|E");
    std::sregex_token_iterator it(value.begin(), value.end(), del, -1);
    std::sregex_token_iterator end;
    
    std::vector<std::string> arr = {"0","0"};
    int count =0;
    while (it != end){
      if (count > 1) CErr("Sci. Notation " + value + " NOT converted properly.");
      arr[count] = std::string(*it);
      ++it;
      ++count;
    }

    double number = std::stod(arr[0]) * std::pow(10,std::stod(arr[1]));
    
    return doubleToString(number);
  }//expandSciNote

  void CQInputFile::parseFreeCQInput (std::string &line){

    ////////////////////////////////////////////////////////////////////////
    // CQ Free Format Input                                               //
    // example, CQ= HF/STO-3G SCF(accuracy=1.e-6)                         //
    // example, ChronuQ: X2C-HF/CD-6-31G RT(time=10fs, stepsize = 1as)    //
    // example, ChronuQ= 4C-HF/ano-rcc hamiltonian(DCB, scalar, atomic)   //
    ////////////////////////////////////////////////////////////////////////
    //Double check identifier is not used (should be done already in parser.cxx but...
    auto const freeCQInput = std::regex("^\\?|CQ[[:blank:]]*=|CQ[[:blank:]]*:|CHRONUSQ[[:blank:]]*=|CHRONUSQ[[:blank:]]*:",std::regex_constants::icase);
    line = std::regex_replace(line, freeCQInput, "");
      
    //Exxtra Options:
    parseFreeCQInputSCF(line);
    parseFreeCQInputSSGuess(line);
    parseFreeCQInputField(line);
    parseFreeCQInputNEO(line);
    //parseFreeCQInputMisc(line);
    
    //Requirements for basic Job:
    if ( not parseFreeCQInputJob(line) ) addData("QM/JOB", "SCF");
    if ( not parseFreeCQInputRef(line) ) addData("QM/REFERENCE", "HF");
    if ( not parseFreeCQInputBas(line) ){
      std::cout << "WARNING, we did not find a valid basis, setting to STO-3G" << std::endl;
      addData("BASIS/BASIS", "STO-3G");
    }

    auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
    line = std::regex_replace(line, freeDividers, " ");
    // the second call remove multiple space left after replace ','
    line = std::regex_replace(line, freeDividers, " ");
    if(line.size()>0) std::cout<<"CQ input ignored: "<< line <<std::endl;

  }; // Free Format Input Parser

  void CQInputFile::parseFreeCQInputNEO (std::string &line){

    /*************************************/
    /* NEO Input                         */
    /* example, NEO(EPC19/prot-pb6-g)    */
    /* example, NEO(EPC19/CD-prot-pb6-g) */
    /* example, NEO(EPC19/ri-prot-pb6-g) */
    /*************************************/
    
    //TODO ability to specify real/cmoplex:
    auto const freeCQInputNEOMETH = std::regex("((2C)|(X2C)|(4C)|(G)|(U)|(R))?(-)?((HF)|(EPC17)|(EPC19))");
    auto const freeCQInputNEOBASIS= std::regex("(CD|RI)?-?(PROT-SP|PROT-PB4-D|PROT-PB4-F1|PROT-PB4-F2|PROT-PB5-F|PROT-PB5-G|PROT-PB6-G)");

    auto const freeCQInputNEO = std::regex("NEO(\\((.*?)\\))?",std::regex_constants::icase);
    std::smatch NEOmatch;

    if( std::regex_search(line, NEOmatch, freeCQInputNEO) ){
      addData("SCF/NEO", "TRUE");
      // str(2) captures what is inside NEO()
      if(NEOmatch.str(2).size()>0) {
        // Parse user-defined input
        // std::cout<<"xsli test NEO Section "<<std::endl;
        std::string NEOInputOptions = NEOmatch.str(2);

        // Methods
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputNEOMETH) ) {
          addData("PROTQM/REFERENCE",NEOmatch.str(0));
        }

        // Basis Sets
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputNEOBASIS) ) {
          if ( ! NEOmatch.str(1).empty()) {
            addData("INTS/ALG","INCORE");
            addData("INTS/RI","DYNAMICERI");
            addData("EPINTS/ALG","INCORE");
            addData("EPINTS/RI","COMBINEAUXBASIS");
          }//CD/RI
          addData("PBASIS/BASIS",NEOmatch.str(2));
        }
      } else{
        //use defaults
        std::cout<<" WARNING: no NEO specifications given, using defaults" << std::endl;
        addData("PROTQM/REFERENCE","UHF");
        addData("PBASIS/BASIS","PROT-SP");
      }

      // We need to delete the NEO section so that we can parse the electronic section properly
      line = std::regex_replace(line, freeCQInputNEO, "");
    } // NEO Input

  };

  bool CQInputFile::parseFreeCQInputJob (std::string &line){
    std::smatch methodMatchCC;
    std::smatch methodMatchLR;
    std::smatch methodMatchCI;
    std::smatch methodMatchRT;
    //Method Choices TODO
    auto const freeCQInputCC  = std::regex("CCSD(T|\\(T\\))?/?",std::regex_constants::icase);
    bool isCC = std::regex_search(line, methodMatchCC, freeCQInputCC);
    auto const freeCQInputLR = std::regex("(LR)(\\((.*?)\\))?", std::regex_constants::icase);
    bool isLR = std::regex_search(line, methodMatchLR, freeCQInputLR);
    auto const freeCQInputCI = std::regex("((2C)|(X2C)|(4C)|(G))?(-)?((CASSCF)|(DASSCF)|(CASCI)|(DASCI)|(CI))\\((([^()]*(\\([^()]*\\))?[^()]*)*)\\)", std::regex_constants::icase);
    bool isCI = std::regex_search(line, methodMatchCI, freeCQInputCI);
    auto const freeCQInputRT = std::regex("(RT|REALTIME)(\\((.*?)\\))?", std::regex_constants::icase);
    bool isRT = std::regex_search(line, methodMatchRT, freeCQInputRT);
    
    //Verify one job type:
    if ( isCC + isLR + isCI + isRT > 1 ){
      CErr("Multiple job types detected, please only choose one");
    }

    if ( isCC ) {
      //CC job
      addData("QM/JOB", "CC");
      std::cout<< "Coupled Cluster JOB IDENTIFIED: " << methodMatchCC[0].str() << std::endl;
      parseFreeCQInputCC(line, freeCQInputCI);
      line = std::regex_replace(line, freeCQInputCC, "");
      return true;
    }//CC Match
    else if ( isLR ) {
      std::smatch lrMatch;
      std::string lrInputOptions = methodMatchLR.str(3);
      //LR job
      addData("QM/JOB", "RESP");
      addData("RESPONSE/TYPE", "RESIDUE");
      std::cout<< "RESPONSE JOB IDENTIFIED: " << methodMatchLR[0].str() << std::endl;
      //LR opts:
      
      // Is the default but in here for completeness
      auto const freeCQInputLRFULL = std::regex("(FULLMAT|FULL|FULLMATRIX|DOFULL)", std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputLRFULL) ) {
        //If we ever do a froz core/virt in lin. resp. we'll need to rethink this
        addData("RESPONSE/DOFULL", "TRUE");
        addData("RESPONSE/FULLMAT", "TRUE");
      }
      
      auto const freeCQInputLRroots = std::regex("(NROOTS|NSTATES)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)",std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputLRroots) ) {
        addData("RESPONSE/DOFULL", "FALSE");
        addData("RESPONSE/FULLMAT", "FALSE");
        addData("RESPONSE/NROOTS", expandSciNote(lrMatch.str(2)));
      }
      
      auto const freeCQInputLRiter = std::regex("(NCYCLES|MAXCYCLES|MAXITER|NITER)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)",std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputLRiter) ) {
        addData("RESPONSE/DOFULL", "FALSE");
        addData("RESPONSE/FULLMAT", "FALSE");
        addData("RESPONSE/MAXITER", expandSciNote(lrMatch.str(2)));
      }
      
      auto const freeCQInputDEMIN = std::regex("(DEMIN)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)", std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputDEMIN) ) {
        //addData("RESPONSE/GPLHR_SIGMA", lrMatch.str(2));
        addData("RESPONSE/DOFULL", "FALSE");
        addData("RESPONSE/FULLMAT", "FALSE");
        addData("RESPONSE/DEMIN", expandSciNote(lrMatch.str(2)));
      }
      
      auto const freeCQInputSIGMA = std::regex("(SIGMA)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)", std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputSIGMA) ) {
        addData("RESPONSE/DOFULL", "FALSE");
        addData("RESPONSE/FULLMAT", "FALSE");
        addData("RESPONSE/GPLHR_SIGMA", expandSciNote(lrMatch.str(2)));
      }
      
      auto const freeCQInputCONV = std::regex("(CONV|ACCURACY)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)", std::regex_constants::icase);
      if ( std::regex_search(lrInputOptions, lrMatch, freeCQInputCONV) ) {
        addData("RESPONSE/CONV", expandSciNote(lrMatch.str(2)));
      }
      //Clear Line
      line = std::regex_replace(line, freeCQInputLR, "");
      return true;  
    }//LR Match
    else if ( isCI ) {
      std::cout<< "CI JOB IDENTIFIED: " << methodMatchCI[0].str() << std::endl;
      parseFreeCQInputCI(line, freeCQInputCI);
      line = std::regex_replace(line, freeCQInputCI, "");
      return true;
    }//CI Match
    else if ( isRT ){
      std::cout<< "RT JOB IDENTIFIED: " << methodMatchRT[0].str() << std::endl;
      std::cout<< "WARNING: for realtime propagator, integrals must be COMPLEX" << std::endl;
      parseFreeCQInputRT(line, freeCQInputRT);
      line = std::regex_replace(line, freeCQInputRT, "");
      return true;
    }//RT Match
    return false;
  }//parseFreeCQInputJob
  
  bool CQInputFile::parseFreeCQInputRef (std::string &line){
    std::smatch methodMatch;
    std::smatch ksMatch;
    
    //////////////////////////////
    // Reference Input          //
    // example, SF-X2C-HF       //
    // example, B3LYP           //
    // example, 4C-B3LYP        //
    // example, X2C-B3LYP       //
    // example, B3LYP           //
    // example, B3LYP           //
    //////////////////////////////
    
    auto const freeCQInputHF = std::regex("(SF)?[-]?(2C|X2C|4C)?[-]?(R|U|G)?[-]?HF/?",std::regex_constants::icase);
    //Make regex key from defined DFT keys:
    std::string allowedDFT = "(";
    for ( int i = 0; i < KSRefs.size(); i++ ){
      allowedDFT += escapeRegEx(KSRefs[i]);
      allowedDFT += "|";
    }
    allowedDFT.pop_back();
    allowedDFT += ")";

    //Test if an allowed DFT is present:
    bool foundKS = false;
    auto const dvar = std::regex(allowedDFT, std::regex_constants::icase);
    line = " " + line + " ";
    if(std::regex_search(line, ksMatch, dvar)){
      foundKS = true;
    }    
    
    if ( std::regex_search(line, methodMatch, freeCQInputHF) ) {
      //Hartree Fock
      std::cout<< "HF REFERENCE IDENTIFIED: " << methodMatch[0].str() << std::endl;
      
      //(SF)?[-]?(2C|X2C|4C)?[-]?(R|U|G)?[-]?HF/?
      //  1           2             3
      
      //Spinfree X2C
      if ( methodMatch[1].str().size() > 0 ) addData("QM/X2CTYPE","SPINFREE");
      
      //Reference (x2c || rug):
      addData("QM/REFERENCE", methodMatch[2].str() + methodMatch[3].str() + "HF" );
      
      //Clear line:
      line = std::regex_replace(line, freeCQInputHF, "");
      //Check for KS and HF refs:
      if ( foundKS ) CErr("Both HF and DFT reference found, please choose one.");
      return true;
    }//HF Reference
    else if ( foundKS ) {
      //DFT (have to do both search for DFT and then catch method)
      auto const freeCQInputDFT = std::regex("(SF)?[-]?(2C|X2C|4C)?[-]?(R|U|G)?[-]?(" + ksMatch[1].str() + ")/?",std::regex_constants::icase);
      std::regex_search(line, methodMatch, freeCQInputDFT);
      std::cout<< "DFT REFERENCE IDENTIFIED: " << methodMatch[0].str() << std::endl;
      
      //(SF)?[-]?(2C|X2C|4C)?[-]?(R|U|G)?[-]?(DFT)/?
      //  1           2             3          4
      
      //Spinfree X2C
      if ( methodMatch[1].str().size() > 0 ) addData("QM/X2CTYPE","SPINFREE");
      
      //Reference (x2c || rug):
      addData("QM/REFERENCE", methodMatch[2].str() + methodMatch[3].str() + methodMatch[4].str() );

      //DEFAULTS:
      addData("DFTINT/INHOUSE","FALSE");
      addData("DFTINT/GAUXC","TRUE");
      
      //Clear line:
      line = std::regex_replace(line, freeCQInputDFT, "");
      return true;
    }//DFT Reference
    return false;
  }//parseFreeCQInputRef
  
  bool CQInputFile::parseFreeCQInputBas (std::string &line){
    std::smatch basisMatch;
   
    //Make regex key from defined keys:
    std::string allowedBasis = "(\\ ";
    for ( const auto& n : basisKeyword ){
      allowedBasis += escapeRegEx(n.first);
      allowedBasis += "\\ |\\ ";
    }
    allowedBasis.pop_back();
    allowedBasis += "\\ )";

    //The Pople sets (at least) have substrings as part so need to put some 
    //flavor of separator so you get the one you request:
    line = " " + line + " ";
    auto const bvar = std::regex(allowedBasis, std::regex_constants::icase);
    if(std::regex_search(line, basisMatch, bvar)){
      std::string match = basisMatch[1].str();
      match.erase(match.begin());
      match.erase(match.end()-1);
      std::cout<<"Basis Specified: "<<match<<std::endl;
      addData("BASIS/BASIS", match);
      line = std::regex_replace(line, bvar, "");
      return true;
    }//if basismatch
    return false;
  }//parseFreeCQInputBas

  void CQInputFile::parseFreeCQInputSCF (std::string &line) {

    ////////////////////////////////////////////////////
    // SCF Input                                      //
    // example, SCF(accuracy=1.e-8)                   //
    // example, SCF(cdiis, maxiteration=100)          //
    // example, SCF(energyonly)                       //
    ////////////////////////////////////////////////////

    auto const freeCQInputSCF = std::regex("(SCF)(\\((.*?)\\))", std::regex_constants::icase);
    std::smatch scfMatch;

    if (std::regex_search(line, scfMatch, freeCQInputSCF)) {
      std::string scfInputOptions = scfMatch.str(3);

      // read in SCF accuracy
      auto const freeCQInputSCFAccuracy = std::regex("accuracy\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFAccuracy) ) {
        addData("SCF/ACCURACY", expandSciNote(scfMatch.str(1)));
        //std::cout<<"xsli test read in accuracy = "<<std::stod(scfMatch.str(1))<<std::endl;
      }

      auto const freeCQInputEnergyOnly = std::regex("energyonly|skip", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputEnergyOnly) ) {
        addData("SCF/ALG", "SKIP");
      }

      auto const freeCQInputSCFMethod = std::regex("(diis)|(nodiis)|(cdiis)|(ediis)|(qc)|(cediis)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFMethod) ) {
        if(!scfMatch.str(1).empty()) addData("SCF/DIISALG","CDIIS");// std::cout<<"xsli test read in SCF method = DIIS "<<std::stod(scfMatch.str(1))<<std::endl;
        if(!scfMatch.str(2).empty()) addData("SCF/DIISALG","NONE"); // std::cout<<"xsli test read in SCF method = NoDIIS "<<std::stod(scfMatch.str(2))<<std::endl;
        if(!scfMatch.str(3).empty()) addData("SCF/DIISALG","CDIIS"); // std::cout<<"xsli test read in SCF method = CDIIS "<<std::stod(scfMatch.str(3))<<std::endl;
        if(!scfMatch.str(4).empty()) addData("SCF/DIISALG","EDIIS"); // std::cout<<"xsli test read in SCF method = EnDIIS "<<std::stod(scfMatch.str(4))<<std::endl;
        if(!scfMatch.str(4).empty()) addData("SCF/ALG","NR"); // std::cout<<"xsli test read in SCF method = QC "<<std::stod(scfMatch.str(5))<<std::endl;
        if(!scfMatch.str(5).empty()) addData("SCF/DIISALG","CEDIIS"); // std::cout<<"xsli test read in SCF method = CEDIIS "<<std::stod(scfMatch.str(5))<<std::endl;
      }

      // Check for SCF maxSteps
      auto const freeCQInputSCFSteps = std::regex("(NITER|NITERATIONS|MAXSTEP|MAXSTEPS|MAXITERATION|MAXITERATIONS|MAXCYCLE|MAXCYCLES)\\s*\\=\\s*(\\d*\\.?\\d*(?:[de][+\\-]?\\d*)?)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFSteps) ) {
        addData("SCF/MAXITER", expandSciNote(scfMatch.str(2)));
        //std::cout<<"xsli test read in MAXITERATIONS = "<<std::stod(scfMatch.str(2))<<std::endl;
      }
    }

  };

#if(0)
  void SCFControls::parseSection(const std::map<std::string,std::string> &dict) {
    if (dict.count("ENERGYONLY")) {
      scfAlg = _CONVENTIONAL_SCF;
      energyOnly = true;
    }
    if (dict.count("MAXITER")) maxSCFIter = std::stoi(dict.at("MAXITER"));
    if (dict.count("ACCURACY")) {
      rmsdPConvTol = std::stod(dict.at("ACCURACY"));
      maxdPConvTol = rmsdPConvTol*100;
      eneConvTol   = rmsdPConvTol;
    }

      if (dict.at("DIISALG") == "CDIIS") diisAlg = CDIIS;
      else if (dict.at("DIISALG") == "EDIIS") diisAlg = EDIIS;
      else if (dict.at("DIISALG") == "DIIS") diisAlg = CDIIS;
      else if (dict.at("DIISALG") == "NODIIS") diisAlg = NONE;

      if (dict.at("ALG") == "QC") scfAlg = _NEWTON_RAPHSON_SCF;

      if (dict.count("MAXITER")) maxSCFIter = std::stoi(dict.at("MAXITER"));
  }
#endif

  void CQInputFile::parseFreeCQInputSSGuess (std::string &line) {

    /*****************************************************************************/
    /* Guess Input                                                               */
    /* example, Guess(core)                                                      */
    /* example, NEOGuess(read)                                                   */
    /* example, Guess(swap=a3-a4,swap=alpha-homo-lumo,swap=b5-b6,swap=11-17)  */
    /*****************************************************************************/

    SingleSlaterGuessOptions ssGuessOptions;

    auto const freeCQInputGuess = std::regex("(GUESS)(\\((.*?)\\))?", std::regex_constants::icase);
    std::smatch Guessmatch;
    if (std::regex_search(line, Guessmatch, freeCQInputGuess)) {
      // str(3) captures the input inside the parentheses
      if (Guessmatch.str(3).size() > 0) {
        std::string GuessInputOptions = Guessmatch.str(3);

        auto const freeCQInputGuessType = std::regex("(SAD)|(CORE)|(READ)|(FCHK)|(ONLY)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputGuessType) ) {
//          if(!Guessmatch.str(1).empty()) ssGuessOptions.electronicGuess = SADGuess;
//          if(!Guessmatch.str(2).empty()) ssGuessOptions.electronicGuess = CoreGuess;
//          if(!Guessmatch.str(3).empty()) ssGuessOptions.electronicGuess = ReadBin;
//          if(!Guessmatch.str(4).empty()) ssGuessOptions.electronicGuess = ReadGaussFCHK;
          if(!Guessmatch.str(1).empty()) addData("SCF/GUESS", "SAD");
          if(!Guessmatch.str(2).empty()) addData("SCF/GUESS", "CORE");
          if(!Guessmatch.str(3).empty()) addData("SCF/GUESS", "READMO");
          if(!Guessmatch.str(4).empty()) addData("SCF/GUESS", "FCHKMO");
          if(!Guessmatch.str(5).empty()) addData("SCF/ALG", "SKIP");
        }
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessType, "");

        auto const freeCQInputGuessSwap = std::regex("(SWAP)\\s*=((\\s*([ab]?)(\\d+)-\\4(\\d+))|(homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
        while( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputGuessSwap) ) {
          if(!Guessmatch.str(3).empty()) {
            // std::cout<<"xsli test guess swap = "<<Guessmatch.str(4)<<Guessmatch.str(5)<<Guessmatch.str(6)<<std::endl;
            if(Guessmatch.str(4)=="a" or Guessmatch.str(4)=="A" or Guessmatch.str(4).size()==0)
              ssGuessOptions.alphaElectronicMOSwap.push_back({(size_t)std::stoi(Guessmatch.str(5)), (size_t)std::stoi(Guessmatch.str(5))});
            else if(Guessmatch.str(4)=="b" or Guessmatch.str(4)=="B")
              ssGuessOptions.betaElectronicMOSwap.push_back({(size_t)std::stoi(Guessmatch.str(5)), (size_t)std::stoi(Guessmatch.str(5))});

            auto const freeCQInputGuessSwap1 = std::regex("(SWAP)\\s*=\\s*([ab]?)(\\d+)-\\2(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessSwap1, "");
          }
          if(!Guessmatch.str(7).empty() or !Guessmatch.str(8).empty()) {
            // std::cout<<"xsli test guess swap = HOMO-LUMO"<<std::endl;
            ssGuessOptions.alphaElectronicMOSwap.push_back({0,0});
            auto const freeCQInputGuessSwap2 = std::regex("(SWAP)\\s*=((homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessSwap2, "");
          }
        }

        for (size_t i = 0; i < ssGuessOptions.alphaElectronicMOSwap.size(); i++) {
          for (size_t j = 0; j < ssGuessOptions.alphaElectronicMOSwap[i].size(); j++) {
            addData("SCF/GUESS/ALPHASWAP[" + std::to_string(i) + "][" + std::to_string(j) + "]",
                    std::to_string(ssGuessOptions.alphaElectronicMOSwap[i][j]));
          }
        }

        for (size_t i = 0; i < ssGuessOptions.betaElectronicMOSwap.size(); i++) {
          for (size_t j = 0; j < ssGuessOptions.betaElectronicMOSwap[i].size(); j++) {
            addData("SCF/GUESS/BETASWAP[" + std::to_string(i) + "][" + std::to_string(j) + "]",
                    std::to_string(ssGuessOptions.betaElectronicMOSwap[i][j]));
          }
        }

        auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeDividers, "");
        if(!GuessInputOptions.empty()) CErr("Unrecognized Guess Input Options: "+GuessInputOptions);
        else line = std::regex_replace(line, freeCQInputGuess, "");
      }

    }

#if(0)
    auto const freeCQInputNEOGuess = std::regex("(NEOGUESS)(\\((.*?)\\))?", std::regex_constants::icase);
    if (std::regex_search(line, Guessmatch, freeCQInputNEOGuess)) {
      // str(3) captures the input inside the parentheses
      if (Guessmatch.str(3).size() > 0) {
        std::string GuessInputOptions = Guessmatch.str(3);

        auto const freeCQInputNEOGuessType = std::regex("(SAD)|(CORE)|(READ)|(FCHK)|(CLASSICAL)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputNEOGuessType) ) {
//          if(!Guessmatch.str(1).empty()) ssGuessOptions.nuclearGuess = SADGuess;
//          if(!Guessmatch.str(2).empty()) ssGuessOptions.nuclearGuess = CoreGuess;
//          if(!Guessmatch.str(3).empty()) ssGuessOptions.nuclearGuess = ReadBin;
//          if(!Guessmatch.str(4).empty()) ssGuessOptions.nuclearGuess = ReadGaussFCHK;
//          if(!Guessmatch.str(5).empty()) ssGuessOptions.nuclearGuess = ClassicalGuess;
          if(!Guessmatch.str(1).empty()) addData("SCF/QP_GUESS", "SAD");
          if(!Guessmatch.str(2).empty()) addData("SCF/QP_GUESS", "CORE");
          if(!Guessmatch.str(3).empty()) addData("SCF/QP_GUESS", "READMO");
          if(!Guessmatch.str(4).empty()) addData("SCF/QP_GUESS", "FCHKMO");
          if(!Guessmatch.str(5).empty()) addData("SCF/QP_GUESS", "CLASSICAL");
        }
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputNEOGuessType, "");

        auto const freeCQInputNEOGuessSwap = std::regex("(SWAP)\\s*=((\\s*([ab]?)(\\d+)-\\4(\\d+))|(homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
        while( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputNEOGuessSwap) ) {
          if(!Guessmatch.str(3).empty()) {
            // std::cout<<"xsli test guess swap = "<<Guessmatch.str(4)<<Guessmatch.str(5)<<Guessmatch.str(6)<<std::endl;
            auto const freeCQInputNEOGuessSwap1 = std::regex("(SWAP)\\s*=\\s*([ab]?)(\\d+)-\\2(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputNEOGuessSwap1, "");
          }
          if(!Guessmatch.str(7).empty() or !Guessmatch.str(8).empty()) {
            // std::cout<<"xsli test guess swap = HOMO-LUMO"<<std::endl;
            auto const freeCQInputGuessNEOSwap2 = std::regex("(SWAP)\\s*=((homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessNEOSwap2, "");
          }
        }

        auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeDividers, "");
        if(!GuessInputOptions.empty()) CErr("Unrecognized NEO Guess Input Options: "+GuessInputOptions);
        else line = std::regex_replace(line, freeCQInputNEOGuess, "");
      }

    }
#endif
  };

  void SingleSlaterGuessOptions::parseSection(const std::map<std::string,std::string> &dict) {
    if (dict.count("GUESS")) {
      if (dict.at("GUESS") == "CORE") electronicGuess = CoreGuess;
      else if (dict.at("GUESS") == "SAD") electronicGuess = SADGuess;
      else if (dict.at("GUESS") == "TIGHT") electronicGuess = TightGuess;
      else if (dict.at("GUESS") == "RANDOM") electronicGuess = RandomGuess;
      else if (dict.at("GUESS") == "READMO") electronicGuess = ReadBin;
      else if (dict.at("GUESS") == "FCHKMO") electronicGuess = ReadGaussFCHK;
      else if (dict.at("GUESS") == "CLASSICAL") electronicGuess = ClassicalGuess;
    }
    if (dict.count("QP_GUESS")) {
      if (dict.at("QP_GUESS") == "CORE") nuclearGuess = CoreGuess;
      else if (dict.at("QP_GUESS") == "SAD") nuclearGuess = SADGuess;
      else if (dict.at("QP_GUESS") == "TIGHT") nuclearGuess = TightGuess;
      else if (dict.at("QP_GUESS") == "RANDOM") nuclearGuess = RandomGuess;
      else if (dict.at("QP_GUESS") == "READMO") nuclearGuess = ReadBin;
      else if (dict.at("QP_GUESS") == "FCHKMO") nuclearGuess = ReadGaussFCHK;
      else if (dict.at("QP_GUESS") == "CLASSICAL") nuclearGuess = ClassicalGuess;
    }

  }

  void CQInputFile::parseFreeCQInputRT (std::string &line, const std::regex &freeCQInputRT) {

    /***********************************************************************/
    /* RT Input                                                            */
    /* example, RT(MMUT, stepsize= 1.0 as, nsteps = 1000, maxtime = 1.0fs) */
    /***********************************************************************/

    addData("QM/JOB", "RT");
    
    TDSCFOptions tdSCFControls;
    std::smatch RTmatch;


    if (std::regex_search(line, RTmatch, freeCQInputRT)) {
      // str(3) captures the input inside the parentheses
      if (RTmatch.str(3).size() > 0) {

        std::string RTInputOptions = RTmatch.str(3);

        // Check for timestep
        auto const freeCQInputStepsize = std::regex("(STEPSIZE|TIMESTEP)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|(au))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputStepsize) ) {
          tdSCFControls.deltaT = std::stod(RTmatch.str(2));
          if(!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) tdSCFControls.deltaT/=(FSPerAUTime()*1.e3);
          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty()) tdSCFControls.deltaT/=FSPerAUTime();
          // std::cout<<"xsli test read in timestep = "<<tdSCFControls.deltaT<<" au"<<std::endl;
          addData("RT/DELTAT", doubleToString(tdSCFControls.deltaT));
          addData("RT/UNITS", "AU");
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputStepsize, "");
        }

        // Check for maxSteps
        auto const freeCQInputNSteps = std::regex("(NSTEPS|MAXSTEPS|MAXITERATIONS|MAXITER)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputNSteps) ) {
          tdSCFControls.maxSteps = std::stoi(RTmatch.str(2));
          tdSCFControls.tMax = tdSCFControls.maxSteps*tdSCFControls.deltaT;
          // std::cout<<"xsli test read in maxSteps = "<<tdSCFControls.maxSteps<<std::endl;
          addData("RT/TMAX", doubleToString(tdSCFControls.tMax));
          addData("RT/MAXSTEPS", std::to_string(tdSCFControls.maxSteps));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputNSteps, "");
        }

        // Check for maxTime - this will override maxSteps
        auto const freeCQInputMaxTime = std::regex("(MAXTIME)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|(au))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputMaxTime) ) {
          tdSCFControls.tMax = std::stod(RTmatch.str(2));
          if(!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) tdSCFControls.tMax/=(FSPerAUTime()*1.e3);
          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty()) tdSCFControls.tMax/=FSPerAUTime();
          tdSCFControls.maxSteps = (tdSCFControls.tMax + tdSCFControls.deltaT / 4) / tdSCFControls.deltaT;
          // std::cout<<"xsli test read in maxTime = "<<tdSCFControls.tMax<<std::endl;
          addData("RT/TMAX", doubleToString(tdSCFControls.tMax));
          addData("RT/MAXSTEPS", std::to_string(tdSCFControls.maxSteps));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputMaxTime, "");
        }

        // Check for restart algorithm
        auto const freeCQInputRestartAlgorithm = std::regex("(RESTARTALGORITHM|RESTARTALG)\\s*=\\s*((MMUT)|(MODIFIEDMIDPOINT)|(FORWARDEULER)|(EULER)|(EXPLICITMAGNUS2)|(EXPLICITMAGNUSTWO)|(MAGNUS2)|(MAGNUSTWO))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRestartAlgorithm) ) {
//          if (!RTmatch.str(3).empty() or !RTmatch.str(4).empty()) tdSCFControls.restartAlgorithm = RTModifiedMidpoint;
//          else if (!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) tdSCFControls.restartAlgorithm = RTForwardEuler;
//          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty() or !RTmatch.str(9).empty() or !RTmatch.str(10).empty()) tdSCFControls.restartAlgorithm = RTExplicitMagnus2;
          std::string str;
          if (!RTmatch.str(3).empty() or !RTmatch.str(4).empty()) str = "MMUT";
          else if (!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) str = "FORWARDEULER";
          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty() or !RTmatch.str(9).empty() or !RTmatch.str(10).empty()) str = "MAGNUS2";
          //std::cout<<"xsli test read in restart algorithm = "<< std::to_string(tdSCFControls.restartAlgorithm) <<std::endl;
          addData("RT/RESTARTSTEP", str);
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRestartAlgorithm, "");
        }

        // Check for subsystem algorithms before the global algorithm.
        auto const freeCQInputSubsystemAlgorithm = std::regex(
          "([A-Z][A-Z0-9_]*)_INTALG\\s*=\\s*(MMUT|MODIFIEDMIDPOINT|FORWARDEULER|EULER|EXPLICITMAGNUS2|EXPLICITMAGNUSTWO|MAGNUS2|MAGNUSTWO|RK4|RUNGEKUTTAFOURTHORDER|BORT)\\s*([,;:]|$)",
          std::regex_constants::icase);
        while(std::regex_search(RTInputOptions, RTmatch, freeCQInputSubsystemAlgorithm)) {
          std::string label = RTmatch.str(1);
          std::transform(label.begin(), label.end(), label.begin(),
            [](unsigned char c) { return std::toupper(c); });
          addData("RT/" + label + "_INTALG",
            canonicalRealTimeAlgorithm(RTmatch.str(2)));
          RTInputOptions = RTmatch.prefix().str() + RTmatch.suffix().str();
        }

        // Check for RT algorithm, this has to be checked after restart algorithm to avoid conflict
        auto const freeCQInputRTAlgorithm = std::regex("(MMUT|MODIFIEDMIDPOINT|FORWARDEULER|EULER|EXPLICITMAGNUS2|EXPLICITMAGNUSTWO|MAGNUS2|MAGNUSTWO|RK4|RUNGEKUTTAFOURTHORDER|BORT)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRTAlgorithm) ) {
          addData("RT/INTALG", canonicalRealTimeAlgorithm(RTmatch.str(1)));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRTAlgorithm, "");
        }

        // Check for the frequency of save
        auto const freeCQInputISave = std::regex("(ISAVE|AUTOSAVE)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputISave) ) {
          tdSCFControls.iSave = std::stoi(RTmatch.str(2));
          // std::cout<<"xsli test read in iSave = "<<tdSCFControls.iSave<<std::endl;
          addData("RT/SAVESTEP", std::to_string(tdSCFControls.iSave));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputISave, "");
        }

        // Check for the frequency of cubegen
        auto const freeCQInputICube = std::regex("(ICUBE|CUBEGEN)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputICube) ) {
          tdSCFControls.iCube = std::stoi(RTmatch.str(2));
          // std::cout<<"xsli test read in iCube = "<<tdSCFControls.iCube<<std::endl;
          addData("RT/SAVECUBE", std::to_string(tdSCFControls.iCube));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputICube, "");
        }

        // Check for the frequency of print
        auto const freeCQInputIPrint = std::regex("(IPRINT|AUTOPRINT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputIPrint) ) {
          tdSCFControls.iPrint = std::stoi(RTmatch.str(2));
          // std::cout<<"xsli test read in iPrint = "<<tdSCFControls.iPrint<<std::endl;
          addData("RT/PRINTSTEP", std::to_string(tdSCFControls.iPrint));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputIPrint, "");
        }

        // Check for the frequency of automatic restart
        auto const freeCQInputIRestart = std::regex("(IRESTART|AUTORESTART)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputIRestart) ) {
          tdSCFControls.iRestart = std::stoi(RTmatch.str(2));
          // std::cout<<"xsli test read in iRestart = "<< std::to_string(tdSCFControls.iRestart) <<std::endl;
          addData("RT/IRSTRT", std::to_string(tdSCFControls.iRestart));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputIRestart, "");
        }

        // Check for restart, this has to be checked after iRESTART to avoid conflict
        auto const freeCQInputRestart = std::regex("RESTART(\\s*=\\s*(\\d+))?\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRestart) ) {
          tdSCFControls.restoreFromStep = -1;
          if(!RTmatch.str(2).empty()) tdSCFControls.restoreFromStep = std::stoi(RTmatch.str(2));
          // std::cout<<"xsli test read in do restart = "<<tdSCFControls.restoreFromStep<<std::endl;
          addData("RT/RESTARTFROM", std::to_string(tdSCFControls.restoreFromStep));
          addData("RT/RESTART", "TRUE");
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRestart, "");
        }

        auto const freeCQInputrtGaunt = std::regex("(RTGAUNT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtGaunt) ) {
          tdSCFControls.rtGaunt = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtGaunt = "<<tdSCFControls.rtGaunt<<std::endl;
          addData("RT/RTGAUNT", std::to_string(tdSCFControls.rtGaunt));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtGaunt, "");
        }

        auto const freeCQInputrtGauge = std::regex("(RTGAUGE)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtGauge) ) {
          tdSCFControls.rtGauge = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtGauge = "<<tdSCFControls.rtGauge<<std::endl;
          addData("RT/RTGAUGE", std::to_string(tdSCFControls.rtGauge));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtGauge, "");
        }

        auto const freeCQInputrtBreit = std::regex("(RTBREIT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtBreit) ) {
          tdSCFControls.rtBreit = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtBreit = "<<tdSCFControls.rtBreit<<std::endl;
          addData("RT/RTBREIT", std::to_string(tdSCFControls.rtBreit));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtBreit, "");
        }

        auto const freeCQInputRtprintden = std::regex("(RTPRINTDEN)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRtprintden) ) {
          tdSCFControls.Rtprintden = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in Rtprintden = "<<tdSCFControls.Rtprintden<<std::endl;
          addData("RT/RTPRINTDEN", std::to_string(tdSCFControls.Rtprintden));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRtprintden, "");
        }

        // Check for BORTPrinting options
        auto const freeCQInputBORTPrint = std::regex("BORTPRINT(\\s*=\\s*(\\d+))?\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputBORTPrint) ) {
          if(!RTmatch.str(2).empty()) tdSCFControls.BORTPrintLevel = std::stoi(RTmatch.str(2));
          addData("RT/BORTPRINTLEVEL", std::to_string(tdSCFControls.BORTPrintLevel));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputBORTPrint, "");
        }

        // Check for BORT SCF Convergence accuracy
        auto const freeCQInputBORTAccuracy = std::regex("(BORTACCURACY)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputBORTAccuracy) ) {
          tdSCFControls.BORTAccuracy = std::stod(RTmatch.str(2));
          addData("RT/BORTACCURACY", doubleToString(tdSCFControls.BORTAccuracy));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputBORTAccuracy, "");
        }

        auto const freeCQInputOrbitalPopFreq = std::regex("(ORBITALPOPFREQ)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputOrbitalPopFreq) ) {
          tdSCFControls.orbitalPopFreq = std::stoi(RTmatch.str(2));
          std::cout<<"test read in OrbitalPopFreq = "<<tdSCFControls.orbitalPopFreq<<std::endl;
          addData("RT/ORBITALPOPFREQ", std::to_string(tdSCFControls.orbitalPopFreq));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputOrbitalPopFreq, "");
        }

        auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
        RTInputOptions = std::regex_replace(RTInputOptions, freeDividers, "");
        if(!RTInputOptions.empty()) CErr("Unrecognized RT Input Options: "+RTInputOptions);
        else line = std::regex_replace(line, freeCQInputRT, "");

        

      }
    };
  }

  void TDSCFOptions::parseSection(const std::map<std::string,std::string> &dict) {
    if (dict.count("DELTAT")) deltaT = std::stod(dict.at("DELTAT"));
    if (dict.count("TMAX")) tMax = std::stod(dict.at("TMAX"));
    if (dict.count("MAXSTEPS")) maxSteps = std::stoi(dict.at("MAXSTEPS"));
    // Initialize maxSteps if not set
    if (maxSteps==0 and (tMax!=0. and deltaT != 0.)) maxSteps = (tMax + deltaT/4) / deltaT;
    if (dict.count("IRSTRT")) iRestart = std::stoi(dict.at("IRSTRT"));
    if (dict.count("SAVESTEP")) iSave = std::stoi(dict.at("SAVESTEP"));
    if (dict.count("SAVECUBE")) iCube = std::stoi(dict.at("SAVECUBE"));
    if (dict.count("PRINTSTEP")) iPrint = std::stoi(dict.at("PRINTSTEP"));
    if (dict.count("RESTARTFROM")) restoreFromStep = std::stoi(dict.at("RESTARTFROM"));
    if (dict.count("INTALG")) integrationAlgorithm = parseRealTimeAlgorithm(dict.at("INTALG"));
    const std::string suffix = "_INTALG";
    for(const auto& [key, value] : dict) {
      if(key.size() <= suffix.size() or key.compare(key.size() - suffix.size(), suffix.size(), suffix) != 0)
        continue;
      auto label = key.substr(0, key.size() - suffix.size());
      subsystemIntegrationAlgorithms[label] = parseRealTimeAlgorithm(value);
    }
    if (dict.count("BORTPRINTLEVEL")) BORTPrintLevel = std::stoi(dict.at("BORTPRINTLEVEL"));
    if (dict.count("BORTACCURACY")) BORTAccuracy = std::stod(dict.at("BORTACCURACY"));
    if (dict.count("RESTARTSTEP")) {
      if (dict.at("RESTARTSTEP") == "MMUT") restartAlgorithm = RestartAlgorithm::ModifiedMidpoint;
      else if (dict.at("RESTARTSTEP") == "FORWARDEULER") restartAlgorithm = RestartAlgorithm::ForwardEuler;
      else if (dict.at("RESTARTSTEP") == "MAGNUS2") restartAlgorithm = RestartAlgorithm::ExplicitMagnus2;
    }
    if (dict.count("RTGAUNT")) rtGaunt = std::stoi(dict.at("RTGAUNT"));
    if (dict.count("RTPRINTDEN")) Rtprintden = std::stoi(dict.at("RTPRINTDEN"));
    if (dict.count("ORBITALPOPFREQ")) orbitalPopFreq = std::stoi(dict.at("ORBITALPOPFREQ"));
    if (dict.count("RTGAUGE")) rtGauge = std::stoi(dict.at("RTGAUGE"));
    if (dict.count("RTBREIT")){ 
      rtBreit = std::stoi(dict.at("RTBREIT"));
      rtGauge = rtBreit;
      rtGaunt = rtBreit;
    }
    if (dict.count("SAVEONEPDM")){ 
      saveOnePDM = true;
    }

  }

  void CQInputFile::parseFreeCQInputField (std::string &line) {

    /**************************************************************************/
    /* Field Input                                                            */
    /* example, Field(delta, start= 1.0 as, end = 10 fs, amplitude = 10 au)   */
    /**************************************************************************/
    auto const freeCQInputField = std::regex("FIELD(\\((.*?)\\))?", std::regex_constants::icase);
    std::smatch Fieldmatch;

    if (std::regex_search(line, Fieldmatch, freeCQInputField)) {
      // str(3) captures the input inside the parentheses
      if (Fieldmatch.str(2).size() > 0) {
        std::string FieldInputOptions = Fieldmatch.str(2);
        // std::cout << "xsli test Field: " << FieldInputOptions << std::endl;

        // step function
        auto const freeCQInputFieldStep = std::regex("(DELTA|STEP|CONSTANT)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldStep) ) {
          // std::cout<<"xsli test read in field type = "<<Fieldmatch.str(0)<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldStep, "");
        }

        // ton or start time
        auto const freeCQInputFieldStart = std::regex("(TON|TIMEON|TSTART|START|STARTTIME)\\s*=\\s*(((-?\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|(au)))|(ALWAYS))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldStart) ) {
          double starttime;
          if (!Fieldmatch.str(11).empty()) starttime = -1.0;
          else {
            starttime = std::stod(Fieldmatch.str(4));
            if (starttime < 0) starttime = -1.0;
            else {
              if (!Fieldmatch.str(7).empty() or !Fieldmatch.str(8).empty()) starttime /= (FSPerAUTime() * 1.e3);
              else if (!Fieldmatch.str(9).empty() or !Fieldmatch.str(10).empty()) starttime /= FSPerAUTime();
            }
          }
          // std::cout<<"xsli test read in tOn = "<<starttime<<" au"<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldStart, "");
        }

        // toff or end time
        auto const freeCQInputFieldEnd = std::regex("(TOFF|TIMEOFF|TEND|END|ENDTIME)\\s*=\\s*(((-?\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|au))|(NEVER))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldEnd) ) {
          double endtime;
          if (!Fieldmatch.str(11).empty()) endtime = -1.0;
          else {
            endtime = std::stod(Fieldmatch.str(4));
            if(endtime<0) endtime = -1.0;
            else {
              if(!Fieldmatch.str(7).empty() or !Fieldmatch.str(8).empty()) endtime/=(FSPerAUTime()*1.e3);
              else if (!Fieldmatch.str(9).empty() or !Fieldmatch.str(10).empty()) endtime/=FSPerAUTime();
            }
          }
          // std::cout<<"xsli test read in tOff = "<<endtime<<" au"<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldEnd, "");
        }


        // field amplitude
        auto const freeCQInputFieldAmp = std::regex("(AMP|AMPLITUDE)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((EV)|(MEV)|(AU))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldAmp) ) {
          double amp;
          amp = std::stod(Fieldmatch.str(2));
          if(!Fieldmatch.str(5).empty()) amp/=EVPerHartree();
          else if (!Fieldmatch.str(6).empty()) amp/=(EVPerHartree()*1.e3);
          // std::cout<<"xsli test read in amp = "<<amp<<" au"<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldAmp, "");
        }

        auto const freeDividers = std::regex("\\s+|,+", std::regex_constants::icase);
        FieldInputOptions = std::regex_replace(FieldInputOptions, freeDividers, "");
        if (!FieldInputOptions.empty()) CErr("Unrecognized Field Input Options: " + FieldInputOptions);
        else line = std::regex_replace(line, freeCQInputField, "");
      }
    }
  }//parseFreeCQInputField

  void CQInputFile::parseFreeCQInputCC (std::string &line, const std::regex &freeCQInputCI){
    /****************************************************************************************/
    /* CC Input                                                                           */
    /* example, CCSDT(skipscf,maxiter=1000,etol=1d-6)                                     */
    /* example, X2C-CCSD                                      *                           */
    /****************************************************************************************/
    
  } //parseFreeCQInputCC

  void CQInputFile::parseFreeCQInputCI (std::string &line, const std::regex &freeCQInputCI){

    /****************************************************************************************/
    /* CI Input                                                                             */
    /* example, CASSCF(10o,5e,accuracy=1.e-4,nstates=5)                                     */
    /* example, X2C-DASCI(40o,20e,nDAS=5,maxexcitation=2)                                   */
    /* example, 4C-RASCI(40o,20e,RAS1(10o,10e,1h),RAS2(20o,10e),RAS3(10o,0e,2p),nstates=10) */
    /****************************************************************************************/

    //"((2C)|(X2C)|(4C)|(G))?(-)?((CASSCF)|(DASSCF)|(CASCI)|(DASCI)|(CI))"
    //            1                                        2
    
    addData("QM/JOB", "CI");
    
    std::smatch CIMatch;
    std::regex_search(line, CIMatch, freeCQInputCI);

    //Non-switched defaults
    addData("CI/PRINTRDMS","2");
    //TODO:addData("CI/PRINTDETOCC","TRUE"); //needs to be put into dev first

    // CIMATCH captures the CI type (e.g., X2C, 4C)
    // auto const meth = std::regex("(SF)?-?(R|U|X2C|G)");
    //if (std::regex_search(line, CIMatch, meth) {
    //  //TODO G/X2C/etc.
    //}

    bool doSCF = false;
    //Redefine maximum number of casscf iterations:
    std::regex nSCFITR("(nscf|ncasscf|nscfiter|ncasscfiter|maxscfiter)\\=([0-9]*)",std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, nSCFITR) ) {
      std::cout << "Max. Number CASSCF Iterations: " + CIMatch.str(2) << std::endl;
      addData("CI/MAXSCFITER",CIMatch.str(2));
      doSCF = true;
    }//Max SCF Iter

    // read in CAS algo
    auto const freeCQInputCASAlgo = std::regex("(scfalg|casalg)\\=([a-z]*)", std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputCASAlgo) ) {
      addData("CI/SCFALG",CIMatch.str(2));
    }//cas alg.

    //state avgeraging (defaults to true)
    auto const freeCQInputSACAS = std::regex("(stateaverage)\\=(true|false)", std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputSACAS) ) {
      addData("CI/STATEAVERAGE",CIMatch.str(2));
    }//bool stavg

    auto const freeCQInputCASSCF = std::regex("scf",std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputCASSCF) or doSCF ) {
      doSCF = true;
      std::cout << "JOBTYPE=SCF" << std::endl;
      addData("CI/JOBTYPE","DASSCF");
    } else {
      std::cout << "JOBTYPE=CI" << std::endl;
      addData("CI/JOBTYPE","DASCI");
    }//parse type of ci

    //Determine active space:
    int nelec = 0;
    int norbt = 0;

    std::regex nEle("\\b([0-9]+)(e|E)\\b"); // an integer number followed by "e" - number of electrons
    std::regex_search(line, CIMatch, nEle);
    std::cout<<"      nElectrons = "+CIMatch.str(1)<< std::endl;
    addData("CI/NACTELEC",CIMatch.str(1));
    nelec = std::stoi(CIMatch.str(1));

    std::regex nOrb("\\b([0-9]+)(o|O)\\b"); // an integer number followed by "o" - number of orbitals
    std::regex_search(line, CIMatch, nOrb);
    std::cout<<"      nOrbitals = "+CIMatch.str(1)<< std::endl;
    addData("CI/NACTORB",CIMatch.str(1));
    norbt = std::stoi(CIMatch.str(1));

    // read in CI eigensolver accuracy
    auto const freeCQInputCIAccuracy = std::regex("accuracy\\s*=\\s*((\\d+\\.?\\d*|\\.\\d+)(e[-+]?\\d+)?)\\s*([,;:]|$)", std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputCIAccuracy) ) {
      addData("CI/CICONV",CIMatch.str(1));
    }//solver acc

    // read in CI eigensolver algo
    auto const freeCQInputCIAlgo = std::regex("(alg|cialg|cidiagalg)\\=([a-z]*)", std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputCIAlgo) ) {
      addData("CI/CIDIAGALG",CIMatch.str(2));
    } else {
      addData("CI/CIDIAGALG","DAVIDSON");
      addData("CI/CISIGMA2EALG","KNOWLESHANDY");
    }//solver algo

    // read in CI eigensolver algo (david above defaults kh, reg. default naive
    auto const freeCQInputCI2eAlgo = std::regex("(cisigma2ealg)\\=([a-z]*)", std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, freeCQInputCI2eAlgo) ) {
      addData("CI/CISIGMA2EALG",CIMatch.str(2));
    }//sigma2ealg

    // nstates = an integer number == number of eigenstates to solve for
    int nroots = 3;
    std::regex nStates("(nstates|nroots)\\=([0-9]*)",std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, nStates) ) {
      nroots = std::stoi(CIMatch.str(2));
    }//Number roots
    //std::cout<<"      Num. Roots = " + std::to_string(nroots) << std::endl;
    addData("CI/NROOTS", std::to_string(nroots));

    // nosc = an integer number == number of eigenstates to calculate f for
    // TODO Raj's osc. order may mess this up
    int nosc = 1;
    std::regex nOsc("(nosc|oscistren)\\=([0-9]*)",std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, nOsc) ) {
      nosc = std::stoi(CIMatch.str(2));
    }//Number Oscistren
    addData("CI/OSCISTREN", "TRUE");
    addData("CI/OSCISTREN_INITSTATES", std::to_string(std::min(nosc,nroots)));

    // nDAS = an integer number == number of DAS TODO self-defined das line?
    int nDAS = 1;
    std::string dasStream = "";
    std::regex nDASspace("ndas\\=([0-9]*)",std::regex_constants::icase);
    if ( std::regex_search(line, CIMatch, nDASspace) ) {
      nDAS = std::stoi(CIMatch.str(1));
      if ( nDAS > norbt ) {
        CErr("You've declared more DAS spaces than you have active orbitals, pull yourself together.");
      }
      // auto: set up das spaces
      int distElec = 0;
      int distOrb  = 0;
      int distDAS  = 0;
      //In case nOrbt doesn't go evenly into nSpace:
      for ( int i = 0; i < norbt % nDAS; i++ ) {
        int orbPerSpace = int(double(norbt)/double(nDAS)+0.5);
        if ( 0 < nelec - orbPerSpace ) {
          distElec = orbPerSpace;
          nelec -= orbPerSpace;
        } else if ( 0 < nelec ) {
          distElec = nelec;
          nelec -= nelec;
        } else {
          distElec = 0;
          nelec = 0;
        }//Reducing number electrons into spaces
        distOrb += orbPerSpace;
        distDAS += 1;
        dasStream += "{[" + std::to_string(orbPerSpace) + "o," \
            + std::to_string(distElec) + "e,\"Space " + std::to_string(distDAS) + "\"]}";
      }//non-equal das spaces
      //remainder and/or if norbt % ndas == 0
      for ( int i = 0; i < nDAS - distDAS; i++) {
        int orbPerSpace = (norbt-distOrb)/(nDAS-distDAS);
        if ( 0 < nelec - orbPerSpace ) {
          distElec = orbPerSpace;
          nelec -= orbPerSpace;
        } else if ( 0 < nelec ) {
          distElec = nelec;
          nelec -= nelec;
        } else {
          distElec = 0;
          nelec = 0;
        }//elect check
        dasStream += "{[" + std::to_string(orbPerSpace) \
            + "o," + std::to_string(distElec) + "e,\"Space " + std::to_string(distDAS + i) + "\"]}";
        }//Remaining das spaces
    }//DAS

    std::cout << "      Splitting into " << nDAS << " DAS spaces" << std::endl;
    if ( nDAS > 1 ) {
      //stream printing
      std::cout << "      DAS STREAM: " << dasStream << std::endl;
    }//stream print
    addData("CI/DAS",std::to_string(nDAS) + "\n" + dasStream);

  };//parseFreeCQInputCI
}; // namespace ChronusQ

