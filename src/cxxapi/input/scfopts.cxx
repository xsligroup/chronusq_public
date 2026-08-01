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

namespace ChronusQ {

  namespace {

    SS_GUESS parseSCFGuess(const std::string& value, const std::string& keyword) {
      if(value == "CORE") return CORE;
      if(value == "SAD") return SAD;
      if(value == "TIGHT") return TIGHT;
      if(value == "RANDOM") return RANDOM;
      if(value == "READMO") return READMO;
      if(value == "READDEN") return READDEN;
      if(value == "SCF") return SCF;
      if(value == "FCHKMO") return FCHKMO;
      if(value == "CLASSICAL") return NEOConvergeClassical;
      CErr("Unrecognized entry for " + keyword);
      return SAD;
    }

  }


  std::set<std::string> CQSCF_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "ENETOL",
      "DENTOL",
      "FDCTOL",
      "MAXITER",
      "INCFOCK",
      "NINCFOCK",
      "GUESS",
      "ALG",
      "EXTRAP",
      "DIIS",
      "DIISALG",
      "NKEEP",
      "DAMP",
      "DAMPPARAM",
      "DAMPERROR",
      "FIELD",
      "PRINTMOS",
      "NEO",
      "PROT_GUESS",
      "SWAPMO",
      "SWITCH",
      "NRAPPROX",
      "NRTRUST",
      "NRLEVELSHIFT",
      "PRINTCONTRACTIONTIMING" ,
      "ACCURACY",
      "CUBE",
      "ORBPROP",
      "NEOOPTIMIZEONLY",
      "NEOOPTIMIZEFIRST",
      "NEOSTEPWISEOPTIMIZE",
      "REMOVELINEARDEP",
      "REMOVEATOMICLINEARDEPONLY",
      "LINEARDEPTOL",
      "DENMODIFIER"
    };

    // <LABEL>_GUESS for Multiparticle SCF guesses
    // Accept all <LABEL>_GUESS here but will error out in resolveSubsystemGuessOptions for invalid "LABEL"
    const std::string suffix = "_GUESS";
    for(const auto& item : inputSection) {
      const auto& key = item.first;
      if(key.size() > suffix.size() and
         key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0)
        allowedKeywords.insert(key);
    }

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }


  void HandlePostSCFRestarts(std::ostream &out, CQInputFile &input,
                             SCFControls &scfControls) {

    bool restartRT = false;
    bool restartMD = false;
    
    if ( input.containsSection("RT") )
      OPTOPT( restartRT |= input.getData<bool>("RT/RESTART"); )
    
    if ( input.containsSection("DYNAMICS") ){
      std::string restart = "FALSE";
      OPTOPT( restart = input.getData<std::string>("DYNAMICS/RESTART");)
      trim(restart);
      restartMD = (restart.compare("FALSE") != 0);
    }
    
    // Add additional checks like the one below for restarting other post SCF
    // Add additional checks like the one below for restarting other post SCF
    // RESP restart is NYI, so this is a commented out placeholder example
    //
    // else if ( input.containsSection("RESP") )
    //   OPTOPT( restart |= input.getData<bool>("RESP.RESTART"); )

    // Skip SCF is the we are doing a restart for RT/MD
    if ( restartRT or restartMD ) {
      // Since scfControls currently is global, for BOMD restart we can't skip SCF
      // TODO: allow BOMD job to have its own scfControls  
      if (restartMD and parseJob(input.getData<std::string>("QM/JOB"))==JobType::BOMD )
        return;
      out << "  *** RESTART requested -- SCF/GUESS set to READMO and SCF/ALG set to SKIP ***";
      out << std::endl;

      scfControls.guess = READMO;
      //scfControls.scfAlg = _SKIP_SCF;
      scfControls.scfAlg = _CONVENTIONAL_SCF;
      scfControls.energyOnly = true;
    }

  }

  std::unordered_map<std::string,int> SpinMap = {
    { "A" , 0  },
    { "B" , 1  }
  };

  void HandleOrbitalSwaps(std::ostream &out, CQInputFile &input,
    SingleSlaterBase &ss, std::string prefix) {

    // MO swapping
    std::string swapMOStrings;
    std::string key = "SCF/" + prefix + "SWAPMO";
    OPTOPT( swapMOStrings = input.getData<std::string>(key));
    if ( not swapMOStrings.empty() ) {
      std::cout << "  * Manually MO Swapping Detected: " << std::endl;

      if( ss.scfControls.guess != READMO and ss.scfControls.guess != FCHKMO )
        CErr("MO swapping only for user-specified guess MOs");

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

        // mo[1] only present for unrestricted calcs
        if( spinDir == "B" and not (ss.nC == 1 and not ss.iCS) ) CErr("Swapping of beta MOs is only valid for open-shell 1c");

        ss.moPairs[SpinMap[spinDir]].emplace_back(std::stoul(moTokens[0]), std::stoul(moTokens[1]));
      }

    }

  }

  SCFControls CQSCFOptions(std::ostream &out, CQInputFile &input, EMPerturbation &pert) {

    SCFControls scfControls;

    // SCF section not required
    if( not input.containsSection("SCF") ) {
     
      // Restart jobs
      HandlePostSCFRestarts(out, input, scfControls);

      return scfControls;
    }

    // Optionally parse guess

    OPTOPT( scfControls.rmsdPConvTol =
      input.getData<double>("SCF/ACCURACY"); )

    scfControls.maxdPConvTol = scfControls.rmsdPConvTol*100;
    scfControls.eneConvTol = scfControls.rmsdPConvTol*100;

    // Energy convergence tolerance
    //OPTOPT( scfControls.eneConvTol =
    //          input.getData<double>("SCF/ENETOL"); )

    // Energy convergence tolerance
    //OPTOPT( scfControls.denConvTol =
    //          input.getData<double>("SCF/DENTOL"); )

    // Energy Gradient convergence tolerance
    //OPTOPT( scfControls.FDCConvTol =
    //          input.getData<double>("SCF/FDCTOL"); )

    // Maximum SCF iterations
    OPTOPT( scfControls.maxSCFIter =
              input.getData<size_t>("SCF/MAXITER"); )


    // Incremental Fock Options
    OPTOPT(
      scfControls.doIncFock = input.getData<bool>("SCF/INCFOCK");
    )
    OPTOPT(
      scfControls.nIncFock = input.getData<size_t>("SCF/NINCFOCK");
    )


    // Guess
    std::string guessString = "SAD";
    OPTOPT( guessString = input.getData<std::string>("SCF/GUESS"); )
    trim(guessString);
    scfControls.guess = parseSCFGuess(guessString, "SCF/GUESS");
    

    // Proton Guess For Legacy NEO Calculations
    const bool hasProtGuess = input.containsData("SCF/PROT_GUESS");
    SS_GUESS protGuess = TIGHT;
    if(hasProtGuess) {
      std::string protGuessString = input.getData<std::string>("SCF/PROT_GUESS");
      trim(protGuessString);
      protGuess = parseSCFGuess(protGuessString, "SCF/PROT_GUESS");
    }

    // Get all LABEL_GUESS keywords specified by the user under SCF section
    // And record them in subsystemGuesses
    const std::string guessSuffix = "_GUESS";
    for(const auto& [key, value] : input.getSection("SCF")) {
      if(key == "GUESS" or key == "PROT_GUESS" or
         key.size() <= guessSuffix.size() or
         key.compare(key.size() - guessSuffix.size(), guessSuffix.size(), guessSuffix) != 0)
        continue;
      std::string label = key.substr(0, key.size() - guessSuffix.size());
      if(label.empty()) continue;
      std::string selection = value;
      trim(selection);
      scfControls.subsystemGuesses[label] = parseSCFGuess(selection, "SCF/" + key);
    }

    // Handle legacy NEO input keyword PROT_GUESS (Transfer to QP_GUESS)
    if(hasProtGuess) {
      if(scfControls.subsystemGuesses.count("QP"))
        CErr("Cannot set both SCF/PROT_GUESS and SCF/QP_GUESS");
      scfControls.subsystemGuesses["QP"] = protGuess;
      scfControls.prot_guess = protGuess == TIGHT ? NEOTightParticle : protGuess;
    }
    
    std::string neoonlyoptstring;
    OPTOPT( neoonlyoptstring = input.getData<std::string>("SCF/NEOOPTIMIZEONLY"); )
    trim(neoonlyoptstring);
    if(!neoonlyoptstring.empty())
    {
      if ( ! neoonlyoptstring.compare("ELECTRONIC"))
        scfControls.NEOSubSystemOpt.push_back("E");   // "Electronic" subsystem is now handled as "E" in the MultiParticleSS
      else if ( ! neoonlyoptstring.compare("PROTONIC") or ! neoonlyoptstring.compare("QP"))
        scfControls.NEOSubSystemOpt.push_back("QP");  // "Protonic" subsystem is now handled as "QP" in the MultiParticleSS
      else
        CErr("Unrecognized entry for SCF/NEOOPTIMIZEONLY");
    }

    OPTOPT(scfControls.NEOStepwiseOpt = input.getData<bool>("SCF/NEOSTEPWISEOPTIMIZE");)
    if(scfControls.NEOStepwiseOpt)
    {
      if(scfControls.NEOSubSystemOpt.size())
        CErr("Cannot set both NEOOPTIMIZEONLY and NEOSTEPWISEOPTIMIZE");
      std::cout << "NEO Stepwise optimization requested" << std::endl;
      std::string neooptfirststring;
      OPTOPT( neooptfirststring = input.getData<std::string>("SCF/NEOOPTIMIZEFIRST"); )
      trim(neooptfirststring);
      if(!neooptfirststring.empty())
      {
        if ( ! neooptfirststring.compare("ELECTRONIC"))
          scfControls.NEOSubSystemOpt.push_back("E");
        else if ( ! neooptfirststring.compare("PROTONIC") or ! neooptfirststring.compare("QP"))
          scfControls.NEOSubSystemOpt.push_back("QP");
        else
          CErr("Unrecognized entry for SCF/NEOOPTIMIZEFIRST");
      }
      else
      {
        scfControls.NEOSubSystemOpt.push_back("QP");
      }
      std::cout << "NEO Stepwise optimization will begin with the " << scfControls.NEOSubSystemOpt[0] << " subsystem" << std::endl;
    }
 
    // ALGORITHM
    std::string algString = "CONVENTIONAL";
    OPTOPT( algString = input.getData<std::string>("SCF/ALG"); )
    if( not algString.compare("CONVENTIONAL") )
      scfControls.scfAlg = _CONVENTIONAL_SCF;
    else if( not algString.compare("NR") )
      scfControls.scfAlg = _NEWTON_RAPHSON_SCF;
    else if( not algString.compare("SKIP") ) {
      scfControls.scfAlg = _CONVENTIONAL_SCF;
      scfControls.energyOnly = true;
    } else 
      CErr("Unrecognized entry for SCF/ALG!");


    // Newton-Raphson SCF Approximation
    std::string nrAlgString = "BFGS";
    OPTOPT( nrAlgString = input.getData<std::string>("SCF/NRAPPROX"); )
    if( not nrAlgString.compare("FULL") )
      scfControls.nrAlg = FULL_NR;
    else if( not nrAlgString.compare("BFGS") )
      scfControls.nrAlg = QUASI_BFGS;
    else if( not nrAlgString.compare("SR1") )
      scfControls.nrAlg = QUASI_SR1;
    else if( not nrAlgString.compare("GRADDESCENT") )
      scfControls.nrAlg = GRAD_DESCENT;
    else 
      CErr("Unrecognized entry for SCF/NRAPPROX");


    // Newton-Raphson SCF Initial trust region and level-shift
    OPTOPT(
      scfControls.nrTrust = input.getData<double>("SCF/NRTRUST");
    )
    OPTOPT(
      scfControls.nrLevelShift = input.getData<double>("SCF/NRLEVELSHIFT");
    )

    // For modified SCF Procedures
    std::string rdmbuilderstr = "";
    OPTOPT( rdmbuilderstr = input.getData<std::string>("SCF/DENMODIFIER"));
    if(!rdmbuilderstr.compare("MOM"))
    {
      if(scfControls.guess != READMO) CErr("MOM requires a guess set of orbitals using READMO");
      scfControls.rdmBuilderType = RDM_BUILDER_TYPE::MOM;
    }
    else if(!rdmbuilderstr.compare("AUFBAU") || rdmbuilderstr.empty())
    {
      scfControls.rdmBuilderType = RDM_BUILDER_TYPE::AUFBAU;
    }
    else
    {
      CErr("Unrecognized option for SCF/DENMODIFIER");
    }

    std::string protrdmbuilderstr = "";
    OPTOPT( protrdmbuilderstr = input.getData<std::string>("SCF/PROT_DENMODIFIER"));
    if(!protrdmbuilderstr.compare("NEOSTATEAVERAGE"))
    {
      scfControls.protrdmBuilderType = RDM_BUILDER_TYPE::NEOSTATEAVERAGE;
      size_t NStates = 0;
      OPTOPT(NStates = input.getData<size_t>("SCF/NEOSTATEAVERAGESTATES"));
      if(!NStates)
        CErr("Requesting a NEO-SCF-StateAveraged Calcualtion requires NEOSTATEAVERAGESTATES > 0!");
      scfControls.NEOStateAverageNStates = NStates;
      //OPTOPT(scfControls.NEOStateAverageNStates = input.getData<size_t>("SCF/NEOSTATEAVERAGESTATES"));
    }
    else if(!protrdmbuilderstr.compare("MOM"))
    {
      if(scfControls.prot_guess != READMO) CErr("Proton MOM requires a guess set of orbitals using READMO");
      scfControls.protrdmBuilderType = RDM_BUILDER_TYPE::MOM;
    }
    else if(!protrdmbuilderstr.compare("NEOFINITETEMP"))
    {
      scfControls.protrdmBuilderType = RDM_BUILDER_TYPE::FINITETEMP;
      OPTOPT(scfControls.finitetemp = input.getData<double>("SCF/NEOFINITET");)
      CErr("NEO Finite Temp NYI!");
    }
    else if(!protrdmbuilderstr.compare("AUFBAU") || protrdmbuilderstr.empty())
    {
      scfControls.protrdmBuilderType = RDM_BUILDER_TYPE::AUFBAU;
    }
    else
    {
      CErr("Unrecognized option for SCF/PROT_DENMODIFIER");
    }


    // Restart jobs
    HandlePostSCFRestarts(out, input, scfControls);


    // Toggle extrapolation in its entireity
    OPTOPT(
      scfControls.doExtrap =
        input.getData<bool>("SCF/EXTRAP");
    )

    // Handle DIIS options
    std::string diisAlgString = "CDIIS"; 
    OPTOPT( diisAlgString = input.getData<std::string>("SCF/DIISALG"); )
    if( not diisAlgString.compare("CEDIIS"))
      scfControls.diisAlg = CEDIIS;
    else if( not diisAlgString.compare("CDIIS"))
      scfControls.diisAlg = CDIIS;
    else if( not diisAlgString.compare("EDIIS"))
      scfControls.diisAlg = EDIIS;
    else
        CErr("Unrecognized entry for SCF/DIISALG!");


    // Check if it specifically says DIIS=FALSE
    OPTOPT(
      bool doDIIS = input.getData<bool>("SCF/DIIS");
      if( not doDIIS ) 
        scfControls.diisAlg = NONE;
    );

    // Number of terms for keep for DIIS
    OPTOPT( scfControls.nKeep = input.getData<size_t>("SCF/NKEEP"); )
    if( scfControls.nKeep < 1 )
      scfControls.nKeep = 1;

    // Point at which to switch from EDIIS to CDIIS
    OPTOPT( scfControls.cediisSwitch = input.getData<double>("SCF/SWITCH"); )
    if( scfControls.cediisSwitch <= 0. )
      CErr("CEDIIS Switch is less than or equal to zero");

    // Parse Damping options
    OPTOPT(
      scfControls.doDamp = input.getData<bool>("SCF/DAMP");
    );

    OPTOPT(
      scfControls.dampStartParam =
        input.getData<double>("SCF/DAMPPARAM");
    );

    OPTOPT(
      scfControls.dampError =
        input.getData<double>("SCF/DAMPERROR");
    );


    // SCF Field
    std::string fieldStr;
    OPTOPT(
      fieldStr = input.getData<std::string>("SCF/FIELD");
    )
    EMPerturbation parsedField;
    if (!fieldStr.empty())
        handleField(fieldStr, parsedField);
    pert.addField(parsedField);

    // Printing Options
    if ( input.containsData("SCF/PRINTMOS") ) {
      try { scfControls.printMOCoeffs = input.getData<size_t>("SCF/PRINTMOS"); }
      catch(...) {
        CErr("Invalid PRINTMOS input. Please use number 0 ~ 9.");
      }
    }
    if (scfControls.printMOCoeffs >= 10 ) CErr("SCF print level is not valid!");


    // Parse whether to print contraction timing during SCF
    OPTOPT( scfControls.printContractionTiming = input.getData<bool>("SCF/PRINTCONTRACTIONTIMING"); )

    // Linear dependency tolerance
    if ( input.containsData("SCF/REMOVELINEARDEP")) {
      scfControls.rmLinearDep = input.getData<bool>("SCF/REMOVELINEARDEP");
    }
    if ( input.containsData("SCF/REMOVEATOMICLINEARDEPONLY")) {
      scfControls.rmAtomicLinDepOnly = input.getData<bool>("SCF/REMOVEATOMICLINEARDEPONLY");
    }
    if ( input.containsData("SCF/LINEARDEPTOL")) {
      scfControls.linearDepTol = input.getData<double>("SCF/LINEARDEPTOL");
    }





    // Handling eqivalences in the input options


    // Setting the damp param to 0. is equivalent to
    // turning damping off
    if( scfControls.dampStartParam == 0. )
      scfControls.doDamp = false;

    // Turning off both damping and DIIS is equivalent
    // to turning off extrapolation entirely
    if( not scfControls.doDamp and
        scfControls.diisAlg == NONE )
      scfControls.doExtrap = false;


    return scfControls;

  }; // CQSCFOptions

  std::set<std::string> CQORBPROP_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "RDF",
      "INITIALRAD",
      "FINALRAD",
      "RADPTS",
      "ANGPTS",
      "ORBEN",
      "NUMMOS"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);

  }

  void ParseOrbitalPropSubsection(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> ss) {

    // check if [SCF/ORBPROP] section
    if( not input.containsSection("SCF/ORBPROP") ) return;

    std::cout << " Found [SCF/ORBPROP] section" << std::endl;

    ss->orbProp = true; //Do orbital properties if subsection present

    std::set<std::string> invalidKeywords = CQORBPROP_VALID(input.getSection("SCF/ORBPROP"));
    printInvalidKeys(invalidKeywords, "SCF/ORBPROP");

    // Turn on RDF analysis
    OPTOPT(
      ss->doRDFs = input.getData<bool>("SCF/ORBPROP/RDF");
    );

    // Turn on orbital energy analysis
    OPTOPT(
      ss->doOrbEne = input.getData<bool>("SCF/ORBPROP/ORBEN");
    );

    // Check variables related to RDF analysis
    if( ss->doRDFs ){

      OPTOPT(
        ss->initialRad =
          input.getData<double>("SCF/ORBPROP/INITIALRAD");
      );

      OPTOPT(
        ss->finalRad =
          input.getData<double>("SCF/ORBPROP/FINALRAD");
      );

      OPTOPT(
        ss->numRadPts =
          input.getData<size_t>("SCF/ORBPROP/RADPTS");
      );

      OPTOPT(
        ss->numAngPts =
          input.getData<size_t>("SCF/ORBPROP/ANGPTS");
      );

      OPTOPT(
        ss->numMOs =
          input.getData<size_t>("SCF/ORBPROP/NUMMOS");
      );

      size_t N = ss->numAngPts;
      if( N!=6 && N!=14 && N!=26 && N!=38 && N!=50 && N!=74 && N!=86 && N!=146 &&
        N!=170 && N!=194 && N!=230 && N!=266 && N!=302 && N!=590 && N!=974 )
          CErr("Number of Angular Points NYI. See src/grid/lebedev.cxx");

    }

  }

}; // namespace ChronusQ
