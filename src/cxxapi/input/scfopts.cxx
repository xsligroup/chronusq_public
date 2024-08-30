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


  void CQSCF_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
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
      "CUBE"
    };

    // Specified keywords
    std::vector<std::string> scfKeywords = input.getDataInSection("SCF");

    // Make sure all of scfKeywords in allowedKeywords

    for( auto &keyword : scfKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword SCF." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }


  void HandlePostSCFRestarts(std::ostream &out, CQInputFile &input,
                             SCFControls &scfControls) {

    bool restartRT = false;
    bool restartMD = false;
    
    if ( input.containsSection("RT") )
      OPTOPT( restartRT |= input.getData<bool>("RT.RESTART"); )
    
    if ( input.containsSection("DYNAMICS") ){
      std::string restart = "FALSE";
      OPTOPT( restart = input.getData<std::string>("DYNAMICS.RESTART");)
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
      if (restartMD and parseJob(input.getData<std::string>("QM.JOB"))==JobType::BOMD )
        return;
      out << "  *** RESTART requested -- SCF.GUESS set to READMO and SCF.ALG set to SKIP ***";
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
    SingleSlaterBase &ss) {

    // MO swapping
    std::string swapMOStrings;
    OPTOPT( swapMOStrings = input.getData<std::string>("SCF.SWAPMO"));
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
      input.getData<double>("SCF.ACCURACY"); )

    scfControls.maxdPConvTol = scfControls.rmsdPConvTol*100;
    scfControls.eneConvTol = scfControls.rmsdPConvTol*100;

    // Energy convergence tolerance
    //OPTOPT( scfControls.eneConvTol =
    //          input.getData<double>("SCF.ENETOL"); )

    // Energy convergence tolerance
    //OPTOPT( scfControls.denConvTol =
    //          input.getData<double>("SCF.DENTOL"); )

    // Energy Gradient convergence tolerance
    //OPTOPT( scfControls.FDCConvTol =
    //          input.getData<double>("SCF.FDCTOL"); )

    // Maximum SCF iterations
    OPTOPT( scfControls.maxSCFIter =
              input.getData<size_t>("SCF.MAXITER"); )


    // Incremental Fock Options
    OPTOPT(
      scfControls.doIncFock = input.getData<bool>("SCF.INCFOCK");
    )
    OPTOPT(
      scfControls.nIncFock = input.getData<size_t>("SCF.NINCFOCK");
    )


    // Guess
    std::string guessString = "SAD";
    OPTOPT( guessString = input.getData<std::string>("SCF.GUESS"); )
    trim(guessString);
    if( not guessString.compare("CORE") )
      scfControls.guess = CORE;
    else if( not guessString.compare("SAD") )
      scfControls.guess = SAD;
    else if( not guessString.compare("TIGHT") )
      scfControls.guess = TIGHT;
    else if( not guessString.compare("RANDOM") )
      scfControls.guess = RANDOM;
    else if( not guessString.compare("READMO") )
      scfControls.guess = READMO;
    else if( not guessString.compare("READDEN") )
      scfControls.guess = READDEN;
    else if( not guessString.compare("FCHKMO") )
      scfControls.guess = FCHKMO;
    else if( not guessString.compare("CLASSICAL") )
      scfControls.guess = NEOConvergeClassical;
    else
      CErr("Unrecognized entry for SCF.GUESS");
    

    // Proton Guess For NEO Calculations
    std::string prot_guessString = "TIGHT";
    OPTOPT( prot_guessString = input.getData<std::string>("SCF.PROT_GUESS"); )
    trim(prot_guessString);
    if( not prot_guessString.compare("CORE") )
      scfControls.prot_guess = CORE;
    else if( not prot_guessString.compare("RANDOM") )
      scfControls.prot_guess = RANDOM;
    else if( not prot_guessString.compare("READMO") )
      scfControls.prot_guess = READMO;
    else if( not prot_guessString.compare("READDEN") )
      scfControls.prot_guess = READDEN;
    else if( not prot_guessString.compare("TIGHT") )
      scfControls.prot_guess = NEOTightProton;
    else
      CErr("Unrecognized entry for SCF.PROT_GUESS");
    

    // ALGORITHM
    std::string algString = "CONVENTIONAL";
    OPTOPT( algString = input.getData<std::string>("SCF.ALG"); )
    if( not algString.compare("CONVENTIONAL") )
      scfControls.scfAlg = _CONVENTIONAL_SCF;
    else if( not algString.compare("NR") )
      scfControls.scfAlg = _NEWTON_RAPHSON_SCF;
    else if( not algString.compare("SKIP") ) {
      scfControls.scfAlg = _CONVENTIONAL_SCF;
      scfControls.energyOnly = true;
    } else 
      CErr("Unrecognized entry for SCF.ALG!");


    // Newton-Raphson SCF Approximation
    std::string nrAlgString = "BFGS";
    OPTOPT( nrAlgString = input.getData<std::string>("SCF.NRAPPROX"); )
    if( not nrAlgString.compare("FULL") )
      scfControls.nrAlg = FULL_NR;
    else if( not nrAlgString.compare("BFGS") )
      scfControls.nrAlg = QUASI_BFGS;
    else if( not nrAlgString.compare("SR1") )
      scfControls.nrAlg = QUASI_SR1;
    else if( not nrAlgString.compare("GRADDESCENT") )
      scfControls.nrAlg = GRAD_DESCENT;
    else 
      CErr("Unrecognized entry for SCF.NRAPPROX");


    // Newton-Raphson SCF Initial trust region and level-shift
    OPTOPT(
      scfControls.nrTrust = input.getData<double>("SCF.NRTRUST");
    )
    OPTOPT(
      scfControls.nrLevelShift = input.getData<double>("SCF.NRLEVELSHIFT");
    )


    // Restart jobs
    HandlePostSCFRestarts(out, input, scfControls);


    // Toggle extrapolation in its entireity
    OPTOPT(
      scfControls.doExtrap =
        input.getData<bool>("SCF.EXTRAP");
    )

    // Handle DIIS options
    std::string diisAlgString = "CDIIS"; 
    OPTOPT( diisAlgString = input.getData<std::string>("SCF.DIISALG"); )
    if( not diisAlgString.compare("CEDIIS"))
      scfControls.diisAlg = CEDIIS;
    else if( not diisAlgString.compare("CDIIS"))
      scfControls.diisAlg = CDIIS;
    else if( not diisAlgString.compare("EDIIS"))
      scfControls.diisAlg = EDIIS;
    else
        CErr("Unrecognized entry for SCF.DIISALG!");


    // Check if it specifically says DIIS=FALSE
    OPTOPT(
      bool doDIIS = input.getData<bool>("SCF.DIIS");
      if( not doDIIS ) 
        scfControls.diisAlg = NONE;
    );

    // Number of terms for keep for DIIS
    OPTOPT( scfControls.nKeep = input.getData<size_t>("SCF.NKEEP"); )
    if( scfControls.nKeep < 1 )
      scfControls.nKeep = 1;

    // Point at which to switch from EDIIS to CDIIS
    OPTOPT( scfControls.cediisSwitch = input.getData<double>("SCF.SWITCH"); )
    if( scfControls.cediisSwitch <= 0. )
      CErr("CEDIIS Switch is less than or equal to zero");

    // Parse Damping options
    OPTOPT(
      scfControls.doDamp = input.getData<bool>("SCF.DAMP");
    );

    OPTOPT(
      scfControls.dampStartParam =
        input.getData<double>("SCF.DAMPPARAM");
    );

    OPTOPT(
      scfControls.dampError =
        input.getData<double>("SCF.DAMPERROR");
    );


    // SCF Field
    std::string fieldStr;
    OPTOPT(
      fieldStr = input.getData<std::string>("SCF.FIELD");
    )
    EMPerturbation parsedField;
    if (!fieldStr.empty())
        handleField(fieldStr, parsedField);
    pert.addField(parsedField);

    // Printing Options
    if ( input.containsData("SCF.PRINTMOS") ) {
      try { scfControls.printMOCoeffs = input.getData<size_t>("SCF.PRINTMOS"); }
      catch(...) {
        CErr("Invalid PRINTMOS input. Please use number 0 ~ 9.");
      }
    }
    if (scfControls.printMOCoeffs >= 10 ) CErr("SCF print level is not valid!");


    // Parse whether to print contraction timing during SCF
    OPTOPT( scfControls.printContractionTiming = input.getData<bool>("SCF.PRINTCONTRACTIONTIMING"); )





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

}; // namespace ChronusQ
