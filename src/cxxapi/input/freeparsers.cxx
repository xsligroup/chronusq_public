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
#include <regex>
#include <orbitalmodifieroptions.hpp>

namespace ChronusQ {

  std::string doubleToString(double value) {
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return out.str();
  }


  void CQInputFile::parseFreeCQInput (std::string &line){

    /********************************************************************************/
    /* CQ Free Format Input                                                         */
    /* example, CQ= HF/STO-3G NEO(EPC17/PROT-BP4-D) SCF(accuracy=1.e-6)             */
    /* example, ChronuQ: X2C-HF/CD-6-31G RT(time=10fs, stepsize = 1as)              */
    /* example, ChronuQ= 4C-HF/ano-rcc hamiltonian(DCB, scalar, atomic)             */
    /********************************************************************************/
//    auto const freeCQInput = std::regex("CQ[[:blank:]]*=|CQ[[:blank:]]*:|CHRONUSQ[[:blank:]]*=|CHRONUSQ[[:blank:]]*:",std::regex_constants::icase);
//    if(!std::regex_search(line, freeCQInput)) return;
//    std::cout<<"xsli test CQ Input"<<std::endl;
//    line = std::regex_replace(line, freeCQInput, "");

    // Parse NEO Section
    parseFreeCQInputNEO(line);
    parseFreeCQInputElectron(line);
    parseFreeCQInputSCF(line);
    parseFreeCQInputSSGuess(line);
    parseFreeCQInputRT(line);
    parseFreeCQInputField(line);
    parseFreeCQInputCI(line);


    auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
    line = std::regex_replace(line, freeDividers, " ");
    // the second call remove multiple space left after replace ','
    line = std::regex_replace(line, freeDividers, " ");
    if(line.size()>0) std::cout<<"CQ input ignored: "<< line <<std::endl;

  }; // Free Format Input Parser



  void CQInputFile::parseFreeCQInputNEO (std::string &line){

    auto const freeCQInputHF    = std::regex("((2C)|(X2C)|(4C)|(G)-?)?HF/?",std::regex_constants::icase);
    auto const freeCQInputCCSD  = std::regex("((2C)|(X2C)|(4C)|(G)-?)?CCSD((T)|(\\(T\\)))?/?",std::regex_constants::icase);

    /*************************************/
    /* NEO Input                         */
    /* example, NEO(EPC19/prot-pb6-g)    */
    /* example, NEO(EPC19/CD-prot-pb6-g) */
    /* example, NEO(EPC19/ri-prot-pb6-g) */
    /*************************************/
    auto const freeCQInputEPC17 = std::regex("EPC17",std::regex_constants::icase);
    auto const freeCQInputEPC19 = std::regex("EPC19",std::regex_constants::icase);

    auto const freeCQInputPROTSP    = std::regex("((CD)|(RI)-?)?PROT-SP",std::regex_constants::icase);
    auto const freeCQInputPROTPB4D  = std::regex("((CD)|(RI)-?)?PROT-PB4-D",std::regex_constants::icase);
    auto const freeCQInputPROTPB4F1 = std::regex("((CD)|(RI)-?)?PROT-PB4-F1",std::regex_constants::icase);
    auto const freeCQInputPROTPB4F2 = std::regex("((CD)|(RI)-?)?PROT-PB4-F2",std::regex_constants::icase);
    auto const freeCQInputPROTPB5G  = std::regex("((CD)|(RI)-?)?PROT-PB5-G",std::regex_constants::icase);
    auto const freeCQInputPROTPB6G  = std::regex("((CD)|(RI)-?)?PROT-PB6-G",std::regex_constants::icase);

    auto const freeCQInputNEO = std::regex("NEO(\\((.*?)\\))?",std::regex_constants::icase);
    std::smatch NEOmatch;

    if( std::regex_search(line, NEOmatch, freeCQInputNEO) ){
      // str(2) captures what is inside NEO()
      if(NEOmatch.str(2).size()>0) {
        // Parse user-defined input
        std::cout<<"xsli test NEO Section "<<std::endl;
        std::string NEOInputOptions = NEOmatch.str(2);

        // Methods
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputHF) ) {
          std::cout<<"xsli test NEO HF"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputEPC17) ) {
          std::cout<<"xsli test NEO EPC17"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputEPC19) ) {
          std::cout<<"xsli test NEO EPC19"<<std::endl;
        }

        // Basis Sets
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTSP) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PORT-SP"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PORT-SP"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PORT-SP"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4D) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-D"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4F1) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-F1"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-F1"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-F1"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4F2) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-F2"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-F2"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-F2"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB5G) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-D"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB6G) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB5-G"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB5-G"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB5-G"<<std::endl;
        }

      } else {
        // Choose default parameters
      }
      // We need to delete the NEO section so that we can parse the electronic section properly
      line = std::regex_replace(line, freeCQInputNEO, "");
    } // NEO Input

  };



  void CQInputFile::parseFreeCQInputElectron (std::string &line){

    auto const freeCQInputHF    = std::regex("((2C)|(X2C)|(4C)|(G)-?)?HF/?",std::regex_constants::icase);
    // Match strings                             1  2      3       4     5         6 7      8
    auto const freeCQInputCCSD  = std::regex("((2C)|(X2C)|(4C)|(G)-?)?CCSD((T)|(\\(T\\)))?/?",std::regex_constants::icase);

    /****************************/
    /* Electronic Input         */
    /* example, B3LYP/6-31G     */
    /* example, 4C-B3LYP/6-31G  */
    /* example, X2C-B3LYP/6-31G */
    /* example, B3LYP/RI-6-31G  */
    /* example, B3LYP/CD-6-31G  */
    /****************************/

    auto const freeCQInputB3LYP = std::regex("((2C)|(X2C)|(4C)|(G)-?)?B3LYP/?",std::regex_constants::icase);
    auto const freeCQInputPBE   = std::regex("((2C)|(X2C)|(4C)|(G)-?)?PBE/?",std::regex_constants::icase);

    std::smatch methodMatch;
    if ( std::regex_search(line, methodMatch, freeCQInputHF) ) {
      if( methodMatch.str(1).size()==0 ) std::cout<< "xsli test HF " <<std::endl;
      if( methodMatch.str(2).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(2)<<std::endl;
      if( methodMatch.str(3).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(3)<<std::endl;
      if( methodMatch.str(4).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(4)<<std::endl;
      if( methodMatch.str(5).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(5)<<std::endl;
      line = std::regex_replace(line, freeCQInputHF, "");
    }
    else if ( std::regex_search(line, methodMatch, freeCQInputCCSD) ) {
      if( methodMatch.str(1).size()==0 ) std::cout<< "xsli test CCSD " <<std::endl;
      if( methodMatch.str(2).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(2)<<std::endl;
      if( methodMatch.str(3).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(3)<<std::endl;
      if( methodMatch.str(4).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(4)<<std::endl;
      if( methodMatch.str(5).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(5)<<std::endl;
      if( methodMatch.str(7).size()>0 ) std::cout<< "xsli test CCSDT "<<std::endl;
      if( methodMatch.str(8).size()>0 ) std::cout<< "xsli test CCSD(T) "<<std::endl;
      line = std::regex_replace(line, freeCQInputCCSD, "");
    }


    auto const freeCQInputSTO3G   = std::regex("((CD)|(RI)-?)?STO-3G",std::regex_constants::icase);
    auto const freeCQInput321G    = std::regex("((CD)|(RI)-?)?3-21G",std::regex_constants::icase);
    auto const freeCQInput631G    = std::regex("((CD)|(RI)-?)?6-31G",std::regex_constants::icase);
    auto const freeCQInput6311G   = std::regex("((CD)|(RI)-?)?6-311G",std::regex_constants::icase);

    std::smatch basisMatch;
    if ( std::regex_search(line, basisMatch, freeCQInputSTO3G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test STO-3G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test CD-STO-3G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test RI-STO-3G"<<std::endl;
      line = std::regex_replace(line, freeCQInputSTO3G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput321G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 3-21G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test CD-3-21G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test RI-3-21G"<<std::endl;
      line = std::regex_replace(line, freeCQInput321G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput631G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      line = std::regex_replace(line, freeCQInput631G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput6311G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 6-311G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test 6-311GG"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test 6-311G"<<std::endl;
      line = std::regex_replace(line, freeCQInput6311G, "");
    }

  };

  void CQInputFile::parseFreeCQInputSCF (std::string &line) {

    /**************************************************/
    /* SCF Input                                      */
    /* example, SCF(accuracy=1.e-8)                   */
    /* example, SCF(cdiis, maxiteration=100)          */
    /* example, SCF(energyonly)                       */
    /**************************************************/

    auto const freeCQInputSCF = std::regex("(SCF)(\\((.*?)\\))", std::regex_constants::icase);
    std::smatch scfMatch;

    if (std::regex_search(line, scfMatch, freeCQInputSCF)) {
      std::string scfInputOptions = scfMatch.str(3);

      // read in SCF accuracy
      auto const freeCQInputSCFAccuracy = std::regex("accuracy\\s*=\\s*((\\d+\\.?\\d*|\\.\\d+)(e[-+]?\\d+)?)\\s*([,;:]|$)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFAccuracy) ) {
        addData("SCF.ACCURACY", scfMatch.str(1));
        std::cout<<"xsli test read in accuracy = "<<std::stod(scfMatch.str(1))<<std::endl;
      }

      auto const freeCQInputEnergyOnly = std::regex("energyonly|skip", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputEnergyOnly) ) {
        addData("SCF.ENERGYONLY", "SKIP");
        std::cout<<"xsli test read in energyonly"<<std::endl;
      }

      auto const freeCQInputSCFMethod = std::regex("(diis)|(nodiis)|(cdiis)|(ediis)|(qc)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFMethod) ) {
        if(!scfMatch.str(1).empty()) addData("SCF.DIISALG","CDIIS");// std::cout<<"xsli test read in SCF method = DIIS "<<std::stod(scfMatch.str(1))<<std::endl;
        if(!scfMatch.str(2).empty()) addData("SCF.DIISALG","NONE"); // std::cout<<"xsli test read in SCF method = NoDIIS "<<std::stod(scfMatch.str(2))<<std::endl;
        if(!scfMatch.str(3).empty()) addData("SCF.DIISALG","CDIIS"); // std::cout<<"xsli test read in SCF method = CDIIS "<<std::stod(scfMatch.str(3))<<std::endl;
        if(!scfMatch.str(4).empty()) addData("SCF.DIISALG","EDIIS"); // std::cout<<"xsli test read in SCF method = EnDIIS "<<std::stod(scfMatch.str(4))<<std::endl;
        if(!scfMatch.str(4).empty()) addData("SCF.ALG","NR"); // std::cout<<"xsli test read in SCF method = QC "<<std::stod(scfMatch.str(5))<<std::endl;
      }

      // Check for SCF maxSteps
      auto const freeCQInputSCFSteps = std::regex("(MAXSTEP|MAXSTEPS|MAXITERATION|MAXITERATIONS|MAXCYCLE|MAXCYCLES)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
      if ( std::regex_search(scfInputOptions, scfMatch, freeCQInputSCFSteps) ) {
        addData("SCF.MAXITER", scfMatch.str(2));
        std::cout<<"xsli test read in MAXITERATIONS = "<<std::stod(scfMatch.str(2))<<std::endl;
      }
    }

  };

  void SCFControls::parseSection(const InputMap &dict) {
    if (dict.count("ENERGYONLY")) {
      scfAlg = _CONVENTIONAL_SCF;
      energyOnly = true;
    }
    if (dict.count("MAXITER")) maxSCFIter = std::stoi(dict.at("MAXITER"));
    if (dict.count("ACCURACY")) {
      rmsdPConvTol = std::stod(dict.at("ACCURACY"));
      maxdPConvTol = rmsdPConvTol*100;
      eneConvTol   = rmsdPConvTol*100;
    }

      if (dict.at("DIISALG") == "CDIIS") diisAlg = CDIIS;
      else if (dict.at("DIISALG") == "EDIIS") diisAlg = EDIIS;
      else if (dict.at("DIISALG") == "DIIS") diisAlg = CDIIS;
      else if (dict.at("DIISALG") == "NODIIS") diisAlg = NONE;

      if (dict.at("ALG") == "QC") scfAlg = _NEWTON_RAPHSON_SCF;

      if (dict.count("MAXITER")) maxSCFIter = std::stoi(dict.at("MAXITER"));
  }

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

        auto const freeCQInputGuessType = std::regex("(SAD)|(CORE)|(READ)|(FCHK)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputGuessType) ) {
//          if(!Guessmatch.str(1).empty()) ssGuessOptions.electronicGuess = SADGuess;
//          if(!Guessmatch.str(2).empty()) ssGuessOptions.electronicGuess = CoreGuess;
//          if(!Guessmatch.str(3).empty()) ssGuessOptions.electronicGuess = ReadBin;
//          if(!Guessmatch.str(4).empty()) ssGuessOptions.electronicGuess = ReadGaussFCHK;
          if(!Guessmatch.str(1).empty()) addData("SCF.GUESS", "SAD");
          if(!Guessmatch.str(2).empty()) addData("SCF.GUESS", "CORE");
          if(!Guessmatch.str(3).empty()) addData("SCF.GUESS", "READMO");
          if(!Guessmatch.str(4).empty()) addData("SCF.GUESS", "FCHKMO");
        }
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessType, "");

        auto const freeCQInputGuessSwap = std::regex("(SWAP)\\s*=((\\s*([ab]?)(\\d+)-\\4(\\d+))|(homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
        while( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputGuessSwap) ) {
          if(!Guessmatch.str(3).empty()) {
            std::cout<<"xsli test guess swap = "<<Guessmatch.str(4)<<Guessmatch.str(5)<<Guessmatch.str(6)<<std::endl;
            if(Guessmatch.str(4)=="a" or Guessmatch.str(4)=="A" or Guessmatch.str(4).size()==0)
              ssGuessOptions.alphaElectronicMOSwap.push_back({(size_t)std::stoi(Guessmatch.str(5)), (size_t)std::stoi(Guessmatch.str(5))});
            else if(Guessmatch.str(4)=="b" or Guessmatch.str(4)=="B")
              ssGuessOptions.betaElectronicMOSwap.push_back({(size_t)std::stoi(Guessmatch.str(5)), (size_t)std::stoi(Guessmatch.str(5))});

            auto const freeCQInputGuessSwap1 = std::regex("(SWAP)\\s*=\\s*([ab]?)(\\d+)-\\2(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessSwap1, "");
          }
          if(!Guessmatch.str(7).empty() or !Guessmatch.str(8).empty()) {
            std::cout<<"xsli test guess swap = HOMO-LUMO"<<std::endl;
            ssGuessOptions.alphaElectronicMOSwap.push_back({0,0});
            auto const freeCQInputGuessSwap2 = std::regex("(SWAP)\\s*=((homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputGuessSwap2, "");
          }
        }

        for (size_t i = 0; i < ssGuessOptions.alphaElectronicMOSwap.size(); i++) {
          for (size_t j = 0; j < ssGuessOptions.alphaElectronicMOSwap[i].size(); j++) {
            addData("SCF.GUESS.ALPHASWAP[" + std::to_string(i) + "][" + std::to_string(j) + "]",
                    std::to_string(ssGuessOptions.alphaElectronicMOSwap[i][j]));
          }
        }

        for (size_t i = 0; i < ssGuessOptions.betaElectronicMOSwap.size(); i++) {
          for (size_t j = 0; j < ssGuessOptions.betaElectronicMOSwap[i].size(); j++) {
            addData("SCF.GUESS.BETASWAP[" + std::to_string(i) + "][" + std::to_string(j) + "]",
                    std::to_string(ssGuessOptions.betaElectronicMOSwap[i][j]));
          }
        }

        auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeDividers, "");
        if(!GuessInputOptions.empty()) CErr("Unrecognized Guess Input Options: "+GuessInputOptions);
        else line = std::regex_replace(line, freeCQInputGuess, "");
      }

    }

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
          if(!Guessmatch.str(1).empty()) addData("SCF.PROT_GUESS", "SAD");
          if(!Guessmatch.str(2).empty()) addData("SCF.PROT_GUESS", "CORE");
          if(!Guessmatch.str(3).empty()) addData("SCF.PROT_GUESS", "READMO");
          if(!Guessmatch.str(4).empty()) addData("SCF.PROT_GUESS", "FCHKMO");
          if(!Guessmatch.str(5).empty()) addData("SCF.PROT_GUESS", "CLASSICAL");
        }
        GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputNEOGuessType, "");

        auto const freeCQInputNEOGuessSwap = std::regex("(SWAP)\\s*=((\\s*([ab]?)(\\d+)-\\4(\\d+))|(homo-lumo)|(lumo-homo))\\s*([,;:]|$)", std::regex_constants::icase);
        while( std::regex_search(GuessInputOptions, Guessmatch, freeCQInputNEOGuessSwap) ) {
          if(!Guessmatch.str(3).empty()) {
            std::cout<<"xsli test guess swap = "<<Guessmatch.str(4)<<Guessmatch.str(5)<<Guessmatch.str(6)<<std::endl;
            auto const freeCQInputNEOGuessSwap1 = std::regex("(SWAP)\\s*=\\s*([ab]?)(\\d+)-\\2(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
            GuessInputOptions = std::regex_replace(GuessInputOptions, freeCQInputNEOGuessSwap1, "");
          }
          if(!Guessmatch.str(7).empty() or !Guessmatch.str(8).empty()) {
            std::cout<<"xsli test guess swap = HOMO-LUMO"<<std::endl;
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
  };

  void SingleSlaterGuessOptions::parseSection(const InputMap &dict) {
    if (dict.count("GUESS")) {
      if (dict.at("GUESS") == "CORE") electronicGuess = CoreGuess;
      else if (dict.at("GUESS") == "SAD") electronicGuess = SADGuess;
      else if (dict.at("GUESS") == "TIGHT") electronicGuess = TightGuess;
      else if (dict.at("GUESS") == "RANDOM") electronicGuess = RandomGuess;
      else if (dict.at("GUESS") == "READMO") electronicGuess = ReadBin;
      else if (dict.at("GUESS") == "FCHKMO") electronicGuess = ReadGaussFCHK;
      else if (dict.at("GUESS") == "CLASSICAL") electronicGuess = ClassicalGuess;
    }
    if (dict.count("PROT_GUESS")) {
      if (dict.at("PROT_GUESS") == "CORE") nuclearGuess = CoreGuess;
      else if (dict.at("PROT_GUESS") == "SAD") nuclearGuess = SADGuess;
      else if (dict.at("PROT_GUESS") == "TIGHT") nuclearGuess = TightGuess;
      else if (dict.at("PROT_GUESS") == "RANDOM") nuclearGuess = RandomGuess;
      else if (dict.at("PROT_GUESS") == "READMO") nuclearGuess = ReadBin;
      else if (dict.at("PROT_GUESS") == "FCHKMO") nuclearGuess = ReadGaussFCHK;
      else if (dict.at("PROT_GUESS") == "CLASSICAL") nuclearGuess = ClassicalGuess;
    }

  }

  void CQInputFile::parseFreeCQInputRT (std::string &line) {

    /***********************************************************************/
    /* RT Input                                                            */
    /* example, RT(MMUT, stepsize= 1.0 as, nsteps = 1000, maxtime = 1.0fs) */
    /***********************************************************************/

    TDSCFOptions tdSCFControls;
    auto const freeCQInputRT = std::regex("(RT|REALTIME)(\\((.*?)\\))?", std::regex_constants::icase);
    std::smatch RTmatch;

    if (std::regex_search(line, RTmatch, freeCQInputRT)) {
      // str(3) captures the input inside the parentheses
      if (RTmatch.str(3).size() > 0) {

        std::string RTInputOptions = RTmatch.str(3);

        // Check for timestep
        auto const freeCQInputStepsize = std::regex("(STEPSIZE|TIMESTEP)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|(au))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputStepsize) ) {
          tdSCFControls.deltaT = std::stod(RTmatch.str(2));
          if(!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) tdSCFControls.deltaT/=(FSPerAUTime*1.e3);
          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty()) tdSCFControls.deltaT/=FSPerAUTime;
          std::cout<<"xsli test read in timestep = "<<tdSCFControls.deltaT<<" au"<<std::endl;
          addData("RT.DELTAT", doubleToString(tdSCFControls.deltaT));
          addData("RT.UNITS", "AU");
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputStepsize, "");
        }

        // Check for maxSteps
        auto const freeCQInputNSteps = std::regex("(NSTEPS|MAXSTEPS|MAXITERATIONS|MAXITER)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputNSteps) ) {
          tdSCFControls.maxSteps = std::stoi(RTmatch.str(2));
          tdSCFControls.tMax = tdSCFControls.maxSteps*tdSCFControls.deltaT;
          std::cout<<"xsli test read in maxSteps = "<<tdSCFControls.maxSteps<<std::endl;
          addData("RT.TMAX", doubleToString(tdSCFControls.tMax));
          addData("RT.MAXSTEPS", std::to_string(tdSCFControls.maxSteps));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputNSteps, "");
        }

        // Check for maxTime - this will override maxSteps
        auto const freeCQInputMaxTime = std::regex("(MAXTIME)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((as)|(attosecond)|(fs)|(femtosecond)|(au))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputMaxTime) ) {
          tdSCFControls.tMax = std::stod(RTmatch.str(2));
          if(!RTmatch.str(5).empty() or !RTmatch.str(6).empty()) tdSCFControls.tMax/=(FSPerAUTime*1.e3);
          else if (!RTmatch.str(7).empty() or !RTmatch.str(8).empty()) tdSCFControls.tMax/=FSPerAUTime;
          tdSCFControls.maxSteps = (tdSCFControls.tMax + tdSCFControls.deltaT / 4) / tdSCFControls.deltaT;
          std::cout<<"xsli test read in maxTime = "<<tdSCFControls.tMax<<std::endl;
          addData("RT.TMAX", doubleToString(tdSCFControls.tMax));
          addData("RT.MAXSTEPS", std::to_string(tdSCFControls.maxSteps));
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
          addData("RT.RESTARTSTEP", str);
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRestartAlgorithm, "");
        }

        // Check for RT algorithm, this has to be checked after restart algorithm to avoid conflict
        auto const freeCQInputRTAlgorithm = std::regex("(MMUT)|(MODIFIEDMIDPOINT)|(FORWARDEULER)|(EULER)|(EXPLICITMAGNUS2)|(EXPLICITMAGNUSTWO)|(MAGNUS2)|(MAGNUSTWO)|(RK4)|(RUNGEKUTTAFOURTHORDER)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRTAlgorithm) ) {
//          if (!RTmatch.str(1).empty() or !RTmatch.str(2).empty()) tdSCFControls.integrationAlgorithm = RTModifiedMidpoint;
//          else if (!RTmatch.str(3).empty() or !RTmatch.str(4).empty()) tdSCFControls.integrationAlgorithm = RTForwardEuler;
//          else if (!RTmatch.str(5).empty() or !RTmatch.str(6).empty() or !RTmatch.str(7).empty() or !RTmatch.str(8).empty()) tdSCFControls.integrationAlgorithm = RTExplicitMagnus2;
          std::string str;
          if (!RTmatch.str(1).empty() or !RTmatch.str(2).empty()) str = "MMUT";
          else if (!RTmatch.str(3).empty() or !RTmatch.str(4).empty()) str = "FORWARDEULER";
          else if (!RTmatch.str(5).empty() or !RTmatch.str(6).empty() or !RTmatch.str(7).empty() or !RTmatch.str(8).empty()) str = "MAGNUS2";
          else if (!RTmatch.str(9).empty() or !RTmatch.str(10).empty()) str = "RK4"; // New check for RK4 or RUNGEKUTTAFOURTHORDER
          //std::cout<<"xsli test read in RT algorithm = "<<tdSCFControls.integrationAlgorithm<<std::endl;
          addData("RT.INTALG", str);
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRTAlgorithm, "");
        }

        // Check for the frequency of save
        auto const freeCQInputISave = std::regex("(ISAVE|AUTOSAVE)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputISave) ) {
          tdSCFControls.iSave = std::stoi(RTmatch.str(2));
          std::cout<<"xsli test read in iSave = "<<tdSCFControls.iSave<<std::endl;
          addData("RT.SAVESTEP", std::to_string(tdSCFControls.iSave));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputISave, "");
        }

        // Check for the frequency of print
        auto const freeCQInputIPrint = std::regex("(IPRINT|AUTOPRINT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputIPrint) ) {
          tdSCFControls.iPrint = std::stoi(RTmatch.str(2));
          std::cout<<"xsli test read in iPrint = "<<tdSCFControls.iPrint<<std::endl;
          addData("RT.PRINTSTEP", std::to_string(tdSCFControls.iPrint));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputIPrint, "");
        }

        // Check for the frequency of automatic restart
        auto const freeCQInputIRestart = std::regex("(IRESTART|AUTORESTART)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputIRestart) ) {
          tdSCFControls.iRestart = std::stoi(RTmatch.str(2));
          std::cout<<"xsli test read in iRestart = "<< std::to_string(tdSCFControls.iRestart) <<std::endl;
          addData("RT.IRSTRT", std::to_string(tdSCFControls.iRestart));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputIRestart, "");
        }

        // Check for restart, this has to be checked after iRESTART to avoid conflict
        auto const freeCQInputRestart = std::regex("RESTART(\\s*=\\s*(\\d+))?\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRestart) ) {
          tdSCFControls.restoreFromStep = -1;
          if(!RTmatch.str(2).empty()) tdSCFControls.restoreFromStep = std::stoi(RTmatch.str(2));
          std::cout<<"xsli test read in do restart = "<<tdSCFControls.restoreFromStep<<std::endl;
          addData("RT.RESTARTFROM", std::to_string(tdSCFControls.restoreFromStep));
          addData("RT.RESTART", "TRUE");
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRestart, "");
        }

        auto const freeCQInputrtGaunt = std::regex("(RTGAUNT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtGaunt) ) {
          tdSCFControls.rtGaunt = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtGaunt = "<<tdSCFControls.rtGaunt<<std::endl;
          addData("RT.RTGAUNT", std::to_string(tdSCFControls.rtGaunt));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtGaunt, "");
        }

        auto const freeCQInputrtGauge = std::regex("(RTGAUGE)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtGauge) ) {
          tdSCFControls.rtGauge = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtGauge = "<<tdSCFControls.rtGauge<<std::endl;
          addData("RT.RTGAUGE", std::to_string(tdSCFControls.rtGauge));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtGauge, "");
        }

        auto const freeCQInputrtBreit = std::regex("(RTBREIT)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputrtBreit) ) {
          tdSCFControls.rtBreit = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in rtBreit = "<<tdSCFControls.rtBreit<<std::endl;
          addData("RT.RTBREIT", std::to_string(tdSCFControls.rtBreit));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputrtBreit, "");
        }

        auto const freeCQInputRtprintden = std::regex("(RTPRINTDEN)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputRtprintden) ) {
          tdSCFControls.Rtprintden = std::stoi(RTmatch.str(2));
          std::cout<<"zxc test read in Rtprintden = "<<tdSCFControls.Rtprintden<<std::endl;
          addData("RT.RTPRINTDEN", std::to_string(tdSCFControls.Rtprintden));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputRtprintden, "");
        }

        auto const freeCQInputOrbitalPopFreq = std::regex("(ORBITALPOPFREQ)\\s*=\\s*(\\d+)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(RTInputOptions, RTmatch, freeCQInputOrbitalPopFreq) ) {
          tdSCFControls.orbitalPopFreq = std::stoi(RTmatch.str(2));
          std::cout<<"test read in OrbitalPopFreq = "<<tdSCFControls.orbitalPopFreq<<std::endl;
          addData("RT.ORBITALPOPFREQ", std::to_string(tdSCFControls.orbitalPopFreq));
          RTInputOptions = std::regex_replace(RTInputOptions, freeCQInputOrbitalPopFreq, "");
        }

        auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
        RTInputOptions = std::regex_replace(RTInputOptions, freeDividers, "");
        if(!RTInputOptions.empty()) CErr("Unrecognized RT Input Options: "+RTInputOptions);
        else line = std::regex_replace(line, freeCQInputRT, "");

        

      }
    };
  }

  void TDSCFOptions::parseSection(const InputMap &dict) {
    if (dict.count("DELTAT")) deltaT = std::stod(dict.at("DELTAT"));
    if (dict.count("TMAX")) tMax = std::stod(dict.at("TMAX"));
    if (dict.count("MAXSTEPS")) maxSteps = std::stoi(dict.at("MAXSTEPS"));
    // Initialize maxSteps if not set
    if (maxSteps==0 and (tMax!=0. and deltaT != 0.)) maxSteps = (tMax + deltaT/4) / deltaT;
    if (dict.count("IRSTRT")) iRestart = std::stoi(dict.at("IRSTRT"));
    if (dict.count("SAVESTEP")) iSave = std::stoi(dict.at("SAVESTEP"));
    if (dict.count("PRINTSTEP")) iPrint = std::stoi(dict.at("PRINTSTEP"));
    if (dict.count("RESTARTFROM")) restoreFromStep = std::stoi(dict.at("RESTARTFROM"));
    if (dict.count("INTALG")) {
      if (dict.at("INTALG") == "MMUT") integrationAlgorithm = RealTimeAlgorithm::RTModifiedMidpoint;
      else if (dict.at("INTALG") == "FORWARDEULER") integrationAlgorithm = RealTimeAlgorithm::RTForwardEuler;
      else if (dict.at("INTALG") == "MAGNUS2") integrationAlgorithm = RealTimeAlgorithm::RTExplicitMagnus2;
      else if (dict.at("INTALG") == "RK4") integrationAlgorithm = RealTimeAlgorithm::RTRungeKuttaOrderFour;
    }
    if (dict.count("PROT_INTALG")) {
      if (dict.at("PROT_INTALG") == "MMUT") protIntegrationAlgorithm = RealTimeAlgorithm::RTModifiedMidpoint;
      else if (dict.at("PROT_INTALG") == "FORWARDEULER") protIntegrationAlgorithm = RealTimeAlgorithm::RTForwardEuler;
      else if (dict.at("PROT_INTALG") == "MAGNUS2") protIntegrationAlgorithm = RealTimeAlgorithm::RTExplicitMagnus2;
      else if (dict.at("PROT_INTALG") == "RK4") protIntegrationAlgorithm = RealTimeAlgorithm::RTRungeKuttaOrderFour;
    }
    // If the user didn't provide a specific protonic integration algorithm, use the same as electronic
    if (protIntegrationAlgorithm == RealTimeAlgorithm::Uninitialized) protIntegrationAlgorithm = integrationAlgorithm;
    // Currently we don't allow mixing of MMUT and RK4
    if( protIntegrationAlgorithm == RealTimeAlgorithm::RTRungeKuttaOrderFour and integrationAlgorithm == RealTimeAlgorithm::RTModifiedMidpoint){
      std::cout << "Warining: INTALG=MMUT and PROT_INTALG=RK4. This creates a mismatch in delta T " << std::endl;
      std::cout << "Forcing electronic subsystem to use magnus2 propagation: INTALG=MAGNUS2 " << std::endl;
      integrationAlgorithm = RealTimeAlgorithm::RTExplicitMagnus2;
    }
    if( integrationAlgorithm == RealTimeAlgorithm::RTRungeKuttaOrderFour and protIntegrationAlgorithm == RealTimeAlgorithm::RTModifiedMidpoint){
      std::cout << "Warining: INTALG=RK4 and PROT_INTALG=MMUT. This creates a mismatch in delta T " << std::endl;
      std::cout << "Forcing protonic subsystem to use magnus2 propagation: PROT_INTALG=MAGNUS2 " << std::endl;
      protIntegrationAlgorithm = RealTimeAlgorithm::RTExplicitMagnus2;
    }
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
    /* example, Filed(delta, start= 1.0 as, end = 10 fs, amplitude = 10 au)   */
    /**************************************************************************/
    auto const freeCQInputField = std::regex("FIELD(\\((.*?)\\))?", std::regex_constants::icase);
    std::smatch Fieldmatch;

    if (std::regex_search(line, Fieldmatch, freeCQInputField)) {
      // str(3) captures the input inside the parentheses
      if (Fieldmatch.str(2).size() > 0) {
        std::string FieldInputOptions = Fieldmatch.str(2);
        std::cout << "xsli test Field: " << FieldInputOptions << std::endl;

        // step function
        auto const freeCQInputFieldStep = std::regex("(DELTA|STEP|CONSTANT)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldStep) ) {
          std::cout<<"xsli test read in field type = "<<Fieldmatch.str(0)<<std::endl;
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
              if (!Fieldmatch.str(7).empty() or !Fieldmatch.str(8).empty()) starttime /= (FSPerAUTime * 1.e3);
              else if (!Fieldmatch.str(9).empty() or !Fieldmatch.str(10).empty()) starttime /= FSPerAUTime;
            }
          }
          std::cout<<"xsli test read in tOn = "<<starttime<<" au"<<std::endl;
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
              if(!Fieldmatch.str(7).empty() or !Fieldmatch.str(8).empty()) endtime/=(FSPerAUTime*1.e3);
              else if (!Fieldmatch.str(9).empty() or !Fieldmatch.str(10).empty()) endtime/=FSPerAUTime;
            }
          }
          std::cout<<"xsli test read in tOff = "<<endtime<<" au"<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldEnd, "");
        }


        // field amplitude
        auto const freeCQInputFieldAmp = std::regex("(AMP|AMPLITUDE)\\s*=\\s*(\\d+(\\.\\d+)?)\\s*((EV)|(MEV)|(AU))\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(FieldInputOptions, Fieldmatch, freeCQInputFieldAmp) ) {
          double amp;
          amp = std::stod(Fieldmatch.str(2));
          if(!Fieldmatch.str(5).empty()) amp/=EVPerHartree;
          else if (!Fieldmatch.str(6).empty()) amp/=(EVPerHartree*1.e3);
          std::cout<<"xsli test read in amp = "<<amp<<" au"<<std::endl;
          FieldInputOptions = std::regex_replace(FieldInputOptions, freeCQInputFieldAmp, "");
        }

        auto const freeDividers = std::regex("\\s+|,+", std::regex_constants::icase);
        FieldInputOptions = std::regex_replace(FieldInputOptions, freeDividers, "");
        if (!FieldInputOptions.empty()) CErr("Unrecognized Field Input Options: " + FieldInputOptions);
        else line = std::regex_replace(line, freeCQInputField, "");
      }
    }
  }

  void CQInputFile::parseFreeCQInputCI (std::string &line) {

    /****************************************************************************************/
    /* CI Input                                                                             */
    /* example, CASSCF(10o,5e,accuracy=1.e-4,nstates=5)                                     */
    /* example, X2C-DASCI(40o,20e,nDAS=5,maxexcitation=2)                                   */
    /* example, 4C-RASCI(40o,20e,RAS1(10o,10e,1h),RAS2(20o,10e),RAS3(10o,0e,2p),nstates=10) */
    /****************************************************************************************/

    auto const freeCQInputCI = std::regex("((2C)|(X2C)|(4C)|(G))?(-)?((CASSCF)|(RASSCF)|(DASSCF)|(CASCI)|(RASCI)|(DASCI))\\((([^()]*(\\([^()]*\\))?[^()]*)*)\\)", std::regex_constants::icase);
    std::smatch CIMatch;

    if (std::regex_search(line, CIMatch, freeCQInputCI)) {

      std::cout<< "xsli test CI-0 type: " <<CIMatch.str(0)<<std::endl;
      // str(1) captures the CI type (e.g., X2C, 4C)
      if (CIMatch.str(1).size() > 0) {
        if( CIMatch.str(2).size()>0 ) std::cout<< "xsli test CI type: " <<CIMatch.str(2)<<std::endl;
        if( CIMatch.str(3).size()>0 ) std::cout<< "xsli test CI type: " <<CIMatch.str(3)<<std::endl;
        if( CIMatch.str(4).size()>0 ) std::cout<< "xsli test CI type: " <<CIMatch.str(4)<<std::endl;
        if( CIMatch.str(5).size()>0 ) std::cout<< "xsli test CI type: " <<CIMatch.str(5)<<std::endl;
      }

      bool CASSCF(false),RASSCF(false),DASSCF(false),CASCI(false),RASCI(false),DASCI(false);
      // str(7) captures the CI methods
      if (CIMatch.str(7).size() > 0) {
        if( CIMatch.str(8).size()>0 ) CASSCF = true; //std::cout<< "xsli test CI method: " <<CIMatch.str(8)<<std::endl;
        if( CIMatch.str(9).size()>0 ) RASSCF = true; //std::cout<< "xsli test CI method: " <<CIMatch.str(9)<<std::endl;
        if( CIMatch.str(10).size()>0 ) DASSCF =  true; //std::cout<< "xsli test CI method: " <<CIMatch.str(10)<<std::endl;
        if( CIMatch.str(11).size()>0 ) CASCI =  true; //std::cout<< "xsli test CI method: " <<CIMatch.str(11)<<std::endl;
        if( CIMatch.str(12).size()>0 ) RASCI = true; //std::cout<< "xsli test CI method: " <<CIMatch.str(12)<<std::endl;
        if( CIMatch.str(13).size()>0 ) DASCI =  true; //std::cout<< "xsli test CI method: " <<CIMatch.str(13)<<std::endl;
      }

      // str(14) captures the input inside the parentheses
      if (CIMatch.str(14).size() > 0) {
        std::string CIInputOptions = CIMatch.str(14);
        std::cout << "xsli test CI Options: " << CIInputOptions << std::endl;

        // read in CI eigensolver accuracy
        auto const freeCQInputCIAccuracy = std::regex("accuracy\\s*=\\s*((\\d+\\.?\\d*|\\.\\d+)(e[-+]?\\d+)?)\\s*([,;:]|$)", std::regex_constants::icase);
        if ( std::regex_search(CIInputOptions, CIMatch, freeCQInputCIAccuracy) ) {
          std::cout<<"xsli test read in accuracy = "<<std::stod(CIMatch.str(1))<<std::endl;
        }

        std::regex nEle("\\b([0-9]+)(e|E)\\b"); // an integer number followed by "e" - number of electrons
        std::regex nOrb("\\b([0-9]+)(o|O)\\b"); // an integer number followed by "o" - number of orbitals
        std::regex_search(CIInputOptions, CIMatch, nEle);
        std::cout<<"      nElectrons = "+CIMatch.str(1)<< std::endl;
        std::regex_search(CIInputOptions, CIMatch, nOrb);
        std::cout<<"      nOrbitals = "+CIMatch.str(1)<< std::endl;

        // nstates = an integer number - number of eigenstates to solve for
        std::regex nStates("(nstates)\\s*(=)\\s*([0-9]+)\\s*([,;:]|$)",std::regex_constants::icase);
        std::regex_search(CIInputOptions, CIMatch, nStates);
        std::cout<<"      nStates = "+CIMatch.str(3)<< std::endl;

        if(DASCI or DASSCF) {
          // nDAS = an integer number - number of DAS
          std::regex nDAS("(ndas)\\s*(=)\\s*([0-9]+)\\s*([,;:]|$)", std::regex_constants::icase);
          std::regex_search(CIInputOptions, CIMatch, nDAS);
          std::cout << "      nDAS = " + CIMatch.str(3) << std::endl;

          // maxexcitation = an integer number - max number of inter-space excitations
          std::regex MaxExcitation("(maxexcitation)\\s*(=)\\s*([0-9]+)\\s*([,;:]|$)", std::regex_constants::icase);
          std::regex_search(CIInputOptions, CIMatch, MaxExcitation);
          std::cout << "      MaxExcitation = " + CIMatch.str(3) << std::endl;
        }

        if(RASCI or RASSCF) {
          auto const freeCQInputRAS1 = std::regex("RAS1(\\((.*?)\\))", std::regex_constants::icase);
          auto const freeCQInputRAS2 = std::regex("RAS2(\\((.*?)\\))", std::regex_constants::icase);
          auto const freeCQInputRAS3 = std::regex("RAS3(\\((.*?)\\))", std::regex_constants::icase);
          std::regex nHole("\\b([0-9]+)(h|H)\\b"); // an integer number followed by "e" - number of electrons
          std::regex nParticle("\\b([0-9]+)(p|P)\\b"); // an integer number followed by "o" - number of orbitals

          if(std::regex_search(CIInputOptions, CIMatch, freeCQInputRAS1)) {
            std::string RAS1InputOptions = CIMatch.str(1);
            std::regex_search(RAS1InputOptions, CIMatch, nEle);
            std::cout<<"      RAS1-nElectrons = "+CIMatch.str(1)<< std::endl;
            std::regex_search(RAS1InputOptions, CIMatch, nOrb);
            std::cout<<"      RAS1-nOrbitals = "+CIMatch.str(1)<< std::endl;
            if(std::regex_search(RAS1InputOptions, CIMatch, nHole)) {
              std::cout<<"      RAS1-nHoles = "+CIMatch.str(1)<< std::endl;
            }
            if(std::regex_search(RAS1InputOptions, CIMatch, nParticle)){
              std::cout<<"      RAS1-nParticles = "+CIMatch.str(1)<< std::endl;
            }
          }

          if(std::regex_search(CIInputOptions, CIMatch, freeCQInputRAS2)) {
            std::string RAS2InputOptions = CIMatch.str(1);
            std::regex_search(RAS2InputOptions, CIMatch, nEle);
            std::cout<<"      RAS2-nElectrons = "+CIMatch.str(1)<< std::endl;
            std::regex_search(RAS2InputOptions, CIMatch, nOrb);
            std::cout<<"      RAS2-nOrbitals = "+CIMatch.str(1)<< std::endl;
            if(std::regex_search(RAS2InputOptions, CIMatch, nHole)) {
              std::cout<<"      RAS2-nHoles = "+CIMatch.str(1)<< std::endl;
            }
            if(std::regex_search(RAS2InputOptions, CIMatch, nParticle)){
              std::cout<<"      RAS2-nParticles = "+CIMatch.str(1)<< std::endl;
            }
          }

          if(std::regex_search(CIInputOptions, CIMatch, freeCQInputRAS3)) {
            std::string RAS3InputOptions = CIMatch.str(1);
            std::regex_search(RAS3InputOptions, CIMatch, nEle);
            std::cout<<"      RAS3-nElectrons = "+CIMatch.str(1)<< std::endl;
            std::regex_search(RAS3InputOptions, CIMatch, nOrb);
            std::cout<<"      RAS3-nOrbitals = "+CIMatch.str(1)<< std::endl;
            if(std::regex_search(RAS3InputOptions, CIMatch, nHole)) {
              std::cout<<"      RAS3-nHoles = "+CIMatch.str(1)<< std::endl;
            }
            if(std::regex_search(RAS3InputOptions, CIMatch, nParticle)){
              std::cout<<"      RAS3-nParticles = "+CIMatch.str(1)<< std::endl;
            }
          }

        }

      }

      line = std::regex_replace(line, freeCQInputCI, "");
    }

  };

}; // namespace ChronusQ

