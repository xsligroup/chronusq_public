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
#include <cxxapi/options.hpp>
#include <cxxapi/output.hpp>
#include <cerr.hpp>
#include <regex>
#include <set>
#include <corehbuilder.hpp>
#include <corehbuilder/nonrel.hpp>
#include <corehbuilder/fourcomp.hpp>
#include <fockbuilder.hpp>
#include <fockbuilder/rofock.hpp>
#include <fockbuilder/fourcompfock.hpp>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/giaodirecteri.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/impl.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>

#include <singleslater/neoss.hpp>
#include <singleslater/multiparticless.hpp>

namespace ChronusQ {

  void resolveSubsystemGuessOptions(std::vector<QuantumSubsystem>& quantumSubsystems, const SCFControls& scfControls) {

    // Track the XXX_GUESS (in the SCF input section) that are used the subsystems,
    // and error out at the end when invalid XXX_GUESS are set by users
    std::set<std::string> recognizedGuessLabels;

    for(auto& sys : quantumSubsystems) {
      SCFControls subsystemSCFControls = scfControls;
      subsystemSCFControls.guessBasis = sys.guessBasis;

      // Electronic subsystem uses scfControls.guess (parsed from SCF/GUESS)
      // Other subsystem uses tight guess by default
      SS_GUESS subsystemGuess = sys.label == "E" ? scfControls.guess : TIGHT;
      // An input label such as QP applies to every runtime subsystem (QP0, QP1, QP2 ... when doing distinguisable)
      if(auto inputLabelGuess = scfControls.subsystemGuesses.find(sys.inputLabel);
          inputLabelGuess != scfControls.subsystemGuesses.end()) {
        subsystemGuess = inputLabelGuess->second;
        recognizedGuessLabels.insert(inputLabelGuess->first);
      }
      // An exact runtime label such as QP0 overrides the input-label selection
      if(auto runtimeLabelGuess = scfControls.subsystemGuesses.find(sys.label);
          runtimeLabelGuess != scfControls.subsystemGuesses.end()) {
        subsystemGuess = runtimeLabelGuess->second;
        recognizedGuessLabels.insert(runtimeLabelGuess->first);
      }
      if(subsystemGuess == NEOConvergeClassical and sys.label != "E")
        CErr("SCF/CLASSICAL guess is only valid for subsystem E; specify " +
             sys.inputLabel + "_GUESS for subsystem " + sys.label);
      // Non-electronic TIGHT guess uses the NEO-specific NEOTightParticle function
      subsystemSCFControls.guess = (subsystemGuess == TIGHT and sys.label != "E") ? NEOTightParticle : subsystemGuess;
      sys.ssOptions.scfControls = subsystemSCFControls;
    }

    // Reject input XXX_GUESS that did not match either an input or runtime subsystem label.
    for(const auto& [label, guess] : scfControls.subsystemGuesses)
      if(not recognizedGuessLabels.count(label))
        CErr("SCF guess specified for unknown quantum subsystem " + label);
  }

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  std::set<std::string> CQQM_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "REFERENCE",
      "JOB",
      "X2CTYPE",
      "SPINORBITSCALING",
      "ATOMICX2C",
      "SNSOTYPE",
      "IGNOREVPP",
      "DKSTYPE",
      "ONECENTERK"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  /**
   *
   *  Check valid keywords in the Proton QM section.
   *
  */
  std::set<std::string> CQPROTQM_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "REFERENCE",
      "IGNOREPROTONTWOBODY",
      "ONECENTERK",
      "ERFOMEGA",
      "DISTINGUISHABLE"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  std::set<std::string> CQQUANTUMSUBSYSTEMQM_VALID(const std::map<std::string, std::string>& inputSection) {
  
    std::set<std::string> allowedKeywords = {
      "REFERENCE",
      "ONECENTERK",
      "DISTINGUISHABLE",
      "PARTICLECHARGE",
      "PARTICLEMASS",
      "IGNOREPROTONTWOBODY"
    };
  
    return CQInvalidKeywords(allowedKeywords, inputSection);
  }


  /**
   *
   *  Check valid keywords in the section.
   *
  */
  std::set<std::string> CQDFTINT_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "EPS",
      "NANG",
      "NRAD",
      "NMACRO",
      "INHOUSE",
      "GAUXC",
      "BASISTOL"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  /**
   * \brief Parse the SingleSlater Referece information using
   * the input file.
   * 
   * \param [in]  out     Output device for data / error output.
   * \param [in]  tokens  Vector of string that contains reference information
   *
   * \returns RefOptions object that stores all reference options.
   *
   */
  RefOptions parseRef(std::ostream &out, 
    Molecule &mol, std::vector<std::string> &tokens) {

    // Initialize return
    RefOptions ref;

    // Determine the Real/Complex flag
    if( tokens.size() == 1 )      ref.RCflag = "AUTO";
    else if( tokens.size() == 2 ) ref.RCflag = tokens[0];
    else CErr("QM.REFERENCE Field not valid",out);


    // Kohn-Sham Keywords
    std::vector<std::string> KSRefs {
      "SLATER", 
      "B88",
      "LSDA",
      "SVWN5",
      "BLYP",
      "PBEXPBEC",
      "B3LYP",
      "B3PW91",
      "PBE0",
      "BHANDHLYP",
      "BHANDH",
      "CAMB3LYP",
      "HSE06",
      "LRCWPBE",
      "LCWPBE",
      "WB97",
      "WB97X",
      "LDA",
      //Gets definition of  Custom Functional from [GAUXC] section of input
      "CUSTOM"
    };

    std::vector<std::string> EPCRefs {
      // CQ allowed EPC functional
      "EPC17",
      "EPC19",
      // Gauxc allowed EPC functional
      "EPC17_1",
      "EPC17_2",
      "EPC18_1",
      "EPC18_2"
    };

    KSRefs.insert(KSRefs.begin(), EPCRefs.begin(), EPCRefs.end());

    // All reference keywords
    std::vector<std::string> rawRefs(KSRefs);
    rawRefs.insert(rawRefs.begin(),"HF");

    // Construct R/U/RO/G/X2C/4C reference keywords
    std::vector<std::string> RRefs, URefs, RORefs, GRefs, X2CRefs, TwoCRefs, FourCRefs;
    for(auto &f : rawRefs) {
      RRefs.emplace_back( "R" + f );
      URefs.emplace_back( "U" + f );
      GRefs.emplace_back( "G" + f );
      X2CRefs.emplace_back( "X2C" + f );
      TwoCRefs.emplace_back( "2C" + f );
      FourCRefs.emplace_back( "4C" + f );
    }
    RORefs.emplace_back( "ROHF" );
    

    // This is the reference string to be parsed
    std::string refString = tokens.back();

    // If dispersion is enabled, parse the d3 model after '-'
#ifdef CQ_HAS_D3
    std::string d3String;
    if( auto pos = refString.find('-'); pos != std::string::npos ) {
      d3String = refString.substr(pos + 1);
      refString  = refString.substr(0, pos);
    }
#endif    
   
    // Boolean for 2cHF specified as GHF in input
    bool isGHF = false;

    // Determine type of reference
    if ( std::find(rawRefs.begin(),rawRefs.end(),refString) != rawRefs.end() )
      ref.refType = isRawRef;
    else if ( std::find(RRefs.begin(),RRefs.end(),refString) != RRefs.end() )
      ref.refType = isRRef;
    else if ( std::find(URefs.begin(),URefs.end(),refString) != URefs.end() )
      ref.refType = isURef;
    else if ( std::find(RORefs.begin(),RORefs.end(),refString) != RORefs.end() )
      ref.refType = isRORef;
    else if ( std::find(GRefs.begin(),GRefs.end(),refString) != GRefs.end() ){
      ref.refType = isTwoCRef;
      isGHF = true;
    }else if ( std::find(TwoCRefs.begin(),TwoCRefs.end(),refString) != TwoCRefs.end() )
      ref.refType = isTwoCRef;
    else if ( std::find(X2CRefs.begin(),X2CRefs.end(),refString) != X2CRefs.end() ){
      ref.refType = isTwoCRef;
      ref.isX2CRef = true;
    }else if ( std::find(FourCRefs.begin(),FourCRefs.end(),refString) != FourCRefs.end() )
      ref.refType = isFourCRef;
    else 
      CErr(refString + " is not a valid QM.REFERENCE",out);


    // Cleanup the reference string
    if( ref.refType != isRawRef ) {
      if( ref.refType == isTwoCRef and ref.isX2CRef ) {
        refString.erase(0,3);
      } else if( ref.refType == isFourCRef or ( ref.refType == isTwoCRef and not isGHF ) ) {
        refString.erase(0,2);
      } else {
        refString.erase(0,1);
      }
    }
    // Handle KS related queries
    ref.isKSRef = 
      std::find(KSRefs.begin(),KSRefs.end(),refString) != KSRefs.end();

    if( ref.isKSRef )
      ref.funcName = refString;

    ref.isEPCRef =
      std::find(EPCRefs.begin(),EPCRefs.end(),refString) != EPCRefs.end();

    // Handle D3 dispersion parsing
#ifdef CQ_HAS_D3
    if( !d3String.empty() ) {
      
      if (ref.isEPCRef) CErr("D3 dispersion not supported for EPC references", out);

      // normalize suffix to upper
      std::string d3ModelString = d3String;
      std::transform(d3ModelString.begin(), d3ModelString.end(), d3ModelString.begin(), [](unsigned char c){ return std::toupper(c); });
      if(d3ModelString == "D3") d3ModelString = "D3BJ"; // For plain -D3, default to -D3BJ

      // Allowed D3 Models
      static const std::vector<std::string> d3Models { "D3BJ","D3ZERO","D3BJM","D3ZEROM","D3OP" };
      if( std::find(d3Models.begin(), d3Models.end(), d3ModelString) == d3Models.end() )
        CErr("QM.REFERENCE: unknown D3 Model '-" + d3String + "'", out);

      ref.useD3   = true;
      ref.d3ModelString = d3ModelString;

      ref.d3RefString = refString;
      if(ref.d3RefString == "PBEXPBEC") ref.d3RefString = "PBE";
      std::transform(ref.d3RefString.begin(), ref.d3RefString.end(), ref.d3RefString.begin(),
                     [](unsigned char c){ return std::tolower(c); });

    }
    std::cout << "      D3 Reference: " << ref.d3RefString   << std::endl;
    std::cout << "      D3 Model:     " << ref.d3ModelString << std::endl;
#endif

    // Raw reference
    if( ref.refType == isRawRef ) {
      out << "  *** Auto-determination of reference: " << refString << " -> ";
      ref.iCS = mol.multip == 1;

      if(ref.iCS){
        out << "R" << refString;
        ref.refType = isRRef;
      }else{
        out << "U" << refString;
        ref.refType = isURef;
      }

      out << " ***" << std::endl;
      
    } else if( ref.refType == isRRef )
      if( mol.multip != 1 )
        CErr("Spin-Restricted Reference only valid for singlet spin multiplicities",out);
      else
        ref.iCS = true;
    else if( ref.refType == isURef or ref.refType == isRORef )
      ref.iCS = false;
    else if( ref.refType == isTwoCRef ) {
      ref.iCS = false; ref.nC = 2;
    }
    else if( ref.refType == isFourCRef ) {
      ref.iCS = false; ref.nC = 4;
    }

    // Determine Real/Complex if need be
    if(not ref.RCflag.compare("AUTO") ) {
      if( ref.nC == 2 or ref.nC == 4 )
        ref.RCflag = "COMPLEX";
      else
        ref.RCflag = "REAL";

      out << "  *** Auto-determination of wave function field: AUTO -> " 
          << ref.RCflag << " ***" << std::endl;
    }

    out << "\n\n";

    return ref;
  }

  /**
   * \brief Construct a list of DFT Functional objects based on input name.
   * 
   * \param [in]  funcName   Input functional name
   * \param [out] funcList   Vector that stores constructed DFT functional object
   *
   */
  void buildFunclist(std::vector<std::shared_ptr<DFTFunctional>> &funcList,
    std::string funcName) {

    if( not funcName.compare("B88") )
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<BEightyEight>()
        )
      );

    if( not funcName.compare("SLATER") )
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<SlaterExchange>()
        )
      );

    if( not funcName.compare("LSDA") or not funcName.compare("LDA") ) {

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<SlaterExchange>()
        )
      );

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<VWNV>()
        )
      );

    }

    if( not funcName.compare("BLYP") ) {

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<BEightyEight>()
        )
      );

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<LYP>()
        )
      );

    }

    if( not funcName.compare("SVWN5") ) {

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<SlaterExchange>()
        )
      );

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<VWNV_G>()
        )
      );

    }

    if( not funcName.compare("PBEXPBEC") ) {

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<PBEX>()
        )
      );

      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<PBEC>()
        )
      );

    }

    if( not funcName.compare("B3LYP") ) 
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<B3LYP>()
        )
      );

    if( not funcName.compare("B3PW91") ) 
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<B3PW91>()
        )
      );

    if( not funcName.compare("PBE0") ) 
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<PBE0>()
        )
      );

    if( not funcName.compare("BHANDH") )
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<BHANDH>()
        )
      );

    if( not funcName.compare("BHANDHLYP") )
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<BHANDHLYP>()
        )
      );

    if (not funcName.compare("EPC17") or not funcName.compare("EPC17_2"))
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<EPC17>("EPC-17")
        )
      );

    if (not funcName.compare("EPC19"))
      funcList.push_back(
        std::dynamic_pointer_cast<DFTFunctional>(
          std::make_shared<EPC19>("EPC-19")
        )
      );
  }

  /**
   * \brief Parse the Integration information using
   * the input file.
   * 
   * \param [in]  out       Output device for data / error output.
   * \param [in]  input     Input file datastructure
   * \param [out] intParam  Object that stores parsed info
   *
   */
  void parseIntParam(std::ostream &out, CQInputFile &input, 
    IntegrationParam &intParam) {

    if( input.containsSection("DFTINT") ) {

      OPTOPT( intParam.epsilon      = input.getData<double>("DFTINT/EPS")  );
      OPTOPT( intParam.nAng         = input.getData<size_t>("DFTINT/NANG") );
      OPTOPT( intParam.nRad         = input.getData<size_t>("DFTINT/NRAD") );
      OPTOPT( intParam.nRadPerBatch = input.getData<size_t>("DFTINT/NMACRO") );
      bool gauFlag1(false), gauFlag2(false);
      OPTOPT( gauFlag1              = not input.getData<bool>("DFTINT/INHOUSE") );
      OPTOPT( gauFlag2              = input.getData<bool>("DFTINT/GAUXC") ); 
      intParam.useGauXC = gauFlag1 or gauFlag2;
      OPTOPT( intParam.basisTol = input.getData<double>("DFTINT/BASISTOL")  );

    }

    if( not intParam.useGauXC ){
    out << "\nDFT Integration Settings:\n" << BannerTop << "\n\n" ;
    out << std::left;

    out << "  " << std::setw(28) << "Screening Tolerance:";
    out << intParam.epsilon << std::endl;

    out << "  " << std::setw(28) << "Angular Grid:";
    out <<  "Lebedev (" << intParam.nAng << ")" << std::endl;
    out << "  " << std::setw(28) << "Radial Grid:";
    out <<  "Euler-Maclaurin (" << intParam.nRad << ")" << std::endl;
    out << "  " << std::setw(28) << "Macro Batch Size:";
    out <<  intParam.nRadPerBatch << " Radial Points" << std::endl;
    out << "  " << std::setw(28) << "DFT Engine:";
    out <<  (intParam.useGauXC ? "GauXC" : "In-house") << " DFT Engine" << std::endl;

    out << std::endl << BannerEnd << std::endl;
    } else {
    out << "\nDFT Integration Settings:\n" << BannerTop << "\n\n" ;
    out << std::left;

    out <<  (intParam.useGauXC ? "GauXC" : "In-house") << " DFT Engine" << std::endl;

    out << std::endl << BannerEnd << std::endl;
    }
  }

  /**
   *  \brief Parse Hamiltonian options using the input 
   *  file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] basis  Basis Set
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *  \param [in]  refOptions Object that stores reference info
   *  \param [out] hamiltonianOptions Objects that get parsed
   *
   *
   */
  void parseHamiltonianOptions(std::ostream &out, CQInputFile &input, 
    BasisSet &basis, RefOptions &refOptions, HamiltonianOptions &hamiltonianOptions, std::string section) {

    // Parse hamiltonianOptions
    hamiltonianOptions.basisType = basis.basisType;

    std::string X;

    // Parse X2C option
    // X2CType = off (default), spinfree, onee, twoe
    X = "DEFAULT";
    OPTOPT( X = input.getData<std::string>(section + "/X2CTYPE")  );
    trim(X);
    if( not X.compare("SPINFREE") ) {

      hamiltonianOptions.x2cType = X2C_TYPE::ONEE;
      hamiltonianOptions.OneEScalarRelativity = true;
      hamiltonianOptions.OneESpinOrbit = false;
      hamiltonianOptions.SNSO = false;
      hamiltonianOptions.AtomicMeanField = false;

    } else if( not X.compare("FOCK")) {

      hamiltonianOptions.x2cType = X2C_TYPE::FOCK;
      hamiltonianOptions.OneEScalarRelativity = true;
      hamiltonianOptions.OneESpinOrbit = true;
      hamiltonianOptions.SNSO = false;
      hamiltonianOptions.AtomicMeanField = false;

    } else if( not X.compare("ONEE") or not X.compare("ONEELECTRON")
               or ( not X.compare("DEFAULT") and refOptions.isX2CRef ) ) {
      // Legacy X2C- reference is equilvalent to 2C- reference + OneE-X2C

      hamiltonianOptions.x2cType = X2C_TYPE::ONEE;
      hamiltonianOptions.OneEScalarRelativity = true;
      hamiltonianOptions.OneESpinOrbit = true;
      hamiltonianOptions.SNSO = true;
      hamiltonianOptions.AtomicMeanField = false;

    } else if( not X.compare("TWOE") or not X.compare("TWOELECTRON")) {

      hamiltonianOptions.x2cType = X2C_TYPE::TWOE;
      CErr(X + " NYI",out);

    } else if( not X.compare("OFF")
               or ( not X.compare("DEFAULT") and not refOptions.isX2CRef ) ) {

      if ( refOptions.refType == isFourCRef ) {

        hamiltonianOptions.OneEScalarRelativity = true;
        hamiltonianOptions.OneESpinOrbit = true;

      } else {

        hamiltonianOptions.OneEScalarRelativity = false;
        hamiltonianOptions.OneESpinOrbit = false;

      }

      hamiltonianOptions.SNSO = false;
      hamiltonianOptions.AtomicMeanField = false;

    } else  {

      CErr(X + " not a valid " + section + "/X2CTYPE",out);

    }



    // Parse one-electron spin-orbie scaling option
    // SpinOrbitScaling  = noscaling, boettger (dafault), atomicmeanfield (amfi)
    X = "DEFAULT"; // Unspecified value
    OPTOPT( X = input.getData<std::string>(section + "/SPINORBITSCALING")  );
    trim(X);
    if( not X.compare("NOSCALING") ) {

      hamiltonianOptions.SNSO = false;
      hamiltonianOptions.AtomicMeanField = false;

    } else if( not X.compare("SNSO") ) {

      if( not hamiltonianOptions.OneESpinOrbit ) 
        CErr("Spin-Orbit Scaling = "+ X + " is not compatible with X2CType = SpinFree",out);
      hamiltonianOptions.SNSO = true;
      hamiltonianOptions.AtomicMeanField = false;

    } else if( not X.compare("AMFI") or not X.compare("ATOMICMEANFIELD")) {

      if( not hamiltonianOptions.OneESpinOrbit ) 
        CErr("Spin-Orbit Scaling = "+ X + " is not compatible with X2CType = SpinFree",out);
      hamiltonianOptions.SNSO = false;
      hamiltonianOptions.AtomicMeanField = true;
      CErr("AMFI NYI!",out);

    } else if( not X.compare("DEFAULT") ) {

      if ( hamiltonianOptions.OneESpinOrbit
           and refOptions.refType != isFourCRef
           and hamiltonianOptions.x2cType != X2C_TYPE::FOCK) {

        hamiltonianOptions.SNSO = true;
        hamiltonianOptions.AtomicMeanField = false;

      } else {

        hamiltonianOptions.SNSO = false;
        hamiltonianOptions.AtomicMeanField = false;

      }

    } else {

      CErr(X + " not a valid " + section + "/SPINORBITSCALING",out);

    }

    // Parse screened nuclear spin–orbit approximation
    X = "BOETTGER"; // Unspecified value
    OPTOPT( X = input.getData<std::string>(section + "/SNSOTYPE")  );
    trim(X);
    if( not X.compare("BOETTGER") ) {

      hamiltonianOptions.snsoType = SNSO_TYPE::BOETTGER;

    } else if( not X.compare("DC") ) {

      hamiltonianOptions.snsoType = SNSO_TYPE::DC;

    } else if( not X.compare("DCB") ) {

      hamiltonianOptions.snsoType = SNSO_TYPE::DCB;

    } else if( not X.compare("ROW_DEP_DCB")) {

      hamiltonianOptions.snsoType = SNSO_TYPE::ROW_DEP_DCB;

    } else {

      CErr(X + " not a valid " + section + "/SNSOTYPE",out);

    }

    // Parse Atomic X2C
    hamiltonianOptions.AtomicX2C = parseAtomicType(out,input,hamiltonianOptions.AtomicX2CType,section);


    // Parse Finite Width Nuclei
    std::string finiteCore = "DEFAULT";
    OPTOPT( finiteCore = input.getData<std::string>("INTS/FINITENUCLEI"); )
    trim(finiteCore);
    if( not finiteCore.compare("TRUE") )
      hamiltonianOptions.finiteWidthNuc = true;
    else if( not finiteCore.compare("FALSE") )
      hamiltonianOptions.finiteWidthNuc = false;
    else if( not finiteCore.compare("DEFAULT") )
      hamiltonianOptions.finiteWidthNuc = refOptions.refType == isFourCRef or hamiltonianOptions.OneEScalarRelativity;
    else
      CErr(finiteCore + " not a valid INTS.ALG",out);


    // Parse Integral library
    OPTOPT( hamiltonianOptions.Libcint = input.getData<bool>("INTS/LIBCINT") )

    if (hamiltonianOptions.Libcint and basis.forceCart)
      CErr("Libcint + cartesian GTO NYI.");


    OPTOPT( hamiltonianOptions.BareCoulomb = input.getData<bool>("INTS/BARECOULOMB") )
    OPTOPT( hamiltonianOptions.BareCoulomb = input.getData<bool>("INTS/LLLL") )
    //OPTOPT( hamiltonianOptions.DiracCoulombSSSS = input.getData<bool>("INTS/SSSS") )
    //OPTOPT( hamiltonianOptions.DiracCoulomb = input.getData<bool>("INTS/DIRACCOULOMB") )
    //OPTOPT( hamiltonianOptions.Gauge = input.getData<bool>("INTS/GAUGE") )
    //OPTOPT( hamiltonianOptions.Gaunt = input.getData<bool>("INTS/GAUNT") )
    //try{
    //  if ( input.getData<bool>("INTS/BREIT") ) {
    //    hamiltonianOptions.DiracCoulomb = true;
    //    hamiltonianOptions.DiracCoulombSSSS = true;
    //    hamiltonianOptions.Gaunt = true;
    //    hamiltonianOptions.Gauge = true;
    //  }
    //} catch(...) {}

    hamiltonianOptions.DiracCoulomb = false;
    hamiltonianOptions.DiracCoulombSSSS = false;
    hamiltonianOptions.Gaunt = false;
    hamiltonianOptions.Gauge = false;

    // Parse 4C options
    // Dirac-Coulomb
    try { 
      std::string DCOptions = "FALSE";
      try {
        DCOptions = input.getData<std::string>("INTS/DIRACCOULOMB");
      } catch (...) {
        DCOptions = input.getData<std::string>("INTS/DC");
      }
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      if( std::regex_search(DCOptions, regexTRUE) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombType = TYPE_4C::All;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(DCOptions, regexFALSE) ){
        hamiltonianOptions.DiracCoulomb = false;
      }

      if( std::regex_search(DCOptions, regexALL) ){
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombType = TYPE_4C::All;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(DCOptions, regexSF) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombType = TYPE_4C::SpinFree;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(DCOptions, regexSD) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombType = TYPE_4C::SpinDependent;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::None;
      }

      if( std::regex_search(DCOptions, regex3C) ){
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
      } else if ( std::regex_search(DCOptions, regex2C) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
      } else if ( std::regex_search(DCOptions, regex1C) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
      } else if ( std::regex_search(DCOptions, regexAMF) ) {
        hamiltonianOptions.DiracCoulomb = true;
        hamiltonianOptions.DiracCoulombApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
      }
  
    } catch(...) {}

    // by default, DC includes SSSS unless "false" is set upon input
    if(hamiltonianOptions.DiracCoulomb) hamiltonianOptions.DiracCoulombSSSS = true;
   
    // SSSS
    try { 
      std::string SSSSOptions = "FALSE";
      SSSSOptions = input.getData<std::string>("INTS/SSSS");
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      if( std::regex_search(SSSSOptions, regexTRUE) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSType = TYPE_4C::All;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(SSSSOptions, regexFALSE) ){
        hamiltonianOptions.DiracCoulombSSSS = false;
      }

      if( std::regex_search(SSSSOptions, regexALL) ){
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSType = TYPE_4C::All;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(SSSSOptions, regexSF) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSType = TYPE_4C::SpinFree;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(SSSSOptions, regexSD) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSType = TYPE_4C::SpinDependent;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::None;
      }

      if( std::regex_search(SSSSOptions, regex3C) ){
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
      } else if ( std::regex_search(SSSSOptions, regex2C) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
      } else if ( std::regex_search(SSSSOptions, regex1C) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
      } else if ( std::regex_search(SSSSOptions, regexAMF) ) {
        hamiltonianOptions.DiracCoulombSSSS = true;
        hamiltonianOptions.SSSSApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
      }
  
    } catch(...) {}

    bool hasGauntInput = false;
    // Gaunt
    try { 
      std::string GauntOptions = "FALSE";
      GauntOptions = input.getData<std::string>("INTS/GAUNT");
      hasGauntInput = true;
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      auto const regexSCALE = std::regex(R"([+-]?(([0-9]+(\.[0-9]*)?)|(\.[0-9]+))([eE][+-]?[0-9]+)?)");

      std::smatch scale;

      if( std::regex_search(GauntOptions, regexTRUE) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
       } else if ( std::regex_search(GauntOptions, scale, regexSCALE)){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;

        double scale_ = std::stod(scale.str(0));
        hamiltonianOptions.GauntScale = scale_;  
      } else if ( std::regex_search(GauntOptions, regexFALSE) ){
        hamiltonianOptions.Gaunt = false;
      }

      if( std::regex_search(GauntOptions, regexALL) ){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(GauntOptions, regexSF) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::SpinFree;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(GauntOptions, regexSD) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::SpinDependent;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
      }

      if( std::regex_search(GauntOptions, regex3C) ){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
      } else if ( std::regex_search(GauntOptions, regex2C) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
      } else if ( std::regex_search(GauntOptions, regex1C) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
      } else if ( std::regex_search(GauntOptions, regexAMF) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
      }
  
    } catch(...) {}


    bool hasGaugeInput = false;
    // Gauge
    try { 
      std::string GaugeOptions = "FALSE";
      GaugeOptions = input.getData<std::string>("INTS/GAUGE");
      hasGaugeInput = true;
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      auto const regexSCALE = std::regex(R"([+-]?(([0-9]+(\.[0-9]*)?)|(\.[0-9]+))([eE][+-]?[0-9]+)?)");

      std::smatch scale;

      if( std::regex_search(GaugeOptions, regexTRUE) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(GaugeOptions, scale, regexSCALE)){
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;

        double scale_ = std::stod(scale.str(0));
        hamiltonianOptions.GaugeScale = scale_;
      } else if ( std::regex_search(GaugeOptions, regexFALSE) ){
        hamiltonianOptions.Gauge = false;
      }

      if( std::regex_search(GaugeOptions, regexALL) ){
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(GaugeOptions, regexSF) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::SpinFree;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(GaugeOptions, regexSD) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::SpinDependent;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      }

      if( std::regex_search(GaugeOptions, regex3C) ){
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
      } else if ( std::regex_search(GaugeOptions, regex2C) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
      } else if ( std::regex_search(GaugeOptions, regex1C) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
      } else if ( std::regex_search(GaugeOptions, regexAMF) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
      }
  
    } catch(...) {}

    // Breit = 1/2 Gaunt + gauge
    try { 
      std::string BreitOptions = "FALSE";
      BreitOptions = input.getData<std::string>("INTS/BREIT");

      if (hasGauntInput) {
        CErr("INT/BREIT and INT/GAUNT options may conflict. Please only set one of them.");
      }
      if (hasGaugeInput) {
        CErr("INT/BREIT and INT/GAUGE options may conflict. Please only set one of them.");
      }

      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      auto const regexSCALE = std::regex(R"([+-]?(([0-9]+(\.[0-9]*)?)|(\.[0-9]+))([eE][+-]?[0-9]+)?)");

      std::smatch scale;

      if( std::regex_search(BreitOptions, regexTRUE) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(BreitOptions, scale, regexSCALE)){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;

        double scale_ = std::stod(scale.str(0));
        hamiltonianOptions.GauntScale = scale_;  
        hamiltonianOptions.GaugeScale = scale_;
        
      } else if ( std::regex_search(BreitOptions, regexFALSE) ){
        hamiltonianOptions.Gaunt = false;
        hamiltonianOptions.Gauge = false;
      }

      if( std::regex_search(BreitOptions, regexALL) ){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(BreitOptions, regexSF) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::SpinFree;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::SpinFree;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      } else if ( std::regex_search(BreitOptions, regexSD) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::SpinDependent;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::SpinDependent;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
      }

      if( std::regex_search(BreitOptions, regex3C) ){
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::ThreeCenter;
      } else if ( std::regex_search(BreitOptions, regex2C) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::TwoCenter;
      } else if ( std::regex_search(BreitOptions, regex1C) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::OneCenter;
      } else if ( std::regex_search(BreitOptions, regexAMF) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::AtomicMeanField;
      }
  
    } catch(...) {}

    if (refOptions.refType != isFourCRef
        and hamiltonianOptions.x2cType != X2C_TYPE::FOCK) {

      hamiltonianOptions.BareCoulomb = false;
      hamiltonianOptions.DiracCoulomb = false;
      hamiltonianOptions.DiracCoulombSSSS = false;
      hamiltonianOptions.Gaunt = false;
      hamiltonianOptions.Gauge = false;

    }

    if ((hamiltonianOptions.Gauge or hamiltonianOptions.DiracCoulombSSSS)
        and not hamiltonianOptions.Libcint)
      CErr("4C Gauge and SSSS terms NYI with libint. "
            "Please use libcint = true instead.", out);

    if(refOptions.refType == isFourCRef){
      hamiltonianOptions.Libcint=true;
      hamiltonianOptions.x2cType = X2C_TYPE::OFF;
    }

    hamiltonianOptions.updateGaunt = hamiltonianOptions.Gaunt;
    hamiltonianOptions.updateGauge = hamiltonianOptions.Gauge;

    // DKS
    X = "DEFAULT";
    OPTOPT( X = input.getData<std::string>(section + "/DKSTYPE")  );
    trim(X);
    if( not X.compare("VLL") ) {
      if( not ( refOptions.isKSRef and (refOptions.refType == isFourCRef ))){
    	CErr("DKS requires 4C + a DFT Functional", out);
      } else {
	  hamiltonianOptions.dksType = DKS_TYPE::VLL;
      } 
    } else if( not X.compare("FULL") ) {
      if( not ( refOptions.isKSRef and (refOptions.refType == isFourCRef ))){
    	CErr("DKS requires 4C + a DFT Functional", out);
      } else {
	  hamiltonianOptions.dksType = DKS_TYPE::FULL;
      } 
    } else if ( not X.compare("DEFAULT") ) {
      if( not ( refOptions.isKSRef and (refOptions.refType == isFourCRef ))){
   	hamiltonianOptions.dksType = DKS_TYPE::OFF;
      } else {
          hamiltonianOptions.dksType = DKS_TYPE::FULL;
      }
    } else { 
	CErr(X + " not a valid " + section + "/DKSTYPE",out);
    } //DKS

    // For NEO (and in particular post-NEO-HF methods)
    OPTOPT(hamiltonianOptions.ignoreProtonTwoBody = input.getData<bool>(section + "/IGNOREPROTONTWOBODY"));

    // For RI J/K contraction with 3-index ERI
    OPTOPT(hamiltonianOptions.oneCenterK = input.getData<bool>(section + "/ONECENTERK"));

  }

  /**
   *  \brief Parse atomic X2C options using the input 
   *  file.
   *
   *  \param [in]  out            Output device for data / error output.
   *  \param [in]  input          Input file datastructure
   *  \param [out] atomicX2CType  Objects that get parsed
   *
   *  \returns boolean that tells whether atomic X2C is used.
   *
   */
  bool parseAtomicType(std::ostream &out, CQInputFile &input, 
    ATOMIC_X2C_TYPE &atomicX2CType, std::string section) {

    // Parse Atomic X2C option
    // AtomicX2C  = ALH, ALU, DLH, DLU, OFF (default)
    bool atomic = false;
    std::string X = "OFF";
    OPTOPT( X = input.getData<std::string>(section + "/ATOMICX2C")  );
    trim(X);
    if( not X.compare("ALH") ) {
      atomic = true;
      atomicX2CType = {true,true};
    } else if( not X.compare("ALU") ) {
      atomic = true;
      atomicX2CType = {true,false};
    } else if( not X.compare("DLH") ) {
      atomic = true;
      atomicX2CType = {false,true};
    } else if( not X.compare("DLU") ) {
      atomic = true;
      atomicX2CType = {false,false};
    } else if( not X.compare("OFF") ){
      atomic = false;
    } else {
      CErr(X + " not a valid " + section + "/ATOMICX2C",out);
    }

    return atomic;

  }


  /**
   *  \brief Construct a SingleSlaterOptions object using the input
   *  file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns a SingleSlaterOptions object
   *    constructed from the input options.
   *
   */
  SingleSlaterOptions getSingleSlaterOptions(
      std::ostream &out, CQInputFile &input,
      Molecule &mol, BasisSet &basis,
      Particle p, std::string section) {

    out << "  *** Parsing " << section << "/REFERENCE options ***\n";

    SingleSlaterOptions options;

    // Attempt to find reference
    std::string reference;
    try { 
      reference = input.getData<std::string>(section + "/REFERENCE");
    } catch(...) {
      CErr(section + "/REFERENCE Keyword not found!",out);
    }

    // Digest reference string
    // Trim Spaces
    trim(reference);

    // Split into tokens
    std::vector<std::string> tokens;
    split(tokens,reference);
    for(auto &X : tokens) trim(X);

    // Parse reference information
    options.refOptions = parseRef(out,mol,tokens);


    // FIXME: Should put this somewhere else
    // Parse KS integration
    if( options.refOptions.isKSRef ){
      
      parseIntParam(out, input, options.intParam);

      if( GauXCUtils::is_range_separated(options.refOptions.funcName) and not options.intParam.useGauXC )
        CErr(options.refOptions.funcName + " requires the GauXC DFT engine. Set DFTINT/GAUXC = TRUE.", out);
      
      // Determine EPC functional is available in CQ/GauXC, and edit keywords in needed
      if(options.refOptions.isEPCRef){
        bool isHMass = true;
        for(const auto & atomQIndex : mol.atomsQ)
        {
          if(mol.atoms[atomQIndex].atomicMass != atomicReference["H-1"].atomicMass)
            bool isHMass = false;
        }
        if(!options.refOptions.funcName.compare("ECP19") && !isHMass)
          CErr("EPC19 currently hardcoded for H-1, while H-2 was requested!");

        if(not options.intParam.useGauXC){
          if(!options.refOptions.funcName.compare("EPC17_1"))
            CErr("EPC17_1 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("EPC18_1"))
            CErr("EPC18_1 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("EPC18_2"))
            CErr("EPC18_2 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("LDA"))
            CErr("LDA Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("EPC17_2"))
            options.refOptions.funcName = "EPC17"; 
        } else{
          if(!options.refOptions.funcName.compare("EPC19"))
            CErr("EPC19 Not Yet Implemented In GauXC");
          else if(!options.refOptions.funcName.compare("EPC17"))
            options.refOptions.funcName = "EPC17_2"; 
        } 
      }
    }



    // Parse hamiltonianOptions
    parseHamiltonianOptions(out,input,basis,options.refOptions,options.hamiltonianOptions,section);

    // Error checking here to have access to the particle set
    if(mol.nTotalP > 1 && options.hamiltonianOptions.ignoreProtonTwoBody)
    {
      //CErr("Requested turning of Proton Two Body Interaction with >1 Quantum Proton is illegal!");
    }

    options.hamiltonianOptions.particle = p;
    options.hamiltonianOptions.particle2 = p;

    out << options.hamiltonianOptions << std::endl;

    return options;

  }


  /**
   *  \brief Construct a SingleSlater object from the options.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns shared_ptr to a SingleSlaterBase object
   *    constructed from the input options.
   *
   */ 
  std::shared_ptr<SingleSlaterBase>
  SingleSlaterOptions::buildSingleSlater(
      std::ostream &out,
      Molecule &mol, BasisSet &basis,
      std::shared_ptr<IntegralsBase> aoints,
      const IntegralOptions &aoints_options) const {


    // Build Functional List
    std::vector<std::shared_ptr<DFTFunctional>> funcList;
    if( refOptions.isKSRef )
      buildFunclist(funcList, refOptions.funcName);

    // Sanity Checks
    bool isGIAO = basis.basisType == COMPLEX_GIAO;

    if (hamiltonianOptions.OneESpinOrbit and not refOptions.RCflag.compare("REAL") )
      CErr("Real + Spin-orbit calculation not valid",out);

    if( isGIAO and not refOptions.RCflag.compare("REAL") )
      CErr("Real + GIAO not valid",out);

    //if( isGIAO and refOptions.isKSRef )
    //  CErr("KS + GIAO  NYI!",out);


    // Override core hamiltoninan type for X2C

    Particle p = hamiltonianOptions.particle;
    std::optional<size_t> nParticleOverride = hamiltonianOptions.nParticleOverride;

    if( p.charge < 0 && refOptions.isEPCRef )
      CErr("EPC functionals only valid on proton references!");

    if( p.charge >= 0 && refOptions.isKSRef && !refOptions.isEPCRef )
      CErr("Proton Kohn Sham references require EPC functionals");

  #define KS_LIST(T) \
    refOptions.funcName,funcList,MPI_COMM_WORLD,intParam,mol,basis,std::dynamic_pointer_cast<Integrals<T>>(aoints),refOptions.nC,refOptions.iCS,p,nParticleOverride

  #define HF_LIST(T) \
    MPI_COMM_WORLD,mol,basis,std::dynamic_pointer_cast<Integrals<T>>(aoints),refOptions.nC,refOptions.iCS,p,nParticleOverride

    // Construct the SS object
    std::shared_ptr<SingleSlaterBase> ss;

    if( not refOptions.RCflag.compare("REAL") )
      if( refOptions.isKSRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<double,double>>( KS_LIST(double) )
          );
      else if(refOptions.refType == isRORef)
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<double,double>>(
            "Real Restricted Open-shell Hartree-Fock", "R-ROHF", HF_LIST(double) )
          );
      else if( not isGIAO )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<double,double>>( HF_LIST(double) )
          );
      else
        CErr("GIAO + REAL is not a valid option.",out);

    else if( not refOptions.RCflag.compare("COMPLEX") and not isGIAO )
      if( refOptions.isKSRef and refOptions.isX2CRef)
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex,double>>(
              "Exact Two Component", "X2C-", KS_LIST(double)
            )
          );
      else if( refOptions.isKSRef and (refOptions.refType == isFourCRef))
         ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex,double>>(
              "Four Component","4C-", KS_LIST(double)
            )
          );
      else if( refOptions.isKSRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex,double>>( KS_LIST(double) )
          );
      else if( refOptions.isX2CRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,double>>(
              "Exact Two Component Hartree-Fock","X2C-HF",HF_LIST(double)
            )
          );
      else if( refOptions.refType == isRORef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,double>>(
              "Complex Restricted Open-shell Hartree-Fock", "C-ROHF", HF_LIST(double)
            )
          );
      else if( refOptions.refType == isFourCRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,double>>(
              "Four Component","4C-HF",HF_LIST(double)
            )
          );
      else // isGRef or isTwoCRef
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,double>>( HF_LIST(double) )
          );
    else
      if( refOptions.isKSRef and refOptions.isX2CRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex,dcomplex>>(
              "Exact Two Component", "X2C-", KS_LIST(dcomplex)
            )
          );
      else if( refOptions.isKSRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<KohnSham<dcomplex,dcomplex>>( KS_LIST(dcomplex) )
          );
      else if( refOptions.isX2CRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,dcomplex>>(
              "Exact Two Component","X2C-HF",HF_LIST(dcomplex)
            )
          );
      else if( refOptions.refType == isRORef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,dcomplex>>(
              "Complex Restricted Open-shell Hartree-Fock", "C-ROHF", HF_LIST(dcomplex)
            )
          );
      else if( refOptions.refType == isFourCRef )
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,dcomplex>>(
              "Four Component","4C-HF",HF_LIST(dcomplex)
            )
          );
      else // isGRef or isTwoCRef
        ss = std::dynamic_pointer_cast<SingleSlaterBase>(
            std::make_shared<HartreeFock<dcomplex,dcomplex>>( HF_LIST(dcomplex) )
          );


    // update IntegralsBase options
    aoints->options_ = hamiltonianOptions;

    // Construct CoreHBuilder
    if( refOptions.refType == isFourCRef ) {

      if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {
        p->coreHBuilder = std::make_shared<FourComponent<double,double>>(
            *std::dynamic_pointer_cast<Integrals<double>>(aoints), hamiltonianOptions);

        p->fockBuilder = std::make_shared<FourCompFock<double,double>>(hamiltonianOptions);

//        CErr("4C + Real WFN is not a valid option",std::cout);
      } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

        std::shared_ptr<Integrals<double>> aoints_double = std::dynamic_pointer_cast<Integrals<double>>(aoints);

        p->coreHBuilder = std::make_shared<FourComponent<dcomplex,double>>(
            *std::dynamic_pointer_cast<Integrals<double>>(aoints), hamiltonianOptions);

        p->fockBuilder = std::make_shared<FourCompFock<dcomplex,double>>(hamiltonianOptions);

      } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

        CErr("4C + GIAO NYI",std::cout);

      } else {

        CErr("Complex INT + Real WFN is not a valid option",std::cout);

      }
    } else if( hamiltonianOptions.OneEScalarRelativity ) {

      if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

        if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<dcomplex,double>>(hamiltonianOptions);
        else p->fockBuilder = std::make_shared<FockBuilder<dcomplex,double>>(hamiltonianOptions);

      } else if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

        if (not hamiltonianOptions.OneESpinOrbit) {

          if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<double,double>>(hamiltonianOptions);
          else p->fockBuilder = std::make_shared<FockBuilder<double,double>>(hamiltonianOptions);

        } else

          CErr("OneE-X2C-SpinOrbit + Real WFN is not a valid option",std::cout);

      } else if (auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

        if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<dcomplex,dcomplex>>(hamiltonianOptions);
        else p->fockBuilder = std::make_shared<FockBuilder<dcomplex,dcomplex>>(hamiltonianOptions);

      } else {

        CErr("Complex INT + Real WFN is not a valid option",std::cout);

      }
    } else {

      if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

        p->coreHBuilder = std::make_shared<NRCoreH<double,double>>(
            *std::dynamic_pointer_cast<Integrals<double>>(aoints), hamiltonianOptions);

        if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<double,double>>(hamiltonianOptions);
        else p->fockBuilder = std::make_shared<FockBuilder<double,double>>(hamiltonianOptions);

      } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

        p->coreHBuilder = std::make_shared<NRCoreH<dcomplex,double>>(
            *std::dynamic_pointer_cast<Integrals<double>>(aoints), hamiltonianOptions);

        if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<dcomplex,double>>(hamiltonianOptions);
        else p->fockBuilder = std::make_shared<FockBuilder<dcomplex,double>>(hamiltonianOptions);

      } else if (auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

        p->coreHBuilder = std::make_shared<NRCoreH<dcomplex,dcomplex>>(
            *std::dynamic_pointer_cast<Integrals<dcomplex>>(aoints), hamiltonianOptions);

        if(refOptions.refType == isRORef) p->fockBuilder = std::make_shared<ROFock<dcomplex,dcomplex>>(hamiltonianOptions);
        else p->fockBuilder = std::make_shared<FockBuilder<dcomplex,dcomplex>>(hamiltonianOptions);

      } else {
        CErr("Complex INT + Real WFN is not a valid option",std::cout);
      }
    }



    // Construct ERIContractions
    if(refOptions.refType == isFourCRef) {


      if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

        std::shared_ptr<TwoPInts<double>> &TPI =
            std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

        if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<double>>(TPI)) {

          TPI = std::make_shared<InCoreRelERI<double>>(basis.nBasis,
              hamiltonianOptions.DiracCoulomb, hamiltonianOptions.Gaunt,
              hamiltonianOptions.DiracCoulombSSSS, hamiltonianOptions.Gauge);

          p->TPI = std::make_shared<InCoreRelERIContraction<double,double>>(TPI);

        } else if (auto tpi_typed = std::dynamic_pointer_cast<InCoreRITPI<double>>(TPI)) {

          auto cdRelTPI = std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(tpi_typed);
          if (not cdRelTPI)
            CErr("Only Cholesky-type 4-component RI implemented.", std::cout);

          std::shared_ptr<InCoreRelERI<double>> relTPI =
              std::make_shared<InCoreRelERI<double>>(cdRelTPI, hamiltonianOptions, aoints_options.cdriintsoptions);

          TPI = relTPI;

          p->TPI = std::make_shared<InCoreRelERIContraction<double,double>>(TPI);

        } else if (auto tpi_typed = std::dynamic_pointer_cast<DirectTPI<double>>(TPI)) {

          p->TPI = std::make_shared<GTODirectRelERIContraction<double,double>>(tpi_typed);

        } else if (TPI) {
          CErr("Invalid TPInts type for Four-component Wavefunction<double,double>",std::cout);
        }

        p->TPI->printContractionTiming = scfControls.printContractionTiming;

      } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

        std::shared_ptr<TwoPInts<double>> &TPI =
            std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

        if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<double>>(TPI)) {

          TPI = std::make_shared<InCoreRelERI<double>>(basis.nBasis,
              hamiltonianOptions.DiracCoulomb, hamiltonianOptions.Gaunt,
              hamiltonianOptions.DiracCoulombSSSS, hamiltonianOptions.Gauge);

          p->TPI = std::make_shared<InCoreRelERIContraction<dcomplex,double>>(TPI);

        } else if (auto tpi_typed = std::dynamic_pointer_cast<InCoreRITPI<double>>(TPI)) {

          auto cdRelTPI = std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(tpi_typed);
          if (not cdRelTPI)
            CErr("Only Cholesky-type 4-component RI implemented.", std::cout);

          std::shared_ptr<InCoreRelERI<double>> relTPI =
              std::make_shared<InCoreRelERI<double>>(cdRelTPI, hamiltonianOptions, aoints_options.cdriintsoptions);

          TPI = relTPI;

          p->TPI = std::make_shared<InCoreRelERIContraction<dcomplex,double>>(TPI);

        } else if (auto tpi_typed = std::dynamic_pointer_cast<DirectTPI<double>>(TPI)) {

          p->TPI = std::make_shared<GTODirectRelERIContraction<dcomplex,double>>(tpi_typed);

        } else if (TPI) {
          CErr("Invalid TPInts type for Four-component Wavefunction<dcomplex,double>",std::cout);
        }
        
        p->TPI->printContractionTiming = scfControls.printContractionTiming;

      } else if (auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

        CErr("Complex INT Four-component Wavefunction method NYI",std::cout);

      } else {

        CErr("Complex INT + Real WFN is not a valid option",std::cout);

      }
    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

      std::shared_ptr<TwoPInts<double>> TPI =
          std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

      p->TPI = makeTPIContraction<double,double>(TPI);

      p->TPI->printContractionTiming = scfControls.printContractionTiming;

      if (hamiltonianOptions.oneCenterK) {
        p->TPI->setOneCenterK(true);
        p->TPI->setMapCen2BfSt(basis.mapCen2BfSt);
      }
    
    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

      std::shared_ptr<TwoPInts<double>> TPI =
          std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

      p->TPI = makeTPIContraction<dcomplex,double>(TPI);

      p->TPI->printContractionTiming = scfControls.printContractionTiming;

      if (hamiltonianOptions.oneCenterK) {
        p->TPI->setOneCenterK(true);
        p->TPI->setMapCen2BfSt(basis.mapCen2BfSt);
      }
    
    } else if (auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

      std::shared_ptr<TwoPInts<dcomplex>> TPI =
          std::dynamic_pointer_cast<Integrals<dcomplex>>(aoints)->TPI;

      if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<dcomplex>>(TPI)) {

        p->TPI = std::make_shared<InCore4indexTPIContraction<dcomplex,dcomplex>>(tpi_typed);

      } else if (auto tpi_typed = std::dynamic_pointer_cast<DirectTPI<dcomplex>>(TPI)) {

        p->TPI = std::make_shared<GIAODirectERIContraction>(tpi_typed);

      } else {

        CErr("Invalid TPInts type for Wavefunction<dcomplex,dcomplex>",std::cout);

      }

      p->TPI->printContractionTiming = scfControls.printContractionTiming;

    } else {

      CErr("Complex INT + Real WFN is not a valid option",std::cout);

    }



    ss->scfControls = scfControls;



    return ss;

  }; // SingleSlaterOptions::buildSingleSlater







  /**
   *  Outputs relevant information for the HamiltonianOptions struct
   *  to a specified output.
   *
   *  \param [in/out] out     Ouput device
   *  \param [in]     options HamiltonianOptions object to output.
   */
  std::ostream& operator<<(std::ostream &out, const HamiltonianOptions &options) {

    out << std::endl << "Hamiltonian Options";
    out << ":" << std::endl << BannerTop << std::endl << std::endl;


    const int fieldNameWidth(40);

    out << "  " << std::setw(fieldNameWidth) << "Integral:" << std::endl;
    out << bannerMid << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Basis Type:";
    switch (options.basisType) {
    case REAL_GTO:
      out << "REAL_GTO";
      break;
    case COMPLEX_GIAO:
      out << "COMPLEX_GIAO";
      break;
    case COMPLEX_GTO:
      out << "COMPLEX_GTO";
      break;
    }
    out << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Finite Width Nuclei:"
        << (options.finiteWidthNuc ? "True" : "False") << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Using Libcint:"
        << (options.Libcint ? "True" : "False") << std::endl;
    out << std::endl;


    out << "  " << std::setw(fieldNameWidth) << "One-Component Options:" << std::endl;
    out << bannerMid << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Perturbative Scalar Relativity:"
        << (options.PerturbativeScalarRelativity ? "On" : "Off") << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Perturbative Spin-orbit Relativity:"
        << (options.PerturbativeSpinOrbit ? "On" : "Off") << std::endl;
    out << std::endl;


    out << "  " << std::setw(fieldNameWidth) << "Two-Component Options:" << std::endl;
    out << bannerMid << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "One-Electron Scalar Relativity:"
        << (options.OneEScalarRelativity ? "On" : "Off") << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "One-Electron Spin-orbit Relativity:"
        << (options.OneESpinOrbit ? "On" : "Off") << std::endl;

    std::function<std::string(SNSO_TYPE)> snsoTypeToString = [](SNSO_TYPE type) {
      switch (type) {
      case SNSO_TYPE::BOETTGER:
        return "Boettger";
      case SNSO_TYPE::DC:
        return "Dirac-Coulomb";
      case SNSO_TYPE::DCB:
        return "Dirac-Coulomb-Breit";
      case SNSO_TYPE::ROW_DEP_DCB:
        return "Row-dependent Dirac-Coulomb-Breit";
      }
      return "Unknown";
    };

    out << "  " << std::setw(fieldNameWidth) << "Screened Nuclear Spin-orbit Approximation:"
        << (options.SNSO ? snsoTypeToString(options.snsoType) : "Off") << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Atomic Mean Field Spin-orbit:"
        << (options.AtomicMeanField ? "On" : "Off") << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Atomic X2C:"
        << (options.AtomicX2C ? "On" : "Off") << std::endl;
    if (options.AtomicX2C)
      out << "  " << std::setw(fieldNameWidth) << "Atomic X2C Type:"
          << options.AtomicX2CType.toString() << std::endl;
    out << std::endl;


    out << "  " << std::setw(fieldNameWidth) << "Four-Component Options:" << std::endl;
    out << bannerMid << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Bare Coulomb (LLLL) Term: "
        << (options.BareCoulomb ? "On" : "Off") << std::endl;

    char TYPE_4C_NAME[3][20] = { "All", "Spin Free Only", "Spin Dependent Only" };
    char TYPE_4C_APPROXIMATION[5][20] = {"None", "Three Center", "Two Center", "One Center", "Atomic Mean Field" };

    out << "  " << std::setw(fieldNameWidth) << "Dirac Coulomb (w/o SSSS) Term: "
        << (options.DiracCoulomb ? "On" : "Off") << std::endl;
    if(options.DiracCoulomb)
    out << "  " << std::setw(fieldNameWidth) << "Contribution---"
        << TYPE_4C_NAME[static_cast<int>(options.DiracCoulombType)] << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Approximation---"
        << TYPE_4C_APPROXIMATION[static_cast<int>(options.DiracCoulombApproximationType)] << std::endl;

    out << "  " << std::setw(fieldNameWidth) << "SSSS Term: "
        << (options.DiracCoulombSSSS ? "On" : "Off") << std::endl;
    if(options.DiracCoulombSSSS)
    out << "  " << std::setw(fieldNameWidth) << "Contribution---"
        << TYPE_4C_NAME[static_cast<int>(options.SSSSType)] << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Approximation---"
        << TYPE_4C_APPROXIMATION[static_cast<int>(options.SSSSApproximationType)] << std::endl;

    out << "  " << std::setw(fieldNameWidth) << "Gaunt Term: "
        << (options.Gaunt ? "On" : "Off") << std::endl;
    if(options.Gaunt)
    out << "  " << std::setw(fieldNameWidth) << "Contribution---"
        << TYPE_4C_NAME[static_cast<int>(options.GauntType)] << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Approximation---"
        << TYPE_4C_APPROXIMATION[static_cast<int>(options.GauntApproximationType)] << std::endl;

    out << "  " << std::setw(fieldNameWidth) << "Gauge Term: "
        << (options.Gauge ? "On" : "Off") << std::endl;
    if(options.Gauge)
    out << "  " << std::setw(fieldNameWidth) << "Contribution---"
        << TYPE_4C_NAME[static_cast<int>(options.GaugeType)] << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Approximation---"
        << TYPE_4C_APPROXIMATION[static_cast<int>(options.GaugeApproximationType)] << std::endl;
    out << std::endl;


    char TYPE_DKS_NAME[3][20] = {"Off","VLL","FULL"};

    out << "  " << std::setw(fieldNameWidth) << "Dirac-Kohn-Sham Options:" << std::endl;
    out << bannerMid << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Approximation---"
	<< TYPE_DKS_NAME[static_cast<int>(options.dksType)] << std::endl;
    if(options.Gaunt)
    out << "  " << std::setw(fieldNameWidth) << "Gaunt Scaling---"
        << options.GauntScale << std::endl;
    if(options.Gauge)
    out << "  " << std::setw(fieldNameWidth) << "Gauge Scaling---"
        << options.GaugeScale << std::endl;

    out << std::endl << BannerEnd << std::endl;

    return out; // Return std::ostream reference

  }


  // Regular SingleSlater wrapper
  SingleSlaterOptions CQSingleSlaterOptions(
    std::ostream &out, CQInputFile &input,
    Molecule &mol, BasisSet &basis) {

    return getSingleSlaterOptions(out, input, mol, basis, {-1., 1.}, "QM");

  }


  // NEO SingleSlater wrapper
  std::tuple<std::shared_ptr<SingleSlaterBase>, SingleSlaterOptions, SingleSlaterOptions> CQNEOSSOptions(
    std::ostream &out, CQInputFile &input,
    Molecule &mol,
    BasisSet &ebasis, BasisSet &pbasis,
    std::shared_ptr<IntegralsBase> eaoints, 
    std::shared_ptr<IntegralsBase> paoints,
    std::shared_ptr<IntegralsBase> epaoints,
    SCFControls scfControls) {

    Particle p{-1., 1.};
#define NEO_LIST(T) \
    MPI_COMM_WORLD,mol,ebasis,std::dynamic_pointer_cast<Integrals<T>>(epaoints),1,false,p

    SingleSlaterOptions essopt = getSingleSlaterOptions(out, input, mol, ebasis, {-1., 1.}, "QM");

    // If doing a NEO calculation with different Deuterium, scale the mass here
    // Note this also checks that all quantum particles are the same mass, as right now
    // this is the only options with how NEOSS is constructed
    double massAMU, mass;
    size_t nelec;
    if(mol.atomsQ.size())
    {
      massAMU = mol.atoms[mol.atomsQ[0]].atomicMass;
      nelec = mol.atoms[mol.atomsQ[0]].atomicNumber;
      for(const auto & atomQIndex : mol.atomsQ)
      {
        if(mol.atoms[atomQIndex].atomicMass != massAMU)
          CErr("All particles for a NEO calculation must have the same mass!");
      }
    }
    // HardCoded masses for H/D/T based on NIST standards:
    if(massAMU == atomicReference["H-1"].atomicMass)
    {
      mass = ProtMassPerE(); // https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
    }
    else if(massAMU == atomicReference["H-2"].atomicMass)
    {
      mass = DeutMassPerE(); // https://physics.nist.gov/cgi-bin/cuu/Value?mdsme
    }
    else if(massAMU == atomicReference["H-3"].atomicMass)
    {
      mass = TritMassPerE(); // https://physics.nist.gov/cgi-bin/cuu/Value?mtsme
    }
    else
    {
      // Atomic masses are mass of nuclei + mass of associated electrons, so we need
      // to subtract out the electron mass to just get the nuclear mass
      mass = massAMU * AUPerAMU() - nelec;
    }
   
    SingleSlaterOptions pssopt = getSingleSlaterOptions(out, input, mol, pbasis, {1., mass}, "PROTQM");

    std::shared_ptr<SingleSlaterBase> ess = essopt.buildSingleSlater(out,  mol, ebasis, eaoints);
    std::shared_ptr<SingleSlaterBase> pss = pssopt.buildSingleSlater(out,  mol, pbasis, paoints);

    std::shared_ptr<SingleSlaterBase> neoss;

    if(auto ess_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ess)) {
      if(auto pss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(pss)) {
        
        if(scfControls.printContractionTiming){
          ess_t->TPI->printContractionTiming = true;
          pss_t->TPI->printContractionTiming = true;
        }
        
        auto neoss_t = std::make_shared<NEOSS<double,double>>(NEO_LIST(double));
        auto epaoints_t = std::dynamic_pointer_cast<Integrals<double>>(epaoints);
        neoss_t->addSubsystem("Electronic", ess_t, {});
        neoss_t->addSubsystem("Protonic", pss_t, {{"Electronic", {true, epaoints_t->TPI}}});
        neoss_t->setOrder({"Protonic", "Electronic"});
        neoss = std::dynamic_pointer_cast<SingleSlaterBase>(neoss_t);

        // Handle the fact that VXC will be formed by the NEOKohnShamBuilder
        if( pssopt.refOptions.isEPCRef ) {
          auto pks_t = std::dynamic_pointer_cast<KohnSham<double,double>>(pss_t);
          pks_t->doVXC_ = false;
          if( auto eks_t = std::dynamic_pointer_cast<KohnSham<double,double>>(ess_t) ) {
            eks_t->doVXC_ = false;
          }
        }
      }
      else
        CErr("Electrons and protons must use the same field (real/real) or (complex/complex)");
    }
    else if(auto ess_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ess)) {
      if(auto pss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(pss)) {
        
        if(scfControls.printContractionTiming){
          ess_t->TPI->printContractionTiming = true;
          pss_t->TPI->printContractionTiming = true;
        }
        
        auto neoss_t = std::make_shared<NEOSS<dcomplex,double>>(NEO_LIST(double));
        auto epaoints_t = std::dynamic_pointer_cast<Integrals<double>>(epaoints);
        neoss_t->addSubsystem("Electronic", ess_t, {});
        neoss_t->addSubsystem("Protonic", pss_t, {{"Electronic", {true, epaoints_t->TPI}}});
        neoss_t->setOrder({"Protonic", "Electronic"});
        neoss = std::dynamic_pointer_cast<SingleSlaterBase>(neoss_t);

        // Handle the fact that VXC will be formed by the NEOKohnShamBuilder
        if( pssopt.refOptions.isEPCRef ) {
          auto pks_t = std::dynamic_pointer_cast<KohnSham<dcomplex,double>>(pss_t);
          pks_t->doVXC_ = false;
          if( auto eks_t = std::dynamic_pointer_cast<KohnSham<dcomplex,double>>(ess_t) ) {
            eks_t->doVXC_ = false;
          }
        }
      }
      else
        CErr("Electrons and protons must use the same field (real/real) or (complex/complex)");
    }
    else if(auto ess_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ess)) {
      if(auto pss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(pss)) {

        if(scfControls.printContractionTiming){
          ess_t->TPI->printContractionTiming = true;
          pss_t->TPI->printContractionTiming = true;
        }

        auto neoss_t = std::make_shared<NEOSS<dcomplex,dcomplex>>(NEO_LIST(dcomplex));
        auto epaoints_t = std::dynamic_pointer_cast<Integrals<dcomplex>>(epaoints);
        neoss_t->addSubsystem("Electronic", ess_t, {});
        neoss_t->addSubsystem("Protonic", pss_t, {{"Electronic", {true, epaoints_t->TPI}}});
        neoss_t->setOrder({"Protonic", "Electronic"});
        neoss = std::dynamic_pointer_cast<SingleSlaterBase>(neoss_t);

        // Handle the fact that VXC will be formed by the NEOKohnShamBuilder
        if( pssopt.refOptions.isEPCRef ) {
          auto pks_t = std::dynamic_pointer_cast<KohnSham<dcomplex,dcomplex>>(pss_t);
          pks_t->doVXC_ = false;
          if( auto eks_t = std::dynamic_pointer_cast<KohnSham<dcomplex,dcomplex>>(ess_t) ) {
            eks_t->doVXC_ = false;
          }
        }
      }
      else
        CErr("Electrons and protons must use the same field (real/real) or (complex/complex)");
    }
    else {
      CErr("NEO w/ mixed GTO/GIAO NYI, \nor Unreconized MatsT/IntsT combination in CQNEOSSOptions");
    }

    // Need to copy this over 
    epaoints->options_.erfOmega=pssopt.hamiltonianOptions.erfOmega;
    
    return {neoss, essopt, pssopt};

  }

  std::shared_ptr<SingleSlaterBase> CQMultiParticleSSOptions(
    std::ostream& out,
    CQInputFile& input,
    Molecule& mol,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    SCFControls scfControls) {
  
    if(quantumSubsystems.empty())              CErr("Cannot construct MultiParticleSS without quantum subsystems.");
    if(quantumSubsystems.front().label != "E") CErr("Electronic subsystem must be the first subsystem.");
  
    std::vector<std::shared_ptr<SingleSlaterBase>> subSS;
    subSS.reserve(quantumSubsystems.size());
  
    // Parse and retain every subsystem's own reference/hamiltonian options
    // before resolving guesses. Guess resolution needs the complete label set
    // and must not be overwritten by a later options parse.
    for(auto& sys : quantumSubsystems) {
      if(!sys.basis)     CErr("Missing basis for quantum subsystem " + sys.label);
      if(!sys.integrals) CErr("Missing integrals for quantum subsystem " + sys.label);
  
      sys.ssOptions = getSingleSlaterOptions(out, input, mol, *sys.basis, sys.particle, sys.qmSection);
      sys.ssOptions.hamiltonianOptions.savFilePrefix = "MULTISS/" + sys.label + "/";
      if(sys.particle.charge >= 0) {
        if(sys.ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::OFF)  CErr("X2C for other particles not implemented yet");
      }
      if(sys.nQuantumParticles > 0) sys.ssOptions.hamiltonianOptions.nParticleOverride = sys.nQuantumParticles;
    }

    resolveSubsystemGuessOptions(quantumSubsystems, scfControls);

    for(auto& sys : quantumSubsystems) {
      subSS.push_back(sys.ssOptions.buildSingleSlater(out,  mol, *sys.basis, sys.integrals));
    }
  
    if(std::dynamic_pointer_cast<SingleSlater<double,double>>(subSS.front())) {
      return buildMultiParticleSS<double,double>(out, mol, quantumSubsystems, quantumPairInteractions, subSS, scfControls);
    }
    else if(std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(subSS.front())) {
      return buildMultiParticleSS<dcomplex,double>(out, mol, quantumSubsystems, quantumPairInteractions, subSS, scfControls);
    }
    else if(std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(subSS.front())) {
      return buildMultiParticleSS<dcomplex,dcomplex>(out, mol, quantumSubsystems, quantumPairInteractions, subSS, scfControls);
    }
    else {
      CErr("Unsupported MatsT/IntsT combination in CQMultiParticleSSOptions");
    }
    return nullptr;
  }


  template <typename MatsT, typename IntsT>
  std::shared_ptr<SingleSlaterBase> buildMultiParticleSS(
    std::ostream& out,
    Molecule& mol,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    const std::vector<std::shared_ptr<SingleSlaterBase>>& subSS,
    SCFControls scfControls) {
  
    if(quantumSubsystems.size() != subSS.size()) CErr("Quantum subsystem descriptor/object count mismatch.");
  
    const auto& firstSys = quantumSubsystems.front();
    auto firstInts = std::dynamic_pointer_cast<Integrals<IntsT>>(firstSys.integrals);
  
    auto multiSS = std::make_shared<MultiParticleSS<MatsT,IntsT>>(MPI_COMM_WORLD, mol, *firstSys.basis, firstInts, 1, false, firstSys.particle);
    multiSS->scfControls = scfControls;
  
    for(size_t i = 0; i < quantumSubsystems.size(); ++i) {
      const auto& sys = quantumSubsystems[i];
      auto ss = std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(subSS[i]);
      if(!ss) CErr("All quantum subsystems must use the same MatsT/IntsT combination. Mismatch for " + sys.label);
      if(scfControls.printContractionTiming && ss->TPI) ss->TPI->printContractionTiming = true;
      if(!sys.atomIndices.empty()) ss->ownedAtomIndices = sys.atomIndices; 
      
      multiSS->addSubsystem(sys.label, ss);
      multiSS->setSubsystemGuessOptions(sys.label, sys.ssOptions);
    }
  
    for(const auto& interaction : quantumPairInteractions) {
      auto pairInts = std::dynamic_pointer_cast<Integrals<IntsT>>(interaction.integrals);
      if(!pairInts) CErr("Invalid pair integral scalar type for " + interaction.labelA + "-" + interaction.labelB);
      if(!pairInts->TPI) CErr("Missing two-particle integrals for " + interaction.labelA + "-" + interaction.labelB);
      
      multiSS->addInteraction(interaction.labelA, interaction.labelB, pairInts->TPI, false);
    }

    // Set up cross-particle correlation 
    // (Ideally this should be general for any pair of subsystems, but currently only supports EPC)
    bool haveEPC = false;
    for(size_t i = 0; i < quantumSubsystems.size(); ++i) {
      const auto& sys = quantumSubsystems[i];
      if(!sys.ssOptions.refOptions.isEPCRef) continue;

      auto ks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(subSS[i]);
      if(!ks) CErr("EPC reference subsystem is not a KohnSham object: " + sys.label);

      multiSS->setInterFunctionals("E", sys.label, ks->functionals);
      // MultiParticleSS owns inter-XC; locally this subsystem remains HF with
      //   full exchange instead of inheriting the electronic GauXC xHFX.
      ks->functionals.clear();
      ks->doVXC_ = false;
      haveEPC = true;
      // This is hacky but makes sure non-electron particles does not use electron xHFX
      // in KohnSham::formFock
      // TODO: This needs to be fixed but probably Alternatives either breaks GauXC CUSTOM functionals 
      // or require adding per-subsystem exchange state/API. 
      ks->intParam.useGauXC = false; 
    }

    // With EPC present, the electron's intra-XC is formed by the unified driver instead its own formVXC
    if( haveEPC )
      if( auto eks = std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(subSS.front()) )
        eks->doVXC_ = false;

    // Use the electronic ("E") subsystem's numerical grid for the shared EPC pass
    multiSS->setIntParam(quantumSubsystems.front().ssOptions.intParam);

    multiSS->setSubSetup();
    multiSS->printSetup(out);
  
    return std::dynamic_pointer_cast<SingleSlaterBase>(multiSS);
  }


}; // namespace ChronusQ
