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
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>

#include <singleslater/neoss.hpp>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQQM_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "REFERENCE",
      "NUCREFERENCE",
      "JOB",
      "X2CTYPE",
      "SPINORBITSCALING",
      "ATOMICX2C",
      "SNSOTYPE"
    };

    // Specified keywords
    std::vector<std::string> qmKeywords = input.getDataInSection("QM");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : qmKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword QM." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQDFTINT_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "EPS",
      "NANG",
      "NRAD",
      "NMACRO",
      "INHOUSE",
      "GAUXC"
    };

    // Specified keywords
    std::vector<std::string> dftintKeywords = input.getDataInSection("DFTINT");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : dftintKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword DFTINT." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
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
      "BHANDH"
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
    }
    RORefs.emplace_back( "ROHF" );
    FourCRefs.emplace_back( "4CHF" );

    // This is the reference string to be parsed
    std::string refString = tokens.back();
   
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

      OPTOPT( intParam.epsilon      = input.getData<double>("DFTINT.EPS")  );
      OPTOPT( intParam.nAng         = input.getData<size_t>("DFTINT.NANG") );
      OPTOPT( intParam.nRad         = input.getData<size_t>("DFTINT.NRAD") );
      OPTOPT( intParam.nRadPerBatch = input.getData<size_t>("DFTINT.NMACRO") );
      bool gauFlag1(false), gauFlag2(false);
      OPTOPT( gauFlag1              = not input.getData<bool>("DFTINT.INHOUSE") );
      OPTOPT( gauFlag2              = input.getData<bool>("DFTINT.GAUXC") ); 
      intParam.useGauXC = gauFlag1 or gauFlag2;

    }


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
    BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
    RefOptions &refOptions, HamiltonianOptions &hamiltonianOptions, std::string section) {

    // Parse hamiltonianOptions
    hamiltonianOptions.basisType = basis.basisType;

    std::string X;

    // Parse X2C option
    // X2CType = off (default), spinfree, onee, twoe
    X = "DEFAULT";
    OPTOPT( X = input.getData<std::string>(section + ".X2CTYPE")  );
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

      CErr(X + " not a valid " + section + ".X2CTYPE",out);

    }



    // Parse one-electron spin-orbie scaling option
    // SpinOrbitScaling  = noscaling, boettger (dafault), atomicmeanfield (amfi)
    X = "DEFAULT"; // Unspecified value
    OPTOPT( X = input.getData<std::string>(section + ".SPINORBITSCALING")  );
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

      CErr(X + " not a valid " + section + ".SPINORBITSCALING",out);

    }

    // Parse screened nuclear spin–orbit approximation
    X = "BOETTGER"; // Unspecified value
    OPTOPT( X = input.getData<std::string>(section + ".SNSOTYPE")  );
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

      CErr(X + " not a valid " + section + ".SNSOTYPE",out);

    }

    // Parse Atomic X2C
    hamiltonianOptions.AtomicX2C = parseAtomicType(out,input,hamiltonianOptions.AtomicX2CType,section);


    // Parse Finite Width Nuclei
    std::string finiteCore = "DEFAULT";
    OPTOPT( finiteCore = input.getData<std::string>("INTS.FINITENUCLEI"); )
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
    OPTOPT( hamiltonianOptions.Libcint = input.getData<bool>("INTS.LIBCINT") )

    if (hamiltonianOptions.Libcint) {
      if (basis.forceCart)
        CErr("Libcint + cartesian GTO NYI.");
      if (auto aoi = std::dynamic_pointer_cast<Integrals<double>>(aoints))
        if (auto rieri = std::dynamic_pointer_cast<InCoreAuxBasisRIERI<double>>(aoi->TPI))
          if (rieri->auxbasisSet()->forceCart)
            CErr("Libcint + cartesian GTO NYI.");
    }


    OPTOPT( hamiltonianOptions.BareCoulomb = input.getData<bool>("INTS.BARECOULOMB") )
    OPTOPT( hamiltonianOptions.BareCoulomb = input.getData<bool>("INTS.LLLL") )
    //OPTOPT( hamiltonianOptions.DiracCoulombSSSS = input.getData<bool>("INTS.SSSS") )
    //OPTOPT( hamiltonianOptions.DiracCoulomb = input.getData<bool>("INTS.DIRACCOULOMB") )
    //OPTOPT( hamiltonianOptions.Gauge = input.getData<bool>("INTS.GAUGE") )
    //OPTOPT( hamiltonianOptions.Gaunt = input.getData<bool>("INTS.GAUNT") )
    //try{
    //  if ( input.getData<bool>("INTS.BREIT") ) {
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
        DCOptions = input.getData<std::string>("INTS.DIRACCOULOMB");
      } catch (...) {
        DCOptions = input.getData<std::string>("INTS.DC");
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
      SSSSOptions = input.getData<std::string>("INTS.SSSS");
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

    // Gaunt
    try { 
      std::string GauntOptions = "FALSE";
      GauntOptions = input.getData<std::string>("INTS.GAUNT");
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      if( std::regex_search(GauntOptions, regexTRUE) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
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


    // Gauge
    try { 
      std::string GaugeOptions = "FALSE";
      GaugeOptions = input.getData<std::string>("INTS.GAUGE");
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      if( std::regex_search(GaugeOptions, regexTRUE) ) {
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
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
      BreitOptions = input.getData<std::string>("INTS.BREIT");
      auto const regexTRUE = std::regex("true|on",std::regex_constants::icase);
      auto const regexALL = std::regex("all|exact",std::regex_constants::icase);
      auto const regexFALSE = std::regex("off|none|false",std::regex_constants::icase);
      auto const regexSF = std::regex("sf|spinfree",std::regex_constants::icase);
      auto const regexSD = std::regex("sd|spindependent",std::regex_constants::icase);
      auto const regex3C = std::regex("3c|3center|threecenter",std::regex_constants::icase);
      auto const regex2C = std::regex("2c|2center|twocenter",std::regex_constants::icase);
      auto const regex1C = std::regex("1c|1center|onecenter",std::regex_constants::icase);
      auto const regexAMF = std::regex("amf|atomicmeanfield",std::regex_constants::icase);

      if( std::regex_search(BreitOptions, regexTRUE) ) {
        hamiltonianOptions.Gaunt = true;
        hamiltonianOptions.GauntType = TYPE_4C::All;
        hamiltonianOptions.GauntApproximationType = APPROXIMATION_TYPE_4C::None;
        hamiltonianOptions.Gauge = true;
        hamiltonianOptions.GaugeType = TYPE_4C::All;
        hamiltonianOptions.GaugeApproximationType = APPROXIMATION_TYPE_4C::None;
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
    OPTOPT( X = input.getData<std::string>(section + ".ATOMICX2C")  );
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
      CErr(X + " not a valid " + section + ".ATOMICX2C",out);
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
      std::shared_ptr<IntegralsBase> aoints,
      Particle p, std::string section) {

    out << "  *** Parsing " << section << ".REFERENCE options ***\n";

    SingleSlaterOptions options;

    // Attempt to find reference
    std::string reference;
    try { 
      reference = input.getData<std::string>(section + ".REFERENCE");
    } catch(...) {
      CErr(section + ".REFERENCE Keyword not found!",out);
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
      
      // Determine EPC functional is available in CQ/GauXC, and edit keywords in needed
      if(options.refOptions.isEPCRef){
        if(not options.intParam.useGauXC){
          if(!options.refOptions.funcName.compare("EPC17_1"))
            CErr("EPC17_1 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("EPC18_1"))
            CErr("EPC18_1 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
          else if(!options.refOptions.funcName.compare("EPC18_2"))
            CErr("EPC18_2 Not Implemented In-House. Turn inhouse flag to false to use GauXC");
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
    parseHamiltonianOptions(out,input,basis,aoints,
        options.refOptions,options.hamiltonianOptions,section);

    options.hamiltonianOptions.particle = p;

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
      std::shared_ptr<IntegralsBase> aoints) const {


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

    if( p.charge < 0 && refOptions.isEPCRef )
      CErr("EPC functionals only valid on proton references!");

    if( p.charge >= 0 && refOptions.isKSRef && !refOptions.isEPCRef )
      CErr("Proton Kohn Sham references require EPC functionals");

  #define KS_LIST(T) \
    refOptions.funcName,funcList,MPI_COMM_WORLD,intParam,mol,basis,std::dynamic_pointer_cast<Integrals<T>>(aoints),refOptions.nC,refOptions.iCS,p

  #define HF_LIST(T) \
    MPI_COMM_WORLD,mol,basis,std::dynamic_pointer_cast<Integrals<T>>(aoints),refOptions.nC,refOptions.iCS,p

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

      size_t nERI4DCB = 0; // Bare-Coulomb

      if( hamiltonianOptions.Gaunt ) nERI4DCB = 23; // Dirac-Coulomb-Gaunt
      else if( hamiltonianOptions.DiracCoulomb ) nERI4DCB = 4; // Dirac-Coulomb

      
      if( hamiltonianOptions.DiracCoulombSSSS ) nERI4DCB += 16; // Dirac-Coulomb-SSSS

      if( hamiltonianOptions.Gauge ) nERI4DCB += 26; // Gauge


      if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

        std::shared_ptr<TwoPInts<double>> &TPI =
            std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

        if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<double>>(TPI)) {

          TPI = std::make_shared<InCore4indexRelERI<double>>(basis.nBasis,nERI4DCB);

          p->TPI = std::make_shared<InCore4indexRelERIContraction<double,double>>(TPI);

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

          TPI = std::make_shared<InCore4indexRelERI<double>>(basis.nBasis,nERI4DCB);

          p->TPI = std::make_shared<InCore4indexRelERIContraction<dcomplex,double>>(TPI);

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

      if (auto tpi_typed = std::dynamic_pointer_cast<DirectTPI<double>>(TPI)) {

        p->TPI = std::make_shared<GTODirectTPIContraction<double,double>>(tpi_typed);

      } else if (auto tpi_typed = std::dynamic_pointer_cast<InCoreRITPI<double>>(TPI)) {

        p->TPI = std::make_shared<InCoreRITPIContraction<double,double>>(tpi_typed);

      } else if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<double>>(TPI)) {

        p->TPI = std::make_shared<InCore4indexTPIContraction<double,double>>(tpi_typed);

      } else {

        CErr("Invalid TPInts type for Wavefunction<double,double>",std::cout);

      }
      
      p->TPI->printContractionTiming = scfControls.printContractionTiming;
    
    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

      std::shared_ptr<TwoPInts<double>> TPI =
          std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI;

      if (auto tpi_typed = std::dynamic_pointer_cast<DirectTPI<double>>(TPI)) {

        p->TPI = std::make_shared<GTODirectTPIContraction<dcomplex,double>>(tpi_typed);

      } else if (auto tpi_typed = std::dynamic_pointer_cast<InCoreRITPI<double>>(TPI)) {

        p->TPI = std::make_shared<InCoreRITPIContraction<dcomplex,double>>(tpi_typed);

      } else if (auto tpi_typed = std::dynamic_pointer_cast<InCore4indexTPI<double>>(TPI)) {

        p->TPI = std::make_shared<InCore4indexTPIContraction<dcomplex,double>>(tpi_typed);

      } else {

        CErr("Invalid TPInts type for Wavefunction<dcomplex,double>",std::cout);

      }
      
      p->TPI->printContractionTiming = scfControls.printContractionTiming;
    
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


    out << std::endl << BannerEnd << std::endl;

    return out; // Return std::ostream reference

  }


  // Regular SingleSlater wrapper
  SingleSlaterOptions CQSingleSlaterOptions(
    std::ostream &out, CQInputFile &input,
    Molecule &mol, BasisSet &basis,
    std::shared_ptr<IntegralsBase> aoints) {

    return getSingleSlaterOptions(out, input, mol, basis, aoints, {-1., 1.}, "QM");

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

    SingleSlaterOptions essopt = getSingleSlaterOptions(out, input, mol, ebasis, eaoints, {-1., 1.}, "QM");
    SingleSlaterOptions pssopt = getSingleSlaterOptions(out, input, mol, pbasis, paoints, {1., ProtMassPerE}, "PROTQM");

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
    else {
      CErr("NEO + GIAO NYI!");
    }

    return {neoss, essopt, pssopt};

  }


}; // namespace ChronusQ
