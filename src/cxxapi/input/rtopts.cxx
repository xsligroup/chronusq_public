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
#include<regex>
#include <cerr.hpp>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQRT_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "TYPE",          // Type of dynamics: BOMD (Default), Ehrenfest, RT
      "TMAX",          // The total time for the whole dynamics: 100 fs (Default)
      "MAXSTEPS",
      "UNITS",         // The units of time: FS (Default), AU //TODO: this option is never used
      "DELTAT",
      "IRSTRT",
      "FIELD",
      "FIELDINDEPENDENTHAMILTONIAN",
      "INTALG",
      "PROT_INTALG",
      "RESTARTALG",
      "RESTARTSTEP",
      "RESTARTFROM",
      "SAVESTEP",
      "SAVEONEPDM",
      "RESTART",
      "SCFFIELD",
      "PRINTLEVEL",
      "CIPOPULATION",
      "INITTYPE",
      "LCWEIGHTS",
      "LCSTATES",
      "COEFFS",
      "DETS",
      "REALTIMECORRELATIONFUNC",
      "REALTIMECORRELATIONFUNCSTART",
      "PRINTDEN",
      "PRINTCONTRACTIONTIMING",
      "PRINTSTEP",
      "RTGAUNT",
      "RTPRINTDEN",
      "RTGAUGE",
      "RTBREIT",
      "ORBITALPOPFREQ"
    };

    // Specified keywords
    std::vector<std::string> rtKeywords = input.getDataInSection("RT");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : rtKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword RT." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  \brief Construct a RealTime object using the input 
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

  std::shared_ptr<RealTimeBase> CQRealTimeOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> &ss, std::shared_ptr<MCWaveFunctionBase> &mcscf,
    std::shared_ptr<TDEMPerturbation>& tdPert,
    EMPerturbation& scfPert ) {

    std::shared_ptr<RealTimeBase> rt;
    if(mcscf) {
        rt = CQRealTimeMultiSlaterOptions(out, input, mcscf, scfPert);
    }
    std::cout << "returning RTMULTISLATER OPTS FROM RT" << std::endl;
    return rt;
  }

  std::shared_ptr<RealTimeBase> CQRealTimeMultiSlaterOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<MCWaveFunctionBase> &mcscf,
    EMPerturbation& scfPert ) {
    if( not input.containsSection("RT") )
      CErr("RT Section must be specified for RT job",out);


    out << "  *** Parsing MRRT options ***\n";

    std::shared_ptr<RealTimeMultiSlaterBase> rt;

    // Determine  reference and construct RT object

    auto inputRTAlg = RealTimeAlgorithm::RTRungeKuttaOrderFour;
    std::shared_ptr<RealTimeMultiSlaterVectorManagerBase> vecManager;
    try {
      auto intAlg = input.getData<std::string>("RT.INTALG");

      if ( not intAlg.compare("SSO") ) { 
        inputRTAlg = RealTimeAlgorithm::RTSymplecticSplitOperator;
        vecManager = std::make_shared<RealTimeMultiSlaterVectorManagerSSO<double*>>();
      }
      else if ( not intAlg.compare("RK4") ) {
      }
      else {
          std::cout << "Could not understand RT.INTALG. Defaulting to RK4.";
          std::cout << std::endl;
      }
    }
    catch(...) {
      std::cout << "Defaulting to RK4 integration algorithm" << std::endl;
    };
    
    if (!vecManager) {
        vecManager = std::make_shared<RealTimeMultiSlaterVectorManagerRK4<dcomplex*>>();
    }
    vecManager->set_vecSize_(mcscf->NDet);

    // Determines the initial CI vector to be propogated 
    HandleRTInitState(out,input,vecManager);

    bool found = false;

    #define CONSTRUCT_RT_MR(_REF,_MT,_IT,_RTALG,_MRWFN)             \
    if( not found ) try {                          \
      if (std::dynamic_pointer_cast<_REF<_MT, _IT> >(_MRWFN)){ \
      rt =     std::make_shared< RealTimeMultiSlater<_MT, _IT> >(  \
              std::dynamic_pointer_cast<_REF<_MT, _IT> >(mcscf), vecManager, _RTALG); \
      found = true;                                \
      } \
    } catch(...) {  }

    // Construct RT object (MatsT, IntsT, PropT)
    if (inputRTAlg == RealTimeAlgorithm::RTSymplecticSplitOperator ) {
        CONSTRUCT_RT_MR( MCWaveFunction, double, double, RealTimeAlgorithm::RTSymplecticSplitOperator, mcscf);
    }
    //CONSTRUCT_RT_MR( MCWaveFunction, double, double, RealTimeAlgorithm::RTRungeKuttaOrderFour,mcscf);
    CONSTRUCT_RT_MR( MCWaveFunction, dcomplex, double, RealTimeAlgorithm::RTRungeKuttaOrderFour,mcscf);
    // no magnetic fields... yet.
    //CONSTRUCT_RT_MR( MCWaveFunction, dcomplex, dcomplex );

    // Parse Options
    try {
      rt->intScheme.tMax = input.getData<double>("RT.TMAX");
    } catch(...) {
      CErr("Must specify RT.TMAX for simulation length");
    }

    try {
      rt->intScheme.deltaT = input.getData<double>("RT.DELTAT");
    } catch(...) {
      CErr("Must specify RT.DELTAT for integration time step");
    }

    // Determine Integration Algorithm
    rt->intScheme.integrationAlgorithm = inputRTAlg;

    // Set SCF perturbation
    rt->setSCFPerturbation( scfPert );
    // Inclusion of SCF perturbation
    OPTOPT(
      rt->intScheme.includeSCFField = input.getData<bool>("RT.SCFFIELD");
    )

    // Save frequency
    OPTOPT(
      rt->intScheme.iSave = input.getData<size_t>("RT.SAVESTEP")
    )

    // Whether we are restarting an RT calculation
    OPTOPT(
      rt->restart = input.getData<bool>("RT.RESTART");
    )

    // Amount of printing in the RT calc
    OPTOPT(
      rt->printLevel = input.getData<size_t>("RT.PRINTLEVEL");
    )

    OPTOPT(
      rt->CIPopFreq = input.getData<size_t>("RT.CIPOPULATION");
    )

    OPTOPT(
      rt->RealTimeCorrelationFunctionFreq = input.getData<size_t>("RT.REALTIMECORRELATIONFUNC");
    )

    OPTOPT(
      rt->RealTimeCorrelationFunctionStart = input.getData<double>("RT.REALTIMECORRELATIONFUNCSTART");
    )

    OPTOPT(
      rt->time_independent_ham = input.getData<bool>("RT.FIELDINDEPENDENTHAMILTONIAN");
    )

    // Whether to print time-dependent density
    OPTOPT(
      rt->printCIVec = input.getData<bool>("RT.PRINTCIVEC");
    )

    return rt;
  }; // CQRealTimeMultiSlaterOptions
  
std::shared_ptr<TDEMFieldBase> parseRTField(std::string& fieldStr, std::ostream& out){
    // Split line on white space
    std::vector<std::string> tokens;
    split(tokens,fieldStr," \t");

    if( tokens.size() == 0 ) return std::shared_ptr<TDEMFieldBase>(nullptr);

    for(auto &X : tokens) trim(X);
    
    // Only Dipole fields for now
    if( tokens.size() != 5 )
      CErr("\"" + fieldStr + "\" not a vaild FIELD specification",out);

    // Determine field type
    std::string fieldTypeStr = tokens[1];

    EMFieldTyp fieldType;
    if( not fieldTypeStr.compare("ELECTRIC") )
      fieldType = Electric;
    else if( not fieldTypeStr.compare("MAGNETIC") )
      CErr("Magnetic Fields NYI");
    else
      CErr(fieldTypeStr + "not a valid Field type");

    // Only DIPOLE implemented
    cart_t DipoleField = {std::stod(tokens[2]), std::stod(tokens[3]), 
                          std::stod(tokens[4])};

    // Handle envelope specification
    std::string envelope = tokens[0];
    // This string will have things in the following form
    // envelopeType[parameters](tOn,tOff) fieldType x y z
    //   - STEPFIELD(tOn,tOff) ELECTRIC x y z
    //   - LINEARRAMP(tOn,tOff) ELECTRIC x y z
    //   - PLANEWAVE[omega](tOn,tOff) ELECTRIC x y z
    //   - GAUSSIAN[alpha](tOn,tOff) ELECTRIC x y z

    // Define all the containers, regex, and lambdas that we need
    std::string timeParameters;
    std::string envelopeParameters;
    std::string envelopeType;
    auto const timeParametersMatcher = std::regex("\\((.*?)\\)");
    auto const envelopeParametersMatcher = std::regex("\\[(.*?)\\]");
    auto const envelopeTypeMatcher = std::regex("^[^([]+");
    auto get_matched_str = [&] (auto const& matcher, std::string& out_str) {
        std::smatch matched_section;
        out_str = "";
        if (std::regex_search(envelope, matched_section, matcher))
            out_str = matched_section.str();
    };
    auto erase_char = [&] (std::string& in_str, const std::vector<char> remove_char_vec) {
        for (char remove_char : remove_char_vec)
            in_str.erase(std::remove(in_str.begin(), in_str.end(), remove_char), in_str.end());
        in_str.erase(std::remove(in_str.begin(), in_str.end(), ' '), in_str.end());
    };
    // Parse
    get_matched_str(envelopeTypeMatcher, envelopeType);
    get_matched_str(envelopeParametersMatcher, envelopeParameters);
    get_matched_str(timeParametersMatcher, timeParameters);
    erase_char(envelopeType, {});
    erase_char(envelopeParameters, {'[', ']'});
    erase_char(timeParameters, {'(',  ')'});

    // Determine if valid specifcation
    if (timeParameters.empty())
      CErr(envelope + " not a valid specification. Missing time on and off specifications.",out);
    
    std::vector<std::string> timeTokens;
    split(timeTokens, timeParameters, ",");
    if( timeTokens.size() != 2 )
      CErr("TDField needs at most 2 time arguments. For example, 'envelopeType[parameters](tOn,tOff) fieldType x y z'", out);
    double stepOn  = std::stod(timeTokens[0]);
    double stepOff = std::stod(timeTokens[1]);
    if( stepOff <= stepOn )
      CErr("STEPOFF must be > STEPON for TDField");

    std::vector<std::string> parameterTokens;
    split(parameterTokens, envelopeParameters, ",");
    if ( parameterTokens.size() > 2)
      CErr("Too many parameters provided. All currently implemented envelopes accept at most 2 parameters.", out);


    //   - PLANEWAVE[omega](tOn,tOff) ELECTRIC x y z
    //   - GAUSSIAN[alpha](tOn,tOff) ELECTRIC x y z
    // STEPFIELD
    if( envelopeType == "STEPFIELD" ) {
      if ( parameterTokens.size() != 0)
        CErr("Too many parameters provided. STEPFIELD envelope accepts no parameters.", out);
      // Append Field
      // XXX: Should store pointer to field base
      // and then append after envelope is determined
      // auto fieldenvelope = FieldEnvelope<FieldEnvelopeType::Step>(stepOn,stepOff);
    //auto a = TDEMField(fieldType,fieldenvelope, DipoleField);
    //auto a2 = std::make_shared<TDEMField<cart_t>>(fieldType,fieldenvelope, DipoleField);
      return 
      std::move(
        std::dynamic_pointer_cast<TDEMFieldBase, TDEMField<cart_t> >(
          std::make_shared<TDEMField<cart_t>>(fieldType,
          FieldEnvelope<FieldEnvelopeType::Step>(stepOn,stepOff)
          , DipoleField)
        )
      );
    } else if ( envelopeType == "LINEARRAMPFIELD" ) {
      if ( parameterTokens.size() != 0)
        CErr("Too many parameters provided. LINEARRAMP envelope accepts no parameters.", out);
      //return std::make_shared<TDEMField>(fieldType, LinRampField(stepOn,stepOff), DipoleField);
      return 
      std::move(
        std::dynamic_pointer_cast<TDEMFieldBase, TDEMField<cart_t>>(
          std::make_shared<TDEMField<cart_t>>(fieldType,
          FieldEnvelope<FieldEnvelopeType::LinRamp>(stepOn,stepOff)
          , DipoleField)
        )
      );
    } else if ( envelopeType == "PLANEWAVEFIELD" ) {
        bool doCos = true;
        if (parameterTokens.size() == 2) {
            if (parameterTokens[1] == "SIN" or parameterTokens[1] == "FALSE")
              doCos = false;
            else if (parameterTokens[1] == "COS" or parameterTokens[1] == "TRUE")
                doCos = true;
            else 
                CErr("PlaneWaveField needs at most 2 parameter arguments in the correct order. For example, 'PlaneWaveField[omega,COS(true) or SIN(false)](tOn,tOff) fieldType x y z'", out);
            parameterTokens.erase(parameterTokens.end() - 1);
        }
        std::vector<double> processedEnvelopeParameters(parameterTokens.size());
        std::transform(parameterTokens.begin(), parameterTokens.end(), processedEnvelopeParameters.begin(), [](const std::string& val){ return std::stod(val); });
        double omega = processedEnvelopeParameters[0];
      //return std::make_shared<TDEMField>(fieldType, PlaneWaveField(stepOn,stepOff,omega, doCos), DipoleField);
      return 
      std::move(
        std::dynamic_pointer_cast<TDEMFieldBase, TDEMField<cart_t>>(
          std::make_shared<TDEMField<cart_t>>(fieldType,
          FieldEnvelope<FieldEnvelopeType::PlaneWave>(stepOn,stepOff, omega, doCos)
          , DipoleField)
        )
      );
    } else if ( envelopeType == "GAUSSIANFIELD" ) {
       std::vector<double> processedEnvelopeParameters(parameterTokens.size());
       std::transform(parameterTokens.begin(), parameterTokens.end(), processedEnvelopeParameters.begin(), [](const std::string& val){ return std::stod(val); });
        double alpha = processedEnvelopeParameters[0];
        if ( alpha < 0.0 )
          CErr("For a GAUSSIAN envelope, alpha should be postive since the implementation includes a negative sign.", out);
      //return std::make_shared<TDEMField>(fieldType, GaussianField(stepOn,stepOff,alpha), DipoleField);
      return 
      std::move(
        std::dynamic_pointer_cast<TDEMFieldBase, TDEMField<cart_t>>(
          std::make_shared<TDEMField<cart_t>>(fieldType,
          FieldEnvelope<FieldEnvelopeType::Gaussian>(stepOn,stepOff, alpha)
          , DipoleField)
        )
      );

    } else CErr("Envelope not recognized or not implemented.");
  return std::shared_ptr<TDEMFieldBase>(nullptr);
  } // parseRTField

  void HandleRTInitState(std::ostream & out, CQInputFile & input, std::shared_ptr<RealTimeMultiSlaterVectorManagerBase> & vecManager)
  {
    std::string inittype;
    OPTOPT(inittype = input.getData<std::string>("RT.INITTYPE"))
    if(inittype.empty())
    {
      out << "** No detailed initial state detected, defaulting to propogating the lowest CIVec **" << std::endl;
      vecManager->initmethod = MSInitialState::LinearCombination;
      vecManager->init_detail.push_back(std::make_pair<double,size_t>(1.0,1));
      return;
    }
    else if(inittype=="LINEARCOMBINATION")
    {
      vecManager->initmethod = MSInitialState::LinearCombination;
      std::string weights,states;
      OPTOPT(weights = input.getData<std::string>("RT.LCWEIGHTS"));
      OPTOPT(states = input.getData<std::string>("RT.LCSTATES"));
      std::vector<std::string> weighttokens,statetokens;
      split(weighttokens,weights," ,;");
      split(statetokens,states," ,;");
      if(weighttokens.size()!=statetokens.size())
      {
        CErr("Number of RT Initial state weights != number of states!");
      }
      if(!weighttokens.size())
      {
        std::cout << "No Linear Combination weights specified, defaulting to lowest CI Vector" << std::endl;
        vecManager->init_detail.push_back(std::make_pair<double,size_t>(1.0,1));
      }
      std::vector<double> dweights;
      double sum;
      for(auto w : weighttokens)
        dweights.push_back(std::stod(w));
      sum = std::accumulate(dweights.begin(),dweights.end(),0.0);
      for(size_t i = 0; i < dweights.size(); i++)
        vecManager->init_detail.push_back(std::make_pair<double,size_t>(dweights[i]/std::sqrt(sum),std::stoi(statetokens[i])));

      out << "The initial state will be constructed of the following linear combination:" << std::endl;
      out << " Coeff, state:" << std::endl;
      for(auto lcstate : vecManager->init_detail)
      {
        out << lcstate.first << " , " << lcstate.second << std::endl;
      }
      out << std::endl;

    }
    else if(inittype=="CUSTOMCI")
    {
      vecManager->initmethod = MSInitialState::CustomCI;
      std::cout << "Creating a custom initial CI Vector" << std::endl;
      std::cout << "This is a very tenuous procedure!" << std::endl;
      std::cout << "You should use very cautiously and probably only" << std::endl
                << "if you REALLY know what you're doing!" << std::endl;
      std::string coeffs,dets;
      OPTOPT(coeffs = input.getData<std::string>("RT.COEFFS"));
      OPTOPT(dets = input.getData<std::string>("RT.DETS"));
      std::vector<std::string> coefftokens,dettokens;
      split(coefftokens,coeffs," ,;");
      split(dettokens,dets," ,;");
      if(coefftokens.size()!=dettokens.size())
      {
        CErr("Number of RT Initial weights != number of determinants!");
      }
      if(!coefftokens.size())
      {
        CErr("CustomCI requested but no parameters provided!");
      }
      std::vector<double> cweights;
      double sum;
      for(size_t i = 0; i < coefftokens.size(); i++)
        cweights.push_back(std::stof(coefftokens[i]));
      sum=std::accumulate(cweights.begin(),cweights.end(),0.0);
      for(size_t i = 0; i < cweights.size(); i++)
        vecManager->init_detail.push_back(std::make_pair<double,size_t>(cweights[i]/std::sqrt(sum),std::stoi(dettokens[i])));

    }
    else
    {
      CErr("Unrecognized option for RT.initmethod: " + inittype);
    }
    return;
  };
}; // namespace ChronusQ
