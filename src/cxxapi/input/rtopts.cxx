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
      "INTALG",
      "RESTARTALG",
      "RESTARTSTEP",
      "RESTARTFROM",
      "SAVESTEP",
      "RESTART",
      "SCFFIELD",
      "PRINTLEVEL",
      "ORBITALPOPFREQ",
      "PRINTDEN",
      "PRINTCONTRACTIONTIMING",
      "PRINTSTEP",
      "RTGAUNT",
      "RTPRINTDEN",
      "RTGAUGE",
      "RTBREIT"
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
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> &ss,
    EMPerturbation& scfPert ) {

    if( not input.containsSection("RT") )
      CErr("RT Section must be specified for RT job",out);


    out << "  *** Parsing RT options ***\n";

    std::shared_ptr<RealTimeBase> rt;

  
    // Determine  reference and construct RT object


    bool found = false;


    #define CONSTRUCT_RT(_REF,_MT,_IT)             \
    if( not found ) try {                          \
      rt = std::dynamic_pointer_cast<RealTimeBase>( \
          std::make_shared< RealTime<_REF,_IT> >(  \
            dynamic_cast< _REF<_MT,_IT>& >(*ss)    \
          )                                        \
        );                                         \
      found = true;                                \
    } catch(...) { }


    // Construct RT object
    CONSTRUCT_RT( NEOSS, double, double     );
    CONSTRUCT_RT( NEOSS, dcomplex, double   );
    CONSTRUCT_RT( NEOSS, dcomplex, dcomplex );

    CONSTRUCT_RT( HartreeFock, double, double     );
    CONSTRUCT_RT( HartreeFock, dcomplex, double   );
    CONSTRUCT_RT( HartreeFock, dcomplex, dcomplex );

    CONSTRUCT_RT( KohnSham, double, double     );
    CONSTRUCT_RT( KohnSham, dcomplex, double   );
    CONSTRUCT_RT( KohnSham, dcomplex, dcomplex );

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
    try {
      auto intAlg = input.getData<std::string>("RT.INTALG");

      if ( not intAlg.compare("MAGNUS2") ) { 
        rt->intScheme.intAlg = ExpMagnus2;
      }
      else if ( not intAlg.compare("MMUT") ) {
      }
      else {
          std::cout << "Could not understand RT.INTALG. Defaulting to MMUT.";
          std::cout << std::endl;
      }
    }
    catch(...) {
      std::cout << "Defaulting to MMUT integration algorithm" << std::endl;
    };

    // Get restart step if explicit leapfrog method
    if ( rt->intScheme.intAlg == MMUT ) {
      try {
        auto intRstrt = input.getData<std::string>("RT.RESTARTSTEP");
        if ( not intRstrt.compare("FORWARDEULER") ) {
          rt->intScheme.rstStep = ForwardEuler;
        }
        else if ( not intRstrt.compare("MAGNUS2") ) {
        }
        else {
          std::cout << "Could not understand RT.RESTARTSTEP. Defaulting to Magnus 2.";
          std::cout << std::endl;
        }
      }
      catch(...) {
        std::cout << "Defaulting to Magnus 2 restart for MMUT" << std::endl;
      }
    }



    // Set SCF perturbation
    rt->setSCFPerturbation( scfPert );
    // Inclusion of SCF perturbation
    OPTOPT(
      rt->intScheme.includeSCFField = input.getData<bool>("RT.SCFFIELD");
    )

    
    // MMUT Restart
    OPTOPT(
      rt->intScheme.iRstrt = input.getData<size_t>("RT.IRSTRT");
    )

    // Handle field specification
    try {

      // Get raw string from input
      std::string fieldSpec = input.getData<std::string>("RT.FIELD");
      std::istringstream fieldStream(fieldSpec);
 
      // Loop over field specification lines
      for(std::string fieldStr; std::getline(fieldStream, fieldStr); ) {
  
        // Split line on white space
        std::vector<std::string> tokens;
        split(tokens,fieldStr," \t");

        if( tokens.size() == 0 ) continue;


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


        // STEPFIELD
        if( envelope.find("STEPFIELD") != std::string::npos ) {

          // Determine if valid specifcation
          auto pStart = envelope.find("(");
          auto pEnd   = envelope.find(")");
          auto pSplit = envelope.find(",");


          if( pStart == std::string::npos or pEnd == std::string::npos
              or pSplit == std::string::npos )
            CErr(envelope + " not a valid STEPFIELD specification",out);

          envelope.erase(envelope.begin() + pEnd,envelope.end());
          envelope.erase(envelope.begin(), envelope.begin() + pStart+1);


          std::vector<std::string> tokens2;
          split(tokens2,envelope,",");

          if( tokens2.size() != 2 )
            CErr("STEPFIELD takes 2 arguements",out);

          double stepOn  = std::stod(tokens2[0]);
          double stepOff = std::stod(tokens2[1]);
   
          if( stepOff <= stepOn )
            CErr("STEPOFF must be > STEPON for STEPFIELD");


          // Append Field
          // XXX: Should store pointer to field base
          // and then append after envelope is determined
          rt->addField(fieldType, 
            StepField(stepOn,stepOff),
            DipoleField);


        } else CErr("Only STEPFIELD Implemented");
       
      }

    } catch( std::runtime_error &e ) {

      throw;

    } catch(...) { 

      out << "  *** Defaulting to Trivial Propagation from SCF Density ***\n";

    }

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

    // Amount of printing in the RT calc
    OPTOPT(
      rt->orbitalPopFreq = input.getData<size_t>("RT.ORBITALPOPFREQ");
    )

    // Whether to print time-dependent density
    OPTOPT(
      rt->printDen = input.getData<bool>("RT.PRINTDEN");
    )

    // Parse whether to print contraction timing during RT propagation
    OPTOPT( 
      rt->printContractionTiming = input.getData<bool>("RT.PRINTCONTRACTIONTIMING"); 
    )

    return rt;

  }; // CQRealTimeOpts

}; // namespace ChronusQ
