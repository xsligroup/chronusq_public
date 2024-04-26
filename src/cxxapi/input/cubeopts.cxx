/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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
#include <regex>

namespace ChronusQ {

  void CQCUBE_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "CUBEDEN",
      "POINTS",
      "STEPS",
      "NAME",
      "RES",
      "PADDING"
    };

    // Specified keywords
    std::vector<std::string> cubeKeywords = input.getDataInSection("CUBE");

    // Make sure all of cubeKeywords in allowedKeywords

    for( auto &keyword : cubeKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword CUBE." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  std::shared_ptr<CubeGen> CQCUBEOptions(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> &ss) {

    // CUBE section not required
    if( not input.containsSection("CUBE") ) return nullptr;

    std::cout << " Found [CUBE] Section " << std::endl;

    // >>Keywords needed for constructor
    std::string CubegenFileName = "";
    std::string resString = "";
    size_t npts = 0;
    double step = 0.0;
    double pad = 0.0;

    // change resolution 
    OPTOPT( resString = input.getData<std::string>("CUBE.RES") );

    // change padding. Used to avoid cube cutoffs
    OPTOPT( pad = input.getData<double>("CUBE.PADDING") );

    // collects naming scheme from user (optional)
    OPTOPT( CubegenFileName = input.getData<std::string>("CUBE.NAME") );

    // Create custom grid with points and stepsize 
    // Assumes cube 
    OPTOPT( npts = input.getData<size_t>("CUBE.POINTS") );
    OPTOPT( step = input.getData<double>("CUBE.STEPS") );

    // Handle strange cases
    if( npts>0 and step==0.0 ) CErr("Number of points in grid also requires step size");
    if( npts==0 and step>0.0 ) CErr("Step size in grid also requires number of points");
    if( step<0.0 ) CErr("Step size needs to be positive");

    if( !resString.empty() and npts>0 and step>0.0 ) std::cout << "   ***WARNING: resolution overwrites custom grid" << std::endl; 
    if( abs(pad) != 0.0 and npts>0 and step>0.0 ) std::cout << "   ***WARNING: padding is not used for custom grid construction" << std::endl;

    std::shared_ptr<CubeGen> cubeptr;

    // Construct cubegen
    if( !resString.empty() and abs(pad) != 0.0 ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(ss, CubegenFileName, resString, pad));

    } else if( !resString.empty() ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(ss, CubegenFileName, resString));

    } else if( npts > 0 ){

      std::array<size_t,3> grid = {npts, npts, npts};
      std::array<double,3> units = {step, step, step};
      cubeptr = std::make_shared<CubeGen>(CubeGen(ss, CubegenFileName, grid, units));

    } else {

      cubeptr = std::make_shared<CubeGen>(CubeGen(ss));

    }

    // >>Keywords not needed for constructor

    // generate cube file for charge density
    OPTOPT( cubeptr->setDenEval(input.getData<bool>("CUBE.CUBEDEN")) );

    // OPTOPT( ss.scfControls.moCube = 
      // input.getData<bool>("CUBE.CUBEMO") );

    return cubeptr;

  }; // CQSCFOptions

}; // namespace ChronusQ
