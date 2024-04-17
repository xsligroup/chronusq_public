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
      "RES"
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

  void CQCUBEOptions(std::ostream &out, CQInputFile &input,
    SingleSlaterBase &ss) {

    // CUBE section not required
    if( not input.containsSection("CUBE") ) return;

    std::cout << " CUBE INITIALIZED " << std::endl;
    // std::string chargeDensitySwitch;
    // OPTOPT(chargeDensitySwitch = input.getData<bool>("CUBE.CUBEDEN"); )

    // auto const regexCubeTrue = std::regex("true|on|1",std::regex_constants::icase);

    // // check for charge density 
    // if (std::regex_search(chargeDensitySwitch, regexCubeTrue)) {
    //   ss.scfControls.denCube = true;
    // }

    // Cube Options
    // generate cube file for charge density
    OPTOPT( ss.scfControls.denCube = input.getData<bool>("CUBE.CUBEDEN") );

    // collects NPTS from user for cube quality
    OPTOPT( ss.scfControls.res = input.getData<std::string>("CUBE.RES") );

    // collects NPTS from user for cube quality
    OPTOPT( ss.scfControls.npts = input.getData<size_t>("CUBE.POINTS") );

    OPTOPT( ss.scfControls.steps = input.getData<double>("CUBE.STEPS") );

    // collects naming scheme from user (optional)
    OPTOPT( ss.scfControls.CubegenFileName = input.getData<std::string>("CUBE.NAME") );

    // OPTOPT( ss.scfControls.moCube = 
      // input.getData<bool>("CUBE.CUBEMO") );


  }; // CQSCFOptions

}; // namespace ChronusQ
