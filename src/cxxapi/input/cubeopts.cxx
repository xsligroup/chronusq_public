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

  void CQCUBE_VALID( std::ostream &out, CQInputFile &input, std::string section ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "DEN",
      "ORB",
      "MOS",
      "POINTS",
      "STEPS",
      "NAME",
      "RES",
      "PADDING",
      "MAGANDPHASE"
    };

    std::string subSection = section + "CUBE";

    // Specified keywords
    std::vector<std::string> cubeKeywords = input.getDataInSection(subSection);

    // Make sure all of cubeKeywords in allowedKeywords

    for( auto &keyword : cubeKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword " + subSection + "." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  std::shared_ptr<CubeGen> CQCUBEOptions(std::ostream &out, CQInputFile &input,
    std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> &basis) {

    // CUBE section not required
    if( not input.containsSection("CUBE") ) return nullptr;

    std::cout << " Found [CUBE] Section " << std::endl;

    // >>Keywords needed for constructor
    std::string resString = "";
    std::string spts;
    std::string ssteps;
    std::vector<std::string> nptstokens;
    std::vector<std::string> stepstokens;
    std::vector<size_t> npts;
    std::vector<double> steps;
    double pad = 0.0;

    // change resolution 
    OPTOPT( resString = input.getData<std::string>("CUBE.RES") );

    // change padding. Used to avoid cube cutoffs
    OPTOPT( pad = input.getData<double>("CUBE.PADDING") );

    // Create custom grid with points and stepsize 
    OPTOPT( spts = input.getData<std::string>("CUBE.POINTS") );
    OPTOPT( ssteps = input.getData<std::string>("CUBE.STEPS") );

    split(nptstokens, spts, " ,;");
    for (auto & npt: nptstokens)
      npts.push_back(std::stoi(npt));

    split(stepstokens, ssteps, " ,;");
    for (auto & nstep: stepstokens)
      steps.push_back(std::stod(nstep));

    // Handle strange cases
    if(!npts.empty()) {
        if(steps.empty()) {
          CErr("Number of points in grid also requires step size");
        } else {
          // cubic
          if (npts.size()==1) {
            npts.push_back(npts[0]);
            npts.push_back(npts[0]); 
          }
          if (steps.size()==1) {
            steps.push_back(steps[0]);
            steps.push_back(steps[0]); 
          }
          // strange cases
          if( npts[0]*npts[1]*npts[2]<= 0 ) CErr("Step needs to be positive"); 
          if( steps[0]*steps[1]*steps[2]<= 0 ) CErr("Step size needs to be positive"); 
        }
    } else {
      if(!steps.empty()) CErr("Step size in grid also requires number of points");
    }

    if( !resString.empty() and !npts.empty() and !steps.empty() ) std::cout << "   ***WARNING: resolution overwrites custom grid" << std::endl; 
    if( abs(pad) != 0.0 and !npts.empty() and !steps.empty() ) std::cout << "   ***WARNING: padding is not used for custom grid construction" << std::endl;

    std::shared_ptr<CubeGen> cubeptr;

    // Construct cubegen
    if( !resString.empty() and abs(pad) != 0.0 ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol,basis, resString, pad));

    } else if( !resString.empty() ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol,basis, resString));

    } else if( !npts.empty() and !steps.empty() ){

      std::array<size_t,3> grid = {npts[0], npts[1], npts[2]};
      std::array<double,3> units = {steps[0], steps[1], steps[2]};
      cubeptr = std::make_shared<CubeGen>(CubeGen(mol, basis, grid, units));

    } else {

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol, basis));

    }

    // >>Keywords not needed for constructor

    auto &cubeOptions = cubeptr->getCubeOptions();
    CQCUBEOptionalKeywords(out,input,cubeOptions,"");

    return cubeptr;

  }; // CQCUBEOptions

  void handle_orbital_requests(std::ostream&out, CQInputFile & input, CubeGenOptions &cubeOpts, std::string subSection)
  {

    std::string OrbRequestString;
    OPTOPT(OrbRequestString = input.getData<std::string>(subSection+"CUBE.MOS"));
    // Default is all orbitals unless requested otherwise
    if(OrbRequestString.empty() || OrbRequestString=="ALL")
    {
      cubeOpts.whichMO = MO_CLASSES::ALL;
      std::cout << "Generating Cubes All Orbitals" << std::endl;
      return;
    }
    // Future options handled here

    // If no string matching, assume user requested a custom list
    cubeOpts.whichMO = MO_CLASSES::CUSTOM;
    std::vector<std::string> OrbRequestTokens;
    split(OrbRequestTokens,OrbRequestString,", ");
    for(auto & mo : OrbRequestTokens)
    {
      try
      {
        // moindex is 1 indexed
        std::vector<std::string> mo2;
        split(mo2,mo,"-");
        if(mo2.size()==1)
        {
          size_t moindex = std::stoul(mo);
          cubeOpts.addMOtoList(moindex);
          std::cout << "Generating Cube for Orbital #" << mo << std::endl;
        }
        else if(mo2.size()==2)
        {
          for(size_t i = std::stoul(mo2[0]); i <= std::stoul(mo2[1]); i++)
          {
            cubeOpts.addMOtoList(i);
            std::cout << "Generating Cube for Orbital #" << i << std::endl;
          }
        }
      }
      catch(...)
      {
        CErr("Unrecognized token in " + subSection+"CUBE.MOS");
      }
    }
  }

  // Handle keywords for an existing cube pointer
  void CQCUBEOptionalKeywords(std::ostream &out, CQInputFile &input, CubeGenOptions &cubeOpts, std::string subSection){


    // collects naming scheme from user (optional)
    OPTOPT( cubeOpts.cubeFileName = input.getData<std::string>(subSection+"CUBE.NAME") );

    // generate cube file for density
    OPTOPT( cubeOpts.denCube = input.getData<bool>(subSection+"CUBE.DEN") );

    // generate cube file for density
    OPTOPT( cubeOpts.orbCube = input.getData<bool>(subSection+"CUBE.ORB") );
    if(cubeOpts.orbCube)
    {
      handle_orbital_requests(out,input,cubeOpts,subSection);
    }

    // If Magnitude & Phase Cubes are requested
    OPTOPT( cubeOpts.MagnitudeAndPhase = input.getData<bool>(subSection+"CUBE.MAGANDPHASE") );

  }; //CQCUBEOptionalKeywords

}; // namespace ChronusQ
