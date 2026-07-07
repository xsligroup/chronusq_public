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

  std::set<std::string> CQCUBE_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "DEN",
      "ORB",
      "MOS",
      "ROOTS",
      "POINTS",
      "STEPS",
      "NAME",
      "RES",
      "PADDING",
      "MAGANDPHASE"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  std::shared_ptr<CubeGen> CQCUBEOptions(std::ostream &out, CQInputFile &input,
    std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> &basis, EMPerturbation &emPert, double particleCharge) {

    std::string substring = "";
    // CUBE section not required
    if( not input.containsSection("CUBE") ) 
    {
      // Check in case [SUBSECTION/CUBE] is specified instead
      if(input.containsSection("SCF/CUBE"))
      {
        std::cout << "[CUBE] settings detected in [SCF/CUBE]" << std::endl;
        std::cout << "These settings are useed globally!" << std::endl;
        substring = "SCF.";
      }
      else if(input.containsSection("MCSCF/CUBE"))
      {
        std::cout << "[CUBE] settings detected in [MCSCF/CUBE]" << std::endl;
        std::cout << "These settings are useed globally!" << std::endl;
        substring = "MCSCF.";
      }
      else 
      {
        return nullptr;
      }
    }

    std::cout << " Found [CUBE] Section " << std::endl;

    // >>Keywords needed for constructor
    std::string resString = "COARSE";
    std::string spts;
    std::string ssteps;
    std::vector<std::string> nptstokens;
    std::vector<std::string> stepstokens;
    std::vector<size_t> npts;
    std::vector<double> steps;
    double pad = 0.0;

    // change resolution 
    OPTOPT( resString = input.getData<std::string>(substring+"CUBE/RES") );

    // change padding. Used to avoid cube cutoffs
    OPTOPT( pad = input.getData<double>(substring+"CUBE/PADDING") );

    // Create custom grid with points and stepsize 
    OPTOPT( spts = input.getData<std::string>(substring+"CUBE/POINTS") );
    OPTOPT( ssteps = input.getData<std::string>(substring+"CUBE/STEPS") );

    split(nptstokens, spts, " ,;");
    for (auto & npt: nptstokens)
      npts.push_back(std::stoi(npt));

    split(stepstokens, ssteps, " ,;");
    for (auto & nstep: stepstokens)
      steps.push_back(std::stod(nstep));

    bool customRes = false;
    // Handle strange cases
    if(!npts.empty()) {
        customRes = true;
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
          if( steps[0]*steps[1]*steps[2] < 0 ) CErr("Step size needs to be positive");
        }
    } else {
      if(!steps.empty()) CErr("Step size in grid also requires number of points");
    }

    if( !resString.empty() and !npts.empty() and !steps.empty() ) std::cout << "   ***WARNING: resolution overwrites custom grid" << std::endl; 
    if( abs(pad) != 0.0 and !npts.empty() and !steps.empty() ) std::cout << "   ***WARNING: padding is not used for custom grid construction" << std::endl;

    std::shared_ptr<CubeGen> cubeptr;

    // Construct cubegen
    if( !customRes and abs(pad) != 0.0 ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol,basis, resString, emPert, particleCharge, pad));

    } else if( !customRes ){

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol,basis, resString, emPert, particleCharge));

    } else if( !npts.empty() and !steps.empty() ){

      std::array<size_t,3> grid = {npts[0], npts[1], npts[2]};
      std::array<double,3> units = {steps[0], steps[1], steps[2]};
      cubeptr = std::make_shared<CubeGen>(CubeGen(mol, basis, grid, units, emPert, particleCharge));

    } else {

      cubeptr = std::make_shared<CubeGen>(CubeGen(mol, basis, emPert, particleCharge));

    }

    // >>Keywords not needed for constructor

    auto &cubeOptions = cubeptr->getCubeOptions();
    CQCUBEOptionalKeywords(out,input,cubeOptions,substring);

    return cubeptr;

  }; // CQCUBEOptions

  void handle_orbital_requests(std::ostream&out, CQInputFile & input, CubeGenOptions &cubeOpts, std::string subSection)
  {

    std::string OrbRequestString;
    OPTOPT(OrbRequestString = input.getData<std::string>(subSection+"CUBE/MOS"));
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
        CErr("Unrecognized token in " + subSection+"CUBE/MOS");
      }
    }
  }

  void handle_root_requests(std::ostream&out, CQInputFile & input, CubeGenOptions &cubeOpts, std::string subSection)
  {

    std::string RootRequestString;
    OPTOPT(RootRequestString = input.getData<std::string>(subSection+"CUBE/ROOTS"));
    if(!RootRequestString.empty() && subSection=="SCF/")
      CErr("Requesting multiple roots to generate cubes from an SCF calculation doesn't make sense!");
    // Default is all orbitals unless requested otherwise
    if(RootRequestString.empty() || RootRequestString=="GS")
    {
      cubeOpts.whichCIRoots = CI_CUBE_ROOT_CLASSES::GS;
      std::cout << "Generating Cubes for ONLY the lowest energy CI Root" << std::endl;
      return;
    }
    if(RootRequestString=="ALL")
    {
      cubeOpts.whichCIRoots = CI_CUBE_ROOT_CLASSES::ALL;
      std::cout << "Generating Cubes for ALL MCSCF Roots" << std::endl;
      return;
    }
    if(RootRequestString=="AVERAGE")
    {
      cubeOpts.whichCIRoots = CI_CUBE_ROOT_CLASSES::AVERAGE;
      std::cout << "Generating Cubes for The State Averaged Density" << std::endl;
      return;
    }

    // If no string matching, assume user requested a custom list
    cubeOpts.whichCIRoots = CI_CUBE_ROOT_CLASSES::CUSTOM;
    std::vector<std::string> RootRequestTokens;
    split(RootRequestTokens,RootRequestString,", ");
    for(auto & Root : RootRequestTokens)
    {
      try
      {
        // roots are 1 indexed
        std::vector<std::string> root2;
        split(root2,Root,"-");
        if(root2.size()==1)
        {
          size_t rootindex = std::stoul(Root);
          cubeOpts.addRoottoList(rootindex);
          std::cout << "Generating Cube for CI Root #" << Root << std::endl;
        }
        else if(root2.size()==2)
        {
          for(size_t i = std::stoul(root2[0]); i <= std::stoul(root2[1]); i++)
          {
            cubeOpts.addRoottoList(i);
            std::cout << "Generating Cube for CI Root #" << i << std::endl;
          }
        }
      }
      catch(...)
      {
        CErr("Unrecognized token in " + subSection+"CUBE/ROOTS");
      }
    }
  }

  // Handle keywords for an existing cube pointer
  void CQCUBEOptionalKeywords(std::ostream &out, CQInputFile &input, CubeGenOptions &cubeOpts, std::string subSection){


    // collects naming scheme from user (optional)
    OPTOPT( cubeOpts.cubeFileName = input.getData<std::string>(subSection+"CUBE/NAME") );

    // generate cube file for density
    OPTOPT( cubeOpts.denCube = input.getData<bool>(subSection+"CUBE/DEN") );
    if(cubeOpts.denCube)
    {
      handle_root_requests(out,input,cubeOpts,subSection);
    }

    // generate cube file for density
    OPTOPT( cubeOpts.orbCube = input.getData<bool>(subSection+"CUBE/ORB") );
    if(cubeOpts.orbCube)
    {
      handle_orbital_requests(out,input,cubeOpts,subSection);
    }

    // Handle Roots for MCSCF
    handle_root_requests(out,input,cubeOpts,subSection);

    // If Magnitude & Phase Cubes are requested
    OPTOPT( cubeOpts.MagnitudeAndPhase = input.getData<bool>(subSection+"CUBE/MAGANDPHASE") );

  }; //CQCUBEOptionalKeywords

  void ParseCubeSubsection(std::ostream &out, CQInputFile &input, std::string subsection,
    CubeGenOptions & cubeopts, std::shared_ptr<CubeGen> cube) {

    if( cube ){

      // Copy the internal cube options to the input cube options
      cubeopts = cube->getCubeOptions();

    }

    // Check if additional options are specified in [subsection/CUBE]
    if( not input.containsSection(subsection+"/CUBE") ) return;

    // Handle the case where [subsection/CUBE] is provided, but [CUBE] is not
    // In this case, cube is not initialized and cubes cannot be generated
    if(input.containsSection(subsection+"/CUBE") && !cube)
      CErr("Must specify Cube Settings in [CUBE] section!");

    std::cout << " Found [" << subsection << "/CUBE] section" << std::endl;
    std::set<std::string> invalidKeywords = CQCUBE_VALID(input.getSection(subsection+"/CUBE"));
    printInvalidKeys(invalidKeywords, subsection+"/CUBE");
    CQCUBEOptionalKeywords(out,input,cubeopts,subsection+"/");

  }; // ParseCubeSubsection


}; // namespace ChronusQ
