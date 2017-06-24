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
#include <util/mpi.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQBASIS_VALID( std::ostream &out, CQInputFile &input, std::string section) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "FORCECART",
      "BASIS",
      "BASISTYPE",
      "BASISDEF",
      "DEFINEBASIS"
    };

    // Specified keywords
    std::vector<std::string> basisKeywords = input.getDataInSection(section);

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : basisKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword " + section + "." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  Construct a basis set object using the input file.
   *
   *  \param [in] out   Output device for data output
   *  \param [in] input Input file datastructure
   *  \param [in] mol   Molecule object from which to construct the BasisSet object
   *
   *  \returns Appropriate BasisSet object for the input parameters
   */
  std::shared_ptr<BasisSet> CQBasisSetOptions(std::ostream &out, CQInputFile &input,
    Molecule &mol, std::string section) {

    // Determine if we're forcing cartesian functions
    bool forceCart(false);
    OPTOPT( forceCart = input.getData<bool>(section+".FORCECART"); );

    // Determine if we're parsing the basis from a basis file or the input file
    bool inputBasis(false);
    OPTOPT( inputBasis = input.getData<bool>(section+".DEFINEBASIS") )

    // Determine if we're parsing a protonic basis
    bool pBasis(false);
    OPTOPT( pBasis = (!section.compare("PBASIS")); ); 

    // Find the Basis Definition
    std::string basisName;
    OPTOPT( basisName = input.getData<std::string>(section+".BASIS"); )
    std::string basisDef;
    if ( inputBasis )
      OPTOPT( basisDef = input.getData<std::string>(section+".BASISDEF"); )

    // Check for consistency
    if ( basisName.empty() and basisDef.empty() ) {
      if ( section == "BASIS" )
        CErr("Basis file or specification not found!");
      else if ( section == "PBASIS" )
        CErr("Proton Basis file or specification not found!");
      else
        return std::make_shared<BasisSet>();
    }


    BASIS_FUNCTION_TYPE bType = REAL_GTO;
    try{
      std::string bTypeString = input.getData<std::string>(section+".BASISTYPE");
      if(not bTypeString.compare("GIAO")) bType = COMPLEX_GIAO;
      else if(not bTypeString.compare("GTO")) bType = REAL_GTO;
      else CErr(section+".BASISTYPE not valid.",out);
    } catch(...) {;}

    // Construct the BasisSet object
    std::shared_ptr<BasisSet> basis =
        std::make_shared<BasisSet>(basisName,basisDef,inputBasis,
                                   mol,bType,forceCart,MPIRank() == 0,
                                   pBasis);

    // Ouput BasisSet information
    out << *basis << std::endl;

    return basis; // Return BasisSet object (no intermediates)

  }; // CQBasisSetOptions


}; // namespace ChronusQ
