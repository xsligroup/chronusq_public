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
   * Check valid MP2 keywords
   * 
   */
  void CQMP2_VALID(std::ostream &out, CQInputFile & input)
  {
    if( not input.containsSection("MP2")) return; 

    std::vector<std::string> allowedKeywords = {
      "NOS",
    };

    // Specified keywords
    std::vector<std::string> MP2Keywords = input.getDataInSection("MP2");

    // Make sure all of the basis keywords in allowed keywords
    for( auto &keyword : MP2Keywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() )
      CErr("Keyword MP2." + keyword + " is not recognized",std::cout);// Error
    }

    // MP2 shouldn't be run on DFT references
    if( input.containsData("QM.REFERENCE") ) {

      std::string ref = input.getData<std::string>("QM.REFERENCE");

      bool isKS  = not (ref.find("HF") != std::string::npos);

      if( isKS )
        CErr("MP2 + KS not allowed");
    }
  }; // CQMP2_VALID


  std::shared_ptr<MP2Base> CQMP2Options(std::ostream & out, CQInputFile & input, 
    std::shared_ptr<SingleSlaterBase> & ss)
    {
      std::shared_ptr<MP2Base> mp;
      
      try{
        mp = std::make_shared<MP2<double,double>>(std::dynamic_pointer_cast<SingleSlater<double,double>>(ss));
      }
      catch(...)
      {
        CErr("MP2 only implemented and tested for real MatsT & IntsT");
      }

      // Attempt to make a NEO-MP2 object
      if(std::dynamic_pointer_cast<NEOSS<double,double>>(ss))
      {
        std::cout << "Running NEOMP2" << std::endl;
        mp = std::make_shared<NEOMP2<double,double>>(std::dynamic_pointer_cast<NEOSS<double,double>>(ss));
      }

      OPTOPT(mp->makeMP2NOs=input.getData<bool>("MP2.NOS"));

      return mp;
    }; // CQMP2Options

}; // namespace ChronusQ
