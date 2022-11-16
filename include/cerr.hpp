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
#pragma once

#include <chronusq_sys.hpp>
#include <libint2/cxxapi.h>
#include <cxxapi/boilerplate.hpp>

#define __CERR_RUNTIMEERR__ // Throw a runtime error on CErr

namespace ChronusQ {

  /**
   *  Standardized error handling.
   *
   *  Prints a message and properly cleans up the ChronusQ runtime
   */ 
  inline void CErr(const std::string &msg = "Die Die Die", 
    std::ostream &out = std::cout) {


    time_t currentTime;
    time(&currentTime); 

    if(MPIRank() == 0)
    out << msg << std::endl << "Job terminated: " << ctime(&currentTime)
        << std::endl;

#ifdef __CERR_RUNTIMEERR__
    throw std::runtime_error("FATAL");
#else
    finalize();
    exit(EXIT_FAILURE);
#endif
  };

  /**
   *  Standardized error handelling.
   *
   *  Prints a message and properly cleans up the ChronusQ runtime
   */ 
  inline void CErr(std::exception_ptr eptr, std::ostream &out = std::cout) {

    try{
      if(eptr) std::rethrow_exception(eptr);
    } catch( const std::exception & e) {
      std::stringstream ss;
      ss << "Caught \"" << e.what() << "\"" << std::endl;
      CErr(ss.str(),out);
    }

  };

};

