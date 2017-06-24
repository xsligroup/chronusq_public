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

#include <util/files.hpp>

#include <response/enums.hpp>
#include <response/settings.hpp>
#include <response/results.hpp>

namespace ChronusQ {
  
		struct ResponseBase {

    SafeFile savFile; ///< Data File

    ResponseSettings         genSettings;
    FDResponseSettings       fdrSettings;
    ResidueResponseSettings  resSettings;

    EMPerturbation           scfPert;

    FDObservables  fdObs;
    ResObservables resObs;

    ResponseBase( ResponseType job = RESIDUE) { genSettings.jobType = job; }

    ResponseBase( const ResponseBase & other ) :
      genSettings(other.genSettings),
      fdrSettings(other.fdrSettings),
      resSettings(other.resSettings) { }



    virtual size_t getNSingleDim( const bool doTDA = false) = 0;


    // Memory allocation

    inline void allocResults() {

      if( genSettings.jobType == RESIDUE ) allocResidueResults();
      else                                 allocFDRResults(); 

    };

    virtual void allocResidueResults() = 0; // Residue memory allocation
    virtual void allocFDRResults()     = 0; // FDR memory allocation



    // Procedural functions
    virtual void run() = 0;




    // Property Evaluation
    virtual void residueProperties() = 0;
    virtual void fdrProperties()     = 0;

    

    // Result Printing
    inline void printResults(std::ostream &out) {

      if( genSettings.jobType == RESIDUE ) printResidueResults(out);
      else                                 printFDRResults(out);

    };

    virtual void printResidueResults(std::ostream &out) = 0;
    virtual void printFDRResults(std::ostream &out)     = 0;



  };

};

