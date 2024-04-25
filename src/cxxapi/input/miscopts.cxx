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

#include <util/mpi.hpp>
#include <util/threads.hpp>
#include <util/timer.hpp>
#include <util/print.hpp>

namespace ChronusQ {

  void CQMISC_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "MEM",
      "MEMBLK",
      "MEMTYPE",
      "NSMP",
      "TIMER",
      "DEBUGTIMING",
      "TIMERUNIT"
    };

    // Specified keywords
    std::vector<std::string> miscKeywords = input.getDataInSection("MISC");

    // Make sure all of miscKeywords in allowedKeywords
    for( auto &keyword : miscKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword MISC." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  void CQMiscOptions(std::ostream &out,
    CQInputFile &input) {

    size_t mem     = 256e6; // Default 256 MB allocation
    size_t blkSize = 2048;  // Default 2KB block size

    // Set max threads
    OPTOPT(
      SetNumThreads(input.getData<size_t>("MISC.NSMP"));
    )

    size_t nThreads = GetNumThreads();
    ProgramTimer::initialize("Chronus Quantum", nThreads);

    // Get timer options
    TimerOpts topts;
    OPTOPT( topts.doSummary = input.getData<bool>("MISC.TIMER"); )
    OPTOPT( topts.dbgPrint = input.getData<size_t>("MISC.DEBUGTIMING"); )
    std::string unit = "";
    OPTOPT( unit = input.getData<std::string>("MISC.TIMERUNIT"); )
    if (unit == "SECONDS" || unit == "SECOND" || unit == "S")
      topts.unit = SECONDS;
    else if ( unit == "MINUTES" || unit == "MINUTE" || unit == "MIN" ||
              unit == "M" )
      topts.unit = MINUTES;
    else if ( unit == "MILLISECONDS" || unit == "MILLISECOND" ||
              unit == "MILLI" || unit == "MS" )
      topts.unit = MILLISECONDS;
    else if ( unit != "" ) {
      std::string err = "Unrecognized unit " + unit + " for MISC.TIMERUNIT";
      CErr(err);
    }

    ProgramTimer::instance()->options = topts;

    // Determine if memory allocation was specified
    OPTOPT(
      std::string memStr = input.getData<std::string>("MISC.MEM");
      trim(memStr);

      size_t posKB = memStr.find("KB");
      size_t posMB = memStr.find("MB");
      size_t posGB = memStr.find("GB");

      if( posKB != std::string::npos ) {

        memStr.erase(posKB,2);
        trim(memStr);
        mem = std::stod(memStr) * 1e3;

      } else if( posMB != std::string::npos ) {

        memStr.erase(posMB,2);
        trim(memStr);
        mem = std::stod(memStr) * 1e6;

      } else if( posGB != std::string::npos ) {

        memStr.erase(posGB,2);
        trim(memStr);
        mem = std::stod(memStr) * 1e9;

      } else 
        mem = std::stod(memStr);

    )

    OPTOPT(blkSize = input.getData<size_t>("MISC.MEMBLK");)

    std::string memTypeStr = "DEFAULT";
    OPTOPT(
        memTypeStr = input.getData<std::string>("MISC.MEMTYPE");
        trim(memTypeStr);
    )
    CQMemBackendType memType = CQMemBackendType::PREALLOCATED;
    if (memTypeStr == "PREALLOCATED" or memTypeStr == "DEFAULT") {
      memType = CQMemBackendType::PREALLOCATED;
    } else if (memTypeStr == "OS"
            or memTypeStr == "OS_DIRECT"
            or memTypeStr == "OSDIRECT") {
      memType = CQMemBackendType::OS_DIRECT;
    } else {
      CErr("Unrecognized memory type " + memTypeStr);
    }


    out << "\n\n";

    std::pair<double, char> mem_postfix = memSize(mem);
    out << "  *** Allocating " << std::fixed << std::setprecision(1)
        << mem_postfix.first << " " << mem_postfix.second << "B *** \n";
    out << "  *** ChronusQ will use " << GetNumThreads() 
        << " OpenMP threads ***\n";
    out << "  *** ChronusQ will use " << MPISize() 
        << " MPI Processes ***\n\n";
    out << "\n\n";

    ProgramTimer::tick("Memory Allocation");
    CQMemManager::get().initialize(memType,mem,blkSize);
    ProgramTimer::tock();

  }; // CQMiscOptions

  std::pair<double, char> memSize(size_t mem) {

    std::string postfixes = " KMGT";
    size_t indx = std::floor(std::log10(mem))/3;
    if (mem == 0) indx = 0;
    else if (indx > 4) indx = 4;
    char postfix = postfixes.c_str()[indx];

    return std::make_pair(static_cast<double>(mem) / std::pow(1e3, indx), postfix);

  }


}; // namespace ChronusQ
