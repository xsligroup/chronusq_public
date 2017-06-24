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

namespace ChronusQ {

  // String constants

  constexpr char  bannerTop[] = "--------------------------------------------------------------------------------";
  constexpr char  bannerMid[] = "  ----------------------------------------------------------------------------  ";
  constexpr char  bannerEnd[] = "--------------------------------------------------------------------------------";
  constexpr char  BannerTop[] = "================================================================================";
  constexpr char  BannerMid[] = "--------------------------------------------------------------------------------";
  constexpr char  BannerEnd[] = "================================================================================";

  constexpr char CQBanner[] = 
"    ______ __                                      ____  \n"
"   / ____// /_   _____ ____   ____   __  __ _____ / __ \\ \n"
"  / /    / __ \\ / ___// __ \\ / __ \\ / / / // ___// / / / \n"
" / /___ / / / // /   / /_/ // / / // /_/ /(__  )/ /_/ /  \n"
" \\____//_/ /_//_/    \\____//_/ /_/ \\__,_//____/ \\___\\_\\  \n";


  /**
   *  Dump the contents of a file to a specified output device.
   *
   *  \param [in]  out    Output device
   *  \param [out] fName  Name of file to output.
   */ 
  inline void CatFile(std::ostream &out, std::string fName) {

    std::ifstream inFile(fName);
    std::string line;
    while(std::getline(inFile,line))
      out << line << std::endl;
    
  }; // CatFile

  /**
   *  Ouput the standard ChronusQ file header to a specified 
   *  output defice.
   *
   *  \param [in]  out    Output device
   */ 
  inline void CQOutputHeader(std::ostream &out) {

    time_t currentClockTime;
    time(&currentClockTime);

    out << "ChronusQ Job Started: " << ctime(&currentClockTime) << std::endl;
    out << CQBanner << std::endl;
    
    out << "Release Version: "  
       << ChronusQ_VERSION_MAJOR << "." 
       << ChronusQ_VERSION_MINOR << "."
       << ChronusQ_VERSION_PATCH << std::endl;

    out << std::endl << std::endl << "Contributors List:" << 
        std::endl << BannerTop << std::endl;
    CatFile(out,ChronusQ_AUTHOR_LIST);

  }; // CQOuputHeader

  /**
   *  Ouput the standard ChronusQ file footer to a specified 
   *  output defice.
   *
   *  \param [in]  out    Output device
   */ 
  inline void CQOutputFooter(std::ostream &out) {

    time_t currentClockTime;
    time(&currentClockTime);

    out << std::endl << std::endl;
    out << "ChronusQ Job Ended: " << ctime(&currentClockTime) << std::endl;

  }; // CQOutputFooter


  // Specific method for printing sections of relevance to CQ
  void printTimerSummary(std::ostream& out);

};


