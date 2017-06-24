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

  void ParseSCFCubeSubsection(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<CubeGen> cube) {

    if( cube ){

      auto &cubeOptions = cube->getCubeOptions();
      ss->cubeOptsSS = cubeOptions;

    }

    // check if [SCF.CUBE] section
    if( not input.containsSection("SCF.CUBE") ) return;

    std::cout << " Found [SCF.CUBE] section" << std::endl;
    CQCUBE_VALID(out,input,"SCF.");
    CQCUBEOptionalKeywords(out,input,ss->cubeOptsSS,"SCF.");

  }

}
