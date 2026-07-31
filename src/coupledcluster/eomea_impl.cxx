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

#include <coupledcluster/EOMEA_1h2p.hpp>

namespace ChronusQ {

  template class EOMEA<dcomplex>;
  template class EOMEA<double>;


  template <typename MatsT>
  std::shared_ptr<EOMCCBase<MatsT>> build_EOMEA(
                                              const SafeFile &savFile,
                                              CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) { 


    std::shared_ptr<EOMCCBase<MatsT>> eomcc = nullptr;
    eomcc = 
      std::make_shared<EOMEA<MatsT>> (
        savFile,
        intermediates, eomSettings, ccSettings);
    return eomcc;
  }

  template std::shared_ptr<EOMCCBase<dcomplex>> build_EOMEA(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings);
  template std::shared_ptr<EOMCCBase<double>> build_EOMEA(
                                              const SafeFile &savFile,
                                              CCIntermediates<double> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings);
}

