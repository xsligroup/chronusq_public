/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/giaodirecteri.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>
#include <typeinfo>
#include <memory>


namespace ChronusQ {

  /**
   *  \brief The pointer convertor. This static function converts
   *  the underlying polymorphism correctly to hold a different
   *  type of matrices. It is called when the corresponding
   *  SingleSlater object is being converted.
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  std::shared_ptr<TPIContractions<MatsU,IntsT>>
  TPIContractions<MatsT,IntsT>::convert(const std::shared_ptr<TPIContractions<MatsT,IntsT>>& ch) {

    if (not ch) return nullptr;

    const std::type_info &tID(typeid(*ch));

    if (tID == typeid(InCore4indexTPIContraction<MatsT,IntsT>)) {
      return std::make_shared<InCore4indexTPIContraction<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<InCore4indexTPIContraction<MatsT,IntsT>>(ch));

    } else if (tID == typeid(GTODirectTPIContraction<MatsT,IntsT>)) {
      return std::make_shared<GTODirectTPIContraction<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<GTODirectTPIContraction<MatsT,IntsT>>(ch));

    } else if (tID == typeid(GTODirectRelERIContraction<MatsT,IntsT>)) {
      return std::make_shared<GTODirectRelERIContraction<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ch));

    } else if (tID == typeid(GIAODirectERIContraction)) {
      return std::dynamic_pointer_cast<TPIContractions<MatsU,IntsT>>(ch);

    } else if (tID == typeid(InCoreRITPIContraction<MatsT,IntsT>)) {
      return std::make_shared<InCoreRITPIContraction<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<InCoreRITPIContraction<MatsT,IntsT>>(ch));

    } else {
      std::stringstream errMsg;
      errMsg << "TPIContractions implementation \"" << tID.name() << "\" not registered in convert." << std::endl;
      CErr(errMsg.str(),std::cout);
    }

    return nullptr;

  }

}
