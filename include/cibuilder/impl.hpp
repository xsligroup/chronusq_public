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

#include <cibuilder.hpp>
#include <cibuilder/casci/impl.hpp>
#include <cibuilder/rasci/impl.hpp>

// #define _DEBUG_CIBuilder_IMPL

namespace ChronusQ {

  // Different type
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  CIBuilder<MatsT,IntsT>::CIBuilder(const CIBuilder<MatsU,IntsT> & other) {};
  
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  CIBuilder<MatsT,IntsT>::CIBuilder(CIBuilder<MatsU,IntsT> && other) {};

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  std::shared_ptr<CIBuilder<MatsU,IntsT>>
  CIBuilder<MatsT,IntsT>::convert(const std::shared_ptr<CIBuilder<MatsT,IntsT>> & ci) {

    if (not ci) return nullptr;

    const std::type_info & tID(typeid(*ci));
    
    std::shared_ptr<CIBuilder<MatsU,IntsT>> output_ptr;

    if (tID == typeid(CASCI<MatsT,IntsT>)) {
      output_ptr = std::dynamic_pointer_cast<CIBuilder<MatsU, IntsT>>(
                     std::make_shared<CASCI<MatsU,IntsT>>(
                      *std::dynamic_pointer_cast<CASCI<MatsT,IntsT>>(ci)));

    } else if (tID == typeid(RASCI<MatsT,IntsT>)) {
      output_ptr = std::dynamic_pointer_cast<CIBuilder<MatsU, IntsT>>(
                     std::make_shared<RASCI<MatsU,IntsT>>(
                       *std::dynamic_pointer_cast<RASCI<MatsT,IntsT>>(ci)));
    }

    return output_ptr;
  } //CIBuilder::convert

}; // namespace ChronusQ
