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
#include <realtime.hpp>


namespace ChronusQ {

  template <template <typename, typename> class _SSTyp, typename IntsT>
  template <typename MatsT>
  void RealTime<_SSTyp,IntsT>::alloc() {

    // XXX: Member functions can't be partially specialized,
    //   so we're just going to cram this all into a single if statement...
    if( std::is_same<NEOSS<dcomplex,IntsT>,_SSTyp<dcomplex,IntsT>>::value ) {
      auto prop_c = dynamic_cast<NEOSS<dcomplex,IntsT>*>(&propagator_);
      auto map = prop_c->getSubsystemMap();
      auto order = prop_c->getOrder();

      assert( !map.empty() );

      // Loop over all subsystems
      for( auto& label: order ) {

        // Easier to read name for the subsystem
        auto& system = map[label];

        // Information for RealTime only
        size_t NB = system->onePDM->dimension();
        bool hasZ = system->onePDM->hasZ();
        bool hasXY= system->onePDM->hasXY();

        systems_.push_back(system.get());

        DOSav.push_back(
          std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(
            NB, hasXY, hasZ
          ));
        UH.push_back(
          std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(
            NB, hasXY, hasZ
          ));

      }

    }
    else {

      size_t NB = propagator_.onePDM->dimension();
      bool hasZ = propagator_.onePDM->hasZ();
      bool hasXY= propagator_.onePDM->hasXY();

      systems_.push_back(&propagator_);

      DOSav.emplace_back(
        std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(
          NB, hasXY, hasZ
        ));
      UH.emplace_back(
        std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(
          NB, hasXY, hasZ
        ));

    }

  };

}; // namespace ChronusQ


