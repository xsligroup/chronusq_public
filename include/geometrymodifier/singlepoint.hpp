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
#include <geometrymodifier.hpp>

namespace ChronusQ {

  class SinglePoint : public GeometryModifier {

    bool updated = false;

    public:

      SinglePoint() = default;

      bool hasNext() {
        return !updated;
      }

      void update(bool, Molecule&, bool, TDSCFOptions&, std::shared_ptr<SingleSlaterBase> ss, 
          EMPerturbation& emPert, std::vector<std::shared_ptr<CubeGen>> cubes) {
        updated = true;
      }
  };

} // namespace ChronusQ
