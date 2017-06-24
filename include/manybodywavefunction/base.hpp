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

#include <util/math.hpp>
#include <detstringmanager.hpp>

namespace ChronusQ {

  /**
   *  \brief The ManyBodyWavefunctionBase class. The abstraction of information
   *  relating to the ManyBodyWavefunction class which are independent of storage
   *  type.
   *
   *
   *  See WaveFunction for further docs.
   */
  class ManyBodyWavefunctionBase {

  public:
    ManyBodyWavefunctionBase()                           = default;
    ManyBodyWavefunctionBase(const ManyBodyWavefunctionBase &) = default;
    ManyBodyWavefunctionBase(ManyBodyWavefunctionBase &&)      = default;

    /**
     *  ManyBodyWavefunctionBase Constructor. Constructs a ManyBodyWaveFunctionBase object
     *
     */
    // ManyBodyWavefunctionBase() {

    // }; // ManyBodyWavefunctionBase Constructor.

    virtual ~ManyBodyWavefunctionBase() = default;

  }; // class ManyBodyWavefunctionBase

}; // namespace ChronusQ
