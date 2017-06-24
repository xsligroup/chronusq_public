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

#include <memory>

namespace ChronusQ {

  // T must have a default constructor
  //   (Concepts _should_ let us enforce this, eventually)
  // Can be inherited into other singletons, but be sure this is really what
  //   you want
  template <typename T>
  class Singleton {

    protected:

    Singleton(){};
    static std::unique_ptr<T> instance_;

    public:

    static T* instance() {
      if ( !instance_ ) instance_ = std::make_unique<T>();
      return instance_.get();
    }

  };
  template <typename T> std::unique_ptr<T> Singleton<T>::instance_;

} // namespace ChronusQ
