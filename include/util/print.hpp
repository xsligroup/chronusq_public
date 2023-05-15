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
#include <utility>
#include <string>
#include <cmath>

namespace ChronusQ {

  inline void FormattedLine(std::ostream &out, std::string s) {
    out << std::setw(45) << "  " + s << std::endl;
  }

  template <typename T>
  void FormattedLine(std::ostream &out, std::string s, T v) {
    out << std::setw(45) << "  " + s << v << std::endl;
  }

  template <typename T, typename U>
  void FormattedLine(std::ostream &out, std::string s, T v, U u) {
    out << std::setw(45) << "  " + s << v << u << std::endl;
  }

  template <typename T, typename U>
  void FormattedLine(std::ostream &out, std::string s1, T v, std::string s2, U u) {
    out << std::setw(45) << "  " + s1 << v << "  " + s2 << u << std::endl;
  }

  std::pair<double, char> memSize(size_t mem);
 
}; // namespace ChronusQ
