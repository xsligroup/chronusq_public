/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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
  
namespace {  
/* Bit Masks in this file:
 *
 *   Type A: 0...0 1 0...0 -> bit at ith  
 *   Type B: 0...0 1 1...1 -> bits previous to ith, not including ith
 *   Type C: 1...1 1 0...0 -> bits after ith, including ith 
 *
 */             

/*
 * \brief generating bit mask type A at compile time
 *    sequence start from 00...01 to 10...00
 */ 
template <typename T, 
          std::size_t N = sizeof(T) * 8, 
          std::size_t ... Is>
constexpr std::array<T, N> genBitAtIthMasks(
    std::index_sequence<Is...>) {
  return { T(1) << Is ... };
}

/*
 * \brief generating bit mask type B at compile time
 *    sequence start from 00...00  to 11...11
 */ 
template <typename T, 
          std::size_t N = sizeof(T) * 8, 
          std::size_t ... Is>
constexpr std::array<T, N+1> genBitsPrevToIthMasks(
    std::index_sequence<Is...>) {
  constexpr T DETST_MAX = std::numeric_limits<T>::max();
  return { T(0),  static_cast<T>(DETST_MAX >> (N - Is - 1)) ... };
}

/*
 * \brief generating bit mask type C at compile time
 *    sequence start from 11...11 to 00...00
 */ 
template <typename T, 
          std::size_t N = sizeof(T) * 8, 
          std::size_t ... Is>
constexpr std::array<T, N+1> genBitsAfterIthBitMasks(
    std::index_sequence<Is...>) {
  constexpr T DETST_MAX = std::numeric_limits<T>::max();
  return { static_cast<T>(DETST_MAX << Is) ..., T(0) };
}

} // namespace 

/*
 * \brief the struct that holds all the masks
 *
 */ 
template <typename T, std::size_t N = sizeof(T) * 8>
struct BITMASKS {
  static constexpr std::array<T, N> bitAt
    = genBitAtIthMasks<T>(std::make_index_sequence<N>{});
  static constexpr std::array<T, N+1> bitsPrevTo
    = genBitsPrevToIthMasks<T>(std::make_index_sequence<N>{});
  static constexpr std::array<T, N+1> bitsAfter
    = genBitsAfterIthBitMasks<T>(std::make_index_sequence<N>{});
}; // struct BITMASKS

template struct BITMASKS<uint8_t>;
template struct BITMASKS<uint16_t>;
template struct BITMASKS<uint32_t>;
template struct BITMASKS<uint64_t>;

}; // namespace ChronusQ
