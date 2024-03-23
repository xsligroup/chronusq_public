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
  
/* 
 * Addressing Array in this file follows 
 * lexicographic order
 *
 * This file is intended to create addressing
 * array at compile time 
 *
 */ 

template <typename T, std::size_t N>
using ArrayMat = std::array<std::array<T, N>, N>;

namespace {

/*
 * Explicit memoization for combination and bitWeight
 */

template <typename T, std::size_t I, std::size_t J>
constexpr T combination();

template <typename T, std::size_t I, std::size_t J>
constexpr T bitWeight();

template <typename T, std::size_t I, std::size_t J>
constexpr T comb = combination<T, I, J>();

template <typename T, std::size_t I, std::size_t J>
constexpr T bw = bitWeight<T, I, J>();

/* 
 *  Compute combinations 
 */ 
template <typename T, std::size_t I, std::size_t J>
constexpr T combination() {
  if constexpr (J == 0 or I == J) {
    return T(1);
  } else if (J > 0 and I > J) {
    constexpr T K = std::min(J, I - J);
    return comb<T, I - 1, K - 1> * I / T(K);
  } 
  // I < J
  return T(0);
}

/* 
 *  Compute bit weights for lexicographic order
 */ 
template <typename T, std::size_t I, std::size_t J> 
constexpr T bitWeight() {
  if constexpr (I < J) {
    return bw<T, I, J - 1> + comb<T, J - 1, I>;
  }
  // I >= J
  return T(0);
}

/* 
 *  Compute a row of weights for lexicographic order
 */ 
template <typename T, std::size_t I,
          std::size_t N,
          std::size_t ... Js>
constexpr std::array<T, N> bitWeightsOfARow(
    std::index_sequence<Js...>) {
  return { bitWeight<T, I, Js>() ... };
}

/* 
 *  Compute addressinng array  
 */ 
template <typename T, 
          std::size_t N, 
          std::size_t ... Is>
constexpr ArrayMat<T, N> computeAddressingArray(
    std::index_sequence<Is...>) {
  return { bitWeightsOfARow<T, Is, N>(std::make_index_sequence<N>{}) ... };
}

} // namespace 

/*
 * initialize addressing array
 */ 
template <typename T, std::size_t N = sizeof(T) * 8> 
struct BITADDRESSING {
  static constexpr ArrayMat<T, N> array 
    = computeAddressingArray<T, N>(std::make_index_sequence<N>{});
}; // struct BITADDRESSING

template struct BITADDRESSING<uint8_t>; 
template struct BITADDRESSING<uint16_t>;
template struct BITADDRESSING<uint32_t>;
template struct BITADDRESSING<uint64_t>;

} // namespace ChronusQ
