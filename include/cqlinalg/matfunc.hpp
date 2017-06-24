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

#include <cqlinalg/cqlinalg_config.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/factorization.hpp>

namespace ChronusQ {

  template <typename F, typename _F1, typename _F2>
  void MatDiagFunc(const F &func, size_t N, _F1 *A, size_t LDA, _F2 *B,
                   size_t LDB);

  template <typename _FExp, typename _F1, typename _F2>
  void MatExp(char ALG, size_t N, _FExp ALPHA, _F1 *A, size_t LDA, 
    _F2 *ExpA, size_t LDEXPA);

  template <typename _F1, typename _F2, typename _FC>
  void MatSeries(size_t NC, size_t N, _F1 *A, size_t LDA, _F2 *B,
    size_t LDB, _FC *C);

    /*
     * Computes the unitary Matrix Exponential (ExpA) Using Taylor Expansion
     * from an antisymmetric matrix (A)
     */
  template <typename MatsU>
  void MatExp( size_t N, MatsU *A, size_t LDA,
    MatsU *ExpA, size_t LDEXPA);


}; // namespace ChronusQ

