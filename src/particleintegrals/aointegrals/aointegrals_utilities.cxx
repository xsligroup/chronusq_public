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

#include <cqlinalg.hpp>
#include <cqlinalg/blasutil.hpp>
#include <lapack.hh>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>

#include <util/threads.hpp>
#include <chrono>
// Debug directives
//#define _DEBUGORTHO
//#define __DEBUGERI__


namespace ChronusQ {



  /**
   *  \brief Allocate and evaluate the Schwarz bounds over the
   *  CGTO shell pairs.
   */ 
//  template <>
//  void AOIntegrals<dcomplex>::computeSchwarz() {
//    CErr("Only real GTOs are allowed",std::cout);
//  };
  template <typename IntsT>
  void DirectTPI<IntsT>::computeSchwarz() {

    if( schwarz() != nullptr ) CQMemManager::get().free(schwarz());
    if( schwarz2() != nullptr ) CQMemManager::get().free(schwarz2());

    // Allocate the schwarz tensor
    size_t nShell = basisSet().nShell;
    schwarz() = CQMemManager::get().malloc<double>(nShell*nShell);
    if (&basisSet() != &basisSet2())
      schwarz2() = CQMemManager::get().malloc<double>(basisSet2().nShell*basisSet2().nShell);

    // Define the libint2 integral engine
    libint2::Engine engine(libint2::Operator::coulomb,
      std::max(basisSet().maxPrim, basisSet2().maxPrim),
      std::max(basisSet().maxL,basisSet2().maxL), 0);

    engine.set_precision(0.); // Don't screen prims during evaluation

    const auto &buf_vec = engine.results();

    auto topSch = std::chrono::high_resolution_clock::now();
  
    size_t n1,n2;
    for(auto s1(0ul); s1 < basisSet().nShell; s1++) {
      n1 = basisSet().shells[s1].size(); // Size shell 1
    for(auto s2(0ul); s2 <= s1; s2++) {
      n2 = basisSet().shells[s2].size(); // Size shell 2



      // Evaluate the shell quartet (s1 s2 | s1 s2)
      engine.compute(
        basisSet().shells[s1],
        basisSet().shells[s2],
        basisSet().shells[s1],
        basisSet().shells[s2]
      );

      if(buf_vec[0] == nullptr) continue;

      // Allocate space to hold the diagonals
      double* diags = CQMemManager::get().malloc<double>(n1*n2);

      for(auto i(0), ij(0); i < n1; i++)
      for(auto j(0); j < n2; j++, ij++)
        diags[i + j*n1] = buf_vec[0][ij*n1*n2 + ij];


      schwarz()[s1 + s2*basisSet().nShell] =
        std::sqrt(lapack::lange(lapack::Norm::Inf,n1,n2,diags,n1));

      // Free up space
      CQMemManager::get().free(diags);

    } // loop s2
    } // loop s1

    if (&basisSet() != &basisSet2()) {
      // compute (rs|rs)
      for(auto s1(0ul); s1 < basisSet2().nShell; s1++) {
        n1 = basisSet2().shells[s1].size(); // Size shell 1
      for(auto s2(0ul); s2 <= s1; s2++) {
        n2 = basisSet2().shells[s2].size(); // Size shell 2



        // Evaluate the shell quartet (s1 s2 | s1 s2)
        engine.compute(
          basisSet2().shells[s1],
          basisSet2().shells[s2],
          basisSet2().shells[s1],
          basisSet2().shells[s2]
        );

        if(buf_vec[0] == nullptr) continue;

        // Allocate space to hold the diagonals
        double* diags = CQMemManager::get().malloc<double>(n1*n2);

        for(auto i(0), ij(0); i < n1; i++)
        for(auto j(0); j < n2; j++, ij++)
          diags[i + j*n1] = buf_vec[0][ij*n1*n2 + ij];


        schwarz2()[s1 + s2*basisSet2().nShell] =
          std::sqrt(lapack::lange(lapack::Norm::Inf,n1,n2,diags,n1));

        // Free up space
        CQMemManager::get().free(diags);

      } // loop s2
      } // loop s1
      // done computing (rs|rs)
    }

    auto botSch = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> durSch = botSch - topSch;

    HerMat('L',basisSet().nShell,schwarz(),basisSet().nShell);
    if (&basisSet() != &basisSet2())
      HerMat('L',basisSet2().nShell,schwarz2(),basisSet2().nShell);

#if 0
    prettyPrintSmart(std::cout,"Schwarz",schwarz,basisSet_.nShell,
      basisSet_.nShell,basisSet_.nShell);
#endif

  }; // DirectERI<double>::computeSchwarz
  template void DirectTPI<double>::computeSchwarz();
  template void DirectTPI<dcomplex>::computeSchwarz();

}; // namespace ChronusQ

