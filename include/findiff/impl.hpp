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

#include <findiff.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>


#define COEFF_DEBUG

namespace ChronusQ {

  template <typename RetT, typename ChangeT>
  void FiniteDifferencer<RetT,ChangeT>::makeTable() {

    // See Fornberg, B. (1988). Mathematics of computation, 51(184), 699-706.

    std::vector<double> alphas{0.};
    for (size_t ai = 1; ai < maxNodes_/2 + 1; ai++) {
      alphas.push_back(double(ai));
      alphas.push_back(-1*double(ai));
    }

#ifdef COEFF_DEBUG
    std::cout << "Alpha size: " << alphas.size() << std::endl;
    for ( size_t i = 0; i < alphas.size(); i++ )
      std::cout << "Alpha_" << i << ": " << alphas[i] << std::endl;
#endif

    // Allocate coefficients
    coeffs_.clear();
    for ( size_t mi = 0; mi <= order_; mi++ ) {
      coeffs_.push_back({});
      for ( size_t ni = 0; ni <= maxNodes_; ni++ ) {
        coeffs_[mi].push_back(std::vector<double>(ni+1, 0.));
      }
    }

    coeffs_[0][0][0] = 1.;
    double c1 = 1.;

    for ( size_t n = 1; n <= maxNodes_; n++ ) {
      
      double c2 = 1.;

      for ( size_t v = 0; v < n; v++ ) {
        
        double c3 = alphas[n] - alphas[v];
        c2 *= c3;

        for ( size_t m = 0; m <= std::min(n,order_); m++ ) {

          double tval = alphas[n] * coeffs_[m][n-1][v];
          if ( m != 0 ) tval -= m * coeffs_[m-1][n-1][v];
          
          coeffs_[m][n][v] = tval / c3;

        }

      }

      for ( size_t m = 0; m <= std::min(n,order_); m++ ) {

        double tval = -1 * alphas[n-1] * coeffs_[m][n-1][n-1];
        if ( m != 0 ) tval += m * coeffs_[m-1][n-1][n-1];

        coeffs_[m][n][n] = c1 / c2 * tval;
      }

      c1 = c2;

    }

#ifdef COEFF_DEBUG
    for ( size_t m = 0; m < coeffs_.size(); m++ )
    for ( size_t n = 0; n < coeffs_[m].size(); n++ )
    for ( size_t v = 0; v < coeffs_[m][n].size(); v++ ) {
      std::cout << "[" << m << "," << n << "," << v << "] :  ";
      std::cout << coeffs_[m][n][v] << std::endl;
    }
#endif

  } // FiniteDifferencer::makeTable

  
  template <typename RetT, typename ChangeT>
  void FiniteDifferencer<RetT, ChangeT>::doDifference() {

    std::cout << " *** Calculating numerical derivative of order " << order_;
    std::cout << " with accuracy of order " << acc_ << " ***" << std::endl;

    // Accuracy is generally n - m + 1, except for even orders on central diff
    //   where you can get away with one n less.
    maxNodes_ = acc_ + order_ + (order_ % 2) - 2;

    std::cout << "   * Required nodes: " << maxNodes_ + 1 - (order_ % 2);
    std::cout << std::endl;

    // Generate coefficients
    makeTable();

    // Temporary result storage
    std::vector<std::vector<RetT>> resVals;

    // Oth node only necessary for evens
    if ( order_ % 2 == 0 ) {
      for ( auto &f: changeFuncs_ ) {
        f(0.*h_);
        resVals.push_back( valFunc_() );
      }
    }
    else {
      for ( auto &f: changeFuncs_ )
        resVals.push_back({});
    }

    // Loop over nodes, then over changeFuncs
    for ( size_t inode = 1; inode <= maxNodes_; inode++ ) {

      ChangeT change = (inode + 1)/2 * h_;
      if ( inode % 2 == 0 ) change *= -1;

      for ( auto &f: changeFuncs_ ) {
        f(change);
        resVals.push_back( valFunc_() );
      }

    }


    // Nodes finished being computed, now combine
    // This can be further generalized (e.g. for matrices) by letting the
    //   user define this reduction function.
    auto& orderC = coeffs_[order_][maxNodes_];
    size_t nFuncs = changeFuncs_.size();

    results_ = std::vector<std::vector<RetT>>(nFuncs);

    for ( auto v = 0; v < orderC.size(); v++ ) {
      for ( auto ifunc = 0; ifunc < nFuncs; ifunc++ ) {

        // Deal with unitialized results
        if ( results_[ifunc].empty() ) {

          results_[ifunc].resize(resVals[v*nFuncs+ifunc].size());

          std::transform(resVals[v*nFuncs+ifunc].begin(),
                         resVals[v*nFuncs+ifunc].end(),
                         results_[ifunc].begin(),
                         [&](RetT val){return orderC[v] * val;});

        }
        else {

          std::transform(resVals[v*nFuncs+ifunc].begin(),
                         resVals[v*nFuncs+ifunc].end(),
                         results_[ifunc].begin(),
                         results_[ifunc].begin(),
                         [&](RetT val1, RetT val2) {
                           return orderC[v] * val1 + val2;
                         });
        }
      }
    }


    // Finally, rescale to the real grid size
    for ( auto& funcRes: results_ ) {
      std::transform(funcRes.begin(), funcRes.end(), funcRes.begin(),
        [&](RetT val){ return val / std::pow( h_, order_ ); });
    }

  };

}
