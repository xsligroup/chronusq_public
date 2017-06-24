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

#include <functional>
#include <vector>

namespace ChronusQ {

  template <typename RetT, typename ChangeT>
  class FiniteDifferencer {


    // Function that returns a vector of values to take the difference over
    std::function<std::vector<RetT>()> valFunc_;

    // Functions that change paramaters necessary for a new call to valFunc_
    std::vector<std::function<void(ChangeT)>> changeFuncs_;

    size_t order_ = 0;  ///< Derivative order
    size_t acc_ = 0;    ///< Accuracy order

    size_t maxNodes_ = 0; ///< Maximum nodes required

    // Results
    std::vector<std::vector<RetT>> results_;

    // Table of coefficents - m, n, v indexed where:
    //   m = derivative level
    //   n = accuracy level
    //   v = node
    std::vector<std::vector<std::vector<double>>> coeffs_;

    ChangeT h_; ///< Size of change to give to changeFuncs_

    // Internal procedural
    void makeTable();

    public:


    FiniteDifferencer() = default;
    FiniteDifferencer(const FiniteDifferencer&) = default;
    FiniteDifferencer(FiniteDifferencer&&) = default;


    // Setters
    void setValFunc(std::function<std::vector<RetT>()> func) {
      valFunc_ = func;
    };

    void setChangeFunc(std::vector<std::function<void(ChangeT)>> func) {
      changeFuncs_ = func;
    };

    void addChangeFunc(std::function<void(ChangeT)> func) {
      changeFuncs_.push_back(func);
    };


    void setDerOrder(size_t num) { order_ = num; };
    void setAccOrder(size_t num) { acc_ = num; };
    void setStepSize(ChangeT size) { h_ = size; };


    // Getters
    size_t getDerOrder() { return order_; };
    size_t getAccOrder() { return acc_; };
    ChangeT getStepSize() { return h_; };
    std::vector<std::vector<RetT>> getResults() { return results_; }


    // Procedural
    void doDifference();


  }; // class FiniteDifferencer

} // namespace ChronusQ 

#include <findiff/impl.hpp>
