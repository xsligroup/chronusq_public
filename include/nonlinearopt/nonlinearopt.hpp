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

#include <chronusq_sys.hpp>

namespace ChronusQ
{

    class NonLinearOpt
    {
    public:
        size_t nStep = 300;        // Number of optimization steps to perform
        bool printOptOutput = false;  // Whether to print output for function optimization
        bool printVariables = false;       // Whether to print variables at each step
        bool printGrad = false;    // Whether to print the gradient at each step

        virtual double objFunc() = 0;  // Function that defines the objective function to be optimized
        virtual void gradFunc() = 0; // Function that computes the gradient of the objective function
        virtual bool optimize() = 0;  // Optimization function
    };
}