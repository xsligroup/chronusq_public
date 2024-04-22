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


#include <nonlinearopt/goldenSection.hpp>

namespace ChronusQ
{

    class GoldenSectionTest : public GoldenSectionSearch<double>
    {
    public:
        double *x;
        double *dx;
        double *initX;
        double yMin;
        double aMin;

        GoldenSectionTest() : GoldenSectionSearch<double>(1, 1E-18)
        {
            // Header
            std::cout << "======================================================" << std::endl;
            std::cout << "         Golden Section Search TEST:" << std::endl;
            std::cout << "  Optimization function = x^2" << std::endl;
            std::cout << "  Solving for x" << std::endl;
            std::cout << "======================================================" << std::endl;
            // Initialize 
            x = CQMemManager::get().malloc<double>(1);
            dx = CQMemManager::get().malloc<double>(1);
            initX = CQMemManager::get().malloc<double>(1);
            GoldenSectionSearch<double>::setXPointer(x);
            x[0] = -3.;
            dx[0] = 1.;
            initX[0] = -3.;
            GoldenSectionSearch<double> :: setXPointer(x);

            // Change optimization options
            this->printLineSearch = true;

            // Optimize
            this->optimize(initX,dx,0.,6.,yMin,aMin);

            // Print out data
            std::cout << "Minimum = " << yMin << std::endl << std::endl;
            std::cout << "X:" << std::endl;
            std::cout << std::scientific << std::setprecision(8) << x[0] << std::endl;
            std::cout << std::endl << std::endl;
        };
        GoldenSectionTest(const BFGSTest &) = delete;
        GoldenSectionTest(BFGSTest &&) = delete;

        // Simple test function for optmization
        // objFunc = x^2 
        inline double objFunc() { double y = x[0]*x[0]; return y; }
        void gradFunc()
        {
            double gx = 2.*x[0];
        }

        ~GoldenSectionTest()
        {
            CQMemManager::get().free(x,initX,dx);
        };
    };
}