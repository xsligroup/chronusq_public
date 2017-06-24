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


#include <nonlinearopt/bfgs.hpp>

namespace ChronusQ
{

    class BFGSTest : public BFGS<double>
    {
    public:
        double *x;
        double *gx;
        double yMin;
        bool converged = false;

        BFGSTest() : BFGS<double>(2, 1E-18)
        {

            // Print Header output
            std::cout << "======================================================" << std::endl;
            std::cout << "               BFGS TEST:" << std::endl;
            std::cout << "  Optimization function = x^2 * y^2 - exp( -(x-1)^2 )" << std::endl;
            std::cout << "  Solving for x,y" << std::endl;
            std::cout << "======================================================" << std::endl;
            // Initialize 
            x = CQMemManager::get().malloc<double>(2);
            gx = CQMemManager::get().malloc<double>(2);
            BFGS<double>::setXPointer(x);
            BFGS<double>::setGPointer(gx);
            x[0] = 3.;
            x[1] = -2.;

            // Change optimization options
            this->printOptOutput = true;
            this->sdMaxDisp = 1.0;
            this->lineAlg = SteapestDescent;

            // Optimize
            converged = this->optimize();
            yMin = this->objFunc();
            this->gradFunc();

            // Print out data
            std::cout << std::scientific << std::setprecision(8)  << "Minimum = " << yMin << std::endl << std::endl;
            std::cout << "X:" << std::endl;
            std::cout << std::scientific << std::setprecision(8)  << x[0] << std::endl;
            std::cout << std::scientific << std::setprecision(8)  << x[1] << std::endl << std::endl;
            std::cout << "Gradient:" << std::endl;
            std::cout << std::scientific << std::setprecision(8) << gx[0] << std::endl;
            std::cout << std::scientific << std::setprecision(8)  << gx[1] << std::endl;

            std::cout << std::endl << std::endl;
        };
        BFGSTest(const BFGSTest &) = delete;
        BFGSTest(BFGSTest &&) = delete;

        // Simple test function for optmization
        // objFunc = x^2 * y^2 - exp(- (x-1)^2 * (y+1)^2 )
        // x = x[0]
        // y = x[1]
        inline double objFunc() { double y = (x[0]*x[0]*x[1]*x[1]) - std::exp( -(x[0]-1.)*(x[0]-1.) ); return y;};
        void gradFunc()
        {
            gx[0] = 2.*x[0]*x[1]*x[1] + 2.*(x[0]-1.)*std::exp(-(x[0]-1.)*(x[0]-1.) );
            gx[1] = 2.*x[1]*x[0]*x[0];
        }

        ~BFGSTest()
        {
            CQMemManager::get().free(x,gx);
        };
    };
}