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


#include <interpolate.hpp>

namespace ChronusQ
{

    class InterpolateTest 
    {
    public:
        size_t nInter = 2;
        std::vector<double> energy;
        double *B;
        double convergence = 1.E-18;
        bool converged = false;
        std::vector<double> coeffs;
        double *gOpt;


        InterpolateTest() : energy(nInter,0.), coeffs(nInter,0.)
        {
            // Print Header output
            std::cout << "======================================================" << std::endl;
            std::cout << "               INTERPOLATE TEST:" << std::endl;
            std::cout << "  Optimization function = E*c - 1/2 c*B*c" << std::endl;
            std::cout << "  Solving for the c vector" << std::endl;
            std::cout << "  E = [0,1]      B = [0., 2.] " << std::endl;
            std::cout << "                     [2., 0.] " << std::endl;
            std::cout << "======================================================" << std::endl;
            // Initialize 
            B = CQMemManager::get().malloc<double>(nInter*nInter);
            energy[1] = 1.;
            for(size_t i=0; i<nInter*nInter; i++ )
                B[i] = 0.;
            B[1] = 2.;
            B[2] = 2.;

            ENERGYDIIS<double> interp(nInter, B, nInter, energy, 0, convergence);
            

            // Change optimization options
            // Use steapestDescent to make sure Gradient is correct
            interp.printOptOutput = true;
            interp.sdMaxDisp = 0.2;
            interp.lineAlg = SteapestDescent;

            // Optimize
            converged = interp.interpolate();
            std::cout << "Interpolate Converged = " << converged << std::endl;

            // Print out data
            std::cout << "Coeffs:" << std::endl;
            std::cout << std::fixed << std::setprecision(8) << interp.coeffs[0] << std::endl;
            std::cout << std::fixed << std::setprecision(8) << interp.coeffs[1] << std::endl << std::endl;
            for( size_t i=0; i<nInter; i++)
                coeffs[i] = interp.coeffs[i];

            gOpt = CQMemManager::get().malloc<double>(2);
            interp.gradFunc();
            std::copy_n(interp.gPointer(),2,gOpt);
            std::cout << "Gradient:" << std::endl;
            std::cout << std::scientific << std::setprecision(8) << gOpt[0] << std::endl;
            std::cout << std::scientific << std::setprecision(8) << gOpt[1] << std::endl;
            std::cout << std::endl << std::endl;
        };
        InterpolateTest(const BFGSTest &) = delete;
        InterpolateTest(BFGSTest &&) = delete;

        ~InterpolateTest()
        {
            CQMemManager::get().free(B,gOpt);
        };
    };
}