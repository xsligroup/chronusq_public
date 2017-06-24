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
#include <nonlinearopt/nonlinearopt.hpp>
#include <cqlinalg/blas3.hpp>

namespace ChronusQ
{

    /*
    *   Brief: This object performs a line search using the Golden
    *          section search algorithm
    * 
    *   Luenberger, D. G.; Ye, Y. Linear and Nonlinear Programming; Hillier, F. S., Ed.; Springer Science+Business Media.
    * 
    *   N = number of positions i.e. dimension of the problem
    *   x = pointer to the current set of positions (should be set using setXPointer)
    */
    template <typename VT>
    class GoldenSectionSearch : public NonLinearOpt
    {
    protected:
        // Function to be Optimized declarations
        size_t N; ///< Dimension of optimization problem
        VT *x;                        ///< Pointer to variables

    public:
        // Optimization parameters
        double convergence = 1E-12;   ///< Convergence criteria for objective function
        bool printLineSearch = false; ///< Whether to print line search results to output
        size_t printPrecision = 10;   ///< Precision for all printing

        // Constructor
        GoldenSectionSearch() = delete;
        GoldenSectionSearch(size_t N, double conv) : N(N), convergence(conv){};

        // Remove copy/move constructors
        GoldenSectionSearch(const GoldenSectionSearch<VT> &) = delete;
        GoldenSectionSearch(GoldenSectionSearch<VT> &&) = delete;

        // Destructor:
        ~GoldenSectionSearch(){};

        // Getter and Setter Functions
        inline VT* xPointer(){return x;};
        inline void setXPointer( VT* xNew) { x = xNew;};

        // Optimization functions
        bool optimize() { return false; };                        ///< Dummy routine for nonlinear interface. This is not used
        void optimize(VT *initX, VT *dx, VT searchMin, VT searchMax, double &ymin, VT &amin); ///< Perform line search on interval

        // These functions should be overloaded by inheriting object
        double objFunc() { return 1.; };
        void gradFunc(){};
    };


    /*
    *   Brief: This function performs a GoldenSection search along the search direction(dx) with the initial set of 
    *          positions at initX.  SearchMin and searchMax determine the interval in which to search.  The results
    *          are the minimum of the objective function (ymin) and the displacement along the search direction (amin).
    */
    template <typename VT>
    void GoldenSectionSearch<VT>::optimize(VT *initX, VT *dx, VT searchMin, VT searchMax, double &ymin, VT &amin)
    {

        ymin = VT(1.E10);
        amin = VT(0.);
        VT a = searchMin;
        VT b = searchMax;
        VT gr = VT((1. + std::sqrt(5.)) / 2.);
        VT c = b - (b - a) / gr;
        VT d = a + (b - a) / gr;
        size_t xSearch = 0;
        std::copy_n(initX,N,x);
        double yPrev = this->objFunc();
        for (size_t iSearch = 0; iSearch < this->nStep; iSearch++)
        {
            for (size_t i = 0; i < N; i++)
                x[i] = initX[i] + c * dx[i];
            double yc = this->objFunc();

            for (size_t i = 0; i < N; i++)
                x[i] = initX[i] + d * dx[i];
            double yd = this->objFunc();

            if( printLineSearch )
                std::cout << std::scientific << std::setprecision(printPrecision) << "Golden Section Search Step " << iSearch+1 << ":  Function value = " << (yc+yd)/VT(2.) << "   Change = " << std::abs(yc-yd) << std::endl;

            // Check convergence
            bool converged = (std::abs(yc-yd)<convergence) and (std::abs(yc-yPrev)<convergence) and (std::abs(yd-yPrev)<convergence);
            if ( converged )
            {
                xSearch = iSearch;
                amin = (a + b) / VT(2.);
                ymin = (yc + yd) / 2.;
                break;
            }
            yPrev = (yc + yd) / 2.;

            // update variables for next iteration
            if (yc < yd)
            {
                b = d;
            }
            else
            {
                a = c;
            }
            c = b - (b - a) / gr;
            d = a + (b - a) / gr;
        }

        if (printLineSearch)
        {
            std::cout << std::endl;
            std::cout << "Number Golden Section Search Steps = " << xSearch << std::endl;
            std::cout << std::scientific << std::setprecision(printPrecision)  << "Optimized Displacement  = " << amin << "    Function Minimum = " << ymin << std::endl;
        }
    }
};