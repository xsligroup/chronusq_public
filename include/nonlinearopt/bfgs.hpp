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
#include <nonlinearopt/goldenSection.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>

namespace ChronusQ
{

    enum LineSearchAlgorithm
    {
        GoldenSection,
        SteapestDescent,
        ParticleSwarm
    };

    /*
    *   Brief: This class computes the BFGS quasi-Newton step from the positions,
    *           and gradients.  The inverse Hession is updated in place for each
    *           step.  At the end of the optimization, xBFGS will contain the new set of
    *           positions. 
    * 
    *   NOTE: The objFunc and gradFunc are called at each step and should be overloaded
    *         by the inheriting object
    * 
    *   N       = number of positions/gradients i.e. dimension of the problem
    *   xBFGS   = pointer to the current set of positions (should be set with setXPointer)
    *   gBFGS   = pointer to the gradient at the current position (should be set with setGPointer)
    */
    template <typename VT>
    class BFGS : public GoldenSectionSearch<VT>
    {
    protected:
        // Function to be Optimized declarations
        size_t N; ///< Dimension of optimization problem
        VT *xBFGS; ///< Pointer to variables
        VT *gBFGS; ///< Pointer to Gradient vector

    public:
        // Optimization parameters
        double convergence = 1E-12;                  ///< Convergence criteria for objective function
        bool updateInvHess = true;                   ///< Whether to update the Inverse Hessian with BFGS update
        bool printSearchDirection = false;           ///< Whether to print the search direction in BFGS opt
        size_t printPrecision = 10;                  ///< Level of Precision for printing
        bool resetInvHess = true;                    ///< Whether to reset the invHess to unit matrix every N steps
        bool scaleInvHess = true;                    ///< Whether to use to rescale the invHess during the optimization
        LineSearchAlgorithm lineAlg = GoldenSection; ///< Choose algorithm to perform line Search.
                                                     ///< NOTE: GoldenSectionSearch requires many calls to the objFunc
                                                     ///<       but zero calls to gradFunc
        double goldenSectionMaxDisp = 5.;            ///< The Maximum displacement for the GoldenSectionSearch

        double sdMaxDisp = 1.0;     ///< The maximum displacement for SteapestDescent along the search direction
        double scaleFailedSD = 0.1; ///< Amount to scale the SteapestDescent step if the objFunc increases by too much
        double percentFailSD = 10.; ///< The percent amount the function can go up by without reducing the search step

        // Constructor
        BFGS() = delete;
        BFGS(size_t N, double conv) : N(N), convergence(conv), GoldenSectionSearch<VT>(N, 0.01 * conv){

                                                                                                 };

        // Remove copy/move constructors
        BFGS(const BFGS<VT> &) = delete;
        BFGS(BFGS<VT> &&) = delete;

        // Destructor:
        ~BFGS(){};

        // Getter and Setter Pointer Functions
        inline VT *gPointer() { return gBFGS; };
        inline void setGPointer(VT *g) { gBFGS = g; };
        inline VT *xPointer() { return xBFGS; };
        inline void setXPointer(VT *xNew)
        {
            xBFGS = xNew;
            if (lineAlg == GoldenSection)
                GoldenSectionSearch<VT>::setXPointer(xNew);
        };

        // Optimization functions
        bool optimize();                                            ///< Will return true for a successful optimization and false for unsucessful
        void lineSearch(VT *initX, VT *dx, double &ymin, VT &amin); ///< Line search algorithms

        // These functions should be overloaded by inheriting object
        double objFunc() { return 1.; };
        void gradFunc(){};

        // Functions used during optimization
        inline void makeUnit(VT *invHess)
        {
            for (size_t i = 0; i < N * N; i++)
                invHess[i] = VT(0.);
            for (size_t i = 0; i < N; i++)
                invHess[i + i * N] = VT(1.);
        };
    };

    template <typename VT>
    bool BFGS<VT>::optimize()
    {
        VT *invHess = CQMemManager::get().malloc<VT>(N * N);
        VT *gPrev = CQMemManager::get().malloc<VT>(N);
        VT *dx = CQMemManager::get().malloc<VT>(N);
        VT *dg = CQMemManager::get().malloc<VT>(N);
        VT *SCR = CQMemManager::get().malloc<VT>(N * N);
        VT *U = CQMemManager::get().malloc<VT>(N * N);
        VT *hessUpdate = CQMemManager::get().malloc<VT>(N * N);
        VT alpha = VT(0.);
        VT yPrev = VT(0.);
        double ymin;
        VT amin;

        bool result = false;

        // Compute initial/guess values before loop
        this->gradFunc();
        ymin = this->objFunc();

        if (this->printVariables)
        {
            std::cout << "BFGS: Initial Parameters:" << std::endl;
            std::cout << std::fixed << std::setprecision(printPrecision);
            for (size_t i = 0; i < N; i++)
                std::cout << xBFGS[i] << std::endl;
            std::cout << std::endl;
        }

        if (this->printGrad)
        {
            std::cout << "BFGS: Initial Gradient:" << std::endl;
            std::cout << std::fixed << std::setprecision(printPrecision);
            for (size_t i = 0; i < N; i++)
                std::cout << gBFGS[i] << std::endl;
            std::cout << std::endl;
        }

        if (this->printOptOutput)
        {
            //std::cout << std::endl;
            std::cout << std::scientific << std::setprecision(printPrecision) << "BFGS: Step 0:  Function = " << ymin << std::endl;
            //std::cout << std::endl;
            yPrev = ymin;
        }

        // BFGS optimization loop
        for (size_t iStep = 1; iStep <= this->nStep; iStep++)
        {
            // Compute inverse Hessian
            if (iStep == 1 or ((iStep + 1) % N == 0 and resetInvHess))
            {
                makeUnit(invHess);
            }
            else if (updateInvHess)
            {
                // Compute dx dg inner product
                VT dxdg = blas::dot(N, dx, 1, dg, 1);

                // Using a scaling for invHessian. See following for details:
                // Luenberger, D. G.; Ye, Y. Linear and Nonlinear Programming; Hillier, F. S., Ed.; Springer Science+Business Media.
                VT scale = VT(1.);
                if (scaleInvHess)
                {
                    VT dgHdg = VT(0.);
                    for (size_t j = 0; j < N; j++)
                        for (size_t i = 0; i < N; i++)
                            dgHdg += dg[i] * invHess[i + j * N] * dg[j];
                    scale = dxdg / dgHdg;
                    if (std::isnan(scale))
                        scale = VT(1.);
                }

                // Compute inverse Hessian transformation ( I-( dg*dx^T / dxdg) )
                makeUnit(U);

                for (size_t j = 0; j < N; j++)
                    for (size_t i = 0; i < N; i++)
                        U[i + j * N] -= dg[i] * dx[j] / dxdg;

                // Transform invHess  =>  U^C * invHess * U
                blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans, blas::Op::NoTrans, N, N, N, VT(scale), U, N, invHess, N, VT(0.), SCR, N);
                blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, N, N, N, VT(1.), SCR, N, U, N, VT(0.), invHess, N);

                // add (dx*dx^T/dotProd) to invHess
                for (size_t j = 0; j < N; j++)
                    for (size_t i = 0; i < N; i++)
                        invHess[i + j * N] += dx[i] * dx[j] / dxdg;

                // Check invHess for nan
                for (size_t i = 0; i < N * N; i++)
                    if (std::isnan(invHess[i]))
                    {
                        std::cout << "WARNING: Found a NaN in BFGS inverse Hessian" << std::endl;
                        makeUnit(invHess);
                        break;
                    }
            }

            // Compute search direction
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, N, 1, N, VT(-1.), invHess, N, gBFGS, N, VT(0.), dx, N);

            if (printSearchDirection)
            {
                std::cout << "DX:" << std::endl;
                for (size_t i = 0; i < N; i++)
                    std::cout << dx[i] << std::endl;
                std::cout << std::endl;
            }

            // Perform Line Search
            std::copy_n(xBFGS, N, SCR);
            BFGS<VT>::lineSearch(SCR, dx, ymin, amin);

            for (size_t i = 0; i < N; i++)
                dx[i] = amin * dx[i];

            // Compute new Gradient and positions
            for (size_t i = 0; i < N; i++)
                xBFGS[i] = SCR[i] + dx[i];

            ymin = this->objFunc();
            std::copy_n(gBFGS, N, gPrev);
            this->gradFunc(); // This overwrites gBFGS
            for (size_t i = 0; i < N; i++)
                dg[i] = gBFGS[i] - gPrev[i];

            // Print output
            if (this->printVariables)
            {
                std::cout << "BFGS: Parameters:" << std::endl;
                std::cout << std::scientific << std::setprecision(printPrecision);
                for (size_t i = 0; i < N; i++)
                    std::cout << xBFGS[i] << std::endl;
                std::cout << std::endl;
            }

            if (this->printGrad)
            {
                std::cout << "BFGS: Gradient:" << std::endl;
                std::cout << std::scientific << std::setprecision(printPrecision);
                for (size_t i = 0; i < N; i++)
                    std::cout << gBFGS[i] << std::endl;
                std::cout << std::endl;
            }

            if (this->printOptOutput)
            {
                //std::cout << std::endl;
                std::cout << std::scientific << std::setprecision(printPrecision) << "BFGS: Step " << iStep << ":  Function = " << ymin << "   Change = " << ymin - yPrev << std::endl;
                //std::cout << std::endl;
            }

            // Evaluate Convergence
            if (std::abs(yPrev - ymin) < convergence)
            {
                result = true;
                break;
            }
            yPrev = ymin;
        }

        CQMemManager::get().free(invHess, gPrev, dx, dg, SCR, hessUpdate, U);
        return result;
    };

    template <typename VT>
    void BFGS<VT>::lineSearch(VT *initX, VT *dx, double &ymin, VT &amin)
    {
        double yPrev = ymin;
        switch (lineAlg)
        {
        case GoldenSection:
        {
            GoldenSectionSearch<VT>::optimize(initX, dx, VT(convergence), VT(goldenSectionMaxDisp), ymin, amin);
            // if step size is too large then decrease
            // the step size. Here we are using percent
            // difference to determine a large increase
            if (100. * (ymin - yPrev) / std::abs(yPrev) > percentFailSD)
            {
                GoldenSectionSearch<VT>::optimize(initX, dx, VT(convergence), VT(scaleFailedSD * goldenSectionMaxDisp), ymin, amin);
            }

            break;
        }
        case SteapestDescent:
        {
            amin = sdMaxDisp;
            for (size_t i = 0; i < N; i++)
                xBFGS[i] = initX[i] + amin * dx[i];
            ymin = this->objFunc();

            // if step size is too large then decrease
            // the step size. Here we are using percent
            // difference to determine a large increase
            if (100. * (ymin - yPrev) / std::abs(yPrev) > percentFailSD)
            {
                amin *= scaleFailedSD;
                for (size_t i = 0; i < N; i++)
                    xBFGS[i] = initX[i] + amin * dx[i];
                ymin = this->objFunc();
            }
            break;
        }
        case ParticleSwarm:
        {
            CErr("Particle Swarmm Line Search has not been implemented yet");
        }
        }
    };
};
