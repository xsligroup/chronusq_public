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

//#define INTERPOLATE_DEBUG
//#define INTERPOLATE_DEBUG_GRADIENT
//#define INTERPOLATE_PRINT_COEFFS
//#define INTERPOLATE_PRINT_WARN
//#define INTERPOLATE_USE_GAUSSIAN_ALG

#include <nonlinearopt/bfgs.hpp>
#include <nonlinearopt/goldenSection.hpp>

namespace ChronusQ
{

    /**
   *   \brief The Energy DIIS class. A class to perform a DIIS interpolation 
   *    based on a series of energies and couplings.
   * 
   *    In general it solves the equation energy*c - 1/2 c^T*B*c where c are a set of
   *    coefficients between 0. and 1.  This object uses a set of optimization
   *    variable that maintain this constraint, which allows for unconstrained
   *    optimization of the function.  
   *
   *    This object inherits the BFGS optimizer to find the  
   *    minimum for the coefficients.
   */

    template <typename T>
    class ENERGYDIIS : public BFGS<T>
    {

    protected:
        // Useful typedefs
        typedef T *oper_t;
        typedef std::vector<oper_t> oper_t_coll;
        typedef std::vector<oper_t_coll> oper_t_coll2;
        T *xOpt; ///< Unconstrained optimization variables
        T *gOpt; ///< gradient of optimization variables

    public:
        // Input parameters
        size_t nInter;         ///< Size of extrapolation space
        std::vector<T> energy; ///< Energies of each SCF interation
        double *B;             ///< Error metric
        size_t nDimB;          ///< The leading dimension of B
        size_t target;         ///< Integer for the target state (guess)

        // Output parameters
        std::vector<T> coeffs; ///< Vector of extrapolation coeficients

        // Constructor

        /*
        *  DIIS Constructor. Constructs a DIIS object
        */
        ENERGYDIIS(size_t nInter, double *B, size_t nD, std::vector<T> en, size_t t, double conv) : 
            nInter(nInter), B(B), nDimB(nD), energy(en), target(t), BFGS<T>(nInter, conv)
        {
            // Allocate data structures
            coeffs.resize(nInter);
            for( size_t i=0; i<nInter; i++ )
                coeffs[i] = 0.;
            xOpt = CQMemManager::get().malloc<T>(nInter);
            gOpt = CQMemManager::get().malloc<T>(nInter);
            BFGS<T> :: setXPointer( xOpt );
            BFGS<T> :: setGPointer( gOpt );
            double eMin = *min_element(energy.begin(),energy.end());
            for( auto &e : energy)
                e -= eMin;
        };

        // Constructors for default, copy, and move
        ENERGYDIIS() = delete;
        ENERGYDIIS(const ENERGYDIIS &) = delete;
        ENERGYDIIS(ENERGYDIIS &&) = delete;

        // Destructor
        ~ENERGYDIIS()
        {
            CQMemManager::get().free(xOpt, gOpt);
        };

        // Getter functions for pointers
        inline T* xPointer() { return xOpt; }
        inline T* gPointer() { return gOpt; }

        // Public Member functions
        bool interpolate();

        // Optimization functions
        void computeCoeffs();
        double objFunc();
        void gradFunc();

    }; // class DIIS


    /*
    *   Brief: The variables are optimized using a set of unconstrained
    *          variables. coeffs[i] = xOpt[i]*xOpt[i] / (\sum_k xOpt[k]*xOpt[k])
    */ 
    template <typename T>
    void ENERGYDIIS<T>::computeCoeffs()
    {
        T sumC = T(0.);
        for (size_t i = 0; i < nInter; i++)
        {
            coeffs[i] = xOpt[i] * xOpt[i];
            sumC += coeffs[i];
        }

        for (size_t i = 0; i < nInter; i++)
            coeffs[i] /= sumC;
    }


    /*
    *   Brief: This function computes the equation we are optimizing
    *          val = E^T*coeffs - 1/2 coeffs^T B coeffs
    */
    template <typename T>
    double ENERGYDIIS<T>::objFunc()
    {
        computeCoeffs();

        double val = T(0.);
        for (size_t i = 0; i < nInter; i++)
        {
            val += coeffs[i] * energy[i];
            for (size_t j = 0; j < nInter; j++)
                val -= T(0.5) * coeffs[i] * B[j + i * nDimB] * coeffs[j];
        }

        return val;
    }

    /*
    *   Brief: This is the gradient for the objective function with respect
    *          to the unconstrained optimization variables(xOpt).
    * 
    */
    template <typename T>
    void ENERGYDIIS<T>::gradFunc()
    {
        computeCoeffs();

        // Renormalize parameters
        // This helps keep the gradient consistent during the optimization
        for (size_t i = 0; i < nInter; i++)
            this->xOpt[i] = std::sqrt(std::abs(coeffs[i]));

        T *lag = CQMemManager::get().malloc<T>(nInter);
        T sumC = T(0.);
        for (size_t i = 0; i < nInter; i++)
        {
            sumC += this->xOpt[i] * this->xOpt[i];
            lag[i] = energy[i];
            for (size_t j = 0; j < nInter; j++)
                lag[i] -= B[j + i * nDimB] * coeffs[j];
        }

        for (size_t k = 0; k < nInter; k++)
        {
            T cpartial;
            this->gOpt[k] = T(0.);
            for (size_t i = 0; i < nInter; i++)
            {
                cpartial = T(0.);
                if (i == k)
                {
                    cpartial = T(2.) * this->xOpt[i] / sumC * (T(1.) - coeffs[i]);
                }
                else
                {
                    cpartial = T(-2.) * coeffs[i] * this->xOpt[k] / sumC;
                }
                this->gOpt[k] += cpartial * lag[i];
            }
        }

#ifdef INTERPOLATE_DEBUG_GRADIENT
        this->printGrad = true;
        std::cout << "Coeffs:" << std::endl;
        for (size_t i = 0; i < nInter; i++)
            std::cout << coeffs[i] << std::endl;
        std::cout << std::endl;
#endif

        CQMemManager::get().free(lag);
    }


    /*
    *   Brief: This routine drives the optimization and returns the result
    */

    template <typename T>
    bool ENERGYDIIS<T>::interpolate()
    {
#ifdef INTERPOLATE_DEBUG
        std::cout << "Energies:" << std::endl;
        for (size_t i = 0; i < nInter; i++)
            std::cout << energy[i] << std::endl;
        std::cout << std::endl;
        prettyPrintSmart(std::cout, "B Matrix", B, nInter, nInter, nDimB);
        std::cout << std::endl;
        this->printOptOutput = true;
        this->printSearchDirection = true;
        this->printLineSearch = true;
        //this->lineAlg = SteapestDescent;  // This gives a good test of the gradient
        //this->sdMaxDisp = 0.05;
        this->updateInvHess = true;
#endif

        // Make the target guess twice as large as the rest unless nInter == 2
        T targetGuess = nInter == 2 ? T(1./double(nInter)) : T( 2./double(nInter) );
        T restGuess = std::sqrt((1.-targetGuess) / T(nInter-1)); 
        targetGuess = std::sqrt(targetGuess);
        for (size_t i = 0; i < nInter; i++)
            this->xOpt[i] = T(restGuess);
        this->xOpt[target] = targetGuess;
        double guessE = this->objFunc();

#ifdef INTERPOLATE_PRINT_COEFFS
        std::cout << "EDIIS Guess Coefficients:  " <<  guessE << std::endl;
        for (size_t i = 0; i < nInter; i++)
        {
            std::cout << coeffs[i];
            if (i == target)
                std::cout << "***";
            std::cout << std::endl;
        }
        std::cout << std::endl;
#endif

        bool conv = BFGS<T> :: optimize();
        computeCoeffs();

#ifdef INTERPOLATE_PRINT_COEFFS
        std::cout << "EDIIS Coefficients:" << std::endl;
        for (size_t i = 0; i < nInter; i++)
        {
            std::cout << coeffs[i];
            if (i == target)
                std::cout << "***";
            std::cout << std::endl;
        }
        std::cout << std::endl;
#endif

        // To force progress, make sure target coeff is greater than 1E-8
        if (coeffs[target] < 1E-4)
        {
#ifdef INTERPOLATE_PRINT_WARN
            std::cout << "EDIIS Warning: found solution with latest coefficient zero." << std::endl;
#endif

#ifndef INTERPOLATE_USE_GAUSSIAN_ALG
            this->xOpt[target] = std::sqrt(1E-4);
            this->objFunc();
#endif

            // NOTE: The following is an implementation of what Kudin et al. describe to do
            //       when the target is zero. However, I found that making the target 1E-4 
            //       does not change the coeffs a lot and they are still close to the optimal
            //       ones. 
            //
            // Mix the other degrees of freedom with the targe state to find local minimum
            // This is Gaussian's method and provided in:
            // Kudin, K. N.; Scuseria, G. E.; Cancès, E. A Black-Box Self-Consistent Field Convergence Algorithm: One Step Closer. J. Chem. Phys. 2002, 116 (19), 8255–8261.

#ifdef INTERPOLATE_USE_GAUSSIAN_ALG
            // Start with target as unity
            std::fill_n(this->xOpt, nInter, T(0.));
            this->xOpt[target] = 1.;

            std::vector<size_t> vopt = {target};
            std::vector<size_t> vmix;
            for (size_t i = 0; i < nInter; i++)
                if (i != target)
                    vmix.push_back(i);

            T *dx = CQMemManager::get().malloc<T>(nInter);
            T *initX = CQMemManager::get().malloc<T>(nInter);
            double ymin;
            T amin;
            for (auto &iMix : vmix)
            {
                // Compute dx
                std::fill_n(dx, nInter, T(0.));
                for (auto &iOpt : vopt)
                    dx[iOpt] = T(-1.);
                dx[iMix] = T(1.);

                // Perform Line Search
                std::copy_n(xOpt, nInter, initX);
                GoldenSectionSearch<T>::optimize(initX, dx, T(1E-8), T(0.9), ymin, amin);

                vopt.push_back(iMix);
            }
            CQMemManager::get().free(dx, initX);
#endif

#ifdef INTERPOLATE_PRINT_COEFFS
            std::cout << "EDIIS Coefficients:" << std::endl;
            for (size_t i = 0; i < nInter; i++)
            {
                std::cout << coeffs[i];
                if (i == target)
                    std::cout << "***";
                std::cout << std::endl;
            }
            std::cout << std::endl;
#endif
        }
        return conv;
    }

}; // namespace ChronusQ
