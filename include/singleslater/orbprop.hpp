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

#include <singleslater.hpp>
#include <grid/integrator.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeOrbitalProps() {

    ROOT_ONLY(comm);

    std::cout << std::endl;
    std::cout << "Orbital Properties" << std::endl;
    std::cout << bannerTop << std::endl << std::endl;

    if(doOrbEne) printOrbitalEnergies();

    if(doRDFs) computeOrbitalRDFs();

  }; // computeOrbitalProps

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printOrbitalEnergies() {

    std::cout << std::endl;
    std::cout << "Orbital Energies" << std::endl;
    std::cout << bannerTop << std::endl << std::endl;

    // loop over MOs
    size_t NB = this->nAlphaOrbital();
    size_t NOrb = NB*this->nC;
    size_t negEShift = (this->nC == 4) ? (2*NB) : 0; //Skip past neg-E solutions when 4c
    std::cout << "Starting at orbital: " << negEShift+1 << std::endl;
    std::cout << std::setprecision(8) << std::fixed;

    // mo[0] orbital energies
    for( size_t iOrb = negEShift; iOrb < NOrb; iOrb++)
      std::cout << "MO #" << iOrb+1 << ": " << this->eps1[iOrb] << std::endl;

    // mo[1] orbital energies
    if( this->nC == 1 and not this->iCS){
      std::cout << std::endl;
      std::cout << "Beta Orbital Energies" << std::endl;
      std::cout << bannerTop << std::endl << std::endl;
      for( size_t iOrb = negEShift; iOrb < NOrb; iOrb++)
        std::cout << "MO #" << iOrb+1 << ": " << this->eps2[iOrb] << std::endl;
    }

  }; // printOrbitalEnergies

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeOrbitalRDFs() {

    std::cout << std::endl;
    std::cout << "Orbital Radial Distribution Functions" << std::endl;
    std::cout << bannerTop << std::endl << std::endl;

    //TODO: implement USCF

    const double normVal = 4.0*M_PI; // Lebedev weights normalized to 1.0, not 4pi

    // Logarithmic grid
    double radSpacing = exp((log(finalRad) - log(initialRad)) / (numRadPts-1));
    Integrator1D<Lebedev> integrator(Lebedev{numAngPts});
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "# of Lebedev angular points: " << integrator.getnPts() << std::endl;
    std::cout << "Radial starting point of exponential grid (in a.u.), final point, and # radial points: "
        << initialRad << ", " << finalRad << ", and " << numRadPts << std::endl;

    // loop over MOs
    size_t NB = this->nAlphaOrbital();
    size_t NOrb = NB*this->nC;
    // TODO: beta orbital in 1c
    MatsT * locMO = this->mo[0].pointer();
    size_t negEShift = (this->nC == 4) ? (2*NB) : 0; //Skip past neg-E solutions when 4c
    std::cout << "Starting at orbital: " << negEShift+1 << std::endl;
    size_t upperOrb = (numMOs == 0) ? NOrb : negEShift+numMOs; // User can define number of MOs

    for( size_t iOrb = negEShift; iOrb < upperOrb; iOrb++){
      std::vector<double> BASIS(numRadPts*NB);
      MatsT * thisMO = locMO + (iOrb) * NOrb; // Doing an MO at a time

      // Radial expectation values are computed with trapezoidal rule
      std::vector<double> radialPt(numRadPts); // The radius
      std::vector<double> radialDistFnc(numRadPts); // Radial distribution function
      double radExpectNum = 0.0; // Numerical radial expectation value numerator
      double radExpectDenom = 0.0; // Numerical radial expectation value denominator
      double radExpectValue = 0.0; // Numerical radial expectation value

      for (auto iRad = 0l; iRad < numRadPts; iRad++) {

        radialPt[iRad] = initialRad * pow(radSpacing,iRad);

        // Scales Lebedev point by rad and calculates MO
        auto integrateFunction = [&, this](const std::array<double, 3>& point) -> MatsT {

          // Scale point on unit sphere
          cart_t frad = {radialPt[iRad]*point[0],radialPt[iRad]*point[1],radialPt[iRad]*point[2]};

          // Evaluate basis function on grid
          evalShellSet(NOGRAD,this->basisSet().shells,&frad[0],1,&BASIS[iRad*NB],false);

          // Dot product of basis function and coefficient
          MatsT val_alpha = blas::dot(NB,thisMO,1,&BASIS[iRad*NB],1);
          MatsT val_beta = MatsT(0.0);
          size_t skipToBeta = (this->nC == 4) ? (2*NB) : NB; // Skip to beta (large)
          if( this->nC>1 ) val_beta = blas::dot(NB,thisMO+skipToBeta,1,&BASIS[iRad*NB],1);

          // For 4c, this is just large component contribution
          return std::norm(val_alpha) + std::norm(val_beta); // Return wave function squared (will be zero for p orbital if not)

        };

        // Integrate the contraction (MO coefficient and basis function) over the grid
        MatsT result = integrator.integrate<MatsT>(integrateFunction);

        radialDistFnc[iRad] = radialPt[iRad]*radialPt[iRad]*std::real(result)*normVal; // r^2 |R(r)|^2 4pi

        // Trapezoidal rule for numerical radial expectation value
        if( iRad > 0 ){
          radExpectNum += 0.5*(radialPt[iRad-1]*radialDistFnc[iRad-1]+radialPt[iRad]*radialDistFnc[iRad])*(radialPt[iRad]-radialPt[iRad-1]);
          radExpectDenom += 0.5*(radialDistFnc[iRad-1]+radialDistFnc[iRad])*(radialPt[iRad]-radialPt[iRad-1]);
        }

        // Print the result of the integration for each radial point
        std::cout << "MO #" << iOrb + 1
                  << " radial distribution function for point " << radialPt[iRad] << ": " << radialDistFnc[iRad] << std::endl;

      } // loop over radial points

      constexpr double radialStabilityThresh = 1e-7;
      if( radExpectNum>radialStabilityThresh and radExpectDenom>radialStabilityThresh ) radExpectValue = radExpectNum/radExpectDenom;
      std::cout << "MO #" << iOrb + 1
        << " numerical radial expectation value = " << radExpectNum << "/" << radExpectDenom << " = " << radExpectValue << std::endl;


      std::cout << std::endl;

    } // loop over MOs

  }; // computeOrbitalRDF

}; // namespace ChronusQ

