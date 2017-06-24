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

#include <grid/quadrature.hpp>

namespace ChronusQ {


  void GaussLegendre::generateQuadrature() {

    Quadrature::generateQuadrature();

    //std::cerr << "GaussLegendre" << std::endl;
    //std::cerr << "Limits " << this->lowBound << " " << this->upBound <<std::endl;
    
    assert(upBound > lowBound);
    this->nPts = this->nPts + this->nPts % 2;
    int m = (this->nPts + 1) / 2;
    for(size_t i = 1; i <= m; i++){

      // Legendre indexing starts and 1
      double z = std::cos(M_PI*(double(i) - 0.25)/(double(this->nPts) + 0.5));
      double pp(0.);
      double z1(0.);
      double eps(3.e-11);
  
      // Iteratively determine the i-th root 
      while(std::abs(z-z1) > eps) {
        double p1(1.), p2(0.);
  
        // Loop over the recurrence relation to evaluate the
        // Legendre polynomial at position z
        for(int j = 1; j <= this->nPts; j++){
          double p3 = p2;
          p2 = p1; 
          p1 = ( (2. * double(j) - 1.) * z * p2 - (double(j) - 1.) * p3) / 
               double(j);
        } // end j for
  
        // p1 is now the desired Legrendre polynomial. We next compute 
        // pp, its derivative, by a standard relation involving also p2,
        // the polynomial of one lower order
        pp = double(this->nPts) * (z * p1 - p2) / (z * z - 1.);
        z1 = z;
        z  = z1 - p1 / pp;
      } // end while
  
      double pt = z;
      double wgt = 2. / (1. - z*z) / pp / pp;

      // Transform points and populate arrays
      std::tie(pts[i-1],weights[i-1]) = unitBoundTransform(lowBound,upBound,pt,wgt);
      //if(this->nPts % 2 == 0 or i != m){
        double reflectedPoint = pts[i-1];
        reflectedPoint -= (this->upBound + this->lowBound);
        reflectedPoint = -reflectedPoint;
        pts[this->nPts-i] = reflectedPoint; 
        weights[this->nPts-i] = weights[i-1]; 
      //}

    }; // Loop over points
    this->printQuadrature();
  };


};
