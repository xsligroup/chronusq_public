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

using cart_t = std::array<double,3>;

namespace ChronusQ {

  /**
   *  \brief The Quadrature Struct.
   *
   *  Contains information pertinant to a quadrature rule
   *  for numerical integration.
   */
  template <typename _PtTyp>
  struct Quadrature {

    typedef _PtTyp point_type; ///< Storage type for the point

    size_t              nPts;    ///< Number of integration points
    std::vector<_PtTyp> pts;     ///< Integration points
    std::vector<double> weights; ///< Quadrature weights

    double upBound;  ///< Upper bound for the numerical integration
    double lowBound; ///< Lower bound for the numerical integration


    // Deleted ctors
    Quadrature() = delete;


    Quadrature(size_t N, double a, double b):
      upBound(b), lowBound(a), nPts(N) {}; // Quadrature


    Quadrature(const Quadrature &other) = default;
    Quadrature(Quadrature &&other)      = default;


    /**
     *  \brief Generic scheme for reseting and allocating
     *  space for the quadrature storage. 
     *
     *  To be specificed in all derived functions.
     */  
    virtual void generateQuadrature() {
      weights.clear(); if(nPts) weights.resize(nPts);
      pts.clear(); if(nPts) pts.resize(nPts);
    };

    void printQuadrature(){
      std::cerr << std::endl;
      std::cerr << "Start Printing Grid " << std::endl;
      std::cerr << std::endl;
      for(auto i = 0; i < this->nPts; i++) 
        std::cerr << i << " "<< this->pts[i] << ", " << this->weights[i] << std::endl;
      double sumW = 0.;
      for(auto i = 0; i < this->nPts; i++) 
        sumW += weights[i];
      std::cerr << "SUM W " << sumW <<std::endl;
      std::cerr << std::endl;
      std::cerr << "End Printing Grid " << std::endl;
      std::cerr << std::endl;
    };

  }; // struct Quadrature

  // Template for Quadrature specialization
  #define QuadratureTemplate(NAME,TYP) \
  struct NAME : public Quadrature<TYP> {\
  \
    NAME(size_t N, \
         double a = -std::numeric_limits<double>::infinity(),\
         double b = std::numeric_limits<double>::infinity()) : \
      Quadrature<TYP>(N,a,b) {}; \
\
    NAME(const NAME &) = default;\
    NAME(NAME &&)      = default;\
\
    void generateQuadrature();  \
\
  };


  // Declarations of quadrature specializations (see src/grid for details)
  QuadratureTemplate(GaussChebFst,double);
  QuadratureTemplate(GaussChebSnd,double);
  QuadratureTemplate(GaussLegendre,double);
  QuadratureTemplate(EulerMac,double);
  QuadratureTemplate(Lebedev,cart_t);


  // Utility Functions
  std::pair<double,double> unitBoundTransform(double lowBound, double upBound, double pt, double wgt);

}; // namespace ChronusQ

