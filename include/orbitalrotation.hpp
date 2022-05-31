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
#include <memmanager.hpp>
#include <mcwavefunction.hpp>
#include <matrix/squarematrix.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

// Orbital Rotation Headers

//#define DEBUG_OrbitalRotation

namespace ChronusQ {
  
  enum OrbitalRotationAlgorithm {
    ORB_ROT_APPROX_QUASI_2ND_ORDER,
    ORB_ROT_QUASI_2ND_ORDER,
    ORB_ROT_2ND_ORDER,
  };
  
  /* 
   * Brief Definition of OrbitalRotationSettings
   */
  
  struct OrbitalRotationSettings {
    
    // blocks for rotation
    bool rotate_within_correlated   = true;
    bool rotate_inact_correlated    = true;
    bool rotate_correlated_virtual  = true;
    bool rotate_inact_virtual       = true;
    bool rotate_negative_positive   = false;

    OrbitalRotationAlgorithm alg = ORB_ROT_APPROX_QUASI_2ND_ORDER;

    // handle hessians 
    double hessianDiagScale        = 2.0;
    double hessianDiagDampTol      = 20.0;
    double hessianDiagDamp         = 10.0;
    double hessianDiagMinTol       = 1.0e-3;
     
    double XDampTol  = 0.5;
    
    OrbitalRotationSettings() = default;
    OrbitalRotationSettings(const OrbitalRotationSettings &) = default;
    OrbitalRotationSettings(OrbitalRotationSettings      &&) = default;
  
    void print(bool fourComp);
  };
  
  /* 
   * \brief the OrbitalRotation class. The class can perform 
   * post-SCF orbital rotation using Newton-Ralphson method
   *
   * \warning:
   *    make sure the dimension of input are correct
   *    make sure the input RDMs are in same definition.
   *
   * Exponential parametrization of the MO coefficient:
   * |\Psi>_new = exp(X) |\Psi>_old 
   *
   */
  template <typename MatsT, typename IntsT>
  class OrbitalRotation { 
    
  protected:
    
    MCWaveFunction<MatsT,IntsT> & mcwfn_;
    std::shared_ptr<SquareMatrix<MatsT>> orbitalGradient_ = nullptr;

  public:
    
    OrbitalRotationSettings & settings;

    // delete default Constructor
    OrbitalRotation() = delete;
    OrbitalRotation(
      MCWaveFunction<MatsT,IntsT> & mcwfn, 
      OrbitalRotationSettings & input_settings):
      mcwfn_(mcwfn), settings(input_settings) {
      
      auto & mopart = mcwfn_.MOPartition;

      if (mopart.nCorrO == 0) CErr("the correlated space cannot be empty in orbital rotation"); 
      
      if (mopart.nInact == 0) {
        std::cout << "  No Inactive in Orbital Rotations" << std::endl;
        settings.rotate_inact_correlated = false;
        settings.rotate_inact_virtual    = false;
      }
      
      if (mopart.nFVirt == 0) {
        std::cout << "  No Frozen Virtual in Orbital Rotations" << std::endl;
        settings.rotate_correlated_virtual = false; 
        settings.rotate_inact_virtual      = false;
      }
       

    
    }; // constructor
    
    // use default copy and move constructor
    OrbitalRotation(const OrbitalRotation<MatsT,IntsT> &) = default;
    OrbitalRotation(OrbitalRotation<MatsT,IntsT> &&)      = default;
	
	~OrbitalRotation() = default;
    
    // Perform one-step Orbital Rotation based on 
    // the underlying oneRDM and twoRDM
    double computeOrbGradient(EMPerturbation &, SquareMatrix<MatsT> &, InCore4indexTPI<MatsT> &);
	void rotateMO(EMPerturbation &, SquareMatrix<MatsT> &, InCore4indexTPI<MatsT> &);
	
    // generate improved virtual orbitals
    void generateIVOs(EMPerturbation &, SquareMatrix<MatsT> &);
    
	void formGeneralizedFock1(EMPerturbation &, SquareMatrix<MatsT> &, MatsT *, const std::string &, bool deltaPQ = false);
	void formGeneralizedFock2(EMPerturbation &, SquareMatrix<MatsT> &, InCore4indexTPI<MatsT> &, 
      MatsT *, const std::string &);
    void computeOrbOrbHessianDiag(EMPerturbation &, SquareMatrix<MatsT> &, InCore4indexTPI<MatsT> &, 
      MatsT *);
    //void computeOrbOrbHessian(MatsT *);
    
  }; // class OrbitalRotation
  

}; // namespace ChronusQ

// Include declaration for specialization of OrbitalRotation

