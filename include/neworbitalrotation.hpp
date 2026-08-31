/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <posthartreefock.hpp>
#include <matrix.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <util/print.hpp>

// Orbital Rotation Headers

//#define DEBUG_OrbitalRotation

namespace ChronusQ {
  
enum OrbitalRotationAlgorithm {
  ORB_ROT_APPROX_QUASI_2ND_ORDER,
  ORB_ROT_QUASI_2ND_ORDER,
  ORB_ROT_2ND_ORDER,
}; // enum OrbitalRotationAlgorithm

/* 
 * Brief Definition of OrbitalRotationSettings
 */

struct OrbitalRotationSettings {
  
  // blocks for rotation
  bool rotate_within_correlated   = false;
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

  void print(bool fourComp) const {
    if (alg == OrbitalRotationAlgorithm::ORB_ROT_APPROX_QUASI_2ND_ORDER) {
      FormattedLine(std::cout,"  Orbital Rotation Algorithm:",  "Approximated Quasi-Second Order");
    } else if (alg == OrbitalRotationAlgorithm::ORB_ROT_QUASI_2ND_ORDER) {
      FormattedLine(std::cout,"  Orbital Rotation Algorithm:",  "Quasi-Second Order");
    } else if (alg == ORB_ROT_2ND_ORDER) {
      FormattedLine(std::cout,"  Orbital Rotation Algorithm:",  "Second Order");
    } else CErr("NYI Orbital Rotation Algorithm");
 
    FormattedLine(std::cout, "  Rotation Blocks:");
    FormattedLine(std::cout, "  - Within Correlated:",               rotate_within_correlated);
    FormattedLine(std::cout, "  - Between Inactive and Correlated:", rotate_inact_correlated);
    FormattedLine(std::cout, "  - Between Inactive and Virtual:",    rotate_inact_virtual);
    FormattedLine(std::cout, "  - Between Correlated and Virtual:",  rotate_correlated_virtual);
    if(fourComp) 
      FormattedLine(std::cout, "  - Between Negative and Positive:",  rotate_negative_positive);
  } // print 
}; // struct OrbitalRotationSettings
  
/* 
 * \brief the NewOrbitalRotation class. The class can perform 
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
class NewOrbitalRotation { 
  
protected:
  
  PostHartreeFock<MatsT,IntsT> & postHF_;
  std::shared_ptr<cqmatrix::Matrix<MatsT>> orbitalGradient_ = nullptr;

public:
  
  OrbitalRotationSettings& settings;

  // delete default Constructor
  NewOrbitalRotation() = delete;
  NewOrbitalRotation(
    PostHartreeFock<MatsT,IntsT>& postHF, 
    OrbitalRotationSettings& input_settings):
    postHF_(postHF), settings(input_settings) {
    
    const auto& corrS = postHF_.corrSpace;

    if (corrS.nCorrO == 0) CErr("the correlated space cannot be empty in orbital rotation"); 
    
    if (corrS.nInact == 0) {
      std::cout << "  No Inactive in Orbital Rotations" << std::endl;
      settings.rotate_inact_correlated = false;
      settings.rotate_inact_virtual    = false;
    }
    
    if (corrS.nSVirt == 0) {
      std::cout << "  No Frozen Virtual in Orbital Rotations" << std::endl;
      settings.rotate_correlated_virtual = false; 
      settings.rotate_inact_virtual      = false;
    }
  
  }; // constructor
  
  // use default copy and move constructor
  NewOrbitalRotation(const NewOrbitalRotation<MatsT,IntsT> &) = default;
  NewOrbitalRotation(NewOrbitalRotation<MatsT,IntsT> &&)      = default;
  
  ~NewOrbitalRotation() = default;
  
  // Perform one-step Orbital Rotation based on 
  // the underlying oneRDM and twoRDM
  double computeOrbGradient(EMPerturbation &, cqmatrix::Matrix<MatsT> &, InCore4indexTPI<MatsT> &);
  void rotateMO(EMPerturbation &, cqmatrix::Matrix<MatsT> &, InCore4indexTPI<MatsT> &);
  
  // generate improved virtual orbitals
  void generateIVOs(EMPerturbation &, cqmatrix::Matrix<MatsT> &);
  
  void formGeneralizedFock1(EMPerturbation &, cqmatrix::Matrix<MatsT> &, MatsT *, const std::string &, bool deltaPQ = false);
  void formGeneralizedFock2(EMPerturbation &, cqmatrix::Matrix<MatsT> &, InCore4indexTPI<MatsT> &, 
    MatsT *, const std::string &);
  void computeOrbOrbHessianDiag(EMPerturbation &, cqmatrix::Matrix<MatsT> &, InCore4indexTPI<MatsT> &, 
    MatsT *);
  //void computeOrbOrbHessian(MatsT *);
  
}; // class NewOrbitalRotation
  

} // namespace ChronusQ
