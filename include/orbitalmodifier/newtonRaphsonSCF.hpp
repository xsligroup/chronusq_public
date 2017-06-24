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

#include <orbitalmodifier/orbitaloptimizer.hpp>
#include <cqlinalg/matfunc.hpp>

//#define _NRSCF_DEBUG
//#define _NRSCF_DEBUG_GRAD
//#define _NRSCF_DEBUG_UNITARY
//#define _NRSCF_DEBUG_BFGS
//#define _NRSCF_DEBUG_SR1
//#define _NRSCF_PRINT_MOS



namespace ChronusQ {

/*
 *    Brief: This is a struct to hold the information about the nonredundant 
 *           rotations in NewtonRaphsonSCF
 */
struct NRRotOptions{
  std::pair<std::set<size_t>,std::set<size_t>> rotIndices;   ///< Which indices are being rotated e.g. first=>occupied second=>virtual
  size_t spaceIndex = 0;                    ///< Which set of fock,den,mo matrices it is associated with
};

/*
 *  Brief: This is a simple function to compute the rotation indices for 
 *         Slater determinant wave functions.  In the std::pair, the first set
 *         is the set of occupied indices, and the second set is the set of 
 *         virtual indices. These are used to determine the nonredundant 
 *         parameters in the NewtonRaphsonSCF object. The implementation
 *         is in include/orbitalmodifier/newtonraphsonscf/impl.hpp
 */ 
std::pair<std::set<size_t>,std::set<size_t>> ssNRRotIndices( size_t nOcc, size_t NB, size_t shift=0 );

template<typename MatsT>
class NewtonRaphsonSCF : public OrbitalOptimizer<MatsT> {
private:
  size_t nParam;                  ///< Total number of independent parameters
  bool refMOsAllocated = false;   ///< Whether the reference MOs were allocated yet
  bool storedRef       = false;   ///< Save the MO's on the first iteration
  MatsT* orbRot;                  ///< Orbital Rotation Parameters
  MatsT* orbGrad;                 ///< Orbital Gradient
  MatsT* orbDiagHess;             ///< Diagonal Hessian Approximation
  std::vector<cqmatrix::Matrix<MatsT>> refMO;   ///< The reference set of MO's that are being rotated
  std::vector<NRRotOptions> rotOpt;         ///< Data structure to generate rotation matrices

  // Quasi-Newton Data Structures
  std::vector<MatsT*> qnOrbRot;             ///< The Set of variable iterates for Quasi-Newton
  std::vector<MatsT*> qnOrbGrad;            ///< The set of gradients for Quasi-Newton

public:
  // Constructor
  NewtonRaphsonSCF() = delete;
  NewtonRaphsonSCF(std::vector<NRRotOptions> nrrot, SCFControls sC, MPI_Comm comm, OrbitalModifierDrivers<MatsT> modOpt):
    OrbitalOptimizer<MatsT>(sC, comm, modOpt), rotOpt(nrrot) {

    sanityChecks();

    // Compute the total number of nonredundant parameters
    nParam = 0;
    for( auto& rO : rotOpt )
      nParam += rO.rotIndices.first.size() * rO.rotIndices.second.size();

    // Alloc Quasi-Newton matrices
    alloc();
  }

  // Destructor
  ~NewtonRaphsonSCF() {
    CQMemManager::get().free(orbRot, orbGrad, orbDiagHess);

    if( this->scfControls.nrAlg == QUASI_BFGS or this->scfControls.nrAlg == QUASI_SR1 ) {
      for( auto* p : qnOrbRot )
        CQMemManager::get().free(p);
      for( auto* p : qnOrbGrad )
        CQMemManager::get().free(p);
    }
  }

  // Functions
  void sanityChecks();
  void alloc();
  void getNewOrbitals(EMPerturbation&, vecMORef<MatsT>&, vecEPtr&);
  void NewtonRaphsonIteration();
  void printRunHeader(std::ostream&, EMPerturbation&) const;

  // Newton-Raphson Approximation Steps
  void fullNRStep();
  void qnBFGSStep();
  void qnSR1Step();
  void gradDescentStep();

  // Common Functions
  void computeGradient(vecMORef<MatsT>&);
  void computeDiagHess(vecEPtr&);
  void rotateMOs(vecMORef<MatsT>&);
  void saveRefMOs(vecMORef<MatsT>& mo);
  std::vector<cqmatrix::Matrix<MatsT>> computeUnitary();


  // Line Search and step functions
  void takeStep(MatsT*);

  // Quasi-Newton Functions
  void qnSetup();
  void computeBFGS(size_t N, const std::vector<MatsT*>& x, const std::vector<MatsT*>& g, MatsT* dx);
  void computeSR1(size_t N, const std::vector<MatsT*>& x, const std::vector<MatsT*>& g, MatsT* dx);
  void computeMatrixInverse(const size_t, MatsT*, const size_t, const double);

  template<typename MatsU>
  void printParamVec(std::string, MatsU*) const;

  // Convergence Testing functions
  double computeFDCConv();
};
};   // namespace ChronusQ
