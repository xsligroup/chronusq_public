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

namespace ChronusQ {

template<typename MatsT>
class ConventionalSCF : public OrbitalOptimizer<MatsT> {

protected:
  // Orthogonal Fock/Density Matrices
  std::vector<cqmatrix::Matrix<MatsT>> fockMatrixOrtho;
  std::vector<cqmatrix::Matrix<MatsT>> onePDMOrtho;
  std::vector<cqmatrix::Matrix<MatsT>> orbGrad;

  // DIIS matrices
  std::vector<std::vector<cqmatrix::Matrix<MatsT>>> diisFock;     ///< List of AO Fock matrices for DIIS extrap
  std::vector<std::vector<cqmatrix::Matrix<MatsT>>> diisOnePDM;   ///< List of AO Density matrices for DIIS extrap
  std::vector<std::vector<cqmatrix::Matrix<MatsT>>> diisError;    ///< List of orthonormal [F,D] for DIIS extrap
  std::vector<double> diisEnergy;                             ///< List of energies for EDIIS
  std::shared_ptr<cqmatrix::Matrix<double>> diisBMat;             ///< Matrix of couplings between EDIIS elements
                                                              ///<   diisBMat needs to be a shared_ptr because there is
                                                              ///<   no resize/default constructor for cqmatrix::Matrix.

  // Damping Matrices
  std::vector<cqmatrix::Matrix<MatsT>> prevFock;     ///< AO Fock from the previous SCF iteration
  std::vector<cqmatrix::Matrix<MatsT>> prevOnePDM;   ///< AO Density from the previous SCF iteration

public:
  // Constructor
  ConventionalSCF() = delete;
  ConventionalSCF(SCFControls sC, MPI_Comm comm, OrbitalModifierDrivers<MatsT> mod):
    OrbitalOptimizer<MatsT>(sC, comm, mod) {

    // Allocate ortho Fock and Den
    vecShrdPtrMat<MatsT> fock = this->orbitalModifierDrivers.getFock();
    for( auto& f : fock )
      fockMatrixOrtho.emplace_back(f->dimension());
    vecShrdPtrMat<MatsT> den = this->orbitalModifierDrivers.getOnePDM();
    for( auto& d : den )
      onePDMOrtho.emplace_back(d->dimension());
    for( auto& d : den )
      orbGrad.emplace_back(d->dimension());

    if( this->scfControls.doExtrap ) allocExtrapStorage();
  };

  // Destructor
  ~ConventionalSCF() {
    diisFock.clear();
    diisOnePDM.clear();
    diisError.clear();
    diisEnergy.clear();
    prevFock.clear();
    prevOnePDM.clear();
    fockMatrixOrtho.clear();
    onePDMOrtho.clear();
    diisBMat = nullptr;
  };

  // ModifyOrbital Functions
  void getNewOrbitals(EMPerturbation&, vecMORef<MatsT>&, vecEPtr&);
  void printRunHeader(std::ostream&, EMPerturbation&) const;

  // Orthogonalization Functions
  void ao2orthoFock(vecShrdPtrMat<MatsT> fock = {});
  void ao2orthoDen(vecShrdPtrMat<MatsT> den = {});
  void ortho2aoMOs(vecMORef<MatsT>&);
  void ao2orthoMOs(vecMORef<MatsT>&);

  // SCF extrapolation functions (see include/singleslater/extrap.hpp for docs)
  void allocExtrapStorage();
  void diagOrthoFock(vecMORef<MatsT>&, vecEPtr&);
  void modifyFock(EMPerturbation&);
  void fockDamping();
  void scfCDIIS(size_t, size_t);
  void scfEDIIS(size_t, size_t, EMPerturbation&);
  void scfCEDIIS(size_t, size_t, EMPerturbation&);
  void ediisErrorMetric(size_t, size_t);
  void diisCombineMat(std::vector<MatsT>, size_t);
  void computeOrbGradient(std::vector<cqmatrix::Matrix<MatsT>>&);
  void FDCommutator(std::vector<cqmatrix::Matrix<MatsT>>&);
  double computeFDCConv();

};   // ConventionalSCF
};   // namespace ChronusQ
