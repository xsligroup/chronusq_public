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

#include <orbitalmodifiernew/orbitaloptimizernew.hpp>

namespace ChronusQ {

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
class ConventionalSCFNew : public OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT> {

protected:

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
  ConventionalSCFNew() = delete;
  ConventionalSCFNew(SCFControls sC, singleSlaterT<MatsT,IntsT> &referenceSS, MPI_Comm comm):
  OrbitalOptimizerNew<singleSlaterT,MatsT,IntsT>(sC, referenceSS, comm) {

    vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
    for( auto& d : onePDM )
      orbGrad.emplace_back(d->dimension());

    if( this->scfControls.doExtrap ) allocExtrapStorage();
  };

  // Destructor
  ~ConventionalSCFNew() {
    diisFock.clear();
    diisOnePDM.clear();
    diisError.clear();
    diisEnergy.clear();
    prevFock.clear();
    prevOnePDM.clear();
    diisBMat = nullptr;
  };

  // ModifyOrbital Functions
  void getNewOrbitals(EMPerturbation&);
  void printRunHeader(std::ostream&, EMPerturbation&) const;


  // SCF extrapolation functions (see include/singleslater/extrap.hpp for docs)
  void allocExtrapStorage();
  void modifyFock(EMPerturbation&);
  void fockDamping();
  void scfCDIIS(size_t, size_t);
  void scfEDIIS(size_t, size_t, EMPerturbation&);
  void scfCEDIIS(size_t, size_t, EMPerturbation&);
  void ediisErrorMetric(size_t, size_t);
  void diisCombineMat(std::vector<MatsT>, size_t);
  //void computeOrbGradient(std::vector<cqmatrix::Matrix<MatsT>>&);
  void FDCommutator(std::vector<cqmatrix::Matrix<MatsT>>&);
  double computeFDCConv();

};   // ConventionalSCF
};   // namespace ChronusQ
