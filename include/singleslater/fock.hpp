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
#include <corehbuilder.hpp>
#include <fockbuilder.hpp>
#include <physcon.hpp>

#include <util/timer.hpp>
#include <cqlinalg/blasext.hpp>

#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>

#include <intermediates.hpp>


//#define _DEBUGORTHO
//#define __DEBUGTPB__

namespace ChronusQ {

  /**
   *  \brief Forms the Fock matrix for a single slater determinant using
   *  the 1PDM.
   *
   *  \param [in] increment Whether or not the Fock matrix is being 
   *  incremented using a previous density
   *
   *  Populates / overwrites fock strorage
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formFock(
    EMPerturbation &pert, bool increment, double xHFX) {

    fockBuilder->formFock(*this, pert, increment, xHFX);

  }; // SingleSlater::fockFock


  /**
   *  \brief Compute the Core Hamiltonian.
   *
   *  \param [in] typ Which Hamiltonian to build
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formCoreH(EMPerturbation& emPert, bool save) {

    size_t NB = basisSet().nBasis;
    if( nC == 4 ) NB = 2 * NB;

    if( coreH != nullptr ) {
      coreH->clear();

    } else {

      if(not iCS and nC == 1 and basisSet().basisType == COMPLEX_GIAO)
        coreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false);
      else if(nC == 2 or nC == 4)
        coreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true);
      else
        coreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false, false);

    }

    if (MPIRank(comm) == 0) {

    ProgramTimer::tick("Form Core H");


    // Make a copy of HamiltonianOptions
    HamiltonianOptions hamiltonianOptions = coreHBuilder->getHamiltonianOptions();

    // Prepare one-electron integrals
    std::vector<std::pair<OPERATOR,size_t>> ops;
    if (hamiltonianOptions.basisType == REAL_GTO)
      ops = {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0},
             {LEN_ELECTRIC_MULTIPOLE,3},
             {VEL_ELECTRIC_MULTIPOLE,3}, {MAGNETIC_MULTIPOLE,2}};
    else
      ops = {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0},
             {LEN_ELECTRIC_MULTIPOLE,3},
             {MAGNETIC_MULTIPOLE,1}};

    // Multipole integrals NYI for 4C
    // if (nC == 4) ops.resize(3);
    // ADL implemeting 4c dipole
    if (nC == 4) ops = {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0}, {LEN_ELECTRIC_MULTIPOLE,1}};


    // In case of X2C coreHBuilder, here we only compute
    // non-relativistic one electron integrals for contracted basis functions.
    // Relativistic integrals will be computed for uncontracted basis functions
    // in X2C CoreHBuilder.
    if (hamiltonianOptions.x2cType != X2C_TYPE::OFF) {
      hamiltonianOptions.OneEScalarRelativity = false;
      hamiltonianOptions.OneESpinOrbit = false;
    }

    this->aoints_->computeAOOneP(this->molecule(),
        this->basisSet(),emPert, ops, hamiltonianOptions); // compute the necessary 1e ints

    // Compute core Hamiltonian
    coreHBuilder->computeCoreH(emPert,coreH);


    // Compute Orthonormalization trasformations
    computeOrtho();

    std::string prefix = this->particle.charge < 0 ? "INTS/" : "PINTS/";

    if (save) CQIntermediates::getInstance().addData("INTS/CORE_HAMILTONIAN", coreH);
    // Save the Core Hamiltonian
    if( savFile.exists() && save) {

      try{
        savFile.safeWriteData("INTS/CORE_HAMILTONIAN", *coreH);
      } catch(...) {
        CErr("Error saving core Hamiltonian. Please use -s. See Running ChronusQ section of wiki.");
      }

    }

    ProgramTimer::tock("Form Core H");
    } // End building CoreH on root

#ifdef CQ_ENABLE_MPI
    // BCast core Hamiltonian to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the core Hamiltonian ***\n";
      for(auto mat : coreH->SZYXPointers())
        MPIBCast(mat,NB*NB,0,comm);
    }
#endif

  }; // SingleSlater<MatsT,IntsT>::computeCoreH


  template <typename MatsT, typename IntsT>
  std::vector<double> SingleSlater<MatsT,IntsT>::getGrad(EMPerturbation& pert,
    bool equil, bool saveInts, double xHFX) {

    // Get constants
    size_t NB = basisSet().nBasis;
    size_t nSQ  = NB*NB;

    size_t nAtoms = this->molecule().nAtoms;
    size_t nGrad = 3*nAtoms;

    size_t nSp = fockMatrix->nComponent();
    bool hasXY = fockMatrix->hasXY();
    bool hasZ = fockMatrix->hasZ();


    // Total gradient
    std::vector<double> gradient(nGrad, 0.);

    HamiltonianOptions opts = this->aoints_->options_;

    auto printGrad = [&](std::string name, std::vector<double>& vecgrad) {
      std::cout << name << std::endl;
      std::cout << std::setprecision(12);
      for( auto iAt = 0; iAt < nAtoms; iAt++ ) {
        std::cout << " Gradient@I = " << iAt << ":";
        for( auto iCart = 0; iCart < 3; iCart++ ) {
          std::cout << "  " << vecgrad[iAt*3 + iCart];
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    };


    // Core H contribution
    this->aoints_->computeGradInts(this->molecule_, basisSet_, pert,
      {{OVERLAP, 1},
       {KINETIC, 1},
       {NUCLEAR_POTENTIAL, 1}
       },
       opts
    );
    std::vector<double> coreGrad = coreHBuilder->getGrad(pert, *this);


    // 2e contribution
    this->aoints_->computeGradInts(this->molecule_, basisSet_, pert,
      {{ELECTRON_REPULSION, 1}},
      opts
    );
    std::vector<double> twoEGrad = fockBuilder->getGDGrad(*this, pert, xHFX);

    
    // Pulay gradient contribution
    bool useW = true;
    std::vector<double> pulayGrad = fockBuilder->getPulayGrad(*this, equil, useW);


    // Add the nuclear gradient and assemble the total gradient
    std::vector<double> nucGrad;
    for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) {
      size_t iAt = iGrad/3;
      size_t iXYZ = iGrad%3;
      nucGrad.push_back(this->molecule().nucRepForce[iAt][iXYZ]);
      gradient[iGrad] = coreGrad[iGrad] + twoEGrad[iGrad] + pulayGrad[iGrad] + nucGrad[iGrad];
    }

    //printGrad("Nuclear Gradient:", nucGrad);
    //printGrad("Core H Gradient:", coreGrad);
    //printGrad("G Gradient:", twoEGrad);
    //printGrad("Pulay Gradient:", pulayGrad);
    //printGrad("Total HF Gradient:", gradient);

    //this->onePDM->output(std::cout, "OnePDM in Gradient Contractions", true);

    return gradient;

  };


  /**
   *  \brief Allocate, compute and store the orthonormalization matricies 
   *  over the CGTO basis.
   *
   *  Computes either the Lowdin or Cholesky transformation matricies based
   *  on orthoType
   */ 
  template <typename MatsT, typename IntsT> 
  void SingleSlater<MatsT,IntsT>::computeOrtho() {

    size_t NB = this->basisSet().nBasis;
    if( nC == 4 ) NB = 2 * NB;
    size_t nSQ  = NB*NB;
    size_t NBC = this->nC*basisSet().nBasis;

    // Allocate scratch
    cqmatrix::Matrix<MatsT> overlapSpinor(NB);
    overlapSpinor.clear();
    cqmatrix::Matrix<MatsT> overlapAB(NBC);
    overlapAB.clear();

    // Copy the overlap over to scratch space
    if ( nC != 4 ) {
      std::copy_n(this->aoints_->overlap->pointer(),nSQ,overlapSpinor.pointer());
    } else if( nC == 4 ) {
      // HBL 4C May need a Ints type check (SetMat) to capture GIAO.
      SetMatRE('N',NB/2,NB/2,1.,
               reinterpret_cast<double*>(this->aoints_->overlap->pointer()),NB/2,
               overlapSpinor.pointer(),NB);
      SetMatRE('N',NB/2,NB/2,1./(2*SpeedOfLight*SpeedOfLight),
               reinterpret_cast<double*>(this->aoints_->kinetic->pointer()),NB/2,
               overlapSpinor.pointer()+NB*NB/2+NB/2,NB);
      //prettyPrintSmart(std::cout,"S Metric Spinor",overlapSpinor.pointer(),NB,NB,NB);
    }
    orthoSpinor = std::make_shared<Orthogonalization<MatsT>>(overlapSpinor, scfControls.linearDepTol);

    // Copy to block diagonal for alpha/beta basis
    if( nC > 1 ){
      SetMat('N',NBC/2,NBC/2,MatsT(1.),overlapSpinor.pointer(), NBC/2, overlapAB.pointer(),NBC);
      size_t disp = NBC/2 + NBC/2*NBC;
      SetMat('N',NBC/2,NBC/2,MatsT(1.),overlapSpinor.pointer(), NBC/2, overlapAB.pointer()+disp,NBC);
      //prettyPrintSmart(std::cout,"S Metric AB",overlapAB.pointer(),NBC,NBC,NBC);
    } else {
      overlapAB = overlapSpinor;
    }
    orthoAB = std::make_shared<Orthogonalization<MatsT>>(overlapAB, scfControls.linearDepTol);

  }; // computeOrtho

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::MOFOCK() {

    const size_t NB   = this->nAlphaOrbital();
    const size_t NBC  = nC * NB;

    if( fockMO.empty() ) {
      fockMO.emplace_back(NBC);
      if( nC == 1 and not iCS )
        fockMO.emplace_back(NBC);
    }

    if( MPIRank(comm) == 0 ) {

      if ( nC == 2 )
        fockMO[0] = fockMatrix->template spinGather<MatsT>();
      else if ( nC == 1 )
        fockMO = fockMatrix->template spinGatherToBlocks<MatsT>(false, not iCS);

      fockMO[0] = fockMO[0].transform('N',this->mo[0].pointer(),NBC,NBC);

      if( nC == 1 and not iCS ) {

        fockMO[1] = fockMO[1].transform('N',this->mo[1].pointer(),NBC,NBC);

      }

    }

  //prettyPrintSmart(std::cout,"MOF",fockMO[0].pointer(),NBC,NBC,NBC);

#ifdef CQ_ENABLE_MPI

    if(MPISize(comm) > 1) {
      std::cerr  << "  *** Scattering MO-FOCK ***\n";
      for(int k = 0; k < fockMO.size(); k++)
        MPIBCast(fockMO[k].pointer(),NBC*NBC,0,comm);
    }

#endif
  };

}; // namespace ChronusQ

