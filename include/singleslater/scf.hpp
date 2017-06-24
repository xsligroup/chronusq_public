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

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <singleslater.hpp>
#include <util/math.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

/**
 *  \brief Saves the current state of wave function
 *
 *  Saves a copy of the current AO 1PDM and orthonormal Fock
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::saveCurrentState(bool saveMO) {

  ROOT_ONLY(comm);

  // Checkpoint if file exists
  if( savFile.exists() ) {

    size_t NB  = this->nAlphaOrbital();
    size_t NBC = this->nC * NB;

    size_t t_hash = std::is_same<MatsT, double>::value ? 1 : 2;

    // Save Field type
    std::string prefix = "SCF/";
    if( this->particle.charge == 1.0 ) prefix = "PROT_" + prefix;

    savFile.safeWriteData(prefix + "FIELD_TYPE", &t_hash, {1});

    savFile.safeWriteData(prefix + "1PDM", *this->onePDM);

    savFile.safeWriteData(prefix + "FOCK", *fockMatrix);

    savFile.safeWriteData(prefix + "1PDM_ORTHO", *onePDMOrtho);

    savFile.safeWriteData(prefix + "FOCK_ORTHO", *fockMatrixOrtho);

    savFile.safeWriteData(prefix + "ORTHO", orthoSpinor->forwardPointer()->pointer(), {NB, NB});

    savFile.safeWriteData(prefix + "ORTHO_INV", orthoSpinor->backwardPointer()->pointer(), {NB, NB});

    // Save MOs
    if (saveMO) {
      savFile.safeWriteData(prefix + "MO1", this->mo[0].pointer(), {NBC, NBC});
      if (this->nC == 1 and not this->iCS) savFile.safeWriteData(prefix + "MO2", this->mo[1].pointer(), {NBC, NBC});
    }

    // Save Energies
    savFile.safeWriteData(prefix + "TOTAL_ENERGY", &this->totalEnergy, {1});
    savFile.safeWriteData(prefix + "ONE_BODY_ENERGY", &this->OBEnergy, {1});
    savFile.safeWriteData(prefix + "MANY_BODY_ENERGY", &this->MBEnergy, {1});

    // Save Multipoles
    savFile.safeWriteData(prefix + "LEN_ELECTRIC_DIPOLE", &this->elecDipole[0], {3});
    savFile.safeWriteData(prefix + "LEN_ELECTRIC_QUADRUPOLE", &this->elecQuadrupole[0][0], {3, 3});
    savFile.safeWriteData(prefix + "LEN_ELECTRIC_OCTUPOLE", &this->elecOctupole[0][0][0], {3, 3, 3});

    // Save Spin
    savFile.safeWriteData(prefix + "S_EXPECT", &this->SExpect[0], {3});
    savFile.safeWriteData(prefix + "S_SQUARED", &this->SSq, {1});
  }

};   // SingleSlater<MatsT>::saveCurrentState()

  template <typename MatsT, typename IntsT>
  bool SingleSlater<MatsT, IntsT>::checkStability() {

    double W;
    MatsT* J;
    std::tie(W,J) = this->getStab();
    std::cout << "  * LOWEST STABILITY EIGENVALUE " << 
      std::scientific << W << "\n";

    if( W < 0. and std::abs(W) > 1e-08 )
      std::cout << "  * LOWEST STABILITY EIGENVALUE NEGATIVE: " 
        << "PERFORMING THOULESS ROTATION\n";
    else { 

      std::cout << "  * LOWEST STABILITY EIGENVALUE POSITIVE: " 
        << "WAVE FUNCTION 2nd ORDER STABLE\n";

      CQMemManager::get().free(J); return true; 

    }

    const size_t NB   = this->nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (nC == 2) ? this->nO : this->nOA;
    const size_t NV    = (nC == 2) ? this->nV : this->nVA;
    const size_t nOAVA = this->nOA * this->nVA;
    const size_t nOBVB = this->nOB * this->nVB;


    MatsT* ROT    = CQMemManager::get().malloc<MatsT>(NBC2);
    MatsT* EXPROT = CQMemManager::get().malloc<MatsT>(NBC2);
    std::fill_n(ROT,NBC2,0.);

    for(auto i = 0ul, ai = 0ul; i < NO;  i++)
    for(auto a = NO;            a < NBC; a++, ai++) {

      ROT[a + i*NBC] =  J[ai];
      ROT[i + a*NBC] = -SmartConj(J[ai]);

    }      

    // FIXME: need to generalize MatExp to take non-hermetian and real
    // matricies
    //MatExp('D',NBC,T(-1.),ROT,NBC,EXPROT,NBC);

    // Taylor
    MatsT s = 1.;
    std::copy_n(ROT,NBC2,EXPROT); // n = 1
    blas::scal(NBC2,-s,EXPROT,1);
    for(auto j = 0; j < NBC; j++) EXPROT[j*(NBC+1)] += 1.; // n = 0

    MatsT* SCR  = CQMemManager::get().malloc<MatsT>(NBC2);
    MatsT* SCR2 = CQMemManager::get().malloc<MatsT>(NBC2);
    std::copy_n(ROT,NBC2,SCR);

    size_t tayMax = 30; 
    for(auto n = 2; n <= tayMax; n++) {

      MatsT* M = nullptr;
      if( n % 2 ) {
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBC,NBC,NBC,MatsT(1.),ROT,NBC,SCR2,NBC,MatsT(0.),SCR,NBC);
        M = SCR;
      } else {
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBC,NBC,NBC,MatsT(1.),ROT,NBC,SCR,NBC,MatsT(0.),SCR2,NBC);
        M = SCR2;
      }

      MatsT fact = std::pow(-s,n)/factorial(n);
      MatAdd('N','N',NBC,NBC,MatsT(1.),EXPROT,NBC,fact,M,NBC,
        EXPROT,NBC);

    }

    /*
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NBC,NBC,NBC,T(1.),EXPROT,NBC,EXPROT,NBC,T(0.),SCR,NBC);
   // prettyPrintSmart(std::cerr,"ROT",ROT,NBC,NBC,NBC);
   // prettyPrintSmart(std::cerr,"EXPROT",EXPROT,NBC,NBC,NBC);
    prettyPrintSmart(std::cout,"SCR",SCR,NBC,NBC,NBC);
   // CErr();
   */


    // MO1 = MO1 * EXPROT
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBC,NBC,NBC,MatsT(1.),this->mo[0].pointer(),NBC,EXPROT,NBC,MatsT(0.),
      ROT,NBC);
    std::copy_n(ROT,NBC2,this->mo[0].pointer());




    orthoAOMO();
    this->formDensity();
    CQMemManager::get().free(J,ROT,EXPROT);

    return false;

  }

  



  /**
   *  \brief Diagonalize the orthonormal fock matrix
   *
   *  General purpose routine which diagonalizes the orthonormal
   *  fock matrix and stores a set of orthonormal MO coefficients
   *  (in WaveFunction::mo1 and possibly WaveFunction::mo2) and
   *  orbital energies. General for both 1 and 2 spin components
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::diagOrthoFock() {

    ROOT_ONLY(comm); 
    size_t NB = this->nAlphaOrbital() * nC;
    size_t NB2 = NB*NB;
    bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT,IntsT>>(fockBuilder) != nullptr);


    // Copy over the fockMatrixOrtho into MO storage
    if(nC == 1 and iCS) 
      this->mo = fockMatrixOrtho->template spinGatherToBlocks<MatsT>(false,false);
    else if(iRO)
      this->mo[0] = fockMatrixOrtho->S();
    else if(nC == 1)
      this->mo = fockMatrixOrtho->template spinGatherToBlocks<MatsT>(false);
    else {
      this->mo[0] = fockMatrixOrtho->template spinGather<MatsT>();
    }

    // Diagonalize the Fock Matrix
    int INFO = HermetianEigen('V', 'L', NB, this->mo[0].pointer(), NB, this->eps1);
    if( INFO != 0 ) CErr("HermetianEigen failed in Fock1",std::cout);

    if(iRO) {
      this->mo[1] = this->mo[0]; // for ROHF
      std::copy_n(this->eps1, NB, this->eps2);
    } else if(nC == 1 and not iCS) {
      INFO = HermetianEigen('V', 'L', NB, this->mo[1].pointer(), NB, this->eps2);
      if( INFO != 0 ) CErr("HermetianEigen failed in Fock2",std::cout);
    }

#if 0
    printMO(std::cout);
#endif

};   // SingleSlater<MatsT>::diagOrthoFock

/**
 *  \brief Diagonalize the AO fock matrix
 *
 *  General purpose routine which diagonalizes the orthonormal
 *  fock matrix and stores a set of orthonormal MO coefficients
 *  (in WaveFunction::mo1 and possibly WaveFunction::mo2) and
 *  orbital energies. General for both 1 and 2 spin components
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::diagAOFock() {

  ROOT_ONLY(comm);
  size_t NB  = this->nAlphaOrbital() * nC;
  size_t NB2 = NB * NB;
  bool iRO   = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr);

  // Copy over the fockMatrix into MO storage
  if( nC == 1 and iCS )
    this->mo = fockMatrix->template spinGatherToBlocks<MatsT>(false, false);
  else if( iRO )
    this->mo[0] = fockMatrix->S();
  else if( nC == 1 )
    this->mo = fockMatrix->template spinGatherToBlocks<MatsT>(false);
  else {
    this->mo[0] = fockMatrix->template spinGather<MatsT>();

    prettyPrintSmart(std::cout, "AO Fock", this->mo[0].pointer(), NB, NB, NB);
  }

  // Diagonalize the Fock Matrix
  int INFO = HermetianEigen('V', 'L', NB, this->mo[0].pointer(), NB, this->eps1);
  if( INFO != 0 ) CErr("HermetianEigen failed in Fock1", std::cout);

  if( iRO ) {
    this->mo[1] = this->mo[0];   // for ROHF
    std::copy_n(this->eps1, NB, this->eps2);
  } else if( nC == 1 and not iCS ) {
    INFO = HermetianEigen('V', 'L', NB, this->mo[1].pointer(), NB, this->eps2);
    if( INFO != 0 ) CErr("HermetianEigen failed in Fock2", std::cout);
  }

#if 0
    printMO(std::cout);
#endif

};   // SingleSlater<MatsT>::diagOrthoFock

/**
 *  \brief Transforms all of the spin components of the AO fock
 *  matrix to the orthonormal basis.
 *
 *  Populates / overwrites fockMatrixOrtho storage
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::ao2orthoFock() {

  ROOT_ONLY(comm);

  *fockMatrixOrtho = orthoSpinor->nonortho2ortho(*fockMatrix);

};   // SingleSlater<MatsT>::ao2orthoFock

/**
 *  \brief Transforms all of the spin components of the AO density
 *  matrix to the orthonormal basis.
 *
 *  Populates / overwrites onePDMOrtho storage
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::ao2orthoDen() {

  ROOT_ONLY(comm);

  // NOTE: Density transforms the opposite way as operators, same as the coeffs
  *onePDMOrtho = orthoSpinor->ortho2nonortho(*this->onePDM);

};   // SingleSlater<MatsT>::ao2orthoFock

/**
 *  \brief Transforms all of the spin compoenents of the orthonormal
 *  1PDM to the AO basis.
 *
 *  Populates / overwrites onePDM storage
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::ortho2aoDen() {

  ROOT_ONLY(comm);

  // NOTE: Density transforms the opposite way as operators, same as the coeffs
  *this->onePDM = orthoSpinor->nonortho2ortho(*onePDMOrtho);

#if 0
    print1PDMOrtho(std::cout);
#endif

};   // SingleSlater<MatsT>::ao2orthoFock

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::ortho2aoMOs() {

  if( MPIRank(comm) == 0 )
    orthoAB->ortho2nonorthoCoeffs(this->mo);

#ifdef CQ_ENABLE_MPI

  // Broadcast the updated MOs to all MPI processes
  if( MPISize(comm) > 1 ) {

      std::cerr  << "  *** Scattering the AO-MOs ***\n";
      size_t Nmo = this->mo[0].dimension();
      MPIBCast(this->mo[0].pointer(),Nmo*Nmo,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->mo[1].pointer(),Nmo*Nmo,0,comm);

      std::cerr  << "  *** Scattering EPS ***\n";
      MPIBCast(this->eps1,Nmo,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->eps2,Nmo,0,comm);

      std::cerr  << "  *** Scattering FOCK ***\n";
      size_t fockDim = fockMatrix->dimension();
      for(MatsT *mat : fockMatrix->SZYXPointers())
        MPIBCast(mat,fockDim*fockDim,0,comm);

    }

#endif

  MOFOCK();   // Form the MO fock matrix

#if 0
    // Check proper orthonormalized wrt overlap

    size_t NB = basisSet().nBasis;
    MatsT* SCR2 = CQMemManager::get().malloc<MatsT>(this->nC*this->nC*NB*NB);
    MatsT* SCR3 = CQMemManager::get().malloc<MatsT>(this->nC*this->nC*NB*NB);


    // MO1 inner product
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,this->nC*NB,NB,T(1.),this->aoints_->overlap,NB,
      this->mo[0].pointer(),this->nC*NB,T(0.),SCR2,this->nC*NB);
    if(this->nC == 2)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,this->nC*NB,NB,T(1.),this->aoints_->overlap,NB,
        this->mo[0].pointer()+NB,this->nC*NB,T(0.),SCR2+NB,this->nC*NB);
   
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,this->nC*NB,this->nC*NB,this->nC*NB,T(1.),this->mo[0].pointer(),
      this->nC*NB,SCR2,this->nC*NB,T(0.),SCR3,this->nC*NB);

    for(auto i = 0; i < this->nC*NB; i++)
      SCR3[i*(this->nC*NB + 1)] -= 1.;


    std::cerr << "Error in orthonormazation of MO1 = " 
      << lapack::lange(lapack::Norm::Fro,this->nC*NB,this->nC*NB,SCR3,this->nC*NB)
      << std::endl;
             


    if(this->nC == 1 and not this->iCS) {
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,T(1.),this->aoints_->overlap,NB,this->mo[1].pointer(),NB,T(0.),SCR2,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,T(1.),this->mo[1].pointer(),NB,SCR2,NB,T(0.),SCR3,NB);

      for(auto i = 0; i < this->nC*NB; i++)
        SCR3[i*(this->nC*NB + 1)] -= 1.;

      std::cerr << "Error in orthonormazation of MO2 = " 
        << lapack::lange(lapack::Norm::Fro,NB,NB,SCR3,NB) << std::endl;
    }

    CQMemManager::get().free(SCR2,SCR3);

#endif

};   // SingleSlater<MatsT>::ortho2aoMOs

template<typename MatsT, typename IntsT>
std::vector<NRRotOptions> SingleSlater<MatsT,IntsT>::buildRotOpt(){

    // Generate NRRotationOptions
    std::vector<NRRotOptions> rotOpt;
    size_t NB = this->nC*this->basisSet().nBasis;
    if( this->nC == 1 and this->iCS ){
      rotOpt.emplace_back();
      rotOpt[0].spaceIndex = 0;
      rotOpt[0].rotIndices = ssNRRotIndices(this->nOA, NB);
    } else if( this->nC == 1 ) {
      rotOpt.emplace_back();
      rotOpt[0].spaceIndex = 0;
      rotOpt[0].rotIndices = ssNRRotIndices(this->nOA, NB);
      rotOpt.emplace_back();
      rotOpt[1].spaceIndex = 1;
      rotOpt[1].rotIndices = ssNRRotIndices(this->nOB, NB);
    } else if( this->nC == 2 ){
      rotOpt.emplace_back();
      rotOpt[0].spaceIndex = 0;
      rotOpt[0].rotIndices = ssNRRotIndices(this->nO, NB);
    } else if( this->nC == 4 ){
      size_t N = NB/2;
      rotOpt.emplace_back();
      rotOpt[0].spaceIndex = 0;
      rotOpt[0].rotIndices = ssNRRotIndices(this->nO, N, N); // only rotate the positive energy orbitals
    }
    return rotOpt;
}

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::initializeSCF() {

  this->moCoefficients.clear();

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(this->fockBuilder) != nullptr);

  // Setup MO reference vector
  if( iRO ) {
    this->moCoefficients.emplace_back(this->mo[0]);
  } else {
    for( auto& m : this->mo )
      this->moCoefficients.emplace_back(m);
  }

  // Setup Eigenvalue vector
  if( this->nC == 1 and not(iCS or iRO) ) {
    this->moEigenvalues = {this->eps1, this->eps2};
  } else {
    this->moEigenvalues = {this->eps1};
  }
}

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::runSCF(EMPerturbation& pert) {

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(this->fockBuilder) != nullptr);

  // Initialize properties
  //if( not std::dynamic_pointer_cast<SkipSCF<MatsT>>(this->orbitalModifier) ) getNewOrbitals();
  //this->computeProperties(pert);

  // Setup MO reference vector
  std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> moRefs;
  if( iRO ) {
    moRefs.emplace_back(this->mo[0]);
  } else {
    for( auto& m : this->mo )
      moRefs.emplace_back(m);
  }

  // Setup Eigenvalue vector
  std::vector<double*> epsVec;
  if( this->nC == 1 and not(iCS or iRO) ) {
    epsVec = {this->eps1, this->eps2};
  } else {
    epsVec = {this->eps1};
  }

  // Run modify orbitals
  this->orbitalModifier->runOrbitalModifier(pert, moRefs, epsVec);

  // Orthogonalize MO's for NewtonRaphsonSCF and check Stability
  if( std::dynamic_pointer_cast<NewtonRaphsonSCF<MatsT>>(this->orbitalModifier) ){
    this->orthoAOMO();
    if( this->scfControls.nrAlg == FULL_NR ){
      bool converged = this->checkStability();
      if( not converged ) {
        this->orbitalModifier->runOrbitalModifier(pert, moRefs, epsVec);
        this->orthoAOMO();
        converged = this->checkStability();
        if( not converged ) CErr("Newton-Raphson SCF failed to converge to a minimum");
      }
    }
  }

  this->ao2orthoFock();   // SCF Does not update the fockMatrixOrtho
  this->MOFOCK();
  saveCurrentState();

#ifdef CQ_ENABLE_MPI

  // Broadcast the updated MOs to all MPI processes
  if( MPISize(comm) > 1 ) {

      std::cerr  << "  *** Scattering the AO-MOs ***\n";
      size_t Nmo = this->mo[0].dimension();
      MPIBCast(this->mo[0].pointer(),Nmo*Nmo,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->mo[1].pointer(),Nmo*Nmo,0,comm);

      std::cerr  << "  *** Scattering EPS ***\n";
      MPIBCast(this->eps1,Nmo,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->eps2,Nmo,0,comm);

      std::cerr  << "  *** Scattering FOCK ***\n";
      size_t fockDim = fockMatrix->dimension();
      for(MatsT *mat : fockMatrix->SZYXPointers())
        MPIBCast(mat,fockDim*fockDim,0,comm);

      std::cerr  << "  *** Scattering the 1PDM ***\n";
      size_t denDim = this->onePDM->dimension();
      for(auto p : this->onePDM->SZYXPointers())
        MPIBCast(p,denDim*denDim,0,comm);
    }

#endif
};   // SingleSlater<MatsT,IntsT> :: runOrbitalModifier

/*
 *     Brief: Function to generate shared pointers to Fock Matrix for modify orbitals
 *            Here, the OrbitalModifier object cannot modify the fockMatrix in SingleSlater
 *            since they are temporary objects. Additionally, the interface assumes the
 *            matrices are spin gathered.
 */
template<typename MatsT, typename IntsT>
std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> SingleSlater<MatsT, IntsT>::getFock() {

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr);
  if( this->nC == 1 and iCS ) {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(MatsT(0.5) * fockMatrix->S())};
  } else if( this->nC == 1 and iRO ) {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(MatsT(1.) * fockMatrix->S())};
  } else if( this->nC == 1 ) {
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fA = std::make_shared<cqmatrix::Matrix<MatsT>>(
                    MatsT(0.5) * fockMatrix->S() + MatsT(0.5) * fockMatrix->Z()
                );
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fB = std::make_shared<cqmatrix::Matrix<MatsT>>(
                    MatsT(0.5) * fockMatrix->S() - MatsT(0.5) * fockMatrix->Z()
                );
    return {fA, fB};
  } else {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(fockMatrix->template spinGather<MatsT>())};
  }
};   // SingleSlater<MatsT,IntsT> :: getFock

/*
 *     Brief: Function to generate shared pointers to Fock Matrix for modify orbitals
 *            Here, the OrbitalModifier object cannot modify the fockMatrix in SingleSlater
 *            since they are copies of the objects. Additionally, the interface assumes the
 *            matrices are spin gathered.
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::setOnePDMOrtho(cqmatrix::Matrix<MatsT> *tempOnePDMOrtho) {

  if(nC == 1) {
    if(iCS) {
      *onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(tempOnePDMOrtho[0]);
    } else {
      *onePDMOrtho = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(tempOnePDMOrtho[0],tempOnePDMOrtho[1]);
    }
  } else {
    *onePDMOrtho = tempOnePDMOrtho[0].template spinScatter<MatsT>();
  }

};   // SingleSlater<MatsT,IntsT> :: setOnePDMOrtho

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::setOnePDMAO(cqmatrix::Matrix<MatsT> *tempOnePDMAO) {

  // Scatter to spin blocks on root process
  if( MPIRank(comm) == 0 ) {
    if(nC == 1) {
      if(iCS) {
        *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(tempOnePDMAO[0]);
      } else {
        *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(tempOnePDMAO[0],tempOnePDMAO[1]);
      }
    } else {
      *this->onePDM = tempOnePDMAO[0].template spinScatter<MatsT>();
    }
  } 

  #ifdef CQ_ENABLE_MPI
    // Broadcast the 1PDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      size_t NB  = this->nAlphaOrbital() * nC;
      std::cerr  << "  *** Scattering the 1PDMAO ***\n";
      for(auto p : this->onePDM->SZYXPointers())
        MPIBCast(p,NB*NB/nC/nC,0,comm);
    }
  #endif

};   // SingleSlater<MatsT,IntsT> :: setOnePDMAO

template<typename MatsT, typename IntsT>
std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> SingleSlater<MatsT, IntsT>::getOnePDM() {

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr);
  if( this->nC == 1 and iCS ) {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(MatsT(0.5) * this->onePDM->S())};
  } else if( this->nC == 1 and iRO ) {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(MatsT(0.5) * this->onePDM->S() + MatsT(0.5)*this->onePDM->Z())};
  } else if( this->nC == 1 ) {
    std::shared_ptr<cqmatrix::Matrix<MatsT>> dA = std::make_shared<cqmatrix::Matrix<MatsT>>(
            MatsT(0.5) * this->onePDM->S() + MatsT(0.5) * this->onePDM->Z()
         );
    std::shared_ptr<cqmatrix::Matrix<MatsT>> dB = std::make_shared<cqmatrix::Matrix<MatsT>>(
            MatsT(0.5) * this->onePDM->S() - MatsT(0.5) * this->onePDM->Z()
         );
    return {dA, dB};
  } else {
    return {std::make_shared<cqmatrix::Matrix<MatsT>>(this->onePDM->template spinGather<MatsT>())};
  }
};   // SingleSlater<MatsT,IntsT> :: getOnePDM

template<typename MatsT, typename IntsT>
std::vector<cqmatrix::Matrix<MatsT>> SingleSlater<MatsT, IntsT>::getOnePDMOrtho() {

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr);
  if( this->nC == 1 and iCS ) {
    return {cqmatrix::Matrix<MatsT>(MatsT(0.5) * onePDMOrtho->S())};
  } else if( this->nC == 1 and iRO ) {
    return {cqmatrix::Matrix<MatsT>(MatsT(0.5) * onePDMOrtho->S() + MatsT(0.5)*onePDMOrtho->Z())};
  } else if( this->nC == 1 ) {
    cqmatrix::Matrix<MatsT> dA = MatsT(0.5) * onePDMOrtho->S() + MatsT(0.5)*onePDMOrtho->Z();
    cqmatrix::Matrix<MatsT> dB = MatsT(0.5) * onePDMOrtho->S() - MatsT(0.5)*onePDMOrtho->Z();
    return {dA, dB};
  } else {
    return {cqmatrix::Matrix<MatsT>(this->onePDMOrtho->template spinGather<MatsT>())};
  }
};   // SingleSlater<MatsT,IntsT> :: getOnePDM

/*
 *     Brief: Function to return the orthogonalization pointers for SCF
 *
 */
template<typename MatsT, typename IntsT>
std::vector<std::shared_ptr<Orthogonalization<MatsT>>> SingleSlater<MatsT, IntsT>::getOrtho() {

  bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr);
  if( this->nC == 1 and not(iCS or iRO) ) {
    return {orthoAB, orthoAB};
  } else {
    return {orthoAB};
  }
};   // SingleSlater<MatsT,INtsT> :: getOrtho

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::printProperties() {
  printMOInfo(std::cout);
  //if( this->nC != 4 ) this->printMultipoles(std::cout);
  this->printMultipoles(std::cout);
  if( this->nC != 4 ) this->printSpin(std::cout);
  printMiscProperties(std::cout);
}

#ifdef TEST_MOINTSTRANSFORMER
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::MOIntsTransformationTest(EMPerturbation &pert) {
   
    // test on MO integral transfromations
    // MOIntsTransformer<MatsT, IntsT> N5TF(*this, INCORE_N5);
    MOIntsTransformer<MatsT, IntsT> N6TF(*this, INCORE_N6);  

    std::cout << "\n --------- Test on MO Ints Transformation----- \n" << std::endl;
    
    size_t NB  = this->nAlphaOrbital() * nC;
    size_t nMO = (this->nC == 4) ? NB / 2: NB;
    size_t nOcc = (this->nC == 1) ? this->nO/2: this->nO; 

    InCore4indexTPI<MatsT> N6MOERI(nOcc); 
    // InCore4indexTPI<MatsT> N5MOERI(nMO); 
    OnePInts<MatsT> hCore(nOcc); 

#if 1
    std::cout << "---- Test: Reconstruct SCF Energy" << std::endl; 
    N6TF.transformHCore(pert, hCore.pointer(), "ij");
    auto timeIdN6 = tick();
    N6TF.transformTPI(pert, N6MOERI.pointer(), "ijkl", false, false);
    auto timeDur = tock(timeIdN6);
    
    MatsT SCFEnergy = MatsT(0.);
    if(this->nC > 1) {
      for (auto i = 0; i < this->nO; i++) {
        SCFEnergy += hCore(i, i);
        for (auto j = 0; j < this->nO; j++)
          SCFEnergy += 0.5 * (N6MOERI(i, i, j, j) - N6MOERI(i, j, j, i)); 
      }
    } else {
      for (auto i = 0; i < this->nO/2; i++) {
        SCFEnergy += hCore(i, i);
        for (auto j = 0; j < this->nO/2; j++)
          SCFEnergy += N6MOERI(i, i, j, j) - 0.5 * N6MOERI(i, j, j, i); 
      }
      SCFEnergy *= 2.0;
    }

    std::cout << "SSFOCK_N6 SCF Energy:" << std::setprecision(16) << SCFEnergy << std::endl;
    // N6MOERI.output(std::cout, "SSFOCK_N6 ERI", true);

    std::cout << " - Time (N6) for transforming ijkl " <<  " = " << timeDur << " s\n"; 

    MatsT * h1e_ii  = CQMemManager::get().malloc<MatsT>(nOcc);
    MatsT * GDjj_ii = CQMemManager::get().malloc<MatsT>(nOcc);

    double fc1C = (this->nC == 1) ? 2.0 : 1.0;

    N6TF.transformHCore(pert, h1e_ii, "ii", true);
    N6TF.transformGD(pert, 'i', GDjj_ii, "jj", true, true, "WithInactive-i");  

    MatsT ECore = 0.;
    // compute core enenrgy
    for (auto i = 0ul; i < nOcc; i++) {
      ECore += h1e_ii[i] + 0.5 * GDjj_ii[i];
    }
    
    SCFEnergy = ECore * fc1C;

    std::cout << "SSFOCK_N6 SCF Energy using formFock exchange:" << std::setprecision(16) << SCFEnergy << std::endl;
#else     
    
    std::vector<std::string> testcases = {"ijkl", "abcd", "pqia", "ipab", "ijab", "pqrs"};
    for (auto & moType: testcases) {
      N6MOERI.clear();
      N5MOERI.clear();

      std::cout << "---- Test: " << moType << std::endl; 
      auto timeIdN6 = tick();
      N6TF.transformTPI(pert, N6MOERI.pointer(), moType, true, false);
      auto timeDur = tock(timeIdN6);
      
      std::cout << " - Time (N6) for transforming " << moType  <<  " = " << timeDur << " s\n"; 

      if (this->aoints_->TPITransAlg == TPI_TRANSFORMATION_ALG::INCORE_N6) {
        auto timeIdN5 = tick();
        N5TF.transformTPI(pert, N5MOERI.pointer(), moType, true, false);
        auto timeDur = tock(timeIdN5);
        std::cout << " - Time (N5) for transforming " << moType  <<  " = " << timeDur << " s\n"; 
#pragma omp parallel for schedule(static) collapse(4) default(shared)       
        for (auto i = 0; i < nMO; i++) 
        for (auto j = 0; j < nMO; j++) 
        for (auto k = 0; k < nMO; k++) 
        for (auto l = 0; l < nMO; l++) 
          N6MOERI(i, j, k, l) -= N5MOERI(i, j, k, l);
        
        N6TF.printOffSizes(N6TF.parseMOType(moType));
        N6MOERI.output(std::cout, "INCORE_N6 ERI - INCORE_N5 ERI", true);
      }
    }

#endif

    std::cout << "\n --------- End of the Test (on MO Ints Transformation)----- \n" << std::endl;
  }; // SingleSlater<MatsT>::MOIntsTransformationTest
  
#endif    
  
  /**
   *  \brief generate MOIntsTranformer using this singleslater as reference
   */
  template <typename MatsT, typename IntsT>
  std::shared_ptr<MOIntsTransformer<MatsT, IntsT>> 
    SingleSlater<MatsT, IntsT>::generateMOIntsTransformer() {
      return std::make_shared<MOIntsTransformer<MatsT, IntsT>>(*this, this->aoints_->TPITransAlg);
  }

/**
 *  \brief Reorthogonalize the MOs wrt overlap
 */
template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::orthoAOMO() {

    const size_t NB   = this->nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = this->nC * NB;
    const size_t NBC2 = NBC * NBC;

    // Compute MO orthogonalization on Root process then broadcast
    if( MPIRank(comm) == 0 ){
      if( this->nC > 1 or this->iCS ){
        size_t nOcc = this->nC == 1 ? nOA : nO;
        size_t disp = this->nC==4 ? NBC/2 : 0;  // displace to positive energy states
        orthoAB->orthogonalizeStates(this->mo[0], nOcc, disp);
      } else {
        orthoAB->orthogonalizeStates(this->mo[0], nOA);
        orthoAB->orthogonalizeStates(this->mo[1], nOB);
      }
    }

#ifdef CQ_ENABLE_MPI

  // Broadcast the updated MOs to all MPI processes
  if( MPISize(comm) > 1 ) {

    std::cerr << "  *** Scattering the AO-MOs ***\n";
    MPIBCast(this->mo[0].pointer(), nC * nC * NB * NB, 0, comm);
    if( nC == 1 and not iCS ) MPIBCast(this->mo[1].pointer(), nC * nC * NB * NB, 0, comm);
  }

#endif
}

template<typename MatsT, typename IntsT>
void SingleSlater<MatsT, IntsT>::setDenEqCoeff(bool val){
    this->denEqCoeff_ = val;
}
};   // namespace ChronusQ
