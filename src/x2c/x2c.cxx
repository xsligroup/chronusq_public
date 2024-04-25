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
#include <corehbuilder/x2c.hpp>
#include <corehbuilder/x2c/atomic.hpp>
#include <corehbuilder/matrixcoreh.hpp>
#include <corehbuilder/nonrel.hpp>
#include <corehbuilder/fourcomp.hpp>
#include <fockbuilder/matrixfock.hpp>
#include <fockbuilder/fourcompfock.hpp>
#include <particleintegrals/onepints/relativisticints.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <matrix.hpp>
#include <physcon.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/matfunc.hpp>

//#define DebugX2Cprint
//#define DebugX2Cprint2 
//#define oldimpl

namespace ChronusQ {

  template <typename MatsT>
  void GatherUSpin(const cqmatrix::Matrix<MatsT> &UL, const cqmatrix::Matrix<MatsT> &US, MatsT *U) {
    size_t NP = UL.dimension() / 2;
    SetMat('N', NP, 2*NP, 1.0, UL.pointer(), 2*NP, U, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, US.pointer(), 2*NP, U + NP, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, UL.pointer() + NP, 2*NP, U + 2*NP, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, US.pointer() + NP, 2*NP, U + 3*NP, 4*NP);
  }

  template void GatherUSpin(const cqmatrix::Matrix<double> &UL, const cqmatrix::Matrix<double> &US, double *U);
  template void GatherUSpin(const cqmatrix::Matrix<dcomplex> &UL, const cqmatrix::Matrix<dcomplex> &US, dcomplex *U);

  template <typename MatsT>
  void ReOrganizeMOSpin(const cqmatrix::Matrix<MatsT> &moSpin, cqmatrix::Matrix<MatsT> &mo) {

    size_t NP = moSpin.dimension() / 4;

    SetMat('N', NP, 4*NP,
           MatsT(1.), moSpin.pointer(), 4 * NP,
           mo.pointer(), 4 * NP);
    SetMat('N', NP, 4*NP,
           MatsT(1.), moSpin.pointer() + 2 * NP, 4 * NP,
           mo.pointer() + NP, 4 * NP);
    SetMat('N', NP, 4*NP,
           MatsT(1.), moSpin.pointer() + NP, 4 * NP,
           mo.pointer() + 2 * NP, 4 * NP);
    SetMat('N', NP, 4*NP,
           MatsT(1.), moSpin.pointer() + 3 * NP, 4 * NP,
           mo.pointer() + 3 * NP, 4 * NP);
  }

  template void ReOrganizeMOSpin(const cqmatrix::Matrix<double> &moSpin, cqmatrix::Matrix<double> &mo);
  template void ReOrganizeMOSpin(const cqmatrix::Matrix<dcomplex> &moSpin, cqmatrix::Matrix<dcomplex> &mo);



  /**
   *  \brief Boettger scaling for spin-orbit operator
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::SNSOScale(
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH,
      SNSO_TYPE snso_type) {

    if( this->basisSet_.maxL > 7 ) CErr("SNSO scaling for L > 7 NYI");

    size_t NB = basisSet_.nBasis;

    size_t n1, n2;
    std::array<double,8> Ql;
    switch (snso_type) {
      case SNSO_TYPE::BOETTGER:
        Ql={0.,2.,10.,28.,60.,110.,182.,280.};
        break;
      case SNSO_TYPE::DC:
        Ql={0.,2.32,10.64,28.38,60.,110.,182.,280.};
        break;
      case SNSO_TYPE::DCB:
        Ql={0.,2.97,11.93,29.84,64.,115.,188.,287.};
        break;
      default:
        CErr("Unknown row-independent SNSO type");
    }


    for(auto s1(0ul), i(0ul); s1 < this->basisSet_.nShell; s1++, i+=n1) {
      n1 = this->basisSet_.shells[s1].size();

      size_t L1 = this->basisSet_.shells[s1].contr[0].l;
      if ( L1 == 0 ) continue;

      auto Z1 = this->molecule_.atoms[this->basisSet_.mapSh2Cen[s1]].nucCharge;


    for(auto s2(0ul), j(0ul); s2 < this->basisSet_.nShell; s2++, j+=n2) {
      n2 = this->basisSet_.shells[s2].size();

      size_t L2 = this->basisSet_.shells[s2].contr[0].l;
      if ( L2 == 0 ) continue;

      auto Z2 = this->molecule_.atoms[this->basisSet_.mapSh2Cen[s2]].nucCharge;

      MatsT fudgeFactor = -1 * std::sqrt(
        Ql[L1] * Ql[L2] /
        Z1 / Z2
      );

      MatAdd('N','N',n1,n2,MatsT(1.),coreH->Z().pointer() + i + j*NB,NB,
          fudgeFactor,coreH->Z().pointer() + i + j*NB,NB,
          coreH->Z().pointer() + i + j*NB,NB);

      MatAdd('N','N',n1,n2,MatsT(1.),coreH->Y().pointer() + i + j*NB,NB,
          fudgeFactor,coreH->Y().pointer() + i + j*NB,NB,
          coreH->Y().pointer() + i + j*NB,NB);

      MatAdd('N','N',n1,n2,MatsT(1.),coreH->X().pointer() + i + j*NB,NB,
          fudgeFactor,coreH->X().pointer() + i + j*NB,NB,
          coreH->X().pointer() + i + j*NB,NB);

    } // loop s2
    } // loop s1
  }


  /**
   *  \brief Row-dependent DCB SNSO approximation for spin-orbit operator (ehrmanj row-dependent)
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::RowDepDCB_SNSO(
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

    size_t NB = basisSet_.nBasis;

    size_t n1, n2;
    // ehrmanj row-dependent Ql values
    std::array<double,8> Ql1;
    std::array<double,8> Ql2;
    // ehrmanj end

    if( this->basisSet_.maxL > 7 ) CErr("Boettger scaling for L > 7 NYI");

    for(auto s1(0ul), i(0ul); s1 < this->basisSet_.nShell; s1++, i+=n1) {
      n1 = this->basisSet_.shells[s1].size();

      size_t L1 = this->basisSet_.shells[s1].contr[0].l;
      if ( L1 == 0 ) continue;

      auto Z1 = this->molecule_.atoms[this->basisSet_.mapSh2Cen[s1]].nucCharge;


      for(auto s2(0ul), j(0ul); s2 < this->basisSet_.nShell; s2++, j+=n2) {
        n2 = this->basisSet_.shells[s2].size();

        size_t L2 = this->basisSet_.shells[s2].contr[0].l;
        if ( L2 == 0 ) continue;

        auto Z2 = this->molecule_.atoms[this->basisSet_.mapSh2Cen[s2]].nucCharge;
        // ehrmanj row-dependent Ql values
        if ( Z1 <= 2 ){
          Ql1={0.,2.97,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z2 <= 2 ){
          Ql2={0.,2.97,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z1 >= 3 && Z1 <= 10){
          Ql1={0.,2.80,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z2 >= 3 && Z2 <= 10){
          Ql2={0.,2.80,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z1 >= 11 && Z1 <= 18){
          Ql1={0.,2.95,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z2 >= 11 && Z2 <= 18){
          Ql2={0.,2.95,11.93,29.84,64.,115.,188.,287.};
        }
        if ( Z1 >= 19 && Z1 <= 36){
          Ql1={0.,3.09,11.49,29.84,64.,115.,188.,287.};
        }
        if ( Z2 >= 19 && Z2 <= 36){
          Ql2={0.,3.09,11.49,29.84,64.,115.,188.,287.};
        }
        if ( Z1 >= 37 && Z1 <= 54){
          Ql1={0.,3.02,11.91,29.84,64.,115.,188.,287.};
        }
        if ( Z2 >= 37 && Z2 <= 54){
          Ql2={0.,3.02,11.91,29.84,64.,115.,188.,287.};
        }
        if ( Z1 >= 55){
          Ql1={0.,2.85,12.31,30.61,64.,115.,188.,287.};
        }
        if ( Z2 >= 55){
          Ql2={0.,2.85,12.31,30.61,64.,115.,188.,287.};
        }//ehrmanj row-dependent Ql values

        MatsT fudgeFactor = -1 * std::sqrt(
            Ql1[L1] * Ql2[L2] /
            Z1 / Z2
        );
        std::cout<<Ql1[L1]<<std::endl;
        MatAdd('N','N',n1,n2,MatsT(1.),coreH->Z().pointer() + i + j*NB,NB,
               fudgeFactor,coreH->Z().pointer() + i + j*NB,NB,
               coreH->Z().pointer() + i + j*NB,NB);

        MatAdd('N','N',n1,n2,MatsT(1.),coreH->Y().pointer() + i + j*NB,NB,
               fudgeFactor,coreH->Y().pointer() + i + j*NB,NB,
               coreH->Y().pointer() + i + j*NB,NB);

        MatAdd('N','N',n1,n2,MatsT(1.),coreH->X().pointer() + i + j*NB,NB,
               fudgeFactor,coreH->X().pointer() + i + j*NB,NB,
               coreH->X().pointer() + i + j*NB,NB);

      } // loop s2
    } // loop s1
  }

  /**
   *  \brief Compute the X2C Core Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C(EMPerturbation &emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

#ifdef REAL_SPACE_X2C_ALGORITHM
    computeX2C_realSpace(emPert, coreH);
    return;
#endif

    IntsT* XXX = reinterpret_cast<IntsT*>(NULL);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    uncontractedInts_.computeAOOneP(
        molecule_, uncontractedBasis_, emPert,
        {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0}},
        ssOptions_.hamiltonianOptions);

    // Make copy of integrals
    IntsT *overlap   = CQMemManager::get().malloc<IntsT>(NP*NP);
    std::copy_n(uncontractedInts_.overlap->pointer(), NP*NP, overlap);

    // Compute the mappings from primitives to CGTOs
    mapPrim2Cont = CQMemManager::get().malloc<IntsT>(NP*NB);
    basisSet_.makeMapPrim2Cont(overlap,mapPrim2Cont);

    // Allocate Scratch Space (enough for 2*NP x 2*NP complex matricies)
    IntsT *SCR1  = CQMemManager::get().malloc<IntsT>(8*NP*NP);
    MatsT *CSCR1 = reinterpret_cast<MatsT*>(SCR1);

    // Singular value storage (initially S then T)
    p = CQMemManager::get().malloc<double>(NP);
    IntsT* SS = p;

    // Get SVD of uncontracted overlap
    // Store the left singular vectors in S
    nPrimUse_ = ORTH(NP,NP,overlap,NP,SS,XXX,NP);

    size_t NPU = nPrimUse_;

    // Form orthonormal transformation matrix in S
    for(auto i = 0ul; i < NPU; i++)
      blas::scal(NP,IntsT(1.)/std::sqrt(SS[i]),
          overlap + i*NP,1);

    // Transform T into the orthonormal basis
    // T -> TO
    std::shared_ptr<OnePInts<IntsT>> kinetic =
        std::dynamic_pointer_cast<OnePInts<IntsT>>(
            ParticleIntegrals::transform(
                *uncontractedInts_.kinetic, 'N', overlap, NPU, NP));

    // Get the SVD of TO
    // Store the left singular vectors in TO
    lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec, 
      NPU,NPU,kinetic->pointer(),NPU,SS,XXX,NPU,XXX,NPU);

    // Transformation matrix
    UK = CQMemManager::get().malloc<IntsT>(NP*NPU);
    std::fill_n(UK,NP*NPU,IntsT(0.0));

    // Form UK = S * T
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NP,NPU,NPU,IntsT(1.),overlap,NP,
      kinetic->pointer(),NPU,IntsT(0.),UK,NP);

    // Allocate and for "P^2" potential
    std::shared_ptr<OnePRelInts<IntsT>> potential =
        std::dynamic_pointer_cast<OnePRelInts<IntsT>>(
            ParticleIntegrals::transform(
                *uncontractedInts_.potential, 'N', UK, NPU, NP));

    // P^2 -> P^-1
    for(auto i = 0; i < NPU; i++) SS[i] = 1./std::sqrt(2*SS[i]);

    // Transform PVP into "P^-1" basis
    for (auto &oei : potential->SZYX())
      for(auto j = 0; j < NPU; j++)
        for(auto i = 0; i < NPU; i++)
          oei(i,j) *= SS[i] * SS[j];

    // Allocate 4C CORE Hamiltonian

    // CH = [ V    cp       ]
    //      [ cp   W - 2mc^2]
    MatsT *CH4C = CQMemManager::get().malloc<MatsT>(16*NPU*NPU);
    std::fill_n(CH4C,16*NPU*NPU,MatsT(0.));

    // Allocate W separately  as it's needed later
    size_t LDW = 2*NPU;
    cqmatrix::Matrix<MatsT> Wp(potential->template formW<MatsT>());

    // Subtract out 2mc^2 from W diagonals
    const double WFact = 2. * SpeedOfLight * SpeedOfLight;
    for(auto j = 0ul; j < 2*NPU; j++) Wp(j,j) -= WFact;

    // Copy W into the 4C CH storage
    MatsT *CHW = CH4C + 8*NPU*NPU + 2*NPU;
    SetMat('N',2*NPU,2*NPU,MatsT(1.),Wp.pointer(),LDW,CHW,4*NPU);

    // P^-1 -> P
    for(auto i = 0; i < NPU; i++) SS[i] = 1./SS[i];

    // V = [ V  0 ]
    //     [ 0  V ]
    MatsT * CHV = CH4C;
    SetMatDiag(NPU,NPU,potential->pointer(),NPU,CHV,4*NPU);

    // Set the diagonal cp blocks of CH
    // CP = [cp 0  ]
    //      [0  cp ]
    MatsT *CP11 = CH4C + 8*NPU*NPU;
    MatsT *CP12 = CP11 + 4*NPU*NPU + NPU;
    MatsT *CP21 = CH4C + 2*NPU;
    MatsT *CP22 = CP21 + 4*NPU*NPU + NPU;

    for(auto j = 0; j < NPU; j++) {
      CP11[j + 4*NPU*j] = SpeedOfLight * SS[j];
      CP12[j + 4*NPU*j] = SpeedOfLight * SS[j];
      CP21[j + 4*NPU*j] = SpeedOfLight * SS[j];
      CP22[j + 4*NPU*j] = SpeedOfLight * SS[j];
    }

    // Diagonalize the 4C CH
    double *CHEV = CQMemManager::get().malloc<double>(4*NPU);
    std::fill_n(CHEV, 4*NPU, 0.0);

    HermetianEigen('V','U',4*NPU,CH4C,4*NPU,CHEV);


    // Get pointers to "L" and "S" components of eigenvectors
    MatsT *L = CH4C + 8*NPU*NPU;
    MatsT *S = L + 2*NPU;


    // Invert "L"; L -> L^-1
    LUInv(2*NPU,L,4*NPU);


    // Reuse the charge conjugated space for X and Y
    X = std::make_shared<cqmatrix::Matrix<MatsT>>(2*NPU);
    X->clear();
    Y = std::make_shared<cqmatrix::Matrix<MatsT>>(2*NPU);
    Y->clear();

    // Form X = S * L^-1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,MatsT(1.),S,4*NPU,L,4*NPU,
               MatsT(0.),X->pointer(),X->dimension());

    // Form Y = sqrt(1 + X**H * X)

    // Y = X**H * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,
               MatsT(1.),X->pointer(),X->dimension(),
               X->pointer(),X->dimension(),
               MatsT(0.),Y->pointer(),Y->dimension());

    // Y = Y + I
    for(auto j = 0; j < 2*NPU; j++) (*Y)(j,j) += 1.0;

    // Y = Y^-0.5
    MatDiagFunc(std::function<double(double)>([](double x){ return std::pow(x, -0.5); }),
                2*NPU, Y->pointer(), Y->dimension(), Y->pointer(), Y->dimension());

    // Build the effective two component CH in "L"
    cqmatrix::Matrix<MatsT> FullCH2C(2*NPU);

    // Copy potential into spin diagonal blocks of 2C CH
    SetMatDiag(NPU,NPU,potential->pointer(),NPU,FullCH2C.pointer(),2*NPU);

    // Construct 2C CH in the uncontracted basis
    // 2C CH = Y * (V' + cp * X + X**H * cp + X**H * W' * X) * Y

    // SCR1 = cp * X
    for(auto j = 0; j < 2*NPU; j++)
    for(auto i = 0; i < NPU; i++) {
      CSCR1[i + 2*NPU*j] = SpeedOfLight * SS[i] * (*X)(i,j);
      CSCR1[i + NPU + 2*NPU*j] = SpeedOfLight * SS[i] * (*X)(i + NPU, j);
    }

    // 2C CH += SCR1 + SCR1**H
    MatAdd('N','N',2*NPU,2*NPU,MatsT(1.),FullCH2C.pointer(),2*NPU,MatsT(1.),
      CSCR1,2*NPU, FullCH2C.pointer(),2*NPU);
    MatAdd('N','C',2*NPU,2*NPU,MatsT(1.),FullCH2C.pointer(),2*NPU,MatsT(1.),
      CSCR1,2*NPU, FullCH2C.pointer(),2*NPU);


    // SCR1 = X**H * W
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,MatsT(1.),X->pointer(),X->dimension(),
               Wp.pointer(),LDW,MatsT(0.),CSCR1,2*NPU);

    // 2C CH += SCR1 * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,MatsT(1.),CSCR1,2*NPU,
               X->pointer(),X->dimension(),MatsT(1.),FullCH2C.pointer(),2*NPU);

    // SCR1 = CH2C * Y
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,MatsT(1.),FullCH2C.pointer(),2*NPU,
               Y->pointer(),Y->dimension(),MatsT(0.),CSCR1,2*NPU);


    // 2C CH = Y * SCR1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NPU,2*NPU,2*NPU,MatsT(1.),Y->pointer(),Y->dimension(),CSCR1,2*NPU,
               MatsT(0.),FullCH2C.pointer(),2*NPU);

    // Allocate memory for the uncontracted spin components
    // of the 2C CH
    cqmatrix::PauliSpinorMatrices<MatsT> HUn(
        FullCH2C.template spinScatter<MatsT>(
            ssOptions_.hamiltonianOptions.OneESpinOrbit,ssOptions_.hamiltonianOptions.OneESpinOrbit));

    // Partition the scratch space into one complex and one real NP x NP
    // matrix
    IntsT *SUK   = SCR1;
    IntsT *CPSUK = SUK + NP*NPU;

    // Store the Product of S and UK
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NP,NPU,NP,IntsT(1.),uncontractedInts_.overlap->pointer(),NP,
         UK,NP,IntsT(0.),SUK,NP);
    // Store the Product of mapPrim2Cont and SUK
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NPU,NP,IntsT(1.),mapPrim2Cont,NB,
         SUK,NP,IntsT(0.),CPSUK,NB);

    // Transform the spin components of the 2C CH into R-space
    *coreH = HUn.transform('C', CPSUK, NB, NB);

    CQMemManager::get().free(overlap, SCR1, CH4C, CHEV);


  }

  template void X2C<dcomplex,double>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>);

  /**
   *  \brief X2C core Hamiltonian with GIAO.
   */
  template <>
  void X2C<dcomplex, dcomplex>::computeOneEX2C(EMPerturbation &emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreH) {

#ifdef REAL_SPACE_X2C_ALGORITHM
    CErr("No Real Space for GIAO!",std::cout);
#endif

    dcomplex* XXX = reinterpret_cast<dcomplex*>(NULL);

    auto magAmp = emPert.getDipoleAmp(Magnetic);
    dcomplex onei = dcomplex(0,1);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    uncontractedInts_.computeAOOneP(
        molecule_, uncontractedBasis_, emPert,
        {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0}, 
         {LEN_ELECTRIC_MULTIPOLE,2}, {MAGNETIC_MULTIPOLE,1},
         {MAGNETIC_4COMP_rVr,2},
         {MAGNETIC_4COMP_PVrprVP,2},
         {MAGNETIC_4COMP_PVrmrVP,2}},
        ssOptions_.hamiltonianOptions);

    std::cout << " Building GIAO core Hamiltonian." << std::endl;
    
    // Make copy of integrals
    dcomplex *overlap   = CQMemManager::get().malloc<dcomplex>(NP*NP);
    std::copy_n(uncontractedInts_.overlap->pointer(), NP*NP, overlap);

    // Compute the mappings from primitives to GIAOs
    mapPrim2Cont = CQMemManager::get().malloc<dcomplex>(NP*NB);
    basisSet_.makeMapPrim2Cont(overlap,mapPrim2Cont);

    // Allocate Scratch Space (enough for 2*NP x 2*NP complex matricies)
    dcomplex *SCR1  = CQMemManager::get().malloc<dcomplex>(8*NP*NP);
    std::fill_n(SCR1,8*NP*NP,dcomplex(0.));
    dcomplex *CSCR1 = reinterpret_cast<dcomplex*>(SCR1);

    // Make a copy of the overlap for later
    dcomplex* SCPY = CQMemManager::get().malloc<dcomplex>(4*NP*NP);
    std::fill_n(SCPY,4*NP*NP,dcomplex(0.));
    dcomplex* M = CQMemManager::get().malloc<dcomplex>(4*NP*NP);
    dcomplex* VCPY = CQMemManager::get().malloc<dcomplex>(4*NP*NP);
    std::fill_n(VCPY,4*NP*NP,dcomplex(0.));
    dcomplex* Ms = CQMemManager::get().malloc<dcomplex>(NP*NP);  // Scalar component of M 

    // Construct M Matrix
    // M = 1/2 * (\sigma\dot\pi) * (\sigma\dot\pi)

    // assemble the scalar part of M (Ms)
    for(auto k = 0ul; k < NP*NP; k++)
      Ms[k] = (uncontractedInts_.kinetic->pointer()[k]);// + this->aoints.potential[k]);

    // this part add the angular momentum term 
    for ( auto index = 0 ; index < 3 ; index++ ) {
      MatAdd('N','N',NP,NP,-0.5*magAmp[index]*onei,
        (*uncontractedInts_.magnetic)[index]->pointer(),NP,dcomplex(1.),Ms,NP,Ms,NP);
    } // for ( auto inde = 0 ; inde < 3 ; inde++ ) 

    // this part add the length gauge electric quadrupole term (diamagnetic term)
    const std::array<std::string,3> diagindex = { "XX","YY","ZZ" };

    double diagcoeff[3];
    diagcoeff[0] = 1.0/8.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]); 
    diagcoeff[1] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);    
    diagcoeff[2] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);    

    // add diagonal part
    for ( auto index = 0 ; index < 3 ; index++ ) { 
      MatAdd('N','N',NP,NP, 
        dcomplex(diagcoeff[index]),
        (*uncontractedInts_.lenElectric)[diagindex[index]]->pointer(),
        NP,dcomplex(1.),Ms,NP,Ms,NP);
    }   

    const std::array<std::string,3> offindex = { "XY","XZ","YZ" };
    
    double offcoeff[3];
    offcoeff[0] = -1.0/4.0*magAmp[0]*magAmp[1];
    offcoeff[1] = -1.0/4.0*magAmp[0]*magAmp[2];
    offcoeff[2] = -1.0/4.0*magAmp[1]*magAmp[2];
   
    // add off diagonal part
    for ( auto index = 0 ; index < 3 ; index++ ) { 
      MatAdd('N','N',NP,NP, 
        dcomplex(offcoeff[index]),
        (*uncontractedInts_.lenElectric)[offindex[index]]->pointer(),
        NP,dcomplex(1.),Ms,NP,Ms,NP);
    } 

    // M = [ M1  M2 ]
    //     [ M3  M4 ]
    dcomplex *M1 = M;
    dcomplex *M2 = M1 + 2*NP*NP;
    dcomplex *M3 = M1 + NP;
    dcomplex *M4 = M2 + NP;

    // M1 = Ms + Mz
    MatAdd('N','N',NP,NP,dcomplex(1.),Ms,NP,0.5*dcomplex(magAmp[2]),
      overlap,NP,M1,2*NP);
    // M4 = Ms - Mz
    MatAdd('N','N',NP,NP,dcomplex(1.),Ms,NP,0.5*dcomplex(-magAmp[2]),
      overlap,NP,M4,2*NP);

    // M2 = Mx - iMy
    MatAdd('N','N',NP,NP,0.5*dcomplex(magAmp[0]),overlap,NP,
      0.5*dcomplex(0.,-1.)*magAmp[1],overlap,NP,M2,2*NP);
    // M3 = Mx + iMy
    MatAdd('N','N',NP,NP,0.5*dcomplex(magAmp[0]),overlap,NP,
      0.5*dcomplex(0.,1.)*magAmp[1],overlap,NP,M3,2*NP);

#ifdef DebugX2Cprint  
    prettyPrintSmart(std::cout,"2c M",M,2*NP,2*NP,2*NP);
#endif 
    
    // need to copy a block instead of the whole matrix
    SetMat('N',NP,NP,dcomplex(1.),overlap,NP,SCPY,2*NP);
    SetMat('N',NP,NP,dcomplex(1.),overlap,NP,SCPY+2*NP*NP+NP,2*NP);
    SetMat('N',NP,NP,dcomplex(1.),uncontractedInts_.potential->pointer(),NP,VCPY,2*NP);
    SetMat('N',NP,NP,dcomplex(1.),uncontractedInts_.potential->pointer(),NP,VCPY+2*NP*NP+NP,2*NP);

    // Singular value storage (initially S then T)
    p = CQMemManager::get().malloc<double>(2*NP);
    std::fill_n(p, 2*NP, 0.0);
    double* SS = p;
    
    // Get SVD of uncontracted overlap
    // Store the left singular vectors in S
    lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec, 
      2*NP,2*NP,SCPY,2*NP,SS,XXX,2*NP,XXX,2*NP);
    double minSS = *std::min_element(SS,SS+2*NP);
    if( minSS < 1e-10 ) CErr("Stop: Uncontracted Overlap is Singular");

    // FIXME: Reducing linear dependency -- something like
    // nPrimUse_ = ORTH(2*NP,2*NP,SCPY,2*NP,SS,XXX,2*NP);
    // size_t NPU = nPrimUse_;
    // if (nPrimUse_ < 2*NP)
    //    std::cout << "Info: Uncontracted Overlap is Singular." << std::endl;  

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"svd of 2c overlap",SCPY,2*NP,2*NP,2*NP);
#endif 

    // Form orthonormal transformation matrix in S
    for(auto i = 0ul; i < 2*NP; i++)
      blas::scal(2*NP,dcomplex(1.)/std::sqrt(SS[i]),
        SCPY + i*2*NP,1);

#ifdef DebugX2Cprint  
    prettyPrintSmart(std::cout,"ortho matrix",SCPY,2*NP,2*NP,2*NP);
    prettyPrintSmart(std::cout,"overlap singular values",SS,2*NP,1,2*NP);
#endif 

    // Transform M into the orthonormal basis (TangDD: M == T if no Magnetic)
    // M -> MO
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),SCPY,2*NP,
      M,2*NP,dcomplex(0.),SCR1,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),SCR1,2*NP,
      SCPY,2*NP,dcomplex(0.),M,2*NP);

#ifdef DebugX2Cprint  
    prettyPrintSmart(std::cout,"M in orthogonal basis",M,2*NP,2*NP,2*NP);
#endif

    // Get the SVD of MO
    // Store the left singular vectors in MO
    lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec, 
      2*NP,2*NP,M,2*NP,SS,XXX,2*NP,XXX,2*NP);

    minSS = *std::min_element(SS,SS+2*NP);
    if( minSS < 1e-10 ) CErr("Uncontracted Kinetic Energy Tensor is Singular");

#ifdef DebugX2Cprint       
    prettyPrintSmart(std::cout,"svd of ortho M",M,2*NP,2*NP,2*NP);
    prettyPrintSmart(std::cout,"M singular values",SS,2*NP,1,2*NP);
#endif

    // Transformation matrix
    UK = CQMemManager::get().malloc<dcomplex>(2*NP*2*NP);
    std::fill_n(UK,2*NP*2*NP,dcomplex(0.0));

    // Form UK = US (Stored in S) * UT (Stored in M)
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),SCPY,2*NP,
      M,2*NP,dcomplex(0.),UK,2*NP);

#ifdef DebugX2Cprint       
    prettyPrintSmart(std::cout,"UK",UK,2*NP,2*NP,2*NP);
#endif

    // Allocate and for "P^2" potential
    dcomplex *P2P = CQMemManager::get().malloc<dcomplex>(2*NP*2*NP);
    std::fill_n(P2P,2*NP*2*NP,dcomplex(0.0));

    // P2P = UK**H * V * UK  -- Potential in P^2 basis
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),UK,2*NP,
      VCPY,2*NP,dcomplex(0.),SCR1,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),SCR1,2*NP,
      UK,2*NP,dcomplex(0.),P2P,2*NP);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"P2P",P2P,2*NP,2*NP,2*NP);
#endif

    // Allocate W separately as it's needed later
    size_t LDW = 2*NP;
    W = std::make_shared<cqmatrix::Matrix<dcomplex>>(
        std::dynamic_pointer_cast<OnePRelInts<dcomplex>>(
            uncontractedInts_.potential)->template formW<dcomplex>());
    auto Wp = W->pointer();

#ifdef DebugX2Cprint       
    prettyPrintSmart(std::cout,"W only pVp",Wp,2*NP,2*NP,2*NP);
#endif

    // need to copy a block instead of the whole matrix
    // Components of 1/(4c^2)*(\sigma\dot)V(\sigma\dot)
    dcomplex *AVA     = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *DOT     = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *CROSSx  = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *CROSSy  = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *CROSSz  = CQMemManager::get().malloc<dcomplex>(NP*NP);

    // calculate AVA part
    memset(AVA,0,NP*NP*sizeof(dcomplex));
    
    diagcoeff[0] = 1.0/4.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]); 
    diagcoeff[1] = 1.0/4.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);    
    diagcoeff[2] = 1.0/4.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);    

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"rVr XX",(*uncontractedInts_.rVr)["XX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"rVr YY",(*uncontractedInts_.rVr)["YY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"rVr ZZ",(*uncontractedInts_.rVr)["ZZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"rVr XY",(*uncontractedInts_.rVr)["XY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"rVr XZ",(*uncontractedInts_.rVr)["XZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"rVr YZ",(*uncontractedInts_.rVr)["YZ"]->pointer(),NP,NP,NP);
#endif

    // add diagonal part:  xVx(By^2+Bz^2)...
    for ( auto index = 0 ; index < 3 ; index++ ) { 
      MatAdd('N','N',NP,NP, 
        dcomplex(diagcoeff[index]),
        (*uncontractedInts_.rVr)[diagindex[index]]->pointer(),
        NP,dcomplex(1.),AVA,NP,AVA,NP);
    }   
    
    offcoeff[0] = -1.0/2.0*magAmp[0]*magAmp[1];
    offcoeff[1] = -1.0/2.0*magAmp[0]*magAmp[2];
    offcoeff[2] = -1.0/2.0*magAmp[1]*magAmp[2];
   
    // add off diagonal part:  -xVy2BxBy...
    for ( auto index = 0 ; index < 3 ; index++ ) { 
      MatAdd('N','N',NP,NP, 
        dcomplex(offcoeff[index]),
        (*uncontractedInts_.rVr)[offindex[index]]->pointer(),
        NP,dcomplex(1.),AVA,NP,AVA,NP);
    }   

    // AVA part end

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"PVrprVP XX",(*uncontractedInts_.PVrprVP)["XX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP XY",(*uncontractedInts_.PVrprVP)["XY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP XZ",(*uncontractedInts_.PVrprVP)["XZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP YX",(*uncontractedInts_.PVrprVP)["YX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP YY",(*uncontractedInts_.PVrprVP)["YY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP YZ",(*uncontractedInts_.PVrprVP)["YZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP ZX",(*uncontractedInts_.PVrprVP)["ZX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP ZY",(*uncontractedInts_.PVrprVP)["ZY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrprVP ZZ",(*uncontractedInts_.PVrprVP)["ZZ"]->pointer(),NP,NP,NP);

    prettyPrintSmart(std::cout,"PVrmrVP XX",(*uncontractedInts_.PVrmrVP)["XX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP XY",(*uncontractedInts_.PVrmrVP)["XY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP XZ",(*uncontractedInts_.PVrmrVP)["XZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP YX",(*uncontractedInts_.PVrmrVP)["YX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP YY",(*uncontractedInts_.PVrmrVP)["YY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP YZ",(*uncontractedInts_.PVrmrVP)["YZ"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP ZX",(*uncontractedInts_.PVrmrVP)["ZX"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP ZY",(*uncontractedInts_.PVrmrVP)["ZY"]->pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"PVrmrVP ZZ",(*uncontractedInts_.PVrmrVP)["ZZ"]->pointer(),NP,NP,NP);
#endif

    // CROSS part start (pVxA + AVxp)

    // notice that PVrprVP(alpha,beta) = PVrprVP(beta,alpha)
    //             PVrmrVP(alpha,beta) =-PVrmrVP(beta,alpha)

    // 0: xx       1: xy      2: xz 
    // 3: yx       4: yy      5: yz
    // 6: zx       7: zy      8: zz   

    memset(CROSSx,0,NP*NP*sizeof(dcomplex));  
    for ( auto index = 0 ; index < 3 ; index++ ) { 
      MatAdd('N','N',NP,NP, 
        dcomplex(1.0),
        (*uncontractedInts_.PVrmrVP)[diagindex[index]]->pointer(),
        NP,dcomplex(1.),CROSSx,NP,CROSSx,NP);
    }   

    memset(CROSSy,0,NP*NP*sizeof(dcomplex));  
    memset(CROSSz,0,NP*NP*sizeof(dcomplex));  
    
    SetMat('N',NP,NP,dcomplex(1.),CROSSx,NP,CROSSy,NP);
    SetMat('N',NP,NP,dcomplex(1.),CROSSx,NP,CROSSz,NP);

    blas::scal(NP*NP,dcomplex(0.5*magAmp[0]),CROSSx,1);
    blas::scal(NP*NP,dcomplex(0.5*magAmp[1]),CROSSy,1);
    blas::scal(NP*NP,dcomplex(0.5*magAmp[2]),CROSSz,1);

    // CROSSx = CROSSx -0.5*Bx*(pxVrx-rxVpx)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[0]),
      (*uncontractedInts_.PVrmrVP)["XX"]->pointer(),
      NP,dcomplex(1.),CROSSx,NP,CROSSx, NP);

    // CROSSx = CROSSx -0.5*By*(pyVrx-rxVpy)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[1]),
      (*uncontractedInts_.PVrmrVP)["YX"]->pointer(),
      NP,dcomplex(1.),CROSSx,NP,CROSSx, NP);

    // CROSSx = CROSSx -0.5*Bz*(pzVrx-rxVpz)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[2]),
      (*uncontractedInts_.PVrmrVP)["ZX"]->pointer(),
      NP,dcomplex(1.),CROSSx,NP,CROSSx, NP);


    // CROSSy = CROSSy -0.5*Bx*(pxVry-ryVpx)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[0]),
      (*uncontractedInts_.PVrmrVP)["XY"]->pointer(),
      NP,dcomplex(1.),CROSSy,NP,CROSSy, NP);

    // CROSSy = CROSSy -0.5*By*(pyVry-ryVpy)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[1]),
      (*uncontractedInts_.PVrmrVP)["YY"]->pointer(),
      NP,dcomplex(1.),CROSSy,NP,CROSSy, NP);

    // CROSSy = CROSSy -0.5*Bz*(pzVry-ryVpz)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[2]),
      (*uncontractedInts_.PVrmrVP)["ZY"]->pointer(),
      NP,dcomplex(1.),CROSSy,NP,CROSSy, NP);


    // CROSSz = CROSSz -0.5*Bx*(pxVrz-rzVpx)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[0]),
      (*uncontractedInts_.PVrmrVP)["XZ"]->pointer(),
      NP,dcomplex(1.),CROSSz,NP,CROSSz, NP);

    // CROSSz = CROSSz -0.5*By*(pyVrz-rzVpy)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[1]),
      (*uncontractedInts_.PVrmrVP)["YZ"]->pointer(),
      NP,dcomplex(1.),CROSSz,NP,CROSSz, NP);

    // CROSSz = CROSSz -0.5*Bz*(pzVrz-rzVpz)        
    MatAdd('N','N',NP,NP, 
      dcomplex(-0.5*magAmp[2]),
      (*uncontractedInts_.PVrmrVP)["ZZ"]->pointer(),
      NP,dcomplex(1.),CROSSz,NP,CROSSz, NP);

    blas::scal(NP*NP,dcomplex(0.0,-1.0),CROSSx,1);
    blas::scal(NP*NP,dcomplex(0.0,-1.0),CROSSy,1);
    blas::scal(NP*NP,dcomplex(0.0,-1.0),CROSSz,1);
    
    // CROSS part end 

    // now calculate DOT (pV.A + AV.p)

    memset(DOT,0,NP*NP*sizeof(dcomplex));  

    // DOT = DOT + Bx*[(DyVrz+rzVDy)-(DzVry+ryVDz)        
    MatAdd('N','N',NP,NP, 
      dcomplex(magAmp[0]),
      (*uncontractedInts_.PVrprVP)["YZ"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);

    MatAdd('N','N',NP,NP, 
      dcomplex(-magAmp[0]),
      (*uncontractedInts_.PVrprVP)["ZY"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);


    // DOT = DOT + By*[(DzVrx+rxVDz)-(DxVrz+rzVDx)        
    MatAdd('N','N',NP,NP, 
      dcomplex(magAmp[1]),
      (*uncontractedInts_.PVrprVP)["ZX"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);

    MatAdd('N','N',NP,NP, 
      dcomplex(-magAmp[1]),
      (*uncontractedInts_.PVrprVP)["XZ"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);


    // DOT = DOT + Bz*[(DxVry+ryVDx)-(DyVrx+rxVDy)        
    MatAdd('N','N',NP,NP, 
      dcomplex(magAmp[2]),
      (*uncontractedInts_.PVrprVP)["XY"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);

    MatAdd('N','N',NP,NP, 
      dcomplex(-magAmp[2]),
      (*uncontractedInts_.PVrprVP)["YX"]->pointer(),
      NP,dcomplex(1.),DOT,NP,DOT, NP);

    blas::scal(NP*NP,dcomplex(0.0,0.5),DOT,1);
   
    // DOT end   

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"AVA",AVA,NP,NP,NP);
    prettyPrintSmart(std::cout,"DOT",DOT,NP,NP,NP);
    prettyPrintSmart(std::cout,"CROSSx one term",CROSSx,NP,NP,NP);
    prettyPrintSmart(std::cout,"CROSSy one term",CROSSy,NP,NP,NP);
    prettyPrintSmart(std::cout,"CROSSz one term",CROSSz,NP,NP,NP);
#endif

    // add DOT and AVA to the diagonal blocks 

    MatAdd('N','N',NP,NP,dcomplex(1.),DOT,NP,dcomplex(1.),Wp,2*NP,Wp,2*NP);
    MatAdd('N','N',NP,NP,dcomplex(1.),AVA,NP,dcomplex(1.),Wp,2*NP,Wp,2*NP);
    MatAdd('N','N',NP,NP,dcomplex(1.),DOT,NP,dcomplex(1.),Wp+2*NP*NP+NP,2*NP,Wp+2*NP*NP+NP,2*NP);
    MatAdd('N','N',NP,NP,dcomplex(1.),AVA,NP,dcomplex(1.),Wp+2*NP*NP+NP,2*NP,Wp+2*NP*NP+NP,2*NP);

    // add CROSS x,y,z component to W 
    // Waa +=  i CROSSz

    MatAdd('N','N',NP,NP,dcomplex(1.),Wp,2*NP,dcomplex(0.,1.),CROSSz,NP,Wp,2*NP);

    // Wbb += -i CROSSz

    MatAdd('N','N',NP,NP,dcomplex(1.),Wp+2*NP*NP+NP,2*NP,dcomplex(0.,-1.),CROSSz,NP,
      Wp+2*NP*NP+NP,2*NP);

    // Wab +=  CROSSy + i CROSSx

    MatAdd('N','N',NP,NP,dcomplex(1.),Wp+2*NP*NP,2*NP,dcomplex(1.0),CROSSy,NP,
      Wp+2*NP*NP,2*NP);
    MatAdd('N','N',NP,NP,dcomplex(1.),Wp+2*NP*NP,2*NP,dcomplex(0.0,1.0),CROSSx,NP,
      Wp+2*NP*NP,2*NP);

    // Wba += -CROSSy + i CROSSx

    MatAdd('N','N',NP,NP,dcomplex(1.),Wp+NP,2*NP,dcomplex(-1.0),CROSSy,NP,
      Wp+NP,2*NP);
    MatAdd('N','N',NP,NP,dcomplex(1.),Wp+NP,2*NP,dcomplex(0.0,1.0),CROSSx,NP,
      Wp+NP,2*NP);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"W end",Wp,2*NP,2*NP,2*NP);
#endif

    // do P^2 transform for W
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),UK,2*NP,
      Wp,2*NP,dcomplex(0.),SCR1,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),SCR1,2*NP,
      UK,2*NP,dcomplex(0.),Wp,2*NP);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"W p2",Wp,2*NP,2*NP,2*NP);
#endif

    // P^2 -> P^-1
    for(auto i = 0; i < 2*NP; i++) SS[i] = 1./std::sqrt(2*SS[i]);

    // Transform W into "P^-1" basis
    for(auto j = 0; j < 2*NP; j++) 
    for(auto i = 0; i < 2*NP; i++){
      Wp[i+j*2*NP] *= SS[i] * SS[j];
    }

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"W p1",Wp,2*NP,2*NP,2*NP);
#endif

    // Subtract out 2mc^2 from W diagonals
    const dcomplex WFact = 2. * SpeedOfLight * SpeedOfLight;
    for(auto j = 0ul; j < 2*NP; j++) Wp[j + LDW*j] -= WFact;

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"W final",Wp,2*NP,2*NP,2*NP);
#endif
   
    // W End

    // Allocate 4C CORE Hamiltonian

    // CH = [ V    cp       ]
    //      [ cp   W - 2mc^2]
    dcomplex *CH4C = CQMemManager::get().malloc<dcomplex>(16*NP*NP);
    memset(CH4C,0,16*NP*NP*sizeof(dcomplex));

    // Copy W into the 4C CH storage
    dcomplex *CHW = CH4C + 8*NP*NP + 2*NP;
    SetMat('N',2*NP,2*NP,dcomplex(1.),Wp,LDW,CHW,4*NP);

    // P^-1 -> P
    for(auto i = 0; i < 2*NP; i++) SS[i] = 1./SS[i];

    // V = [ P2P  0   ]
    //     [ 0    P2P ]
    dcomplex * V1 = CH4C;
    SetMat('N',2*NP,2*NP,dcomplex(1.),P2P,2*NP,V1,4*NP);

    // Set the diagonal cp blocks of CH
    // CP = [ cp  0 ]
    //      [ 0  cp ]
    dcomplex *CP11 = CH4C + 8*NP*NP;
    dcomplex *CP21 = CH4C + 2*NP;
   
    for(auto j = 0; j < 2*NP; j++) {
      CP11[j + 4*NP*j] = SpeedOfLight * SS[j];
      CP21[j + 4*NP*j] = SpeedOfLight * SS[j];
    }

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"4CCH",CH4C,4*NP,4*NP,4*NP);
#endif

    // Diagonalize the 4C CH
    double *CHEV = CQMemManager::get().malloc<double>(4*NP);
    HermetianEigen('V','U',4*NP,CH4C,4*NP,CHEV);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"4C eigen values",CHEV,4*NP,1,4*NP);
#endif

#ifdef oldimpl
    // Old Implementations 
    // Only for debug use!

    // Get pointers to "L" and "S" components of eigenvectors
    dcomplex *L = CH4C + 8*NP*NP;
    dcomplex *S = L + 2*NP;

    // Invert "L"; L -> L^-1
    LUInv(2*NP,L,4*NP);

    // Reuse the charge conjugated space for X and Y
    dcomplex *X = CH4C;
    dcomplex *Y = X + 2*NP;

    // Form X = S * L^-1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),S,4*NP,
      L,4*NP,dcomplex(0.),X,4*NP);

    // Form Y = sqrt(1 + X**H * X)

    // Y = X**H * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),X,4*NP,
      X,4*NP,dcomplex(0.),Y,4*NP);

    // Y = Y + I
    for(auto j = 0; j < 2*NP; j++) Y[j + 4*NP*j] += 1.0;

    // Y -> V * y * V**H 
    // XXX: Store the eigenvalues of Y in CHEV
    HermetianEigen('V','U',2*NP,Y,4*NP,CHEV);

    // SCR1 -> V * y^-0.25
    for(auto j = 0ul; j < 2*NP; j++)
    for(auto i = 0ul; i < 2*NP; i++)
      CSCR1[i + 2*NP*j] = Y[i + 4*NP*j] * std::pow(CHEV[j],-0.25);

    // Y = SCR1 * SCR1**H
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,
      CSCR1,2*NP,dcomplex(0.),Y,4*NP);
    
    // Build the effective two component CH in "L"
    dcomplex *FullCH2C = L;

    // Zero it out
    for(auto j = 0; j < 2*NP; j++)
    for(auto i = 0; i < 2*NP; i++)
      FullCH2C[i + 4*NP*j] = 0.;

    // Copy P2P into spin diagonal blocks of 2C CH
    dcomplex *CH2C1 = FullCH2C;
    dcomplex *CH2C2 = CH2C1 + 4*NP*NP + NP;

    SetMat('N',2*NP,2*NP,dcomplex(1.),P2P,2*NP,CH2C1,4*NP);

    // Construct 2C CH in the uncontracted basis
    // 2C CH = Y * (V' + cp * X + X**H * cp + X**H * W' * X) * Y

    // SCR1 = cp * X
    for(auto j = 0; j < 2*NP; j++)
    for(auto i = 0; i < 2*NP; i++) {
      CSCR1[i + 2*NP*j] = SpeedOfLight * SS[i] * X[i + 4*NP*j];
    }

    // 2C CH += SCR1 + SCR1**H
    MatAdd('N','N',2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C,4*NP);
    MatAdd('N','C',2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C,4*NP);

    // SCR1 = X**H * W
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),X,4*NP,
      Wp,LDW,dcomplex(0.),CSCR1,2*NP);
    
    // 2C CH += SCR1 * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,
      X,4*NP,dcomplex(1.),FullCH2C,4*NP);

    // SCR1 = CH2C * Y
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),FullCH2C,4*NP,
      Y,4*NP,dcomplex(0.),CSCR1,2*NP);

    // 2C CH = Y * SCR1
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),Y,4*NP,
      CSCR1,2*NP,dcomplex(0.),FullCH2C,4*NP);

    // Allocate memory for the uncontracted spin components of the 2C CH
    dcomplex *HUnS = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *HUnZ = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *HUnX = CQMemManager::get().malloc<dcomplex>(NP*NP);
    dcomplex *HUnY = CQMemManager::get().malloc<dcomplex>(NP*NP);

    // first half of the scratch space is SUK(2*NP x 2*NP) 
    // matrix
    dcomplex   * SUK   = SCR1;
    dcomplex * CSCR2 = SUK+4*NP*NP;

    // recopy overlap matrix 2c 
    memset(SCPY,0,4*NP*NP*sizeof(dcomplex));
    SetMat('N',NP,NP,dcomplex(1.),uncontractedInts_.overlap->pointer(),NP,SCPY,2*NP);
    SetMat('N',NP,NP,dcomplex(1.),uncontractedInts_.overlap->pointer(),NP,SCPY+2*NP*NP+NP,2*NP);

    LUInv(2*NP,UK,2*NP);

    // Transform the spin components of the 2C CH into R-space
    //
    // H(k) -> SUK * H(k) * (SUK)**H
    //
    // ** Using the fact that H(k) is hermetian
    // CSCR2 = SUK * H(k) -> CSCR2**H = H(k) * (SUK)**H
    // H(k) -> SUK * CSCR2**H
    //

    // SCR2 = SUK**H * CH2C 
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),UK,2*NP,
      FullCH2C,4*NP,dcomplex(0.),CSCR2,2*NP);

    // 2C CH = SCR2 * SUK
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),CSCR2,2*NP,
      UK,2*NP,dcomplex(0.),FullCH2C,4*NP);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"FullCH2C",FullCH2C,2*NP,2*NP,4*NP);
#endif

    SpinScatter(NP,FullCH2C,4*NP,HUnS,NP,HUnZ,NP,HUnY,NP,HUnX,NP);

    // Transform H(k) into the contracted basis

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnS,
      NP,dcomplex(0.),CSCR1,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnS,NB);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnZ,
      NP,dcomplex(0.),CSCR1,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnZ,NB);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnY,
      NP,dcomplex(0.),CSCR1,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnY,NB);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,HUnX,
      NP,dcomplex(0.),CSCR1,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NP,dcomplex(1.),mapPrim2Cont,NB,CSCR1,
      NB,dcomplex(0.),HUnX,NB);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"CH0 no Boettger",HUnS,NB,NB,NB);
    prettyPrintSmart(std::cout,"CHx no Boettger",HUnX,NB,NB,NB);
    prettyPrintSmart(std::cout,"CHy no Boettger",HUnY,NB,NB,NB);
    prettyPrintSmart(std::cout,"CHz no Boettger",HUnZ,NB,NB,NB);
#endif
      
#else
    // Get pointers to "L" and "S" components of eigenvectors
    dcomplex *L = CH4C + 8*NP*NP;
    dcomplex *S = L + 2*NP;

    // Invert "L"; L -> L^-1
    LUInv(2*NP,L,4*NP);

    // Save X and Y
    X = std::make_shared<cqmatrix::Matrix<dcomplex>>(2*NP);
    X->clear();
    Y = std::make_shared<cqmatrix::Matrix<dcomplex>>(2*NP);
    Y->clear();

    // Form X = S * L^-1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),S,4*NP,
      L,4*NP,dcomplex(0.),X->pointer(),X->dimension());

    // Form Y = sqrt(1 + X**H * X)

    // Y = X**H * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),X->pointer(),X->dimension(),
      X->pointer(),X->dimension(),dcomplex(0.),Y->pointer(),Y->dimension());

    // Y = Y + I
    for(auto j = 0; j < 2*NP; j++) (*Y)(j,j) += 1.0;

    // Y -> V * y * V**H 
    // XXX: Store the eigenvalues of Y in CHEV
    HermetianEigen('V','U',2*NP,Y->pointer(),Y->dimension(),CHEV);

    // SCR1 -> V * y^-0.25
    for(auto j = 0ul; j < 2*NP; j++)
    for(auto i = 0ul; i < 2*NP; i++)
      CSCR1[i + 2*NP*j] = (*Y)(i,j) * std::pow(CHEV[j],-0.25);

    // Y = SCR1 * SCR1**H
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,
      CSCR1,2*NP,dcomplex(0.),Y->pointer(),Y->dimension());

    // Build the effective two component CH 
    cqmatrix::Matrix<dcomplex> FullCH2C(2*NP);
    FullCH2C.clear();

    // Copy P2P into spin diagonal blocks of 2C CH
    SetMat('N',2*NP,2*NP,dcomplex(1.),P2P,2*NP,FullCH2C.pointer(),2*NP);

    // Construct 2C CH in the uncontracted basis
    // 2C CH = Y * (V' + cp * X + X**H * cp + X**H * W' * X) * Y

    // SCR1 = cp * X
    for(auto j = 0; j < 2*NP; j++)
    for(auto i = 0; i < 2*NP; i++) {
      CSCR1[i + 2*NP*j] = SpeedOfLight * SS[i] * (*X)(i,j);
    }

    // 2C CH += SCR1 + SCR1**H
    MatAdd('N','N',2*NP,2*NP,dcomplex(1.),FullCH2C.pointer(),2*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C.pointer(),2*NP);
    MatAdd('N','C',2*NP,2*NP,dcomplex(1.),FullCH2C.pointer(),2*NP,dcomplex(1.),
      CSCR1,2*NP, FullCH2C.pointer(),2*NP);

    // SCR1 = X**H * W
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,dcomplex(1.),X->pointer(),X->dimension(),
               Wp,LDW,dcomplex(0.),CSCR1,2*NP);

    // 2C CH += SCR1 * X
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,dcomplex(1.),CSCR1,2*NP,
               X->pointer(),X->dimension(),dcomplex(1.),FullCH2C.pointer(),2*NP);

    // SCR1 = CH2C * Y
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,dcomplex(1.),FullCH2C.pointer(),2*NP,
               Y->pointer(),Y->dimension(),dcomplex(0.),CSCR1,2*NP);


    // 2C CH = Y * SCR1
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,dcomplex(1.),Y->pointer(),Y->dimension(),CSCR1,2*NP,
               dcomplex(0.),FullCH2C.pointer(),2*NP);


    // Transform the spin components of the 2C CH into R-space
    //
    // H(k) -> SUK * H(k) * (SUK)**H
    //
    // ** Using the fact that H(k) is hermetian
    // CSCR2 = SUK * H(k) -> CSCR2**H = H(k) * (SUK)**H
    // H(k) -> SUK * CSCR2**H
    //

    // first half of the scratch space is SUK(2*NP x 2*NP) 
    // matrix
    dcomplex   * SUK   = SCR1;
    dcomplex * CSCR2 = SUK+4*NP*NP;

    // Recover R-Space
    LUInv(2*NP,UK,2*NP);

    // SCR2 = SUK**H * CH2C 
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),UK,2*NP,
      FullCH2C.pointer(),2*NP,dcomplex(0.),CSCR2,2*NP);

    // 2C CH = SCR2 * SUK
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NP,2*NP,dcomplex(1.),CSCR2,2*NP,
      UK,2*NP,dcomplex(0.),FullCH2C.pointer(),2*NP);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"FullCH2C",FullCH2C.pointer(),2*NP,2*NP,2*NP);
#endif 

    // Allocate memory for the uncontracted spin components
    // of the 2C CH
    cqmatrix::PauliSpinorMatrices<dcomplex> HUn(
        FullCH2C.template spinScatter<dcomplex>(
            ssOptions_.hamiltonianOptions.OneESpinOrbit,ssOptions_.hamiltonianOptions.OneESpinOrbit));

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"CH0 in Primitive",HUn.S().pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"CHx in Primitive",HUn.X().pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"CHy in Primitive",HUn.Y().pointer(),NP,NP,NP);
    prettyPrintSmart(std::cout,"CHz in Primitive",HUn.Z().pointer(),NP,NP,NP);
#endif

    // Transform the spin components of the 2C CH into Contracted Basis
    *coreH = HUn.transform('C', mapPrim2Cont, NB, NB);

#ifdef DebugX2Cprint
    prettyPrintSmart(std::cout,"CH0 no Boettger",coreH->S().pointer(),NB,NB,NB);
    prettyPrintSmart(std::cout,"CHx no Boettger",coreH->X().pointer(),NB,NB,NB);
    prettyPrintSmart(std::cout,"CHy no Boettger",coreH->Y().pointer(),NB,NB,NB);
    prettyPrintSmart(std::cout,"CHz no Boettger",coreH->Z().pointer(),NB,NB,NB);
#endif

    CQMemManager::get().free(overlap, SCR1, SCPY, VCPY, CH4C, CHEV,
                     M, Ms, P2P, AVA, DOT, CROSSx, CROSSy, CROSSz);

#endif
  }


  template void X2C<double,double>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>);

  /**
   *  \brief Compute the picture change matrices UL, US
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C_Umatrix() {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NPU= nPrimUse_;
    size_t NB = basisSet_.nBasis;

    // UL = UK * Y * UK^-1 * UP2C
    // US = 2 * SpeedOfLight * UK * p^-1 * X * Y * UK^-1 * UP2C

    // 1.  UP2CSUK = UP2C * S * UK
    IntsT *UP2CS = CQMemManager::get().malloc<IntsT>(NB*NP);
    std::fill_n(UP2CS,NB*NP,IntsT(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,IntsT(1.),mapPrim2Cont,NB,
      uncontractedInts_.overlap->pointer(),NP,IntsT(0.),UP2CS,NB);
    IntsT *UP2CSUK = CQMemManager::get().malloc<IntsT>(4*NP*NPU);
    std::fill_n(UP2CSUK,4*NP*NPU,IntsT(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NPU,NP,IntsT(1.),UP2CS,NB,UK,NP,IntsT(0.),UP2CSUK,2*NB);
    SetMatDiag(NB,NPU,UP2CSUK,2*NB,UP2CSUK,2*NB);

    // 2. R^T = UP2C * S * UK * Y^T
    MatsT *RT = CQMemManager::get().malloc<MatsT>(4*NB*NPU);
    std::fill_n(RT,4*NB*NPU,MatsT(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NB,2*NPU,2*NPU,MatsT(1.),UP2CSUK,2*NB,
               Y->pointer(),Y->dimension(),MatsT(0.),RT,2*NB);

    // 3. Xp = 2 c p^-1 X
    double twoC = 2 * SpeedOfLight;
    double *twoCPinv = CQMemManager::get().malloc<double>(NPU);
    for(size_t i = 0; i < NPU; i++) twoCPinv[i] = twoC/p[i];
    MatsT *twoCPinvX = CQMemManager::get().malloc<MatsT>(4*NPU*NPU);
    for(size_t j = 0; j < 2*NPU; j++)
    for(size_t i = 0; i < NPU; i++) {
      twoCPinvX[i + 2*NPU*j] = twoCPinv[i] * (*X)(i,j);
      twoCPinvX[i + NPU + 2*NPU*j] = twoCPinv[i] * (*X)(i + NPU, j);
    }

    // 4. UK2c = [ UK  0  ]
    //           [ 0   UK ]
    IntsT *UK2c = UP2CSUK;
    SetMatDiag(NP,NPU,UK,NP,UK2c,2*NP);

    // 5. US = UK2c * Xp * RT^T
    UL = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(UL, 4*NP*NB, MatsT(0.));
    US = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(US, 4*NP*NB, MatsT(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NPU,2*NB,2*NPU,MatsT(1.),twoCPinvX,2*NPU,
      RT,2*NB,MatsT(0.),UL,2*NPU);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NPU,MatsT(1.),UK2c,2*NP,
      UL,2*NPU,MatsT(0.),US,2*NP);

    // 6. UL = UK2c * RT^T
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NB,2*NPU,MatsT(1.),UK2c,2*NP,
      RT,2*NB,MatsT(0.),UL,2*NP);

    CQMemManager::get().free(UP2CS, UP2CSUK, RT, twoCPinv, twoCPinvX);

  }

  template void X2C<dcomplex,double>::computeOneEX2C_Umatrix();

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_Umatrix() {
    
    //CErr("X2C + Complex Ints NYI",std::cout);

    // Mind that UK & p are 2-component for GIAO as T is no longer DiagMat
    // After X2C, stored UK is actually UK^-1

    // NPU is temperaily disabled in GIAO
    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // 1. UP2CSUK = UP2C * S * UK  (in 2 component)
    // Compute UP2CS
    dcomplex *UP2CS = CQMemManager::get().malloc<dcomplex>(4*NB*NP);
    std::fill_n(UP2CS,4*NB*NP,dcomplex(0.));
    dcomplex *UP2CSUK = CQMemManager::get().malloc<dcomplex>(4*NP*NP);
    std::fill_n(UP2CSUK,4*NP*NP,dcomplex(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,dcomplex(1.),mapPrim2Cont,NB,
      uncontractedInts_.overlap->pointer(),NP,dcomplex(0.),UP2CS,2*NB);
    SetMatDiag(NB,NP,UP2CS,2*NB,UP2CS,2*NB);
    // Recover UK
    LUInv(2*NP,UK,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NB,2*NP,2*NP,dcomplex(1.),UP2CS,2*NB,UK,2*NP,dcomplex(0.),UP2CSUK,2*NB);

    // 2. R^T = UP2CSUK * Y^T
    dcomplex *RT = CQMemManager::get().malloc<dcomplex>(4*NB*NP);
    std::fill_n(RT,4*NB*NP,dcomplex(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NB,2*NP,2*NP,dcomplex(1.),UP2CSUK,2*NB,
               Y->pointer(),Y->dimension(),dcomplex(0.),RT,2*NB);  

    // 3. Xp = 2 c p^-1 X
    // Mind that p is reverted by the end of computeOneEX2C 
    double twoC = 2 * SpeedOfLight;
    double *twoCPinv = CQMemManager::get().malloc<double>(2*NP);
    for(size_t i = 0; i < 2*NP; i++) twoCPinv[i] = twoC/p[i];
    dcomplex *twoCPinvX = CQMemManager::get().malloc<dcomplex>(4*NP*NP);
    for(size_t j = 0; j < 2*NP; j++)
    for(size_t i = 0; i < 2*NP; i++) {
      twoCPinvX[i + 2*NP*j] = twoCPinv[i] * (*X)(i,j);
    }

    // 4. US = UK2c * Xp * RT^T
    UL = CQMemManager::get().malloc<dcomplex>(4*NP*NB);
    std::fill_n(UL, 4*NP*NB, dcomplex(0.));
    US = CQMemManager::get().malloc<dcomplex>(4*NP*NB);
    std::fill_n(US, 4*NP*NB, dcomplex(0.));
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NB,2*NP,dcomplex(1.),twoCPinvX,2*NP,
      RT,2*NB,dcomplex(0.),UL,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,dcomplex(1.),UK,2*NP,
      UL,2*NP,dcomplex(0.),US,2*NP);

    // 5. UL = UK2c * RT^T
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NB,2*NP,dcomplex(1.),UK,2*NP,
      RT,2*NB,dcomplex(0.),UL,2*NP);

    CQMemManager::get().free(UP2CS, UP2CSUK, RT, twoCPinv, twoCPinvX);
    
  }

  template void X2C<double,double>::computeOneEX2C_Umatrix();

  /**
   *  \brief Compute the X2C Core Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C_UDU(EMPerturbation& emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // Allocate W separately  as it's needed later
    W = std::make_shared<cqmatrix::Matrix<MatsT>>(
        std::dynamic_pointer_cast<OnePRelInts<IntsT>>(
            uncontractedInts_.potential)->template formW<MatsT>());

    // T2c = [ T  0 ]
    //       [ 0  T ]
    const OnePInts<IntsT> &T2c = uncontractedInts_.kinetic->
        template spatialToSpinBlock<IntsT>();

    // V2c = [ V  0 ]
    //       [ 0  V ]
    const OnePInts<IntsT> &V2c = uncontractedInts_.potential->
        template spatialToSpinBlock<IntsT>();

    cqmatrix::Matrix<MatsT> Hx2c(2*NB);
    MatsT *SCR = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(SCR,4*NP*NB,MatsT(0.));

    // Hx2c = UL^H * T2c * US
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,MatsT(1.),T2c.pointer(),2*NP,
      US,2*NP,MatsT(0.),SCR,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NB,2*NB,2*NP,MatsT(1.),UL,2*NP,
      SCR,2*NP,MatsT(0.),Hx2c.pointer(),2*NB);
    // Hx2c += US^H * T2c * UL
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,MatsT(1.),T2c.pointer(),2*NP,
      UL,2*NP,MatsT(0.),SCR,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NB,2*NB,2*NP,MatsT(1.),US,2*NP,
      SCR,2*NP,MatsT(1.),Hx2c.pointer(),2*NB);
    // Hx2c -= US^H * T2c * US
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,MatsT(1.),T2c.pointer(),2*NP,
      US,2*NP,MatsT(0.),SCR,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NB,2*NB,2*NP,MatsT(-1.),US,2*NP,
      SCR,2*NP,MatsT(1.),Hx2c.pointer(),2*NB);
    // Hx2c += UL^H * V2c * UL
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,MatsT(1.),V2c.pointer(),2*NP,
      UL,2*NP,MatsT(0.),SCR,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NB,2*NB,2*NP,MatsT(1.),UL,2*NP,
      SCR,2*NP,MatsT(1.),Hx2c.pointer(),2*NB);
    // Hx2c += 1/(4*C**2) US^H * W * US
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NP,
      MatsT(0.25/SpeedOfLight/SpeedOfLight),W->pointer(),2*NP,
      US,2*NP,MatsT(0.),SCR,2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*NB,2*NB,2*NP,MatsT(1.),US,2*NP,
      SCR,2*NP,MatsT(1.),Hx2c.pointer(),2*NB);

    *coreH = Hx2c.template spinScatter<MatsT>(
        ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);

    CQMemManager::get().free(SCR);
  }

  template void X2C<dcomplex,double>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>);

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>);

  /**
   *  \brief Compute the X2C Core Hamiltonian correction to NR
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C_corr(EMPerturbation &emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

    computeOneEX2C(emPert, coreH);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> NRcoreH =
        std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NP);
    NRcoreH->clear();

    NRCoreH<MatsT, IntsT>(uncontractedInts_, ssOptions_.hamiltonianOptions)
        .computeNRCH(emPert, NRcoreH);

    *coreH -= NRcoreH->transform('C', mapPrim2Cont, NB, NB);

  }

  template void X2C<dcomplex,double>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>);

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>);

  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeFockX2C(EMPerturbation &emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> fockMatrix,
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> pchgDipole_,
      bool incore, double threshSchwarz) {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;
    nPrimUse_ = NP;

    // Prepare four-component (4C) single slater options
    SingleSlaterOptions fourCoptions(ssOptions_);
    fourCoptions.hamiltonianOptions.x2cType = X2C_TYPE::OFF;
    fourCoptions.refOptions.refType = isFourCRef;
    fourCoptions.refOptions.isKSRef = false;
    fourCoptions.refOptions.nC = 4;
    fourCoptions.refOptions.iCS = false;

    std::shared_ptr<Integrals<IntsT>> fourCInts = std::make_shared<Integrals<IntsT>>(uncontractedInts_);
    if (incore) {
      fourCInts->TPI =
          std::make_shared<InCore4indexTPI<IntsT>>(NP);
    } else {
      fourCInts->TPI =
          std::make_shared<DirectTPI<IntsT>>(
              uncontractedBasis_,uncontractedBasis_,
              molecule_,threshSchwarz);
    }

    // Construct 4C single slater object
    std::shared_ptr<SingleSlaterBase> ss =
        fourCoptions.buildSingleSlater(
            std::cout, molecule_, uncontractedBasis_, fourCInts);

    std::shared_ptr<SingleSlater<MatsT,IntsT>> ptr = std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ss);
    SingleSlater<MatsT,IntsT> &fourCompSS = *ptr;

    // Compute 4C core Hamiltonian
    fourCompSS.formCoreH(emPert, true);
    uncontractedInts_ = *(fourCompSS.aoints_);

    if (ssOptions_.hamiltonianOptions.x2cType == X2C_TYPE::ONEE) {

      // For One-electron X2C, solve core-Hamiltonian-only HC = eSC
      fourCompSS.fockBuilder = std::make_shared<MatrixFock<MatsT, IntsT>>(
          ssOptions_.hamiltonianOptions, *fourCompSS.coreH);

      fourCompSS.formGuess(fourCoptions);

    } else if (ssOptions_.hamiltonianOptions.x2cType == X2C_TYPE::FOCK) {

      fourCompSS.aoints_->computeAOTwoE(uncontractedBasis_, molecule_, emPert);

      // For Fock X2C, solve four-component SCF
      fourCompSS.formGuess(fourCoptions);
      fourCompSS.buildOrbitalModifierOptions();
      fourCompSS.runSCF(emPert);
    }

    ROOT_ONLY(ss->comm);

    computeFockX2C_Umatrix(fourCompSS.mo[0]);

    // Construct U matrix
    MatsT *U = CQMemManager::get().malloc<MatsT>(8*NP*NB);
    SetMat('N', NP, 2*NB, 1.0, UL, 2*NP, U, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, US, 2*NP, U + NP, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, UL + NP, 2*NP, U + 2*NP, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, US + NP, 2*NP, U + 3*NP, 4*NP);

    // Generate X2C core Hamiltonian
    cqmatrix::Matrix<MatsT> fourCompCoreH = fourCompSS.coreH->template spinGather<MatsT>();

    *coreH = fourCompCoreH.transform('N', U, 2*NB, 4*NP).template spinScatter<MatsT>(
        ssOptions_.hamiltonianOptions.OneESpinOrbit,ssOptions_.hamiltonianOptions.OneESpinOrbit);

    if (ssOptions_.hamiltonianOptions.x2cType == X2C_TYPE::FOCK
        and fockMatrix) {
      cqmatrix::Matrix<MatsT> fourCompFock = fourCompSS.fockMatrix->template spinGather<MatsT>();

      *fockMatrix = fourCompFock.transform('N', U, 2 * NB, 4 * NP).template spinScatter<MatsT>(
          ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);


      cqmatrix::Matrix<dcomplex> fourCompDipoleX = (*(fourCompSS.aoints_->lenElectric4C))[0].template spinGather<dcomplex>();
      cqmatrix::Matrix<dcomplex> fourCompDipoleY = (*(fourCompSS.aoints_->lenElectric4C))[1].template spinGather<dcomplex>();
      cqmatrix::Matrix<dcomplex> fourCompDipoleZ = (*(fourCompSS.aoints_->lenElectric4C))[2].template spinGather<dcomplex>();

      *pchgDipole_[0] = fourCompDipoleX.transform('N', U, 2 * NB, 4 * NP).template spinScatter<dcomplex>(
          ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);
      *pchgDipole_[1] = fourCompDipoleY.transform('N', U, 2 * NB, 4 * NP).template spinScatter<dcomplex>(
          ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);
      *pchgDipole_[2] = fourCompDipoleZ.transform('N', U, 2 * NB, 4 * NP).template spinScatter<dcomplex>(
          ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);

    }

    CQMemManager::get().free(U);

  } // X2C::computeFockX2C

  template void X2C<dcomplex,double>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>>, bool, double);

  template<> void X2C<dcomplex,dcomplex>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>>, bool, double) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>, 
      std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>>, bool, double);


  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeFockX2C_Umatrix(const cqmatrix::Matrix<MatsT> &fourCompMOSpin) {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // Compute transformation matrices in primitives
    cqmatrix::Matrix<IntsT> S2c(uncontractedInts_.overlap->matrix().template spatialToSpinBlock<IntsT>());
    cqmatrix::Matrix<IntsT> T2c(uncontractedInts_.kinetic->matrix().template spatialToSpinBlock<IntsT>());

    // Get and reorganize coefficients
    cqmatrix::Matrix<MatsT> fourCompMO(fourCompMOSpin.dimension());
    fourCompMO.clear();
    ReOrganizeMOSpin(fourCompMOSpin, fourCompMO);

    // Get pointers to "L" and "S" components of eigenvectors
    MatsT *coef = fourCompMO.pointer();
    size_t ldCoef = fourCompMO.dimension();
    MatsT *L = coef + 2*NP * ldCoef;
    MatsT *S = L + 2*NP;


    // Invert "L"; L -> L^-1
    LUInv(2*NP, L, ldCoef);

    // Compute X
    X = std::make_shared<cqmatrix::Matrix<MatsT>>(2*NP);

    // Form X = S * L^-1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,MatsT(1.),S,ldCoef,L,ldCoef,
               MatsT(0.),X->pointer(),X->dimension());

    // Compute UL and US

    // UL = S^-1/2 ( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2 S^1/2
    cqmatrix::Matrix<IntsT> Shalf(2*NP);
    cqmatrix::Matrix<IntsT> SinvHalf(2*NP);

    MatDiagFunc(std::function<double(double)>([](double x){ return std::sqrt(x); }),
                S2c.dimension(), S2c.pointer(), S2c.dimension(),
                Shalf.pointer(), Shalf.dimension());
    MatDiagFunc(std::function<double(double)>([](double x){ return 1.0/std::sqrt(x); }),
                S2c.dimension(), S2c.pointer(), S2c.dimension(),
                SinvHalf.pointer(), SinvHalf.dimension());

    // S + 1/2c^2 X^H T X
    const double TFact = 0.5 / (SpeedOfLight * SpeedOfLight);
    Y = std::make_shared<cqmatrix::Matrix<MatsT>>(S2c + TFact * T2c.transform('N', X->pointer(), X->dimension(), X->dimension()));

    // S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2
    *Y = Y->transform('N', SinvHalf.pointer(), SinvHalf.dimension(), SinvHalf.dimension());

    // ( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2
    MatDiagFunc(std::function<double(double)>([](double x){ return 1.0/std::sqrt(x); }),
                Y->dimension(), Y->pointer(), Y->dimension(),
                Y->pointer(), Y->dimension());

    cqmatrix::Matrix<MatsT> SCR(2*NP);

    // CSCR1 = (( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2 S^1/2)^T
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), Shalf.pointer(), Shalf.dimension(),
               Y->pointer(), Y->dimension(),
               MatsT(0.0), SCR.pointer(), SCR.dimension());

    // compute UL
    cqmatrix::Matrix<MatsT> ULsub(2*NP);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), SinvHalf.pointer(), SinvHalf.dimension(),
               SCR.pointer(), SCR.dimension(),
               MatsT(0.0), ULsub.pointer(), ULsub.dimension());

    // compute US = X UL
    cqmatrix::Matrix<MatsT> USsub(2*NP);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), X->pointer(), X->dimension(),
               ULsub.pointer(), ULsub.dimension(),
               MatsT(0.0), USsub.pointer(), USsub.dimension());

    // Compute the mappings from primitives to CGTOs
    mapPrim2Cont = CQMemManager::get().malloc<IntsT>(NB*NP);
    basisSet_.makeMapPrim2Cont(uncontractedInts_.overlap->pointer(), mapPrim2Cont);
    MatsT *P2C2c = CQMemManager::get().malloc<MatsT>(4*NB*NP);
    SetMatDiag(NB, NP, mapPrim2Cont, NB, P2C2c, 2*NB);

    // Contract transformation matrices with P2C mapping
    UL = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(UL, 4*NP*NB, MatsT(0.));
    US = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(US, 4*NP*NB, MatsT(0.));

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NP, 2*NB, 2*NP,
               MatsT(1.), ULsub.pointer(), 2*NP, P2C2c, 2*NB,
               MatsT(0.), UL, 2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NP, 2*NB, 2*NP,
               MatsT(1.), USsub.pointer(), 2*NP, P2C2c, 2*NB,
               MatsT(0.), US, 2*NP);

  }

  template<> void X2C<dcomplex, dcomplex>::computeFockX2C_Umatrix(const cqmatrix::Matrix<dcomplex>&) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<dcomplex, double>::computeFockX2C_Umatrix(const cqmatrix::Matrix<dcomplex>&);

  template void X2C<double, double>::computeFockX2C_Umatrix(const cqmatrix::Matrix<double>&);




  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::compute_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions) {

    if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::ONEE)
      ROOT_ONLY(ss->comm);

    std::shared_ptr<X2C<MatsT, IntsT>> x2c;

    if (ssOptions.hamiltonianOptions.AtomicX2C) {

      if (ssOptions.hamiltonianOptions.x2cType != X2C_TYPE::ONEE)
        CErr("Atomic X2C is only implemented for One-electron X2C.");

      // Compute non-relativistic integrals for ALH and DLH algorithms
      if (ssOptions.hamiltonianOptions.AtomicX2CType.diagonalOnly) {

        // Build a non-relativistic HamiltonianOptions to avoid pVp and pxVp integral evaluation
        HamiltonianOptions nonRelHoption(ssOptions.hamiltonianOptions);
        nonRelHoption.OneEScalarRelativity = false;
        nonRelHoption.OneESpinOrbit = false;

        std::vector<std::pair<OPERATOR,size_t>> ops{{KINETIC,0}, {NUCLEAR_POTENTIAL,0}};
        aoints->computeAOOneP(mol, basis, emPert, ops, nonRelHoption);
      }

      x2c = std::make_shared<AtomicX2C<MatsT,IntsT>>(
          *std::dynamic_pointer_cast<Integrals<IntsT>>(aoints),
          mol, basis, ssOptions);
    } else {
      x2c = std::make_shared<X2C<MatsT,IntsT>>(
          *std::dynamic_pointer_cast<Integrals<IntsT>>(aoints),
          mol, basis, ssOptions);
    }

    SingleSlater<MatsT, IntsT> &ref = *std::dynamic_pointer_cast<SingleSlater<MatsT, IntsT>>(ss);
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH =
        std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(
            basis.nBasis,
            ssOptions.hamiltonianOptions.OneESpinOrbit,
            ssOptions.hamiltonianOptions.OneESpinOrbit);

    if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::ONEE) {
      x2c->computeOneEX2C(emPert, coreH);

      // Added SNSOScale two-electron relativistic effect
      if (ssOptions.hamiltonianOptions.SNSO)
        switch (ssOptions.hamiltonianOptions.snsoType) {
          case SNSO_TYPE::BOETTGER:
          case SNSO_TYPE::DC:
          case SNSO_TYPE::DCB:
            x2c->SNSOScale(coreH, ssOptions.hamiltonianOptions.snsoType);
            break;
          case SNSO_TYPE::ROW_DEP_DCB:
            x2c->RowDepDCB_SNSO(coreH);
            break;
          default:
            CErr("Unknown SNSO type",std::cout);
        }

#ifdef DebugX2Cprint2
    prettyPrintSmart(std::cout,"CH0",coreH->S().pointer(),basis.nBasis,basis.nBasis,basis.nBasis);
    prettyPrintSmart(std::cout,"CH1",coreH->Z().pointer(),basis.nBasis,basis.nBasis,basis.nBasis);
    prettyPrintSmart(std::cout,"CH2",coreH->Y().pointer(),basis.nBasis,basis.nBasis,basis.nBasis);
    prettyPrintSmart(std::cout,"CH3",coreH->X().pointer(),basis.nBasis,basis.nBasis,basis.nBasis);
#endif
        
    }

    if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::FOCK) {

      std::shared_ptr<DirectTPI<IntsT>> tpi =
          std::dynamic_pointer_cast<DirectTPI<IntsT>>(ref.aoints_->TPI);
      double threshSchwarz = 0.0;
      bool incore = tpi == nullptr;

      if (not incore)
        threshSchwarz = tpi->threshSchwarz();

      // Initialize the returned dipole matrix
      ref.pchgDipole_[0] = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(basis.nBasis, true, true);
      ref.pchgDipole_[1] = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(basis.nBasis, true, true);
      ref.pchgDipole_[2] = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(basis.nBasis, true, true);

      x2c->computeFockX2C(emPert, coreH, ref.fockMatrix, ref.pchgDipole_, incore, threshSchwarz);

#ifdef CQ_ENABLE_MPI
      // BCast fockMatrix to all MPI processes
      if( MPISize(ss->comm) > 1 ) {
        std::cerr  << "  *** Scattering the X2C Fock ***\n";
        size_t NB = ref.fockMatrix->dimension();
        for(auto mat : ref.fockMatrix->SZYXPointers())
          MPIBCast(mat,NB*NB,0,ss->comm);
      }
#endif

      ref.fockBuilder = std::make_shared<MatrixFock<MatsT, IntsT>>(
          ssOptions.hamiltonianOptions, *ref.fockMatrix);
    }

    ref.coreHBuilder = std::make_shared<MatrixCoreH<MatsT, IntsT>>(
        *std::dynamic_pointer_cast<Integrals<IntsT>>(aoints),
        ssOptions.hamiltonianOptions, std::move(*coreH));


    // U transformation to bin file
    if(not (ssOptions.hamiltonianOptions.AtomicX2C
            and ssOptions.hamiltonianOptions.AtomicX2CType.diagonalOnly == true))
      x2c->saveX2C(ss);

//    CErr("Requested X2C type NYI.");

  }

// TangDD Add X2C + GIAO
  template void X2C<dcomplex, dcomplex>::compute_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);
    //CErr("X2C + Complex Ints NYI",std::cout);

  template void X2C<dcomplex, double>::compute_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);

  template void X2C<double, double>::compute_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);


  void compute_X2C_CoreH_Fock(Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions) {

    if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

      X2C<double, double>::compute_CoreH_Fock(
          mol, basis, aoints, emPert, ss, ssOptions);

    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

      X2C<dcomplex, double>::compute_CoreH_Fock(
          mol, basis, aoints, emPert, ss, ssOptions);

    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

      X2C<dcomplex, dcomplex>::compute_CoreH_Fock(
          mol, basis, aoints, emPert, ss, ssOptions);

    } else {

      CErr("Real X2C + Complex Ints invalid",std::cout);
    }

  }

  /**
   *  \brief Save the X2C transformation
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::saveX2C(std::shared_ptr<SingleSlaterBase> ss) {
    ROOT_ONLY(ss->comm);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    if( ss->savFile.exists() ){

      if( UL == nullptr or US == nullptr ){

        std::cout << "    * Saving X2C. U matrices not found. Recomputing." << std::endl << std::endl;
        computeOneEX2C_Umatrix();

      }

      // dimensions: row 2*NP, column 2*NB
      size_t Urow = 2*NP;
      size_t Ucol = 2*NB;

      std::string prefix = "X2C/";
      ss->savFile.safeWriteData(prefix + "UL", UL, {Urow, Ucol});
      ss->savFile.safeWriteData(prefix + "US", US, {Urow, Ucol});

    } else CErr("Could not find savFile in saveX2C");

  }

  template void X2C<dcomplex,double>::saveX2C(std::shared_ptr<SingleSlaterBase>);

  template void X2C<dcomplex,dcomplex>::saveX2C(std::shared_ptr<SingleSlaterBase>);

  template void X2C<double,double>::saveX2C(std::shared_ptr<SingleSlaterBase>);

}; // namespace ChronusQ

