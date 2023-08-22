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

namespace ChronusQ {

  template <typename MatsT>
  void GatherUSpin(const SquareMatrix<MatsT> &UL, const SquareMatrix<MatsT> &US, MatsT *U) {
    size_t NP = UL.dimension() / 2;
    SetMat('N', NP, 2*NP, 1.0, UL.pointer(), 2*NP, U, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, US.pointer(), 2*NP, U + NP, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, UL.pointer() + NP, 2*NP, U + 2*NP, 4*NP);
    SetMat('N', NP, 2*NP, 1.0, US.pointer() + NP, 2*NP, U + 3*NP, 4*NP);
  }

  template void GatherUSpin(const SquareMatrix<double> &UL, const SquareMatrix<double> &US, double *U);
  template void GatherUSpin(const SquareMatrix<dcomplex> &UL, const SquareMatrix<dcomplex> &US, dcomplex *U);

  template <typename MatsT>
  void ReOrganizeMOSpin(const SquareMatrix<MatsT> &moSpin, SquareMatrix<MatsT> &mo) {

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

  template void ReOrganizeMOSpin(const SquareMatrix<double> &moSpin, SquareMatrix<double> &mo);
  template void ReOrganizeMOSpin(const SquareMatrix<dcomplex> &moSpin, SquareMatrix<dcomplex> &mo);



  /**
   *  \brief Boettger scaling for spin-orbit operator
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::BoettgerScale(
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH) {

    size_t NB = basisSet_.nBasis;

    size_t n1, n2;
    std::array<double,6> Ql={0.,2.,10.,28.,60.,110.};

    if( this->basisSet_.maxL > 5 ) CErr("Boettger scaling for L > 5 NYI");

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
   *  \brief Compute the X2C Core Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C(EMPerturbation &emPert,
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH) {

#ifdef REAL_SPACE_X2C_ALGORITHM
    computeX2C_realSpace(emPert, coreH);
    return;
#endif

    IntsT* XXX = reinterpret_cast<IntsT*>(NULL);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    uncontractedInts_.computeAOOneP(memManager_,
        molecule_, uncontractedBasis_, emPert,
        {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0}},
        ssOptions_.hamiltonianOptions);

    // Make copy of integrals
    IntsT *overlap   = memManager_.malloc<IntsT>(NP*NP);
    std::copy_n(uncontractedInts_.overlap->pointer(), NP*NP, overlap);

    // Compute the mappings from primitives to CGTOs
    mapPrim2Cont = memManager_.malloc<IntsT>(NP*NB);
    basisSet_.makeMapPrim2Cont(overlap,mapPrim2Cont,memManager_);

    // Allocate Scratch Space (enough for 2*NP x 2*NP complex matricies)
    IntsT *SCR1  = memManager_.malloc<IntsT>(8*NP*NP);
    MatsT *CSCR1 = reinterpret_cast<MatsT*>(SCR1);

    // Singular value storage (initially S then T)
    p = memManager_.malloc<double>(NP);
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
    UK = memManager_.malloc<IntsT>(NP*NPU);

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
    MatsT *CH4C = memManager_.malloc<MatsT>(16*NPU*NPU);
    std::fill_n(CH4C,16*NPU*NPU,MatsT(0.));

    // Allocate W separately  as it's needed later
    size_t LDW = 2*NPU;
    SquareMatrix<MatsT> Wp(potential->template formW<MatsT>());

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
    double *CHEV = memManager_.malloc<double>(4*NPU);

    HermetianEigen('V','U',4*NPU,CH4C,4*NPU,CHEV,memManager_);


    // Get pointers to "L" and "S" components of eigenvectors
    MatsT *L = CH4C + 8*NPU*NPU;
    MatsT *S = L + 2*NPU;


    // Invert "L"; L -> L^-1
    LUInv(2*NPU,L,4*NPU,memManager_);


    // Reuse the charge conjugated space for X and Y
    X = std::make_shared<SquareMatrix<MatsT>>(memManager_, 2*NPU);
    X->clear();
    Y = std::make_shared<SquareMatrix<MatsT>>(memManager_, 2*NPU);
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
                2*NPU, Y->pointer(), Y->dimension(), Y->pointer(), Y->dimension(), memManager_);

    // Build the effective two component CH in "L"
    SquareMatrix<MatsT> FullCH2C(memManager_, 2*NPU);

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
    PauliSpinorSquareMatrices<MatsT> HUn(
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

    memManager_.free(overlap, SCR1, CH4C, CHEV);


  }

  template void X2C<dcomplex,double>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>);

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<double>>);

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
    IntsT *UP2CS = memManager_.malloc<IntsT>(NB*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NP,NP,IntsT(1.),mapPrim2Cont,NB,
      uncontractedInts_.overlap->pointer(),NP,IntsT(0.),UP2CS,NB);
    IntsT *UP2CSUK = memManager_.malloc<IntsT>(4*NP*NPU);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NPU,NP,IntsT(1.),UP2CS,NB,UK,NP,IntsT(0.),UP2CSUK,2*NB);
    SetMatDiag(NB,NPU,UP2CSUK,2*NB,UP2CSUK,2*NB);

    // 2. R^T = UP2C * S * UK * Y^T
    MatsT *RT = memManager_.malloc<MatsT>(4*NB*NPU);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NB,2*NPU,2*NPU,MatsT(1.),UP2CSUK,2*NB,
               Y->pointer(),Y->dimension(),MatsT(0.),RT,2*NB);

    // 3. Xp = 2 c p^-1 X
    double twoC = 2 * SpeedOfLight;
    double *twoCPinv = memManager_.malloc<double>(NPU);
    for(size_t i = 0; i < NPU; i++) twoCPinv[i] = twoC/p[i];
    MatsT *twoCPinvX = memManager_.malloc<MatsT>(4*NPU*NPU);
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
    UL = memManager_.malloc<MatsT>(4*NP*NB);
    US = memManager_.malloc<MatsT>(4*NP*NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NPU,2*NB,2*NPU,MatsT(1.),twoCPinvX,2*NPU,
      RT,2*NB,MatsT(0.),UL,2*NPU);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*NP,2*NB,2*NPU,MatsT(1.),UK2c,2*NP,
      UL,2*NPU,MatsT(0.),US,2*NP);

    // 6. UL = UK2c * RT^T
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,2*NP,2*NB,2*NPU,MatsT(1.),UK2c,2*NP,
      RT,2*NB,MatsT(0.),UL,2*NP);

    memManager_.free(UP2CS, UP2CSUK, RT, twoCPinv, twoCPinvX);

  }

  template void X2C<dcomplex,double>::computeOneEX2C_Umatrix();

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_Umatrix() {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C_Umatrix();

  /**
   *  \brief Compute the X2C Core Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C_UDU(EMPerturbation& emPert,
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH) {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // Allocate W separately  as it's needed later
    W = std::make_shared<SquareMatrix<MatsT>>(
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

    SquareMatrix<MatsT> Hx2c(memManager_, 2*NB);
    MatsT *SCR = memManager_.malloc<MatsT>(4*NP*NB);

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

    memManager_.free(SCR);
  }

  template void X2C<dcomplex,double>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>);

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C_UDU(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<double>>);

  /**
   *  \brief Compute the X2C Core Hamiltonian correction to NR
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeOneEX2C_corr(EMPerturbation &emPert,
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH) {

    computeOneEX2C(emPert, coreH);

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> NRcoreH =
        std::make_shared<PauliSpinorSquareMatrices<MatsT>>(memManager_, NP);
    NRcoreH->clear();

    NRCoreH<MatsT, IntsT>(uncontractedInts_, ssOptions_.hamiltonianOptions)
        .computeNRCH(emPert, NRcoreH);

    *coreH -= NRcoreH->transform('C', mapPrim2Cont, NB, NB);

  }

  template void X2C<dcomplex,double>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>);

  template<> void X2C<dcomplex,dcomplex>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeOneEX2C_corr(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<double>>);

  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeFockX2C(EMPerturbation &emPert,
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH,
      std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> fockMatrix,
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
          std::make_shared<InCore4indexTPI<IntsT>>(memManager_, NP);
    } else {
      fourCInts->TPI =
          std::make_shared<DirectTPI<IntsT>>(
              memManager_,uncontractedBasis_,uncontractedBasis_,
              molecule_,threshSchwarz);
    }

    // Construct 4C single slater object
    std::shared_ptr<SingleSlaterBase> ss =
        fourCoptions.buildSingleSlater(
            std::cout, memManager_, molecule_, uncontractedBasis_, fourCInts);

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
    MatsT *U = memManager_.malloc<MatsT>(8*NP*NB);
    SetMat('N', NP, 2*NB, 1.0, UL, 2*NP, U, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, US, 2*NP, U + NP, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, UL + NP, 2*NP, U + 2*NP, 4*NP);
    SetMat('N', NP, 2*NB, 1.0, US + NP, 2*NP, U + 3*NP, 4*NP);

    // Generate X2C core Hamiltonian
    SquareMatrix<MatsT> fourCompCoreH = fourCompSS.coreH->template spinGather<MatsT>();

    *coreH = fourCompCoreH.transform('N', U, 2*NB, 4*NP).template spinScatter<MatsT>(
        ssOptions_.hamiltonianOptions.OneESpinOrbit,ssOptions_.hamiltonianOptions.OneESpinOrbit);

    if (ssOptions_.hamiltonianOptions.x2cType == X2C_TYPE::FOCK
        and fockMatrix) {
      SquareMatrix<MatsT> fourCompFock = fourCompSS.fockMatrix->template spinGather<MatsT>();

      *fockMatrix = fourCompFock.transform('N', U, 2 * NB, 4 * NP).template spinScatter<MatsT>(
          ssOptions_.hamiltonianOptions.OneESpinOrbit, ssOptions_.hamiltonianOptions.OneESpinOrbit);
    }

    memManager_.free(U);

  } // X2C::computeFockX2C

  template void X2C<dcomplex,double>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>, bool, double);

  template<> void X2C<dcomplex,dcomplex>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>,
      std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>, bool, double) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<double,double>::computeFockX2C(EMPerturbation&,
      std::shared_ptr<PauliSpinorSquareMatrices<double>>,
      std::shared_ptr<PauliSpinorSquareMatrices<double>>, bool, double);


  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::computeFockX2C_Umatrix(const SquareMatrix<MatsT> &fourCompMOSpin) {

    size_t NP = uncontractedBasis_.nPrimitive;
    size_t NB = basisSet_.nBasis;

    // Compute transformation matrices in primitives
    SquareMatrix<IntsT> S2c(uncontractedInts_.overlap->matrix().template spatialToSpinBlock<IntsT>());
    SquareMatrix<IntsT> T2c(uncontractedInts_.kinetic->matrix().template spatialToSpinBlock<IntsT>());

    // Get and reorganize coefficients
    SquareMatrix<MatsT> fourCompMO(memManager_, fourCompMOSpin.dimension());
    fourCompMO.clear();
    ReOrganizeMOSpin(fourCompMOSpin, fourCompMO);

    // Get pointers to "L" and "S" components of eigenvectors
    MatsT *coef = fourCompMO.pointer();
    size_t ldCoef = fourCompMO.dimension();
    MatsT *L = coef + 2*NP * ldCoef;
    MatsT *S = L + 2*NP;


    // Invert "L"; L -> L^-1
    LUInv(2*NP, L, ldCoef, memManager_);

    // Compute X
    X = std::make_shared<SquareMatrix<MatsT>>(memManager_, 2*NP);

    // Form X = S * L^-1
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               2*NP,2*NP,2*NP,MatsT(1.),S,ldCoef,L,ldCoef,
               MatsT(0.),X->pointer(),X->dimension());

    // Compute UL and US

    // UL = S^-1/2 ( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2 S^1/2
    SquareMatrix<IntsT> Shalf(memManager_, 2*NP);
    SquareMatrix<IntsT> SinvHalf(memManager_, 2*NP);

    MatDiagFunc(std::function<double(double)>([](double x){ return std::sqrt(x); }),
                S2c.dimension(), S2c.pointer(), S2c.dimension(),
                Shalf.pointer(), Shalf.dimension(), memManager_);
    MatDiagFunc(std::function<double(double)>([](double x){ return 1.0/std::sqrt(x); }),
                S2c.dimension(), S2c.pointer(), S2c.dimension(),
                SinvHalf.pointer(), SinvHalf.dimension(), memManager_);

    // S + 1/2c^2 X^H T X
    const double TFact = 0.5 / (SpeedOfLight * SpeedOfLight);
    Y = std::make_shared<SquareMatrix<MatsT>>(S2c + TFact * T2c.transform('N', X->pointer(), X->dimension(), X->dimension()));

    // S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2
    *Y = Y->transform('N', SinvHalf.pointer(), SinvHalf.dimension(), SinvHalf.dimension());

    // ( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2
    MatDiagFunc(std::function<double(double)>([](double x){ return 1.0/std::sqrt(x); }),
                Y->dimension(), Y->pointer(), Y->dimension(),
                Y->pointer(), Y->dimension(), memManager_);

    SquareMatrix<MatsT> SCR(memManager_, 2*NP);

    // CSCR1 = (( S^-1/2 (S + 1/2c^2 X^H T X) S^-1/2 )^-1/2 S^1/2)^T
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::Trans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), Shalf.pointer(), Shalf.dimension(),
               Y->pointer(), Y->dimension(),
               MatsT(0.0), SCR.pointer(), SCR.dimension());

    // compute UL
    SquareMatrix<MatsT> ULsub(memManager_, 2*NP);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), SinvHalf.pointer(), SinvHalf.dimension(),
               SCR.pointer(), SCR.dimension(),
               MatsT(0.0), ULsub.pointer(), ULsub.dimension());

    // compute US = X UL
    SquareMatrix<MatsT> USsub(memManager_, 2*NP);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
               2*NP, 2*NP, 2*NP,
               MatsT(1.0), X->pointer(), X->dimension(),
               ULsub.pointer(), ULsub.dimension(),
               MatsT(0.0), USsub.pointer(), USsub.dimension());

    // Compute the mappings from primitives to CGTOs
    mapPrim2Cont = memManager_.malloc<IntsT>(NB*NP);
    basisSet_.makeMapPrim2Cont(uncontractedInts_.overlap->pointer(), mapPrim2Cont,memManager_);
    MatsT *P2C2c = memManager_.malloc<MatsT>(4*NB*NP);
    SetMatDiag(NB, NP, mapPrim2Cont, NB, P2C2c, 2*NB);

    // Contract transformation matrices with P2C mapping
    UL = memManager_.malloc<MatsT>(4*NP*NB);
    US = memManager_.malloc<MatsT>(4*NP*NB);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NP, 2*NB, 2*NP,
               MatsT(1.), ULsub.pointer(), 2*NP, P2C2c, 2*NB,
               MatsT(0.), UL, 2*NP);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,
               2*NP, 2*NB, 2*NP,
               MatsT(1.), USsub.pointer(), 2*NP, P2C2c, 2*NB,
               MatsT(0.), US, 2*NP);

  }

  template<> void X2C<dcomplex, dcomplex>::computeFockX2C_Umatrix(const SquareMatrix<dcomplex>&) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<dcomplex, double>::computeFockX2C_Umatrix(const SquareMatrix<dcomplex>&);

  template void X2C<double, double>::computeFockX2C_Umatrix(const SquareMatrix<double>&);




  /**
   *  \brief Compute the X2C Core Hamiltonian in real space
   */
  template <typename MatsT, typename IntsT>
  void X2C<MatsT, IntsT>::compute_CoreH_Fock(CQMemManager &mem, Molecule &mol,
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
        aoints->computeAOOneP(mem, mol, basis, emPert, ops, nonRelHoption);
      }

      x2c = std::make_shared<AtomicX2C<MatsT,IntsT>>(
          *std::dynamic_pointer_cast<Integrals<IntsT>>(aoints),
          mem, mol, basis, ssOptions);
    } else {
      x2c = std::make_shared<X2C<MatsT,IntsT>>(
          *std::dynamic_pointer_cast<Integrals<IntsT>>(aoints),
          mem, mol, basis, ssOptions);
    }

    SingleSlater<MatsT, IntsT> &ref = *std::dynamic_pointer_cast<SingleSlater<MatsT, IntsT>>(ss);
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH =
        std::make_shared<PauliSpinorSquareMatrices<MatsT>>(
            mem, basis.nBasis,
            ssOptions.hamiltonianOptions.OneESpinOrbit,
            ssOptions.hamiltonianOptions.OneESpinOrbit);

    if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::ONEE) {
      x2c->computeOneEX2C(emPert, coreH);

      // Added BoettgerScale two-electron relativistic effect
      if (ssOptions.hamiltonianOptions.Boettger)
        x2c->BoettgerScale(coreH);
    }

    if (ssOptions.hamiltonianOptions.x2cType == X2C_TYPE::FOCK) {

      std::shared_ptr<DirectTPI<IntsT>> tpi =
          std::dynamic_pointer_cast<DirectTPI<IntsT>>(ref.aoints_->TPI);
      double threshSchwarz = 0.0;
      bool incore = tpi == nullptr;

      if (not incore)
        threshSchwarz = tpi->threshSchwarz();

      x2c->computeFockX2C(emPert, coreH, ref.fockMatrix, incore, threshSchwarz);

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

//    CErr("Requested X2C type NYI.");

  }

  template <> void X2C<dcomplex, dcomplex>::compute_CoreH_Fock(CQMemManager &mem, Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template void X2C<dcomplex, double>::compute_CoreH_Fock(CQMemManager &mem, Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);

  template void X2C<double, double>::compute_CoreH_Fock(CQMemManager &mem, Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions);


  void compute_X2C_CoreH_Fock(CQMemManager &mem, Molecule &mol,
      BasisSet &basis, std::shared_ptr<IntegralsBase> aoints,
      EMPerturbation &emPert,
      std::shared_ptr<SingleSlaterBase> ss, SingleSlaterOptions ssOptions) {

    if(auto p = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss)) {

      X2C<double, double>::compute_CoreH_Fock(
          mem, mol, basis, aoints, emPert, ss, ssOptions);

    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss)) {

      X2C<dcomplex, double>::compute_CoreH_Fock(
          mem, mol, basis, aoints, emPert, ss, ssOptions);

    } else if(auto p = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss)) {

      X2C<dcomplex, dcomplex>::compute_CoreH_Fock(
          mem, mol, basis, aoints, emPert, ss, ssOptions);

    } else {

      CErr("Real X2C + Complex Ints invalid",std::cout);
    }
  }

}; // namespace ChronusQ

