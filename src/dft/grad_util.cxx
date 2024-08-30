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

#include <matrix.hpp>
#include <dft/util.hpp>

namespace ChronusQ {

    /**
   *  \brief evaluates the nuclear gradient of V (Grad_Den, Grad_GDENX, Grad_GDENY and Grad_GDENZ)
   *  varibles for a given density in input in DENMAT 
   *  (scalar or magnetization).
   *
   *
   *  \param [in]  typ        Type of evaluation to perform (gradient, etc)
   *  \param [in]  NPts       Number of points in the batch
   *  \param [in]  NBE        Effective number of basis to be evalauted (only shell actives)
   *  \param [in]  NB         Total Number of basis.
   *  \param [in]  subMatCut  Pair to handle the cut of the shell submatrix to be evaluated.
   *  \param [in]  SCR1       Pointer to an NB*NB scratch.
   *  \param [in]  SCR2       Pointer to an NB*NPts scratch.
   *  \param [in]  DENMAT     Pointer to 1PDM (scalar, Mk).
   *  \param [out] Den        Pointer to the V variable - SCALAR/Mk
   *  \param [out] GDenX      Pointer to the V variable - Gradient X comp of SCALAR/Mk
   *  \param [out] GDenY      Pointer to the V variable - Gradient Y comp of SCALAR/Mk
   *  \param [out] GDenZ      Pointer to the V variable - Gradient Z comp of SCALAR/Mk
   *  \param [in]  BasisScr   Pointer to Basis set evaluated over batch of points.
   */
  void evalDenGrad(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB,
    std::vector<std::pair<size_t,size_t>> &subMatCut, double* SCR1,
    double *SCR2, std::vector<std::vector<double*>>SCR3, std::vector<std::vector<double*>>SCR4, 
    std::vector<double*> SCR5,
    std::vector<std::vector<double*>>GDENMAT, double *DENMAT, 
    std::vector<double*>GDenX, std::vector<double*>GDenY, std::vector<double*>GDenZ,
    std::vector<double*>GGDenxX, std::vector<double*>GGDenxY, std::vector<double*>GGDenxZ,
    std::vector<double*>GGDenyX, std::vector<double*>GGDenyY, std::vector<double*>GGDenyZ, 
    std::vector<double*>GGDenzX, std::vector<double*>GGDenzY, std::vector<double*>GGDenzZ, 
    double *BasisScr, double *BasisGradScr, size_t nAtoms, BasisSet &basisSet){

    size_t IOff = NPts*NBE; 

    // effective orbitals for gradient
    std::vector<std::vector<size_t>> effOrbsForAtom(nAtoms);
    std::vector<size_t> totOrbsForAtom;

    // loop over pairs in subMatCut
    for (size_t pi = 0; pi < subMatCut.size(); pi++) {
      for (size_t oi = subMatCut[pi].first; oi < subMatCut[pi].second; oi++) {
        for (int i = nAtoms-1; i >=0; i--) {
          if (oi >= basisSet.mapCen2BfSt[i]) {
            effOrbsForAtom[i].emplace_back(oi);
            break;
          }
        }
        totOrbsForAtom.emplace_back(oi);
      }
    }

    for (size_t ic = 0; ic < nAtoms; ic++) {
      for (size_t i = 0; i < effOrbsForAtom[ic].size(); i++) {
        auto it = std::find(totOrbsForAtom.begin(), totOrbsForAtom.end(), effOrbsForAtom[ic][i]);
        size_t oi = std::distance(totOrbsForAtom.begin(), it);
        effOrbsForAtom[ic][i] = oi;
      }
    }

    // SCR1: effective part of the total density matrix (NBE * NBE)
    SubMatSet(NB,NB,NBE,NBE,DENMAT,NB,SCR1,NBE,subMatCut);           

    //// SCR3: effective part of the dP / dR (NAtom * 3 * NBE * NBE)
    //for (size_t ic = 0; ic < nAtoms; ic++)
    //  for(int xyz = 0; xyz < 3; xyz++)
    //    SubMatSet(NB,NB,NBE,NBE,GDENMAT[ic][xyz],NB,SCR3[ic][xyz],NBE,subMatCut);

    // Obtain Sum_nu P_mu_nu Phi_nu ( NBE * NPts ) -> SCR2
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,1.,SCR1,NBE,BasisScr,NBE,0.,SCR2,NBE);

    //// Obtain Sum_nu dP_mu_nu/dR Phi_nu ( NAtom * 3 * NBE * NPts ) -> SCR4
    //for (size_t ic = 0; ic < nAtoms; ic++) 
    //  for(int xyz = 0; xyz < 3; xyz++)
    //    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,1.,SCR3[ic][xyz],NBE,BasisScr,NBE,0.,SCR4[ic][xyz],NBE);

    // If only the density gradient is needed
    if( typ != GRADIENT ) {
      for(auto iPt = 0; iPt < NPts; iPt++) {
        // loop over atom center
        for(auto ic = 0; ic < nAtoms; ic++) {
          size_t NBEAtom = effOrbsForAtom[ic].size();
          GDenX[ic][iPt] = 0.;
          GDenY[ic][iPt] = 0.;
          GDenZ[ic][iPt] = 0.;
          const size_t NBEiPt = iPt*NBE;
          const double *SCR_cur   = SCR2 + NBEiPt;
          //const double *SCR2_curX  = SCR4[ic][0] + NBEiPt;
          //const double *SCR2_curY  = SCR4[ic][1] + NBEiPt;
          //const double *SCR2_curZ  = SCR4[ic][2] + NBEiPt;
          const double *B_cur    = BasisScr + NBEiPt;
          const double *B_curX   = BasisGradScr + NBEiPt;
          const double *B_curY   = B_curX + IOff;
          const double *B_curZ   = B_curY + IOff;

          //// If this is the atomic center, subtract ALL orbitals
          //if ( ic == iAtm and this->scfControls.grid_resp ) {
          ////if ( false ) {
          //  for ( auto& atOrbs: effOrbsForAtom ) {
          //    for ( auto& oi: atOrbs ) {
          //      GDenX[ic][iPt] -= SCR_cur[oi] * B_curX[oi];
          //      GDenY[ic][iPt] -= SCR_cur[oi] * B_curY[oi];
          //      GDenZ[ic][iPt] -= SCR_cur[oi] * B_curZ[oi];
          //    }
          //  }
          //}
      
          for(size_t j = 0; j < NBEAtom; j++) { 
            size_t oi = effOrbsForAtom[ic][j];
            GDenX[ic][iPt] += SCR_cur[oi] * B_curX[oi];
            GDenY[ic][iPt] += SCR_cur[oi] * B_curY[oi];
            GDenZ[ic][iPt] += SCR_cur[oi] * B_curZ[oi];
          }

          GDenX[ic][iPt] *= 2.0;
          GDenY[ic][iPt] *= 2.0;
          GDenZ[ic][iPt] *= 2.0;

          //// add contribution from density matrix gradient
          //for(size_t j = 0; j < NBE; j++) {
          //  GDenX[ic][iPt] += SCR2_curX[j] * B_cur[j]; 
          //  GDenY[ic][iPt] += SCR2_curY[j] * B_cur[j]; 
          //  GDenZ[ic][iPt] += SCR2_curZ[j] * B_cur[j]; 
          //}
        }
      }
    } else { // gradient of density gradient is also needed

      // Obtain Sum_nu P_mu_nu dPhi_nu/dr (3 * NBE * NPts)
      for(int xyz = 0; xyz < 3; xyz++)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,1.,SCR1,NBE,BasisScr+IOff+xyz*IOff,NBE,0.,SCR5[xyz],NBE);

      for(auto iPt = 0; iPt < NPts; iPt++) {
        for(auto ic = 0; ic < nAtoms; ic++) {
          GDenX[ic][iPt] = 0.;
          GDenY[ic][iPt] = 0.;
          GDenZ[ic][iPt] = 0.;
          GGDenxX[ic][iPt] = 0.;
          GGDenxY[ic][iPt] = 0.;
          GGDenxZ[ic][iPt] = 0.;
          GGDenyX[ic][iPt] = 0.;
          GGDenyY[ic][iPt] = 0.;
          GGDenyZ[ic][iPt] = 0.;
          GGDenzX[ic][iPt] = 0.;
          GGDenzY[ic][iPt] = 0.;
          GGDenzZ[ic][iPt] = 0.;
          const size_t NBEiPt = iPt*NBE;
          size_t NBEAtom = effOrbsForAtom[ic].size();
          const double *SCR_cur  = SCR2 + NBEiPt;
          //const double *SCR2_curX  = SCR4[ic][0] + NBEiPt;
          //const double *SCR2_curY  = SCR4[ic][1] + NBEiPt;
          //const double *SCR2_curZ  = SCR4[ic][2] + NBEiPt;
          const double *B_cur    = BasisScr + NBEiPt;
          const double *B_curx   = B_cur + IOff;
          const double *B_cury   = B_curx + IOff;
          const double *B_curz   = B_cury + IOff;
          const double *B_curX   = BasisGradScr  + NBEiPt;
          const double *B_curY   = B_curX + IOff;
          const double *B_curZ   = B_curY + IOff;

          // nuclear gradient of orbital gradient
          const double *B_curxX = B_curZ + IOff;
          const double *B_curxY = B_curxX + IOff;
          const double *B_curxZ = B_curxY + IOff;
          const double *B_curyX = B_curxZ + IOff;
          const double *B_curyY = B_curyX + IOff;
          const double *B_curyZ = B_curyY + IOff;
          const double *B_curzX = B_curyZ + IOff;
          const double *B_curzY = B_curzX + IOff;
          const double *B_curzZ = B_curzY + IOff;
          const double *SCR3_curx = SCR5[0] + NBEiPt;
          const double *SCR3_cury = SCR5[1] + NBEiPt;
          const double *SCR3_curz = SCR5[2] + NBEiPt;
      
          for(size_t j = 0; j < NBEAtom; j++) { 
            size_t oi = effOrbsForAtom[ic][j];
            GDenX[ic][iPt] += SCR_cur[oi] * B_curX[oi];
            GDenY[ic][iPt] += SCR_cur[oi] * B_curY[oi];
            GDenZ[ic][iPt] += SCR_cur[oi] * B_curZ[oi];

            GGDenxX[ic][iPt] += SCR_cur[oi] * B_curxX[oi];
            GGDenxX[ic][iPt] += SCR3_curx[oi] * B_curX[oi];

            GGDenxY[ic][iPt] += SCR_cur[oi] * B_curxY[oi];
            GGDenxY[ic][iPt] += SCR3_curx[oi] * B_curY[oi];

            GGDenxZ[ic][iPt] += SCR_cur[oi] * B_curxZ[oi];
            GGDenxZ[ic][iPt] += SCR3_curx[oi] * B_curZ[oi];

            GGDenyX[ic][iPt] += SCR_cur[oi] * B_curyX[oi];
            GGDenyX[ic][iPt] += SCR3_cury[oi] * B_curX[oi];

            GGDenyY[ic][iPt] += SCR_cur[oi] * B_curyY[oi];
            GGDenyY[ic][iPt] += SCR3_cury[oi] * B_curY[oi];

            GGDenyZ[ic][iPt] += SCR_cur[oi] * B_curyZ[oi];
            GGDenyZ[ic][iPt] += SCR3_cury[oi] * B_curZ[oi];

            GGDenzX[ic][iPt] += SCR_cur[oi] * B_curzX[oi];
            GGDenzX[ic][iPt] += SCR3_curz[oi] * B_curX[oi];

            GGDenzY[ic][iPt] += SCR_cur[oi] * B_curzY[oi];
            GGDenzY[ic][iPt] += SCR3_curz[oi] * B_curY[oi];

            GGDenzZ[ic][iPt] += SCR_cur[oi] * B_curzZ[oi];
            GGDenzZ[ic][iPt] += SCR3_curz[oi] * B_curZ[oi];

          }

          //if ( ic == iAtm and this->scfControls.grid_resp ) {
          ////if (false) {
          //  for ( auto& atOrbs: effOrbsForAtom ) {
          //    for ( auto& oi: atOrbs ) {
          //      GDenX[ic][iPt] -= SCR_cur[oi] * B_curX[oi];
          //      GDenY[ic][iPt] -= SCR_cur[oi] * B_curY[oi];
          //      GDenZ[ic][iPt] -= SCR_cur[oi] * B_curZ[oi];

          //      GGDenxX[ic][iPt] -= SCR_cur[oi] * B_curxX[oi];
          //      GGDenxX[ic][iPt] -= SCR3_curx[oi] * B_curX[oi];

          //      GGDenxY[ic][iPt] -= SCR_cur[oi] * B_curxY[oi];
          //      GGDenxY[ic][iPt] -= SCR3_curx[oi] * B_curY[oi];

          //      GGDenxZ[ic][iPt] -= SCR_cur[oi] * B_curxZ[oi];
          //      GGDenxZ[ic][iPt] -= SCR3_curx[oi] * B_curZ[oi];

          //      GGDenyX[ic][iPt] -= SCR_cur[oi] * B_curyX[oi];
          //      GGDenyX[ic][iPt] -= SCR3_cury[oi] * B_curX[oi];

          //      GGDenyY[ic][iPt] -= SCR_cur[oi] * B_curyY[oi];
          //      GGDenyY[ic][iPt] -= SCR3_cury[oi] * B_curY[oi];

          //      GGDenyZ[ic][iPt] -= SCR_cur[oi] * B_curyZ[oi];
          //      GGDenyZ[ic][iPt] -= SCR3_cury[oi] * B_curZ[oi];

          //      GGDenzX[ic][iPt] -= SCR_cur[oi] * B_curzX[oi];
          //      GGDenzX[ic][iPt] -= SCR3_curz[oi] * B_curX[oi];

          //      GGDenzY[ic][iPt] -= SCR_cur[oi] * B_curzY[oi];
          //      GGDenzY[ic][iPt] -= SCR3_curz[oi] * B_curY[oi];

          //      GGDenzZ[ic][iPt] -= SCR_cur[oi] * B_curzZ[oi];
          //      GGDenzZ[ic][iPt] -= SCR3_curz[oi] * B_curZ[oi];

          //    }
          //  }
          //}

          //for(size_t j = 0; j < NBE; j++) { 

          //  GGDenxX[ic][iPt] += SCR2_curX[j] * B_curx[j];
          //  GGDenxY[ic][iPt] += SCR2_curY[j] * B_curx[j];
          //  GGDenxZ[ic][iPt] += SCR2_curZ[j] * B_curx[j];

          //  GGDenyX[ic][iPt] += SCR2_curX[j] * B_cury[j];
          //  GGDenyY[ic][iPt] += SCR2_curY[j] * B_cury[j];
          //  GGDenyZ[ic][iPt] += SCR2_curZ[j] * B_cury[j];

          //  GGDenzX[ic][iPt] += SCR2_curX[j] * B_curz[j];
          //  GGDenzY[ic][iPt] += SCR2_curY[j] * B_curz[j];
          //  GGDenzZ[ic][iPt] += SCR2_curZ[j] * B_curz[j];

          //}
          // Since we are summing over mu and nu
          // Del (mu nu) = 2 * Del(mu) nu 
          GDenX[ic][iPt] *= 2.0;
          GDenY[ic][iPt] *= 2.0;
          GDenZ[ic][iPt] *= 2.0;
          GGDenxX[ic][iPt] *= 2.0;
          GGDenxY[ic][iPt] *= 2.0;
          GGDenxZ[ic][iPt] *= 2.0;
          GGDenyX[ic][iPt] *= 2.0;
          GGDenyY[ic][iPt] *= 2.0;
          GGDenyZ[ic][iPt] *= 2.0;
          GGDenzX[ic][iPt] *= 2.0;
          GGDenzY[ic][iPt] *= 2.0;
          GGDenzZ[ic][iPt] *= 2.0;

          //// add contribution from density matrix gradient
          //for(size_t j = 0; j < NBE; j++) {
          //  GDenX[ic][iPt] += SCR2_curX[j] * B_cur[j]; 
          //  GDenY[ic][iPt] += SCR2_curY[j] * B_cur[j]; 
          //  GDenZ[ic][iPt] += SCR2_curZ[j] * B_cur[j]; 
          //}
        }
      }
    }
  }; // evalDenGrad

  /**
   *  \brief form the grad U variables given the grad V variables.
   *
   *  \param [in]  isGGA                Whether to include GGA contributions
   *  \param [in]  check_aux            Whether the density matrix size check is done for aux system
   *  \param [in]  epsScreen            Screening tolerance
   *  \param [in]  Scalar               Scalar
   *  \param [in]  dndX, dndY, dndZ     Gradient components of n scalar
   *  \param [in]  Mz, My, Mx           Magnetization components
   *  \param [in]  dMkdX, dMkdY, dMkdZ  Gradient components of Mk component of the magnetization
   *  \param [out] ncoll                U variable for density (+,-) for NPts
   *  \param [out] gammaColl            U variable for Gdensity (++,+-,--) for NPts
   *
   */ 
  template <typename MatsT> 
  void mkAuxVarGrad(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM,
    bool isGGA, double epsScreen, size_t NPts,
    std::vector<double*> dScalardX, std::vector<double*> dScalardY, std::vector<double*> dScalardZ,
    std::vector<double*> dMzdX, std::vector<double*> dMzdY, std::vector<double*> dMzdZ,
    double *dScalardx, double *dScalardy, double *dScalardz,
    double *dMzdx, double *dMzdy, double *dMzdz,
    std::vector<double*> dScalardxX, std::vector<double*> dScalardxY, std::vector<double*> dScalardxZ, 
    std::vector<double*> dScalardyX, std::vector<double*> dScalardyY, std::vector<double*> dScalardyZ, 
    std::vector<double*> dScalardzX, std::vector<double*> dScalardzY, std::vector<double*> dScalardzZ, 
    std::vector<double*> dMzdxX, std::vector<double*> dMzdxY, std::vector<double*> dMzdxZ, 
    std::vector<double*> dMzdyX, std::vector<double*> dMzdyY, std::vector<double*> dMzdyZ, 
    std::vector<double*> dMzdzX, std::vector<double*> dMzdzY, std::vector<double*> dMzdzZ, 
    std::vector<double*> nCollGrad_X, std::vector<double*> nCollGrad_Y,std::vector<double*> nCollGrad_Z, 
    std::vector<double*> gammaCollGrad_X, std::vector<double*> gammaCollGrad_Y, std::vector<double*> gammaCollGrad_Z,
    size_t nAtoms){

    double tmp = 0.;
    std::vector<double*> uPlusGrad_X, uPlusGrad_Y, uPlusGrad_Z;
    std::vector<double*> uMinusGrad_X, uMinusGrad_Y, uMinusGrad_Z;
    for(size_t ic = 0; ic < nAtoms; ic++) {
      uPlusGrad_X.emplace_back(nCollGrad_X[ic]);
      uPlusGrad_Y.emplace_back(nCollGrad_Y[ic]);
      uPlusGrad_Z.emplace_back(nCollGrad_Z[ic]);

      uMinusGrad_X.emplace_back(nCollGrad_X[ic]+1);
      uMinusGrad_Y.emplace_back(nCollGrad_Y[ic]+1);
      uMinusGrad_Z.emplace_back(nCollGrad_Z[ic]+1);
    }

    for(size_t ic = 0; ic < nAtoms; ic++) {
      memset(nCollGrad_X[ic],0,2*NPts*sizeof(double));
      memset(nCollGrad_Y[ic],0,2*NPts*sizeof(double));
      memset(nCollGrad_Z[ic],0,2*NPts*sizeof(double));
    }



    // LDA contributions
    // GradU(+) = 0.5 * (Grad_SCALAR + Grad_MZ)
    // GradU(-) = 0.5 * (Grad_SCALAR - Grad_MZ)
         
    for(size_t ic = 0; ic < nAtoms; ic++) {
      blas::axpy(NPts,0.5,dScalardX[ic],1,uPlusGrad_X[ic],2);
      blas::axpy(NPts,0.5,dScalardY[ic],1,uPlusGrad_Y[ic],2);
      blas::axpy(NPts,0.5,dScalardZ[ic],1,uPlusGrad_Z[ic],2);
      blas::axpy(NPts,0.5,dScalardX[ic],1,uMinusGrad_X[ic],2);
      blas::axpy(NPts,0.5,dScalardY[ic],1,uMinusGrad_Y[ic],2);
      blas::axpy(NPts,0.5,dScalardZ[ic],1,uMinusGrad_Z[ic],2);
    }

    bool skipdMzdX = false;
    bool skipdMzdY = false;
    bool skipdMzdZ = false;
    if( onePDM->hasZ() ) {

      for(size_t ic = 0; ic < nAtoms; ic++) {
    // UKS
        if(not skipdMzdX){
          blas::axpy(NPts,0.5,dMzdX[ic],1,uPlusGrad_X[ic],2);
          blas::axpy(NPts,-0.5,dMzdX[ic],1,uMinusGrad_X[ic],2);
        }
        if(not skipdMzdY){
          blas::axpy(NPts,0.5,dMzdY[ic],1,uPlusGrad_Y[ic],2);
          blas::axpy(NPts,-0.5,dMzdY[ic],1,uMinusGrad_Y[ic],2);
        }
        if(not skipdMzdZ){
          blas::axpy(NPts,0.5,dMzdZ[ic],1,uPlusGrad_Z[ic],2);
          blas::axpy(NPts,-0.5,dMzdZ[ic],1,uMinusGrad_Z[ic],2);
        }
      }
    }  else if ( onePDM->hasXY() ) {
      
      CErr("Nuclear Gradients are not implemented for 2-component systems (GHF and X2C)!", std::cout);

    }

    // GGA Contributions
    // GAMMA(++) = 0.25 * (GSCALAR.GSCALAR + GMZ.GMZ) + 0.5 * GSCALAR.GMZ
    // GAMMA(+-) = 0.25 * (GSCALAR.GSCALAR - GMZ.GMZ) 
    // GAMMA(--) = 0.25 * (GSCALAR.GSCALAR + GMZ.GMZ) - 0.5 * GSCALAR.GMZ
    //2C
    // 2C See J. Chem. Theory Comput. 2017, 13, 2591-2603  
    // GAMMA(++) = 0.25 * (GSCALAR.GSCALAR + SUM_K GMK.GMK) + 0.5 * SING * SQRT(SUM_K (GSCALAR.GMK)^2)
    // GAMMA(+-) = 0.25 * (GSCALAR.GSCALAR - SUM_K (GMZ.GMK) ) 
    // GAMMA(--) = 0.25 * (GSCALAR.GSCALAR + SUM_K GMK.GMK) - 0.5 * SIGN * SQRT(SUM_K (GSCALAR.GMK)^2)

    if(isGGA) {
      for(size_t ic = 0; ic < nAtoms; ic++) {
        // RKS part
        for(auto iPt = 0; iPt < NPts; iPt++) {
          gammaCollGrad_X[ic][3*iPt] = 0.5 * (dScalardx[iPt]*dScalardxX[ic][iPt] + dScalardy[iPt]*dScalardyX[ic][iPt] + dScalardz[iPt]*dScalardzX[ic][iPt]);
          gammaCollGrad_X[ic][3*iPt+1] = gammaCollGrad_X[ic][3*iPt]; 
          gammaCollGrad_X[ic][3*iPt+2] = gammaCollGrad_X[ic][3*iPt];

          gammaCollGrad_Y[ic][3*iPt] = 0.5 * (dScalardx[iPt]*dScalardxY[ic][iPt] + dScalardy[iPt]*dScalardyY[ic][iPt] + dScalardz[iPt]*dScalardzY[ic][iPt]);
          gammaCollGrad_Y[ic][3*iPt+1] = gammaCollGrad_Y[ic][3*iPt]; 
          gammaCollGrad_Y[ic][3*iPt+2] = gammaCollGrad_Y[ic][3*iPt];

          gammaCollGrad_Z[ic][3*iPt] = 0.5 * (dScalardx[iPt]*dScalardxZ[ic][iPt] + dScalardy[iPt]*dScalardyZ[ic][iPt] + dScalardz[iPt]*dScalardzZ[ic][iPt]);
          gammaCollGrad_Z[ic][3*iPt+1] = gammaCollGrad_Z[ic][3*iPt]; 
          gammaCollGrad_Z[ic][3*iPt+2] = gammaCollGrad_Z[ic][3*iPt];
        }

        if( onePDM->hasZ() ) {
        // UKS
          for(auto iPt = 0; iPt < NPts; iPt++) {
            if( true ) {

              double innerdX  = 0.5 * 
                (dMzdx[iPt]*dMzdxX[ic][iPt] + dMzdy[iPt]*dMzdyX[ic][iPt] + dMzdz[iPt]*dMzdzX[ic][iPt]);
              double innerdY  = 0.5 * 
                (dMzdx[iPt]*dMzdxY[ic][iPt] + dMzdy[iPt]*dMzdyY[ic][iPt] + dMzdz[iPt]*dMzdzY[ic][iPt]);
              double innerdZ  = 0.5 * 
                (dMzdx[iPt]*dMzdxZ[ic][iPt] + dMzdy[iPt]*dMzdyZ[ic][iPt] + dMzdz[iPt]*dMzdzZ[ic][iPt]);

              double inner2dX = 0.5  * 
                (  dScalardx[iPt]*dMzdxX[ic][iPt] + dMzdx[iPt]*dScalardxX[ic][iPt] 
                 + dScalardy[iPt]*dMzdyX[ic][iPt] + dMzdy[iPt]*dScalardyX[ic][iPt]
                 + dScalardz[iPt]*dMzdzX[ic][iPt] + dMzdz[iPt]*dScalardzX[ic][iPt]              
                );

              double inner2dY = 0.5  * 
                (  dScalardx[iPt]*dMzdxY[ic][iPt] + dMzdx[iPt]*dScalardxY[ic][iPt] 
                 + dScalardy[iPt]*dMzdyY[ic][iPt] + dMzdy[iPt]*dScalardyY[ic][iPt]
                 + dScalardz[iPt]*dMzdzY[ic][iPt] + dMzdz[iPt]*dScalardzY[ic][iPt]              
                );

              double inner2dZ = 0.5  * 
                (  dScalardx[iPt]*dMzdxZ[ic][iPt] + dMzdx[iPt]*dScalardxZ[ic][iPt] 
                 + dScalardy[iPt]*dMzdyZ[ic][iPt] + dMzdy[iPt]*dScalardyZ[ic][iPt]
                 + dScalardz[iPt]*dMzdzZ[ic][iPt] + dMzdz[iPt]*dScalardzZ[ic][iPt]              
                );
        
              gammaCollGrad_X[ic][3*iPt]   += innerdX;
              gammaCollGrad_Y[ic][3*iPt]   += innerdY;
              gammaCollGrad_Z[ic][3*iPt]   += innerdZ;
              gammaCollGrad_X[ic][3*iPt+1] -= innerdX;
              gammaCollGrad_Y[ic][3*iPt+1] -= innerdY;
              gammaCollGrad_Z[ic][3*iPt+1] -= innerdZ;
              gammaCollGrad_X[ic][3*iPt+2] += innerdX;
              gammaCollGrad_Y[ic][3*iPt+2] += innerdY;
              gammaCollGrad_Z[ic][3*iPt+2] += innerdZ;
        
              gammaCollGrad_X[ic][3*iPt]   += inner2dX;
              gammaCollGrad_Y[ic][3*iPt]   += inner2dY;
              gammaCollGrad_Z[ic][3*iPt]   += inner2dZ;
              gammaCollGrad_X[ic][3*iPt+2] -= inner2dX;
              gammaCollGrad_Y[ic][3*iPt+2] -= inner2dY;
              gammaCollGrad_Z[ic][3*iPt+2] -= inner2dZ;
            }
          } // loop pts
        } else if ( onePDM->hasXY()) {

          CErr("Nuclear Gradients are not implemented for 2-component systems (GHF and X2C)!", std::cout);

        } // loop pts
      } // 2C 
    } //GGA 

  }; // mkAuxVarGrad


  template void mkAuxVarGrad(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> onePDM,
    bool isGGA, double epsScreen, size_t NPts,
    std::vector<double*> dScalardX, std::vector<double*> dScalardY, std::vector<double*> dScalardZ,
    std::vector<double*> dMzdX, std::vector<double*> dMzdY, std::vector<double*> dMzdZ,
    double *dScalardx, double *dScalardy, double *dScalardz,
    double *dMzdx, double *dMzdy, double *dMzdz,
    std::vector<double*> dScalardxX, std::vector<double*> dScalardxY, std::vector<double*> dScalardxZ, 
    std::vector<double*> dScalardyX, std::vector<double*> dScalardyY, std::vector<double*> dScalardyZ, 
    std::vector<double*> dScalardzX, std::vector<double*> dScalardzY, std::vector<double*> dScalardzZ, 
    std::vector<double*> dMzdxX, std::vector<double*> dMzdxY, std::vector<double*> dMzdxZ, 
    std::vector<double*> dMzdyX, std::vector<double*> dMzdyY, std::vector<double*> dMzdyZ, 
    std::vector<double*> dMzdzX, std::vector<double*> dMzdzY, std::vector<double*> dMzdzZ, 
    std::vector<double*> nCollGrad_X, std::vector<double*> nCollGrad_Y,std::vector<double*> nCollGrad_Z, 
    std::vector<double*> gammaCollGrad_X, std::vector<double*> gammaCollGrad_Y, std::vector<double*> gammaCollGrad_Z,
    size_t nAtoms);

  template void mkAuxVarGrad(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> onePDM,
    bool isGGA, double epsScreen, size_t NPts,
    std::vector<double*> dScalardX, std::vector<double*> dScalardY, std::vector<double*> dScalardZ,
    std::vector<double*> dMzdX, std::vector<double*> dMzdY, std::vector<double*> dMzdZ,
    double *dScalardx, double *dScalardy, double *dScalardz,
    double *dMzdx, double *dMzdy, double *dMzdz,
    std::vector<double*> dScalardxX, std::vector<double*> dScalardxY, std::vector<double*> dScalardxZ, 
    std::vector<double*> dScalardyX, std::vector<double*> dScalardyY, std::vector<double*> dScalardyZ, 
    std::vector<double*> dScalardzX, std::vector<double*> dScalardzY, std::vector<double*> dScalardzZ, 
    std::vector<double*> dMzdxX, std::vector<double*> dMzdxY, std::vector<double*> dMzdxZ, 
    std::vector<double*> dMzdyX, std::vector<double*> dMzdyY, std::vector<double*> dMzdyZ, 
    std::vector<double*> dMzdzX, std::vector<double*> dMzdzY, std::vector<double*> dMzdzZ, 
    std::vector<double*> nCollGrad_X, std::vector<double*> nCollGrad_Y,std::vector<double*> nCollGrad_Z, 
    std::vector<double*> gammaCollGrad_X, std::vector<double*> gammaCollGrad_Y, std::vector<double*> gammaCollGrad_Z,
    size_t nAtoms);


  double energy_vxc_grad(bool isGGA, size_t NPts, std::vector<double> &weights, double *vrho, double *vsigma, double *GradDen, double *GradGamma){

    double GradXCEnergy = 0.0;

    for(auto iPt = 0; iPt < NPts; iPt++) { 
        GradXCEnergy += weights[iPt] * vrho[2*iPt] * GradDen[2*iPt]; // a
        GradXCEnergy += weights[iPt] * vrho[2*iPt+1] * GradDen[2*iPt+1]; // b
      if (isGGA) {
        GradXCEnergy += weights[iPt] * vsigma[3*iPt] * GradGamma[3*iPt]; // aa
        GradXCEnergy += weights[iPt] * vsigma[3*iPt+1] * GradGamma[3*iPt+1]; // ab
        GradXCEnergy += weights[iPt] * vsigma[3*iPt+2] * GradGamma[3*iPt+2]; // bb
      }
    }
    return GradXCEnergy;

  };// energy_vxc_grad

}
