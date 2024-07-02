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
   *  \brief wrapper for evaluating and loading the DFT 
   *  kernel derivatives from libxc wrt to the U var.
   *
   *  Note. We compute all the quantities for all points in the batch.
   *
   *  \param [in]  NPts       Number of points in the batch
   *  \param [in]  Den        Pointer to the Uvar Density vector[+,-].
   *  \param [in]  Gamma      Pointer to the GGA Uvar vector[++,+-,--]. 
   *  \param [out] epsEval    Pointer to the energy per unit particle. 
   *
   *  \param [out] VRhoEval   Pointer to the first part der of 
   *                          the energy per unit volume in terms of 
   *                          the dens[+,-]. 
   *
   *  \param [out] VgammaEval Pointer to the first part der of 
   *                          the energy per unit volume in terms of 
   *                          the gammas[++,+-,--].
   *
   *  \param [out] EpsSCR     Pointer for multiple functional eval. See epsEval
   *  \param [out] VRhoSCR    Pointer for multiple functional eval. See VRhoEval
   *  \param [out] VgammaSCR  Pointer for multiple functional eval. See VgammaEval
   */  
  void loadVXCder(
    std::vector<std::shared_ptr<DFTFunctional>> functionals,
    size_t NPts, double *Den, double *Gamma,
    double *epsEval, double* VRhoEval, double *VgammaEval,
    double *EpsSCR, double *VRhoSCR, double* VgammaSCR)
  { 

    for(auto iF = 0; iF < functionals.size(); iF++) {
      double *ES,*VR,*VS;
      if( iF == 0 ) {
        ES = epsEval;
        VR = VRhoEval;
        VS = VgammaEval;
      } else {
        ES = EpsSCR;
        VR = VRhoSCR;
        VS = VgammaSCR;
      }

      if( functionals[iF]->isGGA() )
        functionals[iF]->evalEXC_VXC(NPts,Den,Gamma,ES,VR,VS);
      else
        functionals[iF]->evalEXC_VXC(NPts,Den,ES,VR);
      
      if( iF != 0 ) {
        MatAdd('N','N',NPts,1,1.,epsEval,NPts,1.,EpsSCR,NPts,epsEval,NPts);
        MatAdd('N','N',2*NPts,1,1.,VRhoEval,2*NPts,1.,VRhoSCR,2*NPts,VRhoEval,2*NPts);
      if( functionals[iF]->isGGA() )
        MatAdd('N','N',3*NPts,1,1.,VgammaEval,3*NPts,1.,VgammaSCR,3*NPts,VgammaEval,3*NPts);
      }

    }

  };

  /**
   *  \brief form the U variables given the V variables.
   *
   *  \param [in]  isGGA                Whether to include GGA contributions
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
  void mkAuxVar(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM,
    bool isGGA, 
    double epsScreen, size_t NPts,
    double *Scalar, double *Mz, double *My, double *Mx,
    double *dndX, double *dndY, double *dndZ, 
    double *dMzdX, double *dMzdY, double *dMzdZ, 
    double *dMydX, double *dMydY, double *dMydZ, 
    double *dMxdX, double *dMxdY, double *dMxdZ, 
    double *Mnorm, double *Kx, double *Ky, double *Kz, 
    double *Hx, double *Hy, double *Hz,
    double *DSDMnormv, double *signMDv,
    bool *Msmall, double *nColl, double *gammaColl)
  {

    double tmp = 0.;
    double *uPlus  = nColl;
    double *uMinus = nColl + 1;
    memset(nColl,0,2*NPts*sizeof(double));


    // LDA contributions
    // U(+) = 0.5 * (SCALAR + MZ)
    // U(-) = 0.5 * (SCALAR - MZ)
    // 2C See J. Chem. Theory Comput. 2017, 13, 2591-2603  
    // U(+) = 0.5 * (SCALAR + |M| )
    // U(-) = 0.5 * (SCALAR - |M| )
         
    blas::axpy(NPts,0.5,Scalar,1,uPlus,2);
    blas::axpy(NPts,0.5,Scalar,1,uMinus,2);

    bool skipMz = false;
    if( onePDM->hasZ() and not onePDM->hasXY() ) {
    // UKS
#if VXC_DEBUG_LEVEL < 3
      double MaxDenZ = *std::max_element(Mz,Mz+NPts);
      if (MaxDenZ < epsScreen) skipMz = true;
#endif
        if(not skipMz){
          blas::axpy(NPts,0.5,Mz,1,uPlus,2);
          blas::axpy(NPts,-0.5,Mz,1,uMinus,2);
        }
#if VXC_DEBUG_LEVEL >= 3
      if(skipMz) std::cerr << "Skypped Mz " << std::endl;
#endif
    }  else if ( onePDM->hasXY()) {
      std::fill_n(Msmall,NPts,0.);
    // 2C
    // 2C See J. Chem. Theory Comput. 2017, 13, 2591-2603  
     //Compute and store Mtot
      for(auto iPt = 0; iPt < NPts; iPt++) {
        tmp =  Mx[iPt] * Mx[iPt];
        tmp += My[iPt] * My[iPt];
        tmp += Mz[iPt] * Mz[iPt];
        if (tmp > 1.e-24) {

          Mnorm[iPt] = std::sqrt(tmp);
          Kx[iPt]    = Mx[iPt] / Mnorm[iPt];
          Ky[iPt]    = My[iPt] / Mnorm[iPt];
          Kz[iPt]    = Mz[iPt] / Mnorm[iPt];

        } else {

          Msmall[iPt] = true; 
          Mnorm[iPt]  = (1./3.) * (Mx[iPt] + My[iPt] + Mz[iPt]);
          Kx[iPt]    = 1. / 3.;
          Ky[iPt]    = 1. / 3.;
          Kz[iPt]    = 1. / 3.;
        }

      }
          blas::axpy(NPts,0.5,Mnorm,1,uPlus,2);
          blas::axpy(NPts,-0.5,Mnorm,1,uMinus,2);

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
      // RKS part
      for(auto iPt = 0; iPt < NPts; iPt++) {
        gammaColl[3*iPt] = 0.25 *  (dndX[iPt]*dndX[iPt] + dndY[iPt]*dndY[iPt] + dndZ[iPt]*dndZ[iPt]);
        gammaColl[3*iPt+1] = gammaColl[3*iPt]; 
        gammaColl[3*iPt+2] = gammaColl[3*iPt];
      }

      if( onePDM->hasZ() and not onePDM->hasXY() ) {
      // UKS
        for(auto iPt = 0; iPt < NPts; iPt++) {
          if( not skipMz ) {
            double inner  = 0.25 * 
              (dMzdX[iPt]*dMzdX[iPt] + dMzdY[iPt]*dMzdY[iPt] + dMzdZ[iPt]*dMzdZ[iPt]);
            double inner2 = 0.5  * 
              (dMzdX[iPt]*dndX[iPt] + dMzdY[iPt]*dndY[iPt] + dMzdZ[iPt]*dndZ[iPt]);
      
            gammaColl[3*iPt]   += inner;
            gammaColl[3*iPt+1] -= inner;
            gammaColl[3*iPt+2] += inner;
      
            gammaColl[3*iPt]   += inner2;
            gammaColl[3*iPt+2] -= inner2;
          }
        } // loop pts

      }  else if ( onePDM->hasXY()) {
      // 2C See J. Chem. Theory Comput. 2017, 13, 2591-2603  
        double tmpnMx, tmpnMy, tmpnMz;
        double tmpSign, inner, inner2;

        for(auto iPt = 0; iPt < NPts; iPt++) {
          double signMD = 1. ;
          tmpnMx       = dMxdX[iPt]*dndX[iPt];
          tmpnMx      += dMxdY[iPt]*dndY[iPt];
          tmpnMx      += dMxdZ[iPt]*dndZ[iPt];
          
          tmpnMy       = dMydX[iPt]*dndX[iPt];
          tmpnMy      += dMydY[iPt]*dndY[iPt];
          tmpnMy      += dMydZ[iPt]*dndZ[iPt];
          
          tmpnMz       = dMzdX[iPt]*dndX[iPt];
          tmpnMz      += dMzdY[iPt]*dndY[iPt];
          tmpnMz      += dMzdZ[iPt]*dndZ[iPt];
          
          tmpSign   = tmpnMx * Mx[iPt];
          tmpSign  += tmpnMy * My[iPt];
          tmpSign  += tmpnMz * Mz[iPt];
          if ( std::signbit(tmpSign) ) signMD = -1.;

          // Save the value of signMD if signMD is not null pointer
          if (  signMDv ) signMDv[iPt] = signMD;

          inner  = 
               (dMxdX[iPt]*dMxdX[iPt] + dMxdY[iPt]*dMxdY[iPt] + dMxdZ[iPt]*dMxdZ[iPt]);
          inner += 
               (dMydX[iPt]*dMydX[iPt] + dMydY[iPt]*dMydY[iPt] + dMydZ[iPt]*dMydZ[iPt]);
          inner += 
               (dMzdX[iPt]*dMzdX[iPt] + dMzdY[iPt]*dMzdY[iPt] + dMzdZ[iPt]*dMzdZ[iPt]);
          inner *= 0.25;
          
          double DSDMnorm = 0.0;
          if (Msmall[iPt]) {
            DSDMnorm  = (1./3.) * (tmpnMx + tmpnMy +tmpnMz);
            //DSDMnorm[iPt]  = std::sqrt(tmpnMx * tmpnMx + tmpnMy * tmpnMy + tmpnMz * tmpnMz);
            Hx[iPt]        = (1./3.) * signMD ;
            Hy[iPt]        = Hx[iPt] ;
            Hz[iPt]        = Hx[iPt] ;
          } else {
            DSDMnorm  = std::sqrt(tmpnMx * tmpnMx + tmpnMy * tmpnMy + tmpnMz * tmpnMz);
            Hx[iPt]        = (tmpnMx * signMD) / DSDMnorm;
            Hy[iPt]        = (tmpnMy * signMD) / DSDMnorm;
            Hz[iPt]        = (tmpnMz * signMD) / DSDMnorm;
          }
 
          inner2 = 0.5  *  DSDMnorm * signMD ;
          
          gammaColl[3*iPt]   += inner;
          gammaColl[3*iPt+1] -= inner;
          gammaColl[3*iPt+2] += inner;
      
          gammaColl[3*iPt]   += inner2;
          gammaColl[3*iPt+2] -= inner2;

          // If DSDMnormv is not nullptr, save the value 
          if ( DSDMnormv ) {
            DSDMnormv[iPt] = DSDMnorm;
          } 

        } // loop pts
      } // 2C 
    } //GGA 
  };

  template void mkAuxVar(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> onePDM,
    bool isGGA, 
    double epsScreen, size_t NPts,
    double *Scalar, double *Mz, double *My, double *Mx,
    double *dndX, double *dndY, double *dndZ, 
    double *dMzdX, double *dMzdY, double *dMzdZ, 
    double *dMydX, double *dMydY, double *dMydZ, 
    double *dMxdX, double *dMxdY, double *dMxdZ, 
    double *Mnorm, double *Kx, double *Ky, double *Kz, 
    double *Hx, double *Hy, double *Hz,
    double *DSDMnormv, double *signMDv,
    bool *Msmall, double *nColl, double *gammaColl);

  template void mkAuxVar(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> onePDM,
    bool isGGA, 
    double epsScreen, size_t NPts,
    double *Scalar, double *Mz, double *My, double *Mx,
    double *dndX, double *dndY, double *dndZ, 
    double *dMzdX, double *dMzdY, double *dMzdZ, 
    double *dMydX, double *dMydY, double *dMydZ, 
    double *dMxdX, double *dMxdY, double *dMxdZ, 
    double *Mnorm, double *Kx, double *Ky, double *Kz, 
    double *Hx, double *Hy, double *Hz,
    double *DSDMnormv, double *signMDv,
    bool *Msmall, double *nColl, double *gammaColl);


  /**
   *  \brief evaluates the V (Den, GDENX, GDENY and GDENZ)
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
  void evalDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, double *SCR1,
    double *SCR2, double *DENMAT, double *Den, double *GDenX, double *GDenY, double *GDenZ,
    double *BasisScr){

    size_t IOff = NPts*NBE;

    SubMatSet(NB,NB,NBE,NBE,DENMAT,NB,SCR1,NBE,subMatCut);           

    // Obtain Sum_nu P_mu_nu Phi_nu
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,1.,SCR1,NBE,BasisScr,NBE,0.,SCR2,NBE);

    if( typ != GRADIENT )
      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        double *SCR_cur = SCR2 + iPt*NBE;
        double *B_cur   = BasisScr + iPt*NBE;
  
        for (size_t j = 0; j < NBE; j++) {
          Den[iPt] += SCR_cur[j] * B_cur[j];
          //std::cout << "Basis Val " << B_cur[j] << std::endl;
        }

    } else {

      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        GDenX[iPt] = 0.;
        GDenY[iPt] = 0.;
        GDenZ[iPt] = 0.;
        const size_t NBEiPt = iPt*NBE;
        const double *SCR_cur  = SCR2 + NBEiPt;
        const double *B_cur    = BasisScr + NBEiPt;
        const double *B_curX   = B_cur  + IOff;
        const double *B_curY   = B_curX + IOff;
        const double *B_curZ   = B_curY + IOff;
      
        for(size_t j = 0; j < NBE; j++) { 
          Den[iPt]   += SCR_cur[j] * B_cur[j];
          GDenX[iPt] += SCR_cur[j] * B_curX[j];
          GDenY[iPt] += SCR_cur[j] * B_curY[j];
          GDenZ[iPt] += SCR_cur[j] * B_curZ[j];
        }
        // Since we are summing over mu and nu
        // Del (mu nu) = 2 * Del(mu) nu 
        GDenX[iPt] *= 2.0;
        GDenY[iPt] *= 2.0;
        GDenZ[iPt] *= 2.0;
      }
    }
  };
  // SS: evalDen for GIAO 
  void evalDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, dcomplex *SCR1,
    dcomplex *SCR2, dcomplex *DENMAT, double *Den, double *GDenX, double *GDenY, double *GDenZ,
    dcomplex *BasisScr){

    size_t IOff = NPts*NBE;

    SubMatSet(NB,NB,NBE,NBE,DENMAT,NB,SCR1,NBE,subMatCut);           

    // Obtain Sum_nu P^T_mu_nu Phi_nu
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBE,NPts,NBE,dcomplex(1.),SCR1,NBE,BasisScr,NBE,dcomplex(0.),SCR2,NBE);

    if( typ != GRADIENT ) {
      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        dcomplex *SCR_cur = SCR2 + iPt*NBE;
        dcomplex *B_cur   = BasisScr + iPt*NBE;

        dcomplex dentmp;  
        for (size_t j = 0; j < NBE; j++) 
          dentmp += SCR_cur[j] * std::conj(B_cur[j]);

        if (std::abs(dentmp.imag())> 1.0e-13 )
          std::cout<<"imaginary part of density is nonzero "<<std::imag(dentmp)<<std::endl;

        Den[iPt] = std::real(dentmp);

       } // for(auto iPt = 0

    } else {

      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        GDenX[iPt] = 0.;
        GDenY[iPt] = 0.;
        GDenZ[iPt] = 0.;
        const size_t NBEiPt = iPt*NBE;
        const dcomplex *SCR_cur  = SCR2 + NBEiPt;
        const dcomplex *B_cur    = BasisScr + NBEiPt;
        const dcomplex *B_curX   = B_cur  + IOff;
        const dcomplex *B_curY   = B_curX + IOff;
        const dcomplex *B_curZ   = B_curY + IOff;

        dcomplex dentmp=0.0;
        dcomplex denXtmp=0.0;
        dcomplex denYtmp=0.0;
        dcomplex denZtmp=0.0;
              
        for(size_t j = 0; j < NBE; j++) { 
          dentmp  += SCR_cur[j] * std::conj(B_cur[j]);
          denXtmp += SCR_cur[j] * std::conj(B_curX[j]);
          denYtmp += SCR_cur[j] * std::conj(B_curY[j]);
          denZtmp += SCR_cur[j] * std::conj(B_curZ[j]);
        }

 
        dentmp = dentmp + std::conj(dentmp);
        denXtmp = denXtmp + std::conj(denXtmp);
        denYtmp = denYtmp + std::conj(denYtmp);
        denZtmp = denZtmp + std::conj(denZtmp);

        if (std::abs(dentmp.imag())> 1.0e-13 )
          std::cout<<"imaginary part of density is nonzero "<<std::imag(dentmp)<<std::endl;
        if (std::abs(denXtmp.imag())> 1.0e-13 )
          std::cout<<"imaginary part of Nabla x density is nonzero "<<std::imag(denXtmp)<<std::endl;
        if (std::abs(denYtmp.imag())> 1.0e-13 )
          std::cout<<"imaginary part of Nabla y density is nonzero "<<std::imag(denYtmp)<<std::endl;
        if (std::abs(denZtmp.imag())> 1.0e-13 )
          std::cout<<"imaginary part of Nabla z density is nonzero "<<std::imag(denZtmp)<<std::endl;


        Den[iPt]   = 0.5*std::real(dentmp);
        GDenX[iPt] = std::real(denXtmp);
        GDenY[iPt] = std::real(denYtmp);
        GDenZ[iPt] = std::real(denZtmp);

      }
    }
  }; // evalDen
// this is specifically for non-Hermitian density matrix contracted with complex orbitals
  void evalDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, dcomplex *SCR1,
    dcomplex *SCR2, dcomplex *DENMAT, dcomplex *Den, dcomplex *GDenX, dcomplex *GDenY, dcomplex *GDenZ,
    dcomplex *BasisScr){

    size_t IOff = NPts*NBE;

    SubMatSet(NB,NB,NBE,NBE,DENMAT,NB,SCR1,NBE,subMatCut);           

    // Obtain Sum_nu P^T_mu_nu Phi_nu
    blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBE,NPts,NBE,dcomplex(1.),SCR1,NBE,BasisScr,NBE,dcomplex(0.),SCR2,NBE);

    if( typ != GRADIENT ) {
      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        dcomplex *SCR_cur = SCR2 + iPt*NBE;
        dcomplex *B_cur   = BasisScr + iPt*NBE;

        dcomplex dentmp;  
        for (size_t j = 0; j < NBE; j++) 
          dentmp += SCR_cur[j] * std::conj(B_cur[j]);

        //if (std::abs(dentmp.imag())> 1.0e-13 )
        //  std::cout<<"imaginary part of density is nonzero "<<std::imag(dentmp)<<std::endl;

        Den[iPt] = dentmp;

       } // for(auto iPt = 0

    } else {

      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        GDenX[iPt] = 0.;
        GDenY[iPt] = 0.;
        GDenZ[iPt] = 0.;
        const size_t NBEiPt = iPt*NBE;
        const dcomplex *SCR_cur  = SCR2 + NBEiPt;
        const dcomplex *B_cur    = BasisScr + NBEiPt;
        const dcomplex *B_curX   = B_cur  + IOff;
        const dcomplex *B_curY   = B_curX + IOff;
        const dcomplex *B_curZ   = B_curY + IOff;

        dcomplex dentmp=0.0;
        dcomplex denXtmp=0.0;
        dcomplex denYtmp=0.0;
        dcomplex denZtmp=0.0;
              
        for(size_t j = 0; j < NBE; j++) { 
          dentmp  += SCR_cur[j] * std::conj(B_cur[j]);
          denXtmp += SCR_cur[j] * std::conj(B_curX[j]);
          denYtmp += SCR_cur[j] * std::conj(B_curY[j]);
          denZtmp += SCR_cur[j] * std::conj(B_curZ[j]);
        }

        Den[iPt]   = dentmp;
        GDenX[iPt] += denXtmp;
        GDenY[iPt] += denYtmp;
        GDenZ[iPt] += denZtmp;


      } //for(auto iPt = 0; iPt < NPts; iPt++) 

      for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
        BasisScr[icount] = std::conj(BasisScr[icount]);
 
      // Obtain Sum_nu P_mu_nu Phi^*_nu
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,dcomplex(1.),SCR1,NBE,BasisScr,NBE,dcomplex(0.),SCR2,NBE);
 
      for ( int icount = 0 ; icount < NBE*NPts ; icount++ )      
        BasisScr[icount] = std::conj(BasisScr[icount]);

      for(auto iPt = 0; iPt < NPts; iPt++) {
        const size_t NBEiPt = iPt*NBE;
        const dcomplex *SCR_cur  = SCR2 + NBEiPt;
        const dcomplex *B_cur    = BasisScr + NBEiPt;
        const dcomplex *B_curX   = B_cur  + IOff;
        const dcomplex *B_curY   = B_curX + IOff;
        const dcomplex *B_curZ   = B_curY + IOff;

        dcomplex denXtmp=0.0;
        dcomplex denYtmp=0.0;
        dcomplex denZtmp=0.0;
              
        for(size_t j = 0; j < NBE; j++) { 
          denXtmp += SCR_cur[j] * B_curX[j];
          denYtmp += SCR_cur[j] * B_curY[j];
          denZtmp += SCR_cur[j] * B_curZ[j];
        }

        GDenX[iPt] += denXtmp;
        GDenY[iPt] += denYtmp;
        GDenZ[iPt] += denZtmp;

      } //for(auto iPt = 0; iPt < NPts; iPt++) 

    }//if( typ != GRADIENT )
  }; // //KohnSham<dcomplex,dcomplex>::evalDen


  /**
   *  \brief Evaluate the EXC energy 
   *
   *  \param [in] NPts     Number of points in the batch.
   *  \param [in] weights  Quadrature weights.
   *  \param [in] epsEval  Pointer to the energy per unit particle. 
   *  \param [in] DenS     Pointer to the scalar density.
   */  
  double energy_vxc(size_t NPts, std::vector<double> &weights, double *epsEval, double *DenS){

    double XCEnergy = 0.0;

    for(auto iPt = 0; iPt < NPts; iPt++)  
      XCEnergy += weights[iPt] * epsEval[iPt] * DenS[iPt];
     // std::cerr << "XCEnergy " << XCEnergy << std::endl;
    return XCEnergy;

  };

  /**
   *  \brief Construct the required quantities for the formation of the Z vector, 
   *  for a given density component, given the kernel derivatives wrt U variables. 
   *
   *  See J. Chem. Theory Comput. 2011, 7, 3097–3104 (modified e derived for Scalar and Magn).
   *
   *  \param [in] dentyp       Type of 1PDM (SCALAR, {M_k}).
   *  \param [in] isGGA        Whether to include GGA contributions.
   *  \param [in] NPts         Number of points in the batch
   *
   *  \param [in] VRhoEval     Pointer to the first partial derivative of 
   *                           the energy per unit volume in terms of the dens[+,-]. 
   *
   *  \param [in] VgammaEval   Pointer to the first partial derivative of 
   *                           the energy per unit volume in terms of the 
   *                           gammas[++,+-,--].
   *
   *  \param[out]  ZrhoVar1    Factors to multiply the scalar density for particular Z matrix
   *  \param[out]  ZgammaVar1  Factors to multiply the gradient scalar density for particular Z matrix
   *  \param[out]  ZgammaVar2  Factors to multiply the gradient particular mag 
   *                           density for particular Z matrix
   *  
   *  Note we build 2 * X   in Eq 12 and 13 in J. Chem. Theory Comput. 2011, 7, 3097–3104.
   *  Since ZMAT LDA part needs to factor 0.5 for the symmetrization procedure 
   *  (see Eq. 15 1st term on r.h.s) ZrhoVar1  part does not need this factor anymore.
   *
   *  The 0.5 factors come from the chain rules.
   *
   *  On the other hand since we are building 2 * X, we factor already in both
   *  ZgammaVar1 and ZgammaVar2 (since there is 0.5 coming from the chain rules
   *  as well for the GGA terms). Note there is still a factor of 2 that is included
   *  already in the Grad SCALAR/Mz (the one required in Eq 17).
   *
   *  Notes. The ZrhoVar1   multiply the LDA contribution to ZMAT 
   *                        
   *  Notes. The ZgammaVar1 multiply the GGA Del SCALAR
   *                        contribution to ZMAT 
   *  Notes. The ZgammaVar2 multiply the GGA Del Mk 
   *                        contribution to ZMAT 
   *  
   */  
  template <typename MatsT>
  void constructZVars(std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
    double *VrhoEval, double *VgammaEval, double *ZrhoVar1, 
    double *ZgammaVar1, double *ZgammaVar2) {
    // FIXME: Don't zero out, copy / use MKL VAdd
    memset(ZrhoVar1,0,NPts*sizeof(double));
    if (isGGA) { 
      memset(ZgammaVar1,0,NPts*sizeof(double));
      memset(ZgammaVar2,0,NPts*sizeof(double));
    }


    if( denTyp == SCALAR) {
      // (DE/DN+ DN+/DSCALAR + DE/DN- DN-/DSCALAR )
      // Where DN+/DSCALAR = DN+/DSCALAR = 0.5
      blas::axpy(NPts,0.5,VrhoEval,2,ZrhoVar1,1);
      blas::axpy(NPts,0.5 ,VrhoEval+1,2,ZrhoVar1,1);
    } else {
      // (DE/DN+ DN+/DMk + DE/DN- DN-/DMk )
      // Where DN+/DMk =  0.5
      // Where DN-/DMk = -0.5 (only for UKS) 
      blas::axpy(NPts,0.5,VrhoEval,2,ZrhoVar1,1);
      blas::axpy(NPts,-0.5,VrhoEval+1,2,ZrhoVar1,1);
    }


    if(isGGA) {
      if( denTyp == SCALAR ) {

        // ( DE/DGamma++ DGamma++/DSCALAR +
        //   DE/DGamma+- DGamma+-/DSCALAR +
        //   DE/DGamma-- DGamma--/DSCALAR   )

        //   Where DGamma++/DSCALAR = 0.5 * (Del SCAL + Del Mz --- only UKS)
        //   Where DGamma+-/DSCALAR = 0.5 * (Del SCAL)
        //   Where DGamma--/DSCALAR = 0.5 * (Del SCAL - Del Mz --- only UKS)
        //   The Del SCAL and Del Mz will be assembled later in formZ_vxc
        //   NOTE we are building 2 * Z. So 0.5 ---> 1.

        blas::axpy(NPts,1.,VgammaEval,3  ,ZgammaVar1,1);
        blas::axpy(NPts,1.,VgammaEval+1,3,ZgammaVar1,1);
        blas::axpy(NPts,1.,VgammaEval+2,3,ZgammaVar1,1);

        if( onePDM->hasZ() ) {
          blas::axpy(NPts,1.,VgammaEval,3   ,ZgammaVar2,1);
          blas::axpy(NPts,-1.,VgammaEval+2,3,ZgammaVar2,1);
        }
      } else {

        // ( DE/DGamma++ DGamma++/DMz +
        //   DE/DGamma+- DGamma+-/DMz +
        //   DE/DGamma-- DGamma--/DMz   )

        //   Where DGamma++/DMz = 0.5 * (Del Mz + Del SCAL --- only UKS)
        //   Where DGamma+-/DMz = 0.5 * (- Del Mz)
        //   Where DGamma--/DMz = 0.5 * (Del Mz - Del SCAL --- only UKS)
        //   The Del SCAL and Del Mz will be assembled later in formZ_vxc
        //   NOTE we are building 2 * Z. So 0.5 ---> 1.

        blas::axpy(NPts,1.,VgammaEval,3  ,ZgammaVar2,1);
        blas::axpy(NPts,-1.,VgammaEval+1,3,ZgammaVar2,1);
        blas::axpy(NPts,1.,VgammaEval+2,3,ZgammaVar2,1);

        blas::axpy(NPts,1.,VgammaEval,3   ,ZgammaVar1,1);
        blas::axpy(NPts,-1.,VgammaEval+2,3,ZgammaVar1,1);
      }
    }
  };

  template
  void constructZVars(std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
    double *VrhoEval, double *VgammaEval, double *ZrhoVar1, 
    double *ZgammaVar1, double *ZgammaVar2);

  template
  void constructZVars(std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
    double *VrhoEval, double *VgammaEval, double *ZrhoVar1, 
    double *ZgammaVar1, double *ZgammaVar2);


  /**
   *  \brief assemble the final Z vector -> ZMAT (for a given density component),
   *  given the precomputed required U dependent quantities 
   *  (ZUvar# from constructZVars) and the V variables (DenS/Z/X/Y 
   *  and  GDenS/Z/Y/X from evalDen) in input. It requires 
   *  in input also the pointer (BasisScratch) to all basis set (and their gradient)
   *
   *  See J. Chem. Theory Comput. 2011, 7, 3097–3104 (modified e derived for Total and Magn)
   *
   *  \param [in]  isGGA        Whether to include GGA contributions
   *  \param [in]  NPts         Number of points in the batch
   *  \param [in]  NBE          Effective number of basis to be evalauted (only shell actives)
   *  \param [in]  IOff         Offset used for accessing all basis set (and their gradient) 
   *  \param [in]  epsScreen    Screening tollerance
   *  \param [in]  weights      Quadrature weights
   *
   *  \param [in]  ZrhoVar1     Factors to multiply the scalar density for particular Z matrix
   *  \param [in]  ZgammaVar1   Factors to multiply the gradient scalar density for particular Z matrix
   *  \param [in]  ZgammaVar2   Factors to multiply the gradient particular mag 
   *                            density for particular Z matrix
   *
   *  \param [in]  DenS         Pointer to the V variable - SCALAR
   *  \param [in]  DenZ         Pointer to the V variable - Mz 
   *  \param [in]  DenY         Pointer to the V variable - My
   *  \param [in]  DenZ         Pointer to the V variable - Mx
   *  \param [in]  GDenS        Pointer to the V variable - Gradient X,Y,Z comp of SCALAR
   *  \param [in]  GDenZ        Pointer to the V variable - Gradient X,Y,Z comp of Mz 
   *  \param [in]  GDenY        Pointer to the V variable - Gradient X,Y,Z comp of My
   *  \param [in]  GDenZ        Pointer to the V variable - Gradient X,Y,Z comp of Mx
   *
   *  \param [in]  BasisScratch Pointer to Basis set evaluated over batch of points.
   *
   *  \param [out] Pointer to the ZMAT.
   *   
   *  Note. See Documentations of constructZVars.
   */  
  template <typename MatsT, typename IntsT>
  void formZ_vxc(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, size_t NBE, size_t IOff, 
    double epsScreen, std::vector<double> &weights, double *ZrhoVar1,
    double *ZgammaVar1, double *ZgammaVar2, 
    double *DenS, double *DenZ, double *DenY, double *DenX, 
    double *GDenS, double *GDenZ, double *GDenY, double *GDenX, 
    double *Kx, double *Ky, double *Kz,
    double *Hx, double *Hy, double *Hz,
    IntsT *BasisScratch, IntsT *ZMAT) {

    IntsT Fg;
    // Fx,y,z  ^m(all batch) in J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 17 (changed for Total and Magn)
    IntsT FgX;
    IntsT FgY;
    IntsT FgZ;

    memset(ZMAT,0,IOff*sizeof(IntsT));

    if( not onePDM->hasXY() ) {
      for(auto iPt = 0; iPt < NPts; iPt++) { 
      // LDA part -> Eq. 15 and 16 (see constructZVars docs for the missing factor of 0.5)
        Fg = weights[iPt] * ZrhoVar1[iPt];

#if VXC_DEBUG_LEVEL < 3
        if(std::abs(Fg) > epsScreen)
#endif
          blas::axpy(NBE,Fg,BasisScratch + iPt*NBE,1,ZMAT+iPt*NBE,1);

      // GGA part -> Eq. 15 and 17 (see constructZVars docs for the missing factor of 2)
        if( isGGA ) {
          FgX = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt];
          if( onePDM->hasZ() )
            FgX += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt];

          FgY = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + NPts];
          if( onePDM->hasZ() )
            FgY += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + NPts];

          FgZ = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + 2*NPts];
          if( onePDM->hasZ() )
            FgZ += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + 2*NPts];

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(FgX) > epsScreen)
#endif
            blas::axpy(NBE,FgX,BasisScratch + iPt*NBE + IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(FgY) > epsScreen)
#endif
            blas::axpy(NBE,FgY,BasisScratch + iPt*NBE + 2*IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(FgZ) > epsScreen)
#endif
            blas::axpy(NBE,FgZ,BasisScratch + iPt*NBE + 3*IOff,1,ZMAT+iPt*NBE,1);
        }
      }

    } else {
      //2C
      // 2C See J. Chem. Theory Comput. 2017, 13, 2591-2603  
      if( denTyp == SCALAR) {
        // SCALAR

        for(auto iPt = 0; iPt < NPts; iPt++) { 
      // LDA part -> Eq. 15 and 16 (see constructZVars docs for the missing factor of 0.5)
          Fg = weights[iPt] * ZrhoVar1[iPt];

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(Fg) > epsScreen)
#endif
            blas::axpy(NBE,Fg,BasisScratch + iPt*NBE,1,ZMAT+iPt*NBE,1);

        // GGA part -> Eq. 15 and 17 (see constructZVars docs for the missing factor of 2)
          if( isGGA ) {
            FgX  = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * Hx[iPt] * GDenX[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * Hy[iPt] * GDenY[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * Hz[iPt] * GDenZ[iPt];
  
            FgY  = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * Hx[iPt] * GDenX[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * Hy[iPt] * GDenY[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * Hz[iPt] * GDenZ[iPt + NPts];
  
            FgZ  = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * Hx[iPt] * GDenX[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * Hy[iPt] * GDenY[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * Hz[iPt] * GDenZ[iPt + 2*NPts];

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgX) > epsScreen)
#endif
              blas::axpy(NBE,FgX,BasisScratch + iPt*NBE + IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgY) > epsScreen)
#endif
              blas::axpy(NBE,FgY,BasisScratch + iPt*NBE + 2*IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgZ) > epsScreen)
#endif
              blas::axpy(NBE,FgZ,BasisScratch + iPt*NBE + 3*IOff,1,ZMAT+iPt*NBE,1);
          }
        } // loop over Pts

      } else if (denTyp == MX) {

        // MX
        for(auto iPt = 0; iPt < NPts; iPt++) { 
      // LDA part -> Eq. 15 and 16 (see constructZVars docs for the missing factor of 0.5)
          Fg = weights[iPt] * ZrhoVar1[iPt] * Kx[iPt];

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(Fg) > epsScreen)
#endif
            blas::axpy(NBE,Fg,BasisScratch + iPt*NBE,1,ZMAT+iPt*NBE,1);

        // GGA part -> Eq. 15 and 17 (see constructZVars docs for the missing factor of 2)
          if( isGGA ) {
            FgX  = weights[iPt] * ZgammaVar1[iPt] * Hx[iPt] * GDenS[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * GDenX[iPt];
  
            FgY  = weights[iPt] * ZgammaVar1[iPt] * Hx[iPt] * GDenS[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * GDenX[iPt + NPts];
  
            FgZ  = weights[iPt] * ZgammaVar1[iPt] * Hx[iPt] * GDenS[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * GDenX[iPt + 2*NPts];

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgX) > epsScreen)
#endif
              blas::axpy(NBE,FgX,BasisScratch + iPt*NBE + IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgY) > epsScreen)
#endif
              blas::axpy(NBE,FgY,BasisScratch + iPt*NBE + 2*IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgZ) > epsScreen)
#endif
              blas::axpy(NBE,FgZ,BasisScratch + iPt*NBE + 3*IOff,1,ZMAT+iPt*NBE,1);
          }
        } // loop over Pts

      } else if (denTyp == MY) {

        // MY
        for(auto iPt = 0; iPt < NPts; iPt++) { 
      // LDA part -> Eq. 15 and 16 (see constructZVars docs for the missing factor of 0.5)
          Fg = weights[iPt] * ZrhoVar1[iPt] * Ky[iPt];

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(Fg) > epsScreen)
#endif
            blas::axpy(NBE,Fg,BasisScratch + iPt*NBE,1,ZMAT+iPt*NBE,1);

        // GGA part -> Eq. 15 and 17 (see constructZVars docs for the missing factor of 2)
          if( isGGA ) {
            FgX  = weights[iPt] * ZgammaVar1[iPt] * Hy[iPt] * GDenS[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * GDenY[iPt];
  
            FgY  = weights[iPt] * ZgammaVar1[iPt] * Hy[iPt] * GDenS[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * GDenY[iPt + NPts];
  
            FgZ  = weights[iPt] * ZgammaVar1[iPt] * Hy[iPt] * GDenS[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * GDenY[iPt + 2*NPts];

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgX) > epsScreen)
#endif
              blas::axpy(NBE,FgX,BasisScratch + iPt*NBE + IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgY) > epsScreen)
#endif
              blas::axpy(NBE,FgY,BasisScratch + iPt*NBE + 2*IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgZ) > epsScreen)
#endif
              blas::axpy(NBE,FgZ,BasisScratch + iPt*NBE + 3*IOff,1,ZMAT+iPt*NBE,1);
          }
        } // loop over Pts

      } else if (denTyp == MZ) {

        // MZ
        for(auto iPt = 0; iPt < NPts; iPt++) { 
      // LDA part -> Eq. 15 and 16 (see constructZVars docs for the missing factor of 0.5)
          Fg = weights[iPt] * ZrhoVar1[iPt] * Kz[iPt];

#if VXC_DEBUG_LEVEL < 3
          if(std::abs(Fg) > epsScreen)
#endif
            blas::axpy(NBE,Fg,BasisScratch + iPt*NBE,1,ZMAT+iPt*NBE,1);

        // GGA part -> Eq. 15 and 17 (see constructZVars docs for the missing factor of 2)
          if( isGGA ) {
            FgX  = weights[iPt] * ZgammaVar1[iPt] * Hz[iPt] * GDenS[iPt];
            FgX += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt];
  
            FgY  = weights[iPt] * ZgammaVar1[iPt] * Hz[iPt] * GDenS[iPt + NPts];
            FgY += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + NPts];
  
            FgZ  = weights[iPt] * ZgammaVar1[iPt] * Hz[iPt] * GDenS[iPt + 2*NPts];
            FgZ += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + 2*NPts];

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgX) > epsScreen)
#endif
              blas::axpy(NBE,FgX,BasisScratch + iPt*NBE + IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgY) > epsScreen)
#endif
              blas::axpy(NBE,FgY,BasisScratch + iPt*NBE + 2*IOff,1,ZMAT+iPt*NBE,1);

#if VXC_DEBUG_LEVEL < 3
            if(std::abs(FgZ) > epsScreen)
#endif
              blas::axpy(NBE,FgZ,BasisScratch + iPt*NBE + 3*IOff,1,ZMAT+iPt*NBE,1);
          }  
        } // loop over Pts

      } //MZ

    } //end 2C

  }; // KohnSham<MatsT,IntsT>::formZ_vxc

  template
  void formZ_vxc(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, size_t NBE, size_t IOff, 
    double epsScreen, std::vector<double> &weights, double *ZrhoVar1,
    double *ZgammaVar1, double *ZgammaVar2, 
    double *DenS, double *DenZ, double *DenY, double *DenX, 
    double *GDenS, double *GDenZ, double *GDenY, double *GDenX, 
    double *Kx, double *Ky, double *Kz,
    double *Hx, double *Hy, double *Hz,
    double *BasisScratch, double *ZMAT);

  template
  void formZ_vxc(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, size_t NBE, size_t IOff, 
    double epsScreen, std::vector<double> &weights, double *ZrhoVar1,
    double *ZgammaVar1, double *ZgammaVar2, 
    double *DenS, double *DenZ, double *DenY, double *DenX, 
    double *GDenS, double *GDenZ, double *GDenY, double *GDenX, 
    double *Kx, double *Ky, double *Kz,
    double *Hx, double *Hy, double *Hz,
    double *BasisScratch, double *ZMAT);

  template
  void formZ_vxc(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> onePDM,
    DENSITY_TYPE denTyp, bool isGGA, size_t NPts, size_t NBE, size_t IOff, 
    double epsScreen, std::vector<double> &weights, double *ZrhoVar1,
    double *ZgammaVar1, double *ZgammaVar2, 
    double *DenS, double *DenZ, double *DenY, double *DenX, 
    double *GDenS, double *GDenZ, double *GDenY, double *GDenX, 
    double *Kx, double *Ky, double *Kz,
    double *Hx, double *Hy, double *Hz,
    dcomplex *BasisScratch, dcomplex *ZMAT);

}
