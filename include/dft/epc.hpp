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

#include <dft.hpp>
#include <matrix.hpp>
#include <quantum/base.hpp>

namespace ChronusQ {
  /**
   *  \brief EPC-17 functional
   */
  class EPC17 : public DFTFunctional {

  public:
    EPC17(std::string epcName) : DFTFunctional(epcName) 
    { 
      this->isGGA_ = false; 
      this->isEPC_ = true; 
    };

    void evalEXC_VXC(size_t N, double *rho, double *aux_rho, double *eps, double *vxc, bool electron);
    
    void evalEPCGrad(size_t N, double *rho, double *aux_rho, double *eps, double *aux_eps, double* dede, double* dedp);

  }; // epc17-2 functional


  /**
   *  \brief EPC-19 functional
   *
   *  Ref: J. Chem. Phys. 151, 124102 (2019)
   */
  class EPC19 : public DFTFunctional {

  public:
    EPC19(std::string epcName) : DFTFunctional(epcName) 
    { 
      this->isGGA_ = true; 
      this->isEPC_ = true;
    };

    void evalEXC_VXC(size_t N, double *rho, double *aux_rho, double *sigma, double *aux_sigma,
                     double *cross_sigma, double *eps, double *vrho, double *vsigma, 
                     double *vcsigma, double* epc, bool electron);
  }; // class EPC19

  /**
   *  \brief form the U variables given the V variables.
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
  void mkCrossAuxVar(
    bool check_aux, bool electron,
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM1,
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM2,
    double epsScreen, size_t NPts,
    double *dndX1, double *dndY1, double *dndZ1, 
    double *dMxdX1, double *dMxdY1, double *dMxdZ1, 
    double *dMydX1, double *dMydY1, double *dMydZ1, 
    double *dMzdX1, double *dMzdY1, double *dMzdZ1, 
    double *dndX2, double  *dndY2, double *dndZ2, 
    double *dMxdX2, double *dMxdY2, double *dMxdZ2, 
    double *dMydX2, double *dMydY2, double *dMydZ2, 
    double *dMzdX2, double *dMzdY2, double *dMzdZ2, 
    double *gammaColl){

    // GGA Contributions
    // GAMMA(++) = 0.25 * (GSCALAR(e).GSCALAR(p) + GMZ(e).GMZ(p) + GSCALAR(e).GZ(p) + GZ(e).GSCALAR(p))
    // GAMMA(+-) = 0.25 * (GSCALAR.GSCALAR - GMZ.GMZ) 
    // GAMMA(-+) = 
    // GAMMA(--) = 0.25 * (GSCALAR.GSCALAR + GMZ.GMZ) - 0.5 * GSCALAR.GMZ
    //2C : not yet implemented

    // RKS part
    for(auto iPt = 0; iPt < NPts; iPt++) {
      gammaColl[4*iPt] = 0.25 * (dndX1[iPt]*dndX2[iPt] + dndY1[iPt]*dndY2[iPt] + dndZ1[iPt]*dndZ2[iPt]);
      gammaColl[4*iPt+1] = gammaColl[4*iPt];
      gammaColl[4*iPt+2] = gammaColl[4*iPt];
      gammaColl[4*iPt+3] = gammaColl[4*iPt];
    }

    // main system is electron
    if (electron) {

      // restricted
      for(auto iPt = 0; iPt < NPts; iPt++) {
        double gammaTmp = 0.25 * (dndX1[iPt]*dMzdX2[iPt] + dndY1[iPt]*dMzdY2[iPt] + dndZ1[iPt]*dMzdZ2[iPt]);
        gammaColl[4*iPt]   += gammaTmp;
        gammaColl[4*iPt+1] -= gammaTmp;
        gammaColl[4*iPt+2] += gammaTmp;
        gammaColl[4*iPt+3] -= gammaTmp;
      }

      // UKS
      if( onePDM1->hasZ() ) {
        for(auto iPt = 0; iPt < NPts; iPt++) {
          double gammaTmp1 = 0.25 * (dMzdX1[iPt]*dMzdX2[iPt] + dMzdY1[iPt]*dMzdY2[iPt] + dMzdZ1[iPt]*dMzdZ2[iPt]);
          double gammaTmp2 = 0.25 * (dMzdX1[iPt]*dndX2[iPt] + dMzdY1[iPt]*dndY2[iPt] + dMzdZ1[iPt]*dndZ2[iPt]);

          // aa
          gammaColl[4*iPt]  += gammaTmp1;
          gammaColl[4*iPt]  += gammaTmp2;
          // ab
          gammaColl[4*iPt+1]  -= gammaTmp1;
          gammaColl[4*iPt+1]  += gammaTmp2;
          // ba
          gammaColl[4*iPt+2]  -= gammaTmp1;
          gammaColl[4*iPt+2]  -= gammaTmp2;
          // bb
          gammaColl[4*iPt+3]  += gammaTmp1;
          gammaColl[4*iPt+3]  -= gammaTmp2;

        }
      }
    }
    else {  // main system is proton

      // restricted
      for(auto iPt = 0; iPt < NPts; iPt++) {
        double gammaTmp = 0.25 * (dMzdX1[iPt]*dndX2[iPt] + dMzdY1[iPt]*dndY2[iPt] + dMzdZ1[iPt]*dndZ2[iPt]);

        gammaColl[4*iPt]   += gammaTmp;
        gammaColl[4*iPt+1] += gammaTmp;
        gammaColl[4*iPt+2] -= gammaTmp;
        gammaColl[4*iPt+3] -= gammaTmp;

      }

      // UKS
      if( onePDM2->hasZ() ) {
        for(auto iPt = 0; iPt < NPts; iPt++) {
          double gammaTmp1 = 0.25 * (dMzdX1[iPt]*dMzdX2[iPt] + dMzdY1[iPt]*dMzdY2[iPt] + dMzdZ1[iPt]*dMzdZ2[iPt]);
          double gammaTmp2 = 0.25 * (dndX1[iPt]*dMzdX2[iPt] + dndY1[iPt]*dMzdY2[iPt] + dndZ1[iPt]*dMzdZ2[iPt]);
          // aa
          gammaColl[4*iPt]  += gammaTmp1;
          gammaColl[4*iPt]  += gammaTmp2;
          // ab
          gammaColl[4*iPt+1]  -= gammaTmp1;
          gammaColl[4*iPt+1]  -= gammaTmp2;
          // ba
          gammaColl[4*iPt+2]  -= gammaTmp1;
          gammaColl[4*iPt+2]  += gammaTmp2;
          // bb
          gammaColl[4*iPt+3]  += gammaTmp1;
          gammaColl[4*iPt+3]  -= gammaTmp2;
        }
      }
    }
  }; //KohnSham<MatsT,IntsT>::mkCrossAuxVar

  /**
   *  \brief wrapper for evaluating and loading the NEO-DFT 
   *  kernel derivatives from libxc wrt to the U var.
   *
   *  Note. We compute all the quantities for all points in the batch.
   *
   *  \param [in]  NPts       Number of points in the batch
   *  \param [in]  Den        Pointer to the Uvar Density vector[+,-].
   *  \param [in]  Gamma      Pointer to the GGA Uvar vector[++,+-,--]. 
   *  \param [in]  aux_Den    Pointer to the auxiliary Uvar Density vector[+,-].
   *  \param [in]  aux_Gamma  Pointer to the auxiliary GGA Uvar vector[++,+-,--]. 
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
   *  \param [out] CVgammaEval Pointer to the first part der of 
   *                          the energy per unit volume in terms of 
   *                          the gammas[++,+-,--].
   *
   *  \param [out] EpsSCR     Pointer for multiple functional eval. See epsEval
   *  \param [out] VRhoSCR    Pointer for multiple functional eval. See VRhoEval
   *  \param [out] VgammaSCR  Pointer for multiple functional eval. See VgammaEval
   */  
  void loadEPCVXCder(
    bool electron,
    std::vector<std::shared_ptr<DFTFunctional>> functionals,
    size_t NPts, double *Den1, double *Gamma1,
    double *Den2, double *Gamma2, double *cGamma, double *epsEval, double *VRhoEval, 
    double *VgammaEval,double *CVgammaEval, double *EpsSCR, double *VRhoSCR, 
    double *VgammaSCR, double *CVgammaSCR, double* epcEval);

  void loadEPCGradder(
    std::vector<std::shared_ptr<DFTFunctional>> functionals,
    size_t NPts, double *Den1, double *Den2, double *epsEval, double *epsEval2, 
    double* dede, double* dedp);

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
   void constructEPCZVars(bool electron,
    DENSITY_TYPE denTyp, size_t NPts, 
    double *CVgammaEval, double *ZgammaVar3);

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
  template <typename MatsT>
  void formZ_vxc_epc(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM1,
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM2,
    DENSITY_TYPE denTyp, 
    bool isGGA, size_t NPts, size_t NBE, size_t IOff, 
    double epsScreen, std::vector<double> &weights, double *ZrhoVar1,
    double *ZgammaVar1, double *ZgammaVar2, double *ZgammaVar3,
    double *DenS, double *DenZ, double *DenY, double *DenX, 
    double *GDenS, double *GDenZ, double *GDenY, double *GDenX, 
    double *aux_DenS,  double *aux_DenZ,  double *aux_DenY,  double *aux_DenX, 
    double *aux_GDenS, double *aux_GDenZ, double *aux_GDenY, double *aux_GDenX, 
    double *BasisScratch, double *ZMAT){

    double Fg;
    // Fx,y,z  ^m(all batch) in J. Chem. Theory Comput. 2011, 7, 3097–3104 Eq. 17 (changed for Total and Magn)
    double FgX;
    double FgY;
    double FgZ;

    memset(ZMAT,0,IOff*sizeof(double));

    if( not onePDM1->hasXY() ) {
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
          // add contribution of the epc part 
          FgX += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenS[iPt];

          if( onePDM1->hasZ() )
            FgX += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt];

          // add contribution of the epc part 
          if( onePDM2->hasZ() )
            FgX += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenZ[iPt];

          FgY = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + NPts];
          // add contribution of the epc part 
          FgY += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenS[iPt + NPts];
          if( onePDM1->hasZ() )
            FgY += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + NPts];

          // add contribution of the epc part 
          if( onePDM2->hasZ() )
            FgY += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenZ[iPt + NPts];

          FgZ = weights[iPt] * ZgammaVar1[iPt] * GDenS[iPt + 2*NPts];
          // add contribution of the epc part 
          FgZ += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenS[iPt + 2*NPts];
          if( onePDM1->hasZ() )
            FgZ += weights[iPt] * ZgammaVar2[iPt] * GDenZ[iPt + 2*NPts];

          // add contribution of the epc part 
          if( onePDM2->hasZ() )
            FgZ += 0.5 * weights[iPt] * ZgammaVar3[iPt] * aux_GDenZ[iPt + 2*NPts];

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
      CErr("EPC-19 functional is not implemented for 2-component systems (GHF or X2C)!", std::cout);
    } // 2C

  }; // formZ_vxc_epc


} // namespace ChronusQ
