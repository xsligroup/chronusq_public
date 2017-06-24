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
#include <basisset/basisset_util.hpp>


/*
 * This file is for shared _pure functions_ that can be used as building blocks
 *   for the main evaluation of EXC/VXC/FXC for KS calcs
 */

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
    double *EpsSCR, double *VRhoSCR, double* VgammaSCR);

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
    double *BasisScr);
  //SS: for non-hermitian density matrix contract with complex orbital
  void evalDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, dcomplex *SCR1,
    dcomplex *SCR2, dcomplex *DENMAT, double *Den, double *GDenX, double *GDenY, double *GDenZ,
    dcomplex *BasisScr);
  void evalDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, dcomplex *SCR1,
    dcomplex *SCR2, dcomplex *DENMAT, dcomplex *Den, dcomplex *GDenX, dcomplex *GDenY, dcomplex *GDenZ,
    dcomplex *BasisScr);


  /**
   *  \brief Evaluate the EXC energy 
   *
   *  \param [in] NPts     Number of points in the batch.
   *  \param [in] weights  Quadrature weights.
   *  \param [in] epsEval  Pointer to the energy per unit particle. 
   *  \param [in] DenS     Pointer to the scalar density.
   */  
  double energy_vxc(size_t NPts, std::vector<double> &weights, double *epsEval, double *DenS);


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
    IntsT *BasisScratch, IntsT *ZMAT);

  // Evaluate Density Nuclear Gradient
  void evalDenGrad(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB,
    std::vector<std::pair<size_t,size_t>> &subMatCut, double* SCR1,
    double *SCR2, std::vector<std::vector<double*>>SCR3, std::vector<std::vector<double*>>SCR4, 
    std::vector<double*> SCR5,
    std::vector<std::vector<double*>>GDENMAT, double *DENMAT,  
    std::vector<double*>GDenX, std::vector<double*>GDenY, std::vector<double*>GDenZ,
    std::vector<double*>GGDenxX, std::vector<double*>GGDenxY, std::vector<double*>GGDenxZ, 
    std::vector<double*>GGDenyX, std::vector<double*>GGDenyY, std::vector<double*>GGDenyZ, 
    std::vector<double*>GGDenzX, std::vector<double*>GGDenzY, std::vector<double*>GGDenzZ, 
    double *BasisScr, double *BasisGradScr, size_t nAtoms, BasisSet &basisSet);

  template <typename MatsT>
  void mkAuxVarGrad(
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDM, 
    bool isGGA, double epsScreen, size_t NPts_Batch, 
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
    std::vector<double*> nCollGrad_X, std::vector<double*> nCollGrad_Y,
    std::vector<double*> nCollGrad_Z, 
    std::vector<double*> gammaCollGrad_X, std::vector<double*> gammaCollGrad_Y, 
    std::vector<double*> gammaCollGrad_Z, size_t nAtoms);

  double energy_vxc_grad(bool isGGA, size_t NPts, std::vector<double> &weights, double *vrho,
    double *vsigma, double *GradDen, double *GradGamma);
}
