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

#include <fockbuilder/fourcompfock.hpp>
#include <cqlinalg.hpp>
#include <physcon.hpp>
#include <util/matout.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>

//#define _PRINT_MATRICES

namespace ChronusQ {

  /**   
   *  \brief Forms the 4C GD in Bathes
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formRawGDInBatches(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & onePDMs, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & coulombMatrices, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & exchangeMatrices,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & twoeHs) {
    
      // for incore just wrap around the old loop for now
      //
      // coulombMatrices and exchange Matrices won't be used
      //
      // because coulombMatrix is cqmatrix::Matrix instead of cqmatrix::PauliSpinorMatrices
      // so the dividing is not accurate 
      if( std::dynamic_pointer_cast<InCore4indexRelERIContraction<MatsT,IntsT>>(ss.TPI) ) {
        
        // cache the ss pointers
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> ss1PDM 
          = increment ? ss.deltaOnePDM: ss.onePDM;
        std::shared_ptr<cqmatrix::Matrix<MatsT>> ssCoulombMatrix = ss.coulombMatrix;
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> ssExchangeMatrix = ss.exchangeMatrix;
        std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> ssTwoeH = ss.twoeH;
        
        // allocate scratch space for coulombMatrix and exchangeMatrix 
        ss.coulombMatrix = std::make_shared<cqmatrix::Matrix<MatsT>>(ss.coulombMatrix->dimension());
        ss.exchangeMatrix = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(ss.exchangeMatrix->dimension());

        for (auto i = 0ul; i < onePDMs.size(); i++) {
          
          // update pointers
          if (increment) ss.deltaOnePDM = onePDMs[i];
          else ss.onePDM = onePDMs[i];
          
          ss.twoeH = twoeHs[i];
          formGDInCore(ss, pert, increment, xHFX, HerDen);
        
        }

        // revert ss pointers
        if (increment) ss.deltaOnePDM = ss1PDM;
        else ss.onePDM = ss1PDM;
        ss.coulombMatrix = ssCoulombMatrix;
        ss.exchangeMatrix = ssExchangeMatrix;
        ss.twoeH = ssTwoeH;
      
      } else if( std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI) ) {
        
        formRawGDInBatchesDirect(ss, pert, increment, xHFX, HerDen, 
          onePDMs, coulombMatrices, exchangeMatrices, twoeHs);
      
      } else
        CErr("Unsupported ERIContraction type.");

  };


  /**   
   *  \brief Forms the 4C Fock matrix using AO-direct
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formRawGDInBatchesDirect(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & onePDMs, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & coulombMatrices, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & exchangeMatrices,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & twoeHs) {
    
    // disable libint2
    if (not this->hamiltonianOptions_.Libcint) CErr("4C Integrals Needs Libcint");
    
    GTODirectRelERIContraction<MatsT,IntsT> &relERICon =
        *std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI);

    bool computeCoulomb  = coulombMatrices.size() > 0;
    bool computeExchange = (std::abs(xHFX) > 1e-12) and exchangeMatrices.size() > 0;
    bool computeTwoeHs   = twoeHs.size() > 0;
    
    // disable exchange as NYI
    if (not computeCoulomb and not computeExchange and not computeTwoeHs) {
     CErr("Nothing specified to compute in FockBuilder::formRawGDInBatches");
    } else if(not computeCoulomb and not computeTwoeHs) {
     CErr("Only computeExchange is not supported in FockBuilder::formRawGDInBatches");
    } else if (computeExchange) {
     CErr("computeExchange NYI in FockBuilder::formRawGDInBatches");
    }
    
    auto & coulombContainers = computeCoulomb ? coulombMatrices: twoeHs;
    
    size_t mPDM  = onePDMs.size();
    size_t NB1C  = ss.basisSet().nBasis;
    size_t NB2C  = 2 * NB1C; 
    size_t NB4C  = 4 * NB1C;
    size_t NB1C2 = NB1C*NB1C;
    size_t NB1C4 = NB1C*NB1C*NB1C*NB1C;
    size_t NB1C3 = NB1C*NB1C*NB1C;

    size_t SS = NB2C*NB1C+NB1C;
    size_t LS = NB2C*NB1C;
    size_t SL = NB1C;

    size_t mpiRank   = MPIRank(ss.comm);
    bool   isNotRoot = mpiRank != 0;
    
    // allocate scratch spaces for Coulomb-type of contraction part
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>>
      contractSymm1PDMLLMS, contractSymm1PDMSS, contract1PDMLSpmSL,  
      CScrLLMS, CScrSS, CScrLS; 

    //std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>>
    //  contract1PDMLL, contract1PDMSS, contract1PDMLS, contract1PDMSL, 
    //  XScrLL, XScrSS, XScrLS, XScrSL; 
     
    #define ALLOCATE_PAULISPINOR_SCR(SCR, SCRSIZE, hasXYZ) \
       SCR.push_back(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(SCRSIZE, hasXYZ, hasXYZ)); 
       // no need to intailize those matrices as it will be initialize in twoBodyRelContract
       // SCR.back()->clear();
    
    /* 
     * SCR Usage for different hamiltonian Options:
     *   
     *   Coulomb Terms and eXchange Terms
     *
     * 1. Bare Coulomb:
     *    - C: symmetrized contract1PDMLLMS, CScrLLMS 
     *    - X: contract1PDMLL(only MS), XScrLL 
     * 2. Dirac Coulomb without SSSS:
     *    - C: symmetrized contract1PDMSS, CScrLLMS 
     *         symmetrized contract1PDMLLMS, CScrSS
     *    - X: contract1PDMLS, XScrLS  
     *         contract1PDMSL, XScrSL    // non-Hermitian ?? 
     * 3. Dirac Coulomb SSSS:
     *    - C: symmetrized contract1PDMSS, CScrSS   
     *    - X: contract1PDMSS, XScrSS  
     * 4. Gaunt:
     *    - C: contract1PDMLS +/- contract1PDMSL, CScrLS 
     *    - X: contract1PDMLL, XScrLL 
     *         contract1PDMSS, XScrSS  
     *         contract1PDMLS, XScrLS
     *         contract1PDMSL, XScrSL
     * 5. Gauge:
     *    - C: contract1PDMLS +/- contract1PDMSL, CScrLS
     *         contract1PDMSL, CScrSL
     *    - X: contract1PDMLL, XScrLL
     *         contract1PDMSS, XScrSS
     *         contract1PDMLS, XScrLS
     *         contract1PDMSL, XScrSL
     */
    
    auto & HOps = this->hamiltonianOptions_;

    bool allocateLLMS = HOps.BareCoulomb  or HOps.DiracCoulomb; 
    bool allocateSS   = HOps.DiracCoulomb or HOps.DiracCoulombSSSS;
    bool allocateLSSL = HOps.Gaunt or HOps.Gauge;
    
    // TODO: Implement exchange part accordingly 
    
    for (auto i = 0ul; i < mPDM; i++) {
      // Allocate Scattered Density
      if (allocateLLMS) ALLOCATE_PAULISPINOR_SCR(contractSymm1PDMLLMS, NB1C, false); 
      if (allocateSS)   ALLOCATE_PAULISPINOR_SCR(contractSymm1PDMSS, NB1C, true);    
      if (allocateLSSL) ALLOCATE_PAULISPINOR_SCR(contract1PDMLSpmSL, NB1C, true);    
      // allocate Coulomb SCR
      if (allocateLLMS) ALLOCATE_PAULISPINOR_SCR(CScrLLMS, NB1C, false); 
      if (allocateSS)   ALLOCATE_PAULISPINOR_SCR(CScrSS, NB1C, true);    
      if (allocateLSSL) ALLOCATE_PAULISPINOR_SCR(CScrLS, NB1C, true);    
      // allocate Exchange SCR
      // if (computeExchange) {
      //   if (allocateXScrLL)  ALLOCATE_PAULISPINOR_SCR(XScrLL, NB1C, true); 
      //   if (allocateXScrSS)  ALLOCATE_PAULISPINOR_SCR(XScrSS, NB1C, true);    
      //   if (allocateCXScrLS) ALLOCATE_PAULISPINOR_SCR(XScrLS, NB1C, true);    
      //   if (allocateCXScrSL) ALLOCATE_PAULISPINOR_SCR(XScrSL, NB1C, true);    
      // }
    }

    // allocate dummies
    auto dummy_pauli = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(0, false, false);
    
    // Compute 1/(2mc)^2
    MatsT C2 = 1./(4*SpeedOfLight*SpeedOfLight);
    
    // make SCRs
    auto onePDMLLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB1C, false, false);
    auto onePDMSSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB1C, true, true);
    auto onePDMLSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB1C, true, true);
    auto onePDMSLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB1C, true, true);
    
    // Component Scatter Density
    for (auto i = 0ul; i < mPDM; i++) {
      
      auto onePDMLL = allocateLLMS ? onePDMLLSCR: dummy_pauli;
      auto onePDMSS = allocateSS   ? onePDMSSSCR: dummy_pauli;
      auto onePDMLS = allocateLSSL ? onePDMLSSCR: dummy_pauli;
      auto onePDMSL = allocateLSSL ? onePDMSLSCR: dummy_pauli;
      
      onePDMs[i]->componentScatter(*onePDMLL, *onePDMLS, *onePDMSL, *onePDMSS);
      
      // Take full advantage of Integral symmetry for Coulomb-type of terms
      
      if (allocateLLMS) contractSymm1PDMLLMS[i]->S() = onePDMLL->S() + onePDMLL->S().T(); 
      
      if (allocateSS) {
        contractSymm1PDMSS[i]->S() = onePDMSS->S() + onePDMSS->S().T();        
        contractSymm1PDMSS[i]->X() = onePDMSS->X() - onePDMSS->X().T();        
        contractSymm1PDMSS[i]->Y() = onePDMSS->Y() - onePDMSS->Y().T();        
        contractSymm1PDMSS[i]->Z() = onePDMSS->Z() - onePDMSS->Z().T();        
      }
      
      if (allocateLSSL) {
        contract1PDMLSpmSL[i]->S() = onePDMLS->S() - onePDMSL->S().T();        
        contract1PDMLSpmSL[i]->X() = onePDMLS->X() + onePDMSL->X().T();        
        contract1PDMLSpmSL[i]->Y() = onePDMLS->Y() + onePDMSL->Y().T();        
        contract1PDMLSpmSL[i]->Z() = onePDMLS->Z() + onePDMSL->Z().T();        
      }
    } 
    
#ifdef _PRINT_MATRICES
    prettyPrintSmart(std::cout, "1PDM[MS]", contract1PDM[0].S().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MX]", contract1PDM[0].X().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MY]", contract1PDM[0].Y().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MZ]", contract1PDM[0].Z().pointer(), NB2C, NB2C, NB2C);
#endif

    // Initialization 
    if(not increment) {
      for (auto i = 0ul; i < mPDM; i++) {
        if(computeCoulomb)  coulombMatrices[i]->clear(); 
        if(computeExchange) exchangeMatrices[i]->clear(); 
        if(computeTwoeHs)   twoeHs[i]->clear();
      }
    };


    /**********************************************/
    /*                                            */
    /*              DIRECT COULOMB     	          */
    /*                                            */
    /**********************************************/


    if(this->hamiltonianOptions_.BareCoulomb) { // DIRECT_COULOMB

      /*+++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Direct Coulomb (LL|LL) Contraction */
      /*+++++++++++++++++++++++++++++++++++++++++++++*/

      std::vector<TwoBodyRelContraction<MatsT>> contractLL;

      for (auto i = 0ul; i < mPDM; i++) {
        contractLL.push_back({contractSymm1PDMLLMS[i], CScrLLMS[i], HerDen, BARE_COULOMB});
        // if (computeExchange) contractLL.push_back({contract1PDMLL[i], XScrLL[i]});
      } 
      
      // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
      relERICon.twoBodyRelContract(ss.comm, true, contractLL, pert, computeExchange); 
      
      if (computeCoulomb or computeTwoeHs) { 
        for (auto i = 0ul; i < mPDM; i++) {
          coulombContainers[i]->componentAdd('N', MatsT(1.), "LL", *CScrLLMS[i]);
        }  
      }

      // if (computeExchange) {
      //   for (auto i = 0ul; i < mPDM; i++) {
      //     exchangeMatrices[i]->componentAdd('N', MatsT(1.), "LL", *XScrLL[i]);
      //   }  
      // }

#ifdef _PRINT_MATRICES
      std::cout<<"After BARE COULOMB"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S",           twoeHs[0]->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X",           twoeHs[0]->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y",           twoeHs[0]->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z",           twoeHs[0]->Z().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", exchangeMatrices[0]->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", exchangeMatrices[0]->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", exchangeMatrices[0]->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", exchangeMatrices[0]->Z().pointer(), NB2C, NB2C, NB2C);
#endif


      /*---------------------------------------------*/
      /*   End of Direct Coulomb (LL|LL) Contraction */
      /*---------------------------------------------*/

    } // DIRECT_COULOMB

    /**********************************************/
    /*                                            */
    /*              DIRAC-COULOMB                 */
    /*                                            */
    /**********************************************/

    if(this->hamiltonianOptions_.DiracCoulomb) { // DIRAC_COULOMB

  
      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (LL|LL) Contraction */
      /*++++++++++++++++++++++++++++++++++++++++++++*/
  
      if (computeCoulomb or computeTwoeHs) { 
        std::vector<TwoBodyRelContraction<MatsT>> contractDCLL;
        
        for (auto i = 0ul; i < mPDM; i++) {
          contractDCLL.push_back({contractSymm1PDMSS[i], CScrLLMS[i], HerDen, LLLL});
          contractDCLL.push_back({contractSymm1PDMLLMS[i], CScrSS[i]});
        } 

        // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
        relERICon.twoBodyRelContract(ss.comm, true, contractDCLL, pert, computeExchange,
           this->hamiltonianOptions_.DiracCoulombType,
           this->hamiltonianOptions_.DiracCoulombApproximationType);

        // Add Dirac-Coulomb contributions to the LLLL block
        for (auto i = 0ul; i < mPDM; i++) {
          coulombContainers[i]->componentAdd('N', C2, "LL", *CScrLLMS[i]);
          coulombContainers[i]->componentAdd('N', C2, "SS", *CScrSS[i]);
        }  
      } 

      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* End of Dirac-Coulomb (LL|LL) Contraction   */
      /*++++++++++++++++++++++++++++++++++++++++++++*/

#ifdef _PRINT_MATRICES

      std::cout<<"After LLLL"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S",           twoeHs[0]->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X",           twoeHs[0]->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y",           twoeHs[0]->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z",           twoeHs[0]->Z().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", exchangeMatrice[0]->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", exchangeMatrice[0]->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", exchangeMatrice[0]->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", exchangeMatrice[0]->Z().pointer(), NB2C, NB2C, NB2C);

#endif

      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (LS|SL) Contraction */
      /*++++++++++++++++++++++++++++++++++++++++++++*/

      if( computeExchange ) {
#if 0 
      std::vector<TwoBodyContraction<MatsT>> contractDCLS =
        { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, LLSS},
          {contract1PDMLL.S().pointer(), XScrLLMS},
          {contract1PDMLL.X().pointer(), XScrLLMX},
          {contract1PDMLL.Y().pointer(), XScrLLMY},
          {contract1PDMLL.Z().pointer(), XScrLLMZ},
          {contract1PDMSS.S().pointer(), CScrSSMS},
          {contract1PDMSS.X().pointer(), CScrSSMX},
          {contract1PDMSS.Y().pointer(), CScrSSMY},
          {contract1PDMSS.Z().pointer(), CScrSSMZ},
          {contract1PDMSS.S().pointer(), XScrSSMS},
          {contract1PDMSS.X().pointer(), XScrSSMX},
          {contract1PDMSS.Y().pointer(), XScrSSMY},
          {contract1PDMSS.Z().pointer(), XScrSSMZ},
          {contract1PDMLS.S().pointer(), XScrLSMS},
          {contract1PDMLS.X().pointer(), XScrLSMX},
          {contract1PDMLS.Y().pointer(), XScrLSMY},
          {contract1PDMLS.Z().pointer(), XScrLSMZ} };

      // Call the contraction engine to do the assembly of Dirac-Coulomb LLSS
      relERICon.twoBodyContract(ss.comm, true, contractDCLS, pert);

      // Add Dirac-Coulomb contributions to the LLSS block
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMS, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->S().pointer()+LS, NB2C,
		      ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMX, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->X().pointer()+LS, NB2C,
		      ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMY, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->Y().pointer()+LS, NB2C,
		      ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMZ, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->Z().pointer()+LS, NB2C,
		      ss.exchangeMatrix->Z().pointer()+LS, NB2C);
#endif


#ifdef _PRINT_MATRICES

      std::cout<<"After LLSS"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S",           ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X",           ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y",           ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z",           ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
      } 
    
    } //_DIRAC_COULOMB



    /*************************************/
    /*                                   */
    /*              SSSS                 */
    /*                                   */
    /*************************************/

    if(this->hamiltonianOptions_.DiracCoulombSSSS) { // SSSS

      MatsT C4 = 1./(16*SpeedOfLight*SpeedOfLight*SpeedOfLight*SpeedOfLight);
  
      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (SS|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++++++++++*/
  
      std::vector<TwoBodyRelContraction<MatsT>> contractDCSS;
      for (auto i = 0ul; i < mPDM; i++) {
          contractDCSS.push_back({contractSymm1PDMSS[i], CScrSS[i], HerDen, SSSS});
          // if (computeExchange) contractDCSS.push_back({contract1PDMSS[i], XScrSS[i]}); 
      } 

      // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
      relERICon.twoBodyRelContract(ss.comm, true, contractDCSS, pert, computeExchange,
           this->hamiltonianOptions_.SSSSType,
           this->hamiltonianOptions_.SSSSApproximationType);

      // Add (SS|SS) Coulomb contributions to the SSSS block
      for (auto i = 0ul; i < mPDM; i++) {
        coulombContainers[i]->componentAdd('N', C4, "SS", *CScrSS[i]);
        // if (computeExchange) exchangeMatrices[i]->componentAdd('N', -C4, "SS", *XScrSS[i]);
      }  
      
#ifdef _PRINT_MATRICES
      std::cout<<"After SSSS"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S",           ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X",           ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y",           ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z",           ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif
    }


    /*************************************/
    /*                                   */
    /*       GAUNT  and   GAUGE          */
    /*                                   */
    /*************************************/

    // if the gauge term is included, the Gaunt term needs to be scaled by half
    if(this->hamiltonianOptions_.Gauge) C2=C2/2.0;
    
    if (this->hamiltonianOptions_.Gaunt or this->hamiltonianOptions_.Gauge) {
      
      std::vector<TWOBODY_CONTRACTION_TYPE> contractTypes;
      if (this->hamiltonianOptions_.Gaunt) contractTypes.push_back(GAUNT);
      if (this->hamiltonianOptions_.Gauge) contractTypes.push_back(GAUGE);
      
      for (const auto & contT: contractTypes) {  
        
        std::vector<TwoBodyRelContraction<MatsT>> contractDCGau;

        for (auto i = 0ul; i < mPDM; i++) {
          contractDCGau.push_back({contract1PDMLSpmSL[i], CScrLS[i], HerDen, contT});
          // if (computeExchange) {
          //   contractDCGau.push_back({contract1PDMLS[i], XScrLS[i]});
          //   auto XScrSL_i = not HerDen ? XScrSL[i]: nullptr;
          //   contractDCGau.push_back({contract1PDMSL[i], XScrSL_i});
          //   contractDCGau.push_back({contract1PDMLL[i], XScrLL[i]});
          //   contractDCGau.push_back({contract1PDMSS[i], XScrSS[i]});
          // }
        } 

        // Call the contraction engine to do the assembly of Gaunt/Gauge
        TYPE_4C fourCType;
        APPROXIMATION_TYPE_4C approximate4C;
        if (contT == GAUNT) {
          fourCType = this->hamiltonianOptions_.GauntType;
          approximate4C = this->hamiltonianOptions_.GauntApproximationType;
        } else if(contT == GAUGE) {
          fourCType = this->hamiltonianOptions_.GaugeType;
          approximate4C = this->hamiltonianOptions_.GaugeApproximationType;
        }
        
        relERICon.twoBodyRelContract(ss.comm, true, contractDCGau, pert, computeExchange, fourCType, approximate4C);
      
        for (auto i = 0ul; i < mPDM; i++) {
          // Add (LL|SS)  and (SS|LL) Coulomb contributions
          coulombContainers[i]->componentAdd('N', C2, "LS", *CScrLS[i]);
          
          // get the SL contribution using symmetry
          // CSLMS =  - [CLSMS]^T
          // CSLMX =    [CLSMX]^T
          // CSLMY =    [CLSMY]^T
          // CSLMZ =    [CLSMZ]^T
          CScrLS[i]->S() = - CScrLS[i]->S();
          coulombContainers[i]->componentAdd('T', C2, "SL", *CScrLS[i]);

          // if (computeExchange) {
          //   // Add (LL|SS)  and (SS|LL) Exchange contributions
          //   exchangeMatrices[i]->componentAdd('N', -C2, "LS", *XScrLS[i]);
          //   if (not HerDen) exchangeMatrices[i]->componentAdd('N', -C2, "SL", *XScrSL[i]);  
          //   // Add (LL|LL) exchange contributions
          //   exchangeMatrices[i]->componentAdd('N', -C2, "LL", *XScrLL[i]);
          //   // Add (SS|SS) exchange contributions
          //   exchangeMatrices[i]->componentAdd('N', -C2, "SS", *XScrSS[i]);
          // }
        }  

#ifdef _PRINT_MATRICES
        if (contT == GAUNT) std::cout<<"After GAUNT"<<std::endl;
        if (contT == GAUGE) std::cout<<"After GAUGE"<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB-S",           ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "COULOMB-X",           ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "COULOMB-Y",           ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "COULOMB-Z",           ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif
      
      } // contractType
    
    } // Gaunt and Gauge

    /*******************************/
    /* Final Assembly of 4C Matrix */
    /*******************************/
    ROOT_ONLY(ss.comm);
    
    /*************************************************/
    /* Hermitrize Coulomb and Exchange Part if HerDen*/
    /*************************************************/
    
    if (HerDen) {
      if (computeCoulomb or computeTwoeHs) { 
        for (auto i = 0ul; i < mPDM; i++) {
          coulombContainers[i]->symmetrizeLSSL('C'); 
        }  
      }

      // if (computeExchange) { 
      //   for (auto i = 0ul; i < mPDM; i++) {
      //     exchangeMatrices[i]->symmetrizeLSSL('C');
      //   }  
      // }
    }

    /*************************************************/
    /* Sum Coulomb and Exchange to twoeH             */
    /*************************************************/
    
    if (computeTwoeHs) {
      for (auto i = 0ul; i < mPDM; i++) { 
        // G[D] += 2*J[D]
        if (computeCoulomb) {
          *twoeHs[i] += 2.0 * *coulombMatrices[i];
        } else {
          *twoeHs[i] *= 2.0;
        }

        // // Form GD: G[D] = 2.0*J[D] - K[D]
        // if (computeExchange) {
        //   *twoeHs[i] -= xHFX * *exchangeMatrices[i];
        // } 
      }
    }

  }; // FourCompFock<MatsT, IntsT>::formRawGDInBatchesDirect

  /*******************************************************************************/
  /*                                                                             */
  /* Compute memory requirement for build 4C GD in Batches                       */
  /* Returns:                                                                    */
  /*   size_t SCR size needed for one batch                                      */
  /*   IMPORTANT HERE: size are all in MatsT (dcomplex)                          */
  /*******************************************************************************/
  template <typename MatsT, typename IntsT>
  size_t FourCompFock<MatsT,IntsT>::formRawGDSCRSizePerBatch(SingleSlater<MatsT,IntsT> &ss,
    bool computeExchange, bool HerDen) const {
      
      if (computeExchange) CErr("computeExchange NYI in FockBuilder::formRawGDSCRSizePerBatch"); 

      size_t SCRSize  = 0ul;
       
      if( std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI) ) {
        
        GTODirectRelERIContraction<MatsT,IntsT> &relERICon =
            *std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI);
        
        // Update with contraction SCR 
        #define UPDATE_CONTRACTION_SCR_SIZE(CONTTYPE) \
          auto contSCR = relERICon.directRelScaffoldLibcintSCRSize(CONTTYPE, computeExchange); \
          SCRSize  = std::max(SCRSize, contSCR); 

        if (this->hamiltonianOptions_.BareCoulomb) {  
          UPDATE_CONTRACTION_SCR_SIZE(BARE_COULOMB);
        }
        if (this->hamiltonianOptions_.DiracCoulomb) { 
          UPDATE_CONTRACTION_SCR_SIZE(LLLL);
          if (computeExchange) UPDATE_CONTRACTION_SCR_SIZE(LLSS);
        }
        if (this->hamiltonianOptions_.DiracCoulombSSSS) {
          UPDATE_CONTRACTION_SCR_SIZE(SSSS);
        }
        if (this->hamiltonianOptions_.Gaunt) {
          UPDATE_CONTRACTION_SCR_SIZE(GAUNT);
        }
        if (this->hamiltonianOptions_.Gauge) {
          UPDATE_CONTRACTION_SCR_SIZE(GAUGE);
        }
        
        // plus extra SCR to build component scattered X and AX
        auto & HOps = this->hamiltonianOptions_;

        // see line 197 for sepecific allocations
        bool allocateLLMS = HOps.BareCoulomb  or HOps.DiracCoulomb; 
        bool allocateSS   = HOps.DiracCoulomb or HOps.DiracCoulombSSSS;
        bool allocateLSSL = HOps.Gaunt or HOps.Gauge;
        
        size_t NB1C  = ss.basisSet().nBasis;
        size_t NB1C2 = NB1C*NB1C;
        // density + Coulomb SCR requirements
        SCRSize += allocateLLMS ? NB1C2*2: 0;
        SCRSize += allocateSS   ? NB1C2*8: 0;
        SCRSize += allocateLSSL ? NB1C2*8: 0;
        
        //// eXchange SCR
        //if (computeExchange) {
        //  bool allocateXScrLL = this->hamiltonianOptions_.BareCoulomb or 
        //    this->hamiltonianOptions_.Gaunt or 
        //    this->hamiltonianOptions_.Gauge;
        //  bool allocateXScrSS = this->hamiltonianOptions_.DiracCoulombSSSS or 
        //    this->hamiltonianOptions_.Gaunt or 
        //    this->hamiltonianOptions_.Gauge;
        //  
        //  SCRSize += allocateXScrLL ? NB1C2*4: 0;
        //  SCRSize += allocateXScrSS ? NB1C2*4: 0;
        //  SCRSize += allocateCXScrLS ? NB1C2*4: 0;
        //  SCRSize += allocateCXScrSL ? NB1C2*4: 0;
        //}

      }
      
      return SCRSize;

  }; // FourCompFock<MatsT, IntsT>::formRawGDSCRSizePerBatch



}; // namespace ChronusQ
