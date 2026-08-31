/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
//#define _FOCK_CONTRACTION_TIME
//#define _SS_J_TIME
//#define _SS_K_TIME

namespace ChronusQ {

  /**   
   *  \brief Forms the 4C Fock matrix
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formGD(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {

    if( std::dynamic_pointer_cast<InCoreRelERIContraction<MatsT,IntsT>>(ss.TPI) )
      formGDInCore(ss, pert, increment, xHFX, HerDen);
    else if( std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI) )
      formGDDirect(ss, pert, increment, xHFX, HerDen);
      //formGD3Index(ss, pert, increment, xHFX, HerDen);
    else
      CErr("Unsupported ERIContraction type.");

  };


  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formGDInCore(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {

    InCoreRelERI<IntsT> &relERI =
        *std::dynamic_pointer_cast<InCoreRelERI<IntsT>>(ss.aoints_->TPI);

    bool computeExchange = true;
    bool do_VLL = false;
    double xHFX_VLL = 1.; 
    
    // Check if doing DFT
    // DKS w/ VXCLL requires computeExchange for all functionals (even pure)
    if (typeid(ss) == typeid(KohnSham<MatsT,IntsT>)){
      computeExchange = true;
    if(this->hamiltonianOptions_.dksType == DKS_TYPE::VLL){
      computeExchange = true;
      xHFX_VLL = xHFX;
      do_VLL = true;
    }}
    else {
	    computeExchange = std::abs(xHFX) >= 1e-12;
    }

    
    if (not HerDen and computeExchange) CErr("formGDInCore with exchange term NYI for non-Hermitian density ");
    
    bool RESET = false;
    
    // Decide list of onePDMs to use
    cqmatrix::PauliSpinorMatrices<MatsT> &contract1PDM
        = increment ? *ss.deltaOnePDM : *ss.onePDM;

    size_t NB1C  = ss.basisSet().nBasis;
    size_t NB2C  = 2 * ss.basisSet().nBasis;
    size_t NB4C  = 4 * ss.basisSet().nBasis;
    size_t NB1C2 = NB1C*NB1C;
    size_t NB1C4 = NB1C*NB1C*NB1C*NB1C;
    size_t NB2C2 = NB2C*NB2C;
    size_t NB4C2 = NB4C*NB4C;

    size_t SS = NB2C*NB1C+NB1C;
    size_t LS = NB2C*NB1C;
    size_t SL = NB1C;

    size_t mpiRank   = MPIRank(ss.comm);
    bool   isNotRoot = mpiRank != 0;

    cqmatrix::PauliSpinorMatrices<MatsT> exchangeMatrixLL(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLL(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSL(NB1C);

    MatsT* Scr1 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr2 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr3 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr4 = CQMemManager::get().malloc<MatsT>(NB1C2);
    memset(Scr1,0.,NB1C2*sizeof(MatsT));
    memset(Scr2,0.,NB1C2*sizeof(MatsT));
    memset(Scr3,0.,NB1C2*sizeof(MatsT));
    memset(Scr4,0.,NB1C2*sizeof(MatsT));


    // Compute 1/(2mc)^2
    //dcomplex scale = 1.;
    //dcomplex iscale = dcomplex(0.0, 1.0);
    dcomplex scale = 1./(4*SpeedOfLight()*SpeedOfLight());
    dcomplex iscale = dcomplex(0.0, 1./(4*SpeedOfLight()*SpeedOfLight()));



#if 0
    MatsT* DEN_GATHER = CQMemManager::get().malloc<MatsT>(NB4C2);
    MatsT* LS_GATHER = CQMemManager::get().malloc<MatsT>(NB2C2);

    for(auto i = 0; i< NB2C2; i++) {
      LS_GATHER[i] = dcomplex(std::cos(i), std::sin(i));
    }
    SetMat('N', NB2C, NB2C, MatsT(1.), LS_GATHER, NB2C, DEN_GATHER+NB4C*NB2C, NB4C);
    SetMat('C', NB2C, NB2C, MatsT(1.), LS_GATHER, NB2C, DEN_GATHER+NB2C, NB4C);
    SetMat('N', NB2C, NB2C, MatsT(1.), LS_GATHER, NB2C, DEN_GATHER, NB4C);
    SetMat('N', NB2C, NB2C, MatsT(1.), LS_GATHER, NB2C, DEN_GATHER+NB4C*NB2C+NB2C, NB4C);
    MatAdd('C','N', NB2C, NB2C, MatsT(1.0), LS_GATHER, 
		    NB2C, MatsT(1.0), DEN_GATHER, NB4C, DEN_GATHER, NB4C);
    MatAdd('C','N', NB2C, NB2C, MatsT(1.0), LS_GATHER, 
		    NB2C, MatsT(1.0), DEN_GATHER+NB4C*NB2C+NB2C, 
		    NB4C, DEN_GATHER+NB4C*NB2C+NB2C, NB4C);

    prettyPrintSmart(std::cout,"Initial 4C Density",DEN_GATHER,NB4C,NB4C,NB4C);

    std::fill_n(DEN_GATHER,NB4C2,1.0);
    SpinScatter(NB2C,DEN_GATHER, NB4C, contract1PDM.S().pointer(),
		    NB2C,contract1PDM.Z().pointer(), NB2C,contract1PDM.Y().pointer(),
		    NB2C,contract1PDM.X().pointer(), NB2C);

    CQMemManager::get().free(DEN_GATHER);
    CQMemManager::get().free(LS_GATHER);
#endif

    for(size_t i = 0; i < contract1PDM.nComponent(); i++) {
      cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer(),    NB2C,
             contract1PDMLL[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SS, NB2C,
             contract1PDMSS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+LS, NB2C,
             contract1PDMLS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SL, NB2C,
             contract1PDMSL[c].pointer(), NB1C);
    }
    
#ifdef _PRINT_MATRICES
    prettyPrintSmart(std::cout, "1PDM[MS]", contract1PDM.S().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MX]", contract1PDM.X().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MY]", contract1PDM.Y().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MZ]", contract1PDM.Z().pointer(), NB2C, NB2C, NB2C);
#endif

#if 0
    std::fill_n(contract1PDMLL.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.Z().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.Z().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.Z().pointer(),NB1C2,1.0);
#endif

    // Initializtion of resulting matries
    if(not increment) {
      ss.twoeH->clear();
      ss.coulombMatrix->clear();
      ss.exchangeMatrix->clear();
    };
    


    /**********************************************/
    /*                                            */
    /*   NON-RELATIVISTIC DIRECT COULOMB          */
    /*                                            */
    /**********************************************/

    if(this->hamiltonianOptions_.BareCoulomb) { // DIRECT_COULOMB

#ifdef _FOCK_CONTRACTION_TIME
      auto topBareCoulomb = tick();
#endif

      /*+++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Direct Coulomb (LL|LL) Contraction */
      /*+++++++++++++++++++++++++++++++++++++++++++++*/
  
#ifdef _FOCK_CONTRACTION_TIME
      auto topBC_J = tick();
#endif

      std::vector<TwoBodyContraction<MatsT>> contractLL =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, COULOMB} };

      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractLL, pert);


      /* Store LL block into 2C spin scattered matrices */
      // Assemble 4C coulombMatrix
      SetMat('N', NB1C, NB1C, MatsT(1.), Scr1, NB1C, ss.coulombMatrix->pointer(), NB2C);

#ifdef _FOCK_CONTRACTION_TIME
      // Print out BareCoulomb duration 
      auto durBC_J = tock(topBC_J);
      std::cout << "Bare Coulomb J Contraction to F_LL = " << durBC_J << std::endl;
#endif 




#ifdef _FOCK_CONTRACTION_TIME
      auto topBC_K = tick();
#endif

      // Determine how many (if any) exchange terms to calculate
      std::vector<TwoBodyContraction<MatsT>> contractLLK; 
      if (computeExchange)
      for(size_t i = 0; i < ss.exchangeMatrix->nComponent(); i++) {
  
        cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
        contractLLK.push_back(
          {contract1PDMLL[c].pointer(), exchangeMatrixLL[c].pointer(), HerDen, EXCHANGE}
        );
      }

      // Zero out K[i]
      if(not increment) ss.exchangeMatrix->clear();

      ss.TPI->twoBodyContract(ss.comm, contractLLK, pert);

      // Assemble 4C exchangeMatrix 
      if(computeExchange) {
      for(auto i = 0; i < ss.exchangeMatrix->nComponent();i++){
        cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
        SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL, exchangeMatrixLL[c].pointer(), NB1C,
               (*ss.exchangeMatrix)[c].pointer(), NB2C);
      }
      }

#ifdef _FOCK_CONTRACTION_TIME
      // Print out BareCoulomb duration 
      auto durBC_K = tock(topBC_K);
      std::cout << "Bare Coulomb K Contraction to F_LL  = " << durBC_K << std::endl;
#endif 

#ifdef _PRINT_MATRICES

      std::cout<<"After BARE COULOMB"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);

#endif //_PRINT_MATRICES


      /*---------------------------------------------*/
      /*   End of Direct Coulomb (LL|LL) Contraction */
      /*---------------------------------------------*/

#ifdef _FOCK_CONTRACTION_TIME
      // Print out BareCoulomb duration 
      auto durBareCoulomb = tock(topBareCoulomb);
      std::cout << "Bare Coulomb Contraction duration = " << durBareCoulomb << std::endl;
#endif 

    } // DIRECT_COULOMB


    /**********************************************/
    /*                                            */
    /*              DIRAC-COULOMB      	          */
    /*                                            */
    /**********************************************/

    // ERI: (ab|cd)
    // ERIDCB0: ∇A∙∇B(ab|cd)
    // ERIDCB1: ∇Ax∇B(ab|cd)-X
    // ERIDCB2: ∇Ax∇B(ab|cd)-Y
    // ERIDCB3: ∇Ax∇B(ab|cd)-Z

    if(this->hamiltonianOptions_.DiracCoulomb) { // DIRAC_COULOMB

#ifdef _FOCK_CONTRACTION_TIME
      auto topERIDC = tick();
#endif

#ifdef _FOCK_CONTRACTION_TIME
      auto topDC_LLLL = tick();
#endif

      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (LL|LL) Contraction */
      /*++++++++++++++++++++++++++++++++++++++++++++*/
      std::vector<TwoBodyContraction<MatsT>> contractDCLL =  
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, DC_COULOMB, 0, TRANS_MNKL},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, DC_COULOMB, 1, TRANS_MNKL},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, DC_COULOMB, 2, TRANS_MNKL},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, DC_COULOMB, 3, TRANS_MNKL} };

      // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
      ss.TPI->twoBodyContract(ss.comm, contractDCLL);
  
      // Add Dirac-Coulomb contributions  to the LLLL block
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), 
		      ss.coulombMatrix->pointer(), NB2C, ss.coulombMatrix->pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), 
		      ss.coulombMatrix->pointer(), NB2C, ss.coulombMatrix->pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), 
		      ss.coulombMatrix->pointer(), NB2C, ss.coulombMatrix->pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), 
		      ss.coulombMatrix->pointer(), NB2C, ss.coulombMatrix->pointer(), NB2C);
#ifdef _PRINT_MATRICES

      std::cout<<"After LLLL"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);

#endif //_PRINT_MATRICES

      /*------------------------------------------*/
      /* End of Dirac-Coulomb (LL|LL) Contraction */
      /*------------------------------------------*/
  
#ifdef _FOCK_CONTRACTION_TIME
      // Print out DC-LLLL duration 
      auto durDC_LLLL = tock(topDC_LLLL);
      std::cout << "DC Contribution J to F_LL = " << durDC_LLLL << std::endl;
#endif 
  
  
#ifdef _FOCK_CONTRACTION_TIME
      auto topDC_SSSS = tick();
#endif
  
      /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb C(2)-(SS|SS) Contraction */
      /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
  
      std::vector<TwoBodyContraction<MatsT>> contractSS =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, DC_COULOMB, 0},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, DC_COULOMB, 1},
          {contract1PDMLL.S().pointer(), Scr3, HerDen, DC_COULOMB, 2},
          {contract1PDMLL.S().pointer(), Scr4, HerDen, DC_COULOMB, 3} };

      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSS);
  
      // Store SS block into 2C spin scattered matrices 
      // These scaling factors were modified to take into account the issue of storing the 
      // Coulomb portion of C(2)-(SS|SS) are directly put to twoeH 
      MatAdd('N','N', NB1C, NB1C, MatsT(2.0*scale),  Scr1, NB1C, MatsT(1.0), 
                      ss.twoeH->S().pointer()+SS, NB2C,
                      ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, MatsT(2.0*iscale), Scr2, NB1C, MatsT(1.0), 
                      ss.twoeH->X().pointer()+SS, NB2C,
                      ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, MatsT(2.0*iscale), Scr3, NB1C, MatsT(1.0), 
                      ss.twoeH->Y().pointer()+SS, NB2C,
                      ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, MatsT(2.0*iscale), Scr4, NB1C, MatsT(1.0), 
                      ss.twoeH->Z().pointer()+SS, NB2C,
                      ss.twoeH->Z().pointer()+SS, NB2C);
      //SetMat('N', NB1C, NB1C, MatsT(scale),       Scr1, NB1C, ss.coulombMatrix->pointer()+SS,      NB2C);
      //SetMat('N', NB1C, NB1C, MatsT(-2.0*iscale), Scr2, NB1C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      //SetMat('N', NB1C, NB1C, MatsT(-2.0*iscale), Scr3, NB1C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      //SetMat('N', NB1C, NB1C, MatsT(-2.0*iscale), Scr4, NB1C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
#ifdef _PRINT_MATRICES
  
      std::cout<<"After SSSS"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMBSS-S", ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMBSS-X", ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMBSS-Y", ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMBSS-Z", ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif //_PRINT_MATRICES
  
  
      /*-----------------------------------------------*/
      /* End of Dirac-Coulomb C(2)-(SS|SS) Contraction */
      /*-----------------------------------------------*/
#ifdef _FOCK_CONTRACTION_TIME
      // Print out DC-SSSS duration 
      auto durDC_SSSS = tock(topDC_SSSS);
      std::cout << "DC Contribution J to F_SS = " << durDC_SSSS << std::endl;
#endif 
  
  
#if 1 
      if (computeExchange) {
      /*++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (LL|SS) / (SS|LL) */
      /*++++++++++++++++++++++++++++++++++++++++++*/
#ifdef _FOCK_CONTRACTION_TIME
      auto topDC_SSLL_Ex = tick();
#endif
  
      std::vector<TwoBodyContraction<MatsT>> contractLSScalar =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, DC_EXCHANGE, 0, TRANS_MNKL},
          {contract1PDMLS.X().pointer(), Scr2, HerDen, DC_EXCHANGE, 1, TRANS_MNKL},
          {contract1PDMLS.Y().pointer(), Scr3, HerDen, DC_EXCHANGE, 2, TRANS_MNKL},
          {contract1PDMLS.Z().pointer(), Scr4, HerDen, DC_EXCHANGE, 3, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractLSScalar);
  
      // Add to the LS part of the exchangeMatrix[MS]
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0),
             ss.exchangeMatrix->S().pointer()+LS, NB2C,
             ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0),
             ss.exchangeMatrix->S().pointer()+LS, NB2C,
             ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0),
             ss.exchangeMatrix->S().pointer()+LS, NB2C,
             ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0),
             ss.exchangeMatrix->S().pointer()+LS, NB2C,
             ss.exchangeMatrix->S().pointer()+LS, NB2C);
  
  
  
      std::vector<TwoBodyContraction<MatsT>> contractLSMX =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, DC_EXCHANGE, 0, TRANS_MNKL},
          {contract1PDMLS.S().pointer(), Scr2, HerDen, DC_EXCHANGE, 1, TRANS_MNKL},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, DC_EXCHANGE, 2, TRANS_MNKL},
          {contract1PDMLS.Y().pointer(), Scr4, HerDen, DC_EXCHANGE, 3, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractLSMX);
  
      // Add to the LS part of the exchangeMatrix[MX]
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C,
             MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C,
             ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C,
             MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C,
             ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C,
             MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C,
             ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C,
             MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C,
             ss.exchangeMatrix->X().pointer()+LS, NB2C);
  
  
  
  
  
      std::vector<TwoBodyContraction<MatsT>> contractLSMY =
        { {contract1PDMLS.Y().pointer(), Scr1, HerDen, DC_EXCHANGE, 0, TRANS_MNKL},
          {contract1PDMLS.Z().pointer(), Scr2, HerDen, DC_EXCHANGE, 1, TRANS_MNKL},
          {contract1PDMLS.S().pointer(), Scr3, HerDen, DC_EXCHANGE, 2, TRANS_MNKL},
          {contract1PDMLS.X().pointer(), Scr4, HerDen, DC_EXCHANGE, 3, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractLSMY);
  
      // Add to the LS part of the exchangeMatrix[MY]
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C,
             ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr2, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C,
             ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C,
             ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C,
             ss.exchangeMatrix->Y().pointer()+LS, NB2C);
  
  
  
  
  
      std::vector<TwoBodyContraction<MatsT>> contractLSMZ =
        { {contract1PDMLS.Z().pointer(), Scr1, HerDen, DC_EXCHANGE, 0, TRANS_MNKL},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, DC_EXCHANGE, 1, TRANS_MNKL},
          {contract1PDMLS.X().pointer(), Scr3, HerDen, DC_EXCHANGE, 2, TRANS_MNKL},
          {contract1PDMLS.S().pointer(), Scr4, HerDen, DC_EXCHANGE, 3, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractLSMZ);
  
      // Add to the LS part of the exchangeMatrix[MZ]
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C,
             ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr2, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C,
             ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C,
             ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C,
             MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C,
             ss.exchangeMatrix->Z().pointer()+LS, NB2C);
#endif  
  
  
#ifdef _PRINT_MATRICES
  
      std::cout<<"After Dirac-Coulomb"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif //_PRINT_MATRICES
  
      /*------------------------------------------*/
      /*   End of Dirac-Coulomb (LL|SS) / (SS|LL) */
      /*------------------------------------------*/

#ifdef _FOCK_CONTRACTION_TIME
      // Print out DC-SSSS duration 
      auto durDC_SSLL_Ex = tock(topDC_SSLL_Ex);
      std::cout << "DC Contribution K to F_SL = " << durDC_SSLL_Ex << std::endl;
#endif 
    
    }  // computeExchange
      
#ifdef _FOCK_CONTRACTION_TIME 
      auto durERIDC = tock(topERIDC);
      std::cout << "Dirac-Coulomb Contraction duration   = " << durERIDC << std::endl;
#endif    
    } // Dirac Coulomb
  





    /**********************************************/
    /*                                            */
    /*              GAUNT                         */
    /*                                            */
    /**********************************************/


    //ERI4:    ∇B∙∇C(mn|kl)
    //ERI5 :   ∇Bx∇C(mn|kl)-X
    //ERI6 :   ∇Bx∇C(mn|kl)-Y
    //ERI7 :   ∇Bx∇C(mn|kl)-Z
    //ERI8 :   ∇B_x∇C_y(mn|kl) + ∇B_y∇C_x(mn|kl)
    //ERI9 :   ∇B_y∇C_x(mn|kl)
    //ERI10:   ∇B_x∇C_z(mn|kl) + ∇B_z∇C_x(mn|kl)
    //ERI11:   ∇B_z∇C_x(mn|kl)
    //ERI12:   ∇B_y∇C_z(mn|kl) + ∇B_z∇C_y(mn|kl)
    //ERI13:   ∇B_z∇C_y(mn|kl)
    //ERI14: - ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) + ∇B_z∇C_z(mn|kl)
    //ERI15:   ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
    //ERI16: - ∇B_x∇C_x(mn|kl) + ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
    //ERI17:   ∇B_x∇C_x(mn|kl)
    //ERI18:   ∇B_x∇C_y(mn|kl)
    //ERI19:   ∇B_x∇C_z(mn|kl)
    //ERI20:   ∇B_y∇C_y(mn|kl)
    //ERI21:   ∇B_y∇C_z(mn|kl)
    //ERI22:   ∇B_z∇C_z(mn|kl)


    if(this->hamiltonianOptions_.Gaunt) {//GAUNT

// SS start
      // when using Breit interaction, all the gaunt and gauge multiplied by 1/2 
      if(this->hamiltonianOptions_.Gauge) {//Gauge
        // std::cout<<"use Breit interaction, Gaunt part devide by 2"<<std::endl;
        scale = 0.5* 1./(4*SpeedOfLight()*SpeedOfLight());
        iscale = 0.5* dcomplex(0.0, 1./(4*SpeedOfLight()*SpeedOfLight()));
      }
// SS end


#ifdef _FOCK_CONTRACTION_TIME
      auto topERIDG = tick();
#endif

#if 0 // Gaunt LLLL Spin-Free
      /* Gaunt LLLL Spin-Free */
      std::vector<TwoBodyContraction<MatsT>> contractGLLSF94 =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 0},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 0},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 0} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLLSF94);
  
      // Add to the LL part of 4C exchangeMatrix in Pauli matrix form
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
#endif // Gaunt LLLL Spin-Orbit
  
#if 0 // Gaunt LLLL Spin-Orbit
      /* Gaunt LLLL Spin-Orbit */
      /* Equation (103) */
      std::vector<TwoBodyContraction<MatsT>> contractGLLSO103 =
        { {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLLSO103);
  
      // Add to the LL part of 4C exchangeMatrix in Pauli matrix form
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
  
  
      /* Equation (104)-(106) */
      std::vector<TwoBodyContraction<MatsT>> contractGLLSO104106 =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 1},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 2},
          {contract1PDMSS.S().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLLSO104106);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
  
#endif // Gaunt LLLL Spin-Orbit
  
#if 1 //Gaunt LLLL
      
      if (computeExchange) {
      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gaunt (LL|LL) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/
  
      /* Equation (113) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL113 =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLL113);
  
      // Add to the LL part of 4C exchangeMatrix in Pauli matrix form
      MatAdd('N','N', NB1C, NB1C, -3.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 3.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 3.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, 3.0*iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
  
      /* Equation (114) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL114 =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 10},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 6},
          {contract1PDMSS.Y().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 8} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLL114);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
  
      /* Equation (115) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL115 =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 11},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 4},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 6} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLL115);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
  
  
  
      /* Equation (116) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL116 =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 12},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 2},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 4},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 8} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLL116);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);

#endif //Gaunt LLLL
  
#ifdef _PRINT_MATRICES
  
      std::cout<<"After Gaunt LLLL"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "GAUNT_EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "GAUNT_EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "GAUNT_EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "GAUNT_EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif //_PRINT_MATRICES
   
      if(RESET) {
        ss.coulombMatrix->clear();
        ss.exchangeMatrix->clear();
      }
  
  
  
      /*----------------------------------*/
      /* End of Gaunt (LL|LL) Contraction */
      /*----------------------------------*/
  
  
  
      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gaunt (SS|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/
  
#if 0 // Gaunt SSSS Spin-Free
      /* Gaunt SSSS Spin-Free */
      /* Equation (118) */
      std::vector<TwoBodyContraction<MatsT>> contractGSSSF118 =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, 1},
          {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 0, 1},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 0, 1},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 0, 1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSSSF118);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
#endif // Gaunt SSSS Spin-Free
  
#if 0 // Gaunt SSSS Spin-Orbit
      /* Gaunt SSSS Spin-Orbit */
      /* Equation (119) */
      std::vector<TwoBodyContraction<MatsT>> contractGSSSO119 =
        { {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1, 1},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, 1},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3, 1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSSSO119);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
  
      /* Equation (120)-(122) */
      std::vector<TwoBodyContraction<MatsT>> contractGSSSO120122 =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 1, 1},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 2, 1},
          {contract1PDMLL.S().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 3, 1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSSSO120122);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
  
#endif // Gaunt SSSS Spin-Orbit
  
#if 1 //Gaunt SSSS
  
      /* Equation (129) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS129 =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, TRANS_MNKL},
          {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1, TRANS_MNKL},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, TRANS_MNKL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSS129);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-3.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,    iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,    iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,    iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);

      /* Equation (130) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS130 =
        { {contract1PDMLL.Z().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 10, TRANS_MNKL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3,  TRANS_MNKL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 6, TRANS_MNKL},
          {contract1PDMLL.Y().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 8, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSS130);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,      scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  3.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
  
  
      /* Equation (131) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS131 =
        { {contract1PDMLL.X().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 11, TRANS_MNKL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1,  TRANS_MNKL},
          {contract1PDMLL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 6, TRANS_MNKL},
          {contract1PDMLL.Y().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 4,  TRANS_MNKL}};
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSS131);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,      scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  3.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
  
  
      /* Equation (132) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS132 =
        { {contract1PDMLL.Y().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 12, TRANS_MNKL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 2,  TRANS_MNKL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 4,  TRANS_MNKL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 8, TRANS_MNKL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGSS132);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,      scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  3.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,      scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);

#endif // Gaunt SSSS
  
#ifdef _PRINT_MATRICES
  
      std::cout<<"After Gaunt SSSS"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif //_PRINT_MATRICES
  
      if(RESET) {
        ss.coulombMatrix->clear();
        ss.exchangeMatrix->clear();
      }
      /*------------------------------------*/
      /*   End of Gaunt (SS|SS) Contraction */
      /*------------------------------------*/
      
      } // computeExchange, Gaunt (LL|LL) and (SS|SS) are all exchange terms
  
  
      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gaunt (LL|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/
  
#if 0 // Gaunt LLSS Spin-Free
      /* Gaunt LLSS Spin-Free */
      /* First term in Equations (91) and (136) */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSF91136 =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 0},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSF91136);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
  
      /* Gaunt LLSS Spin-Free */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSF140 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, 2},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 0, 2},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 0, 2},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 0, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSF140);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
#endif // Gaunt LLSS Spin-Free
  
#if 0 // Gaunt LLSS Spin-Orbit
      /* Gaunt LLSS Spin-Orbit */
  
      /* Equation (91) second term */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO91 =
        { {contract1PDMLS.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 1},
          {contract1PDMLS.Y().pointer(), Scr3, HerDen, GAUNT_COULOMB, 2},
          {contract1PDMLS.Z().pointer(), Scr4, HerDen, GAUNT_COULOMB, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO91);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
  
  
      /* Equation (92) first term */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO92 =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 1},
          {contract1PDMLS.S().pointer(), Scr2, HerDen, GAUNT_COULOMB, 2},
          {contract1PDMLS.S().pointer(), Scr3, HerDen, GAUNT_COULOMB, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO92);
  
       // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
  
  
      /* Equation (136) second term*/
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO136 =
        { {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 1, 2},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUNT_COULOMB, 2, 2},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUNT_COULOMB, 3, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO136);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
  
      /* Equation (137) first term */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO137 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 1, 2},
          {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 2, 2},
          {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 3, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO137);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
  
  
      /* Equation (150) */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO150 =
        { {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1, 2},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO150);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
  
      
  
      /* Equation (151) */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO151 =
        { {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3, 2},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO151);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
  
  
  
      /* Equation (152) */
      std::vector<TwoBodyContraction<MatsT>> contractGLSSO152 =
        { {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3, 2},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 1, 2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLSSO152);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
  
  
      std::cout<<"After Gaunt LLSS Spin-Orbit "<<std::endl;
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif // Gaunt LLSS Spin-Orbit
  
  
#if 1 // Gaunt LLSS GAUNT_COULOMB
  
#if 1      
      /* Equation (91) */
      std::vector<TwoBodyContraction<MatsT>> contractGLS91 =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 0},
          {contract1PDMLS.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 1},
          {contract1PDMLS.Y().pointer(), Scr3, HerDen, GAUNT_COULOMB, 2},
          {contract1PDMLS.Z().pointer(), Scr4, HerDen, GAUNT_COULOMB, 3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS91);
  
      // Assemble coulomb contributation direct to twoeH 
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);

      /* Equation (92)Z first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92AZ =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 3},
          {contract1PDMLS.Z().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92AZ);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
  
      /* Equation (92)Z last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92BZ =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 15},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 17},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 18} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92BZ);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);

      /* Equation (92)X first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92AX =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 1},
          {contract1PDMLS.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92AX);
  
      // Assemble coulomb contributation direct to twoeH 
      MatAdd('N','N', NB1C, NB1C,  2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  
  
      /* Equation (92)X last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92BX =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 13},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 5},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 7} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92BX);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  
  
      /* Equation (92)Y first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92AY =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 2},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92AY);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
  
      /* Equation (92)Y last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS92BY =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 14},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 16},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 9} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS92BY);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
#endif
  
#ifdef _PRINT_MATRICES
      std::cout<<"After Gaunt 91-92"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S", ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X", ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y", ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z", ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
#endif
  
      if(RESET) {
        ss.coulombMatrix->clear();
        ss.exchangeMatrix->clear();
      }
  
#if 1      
      /* Equation (136) */
      std::vector<TwoBodyContraction<MatsT>> contractGLS136 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 0, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 1, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUNT_COULOMB, 2, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUNT_COULOMB, 3, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS136);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
  
      /* Equation (137)X first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137AX =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 1, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137AX);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  
  
      /* Equation (137)X last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137BX =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 13, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 5,  TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 7, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137BX);
 
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  
      /* Equation (137)Y first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137AY =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 2, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137AY);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
  
  
      /* Equation (137)Y last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137BY =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 14, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 16, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 9, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137BY);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
  
      /* Equation (137)Z first two terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137AZ =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_COULOMB, 3, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr2, HerDen, GAUNT_COULOMB, 0, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137AZ);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
  
  
      /* Equation (137)Z last term */
      std::vector<TwoBodyContraction<MatsT>> contractGLS137BZ =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUNT_COULOMB, 15, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_COULOMB, 17, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_COULOMB, 18, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS137BZ);
  
      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
  
#endif      
  
#ifdef _PRINT_MATRICES
      std::cout<<"After Gaunt 136-137"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S", ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X", ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y", ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z", ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
#endif
      if(RESET) {
        ss.coulombMatrix->clear();
        ss.exchangeMatrix->clear();
      }
  
#endif  // Gaunt LLSS COULOMB
      
      if (computeExchange) {
#if 1 // Gaunt LLSS EXCHANGE
      /* Equation (159) */
      std::vector<TwoBodyContraction<MatsT>> contractGLS159 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS159);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
  
      /* Equation (160) first four terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS160A =
        { {contract1PDMSL.Z().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 1, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 3, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS160A);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,   -iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
  
  
      /* Equation (160) last three terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS160B =
        { {contract1PDMSL.Z().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 10, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 6, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 8, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS160B);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
  
  
  
      /* Equation (161) first four terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS161A =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 1, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS161A);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,   -iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
  
  
      /* Equation (161) last three terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS161B =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 11, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 4,  TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 6, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS161B);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
  
  
      /* Equation (162) first four terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS162A =
        { {contract1PDMSL.Y().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 3, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 1, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr4, HerDen, GAUNT_EXCHANGE, 2, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS162A);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,   -iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
  
  
      /* Equation (162) last three terms */
      std::vector<TwoBodyContraction<MatsT>> contractGLS162B =
        { {contract1PDMSL.Y().pointer(), Scr1, HerDen, GAUNT_EXCHANGE, 12, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUNT_EXCHANGE, 4,  TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUNT_EXCHANGE, 8, TRANS_KL} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLS162B);
  
      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
  
  
#ifdef _PRINT_MATRICES
      std::cout<<"After Gaunt 159-162"<<std::endl;
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif //_PRINT_MATRICES
  
#endif // Gaunt LLSS EXCHANGE
  
#ifdef _PRINT_MATRICES
  
      std::cout<<"After Gaunt LLSS|SSLL"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
  
#endif //_PRINT_MATRICES
  
  
      /*------------------------------------*/
      /*   End of Gaunt (LL|SS) Contraction */
      /*------------------------------------*/
      } // computeExchange, LLSS exchange part 
  
#ifdef _FOCK_CONTRACTION_TIME
      auto durERIDG = tock(topERIDG);
      std::cout << "Gaunt Contraction duration   = " << durERIDG << std::endl;
#endif

    } //GAUNT








    /**********************************************/
    /*                                            */
    /*              SSSS                          */
    /*                                            */
    /**********************************************/


    //ERI00:    (∇A∙∇B)(∇C∙∇D)(ij|kl)
    //ERI01:    (∇Ax∇B)_x(∇C∙∇D)(ijkl)
    //ERI02:    (∇Ax∇B)_y(∇C∙∇D)(ijkl)
    //ERI03:    (∇Ax∇B)_z(∇C∙∇D)(ijkl)
    //ERI04:    (∇A∙∇B)(∇Cx∇D)_x(ijkl)
    //ERI05:    (∇A∙∇B)(∇Cx∇D)_y(ijkl)
    //ERI06:    (∇A∙∇B)(∇Cx∇D)_z(ijkl)
    //ERI07:    (∇Ax∇B)_x(∇Cx∇D)_x(ijkl)
    //ERI08:    (∇Ax∇B)_x(∇Cx∇D)_y(ijkl)
    //ERI09:    (∇Ax∇B)_x(∇Cx∇D)_z(ijkl)
    //ERI10:    (∇Ax∇B)_y(∇Cx∇D)_x(ijkl)
    //ERI11:    (∇Ax∇B)_y(∇Cx∇D)_y(ijkl)
    //ERI12:    (∇Ax∇B)_y(∇Cx∇D)_z(ijkl)
    //ERI13:    (∇Ax∇B)_z(∇Cx∇D)_x(ijkl)
    //ERI14:    (∇Ax∇B)_z(∇Cx∇D)_y(ijkl)
    //ERI15:    (∇Ax∇B)_z(∇Cx∇D)_z(ijkl)





    if(this->hamiltonianOptions_.DiracCoulombSSSS) {//Dirac-Coulomb-SSSS

      double C4 = 1./(16*SpeedOfLight()*SpeedOfLight()*SpeedOfLight()*SpeedOfLight());
      dcomplex scaleC4 = dcomplex(C4,0.0);
      dcomplex iscaleC4 = dcomplex(0.0,C4);

#ifdef _FOCK_CONTRACTION_TIME
      auto topERIDCSSSS = tick();
#endif

#ifdef _SS_J_TIME
      auto topSS_J1 = tick();
#endif

      /* Equation 70 in the paper */
      std::vector<TwoBodyContraction<MatsT>> contractSSSS70A =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_COULOMB,  0},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_COULOMB,  4},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, SSSS_COULOMB,  5},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, SSSS_COULOMB,  6} };


      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70A);
 

      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*scaleC4, Scr1, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+SS, NB2C, ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscaleC4, Scr2, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+SS, NB2C, ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscaleC4, Scr3, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+SS, NB2C, ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*iscaleC4, Scr4, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+SS, NB2C, ss.twoeH->S().pointer()+SS, NB2C);

#ifdef _SS_J_TIME
      // Print out BareCoulomb duration 
      auto durSS_J1 = tock(topSS_J1);
      //std::cout << "J Contribution to F_SS = " << topSS_J1 << std::endl;
#endif 

#ifdef _SS_K_TIME
      auto topSS_K1 = tick();
#endif
      if (computeExchange) {
#if 1
      std::vector<TwoBodyContraction<MatsT>> contractSSSS70B =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  0},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  7},
          {contract1PDMSS.S().pointer(), Scr3, HerDen, SSSS_EXCHANGE, 11},
          {contract1PDMSS.S().pointer(), Scr4, HerDen, SSSS_EXCHANGE, 15} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70B);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS70C1 =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  1},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  2},
          {contract1PDMSS.Z().pointer(), Scr3, HerDen, SSSS_EXCHANGE,  3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70C1);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS70C2 =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  4},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  5},
          {contract1PDMSS.Z().pointer(), Scr3, HerDen, SSSS_EXCHANGE,  6} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70C2);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS70D =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  8},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 10} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70D);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS70E =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 13},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  9} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70E);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS70F =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 12},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 14} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS70F);
 

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);

#endif
      } // computeExchange

#ifdef _SS_K_TIME
      auto durSS_K1 = tock(topSS_K1);
#endif 

#ifdef _SS_J_TIME
      auto topSS_J2 = tick();
#endif
      /* Equation 71 in the paper */
      std::vector<TwoBodyContraction<MatsT>> contractSSSS71A =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_COULOMB,  3},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_COULOMB, 13},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, SSSS_COULOMB, 14},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, SSSS_COULOMB, 15} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71A);
 

      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscaleC4, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+SS, NB2C, ss.twoeH->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+SS, NB2C, ss.twoeH->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr3, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+SS, NB2C, ss.twoeH->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr4, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+SS, NB2C, ss.twoeH->Z().pointer()+SS, NB2C);

#ifdef _SS_J_TIME
      // Print out BareCoulomb duration 
      auto durSS_J2 = tock(topSS_J2);
      //std::cout << "J Contribution to F_SS = " << topSS_J1 << std::endl;
#endif 

#ifdef _SS_K_TIME
      auto topSS_K2 = tick();
#endif
      if (computeExchange) {
#if 1
      std::vector<TwoBodyContraction<MatsT>> contractSSSS71B =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  0},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 15} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71B);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71C =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  7},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 11} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71C);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71D =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  6},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71D);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71G =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  8},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 10} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71G);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71E =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  4},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71E);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71I =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 12},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 14} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71I);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71F =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  5},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71F);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS71H =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  9},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 13} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS71H);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);

#endif
      } // computeExchange

#ifdef _SS_K_TIME
      auto durSS_K2 = tock(topSS_K2);
#endif 


#ifdef _SS_J_TIME
      auto topSS_J3 = tick();
#endif

      /* Equation 72 in the paper */
      std::vector<TwoBodyContraction<MatsT>> contractSSSS72A =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_COULOMB,  1},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_COULOMB,  7},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, SSSS_COULOMB,  8},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, SSSS_COULOMB,  9} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72A);
 

      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C,  2.0*iscaleC4, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+SS, NB2C, ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scaleC4, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+SS, NB2C, ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scaleC4, Scr3, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+SS, NB2C, ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scaleC4, Scr4, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+SS, NB2C, ss.twoeH->X().pointer()+SS, NB2C);


#ifdef _SS_J_TIME
      // Print out BareCoulomb duration 
      auto durSS_J3 = tock(topSS_J3);
      //std::cout << "J Contribution to F_SS = " << topSS_J1 << std::endl;
#endif 


#ifdef _SS_K_TIME
      auto topSS_K3 = tick();
#endif
      if (computeExchange) {
#if 1
      std::vector<TwoBodyContraction<MatsT>> contractSSSS72B =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  0},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  7} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72B);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72C =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 11},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 15} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72C);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72D =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  4},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72D);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72G =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 12},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 14} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72G);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72E =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  6},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72E);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72H =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  8},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 10} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72H);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72F =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  5},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72F);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS72I =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  9},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 13} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS72I);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);

#endif
      } // computeExchange

#ifdef _SS_K_TIME
      auto durSS_K3 = tock(topSS_K3);
#endif 

#ifdef _SS_J_TIME
      auto topSS_J4 = tick();
#endif

      /* Equation 73 in the paper */
      std::vector<TwoBodyContraction<MatsT>> contractSSSS73A =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_COULOMB,  2},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_COULOMB, 10},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, SSSS_COULOMB, 11},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, SSSS_COULOMB, 12} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73A);
 

      // Assemble 4C twoeH 
      MatAdd('N','N', NB1C, NB1C, 2.0*iscaleC4, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+SS, NB2C, ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+SS, NB2C, ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr3, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+SS, NB2C, ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-2.0*scaleC4, Scr4, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+SS, NB2C, ss.twoeH->Y().pointer()+SS, NB2C);

#ifdef _SS_J_TIME
      // Print out BareCoulomb duration 
      auto durSS_J4 = tock(topSS_J4);
      //std::cout << "J Contribution to F_SS = " << topSS_J1 << std::endl;
#endif 


#ifdef _SS_K_TIME
      auto topSS_K4 = tick();
#endif
      if (computeExchange) {
#if 1
      std::vector<TwoBodyContraction<MatsT>> contractSSSS73B =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  0},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 11} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73B);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,  scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73C =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  7},
          {contract1PDMSS.Y().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 15} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73C);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73D =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  5},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  2} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73D);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73G =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  9},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 13} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73G);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, iscaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73E =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  6},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  3} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73E);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73H =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  8},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 10} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73H);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73F =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE,  4},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE,  1} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73F);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C,-scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractSSSS73I =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, SSSS_EXCHANGE, 12},
          {contract1PDMSS.Z().pointer(), Scr2, HerDen, SSSS_EXCHANGE, 14} };
  
      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractSSSS73I);

      // Assemble 4C exchangeMatrix 
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scaleC4, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);

#endif
     } // computeExchange

#ifdef _SS_K_TIME
      // Print out BareCoulomb duration 
      auto durSS_K4 = tock(topSS_K4);
#endif 

#ifdef _SS_J_TIME
      // Print out SS Coulomb duration 
      std::cout << "SSSS J F_SS duration   = " << durSS_J1 + durSS_J2 + durSS_J3 + durSS_J4 << std::endl;
#endif 

#ifdef _SS_K_TIME
      // Print out SS Exchange duration 
      std::cout << "SSSS K F_SS duration   = " << durSS_K1 + durSS_K2 + durSS_K3 + durSS_K4 << std::endl;
#endif 
      
#ifdef _FOCK_CONTRACTION_TIME
      auto durERIDCSSSS = tock(topERIDCSSSS);
      std::cout << "SSSS duration   = " << durERIDCSSSS << std::endl;
#endif


    } ////Dirac-Coulomb-SSSS



    /**********************************************/
    /*                                            */
    /*              GAUGE                         */
    /*                                            */
    /**********************************************/

    //ERI00:    (ss)(ij|kl)
    //ERI01:    (sσ)_x(ijkl)
    //ERI02:    (sσ)_y(ijkl)
    //ERI03:    (sσ)_z(ijkl)
    //ERI04:    (σs)_x(ijkl)
    //ERI05:    (σs)_y(ijkl)
    //ERI06:    (σs)_z(ijkl)
    //ERI07:    (σ∙σ)(ijkl)
    //ERI08:    (σxσ)_x(ijkl)
    //ERI09:    (σxσ)_y(ijkl)
    //ERI10:    (σxσ)_z(ijkl)
    //ERI11:    (σσ)(xx - yy - zz)(ijkl)
    //ERI12:    (σσ)(-xx + yy - zz)(ijkl)
    //ERI13:    (σσ)(-xx - yy + zz)(ijkl)
    //ERI14:    (σ_x σ_y + σ_y σ_x)(ijkl)
    //ERI15:    (σ_z σ_x + σ_x σ_z)(ijkl)
    //ERI16:    (σ_y σ_z + σ_z σ_y)(ijkl)
    //ERI17:    (σ_x σ_x)(ijkl)
    //ERI18:    (σ_x σ_y)(ijkl)
    //ERI19:    (σ_x σ_z)(ijkl)
    //ERI20:    (σ_y σ_x)(ijkl)
    //ERI21:    (σ_y σ_y)(ijkl)
    //ERI22:    (σ_y σ_z)(ijkl)
    //ERI23:    (σ_z σ_x)(ijkl)
    //ERI24:    (σ_z σ_y)(ijkl)
    //ERI25:    (σ_z σ_z)(ijkl)

    // Gauge contribution
    if(this->hamiltonianOptions_.Gauge) {//Gauge


#ifdef _FOCK_CONTRACTION_TIME
      auto topERIGauge = tick();
#endif


      scale = 0.5* 1./(4*SpeedOfLight()*SpeedOfLight());
      iscale = 0.5* dcomplex(0.0, 1./(4*SpeedOfLight()*SpeedOfLight()));

      dcomplex scalef = 0.5* 1./(4*SpeedOfLight()*SpeedOfLight());
      dcomplex iscalef = 0.5* dcomplex(0.0, 1./(4*SpeedOfLight()*SpeedOfLight()));


      if (computeExchange) {
      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gauge (LL|LL) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/


      /* Equation (208) */
      // 1st Line of Equation (208)
      std::vector<TwoBodyContraction<MatsT>> contractGLL208 =
      { {contract1PDMSS.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0},
        {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1},
        {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2},
        {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 3} };

      // Call the contraction engine to do the assembly
      ss.TPI->twoBodyContract(ss.comm, contractGLL208);

      // Add to the LL part of 4C exchangeMatrix in Pauli matrix form
      MatAdd('N','N', NB1C, NB1C,-scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);


      // 2nd Line of Equation (208)
      std::vector<TwoBodyContraction<MatsT>> contractGLL2082 =
        { {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 6} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2082);

      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);

      // 2rd Line of Equation (208)
      std::vector<TwoBodyContraction<MatsT>> contractGLL2083 =
        { {contract1PDMSS.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 7},
          {contract1PDMSS.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 9},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 10} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2083);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer(), NB2C, ss.exchangeMatrix->S().pointer(), NB2C);

     /* Equation (209) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL209 =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 3},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2},
          {contract1PDMSS.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL209);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLL2092 =
        { {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 6},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5},
          {contract1PDMSS.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2092);
  
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLL2093 =
        { {contract1PDMSS.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 13},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 10},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 15},
          {contract1PDMSS.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2093);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer(), NB2C, ss.exchangeMatrix->Z().pointer(), NB2C);


      /* Equation (210) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL210 =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 2} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL210);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLL2102 =
        { {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4},
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 5} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2102);
  
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
  

      std::vector<TwoBodyContraction<MatsT>> contractGLL2103 =
        { {contract1PDMSS.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 11},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8},  // changed here
          {contract1PDMSS.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 15} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2103);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer(), NB2C, ss.exchangeMatrix->X().pointer(), NB2C);
  

      /* Equation (211) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL211 =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 2},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL211);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLL2112 =
        { {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 5},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2112);
  
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
  

      std::vector<TwoBodyContraction<MatsT>> contractGLL2113 =
        { {contract1PDMSS.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 12},
          {contract1PDMSS.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 9},
          {contract1PDMSS.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14},
          {contract1PDMSS.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLL2113);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer(), NB2C, ss.exchangeMatrix->Y().pointer(), NB2C);
 
      /*++++++++++++++++++++++++++++++++++++*/
      /*  End of Gauge (LL|LL) Contraction  */
      /*++++++++++++++++++++++++++++++++++++*/



      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gauge (SS|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/
  

      /* Equation (227) */
    std::vector<TwoBodyContraction<MatsT>> contractGSS227 =
      { {contract1PDMLL.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0,TRANS_MN_TRANS_KL},
        {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 3,TRANS_MN_TRANS_KL} };
  
    ss.TPI->twoBodyContract(ss.comm, contractGSS227);
  
    MatAdd('N','N', NB1C, NB1C,-scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C,iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C,iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C,iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);

 
    std::vector<TwoBodyContraction<MatsT>> contractGSS2272 =
      { {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 6,TRANS_MN_TRANS_KL} };
  
    ss.TPI->twoBodyContract(ss.comm, contractGSS2272);
  
    MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    
    std::vector<TwoBodyContraction<MatsT>> contractGSS2273 =
      { {contract1PDMLL.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 7,TRANS_MN_TRANS_KL},
        {contract1PDMLL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 9,TRANS_MN_TRANS_KL},
        {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 10,TRANS_MN_TRANS_KL} };
  
    ss.TPI->twoBodyContract(ss.comm, contractGSS2273);
  
    MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);
    MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+SS, NB2C, ss.exchangeMatrix->S().pointer()+SS, NB2C);


      /* Equation (228) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS228 =
        { {contract1PDMLL.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 3,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS228);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGSS2282 =
        { {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 6,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2282);
  
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGSS2283 =
        { {contract1PDMLL.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 13,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 10,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 15,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2283);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+SS, NB2C, ss.exchangeMatrix->Z().pointer()+SS, NB2C);



      /* Equation (229) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS229 =
        { {contract1PDMLL.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 2,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS229);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
  

      std::vector<TwoBodyContraction<MatsT>> contractGSS2292 =
        { {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 5,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2292);
  
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
  

      std::vector<TwoBodyContraction<MatsT>> contractGSS2293 =
        { {contract1PDMLL.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 11,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 15,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2293);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+SS, NB2C, ss.exchangeMatrix->X().pointer()+SS, NB2C);
  


      /* Equation (230) */
      std::vector<TwoBodyContraction<MatsT>> contractGSS230 =
        { {contract1PDMLL.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 2,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS230);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGSS2302 =
        { {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 5,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2302);
  
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);

      std::vector<TwoBodyContraction<MatsT>> contractGSS2303 =
        { {contract1PDMLL.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 12,TRANS_MN_TRANS_KL},
          {contract1PDMLL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 9,TRANS_MN_TRANS_KL},
          {contract1PDMLL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14,TRANS_MN_TRANS_KL},
          {contract1PDMLL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16,TRANS_MN_TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGSS2303);
  
      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+SS, NB2C, ss.exchangeMatrix->Y().pointer()+SS, NB2C);

      /*++++++++++++++++++++++++++++++++++++*/
      /*  End of Gauge (SS|SS) Contraction  */
      /*++++++++++++++++++++++++++++++++++++*/

      } // computeExchange

      /*++++++++++++++++++++++++++++++++++++*/
      /* Start of Gauge (LL|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++*/

      /* Equation (232) */
      std::vector<TwoBodyContraction<MatsT>> contractGLS232 =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 0},
          {contract1PDMLS.X().pointer(), Scr2, HerDen, GAUGE_COULOMB, 1},
          {contract1PDMLS.Y().pointer(), Scr3, HerDen, GAUGE_COULOMB, 2},
          {contract1PDMLS.Z().pointer(), Scr4, HerDen, GAUGE_COULOMB, 3} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS232);
  
      MatAdd('N','N', NB1C, NB1C,  -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);

      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLS2322 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 0, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUGE_COULOMB, 1, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_COULOMB, 2, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_COULOMB, 3, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS2322);
  
      MatAdd('N','N', NB1C, NB1C,  2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr2, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr3, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr4, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+LS, NB2C);
 
      /* Equation (233)X */
      std::vector<TwoBodyContraction<MatsT>> contractGLS233AX =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 4} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233AX);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  

      std::vector<TwoBodyContraction<MatsT>> contractGLS233BX =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 17},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 18},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 19} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233BX);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLS233CX =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 4, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233CX);
  
      MatAdd('N','N', NB1C, NB1C,  2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
  
  
      std::vector<TwoBodyContraction<MatsT>> contractGLS233DX =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 17, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 18, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 19, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233DX);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+LS, NB2C);

      /* Equation (233)Y  */
      std::vector<TwoBodyContraction<MatsT>> contractGLS233AY =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 5} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233AY);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLS233BY =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 20},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 21},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 22} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233BY);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLS233CY =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 5, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233CY);
  
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
  
  
      std::vector<TwoBodyContraction<MatsT>> contractGLS233DY =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 20, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 21, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 22, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233DY);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+LS, NB2C);



      /* Equation (233)Z  */
      std::vector<TwoBodyContraction<MatsT>> contractGLS233AZ =
        { {contract1PDMLS.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 6} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233AZ);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
 

      std::vector<TwoBodyContraction<MatsT>> contractGLS233BZ =
        { {contract1PDMLS.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 23},
          {contract1PDMLS.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 24},
          {contract1PDMLS.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 25} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233BZ);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
  
  
      std::vector<TwoBodyContraction<MatsT>> contractGLS233CZ =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_COULOMB, 6, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233CZ);
  
      MatAdd('N','N', NB1C, NB1C, 2.0*iscale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
  
  
      std::vector<TwoBodyContraction<MatsT>> contractGLS233DZ =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUGE_COULOMB, 23, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr2, HerDen, GAUGE_COULOMB, 24, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr3, HerDen, GAUGE_COULOMB, 25, TRANS_KL} };
  
      ss.TPI->twoBodyContract(ss.comm, contractGLS233DZ);
  
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr1, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr2, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -2.0*scale, Scr3, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+LS, NB2C);
      
      /* Gauge (LL|SS) Contraction Exchange part */
      if (computeExchange) {

      /* Equation (248) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL248 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 3, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL248);

      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2482 =
        { {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 6, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2482);

      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2483 =
        { {contract1PDMSL.S().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 7, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 9, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 10, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2483);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, iscale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+LS, NB2C);


      /* Equation (249) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL249 =
        { {contract1PDMSL.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 3, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL249);

      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2492 =
        { {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 6, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 5, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2492);

      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2493 =
        { {contract1PDMSL.Z().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 13, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 10, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 15, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2493);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+LS, NB2C);


      /* Equation (250) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL250 =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 1, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 2, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL250);

      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2502 =
        { {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 4, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 5, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2502);

      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2503 =
        { {contract1PDMSL.X().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 11, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 8, TRANS_KL},
          {contract1PDMSL.Y().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 15, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2503);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+LS, NB2C);


      /* Equation (251) */
      std::vector<TwoBodyContraction<MatsT>> contractGLL251 =
        { {contract1PDMSL.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 0, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 2, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 3, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 1, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL251);

      MatAdd('N','N', NB1C, NB1C,  scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2512 =
        { {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 5, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 6, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 4, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2512);

      MatAdd('N','N', NB1C, NB1C, iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,  scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);


      std::vector<TwoBodyContraction<MatsT>> contractGLL2513 =
        { {contract1PDMSL.Y().pointer(), Scr1, HerDen, GAUGE_EXCHANGE, 12, TRANS_KL},
          {contract1PDMSL.S().pointer(), Scr2, HerDen, GAUGE_EXCHANGE, 9, TRANS_KL},
          {contract1PDMSL.X().pointer(), Scr3, HerDen, GAUGE_EXCHANGE, 14, TRANS_KL},
          {contract1PDMSL.Z().pointer(), Scr4, HerDen, GAUGE_EXCHANGE, 16, TRANS_KL} };

      ss.TPI->twoBodyContract(ss.comm, contractGLL2513);

      MatAdd('N','N', NB1C, NB1C, -scale, Scr1, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C,-iscale, Scr2, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr3, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -scale, Scr4, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      
      } // computeExchange

#ifdef _PRINT_MATRICES
      std::cout<<"After Gauge LLLL"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif


#ifdef _FOCK_CONTRACTION_TIME
      auto durERIGauge = tock(topERIGauge);
      std::cout << "Gauge duration   = " << durERIGauge << std::endl;
#endif



    } // Gauge
 


    /*******************************/
    /* Final Assembly of 4C Matrix */
    /*******************************/
    ROOT_ONLY(ss.comm);
    
    if (computeExchange) {
    // Copy LS to SL part of the exchangeMatrix[MS]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MX]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MY]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MZ]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+SL, NB2C);
    }
    
    if (HerDen) {
      // Hermitrize matrices to avoid accumulate small rounding errors

      // Copy LS to SL part of the twoeH[MS]
      SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MX]
      SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MY]
      SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MZ]
      SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+SL, NB2C);
    } else {
      // only use symmetry of the integrals here

      // Copy LS to SL part of the twoeH[MS]
      SetMat('T', NB1C, NB1C, -MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MX]
      SetMat('T', NB1C, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MY]
      SetMat('T', NB1C, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+SL, NB2C);
      // Copy LS to SL part of the twoeH[MZ]
      SetMat('T', NB1C, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+SL, NB2C);
    } 

    // Form GD: G[D] = 2.0*J[D] - K[D]
    if( computeExchange ){
      if( do_VLL)  {
        *ss.twoeH -= *ss.exchangeMatrix;
      } else {
        *ss.twoeH -= xHFX * *ss.exchangeMatrix;
      } 
    }
    // G[D] += 2*J[D]
    *ss.twoeH += 2.0 * *ss.coulombMatrix;

    CQMemManager::get().free(Scr1);
    CQMemManager::get().free(Scr2);
    CQMemManager::get().free(Scr3);
    CQMemManager::get().free(Scr4);


#ifdef _PRINT_MATRICES

    prettyPrintSmart(std::cout,"twoeH MS",ss.twoeH->S().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MX",ss.twoeH->X().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MY",ss.twoeH->Y().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MZ",ss.twoeH->Z().pointer(),NB2C,NB2C,NB2C);


    MatsT* TEMP_GATHER1 = CQMemManager::get().malloc<MatsT>(NB4C2);
    MatsT* TEMP_GATHER2 = CQMemManager::get().malloc<MatsT>(NB4C2);

    memset(TEMP_GATHER1,0.,NB4C2*sizeof(MatsT));
    memset(TEMP_GATHER2,0.,NB4C2*sizeof(MatsT));

    std::cout << std::scientific << std::setprecision(16);
    SpinGather(NB2C,TEMP_GATHER1,NB4C,contract1PDM.S().pointer(),NB2C,contract1PDM.Z().pointer(),NB2C,contract1PDM.Y().pointer(),NB2C,contract1PDM.X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"density Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);


    SpinGather(NB2C,TEMP_GATHER2,NB4C,ss.twoeH->S().pointer(),NB2C,ss.twoeH->Z().pointer(),NB2C,ss.twoeH->Y().pointer(),NB2C,ss.twoeH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"twoeH Gather",TEMP_GATHER2,NB4C,NB4C,NB4C,1,12,16);

    SpinGather(NB2C,TEMP_GATHER1,NB4C,ss.coreH->S().pointer(),NB2C,ss.coreH->Z().pointer(),NB2C,ss.coreH->Y().pointer(),NB2C,ss.coreH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"coreH Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);
 
    CQMemManager::get().free(TEMP_GATHER1);
    CQMemManager::get().free(TEMP_GATHER2);

#endif //_PRINT_MATRICES


  }; // FourCompFock<MatsT, IntsT>::formGDInCore


  template <>
  void FourCompFock<double,double>::formGDInCore(SingleSlater<double,double> &ss,
                                                 EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {
    CErr("Real number Four-Component NYI.");
  }


  template <>
  void FourCompFock<dcomplex,dcomplex>::formGDInCore(SingleSlater<dcomplex,dcomplex> &ss,
                                                     EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {
    CErr("Complex integral Four-Component NYI.");
  }





#ifdef HAS_FORMGD_3INDEX
  /**   
   *  \brief Forms the 4C Fock matrix using 3 Index ERI
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formGD3Index(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {

    GTODirectRelERIContraction<MatsT,IntsT> &relERICon =
        *std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI);

    // Decide list of onePDMs to use
    cqmatrix::PauliSpinorMatrices<MatsT> &contract1PDM
        = increment ? *ss.deltaOnePDM : *ss.onePDM;

    size_t NB1C  = ss.basisSet().nBasis;
    size_t NB2C  = 2 * ss.basisSet().nBasis;
    size_t NB4C  = 4 * ss.basisSet().nBasis;
    size_t NB1C2 = NB1C*NB1C;
    size_t NB1C4 = NB1C*NB1C*NB1C*NB1C;
    size_t NB1C3 = NB1C*NB1C*NB1C;
    size_t NB2C2 = NB2C*NB2C;
    size_t NB4C2 = NB4C*NB4C;

    size_t SS = NB2C*NB1C+NB1C;
    size_t LS = NB2C*NB1C;
    size_t SL = NB1C;

    auto MS = SCALAR;

    size_t mpiRank   = MPIRank(ss.comm);
    bool   isNotRoot = mpiRank != 0;

    cqmatrix::PauliSpinorMatrices<MatsT> exchangeMatrixLL(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLL(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSL(NB1C);

    MatsT* Scr1 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr2 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr3 = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* Scr4 = CQMemManager::get().malloc<MatsT>(NB1C2);
    memset(Scr1,0.,NB1C2*sizeof(MatsT));
    memset(Scr2,0.,NB1C2*sizeof(MatsT));
    memset(Scr3,0.,NB1C2*sizeof(MatsT));
    memset(Scr4,0.,NB1C2*sizeof(MatsT));


    // Compute 1/(2mc)^2
    //dcomplex scale = 1.;
    //dcomplex iscale = dcomplex(0.0, 1.0);
    dcomplex scale = 1./(4*SpeedOfLight()*SpeedOfLight());
    dcomplex iscale = dcomplex(0.0, 1./(4*SpeedOfLight()*SpeedOfLight()));

    for(size_t i = 0; i < contract1PDM.nComponent(); i++) {
      cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer(),    NB2C,
             contract1PDMLL[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SS, NB2C,
             contract1PDMSS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+LS, NB2C,
             contract1PDMLS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SL, NB2C,
             contract1PDMSL[c].pointer(), NB1C);
    }

#ifdef _PRINT_MATRICES
    prettyPrintSmart(std::cout, "1PDM[MS]", contract1PDM.S().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MX]", contract1PDM.X().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MY]", contract1PDM.Y().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MZ]", contract1PDM.Z().pointer(), NB2C, NB2C, NB2C);
#endif

 
    if(not increment) {
      ss.coulombMatrix->clear();
      ss.exchangeMatrix->clear();
    };


    /**********************************************/
    /*                                            */
    /*              DIRECT COULOMB     	          */
    /*                                            */
    /**********************************************/


    if(this->hamiltonianOptions_.BareCoulomb) { // DIRECT_COULOMB

#ifdef _FOCK_CONTRACTION_TIME
      auto topBareCoulomb = tick();
#endif
      /*+++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Direct Coulomb (LL|LL) Contraction */
      /*+++++++++++++++++++++++++++++++++++++++++++++*/
  
      std::vector<TwoBodyContraction<MatsT>> contractLL =
        { {contract1PDMLL.S().pointer(), Scr1, HerDen, COULOMB} };
  
      // Determine how many (if any) exchange terms to calculate
      if( std::abs(xHFX) > 1e-12 ) {
        exchangeMatrixLL.clear();
        for(size_t i = 0; i < ss.exchangeMatrix->nComponent(); i++) {
  
          cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
          contractLL.push_back(
            {contract1PDMLL[c].pointer(), exchangeMatrixLL[c].pointer(), HerDen, EXCHANGE}
          );
        }
      }
  
      // Zero out K[i]
      if(not increment) ss.exchangeMatrix->clear();
  
      // Call the contraction engine to do the assembly of direct Coulomb LLLL
      GTODirectTPIContraction<MatsT,IntsT>(ss.TPI->ints()).twoBodyContract(ss.comm, HerDen, contractLL, pert);
  
  
      /* Store LL block into 2C spin scattered matrices */
      // Assemble 4C coulombMatrix
      SetMat('N', NB1C, NB1C, MatsT(1.), Scr1, NB1C, ss.coulombMatrix->pointer(), NB2C);
  
      // Assemble 4C exchangeMatrix 
      for(auto i = 0; i < ss.exchangeMatrix->nComponent();i++){
        cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
        SetMat('N', NB1C, NB1C, MatsT(1.), exchangeMatrixLL[c].pointer(), NB1C,
               (*ss.exchangeMatrix)[c].pointer(), NB2C);
      }

      /*---------------------------------------------*/
      /*   End of Direct Coulomb (LL|LL) Contraction */
      /*---------------------------------------------*/

      // Print out BareCoulomb duration 
#ifdef _FOCK_CONTRACTION_TIME
      auto durBareCoulomb = tock(topBareCoulomb);
      std::cout << "Non-relativistic Coulomb duration = " << durBareCoulomb << std::endl;
#endif

    } // DIRECT_COULOMB


    /* using 3-Index ERI */
    // Loop over first shell
    size_t n1, bf1, bf1_s;
    size_t nERI3 = 37;
    for(auto s1(0), bf1_s(0); s1< ss.basisSet().nShell; bf1_s+=n1, s1++) {

      n1 = ss.basisSet().shells[s1].size();

      relERICon.computeERI3Index(s1);

    // Loop over all basis in the first shell
    for(auto ibatch = 0 ; ibatch < n1; ibatch++){
      auto ERI4bf1 = relERICon.ERI4DCB+nERI3*NB1C3*ibatch;
      bf1 = bf1_s + ibatch;


      /**********************************************/
      /*                                            */
      /*              DIRAC-COULOMB    	            */
      /*                                            */
      /**********************************************/

      // ERI: (ab|cd)
      // ERIDCB0: ∇A∙∇B(ab|cd)
      // ERIDCB1: ∇Ax∇B(ab|cd)-X
      // ERIDCB2: ∇Ax∇B(ab|cd)-Y
      // ERIDCB3: ∇Ax∇B(ab|cd)-Z

      if(this->hamiltonianOptions_.DiracCoulomb) { // DIRAC_COULOMB
  
        /*++++++++++++++++++++++++++++++++++++++++++++*/
        /* Start of Dirac-Coulomb (LL|LL) Contraction */
        /*++++++++++++++++++++++++++++++++++++++++++++*/
    
        std::vector<TwoBodyContraction<MatsT>> contractDCLL =
          { {contract1PDMSS.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+4*NB1C3},
            {contract1PDMSS.X().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+5*NB1C3},
            {contract1PDMSS.Y().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+6*NB1C3},
            {contract1PDMSS.Z().pointer(), Scr4, HerDen, COULOMB, .ERI4=ERI4bf1+7*NB1C3} };
    
        // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
        relERICon.twoBodyContract3Index(ss.comm, contractDCLL);
    
        // Add Dirac-Coulomb contributions  to the LLLL block
        for(auto i=0; i<NB1C; i++){
          ss.coulombMatrix->pointer()[bf1+i*NB2C] += scale*Scr1[i] + iscale*Scr2[i] + iscale*Scr3[i] + iscale*Scr4[i];
        }
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After LLLL Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
        /*------------------------------------------*/
        /* End of Dirac-Coulomb (LL|LL) Contraction */
        /*------------------------------------------*/
    
    
    
    
    
        /*++++++++++++++++++++++++++++++++++++++++++++*/
        /* Start of Dirac-Coulomb (SS|SS) Contraction */
        /*++++++++++++++++++++++++++++++++++++++++++++*/
    
        std::vector<TwoBodyContraction<MatsT>> contractSS =
          { {contract1PDMLL.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1},
            {contract1PDMLL.S().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+NB1C3},
            {contract1PDMLL.S().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+2*NB1C3},
            {contract1PDMLL.S().pointer(), Scr4, HerDen, COULOMB, .ERI4=ERI4bf1+3*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractSS);
    
        // Store SS block into 2C spin scattered matrices 
        // These scaling factors were modified to take into account the issue of storing the 
        // Coulomb portion in the exchange matrix, this will be fixed later
        for(auto i=0; i<NB1C; i++){
          ss.coulombMatrix->pointer()[SS+bf1+i*NB2C]      +=       scale*Scr1[i];
          ss.exchangeMatrix->X().pointer()[SS+bf1+i*NB2C] += -2.0*iscale*Scr2[i];
          ss.exchangeMatrix->Y().pointer()[SS+bf1+i*NB2C] += -2.0*iscale*Scr3[i];
          ss.exchangeMatrix->Z().pointer()[SS+bf1+i*NB2C] += -2.0*iscale*Scr4[i];
        }
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After SSSS Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
    
        /*--------------------------------------------*/
        /* End of Dirac-Coulomb (SS|SS) Contraction */
        /*--------------------------------------------*/
    
    
    
    
    
        /*++++++++++++++++++++++++++++++++++++++++++*/
        /* Start of Dirac-Coulomb (LL|SS) / (SS|LL) */
        /*++++++++++++++++++++++++++++++++++++++++++*/
    
        std::vector<TwoBodyContraction<MatsT>> contractLSScalar =
          { {contract1PDMLS.S().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+4*NB1C3},
            {contract1PDMLS.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+5*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+6*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+7*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractLSScalar);
    
        // Add to the LS part of the exchangeMatrix[MS]
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +iscale*Scr2[i] +iscale*Scr3[i] +iscale*Scr4[i];
        }
    
    
    
        std::vector<TwoBodyContraction<MatsT>> contractLSMX =
          { {contract1PDMLS.X().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+4*NB1C3},
            {contract1PDMLS.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+5*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+6*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+7*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractLSMX);
    
        // Add to the LS part of the exchangeMatrix[MX]
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +iscale*Scr2[i] +scale*Scr3[i] -scale*Scr4[i];
        }
    
    
    
    
        std::vector<TwoBodyContraction<MatsT>> contractLSMY =
          { {contract1PDMLS.Y().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+4*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+5*NB1C3},
            {contract1PDMLS.S().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+6*NB1C3},
            {contract1PDMLS.X().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+7*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractLSMY);
    
        // Add to the LS part of the exchangeMatrix[MY]
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] -scale*Scr2[i] +iscale*Scr3[i] +scale*Scr4[i];
        }
    
    
    
    
    
        std::vector<TwoBodyContraction<MatsT>> contractLSMZ =
          { {contract1PDMLS.Z().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+4*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+5*NB1C3},
            {contract1PDMLS.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+6*NB1C3},
            {contract1PDMLS.S().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+7*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractLSMZ);
    
        // Add to the LS part of the exchangeMatrix[MZ]
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +scale*Scr2[i] -scale*Scr3[i] +iscale*Scr4[i];
        }
    
    
        /*------------------------------------------*/
        /*   End of Dirac-Coulomb (LL|SS) / (SS|LL) */
        /*------------------------------------------*/
    
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After LLSS Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
  
      } //_DIRAC_COULOMB
  



      /**********************************************/
      /*                                            */
      /*              GAUNT                         */
      /*                                            */
      /**********************************************/
    
    
      //ERI4:    ∇B∙∇C(mn|kl)
      //ERI5 :   ∇Bx∇C(mn|kl)-X
      //ERI6 :   ∇Bx∇C(mn|kl)-Y
      //ERI7 :   ∇Bx∇C(mn|kl)-Z
      //ERI8 :   ∇B_x∇C_y(mn|kl) + ∇B_y∇C_x(mn|kl)
      //ERI9 :   ∇B_y∇C_x(mn|kl)
      //ERI10:   ∇B_x∇C_z(mn|kl) + ∇B_z∇C_x(mn|kl)
      //ERI11:   ∇B_z∇C_x(mn|kl)
      //ERI12:   ∇B_y∇C_z(mn|kl) + ∇B_z∇C_y(mn|kl)
      //ERI13:   ∇B_z∇C_y(mn|kl)
      //ERI14: - ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) + ∇B_z∇C_z(mn|kl)
      //ERI15:   ∇B_x∇C_x(mn|kl) - ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
      //ERI16: - ∇B_x∇C_x(mn|kl) + ∇B_y∇C_y(mn|kl) - ∇B_z∇C_z(mn|kl)
      //ERI17:   ∇B_x∇C_x(mn|kl)
      //ERI18:   ∇B_x∇C_y(mn|kl)
      //ERI19:   ∇B_x∇C_z(mn|kl)
      //ERI20:   ∇B_y∇C_y(mn|kl)
      //ERI21:   ∇B_y∇C_z(mn|kl)
      //ERI22:   ∇B_z∇C_z(mn|kl)
    
      if(this->hamiltonianOptions_.Gaunt) {//GAUNT
    
#if 1 //Gaunt LLLL
        /*++++++++++++++++++++++++++++++++++++*/
        /* Start of Gaunt (LL|LL) Contraction */
        /*++++++++++++++++++++++++++++++++++++*/
    
        /* Equation (113) */
        std::vector<TwoBodyContraction<MatsT>> contractGLL113 =
          { {contract1PDMSS.S().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 8*NB1C3},
            {contract1PDMSS.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3},
            {contract1PDMSS.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3},
            {contract1PDMSS.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLL113);
    
        // Add to the LL part of 4C exchangeMatrix in Pauli matrix form
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[bf1+i*NB2C] += -3.0*scale*Scr1[i] +3.0*iscale*Scr2[i] +3.0*iscale*Scr3[i] +3.0*iscale*Scr4[i];
        }
    
    
        /* Equation (114) */
        std::vector<TwoBodyContraction<MatsT>> contractGLL114 =
          { {contract1PDMSS.Z().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+18*NB1C3},
            {contract1PDMSS.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3},
            {contract1PDMSS.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+14*NB1C3},
            {contract1PDMSS.Y().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+16*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLL114);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[bf1+i*NB2C] += scale*Scr1[i] +iscale*Scr2[i] +scale*Scr3[i] +scale*Scr4[i];
        }
    
    
    
        /* Equation (115) */
        std::vector<TwoBodyContraction<MatsT>> contractGLL115 =
          { {contract1PDMSS.X().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+19*NB1C3},
            {contract1PDMSS.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3},
            {contract1PDMSS.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+12*NB1C3},
            {contract1PDMSS.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+14*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLL115);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[bf1+i*NB2C] += scale*Scr1[i] +iscale*Scr2[i] +scale*Scr3[i] +scale*Scr4[i];
        }
    
    
    
        /* Equation (116) */
        std::vector<TwoBodyContraction<MatsT>> contractGLL116 =
          { {contract1PDMSS.Y().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+20*NB1C3},
            {contract1PDMSS.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3},
            {contract1PDMSS.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+12*NB1C3},
            {contract1PDMSS.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+16*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLL116);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[bf1+i*NB2C] += scale*Scr1[i] +iscale*Scr2[i] +scale*Scr3[i] +scale*Scr4[i];
        }
    
    
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After Gaunt LLLL Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
#endif //Gaunt LLLL
    
        /*----------------------------------*/
        /* End of Gaunt (LL|LL) Contraction */
        /*----------------------------------*/
    
    
    
    
        /*++++++++++++++++++++++++++++++++++++*/
        /* Start of Gaunt (SS|SS) Contraction */
        /*++++++++++++++++++++++++++++++++++++*/
    
#if 1 //Gaunt SSSS
    
        /* Equation (129) */
        std::vector<TwoBodyContraction<MatsT>> contractGSS129 =
          { {contract1PDMLL.S().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+27*NB1C3},
            {contract1PDMLL.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+28*NB1C3},
            {contract1PDMLL.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+29*NB1C3},
            {contract1PDMLL.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+30*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGSS129);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[SS+bf1+i*NB2C] += -3.0*scale*Scr1[i] - iscale*Scr2[i] - iscale*Scr3[i] - iscale*Scr4[i];
        }
        
    
        /* Equation (130) */
        std::vector<TwoBodyContraction<MatsT>> contractGSS130 =
          { {contract1PDMLL.Z().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+34*NB1C3},
            {contract1PDMLL.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+30*NB1C3},
            {contract1PDMLL.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+33*NB1C3},
            {contract1PDMLL.Y().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+32*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGSS130);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[SS+bf1+i*NB2C] += scale*Scr1[i] - 3.0*iscale*Scr2[i] + scale*Scr3[i] + scale*Scr4[i];
        }
    
    
        /* Equation (131) */
        std::vector<TwoBodyContraction<MatsT>> contractGSS131 =
          { {contract1PDMLL.X().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+35*NB1C3},
            {contract1PDMLL.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+28*NB1C3},
            {contract1PDMLL.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+33*NB1C3},
            {contract1PDMLL.Y().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+31*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGSS131);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[SS+bf1+i*NB2C] += scale*Scr1[i] - 3.0*iscale*Scr2[i] + scale*Scr3[i] + scale*Scr4[i];
        }
    
    
        /* Equation (132) */
        std::vector<TwoBodyContraction<MatsT>> contractGSS132 =
          { {contract1PDMLL.Y().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+36*NB1C3},
            {contract1PDMLL.S().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+29*NB1C3},
            {contract1PDMLL.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+31*NB1C3},
            {contract1PDMLL.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+32*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGSS132);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[SS+bf1+i*NB2C] += scale*Scr1[i] - 3.0*iscale*Scr2[i] + scale*Scr3[i] + scale*Scr4[i];
        }
    
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After Gaunt SSSS Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
    
#endif // Gaunt SSSS
        /*------------------------------------*/
        /*   End of Gaunt (SS|SS) Contraction */
        /*------------------------------------*/
    
    
    
    
    
    
        /*++++++++++++++++++++++++++++++++++++*/
        /* Start of Gaunt (LL|SS) Contraction */
        /*++++++++++++++++++++++++++++++++++++*/
    
    
#if 1 // Gaunt LLSS COULOMB
    
        /* Equation (91) */
        std::vector<TwoBodyContraction<MatsT>> contractGLS91 =
          { {contract1PDMLS.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3},
            {contract1PDMLS.X().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 9*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+10*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr4, HerDen, COULOMB, .ERI4=ERI4bf1+11*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS91);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[LS+bf1+i*NB2C] += 2.0*scale*Scr1[i] -2.0*iscale*Scr2[i] -2.0*iscale*Scr3[i] -2.0*iscale*Scr4[i];
        }
    
    
        /* Equation (92)X first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92AX =
          { {contract1PDMLS.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+9*NB1C3},
            {contract1PDMLS.X().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+8*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92AX);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += -2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
        /* Equation (92)X last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92BX =
          { {contract1PDMLS.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+21*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+13*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+15*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92BX);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
    
        /* Equation (92)Y first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92AY =
          { {contract1PDMLS.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+10*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92AY);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += -2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
    
        /* Equation (92)Y last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92BY =
          { {contract1PDMLS.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+22*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+24*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+17*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92BY);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
    
        /* Equation (92)Z first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92AZ =
          { {contract1PDMLS.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+11*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92AZ);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += -2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
        /* Equation (92)Z last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS92BZ =
          { {contract1PDMLS.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+23*NB1C3},
            {contract1PDMLS.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+25*NB1C3},
            {contract1PDMLS.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+26*NB1C3} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS92BZ);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
    
    
#ifdef _PRINT_MATRICES
        std::cout<<"After Gaunt 91-92 Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif
    
    
        /* Equation (136) */
        std::vector<TwoBodyContraction<MatsT>> contractGLS136 =
          { {contract1PDMSL.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 9*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr4, HerDen, COULOMB, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS136);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*iscale*Scr2[i] -2.0*iscale*Scr3[i] -2.0*iscale*Scr4[i];
        }
    
    
        /* Equation (137)X first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137AX =
          { {contract1PDMSL.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+9*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+8*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137AX);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += 2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
        /* Equation (137)X last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137BX =
          { {contract1PDMSL.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+21*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+13*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+15*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137BX);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
        /* Equation (137)Y first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137AY =
          { {contract1PDMSL.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137AY);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += 2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
        /* Equation (137)Y last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137BY =
          { {contract1PDMSL.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+22*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+24*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+17*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137BY);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
        /* Equation (137)Z first two terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137AZ =
          { {contract1PDMSL.S().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137AZ);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += 2.0*iscale*Scr1[i] +2.0*scale*Scr2[i];
        }
    
    
        /* Equation (137)Z last term */
        std::vector<TwoBodyContraction<MatsT>> contractGLS137BZ =
          { {contract1PDMSL.X().pointer(), Scr1, HerDen, COULOMB, .ERI4=ERI4bf1+23*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, COULOMB, .ERI4=ERI4bf1+25*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, COULOMB, .ERI4=ERI4bf1+26*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS137BZ);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += -2.0*scale*Scr1[i] -2.0*scale*Scr2[i] -2.0*scale*Scr3[i];
        }
    
    
    
    
#ifdef _PRINT_MATRICES
        std::cout<<"After Gaunt 136-137 Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif
    
    
#endif  // Gaunt LLSS COULOMB
    
    
    
    
#if 1 // Gaunt LLSS EXCHANGE
        /* Equation (159) */
        std::vector<TwoBodyContraction<MatsT>> contractGLS159 =
          { {contract1PDMSL.S().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS159);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->S().pointer()[LS+bf1+i*NB2C] += -scale*Scr1[i] +iscale*Scr2[i] +iscale*Scr3[i] +iscale*Scr4[i];
        }
    
        /* Equation (160) first four terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS160A =
          { {contract1PDMSL.Z().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.S().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS160A);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += 2.0*scale*Scr1[i] +2.0*scale*Scr2[i] -2.0*scale*Scr3[i] -iscale*Scr4[i];
        }
    
    
        /* Equation (160) last three terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS160B =
          { {contract1PDMSL.Z().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+18*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+14*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+16*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS160B);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Z().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +scale*Scr2[i] +scale*Scr3[i];
        }
    
    
    
        /* Equation (161) first four terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS161A =
          { {contract1PDMSL.X().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.S().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS161A);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += 2.0*scale*Scr1[i] -2.0*scale*Scr2[i] +2.0*scale*Scr3[i] -iscale*Scr4[i];
        }
    
    
        /* Equation (161) last three terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS161B =
          { {contract1PDMSL.X().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+19*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Y().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+12*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+14*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS161B);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->X().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +scale*Scr2[i] +scale*Scr3[i];
        }
    
    
        /* Equation (162) first four terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS162A =
          { {contract1PDMSL.Y().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 8*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+11*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+ 9*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.S().pointer(), Scr4, HerDen, EXCHANGE, .ERI4=ERI4bf1+10*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS162A);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += 2.0*scale*Scr1[i] +2.0*scale*Scr2[i] -2.0*scale*Scr3[i] -iscale*Scr4[i];
        }
    
    
        /* Equation (162) last three terms */
        std::vector<TwoBodyContraction<MatsT>> contractGLS162B =
          { {contract1PDMSL.Y().pointer(), Scr1, HerDen, EXCHANGE, .ERI4=ERI4bf1+20*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.X().pointer(), Scr2, HerDen, EXCHANGE, .ERI4=ERI4bf1+12*NB1C3, .intTrans=TRANS_KL},
            {contract1PDMSL.Z().pointer(), Scr3, HerDen, EXCHANGE, .ERI4=ERI4bf1+16*NB1C3, .intTrans=TRANS_KL} };
    
        // Call the contraction engine to do the assembly
        relERICon.twoBodyContract3Index(ss.comm, contractGLS162B);
    
        // Assemble 4C exchangeMatrix 
        for(auto i=0; i<NB1C; i++){
          ss.exchangeMatrix->Y().pointer()[LS+bf1+i*NB2C] += scale*Scr1[i] +scale*Scr2[i] +scale*Scr3[i];
        }
    
    
#ifdef _PRINT_MATRICES
        std::cout<<"After Gaunt 159-162 Iteration #"<<bf1<<std::endl;
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif //_PRINT_MATRICES
    
    
#ifdef _PRINT_MATRICES
    
        std::cout<<"After Gaunt LLSS|SSLL"<<std::endl;
        prettyPrintSmart(std::cout, "COULOMB",    ss.coulombMatrix->pointer(),      NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
        prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
    
#endif //_PRINT_MATRICES
    
#endif // Gaunt LLSS EXCHANGE
    
        /*------------------------------------*/
        /*   End of Gaunt (LL|SS) Contraction */
        /*------------------------------------*/
    
      } //GAUNT
  
    } // Loop over bf1 using 3-Index ERI
    } // Loop over s1 using 3-Index ERI



    /*******************************/
    /* Final Assembly of 4C Matrix */
    /*******************************/
    ROOT_ONLY(ss.comm);

    // Copy LS to SL part of the exchangeMatrix[MS]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MX]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MY]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MZ]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+SL, NB2C);

    // Form GD: G[D] = 2.0*J[D] - K[D]
    if( std::abs(xHFX) > 1e-12 ) {
      *ss.twoeH = -xHFX * *ss.exchangeMatrix;
    } else {
      ss.twoeH->clear();
    }
    // G[D] += 2*J[D]
    *ss.twoeH += 2.0 * *ss.coulombMatrix;


    CQMemManager::get().free(Scr1);
    CQMemManager::get().free(Scr2);
    CQMemManager::get().free(Scr3);
    CQMemManager::get().free(Scr4);


#ifdef _PRINT_MATRICES

    prettyPrintSmart(std::cout,"twoeH MS",ss.twoeH->S().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MZ",ss.twoeH->Z().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MY",ss.twoeH->Y().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MX",ss.twoeH->X().pointer(),NB2C,NB2C,NB2C);


    MatsT* TEMP_GATHER1 = CQMemManager::get().malloc<MatsT>(NB4C2);
    MatsT* TEMP_GATHER2 = CQMemManager::get().malloc<MatsT>(NB4C2);

    memset(TEMP_GATHER1,0.,NB4C2*sizeof(MatsT));
    memset(TEMP_GATHER2,0.,NB4C2*sizeof(MatsT));

    std::cout << std::scientific << std::setprecision(16);
    SpinGather(NB2C,TEMP_GATHER1,NB4C,contract1PDM.S().pointer(),NB2C,contract1PDM.Z().pointer(),NB2C,contract1PDM.Y().pointer(),NB2C,contract1PDM.X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"density Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);


    SpinGather(NB2C,TEMP_GATHER2,NB4C,ss.twoeH->S().pointer(),NB2C,ss.twoeH->Z().pointer(),NB2C,ss.twoeH->Y().pointer(),NB2C,ss.twoeH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"twoeH Gather",TEMP_GATHER2,NB4C,NB4C,NB4C,1,12,16);

    SpinGather(NB2C,TEMP_GATHER1,NB4C,ss.coreH->S().pointer(),NB2C,ss.coreH->Z().pointer(),NB2C,ss.coreH->Y().pointer(),NB2C,ss.coreH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"coreH Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);
 
    CQMemManager::get().free(TEMP_GATHER1);
    CQMemManager::get().free(TEMP_GATHER2);

#endif //_PRINT_MATRICES


  }; // FourCompFock<MatsT, IntsT>::formGD3Index
#endif


  /**   
   *  \brief Forms the 4C Fock matrix using AO-direct
   */
  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formGDDirect(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {

    GTODirectRelERIContraction<MatsT,IntsT> &relERICon =
        *std::dynamic_pointer_cast<GTODirectRelERIContraction<MatsT,IntsT>>(ss.TPI);

    // Decide list of onePDMs to use
    cqmatrix::PauliSpinorMatrices<MatsT> &contract1PDM
        = increment ? *ss.deltaOnePDM : *ss.onePDM;

    size_t NB1C  = ss.basisSet().nBasis;
    size_t NB2C  = 2 * ss.basisSet().nBasis;
    size_t NB4C  = 4 * ss.basisSet().nBasis;
    size_t NB1C2 = NB1C*NB1C;
    size_t NB1C4 = NB1C*NB1C*NB1C*NB1C;
    size_t NB1C3 = NB1C*NB1C*NB1C;
    size_t NB2C2 = NB2C*NB2C;
    size_t NB4C2 = NB4C*NB4C;

    size_t SS = NB2C*NB1C+NB1C;
    size_t LS = NB2C*NB1C;
    size_t SL = NB1C;

    auto MS = SCALAR;
    bool computeExchange = true;
    bool do_VLL = false;

    size_t mpiRank   = MPIRank(ss.comm);
    bool   isNotRoot = mpiRank != 0;

    double xHFX_VLL = 1.; 

    double xGaugeX = this->hamiltonianOptions_.GaugeScale;
    double xGauntX = this->hamiltonianOptions_.GauntScale;

    // Check if doing DFT
    // DKS w/ VXCLL requires computeExchange for all functionals (even pure)
    if (typeid(ss) == typeid(KohnSham<MatsT,IntsT>)){
    if(this->hamiltonianOptions_.dksType == DKS_TYPE::VLL){
      computeExchange = true;
      xHFX_VLL = xHFX;
      xHFX = 1.;
      do_VLL = true;
    }}
    else{
      computeExchange = true;
      xHFX_VLL = 1.;
    }

    cqmatrix::PauliSpinorMatrices<MatsT> exchangeMatrixLL(NB1C);

    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLL(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMLS(NB1C);
    cqmatrix::PauliSpinorMatrices<MatsT> contract1PDMSL(NB1C);

    MatsT* CScrLLMS = CQMemManager::get().malloc<MatsT>(NB1C2);

    MatsT* CScrSSMS = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrSSMX = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrSSMY = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrSSMZ = CQMemManager::get().malloc<MatsT>(NB1C2);

    MatsT* CScrLSMS = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrLSMX = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrLSMY = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* CScrLSMZ = CQMemManager::get().malloc<MatsT>(NB1C2);

    MatsT* XScrLLMS = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLLMX = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLLMY = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLLMZ = CQMemManager::get().malloc<MatsT>(NB1C2);

    MatsT* XScrSSMS = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrSSMX = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrSSMY = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrSSMZ = CQMemManager::get().malloc<MatsT>(NB1C2);

    MatsT* XScrLSMS = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLSMX = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLSMY = CQMemManager::get().malloc<MatsT>(NB1C2);
    MatsT* XScrLSMZ = CQMemManager::get().malloc<MatsT>(NB1C2);


    memset(CScrLLMS,0.,NB1C2*sizeof(MatsT));

    memset(CScrSSMS,0.,NB1C2*sizeof(MatsT));
    memset(CScrSSMX,0.,NB1C2*sizeof(MatsT));
    memset(CScrSSMY,0.,NB1C2*sizeof(MatsT));
    memset(CScrSSMZ,0.,NB1C2*sizeof(MatsT));

    memset(CScrLSMS,0.,NB1C2*sizeof(MatsT));
    memset(CScrLSMX,0.,NB1C2*sizeof(MatsT));
    memset(CScrLSMY,0.,NB1C2*sizeof(MatsT));
    memset(CScrLSMZ,0.,NB1C2*sizeof(MatsT));

    memset(XScrLLMS,0.,NB1C2*sizeof(MatsT));
    memset(XScrLLMX,0.,NB1C2*sizeof(MatsT));
    memset(XScrLLMY,0.,NB1C2*sizeof(MatsT));
    memset(XScrLLMZ,0.,NB1C2*sizeof(MatsT));

    memset(XScrSSMS,0.,NB1C2*sizeof(MatsT));
    memset(XScrSSMX,0.,NB1C2*sizeof(MatsT));
    memset(XScrSSMY,0.,NB1C2*sizeof(MatsT));
    memset(XScrSSMZ,0.,NB1C2*sizeof(MatsT));

    memset(XScrLSMS,0.,NB1C2*sizeof(MatsT));
    memset(XScrLSMX,0.,NB1C2*sizeof(MatsT));
    memset(XScrLSMY,0.,NB1C2*sizeof(MatsT));
    memset(XScrLSMZ,0.,NB1C2*sizeof(MatsT));


    // Compute 1/(2mc)^2
    auto C2 = 1./(4*SpeedOfLight()*SpeedOfLight());

    for(size_t i = 0; i < contract1PDM.nComponent(); i++) {
      cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer(),    NB2C,
             contract1PDMLL[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SS, NB2C,
             contract1PDMSS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+LS, NB2C,
             contract1PDMLS[c].pointer(), NB1C);
      SetMat('N', NB1C, NB1C, MatsT(1.), contract1PDM[c].pointer()+SL, NB2C,
             contract1PDMSL[c].pointer(), NB1C);
    }

#ifdef _PRINT_MATRICES
    prettyPrintSmart(std::cout, "1PDM[MS]", contract1PDM.S().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MX]", contract1PDM.X().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MY]", contract1PDM.Y().pointer(), NB2C, NB2C, NB2C);
    prettyPrintSmart(std::cout, "1PDM[MZ]", contract1PDM.Z().pointer(), NB2C, NB2C, NB2C);
#endif

#if 0
    std::fill_n(contract1PDMLL.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLL.Z().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMSS.Z().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.S().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.X().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.Y().pointer(),NB1C2,1.0);
    std::fill_n(contract1PDMLS.Z().pointer(),NB1C2,1.0);
#endif

    ss.twoeH->clear();
    if(not increment) {
      ss.coulombMatrix->clear();
      ss.exchangeMatrix->clear();
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
#ifdef _FOCK_CONTRACTION_TIME
      auto topBC_FLL = tick();
#endif

      if(this->hamiltonianOptions_.Libcint) {

        std::vector<TwoBodyContraction<MatsT>> contractLL =
          { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, BARE_COULOMB},
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
  
        // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
        relERICon.twoBodyContract(ss.comm, true, contractLL, pert, computeExchange);
  
        /* Store LL block into 2C spin scattered matrices */
        SetMat('N', NB1C, NB1C, MatsT(2.), CScrLLMS, NB1C, ss.twoeH->S().pointer(), NB2C);
        // Scaled by percent HF exchange
        if (computeExchange) {
          SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL*xHFX, XScrLLMS, NB1C, ss.exchangeMatrix->S().pointer(), NB2C);
          SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL*xHFX, XScrLLMX, NB1C, ss.exchangeMatrix->X().pointer(), NB2C);
          SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL*xHFX, XScrLLMY, NB1C, ss.exchangeMatrix->Y().pointer(), NB2C);
          SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL*xHFX, XScrLLMZ, NB1C, ss.exchangeMatrix->Z().pointer(), NB2C);
        }

      } else {
        // using Libint  
        std::vector<TwoBodyContraction<MatsT>> contractLL =
          { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, COULOMB} };
    
        // Determine how many (if any) exchange terms to calculate
        if(computeExchange) {
          exchangeMatrixLL.clear();
          for(size_t i = 0; i < ss.exchangeMatrix->nComponent(); i++) {
    
            cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
            contractLL.push_back(
              {contract1PDMLL[c].pointer(), exchangeMatrixLL[c].pointer(), HerDen, EXCHANGE}
            );
          }
        }
    
        // Zero out K[i]
        if(not increment) ss.exchangeMatrix->clear();
    
        // Call the contraction engine to do the assembly of direct Coulomb LLLL
        GTODirectTPIContraction<MatsT,IntsT>(ss.TPI->ints()).twoBodyContract(ss.comm, true, contractLL, pert);
    
        /* Store LL block into 2C spin scattered matrices */
        // Assemble 4C coulombMatrix
        SetMat('N', NB1C, NB1C, MatsT(1.), CScrLLMS, NB1C, ss.twoeH->S().pointer(), NB2C);
    
        // Assemble 4C exchangeMatrix 
        if(computeExchange) {
          for(auto i = 0; i < ss.exchangeMatrix->nComponent();i++){
            cqmatrix::PAULI_SPINOR_COMPS c = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(i);
            SetMat('N', NB1C, NB1C, MatsT(1.)*xHFX_VLL*xHFX, exchangeMatrixLL[c].pointer(), NB1C,
                   (*ss.exchangeMatrix)[c].pointer(), NB2C);
          }
        }
      }

 #ifdef _PRINT_MATRICES
      std::cout<<"After BARE COULOMB"<<std::endl;
      prettyPrintSmart(std::cout, "COULOMB-S",           ss.twoeH->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-X",           ss.twoeH->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Y",           ss.twoeH->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "COULOMB-Z",           ss.twoeH->Z().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-S", ss.exchangeMatrix->S().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-X", ss.exchangeMatrix->X().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Y", ss.exchangeMatrix->Y().pointer(), NB2C, NB2C, NB2C);
      prettyPrintSmart(std::cout, "EXCHANGE-Z", ss.exchangeMatrix->Z().pointer(), NB2C, NB2C, NB2C);
#endif

#ifdef _FOCK_CONTRACTION_TIME
      // Print out BareCoulomb duration 
      auto durBC_FLL_dir = tock(topBC_FLL);
      std::cout << "Bare Coulomb Contribution to F_LL = " << durBC_FLL_dir << std::endl;
#endif 

      /*---------------------------------------------*/
      /*   End of Direct Coulomb (LL|LL) Contraction */
      /*---------------------------------------------*/

    } // DIRECT_COULOMB

    //if (not HerDen) {
    //ss.twoeH->clear();
    //if(not increment) {
    //  ss.coulombMatrix->clear();
    //  ss.exchangeMatrix->clear();
    //};
    //}
    /**********************************************/
    /*                                            */
    /*              DIRAC-COULOMB                 */
    /*                                            */
    /**********************************************/

    if(this->hamiltonianOptions_.DiracCoulomb) { // DIRAC_COULOMB

  
      /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (C^-2 order) Contraction */
      /*+++++++++++++++++++++++++++++++++++++++++++++++++*/
  
      std::vector<TwoBodyContraction<MatsT>> contractDCLL =
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

      // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
      relERICon.twoBodyContract(ss.comm, true, contractDCLL, pert, 
        computeExchange,
        this->hamiltonianOptions_.DiracCoulombType,
        this->hamiltonianOptions_.DiracCoulombApproximationType);

      // Add Dirac-Coulomb contributions to the LLLL block
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLLMS, NB1C, MatsT(1.0), 
                      ss.twoeH->S().pointer(), NB2C,
                      ss.twoeH->S().pointer(), NB2C);

      // Add Dirac-Coulomb contributions to the SSSS block
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrSSMS, NB1C, MatsT(1.0), 
                      ss.twoeH->S().pointer()+SS, NB2C,
                      ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrSSMX, NB1C, MatsT(1.0), 
                      ss.twoeH->X().pointer()+SS, NB2C,
                      ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrSSMY, NB1C, MatsT(1.0), 
                      ss.twoeH->Y().pointer()+SS, NB2C,
                      ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrSSMZ, NB1C, MatsT(1.0), 
                      ss.twoeH->Z().pointer()+SS, NB2C,
                      ss.twoeH->Z().pointer()+SS, NB2C);

      // Add Dirac-Coulomb contributions to the LLSS block
      MatAdd('N','N', NB1C, NB1C, -C2*xHFX, XScrLSMS, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->S().pointer()+LS, NB2C,
		      ss.exchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2*xHFX, XScrLSMX, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->X().pointer()+LS, NB2C,
		      ss.exchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2*xHFX, XScrLSMY, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->Y().pointer()+LS, NB2C,
		      ss.exchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2*xHFX, XScrLSMZ, NB1C, MatsT(1.0), 
		      ss.exchangeMatrix->Z().pointer()+LS, NB2C,
		      ss.exchangeMatrix->Z().pointer()+LS, NB2C);

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

#endif


#if 0 // LLSS contribution is included in the Dirac-Coulomb term
      if(computeExchange) {
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
#endif

    } //_DIRAC_COULOMB



    /*************************************/
    /*                                   */
    /*              SSSS                 */
    /*                                   */
    /*************************************/

    if(this->hamiltonianOptions_.DiracCoulombSSSS) { // SSSS

#ifdef _FOCK_CONTRACTION_TIME
      auto topDCSS_FSS = tick();
#endif
      double C4 = 1./(16*SpeedOfLight()*SpeedOfLight()*SpeedOfLight()*SpeedOfLight());
  
      /*++++++++++++++++++++++++++++++++++++++++++++*/
      /* Start of Dirac-Coulomb (SS|SS) Contraction */
      /*++++++++++++++++++++++++++++++++++++++++++++*/

      std::vector<TwoBodyContraction<MatsT>> contractDCSS =
        { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, SSSS},
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

      // Call the contraction engine to do the assembly of Dirac-Coulomb LLLL
      relERICon.twoBodyContract(ss.comm, true, contractDCSS, pert, 
        computeExchange,
        this->hamiltonianOptions_.SSSSType,
        this->hamiltonianOptions_.SSSSApproximationType);

      // Add (SS|SS) Coulomb contributions to the SSSS block
      MatAdd('N','N', NB1C, NB1C, 2.0*C4, CScrSSMS, NB1C, MatsT(1.0), 
                      ss.twoeH->S().pointer()+SS, NB2C,
                      ss.twoeH->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C4, CScrSSMX, NB1C, MatsT(1.0), 
                      ss.twoeH->X().pointer()+SS, NB2C,
                      ss.twoeH->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C4, CScrSSMY, NB1C, MatsT(1.0), 
                      ss.twoeH->Y().pointer()+SS, NB2C,
                      ss.twoeH->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C4, CScrSSMZ, NB1C, MatsT(1.0), 
                      ss.twoeH->Z().pointer()+SS, NB2C,
                      ss.twoeH->Z().pointer()+SS, NB2C);

      if (computeExchange) {
      // Add (SS|SS) exchange contributions to the SSSS block
      MatAdd('N','N', NB1C, NB1C, -C4*xHFX, XScrSSMS, NB1C, MatsT(1.0), 
                      ss.exchangeMatrix->S().pointer()+SS, NB2C,
                      ss.exchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C4*xHFX, XScrSSMX, NB1C, MatsT(1.0), 
                      ss.exchangeMatrix->X().pointer()+SS, NB2C,
                      ss.exchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C4*xHFX, XScrSSMY, NB1C, MatsT(1.0), 
                      ss.exchangeMatrix->Y().pointer()+SS, NB2C,
                      ss.exchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C4*xHFX, XScrSSMZ, NB1C, MatsT(1.0), 
                      ss.exchangeMatrix->Z().pointer()+SS, NB2C,
                      ss.exchangeMatrix->Z().pointer()+SS, NB2C);
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


#ifdef _FOCK_CONTRACTION_TIME
      // Print out BareCoulomb duration 
      auto durDCSS_FSS_dir = tock(topDCSS_FSS);
      std::cout << "Dirac Coulomb SSSS Contribution to F_SS = " << durDCSS_FSS_dir << std::endl;
#endif 
    }



    /*************************************/
    /*                                   */
    /*              GAUNT                */
    /*                                   */
    /*************************************/

    // if the gauge term is included, the Gaunt term needs to be scaled by half
    if(this->hamiltonianOptions_.Gauge) C2=C2/2.0;

    if(this->hamiltonianOptions_.Gaunt) { // Gaunt

      std::vector<TwoBodyContraction<MatsT>> contractDCGaunt =
        { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, GAUNT},
  	  //
          {contract1PDMLL.S().pointer(), XScrLLMS},
          {contract1PDMLL.X().pointer(), XScrLLMX},
          {contract1PDMLL.Y().pointer(), XScrLLMY},
          {contract1PDMLL.Z().pointer(), XScrLLMZ},
	  //
          {contract1PDMSS.S().pointer(), CScrSSMS},
          {contract1PDMSS.X().pointer(), CScrSSMX},
          {contract1PDMSS.Y().pointer(), CScrSSMY},
          {contract1PDMSS.Z().pointer(), CScrSSMZ},
	  //
          {contract1PDMSS.S().pointer(), XScrSSMS},
          {contract1PDMSS.X().pointer(), XScrSSMX},
          {contract1PDMSS.Y().pointer(), XScrSSMY},
          {contract1PDMSS.Z().pointer(), XScrSSMZ},
	  //
          {contract1PDMLS.S().pointer(), XScrLSMS},
          {contract1PDMLS.X().pointer(), XScrLSMX},
          {contract1PDMLS.Y().pointer(), XScrLSMY},
          {contract1PDMLS.Z().pointer(), XScrLSMZ},
	  //
          {contract1PDMLS.S().pointer(), CScrLSMS},
          {contract1PDMLS.X().pointer(), CScrLSMX},
          {contract1PDMLS.Y().pointer(), CScrLSMY},
          {contract1PDMLS.Z().pointer(), CScrLSMZ},
	  //
          {contract1PDMSL.S().pointer(), CScrLSMS},
          {contract1PDMSL.X().pointer(), CScrLSMX},
          {contract1PDMSL.Y().pointer(), CScrLSMY},
          {contract1PDMSL.Z().pointer(), CScrLSMZ},
	  //
          {contract1PDMSL.S().pointer(), XScrLSMS},
          {contract1PDMSL.X().pointer(), XScrLSMX},
          {contract1PDMSL.Y().pointer(), XScrLSMY},
          {contract1PDMSL.Z().pointer(), XScrLSMZ},
	};


    if (this->hamiltonianOptions_.updateGaunt == true){

        relERICon.twoBodyContract(ss.comm, true, contractDCGaunt, pert, 
        computeExchange,
        this->hamiltonianOptions_.GauntType,
        this->hamiltonianOptions_.GauntApproximationType);

        ss.gaunttwoeH->clear();
        ss.gauntexchangeMatrix->clear();

      // Add (LL|SS) Coulomb contributions
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMS, NB1C, MatsT(1.0), 
                      ss.gaunttwoeH->S().pointer()+LS, NB2C,
                      ss.gaunttwoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMX, NB1C, MatsT(1.0), 
                      ss.gaunttwoeH->X().pointer()+LS, NB2C,
                      ss.gaunttwoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMY, NB1C, MatsT(1.0), 
                      ss.gaunttwoeH->Y().pointer()+LS, NB2C,
                      ss.gaunttwoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMZ, NB1C, MatsT(1.0), 
                      ss.gaunttwoeH->Z().pointer()+LS, NB2C,
                      ss.gaunttwoeH->Z().pointer()+LS, NB2C);

      if (computeExchange) {
      // Add (LL|LL) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMS, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->S().pointer(), NB2C,
                      ss.gauntexchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMX, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->X().pointer(), NB2C,
                      ss.gauntexchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMY, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Y().pointer(), NB2C,
                      ss.gauntexchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMZ, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Z().pointer(), NB2C,
                      ss.gauntexchangeMatrix->Z().pointer(), NB2C);

      // Add (SS|SS) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMS, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->S().pointer()+SS, NB2C,
                      ss.gauntexchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMX, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->X().pointer()+SS, NB2C,
                      ss.gauntexchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMY, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Y().pointer()+SS, NB2C,
                      ss.gauntexchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMZ, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Z().pointer()+SS, NB2C,
                      ss.gauntexchangeMatrix->Z().pointer()+SS, NB2C);

      // Add (LL|SS) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMS, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->S().pointer()+LS, NB2C,
                      ss.gauntexchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMX, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->X().pointer()+LS, NB2C,
                      ss.gauntexchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMY, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Y().pointer()+LS, NB2C,
                      ss.gauntexchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMZ, NB1C, MatsT(1.0), 
                      ss.gauntexchangeMatrix->Z().pointer()+LS, NB2C,
                      ss.gauntexchangeMatrix->Z().pointer()+LS, NB2C);
      }
    } 

      // Add gaunt contribution 
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaunttwoeH->S().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->S().pointer(), NB2C,
                    ss.twoeH->S().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaunttwoeH->X().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->X().pointer(), NB2C,
                    ss.twoeH->X().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaunttwoeH->Y().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->Y().pointer(), NB2C,
                    ss.twoeH->Y().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaunttwoeH->Z().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->Z().pointer(), NB2C,
                    ss.twoeH->Z().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGauntX, ss.gauntexchangeMatrix->S().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->S().pointer(), NB2C,
                    ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGauntX, ss.gauntexchangeMatrix->X().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->X().pointer(), NB2C,
                    ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGauntX, ss.gauntexchangeMatrix->Y().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->Y().pointer(), NB2C,
                    ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGauntX, ss.gauntexchangeMatrix->Z().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->Z().pointer(), NB2C,
                    ss.exchangeMatrix->Z().pointer(), NB2C);





#ifdef _PRINT_MATRICES
      std::cout<<"After GAUNT"<<std::endl;
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
    /*              GAUGE                */
    /*                                   */
    /*************************************/

    if(this->hamiltonianOptions_.Gauge) { // Gauge

      std::vector<TwoBodyContraction<MatsT>> contractDCGauge =
        { {contract1PDMLL.S().pointer(), CScrLLMS, HerDen, GAUGE},
          //
          {contract1PDMLL.S().pointer(), XScrLLMS},
          {contract1PDMLL.X().pointer(), XScrLLMX},
          {contract1PDMLL.Y().pointer(), XScrLLMY},
          {contract1PDMLL.Z().pointer(), XScrLLMZ},
          //
          {contract1PDMSS.S().pointer(), CScrSSMS},
          {contract1PDMSS.X().pointer(), CScrSSMX},
          {contract1PDMSS.Y().pointer(), CScrSSMY},
          {contract1PDMSS.Z().pointer(), CScrSSMZ},
          //
          {contract1PDMSS.S().pointer(), XScrSSMS},
          {contract1PDMSS.X().pointer(), XScrSSMX},
          {contract1PDMSS.Y().pointer(), XScrSSMY},
          {contract1PDMSS.Z().pointer(), XScrSSMZ},
          //
          {contract1PDMLS.S().pointer(), XScrLSMS},
          {contract1PDMLS.X().pointer(), XScrLSMX},
          {contract1PDMLS.Y().pointer(), XScrLSMY},
          {contract1PDMLS.Z().pointer(), XScrLSMZ},
          //
          {contract1PDMLS.S().pointer(), CScrLSMS},
          {contract1PDMLS.X().pointer(), CScrLSMX},
          {contract1PDMLS.Y().pointer(), CScrLSMY},
          {contract1PDMLS.Z().pointer(), CScrLSMZ},
          //
          {contract1PDMSL.S().pointer(), CScrLSMS},
          {contract1PDMSL.X().pointer(), CScrLSMX},
          {contract1PDMSL.Y().pointer(), CScrLSMY},
          {contract1PDMSL.Z().pointer(), CScrLSMZ},
          //
          {contract1PDMSL.S().pointer(), XScrLSMS},
          {contract1PDMSL.X().pointer(), XScrLSMX},
          {contract1PDMSL.Y().pointer(), XScrLSMY},
          {contract1PDMSL.Z().pointer(), XScrLSMZ},
        };

      // Call the contraction engine to do the assembly of Gaunt

    if (this->hamiltonianOptions_.updateGauge == true){
        relERICon.twoBodyContract(ss.comm, true, contractDCGauge, pert,
        computeExchange,
        this->hamiltonianOptions_.GaugeType,
        this->hamiltonianOptions_.GaugeApproximationType);

        ss.gaugetwoeH->clear();
        ss.gaugeexchangeMatrix->clear();

      // Add (LL|SS) Coulomb contributions
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMS, NB1C, MatsT(1.0), 
                      ss.gaugetwoeH->S().pointer()+LS, NB2C,
                      ss.gaugetwoeH->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMX, NB1C, MatsT(1.0), 
                      ss.gaugetwoeH->X().pointer()+LS, NB2C,
                      ss.gaugetwoeH->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMY, NB1C, MatsT(1.0), 
                      ss.gaugetwoeH->Y().pointer()+LS, NB2C,
                      ss.gaugetwoeH->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, 2.0*C2, CScrLSMZ, NB1C, MatsT(1.0), 
                      ss.gaugetwoeH->Z().pointer()+LS, NB2C,
                      ss.gaugetwoeH->Z().pointer()+LS, NB2C);

      if (computeExchange) {
      // Add (LL|LL) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMS, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->S().pointer(), NB2C,
                      ss.gaugeexchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMX, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->X().pointer(), NB2C,
                      ss.gaugeexchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMY, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Y().pointer(), NB2C,
                      ss.gaugeexchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLLMZ, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Z().pointer(), NB2C,
                      ss.gaugeexchangeMatrix->Z().pointer(), NB2C);


      // Add (SS|SS) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMS, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->S().pointer()+SS, NB2C,
                      ss.gaugeexchangeMatrix->S().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMX, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->X().pointer()+SS, NB2C,
                      ss.gaugeexchangeMatrix->X().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMY, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Y().pointer()+SS, NB2C,
                      ss.gaugeexchangeMatrix->Y().pointer()+SS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrSSMZ, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Z().pointer()+SS, NB2C,
                      ss.gaugeexchangeMatrix->Z().pointer()+SS, NB2C);

      // Add (LL|SS) exchange contributions
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMS, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->S().pointer()+LS, NB2C,
                      ss.gaugeexchangeMatrix->S().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMX, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->X().pointer()+LS, NB2C,
                      ss.gaugeexchangeMatrix->X().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMY, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Y().pointer()+LS, NB2C,
                      ss.gaugeexchangeMatrix->Y().pointer()+LS, NB2C);
      MatAdd('N','N', NB1C, NB1C, -C2, XScrLSMZ, NB1C, MatsT(1.0), 
                      ss.gaugeexchangeMatrix->Z().pointer()+LS, NB2C,
                      ss.gaugeexchangeMatrix->Z().pointer()+LS, NB2C);
      }
    } 

      // Add gauge contribution 
      // Add (LL|SS) Coulomb contributions
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaugetwoeH->S().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->S().pointer(), NB2C,
                    ss.twoeH->S().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaugetwoeH->X().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->X().pointer(), NB2C,
                    ss.twoeH->X().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaugetwoeH->Y().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->Y().pointer(), NB2C,
                    ss.twoeH->Y().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0), ss.gaugetwoeH->Z().pointer(), NB2C, MatsT(1.0), 
                    ss.twoeH->Z().pointer(), NB2C,
                    ss.twoeH->Z().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGaugeX, ss.gaugeexchangeMatrix->S().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->S().pointer(), NB2C,
                    ss.exchangeMatrix->S().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGaugeX, ss.gaugeexchangeMatrix->X().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->X().pointer(), NB2C,
                    ss.exchangeMatrix->X().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGaugeX, ss.gaugeexchangeMatrix->Y().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->Y().pointer(), NB2C,
                    ss.exchangeMatrix->Y().pointer(), NB2C);
      MatAdd('N','N', NB2C, NB2C, MatsT(1.0)*xGaugeX, ss.gaugeexchangeMatrix->Z().pointer(), NB2C, MatsT(1.0), 
                    ss.exchangeMatrix->Z().pointer(), NB2C,
                    ss.exchangeMatrix->Z().pointer(), NB2C);
      

#ifdef _PRINT_MATRICES
      std::cout<<"After GAUGE"<<std::endl;
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



    /*******************************/
    /* Final Assembly of 4C Matrix */
    /*******************************/
    ROOT_ONLY(ss.comm);

    if (computeExchange) {
    // Copy LS to SL part of the exchangeMatrix[MS]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->S().pointer()+LS, NB2C, ss.exchangeMatrix->S().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MX]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->X().pointer()+LS, NB2C, ss.exchangeMatrix->X().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MY]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Y().pointer()+LS, NB2C, ss.exchangeMatrix->Y().pointer()+SL, NB2C);
    // Copy LS to SL part of the exchangeMatrix[MZ]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.exchangeMatrix->Z().pointer()+LS, NB2C, ss.exchangeMatrix->Z().pointer()+SL, NB2C);
    }
    
    // Copy LS to SL part of the twoeH[MS]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->S().pointer()+LS, NB2C, ss.twoeH->S().pointer()+SL, NB2C);
    // Copy LS to SL part of the twoeH[MX]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->X().pointer()+LS, NB2C, ss.twoeH->X().pointer()+SL, NB2C);
    // Copy LS to SL part of the twoeH[MY]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->Y().pointer()+LS, NB2C, ss.twoeH->Y().pointer()+SL, NB2C);
    // Copy LS to SL part of the twoeH[MZ]
    SetMat('C', NB1C, NB1C, MatsT(1.0), ss.twoeH->Z().pointer()+LS, NB2C, ss.twoeH->Z().pointer()+SL, NB2C);
    
    if( computeExchange ){
      *ss.twoeH -= *ss.exchangeMatrix;
    } 
     
    CQMemManager::get().free(CScrLLMS);
    CQMemManager::get().free(CScrSSMS);
    CQMemManager::get().free(CScrSSMX);
    CQMemManager::get().free(CScrSSMY);
    CQMemManager::get().free(CScrSSMZ);
    CQMemManager::get().free(CScrLSMS);
    CQMemManager::get().free(CScrLSMX);
    CQMemManager::get().free(CScrLSMY);
    CQMemManager::get().free(CScrLSMZ);

    CQMemManager::get().free(XScrLLMS);
    CQMemManager::get().free(XScrLLMX);
    CQMemManager::get().free(XScrLLMY);
    CQMemManager::get().free(XScrLLMZ);
    CQMemManager::get().free(XScrSSMS);
    CQMemManager::get().free(XScrSSMX);
    CQMemManager::get().free(XScrSSMY);
    CQMemManager::get().free(XScrSSMZ);
    CQMemManager::get().free(XScrLSMS);
    CQMemManager::get().free(XScrLSMX);
    CQMemManager::get().free(XScrLSMY);
    CQMemManager::get().free(XScrLSMZ);


#ifdef _PRINT_MATRICES

    prettyPrintSmart(std::cout,"twoeH MS",ss.twoeH->S().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MX",ss.twoeH->X().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MY",ss.twoeH->Y().pointer(),NB2C,NB2C,NB2C);
    prettyPrintSmart(std::cout,"twoeH MZ",ss.twoeH->Z().pointer(),NB2C,NB2C,NB2C);


    MatsT* TEMP_GATHER1 = CQMemManager::get().malloc<MatsT>(NB4C2);
    MatsT* TEMP_GATHER2 = CQMemManager::get().malloc<MatsT>(NB4C2);

    memset(TEMP_GATHER1,0.,NB4C2*sizeof(MatsT));
    memset(TEMP_GATHER2,0.,NB4C2*sizeof(MatsT));

    std::cout << std::scientific << std::setprecision(16);
    SpinGather(NB2C,TEMP_GATHER1,NB4C,contract1PDM.S().pointer(),NB2C,contract1PDM.Z().pointer(),NB2C,contract1PDM.Y().pointer(),NB2C,contract1PDM.X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"density Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);


    SpinGather(NB2C,TEMP_GATHER2,NB4C,ss.twoeH->S().pointer(),NB2C,ss.twoeH->Z().pointer(),NB2C,ss.twoeH->Y().pointer(),NB2C,ss.twoeH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"twoeH Gather",TEMP_GATHER2,NB4C,NB4C,NB4C,1,12,16);

    SpinGather(NB2C,TEMP_GATHER1,NB4C,ss.coreH->S().pointer(),NB2C,ss.coreH->Z().pointer(),NB2C,ss.coreH->Y().pointer(),NB2C,ss.coreH->X().pointer(),NB2C);
    prettyPrintSmart(std::cout,"coreH Gather",TEMP_GATHER1,NB4C,NB4C,NB4C,1,12,16);
 
    CQMemManager::get().free(TEMP_GATHER1);
    CQMemManager::get().free(TEMP_GATHER2);

#endif //_PRINT_MATRICES


  }; // FourCompFock<MatsT, IntsT>::formGD3Direct


  template <typename MatsT, typename IntsT>
  void FourCompFock<MatsT,IntsT>::formFock(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX) {

    // General fock build
    FockBuilder<MatsT,IntsT>::formFock(ss, pert, increment, xHFX);


  } // ROFock<MatsT,IntsT>::formFock


}; // namespace ChronusQ

