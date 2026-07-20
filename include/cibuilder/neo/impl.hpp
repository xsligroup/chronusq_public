/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <mcwavefunction.hpp>
#include <cibuilder/neo.hpp>
#include <cibuilder.hpp>
#include <cibuilder/casci/helper.hpp>

#define NEOCASCI_LOOP_INIT() \
  auto neowfn = dynamic_cast<NEOMCWaveFunction<MatsT,IntsT>*>(&mcwfn);\
  auto ewfn = dynamic_cast<MCWaveFunction<MatsT,IntsT>*>(neowfn->ewfn_.get());\
  auto pwfn = dynamic_cast<MCWaveFunction<MatsT,IntsT>*>(neowfn->pwfn_.get());\
  size_t nC = mcwfn.reference().nC; \
  size_t NDet = mcwfn.NDet; \
  std::shared_ptr<const ExcitationList> exList_a = \
    std::dynamic_pointer_cast<CASStringManager>( \
      neowfn->ewfn_->detStr)->excitationList(); \
  std::shared_ptr<const ExcitationList> exList_b = \
    (nC == 1) ?  std::dynamic_pointer_cast<CASStringManager>( \
      neowfn->ewfn_->detStrBeta)->excitationList() : nullptr; \
  std::shared_ptr<const ExcitationList> exList_p = \
    (nC == 1) ?  std::dynamic_pointer_cast<CASStringManager>( \
      neowfn->pwfn_->detStr)->excitationList() : nullptr; \
  size_t nStr_a  = exList_a->nString(); \
  size_t nStr_b = (exList_b) ? exList_b->nString(): 1; \
  size_t nStr_p = (exList_p) ? exList_p->nString(): 1; \
  size_t nNZa = exList_a->nNonZero(); \
  size_t nNZb = (exList_b) ? exList_b->nNonZero(): 0; \
  size_t nNZp = (exList_p) ? exList_p->nNonZero(): 0; 

namespace ChronusQ {

    template<typename MatsT, typename IntsT>
    void NEOCASCI<MatsT,IntsT>::buildFullH(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT* fullH)
    {
        // Cast so we can find the appropriate underlying SS objects
        NEOCASCI_LOOP_INIT();

        auto & ehCore = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("hCoreP_Correlated_Space"));
        auto & PhCore = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("PhCoreP_Correlated_Space"));
        auto & eRI    = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));
        auto & PRI    = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("PRI_Correlated_Space"));
        auto & ePRI   = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("eP_Correlated_Space"));

        // Calculate offset sizes we'll need in constructing the full matrix
        size_t nStr_a2 = nStr_a*nStr_a;
        size_t nStr_b2 = nStr_b*nStr_b;
        size_t nStr_p2 = nStr_p*nStr_p;
        size_t nStr_ab = nStr_a*nStr_b;
        size_t nStr_bp = nStr_b*nStr_p;
        size_t nStr_ap = nStr_a*nStr_p;

        // Allocate needed memory
        size_t nSCR = std::max(std::max(nStr_a,nStr_b),nStr_p);
        size_t nThreads = GetNumThreads();

        // Build a temporary matrix which is alp
        MatsT * abSCR = CQMemManager::get().template malloc<MatsT>(nStr_a2*nStr_b2);
        std::fill_n(abSCR,nStr_a2*nStr_b2,MatsT(0.0));

        CASCI<MatsT,IntsT>::buildFullH(*ewfn,abSCR);
        
        // Temporary storage to avoid reaccessing integrals
        MatsT val;

        // Temporary storage for our blockwise calculation of 
        MatsT * tmpH = CQMemManager::get().template malloc<MatsT>(nSCR*nSCR);

        // Zero out an alpha x alpha block of the full Hamiltonian
        // Zero out the full Hamiltonian
        std::fill_n(fullH,NDet*NDet,MatsT(0.0));

        MatsT * Col = nullptr, *SCR_ith = nullptr;

        // Fill in the first block of the full Hamiltonian
        size_t lastDimOff = nStr_a2*nStr_b2*nStr_p;
        Col = fullH;
        for(size_t ab = 0; ab < nStr_ab; ab++)
        {
            std::copy(abSCR+ab*nStr_a*nStr_b,abSCR+(ab+1)*nStr_a*nStr_b,Col+ab*NDet);
        }
        // Duplicate this block per Proton Block
        lastDimOff = nStr_a2*nStr_b2*nStr_p;
        Col = fullH+lastDimOff;
        for(size_t p = 1; p < nStr_p; p++, Col+=lastDimOff)
        {
            std::copy(fullH,fullH+nStr_a2*nStr_b2*nStr_p-nStr_a*nStr_b*(nStr_p-1),Col+p*nStr_a*nStr_b);
        }
        
        // Build the Proton-Proton block
        std::fill_n(tmpH,nStr_p2,MatsT(0.0));
        size_t nStr_p_nThread = nStr_p * nThreads;


        CASHelper<MatsT,IntsT>::buildFullHOneParticle(mcwfn,tmpH,exList_p,"PhCoreP_Correlated_Space","PRI_Correlated_Space");

   
    // Transpose the full Hamiltonian to make the primary index proton
    IMatCopy('T',nStr_ab,nStr_p*NDet,MatsT(1.0),fullH,nStr_ab,nStr_p*NDet);

    // Update the full Hamiltonian by blockwise adding the just calculated beta piece
    lastDimOff = nStr_p2*nStr_ab;
    Col = fullH;
    for(size_t ab = 0; ab < nStr_ab; ab++, Col+=lastDimOff)
    {
        IMatCopy('T',NDet,nStr_p,MatsT(1.),Col,NDet,nStr_p);
        MatAdd('N','N',nStr_p,nStr_p,MatsT(1.),Col+ab*nStr_p2,nStr_p,MatsT(1.),tmpH,nStr_p,Col+ab*nStr_p2,nStr_p);
        IMatCopy('T',nStr_p,NDet,MatsT(1.0),Col,nStr_p,NDet);
    }
    
    IMatCopy('T',nStr_p*NDet,nStr_ab,MatsT(1.),fullH,nStr_p*NDet,nStr_ab);

    MatsT * aPSCR = CQMemManager::get().template malloc<MatsT>(nStr_ap*nStr_ap);
    std::fill_n(aPSCR,nStr_ap*nStr_ap,MatsT(0.0));

    CASHelper<MatsT,IntsT>::buildFullHTwoParticle(mcwfn,aPSCR,exList_a,exList_p,"eP_Correlated_Space");

    // Temporary storage for our blockwise calculation of 
    MatsT * ZeroMat = CQMemManager::get().template malloc<MatsT>(NDet*NDet);
    std::fill_n(ZeroMat,NDet*NDet,MatsT(0.0));

    // Add the alpha-Proton to the main Hamiltonian
    Col = ZeroMat;
    MatsT * ColBeta, * aPTemp;
    aPTemp = aPSCR;
    for(size_t p = 0; p < nStr_p; p++, Col+=NDet*nStr_ab, aPTemp+=nStr_ap*nStr_a)
    {
        ColBeta = Col;
        for(size_t b = 0; b < nStr_b; b++, ColBeta+=nStr_a)
        {
            for(size_t a = 0; a < nStr_a; a++, ColBeta+=NDet)
            {
                // Note the negative here since (ee|PP) integrals do not carry the attractive sign!
                MatAdd('N','N',nStr_a,nStr_p,MatsT(-1.0),aPTemp+a*nStr_ap,nStr_a,MatsT(1.0),ColBeta,nStr_ab,ColBeta,nStr_ab);
            }
        }
    }

    MatAdd('N','N',NDet,NDet,MatsT(1.),ZeroMat,NDet,MatsT(1.),fullH,NDet,fullH,NDet);

    MatsT * bPSCR = CQMemManager::get().template malloc<MatsT>(nStr_bp*nStr_bp);
    std::fill_n(bPSCR,nStr_bp*nStr_bp,MatsT(0.0));

    CASHelper<MatsT,IntsT>::buildFullHTwoParticle(mcwfn,bPSCR,exList_b,exList_p,"eP_Correlated_Space");

    std::fill_n(ZeroMat,NDet*NDet,MatsT(0.0));

    MatsT * ColAlpha, * bPTemp;
    Col = ZeroMat;
    bPTemp = bPSCR;
    ColAlpha = ZeroMat;
    for(size_t p = 0; p < nStr_p; p++, bPTemp+=nStr_bp*nStr_b)
    {
        for(size_t b = 0; b < nStr_b; b++)
        {
            for(size_t a = 0; a < nStr_a; a++, ColAlpha+=NDet)
            {
                MatAdd('N','N',1,nStr_bp,-MatsT(1.0),bPTemp+b*nStr_bp,1,MatsT(1.0),ColAlpha+a,nStr_a,ColAlpha+a,nStr_a);
            }
        }
    }
    
    MatAdd('N','N',NDet,NDet,MatsT(1.),ZeroMat,NDet,MatsT(1.),fullH,NDet,fullH,NDet);

    CQMemManager::get().free(abSCR,aPSCR,bPSCR,tmpH,ZeroMat);

    return;
    } // NEOCASCI::buildFullH

template <typename MatsT, typename IntsT>
void NEOCASCI<MatsT,IntsT>::buildDiagH(MCWaveFunction<MatsT,IntsT>&mcwfn, MatsT * diagH)
{
    NEOCASCI_LOOP_INIT();

    // Calculate offset sizes we'll need in constructing the full matrix
    size_t nStr_ab = nStr_a*nStr_b;
    size_t nStr_bp = nStr_b*nStr_p;
    size_t nStr_ap = nStr_a*nStr_p;
    size_t nActEA = neowfn->ewfn_->MOPartition.nCorrEA;
    size_t nActEB = neowfn->ewfn_->MOPartition.nCorrEB;
    size_t nActP  = neowfn->pwfn_->MOPartition.nCorrEA;

    MatsT * ediag = CQMemManager::get().template malloc<MatsT>(nStr_ab);
    MatsT * ones = CQMemManager::get().template malloc<MatsT>(nStr_ab);
    std::fill_n(ones,nStr_ab,MatsT(1.0));

    // Uses the standard CASCI diagonal H builder
    CASCI<MatsT,IntsT>::buildDiagH(*(neowfn->ewfn_),ediag);

    MatsT * SCR = CQMemManager::get().template malloc<MatsT>(std::max(nStr_ap,nStr_bp));
    CASHelper<MatsT,IntsT>::buildDiagHOneParticle(mcwfn,SCR,exList_p,"PhCoreP_Correlated_Space","PRI_Correlated_Space");

    // Add the proton block onto the electronic block
    // Note what's happening is diagH <- ProtonDiagonalValue * vector_of_1s_size_nStr_ab + ediag
    for(size_t i = 0; i < nStr_p; i++)
    {
        MatAdd('N','N',1,nStr_ab,SCR[i],ones,1,MatsT(1.),ediag,1,diagH+i*nStr_ab,1);
    }

    std::fill_n(SCR,nStr_ap,MatsT(0.0));
    CASHelper<MatsT,IntsT>::buildDiagHTwoParticle(mcwfn,SCR,nActEA,nActP,exList_a,exList_p,"eP_Correlated_Space");
    CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_ab,nStr_p,diagH);
    CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_a,nStr_p,SCR);
    for(size_t b = 0; b < nStr_b; b++)
        blas::axpy(nStr_ap,-MatsT(1.0),SCR,1,diagH+b*nStr_ap,1);
    CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_p,nStr_ab,diagH);

    std::fill_n(SCR,nStr_bp,MatsT(0.0));
    CASHelper<MatsT,IntsT>::buildDiagHTwoParticle(mcwfn,SCR,nActEB,nActP,exList_b,exList_p,"eP_Correlated_Space");
    for(size_t a = 0; a < nStr_a; a++)
        blas::axpy(nStr_bp,-MatsT(1.0),SCR,1,diagH+a,nStr_a);

    CQMemManager::get().free(SCR,ones,ediag);

    return;
}; // NEOCASCI::buildDiagH

  template <typename MatsT, typename IntsT>
  void NEOCASCI<MatsT,IntsT>::buildSigma(MCWaveFunction<MatsT, IntsT> & mcwfn, 
    size_t nVec, MatsT * C, MatsT * Sigma) {

        NEOCASCI_LOOP_INIT();

        // Calculate offset sizes we'll need in constructing the full matrix
        size_t nStr_ab = nStr_a*nStr_b;
        size_t nStr_bp = nStr_b*nStr_p;
        size_t nStr_ap = nStr_a*nStr_p;

        // empty the sigma guess
        std::fill_n(Sigma,NDet*nVec,MatsT(0.0));

        // Alpha-Alpha block, storage is |a,b,p>
        CASHelper<MatsT,IntsT>::buildSigmaOneParticle(mcwfn,C,Sigma,nVec,nStr_bp,exList_a,"hCoreP_Correlated_Space","ERI_Correlated_Space");

        // Do Alpha-Beta block (hard coded to work with 1C ONLY for now)
        CASHelper<MatsT,IntsT>::buildSigmaTwoParticle(mcwfn,C,Sigma,nVec,nStr_p,exList_a,exList_b,"ERI_Correlated_Space");

        // Transpose |a,b,p> -> |b,a,p>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec*nStr_p,nStr_a,nStr_b,C,Sigma);

        // Do the Beta-Beta block
        CASHelper<MatsT,IntsT>::buildSigmaOneParticle(mcwfn,C,Sigma,nVec,nStr_ap,exList_b,"hCoreP_Correlated_Space","ERI_Correlated_Space");

        // Transpose |b,a,p> -> |a,p,b>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec,nStr_b,nStr_ap,C,Sigma);

        // Do the Alpha-Proton block
        CASHelper<MatsT,IntsT>::buildSigmaTwoParticle(mcwfn,C,Sigma,nVec,nStr_b,exList_a,exList_p,"eP_Correlated_Space",true);

        // Transpose |a,p,b> -> |b,a,p>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec,nStr_ap,nStr_b,C,Sigma);

        // Transpose |b,a,p> -> |p,b,a>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec,nStr_ab,nStr_p,C,Sigma);

        // Do the Proton-Proton block
        CASHelper<MatsT,IntsT>::buildSigmaOneParticle(mcwfn,C,Sigma,nVec,nStr_ab,exList_p,"PhCoreP_Correlated_Space","PRI_Correlated_Space");

        // Transpose |p,b,a> -> |b,p,a>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec*nStr_a,nStr_p,nStr_b,C,Sigma);

        // Do the Beta-Proton block
        CASHelper<MatsT,IntsT>::buildSigmaTwoParticle(mcwfn,C,Sigma,nVec,nStr_a,exList_b,exList_p,"eP_Correlated_Space",true);

        // Transpose |b,p,a> -> |a,b,p>
        CASHelper<MatsT,IntsT>::transposeVectors(nVec,nStr_bp,nStr_a,C,Sigma);

        return;

    } // NEOCASCI::buildSigma

    template <typename MatsT, typename IntsT>
    void NEOCASCI<MatsT,IntsT>::computeOneRDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * C, cqmatrix::Matrix<MatsT> & oneRDM)
    {
        computeTDM(mcwfn,C,C,oneRDM);
    } // NEOCASCI::computeOneRDM

    template <typename MatsT, typename IntsT>
    void NEOCASCI<MatsT,IntsT>::computeTDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * Cm, MatsT * Cn, cqmatrix::Matrix<MatsT> & TDM)
    {
        NEOCASCI_LOOP_INIT();

        TDM.clear();

        CASHelper<MatsT,IntsT>::computeTDM(mcwfn,Cm,Cn,NDet/nStr_a,exList_a,TDM);

        if(nC != 1) return;

        // Need to avoid transposing the same vector twice since Cm == Cn for RDM calculation
        if(Cm != Cn) CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_a,nStr_b*nStr_p,Cm,Cn);
        else CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_a,nStr_b*nStr_p,Cm);

        CASHelper<MatsT,IntsT>::computeTDM(mcwfn,Cm,Cn,NDet/nStr_b,exList_b,TDM);

        // Need to avoid transposing the same vector twice since Cm == Cn for RDM calculation
        if(Cm != Cn) CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_b*nStr_p,nStr_a,Cm,Cn);
        else CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_b*nStr_p,nStr_a,Cm);

        return;
    
    } // NEOCASCI::computeTDM

    template <typename MatsT, typename IntsT>
    void NEOCASCI<MatsT,IntsT>::computePOneRDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * C, cqmatrix::Matrix<MatsT> & PoneRDM)
    {
        computePTDM(mcwfn,C,C,PoneRDM);
    } // NEOCASCI::computePOneRDM

    template <typename MatsT, typename IntsT>
    void NEOCASCI<MatsT,IntsT>::computePTDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * Cm, MatsT * Cn, cqmatrix::Matrix<MatsT> & PTDM)
    {
        NEOCASCI_LOOP_INIT();

        PTDM.clear();

        // First transpose vectors since default storage is |a,b,p>
        if(Cm != Cn) CASHelper<MatsT,IntsT>::transposeVectors(1,NDet/nStr_p,nStr_p,Cm,Cn);
        else CASHelper<MatsT,IntsT>::transposeVectors(1,NDet/nStr_p,nStr_p,Cm);

        CASHelper<MatsT,IntsT>::computeTDM(mcwfn,Cm,Cn,NDet/nStr_p,exList_p,PTDM);

        // Transpose back to default storage
        if(Cm != Cn) CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_p,NDet/nStr_p,Cm,Cn);
        else CASHelper<MatsT,IntsT>::transposeVectors(1,nStr_p,NDet/nStr_p,Cm);

        return;

    } // NEOCASCI::computePTDM


}; // namespace ChronusQ