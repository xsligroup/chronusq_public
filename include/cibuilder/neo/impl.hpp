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
        MatsT * SCR = CQMemManager::get().template malloc<MatsT>(nSCR * nThreads);

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

        // indicies for parallel block
        int i,j,k,l,I,J,K,L,Ia,Ja,Ka,La,Ib,Jb,Kb,Lb,Ip,Jp,Kp,Lp;
        double signij,signkl,signIJ,signKL;
        double small_number = std::numeric_limits<double>::epsilon();

        //size_t nStr_a_nThread = nStr_a * nThreads;

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


#pragma omp parallel default(shared) private(Col,SCR_ith,Lp,k,l,Kp,signkl,i,j,Jp,signij)
{
        // Get this threads ID
        auto myThread = GetThreadID();
        // Offset the CIHCol by the thread amount
        Col = tmpH + nStr_p * myThread;

        // Each thread gets a column for an individual determinant
        for(Lp = myThread; Lp < nStr_p; Lp+=nThreads, Col+=nStr_p_nThread)
        {
            // Zero out my column in the SCR space
            std::fill_n(Col, nStr_p, MatsT(0.0));
            // Get the excitation list for my current det
            const int * exList_Lp = exList_p->pointerAtDet(Lp);

            // Loop over the non-zero excitations
            // Incriment both the current Ekl counter and the excitation list pointer
            for(auto Ekl = 0ul; Ekl < nNZp; Ekl++, exList_Lp+=4)
            {
                // Get the indices of the excitation connection
                UNPACK_EXCITATIONLIST_4(exList_Lp,k,l,Kp,signkl);
                // The connection between this det (indexed by La, already offset 
                // by La*nThread in SCR_ith) is the one particle contribution
                Col[Kp] += signkl * PhCore(k,l);

                // Now that we have a singly excited determinant, generate the 
                // doubly excited determinants off of this by looking at the
                // singles from the Ka'th det
                const int * exList_Kp = exList_p->pointerAtDet(Kp);
                for(auto Eij = 0ul; Eij < nNZp; Eij++, exList_Kp+=4)
                {
                    // Grab the now doubly excited contribution
                    UNPACK_EXCITATIONLIST_4(exList_Kp,i,j,Jp,signij);
                    // Add the two electron parts
                    Col[Jp] += 0.5 * signij * signkl * PRI(i,j,k,l);
                }
            }
        }
}
    
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
    
    // Alpha-Proton double excitations
#pragma omp parallel for schedule(static) default(shared) private(Col,Lp,K,L,Kp,signKL,La,i,j,Ka,signij)
    for(Lp = 0; Lp < nStr_p; Lp++)
    {
        Col = aPSCR + nStr_ap*(Lp*nStr_a);
        const int * exList_Lp_head = exList_p->pointerAtDet(Lp);
        for(La = 0; La < nStr_a; La++, Col+=nStr_ap)
        {
            const int * exList_La_head = exList_a->pointerAtDet(La);
            const int * exList_Lp = exList_Lp_head;
            for(size_t Ekl = 0; Ekl < nNZp; Ekl++, exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,K,L,Kp,signKL);
                const int * exList_La = exList_La_head;
                for(size_t Eij = 0; Eij < nNZa; Eij++, exList_La+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_La,i,j,Ka,signij);
                    Col[Ka+Kp*nStr_a]-=signij*signKL*ePRI(i,j,K,L);
                }
            }
        }
    }

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
                MatAdd('N','N',nStr_a,nStr_p,MatsT(1.0),aPTemp+a*nStr_ap,nStr_a,MatsT(1.0),ColBeta,nStr_ab,ColBeta,nStr_ab);
            }
        }
    }

    MatAdd('N','N',NDet,NDet,MatsT(1.),ZeroMat,NDet,MatsT(1.),fullH,NDet,fullH,NDet);

    MatsT * bPSCR = CQMemManager::get().template malloc<MatsT>(nStr_bp*nStr_bp);
    std::fill_n(bPSCR,nStr_bp*nStr_bp,MatsT(0.0));


    // Beta-Proton double excitations
#pragma omp parallel for schedule(static) default(shared) private(Col,Lp,K,L,Kp,signKL,Lb,i,j,Kb,signij)
    for(Lp = 0; Lp < nStr_p; Lp++)
    {
        Col = bPSCR + nStr_bp*(Lp*nStr_b);
        const int * exList_Lp_head = exList_p->pointerAtDet(Lp);
        for(Lb = 0; Lb < nStr_b; Lb++, Col+=nStr_bp)
        {
            const int * exList_Lb_head = exList_b->pointerAtDet(Lb);
            const int * exList_Lp = exList_Lp_head;
            for(size_t Ekl = 0; Ekl < nNZp; Ekl++, exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,K,L,Kp,signKL);
                const int * exList_Lb = exList_Lb_head;
                for(size_t Eij = 0; Eij < nNZb; Eij++, exList_Lb+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_Lb,i,j,Kb,signij);
                    Col[Kb+Kp*nStr_b]-=signij*signKL*ePRI(i,j,K,L);
                }
            }
        }
    }


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
                MatAdd('N','N',1,nStr_bp,MatsT(1.0),bPTemp+b*nStr_bp,1,MatsT(1.0),ColAlpha+a,nStr_a,ColAlpha+a,nStr_a);
            }
        }
    }
    
    MatAdd('N','N',NDet,NDet,MatsT(1.),ZeroMat,NDet,MatsT(1.),fullH,NDet,fullH,NDet);

    CQMemManager::get().free(SCR,abSCR,aPSCR,bPSCR,tmpH,ZeroMat);

    return;
    } // NEOCASCI::buildFullH

template <typename MatsT, typename IntsT>
void NEOCASCI<MatsT,IntsT>::buildDiagH(MCWaveFunction<MatsT,IntsT>&mcwfn, MatsT * diagH)
{
    NEOCASCI_LOOP_INIT();

    // Grab integrals
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
    size_t nActEA = neowfn->ewfn_->MOPartition.nCorrEA;
    size_t nActEB = neowfn->ewfn_->MOPartition.nCorrEB;
    size_t nActP  = neowfn->pwfn_->MOPartition.nCorrEA;

    MatsT * ediag = CQMemManager::get().template malloc<MatsT>(nStr_ab);
    MatsT * ones = CQMemManager::get().template malloc<MatsT>(nStr_ab);
    std::fill_n(ones,nStr_ab,MatsT(1.0));

    // Uses the standard CASCI diagonal H builder
    CASCI<MatsT,IntsT>::buildDiagH(*(neowfn->ewfn_),ediag);

    // Cannot do this until buildDiagH has the option to have a string passsed in
    // which would point to the proper sets of one and two particle integrals
    // CASCI<MatsT,IntsT>::buildDiagH(*(neowfn->pwfn_),diagH+nStr_ab);
    // for(size_t i = 0; i < nStr_p; i++)
        //std::cout << i << " -> " << diagH[i+nStr_ab] << std::endl;

    MatsT * SCR = CQMemManager::get().template malloc<MatsT>(nStr_p);
    MatsT tmp;

    int i,j,k,l,I,J,K,L,Ia,Ja,Ka,La,Ib,Jb,Kb,Lb,Ip,Jp,Kp,Lp;
    double signij,signkl,signIJ,signKL;
// Calculate the proton-proton diagonal block
#pragma omp parallel for schedule(static) default(shared) private(tmp,Lp,k,l,Kp,signkl,i,j,Ja,signij)
    for(Lp = 0; Lp < nStr_p; Lp++)
    {
        tmp = MatsT(0.);
        const int * exList_Lp = exList_p->pointerAtDet(Lp);
        for(size_t Ekl = 0; Ekl < nNZp; Ekl++, exList_Lp+=4)
        {
            UNPACK_EXCITATIONLIST_4(exList_Lp,k,l,Kp,signkl);
            if(Kp==Lp)
            {
                tmp += signkl * PhCore(k,l);
            }
            const int * exList_Kp = exList_p->pointerAtDet(Kp);
            for(size_t Eij = 0; Eij < nNZp; Eij++, exList_Kp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Kp,i,j,Jp,signij);
                if(Jp==Lp)
                {
                    tmp += 0.5 * signij * signkl * PRI(i,j,k,l);
                }
            }
        }
        SCR[Lp]=tmp;
    }

    // Add the proton block onto the electronic block
    for(size_t i = 0; i < nStr_p; i++)
    {
        MatAdd('N','N',1,nStr_ab,SCR[i],ones,1,MatsT(1.),ediag,1,diagH+i*nStr_ab,1);
    }

    MatsT * dH;
    MatsT val;

    // Do the alpha-proton block
#pragma omp parallel for schedule(static) default(shared) private(dH,Lp,k,l,Kp,signkl,La,i,j,Ka,signij)
    for(Lp = 0; Lp < nStr_p; Lp++)
    {
        dH = diagH + Lp*nStr_ab;
        const int * exList_Lp_head = exList_p->pointerAtDet(Lp);
        for(La = 0; La < nStr_a; La++, dH++)
        {
            val = MatsT(0.0);
            const int * exList_La_head=exList_a->pointerAtDet(La);
            const int * exList_Lp = exList_Lp_head;
            for(size_t Ekl = 0; Ekl < nActP; Ekl++, exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,k,l,Kp,signkl);
                const int * exList_La = exList_La_head;
                for(size_t Eij = 0; Eij < nActEA; Eij++, exList_La+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_La,i,j,Ka,signij);
                    val+= signij * signkl * ePRI(i,j,k,l);
                }
            }
            for(size_t b = 0; b < nStr_b; b++)
            {
                dH[b*nStr_a]-=val;
            }
        }
    }

    // Do the beta-prton block
#pragma omp parallel for schedule(static) default(shared) private(dH,Lp,k,l,Kp,signkl,Lb,i,j,Kb,signij)
    for(Lp = 0; Lp < nStr_p; Lp++)
    {
        dH = diagH + Lp*nStr_ab;
        const int * exList_Lp_head = exList_p->pointerAtDet(Lp);
        for(Lb = 0; Lb < nStr_b; Lb++, dH+=nStr_a)
        {
            val = MatsT(0.0);
            const int * exList_Lb_head = exList_b->pointerAtDet(Lb);
            const int * exList_Lp = exList_Lp_head;
            for(size_t Ekl = 0; Ekl < nActP; Ekl++, exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,k,l,Kp,signkl);
                const int * exList_Lb = exList_Lb_head;
                for(size_t Eij = 0; Eij < nActEB; Eij++, exList_Lb+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_Lb,i,j,Kb,signij);
                    val+=signij*signkl*ePRI(i,j,k,l);
                }
            }
            for(size_t a = 0; a < nStr_a; a++)
            {
                dH[a]-=val;
            }
        }
    }

    CQMemManager::get().free(SCR,ones,ediag);

    return;
}; // NEOCASCI::buildDiagH

  template <typename MatsT, typename IntsT>
  void NEOCASCI<MatsT,IntsT>::buildSigma(MCWaveFunction<MatsT, IntsT> & mcwfn, 
    size_t nVec, MatsT * C, MatsT * Sigma) {

        NEOCASCI_LOOP_INIT();

        auto & ehCore = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("hCoreP_Correlated_Space"));
        auto & PhCore = *(mcwfn.moints->template getIntegral<OnePInts,MatsT>("PhCoreP_Correlated_Space"));
        auto & eRI    = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));
        auto & PRI    = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("PRI_Correlated_Space"));
        auto & ePRI   = *(mcwfn.moints->template getIntegral<InCore4indexTPI,MatsT>("eP_Correlated_Space"));

//        CASCI<MatsT,IntsT>::buildSigma()

        // Calculate offset sizes we'll need in constructing the full matrix
        size_t nStr_a2 = nStr_a*nStr_a;
        size_t nStr_b2 = nStr_b*nStr_b;
        size_t nStr_p2 = nStr_p*nStr_p;
        size_t nStr_ab = nStr_a*nStr_b;
        size_t nStr_bp = nStr_b*nStr_p;
        size_t nStr_ap = nStr_a*nStr_p;

        // Memory allocations for temporary storage
        size_t nSCR = std::max(std::max(nStr_a,nStr_b),nStr_p);
        size_t nThreads = GetNumThreads();
        MatsT * SCR = CQMemManager::get().template malloc<MatsT>(nSCR*nThreads);

        // empty the sigma guess
        std::fill_n(Sigma,NDet*nVec,MatsT(0.0));

        // indicies for parallel block
        int i,j,k,l,I,J,K,L,Ia,Ja,Ka,La,Ib,Jb,Kb,Lb,Ip,Jp,Kp,Lp;
        double signij,signkl,signIJ,signKL;
        double small_number = std::numeric_limits<double>::epsilon();

        MatsT * HC, *Ci, *SCR_ith;
        MatsT val;


        // Alpha-Alpha block, storage is |a,b,p>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,La,k,l,Ka,signkl,i,j,Ja,signij,Kb,Kp)
{
        auto iThread = GetThreadID();
        SCR_ith = SCR + nSCR * iThread;
        for(Ka = iThread; Ka < nStr_a; Ka+=nThreads)
        {
            std::fill_n(SCR_ith,nStr_a,MatsT(0.0));
            const int * exList_Ka = exList_a->pointerAtDet(Ka);
            for(size_t Ekl = 0; Ekl < nNZa; Ekl++, exList_Ka+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Ka,k,l,La,signkl);
                SCR_ith[La] += signkl * ehCore(k,l);

                const int * exList_La = exList_a->pointerAtDet(La);
                
                for(size_t Eij = 0; Eij < nNZa; Eij++, exList_La+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_La,i,j,Ja,signij);
                    SCR_ith[Ja]+=0.5*signij*signkl*eRI(i,j,k,l);
                }
            }
            
            // Screen for small elements
            std::vector<int> SCR_nonZero_ith;
            for(La = 0; La < nStr_a; La++)
            {
                if(std::abs(SCR_ith[La])>small_number)
                {
                    SCR_nonZero_ith.push_back(La);
                }
            }

            // Multiply by sigma
            HC = Sigma;
            Ci = C;
            for(size_t iVec = 0; iVec < nVec; iVec++)
            {
                for(Kp = 0; Kp < nStr_p; Kp++)
                {
                    for(Kb = 0; Kb < nStr_b; Kb++, HC+=nStr_a, Ci+=nStr_a)
                    {
                        for(size_t iSCR = 0; iSCR < SCR_nonZero_ith.size(); iSCR++)
                        {
                            La = SCR_nonZero_ith[iSCR];
                            HC[Ka]+=SCR_ith[La]*Ci[La];
                        }
                    }
                }
            }
        } 
}

    // Alpha-Beta block, storage is |a,b,p>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,Lb,k,l,Kb,signkl,i,j,La,signij,Ka,Kp,val)
{
    int KAddr,LAddr;
    auto iThread = GetThreadID();
    for(Kb = iThread; Kb < nStr_b; Kb+=nThreads)
    {
        const int * exList_Kb_head = exList_b->pointerAtDet(Kb);
        for(Ka = 0, KAddr=Kb*nStr_a; Ka < nStr_a; Ka++, KAddr++)
        {
            const int * exList_Ka_head = exList_a->pointerAtDet(Ka);
            const int * exList_Kb = exList_Kb_head;
            for(size_t Ekl = 0; Ekl < nNZb; Ekl++, exList_Kb+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Kb,l,k,Lb,signkl);
                const int * exList_Ka = exList_Ka_head;
                for(size_t Eij = 0; Eij < nNZa; Eij++, exList_Ka+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_Ka,i,j,La,signij);
                    val = signij * signkl * eRI(i,j,k,l);
                    if(std::abs(val)>small_number)
                    {
                        HC = Sigma;
                        Ci = C;
                        LAddr = La + Lb*nStr_a;
                        for(size_t iVec = 0; iVec < nVec; iVec++)
                        {
                            for(size_t p = 0; p < nStr_p; p++, HC+=nStr_ab, Ci+=nStr_ab)
                            {
                                HC[KAddr]+= val * Ci[LAddr];
                            }
                        }
                    }
                }
            }
        }
    }
}

    // Blockwise (on Proton index) to (b,a,p) indexing
    // |a,b,p> --> |b,a,p>
    HC = Sigma;
    Ci = C;
    for(size_t iVec = 0; iVec < nVec; iVec++)
    {
        for(size_t p = 0; p < nStr_p; p++, HC+=nStr_ab, Ci+=nStr_ab)
        {
            IMatCopy('T',nStr_a,nStr_b,MatsT(1.),HC,nStr_a,nStr_b);
            IMatCopy('T',nStr_a,nStr_b,MatsT(1.),Ci,nStr_a,nStr_b);
        }
    }

    // Beta-beta block, storage is |b,a,p>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,Lb,k,l,Kb,signkl,i,j,Jb,signij,Ka,Kp)
{
    auto iThread = GetThreadID();
    SCR_ith = SCR + nSCR * iThread;
    for(Kb = iThread; Kb < nStr_b; Kb+=nThreads)
    {
        std::fill_n(SCR_ith,nStr_b,MatsT(0.0));
        const int * exList_Kb = exList_b->pointerAtDet(Kb);
        for(size_t Elk = 0; Elk < nNZb; Elk++, exList_Kb+=4)
        {
            UNPACK_EXCITATIONLIST_4(exList_Kb,k,l,Lb,signkl);
            SCR_ith[Lb]+=signkl*ehCore(k,l);
            const int * exList_Lb = exList_b ->pointerAtDet(Lb);
            for(size_t Eij = 0; Eij < nNZb; Eij++,exList_Lb+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lb,i,j,Jb,signij);
                SCR_ith[Jb]+=0.5*signij*signkl*eRI(i,j,k,l);
            }
        }
        
        // Screen for small elements
        std::vector<int> SCR_nonZero_ith;
        for(Lb = 0; Lb < nStr_b; Lb++)
        {
            if(std::abs(SCR_ith[Lb])>small_number)
            {
                SCR_nonZero_ith.push_back(Lb);
            }
        }
        
        // Multiply by sigma
        HC = Sigma;
        Ci = C;
        for(size_t iVec = 0; iVec < nVec; iVec++)
        {
            for(Kp = 0; Kp < nStr_p; Kp++)
            {
                for(Ka = 0; Ka < nStr_a; Ka++, HC+=nStr_b, Ci+=nStr_b)
                {
                    for(size_t iSCR = 0; iSCR < SCR_nonZero_ith.size(); iSCR++)
                    {
                        Lb = SCR_nonZero_ith[iSCR];
                        HC[Kb]+=SCR_ith[Lb]*Ci[Lb];
                    }
                }
            }
        }
    }
}

    // Transpose proton index to primary
    // |b,a,p> --> |p,b,a>
    HC = Sigma;
    Ci = C;
    for(size_t iVec = 0; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet)
    {
        IMatCopy('T',nStr_ab,nStr_p,MatsT(1.),HC,nStr_ab,nStr_p);
        IMatCopy('T',nStr_ab,nStr_p,MatsT(1.),Ci,nStr_ab,nStr_p);
    }

// Do the proton-proton block, storage is |p,b,a>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,Lp,k,l,Kp,signkl,i,j,Jp,signij,Ka,Kb)
{
    auto iThread = GetThreadID();
    SCR_ith = SCR + nSCR * iThread;
    for(Kp = iThread; Kp < nStr_p; Kp+=nThreads)
    {
        std::fill_n(SCR_ith,nStr_p,MatsT(0.0));
        const int * exList_Kp = exList_p->pointerAtDet(Kp);
        for(size_t Elk = 0; Elk < nNZp; Elk++, exList_Kp+=4)
        {
            UNPACK_EXCITATIONLIST_4(exList_Kp,k,l,Lp,signkl);
            SCR_ith[Lp]+=signkl*PhCore(k,l);
            const int * exList_Lp = exList_p->pointerAtDet(Lp);
            for(size_t Eij = 0; Eij < nNZp; Eij++,exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,i,j,Jp,signij);
                SCR_ith[Jp]+=0.5*signij*signkl*PRI(i,j,k,l);
            }
        }
        
        // Screen for small elements
        std::vector<int> SCR_nonZero_ith;
        for(Lp = 0; Lp < nStr_p; Lp++)
        {
            if(std::abs(SCR_ith[Lp])>small_number)
            {
                SCR_nonZero_ith.push_back(Lp);
            }
        }
        
        // Multiply by sigma
        HC = Sigma;
        Ci = C;
        for(size_t iVec = 0; iVec < nVec; iVec++)
        {
            for(Ka = 0; Ka < nStr_a; Ka++)
            {
                for(Kb = 0; Kb < nStr_b; Kb++, HC+=nStr_p, Ci+=nStr_p)
                {
                    for(size_t iSCR = 0; iSCR < SCR_nonZero_ith.size(); iSCR++)
                    {
                        Lp = SCR_nonZero_ith[iSCR];
                        HC[Kp]+=SCR_ith[Lp]*Ci[Lp];
                    }
                }
            }
        }
    }
}

    // Proton-beta block, storage is |p,b,a>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,Lb,k,l,Kb,signkl,i,j,Lp,signij,Kp,Ka,val)
{
    int KAddr,LAddr;
    auto iThread = GetThreadID();
    for(Kb = iThread; Kb < nStr_b; Kb+=nThreads)
    {
        const int * exList_Kb_head = exList_b->pointerAtDet(Kb);
        for(Kp = 0, KAddr=Kb*nStr_p; Kp < nStr_p; Kp++, KAddr++)
        {
            const int * exList_Kp_head = exList_p->pointerAtDet(Kp);
            const int * exList_Kb = exList_Kb_head;
            for(size_t Ekl = 0; Ekl < nNZb; Ekl++, exList_Kb+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Kb,l,k,Lb,signkl);
                const int * exList_Kp = exList_Kp_head;
                for(size_t Eij = 0; Eij < nNZp; Eij++, exList_Kp+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_Kp,i,j,Lp,signij);
                    val = signij * signkl * ePRI(k,l,i,j);
                    if(std::abs(val)>small_number)
                    {
                        HC = Sigma;
                        Ci = C;
                        LAddr = Lp + Lb*nStr_p;
                        for(size_t iVec = 0; iVec < nVec; iVec++)
                        {
                            for(size_t a = 0; a < nStr_a; a++, HC+=nStr_bp, Ci+=nStr_bp)
                            {
                                HC[KAddr]-= val * Ci[LAddr];
                            }
                        }
                    }
                }
            }
        }
    }
}
    // TODO: Figure out how to transpose directly from |p,b,a> to |p,a,b>
    // |p,b,a> --> |a,p,b> --> |p,a,b>

    // |p,b,a> --> |a,p,b>
    HC = Sigma;
    Ci = C;
    for(size_t iVec = 0; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet)
    {
        IMatCopy('T',nStr_bp,nStr_a,MatsT(1.0),HC,nStr_bp,nStr_a);
        IMatCopy('T',nStr_bp,nStr_a,MatsT(1.0),Ci,nStr_bp,nStr_a);
    }

    // |a,p,b> --> |p,a,b>
    HC = Sigma;
    Ci = C;
    for(size_t iVec = 0; iVec < nVec; iVec++)
    {
        for(size_t b = 0; b < nStr_b; b++, HC+=nStr_ap, Ci+=nStr_ap)
        {
            IMatCopy('T',nStr_a,nStr_p,MatsT(1.),HC,nStr_a,nStr_p);
            IMatCopy('T',nStr_a,nStr_p,MatsT(1.),Ci,nStr_a,nStr_p);
        }
    } 

    // Proton-alpha block, storage is |p,a,b>
#pragma omp parallel default(shared) private(HC,Ci,SCR_ith,La,k,l,Ka,signkl,i,j,Lp,signij,Kp,Kb,val)
{
    int KAddr,LAddr;
    auto iThread = GetThreadID();
    for(Ka = iThread; Ka < nStr_a; Ka+=nThreads)
    {
        const int * exList_Ka_head = exList_a->pointerAtDet(Ka);
        for(Kp = 0, KAddr=Ka*nStr_p; Kp < nStr_p; Kp++, KAddr++)
        {
            const int * exList_Kp_head = exList_p->pointerAtDet(Kp);
            const int * exList_Ka = exList_Ka_head;
            for(size_t Ekl = 0; Ekl < nNZa; Ekl++, exList_Ka+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Ka,l,k,La,signkl);
                const int * exList_Kp = exList_Kp_head;
                for(size_t Eij = 0; Eij < nNZp; Eij++, exList_Kp+=4)
                {
                    UNPACK_EXCITATIONLIST_4(exList_Kp,i,j,Lp,signij);
                    val = signij * signkl * ePRI(k,l,i,j);
                    if(std::abs(val)>small_number)
                    {
                        HC = Sigma;
                        Ci = C;
                        LAddr = Lp + La*nStr_p;
                        for(size_t iVec = 0; iVec < nVec; iVec++)
                        {
                            for(size_t b = 0; b < nStr_b; b++, HC+=nStr_ap, Ci+=nStr_ap)
                            {
                                HC[KAddr] -= val * Ci[LAddr];
                            }
                        }
                    }
                }
            }
        }
    }
}
    
    // Transpose back to |a,b,p>
    // |p,a,b> --> |a,b,p>
    HC = Sigma;
    Ci = C;
    for(size_t iVec = 0; iVec < nVec; iVec++, HC+=NDet, Ci+=NDet)
    {
        IMatCopy('T',nStr_p,nStr_ab,MatsT(1.),HC,nStr_p,nStr_ab);
        IMatCopy('T',nStr_p,nStr_ab,MatsT(1.),Ci,nStr_p,nStr_ab);
    }

    CQMemManager::get().free(SCR);

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

        size_t nThreads = GetNumThreads();
        std::vector<cqmatrix::Matrix<MatsT>> SCR;
        for(size_t i = 0; i < nThreads; i++)
            SCR.emplace_back(TDM.nRows());

        // alpha electron part
        int k,l,La,Lb,Lp,Ka,Kb,Kp;
        double signkl;
#pragma omp parallel default(shared) private(La,Lb,Lp,k,l,Ka,signkl)
    {
        auto iThread = GetThreadID();
        auto & tmpRDM = SCR[iThread];
        tmpRDM.clear();
        for(La = iThread; La < nStr_a; La+=nThreads)
        {
            const int * exList_La = exList_a->pointerAtDet(La);
            for(size_t Ekl = 0; Ekl < nNZa; Ekl++, exList_La+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_La,k,l,Ka,signkl);
                auto tmp = MatsT(0.0);
                for(Lb = 0; Lb < nStr_b; Lb++)
                {
                    for(Lp = 0; Lp < nStr_p; Lp++)
                    {
                        tmp += SmartConj(Cm[Ka + Lb*nStr_a + Lp*nStr_a*nStr_b])*Cn[La + Lb*nStr_a + Lp*nStr_a*nStr_b];
                    }
                }
                tmpRDM(k,l) += tmp*signkl;
            }
        }
    }

    // Beta part
#pragma omp parallel default(shared) private(La,Lb,Lp,k,l,Kb,signkl)
    {
        auto iThread = GetThreadID();
        auto & tmpRDM = SCR[iThread];
        for(Lb = iThread; Lb < nStr_b; Lb+=nThreads)
        {
            const int * exList_Lb = exList_b->pointerAtDet(Lb);
            for(size_t Ekl = 0; Ekl < nNZb; Ekl++, exList_Lb+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lb,k,l,Kb,signkl);
                auto tmp = MatsT(0.0);
                for(La = 0; La < nStr_a; La++)
                {
                    for(Lp = 0; Lp < nStr_p; Lp++)
                    {
                        tmp += SmartConj(Cm[La + Kb*nStr_a + Lp*nStr_a*nStr_b])*Cn[La + Lb*nStr_a + Lp*nStr_a*nStr_b];
                    }
                }
                tmpRDM(k,l) += tmp*signkl;
            }
        }
    }

    TDM.clear();
    for(size_t i = 0; i < nThreads; i++)
    {
        TDM+=SCR[i];
    }

//    prettyPrintSmart(std::cout,"Electronic TDM",TDM.pointer(),TDM.dimension(),TDM.dimension(),TDM.dimension());

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

        size_t nThreads = GetNumThreads();
        std::vector<cqmatrix::Matrix<MatsT>> SCR;
        for(size_t i = 0; i < nThreads; i++)
            SCR.emplace_back(PTDM.nRows());

        // alpha electron part
        int k,l,La,Lb,Lp,Ka,Kb,Kp;
        double signkl;
#pragma omp parallel default(shared) private(La,Lb,Lp,k,l,Kp,signkl)
    {
        auto iThread = GetThreadID();
        auto & tmpRDM = SCR[iThread];
        tmpRDM.clear();
        for(Lp = iThread; Lp < nStr_p; Lp+=nThreads)
        {
            const int * exList_Lp = exList_p->pointerAtDet(Lp);
            for(size_t Ekl = 0; Ekl < nNZp; Ekl++, exList_Lp+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Lp,k,l,Kp,signkl);
                auto tmp = MatsT(0.0);
                for(Lb = 0; Lb < nStr_b; Lb++)
                {
                    for(La = 0; La < nStr_a; La++)
                    {
                        tmp += SmartConj(Cm[La + Lb*nStr_a + Kp*nStr_a*nStr_b])*Cn[La + Lb*nStr_a + Lp*nStr_a*nStr_b];
                    }
                }
                tmpRDM(k,l) += tmp*signkl;
            }
        }
    }

    PTDM.clear();
    for(size_t i = 0; i < nThreads; i++)
    {
        PTDM+=SCR[i];
    }

//    prettyPrintSmart(std::cout,"Protonic TDM",PTDM.pointer(),PTDM.dimension(),PTDM.dimension(),PTDM.dimension());

    } // NEOCASCI::computePTDM


}; // namespace ChronusQ