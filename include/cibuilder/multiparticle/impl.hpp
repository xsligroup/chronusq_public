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

#include <mcwavefunction.hpp>
#include <cibuilder/multiparticle.hpp>
#include <cibuilder/casci/helper.hpp>

#define MULTIPARTICLECI_LOOP_INIT() \
    auto mpmcwfn = dynamic_cast<MultiParticleMCWaveFunction<MatsT,IntsT>*>(&mcwfn);\
    size_t NDet = mcwfn.NDet; \
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> p1mcwfn; \
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> p2mcwfn; \
    std::shared_ptr<const ExcitationList> exList_p1_a; \
    std::shared_ptr<const ExcitationList> exList_p1_b; \
    std::shared_ptr<const ExcitationList> exList_p2_a; \
    std::shared_ptr<const ExcitationList> exList_p2_b; \
    std::vector<std::string> subsystems = mpmcwfn->getOrder();\
    std::vector<std::string> CIOrder = mpmcwfn->getCIOrder();\
    size_t nsystems  = subsystems.size();

namespace ChronusQ {

  // SMG 08/06/26
  // Notes on how this code works:
  //
  // The assumption is that our CI Vector is always stored in the order of CIOrder,
  // which is stored in the same order as the order_ class member but expanded to
  // included all alpha & beta pieces more or less as independent wavefunctions 
  // (alpha & beta still often share integral sets, however with this structure it
  // might not be too difficult to do an unrestricted CI calculation if one wanted
  // as a MultiParticleSS calculation)
  //
  // As an example, if MultiParticleSS has E, QP1, QP2, and QP3 in that particular
  // order, then CIOrder will be {EA, EB, QP1, QP2, QP3} (assuming QPN is only alpha)
  // 
  // The general strategy is to make it so one particle and pairwise particle 
  // interactions always operate on the primary one / two indices.  To achieve this,
  // We make use of two primary operations: 
  // Stashing of an index at the back of the CIVector (e.g., |1,2,3,4> -> |2,3,4,1>)
  // Shifting indices up to the second position (e.g., |1,2,3,4> -> |1,3,4,2>)
  // 
  // Any CI operation then is composed of the following stepwise operations:
  // 1) Operate on the particle who's index is currently first
  // 2) Operate on the pairwise interaction between the current first particle and
  //    the current second particle
  // 3) Shift the current second index to the back of the queue of particles which have
  //    yet to be operated on
  // 4) Repeat steps 2-3 for all pairwise interactions which haven't been accounted for
  // 5) Shift the current first index to the back of the CI vector, all of it's 
  //    interactions at this point have been accounted for and so we are essentially
  //    just stashing this index
  // 6) Repeat steps 1-5 until the all interactions have been accounted for
  //
  // Example:
  // For a CI vector with 4 particle types |1,2,3,4> (e.g., for a NEO calculation with two 
  // distinguishable quantum protons we will have EA, EB, QP1, QP2), the state of the CI
  // vector and which operations are done in which state proceeds as follows:
  // 
  // State       |   Hamiltonian Operations
  // |1,2,3,4>   |   1's one and two body Hamiltonian
  // |1,2,3,4>   |   Pairwise interaction between 1 & 2
  // |1,2,3,4>   ->  |1,3,4,2>
  // |1,3,4,2>   |   Pairwise interaction between 1 & 3
  // |1,3,4,2>   ->  |1,4,2,3>
  // |1,4,2,3>   |   Pairwise interaction between 1 & 4
  // |1,4,2,3>   ->  |1,2,3,4>
  // |1,2,3,4>   ->  |2,3,4,1>
  // |2,3,4,1>   |   2's one and two body Hamiltonain
  // |2,3,4,1>   |   Pairwise interaction between 2 & 3
  // |2,3,4,1>   ->  |2,4,3,1>
  // |2,4,3,1>   |   Pairwise interaction between 2 & 4
  // |2,4,3,1>   ->  |2,3,4,1>
  // |2,3,4,1>   ->  |3,4,1,2>
  // |3,4,1,2>   |   3's one and two body Hamiltonian
  // |3,4,1,2>   |   Pairwise interaction between 3 & 4
  // |3.4,1,2>   ->  |4,1,2,3>
  // |4,1,2,3>   |   4's one and two body Hamiltonian
  // |4,1,2,3>   ->  |1,2,3,4>
  //
  // At this point, all terms in the Hamiltonian have been accounted for, and the CIVector
  // has been returned to the standard storage state (which is important for calculating RDM's
  // and printing other properties)
  //
  // For the FullH calcluation, documentation is included within for how the row / column indices
  // are handled (although it is essentially the same procedure twice, once for rows and once for columns)


  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::buildFullH(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * fullH)
  {
    MULTIPARTICLECI_LOOP_INIT();

    // SMG 07/30/26
    // The storage of the fullH is by default
    // (Lp1, Lp2, Lp3, ... LpN, Rp1, Rp2, Rp3 ... RpN)
    // Where L & R are the row / column indices, respectively (might be the other way but it doesn't really matter tbh)
    // For one particle, the easiest way to handle things is when the indices we care about are either both primary
    // (e.g., Lp1 Rp1 Lp2 Lp3 ....) at which point we simply directly add (continuously in memory!)
    // the Lp1 Rp1 blockwise on diagonal everything else (note diagonal implies \delta LpN RpN for all N)

    // Zero out the Full Matrix
    std::fill_n(fullH,NDet*NDet,MatsT(0.0));

    // We'll stash indices at the back end once we've taken care of all their interactions,
    // so we need to keep track of how many dets are currently stashed
    size_t NDetStashed = 1;

    // Allocate storage space for the temporary matrix
    MatsT * tmpH = CQMemManager::get().template malloc<MatsT>(NDet*NDet);
    std::fill_n(tmpH,NDet*NDet,MatsT(0.0));
    
    size_t NDetp2;
    std::string pilabel, pjlabel;

    for(size_t pi = 0; pi < CIOrder.size(); pi++)
    {
        pilabel = CIOrder[pi];
        std::fill_n(tmpH,NDet*NDet,MatsT(0.0));    
        // Unpack the 1 particle objects
        auto piCIBuilder = mpmcwfn->getOneParticleCIHelper(pilabel);
        std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pimcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(piCIBuilder.refwfn);
        std::shared_ptr<const ExcitationList> exList_pi = piCIBuilder.exList;
        std::string onePIntsStr = piCIBuilder.oneBodyIntsStr;
        std::string twoPIntsStr = piCIBuilder.twoBodyIntsStr;
        size_t NDetp1 = piCIBuilder.NDet;

        // Build pi's self-Hamiltonian
        CASHelper<MatsT,IntsT>::buildFullHOneParticle(*pimcwfn,tmpH,exList_pi,onePIntsStr,twoPIntsStr);

        // Add this one particle to the full Hamiltonian
        // Note we're assuming the current storage is
        // (Lpi, Lpi+1, ... LpN, Rpi, Rpi+i, ... RpN)
        CASHelper<MatsT,IntsT>::addBlockToMatrixBlockDiagonal(NDetp1,NDet,MatsT(1.0),tmpH,fullH);
        // Now begin building the pairwise interactions
        for(const auto & interaction : mpmcwfn->getTwoParticleCIHelper(pilabel))
        {
            std::fill_n(tmpH,NDet*NDet,MatsT(0.0));
            pjlabel = interaction.labelp2;
            std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pjmcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(interaction.mcwfnp2);
            std::shared_ptr<const ExcitationList> exList_pj = interaction.exList_p2;
            NDetp2 = interaction.NDetp2;
            std::string twoparticleInts = interaction.twoBodyIntsString;
            CASHelper<MatsT,IntsT>::buildFullHTwoParticle(*pimcwfn,tmpH,exList_pi,exList_pj,twoparticleInts);
            MatsT chargeproduct = (MatsT)interaction.chargeproduct;
            // Add this block to the wavefunction
            CASHelper<MatsT,IntsT>::addBlockToMatrixBlockDiagonal(NDetp1*NDetp2,NDet,chargeproduct,tmpH,fullH);

            // Shift the second index to the back of queue
            shiftSecondIndexFullH(NDetp1,NDetp2,NDetStashed,NDet,fullH);

        }
       
        shiftFrontIndexToBackFullH(NDetp1,NDet,fullH);
        NDetStashed *= NDetp1;

    }
    CQMemManager::get().free(tmpH);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::buildDiagH(MCWaveFunction<MatsT,IntsT> & mcwfn,MatsT * diagH)
  {

    MULTIPARTICLECI_LOOP_INIT();

    MatsT * SCR = CQMemManager::get().template malloc<MatsT>(NDet);
    size_t NDetStashed = 1;
    std::string pilabel, pjlabel;
    size_t NDetp2;
    std::fill_n(diagH,NDet,MatsT(0.0));    

    for(size_t pi = 0; pi < CIOrder.size(); pi++)
    {
        std::fill_n(SCR,NDet,MatsT(0.0));    
        
        pilabel = CIOrder[pi];
        // Unpack the 1 particle objects
        auto piCIBuilder = mpmcwfn->getOneParticleCIHelper(pilabel);
        std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pimcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(piCIBuilder.refwfn);
        std::shared_ptr<const ExcitationList> exList_pi = piCIBuilder.exList;
        std::string onePIntsStr = piCIBuilder.oneBodyIntsStr;
        std::string twoPIntsStr = piCIBuilder.twoBodyIntsStr;
        size_t NDetp1 = piCIBuilder.NDet;

        CASHelper<MatsT,IntsT>::buildDiagHOneParticle(*pimcwfn,SCR,exList_pi,onePIntsStr,twoPIntsStr);

        // Add this diagonal to the full Hamiltonian
        CASHelper<MatsT,IntsT>::addVecToBlockedVector(NDetp1,NDet,MatsT(1.0),SCR,diagH);

        // Now begin building the pairwise interactions
        for(const auto & interaction : mpmcwfn->getTwoParticleCIHelper(pilabel))
        {
            std::fill_n(SCR,NDet,MatsT(0.0));
            pjlabel = interaction.labelp2;
            std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pjmcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(interaction.mcwfnp2);
            std::shared_ptr<const ExcitationList> exList_pj = interaction.exList_p2;
            size_t nActpi = interaction.nCorrp1;
            size_t nActpj = interaction.nCorrp2;
            NDetp2 = interaction.NDetp2;
            std::string twoparticleInts = interaction.twoBodyIntsString;
            CASHelper<MatsT,IntsT>::buildDiagHTwoParticle(*pimcwfn,SCR,nActpi,nActpj,exList_pi,exList_pj,twoparticleInts);
            MatsT chargeproduct = (MatsT)interaction.chargeproduct;
            // Add this block to the wavefunction
            CASHelper<MatsT,IntsT>::addVecToBlockedVector(NDetp1*NDetp2,NDet,chargeproduct,SCR,diagH);

            // Shift the second index to the back of queue
            shiftSecondIndex(NDetp1,NDetp2,NDetStashed,NDet,diagH);
        }

        shiftFrontIndexToBack(NDetp1,NDet,diagH);
        NDetStashed *= NDetp1;
    }

    CQMemManager::get().free(SCR);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::buildSigma(MCWaveFunction<MatsT,IntsT> & mcwfn, size_t nVec, MatsT * C, MatsT * Sigma)
  {
    MULTIPARTICLECI_LOOP_INIT();

    std::fill_n(Sigma,NDet*nVec,MatsT(0.0));

    size_t NDetStashed = 1;
    std::string pilabel, pjlabel;
    size_t NDetp2;

    for(size_t pi = 0; pi < CIOrder.size(); pi++)
    {
        pilabel = CIOrder[pi];

        auto piCIBuilder = mpmcwfn->getOneParticleCIHelper(pilabel);
        std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pimcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(piCIBuilder.refwfn);
        std::shared_ptr<const ExcitationList> exList_pi = piCIBuilder.exList;
        std::string onePIntsStr = piCIBuilder.oneBodyIntsStr;
        std::string twoPIntsStr = piCIBuilder.twoBodyIntsStr;
        size_t NDetp1 = piCIBuilder.NDet;

        CASHelper<MatsT,IntsT>::buildSigmaOneParticle(*pimcwfn,C,Sigma,nVec,NDet/NDetp1,exList_pi,onePIntsStr,twoPIntsStr);

        for(const auto & interaction : mpmcwfn->getTwoParticleCIHelper(pilabel))
        {
            pjlabel = interaction.labelp2;
            std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pjmcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(interaction.mcwfnp2);
            std::shared_ptr<const ExcitationList> exList_pj = interaction.exList_p2;
            NDetp2 = interaction.NDetp2;
            std::string twoparticleInts = interaction.twoBodyIntsString;

            CASHelper<MatsT,IntsT>::buildSigmaTwoParticle(*pimcwfn,C,Sigma,nVec,NDet/(NDetp1*NDetp2),exList_pi,exList_pj,twoparticleInts,interaction.chargeproduct);
            
            // Shift the second index to the back of queue
            for(size_t iVec = 0; iVec < nVec; iVec++)
            {
                shiftSecondIndex(NDetp1,NDetp2,NDetStashed,NDet,C+iVec*NDet);
                shiftSecondIndex(NDetp1,NDetp2,NDetStashed,NDet,Sigma+iVec*NDet);
            }
        }
        for(size_t iVec = 0; iVec < nVec; iVec++)
        {
            shiftFrontIndexToBack(NDetp1,NDet,C+iVec*NDet);
            shiftFrontIndexToBack(NDetp1,NDet,Sigma+iVec*NDet);
        }
        NDetStashed *= NDetp1;

    }
  }

  
  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::shiftFrontIndexToBackFullH(size_t NDetFront, size_t NDet, MatsT * Mat)
  {
    MatsT * tmp = Mat;
    // First do the L vectors
    for(size_t R = 0; R < NDet; R++, tmp+=NDet)
        shiftFrontIndexToBack(NDetFront,NDet,tmp);
    // Transpose the full matrix to handle R's
    IMatCopy('T',NDet,NDet,MatsT(1.0),Mat,NDet,NDet);
    // Do the R vectors
    tmp = Mat;
    for(size_t L = 0; L < NDet; L++, tmp+=NDet)
        shiftFrontIndexToBack(NDetFront,NDet,tmp);
    // Transpose L's back to front
    IMatCopy('T',NDet,NDet,MatsT(1.0),Mat,NDet,NDet);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::shiftFrontIndexToBack(size_t NDetFront, size_t NDet, MatsT * vec)
  {
    CASHelper<MatsT,IntsT>::transposeVectors(1,NDetFront,NDet/NDetFront,vec);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::shiftSecondIndexFullH(size_t NDetFront, size_t NDetSecond, size_t NDetAux, size_t NDet, MatsT * Mat)
  {
    MatsT * tmp = Mat;
    // First shift the L second index
    for(size_t R = 0; R < NDet; R++, tmp+=NDet)
    {
        shiftSecondIndex(NDetFront,NDetSecond,NDetAux,NDet,tmp);
    }
    // Transpose the full matrix to handle R's
    IMatCopy('T',NDet,NDet,MatsT(1.0),Mat,NDet,NDet);
    // Do the R vectors
    tmp = Mat;
    for(size_t L = 0; L < NDet; L++, tmp+=NDet)
        shiftSecondIndex(NDetFront,NDetSecond,NDetAux,NDet,tmp);
    // Transpose L's back to front
    IMatCopy('T',NDet,NDet,MatsT(1.0),Mat,NDet,NDet);

  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::shiftSecondIndex(size_t NDetFront, size_t NDetSecond, size_t NDetAux, size_t NDet, MatsT * vec)
  {
    // To shift the second index, we take the follow set of transposes:
    // |1,2,3,4,...,N> -> |2,3,4,...,N,1> -> |3,4,...N,2,1> -> |1,3,4,...,N,2>
    // We also do this loop over nAuxDet which is the number of indicies which 
    // are stashed at the back of the list of indices, e.g.
    // |3,4,...N,1,2> will need to do this over indices 1 & 2, so NDetAux should
    // be NDet1 * NDet2

    size_t NActiveDet = NDet/NDetAux;
    MatsT * auxvec = vec;
    for(size_t nAux = 0; nAux < NDetAux; nAux++, auxvec += NActiveDet)
    {
        CASHelper<MatsT,IntsT>::transposeVectors(1,NDetFront,NActiveDet/NDetFront,auxvec);
        MatsT * tmp = auxvec;
        for(size_t i = 0; i < NDetFront; i++, tmp += NActiveDet/NDetFront)
        {
            CASHelper<MatsT,IntsT>::transposeVectors(1,NDetSecond,NActiveDet/(NDetFront*NDetSecond),tmp);
        }
        CASHelper<MatsT,IntsT>::transposeVectors(1,NActiveDet/NDetFront,NDetFront,auxvec);
    }
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::computeOneRDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * C, std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> rdms)
  {
    computeTDM(mcwfn,C,C,rdms);
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::computeTDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT* Cm, MatsT* Cn, std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> tdms)
  {
    MULTIPARTICLECI_LOOP_INIT();
    std::unordered_map<std::string,size_t> rdmorder = mpmcwfn->getRDMOrder();
    for(auto & tdm : tdms) tdm.get().clear();
    for(size_t i = 0; i < CIOrder.size(); i++)
    {
        std::string pilabel = CIOrder[i];
        auto piCIBuilder = mpmcwfn->getOneParticleCIHelper(pilabel);
        std::shared_ptr<const ExcitationList> exList_pi = piCIBuilder.exList;
        size_t NDetpi = piCIBuilder.NDet;

        size_t whichTDM = rdmorder.at(pilabel);
        cqmatrix::Matrix<MatsT> tempTDM(tdms[whichTDM].get().nRows());
        tempTDM.clear();
        
        CASHelper<MatsT,IntsT>::computeTDM(mcwfn,Cm,Cn,NDet/NDetpi,exList_pi,tempTDM);
        tdms[whichTDM].get() += tempTDM;
        shiftFrontIndexToBack(NDetpi,NDet,Cm);
        if(Cm!=Cn)shiftFrontIndexToBack(NDetpi,NDet,Cn);
    }
 
  }

  template <typename MatsT, typename IntsT>
  void MultiParticleCASCI<MatsT,IntsT>::computeTwoRDM(MCWaveFunction<MatsT,IntsT> & mcwfn, MatsT * C, std::unordered_map<std::string,std::unordered_map<std::string,InCore4indexTPI<MatsT>>> twordms)
  {
    CErr("TwoRDM's for MultiParticleCASCI NYI!");
  }

}; // namespace ChronusQ