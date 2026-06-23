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

#include <chronusq_sys.hpp>

namespace ChronusQ{

  // Forward declare SingleSlater
  template <typename MatsT, typename IntsT>
  class SingleSlater;

  /*
  * A class which is used for construction of the SCF density
  * 
  * The SCF density can be written C^t[:occ] * C[:occ]
  * This however can be written C^t D[:occ] C
  * where D is the identity matrix with entries of 1 on the
  * diagaonal for occupied orbitals
  * 
  * In general, other flavors of SCF (e.g. MOM-Delta SCF, Thermal SCF)
  * can be achieve analogously by simple changes to the form of D 
  * 
  * This adds a small additional cost to the formation of the density,
  * for a conventional SCF calculation (since it is now 2 full NBxNB
  * matrix multiplications as opposed to one 1 NBxnocc), but this 
  * is typically not a bottleneck for SCF.
  */
  template <typename MatsT, typename IntsT>
  class RDMBuilderBase
  {
    public:
      virtual void buildRDM(SingleSlater<MatsT,IntsT> & ss) = 0;
  };

  /*
  * The default RDM builder builds the aufbau occupation matrix
  */
  template <typename MatsT, typename IntsT>
  class AufbauRDMBuilder : public RDMBuilderBase<MatsT,IntsT>
  {
    public:
    void buildRDM(SingleSlater<MatsT,IntsT>& ss);
 
  };

  /*
  * Maximum Overlap Method
  *
  * Computes the overlap of orbitals from the last SCF iteration
  * to current and selects the orbitals which best match those
  * from the previous iteration 
  */
  template <typename MatsT, typename IntsT>
  class MOMRDMBuilder : public RDMBuilderBase<MatsT,IntsT>
  {
    bool init;
    size_t count = 0;
    std::vector<cqmatrix::Matrix<MatsT>> prev_mo;
    public:
    void buildRDM(SingleSlater<MatsT,IntsT>& ss);
    MOMRDMBuilder(SingleSlater<MatsT,IntsT>&ss)
    {
      init = true;
      for(size_t i = 0; i < ss.mo.size(); i++)
        prev_mo.push_back(ss.mo[i]);
    };
  };

  /*
  * Floating Occupation Number SCF
  */
  template <typename MatsT, typename IntsT>
  class FON : public RDMBuilderBase<MatsT,IntsT>
  {
    public:
    void buildRDM(SingleSlater<MatsT,IntsT>& ss)
      {CErr("Floating Occupation Number SCF NYI!");};
 
  };


  /*
  * NEOStateAveragedRDMBuilder
  *
  * Builds an RDM which equally weights some user number
  * of input orbitals
  *
  * Primarily used for NEO: when there is only one particle
  * every orbital can be thought of as an independent
  * state and essentially SCF can be thought of as
  * doing a SA-CASSCF calculation
  */
  template <typename MatsT, typename IntsT>
  class NEOStateAveragedRDMBuilder : public RDMBuilderBase<MatsT,IntsT>
  {
    public:
    void buildRDM(SingleSlater<MatsT,IntsT>& ss);

    NEOStateAveragedRDMBuilder() = delete;
    NEOStateAveragedRDMBuilder(size_t NS) : NStates(NS) {};

    size_t NStates;
  };

  /*
  * Thermal
  *
  * Builds an RDM based on the thermal occupation of 
  * SCF orbitals
  */
  template <typename MatsT, typename IntsT>
  class Thermal : public RDMBuilderBase<MatsT,IntsT>
  {
    public:
    void buildRDM(SingleSlater<MatsT,IntsT>& ss)
      {CErr("Thermal SCF NYI!");};
 
  };


  // Implementation functions

  template <typename MatsT, typename IntsT>
  void AufbauRDMBuilder<MatsT,IntsT>::buildRDM(SingleSlater<MatsT,IntsT> & ss)
  {
    size_t nC = ss.nC;
    size_t NB = ss.nAlphaOrbital() * nC;

    // Only need to build the oneRDM once (it will never change)
    if(!ss.oneRDM) ss.oneRDM = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
    else return;
    
    ss.oneRDM->clear();

    if(nC == 1)
    {
      for(size_t ii = 0; ii < ss.nOA; ii++) ss.oneRDM->pointer()[ii+ii*NB] = 1.0;

      if(!ss.iCS)
      {
        if(!ss.oneRDMB)
          ss.oneRDMB = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
        ss.oneRDMB->clear();
        for(size_t ii = 0; ii < ss.nOB; ii++) ss.oneRDMB->pointer()[ii+ii*NB] = 1.0;
      }

    }
    else
    {
      for(size_t ii = 0; ii < ss.nO; ii++) ss.oneRDM->pointer()[ii+ii*NB] = 1.0;
    }
  }

  template <typename MatsT, typename IntsT>
  void NEOStateAveragedRDMBuilder<MatsT,IntsT>::buildRDM(SingleSlater<MatsT,IntsT> & ss)
  {
    size_t nC = ss.nC;
    size_t NB = ss.nAlphaOrbital() * nC;

    // Only need to build the oneRDM once (it will never change)
    if(!ss.oneRDM) ss.oneRDM = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
    else return;

    ss.oneRDM->clear();

    if(nC == 1)
    {
      for(size_t ii = 0; ii < NStates; ii++) ss.oneRDM->pointer()[ii+ii*NB] = 1.0/NStates;

      if(!ss.iCS)
      {
        if(!ss.oneRDMB)
          ss.oneRDMB = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
        ss.oneRDMB->clear();
      }

    }
    else
    {
      for(size_t ii = 0; ii < NStates; ii++) ss.oneRDM->pointer()[ii+ii*NB] = 1.0/NStates;
    }

  }

  template <typename MatsT, typename IntsT>
  void MOMRDMBuilder<MatsT,IntsT>::buildRDM(SingleSlater<MatsT,IntsT> & ss)
  {

    // Follows original MOM algorithm, see: 
    // https://doi.org/10.1021/jp801738f
    // Other algorithms might be better, see:
    // https://doi.org/10.1002/jcc.26797

    std::array<SpinType,2> spinmap{SpinType::isAlpha,SpinType::isBeta};
    size_t NB  = ss.nAlphaOrbital() * ss.nC;
    size_t nC = ss.nC;
    size_t offset = nC == 4 ? NB * NB / 2 : 0;

    // Tested:  UHF, 2CHF, GHF
    // Untested: RHF, ROHF
    // Needs further development with overlap: 4C
    if(ss.nC == 4) CErr("MOM with 4C NYI!");
    if(ss.nC == 1 && ss.iCS) CErr("MOM Untested with RHF");
    // Note ROHF errors out before this call to have access to the fock builder

    // SMG 05/26/26
    // This code is pretty fagile in it relies on the call to the constructor
    // and subsequent calls to buildRDM.
    // Implemented this way for future ROHF development
    if(init)
    {
      init = false;

      for(size_t i = 0; i < prev_mo.size(); i++)
        std::copy_n(prev_mo[i].pointer(),NB*NB,ss.mo[i].pointer());
      
      // Build now the Aufbau RDM with the rearranged orbitals
      AufbauRDMBuilder<MatsT,IntsT> subbuilder;
      subbuilder.buildRDM(ss);

      return;
    }

    // The way we are going to do Delta-SCF is by swapping orbitals and then calling the Aufbau builder
    std::vector<std::vector<std::pair<size_t,size_t>>> swaps;

    cqmatrix::Matrix<MatsT> S(NB); 
    S.clear();
    //prettyPrintSmart(std::cout,"OVERLAPSUBMAT",ss.aoints_->overlap->pointer(),NB/nC,NB/nC,NB/nC);

    if(nC == 1) S = ss.aoints_->overlap->matrix();
    else if (nC == 2)
    {
      size_t NBh = NB/2;
      MatAdd('N','N',NBh,NBh,MatsT(1.0),ss.aoints_->overlap->pointer(),NBh,MatsT(0.0),S.pointer(),NB,S.pointer(),NB);
      MatAdd('N','N',NBh,NBh,MatsT(1.0),ss.aoints_->overlap->pointer(),NBh,MatsT(0.0),S.pointer()+NB*NBh+NBh,NB,S.pointer()+NB*NBh+NBh,NB);
    }

    //prettyPrintSmart(std::cout,"OVERLAPSUPERMAT",S.pointer(),NB,NB,NB);

    // At most a system can have alpha and beta orbitals to pick from
    std::array<size_t,2> NToPick{0,0};
    NToPick[0] = nC == 1 ? ss.nOA : ss.nO;
    NToPick[1] = nC == 1 ? ss.nOB : 0;

    // Storage for calculating overlaps
    cqmatrix::Matrix<MatsT> SCR(NB),ovl(NB);
    MatsT * ovlp_metric = CQMemManager::get().malloc<MatsT>(NB);

    for(size_t i = 0; i < ss.mo.size(); i++)
    {
      SCR.clear();
      ovl.clear();
      std::fill_n(ovlp_metric,NB,MatsT(0.0));

      // This is just one metric, other metrics may be useful!
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans, blas::Op::NoTrans, NToPick[i], NB, NB, MatsT(1.), prev_mo[i].pointer()+offset, NB,
        S.pointer(), NB, MatsT(0.), SCR.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), SCR.pointer(), NB,
          ss.mo[i].pointer()+offset, NB, MatsT(0.), ovl.pointer(), NB);

      for(size_t j = 0; j < NB; j++)
        ovlp_metric[j] = blas::dot(NToPick[i],ovl.pointer()+j*NB,1,ovl.pointer()+j*NB,1);

      // Get the indices of the largest overlaps
      std::vector<size_t> idx(NB);
      std::iota(idx.begin(),idx.end(),0);
      std::stable_sort(idx.begin(),idx.end(),[&ovlp_metric](size_t i1, size_t i2){return std::abs(ovlp_metric[i1])>std::abs(ovlp_metric[i2]);});
      
      // idx will be sorted based on the overlap metric, so we just need to figure out which swaps from this list
      std::vector<std::pair<size_t,size_t>> this_spin_swaps;
      std::vector<size_t> to_swap_in;
      std::vector<size_t> to_swap_out;
      for(size_t n = 0; n < NToPick[i]; n++)
      {
        if(std::distance(idx.begin(),std::find(idx.begin(),idx.end(),n)) >= NToPick[i])
          to_swap_out.push_back(n);
        if(idx[n] >= NToPick[i])
          to_swap_in.push_back(idx[n]);
      }

      //if(!to_swap_out.size()) break;
      if(to_swap_out.size()!=to_swap_in.size()) CErr("Issue in MOM!");

      for(size_t n = 0; n < to_swap_out.size(); n++)
        this_spin_swaps.push_back(std::make_pair(to_swap_out[n]+1,to_swap_in[n]+1));

      swaps.push_back(this_spin_swaps);
      ss.swapMOs(swaps,spinmap[i]);
    }

    // Save the previous orbitals, this can change too
    for(size_t i = 0; i < ss.mo.size(); i++)
      prev_mo[i] = ss.mo[i];

    // Build now the Aufbau RDM with the rearranged orbitals
    AufbauRDMBuilder<MatsT,IntsT> subbuilder;
    subbuilder.buildRDM(ss);

  }


}; // namespace ChronusQ
