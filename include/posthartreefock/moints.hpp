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

#include <posthartreefock.hpp>
#include <mointstransformer/impl.hpp>
#include <particleintegrals/dasints.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>
#include <util/matout.hpp>

// #define DEBUG_POSTHF_MOINTS

namespace ChronusQ {

template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::setMORanges() {
  mointsTF->setMORanges(dynamic_cast<const PostHartreeFockBase&>(*this));
} // PostHartreeFock::setMORanges

template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::transformInts(EMPerturbation & pert,
  bool cacheHalfTransTPI) {
  
  // TODO: expand to RI
  // splitting comm for root only parts
  auto old_comm = this->comm;
  MPI_Comm new_comm;

#ifdef CQ_ENABLE_MPI
  int color = (MPIRank(old_comm) == 0); // Set color to 1 only for root process
  MPI_Comm_split(old_comm, color, 0, &new_comm); 
#endif

  // build object and allocate memory
  size_t nCorrO  = this->corrSpace.nCorrO;
  size_t nCoreO  = this->corrSpace.nInact + this->corrSpace.nFCore;

  /*
   * compute inactive core energy
   */

  double fc1C = (this->reference()->nC == 1) ? 2.0 : 1.0;
  ProgramTimer::tick("MOINTSTRANSFORM CORE ENERGY");
  if (MPIRank(this->comm) == 0) {
    MatsT * h1e_II  = CQMemManager::get().malloc<MatsT>(nCoreO);
    MatsT * GD_JJII = CQMemManager::get().malloc<MatsT>(nCoreO);
    
    mointsTF->transformHCore(pert, h1e_II, "II", true);
    // just change the communicator for the ss object not the mointsTF object
    mointsTF->ss_.comm = new_comm;
    mointsTF->transformGD(pert, 'I', GD_JJII, "JJ", true, true, "WithInactive-I");  
    // change it back
    mointsTF->ss_.comm = old_comm;
    // destroy new comm
    #ifdef CQ_ENABLE_MPI
      MPI_Comm_free(&new_comm);
    #endif
    
    MatsT ECore = 0.;
    // compute core enenrgy
    for (auto i = 0ul; i < nCoreO; i++) {
      ECore += h1e_II[i] + 0.5 * GD_JJII[i];
    }
  
    CQMemManager::get().free(h1e_II, GD_JJII);

    this->coreEnergy = std::real(ECore) * fc1C;
  }

  MPIBCast(this->coreEnergy, 0, this->comm); 
  ProgramTimer::tock("MOINTSTRANSFORM CORE ENERGY");

  /*
   * compute hCore and ERI in correlated space
   */ 
  OnePInts<MatsT> hCore_tu(nCorrO);
  InCore4indexTPI<MatsT> ERI_tuvw(nCorrO);
  
  // TODO: MPI-Parallel this
  ProgramTimer::tick("MOINTSTRANSFORM OPI TRANS");
  if (MPIRank(this->comm) == 0) {
    mointsTF->transformHCore(pert, hCore_tu.pointer(), "tu", false, 'I');
  }
  size_t nCorrO2 = nCorrO * nCorrO;
#ifdef CQ_ENABLE_MPI
  MPIBCast(hCore_tu.pointer(), nCorrO2, 0, this->comm);
#endif
  ProgramTimer::tock("MOINTSTRANSFORM OPI TRANS");


  ProgramTimer::tick("MOINTSTRANSFORM TPI TRANS");
  // With in-core AO integrals, use in-core SSFock-N6 transform (same path as the old
  // MCSCF code)
  if (std::dynamic_pointer_cast<InCoreTPI<IntsT>>(
        this->reference()->aoints_->TPI))
    mointsTF->transformTPI(pert, ERI_tuvw.pointer(), "tuvw", false);
  else
    mointsTF->directTransformTPI(pert, ERI_tuvw.pointer(), "tuvw");
  ProgramTimer::tock("MOINTSTRANSFORM TPI TRANS");

  // ERI_tuvw.output(std::cout, "TPI", true);
  
  // for diagonal elements 
  // can also form this using direct transformation
  // TODO: What is this for one component?
  OnePInts<MatsT> antiSymmetricERI_ttuu(nCorrO);
  DASOnePInts<MatsT> hCore_tt(nCorrO, 1ul); 
  #pragma omp parallel for schedule(static) default(shared)       
  for (auto u = 0ul; u < nCorrO; u++) {
    hCore_tt(u, 0) = hCore_tu(u, u);
    for (auto t = 0ul; t < nCorrO; t++) {
      antiSymmetricERI_ttuu(t, u) = ERI_tuvw(t,t,u,u) - ERI_tuvw(t,u,u,t);
    }
  }
  
  this->moints->addIntegral("hCore_Correlated_Space", 
    std::make_shared<OnePInts<MatsT>>(hCore_tu));
  this->moints->addIntegral("ERI_Correlated_Space",  
    std::make_shared<InCore4indexTPI<MatsT>>(ERI_tuvw));

  if (this->saveMOInts and MPIRank(this->comm) == 0 and this->savFile.exists()) {
    double inactiveEnergy = this->reference()->molecule().nucRepEnergy + this->coreEnergy;
    this->savFile.safeWriteData("MOINTS/INACTENERGY", &inactiveEnergy, {1});
    this->savFile.safeWriteData("MOINTS/ONEELEC", hCore_tu.pointer(), {nCorrO, nCorrO});
    this->savFile.safeWriteData("MOINTS/ERI", ERI_tuvw.pointer(),
                                {nCorrO, nCorrO, nCorrO, nCorrO});
  }

  this->moints->addIntegral("hCore_tt",  
    std::make_shared<DASOnePInts<MatsT>>(hCore_tt));
  this->moints->addIntegral("antiSymmetricERI_ttuu",  
    std::make_shared<OnePInts<MatsT>>(antiSymmetricERI_ttuu));

} // PostHartreeFock::transformInts

template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::prepareMOIntegrals(
    EMPerturbation & pert, const DeterminantFactory& detFactory,
    bool cacheHalfTransTPI, bool removeFullSpaceMOInts) {
  
  this->transformInts(pert, cacheHalfTransTPI); 
  
  ProgramTimer::tick("MOINTSTRANSFORM DAS INTS");
  std::vector<std::string> termVec;
  std::vector<size_t> span;
  size_t nt, nu, nw, nv, tOff, uOff, wOff, vOff;
  
  const auto& corrS = this->corrSpace;
  const auto& corrSpaceOff = corrS.nNegMO + corrS.nFCore + corrS.nInact;
  const auto& activeSpaces = detFactory.ketCategoricalSpace()->activeSpaces();
  const auto& ERI  = *(moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space"));
  
  for (const auto& term : detFactory.twoEExTerms()) {
    
    // std::cout << "* computing 2e term:" << term << std::endl; 

    parseTermSpan(term, termVec, span);
    
    nt = activeSpaces[span[0]].nOrbitals; 
    nu = activeSpaces[span[1]].nOrbitals; 
    nw = activeSpaces[span[2]].nOrbitals; 
    nv = activeSpaces[span[3]].nOrbitals; 
    tOff = activeSpaces[span[0]].MOOffset - corrSpaceOff; 
    uOff = activeSpaces[span[1]].MOOffset - corrSpaceOff; 
    wOff = activeSpaces[span[2]].MOOffset - corrSpaceOff; 
    vOff = activeSpaces[span[3]].MOOffset - corrSpaceOff; 
     
    // std::cout << "    - tOff = " << tOff << ", nt = " << nt << std::endl;
    // std::cout << "    - uOff = " << uOff << ", nu = " << nu << std::endl;
    // std::cout << "    - wOff = " << wOff << ", nw = " << nw << std::endl;
    // std::cout << "    - vOff = " << vOff << ", nv = " << nv << std::endl;

    DASTwoPInts<MatsT> ERI_sub(nt, nu, nw, nv);

    if (termVec.size() == 1) {
#pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto vv = 0ul; vv < nv; vv++)
      for (auto ww = 0ul; ww < nw; ww++) {
        auto v = vOff + vv;
        auto w = wOff + ww;
        for (auto u = uOff, uu = 0ul; uu < nu; u++, uu++)
        for (auto t = tOff, tt = 0ul; tt < nt; t++, tt++)      
          ERI_sub(tt, uu, ww, vv) = ERI(t, u, w, v);
      }
    } else if (termVec.size() > 1 and termVec[1] == "X") {
#pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto vv = 0ul; vv < nv; vv++)
      for (auto ww = 0ul; ww < nw; ww++) {
        auto v = vOff + vv;
        auto w = wOff + ww;
        for (auto u = uOff, uu = 0ul; uu < nu; u++, uu++)
        for (auto t = tOff, tt = 0ul; tt < nt; t++, tt++)      
          ERI_sub(tt, uu, ww, vv) = ERI(t, u, w, v) - ERI(t, v, w, u);
      }
    } else {
      CErr(term + " NYI"); 
    } 
    
    this->moints->addIntegral(term, std::make_shared<DASTwoPInts<MatsT>>(ERI_sub));
  } // twoEExTerms
  
  const auto& hCore = *(moints->template getIntegral<OnePInts,MatsT>("hCore_Correlated_Space"));

  for (const auto& term : detFactory.oneEExTerms()) {
    
    // std::cout << "* computing 1e term:" << term << std::endl; 
    
    parseTermSpan(term, termVec, span, "+RI[]");

    nt = activeSpaces[span[0]].nOrbitals; 
    nu = activeSpaces[span[1]].nOrbitals; 
    tOff = activeSpaces[span[0]].MOOffset - corrSpaceOff; 
    uOff = activeSpaces[span[1]].MOOffset - corrSpaceOff; 
    
    DASOnePInts<MatsT> h1e_sub(nt, nu);
    #pragma omp parallel for schedule(static) collapse(2) default(shared)       
    for (auto uu = 0ul; uu < nu; uu++)
    for (auto tt = 0ul; tt < nt; tt++) {
      auto t = tOff + tt;
      auto u = uOff + uu;
      h1e_sub(tt, uu) = hCore(t,u);
    }
    
    for (auto i = 1ul; i < termVec.size(); ++i) {
      const auto& ERI_sub  = *(moints->template getIntegral<DASTwoPInts,MatsT>(termVec[i]));
      #pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto u = 0ul; u < nu; u++)
      for (auto t = 0ul; t < nt; t++) {
        for (auto w = 0ul; w < ERI_sub.nBasis2(); ++w)
          h1e_sub(t, u) -= 0.5 * ERI_sub(t, w, w, u);
      }
    } 
    
    this->moints->addIntegral(term, std::make_shared<DASOnePInts<MatsT>>(h1e_sub));
  } // oneEExcitations 
  ProgramTimer::tock("MOINTSTRANSFORM DAS INTS");
  
  if (removeFullSpaceMOInts) {
    this->moints->erase("hCore_Correlated_Space"); 
    this->moints->erase("ERI_Correlated_Space");
  }
  return; 
} // prepareMOIntegrals 

} // namespace ChronusQ  
