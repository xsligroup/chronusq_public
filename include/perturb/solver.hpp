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

#include <perturb.hpp>
#include <lapack.hh>

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

#define _DEBUG_PERTURB_SOLVER

namespace ChronusQ {

  /**
   *  \brief solve the linear equation
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::solve(MatsT * Sol, size_t i) {

    ProgramTimer::tick("Solve Linear System");

#ifdef _DEBUG_PERTURB_SOLVER
    std::cout<<"solve linear system"<<std::endl;
#endif

    // whether to form full matrix or not
    if( PTopts.doFull ) {
      buildFull(*LHS_, E0_[i]);
      if( not PTopts.doIter ) {
        //std::cout<<"PT2 lapack solver."<<std::endl;
        
        int64_t* IPIV = CQMemManager::get().malloc<int64_t>(SDsize);
        lapack::gesv(SDsize,1,LHS_->pointer(),SDsize,IPIV,Sol,SDsize);
        CQMemManager::get().free(IPIV);

      }
    }

    if( PTopts.doIter ) iterSolve(Sol, E0_[i]);

    ProgramTimer::tock("Solve Linear System");

  } // PERTURB::solve

  /**
   *  \brief Preconditioner for GMRES
   *            using diagonal of LHS inverse
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::PT2preCond(size_t nVec, MatsT shift, SolverVectors<MatsT> &V,
                            SolverVectors<MatsT> &AV, MatsT * currdiag) {

    double small = 1e-12;
    double P;
    for (auto i = 0ul; i < nVec; i++) {
      auto Vi = tryGetRawVectorsPointer(V) + i * SDsize;
      auto AVi = tryGetRawVectorsPointer(AV) + i * SDsize;
      #pragma omp parallel for schedule(static) default(shared)
      for (auto j = 0ul; j < SDsize; j++)
        if (std::abs(currdiag[j]) > small) AVi[j] = Vi[j] / currdiag[j];
    }
  } // PERTURB::PT2preCond


  /**
   *  \brief Solve the linear equation iteratively
   *            using GMRES.
   *            LHS * P = RHS
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::iterSolve(MatsT * Sol, MatsT E0) {

    std::cout<<"LL PT2 itersolve."<<std::endl;
    bool isRoot = MPIRank(this->comm) == 0;

    //build diagLHS for preconditioner
    MatsT  * currdiag = CQMemManager::get().malloc<MatsT>(SDsize);
    std::copy_n(diagLHS, SDsize, currdiag);
    for (auto I = 0ul; I < SDsize; I++) currdiag[I] += E0 - computeShift(I, E0);

    // Simple preconditioner
    typename GMRES<MatsT>::Shift_t pc = [&](size_t nVec, MatsT shift, SolverVectors<MatsT> &V,
                                            SolverVectors<MatsT> &AV) {

#ifdef CQ_ENABLE_MPI
      // disable MPI for now
      ROOT_ONLY(this->comm);
#endif
      this->PT2preCond(nVec, shift, V, AV, currdiag);

    };

    typename GMRES<MatsT>::LinearTrans_t lt = [&](size_t nVec, SolverVectors<MatsT> &V,
                                                    SolverVectors<MatsT> &AV) {

#ifdef CQ_ENABLE_MPI
      // disable MPI for now
      ROOT_ONLY(this->comm);
#endif
      size_t N = this->SDsize;
      if (PTopts.doFull)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,
                   N,MatsT(1.),LHS_->pointer(),N,tryGetRawVectorsPointer(V),N,MatsT(0.),tryGetRawVectorsPointer(AV),N);
      else this->buildAV(tryGetRawVectorsPointer(V), tryGetRawVectorsPointer(AV), E0);

    };


    GMRES<MatsT> gmres(this->comm,this->SDsize,
      PTopts.maxIter,PTopts.convCrit,lt,pc);

    gmres.setRHS(1,Sol,this->SDsize);
    gmres.setShifts(this->shifts.size(), &this->shifts[0]);

    gmres.rhsBS = 1;
    gmres.shiftBS = 1;

    std::cout<<"itersolve run."<<std::endl;
    gmres.run();

    if( isRoot )
      std::copy_n(tryGetRawVectorsPointer(*gmres.getSol()), this->SDsize, Sol);

    std::cout<<"itersolve done."<<std::endl;

    CQMemManager::get().free(currdiag);
  } // PERTRUB::solveLinear

} // namespace ChronusQ



