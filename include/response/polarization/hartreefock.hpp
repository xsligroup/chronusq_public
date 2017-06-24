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

#include <response/polarization.hpp>

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/factorization.hpp>

#include <util/threads.hpp>
#include <util/timer.hpp>


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  template <typename U>
  void PolarizationPropagator<HartreeFock<MatsT, IntsT>>::formLinearTrans_direct_impl(
    MPI_Comm c, RC_coll<U> x, SINGLESLATER_POLAR_COPT op, 
    bool noTrans){ 
  

    if( op != FULL ) CErr("Direct + non-Full NYI");

    HartreeFock<MatsT, IntsT> &hf = dynamic_cast<HartreeFock<MatsT, IntsT>&>(*this->ref_);
		//additions to refactor Transforms and Scaling Functions

    const size_t N        = this->getNSingleDim(this->genSettings.doTDA) * (this->doReduced ? 2 : 1);  
    const size_t tdOffSet = N / 2;
    const size_t chunk    = 600;

    std::shared_ptr<TPIContractions<U,IntsT>> TPI =
        TPIContractions<MatsT,IntsT>::template convert<U>(hf.TPI);

    ProgramTimer::tick("Direct Hessian Contract");

    for(auto &X : x) {

      const size_t nVec = X.nVec;

      if( X.AX ) std::fill_n(X.AX,N*nVec,0.);

      for(size_t k = 0; k < nVec; k += chunk) {

        MPI_Barrier(c); // Sync MPI Processes at begining of each batch of 
                        // vectors

        const size_t nDo = std::min( chunk, nVec - k);

        auto *V_c  = X.X  + k*N;
        auto *HV_c = X.AX + k*N;


        bool scatter = not bool(X.X);
#ifdef CQ_ENABLE_MPI
        scatter = MPIAnyOf(scatter, c);
#endif


        // Transform ph vector MO -> AO
        auto cList = 
          this->template phTransitionVecMO2AO<U>(c, scatter, nDo, N, hf, hf,
            true, V_c, V_c + tdOffSet);

        TPI->twoBodyContract(c,cList,this->scfPert); // form G[V]

        // Only finish transformation on root process
        if( MPIRank(c) == 0 ) {

          // Transform ph vector AO -> MO
          this->phTransitionVecAO2MO(nDo,N,cList,hf,true,HV_c,HV_c + tdOffSet);

          // Scale by diagonals
          this->phEpsilonScale(true,false,nDo,N,hf,V_c,HV_c);
          this->phEpsilonScale(true,false,nDo,N,hf,V_c+tdOffSet,
            HV_c+tdOffSet);

        }

        // Free up transformation memory
        CQMemManager::get().free(cList[0].X);

      }


      if( X.AX and this->incMet and not this->doAPB_AMB )
        SetMat('N', N/2, nVec, U(-1.), X.AX + (N/2), N, X.AX + (N/2), N);

    } // loop over groups of vectors

    ProgramTimer::tock("Direct Hessian Contract");

  };
  

  template <>
  void PolarizationPropagator<HartreeFock<double,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x, SINGLESLATER_POLAR_COPT op, 
    bool noTrans){ 
  
    formLinearTrans_direct_impl(c,x,op,noTrans);

  };

  template <>
  void PolarizationPropagator<HartreeFock<dcomplex,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x, SINGLESLATER_POLAR_COPT op, 
    bool noTrans){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };


  template <>
  void PolarizationPropagator<HartreeFock<dcomplex,dcomplex>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x, SINGLESLATER_POLAR_COPT op, 
    bool noTrans){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };



} // namespace ChronusQ

