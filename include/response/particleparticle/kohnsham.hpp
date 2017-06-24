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

#include <response/particleparticle.hpp>
#include <response/particleparticle/singleslater_helper.hpp>

#include <particleintegrals/twopints.hpp>

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/factorization.hpp>

#include <util/threads.hpp>


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator<KohnSham<MatsT,IntsT>>::formLinearTrans_direct_impl(
    MPI_Comm c, RC_coll<U> x){ 
  

    SingleSlater<MatsT,IntsT> &ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);

    bool doA   = this->tdaOp == PP_A;
    bool doOcc = not this->genSettings.doTDA or this->tdaOp == PP_C;

    const size_t tdOffSet = this->genSettings.doTDA ? 0 : this->aDim_;
    const size_t N        = this->nSingleDim_ ;  
    const size_t chunk    = 600;

    std::shared_ptr<TPIContractions<U,IntsT>> TPI =
        TPIContractions<MatsT,IntsT>::template convert<U>(ss.TPI);

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
        scatter = MPIAnyOf(scatter,c);
#endif

        // Transform pp vector MO -> AO
        auto cList = this->genSettings.doTDA ?
          this->template ppTransitionVecMO2AO<U>(c,scatter,nDo,N,V_c) : 
          this->template ppTransitionVecMO2AO<U>(c,scatter,nDo,N,V_c,
            V_c + tdOffSet);

        TPI->twoBodyContract(c,cList); // form G[V]




        // Only finish transformation on root process
        if( MPIRank(c) == 0 ) {

          // Transform pp vector AO -> MO
          if( this->genSettings.doTDA )
            this->ppTransitionVecAO2MO(nDo,N,cList,HV_c);
          else
            this->ppTransitionVecAO2MO(nDo,N,cList,HV_c,HV_c + tdOffSet);




          // Scale by diagonals
          if( not this->genSettings.doTDA or doA )
            this->ppEpsilonScale(true,false,false,nDo,N,V_c,HV_c);

          if( doOcc )
            this->ppEpsilonScale(true,false,true,nDo,N,V_c + tdOffSet, 
                HV_c + tdOffSet);

        }


        // Free up transformation memory
        CQMemManager::get().free(cList[0].X);

      } // loop over vectors

    } // loop over groups of vectors


  };
  

  template <>
  void ParticleParticlePropagator<KohnSham<double,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    formLinearTrans_direct_impl(c,x);

  };

  template <>
  void ParticleParticlePropagator<KohnSham<dcomplex,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };

  template <>
  void ParticleParticlePropagator<KohnSham<dcomplex,dcomplex>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };





} // namespace ChronusQ

