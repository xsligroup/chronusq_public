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

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/factorization.hpp>

#include <util/threads.hpp>


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator<HartreeFock<MatsT, IntsT>>::formLinearTrans_direct_impl(
    MPI_Comm c, RC_coll<U> x){ 
  

    SingleSlater<MatsT,IntsT> &ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);

    bool doA   = this->tdaOp == PP_A;
    bool doOcc = not this->genSettings.doTDA or this->tdaOp == PP_C;

    const size_t tdOffSet = this->genSettings.doTDA ? 0 : this->aDim_;
    const size_t N        = this->nSingleDim_ ;  
    const size_t chunk    = 600;



#if 1
/****************************/

  
    size_t NB = ss.nAlphaOrbital();
    size_t NO = ss.nOA;

    MatsT* MMat = CQMemManager::get().malloc<MatsT>(NB*NB);
    MatsT* JMMat = CQMemManager::get().malloc<MatsT>(NB*NB);
    MatsT* KMMat = CQMemManager::get().malloc<MatsT>(NB*NB);

    MatsT* MO = ss.mo[0].pointer();

    size_t indmo1 = NO;
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),MO + NO*NB,NB,MO + NO*NB,NB,MatsT(0.),MMat,NB);

    std::fill_n(JMMat,NB*NB,0.);
    std::fill_n(KMMat,NB*NB,0.);

    std::vector<TwoBodyContraction<MatsT>> in = 
      { 
        { MMat,JMMat,true, COULOMB  },
        { MMat,KMMat,true, EXCHANGE }
      };

    ss.TPI->twoBodyContract(c,in);

    MatAdd('N','N',NB,NB,MatsT(1.),JMMat, NB, MatsT(-0.5), KMMat, NB, JMMat, NB);

    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),JMMat,NB,MO   ,NB,MatsT(0.),KMMat,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),MO   ,NB,KMMat,NB,MatsT(0.),MMat ,NB);

/****************************/
#endif

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


          // SCF* reference
          if( this->doStarRef ) {

            const size_t NV_t  = this->getNV1(this->doStarRef); 
            const size_t NV2_t = this->getNV2(this->doStarRef);
            const size_t NO_t  = this->getNO1(this->doStarRef); 
            const size_t NO2_t = this->getNO2(this->doStarRef); 

            const bool doLT = (this->ref_->nC == 2) or (this->spinSepProp != PP_AB);
            auto bmax = [&](size_t a){ return doLT ? a : NV2_t; };
            auto jmax = [&](size_t i){ return doLT ? i : NO2_t; };

            for(auto iVec  = 0; iVec < nDo; iVec++) {

              auto * HV_cc = HV_c + iVec*N;
              auto * V_cc  = V_c  + iVec*N;

              for(auto a = 0ul, ab = 0ul; a < NV_t;    a++      )
              for(auto b = 0ul          ; b < bmax(a); b++, ab++){

                const size_t A = a + NO_t;
                const size_t B = b + NO_t;

                for(auto c = 0ul, cd = 0ul; c < NV_t;    c++      )
                for(auto d = 0ul          ; d < bmax(c); d++, cd++){

                  const size_t C = c + NO_t;
                  const size_t D = d + NO_t;

                  if(a == c) HV_cc[ab] -= V_cc[cd] * MMat[B + D*NB];
                  if(b == d) HV_cc[ab] -= V_cc[cd] * MMat[A + C*NB];
                  if(a == d) HV_cc[ab] += V_cc[cd] * MMat[B + C*NB];
                  if(b == c) HV_cc[ab] += V_cc[cd] * MMat[A + D*NB];

                }


              }

              for(size_t i = 0ul, ij = this->aDim_; i < NO_t;    i++      )
              for(size_t j = 0ul          ; j < jmax(i); j++, ij++){

                for(size_t k = 0ul, kl = this->aDim_; k < NO_t;    k++      )
                for(size_t l = 0ul          ; l < jmax(k); l++, kl++){

                  if(i == k) HV_cc[ij] += V_cc[kl] * MMat[j + l*NB];
                  if(j == l) HV_cc[ij] += V_cc[kl] * MMat[i + k*NB];
                  if(i == l) HV_cc[ij] -= V_cc[kl] * MMat[j + k*NB];
                  if(j == k) HV_cc[ij] -= V_cc[kl] * MMat[i + l*NB];

                }


              }





            }

          }
        }


        // Free up transformation memory
        CQMemManager::get().free(cList[0].X);

      } // loop over vectors

    } // loop over groups of vectors

    CQMemManager::get().free(MMat,JMMat,KMMat);

  };
  

  template <>
  void ParticleParticlePropagator<HartreeFock<double,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    formLinearTrans_direct_impl(c,x);

  };

  template <>
  void ParticleParticlePropagator<HartreeFock<dcomplex,double>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };

  template <>
  void ParticleParticlePropagator<HartreeFock<dcomplex,dcomplex>>::formLinearTrans_direct(
    MPI_Comm c, RC_coll<double> x){ 
  
    CErr("SOMETHING HAS GONE HORRIBLY WRONG: formLinearTrans_direct");

  };





} // namespace ChronusQ

