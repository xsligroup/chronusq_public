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

#include <singleslater.hpp>
#include <singleslater/hartreefock.hpp>
#include <singleslater/kohnsham.hpp>

#include <response/tbase.hpp>

namespace ChronusQ {

  // Forward decl of Reference templated ParticleParticlePropagator class
  template <class Reference> class ParticleParticlePropagator;


  struct SingleSlaterParticleParticleBase { 

    ParticleParticleProp_SpinSep  spinSepProp = PP_AB; // Default to AB
    ParticleParticleTDA           tdaOp       = PP_A;  // Default to A-TDA

    double                        mu          = 0;     // Chemical Potential
    
    bool                          doStarRef   = false; // Do "*" reference
    std::pair<long int, long int> starOrb;             // Orbitals for *-ref
  };


  template <typename MatsT, typename IntsT>
  class ParticleParticlePropagator< SingleSlater<MatsT,IntsT> > : 
    public ResponseRefBase< SingleSlater<MatsT,IntsT> >, 
    public SingleSlaterParticleParticleBase {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;

    protected:

      size_t aDim_  = 0;
      size_t cDim_  = 0;
      bool   doVir_ = false;
      bool   doOcc_ = false;


      inline size_t getNO1(bool doStar) const {

        size_t NO = (this->ref_->nC == 2)  ? this->ref_->nO  :
          (spinSepProp == PP_BB) ? this->ref_->nOB : this->ref_->nOA;

        if( doStar ) {

          bool starIsAlpha = starOrb.first > 0;

          if( this->ref_->nC == 2 ) NO -= 2;
          else if( spinSepProp == PP_BB ) {

            if( not starIsAlpha ) NO--;

          } else if( starIsAlpha ) NO--;

        }


        return NO;

      }


      inline size_t getNO2(bool doStar) const {

        size_t NO = (this->ref_->nC == 2)  ? this->ref_->nO  :
          (spinSepProp == PP_AA) ? this->ref_->nOA : this->ref_->nOB;

        if( doStar ) {

          bool starIsAlpha = starOrb.second > 0;

          if( this->ref_->nC == 2 ) NO -= 2;
          else if( spinSepProp == PP_AA ) {

            if( starIsAlpha ) NO--;

          } else if( not starIsAlpha ) NO--;

        }


        return NO;

      }


      inline size_t getNV1(bool doStar) const {

        size_t NV = (this->ref_->nC == 2)  ? this->ref_->nV  :
          (spinSepProp == PP_BB) ? this->ref_->nVB : this->ref_->nVA;

        if( doStar ) {

          bool starIsAlpha = starOrb.first > 0;

          if( this->ref_->nC == 2 ) NV += 2;
          else if( spinSepProp == PP_BB ) {

            if( not starIsAlpha ) NV++;

          } else if( starIsAlpha ) NV++;

        }


        return NV;

      }


      inline size_t getNV2(bool doStar) const {

        size_t NV = (this->ref_->nC == 2)  ? this->ref_->nV  :
          (spinSepProp == PP_AA) ? this->ref_->nVA : this->ref_->nVB;

        if( doStar ) {

          bool starIsAlpha = starOrb.second > 0;

          if( this->ref_->nC == 2 ) NV += 2;
          else if( spinSepProp == PP_AA ) {

            if( starIsAlpha ) NV++;

          } else if( not starIsAlpha ) NV++;

        }


        return NV;

      }



      inline size_t getADim() const {

        /*
        auto &ss = *this->ref_;
        size_t nVAVA_lt = ss.nVA * (ss.nVA-1) / 2;
        size_t nVBVB_lt = ss.nVB * (ss.nVB-1) / 2;
        size_t nVV_lt   = ss.nV  * (ss.nV-1)  / 2;

        size_t nVAVB    = ss.nVA * ss.nVB;

        if( ss.nC == 2 )                return nVV_lt;
        else if( spinSepProp == PP_AA ) return nVAVA_lt;
        else if( spinSepProp == PP_BB ) return nVBVB_lt;
        else                            return nVAVB;
        */

        const size_t NV1 = getNV1(doStarRef);
        const size_t NV2 = getNV2(doStarRef);

        const bool doLT = (this->ref_->nC == 2) or (spinSepProp != PP_AB);


        return doLT ? NV1*(NV1-1) / 2 : NV1*NV2;


      };

      inline size_t getCDim() const {

        /*
        auto &ss = *this->ref_;
        size_t nOAOA_lt = ss.nOA * (ss.nOA-1) / 2;
        size_t nOBOB_lt = ss.nOB * (ss.nOB-1) / 2;
        size_t nOO_lt   = ss.nO  * (ss.nO-1)  / 2;

        size_t nOAOB    = ss.nOA * ss.nOB;

        if( ss.nC == 2 )                return nOO_lt;
        else if( spinSepProp == PP_AA ) return nOAOA_lt;
        else if( spinSepProp == PP_BB ) return nOBOB_lt;
        else                            return nOAOB;
        */

        const size_t NO1 = getNO1(doStarRef);
        const size_t NO2 = getNO2(doStarRef);

        const bool doLT = (this->ref_->nC == 2) or (spinSepProp != PP_AB);


        return doLT ? NO1*(NO1-1) / 2 : NO1*NO2;


      };


    public:

      inline virtual size_t getNSingleDim(const bool doTDA = false) override {

        size_t N = 0;

        size_t ADIM = getADim();
        size_t CDIM = getCDim();

        if( doTDA ) N = (tdaOp == PP_A) ? ADIM : CDIM;
        else        N = ADIM + CDIM;

        assert(N != 0);

        return N;

      }



      ParticleParticlePropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<SingleSlater<MatsT,IntsT>> ref, 
          MatsT* fullMatrix = nullptr ) : 
        ResponseRefBase<SingleSlater<MatsT,IntsT>>(c,job,ref,fullMatrix),
        SingleSlaterParticleParticleBase() { 

          this->hasResGuess_ = false;

          this->PC_ = [&](size_t nVec, MatsT shift, SolverVectors<MatsT> &V, SolverVectors<MatsT> &AV) {

            preConditioner(nVec,shift,V,AV);

          };

          this->nSPC_ = [&](size_t nVec, SolverVectors<MatsT> &V, SolverVectors<MatsT> &AV) {

            preConditioner(nVec,MatsT(0.),V,AV);

          };

          this->cmplxPC_ = [&](size_t nVec, dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV) {

            preConditioner(nVec,shift,V,AV);

          };

        }

      ParticleParticlePropagator( const ParticleParticlePropagator &other ) : 
        ResponseRefBase<SingleSlater<MatsT,IntsT>>(
          dynamic_cast<const ResponseRefBase<SingleSlater<MatsT,IntsT>>&>(other)
        ), 
        SingleSlaterParticleParticleBase(
          dynamic_cast<const SingleSlaterParticleParticleBase&>(other)
        ) { }


      MatsT*                   formFullFromMemory() override;
      MatsT*                   formFullMatrix()  override;
      void                 formRHS() override;
      std::pair<size_t,MatsT*> formPropGrad(ResponseOperator) override;
      void                 configOptions() override;
      void                 eigVecNorm() override;
      void                 constructShifts() override;
      void                 postLinearSolve() override;
      void                 resGuess(size_t, SolverVectors<MatsT> &, size_t) override { CErr(); };


      // Internal implementations for direct linear transformation
      virtual void formLinearTrans_direct(MPI_Comm, RC_coll<double>)   = 0;
      virtual void formLinearTrans_direct(MPI_Comm, RC_coll<dcomplex>) = 0;

      // Incore linear transformation
      template <typename U>
      void formLinearTrans_incore_impl(RC_coll<U>);

      void formLinearTrans_incore(RC_coll<double> x);

      void formLinearTrans_incore(RC_coll<dcomplex> x) {

        formLinearTrans_incore_impl(x);

      }






      // Toggle between direct / incore
      template <typename U>
      void formLinearTrans_impl(RC_coll<U> x) {
       
        if( this->genSettings.formFullMat ) formLinearTrans_incore(x);
        else  formLinearTrans_direct(this->comm_,x);

      };

      // Interface to ResponseRBase double exposure
      inline void formLinearTrans( RC_coll<double> x ) override  {
       
        formLinearTrans_impl(x);

      };

      // Interface to ResponseRBase dcomplex exposure
      inline void formLinearTrans( RC_coll<dcomplex> x ) override {
       
        formLinearTrans_impl(x);

      };


      void printResMO(std::ostream &out) override {

        size_t N      = this->nSingleDim_;
        size_t nRoots = this->resSettings.nRoots;

        auto *VR = this->resResults.VR; 

        printResMO(out,this->resSettings.nRoots,this->resResults.W,
            {{"f     = ",this->resObs.oscStrength}},VR,VR); 

      }

      template <typename U>
      void printResMO_impl(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, U* VL, U* VR);

      virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, double* VL, 
        double* VR) override {
        printResMO_impl(out,nRoots,W,data,VL,VR); 
      };
      virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, dcomplex* VL, 
        dcomplex* VR) override {
        printResMO_impl(out,nRoots,W,data,VL,VR); 
      };


      template <typename U>
      std::vector< std::pair< std::pair<int,int>, U > >
        getMOContributions(U *V, double tol);


      template <typename U>
      void preConditioner(size_t nVec, U shift, SolverVectors<U> &V, SolverVectors<U> &AV);






      // Helper functions
        
      // Transform pp-transition vector MO -> AO
      template <typename U, typename... Args>
      std::vector<TwoBodyContraction<U>> 
        ppTransitionVecMO2AO(MPI_Comm c, bool scatter, size_t nVec, size_t N, 
          Args... Vs);


      // Transform pp-transition vector AO -> MO
      template <typename U, typename... Args>
      void ppTransitionVecAO2MO( size_t nVec, size_t N, 
        std::vector<TwoBodyContraction<U>> &cList, Args... HVs); 

      // Scale pp-transition vector by diagonals
      template <typename U>
      void ppEpsilonScale(bool doInc, bool doInv, bool doOcc, size_t nVec, 
        size_t N, U* V, U* HV);






  }; 


  template <typename MatsT, typename IntsT>
  class ParticleParticlePropagator< HartreeFock<MatsT,IntsT> > : 
    public ParticleParticlePropagator< SingleSlater<MatsT,IntsT> > {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;


    public:


      ParticleParticlePropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<HartreeFock<MatsT,IntsT>> ref, MatsT* fullMatrix = nullptr ) : 
        ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>(c,job, 
            std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ref), fullMatrix) { }

      ParticleParticlePropagator( const ParticleParticlePropagator &other ) : 
        ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>( 
            dynamic_cast<const ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>&>(
              other) 
            ){ }


      // Inherit the ResponseTBase exposures from 
      // ParticleParticlePropagator<SingleSlater>
      using ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::formLinearTrans;

      template <typename U>
      void formLinearTrans_direct_impl(MPI_Comm, RC_coll<U> x);



      // Interface to ParticleParticlePropagator<SingleSlater> double exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<double> x); 

      // Interface to ParticleParticlePropagator<SingleSlater> dcomplex 
      // exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<dcomplex> x){ 
      
        formLinearTrans_direct_impl(c,x);

      };


  }; 



  template <typename MatsT, typename IntsT>
  class ParticleParticlePropagator< KohnSham<MatsT,IntsT> > : 
    public ParticleParticlePropagator< SingleSlater<MatsT,IntsT> > {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;

    public:


      ParticleParticlePropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<KohnSham<MatsT,IntsT>> ref, MatsT* fullMatrix = nullptr ) : 
        ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>(c,job,
            std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ref), fullMatrix) { }

      ParticleParticlePropagator( const ParticleParticlePropagator &other ) : 
        ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>( 
            dynamic_cast<const ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>&>(
              other) 
            ){ }


      // Inherit the ResponseTBase exposures from 
      // ParticleParticlePropagator<SingleSlater>
      using ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::formLinearTrans;

      template <typename U>
      void formLinearTrans_direct_impl(MPI_Comm, RC_coll<U> x);



      // Interface to ParticleParticlePropagator<SingleSlater> double exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<double> x); 

      // Interface to ParticleParticlePropagator<SingleSlater> dcomplex 
      // exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<dcomplex> x){ 
      
        formLinearTrans_direct_impl(c,x);

      };


  }; 



};



