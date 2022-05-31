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
#include <singleslater/neoss.hpp>

#include <response/tbase.hpp>

namespace ChronusQ {

  // Forward decl of Reference templated PolarizationPropagator class
  template <class Reference> class PolarizationPropagator;


  struct SingleSlaterPolarBase { 

      bool doReduced = false;
      bool doAPB_AMB = false;
      bool incMet    = true;
      bool doStab    = false;
      bool doNR      = false;

  };

  struct OrbitalHessianSettings {

    bool useEigenEnergies = false;

  };

  template <typename MatsT, typename IntsT>
  class PolarizationPropagator< SingleSlater<MatsT, IntsT> > : 
    public ResponseRefBase< SingleSlater<MatsT, IntsT> >, 
    public SingleSlaterPolarBase {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;

    protected:

    public:

      // Orbital Hessian Controls
      OrbitalHessianSettings orbHessCtl;





      template <typename U>
      inline void blockTransform(size_t N, size_t nVec, double fact, 
        U* V1, size_t LDV1, U* V2, size_t LDV2) {


        for( auto j = 0ul; j < nVec; j++ )
        for( auto i = 0ul; i < N ; i++ ) {

          U tmp = V1[i + j*LDV1];

          V1[i + j*LDV1] = fact * ( tmp + V2[i + j*LDV2] );
          V2[i + j*LDV2] = fact * ( tmp - V2[i + j*LDV2] );

        }


      }
      inline virtual size_t getNSingleDim(const bool doTDA = false) {

        size_t N = 0;

        auto &ss = *this->ref_;
        size_t nOAVA = ss.nOA * ss.nVA;
        size_t nOBVB = ss.nOB * ss.nVB;
        size_t nOV   = ss.nO  * ss.nV;

        N = (ss.nC == 1) ? nOAVA + nOBVB : nOV;

        if( not doTDA and not doReduced ) N *= 2;

        assert(N != 0);

        return N;

      }



      PolarizationPropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<SingleSlater<MatsT,IntsT>> ref, MatsT* fullMatrix = nullptr ) : 
        ResponseRefBase<SingleSlater<MatsT,IntsT>>(c,job,ref,fullMatrix),
        SingleSlaterPolarBase() { 

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

      PolarizationPropagator( const PolarizationPropagator &other ) : 
        ResponseRefBase<SingleSlater<MatsT, IntsT>>(
          dynamic_cast<const ResponseRefBase<SingleSlater<MatsT, IntsT>>&>(other)
        ), 
        SingleSlaterPolarBase(
          dynamic_cast<const SingleSlaterPolarBase&>(other)
        ) { }

      MatsT*                   formFullFromMemory(); 
      MatsT*                   formFullMatrix(); 
      void                 formRHS();
      std::pair<size_t,MatsT*> formPropGrad(ResponseOperator);
      void                 configOptions();
      void                 eigVecNorm();
      void                 constructShifts();
      void                 postLinearSolve();
      virtual void         resGuess(size_t, SolverVectors<MatsT> &, size_t);

      MatsT getGDiag(size_t, size_t, bool, SingleSlater<MatsT,IntsT>&, MatsT*,
        MatsT*);


      // Internal implementations for direct linear transformation
      virtual void formLinearTrans_direct(MPI_Comm, RC_coll<double>,
          SINGLESLATER_POLAR_COPT, bool noTrans = false) = 0;
      virtual void formLinearTrans_direct(MPI_Comm, RC_coll<dcomplex>,
          SINGLESLATER_POLAR_COPT, bool noTrans = false) = 0;

      // Incore linear transformation
      template <typename U>
      void formLinearTrans_incore_impl(RC_coll<U>, 
          SINGLESLATER_POLAR_COPT);

      void formLinearTrans_incore(RC_coll<double> x, 
          SINGLESLATER_POLAR_COPT op);

      void formLinearTrans_incore(RC_coll<dcomplex> x, 
          SINGLESLATER_POLAR_COPT op) {

        formLinearTrans_incore_impl(x,op);

      }






      // Toggle between direct / incore
      template <typename U>
      void formLinearTrans(RC_coll<U> x, SINGLESLATER_POLAR_COPT op) {

        ProgramTimer::timeOp("Linear Trans", [&]() {
          if( this->genSettings.formFullMat ) formLinearTrans_incore(x,op);
          else  formLinearTrans_direct(this->comm_,x,op);
        });

      };

      // Interface to ResponseRBase double exposure
      inline void formLinearTrans( RC_coll<double> x ) {
       
        formLinearTrans(x,FULL);

      };

      // Interface to ResponseRBase dcomplex exposure
      inline void formLinearTrans( RC_coll<dcomplex> x ) {
       
        formLinearTrans(x,FULL);

      };


      void printResMO(std::ostream &out) {

        size_t N      = this->nSingleDim_;
        size_t nRoots = this->resSettings.nRoots;

        auto *VR = this->resResults.VR; 
        auto *VL = this->resResults.VL; 
        if( not doReduced ) VL = VR + N/2;

        // Transform X+Y/X-Y to X/Y for print
        if( doAPB_AMB ) blockTransform(nRoots,nRoots,std::sqrt(0.5),VR,N,VL,N);

        printResMO(out,this->resSettings.nRoots,this->resResults.W,
            {{"f     = ",this->resObs.oscStrength}},VL,VR); 

        // Transform back to X+Y/X-Y from X/Y 
        if( doAPB_AMB ) blockTransform(nRoots,nRoots,std::sqrt(0.5),VR,N,VL,N);

      }

      template <typename U>
      void printResMO_impl(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, U* VL, U* VR);

      virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, double* VL, 
        double* VR) {
        printResMO_impl(out,nRoots,W,data,VL,VR); 
      };
      virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
        std::vector<std::pair<std::string,double *>> data, dcomplex* VL, 
        dcomplex* VR) {
        printResMO_impl(out,nRoots,W,data,VL,VR); 
      };


      template <typename U>
      std::vector< std::pair< std::pair<int,int>, U > >
        getMOContributions(U *V, double tol);


      template <typename U>
      void preConditioner(size_t nVec, U shift, SolverVectors<U> &V, SolverVectors<U> &AV);







      // Helper functions

      // Transform ph-transition vector MO -> AO
      template <typename U, typename... Args>
      std::vector<TwoBodyContraction<U>> phTransitionVecMO2AO(
        MPI_Comm c,bool scatter, size_t nVec, size_t N,
        SingleSlater<MatsT,IntsT>& ss1, SingleSlater<MatsT,IntsT>& ss2,
        bool doExchange, Args... Vs);

      // Transform ph-transition vector AO -> MO
      template <typename U, typename... Args>
      void phTransitionVecAO2MO(size_t nVec, size_t N, 
        std::vector<TwoBodyContraction<U>> &cList,SingleSlater<MatsT, IntsT>& ss,bool doExchange, Args... HVs); 

      // Scale transition vector by diagonals of orbital Hessian
      template <typename U>
      void phEpsilonScale(bool doInc, bool doInv, size_t nVec, size_t N,SingleSlater<MatsT, IntsT>& ss, 
       U* V, U* HV);
  }; 


  template <typename MatsT, typename IntsT>
  class PolarizationPropagator< HartreeFock<MatsT, IntsT> > : 
    public PolarizationPropagator< SingleSlater<MatsT, IntsT> > {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;

    public:


      PolarizationPropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<HartreeFock<MatsT, IntsT>> ref, MatsT* fullMatrix = nullptr ) : 
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(c,job,
          std::dynamic_pointer_cast<SingleSlater<MatsT, IntsT>>(ref),fullMatrix) { }

      PolarizationPropagator( const PolarizationPropagator &other ) : 
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(
          dynamic_cast<const PolarizationPropagator<SingleSlater<MatsT, IntsT>>&>(other)
        ){ }


      // Inherit the ResponseTBase exposures from 
      // PolarizationPropagator<SingleSlater>
      using PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formLinearTrans;

      template <typename U>
      void formLinearTrans_direct_impl(MPI_Comm, RC_coll<U> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans);



      // Interface to PolarizationPropagator<SingleSlater> double exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<double> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false); 

      // Interface to PolarizationPropagator<SingleSlater> dcomplex exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<dcomplex> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false){ 
      
        formLinearTrans_direct_impl(c,x,op,noTrans);

      };


  }; 



  template <typename MatsT, typename IntsT>
  class PolarizationPropagator< KohnSham<MatsT, IntsT> > : 
    public PolarizationPropagator< SingleSlater<MatsT, IntsT> > {

      template <typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;


    public:


      PolarizationPropagator( MPI_Comm c, ResponseType job, 
        std::shared_ptr<KohnSham<MatsT, IntsT>> ref, MatsT* fullMatrix = nullptr ) : 
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(c,job,
          std::dynamic_pointer_cast<SingleSlater<MatsT, IntsT>>(ref),fullMatrix) { }

      PolarizationPropagator( const PolarizationPropagator &other ) : 
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(
          dynamic_cast<const PolarizationPropagator<SingleSlater<MatsT, IntsT>>&>(other)
        ){ }
      // Inherit the ResponseTBase exposures from 
      // PolarizationPropagator<SingleSlater>
      using PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formLinearTrans;

      template <typename U>
      void formLinearTrans_direct_impl(MPI_Comm, RC_coll<U> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans);



      // Interface to PolarizationPropagator<SingleSlater> double exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<double> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false); 

      // Interface to PolarizationPropagator<SingleSlater> dcomplex exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<dcomplex> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false){ 
      
        formLinearTrans_direct_impl(c,x,op,noTrans);

      };

  };


 
  //Class specialization to deal with NEOSS Objects
  template <typename MatsT, typename IntsT>
  class PolarizationPropagator< NEOSS<MatsT, IntsT> > :
    public PolarizationPropagator< SingleSlater<MatsT, IntsT> >{

      template<typename U>
      using RC_coll = std::vector<RESPONSE_CONTRACTION<U>>;

      std::vector<MatsT> diagonals;

    public:

      PolarizationPropagator( MPI_Comm c, ResponseType job,
        std::shared_ptr<NEOSS<MatsT, IntsT>> ref, MatsT* fullMatrix = nullptr ) :
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(c,job,
          std::dynamic_pointer_cast<SingleSlater<MatsT, IntsT>>(ref),fullMatrix) {
        this->PC_ = [&](size_t nVec, MatsT shift, SolverVectors<MatsT> &V, SolverVectors<MatsT> &AV) {
          neoPreConditioner(nVec,shift,V,AV);
        };
        this->nSPC_ = [&](size_t nVec, SolverVectors<MatsT> &V, SolverVectors<MatsT> &AV) {
          neoPreConditioner(nVec,MatsT(0.),V,AV);
        };
        this->cmplxPC_ = [&](size_t nVec, dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV) {
          neoPreConditioner(nVec,shift,V,AV);
        };
      }

      PolarizationPropagator( const PolarizationPropagator &other ) :
        PolarizationPropagator<SingleSlater<MatsT, IntsT>>(
          dynamic_cast<const PolarizationPropagator<SingleSlater<MatsT, IntsT>>&>(other)
        ){ }

      // Inherit the ResponseTBase exposures from 
      // PolarizationPropagator<SingleSlater>
      using PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formLinearTrans;

      template <typename U>
      void formLinearTrans_direct_impl(MPI_Comm, RC_coll<U> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans);

      virtual size_t getNSingleDim(const bool);

      size_t getNSingleSSDim(SingleSlater<MatsT,IntsT>& , const bool);

      std::pair<size_t,MatsT*> formPropGrad( ResponseOperator );

      template <typename U>
      void neoPreConditioner(size_t nVec, U shift, SolverVectors<U> &V, SolverVectors<U> &AV);

      template <typename U>
      std::vector< std::pair< std::pair<int,int>, U > >
        getMOContributions(U *V, double tol);

      template <typename U>
      void printResMO_impl( std::ostream &out, size_t nRoots, double *W_print,
      std::vector<std::pair<std::string,double *>> data,U* VL, U* VR);

      void resGuess(size_t, MatsT*, size_t);

      // Interface to PolarizationPropagator<SingleSlater> double exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<double> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false); 
      // Interface to PolarizationPropagator<SingleSlater> dcomplex exposure
      void formLinearTrans_direct(MPI_Comm c, RC_coll<dcomplex> x,
          SINGLESLATER_POLAR_COPT op, bool noTrans = false){ 
      
        formLinearTrans_direct_impl(c,x,op,noTrans);

      };
        
      virtual void printResMO(std::ostream &out, size_t nRoots, double *W, 
        std::vector<std::pair<std::string,double *>> data, double* VL, double* VR) {
          printResMO_impl(out,nRoots,W,data,VL,VR);
      };
      virtual void printResMO(std::ostream &out, size_t nRoots, double *W, 
        std::vector<std::pair<std::string,double *>> data, dcomplex* VL, dcomplex* VR) {
          printResMO_impl(out,nRoots,W,data,VL,VR);
      };
 
  };



};

