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

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/factorization.hpp>

#include <util/threads.hpp>

namespace ChronusQ {


  template <typename MatsT, typename IntsT>
  void ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::configOptions() {

    bool isRoot = MPIRank(this->comm_) == 0;


    std::cout << "  * PROPAGATOR = ParticleParticlePropagator" <<
      std::endl << std::endl;

    if( this->genSettings.printLevel > 1  and isRoot )
    std::cout << 
      "  * PERFORMING RESPONSE OPTION CONFIGURATION FOR SINGLE SLATER "
              << "WAVE FUNCTION\n";


    const size_t NO_t  = getNO1(false); 
    const size_t NO2_t = getNO2(false); 

    starOrb.first  = spinSepProp == PP_BB ? -NO_t  : NO_t  ;
    starOrb.second = spinSepProp == PP_AA ?  NO2_t : -NO2_t;


    // Determine if actually TDA
    auto& ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);

    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 


    bool shouldBeATDA = NO < 2;
    bool shouldBeCTDA = NV < 2;

    if( shouldBeATDA and shouldBeCTDA )
      CErr("System Not Suitable for ParticleParticlePropagator: NO < 2 and NV < 2");

    if( shouldBeATDA ) {

      std::cout << "  * Switching to pp-TDA: NO < 2" << std::endl;
      this->genSettings.doTDA = true;
      tdaOp = PP_A;

    }


    if( shouldBeCTDA ) {

      std::cout << "  * Switching to hh-TDA: NV < 2" << std::endl;
      this->genSettings.doTDA = true;
      tdaOp = PP_C;

    }




    // Set the internal NSingleDim and aux dimension
    this->nSingleDim_ = getNSingleDim(this->genSettings.doTDA);
    aDim_       = getADim();
    cDim_       = getCDim();

    doVir_ = tdaOp == PP_A or not this->genSettings.doTDA;
    doOcc_ = tdaOp == PP_C or not this->genSettings.doTDA;




    // Turn off property evaluation (FIXME: for now)
    this->genSettings.evalProp = false;

    // No such thing as FDR for PPRPA/TDA FIXME: techinically yes, but NYI
    this->genSettings.jobType = RESIDUE;

    // TDA -> Hermetian
    this->genSettings.matIsHer = this->genSettings.doTDA;



    // Direct logic
    if( not this->genSettings.formFullMat ) {

      // Turn off distribution
      this->genSettings.formMatDist     = false;
      this->genSettings.distMatFromRoot = false;

    }


    // Logic for full RESIDUE problem
    if( this->genSettings.doFull and this->genSettings.jobType == RESIDUE ) {

      if( this->genSettings.printLevel > 0 and isRoot)
        std::cout << "  * SETTING NROOTS = NSINGLEDIM\n";

      // Nroots = NSingleDim
      this->resSettings.nRoots = this->nSingleDim_;

    }

    // Compute MU
    if( ss.nC == 1 and not ss.iCS ) {
  

      double HOMO = -std::numeric_limits<double>::infinity();
      if( ss.nOA > 0 ) HOMO = ss.eps1[ss.nOA - 1];
      if( ss.nOB > 0 ) HOMO = std::max(HOMO,ss.eps2[ss.nOB - 1]);

      double LUMO = std::min(ss.eps1[ss.nOA],ss.eps2[ss.nOB]);

      mu = 0.5 * ( HOMO + LUMO );

    } else 
      mu = 0.5 * (ss.eps1[NO] + ss.eps1[NO-1]);

    std::cout << "  * ParticleParticlePropagator has chosen MU = " 
              << mu << std::endl;

    // Output matrix beging diagonalized
    if( ss.nC == 1 ) {

      std::cout << "  * PPSPINMAT = ";
      if( spinSepProp == PP_AA )       std::cout << "AA";
      else if( spinSepProp == PP_AB )  std::cout << "AB";
      else if( spinSepProp == PP_BB )  std::cout << "BB";
      std::cout << std::endl;

    }

    // Ouput TDA matrix
    if( this->genSettings.doTDA ) {


      std::cout << "  * DOING TDA!" << std::endl;

      std::cout << "  * PPTDAMAT = ";
      if( doVir_ )       std::cout << "A";
      else if( doOcc_ )  std::cout << "C";
      std::cout << std::endl;

    }

    std::cout << std::endl << std::endl;


    if( this->genSettings.printLevel > 0 and isRoot ) {
      if( not this->genSettings.formFullMat )
        std::cout << "  * TENSOR CONTRACTIONS WILL BE FORMED DIRECTLY" 
          << std::endl;
    }

  };


  template <typename MatsT, typename IntsT>
  void ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::constructShifts() {

    CErr("FDR Not Implemented for ParticleParticlePropagator: no shifts!");

  };

  template <typename MatsT, typename IntsT>
  void ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::postLinearSolve() {

    CErr("FDR Not Implemented for ParticleParticlePropagator: no postLinearSolve!");

  };

  template <typename MatsT, typename IntsT>
  template <typename U>
  std::vector< std::pair< std::pair<int,int>, U > >
    ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::getMOContributions(U *V, 
        double tol) {

    std::vector< std::pair< std::pair<int,int>, U > > moCont;

    auto& ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);


    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 




    bool doLT = (ss.nC == 2) or (spinSepProp != PP_AB);

    size_t hhSt = this->genSettings.doTDA ? 0 : aDim_;

    auto bmax = [&](size_t a){ return doLT ? a : NV2; };
    auto jmax = [&](size_t i){ return doLT ? i : NO2; };

    if( doVir_ )
    for(auto a = 0ul, ab = 0ul; a < NV;      a++      )
    for(auto b = 0ul;           b < bmax(a); b++, ab++) {

      if( std::abs(V[ab]) > tol ) 
        moCont.push_back( { {a+NO, b+NO2}, V[ab] } );

    }


    if( doOcc_ ) 
    for(auto i = 0ul, ij = hhSt; i < NO;      i++      )
    for(auto j = 0ul;            j < jmax(i); j++, ij++) {

      if( std::abs(V[ij]) > tol ) 
        moCont.push_back( { {-i, -j}, V[ij] } );

    }


    return moCont;

  };


  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::printResMO_impl(
    std::ostream &out, size_t nRoots, double *W_print,
    std::vector<std::pair<std::string,double *>> data, U* VL, U* VR) {

    this->nSingleDim_ = getNSingleDim(this->genSettings.doTDA);

    out << "\n\n\n* RESIDUE EIGENMODES\n\n\n";

    for(auto iRt = 0; iRt < nRoots; iRt++) {

      out << "  Root " << std::setw(7) << std::right << iRt+1 << ":";

      // Energy eigenvalues in various unit systems
      out << std::setw(15) << std::right << "W(Eh) = " 
          << std::setprecision(8) << std::fixed << W_print[iRt];

      out << std::setw(15) << std::right << "W(eV) = " 
          << std::setprecision(8) << std::fixed 
          << W_print[iRt]*EVPerHartree;

      out << "\n";


      if( this->genSettings.evalProp ) {
        for(auto &d : data) {
        out << "       " << std::setw(7) << " " << " ";
        out << std::setw(15) << std::right << d.first 
            << std::setprecision(8) << std::fixed 
            << d.second[iRt];

        out << "\n";
        }
      }

      auto vCont = getMOContributions(VR+iRt*this->nSingleDim_,1e-1);

      // MO contributions
      out << "    MO Contributions:\n";
      for(auto &c : vCont) {

        std::string exLabel = (c.first.first > 0) ? "(++)" : "(--)";

        char spinLabel1 = (this->ref_->nC == 2)  ? ' ' :
                          (spinSepProp == PP_BB) ? 'B' : 'A';
        char spinLabel2 = (this->ref_->nC == 2)  ? ' ' :
                          (spinSepProp == PP_AA) ? 'A' : 'B';

      
        out << "      " << exLabel;
        out << std::setw(4) << std::right 
            << std::abs(c.first.first) + 1 << spinLabel1 << " ; ";
        out << std::setw(4) << std::right 
            << std::abs(c.first.second) + 1  << spinLabel2;

        if(std::is_same<U,double>::value)
          out << "  " << std::fixed << std::setprecision(5) 
                      << std::setw(10) << std::right << c.second << "\n";
        else {
          out << "  " << std::fixed << std::setprecision(5) 
                      << std::setw(10) << std::right << std::abs(c.second);
          out << "  " << std::fixed << std::setprecision(5) 
                      << std::setw(10) << std::right << std::arg(c.second) 
                      << "\n";
        }



      }


      out << "\n" << std::endl;

    }

  }
  

  template <typename MatsT, typename IntsT>
  void ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::eigVecNorm() {

    if(this->genSettings.doTDA) return;

    auto N  = this->nSingleDim_;
    size_t nRoots = this->resSettings.nRoots;

    MatsT* V = this->resResults.VR;

    // Fixme does not propery orthonormalize Y vectors
    auto tdInner = [&](MatsT* Vc) { 
        MatsT inX = blas::dot(aDim_,Vc,1,Vc,1);
        MatsT inY = blas::dot(cDim_,Vc+aDim_,1,Vc+aDim_,1) ;

        bool sign = std::signbit(std::abs(inX) - std::abs(inY));

        int fact = sign ? 1 : -1;

        return fact * std::sqrt(std::abs(inX - inY));
      };




    auto tdMatInner = [&](size_t i, size_t j, MatsT* Vi, size_t LDVi,
          MatsT* Vj, size_t LDVj, MatsT* inner){ 

          MatsT    *Xi = Vi,            *Xj = Vj;
          MatsT    *Yi = Vi + aDim_,    *Yj = Vj + aDim_;

          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,aDim_,MatsT(1.) ,Xi,LDVi,Xj,LDVj,MatsT(0.),inner,i);
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,cDim_,MatsT(-1.),Yi,LDVi,Yj,LDVj,MatsT(1.),inner,i);

        };


    /*
    MatsT* inner = CQMemManager::get().malloc<MatsT>(nRoots*nRoots);

    tdMatInner(N,N,V,N,V,N,inner);
    prettyPrintSmart(std::cout,"BEFORE",inner,N,N,N);
    */

    GramSchmidt(N,0,nRoots,V,N,tdInner,tdMatInner,1);

    /*
    tdMatInner(N,N,V,N,V,N,inner);
    prettyPrintSmart(std::cout,"AFTER",inner,N,N,N);
    */


  };



  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator< SingleSlater<MatsT, IntsT> >::formLinearTrans_incore_impl( 
    std::vector<RESPONSE_CONTRACTION<U>> cList){ 
  
    bool isDist = this->genSettings.isDist(); 

    /*
    // FIXME MPI : This is ad hoc logic, assumes that if the cases when
    // matrix is not distributed that only the root process will be
    // calling this function
    bool needFullMat = not bool(this->fullMatrix_);
#ifdef CQ_ENABLE_MPI
    if( this->genSettings.isDist() ) MPIBCast(needFullMat,0,this->comm_);
#endif
    if(needFullMat) formFullMatrix(); 
    */

    size_t N = this->nSingleDim_;


    size_t maxNVec = 
      std::max_element(cList.begin(),cList.end(),
        [](RESPONSE_CONTRACTION<U> l, RESPONSE_CONTRACTION<U> r) {
          return l.nVec < r.nVec;
        }
      )->nVec;


    // (Local) Dims
    int64_t MLoc(N) , NLoc(N);

#ifdef CQ_ENABLE_MPI
    if( isDist )
      std::tie(MLoc,NLoc) = this->fullMatGrid_->get_local_dims(MLoc,NLoc);


    // ScaLAPACK DESC
    scalapackpp::scalapack_desc descMem;
    if( isDist ) descMem = this->fullMatGrid_->descinit_noerror(N,N,MLoc);
#endif

    



    // FIXME MPI: Because there doesn't exist a PDZGEMM, make a copy of the
    // matrix if the contractions are complex... This is horrifying!
    U* FM = nullptr;
#ifdef CQ_ENABLE_MPI
    if( isDist ) {

      if( std::is_same<U,MatsT>::value ) 
        FM = reinterpret_cast<U*>(this->fullMatrix_);
      else {
        FM = CQMemManager::get().malloc<U>(MLoc*NLoc);
        std::copy_n(reinterpret_cast<double*>(this->fullMatrix_),MLoc*NLoc,FM);
      }

    }
#endif

    for(auto &mat : cList) {

      size_t nVec = mat.nVec;
      auto * X    = mat.X;
      auto *AX    = mat.AX;

#ifdef CQ_ENABLE_MPI
      scalapackpp::scalapack_desc DescAX = mat.DescAX;
#endif

#ifdef CQ_ENABLE_MPI
      if( isDist )
        Gemm_MPI('N','N',N,nVec,N,U(1.),FM,1,1,descMem,X,1,1,mat.DescX,
          U(0.),AX,1,1,DescAX);
      else
#endif
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_,N,X,N,U(0.),AX,N);

      //if( not this->genSettings.matIsHer ) {
      //  if( isDist )
      //    CXXBLACS::PLASCL('F',-1.,1.,CDIM,nVec,AX,ADIM+1,1,DescAX);
      //  else
      //    SetMat('N',CDIM,nVec,U(-1.),AX + ADIM,N,AX + ADIM,N);
      //}



    }
  
    if( FM and (reinterpret_cast<MatsT*>(FM) != this->fullMatrix_) ) 
      CQMemManager::get().free(FM);
  
  };

  template <>
  void ParticleParticlePropagator< SingleSlater<double,double> >::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList){ 

    formLinearTrans_incore_impl(cList);

  }

  template <>
  void ParticleParticlePropagator<SingleSlater<dcomplex,double>>::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList){ 

    CErr("SOMETHING HAS GONE HORRIBLY WRONG formLinearTrans_incore");

  }

  template <>
  void ParticleParticlePropagator<SingleSlater<dcomplex,dcomplex>>::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList){ 

    CErr("SOMETHING HAS GONE HORRIBLY WRONG formLinearTrans_incore");

  }




  template <typename MatsT, typename IntsT>
  MatsT * ParticleParticlePropagator<SingleSlater<MatsT, IntsT>>::formFullMatrix() {
  
    bool isRoot = MPIRank(this->comm_) == 0;
    
    if( isRoot ) std::cout << "\n  * FORMING FULL MATRIX\n";
    
    auto topFull = tick();

    this->nSingleDim_ = this->getNSingleDim(this->genSettings.doTDA);
    size_t N = this->nSingleDim_;

    size_t nForm  = N;
    size_t nStore = nForm;
    size_t nFormP = nForm;   // Persistant for matrix distribution
    size_t nStoreP = nStore; // Persistant for matrix distribution

    // Determine if we want to distribute the matrix at all
    bool isDist = this->genSettings.isDist(); 

    // MPI Communicator for matrix formation
    MPI_Comm matComm = this->comm_;


#ifdef CQ_ENABLE_MPI
    std::shared_ptr<scalapackpp::BlockCyclicDist2D> formGrid;
    // Divide up the work if forming the matrix distributed
    if( this->genSettings.formMatDist ) { 

      if( isRoot )
      std::cout << "    * FORMING THE MATRIX DISTRIBUTED\n";

      // Split the communicator
      int color = MPIRank(this->comm_);
      matComm = MPICommSplit(this->comm_,color,0);
      

      // Create a temorary grid to form the matrix in
      // columns
      formGrid = std::make_shared<scalapackpp::BlockCyclicDist2D>(
          blacspp::Grid(this->comm_, 1, MPISize(this->comm_)),
          N, nForm / MPISize(this->comm_));
      std::tie(N,nForm) = formGrid->get_local_dims(N,nForm);
      nStore = nForm;

    }
#endif

    // Print memory requirements (ish)
    if( isRoot ) {
      size_t NB = this->ref_->nAlphaOrbital();
      size_t vcSize = N*(nForm + nStore)*sizeof(double);
      size_t aoSize = NB*NB*(nForm + nStore)*sizeof(double);
      size_t parSize = nForm*NB*NB*sizeof(double)*(GetNumThreads()-1);
      std::cout << "    * TOTAL VECTOR STORAGE REQUIREMENT             "
        << double(vcSize) / 1e9 << " GB\n";
      std::cout << "    * TOTAL AO SCRATCH STORAGE REQUIREMENT         "
        << double(aoSize) / 1e9 << " GB\n";
      if( GetNumThreads() - 1 )
        std::cout << "    * TOTAL PARALLEL SCRATCH STORAGE REQUIREMENT "
          << double(parSize) / 1e9 << " GB\n";
      std::cout << "\n" << std::flush;

    }

    bool isRootMatComm = MPIRank(matComm) == 0;


    // Allocate space for identity for contraction
    MatsT* V  = CQMemManager::get().malloc<MatsT>(N*nForm);
    std::fill_n(V ,N*nForm ,0.);
    
    // Form the identity in V
#ifdef CQ_ENABLE_MPI
    if( this->genSettings.formMatDist ){ 
      size_t j = 0;
      for(auto i = 0ul; i < nFormP; i++){ 

        if( formGrid->i_own(0,i) ) {
          V[i + j*N] = 1.0; 
          j++;
        }
      }
    } else
#endif
      for(auto i = 0ul; i < nForm; i++){ V[i * (N+1)] = 1.0; }

    // Only allocate HV on root process as it won't be touched by
    // other processes
    MatsT* HV = nullptr;
    if( isRootMatComm ) {
      HV = CQMemManager::get().malloc<MatsT>(N*nStore);
      //std::fill_n(HV,N*nStore,0.);
    }
  


    // Perform contraction directly
    auto topContract = tick();
    RC_coll<MatsT> cList(1);
    cList.back().nVec = nForm;
    cList.back().N    = N;
    cList.back().X    = V;
    cList.back().AX   = HV;

    formLinearTrans_direct(matComm,cList);

    auto durContract = tock(topContract);


    CQMemManager::get().free(V); // Free up some memory

    auto topTrans = tick();
    if( not this->genSettings.formMatDist  and isRootMatComm ){ 

      // Negate the bottom half
      if( not this->genSettings.doTDA ) 
        SetMat('N',cDim_,N,MatsT(-1.),HV + aDim_,N,HV + aDim_,N);
      else if( doOcc_ )
        SetMat('N',cDim_,N,MatsT(-1.),HV,N,HV,N);

    }

    double durTrans = tock(topTrans);
  


    auto topDist = tick();
#ifdef CQ_ENABLE_MPI
    if( isDist ) {

      // Create the grid
      int64_t MB = this->genSettings.MB;
      this->fullMatGrid_ = 
        std::make_shared<scalapackpp::BlockCyclicDist2D>(
          blacspp::Grid::square_grid(this->comm_),MB,MB);

      // Get local dims
      int64_t MLoc, NLoc;
      std::tie(MLoc,NLoc) = this->fullMatGrid_->get_local_dims(N,nStoreP);

      // Set to ScaLAPACK descriptors
      this->descFullMat_ = 
        this->fullMatGrid_->descinit_noerror(N,nStoreP,MLoc);


      if( isRoot ) {

        int64_t NPROCROW = this->fullMatGrid_->grid().npr();
        int64_t NPROCCOL = this->fullMatGrid_->grid().npc();
        size_t NLOCAL = MLoc*NLoc*sizeof(MatsT);

        std::cout << "  * DISTRIBUTING FULL MATRIX TO BLACS GRID\n";
        std::cout << "    * MB       = " << MB << "\n";
        std::cout << "    * NPROCROW = " <<  NPROCROW << "\n";
        std::cout << "    * NPROCCOL = " <<  NPROCCOL << "\n";
        std::cout << "    * LOCALBUF = " << double(NLOCAL) / 1e6 << " MB\n";
        
        std::cout << "\n";

      }

      // Allocate local buffers
      this->fullMatrix_ = nullptr;
      if( MLoc and NLoc )
        this->fullMatrix_ = CQMemManager::get().malloc<MatsT>(MLoc*NLoc);

      if( this->genSettings.distMatFromRoot )
        this->fullMatGrid_->scatter(N,nStoreP,HV,N,this->fullMatrix_,MLoc,0,0);
      else {

        // Get the array descriptiors on the form Grid
        auto curDesc = formGrid->descinit_noerror(N,nFormP,N);

        // Scale / conjugate if necessary Conj(B) -> -Conj(B); C -> -C
       if( not this->genSettings.doTDA )
         SetMat('N',cDim_,nForm,MatsT(-1.),HV + aDim_,N,HV + aDim_,N);
       else if( doOcc_ )
         SetMat('N',cDim_,nForm,MatsT(-1.),HV,N,HV,N);

        // Redistribute
        scalapackpp::wrappers::pgemr2d(N,nFormP,HV,1,1,curDesc,this->fullMatrix_,1,1,
          this->descFullMat_,this->fullMatGrid_->grid().context());

      }


      if( HV ) CQMemManager::get().free(HV);

    } else 
#endif
      this->fullMatrix_ = HV;

    double durDist = tock(topDist);

    //if( isDist ) CErr();


    if( matComm != this->comm_ ) MPICommFree(matComm);

    if( this->savFile.exists() and isRoot and not isDist ) {

      this->savFile.safeWriteData("/RESP/FULLMATRIX",HV,{nStore,N});

    }

    double durFull =  tock(topFull);
    if( isRoot ) {

      std::cout << "  * TIMINGS\n";
      std::cout << "    * TOTAL     " << durFull << " s\n";
      std::cout << "    * CONTRACT  " << durContract << " s\n";
      if( not this->genSettings.doTDA )
        std::cout << "    * TRANS     " << durTrans << " s\n";
#ifdef CQ_ENABLE_MPI
      if( isDist )
        std::cout << "    * DIST      " << durDist << " s\n";
#endif

      std::cout << "\n";
    }


    MPI_Barrier(this->comm_); // Sync MPI processes

    /*
    if( isRoot ) {

      prettyPrintSmart(std::cout,"FM",this->fullMatrix_,N,nStore,N);

    }
    CErr();
    */
  
    return this->fullMatrix_;

  };












  template <typename MatsT, typename IntsT>
  MatsT * ParticleParticlePropagator<SingleSlater<MatsT,IntsT>>::formFullFromMemory() {

    bool isDist = this->genSettings.isDist(); 

    // FIXME MPI : This is ad hoc logic, assumes that if the cases when
    // matrix is not distributed that only the root process will be
    // calling this function
    bool needFullMat = not bool(this->fullMatrix_);
#ifdef CQ_ENABLE_MPI
    if( this->genSettings.isDist() ) MPIBCast(needFullMat,0,this->comm_);
#endif
    if(needFullMat) this->formFullMatrix(); 


    return this->fullMatrix_;

  };


  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator<SingleSlater<MatsT, IntsT>>::preConditioner(size_t nVec, 
      U shift, SolverVectors<U> &V, SolverVectors<U> &AV) {


    auto& ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);

    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 



    bool doLT = (ss.nC == 2) or (spinSepProp != PP_AB);
    size_t hhSt = this->genSettings.doTDA ? 0 : aDim_;







    auto bmax = [&](size_t a){ return doLT ? a : NV2; };
    auto jmax = [&](size_t i){ return doLT ? i : NO2; };

    double *eps1 = nullptr; 
    double *eps2 = nullptr;

    if( ss.nC == 2 or ss.iCS ) {

      eps1 = ss.eps1;
      eps2 = ss.eps1;

    } else {

      if( spinSepProp == PP_AA ) {

        eps1 = ss.eps1;
        eps2 = ss.eps1;

      } else if( spinSepProp == PP_AB ) {

        eps1 = ss.eps1;
        eps2 = ss.eps2;

      } else {

        eps1 = ss.eps2;
        eps2 = ss.eps2;

      }

    }
                 

    size_t NS = this->nSingleDim_;

    bool isRaw = V.underlyingType() == typeid(RawVectors<U>)
            and AV.underlyingType() == typeid(RawVectors<U>);
    
    if (isRaw and MPIRank(this->comm_) == 0) {
      U * V_ptr = tryGetRawVectorsPointer(V);
      U * AV_ptr = tryGetRawVectorsPointer(AV);
      for(auto iVec = 0ul; iVec < nVec; iVec++) {

        auto * AVk = AV_ptr + iVec * NS;
        auto * Vk  = V_ptr  + iVec * NS;

        for(auto a = 0ul, ab = 0ul; a < NV;      a++      )
        for(auto b = 0ul;           b < bmax(a); b++, ab++)
          AVk[ab] = Vk[ab] / (eps1[a + NO] + eps2[b + NO2] - 2.*mu);

        if( doOcc_ )
        for(auto i = 0ul, ij = hhSt; i < NO;      i++      )
        for(auto j = 0ul;            j < jmax(i); j++, ij++)
          AVk[ij] = - Vk[ij] / (eps1[i] + eps2[j]- 2.*mu);


      } // loop over vectors
    } else if (not isRaw) {
      for(auto iVec = 0ul; iVec < nVec; iVec++) {


        for(auto a = 0ul, ab = 0ul; a < NV;      a++      )
        for(auto b = 0ul;           b < bmax(a); b++, ab++)
          AV.set(ab, iVec, V.get(ab, iVec) / (eps1[a + NO] + eps2[b + NO2] - 2.*mu));

        if( doOcc_ )
        for(auto i = 0ul, ij = hhSt; i < NO;      i++      )
          for(auto j = 0ul;            j < jmax(i); j++, ij++)
            AV.set(ij, iVec, - V.get(ij, iVec) / (eps1[i] + eps2[j]- 2.*mu));


      } // loop over vectors
    }

  };


  template <typename MatsT, typename IntsT>
  std::pair<size_t,MatsT*> 
    ParticleParticlePropagator<SingleSlater<MatsT, IntsT>>::formPropGrad(
      ResponseOperator op) {


    CErr("FDR Not Implemented for ParticleParticlePropagator: no formPropGrad!");

    return std::pair<size_t,MatsT*>(); 
  };

  template <typename MatsT, typename IntsT>
  void ParticleParticlePropagator<SingleSlater<MatsT, IntsT>>::formRHS() {

    CErr("FDR Not Implemented for ParticleParticlePropagator: no formRHS!");

  }
  












  template <typename MatsT, typename IntsT>
  template <typename U>
  void ParticleParticlePropagator< SingleSlater<MatsT, IntsT> >::ppEpsilonScale(
    bool doInc, bool doInv, bool doOcc, size_t nVec, size_t N, U* V, U* HV) {




    SingleSlater<MatsT,IntsT> &ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);




    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 

    const bool   doLT  = (ss.nC == 2) or (spinSepProp != PP_AB);

    const auto bmax = [&](size_t a){ return doLT ? a : NV2; };
    const auto jmax = [&](size_t i){ return doLT ? i : NO2; };

    auto epsilonScale = doInv ? []( double x ) { return 1./x; } :
                                []( double x ) { return x;    };

    auto increment    = doInc ? []( U x, U y ) { return x + y; } :
                                []( U x, U y ) { return y;  };


    double *eps1 = nullptr; 
    double *eps2 = nullptr;

    if( ss.nC == 2 or ss.iCS ) {

      eps1 = ss.eps1;
      eps2 = ss.eps1;

    } else {

      if( spinSepProp == PP_AA ) {

        eps1 = ss.eps1;
        eps2 = ss.eps1;

      } else if( spinSepProp == PP_AB ) {

        eps1 = ss.eps1;
        eps2 = ss.eps2;

      } else {

        eps1 = ss.eps2;
        eps2 = ss.eps2;

      }

    }

    for(size_t iVec = 0; iVec < nVec; iVec++) {


      U* V_c  = V  + iVec * N;
      U* HV_c = HV + iVec * N;

      if( doOcc )
        for(auto i = 0ul, ij = 0ul; i < NO;      i++      )
        for(auto j = 0ul;           j < jmax(i); j++, ij++) 
          HV_c[ij] = increment(HV_c[ij],
                               epsilonScale(-(eps1[i] + eps2[j] - 2.*mu)) * 
                               V_c[ij]);

      else
        for(auto a = 0ul, ab = 0ul; a < NV;      a++      )
        for(auto b = 0ul;           b < bmax(a); b++, ab++) 
          HV_c[ab] = increment(HV_c[ab],
                               epsilonScale(eps1[a + NO] + eps2[b + NO2] - 
                                 2.*mu) * V_c[ab]);

    }


  };






  template <typename MatsT, typename IntsT>
  template <typename U, typename... Args>
  std::vector<TwoBodyContraction<U>> 
    ParticleParticlePropagator< SingleSlater<MatsT, IntsT> >::ppTransitionVecMO2AO(
      MPI_Comm c, bool scatter, size_t nVec, size_t N, Args... Vs){ 
    
    static_assert(sizeof...(Vs) > 0 and sizeof...(Vs) < 3,
      "Vs must consist of 1 or 2 pointers");

    constexpr size_t nVs = sizeof...(Vs);

    std::array<U*,nVs> V_arr = { Vs... };


    bool isRoot = MPIRank(c) == 0;
    bool trans  = isRoot or not scatter;


    SingleSlater<MatsT,IntsT> &ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);
    const size_t NB   = ss.nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = ss.nC * NB;
    const size_t NBC2 = NBC * NBC;



    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 


    const bool   doLT  = (ss.nC == 2) or (spinSepProp != PP_AB);

    const auto bmax = [&](size_t a){ return doLT ? a : NV2; };
    const auto jmax = [&](size_t i){ return doLT ? i : NO2; };
    
    U* MOT  = CQMemManager::get().malloc<U>(NBC2);
    U* SCR  = trans ? CQMemManager::get().malloc<U>(NBC2) : nullptr;


    std::vector<TwoBodyContraction<U>> cList;

    auto MOTRANS = [&]( MatsT* CMO1, MatsT* CMO2, U* X ) {
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBC,NBC,NBC,U(1.0),CMO1,NBC,X  ,NBC,U(0.0),SCR,NBC); 
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NBC,NBC,NBC,U(1.0),CMO2,NBC,SCR,NBC,U(0.0),X  ,NBC);
      IMatCopy('C',NBC,NBC,U(1.),X,NBC,NBC);
    };



    size_t nMatPVec = (ss.nC == 1) ? 2 : 8;
    size_t nAlloc = nVec * nMatPVec * NB2;
    //std::cerr << "MO2AO " << nAlloc*sizeof(MatsT) / 1e9 << std::endl;

    U * first = CQMemManager::get().malloc<U>(nAlloc);

    MatsT *mo1 = nullptr; 
    MatsT *mo2 = nullptr;

    if( ss.nC == 2 or ss.iCS ) {

      mo1 = ss.mo[0].pointer();
      mo2 = ss.mo[0].pointer();

    } else {

      if( spinSepProp == PP_AA ) {

        mo1 = ss.mo[0].pointer();
        mo2 = ss.mo[0].pointer();

      } else if( spinSepProp == PP_AB ) {

        mo1 = ss.mo[0].pointer();
        mo2 = ss.mo[1].pointer();

      } else {

        mo1 = ss.mo[1].pointer();
        mo2 = ss.mo[1].pointer();

      }

    }

    for(size_t iVec = 0; iVec < nVec; iVec++) {

      //U* VX_c = V_arr[0] + iVec * N;
      //U* VY_c = nullptr; 
      //if( nVs > 1 ) VY_c = V_arr[1] + iVec * N;


      U* VX_c = doVir_ ? V_arr[0] : nullptr;
      U* VY_c = (nVs > 1) ? V_arr[1] : doOcc_ ? V_arr[0] : nullptr;

      VX_c += iVec * N;
      VY_c += iVec * N;



      // Perform MO -> AO transformation
      if( trans ) {

        // Zero out MOT
        std::fill_n(MOT,NBC2,0.);

        // Place V -> MOT and transform into AO basis

        // VV block
        if( doVir_ )
        for(auto a = 0ul, ab = 0ul; a < NV;       a++      )
        for(auto b = 0ul          ; b < bmax(a) ; b++, ab++) {

          auto A = a + NO;
          auto B = b + NO2;

          MOT[A + B*NBC] = VX_c[ab];
          if( doLT ) MOT[B + A*NBC] = -VX_c[ab];

        }

        // OO block
        if( doOcc_ )
        for(auto i = 0ul, ij = 0ul; i < NO;       i++      )
        for(auto j = 0ul          ; j < jmax(i) ; j++, ij++) {

          MOT[i + j*NBC] = VY_c[ij];
          if( doLT ) MOT[j + i*NBC] = -VY_c[ij];

        }


        MOTRANS(mo1,mo2,MOT);

      }

      // Scatter results if requested
      if( scatter ) MPIBCast(MOT,NBC*NBC,0,c);

      U* AOaa(nullptr),  *AOab(nullptr),  *AOba(nullptr),  *AObb(nullptr);
      U* KAOaa(nullptr), *KAOab(nullptr), *KAOba(nullptr), *KAObb(nullptr);

      AOaa  = first; 
      KAOaa = first + NB2;

      cList.push_back( {AOaa, KAOaa, false, EXCHANGE } );

      if( ss.nC == 2 ) {

        AOab  = first + 2*NB2; 
        KAOab = first + 3*NB2;
        AOba  = first + 4*NB2; 
        KAOba = first + 5*NB2;
        AObb  = first + 6*NB2; 
        KAObb = first + 7*NB2;

        cList.push_back( {AOab, KAOab, false, EXCHANGE } );
        cList.push_back( {AOba, KAOba, false, EXCHANGE } );
        cList.push_back( {AObb, KAObb, false, EXCHANGE } );

      }
      
      if( ss.nC == 1 ) {

        std::copy_n(MOT,NB2,AOaa);
        std::fill_n(KAOaa,NB2,0.);

      } else {

        std::fill_n(KAOaa,NB2,0.);
        std::fill_n(KAOab,NB2,0.);
        std::fill_n(KAOba,NB2,0.);
        std::fill_n(KAObb,NB2,0.);

        SetMat('N',NB,NB,U(1.),MOT             ,NBC,AOaa,NB);
        SetMat('N',NB,NB,U(1.),MOT + NB*NBC    ,NBC,AOab,NB);
        SetMat('N',NB,NB,U(1.),MOT + NB        ,NBC,AOba,NB);
        SetMat('N',NB,NB,U(1.),MOT + NB*(NBC+1),NBC,AObb,NB);

	}

      first += nMatPVec * NB2;



    }


    CQMemManager::get().free(MOT);
    if( SCR ) CQMemManager::get().free(SCR);


    return cList;



  };






  template <typename MatsT, typename IntsT>
  template <typename U, typename... Args>
  void ParticleParticlePropagator< SingleSlater<MatsT, IntsT> >::ppTransitionVecAO2MO(
    size_t nVec, size_t N, std::vector<TwoBodyContraction<U>> &cList, 
    Args... HVs){ 
  
    constexpr size_t nHVs = sizeof...(HVs);

    std::array<U*,nHVs> HV_arr = { HVs... };


    SingleSlater<MatsT,IntsT> &ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);
    const size_t NB   = ss.nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = ss.nC * NB;
    const size_t NBC2 = NBC * NBC;



    const size_t NV  = getNV1(doStarRef); 
    const size_t NV2 = getNV2(doStarRef);
    const size_t NO  = getNO1(doStarRef); 
    const size_t NO2 = getNO2(doStarRef); 



    const bool   doLT  = (ss.nC == 2) or (spinSepProp != PP_AB);

    const auto bmax = [&](size_t a){ return doLT ? a : NV2; };
    const auto jmax = [&](size_t i){ return doLT ? i : NO2; };
    
    U* MOT  = CQMemManager::get().malloc<U>(NBC2);
    U* SCR  = CQMemManager::get().malloc<U>(NBC2);
  
    auto MOTRANS = [&]( MatsT* CMO1, MatsT* CMO2, U* X ) {
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NBC,NBC,NBC,U(1.0),CMO1,NBC,X  ,NBC,U(0.0),SCR,NBC); 
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::ConjTrans,NBC,NBC,NBC,U(1.0),CMO2,NBC,SCR,NBC,U(0.0),X  ,NBC);
      IMatCopy('C',NBC,NBC,U(1.),X,NBC,NBC);
    };

    MatsT *mo1 = nullptr; 
    MatsT *mo2 = nullptr;

    if( ss.nC == 2 or ss.iCS ) {

      mo1 = ss.mo[0].pointer();
      mo2 = ss.mo[0].pointer();

    } else {

      if( spinSepProp == PP_AA ) {

        mo1 = ss.mo[0].pointer();
        mo2 = ss.mo[0].pointer();

      } else if( spinSepProp == PP_AB ) {

        mo1 = ss.mo[0].pointer();
        mo2 = ss.mo[1].pointer();

      } else {

        mo1 = ss.mo[1].pointer();
        mo2 = ss.mo[1].pointer();

      }

    }

    size_t nMat = (ss.nC == 1) ? 1 : 4;
    for(size_t iVec = 0; iVec < nVec; iVec++) {

      //U* HVX_c = HV_arr[0] + iVec * N;
      //U* HVY_c = nullptr;
      //if( nHVs > 1 ) HVY_c = HV_arr[1] + iVec * N;

      U* HVX_c = doVir_ ? HV_arr[0] : nullptr;
      U* HVY_c = (nHVs > 1) ? HV_arr[1] : doOcc_ ? HV_arr[0] : nullptr;

      HVX_c += iVec * N;
      HVY_c += iVec * N;

      size_t indx = iVec * nMat;

      U* Kaa = cList[indx].AX;
      U* Kab(nullptr), *Kba(nullptr), *Kbb(nullptr);

      if( ss.nC == 2 ) {

        Kab = cList[indx + 1].AX;
        Kba = cList[indx + 2].AX;
        Kbb = cList[indx + 3].AX;

      }


      if( ss.nC == 1 ) std::copy_n(Kaa,NB2,MOT);
      else {

        SetMat('N',NB,NB,U(1.),Kaa,NB,MOT             ,NBC);
        SetMat('N',NB,NB,U(1.),Kab,NB,MOT + NB*NBC    ,NBC);
        SetMat('N',NB,NB,U(1.),Kba,NB,MOT + NB        ,NBC);
        SetMat('N',NB,NB,U(1.),Kbb,NB,MOT + NB*(NBC+1),NBC);

      }



      MOTRANS(mo1,mo2,MOT);


      if( doVir_ )
      for(auto a = 0ul, ab = 0ul; a < NV;       a++      )
      for(auto b = 0ul          ; b < bmax(a) ; b++, ab++){ 

        auto A = a + NO;
        auto B = b + NO2;

        HVX_c[ab] += MOT[A + B*NBC];

      }

      if( doOcc_ )
      for(auto i = 0ul, ij = 0ul; i < NO;       i++      )
      for(auto j = 0ul          ; j < jmax(i) ; j++, ij++) 
        HVY_c[ij] += MOT[i + j*NBC];


    }

    CQMemManager::get().free(MOT,SCR);

  }; 

} // namespace ChronusQ


