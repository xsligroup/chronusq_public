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

#include <response/base.hpp>

#include <cerr.hpp>
#include <util/mpi.hpp>
#include <util/timer.hpp>
#include <itersolver.hpp>

#ifdef CQ_ENABLE_MPI
#include <scalapackpp/block_cyclic.hpp>
#include <scalapackpp/types.hpp>
#include <scalapackpp/information.hpp>
#endif


namespace ChronusQ {


  template <typename T>
  struct RESPONSE_CONTRACTION {

    size_t nVec;
    size_t N;

    T* X;
    T* AX;

#ifdef CQ_ENABLE_MPI
    scalapackpp::scalapack_desc DescX;
    scalapackpp::scalapack_desc DescAX;
#endif

  };



  /**
   *  \brief Type dependent interface for Response classes
   *
   *  - Defines problem dimension / (optional) matrix storage
   *  - Stores the response (RESIDUE / FDR) results
   *  - Provides access to memory manager
   *  - Defines the preconditioner
   *  - Allocates memory
   *  - General "run" functions for response calculations
   *
   *  - Provides pure virtual interface for method specific 
   *  specializations
   *    - Form full matrix
   *    - Form RHS
   *    - Configuration options (set up)
   *    - Normalizing eigenvectors
   *    - Forming linear transformations
   *
   *
   */
  template <typename T>
  class ResponseTBase : public ResponseBase {

    template <class U> friend class ResponseTBase;

  protected:

    MPI_Comm      comm_;
    MPI_Comm      rcomm_;
    

    std::function< void(size_t,T,SolverVectors<T>&,SolverVectors<T>&) >       PC_;
    std::function< void(size_t,dcomplex,SolverVectors<dcomplex>&,SolverVectors<dcomplex>&) > cmplxPC_;

    std::function< void(size_t,SolverVectors<T>&,SolverVectors<T>&) >     nSPC_;


#ifdef CQ_ENABLE_MPI

    std::shared_ptr<scalapackpp::BlockCyclicDist2D> fullMatGrid_;
    scalapackpp::scalapack_desc    descFullMat_;

#endif


    size_t        nSingleDim_;
    T*            fullMatrix_;

    bool          alreadyRan_  = false;
    bool          hasResGuess_ = true;




    void writeMeta();

    template <typename U> void iterLinearTrans(size_t nVec, SolverVectors<U> &V, SolverVectors<U> &AV);


  public:

    ResidueResponseResults<T>     resResults;
    FDResponseResults<T,T>        fdrResults;
    FDResponseResults<T,dcomplex> dfdrResults;

    ResponseTBase( MPI_Comm c, ResponseType job, 
      T* fullMatrix = nullptr ) : ResponseBase(job),
        comm_(c), rcomm_(CreateRootComm(c)),
        nSingleDim_(0), 
        fullMatrix_(fullMatrix) {

        if(fullMatrix) genSettings.doFull = true; 

    }

    ~ResponseTBase() {

      reset();
      if( rcomm_ != MPI_COMM_NULL and rcomm_ != MPI_COMM_WORLD )
        MPICommFree(rcomm_);

    };

    ResponseTBase( const ResponseTBase &other ) : 
      ResponseBase(dynamic_cast<const ResponseBase&>(other)),
      comm_(other.comm_), rcomm_(CreateRootComm(other.comm_)),
      fullMatrix_(nullptr),
#ifdef CQ_ENABLE_MPI
      fullMatGrid_(other.fullMatGrid_),
      descFullMat_(other.descFullMat_),
#endif
      hasResGuess_(other.hasResGuess_),
      PC_(other.PC_),
      cmplxPC_(other.cmplxPC_) {

      if (other.fullMatrix_) {
        size_t fullMatSize = CQMemManager::get().getSize(other.fullMatrix_);
        fullMatrix_ = CQMemManager::get().malloc<T>(fullMatSize);
        std::copy_n(other.fullMatrix_, fullMatSize, fullMatrix_);
      }
    }

    inline void reset() {

      resResults. dealloc();
      fdrResults. dealloc();
      dfdrResults.dealloc();

      fdObs. dealloc();
      resObs.dealloc();

      alreadyRan_ = false;
      if (fullMatrix_) CQMemManager::get().free(fullMatrix_);

    }


    // Incore Matrix / Vector formation (Method specific)
    virtual size_t               getNSingleDim(const bool doTDA)= 0;
    virtual T*                   formFullMatrix()               = 0;
    virtual void                 formRHS       ()               = 0;
    virtual std::pair<size_t,T*> formPropGrad(ResponseOperator) = 0;
    virtual void                 configOptions()                = 0;
    virtual void                 eigVecNorm()                   = 0;
    virtual void                 resGuess(size_t, SolverVectors<T> &, size_t)   = 0;

    virtual void formLinearTrans( 
      std::vector<RESPONSE_CONTRACTION<double>>
    ) = 0;

    virtual void formLinearTrans( 
      std::vector<RESPONSE_CONTRACTION<dcomplex>>
    ) = 0;








    virtual T*  formFullFromMemory() { 

      // FIXME: This is ad hoc logic, assumes that if the cases when
      // matrix is not distributed that only the root process will be
      // calling this function
      bool needFullMat = not bool(fullMatrix_);
#ifdef CQ_ENABLE_MPI
      if( genSettings.isDist() ) MPIBCast(needFullMat,0,comm_);
#endif
      if(needFullMat) formFullMatrix(); 
      return fullMatrix_;

    };

    virtual void constructShifts() {

      fdrResults. shifts.clear();
      dfdrResults.shifts.clear();

      for(auto &omega : fdrSettings.bFreq) {
        fdrResults. shifts.emplace_back(omega);
        dfdrResults.shifts.emplace_back(omega,fdrSettings.dampFactor);
      }

    };

    // Memory allocation
    void allocResidueResults(); // Residue memory allocation
    void allocFDRResults();     // FDR memory allocation










    // Procedural functions


    inline virtual void run() {

      // Synchronize MPI processes
      MPI_Barrier(comm_);

      if( alreadyRan_ ) CErr("RESPONSE module only designed to be run once!");

      ProgramTimer::tick("Response Total");

      genSettings.matIsHer = genSettings.doTDA; // TDA always implies hermetian matrix
      configOptions();

      bool rootHasFullMat = (MPIRank(comm_) == 0) and fullMatrix_;
#ifdef CQ_ENABLE_MPI
      MPIBCast(rootHasFullMat,0,comm_);
#endif

      writeMeta();

      if( (genSettings.doFull or genSettings.formFullMat) and 
          not rootHasFullMat ) formFullMatrix();


      allocResults();

      bool isRoot = MPIRank(comm_) == 0;

      if( isRoot and genSettings.printLevel >= 0 )
        std::cout << "  ** STARTING RESPONSE CALCULATION\n\n";

      auto topResp = tick();

      if( genSettings.jobType == RESIDUE ) runResidue();
      else                                 runFDR();

      double botResp = tock(topResp);

      if( isRoot and genSettings.printLevel >= 0 )
        std::cout << "\n\n  ** RESPONSE CALCULATION FINISHED (" 
          << std::fixed << botResp << " s)\n";


      // Synchronize MPI processes
      MPI_Barrier(comm_);
      alreadyRan_ = true;

      ProgramTimer::tock("Response Total");

    };



    inline void runResidue() {

      if( genSettings.doFull ) runFullResidue();
      else                     runIterResidue();

      postDiagonalization();

      if( genSettings.evalProp ) residueProperties();

    }

    void runFullResidue();
    void runIterResidue();



    inline void runFDR() { 

      formRHS();
      constructShifts();

      if( genSettings.doFull ) runFullFDR();
      else                     runIterFDR();

      postLinearSolve();

      if( genSettings.evalProp ) fdrProperties();

      MPI_Barrier(comm_);
    }

    virtual void postLinearSolve(){


      ROOT_ONLY(comm_);

      bool noDamp = (fdrSettings.dampFactor == 0.) and 
                    not fdrSettings.forceDamp;

      size_t N = nSingleDim_;
      size_t nOmega = fdrSettings.bFreq.size();
      size_t nRHS   = fdrSettings.nRHS;

      // Write Solutions to disk
      if( savFile.exists() ) {

        savFile.safeWriteData("/RESP/FDR/OMEGA",&fdrSettings.bFreq[0],
          {nOmega});

        
        savFile.safeWriteData("/RESP/FDR/DAMP",&fdrSettings.dampFactor,
          {1});


        if( noDamp ) {

          savFile.safeWriteData("/RESP/FDR/RHS",fdrResults.RHS, {nRHS,N});

          savFile.safeWriteData("/RESP/FDR/SOLUTION",fdrResults.SOL,
            {nRHS*nOmega,N});

        } else {

          savFile.safeWriteData("/RESP/FDR/RHS",dfdrResults.RHS, {nRHS,N});

          savFile.safeWriteData("/RESP/FDR/SOLUTION",dfdrResults.SOL,
            {nRHS*nOmega,N});
        }

      };

    }


    virtual void postDiagonalization() {


      ROOT_ONLY(comm_);

      // Ensure proper orthogonality of eigenvectors
      eigVecNorm();


      // Write the eigensystem to disk
      if( savFile.exists() ) {

        savFile.safeWriteData("/RESP/RESIDUE/EIGENVALUES",resResults.W,
          {resSettings.nRoots});
        
        if( resSettings.needVR )
        savFile.safeWriteData("/RESP/RESIDUE/EIGENVECTORS",resResults.VR,
          {resSettings.nRoots,nSingleDim_});

      }

    }

    inline void runFullFDR() {

      if(not genSettings.isDist()) ROOT_ONLY(comm_);

      bool noDamp = (fdrSettings.dampFactor == 0.) and 
                    not fdrSettings.forceDamp;

      if( noDamp ) runFullFDR(fdrResults);
      else         runFullFDR(dfdrResults);

    }

    template <typename U>
    void runFullFDR(FDResponseResults<T,U> &results);



    inline void runIterFDR() {

      if((not genSettings.isDist()) and genSettings.formFullMat) 
        ROOT_ONLY(comm_);

      bool noDamp = (fdrSettings.dampFactor == 0.) and 
                    not fdrSettings.forceDamp;

      if( noDamp ) runIterFDR(fdrResults,  PC_);
      else         runIterFDR(dfdrResults, cmplxPC_);

    }

    template <typename U>
    void runIterFDR(FDResponseResults<T,U> &results,
                    std::function< void(size_t,U,SolverVectors<U>&,SolverVectors<U>&) > &preCond );



    // Property Evaluation

    template <typename U>
    std::map<ResponseOperator, U*> 
      evalProperties(std::vector<ResponseOperator> ops, size_t nVec, U* V);



    // Residue properties
    inline void residueProperties() {

      ROOT_ONLY(comm_);

      // If we're doing a residue calculation, make sure we evaluate all
      // of the properties
      //genSettings.aOps = AllOps;

      ProgramTimer::timeOp("Property Eval", [&](){
        residueTMoments();
        residueObservables();
      });
    }

    void residueTMoments();
    void residueObservables();


    // FDR properties
    inline void fdrProperties() {


      ROOT_ONLY(comm_);

      bool noDamp = (fdrSettings.dampFactor == 0.) and 
                    not fdrSettings.forceDamp;

      if( noDamp ) {
        fdrRF(fdrResults);
        // fdrObservables(fdrResults); // Can we define for undamped?
      } else {
        fdrRF(dfdrResults);
        fdrObservables(dfdrResults);
      }

    }

    template <typename U>
    void fdrRF(FDResponseResults<T,U> &results);

    template <typename U>
    void fdrObservables(FDResponseResults<T,U> &results);






    // Print results / properties
      
    inline void printResults(std::ostream &out) {

      ROOT_ONLY(comm_);

      if( genSettings.jobType == RESIDUE ) printResidueResults(out);
      else                                 printFDRResults(out);

    };


    // Print Residue results
          
    inline void printResidueResults(std::ostream &out) {

      ROOT_ONLY(comm_);
      out << "\n\n\nRESIDUE RESPONSE RESULTS";

      printResMO(out);          // Print residual MO contributions
      if( genSettings.evalProp ) {
        printTMoments(out);       // Print transition moments at poles
        printResObservables(out); // Print observables
      }

    }

    virtual void printResMO(std::ostream &out) = 0; // Method specific
    virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
      std::vector<std::pair<std::string,double *>> data, double* VL, 
      double* VR) = 0; // Method specific
    virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
      std::vector<std::pair<std::string,double *>> data, dcomplex* VL, 
      dcomplex* VR) = 0; // Method specific

    void         printTMoments(std::ostream &out);
    void         printResObservables(std::ostream &out);


    // Print FDR Results (NYI)
      
    inline void printFDRResults(std::ostream &out) {

      if(not genSettings.evalProp) return;

      ROOT_ONLY(comm_);
      bool noDamp = (fdrSettings.dampFactor == 0.) and 
                    not fdrSettings.forceDamp;

      out << "\n\n\nFREQUENCY DEPENDENT RESPONSE RESULTS";

      if( noDamp ) printRF(fdrResults, out);
      else {         
        printRF(dfdrResults,out);
        printFDObservables(out);
      }

    }

    template <typename U>
    void printRF(FDResponseResults<T,U> &results, std::ostream &out);

    void printFDObservables(std::ostream &out);



    // Getters
    T*                    fullMatrix()  const { return fullMatrix_; }

#ifdef CQ_ENABLE_MPI
    auto* fullMatGrid() const { return fullMatGrid_.get(); }
#endif


  };


  template <typename T>
  void ResponseTBase<T>::writeMeta() {

    // Write out the meta data to disk

    if( not savFile.exists() ) return;

    int isFull = int(genSettings.doFull);
    savFile.safeWriteData( "/RESP/DOFULL", &isFull, {1} );

    if( genSettings.jobType == RESIDUE ) {

      savFile.safeWriteData( "/RESP/RESIDUE/DEMIN",   
          &resSettings.deMin, {1} );
      savFile.safeWriteData( "/RESP/RESIDUE/GPLHR_M", 
          &resSettings.gplhr_m, {1});
      savFile.safeWriteData( "/RESP/RESIDUE/GPLHR_SIGMA", 
          &resSettings.gplhr_sigma, {1});

    }
  }


  /**
   *  \brief General implementation of linear transformation required
   *  for iterative solvers (GMRES, GPLHR, etc)
   */
  template <typename T>
  template <typename U> 
  void ResponseTBase<T>::iterLinearTrans(size_t nVec, SolverVectors<U> &V, SolverVectors<U> &AV) {

    std::vector<RESPONSE_CONTRACTION<U>> cList(1);
    cList.back().nVec = nVec;
    cList.back().N    = nSingleDim_;
    cList.back().X    = tryGetRawVectorsPointer(V);
    cList.back().AX   = tryGetRawVectorsPointer(AV);

#ifdef CQ_ENABLE_MPI
    int64_t MLoc, NLoc;
    if( genSettings.isDist() ) {
      //std::tie(MLoc,NLoc) = fullMatGrid_->getLocalDims(nSingleDim_,nVec);
      std::tie(MLoc,NLoc) = fullMatGrid_->get_local_dims(nSingleDim_, nVec);
      if( MLoc and NLoc ) {

        cList.back().X  = CQMemManager::get().malloc<U>(MLoc*NLoc);
        cList.back().AX = CQMemManager::get().malloc<U>(MLoc*NLoc);

      } else {

        cList.back().X  = nullptr;
        cList.back().AX = nullptr;

      }

      cList.back().DescX  = fullMatGrid_->descinit_noerror(nSingleDim_,nVec,MLoc);
      cList.back().DescAX = cList.back().DescX;

      fullMatGrid_->scatter(nSingleDim_,nVec,tryGetRawVectorsPointer(V),nSingleDim_,cList.back().X,
        MLoc,0,0);
    }
#endif
    
    this->formLinearTrans( cList );

#ifdef CQ_ENABLE_MPI
    if( genSettings.isDist() ) {

      fullMatGrid_->gather(nSingleDim_,nVec,tryGetRawVectorsPointer(AV),nSingleDim_,cList.back().AX,
        MLoc,0,0);

      if( cList.back().X  ) CQMemManager::get().free(cList.back().X );
      if( cList.back().AX ) CQMemManager::get().free(cList.back().AX);

    }
#endif

  }



























  template <class Reference> 
  class ResponseRefBase : 
    public ResponseTBase<typename Reference::value_type> {

  protected:

    typedef typename Reference::value_type T;
    std::shared_ptr<Reference> ref_;

  public:

    ResponseRefBase( MPI_Comm c, ResponseType job, 
      std::shared_ptr<Reference> ref, T* fullMatrix = nullptr ) : 
      ResponseTBase<T>(c,job,fullMatrix), ref_(ref){ }

    ResponseRefBase( const ResponseRefBase &other ) : 
      ResponseTBase<T>(dynamic_cast<const ResponseTBase<T>&>(other)),
      ref_(other.ref_){ }


    std::shared_ptr<Reference> ref() { return ref_; };

  };

}; // namespace ChronusQ

