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
#include <response.hpp>
#include <physcon.hpp>

#include <cqlinalg/blas3.hpp>
#include <cqlinalg/ortho.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/eig.hpp>
#include <cqlinalg/solve.hpp>
#include <cqlinalg/blasutil.hpp>

namespace ChronusQ {


  enum MOR_TARGET {
    
    OPA_CROSS_SECTION_EDA,
    ECD_CROSS_SECTION_MD

  };

  struct MORSettings {

    size_t nModel    = 16;
    size_t nModelMax = 256;
    double convCrit  = 1e-2;
    bool   doRefine  = false;
    bool   getEig    = false;
    bool   doRelErr  = true;

    MOR_TARGET target = OPA_CROSS_SECTION_EDA;

  };


  struct MORSpecBase {

    MORSettings morSettings;

  };

  template <class Reference>
  class MORSpec : public ResponseTBase<dcomplex>, public MORSpecBase {

    typedef dcomplex T;
    typedef typename Reference::value_type RefT;
    typedef typename Reference::ints_type  IntsT;

    std::shared_ptr<Reference> ref_;
    dcomplex*      modelBasis1_    = nullptr;
    dcomplex*      modelBasis2_    = nullptr;
    dcomplex*      modelBasis1_LT_ = nullptr;
    dcomplex*      modelBasis2_LT_ = nullptr;
    dcomplex*      ritzVecR_       = nullptr;
    dcomplex*      ritzVecL_       = nullptr;

    double*        resNorms_       = nullptr;

    PolarizationPropagator<Reference> respFactory_;
    typename Reference::value_type* respFullMatrix_ = nullptr;

    inline size_t getNSingleDim(const bool doTDA = false) override {

      return this->morSettings.nModel;

    }

  public:

    T*                   formFullMatrix() override;
    void                 formRHS       () override;
    std::pair<size_t,T*> formPropGrad(ResponseOperator) override;
    void                 configOptions() override;
    void                 eigVecNorm() override                   {};
    void                 resGuess(size_t, SolverVectors<T> &, size_t) override { CErr(); };



    void printResMO(std::ostream &out)  override {
     

      double * W_print = this->resResults.W ;

      size_t N = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
      size_t nRoots = this->resSettings.nRoots;
      size_t nModel = this->morSettings.nModel;

      dcomplex *VR = ritzVecR_;
      dcomplex *VL = ritzVecL_ ? ritzVecL_ : VR + N/2;

      respFactory_.printResMO(out,nRoots,W_print,
        {{" f    = ",this->resObs.oscStrength},
         {"|R|   = ",resNorms_}}, VL,VR);

    };

    virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
      std::vector<std::pair<std::string,double *>> data, double* VL, 
      double* VR) override { }
    virtual void printResMO(std::ostream &out, size_t nRoots, double *W,
      std::vector<std::pair<std::string,double *>> data, dcomplex* VL, 
      dcomplex* VR) override { }

    void constructShifts() override {

      respFactory_.constructShifts();
      this->fdrResults. shifts.clear();
      this->dfdrResults.shifts.clear();
      std::copy(
        respFactory_.fdrResults.shifts.begin(),
        respFactory_.fdrResults.shifts.end(),
        std::back_inserter(this->fdrResults.shifts)
      );
      std::copy(
        respFactory_.dfdrResults.shifts.begin(),
        respFactory_.dfdrResults.shifts.end(),
        std::back_inserter(this->dfdrResults.shifts)
      );

    };


    std::vector<std::vector<double>> 
      generateShiftLevels(size_t nLevel, size_t nModel, double wMin, 
        double wMax) {

      std::vector<std::vector<double>> shifts;

      size_t nCurMax = nModel;
      for( auto iLevel = 0; iLevel < nLevel; iLevel++ ) {

        std::vector<double> tmp(nCurMax,0.);

        tmp[0] = wMin;
        for(auto k = 1; k < tmp.size(); k++){
          tmp[k] = tmp[k-1] + (wMax - wMin) / (tmp.size() - 1);
        }

        if(nCurMax == nModel) shifts.emplace_back(tmp);
        else {
          std::vector<double> tmp2;

          for(auto j = 1; j < tmp.size(); j+=2)
            tmp2.emplace_back(tmp[j]);

          shifts.emplace_back(tmp2);
          
        }

        nCurMax = 2*nCurMax - 1;
      }


      return shifts;

    }

    std::vector<double> selectShifts(std::vector<double> &candidate,
        std::vector<double> &omega, std::vector<double> &diff, 
        double tol) {

      std::vector<double> shifts;

      size_t iOmegaSt  = 0;
      auto   itOmegaSt = omega.begin();

      for(auto k = 0; k < candidate.size(); k++) {

        double boundary = candidate[k];
        if( k == candidate.size() - 1 ) 
          boundary = omega.back();
        else                            
          boundary += 0.5 * (candidate[k+1] - candidate[k]);

        auto itOmegaEnd = std::find_if(itOmegaSt,omega.end(),
            [&](double a){ return a > boundary; });

        size_t iOmegaEnd = std::distance(omega.begin(),itOmegaEnd);

        bool needShift = 
          std::any_of( diff.begin() + iOmegaSt, diff.begin() + iOmegaEnd,
            [&](double a){ return std::abs(a) > tol; } );

        if( needShift ) shifts.emplace_back(candidate[k]);

        iOmegaSt  = iOmegaEnd;
        itOmegaSt = itOmegaEnd;

      };

      return shifts;

    }


    void formLinearTrans(
      std::vector<RESPONSE_CONTRACTION<double>> x
    ) override { CErr(); }
    void formLinearTrans( 
      std::vector<RESPONSE_CONTRACTION<dcomplex>> x
    ) override { CErr(); }


    MORSpec( MPI_Comm c, std::shared_ptr<Reference> ref ) : 
      ResponseTBase<dcomplex>(c,FDR), ref_(ref),
      respFactory_(c,FDR,ref) { 
    
      this->fdrSettings.dampFactor  = 0.01; // Default the damping factor
      this->genSettings.doFull      = true;
      this->genSettings.printLevel  = -1;

    }

    ~MORSpec() {
      if(modelBasis1_) CQMemManager::get().free(modelBasis1_);
      if(modelBasis1_LT_) CQMemManager::get().free(modelBasis1_LT_);
    }



    PolarizationPropagator<Reference>& respFactory(){ return respFactory_; }



    void formModelBasis(){ 

      bool isRoot = MPIRank(this->comm_) == 0;

      if( isRoot ) {
        std::cout << "\nSTARTING CONSTRUCTION OF MODEL BASIS\n\n";

        // Print settings
        std::cout << "  * MOR SETTINGS\n";
        std::cout << "    * NMODEL START = " << this->morSettings.nModel 
          << "\n";
        if( this->morSettings.doRefine ) {

          std::cout << "    * MOR WILL PROCEDE ADAPTIVELY\n";
          std::cout << "      * NMODELMAX = " 
            << this->morSettings.nModelMax << "\n";

        } else {

          std::cout << 
            "    * MOR WILL PROCEDE WITH A FIXED NUMBER OF SHIFTS\n";

        }

        std::cout << std::endl << std::endl;

      }

      // Form the full matrix if requested 
      if( respFactory_.genSettings.formFullMat ) 
        respFullMatrix_ = respFactory_.formFullMatrix();

      // Get problem dimension
      size_t N = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
      size_t NRHS = respFactory_.fdrSettings.nRHS;
      size_t nModelMax = this->morSettings.nModelMax;
      size_t nEval     = respFactory_.fdrSettings.bFreq.size();
      





      if( isRoot ) {
        // Allocate max amount of space that could be required for MOR
        std::cout << "  * ALLOCATING MAX SPACE FOR MOR BASIS\n";
        modelBasis1_ = 
          CQMemManager::get().malloc<dcomplex>(N*NRHS*nModelMax);
        std::fill_n(modelBasis1_,N*NRHS*nModelMax,dcomplex(0.));
        modelBasis1_LT_ = 
          CQMemManager::get().malloc<dcomplex>(N*NRHS*nModelMax);
        std::fill_n(modelBasis1_LT_,N*NRHS*nModelMax,dcomplex(0.));
      }








      // Generate Shift Levels
      double wMin = this->fdrSettings.bFreq.front();
      double wMax = this->fdrSettings.bFreq.back();

      if( isRoot ) std::cout << "  * GENERATING SHIFT LEVELS\n";

      auto shiftLevels = 
        generateShiftLevels(this->morSettings.doRefine ? 6 : 1,
          this->morSettings.nModel,wMin,wMax);

      if( isRoot )
        std::cout << "    * GENERATED " << shiftLevels.size() 
                  << " SHIFT LEVELS\n\n";








      // Loop over shift levels
      this->morSettings.nModel = 0;
      std::vector<std::vector<double>> targetFunction;
      bool refineConv(false);
      bool exitMOR(false);
      size_t nModelShift(0);

      for(auto iLevel = 0; iLevel < shiftLevels.size(); iLevel++ ) {



        MPI_Barrier(this->comm_); // Sync processes at top of MOR

        auto & shiftLevel = shiftLevels[iLevel];
        size_t nShift = shiftLevel.size();

        // Synchronize selected shifts
        if( MPISize(this->comm_) > 1 ) {

          MPIBCast(nShift,0,this->comm_); // Send out nShift

          // Clear out shifts on non-root process
          if( not isRoot ) {
            shiftLevel.clear();
            shiftLevel.resize(nShift);
          }

          // Send out shifts
          MPIBCast(&shiftLevel[0],nShift,0,this->comm_);

        }


        nModelShift += nShift;

        // Output the currently selected shifts
        if( isRoot ) {

          std::cout << "  * MOR SHIFT LEVEL " << iLevel << " HAS SELECTED " 
            << nShift << " NEW SHIFTS\n";
          for(auto iShift = 0; iShift < nShift; iShift++)
            std::cout << "    SHIFT " << std::setw(5) << iShift << " = "
              << std::scientific << std::setprecision(8) 
              << shiftLevel[iShift] << std::endl;

          std::cout << "\n\n";

        }


        if( isRoot )
          std::cout << "  * RUNNING FULL PROBLEM FOR SELECTED SHIFTS\n\n\n"; 

        // Solve the full problem for the selected shifts
        PolarizationPropagator<Reference> subProblem(respFactory_);

        // Turn off print for sub problem
        subProblem.genSettings.printLevel = 0;

        subProblem.fdrSettings.bFreq = shiftLevel;
        subProblem.run();

        if( isRoot ) std::cout << "\n\n\n";




        // Copy over the solutions to the model basis
        if( isRoot ) {
          std::cout << "  * COPYING FULL SOLUTION INTO MODEL BASIS\n";
          std::copy_n(subProblem.dfdrResults.SOL,N*NRHS*nShift,
              modelBasis1_ + this->morSettings.nModel*N);
        }

      //prettyPrintSmart(std::cout,"SOLUTION",subProblem.dfdrResults.SOL,
      //  N,NRHS*nShift,N);






        // Orthonormalize the basis
        size_t nOrtho = 0;

        if( isRoot ) {
          std::cout << "  * ORTHONORMALIZING MODEL BASIS\n";

          size_t nVec = this->morSettings.nModel + nShift*NRHS;
          size_t nUse = std::min(N,nVec);

          bool doSVD = nVec >= N;
          doSVD = true;
     
          if( doSVD ) {
            std::cout << "    * SENDING " << nVec << " VECTORS TO SVD\n";

            double *SVal = CQMemManager::get().malloc<double>(nUse);
            std::fill_n(SVal,nUse,0.);
            dcomplex *DUMMY = nullptr;

            lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec,
              N,nVec,modelBasis1_,N,SVal,DUMMY,1,DUMMY,1);

            double orthoTol = 1e-12 * SVal[0] * std::max(N,nVec);
            for(auto k = 0; k < nUse; k++) if(SVal[k] > orthoTol) nOrtho++;

       
            //prettyPrintSmart(std::cout,"S",SVal,nUse,1,nUse);

            CQMemManager::get().free(SVal);

            std::cout << "    * SVD PRODUCED " << nOrtho 
              << " ORTHONORMAL VECTORS\n";
          } else {
          
            std::cout << "    * SENDING " << nVec << " VECTORS TO QR" 
              << std::endl;
            nOrtho = nVec;
            QR(N,nVec,modelBasis1_,N);

          }

        }
       //prettyPrintSmart(std::cout,"MODEL BASIS AFTER",modelBasis1_,N,
       //    nOrtho,N);

        // Broadcast the number of orthonormal vectors
        if( MPISize(this->comm_) > 1 ) MPIBCast(nOrtho,0,this->comm_);


        // Solve the reduced space problem
        this->morSettings.nModel = nOrtho;
        this->nSingleDim_  = this->morSettings.nModel;
        this->reset();
        if( iLevel > 0 ) { 
          if(this->fullMatrix_) CQMemManager::get().free(this->fullMatrix_);
          this->fullMatrix_ = nullptr;
        }

        if( isRoot ) {
          std::cout << "\n  * RUNNING REDUCED SPACE RESPONSE PROBLEM\n";
          std::cout << "    * NMODEL = " << this->morSettings.nModel << "\n";
          std::cout << "    * NEVAL  = " << nEval << "\n\n";
        }

        ResponseTBase<dcomplex>::run();


        if( isRoot ) { // Do reduced dimensional stuff on ROOT





          std::cout << "  * MAKING COPY OF TARGET FUNCTION FOR SHIFT LEVEL "
            << iLevel << "\n";

          // Save a copy of the target function
          targetFunction.emplace_back();
          if( this->morSettings.target == OPA_CROSS_SECTION_EDA )
            std::copy_n(this->fdObs.opaCross_eda,nEval,
              std::back_inserter(targetFunction.back()));
          else
            CErr("MOR TARGETS OTHER THAN OPA CROSS SECTION NYI");

          // Check for convergence if requested
          if( iLevel > 0 ) {

            std::cout << "  * COMPARING TARGET FUNCTIONS FOR SHIFT LEVELS "
              << std::setw(5) << iLevel-1 << " AND " 
              << std::setw(5) << iLevel << "\n";



            double maxTarget = 
              std::abs(*std::max_element(targetFunction.back().begin(),
                                         targetFunction.back().end(),
                                         [&](double x, double y) {
                                           return
                                             std::abs(x) < std::abs(y);
                                         }));


            // Evaluate the differences
            std::vector<double> diff(nEval);
            for(auto k = 0; k < nEval; k++) {

              double w = fdrSettings.bFreq[k];

              // Compare Tr[Alpha] as opposed to Sigma
              diff[k] = std::abs(
                          (targetFunction[iLevel]  [k]/w) -
                          (targetFunction[iLevel-1][k]/w)
                        );

              if(this->morSettings.doRelErr)
                diff[k] /= std::abs(targetFunction[iLevel-1][k]/w);
              else
                diff[k] /= maxTarget;

            }

            // Check the maximum relative difference and break if converged
            double maxDiff = *std::max_element(diff.begin(),diff.end());

            std::cout << "    * DIFF METHOD ";
            if( this->morSettings.doRelErr )
              std::cout << "RELATIVE ERROR";
            else
              std::cout << "NORMALIZED ERROR";
            std::cout << "\n";

            std::cout << "    * CONV CRIT                      = " 
              << std::scientific <<  this->morSettings.convCrit << "\n";
            std::cout << "    * MAX TARGET FUNCTION DIFFERENCE = " 
              << maxDiff << "\n\n";

            if( maxDiff < this->morSettings.convCrit ) refineConv = true;

            // If not converged, select new shifts from the next level
            size_t nNew = 0;
            if( (iLevel != shiftLevels.size() - 1) and not refineConv ) {

              std::cout << "  * SELECTING SHIFTS FOR NEXT SHIFT LEVEL\n";
              std::cout << "    * NCANDIDATE = " << shiftLevels[iLevel+1].size()
                << "\n";

              shiftLevels[iLevel+1] = 
                selectShifts(
                  shiftLevels[iLevel+1],
                  this->fdrSettings.bFreq,
                  diff, this->morSettings.convCrit
                );

              std::cout << "    * NSELECT    = " << shiftLevels[iLevel+1].size()
                << "\n";

              nNew = shiftLevels[iLevel+1].size();

              // If no shifts were selected, then I guess we're kinda converged
              if( nNew == 0 ) refineConv = true;

              nNew = NRHS*nNew + this->morSettings.nModel;

            }

      
            if( nNew > nModelMax ) {
              std::cout << "\n  * SELECTED SHIFTS WOULD GROW THE MODEL BASIS "
                << "PAST THE SPECIFIED NMODELMAX!";
              std::cout << "\n\n";
            }

            // Break if converged or over max basis size
            exitMOR = refineConv or (nNew > nModelMax);

          }

          std::cout << "\n\n";

        }

        if( MPISize(this->comm_) > 1 ) MPIBCast(exitMOR,0,this->comm_);
        if( exitMOR ) break;

      }

      if( isRoot ) {
        if( refineConv and this->morSettings.doRefine ) {
          std::cout << "  ******** ITERATIVE MOR CONVERGED ********\n\n";
          std::cout << "    * NMODEL  = " << this->morSettings.nModel << "\n";
          std::cout << "    * NLINEAR = " << nModelShift << "\n\n";
          std::cout << "  *****************************************\n\n";
        } else if (this->morSettings.doRefine) {
          std::cout << "  ******** ITERATIVE MOR FAILED TO CONVERGE ********\n\n";
          std::cout << "    * NMODEL  = " << this->morSettings.nModel << "\n";
          std::cout << "    * NLINEAR = " << nModelShift << "\n\n";
          std::cout << "  **************************************************\n\n";
        }
      }

      if( MPIRank(this->comm_) == 0 ) {


        if( this->savFile.exists() ) {

          this->savFile.safeWriteData("/MOR/MODELBASIS",modelBasis1_,
              { this->morSettings.nModel, N });

        }



      }

      MPI_Barrier(this->comm_);

    };
    
    inline void run() override {

      // Input RESP arguements are used to populate internal options,
      // this copies them over to respFactory_
      respFactory_.fdrSettings = this->fdrSettings;
      respFactory_.resSettings = this->resSettings;
      respFactory_.genSettings = this->genSettings;

      // MOR ALWAYS has these options
      this->genSettings.doFull          = true;
      this->genSettings.formFullMat     = true;
      this->genSettings.distMatFromRoot = false;
      this->genSettings.formMatDist     = false;

      // Turn on print for resp factory
      respFactory_.genSettings.printLevel = 1;

      respFactory_.configOptions();

      formModelBasis();
      //modelBasis1_ = CQMemManager::get().malloc<dcomplex>(N*N);
      //std::fill_n(modelBasis1_,N*N,0.);
      //for(auto k = 0; k < N; k++)
      //  modelBasis1_[k*(N+1)] = 1.;
      //this->morSettings.nModel = N;

      if( this->morSettings.getEig )
//      extractEig1();
        extractEig2();

    };


    void extractEig1() {

      if( MPISize(this->comm_) > 1 ) CErr("Stupid eigenvector extraction doesn't work with MPI");

      size_t N = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
      size_t nModel = this->morSettings.nModel;
      size_t nRoots = nModel;

      this->genSettings.jobType = RESIDUE;
      this->genSettings.evalProp = false;
      this->resSettings.nRoots = nModel;


      this->reset();
      if(this->fullMatrix_) CQMemManager::get().free(this->fullMatrix_);
      this->fullMatrix_ = nullptr;

      ResponseTBase<dcomplex>::run();

      // Determine if roots are actually roots
      dcomplex *AV = CQMemManager::get().malloc<dcomplex>(
          nModel*N);
      std::fill_n(AV,nModel*N,dcomplex(0.));
      dcomplex *RES = CQMemManager::get().malloc<dcomplex>(
          nModel*N);
      std::fill_n(RES,nModel*N,dcomplex(0.));

      // Compute AV
      std::vector<RESPONSE_CONTRACTION<T>> cList(1);
      cList.back().nVec   = nModel;
      cList.back().N      = N;
      cList.back().X      = modelBasis1_;
      cList.back().AX     = AV;

      respFactory_.formLinearTrans(cList);


      // Compute the residual
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nModel,nModel,dcomplex(1.),modelBasis1_,N,
          this->resResults.VR,nModel,dcomplex(0.),RES,N);
      for(auto k = 0; k < nModel; k++)
        blas::scal(N,-dcomplex(this->resResults.W[k]),RES + k*N,1);

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nModel,nModel,dcomplex(1.),AV,N,
          this->resResults.VR,nModel,dcomplex(1.),RES,N);

      double * rNorms = CQMemManager::get().malloc<double>(nModel);
      for(auto k = 0; k < nModel; k++)
        rNorms[k] = blas::nrm2(N,RES+k*N,1);


      // Sort the vectors on residual
      std::vector<int> indx(nModel,0);
      std::iota(indx.begin(),indx.end(),0);

      std::stable_sort(indx.begin(),indx.end(),
        [&](const int &a, const int &b) {
          return rNorms[a] > rNorms[b];
        }
      );


      for(auto k = 0; k < (nModel - 1); k++) {

        int ind = indx[k];
        while( ind < k ) ind = indx[ind];

        std::swap(rNorms[k],rNorms[ind]);
        std::swap(this->resResults.W[k],this->resResults.W[ind]);
        blas::swap(nModel,this->resResults.VR + k*nModel  , 1, 
                    this->resResults.VR + ind*nModel, 1);

      }

      // Find first small residual
      double eigResTol = 1e-6;
      double * erstKlein = 
        std::find_if(rNorms,rNorms + nModel,
          [&](double x){ return x < eigResTol; });

      // Remove vectors with large residual
      if( erstKlein != rNorms ) {

        size_t k = std::distance(rNorms,erstKlein);

        std::cout << "  * FOUND " << k << " RITZ VECTORS WHICH ARE NOT "
          << "TRUE EIGENVECTORS\n";

        nRoots = nModel - k;
        size_t nCopy = nRoots * nModel;
        std::copy_n(this->resResults.VR + k*nModel,nCopy, 
          this->resResults.VR);

        std::copy_n(this->resResults.W  + k, nRoots, this->resResults.W);

      }




      // Remove roots with negative eigenvalues

      indx.resize(nRoots);
      std::iota(indx.begin(),indx.end(),0);

      std::stable_sort(indx.begin(),indx.end(),
        [&](const int &a, const int &b) {
          return this->resResults.W[a] < this->resResults.W[b];
        }
      );

      for(auto k = 0; k < (nRoots - 1); k++) {

        int ind = indx[k];
        while( ind < k ) ind = indx[ind];

        std::swap(this->resResults.W[k],this->resResults.W[ind]);
        blas::swap(nModel,this->resResults.VR + k*nModel  , 1, 
                    this->resResults.VR + ind*nModel, 1);

      }



      // Find first small residual
      double * erstPositiv = 
        std::find_if(this->resResults.W,this->resResults.W + nRoots,
          [&](double x){ return x > 0.; });

      // Remove vectors with large residual
      if( erstPositiv != this->resResults.W ) {

        size_t k = std::distance(this->resResults.W,erstPositiv);

        std::cout << "  * FOUND " << k << " RITZ VECTORS WHICH CORRESPOND "
          << "TO DE-EXCITATIONS\n";

        nRoots = nRoots - k;
        size_t nCopy = nRoots * nModel;
        std::copy_n(this->resResults.VR + k*nModel,nCopy, 
          this->resResults.VR);

        std::copy_n(this->resResults.W  + k, nRoots, this->resResults.W);

      }



      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nRoots,nModel,dcomplex(1.),modelBasis1_,N,
          this->resResults.VR,nModel,dcomplex(0.),RES,N);
      for(auto k = 0; k < nRoots; k++)
        blas::scal(N,-dcomplex(this->resResults.W[k]),RES + k*N,1);

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nRoots,nModel,dcomplex(1.),AV,N,
          this->resResults.VR,nModel,dcomplex(1.),RES,N);


      for(auto k = 0; k < nRoots; k++) {
        rNorms[k] = blas::nrm2(N,RES+k*N,1);
      }

      


      // Orthonormalize the basis
      // Compute the inner product
      dcomplex * SCR = CQMemManager::get().malloc<dcomplex>(
          nModel*nModel);
      std::fill_n(SCR,nModel*nModel,dcomplex(0.));
      dcomplex * SCR2 = CQMemManager::get().malloc<dcomplex>(
          nModel*nModel);
      std::fill_n(SCR2,nModel*nModel,dcomplex(0.));

      // SCR = V* S V
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(1.),
        modelBasis1_,N,modelBasis1_,N,dcomplex(0.),SCR,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(-1.),
        modelBasis1_+N/2,N,modelBasis1_+N/2,N,dcomplex(1.),SCR,nModel);

      // SCR = c* SCR c
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nModel,nRoots,nModel,dcomplex(1.),SCR,nModel,
          this->resResults.VR,nModel,dcomplex(0.),SCR2,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nRoots,nRoots,nModel,dcomplex(1.),this->resResults.VR,
          nModel,SCR2,nModel,dcomplex(0.),SCR,nRoots);

    //prettyPrintSmart(std::cout,"SCR",SCR,nRoots,nRoots,nRoots);

      for(auto k = 0; k < nRoots; k++)
        blas::scal(nModel,dcomplex(1./std::sqrt(SCR[k*(nRoots+1)])),
            this->resResults.VR + k*nModel,1);



      this->resSettings.nRoots = nRoots;
      this->genSettings.evalProp = true;
      residueProperties();



      /*
      // SCR = V* S V
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(1.),
        modelBasis1_,N,modelBasis1_,N,dcomplex(0.),SCR,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(-1.),
        modelBasis1_+N/2,N,modelBasis1_+N/2,N,dcomplex(1.),SCR,nModel);

      // SCR = c* SCR c
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nModel,nRoots,nModel,dcomplex(1.),SCR,nModel,
          this->resResults.VR,nModel,dcomplex(0.),SCR2,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nRoots,nRoots,nModel,dcomplex(1.),this->resResults.VR,
          nModel,SCR2,nModel,dcomplex(0.),SCR,nRoots);

      prettyPrintSmart(std::cout,"SCR After",SCR,nRoots,nRoots,nRoots);
      */

      /*

      double * O = CQMemManager::get().malloc<double>(nRoots);
      HermetianEigen('V','L',nRoots,SCR,nRoots,O);

      prettyPrintSmart(std::cout,"O",O,nRoots,1,nRoots);

      // c w O^-1/2
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nModel,nRoots,nRoots,dcomplex(1.),this->resResults.VR,
        nModel,SCR,nModel,dcomplex(0.),SCR2,nModel);

      prettyPrintSmart(std::cout,"EV",SCR,nRoots,nRoots,nRoots);

      std::copy_n(SCR2,nModel*nRoots,this->resResults.VR);

      for(auto k = 0; k < nRoots; k++)
        blas::scal(nModel,dcomplex(1./std::sqrt(O[k])),
            this->resResults.VR + k*nModel,1);

      // SCR = V* S V
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(1.),
        modelBasis1_,N,modelBasis1_,N,dcomplex(0.),SCR,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N/2,dcomplex(-1.),
        modelBasis1_+N/2,N,modelBasis1_+N/2,N,dcomplex(1.),SCR,nModel);

      // SCR = c* SCR c
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nModel,nRoots,nModel,dcomplex(1.),SCR,nModel,
          this->resResults.VR,nModel,dcomplex(0.),SCR2,nModel);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nRoots,nRoots,nModel,dcomplex(1.),this->resResults.VR,
          nModel,SCR2,nModel,dcomplex(0.),SCR,nRoots);

      prettyPrintSmart(std::cout,"SCR AFTER",SCR,nRoots,nRoots,nRoots);
      */

    }


    void extractEig2() {

      PolarizationPropagator<SingleSlater<RefT,IntsT>> &ss = 
        dynamic_cast<PolarizationPropagator<SingleSlater<RefT,IntsT>>&>(respFactory_);
    
      bool isRoot = MPIRank(this->comm_) == 0;

      if( isRoot )
        std::cout << "  * EXTRACTING EIGENSYSTEM FROM MOR BASIS" << std::endl;

      if( ss.doReduced ) extractEigRed();
      else               extractEigFull();

    }

    void extractEigRed(){ 
    
      bool isRoot = MPIRank(this->comm_) == 0;
      size_t N = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
      size_t nModel = this->morSettings.nModel;
      size_t nRoots = nModel;

      PolarizationPropagator<SingleSlater<RefT,IntsT>> &ss = 
        dynamic_cast<PolarizationPropagator<SingleSlater<RefT,IntsT>>&>(respFactory_);


      // Reset the Response job
      this->reset();
      if(this->fullMatrix_) CQMemManager::get().free(this->fullMatrix_);
    
      // Allocate full matrix
      this->fullMatrix_ = isRoot ?
        CQMemManager::get().malloc<dcomplex>(nModel*nModel) : 
        nullptr ;
      if( isRoot ) std::fill_n(this->fullMatrix_,nModel*nModel,dcomplex(0.));
    

      // Setup contraction
      std::vector<RESPONSE_CONTRACTION<T>> cList(1);
      cList.back().nVec   = nModel;
      cList.back().N      = N;
      cList.back().X      = modelBasis1_;
      cList.back().AX     = modelBasis1_LT_;


      // If distributed full matrix, allocate some scratch space
      int64_t MLoc, NLoc;
#ifdef CQ_ENABLE_MPI
      if( respFactory_.genSettings.isDist() ) {

        auto * grid = respFactory_.fullMatGrid();
        std::tie(MLoc,NLoc) = grid->get_local_dims(N,nModel);

        if( MLoc and NLoc ) {

          cList.back().X  = CQMemManager::get().malloc<T>(MLoc*NLoc);
          cList.back().AX = CQMemManager::get().malloc<T>(MLoc*NLoc);

        } else {

          cList.back().X  = nullptr;
          cList.back().AX = nullptr;

        }

        cList.back().DescX  = grid->descinit_noerror(N,nModel,MLoc);
        cList.back().DescAX = cList.back().DescX;

        // Scatter modelBasis1_ to the BLACS grid
        grid->scatter(N,nModel,modelBasis1_,N,cList.back().X,MLoc,0,0);

      }
#endif

      bool trans = respFactory_.genSettings.isDist() or isRoot or
        not respFactory_.genSettings.formFullMat;

      // Form the linear transformation
      if( trans ) respFactory_.formLinearTrans(cList,M);


      // If distributed, gather distributed results
#ifdef CQ_ENABLE_MPI
      if( respFactory_.genSettings.isDist() ) {

        auto * grid = respFactory_.fullMatGrid();

        // Gather results to root process modelBasis1_LT
        grid->gather(N,nModel,modelBasis1_LT_,N,cList.back().AX,MLoc,0,0);

        // Dealloc memory if need be
        //if( cList.back().X )  CQMemManager::get().free(cList.back().X );
        //if( cList.back().AX ) CQMemManager::get().free(cList.back().AX);

      }
#endif

      if( isRoot ) {

        // Form full M' = V**H M V
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N,T(1.),
          modelBasis1_,N,modelBasis1_LT_,N,T(0.), this->fullMatrix_,nModel);

        // M' = L L**H
        int INFO = lapack::potrf(lapack::Uplo::Lower,nModel,this->fullMatrix_,nModel);
        if( INFO != 0 ) {

          CErr("Cholesky failed: H is not positive semidefinite");

        }

        // V -> VL**-H
        blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Lower,blas::Op::ConjTrans,blas::Diag::NonUnit,
          N,nModel,dcomplex(1.),this->fullMatrix_,nModel,modelBasis1_,N);

        // MV -> MVL**-H
        blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Lower,blas::Op::ConjTrans,blas::Diag::NonUnit,
          N,nModel,dcomplex(1.),this->fullMatrix_,nModel,modelBasis1_LT_,N);
      }


      dcomplex *KV = nullptr;
      // Alloc space to store KMV and set up contraction
      if( isRoot ) {
        KV = CQMemManager::get().malloc<dcomplex>(nModel*N);
        std::fill_n(KV,nModel*N,dcomplex(0.));
      }

      // If distributed...
      if( respFactory_.genSettings.isDist() ) {

#ifdef CQ_ENABLE_MPI
        auto * grid = respFactory_.fullMatGrid();

        // Scatter modelBasis1_LT_ to the BLACS grid
        grid->scatter(N,nModel,modelBasis1_LT_,N,cList.back().X,MLoc,0,0);
#endif


      } else {

        cList.back().X  = modelBasis1_LT_;
        cList.back().AX = KV;

      }

      if( trans ) respFactory_.formLinearTrans(cList,K);


      // If distributed, gather distributed results
#ifdef CQ_ENABLE_MPI
      if( respFactory_.genSettings.isDist() ) {

        auto * grid = respFactory_.fullMatGrid();

        // Gather results to root process modelBasis1_LT
        grid->gather(N,nModel,KV,N,cList.back().AX,MLoc,0,0);

        // Dealloc memory if need be
        if( cList.back().X )  CQMemManager::get().free(cList.back().X );
        if( cList.back().AX ) CQMemManager::get().free(cList.back().AX);

      }
#endif

      // Form V**H M**H K M V in FM
      if( isRoot )
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nModel,nModel,N,T(1.),modelBasis1_LT_,N,KV,N,T(0.),
          this->fullMatrix_,nModel);



      this->genSettings.jobType = RESIDUE;
      this->genSettings.evalProp = false;
      this->resSettings.nRoots = nModel;
      this->resSettings.needVL = true;
      this->morSettings.nModel = nModel;

      ResponseTBase<dcomplex>::run();


      if( isRoot ) {

        for(auto k = 0; k < nModel; k++) {
          this->resResults.W[k] = std::sqrt(this->resResults.W[k]);
          blas::scal(nModel,dcomplex(std::sqrt(this->resResults.W[k]/2.)),
            this->resResults.VR + k*nModel,1);
          blas::scal(nModel,
            dcomplex(std::sqrt(0.5)/std::sqrt(this->resResults.W[k])),
            this->resResults.VL + k*nModel,1);
        }


        ritzVecR_ = CQMemManager::get().malloc<dcomplex>(N*nModel);
        std::fill_n(ritzVecR_,N*nModel,dcomplex(0.));
        ritzVecL_ = CQMemManager::get().malloc<dcomplex>(N*nModel);
        std::fill_n(ritzVecL_,N*nModel,dcomplex(0.));

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nModel,nModel,dcomplex(1.),modelBasis1_,N,
           this->resResults.VR,nModel,dcomplex(0.),ritzVecR_,N);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nModel,nModel,dcomplex(1.),modelBasis1_LT_,N,
           this->resResults.VL,nModel,dcomplex(0.),ritzVecL_,N);



        resNorms_ = CQMemManager::get().malloc<double>(nModel);
        std::fill_n(resNorms_,nModel,0.);

        dcomplex * RES = CQMemManager::get().malloc<dcomplex>(nModel*N);
        std::fill_n(RES,nModel*N,dcomplex(0.));

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nModel,nModel,dcomplex(1.),KV,N,
           this->resResults.VR,nModel,dcomplex(0.),RES,N);

        for(auto k = 0; k < nModel; k++) {

          MatAdd('N','N',N,1,dcomplex(1.),RES + k*N,N,
            -dcomplex(this->resResults.W[k]*this->resResults.W[k]), 
            ritzVecR_ + k*N, N, RES + k*N, N);

          resNorms_[k] = blas::nrm2(N,RES + k*N,1);

        }

        CQMemManager::get().free(RES,KV);




        this->genSettings.evalProp = true;
        residueProperties();

        if( this->savFile.exists() ) {

          std::cout << "    * WRITING EIGENSYSTEM DATA TO DISK" << std::endl;

          this->savFile.safeWriteData("/MOR/EIG/MODELBASIS_RIGHT",modelBasis1_,
              {nModel, N});
          this->savFile.safeWriteData("/MOR/EIG/RITZVAL",this->resResults.W,
              {nModel});
          this->savFile.safeWriteData("/MOR/EIG/RITZCOEFF",this->resResults.VR,
              {nModel, nModel});
          this->savFile.safeWriteData("/MOR/EIG/RITZVEC_RIGHT",ritzVecR_,
              {nModel, N});
          this->savFile.safeWriteData("/MOR/EIG/RITZVEC_LEFT",ritzVecL_,
              {nModel, N});

          this->savFile.safeWriteData("/MOR/EIG/RESNORMS",resNorms_,
              {nModel});
        }

        respFactory_.blockTransform(N,nModel,std::sqrt(0.5),ritzVecR_,N,
            ritzVecL_,N);

      }

      //CErr();




    }

    void extractEigFull() {





      size_t N = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
      size_t nModel = this->morSettings.nModel;
      size_t nRoots = nModel;

      PolarizationPropagator<SingleSlater<RefT,IntsT>> &ss = 
        dynamic_cast<PolarizationPropagator<SingleSlater<RefT,IntsT>>&>(respFactory_);

      this->reset();
      if(this->fullMatrix_) CQMemManager::get().free(this->fullMatrix_);
      this->fullMatrix_ = nullptr;

      dcomplex one(1.);
      dcomplex *VEXP = 
        CQMemManager::get().malloc<dcomplex>( 2*nModel*N);
      std::fill_n(VEXP,2*nModel*N,dcomplex(0.));
      dcomplex *AVEXP = 
        CQMemManager::get().malloc<dcomplex>( 2*nModel*N);
      std::fill_n(AVEXP,2*nModel*N,dcomplex(0.));



      std::cout << "    * CONSTRUCTING PAIRED MODEL BASIS" << std::endl;
      SetMat('N',N,nModel,one,modelBasis1_,N,VEXP,N);
      SetMat('R',N/2,nModel,one,modelBasis1_,N,VEXP + nModel*N + N/2,N);
      SetMat('R',N/2,nModel,one,modelBasis1_+(N/2),N,VEXP + nModel*N,N);

      std::cout << "    * ORTHONORMALIZING PAIRED BASIS via SVD" << std::endl;
      double * SVAL = CQMemManager::get().malloc<double>(2*nModel);
      std::fill_n(SVAL,2*nModel,0.);
      dcomplex * dummy = nullptr;

      lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec,
        N,2*nModel,VEXP,N,SVAL,dummy,1,dummy,1);

      double orthoTol = 1e-12 * SVAL[0] * N;

      size_t nOrtho = 0;
      for(auto k = 0; k < 2*nModel; k++) if(SVAL[k] > orthoTol) nOrtho++;

      std::cout << "      * SVD RECEIVED " << 2*nModel << " VECTORS"
        << " AND RETURNED " << nOrtho << " VECTORS" << std::endl;

      // H-Orthogonalize the V basis and store in VEXP


      std::cout << "    * ORTHONORMALIZING PAIRED BASIS WRT H" << std::endl;
      ss.incMet = false;
      dcomplex *oldModel    = modelBasis1_;
      dcomplex *oldModel_LT = modelBasis1_LT_;
      modelBasis1_    = VEXP;
      modelBasis1_LT_ = AVEXP;
      this->morSettings.nModel = nOrtho;
      formFullMatrix();
      this->morSettings.nModel = nModel;
      modelBasis1_LT_ = oldModel_LT;
      modelBasis1_    = oldModel;
      ss.incMet = true;


      int INFO = lapack::potrf(lapack::Uplo::Lower,nOrtho,this->fullMatrix_,nOrtho);
      if( INFO != 0 ) {

        CErr("Cholesky failed: H is not positive semidefinite");

      }

      blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Lower,blas::Op::ConjTrans,blas::Diag::NonUnit,
        N,nOrtho,dcomplex(1.),this->fullMatrix_,nOrtho,VEXP,N);

      blas::trsm(blas::Layout::ColMajor,blas::Side::Right,blas::Uplo::Lower,blas::Op::ConjTrans,blas::Diag::NonUnit,
        N,nOrtho,dcomplex(1.),this->fullMatrix_,nOrtho,AVEXP,N);


      std::cout << "    * FORMING S INNER PRODUCT FOR EIGENSYSTEM" 
        << std::endl;
      // FM = V* S V
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nOrtho,nOrtho,N/2,dcomplex(1.),
        VEXP,N,VEXP,N,dcomplex(0.),this->fullMatrix_,nOrtho);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nOrtho,nOrtho,N/2,dcomplex(-1.),
        VEXP+N/2,N,VEXP+N/2,N,dcomplex(1.),this->fullMatrix_,nOrtho);


      std::cout << "    * PERFORMING REDUCED SPACE DIAGONALIZATION" 
        << std::endl;
      // Run the RESIDUE calculation
      this->genSettings.jobType = RESIDUE;
      this->genSettings.evalProp = false;
      this->resSettings.nRoots = nOrtho;
      this->morSettings.nModel = nOrtho;

      ResponseTBase<dcomplex>::run();

      std::cout << "    * EXTRACTING RITZ PAIRS" << std::endl;
      // Extract Ritz Values and reorder eigenvectors
      for(auto k = 0; k < nOrtho/2; k++){ 
        this->resResults.W[k] = 1./this->resResults.W[nOrtho-k-1];
        SetMat('N',nOrtho,1,dcomplex(1.),
          this->resResults.VR + (nOrtho-k-1)*nOrtho, nOrtho,
          this->resResults.VR + k*nOrtho,          nOrtho);
      }


      // Orthonormalize Ritz coeffs wrt S
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nOrtho,nOrtho,N/2,dcomplex(1.),
        VEXP,N,VEXP,N,dcomplex(0.),this->fullMatrix_,nOrtho);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nOrtho,nOrtho,N/2,dcomplex(-1.),
        VEXP+N/2,N,VEXP+N/2,N,dcomplex(1.),this->fullMatrix_,nOrtho);

      dcomplex * SCR = 
        CQMemManager::get().malloc<dcomplex>(nOrtho*nOrtho);
      std::fill_n(SCR,nOrtho*nOrtho,dcomplex(0.));
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nOrtho, nOrtho/2, nOrtho, dcomplex(1.),
          this->fullMatrix_, nOrtho, this->resResults.VR, nOrtho,
          dcomplex(0.), SCR,nOrtho);
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans, nOrtho/2, nOrtho/2, nOrtho, dcomplex(1.),
          this->resResults.VR, nOrtho, SCR, nOrtho,
          dcomplex(0.), this->fullMatrix_, nOrtho/2);

      lapack::gesvd(lapack::Job::OverwriteVec,lapack::Job::NoVec,
        nOrtho/2,nOrtho/2,this->fullMatrix_,nOrtho/2,SVAL,dummy,1,dummy,1);


      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans, nOrtho, nOrtho/2, nOrtho/2, dcomplex(1.),
          this->resResults.VR, nOrtho, this->fullMatrix_, nOrtho/2,
          dcomplex(0.), SCR, nOrtho);

      std::copy_n(SCR,nOrtho*nOrtho/2,this->resResults.VR);

      for(auto k = 0; k < nOrtho / 2; k++)
        blas::scal(nOrtho, dcomplex( 1. / std::sqrt(std::abs(SVAL[k])) ), 
            this->resResults.VR + k*nOrtho, 1) ;



      // Compute Ritz Vectors
      ritzVecR_ = CQMemManager::get().malloc<dcomplex>(N*nOrtho/2);
      std::fill_n(ritzVecR_,N*nOrtho/2,dcomplex(0.));
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nOrtho/2,nOrtho,dcomplex(1.),VEXP,N,
         this->resResults.VR,nOrtho,dcomplex(0.),ritzVecR_,N);


      // Compute Residuals
      dcomplex *RES = 
        CQMemManager::get().malloc<dcomplex>(N*nOrtho/2);
      std::copy_n(ritzVecR_,N*nOrtho/2,RES);

      for(auto k = 0; k < nOrtho/2; k++)
        blas::scal(N,-dcomplex(this->resResults.W[k]), RES + k*N, 1);

      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N/2,nOrtho/2,nOrtho,dcomplex(1.),AVEXP,N,
         this->resResults.VR,nOrtho,dcomplex(1.),RES,N);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N/2,nOrtho/2,nOrtho,dcomplex(-1.),AVEXP+N/2,N,
         this->resResults.VR,nOrtho,dcomplex(1.),RES + N/2,N);

      resNorms_ = CQMemManager::get().malloc<double>(nOrtho/2);
      for(auto k = 0; k < nOrtho/2; k++)
        resNorms_[k] = blas::nrm2(N,RES+k*N,1);

      // Compute properties using Ritz vectors
      this->resSettings.nRoots = nOrtho/2;
      this->morSettings.nModel = nOrtho;
      this->genSettings.evalProp = true;
      CQMemManager::get().free(modelBasis1_, modelBasis1_LT_); // Freeup mem
      modelBasis1_    = VEXP;
      modelBasis1_LT_ = AVEXP;
      residueProperties();

      if( this->savFile.exists() ) {

        std::cout << "    * WRITING EIGENSYSTEM DATA TO DISK" << std::endl;

        this->savFile.safeWriteData("/MOR/EIG/MODELBASIS_RIGHT",modelBasis1_,
            {nOrtho, N});
        this->savFile.safeWriteData("/MOR/EIG/RITZVAL",this->resResults.W,
            {nOrtho / 2});
        this->savFile.safeWriteData("/MOR/EIG/RITZCOEFF",this->resResults.VR,
            {nOrtho / 2, nOrtho});
        this->savFile.safeWriteData("/MOR/EIG/RITZVEC_RIGHT",ritzVecR_,
            {nOrtho / 2, N});

        this->savFile.safeWriteData("/MOR/EIG/RESNORMS",resNorms_,
            {nOrtho/2});
      }


    }



  };

  template <class Reference>
  void MORSpec<Reference>::configOptions() {
    

  };

  template <class Reference>
  void MORSpec<Reference>::formRHS() {


    bool isRoot = MPIRank(this->comm_) == 0;
    bool rootHasRHS = isRoot and respFactory_.dfdrResults.RHS;
    MPIBCast(rootHasRHS,0,this->comm_);

    // Form full RHS if not already formed
    if( not rootHasRHS ) respFactory_.formRHS();

    ROOT_ONLY(this->comm_);

    // If subspace RHS already allocated, remove it
    if( this->dfdrResults.RHS ) CQMemManager::get().free(this->dfdrResults.RHS);

    size_t N      = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
    size_t NRHS   = respFactory_.fdrSettings.nRHS;
    size_t nModel = this->morSettings.nModel;

    // Allocate subspace RHS
    T* RHS = CQMemManager::get().malloc<dcomplex>(NRHS*nModel);
    std::fill_n(RHS,NRHS*nModel,dcomplex(0.));



    // Calculate RHS (subspace) = V^* RHS (full)
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NRHS,nModel,N,dcomplex(1.),respFactory_.dfdrResults.RHS,N,
      modelBasis1_,N,dcomplex(0.),RHS,NRHS);
    IMatCopy('C',NRHS,nModel,dcomplex(1.),RHS,NRHS,nModel);

    this->dfdrResults.RHS = RHS;

  };



  template <class Reference>
  typename MORSpec<Reference>::T* MORSpec<Reference>::formFullMatrix() {

    size_t N      = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
    size_t nModel = this->morSettings.nModel;
    bool isRoot = MPIRank(this->comm_) == 0;

    std::vector<RESPONSE_CONTRACTION<T>> cList(1);
    cList.back().nVec   = nModel;
    cList.back().N      = N;
    cList.back().X      = modelBasis1_;
    cList.back().AX     = modelBasis1_LT_;

    int64_t MLoc, NLoc;

#ifdef CQ_ENABLE_MPI
    if( respFactory_.genSettings.isDist() ) {

      auto * grid = respFactory_.fullMatGrid();
      std::tie(MLoc,NLoc) = grid->get_local_dims(N,nModel);

      if( MLoc and NLoc ) {

        cList.back().X  = CQMemManager::get().malloc<T>(MLoc*NLoc);
        std::fill_n(cList.back().X,MLoc*NLoc,T(0.));
        cList.back().AX = CQMemManager::get().malloc<T>(MLoc*NLoc);
        std::fill_n(cList.back().AX,MLoc*NLoc,T(0.));

      } else {

        cList.back().X  = nullptr;
        cList.back().AX = nullptr;

      }

      cList.back().DescX  = grid->descinit_noerror(N,nModel,MLoc);
      cList.back().DescAX = cList.back().DescX;

      grid->scatter(N,nModel,modelBasis1_,N,cList.back().X,MLoc,0,0);

    }
#endif

    bool trans = respFactory_.genSettings.isDist() or isRoot or
      not respFactory_.genSettings.formFullMat;

    if( trans ) respFactory_.formLinearTrans(cList);

#ifdef CQ_ENABLE_MPI
    if( respFactory_.genSettings.isDist() ) {

      auto * grid = respFactory_.fullMatGrid();
      grid->gather(N,nModel,modelBasis1_LT_,N,cList.back().AX,MLoc,0,0);

      if( cList.back().X )  CQMemManager::get().free(cList.back().X );
      if( cList.back().AX ) CQMemManager::get().free(cList.back().AX);

    }
#endif



    this->fullMatrix_ = isRoot ? 
      CQMemManager::get().malloc<T>(N*N) : nullptr;

    if( isRoot ) {
      std::fill_n(this->fullMatrix_,N*N,T(0.));
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, nModel, nModel, N, T(1.), modelBasis1_,
                 N, modelBasis1_LT_, N, T(0.),
                 this->fullMatrix_, nModel);
    }

    return this->fullMatrix_;

  };

  template <class Reference>
  std::pair<size_t,typename MORSpec<Reference>::T*> 
    MORSpec<Reference>::formPropGrad(ResponseOperator op) {

    bool isRoot = MPIRank(this->comm_) == 0;
    size_t N      = respFactory_.getNSingleDim(respFactory_.genSettings.doTDA);
    size_t nModel = this->morSettings.nModel;

    typename Reference::value_type* g = nullptr;
    size_t nProp = 0;

    std::tie(nProp,g) = respFactory_.formPropGrad(op);


    T* newG = isRoot ? 
      CQMemManager::get().malloc<T>(nProp*nModel) : nullptr;

    if( isRoot ) {
      std::fill_n(newG,nProp*nModel,T(0.));
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nProp,nModel,N,dcomplex(1.),g,N,modelBasis1_,N,
        dcomplex(0.),newG,nProp);
      IMatCopy('C',nProp,nModel,dcomplex(1.),newG,nProp,nModel);
    }


    if( g ) CQMemManager::get().free(g);

    return {nProp, newG};

  };


};

