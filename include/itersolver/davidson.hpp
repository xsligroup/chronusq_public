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
#ifndef __INCLUDED_DAVIDSON_HPP__
#define __INCLUDED_DAVIDSON_HPP__

#include <itersolver.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>
#include <cqlinalg/factorization.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/ortho.hpp>
#include <cqlinalg/eig.hpp>
#include <cerr.hpp>

// #define DEBUG_DAVIDSON


namespace ChronusQ {


  template <typename _F>
  bool Davidson<_F>::runMicro() {
    
    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;
    
    if( isRoot ) {
      std::cout  << "\n\n";
      std::cout  << "  * Davidson Settings:\n";
      std::cout  << "    * Right Eigenvector is Requested. \n";
      
      if(this->DoLeftEigVec) {
        std::cout  << "    * Left Eigenvector is Requested.\n";
        CErr("Do Left Eig Vec is not implemented yet");
      }
      
      if (this->EnergySpecific) {          
        std::cout<< "    * Use Energy Specific:           " << (this->EnergySpecific ? "True" : "False") << "\n"
                 << "      * Number of Low  Energy Roots = " << this->nHighERoots    << "\n"
                 << "      * Number of High Energy Roots = " << this-> nLowERoots    << "\n";
        if(this->adaptiveERef) {
          std::cout<< "      * Use Ground State Energy in each iteration as Reference \n";
        } else {
          std::cout<< "      * Energy Referene  = " << this->EnergyRef    << "\n";
        }
      }          
      std::cout << "\n\n" << std::endl;
    }
    
    std::cout << std::setprecision(10) << std::scientific;
    
    const size_t N   = this->N_;
    const size_t MSS = this->mSS_;
    const size_t nR  = this->nRoots_;
    
    const size_t MSS2 = MSS * MSS;
    const size_t NMSS = N * MSS;
    const size_t NNR  = N * nR;
    
    const size_t nG  = this->nGuess_;
    const size_t NNG = N * nG;  
    
    // Initialize pointers
    std::shared_ptr<SolverVectors<_F>> VR  = nullptr;
    _F *XR  = nullptr, *XRPrev = nullptr;
    std::shared_ptr<SolverVectors<_F>> VL  = nullptr;
    _F *XL  = nullptr; // for left search space
    std::shared_ptr<SolverVectors<_F>> AVR = nullptr;
    _F *Ovlp = nullptr; // overlap is using right eigenvectors
    std::shared_ptr<SolverVectors<_F>> R   = nullptr, S = nullptr; // Scratch space for residue and perturbed vector
    _F *SubA = nullptr;
    _F *SCR = nullptr;

    dcomplex *Eig = nullptr, *EPrev = nullptr;

    VR     = this->vecGen_(MSS);
    AVR    = this->vecGen_(MSS);
    R      = this->vecGen_(nG);
    S      = this->vecGen_(nG);
    SubA   = this->memManager_.template malloc<_F>(MSS2);
    if(this->DoLeftEigVec) {
      VL   = this->vecGen_(MSS);
    }

    XR     = this->memManager_.template malloc<_F>(MSS2);
    XRPrev = this->memManager_.template malloc<_F>(MSS2);
    Eig    = this->memManager_.template malloc<dcomplex>(MSS);
    EPrev  = this->memManager_.template malloc<dcomplex>(MSS);
    Ovlp   = this->memManager_.template malloc<_F>(MSS2);
    SCR    = this->memManager_.template malloc<_F>(MSS2);

    if(this->DoLeftEigVec) {
      XL   = this->memManager_.template malloc<_F>(MSS2);
    }

    // variables during iteration
    std::vector<bool> SiConv(nG); // state_i converged?
    size_t iter   = 0;
    size_t nDo    = nG;
    size_t nExam  = nG;
    size_t nVPrev = 0;     // number of vectors at previous iteration
    size_t nVCur  = nG;    // number of vectors at current iteration
    char   JOBVL  = this->DoLeftEigVec ? 'V': 'N';
    double VecNear = 0.1;  // criterion to determine if two vector are similar by their overlap
    double VecConv = this->convCrit_;
    double EConv   = VecConv * 0.01;

    // generate guess
    // Initailize VR as Guess, if no Guess set,
    if( Guess ) {
      VR->set_data(0, nG, *Guess, 0);

#ifndef DEBUG_DAVIDSON
      std::cout.setstate(std::ios_base::failbit);
#endif
      nVCur = VR->GramSchmidt(0,0,nG,this->memManager_,GramSchmidt_NRe,GramSchmidt_eps);
#ifndef DEBUG_DAVIDSON
      std::cout.clear();
#endif
    } else {
      std::cout << "  * use unit vector guess" << std::endl;
      VR->clear();
      for(auto i = 0ul; i < nG; i++) VR->set(i, i, 1.0);
    } // right vector guess

    // left vector guess
    if(this->DoLeftEigVec) {
      // VL and VR should be biothogonalized
      CErr("Do Left Eig Vec is not implemented yet");
    }

    if( isRoot ) {
      std::cout << "\n\n  * Starting Davidson Iterations" << std::endl;
    } // Root Only

    // ****************************
    // ** Begin Davidson iterations **
    // ****************************
    for( iter = 0; iter < this->maxMicroIter_; iter++) {
      
      auto DavidsonSt = tick();
      
      if( isRoot ) {
        std::cout << "\n    DavidsonIter " << std::setw(5) << iter+1  
                  << ": Number of new vectors = " << std::setw(5) << nDo << std::endl;
      } // Root Only

      auto LTst = tick();
      
      // AVR <- A * VR
//      if (isRoot) {
        SolverVectorsView<_F> VRSend(*VR, nVPrev), AVRRecv(*AVR, nVPrev);
        this->linearTrans_(nDo, VRSend, AVRRecv);
//      }
      
      double LTdur = tock(LTst);

#ifdef DEBUG_DAVIDSON
      VR->print(std::cout, "VR");
      AVR->print(std::cout, "AVR");
#endif
        
      // Construct submatrix of A
      if(this->DoLeftEigVec) {
      //SubA <- VL * AVR
        CErr("Do Left Eig Vec is not implemented yet");
      } else {
      // SubA <- VR_\dagger * AVR
        VR->dot_product(0, *AVR, 0, nVCur,nVCur,SubA,nVCur);
      }

//      if( isRoot ) {
    
#ifdef DEBUG_DAVIDSON
        prettyPrintSmart(std::cout,"HH Davidson SubMatix ",SubA,nVCur,nVCur,nVCur);
#endif 

        // Diagonalize SubA 
        GeneralEigen(JOBVL, 'V', nVCur, SubA, nVCur, Eig, XL, nVCur, XR, nVCur);

        // swap high energy roots for energy specific
        if(this->EnergySpecific) {

          if(this->adaptiveERef)
            this->EnergyRef = std::real(Eig[0]);

          std::vector<size_t> indx(nVCur,0);
          std::iota(indx.begin(), indx.end(), 0);

          Eigen::Map<
          Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
          > XRMap(XR,nVCur,nVCur);


          if ( this->sortByDistance ){
            std::stable_sort(indx.begin(), indx.end(),
                             [&] (size_t i , size_t j) {
              return std::abs(Eig[i] - this->EnergyRef) < std::abs(Eig[j] - this->EnergyRef);
            }
            );
          }
          else {
            std::stable_sort(indx.begin(), indx.end(),
                             [&] (size_t i , size_t j) {
              if ((std::real(Eig[i]) > this->EnergyRef  and std::real(Eig[j]) > this->EnergyRef)
              or (std::real(Eig[i]) < this->EnergyRef  and std::real(Eig[j]) < this->EnergyRef))
                return std::real(Eig[i]) < std::real(Eig[j]);
              else if (std::real(Eig[i]) < this->EnergyRef  and std::real(Eig[j]) > this->EnergyRef )
                return false;
              else
                return true;
            }
            );
          }



          for(auto i = 0ul; i < nVCur - 1; i++){
            size_t ind = indx[i];
            while(ind < i) ind = indx[ind];

            if (ind > i) {
              std::swap(Eig[i],Eig[ind]);
              XRMap.col(i).swap(XRMap.col(ind));
            }
          }

//          REMOVING LINEAR DEPENDENCY?
//          for (auto i = 0ul; i < nVCur; i++){
//            if (std::abs(Eig[i]) < 1e-8) {
//              nVCur = i;
//              std::cout << " Reducing nVCur: " << std::endl;
//              std::cout << nVCur << "  " << Eig[i] << std::endl;
//              break;
//            }
//          }
//          if (nVCur == 0) CErr("Not enough eigenpairs to proceed! ");
        }   

        // print eigenvalues at current iteration
        std::cout << "      - Eigenvalues at the current iteration:" << std::endl;
        for(auto i = 0; i < nExam; i++) {
          std::cout << "        Root " << std::setw(5) << std::right << i << ":"
                    << std::right << std::setw(20) << std::real(Eig[i]);
              
          if( std::is_same<dcomplex,_F>::value ) {
            std::cout << " + " << std::setw(20) << std::imag(Eig[i]) << " i";
          }
              
          std::cout << std::endl;
        }
          
        // Exam Eigenvalues and eigenvectors and do mapping if iter > 0
        std::fill_n(SiConv.begin(),nExam,false);
        if( iter > 0) {
            
          // overlap = (VR XR)_old ^\dagger * (VR XR)_new
          std::fill_n(Ovlp, nVCur*nVPrev, _F(0.)); 
              for(auto i = 0ul; i < nVCur; i++) Ovlp[i + i*nVPrev] = _F(1.); 
              
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nExam,nVCur,nVPrev,_F(1.),XRPrev,nVPrev,Ovlp,nVPrev,_F(0.),SCR,nExam);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nExam,nExam,nVCur ,_F(1.),SCR,nExam,XR,nVCur,_F(0.),Ovlp,nExam);
          
          // mapping old vector to new vectors based on overlap
          std::vector<int> StMap(nExam);
          std::fill_n(StMap.begin(),nExam,-1);
          std::copy_n(Ovlp,nExam*nExam,SCR);
              
          // i -> new state, j -> old state
          int i = 0, j = 0;
          for(i = 0; i < nExam; i++) {
            auto iO = Ovlp + i*nExam;
            j = std::distance(iO, std::max_element(iO, iO+nExam, 
                  [&] (_F A, _F B) {return std::norm(A) < std::norm(B); 
                                  }));
            if( std::abs(iO[j]) > VecNear) {
              StMap[i] = j;
              for (auto k = i+1; k < nExam; k++) Ovlp[j + k*nExam] = _F(0.);
            }
          }
        
          std::copy_n(SCR,nExam*nExam,Ovlp);
          
#ifdef DEBUG_DAVIDSON
          prettyPrintSmart(std::cout,"HH Davidson States Overlap",Ovlp,nExam,nExam,nExam);
#endif 
      
          std::cout << "\n      - Comparison to the previous iteration: " << std::endl;  
          // exam eigenvectors
          // R and S as scratch space to hold full vector old and new repectively
          VR->multiply_matrix(0, blas::Op::NoTrans,nExam,nVPrev,_F(1.),XRPrev,nVPrev,_F(0.),*R, 0);
          VR->multiply_matrix(0, blas::Op::NoTrans,nExam,nVCur,_F(1.),XR,nVCur,_F(0.),*S, 0);
              
          _F phase;
          double maxDel;  // maximum differece in vectors
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            if(j < 0) {
              std::cout << "          New root" << std::setw(5) << i << " is brand new" << std::endl;
              continue; 
            } else if (j!=i) { 
              std::cout << "          New root" << std::setw(5) << i
                        << " was old root" << std::setw(5) << j << std::endl; 
            }

            phase = Ovlp[j + i*nExam];
            phase /= std::abs(phase);

            std::shared_ptr<SolverVectors<_F>> iNew = std::make_shared<SolverVectorsView<_F>>(*S, i);
            std::shared_ptr<SolverVectors<_F>> jOld = std::make_shared<SolverVectorsView<_F>>(*R, j);

            jOld->scale(-phase, 0, 1);
            jOld->axpy(0, 1, _F(1.), *iNew, 0);
            maxDel = jOld->maxNormElement(0, 1);
 
            SiConv[i] = maxDel < VecConv;
            if(SiConv[i]) std::cout << "        Root "   << std::setw(5) << std::right << i
                                    << " has converged"  << std::endl;
            else std::cout << "        Root " <<  std::right << std::setw(5) << i
                           << " has not converged, maximum delta is"<<  std::right  
                           << std::setw(20) << maxDel << std::endl;
          }

          // exam eigenvalues
          dcomplex EDiff;
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            EDiff = Eig[i] - EPrev[j];
            SiConv[i] = SiConv[i] and (std::abs(EDiff) < EConv);  
          } 
        }   // Exam eigenvales and eigenvectors
          
        isConverged = std::all_of(SiConv.begin(), SiConv.begin()+nExam, [&] (bool i) { return i; });
        if(isConverged or nVCur >= MSS) {
        
          double DavidsonDur = tock(DavidsonSt);
          double perLT = LTdur * 100 / DavidsonDur;
        
          std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur 
            << " s  ( " << perLT << " % LT )" << std::endl;
          break;
        }

        // Form residue vectors only for unconverged vectors
        std::vector<int> unConvS(nExam);
        nDo = 0ul;
        std::fill_n(unConvS.begin(),nExam,-1);
        std::fill_n(SCR,nVCur*nExam,_F(0.));
        for (auto i = 0ul; i < nExam; i++) {
          if(not SiConv[i]) {
            unConvS[nDo] = i;
            std::copy_n(XR + i*nVCur,nVCur,SCR + nDo*nVCur);
            nDo++;
          }
        }

        if(nVCur+nDo > MSS)  nDo = MSS - nVCur;
        
        // R <- AVR * XR
        AVR->multiply_matrix(0, blas::Op::NoTrans,nDo,nVCur,_F(1.),SCR,nVCur,_F(0.),*R, 0);

        // S as scratch space, <- eig_i * (VR * XR_i)
        VR->multiply_matrix(0, blas::Op::NoTrans,nDo,nVCur,_F(1.),SCR,nVCur,_F(0.),*S, 0);
            
        for (auto i = 0ul; i < nDo; i++)
          S->scale(this->dcomplexTo_F(Eig[unConvS[i]]), i, 1);
            
        // Compute the residue norm and generate perturbbed vectors
        S->scale(_F(-1.), 0, nDo);
        S->axpy(0, nDo, _F(1.), *R, 0);
        std::cout << "\n      - Residues of non-converged roots: " << std::endl;  
        
        if(EigForT) std::fill_n(EigForT,nG,dcomplex(0.));
        std::fill_n(RelRes, nG, 0.);
        for (auto i = 0ul; i < nDo; i++) {
          auto j = unConvS[i];
          RelRes[i] = S->norm2F(i, 1);
          std::cout << "        Root " << std::setw(5) << std::right << j+1 
                    << " 2nd order lowering " << std::right << std::setw(20) << RelRes[i]*RelRes[i] 
                    << " norm " << std::right << std::setw(20) << RelRes[i] << std::endl;
              
          if(EigForT) this->EigForT[i] = Eig[j];
        }

#ifdef DEBUG_DAVIDSON
        S->print(std::cout, "S before preCond");
#endif
        this->preCondNoShift_(nDo,*S,*S);

#ifdef DEBUG_DAVIDSON
        S->print(std::cout, "S after preCond");
#endif
        // Append the new vectors to VR and orthogoalize against existing ones 
        // Also update the dimensions and save the XRPrev
        VR->set_data(nVCur, nDo, *S, 0);
        std::copy_n(XR,nVCur*nVCur,XRPrev);
        std::copy_n(Eig,nVCur,EPrev);
        nVPrev = nVCur;
        if(this->DoLeftEigVec) {
          CErr("Do Left Eig Vec is not implemented yet");
        } else {
          // disable printing from GramSchmidt
          
#ifndef DEBUG_DAVIDSON
          std::cout.setstate(std::ios_base::failbit);
#endif

#ifdef DEBUG_DAVIDSON
          VR->print(std::cout, "VR before GramSchmidt");
#endif
          nVCur = VR->GramSchmidt(0, nVCur,nDo,this->memManager_,GramSchmidt_NRe,GramSchmidt_eps);
#ifdef DEBUG_DAVIDSON
          VR->print(std::cout, "VR after GramSchmidt");
#endif
          
#ifndef DEBUG_DAVIDSON
          std::cout.clear();
#endif
          
          nDo   = nVCur - nVPrev;
        }
         
        if(iter+1 == this->whenSc) { 
          nExam = nR; 
          nDo   = std::min(nDo, nExam);
          nVCur = nVPrev + nDo;
        }
        
        double DavidsonDur = tock(DavidsonSt);
        double perLT = LTdur * 100 / DavidsonDur;
    
        std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur 
          << " s  ( " << perLT << " % LT )" << std::endl;
    
        if(nDo == 0) {
          isConverged = true;
          break;
        }
//      } // Root Only
    
    } // Davidson iteration    

//    if( isRoot ) {

      // move data before exit runMicro      
      std::copy_n(Eig,nR,this->eigVal_);
      VR->multiply_matrix(0, blas::Op::NoTrans,nR,nVCur,_F(1.),XR,nVCur,_F(0.),*this->VR_, 0);
      //size_t nVSave = isConverged ? nR: nG;  
      //std::copy_n(Eig,nVSave,this->eigVal_);
      //blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVSave,nVCur,_F(1.),VR,N,XR,nVCur,_F(0.),this->VR_,N);

    if( isRoot ) {
      std::cout << "\n  * ";
      if( isConverged and nDo !=0)
        std::cout << "Davidson Converged in " << iter+1 << " Iterations" 
                  << std::endl;
      else if (nDo ==0) {
        double maxResDel = *std::max_element(RelRes, RelRes+nR);
        std::cout << "Davidson expansion finished, and wavefunction converged "
          << "below " << std::setw(20) << std::right << maxResDel << std::endl;
      } else
        std::cout << "Davidson Failed to Converged in " << iter+1 << " Iterations" 
                  << std::endl;
    } // Root Only 

    // Free Scratch space

    if(XR)       this->memManager_.template free(XR);    
    if(XRPrev)   this->memManager_.template free(XRPrev);
    if(SubA)     this->memManager_.template free(SubA);
    if(Eig)      this->memManager_.template free(Eig);
    if(EPrev)    this->memManager_.template free(EPrev); 
    if(Ovlp)     this->memManager_.template free(Ovlp);
    if(SCR)      this->memManager_.template free(SCR);
    if(XL)       this->memManager_.template free(XL);    
    
    return isConverged;
  
  } // Davidson::runMicro

  template <typename _F>
  void Davidson<_F>::restart() {
    // copy full vectors as new guess  
    std::cout << "\n  * Restarting Davidson..." << std::endl;
    
    // restart with only nRoots_ of guess
    // as now it's more close to the solution
    this->kG      = 1;
    this->whenSc  = 1;
    this->nGuess_ = this->nRoots_;

    Guess = this->vecGen_(this->nGuess_);

    // NO MPI
    // ROOT_ONLY(this->comm_);
    this->Guess->set_data(0, this->nGuess_, *this->VR_, 0);
  } // Davidson::restart
  
}; // namespace ChronusQ

#endif
