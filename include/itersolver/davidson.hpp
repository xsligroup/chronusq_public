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
#include <list>

// #define DEBUG_DAVIDSON
// #define DAVIDSON_PRINT_TIMING


namespace ChronusQ {


  template <typename _F>
  bool Davidson<_F>::runMicro() {

    std::shared_ptr<SolverVectors<_F>> &VR = vecs;
    std::shared_ptr<SolverVectors<_F>> &AVR = sigmaVecs;
    
    bool isRoot = MPIRank(this->comm_) == 0;
    bool isConverged = false;
    
    if( isRoot ) {
      std::cout << std::setprecision(3) << std::scientific;
      std::cout << std::endl << std::endl;
      std::cout << "  * Davidson Settings:" << std::endl;
      std::cout << "    * Right Eigenvector is Requested.";
      auto printCheckCrit = [](bool check, double crit) {
        if (check)
          std::cout << crit;
        else
          std::cout << "N/A";
      };
      std::cout << std::endl << "    * Residual convergence    : ";
      printCheckCrit(checkResidueConv, this->convCrit_);
      std::cout << std::endl << "    * Eigenvalue convergence  : ";
      printCheckCrit(checkEigenValueConv, eigenValueCrit);
      std::cout << std::endl << "    * Eigenvector convergence : ";
      printCheckCrit(checkEigenVectorConv, eigenVectorCrit);
      std::cout << std::endl;
      
      if(this->DoLeftEigVec) {
        std::cout  << "    * Left Eigenvector is Requested." << std::endl;
        CErr("Do Left Eig Vec is not implemented yet");
      }
      
      if (!this->energyRefs.empty()) {
        std::cout << "    * Use Energy Specific:           True";
        if (this-> AbsoluteES) std::cout<< " with absolute energy threshold." << std::endl;
        else std::cout<< " with relative energy threshold." << std::endl;
        for (auto & pair: this->energyRefs) {
          std::cout<< "    * Energy Threshold = " << pair.first << " a.u., Number of Roots = "
                   << pair.second << std::endl;
        }
        std::cout<< "    * Number of Low  Energy Roots = " << this-> nLowERoots << std::endl;
      }          
      std::cout << std::endl << std::endl;
    }
    
    std::cout << std::setprecision(10) << std::scientific;

    const size_t MSS = this->mSS_;
    const size_t nR  = this->nRoots_;
    
    const size_t MSS2 = MSS * MSS;
    
    const size_t nG  = this->nGuess_;

    // variables during iteration
    size_t iter   = 0;
    size_t nDo    = nG < MSS ? nG : MSS;
    size_t nExam  = (whenSc > 1 ? nDo : nR);
    size_t nVPrev = 0;     // number of vectors at previous iteration
    size_t nVCur  = nDo;    // number of vectors at current iteration
    const char   JOBVL  = this->DoLeftEigVec ? 'V': 'N';
    const double VecNear = 0.1;  // criterion to determine if two vector are similar by their overlap
    const double VecConv = eigenVectorCrit;
    const double EConv   = eigenValueCrit;
    std::vector<IterDiagConvStatus> SiConv(nExam); // state_i converged?
    
    // Initialize pointers
    _F *XR  = nullptr, *XRPrev = nullptr;
    std::shared_ptr<SolverVectors<_F>> VL  = nullptr;
    _F *XL  = nullptr; // for left search space
    _F *Ovlp = nullptr; // overlap is using right eigenvectors
    _F *SubA = nullptr;
    _F *SCR = nullptr;

    dcomplex *Eig = nullptr, *EPrev = nullptr;

    if (not VR or VR->size() < MSS)
      VR   = this->vecGen_(MSS);
    if (not AVR or AVR->size() < MSS)
      AVR  = this->vecGen_(MSS);
    if (not R or R->size() < nExam)
      R    = this->vecGen_(nExam);
    if (not S or S->size() < nExam)
      S    = this->vecGen_(nExam);
    SubA   = CQMemManager::get().malloc<_F>(MSS2);
    if(this->DoLeftEigVec) {
      VL   = this->vecGen_(MSS);
    }

    XR     = CQMemManager::get().malloc<_F>(MSS2);
    XRPrev = CQMemManager::get().malloc<_F>(MSS2);
    Eig    = CQMemManager::get().malloc<dcomplex>(MSS);
    EPrev  = CQMemManager::get().malloc<dcomplex>(MSS);
    Ovlp   = CQMemManager::get().malloc<_F>(MSS2);
    SCR    = CQMemManager::get().malloc<_F>(MSS2);

    if(this->DoLeftEigVec) {
      XL   = CQMemManager::get().malloc<_F>(MSS2);
    }

    // generate guess
    // Initailize VR as Guess, if no Guess set,
    if( Guess ) {
      VR->set_data(0, nDo, *Guess, 0, true);

#ifndef DEBUG_DAVIDSON
      std::cout.setstate(std::ios_base::failbit);
#endif
      nVCur = VR->GramSchmidt(0,0,nDo,GramSchmidt_NRe,GramSchmidt_eps);
#ifndef DEBUG_DAVIDSON
      std::cout.clear();
#endif
    } else {
      std::cout << "  * use unit vector guess" << std::endl;
      VR->clear();
      for(auto i = 0ul; i < nDo; i++) VR->set(i, i, 1.0);
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
                  << ": Number of new vectors = " << std::setw(5) << nDo
                  << ", Subspace dimension = " << std::setw(5) << nVCur << std::endl;
      } // Root Only

      auto LTst = tick();
      
      // AVR <- A * VR
//      if (isRoot) {
        SolverVectorsView<_F> VRSend(*VR, nVPrev), AVRRecv(*AVR, nVPrev);
        this->linearTrans_(nDo, VRSend, AVRRecv);
#ifdef DEBUG_DAVIDSON
        std::cout << "HH Davidson Norm of AVR: " << AVR->norm2F(nVPrev, nDo) << std::endl;
#endif
//      }
      
      double LTdur = tock(LTst);

#ifdef DAVIDSON_PRINT_TIMING
      if( isRoot ) {
        std::cout << "      Linear transformation took "
        << std::fixed << std::setprecision(6) << LTdur << " s." << std::endl;
      } // Root Only
#endif

#ifdef DEBUG_DAVIDSON
      VR->print(std::cout, "VR");
      AVR->print(std::cout, "AVR");
#endif
        
      // Construct submatrix of A
      if(this->DoLeftEigVec) {
      //SubA <- VL * AVR
        CErr("Do Left Eig Vec is not implemented yet");
      } else {
#ifdef DAVIDSON_PRINT_TIMING
        auto DOTst = tick();
#endif
      // SubA <- VR_\dagger * AVR
//        VR->dot_product(0, *AVR, 0, nVCur,nVCur,SubA,MSS);
        VR->dot_product(nVPrev, *AVR, 0, nDo, nVCur, SubA + nVPrev, MSS);
        VR->dot_product(0, *AVR, nVPrev, nVPrev, nDo, SubA + nVPrev * MSS,MSS);

#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "      Dot production took "
          << std::fixed << std::setprecision(6) << tock(DOTst) << " s." << std::endl;
        } // Root Only
#endif
      }

#ifdef CQ_HAS_TA
      if (TA::initialized())
        TA::get_default_world().gop.fence();
#endif
    
#ifdef DEBUG_DAVIDSON
      prettyPrintSmart(std::cout,"HH Davidson SubMatix ",SubA,nVCur,nVCur,MSS);
#endif 

        // Diagonalize SubA
        if( isRoot ) {
          SetMat('N', nVCur, nVCur, 1.0, SubA, MSS, SCR, nVCur);
#ifdef DAVIDSON_PRINT_TIMING
          auto EIGst = tick();
#endif
          if( DoHerm ){
            HermetianEigen('V', 'L', nVCur, SCR, nVCur, Eig);
            std::copy_n(SCR,nVCur*nVCur,XR);
          }else
            GeneralEigen(JOBVL, 'V', nVCur, SCR, nVCur, Eig, XL, nVCur, XR, nVCur);
#ifdef DEBUG_DAVIDSON
          prettyPrintSmart(std::cout,"HH Davidson XR",XR,nVCur,nVCur,nVCur);
#endif

#ifdef DAVIDSON_PRINT_TIMING
          std::cout << "      Subspace eigen took "
          << std::fixed << std::setprecision(6) << tock(EIGst) << " s." << std::endl;
#endif

        // swap high energy roots for energy specific
        if(not this->energyRefs.empty()) {

#ifdef DAVIDSON_PRINT_TIMING
          auto SWAPst = tick();
#endif

          std::vector<size_t> sortedIndices;
          sortedIndices.reserve(kG * nR);

          double Eoffset = this->AbsoluteES? 0. : std::real(Eig[0]); // the lowest eigenvalue

          // Initialize a candidate list
          std::list<size_t> candList(nVCur);
          std::iota(candList.begin(), candList.end(), 0);
          std::list<size_t> selectedList;

          std::vector<size_t> missings(this->energyRefs.size()+1,0);
          double curERef = std::real(Eig[0]);
          size_t curRequest = this->nLowERoots;
          std::list<size_t>::iterator curIt = candList.begin();
          for (size_t i = 0; i <= this->energyRefs.size(); i++) {
            // Find the next energy reference
            double nextERef = 0.0;
            std::list<size_t>::iterator nextIt;
            if (i == this->energyRefs.size()) {
              nextERef = std::numeric_limits<double>::max();
              nextIt = candList.end();
            } else {
              nextERef = this->energyRefs[i].first + Eoffset;
              nextIt = std::lower_bound(candList.begin(), candList.end(), nextERef,
                                                                  [&Eig](size_t j, double x){ return std::real(Eig[j]) < x; });
            }

            size_t count = 0;
            // Moves indices from candList to selectedList until next energy reference begins
            while (curIt != nextIt and count < curRequest) {
              selectedList.push_back(*curIt);
              curIt = candList.erase(curIt);
              count++;
            }
            missings[i] = curRequest - count;

            if (i == this->energyRefs.size()) break; // Jump out of the loop if it is the last energy reference
            curERef = nextERef;
            curRequest = this->energyRefs[i].second;
            curIt = nextIt;
          }

          if (kG > 1) {
            // Add missing counts
            missings[0] += this->nLowERoots * (kG - 1);
            for (size_t i = 1; i <= this->energyRefs.size(); i++) {
              missings[i] += this->energyRefs[i-1].second * (kG - 1);
            }
          }

          // Fill in missing eigenvalues with nearby eigenvalues.
          curIt = candList.begin();
          std::advance(curIt, missings[0]);
          selectedList.splice(selectedList.end(), candList, candList.begin(), curIt);
          missings[0] = 0;
          for (size_t i = 0; i < this->energyRefs.size(); i++) {
            curERef = this->energyRefs[i].first + Eoffset;
            curIt = std::lower_bound(candList.begin(), candList.end(), curERef,
                                     [&Eig](size_t j, double x){ return std::real(Eig[j]) < x; });
            while (missings[i+1] > 0) {
              if (curIt != candList.begin()) {
                std::list<size_t>::iterator preIt = curIt;
                preIt--;
                if (curIt == candList.end() or std::abs(Eig[*curIt] - curERef) > std::abs(Eig[*preIt] - curERef)) {
                  curIt = preIt;
                }
              }
              if (curIt == candList.end()) {
                CErr("Error in energy specific Davidson, not enough eigenvalues to select!");
              }
              selectedList.push_back(*curIt);
              curIt = candList.erase(curIt);
              missings[i+1] -= 1;
            }
          }

          sortedIndices.insert(sortedIndices.end(), selectedList.begin(), selectedList.end());
          sortedIndices.insert(sortedIndices.end(), candList.begin(), candList.end());

          // Reorder eigenvectors
          Eigen::Map<
              Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
          > XRMap(XR,nVCur,nVCur);

          for(auto i = 0ul; i < kG*nR; i++){
            size_t ind = sortedIndices[i];
            while(ind < i) ind = sortedIndices[ind];

            if (ind > i) {
              std::swap(Eig[i],Eig[ind]);
              XRMap.col(i).swap(XRMap.col(ind));
            }
          }

#ifdef DEBUG_DAVIDSON
          prettyPrintSmart(std::cout,"HH Davidson XR after swap",XR,nVCur,nExam,nVCur);
#endif
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
#ifdef DAVIDSON_PRINT_TIMING
          std::cout << "      Swap high energy roots took "
          << std::fixed << std::setprecision(6) << tock(SWAPst) << " s." << std::endl;
#endif
        }   


        // print eigenvalues at current iteration
        std::cout << "      - Eigenvalues at the current iteration:" << std::endl;
        for(auto i = 0; i < nExam; i++) {
          std::cout << "        Root " << std::setw(5) << std::right << i << ":"
              << std::right << std::scientific << std::setprecision(10) << std::setw(20) << std::real(Eig[i]);
              
          if( std::is_same<dcomplex,_F>::value ) {
            std::cout << " + " << std::right << std::scientific << std::setprecision(10)
                << std::setw(20) << std::imag(Eig[i]) << " i";
          }
              
          std::cout << std::endl;
        }

        } // Root only

        MPIBCast(Eig, nVCur, 0, this->comm_);
        MPIBCast(XR, nVCur * nVCur, 0, this->comm_);
          
        // Exam Eigenvalues and eigenvectors and do mapping if iter > 0
        // R and S as scratch space to hold full vector old and new repectively
#ifdef DAVIDSON_PRINT_TIMING
        auto MMst = tick();
#endif
        VR->multiply_matrix(0, blas::Op::NoTrans,nExam,nVCur,_F(1.),XR,nVCur,_F(0.),*S, 0);
#ifdef DEBUG_DAVIDSON
        std::cout << "HH Davidson Norm of S: " << S->norm2F(0, nExam) << std::endl;
#endif
#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "      Multiply matrix took "
          << std::fixed << std::setprecision(6) << tock(MMst) << " s." << std::endl;
        } // Root Only
#endif

        for (size_t i = 0; i < nExam; i++) SiConv[i].clear();
        if( iter > 0) {

#ifdef DAVIDSON_PRINT_TIMING
          auto EXAMst = tick();
#endif
          // mapping old vector to new vectors based on overlap
          std::vector<int> StMap(nExam);
          // i -> new state, j -> old state
          int i = 0, j = 0;

#ifdef CQ_HAS_TA
          if (TA::initialized())
            TA::get_default_world().gop.fence();
#endif
        if (isRoot) {
          // overlap = (VR XR)_old ^\dagger * (VR XR)_new
          std::fill_n(Ovlp, nVCur*nVPrev, _F(0.));
          for(size_t i = 0ul; i < nVCur; i++) Ovlp[i + i*nVPrev] = _F(1.);

          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nExam,nVCur,nVPrev,_F(1.),XRPrev,nVPrev,Ovlp,nVPrev,_F(0.),SCR,nExam);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nExam,nExam,nVCur ,_F(1.),SCR,nExam,XR,nVCur,_F(0.),Ovlp,nExam);

          std::fill_n(StMap.begin(),nExam,-1);
          std::copy_n(Ovlp,nExam*nExam,SCR);

          for(size_t i = 0; i < nExam; i++) {
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
        }
        MPIBCast(StMap.data(), nExam, 0, this->comm_);
        MPIBCast(Ovlp, nExam * nExam, 0, this->comm_);
          
#ifdef DEBUG_DAVIDSON
          prettyPrintSmart(std::cout,"HH Davidson States Overlap",Ovlp,nExam,nExam,nExam);
#endif 
      
          std::cout << "\n      - Comparison to the previous iteration: " << std::endl;

#ifdef DAVIDSON_PRINT_TIMING
          auto DIFFst = tick();
#endif
          _F phase;
          double maxDel;  // maximum differece in vectors
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            SiConv[i].prev_index = j;
            if(j < 0) {
              std::cout << "          New root" << std::setw(5) << i << " is brand new" << std::endl;
              continue; 
            } else if (j!=i) { 
              std::cout << "          New root" << std::setw(5) << i
                        << " was old root" << std::setw(5) << j << std::endl; 
            }

            phase = Ovlp[j + i*nExam];
            phase /= std::abs(phase);

            R->scale(-phase, j, 1);
            R->axpy(j, 1, _F(1.), *S, i);
            maxDel = R->maxNormElement(j, 1);

            SiConv[i].diff_eigen_vector = maxDel;
            SiConv[i].eigenVector = maxDel < VecConv;
            std::cout << "        Root " <<  std::right << std::setw(5) << i
                      << " maximum delta is"<<  std::right
                      << std::setw(20) << maxDel << std::endl;
          }
#ifdef DAVIDSON_PRINT_TIMING
          if( isRoot ) {
            std::cout << "        Compute difference took "
            << std::fixed << std::setprecision(6) << tock(DIFFst) << " s." << std::endl;
          } // Root Only
#endif

          // exam eigenvalues
          dcomplex EDiff;
          for(i = 0; i < nExam; i++) {
            j = StMap[i];
            if(j < 0)
              continue;
            EDiff = Eig[i] - EPrev[j];
            SiConv[i].diff_eigen_value = std::abs(EDiff);
            SiConv[i].eigenValue = std::abs(EDiff) < EConv;
          }
#ifdef DAVIDSON_PRINT_TIMING
          if( isRoot ) {
            std::cout << "      Exam eigen took "
            << std::fixed << std::setprecision(6) << tock(EXAMst) << " s." << std::endl;
          } // Root Only
#endif
        }   // Exam eigenvales and eigenvectors

        // Form residue vectors only for unconverged vectors
#ifdef DAVIDSON_PRINT_TIMING
        auto RESIDUEst = tick();
#endif
        std::vector<int> residualStates(nExam);
        size_t nRes = 0ul;
        std::fill_n(residualStates.begin(), nExam, -1);
        std::fill_n(SCR,nVCur*nExam,_F(0.));
        for (auto i = 0ul; i < nExam; i++) {
          if(checkResidueConv or not SiConv[i].hasConverged(checkEigenVectorConv, checkEigenValueConv, false)) {
            residualStates[nRes] = i;
            std::copy_n(XR + i*nVCur,nVCur,SCR + nRes*nVCur);
            nRes++;
          }
        }

#ifdef DEBUG_DAVIDSON
        prettyPrintSmart(std::cout,"HH Davidson Unconverge Transform",SCR,nVCur,nRes,nVCur);
#endif

#ifdef DAVIDSON_PRINT_TIMING
        MMst = tick();
#endif
        // R <- AVR * XR
        AVR->multiply_matrix(0, blas::Op::NoTrans,nRes,nVCur,_F(1.),SCR,nVCur,_F(0.),*R, 0);
#ifdef DEBUG_DAVIDSON
        std::cout << "HH Davidson Norm of R: " << R->norm2F(0, nRes) << std::endl;
#endif

#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "        Multiply matrix took "
          << std::fixed << std::setprecision(6) << tock(MMst) << " s." << std::endl;
        } // Root Only
#endif
            
        // Compute the residue norm and generate perturbbed vectors
#ifdef DAVIDSON_PRINT_TIMING
        auto AXPYst = tick();
#endif
        // R <- AVR * XR - eig_i * (VR * XR_i)
        for (auto i = 0ul; i < nRes; i++) {
          R->axpy(i, 1, -dcomplexTo_F(Eig[residualStates[i]]), *S, residualStates[i]);
        }
#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "        Axpy took "
          << std::fixed << std::setprecision(6) << tock(AXPYst) << " s." << std::endl;
        } // Root Only
#endif
        std::cout << "\n      - Residues of roots: " << std::endl;
        
        if(EigForT) std::fill_n(EigForT,nRes,dcomplex(0.));
        std::fill_n(RelRes, nRes, 0.);
        for (auto i = 0ul; i < nRes; i++) {
          auto j = residualStates[i];
          RelRes[i] = R->norm2F(i, 1);
          SiConv[j].residual_norm = RelRes[i];
          SiConv[j].residual = RelRes[i] < this->convCrit_;
          std::cout << "        Root " << std::setw(5) << std::right << j
                    << " 2nd order lowering " << std::right << std::setw(20) << RelRes[i]*RelRes[i] 
                    << " norm " << std::right << std::setw(20) << RelRes[i] << std::endl;
              
          if(EigForT) this->EigForT[i] = Eig[j];
        }
#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "      Residue took "
          << std::fixed << std::setprecision(6) << tock(RESIDUEst) << " s." << std::endl;
        } // Root Only
#endif


        std::cout << std::endl << "    Convergence check:" << std::endl;
        std::cout << "    " << BannerMid << std::endl;

        std::cout << "    " << std::setw(5) << std::left <<  " Root" << "  ";
        std::cout << std::setw(34) << std::left << "Eigenvalue";
        std::cout << std::setw(4) << std::left << "Prev" << "  ";
        std::cout << std::setw(10) << std::left << "|\u0394Eval|";
        std::cout << std::setw(10) << std::left << "max(\u0394Evec)";
        std::cout << std::setw(10) << std::left << "  |res.|";
        std::cout << std::setw(4) << std::right << "Conv";
        std::cout << std::endl;
        std::cout << "    " << std::setw(5) << std::left <<  " ----" << "  ";
        std::cout << std::setw(34) << std::left << "-----------------";
        std::cout << std::setw(4) << std::left <<  "----" << "  ";
        std::cout << std::setw(10) << std::left << "--------";
        std::cout << std::setw(10) << std::left << "--------";
        std::cout << std::setw(10) << std::left << "--------";
        std::cout << std::setw(3) << std::left << "---";
        std::cout << std::endl;

        for (size_t i = 0; i < nR ; i++){

          std::cout << std::setprecision(12) << std::fixed;
          std::cout << "    " << std::setw(5) << std::right << i << "  ";
          std::cout << std::setw(34) << std::left << std::fixed << Eig[i];
          std::cout << std::setprecision(2) << std::scientific;
          if (SiConv[i].prev_index < 0) {
            std::cout << std::setw(4) << std::right << "new" << "  ";
            std::cout << std::setw(10) << std::left << "   -";
            std::cout << std::setw(10) << std::left << "   -";
          } else {
            std::cout << std::setw(4) << std::right << SiConv[i].prev_index << "  ";
            std::cout << std::setw(10) << std::left << SiConv[i].diff_eigen_value;
            std::cout << std::setw(10) << std::left << SiConv[i].diff_eigen_vector;
          }
          if(checkResidueConv or not SiConv[i].hasConverged(checkEigenVectorConv, checkEigenValueConv, false)) {
            std::cout << std::setw(10) << std::left << SiConv[i].residual_norm;
          } else {
            std::cout << std::setw(10) << std::left << "   -";
          }
          std::cout << std::setw(3) << std::left
          << (SiConv[i].hasConverged(checkEigenVectorConv, checkEigenValueConv, checkResidueConv) ? "YES" : "NO");
          std::cout << std::endl;

        }
        std::cout << "    " << BannerMid << std::endl;

        isConverged = std::all_of(SiConv.begin(), SiConv.begin()+nR, [&] (const IterDiagConvStatus &i) {
          return i.hasConverged(checkEigenVectorConv, checkEigenValueConv, checkResidueConv); });
        if(isConverged or nVCur >= MSS) {

          double DavidsonDur = tock(DavidsonSt);
          double perLT = LTdur * 100 / DavidsonDur;

          std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur
          << " s  ( " << perLT << " % LT )" << std::endl;
          break;
        }

        nDo = 0;
        for (size_t i = 0ul; i < nRes; i++) {
          if(not SiConv[residualStates[i]].hasConverged(checkEigenVectorConv, checkEigenValueConv, checkResidueConv)) {
            if (nDo < i) {
              R->set_data(nDo, 1, *R, i, true);
              if(EigForT) EigForT[nDo] = EigForT[i];
            }
            nDo++;
          }
        }

        if(nVCur+nDo > MSS)  nDo = MSS - nVCur;

#ifdef DEBUG_DAVIDSON
        S->print(std::cout, "S before preCond");
#endif
#ifdef DAVIDSON_PRINT_TIMING
        auto PRECONDst = tick();
#endif
        this->preCondNoShift_(nDo,*R,*R);
#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "      Precondition took "
          << std::fixed << std::setprecision(6) << tock(PRECONDst) << " s." << std::endl;
        } // Root Only
#endif

#ifdef DEBUG_DAVIDSON
        S->print(std::cout, "S after preCond");
#endif
        // Append the new vectors to VR and orthogoalize against existing ones 
        // Also update the dimensions and save the XRPrev
#ifdef DAVIDSON_PRINT_TIMING
        auto COPYst = tick();
#endif
        VR->set_data(nVCur, nDo, *R, 0, true);
        std::swap(R,S);
        std::copy_n(XR,nVCur*nVCur,XRPrev);
        std::copy_n(Eig,nVCur,EPrev);
        nVPrev = nVCur;
#ifdef DAVIDSON_PRINT_TIMING
        if( isRoot ) {
          std::cout << "      Copy residue took "
          << std::fixed << std::setprecision(6) << tock(COPYst) << " s." << std::endl;
        } // Root Only
#endif
        if(this->DoLeftEigVec) {
          CErr("Do Left Eig Vec is not implemented yet");
        } else {
          // disable printing from GramSchmidt

#ifdef DAVIDSON_PRINT_TIMING
          auto GramSchmidtSt = tick();
#endif

#ifndef DEBUG_DAVIDSON
          std::cout.setstate(std::ios_base::failbit);
#endif

#ifdef DEBUG_DAVIDSON
          VR->print(std::cout, "VR before GramSchmidt");
#endif
          nVCur = VR->GramSchmidt(0, nVCur,nDo,GramSchmidt_NRe,GramSchmidt_eps);
#ifdef DEBUG_DAVIDSON
          VR->print(std::cout, "VR after GramSchmidt");
#endif
          
#ifndef DEBUG_DAVIDSON
          std::cout.clear();
#endif

#ifdef DAVIDSON_PRINT_TIMING
          if( isRoot ) {
            std::cout << "      Gram-Schmidt took "
            << std::fixed << std::setprecision(6) << tock(GramSchmidtSt) << " s." << std::endl;
          } // Root Only
#endif
          
          nDo   = nVCur - nVPrev;
        }
         
        if(iter+1 == this->whenSc) { 
          nExam = nR; 
          nDo   = std::min(nDo, nExam);
          nVCur = nVPrev + nDo;
          this->kG = 1;
        }
        
        double DavidsonDur = tock(DavidsonSt);
        double perLT = LTdur * 100 / DavidsonDur;
    
        std::cout << "\n      - DURATION = " << std::setprecision(8) << DavidsonDur 
          << " s  ( " << perLT << " % LT )" << std::endl;
    
        if(nDo == 0) {
          std::cout << "  * All new vectors in GramSchmidt are linear dependent of existing vectors!" << std::endl;
          isConverged = convOnGramSchmidt;
          break;
        }
//      } // Root Only
    
    } // Davidson iteration    

//    if( isRoot ) {

      // move data before exit runMicro      
      std::copy_n(Eig,nR,this->eigVal_);
#ifdef DAVIDSON_PRINT_TIMING
      auto FINALst = tick();
#endif
      VR->multiply_matrix(0, blas::Op::NoTrans,nR,nVCur,_F(1.),XR,nVCur,_F(0.),*this->VR_, 0);
#ifdef DAVIDSON_PRINT_TIMING
      if( isRoot ) {
        std::cout << "      Final linear combination took "
        << std::fixed << std::setprecision(6) << tock(FINALst) << " s." << std::endl;
      } // Root Only
#endif
      //size_t nVSave = isConverged ? nR: nG;  
      //std::copy_n(Eig,nVSave,this->eigVal_);
      //blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVSave,nVCur,_F(1.),VR,N,XR,nVCur,_F(0.),this->VR_,N);

    if( isRoot ) {
      std::cout << "\n  * ";
      if( isConverged and nDo !=0)
        std::cout << "Davidson converged in " << iter+1 << " iterations." << std::endl;
      else if (isConverged) {
        double maxResDel = *std::max_element(RelRes, RelRes+nR);
        std::cout << "Davidson converged in " << iter+1 << " iterations "
                  << "by the criteria of GramSchmidt threshold. Residual converge below "
                  << std::setw(20) << std::right << maxResDel << std::endl;
      } else
        std::cout << "Davidson Failed to Converged in " << iter+1 << " Iterations." << std::endl;
    } // Root Only 

    // Free Scratch space

    if(XR)       CQMemManager::get().free(XR);    
    if(XRPrev)   CQMemManager::get().free(XRPrev);
    if(SubA)     CQMemManager::get().free(SubA);
    if(Eig)      CQMemManager::get().free(Eig);
    if(EPrev)    CQMemManager::get().free(EPrev); 
    if(Ovlp)     CQMemManager::get().free(Ovlp);
    if(SCR)      CQMemManager::get().free(SCR);
    if(XL)       CQMemManager::get().free(XL);    
    
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

    if (not Guess or Guess->size() < this->nGuess_)
      Guess = this->vecGen_(this->nGuess_);

    // NO MPI
    // ROOT_ONLY(this->comm_);
    this->Guess->set_data(0, this->nGuess_, *this->VR_, 0, true);
  } // Davidson::restart
  
}; // namespace ChronusQ

#endif
