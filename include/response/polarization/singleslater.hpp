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

#ifdef CQ_ENABLE_MPI
#include <scalapackpp/lascl.hpp>
#endif

namespace ChronusQ {


  template <typename MatsT, typename IntsT>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::configOptions() {

    bool isRoot = MPIRank(this->comm_) == 0;

    if( this->genSettings.printLevel > 0  and isRoot ) {
      std::cout << "  * PROPAGATOR = ParticleHolePropagator" <<
        std::endl << std::endl;

      std::cout << 
        "  * PERFORMING RESPONSE OPTION CONFIGURATION FOR SINGLE SLATER "
                << "WAVE FUNCTION\n";
    }
    SingleSlater<MatsT,IntsT>& ss =dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);
    // Set the internal NSingleDim
    this->nSingleDim_ = getNSingleDim(this->genSettings.doTDA);

    // Direct logic
    if( not this->genSettings.formFullMat ) {

      // Turn off doAPB_AMB and doReduced
      doAPB_AMB = false;
      doReduced = false;

      // Turn off distribution
      this->genSettings.formMatDist     = false;
      this->genSettings.distMatFromRoot = false;

    }


    // Stab or NR
    if( doStab or doNR ) {

      incMet = false; // Turn off metric
      this->genSettings.evalProp = false; // No property evaluation

      // Set misc settings to override stupidity 
      if( doStab ) this->genSettings.jobType = RESIDUE;
      else {
        this->genSettings.jobType = FDR;
        this->genSettings.bOps = { Brillouin };
      }

    }



    // Logic handling for reduced
    if( doReduced ) {
      doAPB_AMB = true; // Reduced -> A+B/A-B

      this->resSettings.needVL = true; // We need left vectors (X-Y)

      // Change to hermetian problem because of BSEPACK's scheme
      if( this->genSettings.doFull and not this->genSettings.doTDA )
        this->genSettings.matIsHer = true;
    }



    // Logic for full RESIDUE problem
    if( this->genSettings.doFull and this->genSettings.jobType == RESIDUE ) {

      if( this->genSettings.printLevel > 0 and isRoot)
        std::cout << "  * SETTING NROOTS = NSINGLEDIM\n";

      // Nroots = NSingleDim
      this->resSettings.nRoots = this->nSingleDim_;

      // Account for pairing if not reduced problem
      if( not this->genSettings.doTDA and not doReduced and incMet ) 
        this->resSettings.nRoots /= 2;

    }



    // If the metric is not included, automatically an hermetian problem:
    if( not incMet ) this->genSettings.matIsHer = true;

  
    // GIAO handling
    if( this->ref_->basisSet().basisType == COMPLEX_GIAO ) {

      auto remove_op = []( std::vector<ResponseOperator> &ops,
                           ResponseOperator op ) {

        ops.erase( std::remove(ops.begin(), ops.end(), op), ops.end() );

      };


      // Remove certain properties from eval list
      remove_op( this->genSettings.aOps, VelElectricDipole     );
      remove_op( this->genSettings.aOps, VelElectricQuadrupole );
      remove_op( this->genSettings.aOps, VelElectricOctupole   );
      remove_op( this->genSettings.aOps, MagneticQuadrupole    );

      remove_op( this->genSettings.bOps, VelElectricDipole     );
      remove_op( this->genSettings.bOps, VelElectricQuadrupole );
      remove_op( this->genSettings.bOps, VelElectricOctupole   );
      remove_op( this->genSettings.bOps, MagneticQuadrupole    );

    }



    // FDR problem settings
    if( this->genSettings.jobType == FDR ) {

      // Hermetian operators
      bool needsP = 
        std::any_of(this->genSettings.aOps.begin(),
        this->genSettings.aOps.end(), isHerOp);

      // Anti-Hermetian operators
      bool needsQ = 
        std::any_of(this->genSettings.aOps.begin(),
        this->genSettings.aOps.end(), isAntiHerOp);

      this->fdrSettings.needP = needsP;
      this->fdrSettings.needQ = needsQ;

      // If reduced and needs both P and Q, expand to A+B/A-B
      if( needsP and needsQ and doReduced) {
        if( this->genSettings.printLevel > 0 and isRoot )
          std::cout << "  * EXPANDING FROM REDUCED PROBLEM TO A+B/A-B DUE TO LHS FDR OPERATORS\n";
        doReduced = false;
      }

      // Determine number of RHS
      this->fdrSettings.nRHS = 0;
      for(auto & op : this->genSettings.bOps) 
        this->fdrSettings.nRHS += OperatorSize[op];
        
      if( this->genSettings.printLevel > 0 and isRoot )
        std::cout << "  * FDR MODULE FOUND " << this->fdrSettings.nRHS 
                  << " RIGHT HAND SIDES\n";
    } 

    if( this->genSettings.printLevel > 0 and isRoot ) {
      if( doReduced )      std::cout << "  * DOING REDUCED PROBLEM\n";
      else if( doAPB_AMB ) std::cout << "  * DOING A+B / A-B PROBLEM\n";
      else                 std::cout << "  * DOING FULL DIMENSIONAL PROBLEM\n";

      if( not this->genSettings.formFullMat )
        std::cout << "  * TENSOR CONTRACTIONS WILL BE FORMED DIRECTLY" 
          << std::endl;
    }


    if( (this->genSettings.printLevel > 0 and isRoot) and 
        (this->genSettings.aOps.size() or this->genSettings.bOps.size()) 
        and (not doStab and not doNR) ) {

      std::cout << "RANK " << MPIRank(this->comm_) << std::endl;;

      std::cout << "  * RESPONSE Will Compute the Following Properties:\n";

      auto print_op = [&](ResponseOperator op) {

        std::cout << "      * ";
        switch (op) {

          case LenElectricDipole: 
            std::cout << "Electric Dipole (Length)\n";
            break;

          case LenElectricQuadrupole: 
            std::cout << "Electric Quadrupole (Length)\n";
            break;

          case LenElectricOctupole: 
            std::cout << "Electric Octupole (Length)\n";
            break;

          case VelElectricDipole: 
            std::cout << "Electric Dipole (Velocity)\n";
            break;

          case VelElectricQuadrupole: 
            std::cout << "Electric Quadrupole (Velocity)\n";
            break;

          case VelElectricOctupole: 
            std::cout << "Electric Octupole (Velocity)\n";
            break;

          case MagneticDipole: 
            std::cout << "Magnetic Dipole\n";
            break;

          case MagneticQuadrupole: 
            std::cout << "Magnetic Quadrupole\n";
            break;default:
            break;
        }


      };

      if( this->genSettings.aOps.size() ) {
        std::cout << "    * AOPS:\n";
        for(auto k = 0; k < this->genSettings.aOps.size(); k++)
          print_op(this->genSettings.aOps[k]);

        std::cout << "\n" << std::endl;
      }

      if( this->genSettings.bOps.size() and 
          this->genSettings.jobType == FDR ) {
        std::cout << "    * BOPS:\n";
        for(auto k = 0; k < this->genSettings.bOps.size(); k++)
          print_op(this->genSettings.bOps[k]);

        std::cout << "\n" << std::endl;
      }

    }

  };


  template <typename MatsT, typename IntsT>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::constructShifts() {

    ResponseTBase<MatsT>::constructShifts();

    if( doReduced ) {

      for(auto &omega : this->fdrResults .shifts) omega *= omega;
      for(auto &omega : this->dfdrResults.shifts) omega *= omega;

    }

  };

  template <typename MatsT, typename IntsT>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::postLinearSolve() {

    ResponseTBase<typename SingleSlater<MatsT, IntsT>::value_type>::postLinearSolve();

    ROOT_ONLY(this->comm_);


    if( not doReduced ) return;

    bool noDamp = (this->fdrSettings.dampFactor == 0.) and 
                  not this->fdrSettings.forceDamp;

    size_t N      = this->nSingleDim_;
    size_t nRHS   = this->fdrSettings.nRHS;
    std::vector<double> freq = this->fdrSettings.bFreq;

    size_t iOff = 0;
    for(auto op : this->genSettings.bOps) {

      bool scaleOp = isHerOp(op) and this->fdrSettings.needQ;
      scaleOp = scaleOp or (isAntiHerOp(op) and this->fdrSettings.needP);

      size_t nOp = OperatorSize[op];

      //std::cerr << "OP = " << op << " " << std::boolalpha 
      //  << scaleOp << std::endl;

      if( scaleOp ) {
        for(auto iO = 0; iO < freq.size(); iO++)
          if( noDamp ) {
            blas::scal(N*nOp,MatsT(freq[iO]),this->fdrResults.SOL + (iO*nRHS + iOff)*N,1);
          } else {
            blas::scal(N*nOp,dcomplex(freq[iO],this->fdrSettings.dampFactor),
              this->dfdrResults.SOL + (iO*nRHS + iOff)*N,1);
          }
      }
      iOff += nOp;
    };

  };

  template <typename MatsT, typename IntsT>
  template <typename U>
  std::vector< std::pair< std::pair<int,int>, U > >
    PolarizationPropagator<SingleSlater<MatsT, IntsT>>::getMOContributions(U *V, 
        double tol) {

    std::vector< std::pair< std::pair<int,int>, U > > moCont;

    auto& ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);
    int nOAVA = ss.nOA * ss.nVA;
    int nOBVB = ss.nOB * ss.nVB;
    int nOV   = ss.nO  * ss.nV;

    int N  = (ss.nC == 1) ? nOAVA  : nOV;
    int NV = (ss.nC == 1) ? ss.nVA : ss.nV;
    int NO = (ss.nC == 1) ? ss.nOA : ss.nO;

    for(auto j = 0; j < N; j++) {

      int i = j / NV;
      int a = (j % NV) + NO;

      if( std::abs(V[j]) > tol ) moCont.push_back( { {a,i}, V[j] } );

    }

    if( ss.nC == 1 ) {

      N  = nOBVB; 
      NV = ss.nVB;
      NO = ss.nOB; 

      for(auto j = 0; j < N; j++) {

        int i = j / NV;
        int a = (j % NV) + NO;

        if( std::abs(V[j+nOAVA]) > tol ) 
          moCont.push_back( { {-a,-i}, V[j+nOAVA] } );

      }

    }

    return moCont;

  };


  template <typename MatsT, typename IntsT>
  template <typename U>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::printResMO_impl(
    std::ostream &out, size_t nRoots, double *W_print,
    std::vector<std::pair<std::string,double *>> data, U* VL, U* VR) {

    SingleSlater<MatsT,IntsT>& ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);
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

      auto xCont = getMOContributions(VR+iRt*this->nSingleDim_,1e-1);
      decltype(xCont) yCont;
      if( not this->genSettings.doTDA ) 
        yCont = getMOContributions(VL+iRt*this->nSingleDim_,1e-1);
      
      // MO contributions
      out << "    MO Contributions:\n";
      for(auto &c : xCont) {

        char spinLabel = (c.first.first > 0) ? 
                           ((this->ref_->nC == 1) ? 'A' : ' ') : 'B';
      
        out << "      ";
        out << std::setw(4) << std::right 
            << std::abs(c.first.second) + 1 << spinLabel << " -> ";
        out << std::setw(4) << std::right 
            << std::abs(c.first.first) + 1<< spinLabel;

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
      for(auto &c : yCont) {

        char spinLabel = (c.first.first > 0) ? 
                           ((this->ref_->nC == 1) ? 'A' : ' ') : 'B';
      
        out << "      ";
        out << std::setw(4) << std::right 
            << std::abs(c.first.second) + 1 << spinLabel << " <- ";
        out << std::setw(4) << std::right 
            << std::abs(c.first.first) + 1 << spinLabel;

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


      out << "\n\n";

    }

  }
  

  template <typename MatsT, typename IntsT>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::eigVecNorm() {

    if(this->genSettings.doTDA) return;

    auto N  = this->nSingleDim_;
    size_t nRoots = this->resSettings.nRoots;

    MatsT* V = this->resResults.VR;

    if( doStab ) {

      for(auto k = 0; k < nRoots; k++) {

        MatsT tnorm = blas::nrm2(N/2, V + k*N, 1);
        if( std::abs(tnorm) < 1e-08 )
          tnorm = blas::nrm2(N/2, V + k*N + N/2, 1);

        blas::scal(N,1./tnorm,V + k*N,1);

      }

      return;

    } else if(not incMet) return;
   


    if( this->genSettings.doFull and not doReduced ) 
      V += N*N/2;


    std::function<MatsT(MatsT*)> tdInner = 
      [&](MatsT* Vc) { 
        return std::sqrt(std::abs(blas::dot(N,Vc,1,Vc,1)));
      };


    if( doAPB_AMB )
      tdInner = 
        [&](MatsT* Vc) { 
          return std::sqrt(std::abs(blas::dot(N/2,Vc,1,Vc+N/2,1) + 
                 blas::dot(N/2,Vc+N/2,1,Vc,1)));
        };
    else 
      tdInner = 
        [&](MatsT* Vc) { 
          return std::sqrt(std::abs(blas::dot(N/2,Vc,1,Vc,1) - 
                 blas::dot(N/2,Vc+N/2,1,Vc+N/2,1)));
        };









    
    std::function<void(size_t,size_t,MatsT*,size_t,MatsT*,size_t,MatsT*)> tdMatInner =
      [&](size_t i, size_t j, MatsT* Vi, size_t LDVi, MatsT* Vj, size_t LDVj, 
        MatsT* inner){ 

        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N,MatsT(1.),Vi,LDVi,Vj,LDVj,MatsT(0.),inner,i);

      };


    if( doAPB_AMB )
      tdMatInner = 
        [&](size_t i, size_t j, MatsT* Vi, size_t LDVi,
          MatsT* Vj, size_t LDVj, MatsT* inner){ 
  
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N/2,MatsT(1.),Vi    ,LDVi,Vj+N/2,LDVj,MatsT(0.),inner,i);
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N/2,MatsT(1.),Vi+N/2,LDVi,Vj    ,LDVj,MatsT(1.),inner,i);
  
        };
    else
      tdMatInner = 
        [&](size_t i, size_t j, MatsT* Vi, size_t LDVi,
          MatsT* Vj, size_t LDVj, MatsT* inner){ 

          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N/2,MatsT(1.) ,Vi    ,LDVi,Vj    ,LDVj,MatsT(0.),inner,i);
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,i,j,N/2,MatsT(-1.),Vi+N/2,LDVi,Vj+N/2,LDVj,MatsT(1.),inner,i);

        };


    if( doReduced ) {

      for(auto k = 0; k < N; k++) 
        this->resResults.W[k] = std::sqrt(this->resResults.W[k]);

    }



    if( not this->genSettings.matIsHer ) {
      if( not doReduced )
        GramSchmidt(N,0,nRoots,V,N,tdInner,tdMatInner,1);
      else
        GramSchmidt(N,0,nRoots,this->resResults.VL,N,V,N,tdMatInner,
          1);
    } else {

      if( not this->fullMatrix_ ) 
        CErr("Reduced space normalization NYI for non-FULL");

      MatsT* M = this->fullMatrix_;
      MatsT* K = this->fullMatrix_ + N;
      // Use BSEPACK's odd normalization scheme...
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,MatsT(1.),M,2*N,
        this->resResults.VR,N,MatsT(0.),this->resResults.VL,N);

      // Triangular linear solve
      blas::trsm(blas::Layout::ColMajor,blas::Side::Left,blas::Uplo::Lower,blas::Op::ConjTrans,blas::Diag::NonUnit,
        N,N,MatsT(1.),M,2*N,this->resResults.VR,N);

      for(auto k = 0; k < N; k++)
      for(auto j = 0; j < N; j++) {

        auto tmp = this->resResults.VR[j + k*N];
        auto omega = this->resResults.W[k];

        this->resResults.VR[j + k*N] = 
          0.5/std::sqrt(omega) * (tmp*omega + this->resResults.VL[j + k*N]);
        this->resResults.VL[j + k*N] = 
          0.5/std::sqrt(omega) * (tmp*omega - this->resResults.VL[j + k*N]);

      }

      this->blockTransform(N,N,std::sqrt(0.5),
        this->resResults.VR,N,this->resResults.VL,N);

    }


    
    if( this->genSettings.doFull and not doReduced ) { 
      // Move over the "EXCITATION" roots to the begning of the storage
      SetMat('N',N/2,1,1.,this->resResults.W + (N/2),N,this->resResults.W ,N);
      SetMat('N',N,N/2,MatsT(1.),V                ,N,this->resResults.VR,N);
    }

  };



  template <typename MatsT, typename IntsT>
  template <typename U>
  void PolarizationPropagator< SingleSlater<MatsT, IntsT> >::formLinearTrans_incore_impl( 
    std::vector<RESPONSE_CONTRACTION<U>> cList,
    SINGLESLATER_POLAR_COPT op ){ 
  
    bool isDist = this->genSettings.isDist(); 

    /*
    // FIXME MPI : This is ad hoc logic, assumes that if the cases when
    // matrix is not distributed that only the root process will be
    // calling this function
    bool needFullMat = not bool(this->fullMatrix_);
#ifdef CQ_ENABLE_MPI
    if( this->genSettings.isDist() ) MPIBCast(&needFullMat,1,0,this->comm_);
#endif
    if(needFullMat) formFullMatrix(); 
    */

    size_t N = this->nSingleDim_;


    if( this->doReduced and op == SINGLESLATER_POLAR_COPT::FULL ) {

      if( this->fdrSettings.needP ) op = KM;
      else                          op = MK;

    }

    bool mContract  = op == SINGLESLATER_POLAR_COPT::M;
    bool kContract  = op == SINGLESLATER_POLAR_COPT::K;
    bool mkContract = op == SINGLESLATER_POLAR_COPT::MK;
    bool kmContract = op == SINGLESLATER_POLAR_COPT::KM;

    size_t maxNVec = 
      std::max_element(cList.begin(),cList.end(),
        [](RESPONSE_CONTRACTION<U> l, RESPONSE_CONTRACTION<U> r) {
          return l.nVec < r.nVec;
        }
      )->nVec;


    // (Local) Dims
    int64_t MLoc,NLoc;
    if( this->doReduced ) {
      MLoc = 2*N; NLoc = N;
    } else if(this->doAPB_AMB) {
      MLoc = N; NLoc = N/2;
    } else {
      MLoc = N; NLoc = N;
    }

#ifdef CQ_ENABLE_MPI
    if( isDist )
      std::tie(MLoc,NLoc) = this->fullMatGrid_->get_local_dims(MLoc,NLoc);
#endif


    size_t nVecLoc(0), SCRMLoc(0);
    if( mkContract or kmContract ) {
      nVecLoc = maxNVec, SCRMLoc = N;
#ifdef CQ_ENABLE_MPI
      if( isDist )
        std::tie(SCRMLoc,nVecLoc) = 
          this->fullMatGrid_->get_local_dims(SCRMLoc,nVecLoc);
#endif
    }



    U* SCR = nullptr;
    if( nVecLoc and SCRMLoc ) {
      SCR = CQMemManager::get().malloc<U>(maxNVec * N);
      std::fill_n(SCR,maxNVec*N,U(0.));
    }


    // ScaLAPACK DESC
#ifdef CQ_ENABLE_MPI
    scalapackpp::scalapack_desc descMem, descSCR;
    if( isDist ) {

      if( this->doReduced )
        descMem = this->fullMatGrid_->descinit_noerror(2*N,N,MLoc);
      else if( this->doAPB_AMB )
        descMem = this->fullMatGrid_->descinit_noerror(N,N/2,MLoc);
      else
        descMem = this->fullMatGrid_->descinit_noerror(N,N,MLoc);

      if( mkContract or kmContract ) 
        descSCR = this->fullMatGrid_->descinit_noerror(N,maxNVec,SCRMLoc);

    }
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
      auto *AX    = SCR ? SCR : mat.AX;

#ifdef CQ_ENABLE_MPI
      scalapackpp::scalapack_desc DescAX = SCR ? descSCR : mat.DescAX;
#endif

      if( this->doReduced ) {

        if( mContract or kmContract ) 
#ifdef CQ_ENABLE_MPI
          if( isDist )
            Gemm_MPI('N','N',N,nVec,N,U(1.),FM,1,1,descMem,X,1,1,mat.DescX,
              U(0.),AX,1,1,DescAX);
          else
#endif
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_,2*N,X,N,U(0.),AX,N);

        if( kContract or mkContract ) 
#ifdef CQ_ENABLE_MPI
          if( isDist )
            Gemm_MPI('N','N',N,nVec,N,U(1.),FM,N+1,1,descMem,X,1,1,mat.DescX,
              U(0.),AX,1,1,DescAX);
          else
#endif
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_+N,2*N,X,N,
              U(0.),AX,N);

        AX     = mat.AX;
#ifdef CQ_ENABLE_MPI
        DescAX = mat.DescAX;
#endif
        

        if( mkContract )
#ifdef CQ_ENABLE_MPI
          if( isDist )
            Gemm_MPI('N','N',N,nVec,N,U(1.),FM,1,1,descMem,SCR,1,1,descSCR,
              U(0.),AX,1,1,DescAX);
          else
#endif
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_,2*N,SCR,N,
              U(0.),AX,N);

        if( kmContract )
#ifdef CQ_ENABLE_MPI
          if( isDist )
            Gemm_MPI('N','N',N,nVec,N,U(1.),FM,N+1,1,descMem,SCR,1,1,descSCR,
              U(0.),AX,1,1,DescAX);
          else
#endif
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_+N,2*N,SCR,N,
              U(0.),AX,N);

      } else if( this->doAPB_AMB ) {

        if( isDist ) {

#ifdef CQ_ENABLE_MPI
          Gemm_MPI('N','N',N/2,nVec,N/2,U(1.),FM,(N/2)+1,1,descMem,
            X,(N/2)+1,1,mat.DescX, U(0.),AX,1,1,DescAX);
          Gemm_MPI('N','N',N/2,nVec,N/2,U(1.),FM,1,1,descMem,
            X,1,1,mat.DescX, U(0.),AX,(N/2)+1,1,DescAX);
#endif

        } else {

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N/2,nVec,N/2,U(1.),this->fullMatrix_+(N/2),N,
            X + (N/2),N, U(0.),AX,N);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N/2,nVec,N/2,U(1.),this->fullMatrix_,N,X,N,
            U(0.),AX+(N/2),N);

        }

      } else {

#ifdef CQ_ENABLE_MPI
        if( isDist )
          Gemm_MPI('N','N',N,nVec,N,U(1.),FM,1,1,descMem,X,1,1,mat.DescX,
            U(0.),AX,1,1,DescAX);
        else
#endif
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,nVec,N,U(1.),this->fullMatrix_,N,X,N,U(0.),AX,N);

        // Undo the metric scaling if need be if the matrix is
        // non hermetian and the incMet flag is turned off
        // i.e. if the TD problem and no metric
        if( not this->genSettings.matIsHer and not this->incMet ) {
#ifdef CQ_ENABLE_MPI
          if( isDist )
            scalapackpp::plascl(scalapackpp::MatrixType::Full,-1.,1.,N/2,nVec,AX,N/2+1,1,DescAX);
          else
#endif
            SetMat('N',N/2,nVec,U(-1.),AX + N/2,N,AX + N/2,N);
        }


      }

    }
  
    if( SCR ) CQMemManager::get().free(SCR);
    if( FM and (reinterpret_cast<MatsT*>(FM) != this->fullMatrix_) ) 
      CQMemManager::get().free(FM);
  
  };

  template <>
  void PolarizationPropagator< SingleSlater<double,double> >::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList,
    SINGLESLATER_POLAR_COPT op ){ 

    formLinearTrans_incore_impl(cList,op);

  }

  template <>
  void PolarizationPropagator<SingleSlater<dcomplex,double>>::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList,
    SINGLESLATER_POLAR_COPT op ){ 

    CErr("SOMETHING HAS GONE HORRIBLY WRONG formLinearTrans_incore");

  }

  template <>
  void PolarizationPropagator<SingleSlater<dcomplex,dcomplex>>::formLinearTrans_incore( 
    std::vector<RESPONSE_CONTRACTION<double>> cList,
    SINGLESLATER_POLAR_COPT op ){ 

    CErr("SOMETHING HAS GONE HORRIBLY WRONG formLinearTrans_incore");

  }



  template <typename MatsT, typename IntsT>
  MatsT * PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formFullMatrix() {
  
    bool isRoot = MPIRank(this->comm_) == 0;
    
    if( isRoot ) std::cout << "\n  * FORMING FULL MATRIX\n";
    
    ProgramTimer::tick("Full Hessian Form");

    this->nSingleDim_ = this->getNSingleDim(this->genSettings.doTDA);
    size_t N = this->nSingleDim_;

    if( not this->genSettings.doTDA and this->doReduced ) N *= 2;
  
    if( std::is_same<MatsT,dcomplex>::value and this->doAPB_AMB )
      CErr("COMPLEX + (A+B)/(A-B) NOT VALID!!!");

    size_t nForm  = this->genSettings.doTDA ? N : N/2;
    size_t nStore = (this->doAPB_AMB and not this->genSettings.doTDA) ? N/2 : N;
    size_t nFormP = nForm; // Persistant for matrix distribution
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

        //auto lc = formGrid->localFromGlobal(0,i);
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
    RC_coll<MatsT> cList(1);
    cList.back().nVec = nForm;
    cList.back().N    = N;
    cList.back().X    = V;
    cList.back().AX   = HV;

    ProgramTimer::timeOp("Hessian Contraction", [&](){
      formLinearTrans_direct(matComm,cList,FULL);
    });


    CQMemManager::get().free(V); // Free up some memory

    ProgramTimer::timeOp("Hessian Transform", [&](){
      if( not this->genSettings.doTDA ) {


        if( this->doAPB_AMB ) {

          if( HV ) this->blockTransform(N/2 ,nForm, 1., HV, N, HV+N/2, N);

        } else if( not this->genSettings.formMatDist  and isRootMatComm ) {

          MatsT fact = incMet ? -1. : 1.;
          // Place upper right "B"
          SetMat('R',N/2,N/2,fact,HV + (N/2),N,HV + N*(N/2),N);

          // Place lower right "-Conj(A)"
          SetMat('R',N/2,N/2,fact,HV,N,HV + (N+1)*(N/2),N);

          // Negate the bottom half
          // SetMat('N',N/2,N,T(-1.),HV + (N/2),N,HV + (N/2),N);

        }

      }

    });
  


    ProgramTimer::tick("Hessian Distribute");
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

        // Scale / conjugate if necessary Conj(B) -> -Conj(B)
      //if( not this->doAPB_AMB )
      //  SetMat('N',N/2,nForm,MatsT(-1.),HV + (N/2),N,HV + (N/2),N);

        // Redistribute
        scalapackpp::wrappers::pgemr2d(N,nFormP,HV,1,1,curDesc,this->fullMatrix_,1,1,
          this->descFullMat_,this->fullMatGrid_->grid().context());

        // Place other blocks if need be
        if( not this->doAPB_AMB ) {

          // A -> -Conj(A) / -Conj(B) -> B
          SetMat('R',N,nForm,MatsT(-1.),HV,N,HV,N);

          // Lower right A
          scalapackpp::wrappers::pgemr2d(N/2,N/2,HV,1,1,curDesc,
            this->fullMatrix_,N/2+1,N/2+1, this->descFullMat_,
            this->fullMatGrid_->grid().context());

          // Upper right B
          scalapackpp::wrappers::pgemr2d(N/2,N/2,HV,N/2+1,1,curDesc,
            this->fullMatrix_,1,N/2+1, this->descFullMat_,
            this->fullMatGrid_->grid().context());

        }

      }


      if( HV ) CQMemManager::get().free(HV);

    } else 
#endif
      this->fullMatrix_ = HV;
    
		ProgramTimer::tock("Hessian Distribute");

    //if( isDist ) CErr();


    if( matComm != this->comm_ ) MPICommFree(matComm);

    if( this->savFile.exists() and isRoot and not isDist ) {

      this->savFile.safeWriteData("/RESP/FULLMATRIX",HV,{nStore,N});

    }

    ProgramTimer::tock("Full Hessian Form");
    if( isRoot ) {

      CQSecond durFull = ProgramTimer::getDurationTotal<CQSecond>(
        "Full Hessian Form");
      CQSecond durContract = ProgramTimer::getDurationTotal<CQSecond>(
        "Hessian Contraction");
      CQSecond durTrans = ProgramTimer::getDurationTotal<CQSecond>(
        "Hessian Transform");
      CQSecond durDist = ProgramTimer::getDurationTotal<CQSecond>(
        "Hessian Distribute");

      std::cout << "  * TIMINGS\n";
      std::cout << "    * TOTAL     " << durFull.count() << " s\n";
      std::cout << "    * CONTRACT  " << durContract.count() << " s\n";
      if( not this->genSettings.doTDA )
        std::cout << "    * TRANS     " << durTrans.count() << " s\n";
#ifdef CQ_ENABLE_MPI
      if( isDist )
        std::cout << "    * DIST      " << durDist.count() << " s\n";
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
  MatsT * PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formFullFromMemory() {

    bool isDist = this->genSettings.isDist(); 

    // FIXME MPI : This is ad hoc logic, assumes that if the cases when
    // matrix is not distributed that only the root process will be
    // calling this function
    bool needFullMat = not bool(this->fullMatrix_);
#ifdef CQ_ENABLE_MPI
    if( this->genSettings.isDist() ) MPIBCast(needFullMat,0,this->comm_);
#endif
    if(needFullMat) this->formFullMatrix(); 

    if( doReduced ) {


      size_t N = this->nSingleDim_;
      int64_t NLoc = N, MLoc = N;

#ifdef CQ_ENABLE_MPI
      if( isDist )
        std::tie(MLoc,NLoc) = this->fullMatGrid_->get_local_dims(N,N);
#endif


      MatsT* full = CQMemManager::get().malloc<MatsT>(NLoc*MLoc);

      // Serial offsets
      MatsT* M = this->fullMatrix_;
      MatsT* K = this->fullMatrix_ + N;

      // ScaLAPACK offsets
      MatsT* FM = this->fullMatrix_;
      int64_t IM = 1  , JM = 1;
      int64_t IK = N+1, JK = 1;

#ifdef CQ_ENABLE_MPI
      // ScaLAPACK DESC
      scalapackpp::scalapack_desc descMem, descFull;
      if( isDist ) {

        int64_t NLocMem, MLocMem;
        std::tie(MLocMem,NLocMem) = this->fullMatGrid_->get_local_dims(2*N,N);

        descMem  = this->fullMatGrid_->descinit_noerror(2*N,N,MLocMem);
        descFull = this->fullMatGrid_->descinit_noerror(N,N  ,MLoc   );

      }
#endif

      if( this->genSettings.jobType == RESIDUE ) {

        if( this->genSettings.isDist() )
          CErr("Reduced Residue Response + ScaLAPACK NYI");

        lapack::potrf(lapack::Uplo::Lower,N,M,2*N); // Cholesky of M
        for(auto k = 0; k < N; k++)
        for(auto j = 0; j < k; j++)
          M[j + 2*k*N] = 0.;

        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,N,N,N,MatsT(1.),M   ,2*N,K,2*N,MatsT(0.),full,N  );
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,MatsT(1.),full,N  ,M,2*N,MatsT(0.),K   ,2*N);

        MatAdd('N','C',N,N,MatsT(0.5),K,2*N,MatsT(0.5),K,2*N,full,N);


      } else if( this->fdrSettings.needP ) {

#ifdef CQ_ENABLE_MPI
        if( isDist )
          Gemm_MPI('N','N',N,N,N,MatsT(1.),FM,IK,JK,descMem,FM,IM,JM,descMem,
            MatsT(0.),full,1,1,descFull);
        else
#endif
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,MatsT(1.),K,2*N,M,2*N,MatsT(0.),full,N);

      } else {

#ifdef CQ_ENABLE_MPI
        if( isDist )
          Gemm_MPI('N','N',N,N,N,MatsT(1.),FM,IM,JM,descMem,FM,IK,JK,descMem,
            MatsT(0.),full,1,1,descFull);
        else
#endif
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,N,N,N,MatsT(1.),M,2*N,K,2*N,MatsT(0.),full,N);

      }

      return full;

    } else if( doAPB_AMB ) {


      int64_t NLoc = this->nSingleDim_, MLoc = NLoc;
#ifdef CQ_ENABLE_MPI
      if( isDist )
        std::tie(MLoc,NLoc) = this->fullMatGrid_->get_local_dims(NLoc,NLoc);
#endif


      MatsT* full = CQMemManager::get().malloc<MatsT>(MLoc*NLoc);

      std::fill_n(full,MLoc*NLoc,0.);


      if( isDist ) {

#ifdef CQ_ENABLE_MPI
        size_t N = this->nSingleDim_;
        int64_t ICTXT = this->fullMatGrid_->grid().context();

        int64_t NLocMem, MLocMem;
        std::tie(MLocMem,NLocMem) = this->fullMatGrid_->get_local_dims(N,N/2);

        auto descMem  = this->fullMatGrid_->descinit_noerror(N,N/2,MLocMem);
        auto descFull = this->fullMatGrid_->descinit_noerror(N,N  ,MLoc   );

        scalapackpp::wrappers::pgemr2d(N/2,N/2,this->fullMatrix_,1,1,descMem,
          full,(N/2)+1,1,descFull,ICTXT);
        scalapackpp::wrappers::pgemr2d(N/2,N/2,this->fullMatrix_,(N/2)+1,1,descMem,
          full,1,(N/2)+1,descFull,ICTXT);
#endif

      } else {

        size_t N = NLoc / 2;

        SetMat('N',N,N,MatsT(1.),this->fullMatrix_    ,MLoc,full + N     ,MLoc);
        SetMat('N',N,N,MatsT(1.),this->fullMatrix_ + N,MLoc,full + MLoc*N,MLoc);

      }

      return full;

    } else return this->fullMatrix_;

  };


  template <typename MatsT, typename IntsT>
  template <typename U>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::preConditioner(size_t nVec,
      U shift, SolverVectors<U> &V, SolverVectors<U> &AV) {


    auto& ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);

    size_t nOAVA = ss.nOA * ss.nVA;
    size_t nOBVB = ss.nOB * ss.nVB;
    size_t nOV   = ss.nO  * ss.nV;

    size_t N  = (ss.nC == 1) ? nOAVA  : nOV;
    size_t NV = (ss.nC == 1) ? ss.nVA : ss.nV;
    size_t NO = (ss.nC == 1) ? ss.nOA : ss.nO;

    std::function< U(double,double) > diag = 
      [&]( double eA, double eI ) {
        return (eA - eI) - shift;
      };


    if( doReduced ) 
      diag = [&]( double eA, double eI ) {
        return (eA - eI)*(eA - eI) - shift*shift;
      };


    double *eps = ss.eps1;

    size_t NS = this->nSingleDim_;
    size_t hNS = NS / 2;
    
    bool isRaw = V.underlyingType() == typeid(RawVectors<U>)
            and AV.underlyingType() == typeid(RawVectors<U>);

    if (isRaw and MPIRank(this->comm_) == 0) {
      U * V_ptr = tryGetRawVectorsPointer(V);
      U * AV_ptr = tryGetRawVectorsPointer(AV);
      for(auto iVec = 0ul; iVec < nVec; iVec++) {
        auto * AVk = AV_ptr + iVec * NS;
        auto * Vk  = V_ptr  + iVec * NS;
        for(auto k = 0ul; k < N; k++) {

          size_t i = k / NV;
          size_t a = (k % NV) + NO;

          AVk[k] = Vk[k] / diag(eps[a],eps[i]);

        } // X update


        if( not doReduced )
        for(auto k = 0ul; k < N; k++) {

          size_t i = k / NV;
          size_t a = (k % NV) + NO;

          AVk[k + hNS] = Vk[k + hNS] / diag(eps[a],eps[i]);

        } // Y update


        if( ss.nC == 1 ) {

          eps = ss.iCS ? eps : ss.eps2;

          N  = nOBVB;
          NV = ss.nVB;
          NO = ss.nOB;

          for(auto k = 0ul; k < N; k++) {

            size_t i = k / NV;
            size_t a = (k % NV) + NO;

            AVk[k + nOAVA] = Vk[k + nOAVA] / diag(eps[a],eps[i]);

          } // X update

          if( not doReduced )
          for(auto k = 0ul; k < N; k++) {

            size_t i = k / NV;
            size_t a = (k % NV) + NO;

            AVk[k + nOAVA + hNS] = Vk[k + nOAVA + hNS] / diag(eps[a],eps[i]);

          } // Y update

        } // Beta update


      } // loop over vectors

    } else if (not isRaw) {
      for(auto iVec = 0ul; iVec < nVec; iVec++) {
        for(auto k = 0ul; k < N; k++) {

          size_t i = k / NV;
          size_t a = (k % NV) + NO;

          AV.set(k, iVec, V.get(k, iVec) / diag(eps[a],eps[i]));

        } // X update


        if( not doReduced )
          for(auto k = 0ul; k < N; k++) {

            size_t i = k / NV;
            size_t a = (k % NV) + NO;

            AV.set(k + hNS, iVec, V.get(k + hNS, iVec) / diag(eps[a],eps[i]));

          } // Y update


        if( ss.nC == 1 ) {

          eps = ss.iCS ? eps : ss.eps2;

          N  = nOBVB;
          NV = ss.nVB;
          NO = ss.nOB;

          for(auto k = 0ul; k < N; k++) {

            size_t i = k / NV;
            size_t a = (k % NV) + NO;

            AV.set(k + nOAVA, iVec, V.get(k + nOAVA, iVec) / diag(eps[a],eps[i]));

          } // X update

          if( not doReduced )
            for(auto k = 0ul; k < N; k++) {

              size_t i = k / NV;
              size_t a = (k % NV) + NO;

              AV.set(k + nOAVA + hNS, iVec, V.get(k + nOAVA + hNS, iVec) / diag(eps[a],eps[i]));

            } // Y update

        } // Beta update

      } // loop over vectors
    }
  };


  template <typename MatsT, typename IntsT>
  std::pair<size_t,MatsT*> 
    PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formPropGrad(
      ResponseOperator op) {
    std::vector<IntsT*> opS;
    std::vector<MatsT*>      opT;

    Integrals<IntsT> &aoi    = *(this->ref_->aoints_);
    SingleSlater<MatsT, IntsT>& ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);

    bool needTrans = true;
    switch (op) {

      case LenElectricDipole: 
        opS = aoi.lenElectric->dipolePointers();
        break;

      case LenElectricQuadrupole: 
        opS = aoi.lenElectric->quadrupolePointers();
        break;

      case LenElectricOctupole: 
        opS = aoi.lenElectric->octupolePointers();
        break;

      case VelElectricDipole: 
        opS = aoi.velElectric->dipolePointers();
        break;

      case VelElectricQuadrupole: 
        opS = aoi.velElectric->quadrupolePointers();
        break;

      case VelElectricOctupole: 
        opS = aoi.velElectric->octupolePointers();
        break;

      case MagneticDipole: 
        opS = aoi.magnetic->dipolePointers();
        break;

      case MagneticQuadrupole: 
        opS = aoi.magnetic->quadrupolePointers();
        break;

      case Brillouin:
        needTrans = false;
        opT = ss.fockMO.size() > 1 ?
              std::vector<MatsT*>{ss.fockMO[0].pointer(), ss.fockMO[1].pointer()}
              : std::vector<MatsT*>{ss.fockMO[0].pointer()};
        break;

    }

    /*
    if(needTrans and (this->ref_->nC == 2 or not this->ref_->iCS)) 
      CErr("NO 2C or U");
      */

    //std::cerr << opS.size() << ", " << opT.size() << std::endl;

    size_t nVec = OperatorSize[op];
    size_t NB = ss.nAlphaOrbital();
    size_t NBC = ss.nC * NB;

    int nOAVA = ss.nOA * ss.nVA;
    int nOBVB = ss.nOB * ss.nVB;
    int nOV   = ss.nO  * ss.nV;

    int N  = (ss.nC == 1) ? nOAVA  : nOV;
    int NV = (ss.nC == 1) ? ss.nVA : ss.nV;
    int NO = (ss.nC == 1) ? ss.nOA : ss.nO;

    MatsT* grad  = CQMemManager::get().malloc<MatsT>(this->nSingleDim_*nVec);

    MatsT* SCR(nullptr);
    if( needTrans ) {
      SCR   = CQMemManager::get().malloc<MatsT>(NBC*NBC);
      opT.emplace_back(CQMemManager::get().malloc<MatsT>(NBC*NBC));
      if( ss.nC == 1 and not ss.iCS)
        opT.emplace_back(CQMemManager::get().malloc<MatsT>(NBC*NBC));
    }

    for(auto iVec = 0; iVec < nVec; iVec++) {
      MatsT* CMO = this->ref_->mo[0].pointer();
      MatsT* CMOB = (ss.nC == 1) ? this->ref_->mo[1].pointer() : CMO + NB;

      MatsT* V = grad + iVec*this->nSingleDim_;

      if( needTrans ) {

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NBC,NB,MatsT(1.),opS[iVec],NB ,CMO,NBC,MatsT(0.),SCR   ,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NBC,NBC,NB,MatsT(1.),CMO     ,NBC,SCR,NB ,MatsT(0.),opT[0],NBC);

        if( ss.nC == 1 and not ss.iCS ) {

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),opS[iVec],NB,CMOB,NB,MatsT(0.),SCR ,NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),CMOB     ,NB,SCR,NB,MatsT(0.),opT[1],NB);

        } else if( ss.nC == 2 ) {

          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NBC,NB,MatsT(1.),opS[iVec],NB ,CMOB,NBC,MatsT(0.),SCR  ,NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NBC,NBC,NB,MatsT(1.),CMOB    ,NBC,SCR,NB ,MatsT(1.),opT[0],NBC);

        }

      }

      MatsT* BOP = (ss.nC == 1 and not ss.iCS) ? opT[1] : opT[0];

      SetMat('N',NV,NO,MatsT(1.),opT[0] + NO,NBC,V,NV);
      if( ss.nC == 1 )
        SetMat('N',ss.nVB,ss.nOB,MatsT(1.),BOP + ss.nOB,NB,V + nOAVA,ss.nVB);

      //for(auto i = 0 , ai = 0; i < NO; i++)
      //for(auto a = NO        ; a < NB; a++, ai++) {

      //  V[ai]                               = opT[0][a + i*NB];
      //  V[ai+nOAVA]                         = opT[0][a + i*NB];

      //}

      if( doReduced ) 


        if( isHerOp(op) )
          for(auto i = 0 , ai = 0; i < NO; i++)
          for(auto a = NO        ; a < NB; a++, ai++) {

            V[ai]         = std::sqrt(0.5)*(V[ai]         + opT[0][i + a*NB]);
            V[ai + nOAVA] = std::sqrt(0.5)*(V[ai + nOAVA] + opT[0][i + a*NB]);

          }
        else
          for(auto i = 0 , ai = 0; i < NO; i++)
          for(auto a = NO        ; a < NB; a++, ai++) {

            V[ai]         = std::sqrt(0.5)*(V[ai]         - opT[0][i + a*NB]);
            V[ai + nOAVA] = std::sqrt(0.5)*(V[ai + nOAVA] - opT[0][i + a*NB]);

          }

      else {

        SetMat('N',NO,NV,MatsT(1.),opT[0] + NO*NBC,NBC,V + this->nSingleDim_/2,NO);
        IMatCopy('T',NO,NV,MatsT(1.),V + this->nSingleDim_/2,NO,NV);
        if( ss.nC == 1 ) {
          SetMat('N',ss.nOB,ss.nVB,MatsT(1.),BOP + ss.nOB*NB,NB,
            V + nOAVA + this->nSingleDim_/2,ss.nOB);
          IMatCopy('T',ss.nOB,ss.nVB,MatsT(1.),
            V + this->nSingleDim_/2 + nOAVA,ss.nOB,ss.nVB);
        }
        //for(auto i = 0 , ai = 0; i < NO; i++)
        //for(auto a = NO        ; a < NB; a++, ai++) {

        //  V[ai + this->nSingleDim_/2]         = opT[0][i + a*NB];
        //  V[ai + nOAVA + this->nSingleDim_/2] = opT[0][i + a*NB];

        //}

      }


    }

    if(SCR) CQMemManager::get().free(SCR);
    if( needTrans )
      for(auto &X : opT ) CQMemManager::get().free(X);


    // Transform to proper form form
    if( doAPB_AMB and not doReduced ) 
      blockTransform(this->nSingleDim_/2,nVec,std::sqrt(0.5),
        grad,this->nSingleDim_,grad+this->nSingleDim_/2,this->nSingleDim_);

#if 0
    if( this->savFile.exists() ) {

      if( op == LenElectricDipole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/LEN_ELEC_DIPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == LenElectricQuadrupole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/LEN_ELEC_QUADRUPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == LenElectricOctupole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/LEN_ELEC_OCTUPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == VelElectricDipole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/VEL_ELEC_DIPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == VelElectricQuadrupole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/VEL_ELEC_QUADRUPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == VelElectricOctupole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/VEL_ELEC_OCTUPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == MagneticDipole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/MAG_DIPOLE",grad,
          { nVec, this->nSingleDim_ } );

      if( op == MagneticQuadrupole )
        this->savFile.safeWriteData( "/RESP/PROPGRAD/MAG_QUADRUPOLE",grad,
          { nVec, this->nSingleDim_ } );
    }
#endif

    return {nVec, grad};

  };

  template <typename MatsT, typename IntsT>
  void PolarizationPropagator<SingleSlater<MatsT, IntsT>>::formRHS() {

    bool isDist = this->genSettings.isDist(); 
    bool isRoot = MPIRank(this->comm_) == 0;
    if( not (isDist and doReduced) ) ROOT_ONLY(this->comm_);

    isDist = isDist and doReduced;
    
    //std::cerr << "Top of formRHS\n";
    std::vector<ResponseOperator> ops = this->genSettings.bOps;

    size_t N      = this->nSingleDim_;
    
    MatsT* g = nullptr; size_t nProp;

    if( this->fdrSettings.nRHS == 0 ) 
      for(auto &op : ops) this->fdrSettings.nRHS += OperatorSize[op];

    if( this->fdrSettings.nRHS == 0 ) CErr("NO RHS GIVEN IN RESPONSE MODULE");

    MatsT* RHSa = nullptr;
    if( isRoot )
      RHSa = CQMemManager::get().malloc<MatsT>(this->fdrSettings.nRHS * N);

    // Space for distributed RHS and property gradient
    int64_t NLoc = N, NRHSLoc = this->fdrSettings.nRHS;


    MatsT* distRHS = nullptr;
    MatsT* distG   = nullptr;
#ifdef CQ_ENABLE_MPI
    scalapackpp::scalapack_desc DescRHS, DescG;
    if( isDist ) {

      std::tie(NLoc,NRHSLoc) = this->fullMatGrid_->get_local_dims(NLoc,NRHSLoc);

      if( NLoc and NRHSLoc )
        distRHS = CQMemManager::get().malloc<MatsT>(NLoc * NRHSLoc);

      DescRHS = 
        this->fullMatGrid_->descinit_noerror(N,this->fdrSettings.nRHS,NLoc);
    }
#endif

    MatsT* RHS = RHSa;

    for(auto &op : ops) {

      if( isRoot ) std::tie(nProp,g) = this->formPropGrad(op);    

      int64_t NPropLoc = nProp;
#ifdef CQ_ENABLE_MPI
      if( isDist ) {

        // Distribute property gradient
        MPIBCast(nProp,0,this->comm_);
        std::tie(NLoc,NPropLoc) = this->fullMatGrid_->get_local_dims(N,nProp);

        if( NLoc and NPropLoc )
          distG = CQMemManager::get().malloc<MatsT>(NLoc*NPropLoc);
        
        DescG = this->fullMatGrid_->descinit_noerror(N,nProp,NLoc);

        this->fullMatGrid_->scatter(N,nProp,g,N,distG,NLoc,0,0);

      } 
#endif

      // Copy G -> RHS
      if( isDist ) std::copy_n(distG,NPropLoc*NLoc,distRHS);
      else         std::copy_n(g,nProp*N,RHS);



      if( doReduced ) {

        if( not this->fullMatrix_ ) CErr("Reduced for non-full NYI");


        /*
        std::vector<std::tuple<size_t,T*,T*>> cList;
        cList.push_back(std::make_tuple(nProp,g,RHS));
        */


        std::vector<RESPONSE_CONTRACTION<MatsT>> cList(1);
        cList.back().nVec   = nProp;
        cList.back().N      = N;
        cList.back().X      = isDist ? distG : g;
        cList.back().AX     = isDist ? distRHS : RHS;
#ifdef CQ_ENABLE_MPI
        cList.back().DescX  = DescG;
        cList.back().DescAX = DescRHS;
#endif


        if( isHerOp(op) and this->fdrSettings.needP )
          formLinearTrans(cList, K);

        else if( isAntiHerOp(op) and this->fdrSettings.needQ ) 
          formLinearTrans(cList, M);


      } else if( doAPB_AMB )

        for(auto j = 0ul; j < nProp; j++)
          blas::swap(N/2,RHS + j*N,1,RHS + (N/2) + j*N,1);

      else if( incMet )

        SetMat('N',N/2,nProp,MatsT(-1.), RHS + (N/2), N, RHS + (N/2), N);
      

#ifdef CQ_ENABLE_MPI
      if( isDist ) {

        // Gather the RHS to root process
        this->fullMatGrid_->gather(N,nProp,RHS,N,distRHS,NLoc,0,0);

      }
#endif


      RHS += nProp*N;

      if( g ) CQMemManager::get().free(g);
      if( distG ) CQMemManager::get().free(distG);
 
    }

    if(this->fdrSettings.dampFactor == 0. and not this->fdrSettings.forceDamp)
      this->fdrResults.RHS = RHSa;
    else 
      this->dfdrResults.RHS = RHSa;


    if( distRHS ) CQMemManager::get().free(distRHS);

  }
  











  template <typename MatsT, typename IntsT>
  template <typename U, typename... Args>
  std::vector<TwoBodyContraction<U>> 
    PolarizationPropagator< SingleSlater<MatsT, IntsT> >::phTransitionVecMO2AO(
      MPI_Comm c, bool scatter, size_t nVec, size_t N,
      SingleSlater<MatsT,IntsT>& ss1, SingleSlater<MatsT,IntsT>& ss2, bool doExchange, Args... Vs) {

    static_assert(sizeof...(Vs) > 0 and sizeof...(Vs) < 3,
      "Vs must consist of 1 or 2 pointers");

    constexpr size_t nVs = sizeof...(Vs);

    std::array<U*,nVs> V_arr = { Vs... };


    bool isRoot = MPIRank(c) == 0;
    bool trans  = isRoot or not scatter;

    //SingleSlater<MatsT, IntsT>& sshold = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);

    const size_t NB   = ss1.nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = ss1.nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (ss1.nC == 2) ? ss1.nO : ss1.nOA;
    const size_t nOAVA = ss1.nOA * ss1.nVA;
    const size_t nOBVB = ss1.nOB * ss1.nVB;

    const size_t NB_AX = ss2.nAlphaOrbital();
    const size_t NB2_AX = NB_AX * NB_AX;
    const size_t NBC_AX = ss2.nC * NB_AX;
    const size_t NBC2_AX = NBC_AX * NBC_AX;

    U* MOT  = CQMemManager::get().malloc<U>(NBC2);
    U* SCR  = trans ? CQMemManager::get().malloc<U>(NBC2) : nullptr;


    std::vector<TwoBodyContraction<U>> cList;

    auto MOTRANS = [&]( MatsT* CMO, U* X ) {
      std::fill_n(SCR,NBC2,U(0.0));
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBC,NBC,NBC,U(1.0),CMO,NBC,X  ,NBC,U(0.0),SCR,NBC); 
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NBC,NBC,NBC,U(1.0),CMO,NBC,SCR,NBC,U(0.0),X  ,NBC);
      IMatCopy('C',NBC,NBC,U(1.),X,NBC,NBC);
    };

    size_t nAXMatPVec = 1;
    if( doExchange )
      nAXMatPVec += 2 * ss2.nC;
    size_t nXMatPVec = 2 * ss1.nC;
    size_t nSpacePVec = (nAXMatPVec * NB2_AX + nXMatPVec * NB2);
    size_t nAlloc = nVec * nSpacePVec;
    //std::cerr << "MO2AO " << nAlloc*sizeof(MatsT) / 1e9 << std::endl;

    U * first = CQMemManager::get().malloc<U>(nAlloc);

    for(size_t iVec = 0; iVec < nVec; iVec++) {

      U* VX_c = V_arr[0] + iVec * N;
      U* VY_c = nullptr; 
      if( nVs > 1 ) VY_c = V_arr[1] + iVec * N;

      MatsT* CMO = ss1.mo[0].pointer();

      // Perform MO -> AO transformation
      if( trans ) {

        // Zero out MOT
        std::fill_n(MOT,NBC2,0.);

        // Place V -> MOT and transform into AO basis (w Scatter)

        // Alpha for RHF and UHF, full for GHF
        for(size_t i = 0, ai = 0;  i < NO; i++) 
        for(size_t a = NO; a < NBC; a++, ai++) {
          MOT[i*NBC + a] = VX_c[ai];
          if( nVs > 1 ) MOT[a*NBC + i] = VY_c[ai];
        }
        // Transform Alpha (full) MOT -> AO basis
        MOTRANS(CMO, MOT);

      }

      // Scatter results if requested
      if( scatter ) MPIBCast(MOT,NBC*NBC,0,c);
 
      U *AOS  = nullptr, *AOZ  = nullptr, *AOY  = nullptr, *AOX  = nullptr;
      U *KAOS = nullptr, *KAOZ = nullptr, *KAOY = nullptr, *KAOX = nullptr;
      U *JAOS = nullptr;

      // Allocate AO storage

      AOS = first;
      AOZ = AOS + NB2;

      JAOS = AOZ + NB2;
      std::fill_n(JAOS,NB2_AX,0.);


      if( doExchange ) {
        KAOS = JAOS + NB2_AX;
        KAOZ = KAOS + NB2_AX;  // FIXME: not for SA

        std::fill_n(KAOS,NB2_AX,0.);
        std::fill_n(KAOZ,NB2_AX,0.);
      }
      
      if( ss1.nC == 2 ) {

        AOY = KAOZ + NB2_AX;
        AOX = AOY + NB2;

      }

      if( ss2.nC == 2 ) {
        if( doExchange ) {
          KAOY = AOX + NB2;
          KAOX = KAOY + NB2_AX;

          std::fill_n(KAOY,NB2_AX,0.);
          std::fill_n(KAOX,NB2_AX,0.);
        }
      }

      first += nSpacePVec;

      cList.push_back( { AOS, JAOS, false, COULOMB  } );
      if (doExchange){
        cList.push_back( { AOS, KAOS, false, EXCHANGE } );
        cList.push_back( { AOZ, KAOZ, false, EXCHANGE } );
        // FIXME: not for SA
  
        if( ss1.nC == 2 && ss2.nC == 2 ) {
          cList.push_back( { AOY, KAOY, false, EXCHANGE } );
          cList.push_back( { AOX, KAOX, false, EXCHANGE } );
        }
      }

      // TODO: what if ss1.nC != ss2.nC?
      if( ss1.nC == 1 && ss2.nC == 1 ) { // RHF / UHF case

        std::copy_n(MOT,NB2,AOS);
        std::copy_n(MOT,NB2,AOZ);

        // Perform transformation
        if( trans ) {

          std::fill_n(MOT,NB2,0.);

          double *eps = ss1.iCS ? ss1.eps1 : ss1.eps2;

          // Beta for RHF / UHF FIXME: not for spin adapted
          for(size_t i = 0, ai = nOAVA;  i < ss1.nOB; i++)
          for(size_t a = ss1.nOB; a < NB; a++, ai++) {
            MOT[i*NB + a] = VX_c[ai];
            if( nVs > 1 ) MOT[a*NB + i] = VY_c[ai];
          }

          CMO = ss1.iCS ? ss1.mo[0].pointer() : ss1.mo[1].pointer();
          // Transform BETA (full) MOT -> AO basis
          MOTRANS(CMO,MOT);
        }


        // Scatter results if requested
        if( scatter ) MPIBCast(MOT,NBC*NBC,0,c);

        MatAdd('N','N',NB,NB,U(1.),AOS,NB,U(1.) ,MOT,NB,AOS,NB);
        MatAdd('N','N',NB,NB,U(1.),AOZ,NB,U(-1.),MOT,NB,AOZ,NB);


      } else { // GHF

        SpinScatter(NB,MOT,NBC,AOS,NB,AOZ,NB,AOY,NB,AOX,NB);

      }

    }


    CQMemManager::get().free(MOT);
    if( SCR ) CQMemManager::get().free(SCR);


    return cList;

  };



  template <typename MatsT, typename IntsT>
  template <typename U, typename... Args>
  void PolarizationPropagator< SingleSlater<MatsT, IntsT> >::phTransitionVecAO2MO(
    size_t nVec, size_t N, std::vector<TwoBodyContraction<U>> &cList, SingleSlater<MatsT, IntsT>& ss,bool doExchange,
    Args... HVs) {

    constexpr size_t nHVs = sizeof...(HVs);

    std::array<U*,nHVs> HV_arr = { HVs... };


   // SingleSlater<MatsT, IntsT> &ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);
    const size_t NB   = ss.nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = ss.nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (ss.nC == 2) ? ss.nO : ss.nOA;
    const size_t nOAVA = ss.nOA * ss.nVA;
    const size_t nOBVB = ss.nOB * ss.nVB;

    U* MOT  = CQMemManager::get().malloc<U>(NBC2);
    U* SCR  = CQMemManager::get().malloc<U>(NBC2);

    const size_t iOff = (ss.nC == 2) ? 5 : 3;

    auto MOTRANS = [&]( MatsT* CMO, U* X ) {
      std::fill_n(SCR,NBC2,U(0.0));
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NBC,NBC,NBC,U(1.0),CMO,NBC,X  ,NBC,U(0.0),SCR,NBC); 
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::ConjTrans,NBC,NBC,NBC,U(1.0),CMO,NBC,SCR,NBC,U(0.0),X  ,NBC);
      IMatCopy('C',NBC,NBC,U(1.),X,NBC,NBC);
    };


    for(size_t iVec = 0; iVec < nVec; iVec++) {

      U* HVX_c = HV_arr[0] + iVec * N;
      U* HVY_c = nullptr;
      if( nHVs > 1 ) HVY_c = HV_arr[1] + iVec * N;

      // Get the index of first G[T(iVec)]

      size_t indx = iOff * iVec; // FIXME: be careful for SA
      U* J_S=nullptr, *K_S=nullptr, *K_Z=nullptr, *K_X=nullptr, *K_Y=nullptr;

      if (doExchange){
        J_S = cList[indx  ].AX;
        K_S = cList[indx+1].AX;
        K_Z = cList[indx+2].AX;
  
        if( ss.nC == 2 ) {
          K_Y = cList[indx+3].AX;
          K_X = cList[indx+4].AX;
        }
        // Form G(S) in K(S)
        MatAdd('N','N',NB,NB,U(2.),J_S,NB,U(-1.),K_S,NB,K_S,NB);
  
        // Negate Ks for Z,Y and X
        blas::scal(NB2,U(-1.),K_Z,1);
      
        if( ss.nC == 2 ) {
          blas::scal(NB2,U(-1.),K_Y,1);
          blas::scal(NB2,U(-1.),K_X,1);
        }
      }
      else{
        J_S = cList[iVec].AX;
      }
      // Transform G[T] into MO basis and extract into HV
            
      // Alpha for RHF/ UHF, full for GHF
      if (doExchange){
        if( ss.nC == 1 )
          MatAdd('N','N',NB,NB,U(0.5),K_S,NB,U(0.5),K_Z,NB,MOT,NB);
        else
          SpinGather(NB,MOT,NBC,K_S,NB,K_Z,NB,K_Y,NB,K_X,NB);
      }
      else{
        SetMat('N',NB,NB,U(1.),J_S,NB,MOT,NB);  
      }
      MatsT* CMO = ss.mo[0].pointer();

      // Transform -> AO basis
      MOTRANS(CMO,MOT);

      for(size_t i = 0, ai = 0;  i < NO; i++) 
      for(size_t a = NO; a < NBC; a++, ai++) {
        HVX_c[ai] += MOT[i*NBC + a]; 
        if( nHVs > 1 ) HVY_c[ai] += MOT[a*NBC + i];
      }

      // Do Beta for RHF / UHF (FIXME: not for SA)
      if( ss.nC == 1 ) {
        if (doExchange){
          MatAdd('N','N',NB,NB,U(0.5),K_S,NB,U(-0.5),K_Z,NB,MOT,NB);
        }
        else{
          SetMat('N',NB,NB,U(1.),J_S,NB,MOT,NB);
        }
        CMO = ss.iCS ? ss.mo[0].pointer() : ss.mo[1].pointer();
        // Transform -> AO basis
        MOTRANS(CMO,MOT);
        
        for(size_t i = 0, ai = nOAVA;  i < ss.nOB; i++)
        for(size_t a = ss.nOB; a < NB; a++, ai++) {
          HVX_c[ai] += MOT[i*NB + a];
           
          if( nHVs > 1 ) HVY_c[ai] += MOT[a*NB + i]; 
        }

      }

    }

    CQMemManager::get().free(MOT,SCR);

  }; 




  template <typename MatsT, typename IntsT>
  template <typename U>
  void PolarizationPropagator< SingleSlater<MatsT, IntsT> >::phEpsilonScale(bool doInc, 
    bool doInv, size_t nVec, size_t N,SingleSlater<MatsT,IntsT>& ss, U* V, U* HV) {
   // SingleSlater<MatsT,IntsT>& ss = dynamic_cast<SingleSlater<MatsT, IntsT>&>(*this->ref_);

    const size_t NB   = ss.nAlphaOrbital();
    const size_t NB2  = NB * NB;
    const size_t NBC  = ss.nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (ss.nC == 2) ? ss.nO : ss.nOA;
    const size_t NV    = (ss.nC == 2) ? ss.nV : ss.nVA;
    const size_t nOAVA = ss.nOA * ss.nVA;
    const size_t nOBVB = ss.nOB * ss.nVB;


    auto epsilonScale = doInv ? []( double x ) { return 1./x; } :
                                []( double x ) { return x;    };

    auto increment    = doInc ? []( U x, U y ) { return x + y; } :
                                []( U x, U y ) { return y;  };


    if( orbHessCtl.useEigenEnergies )
    for(size_t iVec = 0; iVec < nVec; iVec++) {

      U* V_c  = V  + iVec * N;
      U* HV_c = HV + iVec * N;

      double *eps = ss.eps1; 

      for(size_t i = 0, ai = 0;  i < NO; i++) 
      for(size_t a = NO; a < NBC; a++, ai++){ 
        HV_c[ai] = increment(HV_c[ai],
                             epsilonScale(eps[a] - eps[i]) * V_c[ai]);
      }

      if( ss.nC == 1 ) { // RHF / UHF case

        eps = ss.iCS ? ss.eps1 : ss.eps2;

        for(size_t i = 0, ai = nOAVA;  i < ss.nOB; i++) 
        for(size_t a = ss.nOB; a < NB; a++, ai++) 
          HV_c[ai] = increment(HV_c[ai],
                               epsilonScale(eps[a] - eps[i]) * V_c[ai]);

      }

    }
    
    else {

      dcomplex* FCMPLX = nullptr, *FBCMPLX = nullptr;
      if( std::is_same<U,dcomplex>::value and std::is_same<MatsT,double>::value) {

        FCMPLX = CQMemManager::get().malloc<dcomplex>(NBC2);
        std::copy_n(ss.fockMO[0].pointer(),NBC2,FCMPLX);

        if(ss.nC == 1 and not ss.iCS) {
          FBCMPLX = CQMemManager::get().malloc<dcomplex>(NBC2);
          std::copy_n(ss.fockMO[1].pointer(),NBC2,FBCMPLX);
        }
      }

      U* Foo = bool(FCMPLX) ? 
        reinterpret_cast<U*>(FCMPLX) : reinterpret_cast<U*>(ss.fockMO[0].pointer());
      U* Fvv = bool(FCMPLX) ? 
        reinterpret_cast<U*>(FCMPLX + NO*(NBC+1)) : 
        reinterpret_cast<U*>(ss.fockMO[0].pointer() + NO*(NBC+1));

      U* Foob = nullptr, *Fvvb = nullptr;

      if(ss.nC == 1 and not ss.iCS) {
        Foob = bool(FBCMPLX) ? 
          reinterpret_cast<U*>(FBCMPLX) : reinterpret_cast<U*>(ss.fockMO[1].pointer());
        Fvvb = bool(FBCMPLX) ? 
          reinterpret_cast<U*>(FBCMPLX + ss.nOB*(NBC+1)) : 
          reinterpret_cast<U*>(ss.fockMO[1].pointer() + ss.nOB*(NBC+1));
      }


      U fact = doInc ? 1. : 0.;

      for(size_t iVec = 0; iVec < nVec; iVec++) {

        U* V_c  = V  + iVec * N;
        U* HV_c = HV + iVec * N;


        // HV(a,i) = \sum_b F(a,b) V(b,i)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NV,NO,NV,U(1.) ,Fvv,NBC,V_c,NV,fact ,HV_c,NV);
        // HV(a,i) -= \sum_j V(a,j) F(i,j)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NV,NO,NO,U(-1.),V_c,NV,Foo,NBC,U(1.),HV_c,NV);


        U* V_cb  = V_c  + nOAVA;
        U* HV_cb = HV_c + nOAVA;

        if( ss.nC == 1 ) {
          if( ss.iCS ) {

            // HV(a,i) = \sum_b F(a,b) V(b,i)
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NV,NO,NV,U(1.) ,Fvv,NB,V_cb,NV,fact ,HV_cb,NV);
            // HV(a,i) -= \sum_j V(a,j) F(i,j)
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NV,NO,NO,U(-1.),V_cb,NV,Foo,NB,U(1.),HV_cb,NV);

          } else {

            // HV(a,i) = \sum_b F(a,b) V(b,i)
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,ss.nVB,ss.nOB,ss.nVB,U(1.),Fvvb,NB,
              V_cb,ss.nVB,fact ,HV_cb,ss.nVB);
            // HV(a,i) -= \sum_j V(a,j) F(i,j)
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,ss.nVB,ss.nOB,ss.nOB,U(-1.),V_cb,ss.nVB,
              Foob,NB,U(1.),HV_cb,ss.nVB);

          }
        }

      }

      if( FCMPLX ) CQMemManager::get().free(FCMPLX);
      if( FBCMPLX ) CQMemManager::get().free(FBCMPLX);

    }


  };


  template <typename MatsT, typename IntsT>
  void PolarizationPropagator< SingleSlater<MatsT, IntsT> >::resGuess(
      size_t nGuess, SolverVectors<MatsT> &G, size_t LDG) {

    SingleSlater<MatsT, IntsT>& ss = dynamic_cast<SingleSlater<MatsT,IntsT>&>(*this->ref_);
    size_t N = getNSingleDim(true);

    // Initialize an NOV index array
    std::vector<size_t> indx(N,0);
    std::iota(indx.begin(),indx.end(),0);

    size_t NO1 = (ss.nC == 1) ? ss.nOA : ss.nO;
    size_t NV1 = (ss.nC == 1) ? ss.nVA : ss.nV;
    size_t NO2 = ss.nOB;
    size_t NV2 = ss.nVB;

    size_t NOV1 = NO1 * NV1;
    size_t NOV2 = NO2 * NV2;

    double *EPS1 = ss.eps1;
    double *EPS2 = (ss.nC == 1 and not ss.iCS) ? ss.eps2 : ss.eps1;


    auto conv_indx = [&](size_t INDX) -> double {


      // If INDX > NOV1 then BETA else ALPHA/2C
      size_t nv   = (INDX >= NOV1) ? NV2 : NV1;
      size_t no   = (INDX >= NOV1) ? NO2 : NO1;
      size_t indx = (INDX >= NOV1) ? (INDX - NOV1) : INDX;

      double *eps = (INDX >= NOV1) ? EPS2 : EPS1;



      size_t i = indx / nv;
      size_t a = (indx % nv) + no;

      return eps[a] - eps[i];

    };

    std::function< bool(size_t, size_t) > compare = 
      [&](size_t INDX1, size_t INDX2) -> bool {
       
        double delta1 = conv_indx(INDX1); 
        double delta2 = conv_indx(INDX2); 

        return std::abs(delta1 - this->resSettings.deMin) < 
               std::abs(delta2 - this->resSettings.deMin);

      };

    // Sort
    std::sort(indx.begin(),indx.end(),compare);

    // Set guess
    G.clear();

    for(size_t iG = 0; iG < nGuess; iG++) G.set(indx[iG], iG, 1.);

    /*
    // Confirm
    std::vector<double> delta(N); 
    for(auto k=0; k<N; k++) delta[k] = conv_indx(k);

    std::sort(delta.begin(),delta.end());

    //for(auto k = 0; k < nGuess; k++)
    //  std::cerr << delta[k] << ", " << conv_indx(indx[k]) << std::endl;

    */
  }


  /**
   * Get the diagonal double bar integral (ai||ai)
   * - XSCR needs to be 2*nC*NB2
   * - AXSCR needs to be (2*nC+1)*NB2
   * - Both scratch must not contain NaNs
   **/
  template <typename MatsT, typename IntsT>
  MatsT PolarizationPropagator< SingleSlater<MatsT, IntsT> >::getGDiag(
    size_t i, size_t a, bool beta, SingleSlater<MatsT,IntsT>& ss,
    MatsT* XSCR, MatsT* AXSCR) 
  {

    MatsT* mo = beta && !ss.iCS ? ss.mo[1].pointer() : ss.mo[0].pointer();
    size_t NB = ss.nAlphaOrbital();
    size_t NB2 = NB*NB;
    size_t NBC = NB*ss.nC;
    size_t NBC2 = NBC*NBC;

    size_t NO = ss.nC == 2 ?  ss.nO : beta ?  ss.nOB : ss.nOA;

    // D_{m,n} := C_{a,m} \otimes C_{i,n}^T
    // Form in AXSCR, then scatter to XSCR
    // TODO: Move this to blas::ger when this branch is updated
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NBC,NBC,1,MatsT(1.),mo+a*NBC,NBC,mo+i*NBC,NBC,
         MatsT(0.),AXSCR,NBC);

    // Put D into pauli spinor form
    MatsT *AOS, *AOZ, *AOY, *AOX;
    AOS = XSCR;
    AOZ = XSCR + NB2;
    if( ss.nC == 2 ) {
      AOY = XSCR + 2*NB2;
      AOZ = XSCR + 3*NB2;
    }

    if( ss.nC == 1 ) {
      SetMat('N', NB, NB, MatsT(1.), AXSCR, NB, AOS, NB);
      MatsT factor = beta ? MatsT(-1.) : MatsT(1.);
      SetMat('N', NB, NB, factor, AXSCR, NB, AOZ, NB);
    }
    else {
      SpinScatter(NB,AXSCR,NBC,AOS,NB,AOZ,NB,AOY,NB,AOX,NB);
    }

    // G[D] = Contract all ingredients with D
    MatsT *JS, *KS, *KZ, *KY, *KX;

    JS = AXSCR;
    KS = AXSCR + NB2;
    KZ = AXSCR + 2*NB2;
    if( ss.nC == 2 ) {
      KY = AXSCR + 3*NB2;
      KX = AXSCR + 4*NB2;
    }

    std::vector<TwoBodyContraction<MatsT>> cList;
    cList.push_back( { AOS, JS, false, COULOMB  } );
    cList.push_back( { AOS, KS, false, EXCHANGE } );
    cList.push_back( { AOZ, KZ, false, EXCHANGE } );
    if( ss.nC == 2 ) {
      cList.push_back( { AOY, KY, false, EXCHANGE } );
      cList.push_back( { AOX, KX, false, EXCHANGE } );
    }

    ss.TPI->twoBodyContract(MPI_COMM_WORLD,cList);

    // Put G into spin blocked form in XSCR
    MatAdd('N', 'N', NB, NB, MatsT(2.), JS, NB, MatsT(-1.), KS, NB, KS, NB);
    blas::scal(NB2,MatsT(-1.),KZ,1);
    if( ss.nC == 2 ) {
      blas::scal(NB2,MatsT(-1.),KY,1);
      blas::scal(NB2,MatsT(-1.),KX,1);
    }

    if( ss.nC == 1 ) {
      MatsT factor = beta ? -0.5 : 0.5;
      MatAdd('N','N',NB,NB,MatsT(0.5),KS,NB,factor,KZ,NB,XSCR,NB);
    }
    else {
      SpinGather(NB,XSCR,NBC,KS,NB,KZ,NB,KY,NB,KX,NB);
    }

    // Do the final AO->MO contraction
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,1,NBC,NBC,MatsT(1.),mo+a*NBC,NBC,XSCR,NBC,MatsT(0.),AXSCR,NBC);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,1,1,NBC,MatsT(1.),AXSCR,NBC,mo+i*NBC,NBC,MatsT(0.),XSCR,NBC);

    return *XSCR;

  }
} // namespace ChronusQ


