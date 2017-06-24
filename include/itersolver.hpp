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

#include <cqlinalg.hpp>
#include <cerr.hpp>

#include <util/mpi.hpp>
#include <util/math.hpp>

#include <itersolver/solvervectors.hpp>

namespace ChronusQ {

  template <typename _F>
  class IterSolver {

  protected:

    typedef std::function< std::shared_ptr<SolverVectors<_F>>(size_t) >            VecsGen_t;
    typedef std::function< void(size_t,SolverVectors<_F>&,SolverVectors<_F>&) >    LinearTrans_t;
    typedef std::function< void(size_t,_F,SolverVectors<_F>&,SolverVectors<_F>&) > Shift_t;

    MPI_Comm     comm_;

    size_t N_;       ///< Problem dimension
    size_t mSS_;     ///< Max subspace dimension
    size_t maxMacroIter_; ///< Maximum number of Macro iterations
    size_t maxMicroIter_; ///< Maximum number of Micro iterations


    double convCrit_; ///< Convergence criteria
    bool converged_ = false; ///< Iteration has converged

    VecsGen_t vecGen_; ///< Vector generator

    LinearTrans_t linearTrans_;    ///< AX Product
    LinearTrans_t preCondNoShift_; ///< Unshifted preconditioner

    Shift_t shiftVec_;       ///< (A - sB)X given AX
    Shift_t preCondWShift_;  ///< Shifted preconditioner

    inline std::shared_ptr<SolverVectors<_F>> rawVecsGenerator_(size_t nVec) {

      return std::make_shared<RawVectors<_F>>(comm_, this->N_, nVec);

    };

    inline void defaultShiftVec_(size_t nVec, _F shift, SolverVectors<_F> &V, SolverVectors<_F> &AV) {

      AV.axpy(0, nVec, shift, V, 0);

    };

    inline void defaultPreCondShift_(size_t nVec, _F shift, SolverVectors<_F> &V, SolverVectors<_F> &AV) {

      preCondNoShift_(nVec,V,AV);
      if(std::abs(shift) > 1e-15) shiftVec_(nVec,-1./shift,V,AV);

    };

  public:

    // Constructor (unshifted preconditioner)
    IterSolver(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
        comm_(c), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), vecGen_(vecGen),
        linearTrans_(linearTrans), shiftVec_(shiftVec),
        preCondNoShift_(preNoShift) {

      using namespace std::placeholders;
      if( not vecGen_ ) {
        vecGen_ =
            std::bind(&IterSolver<_F>::rawVecsGenerator_,this,_1);
      }

      if( not shiftVec_ ) {
        shiftVec_ =
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4);
      }

      preCondWShift_ =
        std::bind(&IterSolver<_F>::defaultPreCondShift_,this,_1,_2,_3,_4);

    }

    // Constructor (shifted preconditioner)
    IterSolver(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
        comm_(c), N_(N), mSS_(MSS), maxMacroIter_(MAXMACRO),
        maxMicroIter_(MAXMICRO), convCrit_(conv), vecGen_(vecGen),
        linearTrans_(linearTrans), shiftVec_(shiftVec),
        preCondWShift_(preShift) {

      using namespace std::placeholders;
      if( not vecGen_ ) {
        vecGen_ =
            std::bind(&IterSolver<_F>::rawVecsGenerator_,this,_1);
      }

      if( not shiftVec_ ) {
        shiftVec_ =
          std::bind(&IterSolver<_F>::defaultShiftVec_,this,_1,_2,_3,_4);
      }

    }

    bool hasConverged() const { return converged_; }

  };




  template <typename _F>
  class IterLinearSolver : public IterSolver<_F> {

  protected:

    size_t nRHS_;
    std::shared_ptr<SolverVectors<_F>> RHS_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> SOL_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> V_   = nullptr;
    std::shared_ptr<SolverVectors<_F>> AV_  = nullptr;
    std::shared_ptr<SolverVectors<_F>> RES_ = nullptr;

    std::vector<_F> shifts_;

    std::vector<double>              rhsNorm_;
    std::vector<std::vector<double>> resNorm_;

  public:


    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    size_t shiftBS = 1;
    size_t rhsBS   = 1;

    IterLinearSolver(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preNoShift,
                    vecGen,shiftVec){ }

    IterLinearSolver(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,preShift,
                    vecGen,shiftVec){ }


    ~IterLinearSolver() {};




    template <typename T>
    void setRHS(size_t nRHS, T*RHS, size_t LDRHS);

    template <typename T>
    void setShifts(size_t nShift, T* shifts);


    std::shared_ptr<SolverVectors<_F>> getSol() const {
      return SOL_;
    };




    virtual void alloc() {

//      size_t nRHSnSN       = nRHS_ * shifts_.size() * this->N_;
//      size_t nRHSnSN_batch = rhsBS * shiftBS * this->N_;

      SOL_ = this->vecGen_(nRHS_ * shifts_.size());
      V_   = this->vecGen_(rhsBS * shiftBS * this->mSS_);
      AV_  = this->vecGen_(rhsBS * shiftBS * this->mSS_);
      RES_  = this->vecGen_(rhsBS * shiftBS);

    }



    void run();

    virtual void runBatch(size_t nRHS, size_t nShift,
                          std::shared_ptr<SolverVectors<_F>> RHS, _F *shifts,
                          std::shared_ptr<SolverVectors<_F>> SOL, double *RHSNorm);

  };




  struct IterDiagConvStatus {
    bool residual = false;
    bool eigenValue = false;
    bool eigenVector = false;
    double residual_norm = 0.0;
    double diff_eigen_value = 0.0;
    double diff_eigen_vector = 0.0;
    int prev_index = -1;

    bool hasConverged(bool checkEVector = true, bool checkEValue = true, bool checkResidual = true) const {
      return (checkResidual ? residual : true)
          and (checkEValue ? eigenValue : true)
          and (checkEVector ? eigenVector : true);
    }

    void clear() {
      residual = false;
      eigenValue = false;
      eigenVector = false;
      residual_norm = 0.0;
      diff_eigen_value = 0.0;
      diff_eigen_vector = 0.0;
      prev_index = -1;
    }
  };



  template <typename _F>
  class IterDiagonalizer : public IterSolver<_F> {


  protected:

    size_t nRoots_;
    size_t nGuess_;
    std::shared_ptr<SolverVectors<_F>> VR_   = nullptr;
    _F     *VL_   = nullptr;
    std::shared_ptr<SolverVectors<_F>> AVR_  = nullptr;
    _F     *AVL_  = nullptr;
    _F     *RESR_ = nullptr;
    _F     *RESL_ = nullptr;


    dcomplex * eigVal_ = nullptr;

    void RayleighRitz();



  public:

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterSolver<_F>::Shift_t;

    IterDiagonalizer(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      size_t nR,
      size_t nG,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
      IterSolver<_F>(c,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                     preNoShift,vecGen,shiftVec),
        nRoots_(nR), nGuess_(nG){ }

    IterDiagonalizer(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      const size_t MAXMACRO,
      const size_t MAXMICRO,
      double conv,
      size_t nR,
      size_t nG,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterSolver<_F>(c,N,MSS,MAXMACRO,MAXMICRO,conv,linearTrans,
                    preShift,vecGen,shiftVec),
       nRoots_(nR), nGuess_(nG){ }

    ~IterDiagonalizer() {

      if(VL_)   CQMemManager::get().free(VL_);
      if(AVL_)  CQMemManager::get().free(AVL_);
      if(RESR_) CQMemManager::get().free(RESR_);
      if(RESL_) CQMemManager::get().free(RESL_);
      if(eigVal_) CQMemManager::get().free(eigVal_);

    };




    virtual void alloc() {
      
      VR_     = this->vecGen_(this->nRoots_);

      // NO MPI
      // ROOT_ONLY(this->comm_);

//      size_t NNR = this->N_ * this->mSS_;

      eigVal_ = CQMemManager::get().malloc<dcomplex>(this->nRoots_);
      //AVR_    = CQMemManager::get().malloc<_F>(NNR);
      //RESR_   = CQMemManager::get().malloc<_F>(NNR);

    }


    dcomplex* eigVal() { return eigVal_; }
    const dcomplex* eigVal() const { return eigVal_; }
    std::shared_ptr<SolverVectors<_F>> VR() { return VR_; }
    std::shared_ptr<const SolverVectors<_F>> VR() const { return VR_; }


    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)>) = 0;

    void run();

    virtual bool runMicro() = 0;
    virtual void restart() = 0;

  };








  // GPLHR
  template <typename _F>
  class GPLHR : public IterDiagonalizer<_F> {

    double *RelRes = nullptr;
    std::shared_ptr<SolverVectors<_F>> Guess = nullptr;


  protected:


    bool useAdaptiveSigma() const { return false; }


    void getResidualNorms(size_t N, size_t nR, const SolverVectors<_F> &WMAT,
                          double *RelRes, dcomplex *LAMBDA, double nrmA);

    void getTriU(size_t N, _F *TRIUA, size_t LDTRIUA, _F *TRIUB,
      size_t LDTRIUB, _F *MA, size_t LDMA, _F *MB, size_t LDMB);

    bool checkConv(size_t nR, double *RelRes);




    void halfProj(size_t nV, size_t nS, const SolverVectors<_F> &V,
                  SolverVectors<_F> &S, _F *smSCR, size_t LDSCR);

    /**
     *  \brief Overload of halfProj where nV == nS
     */
    inline void halfProj(size_t nV, const SolverVectors<_F> &V,
                         SolverVectors<_F> &S, _F *smSCR, size_t LDSCR) {

      halfProj(nV,nV,V,S,smSCR,LDSCR);

    }

    void halfProj2(size_t nV, size_t nS, const SolverVectors<_F> &V, const SolverVectors<_F> &AV,
                   SolverVectors<_F> &S, SolverVectors<_F> &AS, _F *smSCR, size_t LDSCR);

    /**
     *  \brief Overload of halfProj2 where nV == nS
     */
    inline void halfProj2(size_t nV, const SolverVectors<_F> &V, const SolverVectors<_F> &AV,
                          SolverVectors<_F> &S, SolverVectors<_F> &AS, _F *smSCR, size_t LDSCR) {

      halfProj2(nV,nV,V,AV,S,AS,smSCR,LDSCR);

    }


    void newSMatrix(size_t nR, const SolverVectors<_F> &V, const SolverVectors<_F> &Q,
                    SolverVectors<_F> &S, _F *smSCR, size_t LDSCR);


  public:


    size_t m = 1;
    _F sigma = 0.;
    _F hardLim = -std::numeric_limits<double>::infinity();

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    GPLHR(
      MPI_Comm c,
      const size_t N,
      const size_t MAXITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t(),
      const size_t m = 1 ) :
      IterDiagonalizer<_F>(c,N,2 + m*nR,1,MAXITER,conv,nR,nR,
                           linearTrans,preNoShift,vecGen,shiftVec), m(m) { }

    GPLHR(
      MPI_Comm c,
      const size_t N,
      const size_t MAXITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t(),
      const size_t m = 1) :
      IterDiagonalizer<_F>(c,N,2 + m*nR,1,MAXITER,conv,nR,nR,
                           linearTrans,preShift,vecGen,shiftVec), m(m) { }

    ~GPLHR() {

      if( RelRes ) CQMemManager::get().free(RelRes);

    }

    void setM(size_t _m) {

      m = _m;
      this->mSS_ = (3 + m)*this->nRoots_;

    }

    void alloc() override {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      // ROOT_ONLY(this->comm_);


      // Allocate GPLHR specific Memory

      this->RelRes = CQMemManager::get().malloc<double>(this->nRoots_);

    }

    bool runMicro() override;

    void restart() override;

    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)> func) override {

      if( nGuess != this->nRoots_ )
        CErr("GPLHR Requires nGuess = nRoots",std::cout);

      Guess = this->vecGen_(nGuess);

      // NO MPI
      // ROOT_ONLY(this->comm_);

      func(nGuess, *Guess, this->N_);

    }

  }; // GPLHR





  // Davidson
  template <typename _F>
  class Davidson : public IterDiagonalizer<_F> {

    bool DoLeftEigVec   = false;
//    bool sortByDistance = false;
    size_t GramSchmidt_NRe = 1;
    double GramSchmidt_eps = 1e-12;
    bool DoHerm = false;

    // Convergence options
    bool checkEigenValueConv = true;
    bool checkEigenVectorConv = true;
    bool checkResidueConv = false;
    bool convOnGramSchmidt = true;
    double eigenValueCrit = 1e-7;
    double eigenVectorCrit = 1e-5;

    // Energy specific options
    std::vector<std::pair<double, size_t>> energyRefs;
    bool AbsoluteES = false; // use absolute energy threshold or relative energy threshold
//    size_t nHighERoots  = 0;
    size_t  nLowERoots  = 0;
//    double EnergyRef    = 0.;
//    bool adaptiveERef   = false; // use loweset eigenvalues at current iteration

    double   *RelRes  = nullptr;
    std::shared_ptr<SolverVectors<_F>> Guess = nullptr;
    std::shared_ptr<SolverVectors<_F>> vecs  = nullptr;
    std::shared_ptr<SolverVectors<_F>> sigmaVecs = nullptr;
    std::shared_ptr<SolverVectors<_F>> R   = nullptr, S = nullptr; // Scratch space for residue and perturbed vector
    dcomplex *EigForT = nullptr; // Eigenvalues can be used for preconditioner

    // creating double(dcomplex) operator to accomodate the
    // complex eigenvalue to real arithmetric
    _F dcomplexTo_F(dcomplex &a) { return * reinterpret_cast<_F*>(&a);}

  protected:

  public:

    size_t m =  50;
    size_t whenSc = 2; // Iteration at which to scale back the number of vectors added, or 0 if no change.
    size_t kG = 3;

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterDiagonalizer<_F>::LinearTrans_t;
    using Shift_t       = typename IterDiagonalizer<_F>::Shift_t;

    Davidson(
      MPI_Comm c,
      const size_t N,
      const size_t MAXMACROITER,
      const size_t MAXMICROITER,
      double conv,
      size_t nR,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t(),
      const size_t m = 50,
      const size_t kG = 3):
      IterDiagonalizer<_F>(c,N,m*nR,MAXMACROITER,MAXMICROITER,conv,nR,nR*kG,
                           linearTrans,preNoShift,vecGen,shiftVec), m(m), kG(kG) {
      eigenVectorCrit = conv;
      eigenValueCrit = 1e-2 * conv;
    }

    ~Davidson() {

      if( RelRes  ) CQMemManager::get().free(RelRes);

    }

    void setM(size_t _m) {
      m = _m;
      this->mSS_ = m*this->nRoots_;
    }

    void setkG(size_t _kG) {
      kG = _kG;
      this->nGuess_ = kG*this->nRoots_;
    }

    void setHerm(bool sSym) {
      DoHerm = sSym;
    }

    void setWhenSc(size_t _WhenSc) { whenSc = _WhenSc;}

    void setGramSchmidtRepeat(size_t nRe) { GramSchmidt_NRe = nRe;}

    void setGramSchmidtEps(double eps) { GramSchmidt_eps = eps;}

    void setConvOnGramSchmidt(bool flag) { convOnGramSchmidt = flag; }

    void setEigenVectorConvCheck(bool flag) { checkEigenVectorConv = flag;}

    void setEigenValueConvCheck(bool flag) { checkEigenValueConv = flag;}

    void setResidueConvCheck(bool flag) { checkResidueConv = flag;}

    void setEigenVectorConvCriteria(double crit) { eigenVectorCrit = crit;}

    void setEigenValueConvCriteria(double crit) { eigenValueCrit = crit;}

    void setEigForT(dcomplex * _Eig) {

      if( CQMemManager::get().getSize(_Eig) < this->nGuess_ )
        CErr("Davison EigForT requires a memory block with size at least nGuess ",std::cout);

      EigForT = _Eig;
    }

    std::shared_ptr<SolverVectors<_F>> getGuessScratch() const {
      return Guess;
    }
    std::shared_ptr<SolverVectors<_F>> getSubspaceScratch() const {
      return vecs;
    }
    std::shared_ptr<SolverVectors<_F>> getSigmaVecScratch() const {
      return sigmaVecs;
    }
    std::shared_ptr<SolverVectors<_F>> getScratchR() const {
      return R;
    }
    std::shared_ptr<SolverVectors<_F>> getScratchS() const {
      return S;
    }
    void setGuessScratch(std::shared_ptr<SolverVectors<_F>> scr) {
      Guess = scr;
    }
    void setSubspaceScratch(std::shared_ptr<SolverVectors<_F>> scr) {
      vecs = scr;
    }
    void setSigmaVecScratch(std::shared_ptr<SolverVectors<_F>> scr) {
      sigmaVecs = scr;
    }
    void setScratchR(std::shared_ptr<SolverVectors<_F>> scr) {
      R = scr;
    }
    void setScratchS(std::shared_ptr<SolverVectors<_F>> scr) {
      S = scr;
    }

    void alloc() override {

      IterDiagonalizer<_F>::alloc();

      // NO MPI
      // ROOT_ONLY(this->comm_);

      // Allocate Davidson specific Memory
      this->RelRes  = CQMemManager::get().malloc<double>(this->nGuess_);
    }

    bool runMicro() override;

    void restart() override;

    void setEnergySpecific(std::vector<std::pair<double, size_t>> eRefs,
                           bool AbsoluteES = false, double ABSshift = 0.) {

      this->AbsoluteES = AbsoluteES;
      // double check nRoots is not smaller than total number of high energy roots
      size_t nHighR = 0;
      for (auto & pair: eRefs) {
        nHighR += pair.second;
        if (AbsoluteES) pair.first = pair.first - ABSshift;
      }
      if (this->nRoots_ < nHighR) CErr("Total #Roots required by energy specific is more than nRoots.");
      this->nLowERoots     = this->nRoots_ - nHighR;
      this->energyRefs = eRefs;

    }

    void doLeftEigenvector() {this->DoLeftEigVec = true; }

//    void setSortByDistance(){
//      this->sortByDistance = true;
//    }

    virtual void setGuess(size_t nGuess,
        std::function<void(size_t, SolverVectors<_F> &, size_t)> func) override {

      if( nGuess != this->nGuess_ )
        CErr("Davison Requires nGuess = nGuess_",std::cout);

      if (not Guess or Guess->size() < this->nGuess_)
        Guess = this->vecGen_(nGuess);

      // NO MPI
      // ROOT_ONLY(this->comm_);

      func(nGuess, *Guess, this->N_);

    }

    void clear_scratch() {
      Guess     = nullptr;
      vecs      = nullptr;
      sigmaVecs = nullptr;
      R         = nullptr;
      S         = nullptr;
    }


  }; // Davidson


  template <typename _F>
  class GMRES : public IterLinearSolver<_F> {


    std::shared_ptr<SolverVectors<_F>> HHR_ = nullptr;
    std::shared_ptr<SolverVectors<_F>> U_   = nullptr;
    _F * W_   = nullptr;
    _F * J_   = nullptr;
    _F * R_   = nullptr;

  public:

    using VecsGen_t = typename IterSolver<_F>::VecsGen_t;
    using LinearTrans_t = typename IterLinearSolver<_F>::LinearTrans_t;
    using Shift_t       = typename IterLinearSolver<_F>::Shift_t;


    GMRES(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      double conv,
      const LinearTrans_t &linearTrans,
      const LinearTrans_t &preNoShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterLinearSolver<_F>(c,N,MSS,1,MSS,conv,linearTrans,preNoShift,
                          vecGen,shiftVec){ }

    GMRES(
      MPI_Comm c,
      const size_t N,
      const size_t MSS,
      double conv,
      const LinearTrans_t &linearTrans,
      const Shift_t &preShift,
      const VecsGen_t &vecGen = VecsGen_t(),
      const Shift_t &shiftVec = Shift_t()) :
     IterLinearSolver<_F>(c,N,MSS,1,MSS,conv,linearTrans,preShift,
                          vecGen,shiftVec){ }


    ~GMRES() {

      if(W_  )  CQMemManager::get().free(W_);
      if(J_  )  CQMemManager::get().free(J_);
      if(R_  )  CQMemManager::get().free(R_);

    }


    void alloc() {

      // Standard linear solver allocations
      IterLinearSolver<_F>::alloc();

      size_t nBatch = this->rhsBS * this->shiftBS;
      size_t MSSnBatch = this->mSS_ * nBatch;

      HHR_ = this->vecGen_(nBatch);
      U_   = this->vecGen_(MSSnBatch);

      // No MPI for GMRES
      // ROOT_ONLY(this->comm_);

      J_   = CQMemManager::get().malloc<_F>(2 * MSSnBatch);
      R_   = CQMemManager::get().malloc<_F>(MSSnBatch * this->mSS_);
      W_   = CQMemManager::get().malloc<_F>(MSSnBatch + nBatch);

    }

    void runBatch(size_t nRHS, size_t nShift,
                  std::shared_ptr<SolverVectors<_F>> RHS, _F *shifts,
                  std::shared_ptr<SolverVectors<_F>> SOL, double *RHSNorm );

  };

  

};

