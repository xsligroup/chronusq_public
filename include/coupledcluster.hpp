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
#ifdef CQ_HAS_TA

#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <singleslater.hpp>
#include <tiledarray.h>
#include <coupledcluster/TADIIS.hpp>
#include <itersolver/solvervectors.hpp>
#include <coupledcluster/TAManager.hpp>

//#define DEBUG_CCSD

namespace ChronusQ {

  std::pair<double, char> memSize(size_t mem);

  enum class JobType;

  enum class CC_TYPE { CCSD};

  struct CoupledClusterSettings {
    CC_TYPE cctype = CC_TYPE::CCSD;
    double eConv = 1e-8;  // Convergence criteria of energy
    double tConv = 1e-6;  // Convergence criteria of T amplitudes
    int maxiter = 1000;   // Maximum # of iteration
    bool useDIIS = true;  // Flag for DIIS
    size_t nDIIS = 8;     // # of vectors to keep in DIIS
    size_t blksize = 32;   // Block size of TildeArray
    size_t nEvariation = 0;// Variation of number of electrons
    bool rebuildFock = false; // Rebuild Fock matrix from Core Hamiltonian
    std::vector<size_t> frozen_occupied;
    std::vector<size_t> frozen_virtual; 
  };

  enum class EOM_HBAR_TYPE { EXPLICIT, IMPLICIT, DEBUG };
  enum class EOM_DIAG_METHOD { FULL, DAVIDSON, GPLHR };

  struct EOMSettings{
    size_t nroots = 1;
    EOM_HBAR_TYPE hbar_type = EOM_HBAR_TYPE::IMPLICIT;
    EOM_DIAG_METHOD diag_method = EOM_DIAG_METHOD::DAVIDSON;
    std::vector<size_t> cvs_core;
    std::vector<size_t> cvs_virtual;
    std::vector<size_t> frozen_occupied;
    std::vector<size_t> frozen_virtual;
    size_t davidson_whenSc = 1;
    size_t davidson_max_macro_iter = 20;
    size_t davidson_max_micro_iter = 128;
    double davidson_residual_conv = 1e-5;
    double davidson_eigen_vector_conv = 1e-6;
    double davidson_eigen_value_conv = 1e-7;
    bool davidson_check_residual = true;
    bool davidson_check_eigen_vector = false;
    bool davidson_check_eigen_value = true;
    bool davidson_conv_on_GramSchmidt = true;
    size_t davidson_subspace_multiplier = 8;
    size_t davidson_guess_multiplier = 3;
//    double davidson_Eref = 0.0;
    size_t davidson_nLowRoots = 1;
    std::vector<std::pair<double, size_t>> davidson_Eref;
    bool davidson_ErefAbs = true;
    double davidson_preCond_small = 1e-12;
//    bool davidson_sort_by_distance = false;
    bool davidson_biortho = true;
    size_t GramSchmidt_NRe = 1;
    double GramSchmidt_eps = 1e-12;
    bool oscillator_strength = false;
    bool save_hamiltonian = false;

    bool doCVS() const {
      return cvs_core.size() > 0 or
      cvs_virtual.size() > 0 or
      frozen_occupied.size() > 0 or
      frozen_virtual.size() > 0;
    }

    void printEOMCCSettings(std::ostream &out);

    size_t estimate_mem_peak() const;
  };

  template <typename MatsT>
  struct CCIntermediates {
    using TArray = TA::TArray<MatsT>;

    // Auxiliary variables and functions to help initialization of TA tensors
    size_t nVir;
    size_t nOcc;
    char aoLabel = 'a';
    char vLabel = 'v';  // non-frozen virtual space, include LUMO, Rydberg
    char oLabel = 'o';  // non-frozen occupied space, include HOMO, core
    char VLabel = 'V';  // all virtual space, include LUMO, Rydberg, Free elec
    char OLabel = 'O';  // all occupied space, include HOMO, core, deep core
    char hLabel = 'h';  // HOMO space, active occupied space, CVS Valence
    char lLabel = 'l';  // LUMO space, active virtual space, CVS virtual valence
    char cLabel = 'c';  // core space, external occupied space, CVS core
    char rLabel = 'r';  // Rydberg space, external virtual space, CVS continuum
    char dLabel = 'd';  // deep core space, frozen core space
    char fLabel = 'f';  // free electron space, frozen virtual space

    // Integrals
    double E_ref;
    double E_fzc;

    std::vector<double> eps;
    std::map<std::string,TArray> fockMatrix;
    std::map<std::string,TArray> antiSymMoInts;
    std::map<std::string,TArray> muMatrix;

    // Intermediates
    TArray D_ai;
    TArray D_abij;
    TArray tau;
    TArray tilde_tau;
    TArray F_ae;
    TArray F_mi;
    TArray F_me;
    TArray W_mnij;
    TArray W_abef;
    TArray W_mbej;
    TArray W_mnie;
    TArray W_amef;
    TArray W_mbij;
    TArray W_abei;
    TArray G_ae;
    TArray G_mi;

    // Amplitudes
    std::shared_ptr<EOMCCSDVector<MatsT>> T, Lg;

    template <typename IntsT>
    void initializeIntegrals(const cqmatrix::PauliSpinorMatrices<MatsT> &aoCoreH,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoFock,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoTwoeH,
                             const TwoPInts<IntsT> &aoTPI,
                             const MultipoleInts<IntsT> &lenElectric,
                             CoupledClusterSettings& ccSettings,
                             EOMSettings& eomSettings,
                             MatsT *mo, size_t nO_, size_t nV_,
                             size_t blksize, double nucRepEnergy,
                             bool rebuildFock = false);

    size_t estimate_mem_peak(size_t nDIIS) const;

    ~CCIntermediates();
  };

  template <typename MatsT, typename IntsT>
  class EOMCCSD;

  template <typename MatsT, typename IntsT>
  class CCSD
  {

    using TArray = TA::TArray<MatsT>;
  protected:
    // Input parameters
    SafeFile savFile_;
    char vLabel_;
    char oLabel_;
    CoupledClusterSettings ccSettings_;
    CCIntermediates<MatsT> &intermediates_;

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta;
    std::map<std::string,TArray> &antiSymMoints;

    // Amplitudes
    EOMCCSDVector<MatsT> &T_;
    TArray &T1_;
    TArray &T2_;

    // Result correlation energy
    MatsT CorrE = 0.0;

    // Intermediates
    TArray &tau_;
    TArray &tilde_tau_;
    TArray Fae_;
    TArray &Fae_wDiag_;
    TArray Fmi_;
    TArray &Fmi_wDiag_;
    TArray &Fme_;
    TArray &Wmnij_;
    TArray &Wabef_;
    TArray &Wmbej_;
    TArray &Dai_;
    TArray &Dabij_;

    // Intermediates of Xiaolin Liu
    TArray A_1;
    TArray A_2;
    TArray A_3;
    TArray A_4;
    TArray B_1;
    TArray B_2;
    TArray B_3;
    TArray B_5;
    TArray B_6;

  public:    
    CCSD(const SafeFile &savFile,
         CCIntermediates<MatsT> &intermediates, const CoupledClusterSettings &ccSettings);

    double computeReferenceEnergy(const InCore4indexTPI<MatsT> &moTPI, double nucRepEnergy);
    void getCorrEnergy();
    void printAnalysis();
    void run();
    void runConventional();
    void printBanner(double Eref) const;

    // Reference Gauss, Stanton, J. Chem. Phys. 103, 3561 (1995) DOI:10.1063/1.470240
    void build_tau_and_tilde_tau();
    // One-body intermediates
    void build_tilde_Fae();
    void build_tilde_Fmi();
    void build_tilde_Fme();
    // Two-body intermediates
    void build_tilde_Wmnij();
    void build_tilde_Wabef();
    void build_tilde_Wmbej();
    // Build all intermediates
    void buildIntermediates();
    // T-amplitude equations
    void updateT1(const TArray T1_old, const TArray T2_old);
    void updateT2(const TArray T1_old, const TArray T2_old);

    //One-body intermediates and tensors of Xiaolin Liu
    void formA_1();
    void formA_2();
    void formA_3();
    void formA_4();
    //Two-body intermediates and tensors of Xiaolin Liu
    void formB_1();
    void formB_2();
    void formB_3();
    void formB_5();
    void formB_6();
    //T-amplitude equations of Xiaolin Liu
    void updateT1_xxl(const TArray T1_old, const TArray T2_old);
    void updateT2_xxl(const TArray T1_old, TArray T2_old);

    void initRanges();
    void initIntermediates();
    void initIntermediates_xll();
    void initAmplitudes();
    void doDIIS(EOMCCSDVector<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis );

    ~CCSD();

  };// class CCSD

  enum class EOMCCEigenVecType { RIGHT, LEFT};

  template <typename MatsT, typename IntsT>
  class EOMCCSD {
  protected:
    using TArray = TA::TArray<MatsT>;

    SafeFile savFile_;
    CCIntermediates<MatsT> &intermediates_;
    char vLabel_;
    char oLabel_;
    size_t nOVshift_;
    size_t nO2shift_;
    size_t nV2shift_;

    size_t Hbar_dim = 0;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;
    EOMCCSDVector<MatsT> &Lg_;

    //Intermediates
    // [Gauss:1995:3561] Table III
    // [Asthana:2019:4102] Appendix
    TArray &tau;
    TArray &F_ae;
    TArray &F_mi;
    TArray &F_me;
    TArray &W_mnij;
    TArray &W_abef;
    TArray &W_mbej;
    TArray &W_mnie;
    TArray &W_amef;
    TArray &W_mbij;
    TArray &W_abei;

    // Lambda intermediates:[Gauss:1995:3561] Table III (c)
    TArray &G_ae;
    TArray &G_mi;
    TArray &D_ai;
    TArray &D_abij;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

    MatsT* theta = nullptr;
    std::shared_ptr<SolverVectors<MatsT>> R_;
    std::shared_ptr<SolverVectors<MatsT>> L_;

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta;
    std::map<std::string,TArray> &antiSymMoints;
    std::map<std::string,TArray> &muMatrix;

    CoupledClusterSettings ccSettings_;
    EOMSettings eomSettings;

    // CVS
//    bool doCVS_ = false;
    std::vector<size_t> CVSOindices_;
    std::vector<size_t> CVSVindices_;
    size_t nCVSOCore_;
    size_t nCVSOValance_;
    size_t nCVSOActive_;
    size_t nCVSVContinuum_;
    size_t nCVSVValance_;
    size_t nCVSVActive_;
    size_t CVSoutOfBound_;

    std::vector<std::vector<size_t>> CVSabIndices_;
    std::vector<std::vector<size_t>> CVSijIndices_;

    void assignCVSIndices();

  public:

    double signD(size_t &a, size_t &b, size_t &i, size_t &j) const;
    size_t CVStoCompoundS(size_t a, size_t i) const;
    size_t CVStoCompoundD(size_t a, size_t b, size_t i, size_t j) const;
    size_t CVStoCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const;
    std::pair<size_t, double> CVStoCompoundSD(size_t e, size_t m,
                                              size_t a, size_t b, size_t i, size_t j, size_t ldH) const;
    std::pair<size_t, double> CVStoCompoundDS(size_t a, size_t b, size_t i, size_t j,
                                              size_t e, size_t m, size_t ldH) const;
    std::pair<size_t, double> CVStoCompoundDD(size_t a, size_t b, size_t i, size_t j,
                                              size_t c, size_t d, size_t k, size_t l, size_t ldH) const;
    bool CVSisInBound(size_t idx) const { return idx < CVSoutOfBound_; }
    size_t CVSoneBodySize() const { return nOVshift_; }

    EOMCCSD(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMCCSD();

    void formF_ae();

    void formF_mi();

    void formW_mnij();

    void formW_abef();

    void formW_mbej();

    void formW_mnie();

    void formW_amef();

    void formW_mbij();

    void formW_abei();

    void formR1_tilde(const TArray &R1, const TArray &R2, TArray &tildeR1) const;

    void formR2_tilde(const TArray &R1, const TArray &R2, TArray &tildeR2) const;

    //Lambda iterations
    void initilizeLambda();
    void updateG_ae(const TArray &L2, TArray &G_ae) const;
    void updateG_mi(const TArray &L2, TArray &G_mi) const;
    void formL1_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL1, bool subtractDiagonal = false) const;
    void formL2_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL2, bool subtractDiagonal = false) const;
    void runLambda();

    size_t getHbarDim(bool includeZeroBody = false) const { return Hbar_dim + (includeZeroBody ? 1 : 0); }

    void run();
    void initilizeEOMCCSD();
    void formEOMIntermediates();

    //Helper functions for diagonalization. Davidson/GPLHR
    void buildSigma(const TArray &V1, const TArray &V2, TArray &HV1, TArray &HV2, EOMCCEigenVecType vecType) const;

    cqmatrix::Matrix<MatsT> buildHbarCVS(bool includeGroundState) const;
    void buildDiag(MatsT * diag) const;
    void buildHbar_sigma(MatsT * out, bool diagOnly) const;

    void full_diagonalization();

    void buildRightZeroBody(size_t nVec);

    // Density functions
    void initilizeDensity();
    MatsT calcOscillatorStrength(size_t i);
    void buildDensity(const TArray& t1, const TArray& t2,  const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ij(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ab(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);
    void formRho_ia(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ai(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);

    // Get and Set results
    MatsT* getTheta() const { return theta; }
    void setTheta(MatsT *eVals, size_t n) {
      if (theta) CQMemManager::get().free(theta);
      theta = CQMemManager::get().malloc<MatsT>(n);
      std::copy_n(eVals, n, theta);
    }
    std::shared_ptr<SolverVectors<MatsT>> getR() const { return R_; }
    void setR(std::shared_ptr<SolverVectors<MatsT>> R) { R_ = R; }
    std::shared_ptr<SolverVectors<MatsT>> getL() const { return L_; }
    void setL(std::shared_ptr<SolverVectors<MatsT>> L) { L_ = L; }
    EOMCCSDVector<MatsT>& getLg() const { return Lg_; }

  };

  void runCoupledCluster(JobType jobType, Molecule &mol, std::shared_ptr<SingleSlaterBase> ss,
                         std::shared_ptr<IntegralsBase> aoints,
                         SafeFile &rstFile, CQInputFile &input, std::ostream &output);

}; // namespace ChronusQ
#endif
