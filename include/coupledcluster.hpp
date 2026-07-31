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
#include <itersolver.hpp>
#include <itersolver/solvervectors.hpp>
#include <coupledcluster/TAManager.hpp>

//#define DEBUG_CCSD

namespace ChronusQ {

  std::pair<double, char> memSize(size_t mem);
  enum class JobType;

  enum class CC_TYPE { CCSD, DFCCSD, CCSDT};

  template <typename MatsT>  class CCBase;
  template <typename MatsT>  class CCSD;
  template <typename MatsT>  class RCCSD;
  template <typename MatsT>  class DFCCSD;
  template <typename MatsT>  class CCSDT;
  template <typename MatsT> struct CCIntermediates;


  struct CoupledClusterSettings {
    CC_TYPE cctype = CC_TYPE::CCSD;
    double eConv = 1e-8;  // Convergence criteria of energy
    double tConv = 1e-6;  // Convergence criteria of T amplitudes
    int maxiter = 1000;   // Maximum # of iteration
    bool useDIIS = true;  // Flag for DIIS
    size_t nDIIS = 8;     // # of vectors to keep in DIIS
    size_t blksize = 32;   // Block size of TildeArray
    size_t nEvariation = 0;// Variation of number of electrons
    double denomshift = 0.0;    // Denominator shift to help convergence
    bool pertT3 = false;   // Flag for CCSD(T)
    bool crcc = false; // Flag for CR-CC(2,3)
    bool triplesMPI = false; // Flag for MPI-parallel (T) / CR-CC(2,3) correction
    size_t triples_begin = 0; // Indices to start loop
    size_t triples_end   = 0; // Indices to end loop
    bool loop_abc = true; // Use loop over abc for triples correction by default
    bool CCSDinit = false; // Flag to run CCSD as initial guess for CCSDT
    bool restart = false; // read T from previous calculations
    bool save = false; // save T for next calculation
    bool rebuildFock = false; // Rebuild Fock matrix from Core Hamiltonian
    std::vector<size_t> frozen_occupied;
    std::vector<size_t> frozen_virtual;
    bool skipSCF = false; // Skip SCF calculation
    bool skipCC = false; // Skip CC calculation, just read T amp from bin file
    bool computeDipole = false; // Flag for compute ground state density and dipole
  };

  enum class EOMCCEigenVecType { RIGHT, LEFT};


  enum class EOM_TYPE { EE, IP, EA, DIP };
  enum class EOM_IMPLEMENTATION { EOMCCSD, EOMRCCSD, CVSEOMCCSD, EOMIP_2h1p, EOMIP_3h2p, EOMEA, EOMDIP_3h1p, EOMDIP_4h2p, SOMETHING_WRONG };
  enum class EOM_HBAR_TYPE { EXPLICIT, IMPLICIT, DEBUG };
  enum class EOM_DIAG_METHOD { FULL, DAVIDSON, GPLHR };

  struct EOMSettings{
    size_t nroots = 1;
    size_t nO = 0;
    EOM_TYPE eom_type = EOM_TYPE::EE;
    EOM_HBAR_TYPE hbar_type = EOM_HBAR_TYPE::IMPLICIT;
    EOM_DIAG_METHOD diag_method = EOM_DIAG_METHOD::DAVIDSON;
    EOM_IMPLEMENTATION eom_implementation = EOM_IMPLEMENTATION::SOMETHING_WRONG;
    std::vector<size_t> cvs_core;
    std::vector<size_t> active_occupied;
    std::vector<size_t> active_virtual;
    std::vector<size_t> external_virtual;
    size_t ip_level = 2;
    size_t dip_max_external_hole = 1;
    size_t dip_max_external_particle = 1;
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
    bool all_excited_dipole = false;
    bool save_hamiltonian = false; //currently only work for full_diag
    bool save_r = false;
    bool save_l = false;
    bool print_large_amplitude = false;
    bool ccs_guess = false;
    bool skip_r = false;
    bool restart_r = false;
    bool restart_l = false;
    bool singlet_only = false;

    EOMSettings() {
      find_eom_implementation();
    }

    bool containActive() const {
      size_t n_core = cvs_core.size();
      size_t n_extvir = external_virtual.size();
      if (nO == 0) 
        return  n_core > 0 or n_extvir > 0;
      return (n_core > 0 && n_core < nO) or (n_extvir > 0 && n_extvir < nO);
    }

    void printEOMCCSettings(std::ostream &out, const CoupledClusterSettings& ccSettings);

    void validateOrbitalSpaces(
      size_t nO, size_t nV, const CoupledClusterSettings& ccSettings) const;

    void assignUnspecifiedOrbitalToSpaces(
      size_t nO, size_t nV, const CoupledClusterSettings& ccSettings);


    template<typename MatsT>
    size_t estimate_mem_peak(const CoupledClusterSettings& ccSettings) const;
    EOM_IMPLEMENTATION find_eom_implementation();
      private: 
    size_t n_MBExpansion() const;
    size_t MBExpansionSize() const;
    size_t intermediate_mem(const CoupledClusterSettings& ccSettings) const;
  };

  template <typename MatsT>
  struct CCIntermediates {
    using TArray = TA::TArray<MatsT>;

    // Auxiliary variables and functions to help initialization of TA tensors
    size_t nOcc;
    size_t nVir;
    size_t nRI;
    char aoLabel = 'a';
    char auxLabel = 'b';
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
    char tLabel = 't';  // TiledArray block size

    // Integrals
    double  E_ref;
    double  E_fzc;
    double  E_cc;
    std::vector<double> eps;
    std::map<std::string,TArray> fockMatrix;
    std::map<std::string,TArray> antiSymMoInts;
    std::map<std::string,TArray> moInts;
    std::map<std::string,TArray> riMoInts;
    std::map<std::string,TArray> muMatrix;

    // Dipoles
    std::array<double, 3> Mu_fzc;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;
    
    // Amplitudes
    std::shared_ptr<MBExpansion<MatsT>> T;

    // Intermediates (CCSD)
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

    // Intermediates (RHF reference CCSD)
    TArray W_mbej_baab;
    TArray W_mbej_baba;

    // Intermediates (CCSDT)
    TArray Pabcijk;

    // Intermediates (EOMDIP)
    std::map<std::string, TArray> tempOps;
    std::map<std::string, TArray> sigmaOps;
    TArray tempPerm_oo;
    TArray Id_oo;
    TArray Id_oooo;

    void reorderMOs(MatsT *mo, size_t nO, size_t nV,
                    CoupledClusterSettings& ccSettings,
                    EOMSettings& eomSettings);

    template <typename IntsT>
    void initializeIntegrals(const cqmatrix::PauliSpinorMatrices<MatsT> &aoCoreH,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoFock,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoTwoeH,
                             const TwoPInts<IntsT> &aoTPI,
                             const MultipoleInts<IntsT> &lenElectric,
                             CoupledClusterSettings& ccSettings,
                             EOMSettings& eomSettings,
                             MatsT *mo, size_t nO_, size_t nV_,
                             size_t blksize, double nucRepEnergy, CC_TYPE cctype,
                             double denomshift_, bool pertT3,
                             bool rebuildFock);

    template <typename IntsT>
    void initializeIntegrals_rhfRef(const cqmatrix::PauliSpinorMatrices<MatsT> &aoCoreH,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoFock,
                             const cqmatrix::PauliSpinorMatrices<MatsT> &aoTwoeH,
                             const TwoPInts<IntsT> &aoTPI,
                             const MultipoleInts<IntsT> &lenElectric,
                             CoupledClusterSettings& ccSettings,
                             EOMSettings& eomSettings,
                             MatsT *mo, size_t nO_, size_t nV_,
                             size_t blksize, double nucRepEnergy, CC_TYPE cctype,
                             double denomshift_, bool pertT3,
                             bool rebuildFock);
    
    std::shared_ptr<CCBase<MatsT>> build_cc(
                                                   const SafeFile &savFile,
                                                   CoupledClusterSettings & ccSettings, bool rhfRef = false) {
      std::shared_ptr<CCBase<MatsT>> cc = nullptr; 
      if (ccSettings.cctype == CC_TYPE::CCSD) {
        if (rhfRef) {
          cc = std::dynamic_pointer_cast<CCBase<MatsT>>(
             std::make_shared<RCCSD<MatsT>>(savFile,*this, ccSettings));

        } else {
          cc = std::dynamic_pointer_cast<CCBase<MatsT>>(
             std::make_shared<CCSD<MatsT>>(savFile,*this,ccSettings));
        }
      }
      else if (ccSettings.cctype == CC_TYPE::DFCCSD) {
        cc = std::dynamic_pointer_cast<CCBase<MatsT>>(
           std::make_shared<DFCCSD<MatsT>>(savFile,*this,ccSettings));
      }
      else if (ccSettings.cctype == CC_TYPE::CCSDT) {
        cc = std::dynamic_pointer_cast<CCBase<MatsT>>(
           std::make_shared<CCSDT<MatsT>>(savFile,*this,ccSettings)); 
      }
      return cc;
    }
    
    ~CCIntermediates();
 
  };

  template <typename MatsT>
  class EOMCCSD;
  
  template <typename MatsT>
  class CCBase
  {

    using TArray = TA::TArray<MatsT>;
    protected:
    // Input parameters
    SafeFile savFile_;
    size_t nO_;
    size_t nV_;
    CoupledClusterSettings ccSettings_;
    CCIntermediates<MatsT> &intermediates_;

    // Auxiliary variables and functions to help initialization of TA tensors
    char vLabel_;
    char oLabel_;

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta;
    std::map<std::string,TArray> &antiSymMoints;
    std::map<std::string,TArray> &riMoInts;

    // Amplitudes
    MBExpansion<MatsT> &T_;
    TArray &T1_;
    TArray &T2_;

    // Intermediates
    TArray &Dai_;
    TArray &Dabij_;

   public:
    CCBase(const SafeFile &savFile,
         CCIntermediates<MatsT> &intermediates, const CoupledClusterSettings &ccSettings):
      savFile_(savFile),
      nO_(intermediates.nOcc), nV_(intermediates.nVir),
      intermediates_(intermediates), ccSettings_(ccSettings),
      oLabel_(intermediates.oLabel),
      vLabel_(intermediates.vLabel),
      fockMatrix_ta(intermediates.fockMatrix),
      antiSymMoints(intermediates.antiSymMoInts),
      riMoInts(intermediates.riMoInts),
      T_(*intermediates.T),
      T1_(intermediates.T->get_tensor("OneBody")),
      T2_(intermediates.T->get_tensor("TwoBody")),
      Dai_(intermediates.D_ai),
      Dabij_(intermediates.D_abij){}

    // Result correlation energy
    MatsT CorrE;
    virtual size_t estimate_mem_peak() const = 0;


    virtual void getCorrEnergy();
    virtual void runConventional() = 0;
    void printBanner(double Eref) const;
    
    virtual void buildIntermediates(){}
    virtual void initIntermediates(){}
    virtual void initAmplitudes() = 0;
    virtual void doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis );

    virtual void cleanMemory(){}

    virtual void run(){}
  };
  
  template <typename MatsT>
  class CCSD : public CCBase<MatsT>
  {

    using TArray = TA::TArray<MatsT>;
   protected:
    // Intermediates
    TArray &tau_;
    TArray &tilde_tau_;
    TArray &Fae_;
    TArray &Fmi_;
    TArray &Fme_;
    TArray &Wmnij_;
    TArray &Wabef_;
    TArray &Wmbej_;

    // CCSD(T) objects
    std::vector<double> &eps;
    MatsT PertT3Energy = 0.0;

  public:    
    CCSD(const SafeFile &savFile,
         CCIntermediates<MatsT> &intermediates, const CoupledClusterSettings &ccSettings);

    virtual void initAmplitudes();
    virtual void initIntermediates();
    //virtual void doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis );
    virtual void buildIntermediates();
    //virtual void printBanner(double Eref) const;
    virtual void runConventional();
    virtual void run();
    void cleanMemory();
    void printAnalysis();
    size_t estimate_mem_peak() const override;

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
    // T-amplitude equations
    void updateT1(const TArray T1_old, const TArray T2_old);
    void updateT2(const TArray T1_old, const TArray T2_old);

    // CCSD(T) run
    void runPertT3ijk(const TArray &T1_, const TArray &T2_);
    void runPertT3abc(const TArray &T1_, const TArray &T2_);

    ~CCSD();

  };// class CCSD

  template <typename MatsT>
  class RCCSD : public CCBase<MatsT>
  {

    using TArray = TA::TArray<MatsT>;
  protected:
    // Intermediates
    std::map<std::string,TArray> &moInts;

    // Intermediates
    TArray &tau_RHF;
    TArray &tilde_tau_RHF;
    TArray &Fae_RHF;
    TArray &Fmi_RHF;
    TArray &Fme_RHF;
    TArray &Wmnij_RHF;
    TArray &Wabef_RHF;
    TArray &Wmbej_RHF_baab;
    TArray &Wmbej_RHF_baba;

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta_RHF;

    // Amplitudes
    MBExpansion<MatsT> &T_RHF;
    TArray &T1_RHF;
    TArray &T2_RHF;
    MBExpansion<MatsT> T_RHF_old;

    // Intermediates
    TArray &Dai_RHF;
    TArray &Dabij_RHF;

  public:

    static void convertTto2C(const MBExpansion<MatsT> &T_RHF, MBExpansion<MatsT> &T_2C);
    static void convertTtoRHF(const MBExpansion<MatsT> &T_2C, MBExpansion<MatsT> &T_RHF);

  public:
    RCCSD(const SafeFile &savFile,
         CCIntermediates<MatsT> &intermediatesRref,
         const CoupledClusterSettings &ccSettings);

    virtual void initAmplitudes();
    virtual void initIntermediates();
    virtual void doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis ) override;
    void doDIIS(MBExpansion<MatsT> &T_old, MBExpansion<MatsT> &T_new, std::shared_ptr<DIISTA<MatsT>> diis );
    virtual void buildIntermediates();
    //virtual void printBanner(double Eref) const;
    virtual void runConventional();
    virtual void run();
    virtual void getCorrEnergy() override;
    void cleanMemory();
    size_t estimate_mem_peak() const override;

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
    // T-amplitude equations
    void updateT1(const TArray T1_old, const TArray T2_old);
    void updateT2(const TArray T1_old, const TArray T2_old);

    ~RCCSD();

  };// class RCCSD

  template <typename MatsT>
  class DFCCSD : public CCBase<MatsT>
  {

  using TArray = TA::TArray<MatsT>;
  protected:
    // Intermediates
    // T1-transformed F_N
    TArray F_ov;
    TArray F_oo;
    TArray F_vv;
    TArray F_vo;

    // T1-transformed B(Q)_pq
    TArray B_oo;
    TArray B_vv;
    TArray B_vo;

    // Still need for EOM or CR correction after convergence
    // Consider moving to EOM later? Eh
    TArray &tau_;
    TArray &tilde_tau_;
    TArray &Fae_;
    TArray &Fmi_;
    TArray &Fme_;
    TArray &Wmnij_;
    TArray &Wabef_;
    TArray &Wmbej_;

  public:
    DFCCSD(const SafeFile &savFile,
         CCIntermediates<MatsT> &intermediates, const CoupledClusterSettings &ccSettings);

    virtual void initAmplitudes();
    virtual void initIntermediates();
    void t1TransformIntegrals();
    //virtual void doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis );
    virtual void buildIntermediates();
    //virtual void printBanner(double Eref) const;
    virtual void runConventional();
    virtual void run();
    size_t estimate_mem_peak() const override;

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
    // T-amplitude equations
    void updateT1(const TArray &T1_old, const TArray &T2_old);
    void updateT2(const TArray &T2_old);

    void cleanMemory();
    
    ~DFCCSD();

  };// class DFCCSD

  template <typename MatsT>
  class CCSDT : public CCBase <MatsT>
  {

    using TArray = TA::TArray<MatsT>;
  protected:

    // Amplitudes
    TArray &T3_;

    // Intermediates
    TArray &tau_;
    TArray &tilde_tau_;
    TArray &Fae_;
    TArray &Fmi_;
    TArray &Fme_;
    TArray &Wmnie_;
    TArray &Wamef_;
    TArray &Wmnij_;
    TArray &Wabef_;
    TArray &Wmbej_;
    TArray &Wmbij_;
    TArray &Wabei_;

    // Misc
    bool include_t3_;
    std::vector<double> &eps;

  public:
    CCSDT(const SafeFile &savFile,
          CCIntermediates<MatsT> &intermediates, const CoupledClusterSettings &ccSettings);

    virtual void initAmplitudes();
    virtual void initIntermediates();
    //virtual void doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis );
    virtual void buildIntermediates();
    //virtual void printBanner(double Eref) const;
    virtual void runConventional();
    virtual void run();
    void cleanMemory();
    void printAnalysis();
    size_t estimate_mem_peak() const override;

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
    // more intermediates that are not used in ground-state CCSD
    void formW_mnie();
    void formW_amef();
    void formW_mbij();
    void formW_abei();
    // T-amplitude equations
    void updateT1(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old);
    void updateT2(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old);
    void updateT3(const TArray &T1_old, const TArray &T2_old, const TArray &T3_old);

    ~CCSDT();

  };// namespace CCSDT

  template <typename MatsT>
  class EOMCCBase {

  protected:
    using TArray = TA::TArray<MatsT>;

    SafeFile savFile_;
    CCIntermediates<MatsT> &intermediates_;
    std::vector<std::string> tensor_builder_;

    size_t Hbar_dim = 0;
    std::map<std::string, size_t> Hbar_dimension_offsets;

    // Amplitudes
    std::shared_ptr<MBExpansion<MatsT>> Lg_;

    dcomplex* theta = nullptr;
    std::shared_ptr<SolverVectors<MatsT>> R_;
    std::shared_ptr<SolverVectors<MatsT>> L_;

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta;
    std::map<std::string,TArray> &antiSymMoints;
    std::map<std::string,TArray> &riMoints;
    std::map<std::string,TArray> &muMatrix;

    CoupledClusterSettings ccSettings_;
    EOMSettings eomSettings;

    typename Davidson<MatsT>::LinearTrans_t funcEOM;
    typename Davidson<MatsT>::LinearTrans_t funcRaw;
    typename Davidson<MatsT>::LinearTrans_t funcDebug;
    typename Davidson<MatsT>::LinearTrans_t PCEOM;
    typename Davidson<MatsT>::LinearTrans_t PCRaw;
    typename Davidson<MatsT>::LinearTrans_t PCDebug;

  public:

    void davidsonSolve();
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder(){return typename Davidson<MatsT>::VecsGen_t();}
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){return typename Davidson<MatsT>::LinearTrans_t(); }
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){return typename Davidson<MatsT>::LinearTrans_t();}

    virtual size_t nOVShift() const {return 0;}
    virtual bool isInBound(size_t idx) const {return false;} 
    virtual size_t oneBodySize() const {return 0;}
    void print_largest_values_and_position(size_t i, MatsT* r_vector, size_t Hbar_dim, const EOMSettings& eomSettings, CCIntermediates<MatsT> & intermediates);


    EOMCCBase(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMCCBase();

    //Lambda iterations
    void initializeGroundStateLambda() {
      Lg_ = std::make_shared<MBExpansion<MatsT>>(tensor_builder_);
    }
    virtual void initializeLambda(){}
    virtual void runLambda(){}
    virtual void fillGuess(MatsT *guess_vec, size_t n_vec)const{}

    size_t getHbarDim(bool includeZeroBody = false) const { return Hbar_dim + (includeZeroBody ? 1 : 0); }

    void prepEOMCC();
    virtual void initializeEOMCC() {}
    virtual void formEOMIntermediates(){}

    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const{}

    void full_diagonalization();
    virtual cqmatrix::Matrix<MatsT> buildHbar(bool includeGroundState) const {return cqmatrix::Matrix<MatsT>(Hbar_dim);}
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const {}
    void buildHbar_sigma(MatsT * out, bool diagOnly) const;


    virtual void buildRightZeroBody(size_t nVec){}

    // Density functions
    virtual void initializeDensity(){}
    virtual dcomplex calcOscillatorStrength(size_t i){return 0.0;}

    // Get and Set results
    dcomplex* getTheta() const { return theta; }
    void setTheta(dcomplex *eVals, size_t n) {
      if (theta) CQMemManager::get().free(theta);
      theta = CQMemManager::get().malloc<dcomplex>(n);
      std::copy_n(eVals, n, theta);
    }
    std::shared_ptr<SolverVectors<MatsT>> getR() const { return R_; }
    void setR(std::shared_ptr<SolverVectors<MatsT>> R) { R_ = R; }
    std::shared_ptr<SolverVectors<MatsT>> getL() const { return L_; }
    void setL(std::shared_ptr<SolverVectors<MatsT>> L) { L_ = L; }
    MBExpansion<MatsT>& getLg() const { return *Lg_; }

    // Get EOM tensor dimension infomration
    std::vector<std::string> tensor_dimensions() const {return tensor_builder_;}

    virtual std::array<MatsT, 3> calcGroundDipole() {
      CErr("calcGroundDipole not implemented for this EOMCC type");
      return {0.0, 0.0, 0.0};
    }
    
  };

  template <typename MatsT>
  class EOMCCSD : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;

    // complete EOM Hamiltonian for full diagonalization and debugging
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nOVshift_;
    size_t nO2shift_;
    size_t nV2shift_;
    size_t nV3shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;

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
    TArray &Rho_ij;
    TArray &Rho_ab;
    TArray &Rho_ia;
    TArray &Rho_ai;

  public:


    EOMCCSD(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMCCSD();

  protected:
    std::vector<std::vector<size_t>> abIndices_;
    std::vector<std::vector<size_t>> ijIndices_;
    size_t outOfBound_;

  public:
    virtual size_t nOVShift() const {return nOVshift_;}
    inline double signD(size_t a, size_t b, size_t i, size_t j) const;
    size_t toCompoundS(size_t a, size_t i) const;
    size_t toCompoundD(size_t a, size_t b, size_t i, size_t j) const;
    size_t toCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const;
    std::pair<size_t, double> toCompoundSD(size_t e, size_t m,
                                              size_t a, size_t b, size_t i, size_t j, size_t ldH) const;
    std::pair<size_t, double> toCompoundDS(size_t a, size_t b, size_t i, size_t j,
                                              size_t e, size_t m, size_t ldH) const;
    std::pair<size_t, double> toCompoundDD(size_t a, size_t b, size_t i, size_t j,
                                                     size_t c, size_t d, size_t k, size_t l, size_t ldH) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nOVshift_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);
    
    void fillGuess(MatsT *guess_vec, size_t n_vec)const;
    
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
    virtual void initializeLambda();
    void updateG_ae(const TArray &L2, TArray &G_ae) const;
    void updateG_mi(const TArray &L2, TArray &G_mi) const;
    void formL1_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL1) const;
    void formL2_tilde(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL2) const;
    virtual void runLambda();

    virtual void initializeEOMCC();
    virtual void formEOMIntermediates();

    // for full_diagonalize  algorithm
    virtual cqmatrix::Matrix<MatsT> buildHbar(bool includeGroundState) const;
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;

    // For building CRCC denominator
    std::vector<MatsT> hbar1;
    std::vector<MatsT> hbar2;
    std::vector<MatsT> hbar3;

    MatsT CRCC23Energy_A = 0.0;
    MatsT CRCC23Energy_B = 0.0;
    MatsT CRCC23Energy_C = 0.0;
    MatsT CRCC23Energy_D = 0.0;

    void runCR(const MatsT &CorrE);
    void builddenom(const TArray &T2_);
    void runCRbatchijk(const TArray &T1_, const TArray &T2_, const TArray &L1, const TArray &L2);
    void runCRbatchabc(const TArray &T1_, const TArray &T2_, const TArray &L1, const TArray &L2);

    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec);


    // Density functions
    virtual void initializeDensity();
    virtual dcomplex calcOscillatorStrength(size_t i) override;

    void buildDensity(const TArray& t1, const TArray& t2,  const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    std::array<MatsT, 3> calcTransitionDipole(const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0,const TArray& l1, const TArray& l2, bool isSame = false);
    std::array<MatsT, 3> calcGround2ExcitedTransitionDipole(size_t i);
    std::array<MatsT, 3> calcExcited2GroundTransitionDipole(size_t i);
    std::array<MatsT, 3> calcExcited2ExcitedTransitionDipole(size_t i, size_t j);
    virtual std::array<MatsT, 3> calcGroundDipole() override;

    void formRho_ij(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ab(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);
    void formRho_ia(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ai(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);

  };

  template <typename MatsT>
  class EOMRCCSD : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;

    // RHF-based intermediates for EOM-RCCSD

    // Hamiltonian in TA implementation
    std::map<std::string,TArray> &fockMatrix_ta_RHF;
    std::map<std::string,TArray> &Moints;
    std::map<std::string,TArray> &muMatrix_RHF;

    char vLabel_RHF;
    char oLabel_RHF;
    size_t nO_RHF;
    size_t nV_RHF;
    size_t nOVshift_RHF;

    // Amplitudes
    TArray &T1_RHF;
    TArray &T2_RHF;

    //Intermediates
    // [Gauss:1995:3561] Table III
    // [Asthana:2019:4102] Appendix
    TArray &tau_RHF;
    TArray &F_ae_RHF;
    TArray &F_mi_RHF;
    TArray &F_me_RHF;
    TArray &W_mnij_RHF;
    TArray &W_abef_RHF;
    TArray &W_mbej_RHF_baab;
    TArray &W_mbej_RHF_baba;
    TArray &W_mnie_RHF;
    TArray &W_amef_RHF;
    TArray &W_mbij_RHF;
    TArray &W_abei_RHF;

    // Lambda intermediates:[Gauss:1995:3561] Table III (c)
    TArray &G_ae_RHF;
    TArray &G_mi_RHF;
    TArray &D_ai_RHF;
    TArray &D_abij_RHF;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij_RHF;
    TArray Rho_ab_RHF;
    TArray Rho_ia_RHF;
    TArray Rho_ai_RHF;

  public:


    EOMRCCSD(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediatesRref,
            const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMRCCSD();

  public:
    size_t toCompoundS_RHF(size_t a, size_t i) const;
    size_t toCompoundD_RHF(size_t a, size_t b, size_t i, size_t j) const;

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);

    void fillGuess(MatsT *guess_vec, size_t n_vec)const;

    void formF_ae();

    void formF_mi();

    void formW_mnij();

    void formW_abef();

    void formW_mbej();

    void formW_mnie();

    void formW_amef();

    void formW_mbij();

    void formW_abei();

    void formR1_tilde_RHF(const TArray &R1, const TArray &R2, TArray &tildeR1) const;

    void formR2_tilde_RHF(const TArray &R1, const TArray &R2, TArray &tildeR2) const;

    //Lambda iterations
    virtual void initializeLambda();
    void updateG_ae_RHF(const TArray &L2, TArray &G_ae) const;
    void updateG_mi_RHF(const TArray &L2, TArray &G_mi) const;
    void formL1_tilde_RHF(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL1) const;
    void formL2_tilde_RHF(const TArray &L1, const TArray &L2, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL2) const;
    virtual void runLambda();

    virtual void initializeEOMCC();
    virtual void formEOMIntermediates();

    // for full_diagonalize  algorithm
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;

    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec);


    // Density functions
    virtual void initializeDensity();
    virtual dcomplex calcOscillatorStrength(size_t i) override;

    std::array<MatsT, 3> calcTransitionDipole(const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0,const TArray& l1, const TArray& l2, bool isSame = false);
    std::array<MatsT, 3> calcGround2ExcitedTransitionDipole(size_t i);
    std::array<MatsT, 3> calcExcited2GroundTransitionDipole(size_t i);
    std::array<MatsT, 3> calcExcited2ExcitedTransitionDipole(size_t i, size_t j);
    virtual std::array<MatsT, 3> calcGroundDipole() override;

    // RHF implementation
    void buildDensity_RHF(const TArray& t1, const TArray& t2,  const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ij_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ab_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);
    void formRho_ia_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2, bool isSame = false);
    void formRho_ai_RHF(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const MatsT l0, const TArray& l1, const TArray& l2);

  };

  template <typename MatsT>
  class CVSEOMCCSD : public EOMCCBase <MatsT>{
    using TArray = TA::TArray<MatsT>;

  protected:
    char vLabel_;
    char oLabel_;
    char cLabel_;
    char hLabel_;
    char lLabel_;
    char rLabel_;
    size_t nOVshift_;
    size_t nO2shift_;
    size_t nV2shift_;
    size_t nCVSOCore_;
    size_t nCVSOValance_;
    size_t nCVSOActive_;
    size_t nCVSVContinuum_;
    size_t nCVSVValance_;
    size_t nCVSVActive_;
    size_t CVSoutOfBound_;

    std::vector<std::vector<size_t>> CVSabIndices_;
    std::vector<std::vector<size_t>> CVSijIndices_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;
    
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

    struct block_labels {
      std::vector<std::pair<size_t, size_t>> cc;
      std::vector<std::pair<size_t, size_t>> ch;
      std::vector<std::pair<size_t, size_t>> cr;
      std::vector<std::pair<size_t, size_t>> hc;
      std::vector<std::pair<size_t, size_t>> hh;
      std::vector<std::pair<size_t, size_t>> hr;
      std::vector<std::pair<size_t, size_t>> oc;
      std::vector<std::pair<size_t, size_t>> oh;
      std::vector<std::pair<size_t, size_t>> rc;
      std::vector<std::pair<size_t, size_t>> rr;
      std::vector<std::pair<size_t, size_t>> co;
      std::vector<std::pair<size_t, size_t>> ho;
      std::vector<std::pair<size_t, size_t>> rh;
      std::vector<std::pair<size_t, size_t>> cccc;
      std::vector<std::pair<size_t, size_t>> ccch;
      std::vector<std::pair<size_t, size_t>> cchh;
      std::vector<std::pair<size_t, size_t>> cccr;
      std::vector<std::pair<size_t, size_t>> chcc;
      std::vector<std::pair<size_t, size_t>> chch;
      std::vector<std::pair<size_t, size_t>> chcr;
      std::vector<std::pair<size_t, size_t>> crcc;
      std::vector<std::pair<size_t, size_t>> crch;
      std::vector<std::pair<size_t, size_t>> crrc;
      std::vector<std::pair<size_t, size_t>> crrh;
      std::vector<std::pair<size_t, size_t>> hrrc;
      std::vector<std::pair<size_t, size_t>> hrrh;
      std::vector<std::pair<size_t, size_t>> occr;
      std::vector<std::pair<size_t, size_t>> ochr;
      std::vector<std::pair<size_t, size_t>> rcrr;
      std::vector<std::pair<size_t, size_t>> rhrr;
      std::vector<std::pair<size_t, size_t>> rrcc;
      std::vector<std::pair<size_t, size_t>> rrch;
      std::vector<std::pair<size_t, size_t>> rrco;
      std::vector<std::pair<size_t, size_t>> rrhh;
      std::vector<std::pair<size_t, size_t>> rrho;
      std::vector<std::pair<size_t, size_t>> rroc;
      std::vector<std::pair<size_t, size_t>> rroh;
      std::vector<std::pair<size_t, size_t>> rrrc;
      std::vector<std::pair<size_t, size_t>> rrrh;
      std::vector<std::pair<size_t, size_t>> rrrr;
    };
    block_labels b;
    std::vector<MatsT> scalars_;
    TArray Id_oo;
    TArray Id_oooo;
    std::map<std::string, TArray> &reuse_tmps_;
    std::map<std::string, TArray> &tmps_;


  public:

    CVSEOMCCSD(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);
    ~CVSEOMCCSD();
    void assignCVSIndices();
     
    virtual size_t nOVShift() const {return nOVshift_;}
    size_t toCompoundS(size_t a, size_t i) const;
    size_t toCompoundD(size_t a, size_t b, size_t i, size_t j) const;
    size_t toCompoundSS(size_t a, size_t i, size_t b, size_t j, size_t ldH) const;
    //std::pair<size_t, double> toCompoundSD(size_t e, size_t m,
    //                                          size_t a, size_t b, size_t i, size_t j, size_t ldH) const;
    //std::pair<size_t, double> toCompoundDS(size_t a, size_t b, size_t i, size_t j,
    //                                          size_t e, size_t m, size_t ldH) const;
    //std::pair<size_t, double> toCompoundDD(size_t a, size_t b, size_t i, size_t j,
    //                                          size_t c, size_t d, size_t k, size_t l, size_t ldH) const;

    // for full_diagonalize  algorithm
    //virtual cqmatrix::Matrix<MatsT> buildHbar(bool includeGroundState) const;
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    virtual bool isInBound(size_t idx) const { return idx < CVSoutOfBound_; }
    virtual size_t oneBodySize() const { return nOVshift_; }

    //Lambda iterations
    virtual void initializeLambda();
    void updateG_ae(const TArray &L2, const TArray &L2Val, TArray &G_ae) const;
    void updateG_mi(const TArray &L2, const TArray &L2Val, TArray &G_mi) const;
    void formL1_tilde(const TArray &L1, const TArray &L2, const TArray &L2Val, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL1) const;
    void formL2_tilde(const TArray &L1, const TArray &L2, const TArray &L2Val, const TArray &G_ae, const TArray &G_mi,
                      TArray &tildeL2) const;
    void formL2Val_tilde(const TArray &L1, const TArray &L2, const TArray &L2Val, const TArray &G_ae, const TArray &G_mi, 
                                                TArray &tildeL2Val) const;
    virtual void runLambda();
    void finalizeLambda();
    void cleanNonCVSIntermediates();
//
    virtual void initializeEOMCC();
    virtual void formEOMIntermediates();

    void buildRightZeroBody(size_t nVec);
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;
  protected:
    void formF_ae();
    void formF_mi();
    void formW_mnij();
    void formW_abef();
    void formW_mbej();
    void formW_mnie();
    void formW_amef();
    void formW_mbij();
    void formW_abei();
    void formR1_tilde(const TArray &R1, const TArray &R2, const TArray &R2val, TArray &tildeR1) const;
    void formR2_tilde(const TArray &R1, const TArray &R2, const TArray &R2val, TArray &tildeR2) const;
    void formR2Val_tilde(const TArray &R1, const TArray &R2, const TArray &R2val,TArray &tildeR2val) const;
    void buildSigmaRight(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV) const;
    void buildSigmaLeft(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV) const;
    void buildCVSLeftIntermediates();
  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder() ;
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType) ;
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);

    void fillGuess(MatsT *guess_vec, size_t n_vec)const;
    protected: 
    
    //  
//    virtual void buildRightZeroBody(size_t nVec){}

    // Density functions
    virtual void initializeDensity();
    virtual dcomplex calcOscillatorStrength(size_t i) override;
    void buildDensity(const TArray& t1, const TArray& t2,  const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val, bool isSame = false);
    void formRho_ij(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val, bool isSame = false);
    void formRho_ab(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val);
    void formRho_ia(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val, bool isSame = false);
    void formRho_ai(const TArray& t1, const TArray& t2, const MatsT r0, const TArray& r1, const TArray& r2, const TArray& r2val, const MatsT l0, const TArray& l1, const TArray& l2, const TArray& l2val);
  
  };
  template <typename MatsT>
  class EOMEA : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nV2shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

    std::map<std::string, TArray> &reuse_tmps_;
    std::map<std::string, TArray> &tmps_;

    MatsT ccsd_energy = 0.0;
  public:


    EOMEA(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMEA();

  protected:
    std::vector<std::vector<size_t>> abIndices_;
    size_t outOfBound_;
    
  public:
    inline double signD(size_t a, size_t b, size_t i) const {
      if (a == b) return 0.0;
      if (a > b) {
        return -1.0;
      } else {
        return 1.0;
      }
    }
    size_t toCompoundS(size_t i) const;
    size_t toCompoundD(size_t a, size_t i, size_t j) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nV_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);


    //Lambda iterations
    void formR_tilde(const TArray &R2, const TArray &R3, TArray &tildeR2, TArray &tildeR3) const;
    void formL_tilde(const TArray &L2, const TArray &L3, TArray &tildeL2, TArray &tildeL3) const{}
    virtual void runLambda();

    virtual void initializeEOMCC();
       
    virtual void formEOMIntermediates();

    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec) {}

    // Density functions
    virtual void initializeDensity(){}
    virtual dcomplex calcOscillatorStrength(size_t i) override { return 0;}
    void buildDensity(const TArray& t1, const TArray& t2,  const TArray& r2, const TArray& r3, const TArray& l2, const TArray& l3){}

  };
  template <typename MatsT>
  class EOMIP_2h1p : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nO2shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

//    std::map<std::string, TArray> &reuse_tmps_;
//    std::map<std::string, TArray> &tmps_;
    // Wloch05_134113: Eqs. (A7)--(A12)
    TArray &tau;
    TArray &F_mi;
    TArray &F_ae;
    TArray &F_me;
    TArray &W_mnie;
    TArray &W_mbij;
    TArray &W_mnij;
    TArray &W_mbej;
    TArray &W_amef;

    MatsT ccsd_energy = 0.0;
  public:


    EOMIP_2h1p(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMIP_2h1p();

  protected:
    std::vector<std::vector<size_t>> ijIndices_;
    size_t outOfBound_;
    
  public:
    size_t toCompoundS(size_t i) const;
    size_t toCompoundD(size_t a, size_t i, size_t j) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nO_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);

    // EOMCCSD intermediates as usual
    // Wloch05_134113
    void formF_ae();
    void formF_mi();
    void formW_mnij();
//    void formW_abef();
    void formW_mbej();
    void formW_mnie();
    void formW_amef();
    void formW_mbij();
//    void formW_abei();

    //Lambda iterations
    void formR_tilde(const TArray &R2, const TArray &R3, TArray &tildeR2, TArray &tildeR3) const;
    void formL_tilde(const TArray &L2, const TArray &L3, TArray &tildeL2, TArray &tildeL3) const{}
    virtual void runLambda();

    virtual void initializeEOMCC();
       
    virtual void formEOMIntermediates();

    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec) {}

    // Density functions
    virtual void initializeDensity(){}
    virtual dcomplex calcOscillatorStrength(size_t i) override { return 0;}
    void buildDensity(const TArray& t1, const TArray& t2,  const TArray& r2, const TArray& r3, const TArray& l2, const TArray& l3){}

  };
  template <typename MatsT>
  class EOMIP_3h2p : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nO2shift_;
    size_t nO3shift_;
    size_t nV2shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

//    std::map<std::string, TArray> &reuse_tmps_;
//    std::map<std::string, TArray> &tmps_;
    // Wloch05_134113: Eqs. (A7)--(A12)
    TArray &tau;
    TArray &F_mi;
    TArray &F_ae;
    TArray &F_me;
    TArray &W_mnie;
    TArray &W_mbij;
    TArray &W_mnij;
    TArray &W_mbej;
    TArray &W_amef;
    TArray &W_abef;
    TArray &W_abei;

    MatsT ccsd_energy = 0.0;
  public:


    EOMIP_3h2p(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMIP_3h2p();

  protected:
    std::vector<std::vector<size_t>> ijIndices_;
    std::vector<std::vector<std::vector<size_t>>> ijkIndices_;
    std::vector<std::vector<size_t>> abIndices_;
    size_t outOfBound_;
    
  public:
    size_t toCompoundS(size_t i) const;
    size_t toCompoundD(size_t a, size_t i, size_t j) const;
    size_t toCompoundT(size_t a, size_t b, size_t i, size_t j, size_t k) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nO_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);

    // EOMCCSD intermediates as usual
    // Wloch05_134113
    void formF_ae();
    void formF_mi();
    void formW_mnij();
    void formW_abef();
    void formW_mbej();
    void formW_mnie();
    void formW_amef();
    void formW_mbij();
    void formW_abei();

    //Lambda iterations
    void formR_tilde(const TArray &R1, const TArray &R2, const TArray &R3, TArray &tildeR1, TArray &tildeR2, TArray &tildeR3) const;
    void formL_tilde(const TArray &L1, const TArray &L2, const TArray &L3, TArray &tildeL1, TArray &tildeL2, TArray &tildeL3) const{}
    virtual void runLambda();

    virtual void initializeEOMCC();
       
    virtual void formEOMIntermediates();

    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec) {}

    // Density functions
    virtual void initializeDensity(){}
    virtual dcomplex calcOscillatorStrength(size_t i) override { return 0;}
    void buildDensity(const TArray& t1, const TArray& t2,  const TArray& r1, const TArray& r2, const TArray& r3, const TArray& l1, const TArray& l2, const TArray& l3){}

  };
  template <typename MatsT>
  class EOMDIP_3h1p : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nO2shift_;
    size_t nO3shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;

    TArray &F_mi;
    TArray &F_ae;
    TArray &F_me;
    TArray &tau;
    TArray &W_mnij;
    TArray &W_mbej;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

    MatsT scalar_0, scalar_1, scalar_2, scalar_3, scalar_4;
    TArray & tempPerm_oo;
    TArray Id_oo;
    TArray Id_oooo;
    std::map<std::string, TArray> &tempOps;
    std::map<std::string, TArray> &sigmaOps;

    MatsT ccsd_energy = 0.0;
  public:


    EOMDIP_3h1p(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMDIP_3h1p();

  protected:
    std::vector<std::vector<size_t>> ijIndices_;
    std::vector<std::vector<std::vector<size_t>>> ijkIndices_;
    size_t outOfBound_;
    
    void formF_ae();
    void formF_mi();
    void formW_mnij();
    void formW_mbej();

  public:
    inline double signT(size_t i, size_t j, size_t k) const ;
    inline double signD(size_t i, size_t j) const ;
    size_t toCompoundS(size_t i, size_t j) const;
    size_t toCompoundD(size_t a, size_t i, size_t j, size_t k) const;
    std::pair<size_t, double> toCompoundSS(size_t i, size_t j, size_t k, size_t l, size_t ldH) const;
    std::pair<size_t, double> toCompoundSD(size_t i, size_t j,
                                              size_t k, size_t l, size_t m, size_t a, size_t ldH) const;
    std::pair<size_t, double> toCompoundDS(size_t i, size_t j, size_t k, size_t a,
                                              size_t l, size_t m, size_t ldH) const;
    std::pair<size_t, double> toCompoundDD(size_t i, size_t j, size_t k, size_t a,
                                                      size_t l, size_t m, size_t n, size_t b, size_t ldH) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nO2shift_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);

    void fillGuess(MatsT *guess_vec, size_t n_vec)const;

    //Lambda iterations
    void formR_tilde(const TArray &R2, const TArray &R3, TArray &tildeR2, TArray &tildeR3) const;
    void formL_tilde(const TArray &L2, const TArray &L3, TArray &tildeL2, TArray &tildeL3) const;
    virtual void runLambda();

    virtual void initializeEOMCC();
  protected:
    void buildCCSDEnergy();
       
  public:
    virtual void formEOMIntermediates();

    // for full_diagonalize  algorithm
    virtual cqmatrix::Matrix<MatsT> buildHbar(bool includeGroundState) const;
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec);

    // Density functions
    virtual void initializeDensity();
    virtual dcomplex calcOscillatorStrength(size_t i) override;
    void buildDensity(const TArray& t1, const TArray& t2,  const TArray& r2, const TArray& r3, const TArray& l2, const TArray& l3);
  protected:
    void buildHbarSS(TArray & H_oooo) const;
    void buildHbarTA(TArray & H_oooo, TArray & H_ooooov,TArray & H_ooovoo,TArray & H_ooovooov) const;

  };
  template <typename MatsT>
  class EOMDIP_4h2pCCSDT : public EOMCCBase <MatsT>{
  protected:
    using TArray = TA::TArray<MatsT>;
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fullMat;

    char vLabel_;
    char oLabel_;
    size_t nO_;
    size_t nV_;
    size_t nO2shift_;
    size_t nO3shift_;
    size_t nO4shift_;
    size_t nV2shift_;

    // Amplitudes
    TArray &T1_;
    TArray &T2_;
    TArray &T3_;

    //TArray &F_mi;
    //TArray &F_ae;
    //TArray &F_me;
    //TArray &tau;
    //TArray &W_mnij;
    //TArray &W_mbej;

    //Reduced one-particle EOM-CCSD density matrices
    TArray Rho_ij;
    TArray Rho_ab;
    TArray Rho_ia;
    TArray Rho_ai;

    MatsT scalar_0, scalar_1, scalar_2, scalar_3, scalar_4;
    //TArray & tempPerm_oo;
    TArray Id_oo;
    std::map<std::string, TArray> &reused_;
    std::map<std::string, TArray> &tmps_;
    std::map<std::string, TArray> perm_tmps;

    MatsT ccsd_energy = 0.0;
  public:


    EOMDIP_4h2pCCSDT(const SafeFile &savFile,
            CCIntermediates<MatsT> &intermediates, const EOMSettings &eomSettings,
            const CoupledClusterSettings &ccSettings);

    ~EOMDIP_4h2pCCSDT();

  protected:
    std::vector<std::vector<size_t>> abIndices_;
    std::vector<std::vector<size_t>> ijIndices_;
    std::vector<std::vector<std::vector<size_t>>> ijkIndices_;
    std::vector<std::vector<std::vector<std::vector<size_t>>>> ijklIndices_;
    size_t outOfBound_;
    
    void formF_ae();
    void formF_mi();
    void formW_mnij();
    void formW_mbej();

  public:
    size_t toCompoundS(size_t i, size_t j) const;
    size_t toCompoundD(size_t a, size_t i, size_t j, size_t k) const;
    size_t toCompoundT(size_t a, size_t b, size_t i, size_t j, size_t k, size_t l) const;

    virtual bool isInBound(size_t idx) const { return idx < outOfBound_; }
    virtual size_t oneBodySize() const { return nO2shift_; }

  public:
    virtual typename Davidson<MatsT>::VecsGen_t EmptyDavidsonVectorBuilder();
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType);
    virtual typename Davidson<MatsT>::LinearTrans_t DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag);


    //Lambda iterations
    void formR_tilde(const TArray &R2, const TArray &R3, const TArray &R4, TArray &tildeR2, TArray &tildeR3, TArray &tildeR4) const;
    void formL_tilde(const TArray &L2, const TArray &L3, const TArray &L4, TArray &tildeL2, TArray &tildeL3, TArray &tildeL4) const{}
    virtual void runLambda();

    virtual void initializeEOMCC();
  protected:
    void buildCCSDEnergy();
       
  public:
    virtual void formEOMIntermediates();

    // for full_diagonalize  algorithm
    virtual cqmatrix::Matrix<MatsT> buildHbar(bool includeGroundState) const { return cqmatrix::Matrix<MatsT>(1);}
    virtual void buildDiag(MatsT * diag, const std::vector<double> &eps) const override;


    //Helper functions for diagonalization. Davidson/GPLHR
    virtual void buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const;

    virtual void buildRightZeroBody(size_t nVec) {}

    // Density functions
    virtual void initializeDensity() {}
    virtual dcomplex calcOscillatorStrength(size_t i) override { return 0;}
    void buildDensity(const TArray& t1, const TArray& t2,  const TArray& r2, const TArray& r3, const TArray& l2, const TArray& l3) {}
  protected:
    void buildHbarSS(TArray & H_oooo) const {}
    void buildHbarTA(TArray & H_oooo, TArray & H_ooooov,TArray & H_ooovoo,TArray & H_ooovooov) const {}

  };

  template<typename MatsT, typename IntsT>
  void runCoupledCluster(JobType jobType, Molecule &mol, std::shared_ptr<SingleSlater<MatsT,IntsT>> ccref,
                         std::shared_ptr<IntegralsBase> aoints,
                         SafeFile &rstFile, CQInputFile &input, std::ostream &output);



}; // namespace ChronusQ
#endif
