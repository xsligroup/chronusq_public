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
#include <wavefunction.hpp>
#include <singleslater/base.hpp>
#include <singleslater/rdm.hpp>
#include <particleintegrals/twopints.hpp>
#include <matrix.hpp>
#include <cubegen.hpp>
#include <orthogonalization.hpp>
#include <orbitalmodifier.hpp>
#include <functional>

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {


  // This struct is for each atom in population analysis
  struct popAtom{

    std::array<double,4> total; // length 4 because S, Z, Y, X
    std::array<std::vector<double>,4> angMom; // s, p, d, ....

    void init(size_t maxL) {
        total.fill(0.0);
        for (auto &v : angMom)
            v.assign(maxL, 0.0);
    }

    void add(DENSITY_TYPE comp, double val) {
        total[comp]   += val;
    }

    // Different prefacs are needed since Tot is sometimes Z-q
    void addL(DENSITY_TYPE comp, size_t lVal, double prefacTot, double prefacL, double val) {
        total[comp]   += prefacTot*val;
        angMom[comp][lVal] += prefacL*val;
    }

    double get(DENSITY_TYPE comp) const {
        return total[comp];
    }

    double getL(DENSITY_TYPE comp, size_t lVal) const {
        return angMom[comp][lVal];
    }

  };

  // Declaration of CoreH and Fock builders.
  template <typename MatsT, typename IntsT>
  class CoreHBuilder;
  template <typename MatsT, typename IntsT>
  class FockBuilder;
  template <typename MatsT, typename IntsT>
  class MOIntsTransformer;

  // Use SFINAE to ensure RT functions only get instantiated when MatsT is dcomplex
  template <typename M>
  using enable_if_dcomplex = typename std::enable_if<std::is_same<M, dcomplex>::value, int>::type;
  using RTFockFormation = std::function<void(EMPerturbation&, bool)>;

  /**
   *  \brief The SingleSlater class. The typed abstract interface for all
   *  classes for which the wave function is described by a single slater
   *  determinant (HF, KS, PHF, etc).
   *
   *  Adds knowledge of storage type to SingleSlaterBase
   *
   *  Specializes the WaveFunction class of the same type
   */ 
  template <typename MatsT, typename IntsT>
  class SingleSlater : public SingleSlaterBase, public WaveFunction<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class SingleSlater;

  protected:

    // Useful typedefs
    typedef MatsT*                    oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

    //BasisSet &basisSet_; ///< BasisSet for the GTO basis definition

  private:
  public:

    typedef MatsT value_type;
    typedef IntsT ints_type;

    bool isX2CReference() const override {
      return this->nC == 2 && this->aoints_ != nullptr &&
             this->aoints_->options_.x2cType != X2C_TYPE::OFF;
    }

    //ORTHO_TYPE            orthoType; ///< Orthogonalization scheme

    // Operator storage
    std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> moCoefficients; ///< List of populated MO coefficient matrices
    std::vector<double*> moEigenvalues; ///< List of populated MO eigenvalues
    virtual void initializeSCF() override; ///< Initialize SCF, populate MO coefficients and eigenvalues

    // AO Fock Matrix
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> fockMatrix; ///< List of populated AO Fock matrices
    std::vector<cqmatrix::Matrix<MatsT>> fockMO;     ///< Fock matrix in the MO basis

    // Orthonormal Fock
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> fockMatrixOrtho; ///< List of populated orthonormal Fock matrices

    // Coulomb (J[D])
    std::shared_ptr<cqmatrix::Matrix<MatsT>> coulombMatrix; ///< scalar Coulomb Matrix

    // Exchange (K[D])
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> exchangeMatrix; ///< List of populated exact (HF) exchange matrices

    // Two-electron Hamiltonian (G[D])
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> twoeH; ///< List of populated HF perturbation tensors

    // Orthonormal density
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDMOrtho; ///< List of populated orthonormal 1PDM matrices
    std::vector<cqmatrix::Matrix<MatsT>> onePDMAlphaBetaOrtho; ///< List of populated orthonormal 1PDM matrices
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> deltaOnePDM; ///< Change in density for incremental Fock Build

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH; ///< Core Hamiltonian (scalar and magnetization)
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreHPerturbed; ///< Perturbed Core Hamiltonian (scalar and magnetization)

    // Algorithm Abstractions
    std::shared_ptr<TPIContractions<MatsT,IntsT>> TPI; ///< TPIContractions
    std::shared_ptr<CoreHBuilder<MatsT,IntsT>> coreHBuilder; ///< Builder for CoreH
    std::shared_ptr<FockBuilder<MatsT,IntsT>> fockBuilder;  ///< Builder for Fock
    std::shared_ptr<Orthogonalization<MatsT>> orthoSpinor;  ///< Orthogonalization functions for spinor basis
    std::shared_ptr<Orthogonalization<MatsT>> orthoAB;      ///< Orthogonalization functions alpha/beta basis
    std::shared_ptr<OrbitalModifier<MatsT>> orbitalModifier;  ///< SCF/RT Abstraction Object
    
    std::shared_ptr<cqmatrix::Matrix<MatsT>> tau;      ///< tau matrix for traveling proton basis
    std::shared_ptr<cqmatrix::Matrix<MatsT>> tauOrtho; ///< orthonormal tau matrix for traveling proton basis

    // Method specific property storage
    std::vector<popAtom> mullikenCharges;
    std::vector<popAtom> lowdinCharges;

    // Whether the current Density and Coefficients represent the same wavefunction
    // If so, then RI-K contractions can be done using Coefficients for better performance
    bool denEqCoeff_ = false;

    // Temporary structure to hold x2c picture-changed dipole matrices
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>> pchgDipole_;

    // 4C gaunt
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> gaunttwoeH; ///< contributions from gaunt to twoeH
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> gauntexchangeMatrix; ///< contributions from gaunt to exchangeMatrix
    // 4C gauge 
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> gaugetwoeH; ///< contributions from gaunt to twoeH
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> gaugeexchangeMatrix; ///< contributions from gaunt to exchangeMatrix

    // Density Gradient (vector of length NAtoms*3 for each R)
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> onePDMGrad;
    // Energy weighted density matrix W
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> W;

    std::vector<size_t> ownedAtomIndices; ///< Indices of atoms owned by this SingleSlater TODOAL: maybe move this somewhere else?

    // For modified of SCF calculations, different ways to build the density matrix
    std::shared_ptr<RDMBuilderBase<MatsT,IntsT>> RDMBuilder = nullptr;

    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handles the BasisSet, etc)
     *  \param [in] args Parameter pack for the remaining parameters of the
     *                   WaveFunction constructor. See include/wavefunction.hpp
     *                   for details. 
     */ 
    template <typename... Args>
    SingleSlater(MPI_Comm c, Molecule &mol, BasisSet &basis,
                 std::shared_ptr<Integrals<IntsT>> aoi, Args... args) :
      SingleSlaterBase(c,mol,basis,args...),
      WaveFunctionBase(c,mol,basis,args...),
      QuantumBase(c,args...),
      WaveFunction<MatsT,IntsT>(c,mol,basis,aoi,args...)
      //, basisSet_(basis)
      //, coreType(NON_RELATIVISTIC), orthoType(LOWDIN)
    {
      // Allocate SingleSlater Object
      alloc();

      // Determine Real/Complex part of method string
      if(std::is_same<MatsT,double>::value) {
        refLongName_  = "Real ";
        refShortName_ = "R-";
      } else {
        refLongName_  = "Complex ";
        refShortName_ = "C-";
      }

      // Initialize temporary container for x2c picture change dipole
      pchgDipole_.reserve(3);
      for (size_t i = 1; i <= 3; i++) {
        pchgDipole_.emplace_back(nullptr);
      }

    }; // SingleSlater constructor

    // See include/singleslater/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename MatsU> 
      SingleSlater(const SingleSlater<MatsU,IntsT> &, int dummy = 0);
    template <typename MatsU> 
      SingleSlater(SingleSlater<MatsU,IntsT> &&     , int dummy = 0);

    // Same type
    SingleSlater(const SingleSlater<MatsT,IntsT> &);
    SingleSlater(SingleSlater<MatsT,IntsT> &&);     

    /**
     *  Destructor.
     *
     *  Destructs a SingleSlater object
     */ 
    ~SingleSlater() { dealloc(); }



    // Public Member functions

    //BasisSet& basisSet() { return basisSet_; }
      
      

    // Deallocation (see include/singleslater/impl.hpp for docs)
    void alloc();
    void dealloc();


    // Declarations from QuantumBase 
    // (see include/singleslater/quantum.hpp for docs)
    void constructRDMBuilder();
    void formDensity() override;

    using QuantumBase::computeEnergy;
    void computeEnergy() override;
    void computeMultipole(EMPerturbation &, const std::vector<PROPERTY> &properties = {}) override;
    void compute4CDipole(EMPerturbation &);
    void computeFockX2CDipole(EMPerturbation &);
    void computeSpin() override;
    void computeSpinAndAngularProperties() override;
    void computeSpinAndAngularProperties(const InCore4indexTPI<MatsT>* twoRDM);
    void computeSpinAndAngularProperties(const InCore4indexTPI<MatsT>* twoRDM, size_t activeEndOff);
    void printAngularProperties(std::ostream&, bool withBanner = true) override;

    // Types and helpers for angular orbital expectation evaluation
    struct AngularOrbitalRow {
      double energy;
      double occupation;
      double sz;
      double lz;
      double jz;
      double s2;
      double l2;
      double j2;
    };

    struct AngularResults {
      std::vector<AngularOrbitalRow> rows;
      double operatorBuildTime = 0.0;
      double expectationTime = 0.0;
      double totalTime = 0.0;
    };

    // Compute orbital-level angular expectation rows and timing metrics
    AngularResults computeAngularOrbitalRows(bool includeBeta = true);
    void computeOrbitalProps() override;
    void computeOrbitalRDFs() override;
    void printOrbitalEnergies() override;
    virtual std::vector<double> getEnergySummary() override;

    // Compute various core Hamitlonian
    void formCoreH(EMPerturbation&, bool) override; // Compute the CH
    void computeOrtho();  // Evaluate orthonormalization transformations
    void computeOrthoGrad(); // Evaluate gradient of orthonormalization
    
    // RT functions
    void computeTau();
    void checkIdempotency(std::string system="Electronic");
    void mcWeenyPurification();
    // RT functions that require complex matrix types
    template <typename M = MatsT, enable_if_dcomplex<M> = 0>
    void addTauToFock();
    template <typename M = MatsT, enable_if_dcomplex<M> = 0>
    void RK4Propagation(bool, double, bool, EMPerturbation&, EMPerturbation&, const RTFockFormation& = {});
    template <typename M = MatsT, enable_if_dcomplex<M> = 0>
    void unitaryPropagation(bool, double, bool, EMPerturbation&, const RTFockFormation& = {});
    template <typename M = MatsT, enable_if_dcomplex<M> = 0>
    cqmatrix::PauliSpinorMatrices<MatsT> getTimeDerDen(bool);

    // Method specific properties
    void populationAnalysis();
    void methodSpecificProperties() override {
      populationAnalysis();
    }

    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock (EMPerturbation &, bool increment = false, double xHFX = 1.) override;
    void formFock(EMPerturbation& pert) { formFock(pert,false,1.);};

    // Get the total gradient
    virtual std::vector<double> getGrad(EMPerturbation&, bool equil,
      bool saveInts, double xHFX = 1.) override;

    void formEWDM(bool equil = false);
    void formEWDM_impl(bool equil, const cqmatrix::Matrix<MatsT>& F, const cqmatrix::Matrix<MatsT>& P, cqmatrix::Matrix<MatsT>& W);

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess(EMPerturbation &, const SingleSlaterOptions&);
    void CoreGuess();
    void SADGuess(SingleSlaterOptions);
    void SCFGuess(SingleSlaterOptions);
    void TightGuess();
    void NEOTightParticleGuess();
    void RandomGuess();
    void ReadGuessMO( const std::shared_ptr<BasisSet> guessBasis, std::string prefix );
    void ReadGuess1PDM( const std::shared_ptr<BasisSet> guessBasis, std::string prefix = "" );
    void FchkGuessMO();
    void NEOConvergeClassicalGuess(EMPerturbation &, const SingleSlaterOptions&);
    void computeNaturalOrbitals();
    void getNewOrbitals();

    // ReadGuess1PDM functions
    void readSameTypeDenBin(std::string prefix = "");
    void readDiffTypeDenBin(std::string binName, const std::shared_ptr<BasisSet> guessBasisSet, std::string prefix = "");
    template <typename ScrMatsT>
    void getScr1PDM(SafeFile &, const std::shared_ptr<BasisSet>, std::string prefix = "");
    template <typename ScrMatsT>
    void getScr1PDM(SafeFile &);
    

    // ReadGuessMO functions
    void readSameTypeMOBin(std::string prefix = "");
    void readDiffTypeMOBin(std::string binName, const std::shared_ptr<BasisSet> guessBasisSet, std::string prefix = "" );
    template <typename ScrMatsT>
    void getScrMO(SafeFile &, const std::shared_ptr<BasisSet>, std::string prefix );
    template <typename ScrMatsT>
    void convert1CRto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&);
    template <typename ScrMatsT>
    void convert1CUto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&, size_t nA, size_t nB);
    template <typename ScrMatsT>
    void convert1CUto2CU_sameType(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<ScrMatsT>>&, size_t nA, size_t nB);
    template <typename ScrMatsT>
    void convert1CRto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&, SafeFile &);
    template <typename ScrMatsT>
    void convert1CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&, SafeFile &);
    template <typename ScrMatsT>
    void convert2CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&, SafeFile &);

    // MO Projection functions
    std::tuple<std::shared_ptr<cqmatrix::Matrix<MatsT>>, std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>> projectMatrix( const Molecule& mol, const BasisSet& fromBasis,
      const BasisSet& toBasis, const std::shared_ptr<cqmatrix::Matrix<MatsT>>& fromMatrix, bool doFull = true );
    std::shared_ptr<cqmatrix::Matrix<MatsT>> getProjectionMatrix( const cqmatrix::Matrix<MatsT>& overlap21, const cqmatrix::Matrix<MatsT>& overlap22);
    std::shared_ptr<cqmatrix::Matrix<MatsT>> projectFullMatrix( const std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>& projMat, const std::shared_ptr<cqmatrix::Matrix<MatsT>>&fromMatrix );
    std::shared_ptr<cqmatrix::Matrix<MatsT>> projectHalfMatrix( const std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>& projMat, const std::shared_ptr<cqmatrix::Matrix<MatsT>>&fromMatrix );


    // Fchk-related functions
    std::vector<int> fchkToCQMO();
    std::unordered_map<int,std::vector<int>> returnAngReorder();
    void reorderAngMO(std::vector<int> sl, MatsT* tmo, int sp);
    void reorderSpinMO();

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ao2orthoMOs();
    void ao2orthoDen();
    void ortho2aoDen();
    void ortho2aoMOs();
    void orthoAOMO();

    // Post-processing functions
    virtual bool secondSCF();
    void runCube(std::vector<std::shared_ptr<CubeGen>> cu, std::string prefix, std::shared_ptr<Molecule> mol) override;

    // SCF Specific Functions
    inline virtual double getTotalEnergy() { return this->totalEnergy; };
    virtual void printProperties();
    virtual std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getOnePDM();
    virtual std::vector<cqmatrix::Matrix<MatsT>> getOnePDMOrtho();
    virtual std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getFock();
    virtual void setOnePDMOrtho(cqmatrix::Matrix<MatsT>*);
    virtual void setOnePDMAO(cqmatrix::Matrix<MatsT>*);
    virtual std::vector<std::shared_ptr<Orthogonalization<MatsT>>> getOrtho();
    virtual void runSCF(EMPerturbation&) override;
    virtual std::vector<NRRotOptions> buildRotOpt();

    // Converts 1C SSbase into 2C GHF SSbase
    std::shared_ptr<SingleSlaterBase> convert1CSSToGHFSS(
    SingleSlaterOptions &ssOptions,
    std::ostream &output,
    CQInputFile &input,
    EMPerturbation &emPert) override;
    void MOSpinBlockBySpace(size_t nActEA, size_t nActOA, 
      size_t nActEB, size_t nActOB) override;

    // Misc procedural
    void diagOrthoFock();
    void diagAOFock();
    virtual void saveCurrentState(bool saveMO = true, std::string prefix = "") override;

    // Stability and reopt
    virtual std::pair<double,MatsT*> getStab() = 0;
    bool checkStability();
    virtual MatsT* getNRCoeffs() { return nullptr;};

    // Print functions
    void printFock(std::ostream& ) override   ;
    void print1PDMOrtho(std::ostream&) override ;
    void printGD(std::ostream&)   override    ;
    void printJ(std::ostream&)  override      ;
    void printK(std::ostream&)   override     ;
    void printMiscProperties(std::ostream&) override;
    void printEPS(std::ostream&) override;
    void printMOInfo(std::ostream&, size_t a = 0) override;
    virtual void printFockTimings(std::ostream&) override;

    // Method to produce a test on integral transformation 
#ifdef TEST_MOINTSTRANSFORMER
    void MOIntsTransformationTest(EMPerturbation &);
#endif    
    std::shared_ptr<MOIntsTransformer<MatsT, IntsT>> generateMOIntsTransformer(TPI_TRANSFORMATION_ALG alg = TPI_TRANSFORMATION_ALG::DIRECT_N6);

    // MO Transformations
    void MOFOCK();
  
    // Set the flag that indicates whether 
    // the current Density and Coefficients represent the same wavefucntion
    void setDenEqCoeff(bool val);

    // Obtain localized MOs using Cholesky decomposition on the density matrix
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getCholeskyMOs();

    // Project a AO density onto a MO basis
    cqmatrix::Matrix<MatsT> generateMODensity(const cqmatrix::Matrix<MatsT>&, const cqmatrix::Matrix<MatsT>&);
    // Print the occupation of orbitals based on MO density 
    void printOrbitalPopulation(std::ostream&);

    // Pointer converter
    template <typename MatsU>
    static std::shared_ptr<SingleSlater<MatsU,IntsT>>
    convert(const std::shared_ptr<SingleSlater<MatsT,IntsT>>&);

  }; // class SingleSlater

}; // namespace ChronusQ

// Include declaration of CoreHBuilder and FockBuilder
#include <corehbuilder.hpp>
#include <fockbuilder.hpp>
#include <mointstransformer.hpp>

// Include headers for specializations of SingleSlater
#include <singleslater/hartreefock.hpp> // HF specialization
#include <singleslater/kohnsham.hpp>    // KS specialization

