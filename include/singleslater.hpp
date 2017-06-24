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
#include <particleintegrals/twopints.hpp>
#include <matrix.hpp>
#include <cubegen.hpp>
#include <orthogonalization.hpp>
#include <orbitalmodifier.hpp>

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {

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

    //BasisSet &basisSet_; ///< BasisSet for the GTO basis defintion

  private:
  public:

    typedef MatsT value_type;
    typedef IntsT ints_type;

    //ORTHO_TYPE            orthoType; ///< Orthogonalization scheme

    // Operator storage
    std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>> moCoefficients; ///< List of populated MO coefficient matricies
    std::vector<double*> moEigenvalues; ///< List of populated MO eigenvalues
    virtual void initializeSCF() override; ///< Initialize SCF, populate MO coefficients and eigenvalues

    // AO Fock Matrix
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> fockMatrix; ///< List of populated AO Fock matricies
    std::vector<cqmatrix::Matrix<MatsT>> fockMO;     ///< Fock matrix in the MO basis

    // Orthonormal Fock
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> fockMatrixOrtho; ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    std::shared_ptr<cqmatrix::Matrix<MatsT>> coulombMatrix; ///< scalar Coulomb Matrix

    // Exchange (K[D])
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> exchangeMatrix; ///< List of populated exact (HF) exchange matricies

    // Two-electron Hamiltonian (G[D])
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> twoeH; ///< List of populated HF perturbation tensors

    // Orthonormal density
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> onePDMOrtho; ///< List of populated orthonormal 1PDM matricies
    std::vector<cqmatrix::Matrix<MatsT>> onePDMAlphaBetaOrtho; ///< List of populated orthonormal 1PDM matricies
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

    // Method specific propery storage
    std::vector<double> mullikenCharges;
    std::vector<double> lowdinCharges;

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

    // Density Gradient (vector of N*atoms*3 for each R)
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> onePDMGrad;

    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
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
    void formDensity() override;

    using QuantumBase::computeEnergy;
    void computeEnergy() override;
    void computeMultipole(EMPerturbation &) override;
    void compute4CDipole(EMPerturbation &);
    void computeFockX2CDipole(EMPerturbation &);
    void computeSpin() override;
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
    void RK4Propagation(bool, double, bool, EMPerturbation&, EMPerturbation&);
    template <typename M = MatsT, enable_if_dcomplex<M> = 0>
    void unitaryPropagation(bool, double, bool, EMPerturbation&);
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

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess(const SingleSlaterOptions&) override;
    void CoreGuess();
    void SADGuess(const SingleSlaterOptions&);
    void TightGuess();
    void RandomGuess();
    void ReadGuessMO();
    void ReadGuess1PDM();
    void FchkGuessMO();
    void NEOTightProtonGuess();
    void NEOConvergeClassicalGuess(const SingleSlaterOptions&);
    void computeNaturalOrbitals();
    void getNewOrbitals();

    // ReadGuess1PDM functions
    void readSameTypeDenBin();
    void readDiffTypeDenBin(std::string binName);
    template <typename ScrMatsT>
    void getScr1PDM(SafeFile &);

    // ReadGuessMO functions
    void readSameTypeMOBin();
    void readDiffTypeMOBin(std::string binName);
    template <typename ScrMatsT>
    void getScrMO(SafeFile &);
    template <typename ScrMatsT>
    void convert1CRto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&);
    template <typename ScrMatsT>
    void convert1CUto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&);
    template <typename ScrMatsT>
    void convert1CRto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&);
    template <typename ScrMatsT>
    void convert1CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&);
    template <typename ScrMatsT>
    void convert2CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>&, std::vector<cqmatrix::Matrix<MatsT>>&, SafeFile &);

    // Fchk-related functions
    std::vector<int> fchkToCQMO();
    std::unordered_map<int,std::vector<int>> returnAngReorder();
    void reorderAngMO(std::vector<int> sl, MatsT* tmo, int sp);
    void reorderSpinMO();

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ao2orthoMOs();
    void ao2orthoDen();
    virtual void ortho2aoDen();
    void ortho2aoMOs();
    void orthoAOMO();

    // Post-processing functions
    void runCube(std::vector<std::shared_ptr<CubeGen>> cu, EMPerturbation &emPert, std::string prefix, std::shared_ptr<Molecule> mol) override;

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

    // Misc procedural
    void diagOrthoFock();
    void diagAOFock();
    virtual void saveCurrentState(bool saveMO = true) override;

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
    std::shared_ptr<MOIntsTransformer<MatsT, IntsT>> generateMOIntsTransformer();

    // MO Transformations
    void MOFOCK();
  
    // Set the flag that indicates whether 
    // the current Density and Coefficients represent the same wavefucntion
    void setDenEqCoeff(bool val);

    // Project a AO density onto a MO basis
    cqmatrix::Matrix<MatsT> generateMODensity(const cqmatrix::Matrix<MatsT>&, const cqmatrix::Matrix<MatsT>&);
    // Print the occupation of orbitals based on MO density 
    void printOrbitalPopulation(std::ostream&);

    // Pointer convertor
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

