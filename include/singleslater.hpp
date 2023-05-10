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
#include <orthogonalization.hpp>
#include <modifyorbitals.hpp>

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

    // AO Fock Matrix
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> fockMatrix; ///< List of populated AO Fock matricies
    std::vector<SquareMatrix<MatsT>> fockMO;     ///< Fock matrix in the MO basis

    // Orthonormal Fock
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> fockMatrixOrtho; ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    std::shared_ptr<SquareMatrix<MatsT>> coulombMatrix; ///< scalar Coulomb Matrix

    // Exchange (K[D])
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> exchangeMatrix; ///< List of populated exact (HF) exchange matricies

    // Two-electron Hamiltonian (G[D])
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> twoeH; ///< List of populated HF perturbation tensors

    // Orthonormal density
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> onePDMOrtho; ///< List of populated orthonormal 1PDM matricies
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> deltaOnePDM; ///< Change in density for incremental Fock Build

    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreH; ///< Core Hamiltonian (scalar and magnetization)
    std::shared_ptr<PauliSpinorSquareMatrices<MatsT>> coreHPerturbed; ///< Perturbed Core Hamiltonian (scalar and magnetization)

    // Algorithm Abstractions
    std::shared_ptr<TPIContractions<MatsT,IntsT>> TPI; ///< TPIContractions
    std::shared_ptr<CoreHBuilder<MatsT,IntsT>> coreHBuilder; ///< Builder for CoreH
    std::shared_ptr<FockBuilder<MatsT,IntsT>> fockBuilder;  ///< Builder for Fock
    std::shared_ptr<Orthogonalization<MatsT>> orthoSpinor;  ///< Orthogonalization functions for spinor basis
    std::shared_ptr<Orthogonalization<MatsT>> orthoAB;      ///< Orthogonalization functions alpha/beta basis
    std::shared_ptr<ModifyOrbitals<MatsT>> modifyOrbitals;  ///< SCF/RT Abstraction Object

    // Method specific propery storage
    std::vector<double> mullikenCharges;
    std::vector<double> lowdinCharges;


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
    SingleSlater(MPI_Comm c, CQMemManager &mem, Molecule &mol, BasisSet &basis,
                 std::shared_ptr<Integrals<IntsT>> aoi, Args... args) :
      SingleSlaterBase(c,mem,mol,basis,args...),
      WaveFunctionBase(c,mem,mol,basis,args...),
      QuantumBase(c,mem,args...),
      WaveFunction<MatsT,IntsT>(c,mem,mol,basis,aoi,args...)
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
    void formDensity();
    void computeEnergy();
    void computeMultipole(EMPerturbation &);
    void computeSpin();
    virtual std::vector<double> getEnergySummary();

    // Compute various core Hamitlonian
    void formCoreH(EMPerturbation&, bool); // Compute the CH
    void computeOrtho();  // Evaluate orthonormalization transformations
    void computeOrthoGrad(); // Evaluate gradient of orthonormalization

    // Method specific properties
    void populationAnalysis();
    void methodSpecificProperties() {
      populationAnalysis();
    }

    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock(EMPerturbation &, bool increment = false, double xHFX = 1.);
    void formFock(EMPerturbation& pert) { formFock(pert,false,1.);};

    // Get the total gradient
    virtual std::vector<double> getGrad(EMPerturbation&, bool equil,
      bool saveInts);

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess(const SingleSlaterOptions&);
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

    // SCF Specific Functions
    inline virtual double getTotalEnergy() { return this->totalEnergy; };
    virtual void printProperties();
    virtual std::vector<std::shared_ptr<SquareMatrix<MatsT>>> getOnePDM();
    virtual std::vector<std::shared_ptr<SquareMatrix<MatsT>>> getFock();
    virtual std::vector<std::shared_ptr<Orthogonalization<MatsT>>> getOrtho();
    virtual void runModifyOrbitals(EMPerturbation&);
    virtual std::vector<NRRotOptions> buildRotOpt();

    // Misc procedural
    void diagOrthoFock();
    void diagAOFock();
    virtual void saveCurrentState();

    // Stability and reopt
    virtual std::pair<double,MatsT*> getStab() = 0;
    bool checkStability();
    virtual MatsT* getNRCoeffs() { return nullptr;};

    // Print functions
    void printFock(std::ostream& )    ;
    void print1PDMOrtho(std::ostream&);
    void printGD(std::ostream&)       ;
    void printJ(std::ostream&)        ;
    void printK(std::ostream&)        ;
    void printMiscProperties(std::ostream&);
    void printEPS(std::ostream&);
    void printMOInfo(std::ostream&, size_t a = 0); 
    virtual void printFockTimings(std::ostream&);

    // Method to produce a test on integral transformation 
#ifdef TEST_MOINTSTRANSFORMER
    void MOIntsTransformationTest(EMPerturbation &);
#endif    
    std::shared_ptr<MOIntsTransformer<MatsT, IntsT>> generateMOIntsTransformer();

    // MO Transformations
    void MOFOCK();

  }; // class SingleSlater

}; // namespace ChronusQ

// Include declaration of CoreHBuilder and FockBuilder
#include <corehbuilder.hpp>
#include <fockbuilder.hpp>
#include <mointstransformer.hpp>

// Include headers for specializations of SingleSlater
#include <singleslater/hartreefock.hpp> // HF specialization
#include <singleslater/kohnsham.hpp>    // KS specialization

