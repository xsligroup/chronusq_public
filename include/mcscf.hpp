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
#include <cerr.hpp>
#include <singleslater.hpp>
#include <singleslater/neoss.hpp>
#include <mcwavefunction.hpp>
#include <orbitalrotation.hpp>
//#include <particleintegrals/twopints/incore4indextpi.hpp>

//#define _DEBUG_MCSCF

namespace ChronusQ {
  
  enum CIDiagonalizationAlgorithm {
    CI_FULL_MATRIX,
    CI_DAVIDSON,
    CI_GPLHR,
    SKIP,
  }; // struct CIDiagonalizationAlgorithm
  
  // Settings 
  struct MCSCFSettings {
     
     // CISettings
     CIDiagonalizationAlgorithm ciAlg = CI_FULL_MATRIX; 
     
     // for davidson and gplhr
	 size_t maxCIIter        = 128;        
     double ciVectorConv     = 1.0e-6;    
     size_t maxDavidsonSpace = 50;
     size_t nDavidsonGuess   = 3;
     std::vector<std::pair<double, size_t>> energyRefs;

     // SCF Settings 
     bool doSCF           = false;
     bool doIVOs          = false;
     
     size_t maxSCFIter         = 0;
     double scfEnergyConv      = 1.0e-8; 
     double scfGradientConv    = 1.0e-4;
     OrbitalRotationSettings ORSettings;
  
     MCSCFSettings() {
       ORSettings.rotate_within_correlated = false;
     }
     
     MCSCFSettings(const MCSCFSettings &) = default;
     MCSCFSettings(MCSCFSettings &&) = default;

     void print(bool, size_t);
  }; // struct MCSCFSettings 
  
  template <typename MatsT, typename IntsT>
  class CISolver { 
     
  protected:
    
    CIDiagonalizationAlgorithm alg_ = CI_FULL_MATRIX;

    // for davidson and gplhr
    size_t maxIter_    = 128;       /// < Max Number of CI iteration 
    double vectorConv_ = 1.0e-6;    /// < Convergence criteria in terms of vector residue norm
    size_t maxDavidsonSpace_ = 50;  /// < Max davidson space in terms of n times of NRoots  
    size_t nDavidsonGuess_ = 3;     /// < number of guess in the intial davidson first a few iterations 
    std::vector<std::pair<double, size_t>> energyRefs_;

    void davidsonGS(size_t, size_t, MatsT *, MatsT *);
    void davidsonPC(size_t, size_t, MatsT *, MatsT *, MatsT *, dcomplex *);
  
  public:
    
    // default Constructor
    CISolver() = default;
    CISolver(CIDiagonalizationAlgorithm alg, size_t maxIter = 128,
      double vectorConv =  1.0e-6, size_t maxDSpace = 50, size_t nDGuess = 3,
      std::vector<std::pair<double, size_t>> eRefs = {}) {
      switchAlgorithm(alg, maxIter, vectorConv, maxDSpace, nDGuess, eRefs);
    };

    // typeconversion
	template <typename MatsU>
    CISolver(const CISolver<MatsU,IntsT> &);
	
	template <typename MatsU>
    CISolver(CISolver<MatsU,IntsT> &&);
    
	~CISolver() = default;
  
    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<CISolver<MatsU,IntsT>>
    convert(const std::shared_ptr<CISolver<MatsT,IntsT>>&);

    // get Algorithm
    CIDiagonalizationAlgorithm getAlg(){ return alg_; }

    // switch Algrithm;
    void switchAlgorithm(CIDiagonalizationAlgorithm alg, 
      size_t maxIter = 128, double vectorConv = 1.0e-6, 
      size_t maxDSpace = 50, size_t nDGuess = 3,
      std::vector<std::pair<double, size_t>> eRefs = {}) {
        alg_ = alg;
        maxIter_ = maxIter;
        vectorConv_ = vectorConv;
        maxDavidsonSpace_ = maxDSpace;
        nDavidsonGuess_ = nDGuess;
        if (!eRefs.empty()) energyRefs_ = eRefs;
    }

    // solve CI
	virtual void solveCI(MCWaveFunction<MatsT,IntsT> &,EMPerturbation &);
  
  }; // class CISolver

  template <typename MatsT, typename IntsT>
  class MCSCF : public MCWaveFunction<MatsT,IntsT> {
  
  protected:
    // Useful Typedefs
    typedef MatsT *                   oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;
  
  public:
    
    MCSCFSettings settings;
    std::shared_ptr<CISolver<MatsT,IntsT>>        ciSolver  = nullptr;
    std::shared_ptr<OrbitalRotation<MatsT,IntsT>> moRotator = nullptr;  

    // Reduced density Matrices (RDMs) only span over correlated space
    // SOI: state of interest, for orbital rotation
    // it's either state specific or state averaged RDM
    std::shared_ptr<cqmatrix::Matrix<MatsT>>    oneRDMSOI = nullptr;
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI = nullptr;
    // Disable default, copy and move constructors
    MCSCF()              = delete;
    MCSCF(const MCSCF &) = delete;
    MCSCF(MCSCF &&)      = delete;

    // Constructors
    
    /**
     *  \brief MCSCF Constructor.
     *
     *  Stores references to a "reference" SingleSlater object
     *  and makes a copy of the reference into a complex
     *  SingleSlater object for the propagation.
     */ 
    template <typename MatsU>
    MCSCF(SingleSlater<MatsU,IntsT> & ref, size_t NS) : 
        MCWaveFunction<MatsT,IntsT>(ref, NS) { };  // MCSCF constructor
  
    ~MCSCF(){ dealloc(); }

    // MCSCF procedural functions
    virtual void run(EMPerturbation &);       // From MCWaveFunctionBase

    // Functions that generally are wrappers around MCWavefunction, 
    // but might be overwritten by classes which derive from MCSCF
    virtual void transformInts(EMPerturbation & pert);
    virtual void computeMultipole();
    virtual void populationAnalysis();
    virtual void spinAnalysis();
    virtual double oscillator_strength(size_t, size_t);
    virtual void formNaturalOrbitals();

    
    // compute all RDM (and state average) if no inputs 
    virtual void computeOneRDM();
    virtual void computeOneRDM(size_t);
    void computeTwoRDM();
    void computeTwoRDM(size_t);

    void saveCurrentStates(bool prop = false);

    void printStateEnergy();
    void printMCSCFHeader(EMPerturbation &);
    void printMCSCFFooter();

    // Memory functions
    virtual void alloc();
    void dealloc();

  }; // class MCSCF

  // For NEO
  template <typename MatsT, typename IntsT>
  class NEOMCSCF : public MCSCF<MatsT,IntsT>
  {
    private:

    protected:

    public:
    
    // Reference to the NEOSS object so that we can access any NEO funcationality as needed
    NEOSS<MatsT, IntsT> & neoref_;

    // SMG 07/08/24
    // Eventually the block below will be replaced with just the ePTF, and the individual
    // MCWaveFunction objects will be stored in a vector (or unordered map akin to NEOSS)
    // which can then use similar ApplyToEach methodology for the "SameParticle" 
    // terms while NEOMCSCF handles the "CrossParticle" terms

    // Transformed integrals for the correlating space
    // Storing for now just three separate sets of integrals for the
    // (ee|ee), (PP|PP), and (ee|PP) sets of integrals
    std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> eeTF;
    std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> PPTF;
    std::shared_ptr<MixedMOIntsTransformer<MatsT,IntsT>> ePTF;

    // Following the structure of NEOSS, hold onto two MCWavefunction objects
    // For now, explicit storage of the two, later can take the NEOSS approach
    // of holding onto an arbitrary sized vector of subsystems
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> ewfn_;
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pwfn_;

    // TODO: The individual MCWaveFunctions really should hold onto their own
    // 1RDM's.  Still required to be calculated at the NEOMCSCF level since the
    // CI vector is defined in the |a>|b>|p> basis, but should then pass these
    // down to the individual wavefunctions
    std::vector<cqmatrix::Matrix<MatsT>> PoneRDM;

    // An additional pointer to the base ciBuilder but which now holds onto
    // the specific NEOCIBuilder
    std::shared_ptr<NEOCASCI<MatsT,IntsT>> NEOCIBuilder;

    // Explicit storage for the cross two electron integral terms
    std::shared_ptr<TwoPInts<IntsT>> interIntegrals;
    bool contractfirst;

    // Important things inhereted from MCWaveFunction:
    // -> CIVecs: can just be inhereted as our CI Vector is just a column vector of 
    //            MatsT * type.
    // -> ciBuilder: This should just be inheritable so long as I'm able to inheret
    //               the ciBuilder class into a NEOciBuilder child class?
    // -> oneRDM,TDM's: Any sort of analysis probably will need overwritten, but 
    //                  for now just let it be the base class

    NEOMCSCF(NEOSS<MatsT,IntsT> & neoref, size_t NS)
    : neoref_(neoref),
      MCSCF<MatsT,IntsT>(neoref,NS)
    {
      // Generate the MCWavefunction objects
      ewfn_ = std::make_shared<MCSCF<MatsT,IntsT>>(*(neoref_.getSubSS(std::string("Electronic"))),NS);
      pwfn_ = std::make_shared<MCSCF<MatsT,IntsT>>(*(neoref_.getSubSS(std::string("Protonic"))),NS);

      // Make the MOIntTransformers for the two single slater classes
      eeTF = ewfn_->mointsTF;
      PPTF = pwfn_->mointsTF;
      //eeTF = neoref_.getSubSS(std::string("Electronic"))->generateMOIntsTransformer(TPI_TRANSFORMATION_ALG::INCORE_N5);
      //PPTF = neoref_.getSubSS(std::string("Protonic"))->generateMOIntsTransformer(TPI_TRANSFORMATION_ALG::INCORE_N5);

      // All wavefunction components should have access to the integrals
      // (useful for resusing CI Builder code)
      ewfn_->moints = this->moints;
      pwfn_->moints = this->moints;

      // Grab the crossed integrals between the two 
      interIntegrals = neoref_.getCrossTPIs(std::string("Electronic"),std::string("Protonic")).second;
      contractfirst = neoref_.getCrossTPIs(std::string("Electronic"),std::string("Protonic")).first;
      ePTF = std::make_shared<MixedMOIntsTransformer<MatsT,IntsT>>(*(neoref_.getSubSS(std::string("Electronic"))),
                                                                   *(neoref_.getSubSS(std::string("Protonic"))),
                                                                   interIntegrals,
                                                                   contractfirst); 


      // Make the CI Solver of custom type for our class
      // alloc();
    };

      ~NEOMCSCF(){};
      void transformInts(EMPerturbation & pert) override
      {
        transformMultipleInts(pert);
      }
      void transformMultipleInts(EMPerturbation &);

      void printMCSCFHeader(EMPerturbation &);
      void printMCSCFFooter();
      void printNEOMCSCFState(std::ostream & out, size_t i, double energy, MatsT * C, std::vector<size_t> & sorted_CAddr, size_t N);
      void printDetOrder(std::ostream&,std::shared_ptr<DetStringManager>&);

      // Function call which creates the CI solver and allocates requried memory
      void alloc() override;

      // RDM functionality
      void computeOneRDM(size_t) override;
      void computeOneRDM() override;

      void formNaturalOrbitals() override;

      // Functionality for properties of 1 proton systems
      void ProtonExpectationValue(size_t);
      void ProtonVariance(size_t);
      void ProtonKE(size_t);
      std::array<double,3> ProtonExpectationValue(cqmatrix::Matrix<MatsT>);
      std::array<double,3> ProtonVariance(cqmatrix::Matrix<MatsT>);
      double ProtonKE(cqmatrix::Matrix<MatsT>);

      // Moment calculator
      void computeMultipole() override;
      void computeMultipole(size_t) override;

      // Overwrite MCWaveFunction cubegen functionality
      void runCube(std::vector<std::shared_ptr<CubeGen>>) override;

      // Transform RDM's into PDMs
      std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getOnePDM();
      std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getPOnePDM();

      std::shared_ptr<MCWaveFunctionBase> get_ewfn(){return std::dynamic_pointer_cast<MCWaveFunctionBase>(ewfn_);};
      std::shared_ptr<MCWaveFunctionBase> get_pwfn(){return std::dynamic_pointer_cast<MCWaveFunctionBase>(pwfn_);};

      // Get oscillator strengths for the
      double oscillator_strength(size_t, size_t s1 = 0) override;
      
  }; // class NEOMCSCF

}; // namespace ChronusQ
