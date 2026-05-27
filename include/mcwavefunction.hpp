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
#include <singleslater.hpp>
#include <singleslater/neoss.hpp>
#include <integrals.hpp>
#include <mcwavefunction/base.hpp>

namespace ChronusQ {

  // Declaration of CI Engine.
  template <typename MatsT, typename IntsT>
  class CIBuilder;

  template <typename MatsT, typename IntsT>
  class CASCI;

  template <typename MatsT, typename IntsT>
  class RASCI;

  template <typename MatsT, typename IntsT>
  class NEOCASCI;

  // Forward declaration of MCSCF class
  template <typename MatsT, typename IntsT>
  class MCSCF;


  /**
   *  \brief The MCWaveFunction class. The typed abstract interface for all
   *  classes for which the wave function is described by a single slater
   *  determinant (HF, KS, PHF, etc).
   *
   *  Adds knowledge of storage type to MCWaveFunctionBase
   *
   */
 
  template <typename MatsT, typename IntsT>
  class MCWaveFunction : public MCWaveFunctionBase, 
    public std::enable_shared_from_this<MCWaveFunction<MatsT,IntsT>> {

  template <typename MatsU, typename IntsU>
  friend class MCSCF;

  protected:  
    
    // Useful Typedefs
    typedef MatsT *                oper_t;
    typedef std::vector<oper_t>    oper_t_coll;
    
    
    
  public:

    // Hacky for now to make this public while calculating 1PDMS from 1RDMS in NEO
    SingleSlater<MatsT, IntsT> &   ref_; 
    
	// Integrals here are computed and stored in correalted space
	// Only one set of integrals, means not working for UHF reference
    std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> mointsTF; 
    std::shared_ptr<IntegralsCollection> moints = std::make_shared<IntegralsCollection>();
	// fold the following intgrals into moints
	//oper_t moERI;   // Transformed MO 2e integral in correalted space
    //oper_t hCore;   // 1e integral with frozen core contribution
    //oper_t hCoreP;  // hCore with 2e-contribution folded in.
    
    // CI Vectors
    oper_t_coll CIVecs; 
    std::shared_ptr<CIBuilder<MatsT,IntsT>>  ciBuilder  = nullptr;
    
    // Reduced density Matrices (RDMs) only span over correlated space
    std::vector<cqmatrix::Matrix<MatsT>> oneRDM;
    // oneRDMSOI belongs to MCSCF, but we will maintain a copy here for potential use
    std::shared_ptr<cqmatrix::Matrix<MatsT>> oneRDMSOI;
    std::shared_ptr<InCore4indexTPI<MatsT>> twoRDMSOI;
    // Transition density matrix
    std::vector<std::vector<cqmatrix::Matrix<MatsT>>> TDMs;

    // TODO: Save vectors 
    //oper_t_coll DOSav;
    //oper_t_coll UH;

    /**
     *  \brief MCWaveFunction Constructor by MCWaveFunction
     *
     *  Stores references to a "reference" MCWaveFunction object and
     *  makes a copy of the reference into a complex
     *  MCWaveFunction object for the propagation.
     */ 
    template <typename MatsU>
    MCWaveFunction(SingleSlater<MatsU,IntsT> &ref, size_t NS):
      MCWaveFunctionBase(ref.comm, NS),
      ref_(ref) {
      
      //if (std::is_same<IntsT, dcomplex>::value) 
      //   CErr("MCWaveFunction with dcomplex IntsT is not tested yet!");
      
      mointsTF = ref_.generateMOIntsTransformer();
	
    };  // MCWaveFunction constructor
  
    // See include/mcwavefunction/impl.hpp for documentation 
    MCWaveFunction() = delete;
	// on the following constructors

    // Different type
    template <typename MatsU> 
      MCWaveFunction(const MCWaveFunction<MatsU,IntsT> &, int dummy = 0);
    template <typename MatsU> 
      MCWaveFunction(MCWaveFunction<MatsU,IntsT> &&     , int dummy = 0);

    // Same type
    MCWaveFunction(const MCWaveFunction<MatsT,IntsT> &);
    MCWaveFunction(MCWaveFunction<MatsT,IntsT> &&);     
    
    ~MCWaveFunction(){ dealloc(); }

    // MCWaveFunction procedural functions
    //virtual void run(EMPerturbation &) = 0;  // From MCWaveFunctionBase

    void computeOverlaps(oper_t, std::vector<MatsT>&); // Calculate the overlaps of an arbitrary CI vector with the CI vectors.
    virtual void computeOneRDM(size_t);    
    virtual void computeOneRDM();
    virtual void computeTDMs(); // compute TDMs
    void rdm2pdm(cqmatrix::Matrix<MatsT> &, double scale = 1., bool isTDM = false);
    void pdm2rdm(cqmatrix::Matrix<MatsT> &);

    WaveFunctionBase & referenceWaveFunction() { return dynamic_cast<WaveFunctionBase&>(ref_); }
    size_t getnC() const { return dynamic_cast<WaveFunctionBase&>(ref_).nC; }
    BASIS_FUNCTION_TYPE getBasisType() const { return dynamic_cast<WaveFunctionBase&>(ref_).basisSet_.basisType; }
    SingleSlater<MatsT,IntsT> & reference() const  { return ref_;} 
    void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>& moPairs, SpinType sp) {
      this->reference().swapMOs(moPairs,sp);
    };
    
    void ReadGuessCIVector();
    virtual void saveCurrentStates(bool);
    void setMORanges();
    virtual void transformInts(EMPerturbation &);
    virtual void printMOSpacePartition();
    void printMCStates();
    virtual void printMCState(std::ostream&,size_t,double,MatsT*,std::vector<size_t>&,size_t,const size_t n_item_per_row = 5);
    void print1RDMs();
    void printMOInfo(std::ostream&, size_t a = 0);

    // Properties
    void populationAnalysis(size_t);
    void populationAnalysis();
    virtual double oscillator_strength(size_t, size_t s1 = 0);
    virtual void computeMultipole(size_t);
    virtual void computeMultipole();
    virtual void formNaturalOrbs(size_t);
    virtual void spinAnalysis(size_t);
    virtual void spinAnalysis();

    // For dealing with electric fields
    void precompute_NucEField(EMPerturbation &);

    virtual void run(EMPerturbation & )
      {CErr("Run Called from MCWaveFunction Invalid!");};

    // Post-processing functions
    void runCube(std::vector<std::shared_ptr<CubeGen>> cu);

    // Memory functions
    virtual void alloc();
    void dealloc();

    // Get shared_ptr
    //template<typename MatsU, typename IntsU>
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> getPtr(){return this->shared_from_this();};

    void addMCWaveFunction(std::shared_ptr<MCWaveFunction>,std::string){CErr("");};

  }; // class MCWaveFunction

  /**
   *  \brief The MultiComponentMCWaveFunction class. 
   * 
   *  \todo The goal is to eventually run MCSCF & other Multiconfigurational calculations
   *        exclusively through a class which looks like MultiComponentMCWaveFunction, 
   *        i.e. a container which holds a vector of underlying particle wavefunctions,
   *        and function calls typically follow the ApplyToEach construct of NEOSS
   *
   */
 
  template <typename MatsT, typename IntsT>
  class MultiComponentMCWaveFunction : public MCWaveFunction<MatsT,IntsT>, 
    public std::enable_shared_from_this<MultiComponentMCWaveFunction<MatsT,IntsT>> {

    private:

    protected:

    public:

      // Data members

      // Akin to NEOSS main storage
      std::unordered_map<std::string,std::shared_ptr<MCWaveFunction<MatsT,IntsT>>> subsystems;
      std::vector<std::string> order_;

      // Storage for the cross-particle integrals
      std::unordered_map<std::string,std::unordered_map<std::string,std::shared_ptr<TwoPInts<IntsT>>>> interIntegrals;


      // Constructors
      MultiComponentMCWaveFunction(SingleSlater<MatsT,IntsT> & ref, size_t NS)
      : MCWaveFunction<MatsT,IntsT>(ref,NS)
      {
        
      };
      

      // Add an additional wavefunction to the existing list
      //virtual void addMCWaveFunction(std::shared_ptr<MCWaveFunctionBase>,std::string label);

  
  };

  /**
   *  \brief The NEOMCWaveFunction class for NEO-CI calculations
   * 
   *  \todo For now this is a specialization of the MultiComponentMCWaveFunction class, but eventually
   *        that class will be generalized and this class can likely be removed
   * 
   *
   */
 
  template <typename MatsT, typename IntsT>
  class NEOMCWaveFunction : public MultiComponentMCWaveFunction<MatsT,IntsT>, 
    public std::enable_shared_from_this<NEOMCWaveFunction<MatsT,IntsT>> {

    private:

    protected:
    
    public:

      // NEO-MCWaveFunction Data Members

      // Reference to the NEOSS object so that we can access any NEO funcationality as needed
      NEOSS<MatsT,IntsT> & neoref_;

      // Transformed integrals for the correlating space
      // Storing for now just three separate sets of integrals for the
      // (ee|ee), (PP|PP), and (ee|PP) sets of integrals
      std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> eeTF;
      std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> PPTF;
      std::shared_ptr<MixedMOIntsTransformer<MatsT,IntsT>> ePTF;

      // For convenience (for now), store individual pointers to the electronic
      // and nuclear wavefunctions in the underlying MultiComponentMCWaveFunction
      // storage vector
      std::shared_ptr<MCWaveFunction<MatsT,IntsT>> ewfn_;
      std::shared_ptr<MCWaveFunction<MatsT,IntsT>> pwfn_;

      // An additional pointer to the base ciBuilder but which now holds onto
      // the specific NEOCIBuilder
      // This normally belongs to MCSCF but we'll stash a copy here for convenience (For now)
      std::shared_ptr<NEOCASCI<MatsT,IntsT>> NEOCIBuilder;

      // Explicit storage for the cross two electron integral terms
      std::shared_ptr<TwoPInts<IntsT>> interIntegrals;
      bool contractfirst;

      // Constructors
      NEOMCWaveFunction(NEOSS<MatsT,IntsT> & neoref, size_t NS):
      neoref_(neoref),
      MultiComponentMCWaveFunction<MatsT,IntsT>(neoref,NS)
      {

        // Grab the crossed integrals between the two 
        interIntegrals = neoref_.getCrossTPIs(std::string("Electronic"),std::string("Protonic")).second;
        contractfirst = neoref_.getCrossTPIs(std::string("Electronic"),std::string("Protonic")).first;
        ePTF = std::make_shared<MixedMOIntsTransformer<MatsT,IntsT>>(*(neoref_.getSubSS(std::string("Electronic"))),
                                                                    *(neoref_.getSubSS(std::string("Protonic"))),
                                                                    interIntegrals,
                                                                    contractfirst); 

      };
      // Might actually be some cleaning up to do I'm missing....
      ~NEOMCWaveFunction(){};
      
      void addMCWaveFunction(std::shared_ptr<MCWaveFunctionBase>,std::string label) override;

      // Integral transform
      void transformInts(EMPerturbation & pert) override
        {transformMultipleInts(pert);}; 
      void transformMultipleInts(EMPerturbation &);

      // Printing Functionality
      void printMOSpacePartition() override;
      void printMCState(std::ostream&,size_t,double,MatsT*,std::vector<size_t>&,size_t,const size_t n_item_per_row = 5)override;
      void printNEOMCState(std::ostream & out, size_t i, double energy, MatsT * C, std::vector<size_t> & sorted_CAddr, size_t N);
      void printDetOrder(std::ostream&,std::shared_ptr<DetStringManager>&);
      void printMOInfo(std::ostream &,size_t);
      
      void formNaturalOrbs(size_t) override;
      
      // Function call which creates the CI solver and allocates requried memory
      void alloc() override;

      // RDM functionality
      void computeOneRDM(size_t) override;
      void computeOneRDM() override;

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
      std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getOnePDM();
      std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> getPOnePDM();
      void runCube(std::vector<std::shared_ptr<CubeGen>>) override;

      // Get oscillator strengths for the
      double oscillator_strength(size_t, size_t s1 = 0) override;
 
  };



}; // namespace ChronusQ

// include declaration of CIBuilder
#include <cibuilder.hpp>


