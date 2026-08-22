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
#include <singleslater/multiparticless.hpp>
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
  class MultiParticleCASCI;

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
    std::shared_ptr<SingleSlater<MatsT, IntsT>>   ref_; 
    
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
    MCWaveFunction(std::shared_ptr<SingleSlater<MatsU,IntsT>> ref, size_t NS):
      MCWaveFunctionBase(ref->comm, NS),
      ref_(ref)
      {
      
      //if (std::is_same<IntsT, dcomplex>::value) 
      //   CErr("MCWaveFunction with dcomplex IntsT is not tested yet!");
      
      mointsTF = ref_->generateMOIntsTransformer();
	
    };  // MCWaveFunction constructor

    //MCWaveFunction(std::shared_ptr<SingleSlater<MatsT,IntsT>> ref, size_t NS):
    //  MCWaveFunctionBase(ref->comm, NS)
    //  {
    //  ref_ = ref;
    //  
    //  //if (std::is_same<IntsT, dcomplex>::value) 
    //  //   CErr("MCWaveFunction with dcomplex IntsT is not tested yet!");
    //  
    //  mointsTF = ref_->generateMOIntsTransformer();
	
    //};  // MCWaveFunction constructor
   
   
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
    virtual void computeTwoRDM(size_t, InCore4indexTPI<MatsT> &);
    std::shared_ptr<InCore4indexTPI<MatsT>> computeFull2RDM(size_t);
    virtual void computeTDMs(); // compute TDMs
    void rdm2pdm(cqmatrix::Matrix<MatsT> &, double scale = 1., bool isTDM = false);
    void rdm2pdm(size_t index, double scale = 1., bool isTDM = false);
    void pdm2rdm(cqmatrix::Matrix<MatsT> &);

    WaveFunctionBase & referenceWaveFunction() { return dynamic_cast<WaveFunctionBase&>(*ref_); }
    size_t getnC() const { return dynamic_cast<WaveFunctionBase&>(*ref_).nC; }
    BASIS_FUNCTION_TYPE getBasisType() const { return dynamic_cast<WaveFunctionBase&>(*ref_).basisSet_.basisType; }
    SingleSlater<MatsT,IntsT> & reference() const  { return *ref_;} 
    std::shared_ptr<SingleSlater<MatsT,IntsT>>  ptr_reference() const  {return ref_;} 
    void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>& moPairs, SpinType sp) {
      this->reference().swapMOs(moPairs,sp);
    };
    
    void ReadGuessCIVector();
    virtual void saveCurrentStates(bool);
    void setMORanges();
    virtual void transformInts(EMPerturbation &, std::string label = "");
    virtual void printMOSpacePartition(std::string label = "");
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
    virtual void spinAndAngularAnalysis(size_t);
    virtual void spinAndAngularAnalysis();

    // For dealing with electric fields
    void precompute_NucEField(EMPerturbation &);

    virtual void run(EMPerturbation & )
      {CErr("Run Called from MCWaveFunction Invalid!");};

    // Post-processing functions
    void runCube(std::vector<std::shared_ptr<CubeGen>> cu);

    // Memory functions
    virtual void alloc(bool owner = true);
    void dealloc();

    // Get shared_ptr
    //template<typename MatsU, typename IntsU>
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> getPtr(){return this->shared_from_this();};

    void addMCWaveFunction(std::shared_ptr<MCWaveFunction>,std::string){CErr("");};

    double getParticleCharge()const{return this->reference().particle.charge;};

  }; // class MCWaveFunction

  class MultiParticleMCWaveFunctionBase
  {
    public:
      virtual std::shared_ptr<MCWaveFunctionBase> getSubMCWaveFunctionBase(std::string label) = 0;
      virtual std::vector<std::string> getLabels() = 0;
      virtual void addMCWaveFunction(std::string,std::shared_ptr<MCWaveFunctionBase>) = 0;
      virtual void addInteraction(const std::string label1, const std::string label2, const std::shared_ptr<IntegralsBase>&) = 0;
  };

  /**
   *  \brief The MultiParticleMCWaveFunction class. 
   * 
   *  \todo The goal is to eventually run MCSCF & other Multiconfigurational calculations
   *        exclusively through a class which looks like MultiParticleMCWaveFunction,
   *        i.e. a container which holds a vector of underlying particle wavefunctions,
   *        and function calls typically follow the ApplyToEach construct
   *
   */

  struct OneParticleCIBuilderHelper{
    size_t NDet;
    std::string label;
    std::string oneBodyIntsStr;
    std::string twoBodyIntsStr;
    std::shared_ptr<const ExcitationList> exList;
    std::shared_ptr<MCWaveFunctionBase> refwfn;
    OneParticleCIBuilderHelper(std::string label_,
                               size_t NDet_,
                               std::string oneBodyIntsStr_,
                               std::string twoBodyIntsStr_,
                               std::shared_ptr<const ExcitationList> exList_,
                               std::shared_ptr<MCWaveFunctionBase>refwfn_)
                               :label(label_),
                                NDet(NDet_),
                                oneBodyIntsStr(oneBodyIntsStr_),
                                twoBodyIntsStr(twoBodyIntsStr_),
                                exList(exList_),
                                refwfn(refwfn_) {};
    void print(std::ostream &) const;
  };

  struct TwoParticleCIBuilderHelper{
    size_t NDetp1;
    size_t NDetp2;
    size_t NDet;
    std::string labelp1;
    std::string labelp2;
    std::shared_ptr<const ExcitationList> exList_p1;
    std::shared_ptr<const ExcitationList> exList_p2;
    std::shared_ptr<MCWaveFunctionBase> mcwfnp1;
    std::shared_ptr<MCWaveFunctionBase> mcwfnp2;
    size_t nCorrp1;
    size_t nCorrp2;
    std::string twoBodyIntsString;
    double chargeproduct;
    TwoParticleCIBuilderHelper(size_t NDetp1_,
                               size_t NDetp2_,
                               size_t nCorrp1_,
                               size_t nCorrp2_,
                               std::string labelp1_,
                               std::string labelp2_,
                               std::string twoBodyIntsString_,
                               std::shared_ptr<const ExcitationList> exList_p1_,
                               std::shared_ptr<const ExcitationList> exList_p2_,
                               std::shared_ptr<MCWaveFunctionBase> mcwfnp1_,
                               std::shared_ptr<MCWaveFunctionBase> mcwfnp2_):
                               NDetp1(NDetp1_),
                               NDetp2(NDetp2_),
                               nCorrp1(nCorrp1_),
                               nCorrp2(nCorrp2_),
                               labelp1(labelp1_),
                               labelp2(labelp2_),
                               twoBodyIntsString(twoBodyIntsString_),
                               exList_p1(exList_p1_),
                               exList_p2(exList_p2_),
                               mcwfnp1(mcwfnp1_),
                               mcwfnp2(mcwfnp2_) 
                               {NDet = NDetp1*NDetp2;
                                chargeproduct = mcwfnp1->getParticleCharge() * mcwfnp2->getParticleCharge();};
    void print(std::ostream&) const;
  };

  template <typename MatsT, typename IntsT>
  class MultiParticleMCWaveFunction : public MCWaveFunction<MatsT,IntsT>, public MultiParticleMCWaveFunctionBase,
    public std::enable_shared_from_this<MultiParticleMCWaveFunction<MatsT,IntsT>> {

    private:
      using SubMCWfnPtr = std::shared_ptr<MCWaveFunction<MatsT,IntsT>>;
      // Data members
      std::shared_ptr<MultiParticleSS<MatsT,IntsT>> mcSSref_;
      std::shared_ptr<MultiParticleCASCI<MatsT,IntsT>> multiparticleciBuilder;

      // Akin to MultiParticleSS main storage
      std::unordered_map<std::string,std::shared_ptr<MCWaveFunction<MatsT,IntsT>>> subsystems;
      std::vector<std::string> order_;
      std::vector<std::string> ciBuilderorder_;
      std::unordered_map<std::string,size_t> rdmorder_;

      // Storage for the cross-particle integrals
      std::unordered_map<std::string,std::unordered_map<std::string,std::shared_ptr<TwoPInts<IntsT>>>> interIntegrals;
      std::unordered_map<std::string,std::unordered_map<std::string,std::shared_ptr<MixedMOIntsTransformer<MatsT,IntsT>>>>interparticleMOTransformers;

      // Helpers for MultiParticleCIBuilder
      // (Essentially, unpack all alpha / beta into individual components)
      std::unordered_map<std::string,OneParticleCIBuilderHelper> oneparticlebuilders;
      std::unordered_map<std::string,std::vector<TwoParticleCIBuilderHelper>> twoparticlebuilders;


      template <typename F>
      void ApplyToEach(F func){
        for(auto & label: order_)
          func(subsystems.at(label));
      }
      template <typename F>
      void ApplyToEachLabeled(F func){
        for(auto & label: order_)
          func(subsystems.at(label),label);
      }
    
    protected:

    public:

      // Constructors
      MultiParticleMCWaveFunction(std::shared_ptr<MultiParticleSS<MatsT,IntsT>> ref, size_t NS)
      : mcSSref_(ref),
        MCWaveFunction<MatsT,IntsT>(std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ref),NS),
        MultiParticleMCWaveFunctionBase()
      {
        
      };

      std::vector<std::string> getOrder()const{return order_;};
      std::vector<std::string> getCIOrder()const{return ciBuilderorder_;};
      std::unordered_map<std::string,size_t> getRDMOrder()const{return rdmorder_;};
      OneParticleCIBuilderHelper getOneParticleCIHelper(std::string sys){return oneparticlebuilders.at(sys);};
      std::vector<TwoParticleCIBuilderHelper> getTwoParticleCIHelper(std::string sys){return twoparticlebuilders.at(sys);};
      // Add an additional wavefunction to the existing list
      std::shared_ptr<MCWaveFunctionBase> getSubMCWaveFunctionBase(std::string label) override;
      std::shared_ptr<MCWaveFunction<MatsT,IntsT>> getSubMCWaveFunction(std::string label);
      std::vector<std::string> getLabels() override;
      void addMCWaveFunction(std::string,std::shared_ptr<MCWaveFunctionBase>) override;
      void addInteraction(const std::string label1, const std::string label2, const std::shared_ptr<IntegralsBase>&) override;


      // Functions which override MCWaveFunction functionality
      void alloc(bool owner) override;
      void printMOSpacePartition(std::string label) override;
      void transformInts(EMPerturbation&,std::string) override;

      // RDM Functionality
      void computeOneRDM() override;
      void computeOneRDM(size_t) override;
      void computeTDMs() override;
      void formNaturalOrbs(size_t);

      // Oscillator strength
      double oscillator_strength(size_t, size_t s1 = 0) override {CErr("Oscillator strengths NYI for MultiParticleCI");};

      // Mulitpole calculation
      void computeMultipole() override;
      void computeMultipole(size_t) override;

      void runCube(std::vector<std::shared_ptr<CubeGen>>) override {CErr("CubeGen NYI for MultiParticleCI!");};
  };



}; // namespace ChronusQ

// include declaration of CIBuilder
#include <cibuilder.hpp>


