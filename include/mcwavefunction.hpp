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


  /**
   *  \brief The MCWaveFunction class. The typed abstract interface for all
   *  classes for which the wave function is described by a single slater
   *  determinant (HF, KS, PHF, etc).
   *
   *  Adds knowledge of storage type to MCWaveFunctionBase
   *
   */
 
  template <typename MatsT, typename IntsT>
  class MCWaveFunction : public MCWaveFunctionBase {

  protected:  
    
    // Useful Typedefs
    typedef MatsT *                oper_t;
    typedef std::vector<oper_t>    oper_t_coll;
    
    SingleSlater<MatsT, IntsT> &   ref_; 
    
    bool cacheHalfTransTPI_ = false;
    
  public:
    
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
    virtual void run(EMPerturbation &) = 0;  // From MCWaveFunctionBase

    void computeOverlaps(oper_t, std::vector<MatsT>&); // Calculate the overlaps of an arbitrary CI vector with the CI vectors.
    virtual void computeOneRDM(size_t);    
    virtual void computeOneRDM();
    virtual void computeTDMs(); // compute TDMs
    void rdm2pdm(cqmatrix::Matrix<MatsT> &, double scale = 1.);
    
    WaveFunctionBase & referenceWaveFunction() { return dynamic_cast<WaveFunctionBase&>(ref_); }
    SingleSlater<MatsT,IntsT> & reference() const  { return ref_;} 
    void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>& moPairs, SpinType sp) {
      this->reference().swapMOs(moPairs,sp);
    };
    
    void ReadGuessCIVector();
    virtual void saveCurrentStates(bool);
    void setMORanges();
    void transformInts(EMPerturbation &);
    void printMOSpacePatition();
    void print1RDMs();
    void printMOInfo(std::ostream&, size_t a = 0);

    // Properties
    void populationAnalysis(size_t);
    void populationAnalysis();
    void spinAnalysis(size_t);
    void spinAnalysis();
    double oscillator_strength(size_t, size_t s1 = 0);
    void computeMultipole(size_t);
    void computeMultipole();
    // For dealing with electric fields
    void precompute_NucEField(EMPerturbation &);

    // Post-processing functions
    void runCube(std::vector<std::shared_ptr<CubeGen>> cu, EMPerturbation &emPert);

    // Memory functions
    void alloc();
    void dealloc();

  }; // class MCWaveFunction

}; // namespace ChronusQ

// include declaration of CIBuilder
#include <cibuilder.hpp>


