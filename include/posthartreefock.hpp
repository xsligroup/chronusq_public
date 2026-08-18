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
#include <detfactory.hpp>
#include <integrals.hpp>
#include <posthartreefock/base.hpp>

namespace ChronusQ {

/**
 *  \brief The PostHartreeFock class. The typed abstract interface for all
 *  classes with correlated wave functions (CI, DMRG, PT, CC, etc).
 */

template <typename MatsT, typename IntsT>
class PostHartreeFock : public PostHartreeFockBase {

protected:  
  
public:
  
  std::shared_ptr<SingleSlater<MatsT, IntsT>> ref_;
  
  // 4CAO Dipole Integral:
  std::shared_ptr<std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>>> AODipole4C_ = nullptr;

  // Integrals here are computed and stored in correalted space by default
  // Only one set of integrals, means not working for UHF reference
  std::shared_ptr<MOIntsTransformer<MatsT, IntsT>> mointsTF; 
  std::shared_ptr<IntegralsCollection> moints = std::make_shared<IntegralsCollection>(); ///< MOIntegrals for the storage of integrals
  // fold the following intgrals into moints
  //oper_t moERI;   // Transformed MO 2e integral in correalted space
  //oper_t hCore;   // 1e integral with frozen core contribution
  //oper_t hCoreP;  // hCore with 2e-contribution folded in.
  
  // Reduced density Matrices (RDMs) only span over correlated space
  std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> oneRDM;
  // TODO: Save vectors 
  //oper_t_coll DOSav;
  //oper_t_coll UH;

  /**
   *  \brief PostHartreeFock Constructor by PostHartreeFock
   *
   *  Stores references to a "reference" PostHartreeFock object and
   *  makes a copy of the reference into a complex
   *  PostHartreeFock object for the propagation.
   */ 
  template <typename MatsU>
  PostHartreeFock(std::shared_ptr<SingleSlater<MatsU,IntsT>> ref, size_t NS):
    PostHartreeFockBase(ref->comm, ref, NS),
    ref_(ref) {
    
    if (std::is_same<IntsT, dcomplex>::value) 
       CErr("PostHartreeFock with dcomplex IntsT is not tested yet!");
    
    mointsTF = ref_->generateMOIntsTransformer();
  
  }  // PostHartreeFock constructor

  PostHartreeFock() = delete;
  // Different type
  template <typename MatsU> 
    PostHartreeFock(const PostHartreeFock<MatsU,IntsT> &, int dummy = 0);
  template <typename MatsU> 
    PostHartreeFock(PostHartreeFock<MatsU,IntsT> &&     , int dummy = 0);

  // Same type
  PostHartreeFock(const PostHartreeFock<MatsT,IntsT> &);
  PostHartreeFock(PostHartreeFock<MatsT,IntsT> &&);     
  
  ~PostHartreeFock(){ dealloc(); }

  // PostHartreeFock procedural functions
  virtual void run(EMPerturbation &)      = 0;  // From PostHartreeFockBase
  virtual void computeTDM(size_t, size_t, std::shared_ptr<cqmatrix::Matrix<MatsT>>) = 0;
  virtual void compute2TDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>) = 0;  

  void computeOneRDM(size_t i) { computeTDM(i, i, oneRDM[i]);};    
  void computeOneRDM();
  void rdm2pdm(cqmatrix::Matrix<MatsT> &, double scale = 1., bool isTDM = false);
  void pdm2rdm(cqmatrix::Matrix<MatsT> &);  

  std::shared_ptr<SingleSlater<MatsT,IntsT>> reference() const { return ref_;}
  
  void transformInts(EMPerturbation &, bool);
  void prepareMOIntegrals(EMPerturbation &, const DeterminantFactory&, bool, bool);

  void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>& moPairs, SpinType sp) {
    this->reference()->swapMOs(moPairs,sp);
  }
  
  virtual void saveCurrentStates();
  void setMORanges();
  void printCorrMOSpace();
  void print1RDMs();
  void printMOInfo(std::ostream&, size_t a = 0);

  // Properties
  void populationAnalysis(size_t);
  void populationAnalysis();
  void compute4CAODipole();
  double oscillator_strength(size_t, size_t s1 = 0);
  double secondorder_oscillator_strength(size_t, size_t s1 = 0);
  std::vector<MatsT> computeStateSpecificDipoleMom(size_t);
  void printAllStateSpecificDipoleMom();
  std::vector<MatsT> computeTransitionDipoleMom(size_t, size_t);
  void printAllTransitionDipoleMom();
  double oscillator_strength4C(size_t, size_t s1 = 0);
  void OneRDMDiff(); 
  std::vector<cqmatrix::Matrix<MatsT>> spin_overlap;
  std::vector<cqmatrix::Matrix<MatsT>> spinOverlap();
  void spinAnalysis(size_t, std::vector<cqmatrix::Matrix<MatsT>>*);
  void spinAnalysis();
  void saveOnePDMs(size_t);
  void saveOnePDMs();


  // Post-processing functions
  void runCube(std::vector<std::shared_ptr<CubeGen>>);

  // Memory functions
  void alloc();
  void dealloc();

}; // class PostHartreeFock

} // namespace ChronusQ

// include declaration of CIBuilder
#include <cibuilder.hpp>
