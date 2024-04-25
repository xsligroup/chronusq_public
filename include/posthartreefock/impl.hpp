/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you ca redistribute it and/or modify
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

#include <posthartreefock.hpp>
#include <quantum/preprocessor.hpp>
#include <util/preprocessor.hpp>
#include <util/print.hpp>


// Template for a collective operation on the members of a 
// PostHartreeFock object
// #define DEBUG_MULTISTATEWFN_IMPL
  
namespace ChronusQ {

/**
 *  Constructs a PostHartreeFock object from another of a another (possibly the 
 *  same) type by copy.
 *
 *  \param [in] other PostHartreeFock object to copy
 *  \param [in] dummy Dummy argument to fix calling signature for delegation 
 *    to copy constructor
 */ 
template <typename MatsT, typename IntsT>
template <typename MatsU> 
PostHartreeFock<MatsT,IntsT>::PostHartreeFock(const PostHartreeFock<MatsU,IntsT> & other,int dummy) :
  moints(other.moints),
  PostHartreeFockBase(dynamic_cast<const PostHartreeFockBase &>(other)),
  ref_(SingleSlater<MatsU, IntsT>::template convert<MatsT>(other.reference())) {

#ifdef DEBUG_MULTISTATEWFN_IMPL
  std::cout << "PostHartreeFock<T>::PostHartreeFock(const PostHartreeFock<U>&) "
            << "(this = " << this << ", other = " << &other << ")" 
            << std::endl;
#endif
  
  mointsTF = ref_->generateMOIntsTransformer();
  
  alloc();

  for (auto i = 0ul; i < this->NStates; ++i) {
    *oneRDM[i] = *other.oneRDM[i];
  }
  
} // PostHartreeFock<T>::PostHartreeFock(const PostHartreeFock<U> &)

/**
 *  Constructs a PostHartreeFock object from another of a another (possibly the 
 *  same) type by move.
 *
 *  \warning Deallocates the passed PostHartreeFock object
 *
 *  \param [in] other PostHartreeFock object to move
 *  \param [in] dummy Dummy argument to fix calling signature for delegation 
 *    to move constructor
 */ 
template <typename MatsT, typename IntsT>
template <typename MatsU> 
PostHartreeFock<MatsT,IntsT>::PostHartreeFock(PostHartreeFock<MatsU,IntsT> &&other, int dummy) : 
  moints(other.moints), 
  PostHartreeFockBase(dynamic_cast<PostHartreeFockBase &&>(std::move(other))),
  ref_(SingleSlater<MatsU, IntsT>::template convert<MatsT>(other.reference())) {

#ifdef DEBUG_MULTISTATEWFN_IMPL
  std::cout << "PostHartreeFock<T>::PostHartreeFock(PostHartreeFock<U>&&) "
            << "(this = " << this << ", other = " << &other << ")" 
            << std::endl;
#endif
  
  mointsTF = ref_->generateMOIntsTransformer();
  
  alloc();
  for (auto i = 0ul; i < this->NStates; ++i) {
    *oneRDM[i] = *other.oneRDM[i];
    other.oneRDM[i] = nullptr;
  }

} // PostHartreeFock<T>::PostHartreeFock(PostHartreeFock<U> &&)

// Delagate the copy constructor to the conversion constructors
template <typename MatsT, typename IntsT>
PostHartreeFock<MatsT,IntsT>::PostHartreeFock(const PostHartreeFock<MatsT, IntsT> &other) : 
  PostHartreeFock(other,0){ };
template <typename MatsT, typename IntsT>
PostHartreeFock<MatsT,IntsT>::PostHartreeFock(PostHartreeFock<MatsT, IntsT> &&other) : 
  PostHartreeFock(std::move(other),0){ };

/**
 *  Allocates the internal memory a PostHartreeFock object
 */ 
template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::alloc() {

#ifdef DEBUG_MULTISTATEWFN_IMPL
  std::cout << "PostHartreeFock::alloc (this = " << this << ")" << std::endl;
#endif

  if (reference()->nC == 4 and not this->FourCompNoPair) 
    CErr("NYI for Four Component without NoPair Approximation");
  
  size_t NS     = this->NStates;
  size_t nCorrO = this->corrSpace.nCorrO;

  oneRDM.reserve(NS);
  for (auto i = 0ul; i < NS; i++) {
    oneRDM.emplace_back(std::make_shared<cqmatrix::Matrix<MatsT>>(nCorrO)); 
  }

} // PostHartreeFock<T>::alloc

/**
 *  Deallocates the internal memory a PostHartreeFock object
 */ 
template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::dealloc() {

#ifdef DEBUG_MULTISTATEWFN_IMPL
  std::cout << "PostHartreeFock::dealloc (this = " << this << ")" << std::endl;
#endif

  oneRDM.clear();

} // PostHartreeFock<T>::dealloc

template <typename MatsT, typename IntsT>
void PostHartreeFock<MatsT,IntsT>::saveCurrentStates() {

  ROOT_ONLY(comm); 

  // Checkpoint if file exists
  if( savFile.exists() ) {
    
    size_t NS = this->NStates;
    savFile.safeWriteData("POSTHF/NSTATES", &NS, {1});
    savFile.safeWriteData("POSTHF/Core_ENERGY", &(this->coreEnergy), {1});
    savFile.safeWriteData("POSTHF/STATE_ENERGY", this->StateEnergy.data(), {NS}); 

    auto & corrS = this->corrSpace;
    savFile.safeWriteData("POSTHF/ORB_INDEX", & (corrS.orbIndices[0]),{corrS.nMO});

    // Save oscillator strength
    if(NosS1) {
      savFile.safeWriteData("POSTHF/OSC_STR", osc_str, {NosS1, NS});      
    }
  }  

} // PostHartreeFock<T>::saveCuurentStates

} // namespace ChronusQ

// Other headers
#include <posthartreefock/base/impl.hpp> // base implementation
#include <posthartreefock/moints.hpp>    // MO integral transformation
#include <posthartreefock/print.hpp>     // print implementaion
#include <posthartreefock/property.hpp>  // property implementation
#include <posthartreefock/rdm.hpp>       // density matrix

