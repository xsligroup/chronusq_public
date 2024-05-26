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

#include <mcwavefunction.hpp>
#include <quantum/preprocessor.hpp>
#include <util/preprocessor.hpp>
#include <util/print.hpp>


// Template for a collective operation on the members of a 
// MCWaveFunction object
// #define DEBUG_MULTISTATEWFN_IMPL
  
#define MCWaveFunction_COLLECTIVE_OP(OP_OP,OP_VEC_OP) \
  /* Handle densities and vectors*/\
  OP_VEC_OP(MatsT,this,other,CIVecs); 

namespace ChronusQ {

  /**
   *  Constructs a MCWaveFunction object from another of a another (possibly the 
   *  same) type by copy.
   *
   *  \param [in] other MCWaveFunction object to copy
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to copy constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  MCWaveFunction<MatsT,IntsT>::MCWaveFunction(const MCWaveFunction<MatsU,IntsT> &other,int dummy) :
    moints(other.moints),
    ciBuilder(CIBuilder<MatsU,IntsT>::template convert<MatsT>(other.ciBuilder)),
    MCWaveFunctionBase(dynamic_cast<const MCWaveFunctionBase &>(other)),
    ref_(dynamic_cast<SingleSlater<MatsT,IntsT>&>(other.reference())) { 

#ifdef DEBUG_MULTISTATEWFN_IMPL
    std::cout << "MCWaveFunction<T>::MCWaveFunction(const MCWaveFunction<U>&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif
    
    mointsTF = ref_.generateMOIntsTransformer();

    oneRDM.reserve(NStates);
    for (const cqmatrix::Matrix<MatsU> &mat : other.oneRDM) 
       oneRDM.emplace_back(mat);
    
    MCWaveFunction_COLLECTIVE_OP(COPY_OTHER_MEMBER_OP, COPY_OTHER_MEMBER_VEC_OP);

  }; // MCWaveFunction<T>::MCWaveFunction(const MCWaveFunction<U> &)

  /**
   *  Constructs a MCWaveFunction object from another of a another (possibly the 
   *  same) type by move.
   *
   *  \warning Deallocates the passed MCWaveFunction object
   *
   *  \param [in] other MCWaveFunction object to move
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to move constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  MCWaveFunction<MatsT,IntsT>::MCWaveFunction(MCWaveFunction<MatsU,IntsT> &&other, int dummy) : 
    moints(other.moints), 
    ciBuilder(CIBuilder<MatsU,IntsT>::template convert<MatsT>(other.ciBuilder)),
    MCWaveFunctionBase(dynamic_cast<MCWaveFunctionBase &&>(std::move(other))),
    ref_(dynamic_cast<SingleSlater<MatsT,IntsT>&>(other.reference())) {

#ifdef DEBUG_MULTISTATEWFN_IMPL
    std::cout << "MCWaveFunction<T>::MCWaveFunction(MCWaveFunction<U>&&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif
    
    mointsTF = ref_.generateMOIntsTransformer();
    
    oneRDM.reserve(NStates);
    for (cqmatrix::Matrix<MatsU> &mat : other.oneRDM) 
       oneRDM.emplace_back(std::move(mat));

    MCWaveFunction_COLLECTIVE_OP(MOVE_OTHER_MEMBER_OP, MOVE_OTHER_MEMBER_VEC_OP);

  }; // MCWaveFunction<T>::MCWaveFunction(MCWaveFunction<U> &&)

  // Delagate the copy constructor to the conversion constructors
  template <typename MatsT, typename IntsT>
  MCWaveFunction<MatsT,IntsT>::MCWaveFunction(const MCWaveFunction<MatsT, IntsT> &other) : 
    MCWaveFunction(other,0){ };
  template <typename MatsT, typename IntsT>
  MCWaveFunction<MatsT,IntsT>::MCWaveFunction(MCWaveFunction<MatsT, IntsT> &&other) : 
    MCWaveFunction(std::move(other),0){ };

  /**
   *  Allocates the internal memory a MCWaveFunction object
   */ 
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::alloc() {

#ifdef DEBUG_MULTISTATEWFN_IMPL
    std::cout << "MCWaveFunction::alloc (this = " << this << ")" << std::endl;
#endif

    if (reference().nC == 4 and not this->FourCompNoPair) 
      CErr("NYI for Four Component without NoPair Approximation");
    
    size_t NS     = this->NStates;
    size_t NDet   = this->NDet;
    size_t nCorrO = this->MOPartition.nCorrO;

    CIVecs = std::vector<MatsT*>(NS);
    oneRDM.reserve(NS);

    if (MOPartition.scheme == CAS) {
      ciBuilder = std::make_shared<CASCI<MatsT,IntsT>>();
    } else if (MOPartition.scheme == RAS) {
      ciBuilder = std::make_shared<RASCI<MatsT,IntsT>>();
    } else {
      CErr();
    }

    try {
      for (auto i = 0ul; i < NS; i++) {
        CIVecs[i] = CQMemManager::get().malloc<MatsT>(NDet);
        oneRDM.emplace_back(cqmatrix::Matrix<MatsT>(nCorrO)); 
      }
    } catch (...) {
      CErr("Not enough Memory to allocate CIVector for the specified number of determiants");
    }
  
    // computing list for detstring 
    std::cout << std::endl;
    FormattedLine(std::cout, "Compute Excitation List(s) ...");
    if (detStr)     detStr->computeList();  
    if (detStrBeta) detStrBeta->computeList();  
    std::cout << std::endl;

    if (this->readCI) ReadGuessCIVector();
 
  }; // MCWaveFunction<T>::alloc

  /**
   *  Deallocates the internal memory a MCWaveFunction object
   */ 
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::dealloc() {

#ifdef DEBUG_MULTISTATEWFN_IMPL
    std::cout << "MCWaveFunction::dealloc (this = " << this << ")" << std::endl;
#endif

    oneRDM.clear();
    ciBuilder  = nullptr;    
    MCWaveFunction_COLLECTIVE_OP(DEALLOC_OP_5, DEALLOC_VEC_OP_5);
    
  }; // MCWaveFunction<T>::dealloc


  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::saveCurrentStates(bool sProp) {
  
    ROOT_ONLY(comm); 
  
    // Checkpoint if file exists
    if( savFile.exists() ) {

      std::string prefix = "MCWFN/";

      size_t t_hash = std::is_same<MatsT, double>::value ? 1 : 2;
      
      size_t NS = this->NStates;
      savFile.safeWriteData(prefix + "FIELD_TYPE", &t_hash, {1});
      savFile.safeWriteData(prefix + "NSTATES", &NS, {1});
      savFile.safeWriteData(prefix + "INACT_ENERGY", &(this->InactEnergy), {1});
      savFile.safeWriteData(prefix + "STATE_ENERGY", this->StateEnergy.data(), {NS});

      auto & mopart = this->MOPartition;
      savFile.safeWriteData(prefix + "ORB_INDEX", & (mopart.orbIndices[0]),{mopart.nMO});

      // Save CI coefficients
      for (auto i = 0; i < NS; i++)
        savFile.safeWriteData(prefix + "CIVec_"+std::to_string(i+1), CIVecs[i], {this->NDet});

      // Save properties after SCF
      if( sProp ){

        // Save oscillator strength
        if(NosS1) {
          savFile.safeWriteData(prefix + "OSC_STR", osc_str, {NosS1, NS});
        }

        // Save Multipoles
        if( multipoleMoment ){
          savFile.safeWriteData(prefix + "LEN_ELECTRIC_DIPOLE", & elecDipoles[0][0], {NS, 3});
          savFile.safeWriteData(prefix + "LEN_ELECTRIC_QUADRUPOLE", & elecQuadrupoles[0][0][0], {NS, 3, 3});
          savFile.safeWriteData(prefix + "LEN_ELECTRIC_OCTUPOLE", & elecOctupoles[0][0][0][0], {NS, 3, 3, 3});
        }

      }
    }

  }; // MCWaveFunction<T>::saveCurrentStates

}; // namespace ChronusQ

// Other headers
#include <mcwavefunction/base/impl.hpp> // base implementation
#include <detstringmanager/impl.hpp>    // detstringmanager implementaion
#include <mcwavefunction/moints.hpp>    // MO integral transformation
#include <mcwavefunction/print.hpp>     // print implementaion
#include <mcwavefunction/property.hpp>  // property implementation
#include <mcwavefunction/rdm.hpp>       // density matrix
#include <mcwavefunction/ciguess.hpp>   // Read in CI Vectors
#include <mcwavefunction/cube.hpp>     // Cubegen header
