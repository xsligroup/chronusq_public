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

#include <wavefunction.hpp>
#include <util/preprocessor.hpp>

// Template for a collective operation on the members of a 
// WaveFunction object
  
#define WaveFunction_COLLECTIVE_OP(OP_MEMBER,OP_OP) \
  /* Handle densities */\
  OP_OP(double,this,other,eps1); \
  OP_OP(double,this,other,eps2); \



namespace ChronusQ {

  /**
   *  Constructs a WaveFunction object from another of a another (possibly the 
   *  same) type by copy.
   *
   *  \param [in] other WaveFunction object to copy
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to copy constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  WaveFunction<MatsT,IntsT>::WaveFunction(const WaveFunction<MatsU,IntsT> &other,int dummy) :
    aoints_(other.aoints_),
    QuantumBase(dynamic_cast<const QuantumBase &>(other)),
    WaveFunctionBase(dynamic_cast<const WaveFunctionBase &>(other)),
    Quantum<MatsT>(dynamic_cast<const Quantum<MatsU>&>(other)) {
    //molecule_(other.molecule_), basisSet_(other.basisSet_) {

    mo.reserve(2);
    for (const cqmatrix::Matrix<MatsU> &mat : other.mo)
      mo.emplace_back(mat);

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction<T>::WaveFunction(const WaveFunction<U>&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif
    WaveFunction_COLLECTIVE_OP(COPY_OTHER_MEMBER,COPY_OTHER_MEMBER_OP);

  }; // WaveFunction<T>::WaveFunction(const WaveFunction<U> &)

  /**
   *  Constructs a WaveFunction object from another of a another (possibly the 
   *  same) type by move.
   *
   *  \warning Deallocates the passed WaveFunction object
   *
   *  \param [in] other WaveFunction object to move
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to move constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  WaveFunction<MatsT,IntsT>::WaveFunction(WaveFunction<MatsU,IntsT> &&other, int dummy) : 
    aoints_(other.aoints_),
    QuantumBase(dynamic_cast<QuantumBase &&>(std::move(other))),
    WaveFunctionBase(dynamic_cast<WaveFunctionBase &&>(std::move(other))),
    Quantum<MatsT>(dynamic_cast<Quantum<MatsU>&&>(std::move(other))) {
    //molecule_(other.molecule_), basisSet_(other.basisSet_) {

    mo.reserve(2);
    for (cqmatrix::Matrix<MatsU> &mat : other.mo)
      mo.emplace_back(std::move(mat));

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction<T>::WaveFunction(WaveFunction<U>&&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    WaveFunction_COLLECTIVE_OP(MOVE_OTHER_MEMBER,MOVE_OTHER_MEMBER_OP);

  }; // WaveFunction<T>::WaveFunction(WaveFunction<U> &&)

  // Delagate the copy constructor to the conversion constructors
  template <typename MatsT, typename IntsT>
  WaveFunction<MatsT,IntsT>::WaveFunction(const WaveFunction<MatsT,IntsT> &other) : 
    WaveFunction(other,0){ };
  template <typename MatsT, typename IntsT>
  WaveFunction<MatsT,IntsT>::WaveFunction(WaveFunction<MatsT,IntsT> &&other) : 
    WaveFunction(std::move(other),0){ };





  /**
   *  Allocates the internal memory a WaveFunction object
   */ 
  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::alloc() {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction::alloc (this = " << this << ")" << std::endl;
#endif

    size_t NB = this->nC * this->nAlphaOrbital();

    mo.reserve(2);
    mo.emplace_back(NB);
    eps1 = CQMemManager::get().malloc<double>(NB);

    if( this->nC == 1 and (not this->iCS) ) {
      mo.emplace_back(NB);
      eps2 = CQMemManager::get().malloc<double>(NB);
    }

  }; // WaveFunction<T>::alloc



  /**
   *  Deallocates the internal memory a WaveFunction object
   */ 
  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::dealloc() {

#ifdef _WaveFunctionDebug
    std::cout << "WaveFunction::dealloc (this = " << this << ")" << std::endl;
#endif

    mo.clear();
    WaveFunction_COLLECTIVE_OP(DUMMY3,DEALLOC_OP_5);

  }; // WaveFunction<T>::dealloc


}; // namespace ChronusQ



// Other headers
#include <wavefunction/print.hpp> // Print header


