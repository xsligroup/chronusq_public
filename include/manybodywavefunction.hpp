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
#include <manybodywavefunction/base.hpp>

namespace ChronusQ {

  /**
   *  \brief The ManyBodyWavefunction class. The typed abstract interface for all
   *  many body wave functions. This serves as a base class for all many body 
   *  wave function classes (CI, CC, DMRG, etc.).
   *
   *
   */
 
  template <typename MatsT, typename IntsT>
  class ManyBodyWavefunction : public ManyBodyWavefunctionBase {

  protected:  
    
    
  public:
    
    template <typename MatsU>
    ManyBodyWavefunction()
      {
      if (std::is_same<IntsT, dcomplex>::value) 
         CErr("ManyBodyWavefunction with dcomplex IntsT is not tested yet!");
	
    };  // ManyBodyWavefunction constructor
  
    //ManyBodyWavefunction() = delete;
    // Different type
     //template <typename MatsU> 
    // ManyBodyWavefunction(const ManyBodyWavefunction<MatsU,IntsT> & );
    // template <typename MatsU> 
    // ManyBodyWavefunction(ManyBodyWavefunction<MatsU,IntsT> &&     );

    // Same type
    ManyBodyWavefunction(const ManyBodyWavefunction<MatsT,IntsT> &);
    ManyBodyWavefunction(ManyBodyWavefunction<MatsT,IntsT> &&);     
    
    virtual ~ManyBodyWavefunction() = default;

    
    //WaveFunctionBase & referenceWaveFunction() { return dynamic_cast<WaveFunctionBase&>(ref_); }
    //SingleSlater<MatsT,IntsT> & reference() const  { return ref_;} 

  }; // class ManyBodyWavefunction

}; // namespace ChronusQ



