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

namespace ChronusQ {

  /**
   *  \breif The Hartree--Fock class.
   *
   *  Trivially specializes the SingleSlater class for a Hartree--Fock description of the
   *  many-body wave function
   */ 
  template <typename MatsT, typename IntsT>
  class HartreeFock : virtual public SingleSlater<MatsT,IntsT>,
    public std::enable_shared_from_this<HartreeFock<MatsT,IntsT>> {

    public:

    std::shared_ptr<HartreeFock<MatsT,IntsT>> getPtr(){ return this->shared_from_this(); }

    // Trivially inherit ctors from SingleSlater<T>

    template <typename... Args>
    HartreeFock(MPI_Comm c, Molecule &mol, BasisSet &basis,
                std::shared_ptr<Integrals<IntsT>> aoi, Args... args) :
      SingleSlater<MatsT,IntsT>(c,mol,basis,aoi,args...),
      WaveFunctionBase(c,mol,basis,args...),
      QuantumBase(c,args...) {

      // Append HF tags to reference names
      if(this->nC == 1) {
        if(this->iCS) {
          this->refLongName_  += "Restricted Hartree-Fock";
          this->refShortName_ += "RHF";
        } else {
          this->refLongName_  += "Unrestricted Hartree-Fock";
          this->refShortName_ += "UHF";
        }
      } else {
        this->refLongName_  += "Generalized Hartree-Fock";
        this->refShortName_ += "GHF";
      }

    }; // HartreeFock constructor

    // Allow for reference name specification
    template <typename... Args>
    HartreeFock(std::string rL, std::string rS, MPI_Comm c,
                Molecule &mol, BasisSet &basis,
                std::shared_ptr<Integrals<IntsT>> aoi, Args... args) :
      SingleSlater<MatsT,IntsT>(c,mol,basis,aoi,args...),
      WaveFunctionBase(c,mol,basis,args...),
      QuantumBase(c,args...) {

      this->refLongName_  = rL;
      this->refShortName_ = rS;
      

    }; // HartreeFock constructor (with strings)


    // Copy and Move ctors

    template <typename MatsU, typename IntsU> 
    HartreeFock(const HartreeFock<MatsU,IntsU> &other, int dummy = 0) :
      SingleSlater<MatsT,IntsT>(dynamic_cast<const SingleSlater<MatsU,IntsU>&>(other),dummy),
      QuantumBase(dynamic_cast<const QuantumBase&>(other)),
      WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other))
      { };

    template <typename MatsU, typename IntsU> 
    HartreeFock(HartreeFock<MatsU,IntsU> &&other, int dummy = 0) :
      SingleSlater<MatsT,IntsT>(dynamic_cast<SingleSlater<MatsU,IntsU>&&>(other),dummy),
      QuantumBase(dynamic_cast<QuantumBase&&>(other)),
      WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(other))
      { };

    HartreeFock(const HartreeFock<MatsT,IntsT> &other) :
      SingleSlater<MatsT,IntsT>(dynamic_cast<const SingleSlater<MatsT,IntsT>&>(other),0),
      QuantumBase(dynamic_cast<const QuantumBase&>(other)),
      WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other))
      { }

    HartreeFock(HartreeFock<MatsT,IntsT> &&other) :
      SingleSlater<MatsT,IntsT>(dynamic_cast<SingleSlater<MatsT,IntsT>&&>(other),0),
      QuantumBase(dynamic_cast<QuantumBase&&>(other)),
      WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(other))
      { }


    using QuantumBase::computeEnergy;
    //void getNRCoeffs(MatsT*);
    void computeFullNRStep(MatsT*);
    std::pair<double,MatsT*> getStab();
    void buildOrbitalModifierOptions();


  }; // class HartreeFock

}; // namespace ChronusQ

