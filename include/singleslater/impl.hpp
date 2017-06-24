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

#include <singleslater.hpp>
#include <util/preprocessor.hpp>
#include <quantum/preprocessor.hpp>
#include <corehbuilder/impl.hpp>
#include <particleintegrals/twopints/impl.hpp>
#include <fockbuilder/rofock.hpp>
#include <fockbuilder/fourcompfock.hpp>

namespace ChronusQ {

  /**
   *  Constructs a SingleSlater object from another of a another (possibly the 
   *  same) type by copy.
   *
   *  \param [in] other SingleSlater object to copy
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to copy constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  SingleSlater<MatsT,IntsT>::SingleSlater(const SingleSlater<MatsU,IntsT> &other,int dummy) : 
    //orthoType(other.orthoType),coreType(other.coreType),
    QuantumBase(dynamic_cast<const QuantumBase&>(other)),
    WaveFunctionBase(dynamic_cast<const WaveFunctionBase&>(other)),
    SingleSlaterBase(dynamic_cast<const SingleSlaterBase&>(other)),
    WaveFunction<MatsT,IntsT>(dynamic_cast<const WaveFunction<MatsU,IntsT>&>(other)),
    //basisSet_(other.basisSet_),
    fockMatrix(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.fockMatrix)),
    fockMatrixOrtho(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.fockMatrixOrtho)),
    coulombMatrix(std::make_shared<cqmatrix::Matrix<MatsT>>(*other.coulombMatrix)),
    exchangeMatrix(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.exchangeMatrix)),
    twoeH(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.twoeH)),
    onePDMOrtho(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.onePDMOrtho)),
    deltaOnePDM(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.deltaOnePDM)),
    //coreH(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.coreH)),
    coreH(other.coreH ? std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.coreH) : nullptr),
    coreHPerturbed(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*other.coreHPerturbed)),
    TPI(TPIContractions<MatsU,IntsT>::template convert<MatsT>(other.TPI)),
    coreHBuilder(CoreHBuilder<MatsU,IntsT>::template convert<MatsT>(other.coreHBuilder)),
    fockBuilder(FockBuilder<MatsU,IntsT>::template convert<MatsT>(other.fockBuilder)) {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater<MatsT>::SingleSlater(const SingleSlater<U>&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif
    if( other.orthoSpinor )
      orthoSpinor = std::make_shared<Orthogonalization<MatsT>>(*other.orthoSpinor);
    if( other.orthoAB )
      orthoAB = std::make_shared<Orthogonalization<MatsT>>(*other.orthoAB);

  }; // SingleSlater<MatsT>::SingleSlater(const SingleSlater<U> &)



  /**
   *  Constructs a SingleSlater object from another of a another (possibly the 
   *  same) by move.
   *
   *  \warning Deallocates the passed SingleSlater object
   *
   *  \param [in] other SingleSlater object to move
   *  \param [in] dummy Dummy argument to fix calling signature for delegation 
   *    to move constructor
   */ 
  template <typename MatsT, typename IntsT>
  template <typename MatsU> 
  SingleSlater<MatsT,IntsT>::SingleSlater(SingleSlater<MatsU,IntsT> &&other,int dummy) : 
    //orthoType(other.orthoType),coreType(other.coreType),
    QuantumBase(dynamic_cast<QuantumBase&&>(std::move(other))),
    WaveFunctionBase(dynamic_cast<WaveFunctionBase&&>(std::move(other))),
    SingleSlaterBase(dynamic_cast<SingleSlaterBase&&>(std::move(other))),
    WaveFunction<MatsT,IntsT>(dynamic_cast<WaveFunction<MatsU,IntsT>&&>(std::move(other))),
    //basisSet_(other.basisSet_),
    fockMatrix(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.fockMatrix))),
    fockMatrixOrtho(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.fockMatrixOrtho))),
    coulombMatrix(std::make_shared<cqmatrix::Matrix<MatsT>>(std::move(*other.coulombMatrix))),
    exchangeMatrix(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.exchangeMatrix))),
    twoeH(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.twoeH))),
    onePDMOrtho(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.onePDMOrtho))),
    deltaOnePDM(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.deltaOnePDM))),
    coreH(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.coreH))),
    coreHPerturbed(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(std::move(*other.coreHPerturbed))),
    TPI(TPIContractions<MatsU,IntsT>::template convert<MatsT>(other.TPI)),
    coreHBuilder(CoreHBuilder<MatsU,IntsT>::template convert<MatsT>(other.coreHBuilder)),
    fockBuilder(FockBuilder<MatsU,IntsT>::template convert<MatsT>(other.fockBuilder)) {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater<MatsT>::SingleSlater(SingleSlater<U>&&) "
              << "(this = " << this << ", other = " << &other << ")" 
              << std::endl;
#endif

    if( other.orthoSpinor )
      orthoSpinor = std::make_shared<Orthogonalization<MatsT>>(std::move(*other.orthoSpinor));
    if( other.orthoAB )
      orthoAB = std::make_shared<Orthogonalization<MatsT>>(std::move(*other.orthoAB));

  }; // SingleSlater<MatsT>::SingleSlater(SingleSlater<U> &&)


  // Delagate the copy constructor to the conversion constructors
  template <typename MatsT, typename IntsT>
  SingleSlater<MatsT,IntsT>::SingleSlater(const SingleSlater<MatsT,IntsT> &other) : 
    SingleSlater(other,0){ };
  template <typename MatsT, typename IntsT>
  SingleSlater<MatsT,IntsT>::SingleSlater(SingleSlater<MatsT,IntsT> &&other) : 
    SingleSlater(std::move(other),0){ };


  /**
   *  Allocates the internal memory a SingleSlater object
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::alloc() {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater::alloc (this = " << this << ")" << std::endl;
#endif

    size_t NB = this->basisSet().nBasis;
    moPairs.resize(2, {});
    
    if( nC != 4 ) {
      SPIN_OPERATOR_ALLOC(NB,fockMatrix);
      SPIN_OPERATOR_ALLOC(NB,fockMatrixOrtho);
      SPIN_OPERATOR_ALLOC(NB,onePDMOrtho);
      SPIN_OPERATOR_ALLOC(NB,deltaOnePDM);
      SPIN_OPERATOR_ALLOC(NB,coreHPerturbed);


      SPIN_OPERATOR_ALLOC(NB,exchangeMatrix);
      SPIN_OPERATOR_ALLOC(NB,twoeH);

      SPIN_OPERATOR_ALLOC(NB,gauntexchangeMatrix);
      SPIN_OPERATOR_ALLOC(NB,gaunttwoeH);
      SPIN_OPERATOR_ALLOC(NB,gaugeexchangeMatrix);
      SPIN_OPERATOR_ALLOC(NB,gaugetwoeH);

      coulombMatrix = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);
    } else {

      SPIN_OPERATOR_ALLOC(2*NB,fockMatrix);
      SPIN_OPERATOR_ALLOC(2*NB,fockMatrixOrtho);
      SPIN_OPERATOR_ALLOC(2*NB,onePDMOrtho);
      SPIN_OPERATOR_ALLOC(2*NB,deltaOnePDM);
      SPIN_OPERATOR_ALLOC(2*NB,coreHPerturbed);

      SPIN_OPERATOR_ALLOC(2*NB,exchangeMatrix);
      SPIN_OPERATOR_ALLOC(2*NB,twoeH);
      SPIN_OPERATOR_ALLOC(2*NB,gauntexchangeMatrix);
      SPIN_OPERATOR_ALLOC(2*NB,gaunttwoeH);
      SPIN_OPERATOR_ALLOC(2*NB,gaugeexchangeMatrix);
      SPIN_OPERATOR_ALLOC(2*NB,gaugetwoeH);

      coulombMatrix = std::make_shared<cqmatrix::Matrix<MatsT>>(2*NB);
    }

  }; // SingleSlater<MatsT>::alloc


  /**
   *  Deallocates the internal memory a SingleSlater object
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::dealloc() {

#ifdef _SingleSlaterDebug
    std::cout << "SingleSlater::dealloc (this = " << this << ")" << std::endl;
#endif

    fockMatrix = nullptr;
    fockMatrixOrtho = nullptr;
    onePDMOrtho = nullptr;
    coreHPerturbed = nullptr;

    exchangeMatrix = nullptr;
    twoeH = nullptr;

   coulombMatrix = nullptr;

  }; // SingleSlater<MatsT>::dealloc

  template <typename MatsT, typename IntsT>
  class NEOSS;

  /**
   *  \brief The pointer convertor. This static function converts
   *  the underlying polymorphism correctly to hold a different
   *  type of matrices.
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  std::shared_ptr<SingleSlater<MatsU,IntsT>>
  SingleSlater<MatsT,IntsT>::convert(const std::shared_ptr<SingleSlater<MatsT,IntsT>>& ss) {

    if (not ss) return nullptr;

    const std::type_info &tID(typeid(*ss));

    if (tID == typeid(NEOSS<MatsT,IntsT>)) {
      return std::make_shared<NEOSS<MatsU,IntsT>>(
          *std::dynamic_pointer_cast<NEOSS<MatsT,IntsT>>(ss));

    } else if (tID == typeid(HartreeFock<MatsT,IntsT>)) {
      return std::make_shared<HartreeFock<MatsU,IntsT>>(
          *std::dynamic_pointer_cast<HartreeFock<MatsT,IntsT>>(ss));

    } else if (tID == typeid(KohnSham<MatsT,IntsT>)) {
      return std::make_shared<KohnSham<MatsU,IntsT>>(
          *std::dynamic_pointer_cast<KohnSham<MatsT,IntsT>>(ss));

    } else {
      std::stringstream errMsg;
      errMsg << "SingleSlater implementation \"" << tID.name() << "\" not registered in convert." << std::endl;
      CErr(errMsg.str(),std::cout);
    }

    return nullptr;

  }
 
}; // namespace ChronusQ


// Other implementation files
#include <singleslater/quantum.hpp>   // Quantum declarations
#include <singleslater/fock.hpp>      // Fock matrix header
#include <singleslater/rt.hpp>        // RT header
#include <singleslater/guess.hpp>     // Guess header
#include <singleslater/scf.hpp>       // SCF header
#include <singleslater/print.hpp>     // Print header
#include <singleslater/pop.hpp>       // Population analysis
#include <singleslater/fchk.hpp>      // Fchk-specific header
#include <singleslater/cube.hpp>     // Cubegen header


#include <singleslater/kohnsham/impl.hpp> // KS headers
#include <singleslater/kohnsham/fxc.hpp> // KS headers
#include <singleslater/kohnsham/scf.hpp> // Newton-Raphson functions

#include <singleslater/hartreefock/scf.hpp> // Newton-Raphson Functions

#include <orbitalmodifier/impl.hpp> // OrbitalModifier Implementation headers
