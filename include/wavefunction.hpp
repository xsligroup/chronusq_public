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
#include <quantum.hpp>
#include <wavefunction/base.hpp>
#include <integrals.hpp>
#include <matrix.hpp>

// Debug print triggered by Quantum
  
#ifdef _QuantumDebug
  #define _WaveFunctionDebug
#endif

namespace ChronusQ {


  /**
   *  \brief The WaveFunction class. The typed abstract interface for
   *  all classes which admit a well defined wave function (HF, KS, etc).
   *
   *  Adds knowledge of storage type to WaveFunctionBase
   *
   *  Specializes the Quantum class of the same type; 
   */
  template <typename MatsT, typename IntsT>
  class WaveFunction : virtual public WaveFunctionBase, public Quantum<MatsT> {

    template <typename MatsU, typename IntsU>
    friend class WaveFunction;

  protected:

    // Useful typedefs
    typedef MatsT*               oper_t;
    typedef std::vector<oper_t>  oper_t_coll;

  public:

    std::shared_ptr<Integrals<IntsT>> aoints_; ///< AOIntegrals for the storage of integrals

    // Operator storage

    ///< mo[0] : Full (nC > 1) / ALPHA (nC == 1) MO coefficient matrix
    ///< mo[1] : BETA (nC == 1) MO coefficient matrix
    /// The order of AO Basis in mo[0]:
    ///    * 2C: [Alpha, Beta]
    ///    * 4C: [Alpha Large, Alpha Small, Beta Large, Beta Small]
    std::vector<cqmatrix::Matrix<MatsT>> mo;
    double* eps1; ///< Full (nC > 1) / ALPHA (nC == 1) Fock eigenvalues
    double* eps2; ///< BETA (nC == 1) Fock eigenvalues

    // Constructors

    // Disable default constructor
    WaveFunction() = delete;

    /**
     *  WaveFunction Constructor. Constructs a WaveFunction object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] _nC  Number of spin components (1 and 2 are supported)
     *  \param [in] iCS  Whether or not to treat as closed shell
     */ 
    WaveFunction(MPI_Comm c, Molecule &mol, BasisSet &basis,
                 std::shared_ptr<Integrals<IntsT>> aoi, size_t _nC, bool iCS, Particle p = {-1.0,1.0}) :
      QuantumBase(c,_nC,iCS,p),
      WaveFunctionBase(c, mol, basis,_nC,iCS,p),
      Quantum<MatsT>(c,_nC,iCS,p,basis.nBasis),
      //molecule_(mol), basisSet_(basis), 
      eps1(nullptr), eps2(nullptr), aoints_(aoi) {

      // Compute meta data
      size_t nBasis = basis.nBasis;

      if (p.charge < 0.) {
        this->nO = molecule_.nTotalE;
        this->nV = 2*nBasis - nO;

        if( this->iCS ) {
          this->nOA = this->nO / 2; this->nOB = this->nO / 2;
          this->nVA = this->nV / 2; this->nVB = this->nV / 2;
        } else {
          size_t nSingleE = molecule_.multip - 1;
          this->nOB = (this->nO - nSingleE) / 2;
          this->nOA = this->nOB + nSingleE;
          this->nVA = nBasis - this->nOA;
          this->nVB = nBasis - this->nOB;
        }
      } else {
         this->nO = molecule_.nTotalP;
         this->nV = 2*nBasis - nO;
 
         if( this->iCS ) {
           this->nOA = this->nO / 2; this->nOB = this->nO / 2;
           this->nVA = this->nV / 2; this->nVB = this->nV / 2;
         } else {
           size_t nSingleP = molecule_.multip_proton - 1;
           this->nOB = (this->nO - nSingleP) / 2;
           this->nOA = this->nOB + nSingleP;
           this->nVA = nBasis - this->nOA;
           this->nVB = nBasis - this->nOB;
         }

      }

      // Allocate internal memory
      if(nBasis != 0) alloc();

    };

    // See include/wavefunction/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename MatsU> 
      WaveFunction(const WaveFunction<MatsU,IntsT> &, int dummy = 0);
    template <typename MatsU> 
      WaveFunction(WaveFunction<MatsU,IntsT> &&     , int dummy = 0);

    // Same type
    WaveFunction(const WaveFunction &);
    WaveFunction(WaveFunction &&);     
    

    /**
     *  Deconstructor
     */ 
    ~WaveFunction(){ dealloc(); }


    // Deallocation (see include/wavefunction/impl.hpp for docs)
    void alloc();
    void dealloc();

    // Print Functions
    void printMO(std::ostream&) ;
    virtual void printEPS(std::ostream&);
    virtual void printMOInfo(std::ostream&, size_t a = 0);

    // Swap Function
    void swapMOs(std::vector<std::vector<std::pair<size_t, size_t>>>&, SpinType sp);

  }; // class WaveFunction

}; // namespace ChronusQ

