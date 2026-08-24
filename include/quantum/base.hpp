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

#include <util/mpi.hpp>
#include <util/timer.hpp>
#include <util/typedefs.hpp>

#include <cqlinalg/blas1.hpp>

#include <fields.hpp>

#include <particleintegrals.hpp>

namespace ChronusQ {

  enum DENSITY_TYPE {
    SCALAR=0,MZ=1,MY=2,MX=3
  }; ///< Enumerate the types of densities for contraction

  enum PROPERTY {
    ELECTRIC_DIPOLE,
    ELECTRIC_QUADRUPOLE, 
    ELECTRIC_OCTUPOLE,
    SPIN,
    METHOD_SPECIFIC
  };

  // Helper function for operator traces
  template <typename Left, typename Right>
  static inline double OperatorTrace(size_t N, const Left& op1 , 
    const Right& op2) {
    
    return std::real(blas::dot(N,op1,1,op2,1));
  } 

  /**
   *  \brief The QuantumBase class. The abstraction of information
   *  relating to the Quantum class which are independent of storage
   *  type.
   *
   *  See Quantum for further docs
   *
   */ 
  class QuantumBase {

  protected:
  private:
  public:

    MPI_Comm      comm; ///< MPI Communicator
    
    int   nC;   ///< Number of spin components
    bool  iCS;  ///< is closed shell?
    Particle particle; ///< Particle Type


    // Property storage

    // Length gauge electric multipoles
    cart_t elecDipole;        ///< Electric Dipole in the length gauge
    cartmat_t elecQuadrupole; ///< Electric Quadrupole in the length gauge
    cartrk3_t elecOctupole;   ///< Electric Octupole in the length gauge
    
    // Spin expectation values
    cart_t SExpect; ///< Expectation values of Sx, Sy and Sz
    double    SSq;  ///< Expectation value of S^2
    double    SQuantNum;  ///< spin quantum number(S) from (S(S+1))
    cart_t SExpectLL; ///< 4C-only: LL contribution to Sx, Sy and Sz
    cart_t SExpectSS; ///< 4C-only: SS contribution to Sx, Sy and Sz
    double SSqLL;     ///< 4C-only: LL contribution to S^2
    double SSqSS;     ///< 4C-only: SS contribution to S^2
    double SSqCross;  ///< 4C-only: LL-SS cross contribution to S^2

    // Angular momentum expectation values
    cart_t LExpect; ///< Expectation values of Lx, Ly and Lz
    double    LSq;  ///< Expectation value of L^2
    double    LQuantNum;  ///< orbital angular momentum quantum number(L) from (L(L+1))
    cart_t JExpect; ///< Expectation values of Jx, Jy and Jz
    double    JSq;  ///< Expectation value of J^2
    double    JQuantNum;  ///< total angular momentum quantum number(J) from (J(J+1))
    double    SL; ///< Expectation value of S·L
    double    LS; ///< Expectation value of L·S

    // Energy expectation values
    double OBEnergy;   ///< 1-Body operator contribution to the energy
    double MBEnergy;   ///< Many(2)-Body operator contribution to the energy
    double PPEnergy = 0.;   //<  Many(2)-Body proton-proton repulsion energy
    double extraEnergy = 0;
    double totalEnergy;///< The total energy

    // Options for orbital property analysis
    bool orbProp = false;


    // Constructors
      
    // Disable default constructor
    QuantumBase() = delete;

    // Default the copy and move constructors
    QuantumBase(const QuantumBase&) = default;
    QuantumBase(QuantumBase&&)      = default;

    /**
     *  QuantumBase Constructor. Constructs a QuantumBase object.
     *
     *  \param [in] _nC   Number of spin components (1 and 2 are supported)
     *  \param [in] _iCS  Whether or not system is closed shell
     *                    (only used when _nC == 1)
     */ 
    QuantumBase(MPI_Comm c, size_t _nC, bool _iCS, Particle p, std::optional<size_t> = std::nullopt): 
      nC(_nC), iCS(_iCS), particle(p), comm(c),
      elecDipole({0.,0.,0.}),
      elecQuadrupole{{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
      elecOctupole{
        {
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
          {{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
        }}, SExpect({0.,0.,0.}), SSq(0.),
      SExpectLL({0.,0.,0.}), SExpectSS({0.,0.,0.}), SSqLL(0.), SSqSS(0.), SSqCross(0.),
      LExpect({0.,0.,0.}), LSq(0.), JExpect({0.,0.,0.}), JSq(0.),
      SL(0.), LS(0.),
      OBEnergy(0.), MBEnergy(0.), totalEnergy(0.) { }; // Quantum::Quantum



    // Procedural
      
    /**
     *  Function to form the density. Pure virtual
     */ 
    virtual void formDensity() = 0;

    /**
     *  Function to compute the field free energy expectation value(s)
     */ 
    virtual void computeEnergy() = 0;


    /**
     *  Function to compute the energy expectation values including
     *  field terms
     */ 
    void computeEnergy(EMPerturbation &pert){

      ROOT_ONLY(comm);

      computeEnergy();


      // Field contributions to the energy
      double field_delta(0.);

      // Purely magnetic field contributions here
      auto magDipoleField = pert.getDipoleAmp(Magnetic);


      // Purely electric field contributions here
      // Only calculate the dipole during SCF iterations if there's an applied
      // electric field
      if(pert_has_type(pert,Electric))
      {
        computeMultipole(pert, {ELECTRIC_DIPOLE});

         auto elecDipoleField = pert.getDipoleAmp(Electric);
         field_delta +=
           elecDipoleField[0] * elecDipole[0] +
           elecDipoleField[1] * elecDipole[1] +
           elecDipoleField[2] * elecDipole[2];
      }

      totalEnergy += field_delta; // Increment total energy

    };


    virtual std::vector<double> getEnergySummary() {
      return {totalEnergy, OBEnergy, MBEnergy};
    }
   
    virtual void computeMultipole(EMPerturbation &, const std::vector<PROPERTY> &properties = {}) = 0;
    virtual void computeSpin() = 0;
    virtual void computeSpinAndAngularProperties() = 0;
    virtual void methodSpecificProperties() = 0;
    virtual void computeOrbitalProps() = 0;
    virtual void computeOrbitalRDFs() = 0;
    virtual void printOrbitalEnergies() = 0;
    virtual void printAngularProperties(std::ostream&, bool withBanner = true) = 0;
    virtual bool isX2CReference() const { return false; }


    inline void computeProperties(EMPerturbation &pert, const std::vector<PROPERTY> &properties = {}) {

      ROOT_ONLY(comm);

      ProgramTimer::tick("Compute Properties");
      
      // By default, compute all properties
      std::vector<PROPERTY> propsToCompute = properties;
      if (propsToCompute.empty()) 
        propsToCompute = {ELECTRIC_DIPOLE, ELECTRIC_QUADRUPOLE, ELECTRIC_OCTUPOLE, SPIN, METHOD_SPECIFIC};
      
      auto hasProperty = [&](PROPERTY p){ return std::find(propsToCompute.begin(), propsToCompute.end(), p) != propsToCompute.end(); };

      // Compute multipole properties
      std::vector<PROPERTY> multipoleProps;
      for (const auto& prop : propsToCompute) 
        if (prop == ELECTRIC_DIPOLE || prop == ELECTRIC_QUADRUPOLE || prop == ELECTRIC_OCTUPOLE) 
          multipoleProps.push_back(prop);
      
      if (multipoleProps.size() > 0)
        computeMultipole(pert, multipoleProps);

      // Compute spin information for all component spaces, including 4C.
      // 4C angular-property evaluation uses SSq from computeSpin().
      if(hasProperty(SPIN)) 
        computeSpin();


      // Compute angular momentum properties
      if (nC != 4)
        computeSpinAndAngularProperties();

      // Compute method specific properties (currently only include Mulliken population analysis)
      if(hasProperty(METHOD_SPECIFIC)) 
        methodSpecificProperties();

      if(orbProp) computeOrbitalProps();

      ProgramTimer::tock("Compute Properties");

    };
    

    // Print functions
    virtual void print1PDM(std::ostream&) = 0;
    virtual void printMultipoles(std::ostream&);
    virtual void printSpin(std::ostream&, bool withBanner = true);
    virtual void printMiscProperties(std::ostream&) = 0;
  }; // class QuantumBase

}; // namespace ChronusQ

