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
#include <util/files.hpp>
#include <fields.hpp>
#include <particleintegrals/gradints.hpp>
#include <particleintegrals/onepints.hpp>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/multipoleints.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/incoreasymmritpi.hpp>

namespace ChronusQ {

  enum TPI_TRANSFORMATION_ALG {
      INCORE_N6 = 0, // hack thru ss.formfock
      DIRECT_N6 = 1, // hack thru ss.formfock
      INCORE_N5 = 2,
      DIRECT_N5 = 3, // NYI
  };

  enum CONTRACTION_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  }; ///< 2-e Integral Contraction Algorithm


  /**
   *  \brief class to store a collection of integrals 
   *
   */
  struct IntegralsCollection {

    std::map<std::string, std::shared_ptr<ParticleIntegrals>> integrals;
    
    // Default copy and move ctors
    IntegralsCollection( const IntegralsCollection & ) = default;
    IntegralsCollection( IntegralsCollection && )      = default;
    IntegralsCollection& operator=(const IntegralsCollection&) = default;
    IntegralsCollection() = default;
    
    template <template <typename> class AddT, typename IntsT>
    void addIntegral(const std::string & name, std::shared_ptr<AddT<IntsT>> integral) {
  
      if (this->integrals.find(name) != integrals.end()) { 
        this->integrals.erase(name);
      }
     
      try {
        this->integrals.emplace(name, std::dynamic_pointer_cast<ParticleIntegrals>(integral));
      } catch (...) {
        CErr(  name + " is not an ParticleIntegrals"); 
      }
     
    }; // add to ints
  
    template <template <typename> class GetT, typename IntsT>
    std::shared_ptr<GetT<IntsT>> getIntegral(const std::string & name) const {
      
      if(this->integrals.find(name) == integrals.end()) return nullptr;
      
	  std::shared_ptr<GetT<IntsT>> integral;
    
      try {
        integral = std::dynamic_pointer_cast<GetT<IntsT>>(this->integrals.at(name));
      } catch (...) {
        return nullptr;
      }
    
      return integral; 
    }; // get from ints
    
    void erase(const std::string & name) { 
      try { integrals.erase(name); } catch (...) {}
    };
    void clear() { integrals.clear(); } 
    
    void output(std::ostream &out, const std::string &s = "",
      bool printFull = false) const {
      
      std::string output_str;
      if (s == "") output_str = "* Integrals Collection: ";
      else output_str = "* Integrals Collection in " + s + ": ";

      out << output_str << std::endl;
      for(auto const & iter: integrals) {
        iter.second->output(out, iter.first, printFull); 
      }
    };
  }; // struct IntegralsCollection
  
  /**
   *  \brief Abstract Base class for AOIntegrals
   *
   *  Stores type independent members and interfaces for templated the
   *  AOIntegrals class
   *
   */
  struct IntegralsBase {
    
    SafeFile savFile; ///< Hard storage of integrals
    HamiltonianOptions options_;
    TPI_TRANSFORMATION_ALG TPITransAlg = TPI_TRANSFORMATION_ALG::DIRECT_N6;

    // Default copy and move ctors
    IntegralsBase( const IntegralsBase & ) = default;
    IntegralsBase( IntegralsBase && )      = default;
    IntegralsBase& operator=(const IntegralsBase&) = default;
    IntegralsBase() = default;

    // Interfaces
    virtual void computeAOOneP(Molecule &mol,
        BasisSet &basis, EMPerturbation&,
        const std::vector<std::pair<OPERATOR,size_t>>&,
        const HamiltonianOptions&) = 0;

    virtual void computeAOTwoE(BasisSet&, Molecule&, EMPerturbation&) = 0;
    virtual void computeAOTwoE(BasisSet&, BasisSet&, Molecule&, EMPerturbation&) = 0;

    virtual void computeGradInts(Molecule&, BasisSet&,
        EMPerturbation&, const std::vector<std::pair<OPERATOR,size_t>>&, const
        HamiltonianOptions&) = 0;

    // Print (see src/aointegrals/print.cxx for docs)
    template <typename G> 
      friend std::ostream & operator<<(std::ostream &, const IntegralsBase& );

    virtual ~IntegralsBase() {}

  };

  /**
   *  \brief Templated class to handle the evaluation and storage of 
   *  integral matrices representing quantum mechanical operators in
   *  a finite basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both real- and complex-valued basis sets 
   *  (e.g., GTO and GIAO)
   *
   *  Real-valued arithmetics are kept in AOIntegrals
   */
  template <typename IntsT>
  class Integrals : public IntegralsBase {

  public:

    // 1-particle storage
    std::shared_ptr<OnePInts<IntsT>> overlap   = nullptr;   ///< Overlap matrix
    std::shared_ptr<OnePInts<IntsT>> kinetic   = nullptr;   ///< Kinetic matrix
    std::shared_ptr<OnePInts<IntsT>> potential = nullptr;   ///< Nuclear potential matrix

    std::shared_ptr<MultipoleInts<IntsT>> lenElectric = nullptr;
    std::shared_ptr<MultipoleInts<IntsT>> velElectric = nullptr;
    std::shared_ptr<MultipoleInts<IntsT>> magnetic = nullptr;

    // 4-component dipole storage
    std::shared_ptr<std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>>> lenElectric4C = nullptr;

    // 2-particle storage
    std::shared_ptr<TwoPInts<IntsT>> TPI = nullptr;

    // Gradient storage
    std::shared_ptr<GradInts<OnePInts,IntsT>> gradOverlap = nullptr;
    std::shared_ptr<GradInts<OnePInts,IntsT>> gradKinetic = nullptr;
    std::shared_ptr<GradInts<OnePInts,IntsT>> gradPotential = nullptr;
    std::shared_ptr<GradInts<OnePInts,IntsT>> S0a = nullptr;

    std::shared_ptr<GradInts<TwoPInts,IntsT>> gradERI = nullptr;

    // GIAO integrals
    std::shared_ptr<VectorInts<IntsT>> PVrprVP = nullptr;
    std::shared_ptr<VectorInts<IntsT>> PVrmrVP = nullptr;
    std::shared_ptr<VectorInts<IntsT>> rVr = nullptr; 

    // miscellaneous storage
    IntegralsCollection misc;

    // Constructors
    Integrals() = default;

    Integrals(const Integrals &) = default; // Copy constructor
    Integrals(Integrals &&)      = default; // Move constructor
    Integrals& operator=(const Integrals&) = default;

    // Destructor.
    ~Integrals() {}

    // Integral evaluation
    // Evaluate the 1-particle ints (general)
    virtual void computeAOOneP(Molecule &mol,
        BasisSet &basis, EMPerturbation&,
        const std::vector<std::pair<OPERATOR,size_t>>&,
        const HamiltonianOptions&);

    virtual void computeAOTwoE(BasisSet& basis, Molecule& mol,
      EMPerturbation& emPert) {
      TPI->computeAOInts(basis, mol, emPert, ELECTRON_REPULSION,
                         options_);
    }

    virtual void computeAOTwoE(BasisSet& basis, BasisSet& basis2, 
                               Molecule& mol, EMPerturbation& emPert) {
      TPI->computeAOInts(basis, basis2, mol, emPert, EP_ATTRACTION,
                         options_);
    }

    virtual void computeGradInts(Molecule&, BasisSet&,
        EMPerturbation&, const std::vector<std::pair<OPERATOR,size_t>>&, const
        HamiltonianOptions&);

    template <typename MatsT>
    Integrals<typename std::conditional<
    (std::is_same<IntsT, dcomplex>::value or
     std::is_same<MatsT, dcomplex>::value),
    dcomplex, double>::type> transform(
        const std::vector<OPERATOR>&, const std::vector<std::string>&,
        char TRANS, const MatsT* T, int NT, int LDT) const;

  }; // class Integrals


  // A struct that stores basic integral options, common for incore/direct/RI/CD)
  struct BasicIntsOptions {
    std::string ALG = "DIRECT";  ///< Parse integral algorithm
    std::string TPITRANSALG = "N6"; ///< TPI AO to MO transformation algorithm
    CONTRACTION_ALGORITHM contrAlg = CONTRACTION_ALGORITHM::DIRECT; ///< Alg for 2-body contraction
    double threshSchwarz = 1e-12; ///< Schwarz screening threshold
    std::string RI = "FALSE"; ///< RI algorithm
  };


  // A struct that stores integral options that are specific for RI/CD)
  struct CDRIIntsOptions {
    CHOLESKY_ALG CDalg = CHOLESKY_ALG::DYNAMIC_ERI; ///< Cholesky algorithm
    double CDRI_thresh = 1e-4; ///< Cholesky RI threshold
    double CDRI_sigma = 1e-2; ///< Cholesky RI sigma for span factor algorithm
    bool CDRI_genContr = true; ///< Cholesky RI uncontract basis functions to primitives
    size_t CDRI_max_qual = 1000; ///< Cholesky RI max # of qualified candidates per iteration for span-factor algorithm
    size_t CDRI_minShrinkCycle = 10; ///< Cholesky RI min # of iterations between shrinks for dynamic-all algorithm
    bool CDRI_build4I = false; ///< Cholesky RI explicitly build 4-index
    ASYMM_CD_ALG CDRI_asymmCDalg = ASYMM_CD_ALG::AUTO; ///< NEO Cholesky algorithm for approximating (ee|pp) integral
    bool CDRI_reportError = false; /// < Whether to report the error of approximate (ee|pp) comparing with exact 4-index (ee|pp) 
    bool CDRI_combineBasisTruncate = false;
    double CDRI_combineBasisThresh = 0.0;
  };


  struct IntegralOptions {
    BasicIntsOptions basicintsoptions;
    CDRIIntsOptions cdriintsoptions;
    
    // Set TPITRANSALG for the integral object
    void setTPITransAlg(std::shared_ptr<IntegralsBase> ints) const;
    
    // Build either an symmetric IntegralBase object ( either (ee|ee) or a (pp|pp) )
    std::shared_ptr<IntegralsBase> buildSymmIntegral(std::ostream &out, Molecule &mol, std::shared_ptr<BasisSet> basis,  
        std::shared_ptr<BasisSet> dfbasis, std::string s) const;

    // Build either an asymmetric IntegralBase object ( (ee|pp) )
    std::shared_ptr<IntegralsBase> buildAsymmIntegral(std::ostream &out, Molecule &mol, std::shared_ptr<BasisSet> basis,  
        std::shared_ptr<BasisSet> dfbasis, std::shared_ptr<BasisSet> basis2, IntegralOptions eopts, IntegralOptions popts, 
        std::shared_ptr<IntegralsBase> aoi, std::shared_ptr<IntegralsBase> paoi) const;
    
    // Build all (ee|ee), (pp|pp), (ee|pp) objects (if needed), and return them in a tuple
    static std::tuple<std::shared_ptr<IntegralsBase>, std::shared_ptr<IntegralsBase>, std::shared_ptr<IntegralsBase>> buildAllIntegrals(
        std::ostream &out, Molecule &mol, std::shared_ptr<BasisSet> basis,  std::shared_ptr<BasisSet> dfbasis, 
        std::shared_ptr<BasisSet> basis2, IntegralOptions eopts, IntegralOptions popts, IntegralOptions epopts);
  };


  // Declear function to parse the user-specific integral options
  IntegralOptions getIntegralOptions(std::ostream &, CQInputFile &, 
      std::shared_ptr<BasisSet>,  std::shared_ptr<BasisSet>, 
      std::shared_ptr<BasisSet>, std::string);

}; // namespace ChronusQ
