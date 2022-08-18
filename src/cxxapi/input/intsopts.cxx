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
#include <cxxapi/options.hpp>
#include <cerr.hpp>

#include <particleintegrals/print.hpp>
#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/giaodirecteri.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>

namespace ChronusQ {

  enum class CONTRACTION_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  }; ///< 2-e Integral Contraction Algorithm

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQINTS_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "ALG",          // Direct or Incore?
      "GRADALG",      // Direct or Incore for gradients?
      "TPITRANSALG",   // N5 or N6
      "SCHWARZ",     // double
      "RI",           // AUXBASIS or CHOLESKY or False
      "RITHRESHOLD",  // double
      "RISIGMA",      // double
      "RIMINSHRINK",  // size_t
      "RIMAXQUAL",    // size_t
      "RIGENCONTR",   // True or False
      "RIBUILD4INDEX",// True or False
      "FINITENUCLEI", // True or False
      "LIBCINT"       // Ture or False

      "BARECOULOMB",  // True or False
      "DC",           // True or False, SF, SD, 3C, 2C, 1C, AMF
      "DIRACCOULOMB", // True or False
      "BREIT",        // True or False
      "GAUNT",        // True or False
      "SSSS",         // True or False
      "GAUGE",        // True or False

   };

    // Specified keywords
    std::vector<std::string> intsKeywords = input.getDataInSection("INTS");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : intsKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword INTS." + keyword + " is not recognized",std::cout);// Error
    }

    // Specified NEO keywords
    intsKeywords = input.getDataInSection("PINTS");

    // Make sure all of the basisKeywords in allowedKeywords
    for( auto &keyword : intsKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if ( ipos == allowedKeywords.end() )
        CErr("Keyword PINTS." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  \brief Optionally set the control parameters for an
   *  AOIntegrals object
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object 
   *
   *
   */ 
  std::shared_ptr<IntegralsBase> CQIntsOptions(std::ostream &out, 
      CQInputFile &input, CQMemManager &mem, Molecule &mol,
      std::shared_ptr<BasisSet> basis,  std::shared_ptr<BasisSet> dfbasis, 
      std::shared_ptr<BasisSet> basis2, std::string int_sec) {

    // check the validity of the integral option section 
    if (not int_sec.compare("INTS") and not int_sec.compare("PINTS") and not int_sec.compare("EPINTS"))
      CErr("Found invalue integral section");

    // Parse integral algorithm
    std::string ALG = "DIRECT";
    OPTOPT( ALG = input.getData<std::string>(int_sec+".ALG"); )
    trim(ALG);
    
    std::string TPITRANSALG = "N6"; 
    OPTOPT( TPITRANSALG = input.getData<std::string>(int_sec+".TPITRANSALG"); )
    trim(TPITRANSALG);

    // Control Variables
    CONTRACTION_ALGORITHM contrAlg = CONTRACTION_ALGORITHM::DIRECT; ///< Alg for 2-body contraction
    double threshSchwarz = 1e-12; ///< Schwarz screening threshold
    std::string RI = "FALSE"; ///< RI algorithm
    CHOLESKY_ALG CDalg = CHOLESKY_ALG::DYNAMIC_ERI; ///< Cholesky algorithm
    double CDRI_thresh = 1e-4; ///< Cholesky RI threshold
    double CDRI_sigma = 1e-2; ///< Cholesky RI sigma for span factor algorithm
    bool CDRI_genContr = true; ///< Cholesky RI uncontract basis functions to primitives
    size_t CDRI_max_qual = 1000; ///< Cholesky RI max # of qualified candidates per iteration for span-factor algorithm
    size_t CDRI_minShrinkCycle = 10; ///< Cholesky RI min # of iterations between shrinks for dynamic-all algorithm
    bool CDRI_build4I = false; ///< Cholesky RI explicitly build 4-index

    if( not ALG.compare("DIRECT") )
      contrAlg = CONTRACTION_ALGORITHM::DIRECT;
    else if( not ALG.compare("INCORE") )
      contrAlg = CONTRACTION_ALGORITHM::INCORE;
    else
      CErr(ALG + " not a valid INTS.ALG",out);

    // Parse Schwarz threshold
    OPTOPT( threshSchwarz = input.getData<double>(int_sec+".SCHWARZ"); )

    // Parse RI option
    OPTOPT( RI = input.getData<std::string>(int_sec+".RI");)
    trim(RI);

    if(RI.compare("FALSE")) {
      if (contrAlg != CONTRACTION_ALGORITHM::INCORE) {
        contrAlg = CONTRACTION_ALGORITHM::INCORE;
        std::cout << "Incore ERI algorithm enforced by RI." << std::endl;
      }
      if(not RI.compare("AUXBASIS") ) {
        if (dfbasis->nBasis < 1)
          CErr("Keyword INTS.RI requires a non-empty DFbasis->",std::cout);
      } else if (not RI.compare("TRADITIONAL")) {
        CDalg = CHOLESKY_ALG::TRADITIONAL;
      } else if (not RI.compare("DYNAMICALL")) {
        CDalg = CHOLESKY_ALG::DYNAMIC_ALL;
      } else if (not RI.compare("SPANFACTOR")) {
        CDalg = CHOLESKY_ALG::SPAN_FACTOR;
      } else if (not RI.compare("CHOLESKY") or not RI.compare("DYNAMICERI")) {
        CDalg = CHOLESKY_ALG::DYNAMIC_ERI;
      } else if (not RI.compare("SPANFACTORREUSE")) {
        CDalg = CHOLESKY_ALG::SPAN_FACTOR_REUSE;
      } else {
        CErr(RI + " not a valid INTS.RI",out);
      }
    }

    OPTOPT( CDRI_genContr = input.getData<bool>("INTS.RIGENCONTR"); )
    OPTOPT( CDRI_thresh = input.getData<double>("INTS.RITHRESHOLD"); )
    OPTOPT( CDRI_sigma = input.getData<double>("INTS.RISIGMA"); )
    OPTOPT( CDRI_max_qual = input.getData<size_t>("INTS.RIMAXQUAL"); )
    OPTOPT( CDRI_minShrinkCycle = input.getData<size_t>("INTS.RIMINSHRINK"); )
    OPTOPT( CDRI_build4I = input.getData<bool>("INTS.RIBUILD4INDEX"); )

    std::shared_ptr<IntegralsBase> aoi = nullptr;

    if(basis->basisType == REAL_GTO) {
      std::shared_ptr<Integrals<double>> aoint =
          std::make_shared<Integrals<double>>();

      if(RI.compare("FALSE")) {
        if (basis2)
          CErr("AUXBASIS or CHOLESKY with NEO NYI");
        if(not RI.compare("AUXBASIS"))
          aoint->TPI =
              std::make_shared<InCoreAuxBasisRIERI<double>>(mem,basis->nBasis,dfbasis);
        else
          aoint->TPI =
              std::make_shared<InCoreCholeskyRIERI<double>>(
                  mem, basis->nBasis, CDRI_thresh, CDalg, CDRI_genContr,
                  CDRI_sigma, CDRI_max_qual, CDRI_minShrinkCycle, CDRI_build4I);
      } else if (contrAlg == CONTRACTION_ALGORITHM::INCORE) {
        if (not basis2)
          aoint->TPI =
              std::make_shared<InCore4indexTPI<double>>(mem,basis->nBasis);
        else
          aoint->TPI = 
              std::make_shared<InCore4indexTPI<double>>(mem,basis->nBasis,basis2->nBasis);
      }
      else {
        if (not basis2)
          aoint->TPI =
              std::make_shared<DirectTPI<double>>(mem,*basis,*basis,mol,threshSchwarz);
        else
          aoint->TPI = 
              std::make_shared<DirectTPI<double>>(mem,*basis,*basis2,mol,threshSchwarz);
      }

      aoi = std::dynamic_pointer_cast<IntegralsBase>(aoint);
    } else if(basis->basisType == COMPLEX_GIAO) {
      std::shared_ptr<Integrals<dcomplex>> giaoint =
          std::make_shared<Integrals<dcomplex>>();
      if(RI.compare("FALSE"))
        CErr("GIAO resolution of identity ERI NYI",std::cout);
      else if (contrAlg == CONTRACTION_ALGORITHM::INCORE) {
        if (basis2)
          CErr("GIAO with NEO NYI",std::cout);
        giaoint->TPI =
            std::make_shared<InCore4indexTPI<dcomplex>>(mem,basis->nBasis);
      }
      else {
        if (basis2)
          CErr("GIAO with NEO NYI",std::cout);
        giaoint->TPI =
            std::make_shared<DirectTPI<dcomplex>>(mem,*basis,*basis,mol,threshSchwarz);
      }

      aoi = std::dynamic_pointer_cast<IntegralsBase>(giaoint);
    }
    
    if (not TPITRANSALG.compare("N5")) {
      if (contrAlg == CONTRACTION_ALGORITHM::DIRECT) {
        aoi->TPITransAlg = TPI_TRANSFORMATION_ALG::DIRECT_N5;
        CErr("DIRECT_N5 TPI Transformation is NYI");
      } else
        aoi->TPITransAlg = TPI_TRANSFORMATION_ALG::INCORE_N5;
    } else if (not TPITRANSALG.compare("N6")){
      if (contrAlg == CONTRACTION_ALGORITHM::DIRECT)
        aoi->TPITransAlg = TPI_TRANSFORMATION_ALG::DIRECT_N6;
      else
        aoi->TPITransAlg = TPI_TRANSFORMATION_ALG::INCORE_N6;
    } else {
      CErr(TPITRANSALG + " not a valid INTS.TPITRANSALG",out);
    }

    // Print
    out <<  *aoi << std::endl;

    
    return aoi;

  }; // CQIntsOptions

}; // namespace ChronusQ
