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
#include <particleintegrals/twopints/incoreasymmritpi.hpp>

namespace ChronusQ {

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
      "TPITRANSALG",  // N5 or N6
      "SCHWARZ",      // double
      "RI",           // String, determines which algorithm to use for RI/CD
                      // For INTS/PINTS section: "AUXBASIS" or "TRADITIONAL" or "DYNAMICALL" or "SPANFACTOR" or "DYNAMICERI" or "CHOLESKY" or "SPANFACTOREUSE"
                      // For EPINTS section:     "ELEC_AUX" or "PROT_AUX" or "ELEC_AND_PROT_AUX" or "AUTO"
      "RITHRESHOLD",  // double
      "RISIGMA",      // double
      "RIMINSHRINK",  // size_t
      "RIMAXQUAL",    // size_t
      "RIGENCONTR",   // True or False
      "RIBUILD4INDEX",// True or False
      "RIREPORTERROR",// True of False, keyword only for EPINTS Section and only applies if user chooeses to approximate (ee|pp)
                      // Determine whether to explicitly build full precision (ee|pp) and calculate RMSD against approximate (ee|pp)
      "FINITENUCLEI", // True or False
      "LIBCINT",      // Ture or False

      "BARECOULOMB",  // True or False
      "LLLL",         // True or False
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
   *  \brief Construct a IntegralOptions object using the input file
   * 
   * 
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object
   *  
   *  \returns a IntegralOptions object
   *    constructed from the input options.
   */ 

  IntegralOptions getIntegralOptions(std::ostream &out, CQInputFile &input, 
      std::shared_ptr<BasisSet> basis,  std::shared_ptr<BasisSet> dfbasis, 
      std::shared_ptr<BasisSet> basis2, std::string int_sec){

    // check the validity of the integral option section 
    if (not int_sec.compare("INTS") and not int_sec.compare("PINTS") and not int_sec.compare("EPINTS"))
      CErr("Found invalue integral section");

    
    out << "  *** Parsing " << int_sec << ".REFERENCE options ***\n\n";

    IntegralOptions options;


    // Parse integral algorithm
    std::string ALG = "DIRECT";
    OPTOPT( ALG = input.getData<std::string>(int_sec+".ALG"); )
    trim(ALG);
    // Set integral algorithm
    if( not ALG.compare("DIRECT") )
      options.basicintsoptions.contrAlg = CONTRACTION_ALGORITHM::DIRECT;
    else if( not ALG.compare("INCORE") )
      options.basicintsoptions.contrAlg = CONTRACTION_ALGORITHM::INCORE;
    else
      CErr(ALG + " not a valid "+ int_sec + ".ALG",out);


    std::string TPITRANSALG = "N6"; 
    // Parse TPI AO to MO transformation algorithm
    OPTOPT( TPITRANSALG = input.getData<std::string>(int_sec+".TPITRANSALG"); )
    trim(TPITRANSALG);
    // Set TPI AO to MO transformation algorithm
    options.basicintsoptions.TPITRANSALG = TPITRANSALG;
    
    
    // Parse Schwarz threshold
    OPTOPT( options.basicintsoptions.threshSchwarz = input.getData<double>(int_sec+".SCHWARZ"); )


    // Set RI to be False by Default
    std::string RI = "FALSE";
    // For EPINTS, the default should be set to auto
    if(not int_sec.compare("EPINTS")) RI = "AUTO";
    // Parse RI option
    OPTOPT( RI = input.getData<std::string>(int_sec+".RI");)
    trim(RI);



    // Parse RI options
    if(RI.compare("FALSE")) {
      // If RI keyword anything other than "FALSE", then force contraction alg to be INCORE
      if (options.basicintsoptions.contrAlg != CONTRACTION_ALGORITHM::INCORE) {
        options.basicintsoptions.contrAlg = CONTRACTION_ALGORITHM::INCORE;
        std::cout << "Incore ERI algorithm enforced by RI." << std::endl;
      }

      // RI Keyword for EPINTS will be decoded differently than INTS/PINTS
      if(not int_sec.compare("EPINTS")){
        // Decode RI keywrod for EPINTS sections
        if (not RI.compare("AUTO")){
          options.cdriintsoptions.NEOCDalg = NEO_CD_ALG::AUTO;
        } else if (not RI.compare("ELEC_AUX")){
          options.cdriintsoptions.NEOCDalg = NEO_CD_ALG::ELEC_AUX;
        } else if (not RI.compare("PROT_AUX")){
          options.cdriintsoptions.NEOCDalg = NEO_CD_ALG::PROT_AUX;
        } else if (not RI.compare("ELEC_AND_PROT_AUX")){
          options.cdriintsoptions.NEOCDalg = NEO_CD_ALG::ELEC_AND_PROT_AUX;
        } else {
          CErr(RI + " is not a valid "+ int_sec + ".RI keyword",out);
        }

        std::string printError = "FALSE";
        OPTOPT( printError = input.getData<std::string>(int_sec+".RIPRINTERROR");)
        if(not printError.compare("FALSE")){
          options.cdriintsoptions.CDRI_printError = false;
        } else if(not printError.compare("TRUE")) {
          options.cdriintsoptions.CDRI_printError = true;
        } else {
          CErr(printError + " is not a valid " + int_sec + ".RIPRINTERROR keyword", out);
        }

      } else{
        // Decode RI keywrod for INTS / PINTS sections
          if(not RI.compare("AUXBASIS") ) {
            if (dfbasis->nBasis < 1)
              CErr("Keyword "+ int_sec + ".RI requires a non-empty DFbasis->",std::cout);
          } else if (not RI.compare("TRADITIONAL")) {
            options.cdriintsoptions.CDalg = CHOLESKY_ALG::TRADITIONAL;
          } else if (not RI.compare("DYNAMICALL")) {
            options.cdriintsoptions.CDalg = CHOLESKY_ALG::DYNAMIC_ALL;
          } else if (not RI.compare("SPANFACTOR")) {
            options.cdriintsoptions.CDalg = CHOLESKY_ALG::SPAN_FACTOR;
          } else if (not RI.compare("CHOLESKY") or not RI.compare("DYNAMICERI")) {
            options.cdriintsoptions.CDalg = CHOLESKY_ALG::DYNAMIC_ERI;
          } else if (not RI.compare("SPANFACTORREUSE")) {
            options.cdriintsoptions.CDalg = CHOLESKY_ALG::SPAN_FACTOR_REUSE;
          } else {
            CErr(RI + "is not a valid "+ int_sec + ".RI keyword",out);
          }
      }
    }

    // Set RI option
    options.basicintsoptions.RI = RI;

    // Parse more RI options
    OPTOPT( options.cdriintsoptions.CDRI_genContr = input.getData<bool>(int_sec+".RIGENCONTR"); )
    OPTOPT( options.cdriintsoptions.CDRI_thresh = input.getData<double>(int_sec+".RITHRESHOLD"); )
    OPTOPT( options.cdriintsoptions.CDRI_sigma = input.getData<double>(int_sec+".RISIGMA"); )
    OPTOPT( options.cdriintsoptions.CDRI_max_qual = input.getData<size_t>(int_sec+".RIMAXQUAL"); )
    OPTOPT( options.cdriintsoptions.CDRI_minShrinkCycle = input.getData<size_t>(int_sec+".RIMINSHRINK"); )
    OPTOPT( options.cdriintsoptions.CDRI_build4I = input.getData<bool>(int_sec+".RIBUILD4INDEX"); )

    return options;
  }

  /**
   *  \brief Set TPITransAlg for the IntegralBase Object.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns N/A
   *
   */
  void IntegralOptions::setTPITransAlg(std::shared_ptr<IntegralsBase> ints) const{
    if (not basicintsoptions.TPITRANSALG.compare("N5")) {
      if (basicintsoptions.contrAlg == CONTRACTION_ALGORITHM::DIRECT) {
        ints->TPITransAlg = TPI_TRANSFORMATION_ALG::DIRECT_N5;
        CErr("DIRECT_N5 TPI Transformation is NYI");
      } else
        ints->TPITransAlg = TPI_TRANSFORMATION_ALG::INCORE_N5;
    } else if (not basicintsoptions.TPITRANSALG.compare("N6")){
      if (basicintsoptions.contrAlg == CONTRACTION_ALGORITHM::DIRECT)
        ints->TPITransAlg = TPI_TRANSFORMATION_ALG::DIRECT_N6;
      else
        ints->TPITransAlg = TPI_TRANSFORMATION_ALG::INCORE_N6;
    } else {
      CErr(basicintsoptions.TPITRANSALG + " is not a valid TPITRANSALG",std::cout);
    }
  }



  /**
   *  \brief Construct a symmetric IntegralBase pointer from the IntegralOptions.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns A shared pointer to an IntegralBase object
   *
   */  
  std::shared_ptr<IntegralsBase> 
  IntegralOptions::buildSymmIntegral(std::ostream &out, CQMemManager &mem, Molecule &mol, 
      std::shared_ptr<BasisSet> basis, std::shared_ptr<BasisSet> dfbasis, std::string s) const{
    
    out << "Building integral object for the " << s << " subsystem:\n";

    std::shared_ptr<IntegralsBase> aoi = nullptr;

    if(basis->basisType == REAL_GTO) {
      std::shared_ptr<Integrals<double>> aoint =
          std::make_shared<Integrals<double>>();

      if(basicintsoptions.RI.compare("FALSE")) {
        if(not basicintsoptions.RI.compare("AUXBASIS"))
          aoint->TPI =
              std::make_shared<InCoreAuxBasisRIERI<double>>(mem,basis->nBasis,dfbasis);
        else
          aoint->TPI =
              std::make_shared<InCoreCholeskyRIERI<double>>(
                  mem, basis->nBasis, cdriintsoptions.CDRI_thresh, cdriintsoptions.CDalg, 
                  cdriintsoptions.CDRI_genContr, cdriintsoptions.CDRI_sigma, 
                  cdriintsoptions.CDRI_max_qual, cdriintsoptions.CDRI_minShrinkCycle, 
                  cdriintsoptions.CDRI_build4I);
      } else if (basicintsoptions.contrAlg == CONTRACTION_ALGORITHM::INCORE) {
        aoint->TPI =
            std::make_shared<InCore4indexTPI<double>>(mem,basis->nBasis);
      }
      else {
        aoint->TPI =
            std::make_shared<DirectTPI<double>>(mem,*basis,*basis,mol,basicintsoptions.threshSchwarz);
      }
      aoi = std::dynamic_pointer_cast<IntegralsBase>(aoint);
    } else if(basis->basisType == COMPLEX_GIAO) {
      std::shared_ptr<Integrals<dcomplex>> giaoint =
          std::make_shared<Integrals<dcomplex>>();
      if(basicintsoptions.RI.compare("FALSE"))
        CErr("GIAO resolution of identity ERI NYI",std::cout);
      else if (basicintsoptions.contrAlg == CONTRACTION_ALGORITHM::INCORE) {
        giaoint->TPI =
            std::make_shared<InCore4indexTPI<dcomplex>>(mem,basis->nBasis);
      }
      else {
        giaoint->TPI =
            std::make_shared<DirectTPI<dcomplex>>(mem,*basis,*basis,mol,basicintsoptions.threshSchwarz);
      }
      aoi = std::dynamic_pointer_cast<IntegralsBase>(giaoint);
    }
    
    setTPITransAlg(aoi);

    // Print
    out <<  *aoi << std::endl;
    
    return aoi;
  
  }


  /**
   *  \brief Construct an asymmetric IntegralBase pointer from the IntegralOptions.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns A shared pointer to an IntegralBase object
   *
   */  
  std::shared_ptr<IntegralsBase> 
  IntegralOptions::buildAsymmIntegral(std::ostream &out, CQMemManager &mem, Molecule &mol, std::shared_ptr<BasisSet> basis,  
      std::shared_ptr<BasisSet> dfbasis, std::shared_ptr<BasisSet> basis2, IntegralOptions eopts, IntegralOptions popts,
      std::shared_ptr<IntegralsBase> aoi, std::shared_ptr<IntegralsBase> paoi) const{
    
    
    std::shared_ptr<IntegralsBase> epaoi = nullptr;
    std::shared_ptr<Integrals<double>> epaoint = std::make_shared<Integrals<double>>();

    if (!basis2){
      CErr("Asymmetric integral object requires two different basis");
    } else {
      out << "Building integral object for the electron-quantum proton Coulomb term:\n\n";
      if(basis->basisType == REAL_GTO and basis2->basisType == REAL_GTO) {

        if(basicintsoptions.contrAlg == CONTRACTION_ALGORITHM::DIRECT) {
          epaoint->TPI = std::make_shared<DirectTPI<double>>(mem,*basis,*basis2,mol,basicintsoptions.threshSchwarz);
        } else {
          // By Default (if EPINTS.RI is not specified), Incore algorithm uses 4-index for (ee|pp)
          if(not basicintsoptions.RI.compare("FALSE")){
            epaoint->TPI = std::make_shared<InCore4indexTPI<double>>(mem,basis->nBasis,basis2->nBasis);
          } else{
            // If user specified algorithm an algorithm for NEO CD, then need tocreate corresponding IncoreAsymmRITPI object 
            
            // First detect if there are aux basis available:
            std::shared_ptr<InCoreRITPI<double>> aux1 = std::dynamic_pointer_cast<InCoreRITPI<double>> ((std::dynamic_pointer_cast<Integrals<double>>(aoi)->TPI));
            std::shared_ptr<InCoreRITPI<double>> aux2 = std::dynamic_pointer_cast<InCoreRITPI<double>> ((std::dynamic_pointer_cast<Integrals<double>>(paoi)->TPI));

            // If user choose alg to be auto, use flags to determine which to build:
            bool auto_elec = false, auto_prot = false, auto_elec_and_prot = false, auto_4I = false;

            // Execute user-specified algorithm to approximate (ee|pp)
            if(cdriintsoptions.NEOCDalg != NEO_CD_ALG::AUTO){
              // If user needs to use electronic aux basis, then will check and see if elec aux is available.
              // If not, we build aux from on the fly CD and save in IncoreAsymmRITPI class
              if(cdriintsoptions.NEOCDalg == NEO_CD_ALG::ELEC_AUX or cdriintsoptions.NEOCDalg == NEO_CD_ALG::ELEC_AND_PROT_AUX){
                out << bannerMid << std::endl;
                out << "   Will use (ee|ee) to approxiamate (ee|pp) " << std::endl;
                if(aux1){
                  out << "     * Found existing aux basis from (ee|ee)!" << std::endl;
                } else{
                  out << "     * Can't find existing aux basis from (ee|ee) " << std::endl;
                  out << "       Will use CD to build aux basis on the fly at " << eopts.cdriintsoptions.CDRI_thresh << " threshold" << std::endl; 
                  aux1 = 
                      std::make_shared<InCoreCholeskyRIERI<double>>(
                      mem, basis->nBasis, eopts.cdriintsoptions.CDRI_thresh, eopts.cdriintsoptions.CDalg, 
                      eopts.cdriintsoptions.CDRI_genContr, eopts.cdriintsoptions.CDRI_sigma, 
                      eopts.cdriintsoptions.CDRI_max_qual, eopts.cdriintsoptions.CDRI_minShrinkCycle, 
                      eopts.cdriintsoptions.CDRI_build4I);
                  // If TPI set to be incore 4-index, we can use that to do CD, which brings some saving
                  if (auto eri4I = std::dynamic_pointer_cast<InCore4indexTPI<double>> ((std::dynamic_pointer_cast<Integrals<double>>(aoi)->TPI)))
                    std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux1)->setFourIndexERI(eri4I);
                }
                out << bannerMid << "\n" << std::endl;
              }


              // If user choose to use protnonic aux basis , then will check and see if prot aux is available.
              // If not, we build aux from on the fly CD and save in IncoreAsymmRITPI class
              if(cdriintsoptions.NEOCDalg == NEO_CD_ALG::PROT_AUX or cdriintsoptions.NEOCDalg == NEO_CD_ALG::ELEC_AND_PROT_AUX){
                out << bannerMid << std::endl;
                std::cout << "   Will use (pp|pp) to approxiamate (ee|pp) " << std::endl;
                if(aux2){
                  out << "     * Found existing aux basis from (pp|pp)!" << std::endl;
                } else{
                  out << "     * Can't find existing aux basis from (pp|pp)" << std::endl;
                  out << "       Will use CD to build aux basis on the fly at " << popts.cdriintsoptions.CDRI_thresh << " threshold" << std::endl;
                  aux2 = 
                      std::make_shared<InCoreCholeskyRIERI<double>>(
                      mem, basis2->nBasis, popts.cdriintsoptions.CDRI_thresh, popts.cdriintsoptions.CDalg, 
                      popts.cdriintsoptions.CDRI_genContr, popts.cdriintsoptions.CDRI_sigma, 
                      popts.cdriintsoptions.CDRI_max_qual, popts.cdriintsoptions.CDRI_minShrinkCycle, 
                      popts.cdriintsoptions.CDRI_build4I);
                  // If TPI set to be incore 4-index, we can do dynamicERI CD algorithm on exisiting 4-index ERI, which can bring some savings
                  if (auto eri4I = std::dynamic_pointer_cast<InCore4indexTPI<double>> ((std::dynamic_pointer_cast<Integrals<double>>(paoi)->TPI)))
                    std::dynamic_pointer_cast<InCoreCholeskyRIERI<double>>(aux2)->setFourIndexERI(eri4I);
                }
                out << bannerMid << "\n" << std::endl;
              }


            // Automatically determine what aux basis is available. If none is availale, do 4 index  
            } else {
              if(aux1){
                out << "     * Detected existing aux basis from (ee|ee)!" << std::endl;
                auto_elec = true;
                if(aux2){
                  out << "     * Detected existing aux basis from (pp|pp)!" << std::endl;
                  auto_elec_and_prot = true;
                  auto_elec = false;
                } 
              } else {
                if(aux2){
                  out << "     * Detected existing aux basis from (pp|pp)!" << std::endl;
                  auto_prot = true;
                } else{
                  out << "     * Can't find existing aux basis to approximate (ee|pp)" << std::endl; 
                  out << "       Will use incore 4-index" << std::endl;
                  auto_4I = true;
                }
              } 
            }

            if(cdriintsoptions.NEOCDalg == NEO_CD_ALG::ELEC_AUX or auto_elec){
              epaoint->TPI = std::make_shared<InCoreAsymmRITPI<double>>(mem, aux1, basis2->nBasis, cdriintsoptions.CDRI_build4I);
              std::dynamic_pointer_cast<InCoreAsymmRITPI<double>>(epaoint->TPI)->setPrintError(cdriintsoptions.CDRI_printError);
              out << "Built (ee|pp) object that will use electronic aux basis. " << std::endl;
            } else if (cdriintsoptions.NEOCDalg == NEO_CD_ALG::PROT_AUX or auto_prot){
              epaoint->TPI = std::make_shared<InCoreAsymmRITPI<double>>(mem, basis->nBasis, aux2, cdriintsoptions.CDRI_build4I);
              std::dynamic_pointer_cast<InCoreAsymmRITPI<double>>(epaoint->TPI)->setPrintError(cdriintsoptions.CDRI_printError);
              out << "Built (ee|pp) object that will use protonic aux basis. " << std::endl;
            } else if (cdriintsoptions.NEOCDalg == NEO_CD_ALG::ELEC_AND_PROT_AUX or auto_elec_and_prot){
              epaoint->TPI = std::make_shared<InCoreAsymmRITPI<double>>(mem, aux1, aux2, cdriintsoptions.CDRI_build4I);
              std::dynamic_pointer_cast<InCoreAsymmRITPI<double>>(epaoint->TPI)->setPrintError(cdriintsoptions.CDRI_printError);
              out << "Built (ee|pp) object that will use both electronic and protonic aux basis. " << std::endl;
            } else if (auto_4I) {
              epaoint->TPI = std::make_shared<InCore4indexTPI<double>>(mem, basis->nBasis, basis2->nBasis);
              out << "Built 4-index (ee|pp) object. " << std::endl;
            } else {
              CErr ("aux basis for (ee|pp) is set up wrong. Can't build IncoreAsymmRITPI object!! ");
            }
          }  
        }
      } else{
        CErr("GIAO for (ee|pp) NYI",std::cout);
      }  
    }
    
    epaoi = std::dynamic_pointer_cast<IntegralsBase>(epaoint);
    setTPITransAlg(epaoi);
    
    // Print
    out <<  *epaoi << std::endl;
    
    return epaoi;
  }

  /**
   *  \brief Construct all IntegralBase pointers from the IntegralOptions.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] aoints AOIntegrals object for SingleSlater
   *                     construction
   *
   *  \returns All (ee|ee), (pp|pp), (ee|pp) IntegralBase objects (last two can be nullptrs if not needed), 
   * organized in a stuple of shared pointers. 
   *
   */  
    std::tuple<std::shared_ptr<IntegralsBase>, std::shared_ptr<IntegralsBase>, std::shared_ptr<IntegralsBase>> 
    IntegralOptions::buildAllIntegrals(
        std::ostream &out, CQMemManager &mem, Molecule &mol, std::shared_ptr<BasisSet> basis,  std::shared_ptr<BasisSet> dfbasis, 
        std::shared_ptr<BasisSet> basis2, IntegralOptions eopts, IntegralOptions popts, IntegralOptions epopts){

      out << BannerTop << std::endl;
      out << std::endl;

      // Build Electronic Integrals (ee|ee):
      std::shared_ptr<IntegralsBase> aoi = eopts.buildSymmIntegral(out, mem, mol, basis, dfbasis, "electronic");

      // Build Protonic Integrals (pp|pp) if we have protonic basis:
      std::shared_ptr<IntegralsBase> paoi = basis2 ? popts.buildSymmIntegral(out, mem, mol, basis2, dfbasis, "protonic") : nullptr;

      // Build Electron/Proton Coulumb Integrals (ee|pp) if we have protonic basis:
      std::shared_ptr<IntegralsBase> epaoi = basis2 ? epopts.buildAsymmIntegral(out, mem, mol, basis, dfbasis,basis2, eopts, popts, aoi, paoi) : nullptr;

      return  std::make_tuple(aoi, paoi, epaoi);
    }




  /**
   *  \brief Optionally set the control parameters for an
   *  AOIntegrals object
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] aoints AOIntegrals object 
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
      CErr(ALG + " not a valid "+ int_sec + ".ALG",out);

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
          CErr("Keyword "+ int_sec + ".RI requires a non-empty DFbasis->",std::cout);
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
        CErr(RI + " not a valid "+ int_sec + ".RI",out);
      }
    }

    OPTOPT( CDRI_genContr = input.getData<bool>(int_sec+".RIGENCONTR"); )
    OPTOPT( CDRI_thresh = input.getData<double>(int_sec+".RITHRESHOLD"); )
    OPTOPT( CDRI_sigma = input.getData<double>(int_sec+".RISIGMA"); )
    OPTOPT( CDRI_max_qual = input.getData<size_t>(int_sec+".RIMAXQUAL"); )
    OPTOPT( CDRI_minShrinkCycle = input.getData<size_t>(int_sec+".RIMINSHRINK"); )
    OPTOPT( CDRI_build4I = input.getData<bool>(int_sec+".RIBUILD4INDEX"); )

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
      CErr(TPITRANSALG + " not a valid "+ int_sec + ".TPITRANSALG",out);
    }

    // Print
    out <<  *aoi << std::endl;

    
    return aoi;

  }; // CQIntsOptions

}; // namespace ChronusQ
