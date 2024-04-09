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

namespace ChronusQ {

  /**
   * 
   *  Check valid keywords in the section.
   *
   */
  void CQPERTURB_VALID( std::ostream &out, CQInputFile &input ) {

    if( not input.containsSection("PERTURB") ) return;

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "DOFULL",
      "DOITER",
      "ITERCONV",
      "MAXITER",
      "STATEAVERAGE",
      "EXTENDMS",
      "DOGVVPT",
      "LEVELSHIFT",
      "IMAGINARYSHIFT",
      "FROZENCORE",
      "FROZENVIRTUAL",
      "SELECTVIRTUAL",
      "SOI"
    };

    // Specified keywords
    std::vector<std::string> perturbKeywords = input.getDataInSection("PERTURB");

    // Make sure all of the basis keywords in allowed keywords
    for( auto &keyword : perturbKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() )
      CErr("Keyword PERTURB." + keyword + " is not recognized",std::cout);// Error
    }

    // Check for disallowed combinations (if any)
    // No GIAO + PERTURB
    if( input.containsData("BASIS.BASISTYPE") ) {

      std::string btype =
        input.getData<std::string>("BASIS.BASISTYPE");

      if( not btype.compare("GIAO") )
        CErr("GIAO + PERTURB not allowed");

    }

    // No TDDFT
    if( input.containsData("QM.REFERENCE") ) {

      std::string ref = input.getData<std::string>("QM.REFERENCE");

      bool isKS  = not (ref.find("HF") != std::string::npos);

      if( isKS )
        CErr("PERTURB + KS not allowed");

    }
  } // CQPERTURB_VALID

  /**
   * 
   *  \brief Construct a PERTURB object using the input
   *
   *  \param [in] out Output device for data/error output.
   *  \param [in] input Input file datastructure
   *  \param [in] ref reference: SingleSlater or MCSCF
   *
   *  \returns shared_ptr to a PerturbBase object
   *
   */
  std::shared_ptr<MCWaveFunctionBase> CQPerturbOptions(std::ostream &out,
    CQInputFile &input, std::shared_ptr<MCWaveFunctionBase> &ref) {

    if( not input.containsSection("PERTURB") )
      CErr("PERTURB Section must be specified for PERTURB job",out);

    out << "  *** Parsing PERTURB options ***\n";

    // Parse state of interest
    std::string sSoI;
    std::vector<size_t> SoI;
    OPTOPT( sSoI = input.getData<std::string>("PERTURB.SOI"); )
    if (sSoI.empty()) {
      std::cout << "Getting MCSCF default state of interests." << std::endl;
      for (auto i = 0ul; i < ref->NStates; i++)
        SoI.push_back(i);
    }
    else {
      std::vector<std::string> sSoITokens;
      split(sSoITokens, sSoI, " ,;");
      for (auto & soi_i: sSoITokens)
        SoI.push_back(std::stoul(soi_i) - 1);
    }
    std::cout <<"number of states: "<<SoI.size()<<std::endl;

//    out << "  *** Parsing PERTURB options ***\n";
    std::shared_ptr<MCWaveFunctionBase> perturb = nullptr;
    PerturbOptions *PTopts;
//    MOSpacePartition refMOpart;
    // Determine reference and construct PERTURB object

    bool found = false;

    #define CONSTRUCT_PERTURB(_MT,_IT)             \
    if( not found ) try {                          \
      auto &tmp = dynamic_cast<MCSCF<_MT,_IT>&>(*ref);\
      perturb = std::dynamic_pointer_cast<MCWaveFunctionBase>( \
            std::make_shared<PERTURB<_MT,_IT>>( \
            std::dynamic_pointer_cast<MCSCF<_MT,_IT>>(ref),SoI)); \
      PTopts = &(std::dynamic_pointer_cast<PERTURB<_MT,_IT>>(perturb)->PTopts); \
      std::cout<<"Perturb constructed."<<std::endl; \
      found = true;\
    } catch(...) { } 

    // Construct PERTURB object

    CONSTRUCT_PERTURB( double,   double   );
    CONSTRUCT_PERTURB( dcomplex, double   );
    CONSTRUCT_PERTURB( dcomplex, dcomplex );

    // Parse options


    OPTOPT( PTopts->saFock = input.getData<bool>("PERTURB.STATEAVERAGE"); )
    if (SoI.size() == 1) PTopts->saFock = false;

    OPTOPT( PTopts->extendMS = input.getData<bool>("PERTURB.EXTENDMS"); )
    if( PTopts->extendMS ) PTopts->saFock = true;
    std::cout<<"extended multi-state: "<<PTopts->extendMS<<std::endl;

    OPTOPT( PTopts->doFull = input.getData<bool>("PERTURB.DOFULL"); )
    std::cout<<"DOFULL? "<<PTopts->doFull<<std::endl;

    OPTOPT( PTopts->doIter = input.getData<bool>("PERTURB.DOITER"); )
    std::cout<<"DOITER? "<<PTopts->doIter<<std::endl;

    // TODO
//    OPTOPT( PTopts->doGVV = input.getData<bool>("PERTURB.DOGVVPT"); )
//    if (PTopts->doGVV) PTopts->saFock = false;
//    std::cout<<"GVVPT type: "<<PTopts->doGVV<<std::endl;
//    std::cout<<"stateaverage (saFock): "<<PTopts->saFock<<std::endl;
      
    OPTOPT( PTopts->maxIter = input.getData<size_t>("PERTURB.MAXITER"); )
    std::cout<<"Max num of GMRES iterations: "<<PTopts->maxIter<<std::endl;

    OPTOPT( PTopts->convCrit = input.getData<double>("PERTURB.ITERCONV"); )
    std::cout<<"convergence criteria: "<<PTopts->convCrit<<std::endl;

    OPTOPT( PTopts->levelShift = input.getData<double>("PERTURB.LEVELSHIFT"); )
    std::cout<<"level shift: "<<PTopts->levelShift<<std::endl;

    OPTOPT( PTopts->imaginaryShift = input.getData<double>("PERTURB.IMAGINARYSHIFT"); )
    std::cout<<"imaginary shift: "<<PTopts->imaginaryShift<<std::endl;

    OPTOPT( PTopts->frozenCore = input.getData<size_t>("PERTURB.FROZENCORE"); )
    std::cout<<"frozen core: "<<PTopts->frozenCore<<std::endl;

    OPTOPT( PTopts->frozenVirtual = input.getData<size_t>("PERTURB.FROZENVIRTUAL"); )
    std::cout<<"frozen virtual: "<<PTopts->frozenVirtual<<std::endl;

    OPTOPT( PTopts->selectVirtual = input.getData<std::string>("PERTURB.SELECTVIRTUAL"); )
    std::cout<<"select virtual: "<<PTopts->selectVirtual<<std::endl;

    // reset frozenvirtual based on selectvirtual input
    size_t svirtual = 0;
    size_t fvirtual = PTopts->frozenVirtual;
    if (not PTopts->selectVirtual.empty()) {

      std::vector<std::string> moTokens;
      split(moTokens, PTopts->selectVirtual, ", ");
      for (auto & mo: moTokens) {
        std::vector<std::string> mo2;
        split(mo2, mo, "-");
        if (mo2.size() == 1) {
          svirtual += 1;
        } else if (mo2.size() == 2) {
          for (auto i = std::stoul(mo2[0]); i <= std::stoul(mo2[1]); i++)
            svirtual += 1;
        } else CErr("Unrecogonized pattern in orbital selection");
      }
      fvirtual = ref->MOPartition.nElecMO - ref->MOPartition.nInact - ref->MOPartition.nCorrO - svirtual;
    }

    if ( fvirtual != PTopts->frozenVirtual ) {
      std::cout << "Update FrozenVirtual based on input selectVirtual."<< std::endl;
      PTopts->frozenVirtual = fvirtual;
      std::cout <<"frozen virtual: "<<PTopts->frozenVirtual<<std::endl;
    }


    return perturb;
  }; // CQPerturbOptions

}; // namespace ChronusQ













