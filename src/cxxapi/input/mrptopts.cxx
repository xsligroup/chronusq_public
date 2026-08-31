/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
  std::set<std::string> CQMRPT_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
      "STATEAVERAGE",
      "GVVPT",
      "ENPT",
      "FROZENCORE",
      "FROZENVIRTUAL",
      "SELECTVIRTUAL",
      "LEVELSHIFT",
      "SPINFREE",
      "EPS",
      "NROOTS",
      "SECONDARYROOTS",
      "SOI"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  } // CQMRPT_VALID

  /**
   * 
   *  \brief Construct a MRPT object using the input
   *
   *  \param [in] out Output device for data/error output.
   *  \param [in] input Input file datastructure
   *  \param [in] ref reference: SingleSlater or MCSCF
   *
   *  \returns shared_ptr to a Posthartreefock base object
   *
   */
  std::shared_ptr<PostHartreeFockBase> CQMRPTSettings(std::ostream &out,
    CQInputFile &input, std::shared_ptr<PostHartreeFockBase> &ref) {

    if( not input.containsSection("MRPT") )
      CErr("MRPT Section must be specified for PT job",out);
    out << "\n*** Parsing MRPT Settings ***\n";

    // Parse Target States
    size_t nTargetStates = static_cast<size_t>(-1);
    std::string sTargetStates;
    std::vector<size_t> TargetStates;
    OPTOPT(sTargetStates = input.getData<std::string>("MRPT/SOI");)
    OPTOPT(nTargetStates = input.getData<size_t>("MRPT/NROOTS");)
    if (sTargetStates.empty() && nTargetStates == static_cast<size_t>(-1)) {
      std::cout << "Retrieve CAS Default Number of Roots!" << std::endl;
      for (auto i = 0ul; i < ref->NStates; i++)
        TargetStates.push_back(i);
    } else if (sTargetStates.empty() && nTargetStates != static_cast<size_t>(-1)) {
        for (auto i = 0ul; i < nTargetStates; i++)
        TargetStates.push_back(i);
    }
    else { 
      std::vector<std::string> sTargetStatesTok;
      split(sTargetStatesTok, sTargetStates, " ,;");
      for (auto & i: sTargetStatesTok)
        TargetStates.push_back(std::stoul(i)-1);
    }
    std::shared_ptr<PostHartreeFockBase> mrpt = nullptr;
    MRPTSettings *PTopts;

    bool found = false;
    #define CONSTRUCT_MRPT(_MT,_IT) \
    if( not found ) try {\
      auto &tmp = dynamic_cast<ConfigurationInteraction<_MT,_IT>&>(*ref);\
            mrpt = std::dynamic_pointer_cast<PostHartreeFockBase>(\
            std::make_shared<DasPerturb<_MT,_IT>>(\
            std::dynamic_pointer_cast<ConfigurationInteraction<_MT,_IT>>(ref),TargetStates));\
      PTopts = &(std::dynamic_pointer_cast<DasPerturb<_MT,_IT>>(mrpt)->PTopts);\
      found = true;\
    } catch(...) { } 

    // Construct MRPT object
    CONSTRUCT_MRPT( double,   double   );
    CONSTRUCT_MRPT( dcomplex, double   );
    CONSTRUCT_MRPT( dcomplex, dcomplex );

    // Parse options
    OPTOPT( PTopts->SPINFREE = input.getData<bool>("MRPT/SPINFREE");)
    OPTOPT( PTopts->EPS = input.getData<double>("MRPT/EPS");)
    OPTOPT( PTopts->ENPT = input.getData<bool>("MRPT/ENPT"); )
    OPTOPT( PTopts->GVVPT = input.getData<bool>("MRPT/GVVPT"); )

    if (PTopts->ENPT) {
      PTopts->STATEAVERAGE = false;
      PTopts->GVVPT = false;
      OPTOPT( PTopts->FROZENCORE = input.getData<size_t>("MRPT/FROZENCORE"); )
      OPTOPT( PTopts->FROZENVIRTUAL = input.getData<size_t>("MRPT/FROZENVIRTUAL"); )
      OPTOPT( PTopts->SELECTVIRTUAL = input.getData<std::string>("MRPT/SELECTVIRTUAL"); )
      OPTOPT( PTopts->LEVELSHIFT = input.getData<double>("MRPT/LEVELSHIFT"); )
      size_t svirtual = 0;
      size_t fvirtual = PTopts->FROZENVIRTUAL;
      if (not PTopts->SELECTVIRTUAL.empty()) {
        std::vector<std::string> moTokens;
        split(moTokens, PTopts->SELECTVIRTUAL, ", ");
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
        fvirtual = ref->corrSpace.nElecMO - ref->corrSpace.nInact - 
          ref->corrSpace.nCorrO - svirtual;
      }
    }
    else if (PTopts->GVVPT) {
      PTopts->ENPT = false;
      OPTOPT( PTopts->STATEAVERAGE = input.getData<bool>("MRPT/STATEAVERAGE"); )
      OPTOPT( PTopts->FROZENCORE = input.getData<size_t>("MRPT/FROZENCORE"); )
      OPTOPT( PTopts->FROZENVIRTUAL = input.getData<size_t>("MRPT/FROZENVIRTUAL"); )
      OPTOPT( PTopts->SECONDARYROOTS = input.getData<size_t>("MRPT/SECONDARYROOTS"); )
    }
    else { CErr("Exception: No MRPT2 flavor specified / Incompartible MRPT2 method"); }

    return mrpt;
  }; // CQMRPTSettings

}; // namespace ChronusQ
