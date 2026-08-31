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

namespace ChronusQ {

  void printInvalidKeys(const std::set<std::string> &invalidKeywords,
                        const std::string &prefix) {
    if (invalidKeywords.empty()) return;
    std::cout << "Keywords in section " << prefix << " are not recognized: ";
    for (const auto &key : invalidKeywords) {
      std::cout << key << " ";
    }
    CErr("Keywords in section " + prefix + " are not recognized.",std::cout);
  }

  std::set<std::string> CQInvalidKeywords(
      const std::set<std::string> &allowedKeywords,
      const std::map<std::string, std::string>& inputSection) {
    std::set<std::string> invalidKeywords;

    // Make sure all of basisKeywords in allowedKeywords
    for( const auto &[keyword, value] : inputSection ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() )
        invalidKeywords.insert(keyword);
    }

    return invalidKeywords;
  }

  void CQINPUT_VALID(std::ostream &out, CQInputFile &input) {

    std::set<std::string> invalidKeywords;

    if (input.containsSection("MOLECULE")) {
      invalidKeywords = CQMOLECULE_VALID(input.getSection("MOLECULE"));
      printInvalidKeys(invalidKeywords, "MOLECULE");
    }
    if (input.containsSection("BASIS")) {
      invalidKeywords = CQBASIS_VALID(input.getSection("BASIS"));
      printInvalidKeys(invalidKeywords, "BASIS");
    }
    if (input.containsSection("DFBASIS")) {
      invalidKeywords = CQBASIS_VALID(input.getSection("DFBASIS"));
      printInvalidKeys(invalidKeywords, "DFBASIS");
    }
    if (input.containsSection("INTS")) {
      invalidKeywords = CQINTS_VALID(input.getSection("INTS"));
      printInvalidKeys(invalidKeywords, "INTS");
    }
    if (input.containsSection("QM")) {
      invalidKeywords = CQQM_VALID(input.getSection("QM"));
      printInvalidKeys(invalidKeywords, "QM");
    }
    if (input.containsSection("DFTINT")) {
      invalidKeywords = CQDFTINT_VALID(input.getSection("DFTINT"));
      printInvalidKeys(invalidKeywords, "DFTINT");
    }
    if (input.containsSection("SCF")) {
      invalidKeywords = CQSCF_VALID(input.getSection("SCF"));
      printInvalidKeys(invalidKeywords, "SCF");
    }
    if (input.containsSection("RT")) {
      invalidKeywords = CQRT_VALID(input.getSection("RT"));
      printInvalidKeys(invalidKeywords, "RT");
    }
    if (input.containsSection("RESPONSE")) {
      invalidKeywords = CQRESPONSE_VALID(input.getSection("RESPONSE"));
      printInvalidKeys(invalidKeywords, "RESPONSE");
    }
    if (input.containsSection("MOR")) {
      invalidKeywords = CQMOR_VALID(input.getSection("MOR"));
      printInvalidKeys(invalidKeywords, "MOR");
    }
    if (input.containsSection("CUBE")) {
      invalidKeywords = CQCUBE_VALID(input.getSection("CUBE"));
      printInvalidKeys(invalidKeywords, "CUBE");
    }
    if (input.containsSection("ORBPROP")) {
      invalidKeywords = CQORBPROP_VALID(input.getSection("ORBPROP"));
      printInvalidKeys(invalidKeywords, "ORBPROP");
    }
    if (input.containsSection("MISC")) {
      invalidKeywords = CQMISC_VALID(input.getSection("MISC"));
      printInvalidKeys(invalidKeywords, "MISC");
    }
    if (input.containsSection("CC")) {
      invalidKeywords = CQCC_VALID(input.getSection("CC"));
      printInvalidKeys(invalidKeywords, "CC");
    }
    if (input.containsSection("DYNAMICS")) {
      invalidKeywords = CQDYNAMICS_VALID(input.getSection("DYNAMICS"));
      printInvalidKeys(invalidKeywords, "DYNAMICS");
    }
    if (input.containsSection("MCSCF")) {
      invalidKeywords = CQMCSCF_VALID(input.getSection("MCSCF"));
      printInvalidKeys(invalidKeywords, "MCSCF");
    }
    if (input.containsSection("EOMCC")) {
      invalidKeywords = CQEOMCC_VALID(input.getSection("EOMCC"));
      printInvalidKeys(invalidKeywords, "EOMCC");
    }
    if (input.containsSection("PERTURB")) {
      invalidKeywords = CQPERTURB_VALID(input.getSection("PERTURB"));
      printInvalidKeys(invalidKeywords, "PERTURB");

      // No GIAO + PERTURB
      if( input.containsData("BASIS/BASISTYPE") ) {

        std::string btype =
            input.getData<std::string>("BASIS/BASISTYPE");

        if( not btype.compare("GIAO") )
          CErr("GIAO + PERTURB not allowed");

      }

      // No TDDFT
      if( input.containsData("QM/REFERENCE") ) {

        std::string ref = input.getData<std::string>("QM/REFERENCE");

        bool isKS  = not (ref.find("HF") != std::string::npos);

        if( isKS )
          CErr("PERTURB + KS not allowed");

      }
    }
    if (input.containsSection("MP2")) {
      invalidKeywords = CQMP2_VALID(input.getSection("MP2"));
      printInvalidKeys(invalidKeywords, "MP2");

      if( input.containsData("QM/REFERENCE") ) {

        std::string ref = input.getData<std::string>("QM/REFERENCE");

        bool isKS  = not (ref.find("HF") != std::string::npos);

        if( isKS )
          CErr("MP2 + KS not allowed");
      }
    }
    if (input.containsSection("CI")) {
      invalidKeywords = CQCI_VALID(input.getSection("CI"));
      printInvalidKeys(invalidKeywords, "CI");
    }
    if (input.containsSection("MRPT")) {
      invalidKeywords = CQMRPT_VALID(input.getSection("MRPT"));
      printInvalidKeys(invalidKeywords, "MRPT");

      // No GIAO + MRPT
      if( input.containsData("BASIS/BASISTYPE") ) {

        std::string btype =
            input.getData<std::string>("BASIS/BASISTYPE");

        if( not btype.compare("GIAO") )
          CErr("GIAO + MRPT not allowed");

      }
    }
    if (input.containsSection("GAUXC")) {
      invalidKeywords = CQGAUXC_VALID(input.getSection("GAUXC"));
      printInvalidKeys(invalidKeywords, "GAUXC");
    }
    if (input.containsSection("PHYSCON")) {
      invalidKeywords = CQPHYSCON_VALID(input.getSection("PHYSCON"));
      printInvalidKeys(invalidKeywords, "PHYSCON");
    }

  }

}; // namespace ChronusQ