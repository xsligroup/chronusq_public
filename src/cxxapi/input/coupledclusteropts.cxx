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
  /**
   *
   *  Check valid keywords in the CC section.
   *
  */
  std::set<std::string> CQCC_VALID(const std::map<std::string, std::string>& inputSection) {
   //Allowed keywords
    std::set<std::string> allowedKeywords = {
      "TYPE",
      "USEDIIS",
      "NDIIS",
      "ETOL",
      "TTOL",
      "MAXITER",
      "TABLKSIZE",
      "NEVARIATION",
      "DENOMSHIFT",
      "PERTT3",
      "CRCC",
      "TRIPLESMPI",
      "TRIPLESBEGIN",
      "TRIPLESEND",
      "LOOPABC",
      "CCSDINIT",
      "RESTART",
      "SAVE",
      "FROZENOCCUPIED",
      "FROZENVIRTUAL",
      "REBUILDFOCK",
      "SKIPSCF",
      "SKIPCC",
      "DIPOLE"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

  /**
   *
   *  Check valid keywords in the EOMCC section.
   *
  */
  std::set<std::string> CQEOMCC_VALID(const std::map<std::string, std::string>& inputSection) {
    //Allowed keywords
    std::set<std::string> allowedKeywords = {
        "EOMTYPE",
        "NROOTS",
        "HBARTYPE",
        "DIAGMETHOD",
        "CVSCORE",
        "ACTIVEOCCUPIED",
        "ACTIVEVIRTUAL",
        //"CVSCONTINUUM",
        //"FROZENOCCUPIED",
        //"FROZENVIRTUAL",
        "IPLEVEL",
        "DIPMAXEXTERNALHOLE",
        "DIPMAXEXTERNALPARTICLE",
        "DAVIDSONWHENSC",
        "DAVIDSONMAXMACROITER",
        "DAVIDSONMAXMICROITER",
        "DAVIDSONRCHECKESIDUAL",
        "DAVIDSONCHECKEVEC",
        "DAVIDSONCHECKEVAL",
        "DAVIDSONRESIDUALCONV",
        "DAVIDSONEVECCONV",
        "DAVIDSONEVALCONV",
        "DAVIDSONCONVONGRAMSCHMIDT",
        "DAVIDSONSUBSPACEMULTIPLIER",
        "DAVIDSONGUESSMULTIPLIER",
        "DAVIDSONENERGYSPECIFICABS",
        "DAVIDSONPRECONDSMALL",
//        "DAVIDSONSORTBYDISTANCE",
        "DAVIDSONBIORTHO",
        "GRAMSCHMIDTREPEAT",
        "GRAMSCHMIDTEPS",
        "OSCILLATORSTRENGTH",
        "ALLEXCITEDSTATEDIPOLE",
        "SAVEHAMILTONIAN",
        "PRINTLARGEAMPLITUDE",
        "CCSGUESS",
        "RESTARTL",
        "RESTARTR",
        "SKIPR",
        "SAVEL",
        "SAVER",
        "SINGLETONLY",
        "ENERGYIMAGTOL"
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  }

#ifdef CQ_HAS_TA

  std::vector<size_t> parseOrbitalSelectionInput(std::string input_str);

  //Construct a CoupledClusterSettings object from input file
  CoupledClusterSettings CQCCOptions(std::ostream &out, CQInputFile &input) {

    if( not input.containsSection("CC") )
      CErr("CC Section must be specified for CC job",out);

    out << "\n  *** Parsing CC options ***\n";


    // Determine  reference and construct CC object

    int cc_level = 0;

    CoupledClusterSettings ccSettings;
    OPTOPT(
      std::string ccopts = input.getData<std::string>("CC/TYPE");

      if( not ccopts.compare("CCSD") ) {
        cc_level = 2;
        ccSettings.cctype = CC_TYPE::CCSD;
      }else if( not ccopts.compare("DFCCSD") ) {
        cc_level = 2;
        ccSettings.cctype = CC_TYPE::DFCCSD;
      }else if( not ccopts.compare("CCSDT") ) {
        cc_level = 3; 
        ccSettings.cctype = CC_TYPE::CCSDT;
      }else 
        CErr(ccopts + " NOT RECOGNIZED CC/TYPE");
    );


    if (cc_level){
      if(input.containsData("CC/SKIPCC")){
        OPTOPT(ccSettings.skipCC = input.getData<bool>("CC/SKIPCC");)
        if (ccSettings.skipCC) { // must read T when skipping CC
            ccSettings.restart = true;
            ccSettings.skipSCF = true;
        }
      }

      if(input.containsData("CC/SKIPSCF")){
        OPTOPT(ccSettings.skipSCF = input.getData<bool>("CC/SKIPSCF");)
      }

      if(input.containsData("CC/USEDIIS")){
        OPTOPT(ccSettings.useDIIS = input.getData<bool>("CC/USEDIIS");)
      }

      if(input.containsData("CC/NDIIS")){
        if(not ccSettings.useDIIS){
          std::cout << "You must enable DIIS to specify nDIIS. " << std::endl;
        }
        else{
          OPTOPT(ccSettings.nDIIS = input.getData<size_t>("CC/NDIIS");)
        }  
      }

      if(input.containsData("CC/ETOL")){
        OPTOPT(ccSettings.eConv = input.getData<double>("CC/ETOL");)
      }

      if(input.containsData("CC/TTOL")){
        OPTOPT(ccSettings.tConv = input.getData<double>("CC/TTOL");)
      }

      if(input.containsData("CC/MAXITER")){
        OPTOPT(ccSettings.maxiter = input.getData<int>("CC/MAXITER");)
      }

      if(input.containsData("CC/TABLKSIZE")){
        OPTOPT(ccSettings.blksize = input.getData<int>("CC/TABLKSIZE");)
      }

      if(input.containsData("CC/NEVARIATION")){
        OPTOPT(ccSettings.nEvariation = input.getData<int>("CC/NEVARIATION");)
      }

      if(input.containsData("CC/DENOMSHIFT")){
        OPTOPT(ccSettings.denomshift = input.getData<double>("CC/DENOMSHIFT");)
      }

      if(input.containsData("CC/PERTT3")){
        OPTOPT(ccSettings.pertT3 = input.getData<bool>("CC/PERTT3");)
      }

      if(input.containsData("CC/CRCC")){
        OPTOPT(ccSettings.crcc = input.getData<bool>("CC/CRCC");)
      }

      if(input.containsData("CC/TRIPLESMPI")){
        OPTOPT(ccSettings.triplesMPI = input.getData<bool>("CC/TRIPLESMPI");)
      }

      if(input.containsData("CC/TRIPLESBEGIN")){
        OPTOPT(ccSettings.triples_begin = input.getData<int>("CC/TRIPLESBEGIN");)
      }

      if(input.containsData("CC/TRIPLESEND")){
        OPTOPT(ccSettings.triples_end = input.getData<int>("CC/TRIPLESEND");)
      }

      if(input.containsData("CC/LOOPABC")){
        OPTOPT(ccSettings.loop_abc = input.getData<bool>("CC/LOOPABC");)
      }

      if(input.containsData("CC/CCSDINIT")){
        OPTOPT(ccSettings.CCSDinit = input.getData<bool>("CC/CCSDINIT");)
      }

      if(input.containsData("CC/RESTART")){
        OPTOPT(ccSettings.restart = input.getData<bool>("CC/RESTART");)
      }

      if(input.containsData("CC/SAVE")){
        OPTOPT(ccSettings.save = input.getData<bool>("CC/SAVE");)
      }

      if(input.containsData("CC/FROZENOCCUPIED")){
        std::string frozen_occ_str;
        OPTOPT(frozen_occ_str = input.getData<std::string>("CC/FROZENOCCUPIED"););
        ccSettings.frozen_occupied = parseOrbitalSelectionInput(frozen_occ_str);
        // Frozen occupied requires rebuild Fock!
        // DO NOT USE REBUILD FOCK FOR MMFX2C! YOU WILL HAVE THE WRONG FOCK MATRIX
        // ccSettings.rebuildFock = true;
      }

      if(input.containsData("CC/FROZENVIRTUAL")){
        std::string frozen_vir_str;
        OPTOPT(frozen_vir_str = input.getData<std::string>("CC/FROZENVIRTUAL"););
        ccSettings.frozen_virtual = parseOrbitalSelectionInput(frozen_vir_str);
      }

      if(input.containsData("CC/REBUILDFOCK")){
        /*
         * *** REBUILDFOCK SHOULD NOT BE USED WITH mmfX2C ***
         * because it will build 2c Fock matrix from 2c coreH and 2c TPI using the 2c mmfX2C MO coefficients
         * and thus will not recover the correct reference energy computed using mmfX2C fock matrix
         */
        OPTOPT(ccSettings.rebuildFock = input.getData<bool>("CC/REBUILDFOCK");)
        if (not ccSettings.rebuildFock and ccSettings.nEvariation != 0) {
          CErr("CC/NEVARIATION being non-zero requires CC/REBUILDFOCK = True");
        }
        if (ccSettings.rebuildFock) {
          std::string X = "DEFAULT";
          OPTOPT( X = input.getData<std::string>("QM/X2CTYPE")  );
          trim(X);
          if (not X.compare("FOCK")) {
            CErr("CC/REBUILDFOCK should not be used with MMF-X2C");
          }
        }
      } else if (ccSettings.nEvariation != 0) {
        ccSettings.rebuildFock = true;
        std::cout << "      ccSettings.rebuildFock default to True for inequal number "
                  << "of electrons between reference and CCSD calculation." << std::endl;
      }

      if(input.containsData("CC/DIPOLE")){
        OPTOPT(ccSettings.computeDipole = input.getData<bool>("CC/DIPOLE");)
      }

    }
    else {
      CErr("NYI");
    }

    return ccSettings;
  }



  //Construct a EOMCCSettings object from input file
  EOMSettings CQEOMCCOptions(std::ostream &out, CQInputFile &input) {

    if( not input.containsSection("EOMCC") )
      CErr("EOMCC Section must be specified for EOMCC job",out);

    out << "\n  *** Parsing EOMCC options ***\n";

    EOMSettings eomSettings;
    if(input.containsData("EOMCC/EOMTYPE")){
      std::string eom_type_str = input.getData<std::string>("EOMCC/EOMTYPE");
      if( not eom_type_str.compare("EE") )
        eomSettings.eom_type = EOM_TYPE::EE;
      else if( not eom_type_str.compare("IP") ){
        eomSettings.eom_type = EOM_TYPE::IP;
        eomSettings.ip_level = 2;
      }
      else if( not eom_type_str.compare("EA") )
        eomSettings.eom_type = EOM_TYPE::EA;
      else if( not eom_type_str.compare("DIP") ){
        eomSettings.eom_type = EOM_TYPE::DIP;
        eomSettings.ip_level = 3;
      }
      else
        CErr(eom_type_str + " NOT RECOGNIZED EOMCC/EOMTYPE");
    }

    if(input.containsData("EOMCC/HBARTYPE")){
      std::string hbar_type_str = input.getData<std::string>("EOMCC/HBARTYPE");
      if( not hbar_type_str.compare("EXPLICIT") )
        eomSettings.hbar_type = EOM_HBAR_TYPE::EXPLICIT;
      else if( not hbar_type_str.compare("IMPLICIT") )
        eomSettings.hbar_type = EOM_HBAR_TYPE::IMPLICIT;
      else if( not hbar_type_str.compare("DEBUG") )
        eomSettings.hbar_type = EOM_HBAR_TYPE::DEBUG;
      else
        CErr(hbar_type_str + " NOT RECOGNIZED EOMCC/HBARTYPE");
    }

    if(input.containsData("EOMCC/DIAGMETHOD")){
      std::string diag_method_str = input.getData<std::string>("EOMCC/DIAGMETHOD");
      if( not diag_method_str.compare("FULL") )
        eomSettings.diag_method = EOM_DIAG_METHOD::FULL;
      else if( not diag_method_str.compare("DAVIDSON") )
        eomSettings.diag_method = EOM_DIAG_METHOD::DAVIDSON;
      else if( not diag_method_str.compare("GPLHR") ) {
        eomSettings.diag_method = EOM_DIAG_METHOD::GPLHR;
        CErr("GPLHR-EOMCC NYI");
      } else
        CErr(diag_method_str + " NOT RECOGNIZED EOMCC/DIAGMETHOD");
    }

    if(input.containsData("EOMCC/CVSCORE")){
      std::string cvs_core_str;
      OPTOPT(cvs_core_str = input.getData<std::string>("EOMCC/CVSCORE"););
      eomSettings.cvs_core = parseOrbitalSelectionInput(cvs_core_str);
    }

    if(input.containsData("EOMCC/ACTIVEOCCUPIED")){
      std::string active_occupied_str;
      OPTOPT(active_occupied_str = input.getData<std::string>("EOMCC/ACTIVEOCCUPIED"););
      eomSettings.active_occupied = parseOrbitalSelectionInput(active_occupied_str);
    }

    if(input.containsData("EOMCC/ACTIVEVIRTUAL")){
      std::string active_virtual_str;
      OPTOPT(active_virtual_str = input.getData<std::string>("EOMCC/ACTIVEVIRTUAL"););
      eomSettings.active_virtual = parseOrbitalSelectionInput(active_virtual_str);
    }

    //if(input.containsData("EOMCC/CVSCONTINUUM")){
    //  std::string cvs_vir_str;
    //  OPTOPT(cvs_vir_str = input.getData<std::string>("EOMCC/CVSCONTINUUM"););
    //  eomSettings.cvs_virtual = parseOrbitalSelectionInput(cvs_vir_str);
    //}

    if(input.containsData("EOMCC/IPLEVEL")){
      OPTOPT(eomSettings.ip_level = input.getData<int>("EOMCC/IPLEVEL");)
    }

    if(input.containsData("EOMCC/DIPMAXEXTERNALHOLE")){
      OPTOPT(eomSettings.dip_max_external_hole = input.getData<int>("EOMCC/DIPMAXEXTERNALHOLE");)
    }

    if(input.containsData("EOMCC/DIPMAXEXTERNALPARTICLE")){
      OPTOPT(eomSettings.dip_max_external_particle = input.getData<int>("EOMCC/DIPMAXEXTERNALparticle");)
    }


    //if(input.containsData("EOMCC/FROZENOCCUPIED")){
    //  std::string frozen_occ_str;
    //  OPTOPT(frozen_occ_str = input.getData<std::string>("EOMCC/FROZENOCCUPIED"););
    //  eomSettings.frozen_occupied = parseOrbitalSelectionInput(frozen_occ_str);
    //}

    //if(input.containsData("EOMCC/FROZENVIRTUAL")){
    //  std::string frozen_vir_str;
    //  OPTOPT(frozen_vir_str = input.getData<std::string>("EOMCC/FROZENVIRTUAL"););
    //  eomSettings.frozen_virtual = parseOrbitalSelectionInput(frozen_vir_str);
    //}

    if(input.containsData("EOMCC/NROOTS")){

      std::string nRoots;
      nRoots = input.getData<std::string>("EOMCC/NROOTS");
      if ( not nRoots.empty() ) {
        eomSettings.nroots = HandleNRootsInput(nRoots, eomSettings.davidson_Eref);
        if (eomSettings.davidson_Eref.size() == 1 and
            eomSettings.davidson_Eref[0].first == -std::numeric_limits<double>::infinity()) {
          eomSettings.davidson_Eref.clear();
        }
        size_t lowR = eomSettings.nroots;
        for (auto & pair: eomSettings.davidson_Eref) lowR -= pair.second;
        eomSettings.davidson_nLowRoots = lowR;
      }
    }

    if(input.containsData("EOMCC/DAVIDSONWHENSC")){
      OPTOPT(eomSettings.davidson_whenSc = input.getData<int>("EOMCC/DAVIDSONWHENSC");)
    }

    if(input.containsData("EOMCC/DAVIDSONMAXMACROITER")){
      OPTOPT(eomSettings.davidson_max_macro_iter = input.getData<int>("EOMCC/DAVIDSONMAXMACROITER");)
    }

    if(input.containsData("EOMCC/DAVIDSONMAXMICROITER")){
      OPTOPT(eomSettings.davidson_max_micro_iter = input.getData<int>("EOMCC/DAVIDSONMAXMICROITER");)
    }

    if(input.containsData("EOMCC/DAVIDSONRESIDUALCONV")){
      OPTOPT(eomSettings.davidson_residual_conv = input.getData<double>("EOMCC/DAVIDSONRESIDUALCONV");)
    }

    if(input.containsData("EOMCC/DAVIDSONEVECCONV")){
      OPTOPT(eomSettings.davidson_eigen_vector_conv = input.getData<double>("EOMCC/DAVIDSONEVECCONV");)
    }

    if(input.containsData("EOMCC/DAVIDSONEVALCONV")){
      OPTOPT(eomSettings.davidson_eigen_value_conv = input.getData<double>("EOMCC/DAVIDSONEVALCONV");)
    }

    if(input.containsData("EOMCC/DAVIDSONRCHECKESIDUAL")){
      OPTOPT(eomSettings.davidson_check_residual = input.getData<bool>("EOMCC/DAVIDSONRCHECKESIDUAL");)
    }

    if(input.containsData("EOMCC/DAVIDSONCHECKEVEC")){
      OPTOPT(eomSettings.davidson_check_eigen_vector = input.getData<bool>("EOMCC/DAVIDSONCHECKEVEC");)
    }

    if(input.containsData("EOMCC/DAVIDSONCHECKEVAL")){
      OPTOPT(eomSettings.davidson_check_eigen_value = input.getData<bool>("EOMCC/DAVIDSONCHECKEVAL");)
    }

    if(input.containsData("EOMCC/DAVIDSONCONVONGRAMSCHMIDT")){
      OPTOPT(eomSettings.davidson_conv_on_GramSchmidt = input.getData<bool>("EOMCC/DAVIDSONCONVONGRAMSCHMIDT");)
    }

    if(input.containsData("EOMCC/DAVIDSONSUBSPACEMULTIPLIER")){
      OPTOPT(eomSettings.davidson_subspace_multiplier = input.getData<int>("EOMCC/DAVIDSONSUBSPACEMULTIPLIER");)
    }

    if(input.containsData("EOMCC/DAVIDSONGUESSMULTIPLIER")){
      OPTOPT(eomSettings.davidson_guess_multiplier = input.getData<int>("EOMCC/DAVIDSONGUESSMULTIPLIER");)
    }

    if(input.containsData("EOMCC/DAVIDSONENERGYSPECIFICABS")){
      OPTOPT(eomSettings.davidson_ErefAbs = input.getData<bool>("EOMCC/DAVIDSONENERGYSPECIFICABS");)
    }

    if(input.containsData("EOMCC/DAVIDSONPRECONDSMALL")){
      OPTOPT(eomSettings.davidson_preCond_small = input.getData<double>("EOMCC/DAVIDSONPRECONDSMALL");)
    }

//    if(input.containsData("EOMCC/DAVIDSONSORTBYDISTANCE")){
//      OPTOPT(eomSettings.davidson_sort_by_distance = input.getData<bool>("EOMCC/DAVIDSONSORTBYDISTANCE");)
//    }

    if(input.containsData("EOMCC/DAVIDSONBIORTHO")){
      OPTOPT(eomSettings.davidson_biortho = input.getData<bool>("EOMCC/DAVIDSONBIORTHO");)
    }

    if(input.containsData("EOMCC/GRAMSCHMIDTREPEAT")){
      OPTOPT(eomSettings.GramSchmidt_NRe = input.getData<int>("EOMCC/GRAMSCHMIDTREPEAT");)
    }

    if(input.containsData("EOMCC/GRAMSCHMIDTEPS")){
      OPTOPT(eomSettings.GramSchmidt_eps = input.getData<double>("EOMCC/GRAMSCHMIDTEPS");)
    }

    if(input.containsData("EOMCC/OSCILLATORSTRENGTH")){
      OPTOPT(eomSettings.oscillator_strength = input.getData<bool>("EOMCC/OSCILLATORSTRENGTH");)
    }

    if(input.containsData("EOMCC/ALLEXCITEDSTATEDIPOLE")){
      OPTOPT(eomSettings.all_excited_dipole = input.getData<bool>("EOMCC/ALLEXCITEDSTATEDIPOLE");)
    }

    if(input.containsData("EOMCC/SAVEHAMILTONIAN")){
      OPTOPT(eomSettings.save_hamiltonian = input.getData<bool>("EOMCC/SAVEHAMILTONIAN");)
    }

    if(input.containsData("EOMCC/PRINTLARGEAMPLITUDE")){
      if (eomSettings.eom_type != EOM_TYPE::EE) {
        CErr("PrintLargeAmplitude only implemented with EOMCCSD and CVS-EOMCCSD");
      }
      OPTOPT(eomSettings.print_large_amplitude = input.getData<bool>("EOMCC/PRINTLARGEAMPLITUDE");)
    }

    if(input.containsData("EOMCC/CCSGUESS")){
      OPTOPT(eomSettings.ccs_guess = input.getData<bool>("EOMCC/CCSGUESS");)
    }

    if(input.containsData("EOMCC/SAVER")){
      OPTOPT(eomSettings.save_r = input.getData<bool>("EOMCC/SAVER");)
    }

    if(input.containsData("EOMCC/SAVEL")){
      OPTOPT(eomSettings.save_l = input.getData<bool>("EOMCC/SAVEL");)
    }

    if(input.containsData("EOMCC/RESTARTR")){
      OPTOPT(eomSettings.restart_r = input.getData<bool>("EOMCC/RESTARTR");)
      if (eomSettings.restart_r)
        eomSettings.davidson_guess_multiplier = 1;
    }

    if(input.containsData("EOMCC/RESTARTL")){
      OPTOPT(eomSettings.restart_l = input.getData<bool>("EOMCC/RESTARTL");)
      if (eomSettings.restart_l) {
        if (not eomSettings.restart_r) {
          eomSettings.skip_r = true;
        }
        eomSettings.davidson_guess_multiplier = 1;
      }
    }
    if(input.containsData("EOMCC/SKIPR")){
      OPTOPT(eomSettings.skip_r = input.getData<bool>("EOMCC/SKIPR");)
    }
    if(input.containsData("EOMCC/SINGLETONLY")){
      OPTOPT(eomSettings.singlet_only = input.getData<bool>("EOMCC/SINGLETONLY");)
    }
    if(input.containsData("EOMCC/ENERGYIMAGTOL")){
      OPTOPT(eomSettings.energy_imag_tol = input.getData<double>("EOMCC/ENERGYIMAGTOL");)
    }

    //if (eomSettings.doCVS() && eomSettings.eom_type != EOM_TYPE::DIP) {
    //  if (eomSettings.hbar_type != EOM_HBAR_TYPE::EXPLICIT
    //      or eomSettings.diag_method != EOM_DIAG_METHOD::FULL) {
    //    eomSettings.hbar_type = EOM_HBAR_TYPE::EXPLICIT;
    //    eomSettings.diag_method = EOM_DIAG_METHOD::FULL;
    //    std::cout << "CVS-EOMCC only implemented with full diagonalization of explicit Hbar."
    //                 "Settings changed accordingly." << std::endl;
    //  }
    //}

    eomSettings.eom_implementation = eomSettings.find_eom_implementation();
    return eomSettings;
  }
#endif
};

