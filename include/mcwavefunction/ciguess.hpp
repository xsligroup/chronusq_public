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

#include <mcwavefunction.hpp>
//#include <cqlinalg/blas1.hpp>
//#include <cqlinalg/blas3.hpp>
//#include <cqlinalg/blasutil.hpp>
//#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

namespace ChronusQ {

  /*
   * \brief Compute 1rdm
   *         
   */ 
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::ReadGuessCIVector() {
  // LL: I have assumed it to work with a case which reads CI from multiple binfiles,
  //     such as binfiles from energy specific case with each file covers a subset of CI roots.
  //     but currently reading from multiple scr files or restart files are not implemented in
  //     RunChronusQ.
  //void MCWaveFunction<MatsT,IntsT>::ReadGuessCIVector(std::vector<std::string> binNames) {

    std::vector<std::string> binNames;
    std::cout << "    * Reading in guess CI vectors from file(s): ";


    if ( ref_.scrBinFileName.empty() ) {
      binNames = {this->savFile.fName()};
      std::cout << this->savFile.fName() << std::endl;
    } else {
      binNames = {ref_.scrBinFileName};
      std::cout << ref_.scrBinFileName << std::endl;
    }

    std::cout << "    * Please make sure the space partitions are the same!" << std::endl;

    size_t t_hash = std::is_same<MatsT,double>::value ? 1 : 2;
    size_t d_hash = 1;
    size_t c_hash = 2;
    size_t savHash;
    bool BinExists;

    std::string prefix = "/MCWFN/";

    // Check number of states, data type and number of determinants
    size_t fileNR = 0;
    size_t tempNR = 0;
    for (auto & binName: binNames) {
      SafeFile binFile(binName, BinExists);

      try{
        binFile.readData(prefix + "FIELD_TYPE", &savHash);
      } catch (...) {
        CErr("Cannot find /MCWFN/FIELD_TYPE on rstFile!",std::cout);
      }

      if ( t_hash != savHash ) {
        std::string t_field = (t_hash == d_hash)? "REAL" : "COMPLEX";
        std::string s_field = (savHash == d_hash)? "REAL" : "COMPLEX";

        std::string message = prefix + "FIELD_TYPE on disk (" + s_field +
            ") is incompatible with current FIELD_TYPE (" + t_field + ")";
        CErr(message,std::cout);
      }

      try{
        binFile.readData(prefix + "NSTATES", &tempNR);
      } catch (...) {
        CErr("Cannot find /MCWFN/NSTATES on rstFile!",std::cout);
      }
      fileNR += tempNR;

      auto nDet = binFile.getDims( prefix + "CIVec_1" );
      if ( nDet[0] != this->NDet )
        CErr("CI Vectors dimensions in rstFile is incompatible with current NDet!");
    }
    if (fileNR != this->NStates)
      CErr("Total number of states in rstFiles does not match with current input!");

    // Read in CI vectors
    fileNR = 0;
    for (auto & binName: binNames) {

      SafeFile binFile(binName, BinExists);

      binFile.readData(prefix + "NSTATES", &tempNR);
      for (auto i = 0; i < tempNR; i++) {
        std::cout << "    * Found MCWFN/CIVec_" <<i+1<< std::endl;
        binFile.readData(prefix + "CIVec_"+std::to_string(i+1), CIVecs[fileNR+i]);
      }
      // read in state energies
      binFile.readData(prefix + "STATE_ENERGY", this->StateEnergy.data() + fileNR);
      fileNR += tempNR;
    }

  } // MCWaveFunction::ReadGuessCIVector

}; // namespace ChronusQ


