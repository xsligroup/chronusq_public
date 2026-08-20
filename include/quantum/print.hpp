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

#include <quantum.hpp>
#include <util/matout.hpp>
#include <physcon.hpp>


namespace ChronusQ {

  template <typename MatsT>
  void Quantum<MatsT>::print1PDM(std::ostream &out) {

    size_t NB = onePDM->nRows();

    prettyPrintSmart(out,"1PDM (AO) Scalar",onePDM->S().pointer(),NB,NB,NB);

    if( onePDM->hasZ() )
      prettyPrintSmart(out,"1PDM (AO) MZ",onePDM->Z().pointer(),NB,NB,NB);

    if( onePDM->hasXY() ) {
      prettyPrintSmart(out,"1PDM (AO) MY",onePDM->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"1PDM (AO) MX",onePDM->X().pointer(),NB,NB,NB);
    }

  }; // Quantum<T>::print1PDM

  // QuantumBase::printMultipoles is defined in src/quantum/print.cxx

  void QuantumBase::printSpin(std::ostream &out, bool withBanner) {
    
    if (withBanner)
    {
      out << "Spin Information:" << std::endl;
      out << bannerTop << std::endl << std::endl;
    }

    if (this->nC == 4) {
        out << "  NOTE: Spin expectation values are not computed for 4C." << std::endl;
        return;
    }

    // Print spin expectation values
    out << std::setprecision(12) << std::fixed;
    out << "  <Sx> = " << std::setw(10) << std::right 
                       << this->SExpect[0] << std::endl;
    out << "  <Sy> = " << std::setw(10) << std::right 
                       << this->SExpect[1] << std::endl;
    out << "  <Sz> = " << std::setw(10) << std::right 
                       << this->SExpect[2] << std::endl;
    out << "  <S^2> = " << std::setw(10) << std::right 
                       << this->SSq << std::endl;
    // Print the eigenvalues of the spin operator S^2 
    out << "  spin quantum number = " << std::setw(10) << std::right << this->SQuantNum << "\n"<< std::endl;

    if (this->nC == 4) {
      out << "\n  4C Spin Decomposition (LL/SS):" << std::endl;
      out << "    <Sx>_LL = " << std::setw(12) << std::right << this->SExpectLL[0]
        << "    <Sx>_SS = " << std::setw(12) << std::right << this->SExpectSS[0] << std::endl;
      out << "    <Sy>_LL = " << std::setw(12) << std::right << this->SExpectLL[1]
        << "    <Sy>_SS = " << std::setw(12) << std::right << this->SExpectSS[1] << std::endl;
      out << "    <Sz>_LL = " << std::setw(12) << std::right << this->SExpectLL[2]
        << "    <Sz>_SS = " << std::setw(12) << std::right << this->SExpectSS[2] << std::endl;
      out << "    <S^2>_LL = " << std::setw(12) << std::right << this->SSqLL
        << "    <S^2>_SS = " << std::setw(12) << std::right << this->SSqSS << std::endl;
      out << "    <S^2>_cross = " << std::setw(12) << std::right << this->SSqCross << std::endl;
    }

    if (withBanner)
      out << std::endl << bannerEnd << std::endl << std::endl;

  }; // QuantumBase::printSpin

  
};

