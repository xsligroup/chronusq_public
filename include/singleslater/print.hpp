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

#include <singleslater.hpp>
#include <util/matout.hpp>


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printFock(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;

    prettyPrintSmart(out,"Fock (AO) Scalar",fockMatrix->S().pointer(),NB,NB,NB);

    if( fockMatrix->hasZ() )
      prettyPrintSmart(out,"Fock (AO) MZ",fockMatrix->Z().pointer(),NB,NB,NB);

    if( fockMatrix->hasXY() ) {
      prettyPrintSmart(out,"Fock (AO) MY",fockMatrix->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"Fock (AO) MX",fockMatrix->X().pointer(),NB,NB,NB);
    }


  }; // SingleSlater<T>::printFock

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::print1PDMOrtho(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;

    prettyPrintSmart(out,"1PDM (Ortho) Scalar",onePDMOrtho->S().pointer(),NB,NB,NB);

    if( onePDMOrtho->hasZ() )
      prettyPrintSmart(out,"1PDM (Ortho) MZ",onePDMOrtho->Z().pointer(),NB,NB,NB);

    if( onePDMOrtho->hasXY() ) {
      prettyPrintSmart(out,"1PDM (Ortho) MY",onePDMOrtho->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"1PDM (Ortho) MX",onePDMOrtho->X().pointer(),NB,NB,NB);
    }


  }; // SingleSlater<T>::print1PDMOrtho

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printGD(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;

    prettyPrintSmart(out,"GD (AO) Scalar",twoeH->S().pointer(),NB,NB,NB);
    if (twoeH->hasZ())
      prettyPrintSmart(out,"GD (AO) MZ",twoeH->Z().pointer(),NB,NB,NB);
    if (twoeH->hasXY()) {
      prettyPrintSmart(out,"GD (AO) MY",twoeH->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"GD (AO) MX",twoeH->X().pointer(),NB,NB,NB);
    }

  }; // SingleSlater<T>::printGD


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printJ(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;

    prettyPrintSmart(out,"J (AO) Scalar",coulombMatrix->pointer(),NB,NB,NB);


  }; // SingleSlater<T>::printJ


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printK(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;

    prettyPrintSmart(out,"K (AO) Scalar",exchangeMatrix->S().pointer(),NB,NB,NB);
    if (exchangeMatrix->hasZ())
      prettyPrintSmart(out,"K (AO) MZ",exchangeMatrix->Z().pointer(),NB,NB,NB);
    if (exchangeMatrix->hasXY()) {
      prettyPrintSmart(out,"K (AO) MY",exchangeMatrix->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"K (AO) MX",exchangeMatrix->X().pointer(),NB,NB,NB);
    }

  }; // SingleSlater<T>::printK


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printMiscProperties(std::ostream &out) {

    out << "\nCharge Analysis:\n" << bannerTop << "\n\n";

    out << std::setw(20) << std::left << "  Atom";
    out << std::setw(20) << std::right << "Mulliken Charges";
    out << std::setw(20) << std::right << "Lowdin Charges" << std::endl;

    out << std::right << bannerMid << std::endl;

    Molecule mol = this->molecule();

    // For protonic SS, charge analysis are done for only proton atoms 
    if (this->particle.charge > 0)  mol = mol.retainQNuc();

    for(auto iAtm = 0; iAtm < mol.nAtoms; iAtm++) {

      // Get symbol
      std::map<std::string,Atom>::const_iterator it = 
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){ 
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      out << "  " << std::setw(18) << std::left <<
        (it == atomicReference.end() ? "X" : it->first);

      out << std::setprecision(5) << std::right;

      out << std::setw(20) << mullikenCharges[iAtm];
      out << std::setw(20) << lowdinCharges[iAtm];

      out << std::endl;
    }

    out << std::endl << bannerEnd << std::endl;

  }; // SingleSlater<T>::printMiscProperties


  /**
   * \brief Print out MO eigenvalues
   *
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printEPS(std::ostream &out) {

    // List MO eigenenergies

    out << std::scientific << std::setprecision(4);
    out << "Orbital Eigenenergies " << (this->nC == 1 ? "(Alpha) " : "" )
        << "/ Eh\n" << bannerTop << "\n";

    // Set number of occupied orbitals
    size_t NO = (this->nC == 1 ? this->nOA : this->nO);
    if ( this->nC == 4 ) NO = this->nO + this->nC * this->nAlphaOrbital()/2;

    //for(auto i = 0ul; i < this->nC * this->nAlphaOrbital(); i++) {

    //  if( i == 0 )
    //    out << "Occupied:\n";
    //  else if( i == NO )
    //    out << "\n\nVirtual:\n";

    //  out << std::setw(13) << this->eps1[i];

    //  if( i < NO and (i + 1) % 5 == 0 )  out << "\n";
    //  else if( i >= NO and ((i - NO) + 1) % 5 == 0 ) out << "\n";
    //}

    if( nC != 4 ) {
      for(auto i = 0ul; i < this->nC * this->nAlphaOrbital(); i++) {
  
        if( i == 0 )
          out << "Occupied:\n";
        else if( i == NO )
          out << "\n\nVirtual:\n";
  
        out << std::setw(13) << this->eps1[i];
  
        if( i < NO and (i + 1) % 5 == 0 )  out << "\n";
        else if( i >= NO and ((i - NO) + 1) % 5 == 0 ) out << "\n";
      }
    } else {
      for(auto i = this->nC * this->nAlphaOrbital()/2; i < this->nC * this->nAlphaOrbital(); i++) {
  
        if( i == this->nC * this->nAlphaOrbital()/2 )
          out << "Occupied:\n";
        else if( i == NO )
          out << "\n\nVirtual:\n";
  
        out << std::setw(13) << this->eps1[i];
  
        if( i < NO and (i + 1 - this->nC * this->nAlphaOrbital()/2) % 5 == 0 )  out << "\n";
        else if( i >= NO and ((i - NO) + 1) % 5 == 0 ) out << "\n";
      }

      // Negative Energy States
      for(auto i = 0ul; i < this->nC * this->nAlphaOrbital()/2; i++) {
  
        if( i == 0 )
            out << "\n\nNegative Energy States:\n";

        out << std::setw(13) << this->eps1[i];
  
        if( i < NO and (i + 1) % 5 == 0 )  out << "\n";
        else if( i >= NO and ((i - NO) + 1) % 5 == 0 ) out << "\n";
      }
    }
     
    out << "\n" << bannerEnd << "\n";

    if( this->nC == 1 and not this->iCS ) {
      out << "\n\nOrbital Eigenenergies (Beta) / Eh\n" << bannerTop << "\n";

      for(auto i = 0ul; i < this->nBetaOrbital(); i++) {

        if( i == 0 )
          out << "Occupied:\n";
        if( i == this->nOB )
          out << "\n\nVirtual:\n";

        out << std::setw(13) << this->eps2[i];

        if( i < this->nOB and (i + 1) % 5 == 0 )  out << "\n";
        else if( i >= this->nOB and ((i - this->nOB) + 1) % 5 == 0 ) out << "\n";
      }

      out << "\n" << bannerEnd << "\n";
    }

  } // SingleSlater::printEPS


  /**
   * \brief Print out general MO information
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printMOInfo(std::ostream &out, size_t printMOLevel) {

    out << "\n\n" << "SCF Results:\n" << BannerTop << "\n\n";

    // print MO eigenvalues
    this->printEPS(out);

    if (not printMOLevel)
      printMOLevel = scfControls.printMOCoeffs;
    // print MO coefficients
    WaveFunction<MatsT,IntsT>::printMOInfo(out, printMOLevel);

  }; // SingleSlater<T>::printMOInfo


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT, IntsT>::printFockTimings(std::ostream &out) {

    out << "    Fock Timings:\n";
    out << "      Wall time G[D] = " << std::setw(8)
        << std::setprecision(5)  << std::scientific
        << GDDur << " s\n\n";


  }; // SingleSlater<T>::printFockTimings

}; // namespace ChronusQ

