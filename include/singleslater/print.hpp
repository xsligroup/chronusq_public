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
#include <fockbuilder/rofock.hpp>


namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printFock(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;
    if(this->nC == 4){ NB *= 2;}

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
    if(this->nC == 4){ NB *= 2;}

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
    if(this->nC == 4){ NB *= 2;}

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
    if(this->nC == 4){ NB *= 2;}

    prettyPrintSmart(out,"J (AO) Scalar",coulombMatrix->pointer(),NB,NB,NB);
    if (twoeH->hasZ())
      prettyPrintSmart(out,"J (AO) MZ",twoeH->Z().pointer(),NB,NB,NB);
    if (twoeH->hasXY()) {
      prettyPrintSmart(out,"J (AO) MY",twoeH->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"J (AO) MX",twoeH->X().pointer(),NB,NB,NB);
    }

  }; // SingleSlater<T>::printJ


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printK(std::ostream &out) {

    size_t NB = this->basisSet().nBasis;
    if(this->nC == 4){ NB *= 2;}

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

    out << std::endl << "Mulliken Charge Analysis:" << std::endl << bannerTop << std::endl;

    out << std::setw(15) << std::left << "  Atom";
    out << std::setw(10) << std::right << "RHO";
    if( this->onePDM->hasZ() ) out << std::setw(15) << std::right << "MZ";
    if( this->onePDM->hasXY() ){
      out << std::setw(15) << std::right << "MY";
      out << std::setw(15) << std::right << "MX";
    }
    out << std::endl;

    out << std::right << bannerMid << std::endl;

    Molecule mol = this->molecule();

    if (this->particle.charge > 0) {
      if(!this->ownedAtomIndices.empty())
        mol = mol.retainAtoms(this->ownedAtomIndices);   // atoms carrying this subsystem's basis
      else
        mol = mol.retainQNuc();                          // legacy fallback: all quantum nuclei
    }

    // For angular momentum printing
    constexpr size_t maxLPrint = 6 + 1;
    size_t maxL = this->basisSet().maxL;
    std::array< std::string, maxLPrint > angLabel =
      { "S", "P", "D", "F", "G", "H", "I" };
    if( maxL + 1 > maxLPrint ) std::cout << "***** WARNING: Printing with L>I basis functions is NYI" << std::endl;

    for(auto iAtm = 0; iAtm < mullikenCharges.size(); iAtm++) {

      // Get symbol
      std::map<std::string,Atom>::const_iterator it = 
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){ 
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      out << "  " << std::setw(7) << std::left <<
        (it == atomicReference.end() ? "X" : it->first);
      out << std::setw(5) << std::right << "TOT";

      out << std::fixed << std::setprecision(5) << std::right;

      out << std::setw(13) << std::right << mullikenCharges[iAtm].get(DENSITY_TYPE::SCALAR);
      if( this->onePDM->hasZ() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].get(DENSITY_TYPE::MZ);
      if( this->onePDM->hasXY() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].get(DENSITY_TYPE::MY);
      if( this->onePDM->hasXY() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].get(DENSITY_TYPE::MX);
      out << std::endl;

      for( size_t iAng=0; iAng<maxL+1; iAng++ ){
        if( iAng >= maxLPrint ) out << std::setw(14) << std::right << "X";
        else                   out << std::setw(14) << std::right << angLabel[iAng];
        out << std::setw(13) << std::right << mullikenCharges[iAtm].getL(DENSITY_TYPE::SCALAR,iAng);
        if( this->onePDM->hasZ() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].getL(DENSITY_TYPE::MZ,iAng);
        if( this->onePDM->hasXY() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].getL(DENSITY_TYPE::MY,iAng);
        if( this->onePDM->hasXY() ) out << std::setw(15) << std::right << mullikenCharges[iAtm].getL(DENSITY_TYPE::MX,iAng);
        out << std::endl;
      }

    }

    out << std::endl << bannerEnd << std::endl;

    out << std::endl << "Lowdin Charge Analysis:" << std::endl << bannerTop << std::endl;

    out << std::setw(15) << std::left << "  Atom";
    out << std::setw(5) << std::right << "RHO";
    out << std::endl;

    out << std::right << bannerMid << std::endl;

    for(auto iAtm = 0; iAtm < lowdinCharges.size(); iAtm++) {

      // Get symbol
      std::map<std::string,Atom>::const_iterator it =
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      out << "  " << std::setw(12) << std::left <<
        (it == atomicReference.end() ? "X" : it->first);

      out << std::fixed << std::setprecision(5) << std::right;

      out << std::setw(8) << std::right << lowdinCharges[iAtm].get(DENSITY_TYPE::SCALAR);
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
    const bool iRO = std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(fockBuilder) != nullptr; //ROHF calculation
    const bool iU = this->nC == 1 and not this->iCS and not iRO; // one-component unrestricted calculation such as UHF, UKS, etc.

    auto printSection = [&out](const std::string& label, size_t count, size_t leadingNewlines = 0) {
      for(size_t i = 0; i < leadingNewlines; i++) out << std::endl;
      out << label << ": (" << count << ")" << std::endl;
    };

    auto printOrbitalRange = [&out](size_t start, size_t end) {
      out << "(" << std::setw(3) << start;
      if (start == end)
        out << "    )";
      else
        out << " -" << std::setw(3) << end << ")";
    };

    auto printEnergyLines = [&out, &printOrbitalRange](size_t begin, size_t end, double* eps) {
      for(size_t i = begin; i < end; i++) {
        const size_t local = i - begin;
        if(local % 5 == 0)
          printOrbitalRange(i + 1, std::min(i + 5, end));

        out << std::setw(13) << eps[i];

        if((local + 1) % 5 == 0 or i + 1 == end)
          out << std::endl;
      }
    };

    // List MO eigenenergies

    out << std::scientific << std::setprecision(4);
    out << "Orbital Eigenenergies " << (iU ? "(Alpha) " : "" )
        << "/ Eh" << std::endl << bannerTop << std::endl;

    if(this->nC == 1 and not iU) { // closed shell case
      printSection("Doubly Occupied", this->nOB);
      printEnergyLines(0, this->nOB, this->eps1);

      if(this->nOA > this->nOB) { // ROHF
        printSection("Singly Occupied", this->nOA - this->nOB, 1);
        printEnergyLines(this->nOB, this->nOA, this->eps1);
      }

      printSection("Virtual", this->nVA, 1);
      printEnergyLines(this->nOA, this->nAlphaOrbital(), this->eps1);

    } else {
      const size_t nTotal = this->nC * this->nAlphaOrbital();
      const size_t nNegative = this->nC == 4 ? nTotal / 2 : 0;
      const size_t NO = iU? this->nOA : this->nO;
      const size_t NV = iU? this->nVA : this->nV;

      printSection("Occupied", NO);
      printEnergyLines(nNegative,  nNegative + NO, this->eps1);

      printSection("Virtual", NV, 1);
      printEnergyLines(nNegative + NO, nTotal, this->eps1);

      if (this->nC == 4) {
        printSection("Negative Energy States", nNegative, 1);
        printEnergyLines(0, nNegative, this->eps1);
      }
    }

    out << bannerEnd << std::endl;

    if (iU) {
      out << std::endl << std::endl << "Orbital Eigenenergies (Beta) / Eh"<< std::endl << bannerTop << std::endl;

      printSection("Occupied", this->nOB);
      printEnergyLines(0, this->nOB, this->eps2);

      printSection("Virtual", this->nVB, 1);
      printEnergyLines(this->nOB, this->nBetaOrbital(), this->eps2);

      out << bannerEnd << std::endl;
    }

  } // SingleSlater::printEPS


  /**
   * \brief Print out general MO information
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printMOInfo(std::ostream &out, size_t printMOLevel) {

    out << std::endl << "SCF Results:" << std::endl << BannerTop << std::endl;

    // print MO eigenvalues
    this->printEPS(out);

    if (not printMOLevel)
      printMOLevel = scfControls.printMOCoeffs;
    // print MO coefficients
    WaveFunction<MatsT,IntsT>::printMOInfo(out, printMOLevel);

  }; // SingleSlater<T>::printMOInfo


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT, IntsT>::printFockTimings(std::ostream &out) {

    out << "    Fock Timings:" << std::endl;
    out << "      Wall time G[D] = " << std::setw(8)
        << std::setprecision(5)  << std::scientific
        << GDDur << " s" << std::endl;


  }; // SingleSlater<T>::printFockTimings

}; // namespace ChronusQ

