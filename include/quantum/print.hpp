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

    size_t NB = onePDM->dimension();

    prettyPrintSmart(out,"1PDM (AO) Scalar",onePDM->S().pointer(),NB,NB,NB);

    if( onePDM->hasZ() )
      prettyPrintSmart(out,"1PDM (AO) MZ",onePDM->Z().pointer(),NB,NB,NB);

    if( onePDM->hasXY() ) {
      prettyPrintSmart(out,"1PDM (AO) MY",onePDM->Y().pointer(),NB,NB,NB);
      prettyPrintSmart(out,"1PDM (AO) MX",onePDM->X().pointer(),NB,NB,NB);
    }

  }; // Quantum<T>::print1PDM

  void QuantumBase::printMultipoles(std::ostream &out) {

    out << std::left << std::setprecision(10);
    double xSmall = 1e-10;
    for(auto i = 0; i < 3; i++) {
      if(std::abs(this->elecDipole[i]) < xSmall)
        this->elecDipole[i] = std::abs(this->elecDipole[i]);

    for(auto j = 0; j < 3; j++) {
      if(std::abs(this->elecQuadrupole[i][j]) < xSmall)
        this->elecQuadrupole[i][j] = std::abs(this->elecQuadrupole[i][j]);

    for(auto k = 0; k < 3; k++)
      if(std::abs(this->elecOctupole[i][j][k]) < xSmall)
        this->elecOctupole[i][j][k] = 
          std::abs(this->elecOctupole[i][j][k]);

    }
    }

    out << "\nMultipole Information:" << std::endl;
    out << bannerTop << std::endl << std::endl;;



    // Dipoles
    out << std::setw(50) << std::left <<"Electric Dipole Moment"
                          << "(Debye)" << std::endl;

    out << std::fixed; 
    out <<  std::setw(5) << std::left << "X=" << std::setw(20) << std::right 
           << this->elecDipole[0] / EBohrPerDebye;

    out <<  std::setw(5) << std::left << " Y=" << std::setw(20) << std::right 
           << this->elecDipole[1] / EBohrPerDebye;

    out <<  std::setw(5) << std::left << " Z=" << std::setw(20) << std::right 
           << this->elecDipole[2] / EBohrPerDebye;


    // Quadrupoles
    out << std::endl << std::endl;
    out << std::setw(50) << std::left << "Electric Quadrupole Moment" 
                         <<  "(Debye-\u212B)" << std::endl;

    out << std::left << std::setw(5) <<"XX=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[0][0] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XY=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[0][1] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XZ=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[0][2] * AngPerBohr/EBohrPerDebye 
                     << std::endl;

    out << std::left << std::setw(5) <<"YX=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[1][0] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YY=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[1][1] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YZ=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[1][2] * AngPerBohr/EBohrPerDebye 
                     << std::endl;

    out << std::left << std::setw(5) <<"ZX=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[2][0] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZY=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[2][1] * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZZ=" << std::right << std::setw(20) 
                     << this->elecQuadrupole[2][2] * AngPerBohr/EBohrPerDebye 
                     << std::endl;

    out << std::endl << std::endl;


    // Octupoles
    out << std::endl << std::endl;
    out << std::setw(50) << std::left 
                       << "Electric Octupole Moment" 
                       << "(Debye-\u212B\u00B2)" << std::endl;

    out << std::left << std::setw(5) <<"XXX=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][0][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XXY=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][0][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XXZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][0][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"XYX=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][1][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XYY=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][1][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XYZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][1][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"XZX=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][2][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XZY=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][2][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" XZZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[0][2][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"YXX=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][0][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YXY=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][0][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YXZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][0][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"YYX=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][1][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YYY=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][1][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YYZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][1][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"YZX=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][2][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YZY=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][2][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" YZZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[1][2][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"ZXX=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][0][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZXY=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][0][1] *
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZXZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][0][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"ZYX=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][1][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZYY=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][1][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZYZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][1][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;

    out << std::left << std::setw(5) <<"ZZX=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][2][0] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZZY=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][2][1] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye;

    out << std::left << std::setw(5) <<" ZZZ=" << std::right << std::setw(20) 
                       << this->elecOctupole[2][2][2] * 
                         AngPerBohr * AngPerBohr/EBohrPerDebye << std::endl;




    out << bannerEnd << std::endl << std::endl;
    

  }; // QuantumBase::printMultipoles


  void QuantumBase::printSpin(std::ostream &out) {
    out << "Spin Information:" << std::endl;
    out << bannerTop << std::endl << std::endl;
    out << std::setprecision(5) << std::fixed;
    out << "  <Sx> = " << std::setw(10) << std::right 
                       << this->SExpect[0] << std::endl;
    out << "  <Sy> = " << std::setw(10) << std::right 
                       << this->SExpect[1] << std::endl;
    out << "  <Sz> = " << std::setw(10) << std::right 
                       << this->SExpect[2] << std::endl;
    out << "  <S\u00B2> = " << std::setw(10) << std::right 
                       << this->SSq << std::endl;
    out << std::endl << bannerEnd << std::endl << std::endl;

  }; // QuantumBase::printSpin
};


