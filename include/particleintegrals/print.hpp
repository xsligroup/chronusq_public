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
#include <integrals.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <cxxapi/output.hpp>


namespace ChronusQ {

  std::ostream& operator<<(std::ostream &out,
                           const ParticleIntegrals &ints) {
    ints.output(out);
    return out;
  }

  std::ostream& operator<<(std::ostream &out, const IntegralsBase &aoints) {

/*
    out << "\nIntegral Engine Settings:\n" << BannerTop << "\n\n" ;
    out << std::left;


    XXX: Hard code this for now
    out << "  " << std::setw(28) << "ERI Engine:" << "Libint2" << std::endl;
    out << "  " << std::setw(28) << "One-Body Engine:" 
        << "Libint2 + In-House" << std::endl;


    out << std::endl;
    out << "  " << std::setw(28) << "Core Hamiltonian:";
    if(aoints.coreType == NON_RELATIVISTIC) 
      out << "Non-Relativistic";
    else                      
      out << "Relativistic (X2C)";
    out << std::endl;
    
    if(aoints.coreType == RELATIVISTIC_X2C_1E)
      out << "    * Using Finite Width Gaussian Nuclei\n\n";
*/


    out << std::endl;
    out << "  Property Integrals:\n";
    out << "    * Will Compute Length Gauge Electric Multipoles up to Octupole"
        << std::endl;
    out << "    * Will Compute Velocity Gauge Electric Multipoles up to Octupole"
        << std::endl;
    out << "    * Will Compute Magnetic Multipoles up to Quadrupole"
        << std::endl;
    out << std::endl;


    try{
      out << *dynamic_cast<const Integrals<double>&>(aoints).TPI;
    } catch(const std::bad_cast& e) {
      try{
        out << *dynamic_cast<const Integrals<dcomplex>&>(aoints).TPI;
      } catch(const std::bad_cast& e) {
        CErr("IntsT type of AOIntegrals is neither double nor dcomplex");
      }
    }
    
    out << "    * AO to MO Transformation Algorithm (if used): "; 
    
    if (aoints.TPITransAlg == TPI_TRANSFORMATION_ALG::INCORE_N6)
       out << "INCORE N6";
    else if (aoints.TPITransAlg == TPI_TRANSFORMATION_ALG::DIRECT_N6)
       out << "DIRECT N6";
    else if (aoints.TPITransAlg == TPI_TRANSFORMATION_ALG::INCORE_N5)
       out << "INCORE N5";
    else if (aoints.TPITransAlg == TPI_TRANSFORMATION_ALG::DIRECT_N5)
       out << "DIRECT N5";
    else
       CErr("Unrecognized TPI_TRANSFORMATION_ALG in aoints");

    out << std::endl;

    out << std::endl << BannerEnd << std::endl;

    return out; // return std::ostream reference

  }; // operator<<(AOIntegrals)

}; // namespace ChronusQ

