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

namespace ChronusQ {

  enum ResponseOperator {

    LenElectricDipole,
    LenElectricQuadrupole,
    LenElectricOctupole,
    VelElectricDipole,
    VelElectricQuadrupole,
    VelElectricOctupole,
    MagneticDipole,
    MagneticQuadrupole,
    Brillouin

  };

  static std::map<ResponseOperator,size_t> OperatorSize = {
    { LenElectricDipole,      3  },
    { LenElectricQuadrupole,  6  },
    { LenElectricOctupole,    10 },
    { VelElectricDipole,      3  },
    { VelElectricQuadrupole,  6  },
    { VelElectricOctupole,    10 },
    { MagneticDipole,         3  },
    { MagneticQuadrupole,     6  },
    { Brillouin,              1  }
  };

  static std::vector<ResponseOperator> AllOps = {
    LenElectricDipole,
    LenElectricQuadrupole,
    LenElectricOctupole,
    VelElectricDipole,
    VelElectricQuadrupole,
    VelElectricOctupole,
    MagneticDipole,
    MagneticQuadrupole
  };

  static std::vector<ResponseOperator> HerOps = {
    LenElectricDipole,
    LenElectricQuadrupole,
    LenElectricOctupole
  };

  static std::vector<ResponseOperator> AntiHerOps = {
    VelElectricDipole,
    VelElectricQuadrupole,
    VelElectricOctupole,
    MagneticDipole,
    MagneticQuadrupole
  };


  static inline bool isHerOp(ResponseOperator op){

    return (std::find(HerOps.begin(),HerOps.end(),op) != HerOps.end() );

  }


  static inline bool isAntiHerOp(ResponseOperator op) {

    return not isHerOp(op);

  };





  enum ResponseType {
    RESIDUE,
    FDR
  };


  enum SINGLESLATER_POLAR_COPT {

    FULL,
    M,K,MK,KM

  };


  enum ParticleParticleProp_SpinSep {
    PP_AA,
    PP_AB,
    PP_BB
  };

  enum ParticleParticleTDA {
    PP_A,
    PP_C
  };


} // namespace ChronusQ

