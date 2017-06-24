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

#include <memmanager.hpp>

namespace ChronusQ {

  template <typename T>
  struct ResidueResponseResults {

    // Solution Storage
    double *W  = nullptr;
    T      *VR = nullptr;
    T      *VL = nullptr;


    // Property Storage
    T * tLenElecDipole_ge      = nullptr;
    T * tLenElecQuadrupole_ge  = nullptr;
    T * tLenElecOctupole_ge    = nullptr; 
    T * tVelElecDipole_ge      = nullptr;
    T * tVelElecQuadrupole_ge  = nullptr;
    T * tVelElecOctupole_ge    = nullptr; 

    T * tMagDipole_ge      = nullptr;
    T * tMagQuadrupole_ge  = nullptr;
    T * tMagOctupole_ge    = nullptr; 

    inline void dealloc() {

      if(W) CQMemManager::get().free(W);
      if(VR) CQMemManager::get().free(VR);
      if(VL) CQMemManager::get().free(VL);

      if(tLenElecDipole_ge    ) CQMemManager::get().free(tLenElecDipole_ge    );
      if(tLenElecQuadrupole_ge) CQMemManager::get().free(tLenElecQuadrupole_ge);
      if(tLenElecOctupole_ge  ) CQMemManager::get().free(tLenElecOctupole_ge  ) ;
      if(tVelElecDipole_ge    ) CQMemManager::get().free(tVelElecDipole_ge    );
      if(tVelElecQuadrupole_ge) CQMemManager::get().free(tVelElecQuadrupole_ge);
      if(tVelElecOctupole_ge  ) CQMemManager::get().free(tVelElecOctupole_ge  ) ;


      if(tMagDipole_ge    ) CQMemManager::get().free(tMagDipole_ge    );
      if(tMagQuadrupole_ge) CQMemManager::get().free(tMagQuadrupole_ge);
    }

  };

  template <typename T, typename U>
  struct FDResponseResults {

    std::vector<U> shifts;


    T* RHS = nullptr;
    U* SOL = nullptr;


    // Polarizabilities
    U* ed_ed_Polar = nullptr;
    U* eq_ed_Polar = nullptr;
    U* eo_ed_Polar = nullptr;
    U* eq_eq_Polar = nullptr;
    U* eo_eq_Polar = nullptr;
    U* eo_eo_Polar = nullptr;

    U* md_ed_Polar = nullptr;
    U* md_md_Polar = nullptr;

    inline void dealloc() {

      if(RHS) CQMemManager::get().free(RHS);
      if(SOL) CQMemManager::get().free(SOL);

      if(ed_ed_Polar) CQMemManager::get().free(ed_ed_Polar);
      if(eq_ed_Polar) CQMemManager::get().free(eq_ed_Polar);
      if(eo_ed_Polar) CQMemManager::get().free(eo_ed_Polar);
      if(eq_eq_Polar) CQMemManager::get().free(eq_eq_Polar);
      if(eo_eq_Polar) CQMemManager::get().free(eo_eq_Polar);
      if(eo_eo_Polar) CQMemManager::get().free(eo_eo_Polar);

      if(md_ed_Polar) CQMemManager::get().free(md_ed_Polar);
      if(md_md_Polar) CQMemManager::get().free(md_md_Polar);

    };

  };

  struct FDObservables {

    double * edStrength    = nullptr;
    double * opaCross_eda  = nullptr;
    double * ecd_len_RM    = nullptr;
    double * ecd_vel_PMQ   = nullptr;

    inline void dealloc() {

      if(edStrength)   CQMemManager::get().free(edStrength);
      if(opaCross_eda) CQMemManager::get().free(opaCross_eda);
      if(ecd_len_RM)   CQMemManager::get().free(ecd_len_RM);
      if(ecd_vel_PMQ)  CQMemManager::get().free(ecd_vel_PMQ);

    }
  };

  struct ResObservables {

    double * oscStrength      = nullptr;
    double * rotatory_len_RM  = nullptr;
    double * rotatory_vel_PMQ = nullptr;

    inline void dealloc() {

      if(oscStrength)      CQMemManager::get().free(oscStrength);
      if(rotatory_len_RM)  CQMemManager::get().free(rotatory_len_RM);
      if(rotatory_vel_PMQ) CQMemManager::get().free(rotatory_vel_PMQ);

    }

  };

}

